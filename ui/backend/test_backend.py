"""Self-check: python3 test_backend.py

No network and no Jenkins needed. The important assertion is the last one --
params.yml must match the parameters{} block of make_config_file.groovy,
because Jenkins rejects a whole build if it is sent an undeclared parameter.
"""

import json
import os
import re
import tempfile
from pathlib import Path

os.environ.setdefault("PRISMSEQ_ROOT", "/cmap/obelix/pod/prismSeq")
os.environ.setdefault("SUSHI_UI_DB", str(Path(tempfile.mkdtemp()) / "test.db"))
# Point the screen-type store at a path that does not exist, so the tests read
# the shipped presets.yml. Without this they assert against whatever defaults
# the deployment happens to have saved, and edits made in the UI break them.
os.environ.setdefault("SUSHI_UI_PRESETS", str(Path(tempfile.mkdtemp()) / "presets.yml"))

import catalog  # noqa: E402
import db  # noqa: E402
from app import _stale_config  # noqa: E402

GROOVY = Path(__file__).parents[2] / "scripts" / "make_config_file.groovy"


def test_catalog_loads():
    data = catalog.catalog()
    assert data["groups"] and data["params"] and data["presets"]
    names = [p["name"] for p in data["params"]]
    assert len(names) == len(set(names)), "duplicate param names"


def test_presets_resolve():
    values = catalog.defaults_for("EPS_SEQ", "EPS009")
    assert values["COUNT_THRESHOLD"] == "40", "EPS preset must override the count threshold"
    assert values["nc_variability_threshold"] == "1.5"
    assert values["BUILD_NAME"] == "EPS009"
    assert values["SCREEN"] == "EPS009"
    assert values["BUILD_DIR"].endswith("/EPS_SEQ/EPS009")

    # A preset must not leak into the next one.
    assert catalog.defaults_for("MTS_SEQ", "MTS033")["COUNT_THRESHOLD"] == "100"
    assert catalog.defaults_for("CPS_SEQ", "CPS017")["SYNERGY"] is True
    assert catalog.defaults_for("MTS_SEQ", "MTS033")["SYNERGY"] is False


def test_coerce_drops_unknown_and_normalizes_bools():
    out = catalog.coerce({"COUNT_THRESHOLD": 100, "DRC": "true", "COLLAPSE": False,
                          "TOTALLY_MADE_UP": "x"})
    assert out == {"COUNT_THRESHOLD": "100", "DRC": True, "COLLAPSE": False}, out


def test_validate_catches_module_dependencies():
    values = catalog.defaults_for("MTS_SEQ", "MTS033")
    assert catalog.validate(values) == [], catalog.validate(values)

    assert any("required" in p for p in catalog.validate({**values, "BUILD_NAME": ""}))
    assert any("number" in p for p in catalog.validate({**values, "COUNT_THRESHOLD": "abc"}))
    assert any("post-LFC" in p for p in
               catalog.validate({**values, "GENERATE_QC_TABLES_2": True, "COMPUTE_LFC": False}))
    assert any("Dose response" in p for p in catalog.validate({**values, "AUC_BIOMARKER": True}))
    assert any("Pin to commit" in p for p in catalog.validate({**values, "USE_LATEST": False}))
    assert any("one of" in p for p in catalog.validate({**values, "SEQ_TYPE": "Sanger"}))


def test_module_names_cover_the_toggles_config_json_loses():
    modules = set(catalog.module_names())
    # These decide which scripts run and appear nowhere in config.json.
    for name in ("COLLATE_FASTQ_READS", "CBNORMALIZE", "SYNERGY", "GENERATE_WELL_METRICS",
                 "QC_IMAGES", "UPLOAD_TO_BQ", "CREATE_CELLDB_METADATA"):
        assert name in modules, name


def test_stale_config_detects_silently_ignored_edits():
    with tempfile.TemporaryDirectory() as tmp:
        build_dir = Path(tmp)
        (build_dir / "config.json").write_text(json.dumps({
            "COUNT_THRESHOLD": "40", "CTL_TYPES": "ctl_vehicle"}))
        (build_dir / "qc_params.json").write_text(json.dumps({
            "well_reads_threshold": "40", "cb_mae_threshold": "1"}))

        stale = _stale_config(build_dir, {"COUNT_THRESHOLD": "100", "CTL_TYPES": "ctl_vehicle",
                                          "cb_mae_threshold": "1"})
        conflicts = {c["name"] for c in stale["conflicts"]}
        assert conflicts == {"COUNT_THRESHOLD"}, conflicts

        # well_reads_threshold is recomputed from COUNT_THRESHOLD every run, so a
        # difference there is not a conflict.
        stale = _stale_config(build_dir, {"well_reads_threshold": "100"})
        assert stale is None, stale

        assert _stale_config(build_dir, {"CTL_TYPES": "ctl_vehicle"}) is None


def test_run_records_survive_a_round_trip():
    db.init()
    values = catalog.defaults_for("CPS_SEQ", "CPS017")
    run_id = db.create(launched_by="tester", preset="CPS_SEQ", build_name="CPS017",
                       screen="CPS017", screen_type="CPS_SEQ", build_dir="/tmp/CPS017",
                       params=catalog.coerce(values),
                       modules={n: values[n] for n in catalog.module_names()},
                       git_branch="main")
    assert run_id in db.unfinished_ids()

    run = db.get(run_id)
    assert run["status"] == "QUEUED"
    assert run["modules"]["SYNERGY"] is True
    assert run["params"]["COUNT_THRESHOLD"] == "100"

    db.update(run_id, status="SUCCESS", build_number=42)
    assert db.unfinished_ids() == [i for i in db.unfinished_ids() if i != run_id]
    runs, total = db.list_runs(build_name="CPS")
    assert total == 1 and runs[0]["build_number"] == 42


def test_launch_and_status_round_trip_against_a_stub_jenkins():
    """Covers the money path: preflight -> launch -> queued -> running -> done."""
    from fastapi.testclient import TestClient

    import app as app_module
    import jenkins

    import session

    state = {"number": None, "building": True, "result": None, "sent": None, "auth": None}
    jenkins.trigger = lambda params, auth=None: (
        state.__setitem__("sent", params), state.__setitem__("auth", auth),
        "http://j/queue/item/9")[-1]
    jenkins.queue_build_number = lambda url: (state["number"], None)
    jenkins.build = lambda n: {"building": state["building"], "result": state["result"],
                               "timestamp": 1_700_000_000_000, "duration": 0 if state["building"] else 90_000}
    jenkins.console = lambda n, start: {"text": "hello", "start": start, "next": start + 5, "more": False}
    jenkins.build_url = lambda n: f"http://j/job/sushi/{n}/"

    client = TestClient(app_module.app)

    with tempfile.TemporaryDirectory() as tmp:
        build_dir = Path(tmp) / "MTS099"
        build_dir.mkdir()
        (build_dir / "raw_counts_uncollapsed.csv.gz").write_bytes(b"")
        (build_dir / "sample_meta.csv").write_text("x\n")

        values = catalog.defaults_for("MTS_SEQ", "MTS099")
        values["BUILD_DIR"] = str(build_dir)
        body = {"preset": "MTS_SEQ", "values": values}

        # Preflight is read-only and needs no identity.
        pre = client.post("/api/preflight", json=body).json()
        assert pre["problems"] == [], pre["problems"]
        assert pre["stale_config"] is None

        # Launching does: Jenkins will not queue a build for anonymous.
        assert client.post("/api/runs", json=body).status_code == 401
        assert client.get("/api/session").json()["identity"] is None

        session.verify = lambda user, token: {"id": user, "full_name": "Test Er"}
        signed_in = client.post("/api/session", json={"user": "tester", "token": "tok"})
        assert signed_in.status_code == 200, signed_in.text
        assert signed_in.json()["full_name"] == "Test Er"
        assert client.get("/api/session").json()["identity"]["id"] == "tester"

        # A missing input file must block the launch, not fail an hour in.
        (build_dir / "sample_meta.csv").unlink()
        blocked = client.post("/api/preflight", json=body).json()
        assert any("sample_meta.csv" in p for p in blocked["problems"]), blocked
        assert client.post("/api/runs", json=body).status_code == 422
        (build_dir / "sample_meta.csv").write_text("x\n")

        run_id = client.post("/api/runs", json=body).json()["id"]
        assert state["sent"]["SCREEN_TYPE"] == "MTS_SEQ"
        assert state["sent"]["COUNT_THRESHOLD"] == "100"
        # The build is queued as the signed-in user, and the run log records
        # the name Jenkins gave rather than anything the client claimed.
        assert state["auth"] == ("tester", "tok")
        assert client.get(f"/api/runs/{run_id}").json()["launched_by"] == "Test Er"

        # Still in the Jenkins queue: no build number yet.
        run = client.get(f"/api/runs/{run_id}").json()
        assert run["status"] == "QUEUED" and run["build_number"] is None

        state["number"] = 7
        run = client.get(f"/api/runs/{run_id}").json()
        assert run["status"] == "RUNNING" and run["build_number"] == 7
        assert run["jenkins_url"].endswith("/7/")
        assert "COLLATE_FASTQ_READS" in run["enabled_modules"]

        state["building"], state["result"] = False, "SUCCESS"
        run = client.get(f"/api/runs/{run_id}").json()
        assert run["status"] == "SUCCESS" and run["duration_ms"] == 90_000
        assert run["finished_at"]

        assert client.get(f"/api/runs/{run_id}/log").json()["text"] == "hello"
        assert client.patch(f"/api/runs/{run_id}", json={"notes": "looked fine"}).json()["notes"] == "looked fine"

        outputs = client.get(f"/api/runs/{run_id}/outputs").json()
        assert outputs["exists"]
        assert {f["name"] for f in outputs["files"]} >= {"l2fc.csv", "collapsed_l2fc.csv"}
        assert all(f["exists"] is False for f in outputs["files"]), "nothing was actually produced"

        # An existing config.json silently wins over the form; preflight must say so.
        (build_dir / "config.json").write_text(json.dumps({"COUNT_THRESHOLD": "40"}))
        stale = client.post("/api/preflight", json=body).json()["stale_config"]
        assert [c["name"] for c in stale["conflicts"]] == ["COUNT_THRESHOLD"]

        archived = client.post("/api/runs", json={**body, "archive_existing_config": True}).json()
        assert archived["archived"] and not (build_dir / "config.json").exists()
        assert list(build_dir.glob("config.json.*.bak"))

    assert client.get("/api/catalog").json()["prismseq_root"]
    assert client.get("/api/runs").json()["total"] >= 2


class temp_root:
    """Point catalog at a throwaway PRISMSEQ_ROOT, and drop the scan cache."""

    def __init__(self, path):
        self.path = str(path)

    def __enter__(self):
        self.original = catalog.PRISMSEQ_ROOT
        catalog.PRISMSEQ_ROOT = self.path
        catalog._scan = None
        return self.path

    def __exit__(self, *exc):
        catalog.PRISMSEQ_ROOT = self.original
        catalog._scan = None


def test_build_paths_finds_builds_one_and_two_levels_down():
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        (root / "MTS_SEQ" / "MTS040").mkdir(parents=True)
        (root / "MTS_SEQ" / "MTS040" / "config.json").write_text("{}")
        (root / "MTS_SEQ" / "not_a_build").mkdir()
        (root / "PDEV099").mkdir()
        (root / "PDEV099" / "raw_counts_uncollapsed.csv.gz").write_bytes(b"")
        (root / "empty_parent").mkdir()

        with temp_root(root):
            found = {p["path"]: p for p in catalog.build_paths(refresh=True)}

        assert set(found) == {"MTS_SEQ/MTS040", "PDEV099"}, sorted(found)
        assert found["MTS_SEQ/MTS040"]["has_config"]
        assert found["MTS_SEQ/MTS040"]["parent"] == "MTS_SEQ"
        assert found["PDEV099"]["has_raw_counts"]
        assert found["PDEV099"]["parent"] == "", "a depth-1 build has no parent segment"


def test_resolve_rel_refuses_paths_outside_the_root():
    # BUILD_DIR drives file writes when archiving a config, so traversal out of
    # the root must not survive normalization.
    for bad in ("..", "../etc", "MTS_SEQ/../../etc"):
        try:
            catalog.resolve_rel(bad)
        except ValueError:
            continue
        raise AssertionError(f"resolve_rel accepted {bad!r}")

    assert catalog.resolve_rel("/MTS_SEQ/MTS040/") == "MTS_SEQ/MTS040"
    assert catalog.resolve_rel("MTS_SEQ/./MTS040") == "MTS_SEQ/MTS040"
    assert catalog.resolve_rel("") == ""


def test_infer_preset_prefers_the_directory_over_a_disagreeing_config():
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        # APS007 really does carry SCREEN_TYPE=MTS_SEQ on disk; presets.yml says
        # the directory is what locates the data, so the directory wins.
        (root / "APS_SEQ" / "APS007").mkdir(parents=True)
        (root / "APS_SEQ" / "APS007" / "config.json").write_text(
            json.dumps({"SCREEN_TYPE": "MTS_SEQ"}))
        (root / "CPS021").mkdir()
        (root / "CPS021" / "config.json").write_text(json.dumps({"SCREEN_TYPE": "CPS_SEQ"}))
        (root / "WHAT001").mkdir()

        with temp_root(root):
            assert catalog.infer_preset("APS_SEQ/APS007")[0] == "APS_SEQ"
            # No parent directory to go on, so fall back to the config.
            assert catalog.infer_preset("CPS021")[0] == "CPS_SEQ"
            # No config either, so the build name prefix.
            assert catalog.infer_preset("EPS123")[0] == "EPS_SEQ"
            # Nothing identifies it: default, and say so rather than guess.
            default_id = catalog.catalog()["presets"][0]["id"]
            assert catalog.infer_preset("WHAT001") == (default_id, "nothing in the path identifies a screen type")


def test_defaults_for_path_fills_name_screen_and_dir():
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        build = root / "MTS_SEQ" / "MTS041"
        build.mkdir(parents=True)

        with temp_root(root):
            got = catalog.defaults_for_path("MTS_SEQ/MTS041")
            assert got["preset"] == "MTS_SEQ" and got["inferred"]
            assert got["values"]["BUILD_NAME"] == "MTS041"
            assert got["values"]["SCREEN"] == "MTS041"
            assert got["values"]["BUILD_DIR"] == str(build)
            assert got["exists"] and not got["has_config"]

            # An explicit screen type wins and is reported as manual.
            forced = catalog.defaults_for_path("MTS_SEQ/MTS041", "CPS_SEQ")
            assert forced["preset"] == "CPS_SEQ" and not forced["inferred"]
            assert forced["values"]["SYNERGY"] is True

            # An empty path still yields a usable form, just no build identity.
            empty = catalog.defaults_for_path("")
            assert empty["path"] == "" and not empty["values"]["BUILD_NAME"]


def test_archiving_reports_an_unwritable_build_dir_instead_of_crashing():
    """Most build directories are owned by another account and mode 755."""
    from app import _archive_blockers

    with tempfile.TemporaryDirectory() as tmp:
        build_dir = Path(tmp) / "MTS100"
        build_dir.mkdir()
        assert _archive_blockers(build_dir) == [], "nothing to archive is not a blocker"

        (build_dir / "config.json").write_text("{}")
        assert _archive_blockers(build_dir) == [], "writable directory is fine"

        build_dir.chmod(0o555)
        try:
            blockers = _archive_blockers(build_dir)
            assert len(blockers) == 1 and "config.json" in blockers[0], blockers
            assert "not writable" in blockers[0]
        finally:
            build_dir.chmod(0o755)


def test_screen_type_store_round_trip_and_validation():
    """A bad store would break catalog() for every request, so validate first."""
    with tempfile.TemporaryDirectory() as tmp:
        store = Path(tmp) / "presets.local.yml"
        original_store, original_root = catalog.PRESETS_STORE, catalog.PRISMSEQ_ROOT
        catalog.PRESETS_STORE = store
        catalog.catalog.cache_clear()
        try:
            assert catalog.presets_source() == "shipped"
            shipped = catalog.catalog()["presets"]

            bad = [{"id": "MTS_SEQ", "label": "x", "dir": "MTS_SEQ",
                    "params": {"NOT_A_PARAM": "1"}}]
            assert any("NOT_A_PARAM" in p for p in catalog.validate_presets(bad))

            for case, expect in (
                ([{"id": "", "label": "x", "dir": "d", "params": {}}], "id is required"),
                ([{"id": "a b", "label": "x", "dir": "d", "params": {}}], "may only contain"),
                ([{"id": "A", "label": "x", "dir": "../etc", "params": {}}], "relative path"),
                ([{"id": "A", "label": "x", "dir": "d", "params": {"SCREEN_TYPE": "nope"}}], "must be one of"),
                ([{"id": "A", "label": "x", "dir": "d", "params": {"SYNERGY": "yes"}}], "true or false"),
                ([], "At least one"),
            ):
                assert any(expect in p for p in catalog.validate_presets(case)), (case, expect)

            dup = [dict(shipped[0]), dict(shipped[0])]
            assert any("duplicate" in p for p in catalog.validate_presets(dup))

            # A real edit survives the round trip and takes over from shipped.
            edited = structured_copy(shipped)
            edited[0]["params"]["COUNT_THRESHOLD"] = "77"
            edited.append({"id": "NEW_SEQ", "label": "New screen", "dir": "NEW_SEQ",
                           "params": {"SYNERGY": True}})
            catalog.save_presets(edited)

            assert catalog.presets_source() == "store" and store.exists()
            assert catalog.defaults_for(shipped[0]["id"])["COUNT_THRESHOLD"] == "77"
            assert catalog.defaults_for("NEW_SEQ")["SYNERGY"] is True

            # A new screen type is immediately usable for path inference.
            catalog.PRISMSEQ_ROOT = tmp
            assert catalog.infer_preset("NEW_SEQ/NEW001")[0] == "NEW_SEQ"

            catalog.reset_presets()
            assert catalog.presets_source() == "shipped"
            assert catalog.defaults_for(shipped[0]["id"])["COUNT_THRESHOLD"] != "77"
        finally:
            catalog.PRESETS_STORE, catalog.PRISMSEQ_ROOT = original_store, original_root
            catalog.catalog.cache_clear()


def structured_copy(value):
    return json.loads(json.dumps(value))


def test_separator_pseudo_params_are_not_reported_as_drift():
    """The parameter-separator plugin declares section headings as parameters."""
    import jenkins

    props = [{"parameterDefinitions": [
        {"name": "BUILD_NAME", "_class": "hudson.model.StringParameterDefinition"},
        {"name": "core_modules",
         "_class": "jenkins.plugins.parameter_separator.ParameterSeparatorDefinition",
         "type": "ParameterSeparatorDefinition"},
        {"name": "SCREEN_TYPE", "type": "ChoiceParameterDefinition"},
    ]}]
    assert jenkins.real_parameters(props) == ["BUILD_NAME", "SCREEN_TYPE"]
    assert jenkins.real_parameters([]) == []
    assert jenkins.real_parameters([{}]) == []


def test_params_yml_matches_the_groovy_job():
    """Drift here is the one failure mode that breaks every launch at once."""
    source = GROOVY.read_text()
    block = source[source.index("parameters {"):source.index("environment {")]
    declared = set(re.findall(r"(?:string|booleanParam|choice)\(\s*name:\s*'([^']+)'", block))
    ours = {p["name"] for p in catalog.catalog()["params"]}

    assert not (ours - declared), f"params.yml sends parameters the job does not declare: {sorted(ours - declared)}"
    assert not (declared - ours), f"job declares parameters the UI never sets: {sorted(declared - ours)}"


if __name__ == "__main__":
    for name, fn in sorted(globals().items()):
        if name.startswith("test_"):
            fn()
            print(f"ok  {name}")
    print("\nall checks passed")
