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

    state = {"number": None, "building": True, "result": None, "sent": None}
    jenkins.trigger = lambda params: (state.__setitem__("sent", params), "http://j/queue/item/9")[1]
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
        body = {"preset": "MTS_SEQ", "launched_by": "tester", "values": values}

        pre = client.post("/api/preflight", json=body).json()
        assert pre["problems"] == [], pre["problems"]
        assert pre["stale_config"] is None

        # A missing input file must block the launch, not fail an hour in.
        (build_dir / "sample_meta.csv").unlink()
        blocked = client.post("/api/preflight", json=body).json()
        assert any("sample_meta.csv" in p for p in blocked["problems"]), blocked
        assert client.post("/api/runs", json=body).status_code == 422
        (build_dir / "sample_meta.csv").write_text("x\n")

        run_id = client.post("/api/runs", json=body).json()["id"]
        assert state["sent"]["SCREEN_TYPE"] == "MTS_SEQ"
        assert state["sent"]["COUNT_THRESHOLD"] == "100"

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
