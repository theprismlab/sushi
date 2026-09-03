"""FastAPI backend for the sushi pipeline UI.

Jenkins remains the executor. This service is a launcher, a status mirror and
the durable run log.
"""

import json
import os
import shutil
from datetime import datetime, timezone
from pathlib import Path

from fastapi import FastAPI, HTTPException, Query
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import FileResponse
from fastapi.staticfiles import StaticFiles
from pydantic import BaseModel, Field

import catalog
import db
import jenkins

app = FastAPI(title="sushi pipeline UI")

# No auth by design: this binds to the internal network only, same trust model
# as the Jenkins instance it fronts. Do not expose it publicly -- anyone who
# can reach it can launch a build and write into the build directory.
app.add_middleware(
    CORSMiddleware,
    allow_origins=os.environ.get("CORS_ORIGINS", "http://localhost:5173").split(","),
    allow_methods=["*"],
    allow_headers=["*"],
)

db.init()

# make_config_file.groovy loads an existing config.json and only fills in keys
# that are absent -- the file wins over the parameters we send. These are the
# exceptions it recomputes on every run.
RECOMPUTED_QC_KEYS = {"nc_raw_count_threshold", "well_reads_threshold"}


# --------------------------------------------------------------- catalog


@app.get("/api/catalog")
def get_catalog():
    data = catalog.catalog()
    return {
        "groups": data["groups"],
        "params": data["params"],
        "presets": [
            {"id": p["id"], "label": p["label"], "dir": p["dir"], "reference": p.get("reference")}
            for p in data["presets"]
        ],
        "module_params": catalog.module_names(),
        "prismseq_root": catalog.PRISMSEQ_ROOT,
    }


@app.get("/api/presets/{preset_id}/defaults")
def preset_defaults(preset_id: str, build_name: str = ""):
    try:
        values = catalog.defaults_for(preset_id, build_name.strip())
    except KeyError:
        raise HTTPException(404, f"unknown screen type {preset_id!r}")
    return {
        "values": values,
        # So the form can badge which fields the preset moved off the stock default.
        "from_preset": sorted(catalog.preset(preset_id)["params"]),
    }


@app.get("/api/builds")
def list_build_dirs(preset: str):
    """Existing build directories for a screen type, so the user picks instead of types."""
    try:
        parent = Path(catalog.build_dir_for(preset, ""))
    except KeyError:
        raise HTTPException(404, f"unknown screen type {preset!r}")
    if not parent.is_dir():
        return {"root": str(parent), "exists": False, "builds": []}

    builds = []
    for entry in parent.iterdir():
        if not entry.is_dir() or entry.name.startswith("."):
            continue
        builds.append({
            "name": entry.name,
            "has_raw_counts": (entry / "raw_counts_uncollapsed.csv.gz").exists(),
            "has_config": (entry / "config.json").exists(),
            "modified": _iso(entry.stat().st_mtime),
        })
    # Most recently touched first: the build someone is working on now is the
    # one they are most likely about to launch.
    builds.sort(key=lambda b: b["modified"], reverse=True)
    return {"root": str(parent), "exists": True, "builds": builds}


# -------------------------------------------------------------- preflight


class LaunchRequest(BaseModel):
    preset: str
    launched_by: str = Field(min_length=1)
    values: dict
    archive_existing_config: bool = False


@app.post("/api/preflight")
def preflight(req: LaunchRequest):
    """Everything we can check before spending an hour of compute."""
    try:
        catalog.preset(req.preset)
    except KeyError:
        raise HTTPException(404, f"unknown screen type {req.preset!r}")

    values = catalog.coerce(req.values)
    problems = catalog.validate(values)
    warnings = []

    build_dir = Path(str(values.get("BUILD_DIR", "")).rstrip("/"))
    if not str(build_dir):
        problems.append("Build directory is required.")
    elif not build_dir.is_dir():
        problems.append(f"Build directory does not exist: {build_dir}")
    else:
        raw = build_dir / str(values.get("RAW_COUNTS_UNCOLLAPSED", ""))
        if values.get("COLLATE_FASTQ_READS") and not raw.exists():
            problems.append(f"Collate FASTQ reads is on but {raw.name} is missing from the build directory.")
        if not values.get("CREATE_SAMPLE_META"):
            sample_meta = build_dir / str(values.get("SAMPLE_META", ""))
            if not sample_meta.exists():
                problems.append(
                    f"{sample_meta.name} is missing and 'Pull sample metadata from COMET' is off."
                )

    return {
        "problems": problems,
        "warnings": warnings,
        "stale_config": _stale_config(build_dir, values) if build_dir.is_dir() else None,
    }


def _stale_config(build_dir: Path, values: dict) -> dict | None:
    """Which submitted values an existing config.json would silently override.

    make_config_file.groovy treats an existing config.json as authoritative and
    only fills in missing keys. Without this check, editing a parameter in the
    UI for a re-run appears to work and then does nothing.
    """
    conflicts = []
    for filename, recomputed in (("config.json", set()), ("qc_params.json", RECOMPUTED_QC_KEYS)):
        path = build_dir / filename
        if not path.exists():
            continue
        try:
            existing = json.loads(path.read_text())
        except (OSError, json.JSONDecodeError) as exc:
            conflicts.append({"file": filename, "name": "-", "existing": f"unreadable: {exc}", "submitted": "-"})
            continue
        for name, submitted in values.items():
            if name not in existing or name in recomputed:
                continue
            if str(existing[name]) != str(submitted):
                conflicts.append({
                    "file": filename,
                    "name": name,
                    "existing": str(existing[name]),
                    "submitted": str(submitted),
                })
    if not conflicts:
        return None
    return {
        "conflicts": sorted(conflicts, key=lambda c: (c["file"], c["name"])),
        "note": (
            "An existing config.json takes precedence over the form: the pipeline only "
            "fills in keys it is missing. Archive it to make these values take effect. "
            "Module toggles are unaffected -- they are never stored in config.json."
        ),
    }


# ------------------------------------------------------------------ runs


@app.post("/api/runs")
def launch(req: LaunchRequest):
    check = preflight(req)
    if check["problems"]:
        raise HTTPException(422, {"problems": check["problems"]})

    values = catalog.coerce(req.values)
    build_dir = Path(str(values["BUILD_DIR"]).rstrip("/"))

    archived = []
    if req.archive_existing_config:
        stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
        for filename in ("config.json", "qc_params.json"):
            path = build_dir / filename
            if path.exists():
                backup = path.with_name(f"{filename}.{stamp}.bak")
                shutil.copy2(path, backup)
                path.unlink()
                archived.append(backup.name)

    modules = {name: values.get(name) for name in catalog.module_names()}
    run_id = db.create(
        launched_by=req.launched_by.strip(),
        preset=req.preset,
        build_name=str(values["BUILD_NAME"]),
        screen=str(values.get("SCREEN", "")),
        screen_type=str(values.get("SCREEN_TYPE", "")),
        build_dir=str(build_dir),
        params=values,
        modules=modules,
        git_branch=str(values.get("GIT_BRANCH", "")),
    )

    try:
        queue_url = jenkins.trigger(values)
    except jenkins.JenkinsError as exc:
        db.update(run_id, status="ERROR", error=str(exc), finished_at=db.now())
        raise HTTPException(502, f"Jenkins refused the build: {exc}")

    db.update(run_id, queue_url=queue_url)
    return {"id": run_id, "queue_url": queue_url, "archived": archived}


@app.get("/api/runs")
def get_runs(limit: int = Query(50, le=200), offset: int = 0, preset: str | None = None,
             status: str | None = None, build_name: str | None = None,
             launched_by: str | None = None):
    # Refresh on read rather than from a background poller: a run that finishes
    # while nobody is looking is correct the moment somebody looks.
    for run_id in db.unfinished_ids():
        _refresh(run_id)
    runs, total = db.list_runs(limit=limit, offset=offset, preset=preset, status=status,
                               build_name=build_name, launched_by=launched_by)
    return {"runs": [_public(r) for r in runs], "total": total}


@app.get("/api/runs/{run_id}")
def get_run(run_id: int):
    run = _refresh(run_id) or db.get(run_id)
    if not run:
        raise HTTPException(404, f"no run {run_id}")
    return _public(run)


class Notes(BaseModel):
    notes: str


@app.patch("/api/runs/{run_id}")
def set_notes(run_id: int, body: Notes):
    if not db.get(run_id):
        raise HTTPException(404, f"no run {run_id}")
    db.update(run_id, notes=body.notes)
    return _public(db.get(run_id))


@app.post("/api/runs/{run_id}/stop")
def stop_run(run_id: int):
    run = db.get(run_id)
    if not run:
        raise HTTPException(404, f"no run {run_id}")
    if not run["build_number"]:
        raise HTTPException(409, "build has no number yet; it is still queued")
    try:
        jenkins.stop(run["build_number"])
    except jenkins.JenkinsError as exc:
        raise HTTPException(502, str(exc))
    return _public(_refresh(run_id) or db.get(run_id))


@app.get("/api/runs/{run_id}/log")
def get_log(run_id: int, start: int = 0):
    run = db.get(run_id)
    if not run:
        raise HTTPException(404, f"no run {run_id}")
    if not run["build_number"]:
        return {"text": "", "start": start, "next": start, "more": True,
                "message": "Waiting for Jenkins to assign an executor."}
    try:
        return jenkins.console(run["build_number"], start)
    except jenkins.JenkinsError as exc:
        # Jenkins rotates builds out; the pipeline also copies the console log
        # into the build directory, so fall back to that.
        archived = Path(run["build_dir"]) / "logs" / "console_output.log"
        if start == 0 and archived.exists():
            text = archived.read_text(errors="replace")
            return {"text": text, "start": 0, "next": len(text), "more": False,
                    "message": f"Jenkins log unavailable ({exc}); showing logs/console_output.log."}
        raise HTTPException(502, str(exc))


@app.get("/api/runs/{run_id}/outputs")
def get_outputs(run_id: int):
    """What the run actually produced, so 'outcome' is more than a status word."""
    run = db.get(run_id)
    if not run:
        raise HTTPException(404, f"no run {run_id}")
    build_dir = Path(run["build_dir"])
    if not build_dir.is_dir():
        return {"build_dir": str(build_dir), "exists": False, "files": []}

    expected = [
        (name, run["params"].get(name))
        for name in ("PRISM_BARCODE_COUNTS", "UNKNOWN_BARCODE_COUNTS", "ANNOTATED_COUNTS",
                     "FILTERED_COUNTS", "NORMALIZED_COUNTS", "LFC", "COLLAPSED_LFC",
                     "SYNERGY_FILE", "SKIPPED_WELLS")
    ]
    files = []
    for param, filename in expected:
        if not filename:
            continue
        path = build_dir / filename
        stat = path.stat() if path.exists() else None
        files.append({
            "param": param,
            "name": filename,
            "exists": stat is not None,
            "size": stat.st_size if stat else None,
            "modified": _iso(stat.st_mtime) if stat else None,
            "stale": bool(stat and run["started_at"] and _iso(stat.st_mtime) < run["started_at"]),
        })
    qc = build_dir / "qc_tables"
    return {
        "build_dir": str(build_dir),
        "exists": True,
        "files": files,
        "qc_tables": sorted(p.name for p in qc.iterdir()) if qc.is_dir() else [],
    }


# ----------------------------------------------------------------- health


@app.get("/api/health")
def health():
    """Also surfaces drift between params.yml and the pipeline's parameters block."""
    result = {"jenkins_url": jenkins.BASE, "job_path": jenkins.JOB_PATH,
              "authenticated": jenkins.AUTH is not None, "db": str(db.DB_PATH)}
    try:
        declared = set(jenkins.declared_params())
        ours = {p["name"] for p in catalog.catalog()["params"]}
        result["jenkins_reachable"] = True
        # Sending an undeclared parameter makes Jenkins reject the whole build.
        result["params_we_send_that_job_rejects"] = sorted(ours - declared)
        result["job_params_not_in_ui"] = sorted(declared - ours)
    except (jenkins.JenkinsError, Exception) as exc:  # noqa: BLE001 - network shapes vary
        result["jenkins_reachable"] = False
        result["jenkins_error"] = str(exc)
    return result


# ---------------------------------------------------------------- helpers


def _refresh(run_id: int) -> dict | None:
    """Pull current state for a run from Jenkins. Cheap and idempotent."""
    run = db.get(run_id)
    if not run or run["status"] in db.TERMINAL:
        return run

    number = run["build_number"]
    try:
        if not number:
            if not run["queue_url"]:
                return run
            number, cancel_reason = jenkins.queue_build_number(run["queue_url"])
            if cancel_reason:
                db.update(run_id, status="ABORTED", error=cancel_reason, finished_at=db.now())
                return db.get(run_id)
            if not number:
                return run  # still queued
            db.update(run_id, build_number=number)

        info = jenkins.build(number)
    except jenkins.JenkinsError as exc:
        db.update(run_id, error=str(exc))
        return db.get(run_id)

    fields = {"error": None}
    if info.get("timestamp"):
        fields["started_at"] = _iso(info["timestamp"] / 1000)
    if info.get("building"):
        fields["status"] = "RUNNING"
    else:
        fields["status"] = info.get("result") or "ERROR"
        fields["duration_ms"] = info.get("duration") or None
        if info.get("timestamp") and info.get("duration"):
            fields["finished_at"] = _iso((info["timestamp"] + info["duration"]) / 1000)
    db.update(run_id, **fields)
    return db.get(run_id)


def _public(run: dict) -> dict:
    run = dict(run)
    # API_KEY is not a UI parameter, but the pipeline writes it into config.json;
    # scrub defensively so a future catalog change cannot leak it through here.
    run["params"] = {k: v for k, v in run["params"].items() if "KEY" not in k and "TOKEN" not in k}
    run["jenkins_url"] = jenkins.build_url(run["build_number"]) if run["build_number"] else None
    run["enabled_modules"] = sorted(k for k, v in run["modules"].items() if v)
    return run


def _iso(epoch_seconds: float) -> str:
    return datetime.fromtimestamp(epoch_seconds, timezone.utc).isoformat(timespec="seconds")


# Serve the built SPA if it is there; in dev the Vite server handles this.
_dist = Path(__file__).parent.parent / "frontend" / "dist"
if _dist.is_dir():
    app.mount("/assets", StaticFiles(directory=_dist / "assets"), name="assets")

    @app.get("/{path:path}")
    def spa(path: str):
        candidate = _dist / path
        if path and candidate.is_file():
            return FileResponse(candidate)
        return FileResponse(_dist / "index.html")
