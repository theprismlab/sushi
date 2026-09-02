# sushi pipeline UI

A launcher, status mirror and run log for the sushi pipeline, replacing the
Jenkins parameter form as the day-to-day entry point.

Jenkins is still the executor. This service does not run modules, check out
code or write `config.json` — it triggers the existing `make_config_file`
pipeline over the Jenkins REST API and reads status and console output back.
That keeps one source of truth for orchestration and means the Jenkins form
remains a working fallback.

## What it adds over the form

- **Screen-type presets.** Pick `MTS_SEQ` / `EPS_SEQ` / `APS_SEQ` / `AIR_SEQ` /
  `CPS_SEQ` and every parameter is prefilled from the most recent build of that
  type. Values the preset moved off the stock default are badged, so the diff
  against "normal" is visible instead of buried in 90 fields.
- **Three tiers instead of one flat list.** Build identity, modules and
  controls are on screen; the ~60 filenames, column names and QC thresholds sit
  behind one disclosure.
- **Preflight.** Before a build is queued: required fields, numeric fields,
  module dependencies (`AUC_BIOMARKER` without `DRC`, `GENERATE_QC_TABLES_2`
  without `COMPUTE_LFC`, …), whether the build directory exists, and whether
  the input files those modules need are actually present.
- **The stale-config warning.** `make_config_file.groovy` treats an existing
  `config.json` as authoritative and only fills in keys it is missing. Editing
  a parameter for a re-run therefore appears to work and silently does nothing.
  Preflight diffs the form against the file on disk, names every parameter that
  would be overridden, and offers to archive the file.
- **Run history.** Who launched it, the full parameter set, **which modules
  ran**, status, duration, console log, which output files exist, and free-text
  notes. Module toggles are the important part: Jenkins never writes them into
  `config.json`, so before this there was no record of which steps a build
  actually performed.

## Layout

    ui/backend/params.yml     every parameter the Jenkins job declares: type, default, tier, help
    ui/backend/presets.yml    per-screen-type overrides (only what differs from stock)
    ui/backend/catalog.py     loads the yml, resolves presets, validates
    ui/backend/jenkins.py     Jenkins REST client (trigger, queue, status, console)
    ui/backend/db.py          SQLite run records
    ui/backend/app.py         FastAPI routes
    ui/frontend/              React + Vite SPA

`params.yml` must stay in sync with the `parameters { ... }` block of
`scripts/make_config_file.groovy`: Jenkins rejects an entire build if it is
sent a parameter the job does not declare. `test_backend.py` parses the groovy
and asserts both directions, and `/api/health` reports drift at runtime.

## Configuration

All via environment variables; every one has a working default.

| Variable | Default | Notes |
| --- | --- | --- |
| `JENKINS_URL` | `http://localhost:8080` | |
| `JENKINS_JOB_PATH` | `job/sushi` | URL path of the job. Nested folders: `job/prism/job/sushi` |
| `JENKINS_USER` / `JENKINS_TOKEN` | unset | Needed if the Jenkins instance requires auth to trigger builds |
| `PRISMSEQ_ROOT` | `/cmap/obelix/pod/prismSeq` | Used to derive build directories and list existing builds |
| `SUSHI_UI_DB` | `ui/backend/sushi_ui.db` | Put this somewhere backed up |
| `CORS_ORIGINS` | `http://localhost:5173` | Only needed for local dev against a remote backend |

## Development

    cd ui/backend
    python3 -m venv .venv && .venv/bin/pip install -r requirements.txt
    .venv/bin/python test_backend.py          # self-check, no Jenkins needed
    .venv/bin/uvicorn app:app --reload --port 8000

    cd ui/frontend
    npm install && npm run dev                # http://localhost:5173, proxies /api to :8000

## Deploying on vercingetorix-r8

Prerequisites, neither of which is RHEL 8's default:

- **Python 3.10+.** The stock `python3` on RHEL 8 is 3.6, which pydantic 2 does
  not support. Install the AppStream module (`dnf module install python:3.11`)
  and build the venv with that interpreter.
- **Node 18+**, only to build the SPA. It is not needed at runtime — if you'd
  rather not put node on the VM, run `npm run build` anywhere and rsync
  `ui/frontend/dist/` across.

Build the SPA once; the backend serves `frontend/dist` directly, so there is no
second web server to run.

    cd /opt/sushi-ui/ui/frontend && npm ci && npm run build
    cd /opt/sushi-ui/ui/backend && python3 -m venv .venv && .venv/bin/pip install -r requirements.txt

`/etc/systemd/system/sushi-ui.service`:

    [Unit]
    Description=sushi pipeline UI
    After=network.target jenkins.service

    [Service]
    User=jenkins
    WorkingDirectory=/opt/sushi-ui/ui/backend
    Environment=JENKINS_URL=http://localhost:8080
    Environment=JENKINS_JOB_PATH=job/sushi
    Environment=PRISMSEQ_ROOT=/cmap/obelix/pod/prismSeq
    Environment=SUSHI_UI_DB=/local/jenkins/sushi_ui.db
    ExecStart=/opt/sushi-ui/ui/backend/.venv/bin/uvicorn app:app --host 0.0.0.0 --port 8100
    Restart=on-failure

    [Install]
    WantedBy=multi-user.target

Run as `jenkins` (or a user in its group): archiving a stale `config.json`
writes into the build directory, and the service needs to read
`PRISMSEQ_ROOT`.

## Security

**There is no authentication.** Anyone who can reach the port can launch a
build and archive files in a build directory. This matches the trust model of
the Jenkins instance it fronts, but it means the port must not be exposed
beyond the internal network — bind to an internal interface or firewall 8100.
Identity ("launched by") is self-reported and stored in the browser; it is a
record-keeping convenience, not an access control.

The UI never handles `API_KEY`; the pipeline reads it from
`/local/jenkins/.clue_api_key` as before. Responses strip any parameter whose
name contains `KEY` or `TOKEN` as a defence against a future catalog change.

## Known issues in the pipeline this UI surfaced

These are pre-existing and deliberately **not** changed here:

1. `make_config_file.groovy:452` branches on `params.JOIN_METADATA`, which is
   never declared as a parameter. It is always null, so
   `join_metadata/join_metadata.sh` can never run. The UI does not expose it,
   because sending an undeclared parameter would make Jenkins reject the build.
2. `config.json` in every build directory contains the plaintext clue API key,
   on a shared mount.
3. `APS007/config.json` has `SCREEN_TYPE=MTS_SEQ` and `AIR002/config.json` has
   `SCREEN_TYPE=APS_SEQ` — neither matches its directory. `presets.yml` assumes
   the APS one was a mistake and uses `APS_SEQ`; see the notes in that file.
4. There is no `AIR_SEQ` value in the `SCREEN_TYPE` choice list, so the AIR
   preset sends `APS_SEQ` as AIR002 did. Add the choice to the groovy and to
   `params.yml` to make it first class.
5. `LOW_ABUNDANCE_THRESHOLD` and `READ_DETECTION_LIMIT` are derived from
   `config.COUNT_THRESHOLD ?: params.COUNT_THRESHOLD`, so on a re-run against
   an existing `config.json` they follow the old value while
   `well_reads_threshold` and `nc_raw_count_threshold` follow the new one.
   Archiving the config avoids the split.
