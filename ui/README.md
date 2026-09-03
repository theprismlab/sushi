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

    ui/backend/pyproject.toml  its own uv project -- see below
    ui/backend/params.yml      every parameter the Jenkins job declares: type, default, tier, help
    ui/backend/presets.yml     per-screen-type overrides (only what differs from stock)
    ui/backend/catalog.py      loads the yml, resolves presets, validates
    ui/backend/jenkins.py      Jenkins REST client (trigger, queue, status, console)
    ui/backend/db.py           SQLite run records
    ui/backend/app.py          FastAPI routes
    ui/frontend/               React + Vite SPA

The backend is a **separate uv project**, not part of the root `sushi` package.
The root `pyproject.toml` is what the pipeline image installs
(`pip install -e .` in the Dockerfile), so putting fastapi and uvicorn there
would pull them into every pipeline container — and conversely, syncing the
root project onto the web host would drag in polars, boto3, bigquery and a
`pyfarmhash` build the UI never touches. Two projects, two lockfiles, no
overlap in dependencies.

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
    uv sync
    uv run python test_backend.py             # self-check, no Jenkins needed
    uv run uvicorn app:app --reload --port 8000

    cd ui/frontend
    npm install && npm run dev                # http://localhost:5173, proxies /api to :8000

## Deploying

Two environments on `vercingetorix-r8`, one per branch, auto-deployed on push:

| Branch | Environment | Port | Run database |
| --- | --- | --- | --- |
| `main` | production | 8100 | `/local/jenkins/sushi_ui.db` |
| `develop` | develop | 8101 | `/local/jenkins/sushi_ui_develop.db` |

Push to `develop` and it is live within two minutes. Merge to `main` and
production follows. Nothing to run by hand.

**`develop` is not a sandbox.** There is one Jenkins job, so a build launched
from the develop instance is a real pipeline run writing into the real build
directory. Only the run history is separate. The UI shows an amber banner on
any non-production instance, and the header always shows the environment and
the deployed commit, so "did my change land?" is answerable from the page.

### How it works

`ui/deploy/deploy.sh <branch>` fetches and exits immediately if the branch has
not moved. Otherwise it resets the checkout, builds the SPA, syncs the venv and
restarts the unit — then polls `/api/health` for 30s. **If the new commit does
not come up healthy it resets to the previous commit, rebuilds and restarts**,
so an unattended deploy cannot leave an environment dead.

    ui/deploy/deploy.sh            per-environment deploy; idempotent, cron-safe
    ui/deploy/install.sh           one-time root bootstrap
    ui/deploy/sushi-ui@.service    systemd template, %i is the branch
    ui/deploy/env.main             production config
    ui/deploy/env.develop          develop config

Three decisions worth knowing:

- **Cron polling, not a webhook or hosted CI.** The VM is on an internal
  network at `10.200.96.63`; github.com cannot reach it, so a webhook or a
  GitHub-hosted runner has nothing to talk to. A self-hosted runner would work
  but is another daemon to keep alive for what `git fetch` answers in 200ms.
  Both cron entries exit silently when the branch has not moved.
- **The SPA is built in a `node:20` podman container.** The VM has no node and
  needs none at runtime; this avoids installing a toolchain on the host and
  avoids committing `dist/` to the repo. `node_modules` is gitignored, so
  `git reset --hard` leaves the npm cache intact between deploys.
- **`deploy.sh` is wrapped in a `{ ... }` block.** A deploy can update
  `deploy.sh` itself, and the brace forces bash to read the whole file before
  executing, so a self-update cannot corrupt the run in progress.

### First-time setup

    sudo /path/to/sushi/ui/deploy/install.sh

Clones `/opt/sushi-ui/<branch>` for each environment, installs the systemd
template and `/etc/cron.d/sushi-ui`, does the first deploy, and enables the
units at boot.

A branch that has not merged the UI yet is **skipped** rather than provisioned,
and gets no cron entry — so bootstrapping while only `develop` has the code
works, and you re-run `install.sh` once `main` catches up. Idempotent: re-run
it after changing the unit file, the schedule, or the set of environments.

### Operating it

    systemctl status sushi-ui@main sushi-ui@develop
    journalctl -u sushi-ui@develop -f            # application log
    tail -f /var/log/sushi-ui-deploy.log         # deploy log
    curl -s localhost:8100/api/health | jq       # environment, commit, jenkins, param drift

    sudo /opt/sushi-ui/main/ui/deploy/deploy.sh main --force    # redeploy unchanged
    sudo systemctl stop sushi-ui@develop                        # take an env down

To roll production back, revert on `main` and push — cron redeploys within two
minutes. The run database lives outside the checkout, so no deploy or rollback
can touch run history.

Ports and database paths live in `ui/deploy/env.<branch>`, which is read by
systemd rather than the app, so changing one needs a restart (a `--force`
deploy does that).

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
