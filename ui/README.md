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

Two Jenkins jobs, one per branch, each running the UI as a detached podman
container:

| Job | Branch | Environment | Port | Run database |
| --- | --- | --- | --- | --- |
| `sushi-ui-main` | `main` | production | 8100 | volume `sushi-ui-db-main` |
| `sushi-ui-develop` | `develop` | develop | 8101 | volume `sushi-ui-db-develop` |

Push to `develop` and it is live within two minutes. Merge to `main` and
production follows. Nothing to run by hand.

**`develop` is not a sandbox.** There is one pipeline Jenkins job, so a build
launched from the develop instance is a real pipeline run writing into the real
build directory. Only the run history is separate. The UI shows an amber banner
on any non-production instance, and the header always shows the environment and
deployed commit, so "did my change land?" is answerable from the page.

### Why Jenkins and not systemd

We have no root on that VM, which rules out systemd units, `/etc/cron.d`,
`/opt` and `/var/log`. Jenkins already runs there as a user with podman and the
obelix mount, already polls this repo, and already gives us build history and
console logs — so the deploy is a Jenkins job and needs no privileges we do not
have.

    ui/Dockerfile            multi-stage: node builds the SPA, python serves it
    ui/.dockerignore          keeps the host .venv and dist/ out of the image
    ui/deploy/Jenkinsfile     the deploy pipeline, shared by both jobs

### How it works

[Jenkinsfile](deploy/Jenkinsfile) derives the environment from the branch it
was checked out on, so the two jobs differ only in SCM configuration. It labels
each container with the commit it is running, so the common case — nothing
changed — skips straight to a health check and finishes in seconds.

On a real change: build the image, replace the container, then poll
`/api/health` for 30s. **If the new commit is unhealthy it restarts the
previous image and fails the build**, so an unattended deploy cannot leave an
environment dead. Per-commit image tags are the rollback path, which is why
`post` only prunes dangling images.

Three things about running a server from a Jenkins job:

- **`JENKINS_NODE_COOKIE=dontKillMe` on `podman run`.** Jenkins'
  ProcessTreeKiller reaps the build's descendants when it finishes, which
  includes `conmon` and would kill the container.
- **A `cron('H/15 * * * *')` trigger alongside `pollSCM`.** A detached
  container does not come back after a reboot, and nothing else on the box will
  restart it. The periodic run is the safety net; it is cheap because of the
  commit-label check.
- **`--network=host`.** The UI has to reach Jenkins on `localhost:8080`, and
  this also publishes the port without a rootless port mapping.

The image carries no `.git`, so the job passes `SUSHI_COMMIT` and
`SUSHI_BRANCH` in as environment variables; `/api/health` prefers those and
falls back to `git` when run from a checkout in dev.

### First-time setup

Create two pipeline jobs in Jenkins — no shell access needed:

1. **New Item → Pipeline**, named `sushi-ui-develop`.
2. **Pipeline → Definition: Pipeline script from SCM**, SCM `git`, repo
   `git@github.com:theprismlab/sushi.git`, **Branch Specifier** `*/develop`,
   **Script Path** `ui/deploy/Jenkinsfile`.
3. Save, then **Build Now** once. The triggers in the Jenkinsfile take over
   after that first run.
4. Repeat as `sushi-ui-main` with Branch Specifier `*/main`, once `main` has
   the UI code.

### Operating it

    curl -s localhost:8101/api/health | jq    # environment, commit, jenkins, param drift
    podman ps --filter name=sushi-ui-
    podman logs -f sushi-ui-develop           # application log

Deploy history and console output are in the Jenkins job. To redeploy without a
code change, build with `FORCE_REDEPLOY`. To take an environment down,
`podman rm -f sushi-ui-develop` — but the 15-minute trigger will bring it back,
so disable the job if you want it to stay down.

To roll production back, revert on `main` and push; the job redeploys within two
minutes. Run history lives in a named podman volume, so no deploy, rollback or
`podman rm` touches it.

Ports and volumes are the `ENVIRONMENTS` map at the top of the Jenkinsfile.

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
