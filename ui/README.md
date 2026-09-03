# sushi pipeline UI

A launcher, status mirror and run log for the sushi pipeline, replacing the
Jenkins parameter form as the day-to-day entry point.

Jenkins is still the executor. This service does not run modules, check out
code or write `config.json` — it triggers the existing `make_config_file`
pipeline over the Jenkins REST API and reads status and console output back.
That keeps one source of truth for orchestration and means the Jenkins form
remains a working fallback.

## What it adds over the form

- **The build directory comes first.** You browse the real `prismSeq` listings
  to a directory that exists, and the screen type, build name and screen are
  derived from it. The browser is keyboard-driven: type to narrow the current
  level, `↑`/`↓` to move, `Enter` to open a folder or select a build,
  `Backspace` on an empty filter to go up. Two or more characters also search
  every build path, so you can jump to a build without knowing which
  screen-type folder holds it.
- **Screen-type defaults you can edit.** The Screen types page owns what each
  type prefills. Values the type moved off the stock default are badged on the
  form, so the diff against "normal" is visible instead of buried in 90 fields.
- **Three tiers instead of one flat list.** Build identity, modules and
  controls are on screen; the ~60 filenames, column names and QC thresholds sit
  behind one disclosure.
- **Real values where they exist.** `CONTROL_BARCODE_META` offers the control
  barcode ladders actually registered in cellDB. It stays a text field, because
  the pipeline also accepts a CSV filename in the build directory.
- **Pipeline version is its own section, populated from GitHub.** `GIT_BRANCH`
  is a list of real branches and `COMMIT_ID` a list of real commits on the one
  you picked, subject and date included; with "use latest" on it names the head
  commit you are about to run. These decide which code executes, so they are
  not buried behind the advanced disclosure, and typing them freehand meant a
  typo surfaced only when Jenkins failed to check out. Both degrade to free
  text if GitHub cannot be reached.
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

## Identity

Sign-in is with your **Jenkins user ID and password**. This is not decoration:

- This Jenkins gives anonymous *read* access but not `Job/Build`, so a build
  cannot be queued without credentials at all.
- It removes the self-reported name box. "Launched by" is whatever Jenkins says
  the account's full name is, so the run log cannot be filled in wrongly.

Whatever is typed is passed straight through as HTTP basic auth to
`/me/api/json`, so an **API token** from `/user/<you>/configure` works in the
password field too. Some Jenkins configurations only accept tokens over the
API rather than passwords; if a correct password is rejected, that is why, and
the error message says so.

The credential is verified against `/me/api/json` specifically because that
endpoint answers 403 for anonymous and 401 for bad credentials — so a 200 means
genuinely authenticated. Since this Jenkins grants anonymous read, a check
against almost any other endpoint would let a wrong password "succeed" as
anonymous; there is an explicit guard against that too.

Credentials are exchanged for an opaque session id in an `HttpOnly` cookie and
held **in memory**, because the credential is not just checked at sign-in — it
is replayed to Jenkins to queue the build as you. A restart therefore signs
everyone out, and nothing is written to disk. It also means **the service must
run as a single worker**; with several uvicorn workers a request can land on a
process that has never seen the session.

## Layout

    ui/backend/pyproject.toml  its own uv project -- see below
    ui/backend/params.yml      every parameter the Jenkins job declares: type, default, tier, help
    ui/backend/presets.yml     seed for the screen-type store; the UI edits the store, not this
    ui/backend/catalog.py      loads the yml, resolves screen types, walks the build tree, validates
    ui/backend/jenkins.py      Jenkins REST client (trigger, queue, status, console)
    ui/backend/session.py      Jenkins-credential sessions
    ui/backend/celldb.py       cellDB lookups for form suggestions
    ui/backend/github.py       branch and commit lists for the version section
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

Only one direction of drift is fatal. A parameter we send that the job does not
declare kills the build; a parameter the job declares that we never send just
keeps its job default. The health check reports them separately. Note that the
parameter-separator plugin declares section headings (`core_modules`,
`do_not_edit`, …) as parameter definitions — those are filtered out, or they
show up as permanent phantom drift.

## Screen-type defaults

`presets.yml` in git is a **seed**. The first save on the Screen types page
writes `SUSHI_UI_PRESETS` (`presets.local.yml` beside the backend by default)
and that store is authoritative from then on; "Revert to shipped defaults"
deletes it and falls back to the seed. Put the store somewhere backed up.

Each type carries a `dir`, which is how a build path is mapped back to a screen
type. Inference order:

1. the first path segment matches a type's `dir` (`MTS_SEQ/MTS032` → `MTS_SEQ`)
2. otherwise `SCREEN_TYPE` in the build's existing `config.json`
3. otherwise the build name's first three characters
4. otherwise the first type, and the form says so

Directory before `config.json` is deliberate — see known issue 3 below.
Whatever it picks is a starting point; the form lets you override it.

Saves are validated before they are written (unknown parameter names, choice
and boolean values, duplicate ids, `dir` escaping the root). A bad store would
otherwise make every request fail until someone fixed a file by hand.

## Configuration

All via environment variables. The values below are the **defaults in code**;
what this VM actually needs is in the deploy section, because two of them
differ.

| Variable | Default | Notes |
| --- | --- | --- |
| `JENKINS_URL` | `http://localhost:8080` | this VM runs Jenkins on **8889**, from a bare .war as `espresso` with no systemd unit |
| `JENKINS_JOB_PATH` | `job/sushi` | this VM's job is **`job/run_sushi`**. Nested folders: `job/prism/job/sushi` |
| `JENKINS_PUBLIC_URL` | same as `JENKINS_URL` | where a **browser** reaches Jenkins. `JENKINS_URL` is localhost on the deploy host, which is the client's own laptop when followed as a link |
| `JENKINS_USER` / `JENKINS_TOKEN` | unset | optional service account, used only when no user session applies |
| `PRISMSEQ_ROOT` | `/cmap/obelix/pod/prismSeq` | root of the build tree the browser walks |
| `SUSHI_UI_DB` | `ui/backend/sushi_ui.db` | put this somewhere backed up |
| `SUSHI_UI_PRESETS` | `ui/backend/presets.local.yml` | screen-type store; back this up too |
| `SUSHI_UI_SESSION_TTL` | `43200` (12h) | how long a sign-in lasts |
| `SUSHI_UI_SCAN_TTL` | `600` | build-tree scan cache. Typing an unlisted path still works, so staleness only costs suggestions |
| `CLUE_API_KEY_FILE` | `/local/jenkins/.clue_api_key` | for cellDB lookups; `CLUE_API_KEY` overrides |
| `CLUE_API_URL` | `https://api.clue.io/api/` | |
| `SUSHI_UI_CELLDB_TTL` | `3600` | cellDB lookups are cached this long |
| `GITHUB_REPO` | `theprismlab/sushi` | source of the branch and commit lists |
| `GITHUB_TOKEN` | unset | only needed if the repo goes private or 60 req/hr unauthenticated is not enough |
| `SUSHI_UI_GITHUB_TTL` | `300` | branch and commit lookups are cached this long |
| `CORS_ORIGINS` | `http://localhost:5173` | only needed for local dev against a remote backend |
| `SUSHI_ENV` | `local` | `production` or `develop`; anything else banners itself as non-production |
| `SUSHI_COMMIT` / `SUSHI_BRANCH` | unset | set by the deploy job, since the image carries no `.git` |
| `SUSHI_SERVICE_ACCOUNT` | the process owner | set by the deploy job; in a container the process owner is `root`, not the host account permissions actually depend on |

The UI reads the clue API key for cellDB lookups. It is used server-side only
and never appears in a response — but it is a change from the original design,
in which the UI touched no key at all.

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

| Job | Branch | Environment | Port | State volume |
| --- | --- | --- | --- | --- |
| `sushi-ui-main` | `main` | production | 8100 | `sushi-ui-db-main` |
| `sushi-ui-develop` | `develop` | develop | 8101 | `sushi-ui-db-develop` |

Push to `develop` and it is live within two minutes. Merge to `main` and
production follows. Nothing to run by hand.

**`develop` is not a sandbox.** There is one pipeline Jenkins job, so a build
launched from the develop instance is a real pipeline run writing into the real
build directory. Only the run history and the screen-type store are separate.
The UI shows an amber banner on any non-production instance, and the header
shows the environment and deployed commit, so "did my change land?" is
answerable from the page.

### Why Jenkins and not systemd

We have no root on that VM, which rules out systemd units, `/etc/cron.d`,
`/opt` and `/var/log`. Jenkins already runs there — as `espresso`, from a bare
`jenkins.war` on port 8889 — with podman and the obelix mount, already polls
this repo, and already gives us build history and console logs. So the deploy
is a Jenkins job and needs no privileges we do not have.

It also solves the ownership problem for free. Archiving a stale `config.json`
writes into the build directory, and those are espresso-owned and mode 755, so
only espresso can do it. Jenkins runs as espresso and rootless podman maps
container root to the invoking user, so the container's writes land as espresso
without anyone needing `sudo -u espresso`.

    ui/Dockerfile            multi-stage: node builds the SPA, python serves it
    ui/.dockerignore         keeps the host .venv and dist/ out of the image
    ui/deploy/Jenkinsfile    the deploy pipeline, shared by both jobs

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

Four things about running a server this way:

- **`JENKINS_NODE_COOKIE=dontKillMe` on `podman run`.** Jenkins'
  ProcessTreeKiller reaps the build's descendants when it finishes, which
  includes `conmon` and would kill the container.
- **A `cron('H/15 * * * *')` trigger alongside `pollSCM`.** A detached
  container does not come back after a reboot, and nothing else on the box will
  restart it. The periodic run is the safety net; it is cheap because of the
  commit-label check.
- **`--network=host`.** The UI has to reach Jenkins on `localhost:8889`, and
  this also publishes the port without a rootless port mapping.
- **One uvicorn worker.** Sessions are in memory (see Identity), so the
  `CMD` deliberately has no `--workers`. A deploy also signs everyone out.

Everything that must survive a deploy is on the container's `/data` volume:
`SUSHI_UI_DB` and `SUSHI_UI_PRESETS`. Both would otherwise live inside the
image and silently reset on every deploy. The image carries no `.git`, so the
job passes `SUSHI_COMMIT` and `SUSHI_BRANCH` in; `/api/health` prefers those
and falls back to `git` when run from a checkout in dev.

### First-time setup

Create two pipeline jobs in Jenkins — no shell access needed:

1. **New Item → Pipeline**, named `sushi-ui-develop`.
2. **Pipeline → Definition: Pipeline script from SCM**, SCM `git`, repo
   `git@github.com:theprismlab/sushi.git`, **Branch Specifier** `*/develop`,
   **Script Path** `ui/deploy/Jenkinsfile`.
3. Save, then **Build Now** once. The triggers take over after that first run.
4. Repeat as `sushi-ui-main` with Branch Specifier `*/main`, once `main` has
   the UI code — and read the next section first, because something is already
   using port 8100.

### Operating it

    curl -s localhost:8101/api/health | jq    # environment, commit, jenkins, param drift
    podman ps --filter name=sushi-ui-
    podman logs -f sushi-ui-develop           # application log

Deploy history and console output are in the Jenkins job. To redeploy without a
code change, build with `FORCE_REDEPLOY`. To take an environment down,
`podman rm -f sushi-ui-develop` — but the 15-minute trigger brings it back, so
disable the job if you want it to stay down.

To roll production back, revert on `main` and push; the job redeploys within two
minutes. Run history and screen-type defaults live in a named podman volume, so
no deploy, rollback or `podman rm` touches them.

Ports and volumes are the `ENVIRONMENTS` map at the top of the Jenkinsfile.

### The hand-started instance on :8100

Before the Jenkins jobs existed, the UI was deployed by hand at
`/local/jenkins/sushi-ui` (group `xchipgrp`, setgid) and started with
`bin/sushi-ui {start|stop|restart|status|log}`, running as espresso with its
state at `/local/jenkins/sushi_ui.db` and `/local/jenkins/sushi_ui_presets.yml`
— outside the checkout, so a reclone could not take them with it. That script
uses a pidfile rather than `pkill`, because the pattern `uvicorn app:app` also
matches the shell running the `pkill`.

That instance is still what serves 8100, and it is the fallback if the Jenkins
jobs are not working. Two consequences:

- **`sushi-ui-main` cannot start until it is stopped**, because
  `--network=host` cannot bind a port something else holds. The deploy would
  fail its health check and roll back. `develop` on 8101 is unaffected.
- **Its state is not the volume's state.** To carry the existing run history
  and screen-type edits into the container, copy those two files into the
  `sushi-ui-db-main` volume before the first production deploy.

Node 22 LTS is unpacked at `~/opt/node` on the VM, which is how the SPA used to
be built in place. Deploys no longer need it — the image's node stage does that
— but it is still there for working on a checkout directly on the host.

## Security

Actions are authenticated; reading is not.

- **Launching and stopping a build require a Jenkins session.** That is
  enforced server-side, and Jenkins independently requires `Job/Build`.
- **The SPA is behind a sign-in wall**, but that is a UI gate, not access
  control. Every read endpoint — run history, parameters, console logs, the
  build tree — still answers an unauthenticated request. Anyone who can reach
  the port can read all of it with `curl`.

So the port must not be exposed beyond the internal network. It currently binds
`0.0.0.0:8100`; bind an internal interface or firewall the port. Locking the
read endpoints behind a session too is a small change if that is wanted.

## Known issues in the pipeline this UI surfaced

These are pre-existing and deliberately **not** changed here:

1. `make_config_file.groovy:452` branches on `params.JOIN_METADATA`, which is
   never declared as a parameter. It is always null, so
   `join_metadata/join_metadata.sh` can never run. The UI does not expose it,
   because sending an undeclared parameter would make Jenkins reject the build.
2. `config.json` in every build directory contains the plaintext clue API key,
   on a shared mount.
3. `APS007/config.json` has `SCREEN_TYPE=MTS_SEQ` and `AIR002/config.json` has
   `SCREEN_TYPE=APS_SEQ` — neither matches its directory. This is why screen
   type is inferred from the directory before the config file.
4. There is no `AIR_SEQ` value in the `SCREEN_TYPE` choice list, so the AIR
   type sends `APS_SEQ` as AIR002 did. Add the choice to the groovy and to
   `params.yml` to make it first class.
5. `LOW_ABUNDANCE_THRESHOLD` and `READ_DETECTION_LIMIT` are derived from
   `config.COUNT_THRESHOLD ?: params.COUNT_THRESHOLD`, so on a re-run against
   an existing `config.json` they follow the old value while
   `well_reads_threshold` and `nc_raw_count_threshold` follow the new one.
   Archiving the config avoids the split.
