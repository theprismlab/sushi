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

Sign-in is with a **Jenkins user ID and API token** (from
`/user/<you>/configure` in Jenkins — not the web password). This is not
decoration:

- This Jenkins gives anonymous *read* access but not `Job/Build`, so a build
  cannot be queued without credentials at all.
- It removes the self-reported name box. "Launched by" is whatever Jenkins says
  the account's full name is, so the run log cannot be filled in wrongly.

Credentials are exchanged for an opaque session id in an `HttpOnly` cookie and
held **in memory** — a restart signs everyone out, and nothing writes API
tokens to disk. That means **the service must run as a single worker**; with
several uvicorn workers a request can land on a process that has never seen the
session.

## Layout

    ui/backend/pyproject.toml  its own uv project -- see below
    ui/backend/params.yml      every parameter the Jenkins job declares: type, default, tier, help
    ui/backend/presets.yml     seed for the screen-type store; the UI edits the store, not this
    ui/backend/catalog.py      loads the yml, resolves screen types, walks the build tree, validates
    ui/backend/jenkins.py      Jenkins REST client (trigger, queue, status, console)
    ui/backend/session.py      Jenkins-credential sessions
    ui/backend/celldb.py       cellDB lookups for form suggestions
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
| `JENKINS_URL` | `http://localhost:8080` | this VM runs Jenkins on **8889** |
| `JENKINS_JOB_PATH` | `job/sushi` | this VM's job is **`job/run_sushi`**. Nested folders: `job/prism/job/sushi` |
| `JENKINS_USER` / `JENKINS_TOKEN` | unset | optional service account, used only when no user session applies |
| `PRISMSEQ_ROOT` | `/cmap/obelix/pod/prismSeq` | root of the build tree the browser walks |
| `SUSHI_UI_DB` | `ui/backend/sushi_ui.db` | put this somewhere backed up |
| `SUSHI_UI_PRESETS` | `ui/backend/presets.local.yml` | screen-type store; back this up too |
| `SUSHI_UI_SESSION_TTL` | `43200` (12h) | how long a sign-in lasts |
| `SUSHI_UI_SCAN_TTL` | `600` | build-tree scan cache. Typing an unlisted path still works, so staleness only costs suggestions |
| `CLUE_API_KEY_FILE` | `/local/jenkins/.clue_api_key` | for cellDB lookups; `CLUE_API_KEY` overrides |
| `CLUE_API_URL` | `https://api.clue.io/api/` | |
| `SUSHI_UI_CELLDB_TTL` | `3600` | cellDB lookups are cached this long |
| `CORS_ORIGINS` | `http://localhost:5173` | only needed for local dev against a remote backend |

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

## Deploying on vercingetorix-r8

The VM has Python 3.11.5, which satisfies `requires-python`, and `uv` on PATH.
It had **no node**; Node 22 LTS is now unpacked at `~/opt/node` (no root
needed), so the SPA can be built in place:

    export PATH="$HOME/opt/node/bin:$PATH"
    cd ui/frontend && npm ci && npm run build

Node is a build-time dependency only — the backend serves `frontend/dist`
directly, so there is no second web server to run.

    cd ui/backend && uv sync --no-dev --frozen

`--frozen` installs exactly what `uv.lock` pins and fails rather than
re-resolving; `--no-dev` leaves out httpx, which only `test_backend.py` needs.
`uv sync` creates `.venv` in that directory.

### Where it lives, and who must run it

The deployment is `/local/jenkins/sushi-ui` (group `xchipgrp`, setgid so new
files keep the group). The database and screen-type store sit **outside** the
checkout, at `/local/jenkins/sushi_ui.db` and
`/local/jenkins/sushi_ui_presets.yml`, so a reclone or `git clean` cannot take
them with it.

**Run it as `espresso`.** Archiving a stale `config.json` writes into the build
directory; those are owned by espresso and many are mode 755, so as any other
account the archive step fails and the launch is refused with a permissions
message. espresso is also the account Jenkins itself runs as.

    sudo -u espresso -i
    /local/jenkins/sushi-ui/bin/sushi-ui restart

`sushi-ui {start|stop|restart|status|log}` runs it detached with the right
environment for this host. It uses a pidfile rather than `pkill`, because the
pattern `uvicorn app:app` also matches the shell running the `pkill`; it also
falls back to whoever holds the port, so it can stop an instance the other
account started.

There is **no systemd unit yet** — installing one needs root, so a reboot
currently loses the service. The unit, with this host's real values:

    [Unit]
    Description=sushi pipeline UI
    After=network.target jenkins.service

    [Service]
    User=espresso
    WorkingDirectory=/local/jenkins/sushi-ui/ui/backend
    Environment=JENKINS_URL=http://localhost:8889
    Environment=JENKINS_JOB_PATH=job/run_sushi
    Environment=PRISMSEQ_ROOT=/cmap/obelix/pod/prismSeq
    Environment=SUSHI_UI_DB=/local/jenkins/sushi_ui.db
    Environment=SUSHI_UI_PRESETS=/local/jenkins/sushi_ui_presets.yml
    ExecStart=/local/jenkins/sushi-ui/ui/backend/.venv/bin/uvicorn app:app --host 0.0.0.0 --port 8100
    Restart=on-failure

    [Install]
    WantedBy=multi-user.target

Do **not** add `--workers` — sessions are in memory (see Identity).

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
