"""Run records.

Jenkins already stores build parameters and console logs, but it rotates old
builds out and -- more importantly -- make_config_file.groovy never writes the
module toggles into config.json. So today you cannot tell from a build
directory which modules actually ran. This table is the durable record.
"""

import json
import os
import sqlite3
from contextlib import contextmanager
from datetime import datetime, timezone
from pathlib import Path

DB_PATH = Path(os.environ.get("SUSHI_UI_DB", Path(__file__).parent / "sushi_ui.db"))

SCHEMA = """
CREATE TABLE IF NOT EXISTS runs (
    id            INTEGER PRIMARY KEY AUTOINCREMENT,
    created_at    TEXT    NOT NULL,
    launched_by   TEXT    NOT NULL,
    preset        TEXT    NOT NULL,
    build_name    TEXT    NOT NULL,
    screen        TEXT,
    screen_type   TEXT,
    build_dir     TEXT    NOT NULL,
    params_json   TEXT    NOT NULL,
    modules_json  TEXT    NOT NULL,
    git_branch    TEXT,
    queue_url     TEXT,
    build_number  INTEGER,
    status        TEXT    NOT NULL,
    started_at    TEXT,
    finished_at   TEXT,
    duration_ms   INTEGER,
    notes         TEXT    NOT NULL DEFAULT '',
    error         TEXT
);
CREATE INDEX IF NOT EXISTS runs_created_idx ON runs (created_at DESC);
CREATE INDEX IF NOT EXISTS runs_build_idx   ON runs (build_name);
"""

TERMINAL = ("SUCCESS", "FAILURE", "ABORTED", "UNSTABLE", "NOT_BUILT", "ERROR")


def now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


@contextmanager
def connect():
    conn = sqlite3.connect(DB_PATH, timeout=15)
    conn.row_factory = sqlite3.Row
    try:
        yield conn
        conn.commit()
    finally:
        conn.close()


def init() -> None:
    with connect() as conn:
        conn.executescript(SCHEMA)


def create(*, launched_by, preset, build_name, screen, screen_type, build_dir,
           params, modules, git_branch) -> int:
    with connect() as conn:
        cur = conn.execute(
            """INSERT INTO runs (created_at, launched_by, preset, build_name, screen,
                                 screen_type, build_dir, params_json, modules_json,
                                 git_branch, status)
               VALUES (?,?,?,?,?,?,?,?,?,?, 'QUEUED')""",
            (now(), launched_by, preset, build_name, screen, screen_type, build_dir,
             json.dumps(params, sort_keys=True), json.dumps(modules, sort_keys=True), git_branch),
        )
        return cur.lastrowid


def update(run_id: int, **fields) -> None:
    if not fields:
        return
    assignments = ", ".join(f"{k} = ?" for k in fields)
    with connect() as conn:
        conn.execute(f"UPDATE runs SET {assignments} WHERE id = ?", (*fields.values(), run_id))


def get(run_id: int) -> dict | None:
    with connect() as conn:
        row = conn.execute("SELECT * FROM runs WHERE id = ?", (run_id,)).fetchone()
    return _row(row) if row else None


def list_runs(*, limit=50, offset=0, preset=None, status=None, build_name=None,
              launched_by=None) -> tuple[list[dict], int]:
    clauses, args = [], []
    for column, value in (("preset", preset), ("status", status), ("launched_by", launched_by)):
        if value:
            clauses.append(f"{column} = ?")
            args.append(value)
    if build_name:
        clauses.append("build_name LIKE ?")
        args.append(f"%{build_name}%")
    where = f"WHERE {' AND '.join(clauses)}" if clauses else ""

    with connect() as conn:
        total = conn.execute(f"SELECT COUNT(*) FROM runs {where}", args).fetchone()[0]
        rows = conn.execute(
            f"SELECT * FROM runs {where} ORDER BY id DESC LIMIT ? OFFSET ?",
            (*args, limit, offset),
        ).fetchall()
    return [_row(r) for r in rows], total


def unfinished_ids() -> list[int]:
    placeholders = ",".join("?" * len(TERMINAL))
    with connect() as conn:
        rows = conn.execute(
            f"SELECT id FROM runs WHERE status NOT IN ({placeholders})", TERMINAL
        ).fetchall()
    return [r["id"] for r in rows]


def _row(row: sqlite3.Row) -> dict:
    run = dict(row)
    run["params"] = json.loads(run.pop("params_json"))
    run["modules"] = json.loads(run.pop("modules_json"))
    return run
