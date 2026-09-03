"""Lookups from cellDB, so the form can offer real values instead of free text.

This is the first thing in the UI that needs the clue API key. The key is read
from disk on this host and used only for these server-side lookups -- it is
never included in a response, and the pipeline still does its own cellDB reads
exactly as before.

Every lookup degrades to an empty list plus an error string: cellDB being down
must not stop someone launching a build with a value they already know.
"""

import os
import time
from pathlib import Path

import requests

API_URL = os.environ.get("CLUE_API_URL", "https://api.clue.io/api/").rstrip("/")
KEY_FILE = os.environ.get("CLUE_API_KEY_FILE", "/local/jenkins/.clue_api_key")
TIMEOUT = float(os.environ.get("CLUE_API_TIMEOUT", "20"))

# Ladders are added rarely; an hour of staleness is invisible in practice.
TTL = float(os.environ.get("SUSHI_UI_CELLDB_TTL", "3600"))

_cache: dict[str, tuple[float, list[str], str]] = {}


def _api_key() -> str:
    key = os.environ.get("CLUE_API_KEY", "")
    if key:
        return key.strip()
    try:
        return Path(KEY_FILE).read_text().strip()
    except OSError as exc:
        raise LookupError(f"no clue API key: {exc}") from exc


def _fetch_distinct(resource: str, field: str) -> list[str]:
    r = requests.get(
        f"{API_URL}/{resource}",
        params={"filter": '{"fields":["%s"],"limit":5000}' % field},
        headers={"user_key": _api_key()},
        timeout=TIMEOUT,
    )
    if not r.ok:
        raise LookupError(f"{resource} returned {r.status_code}")
    values = {str(row.get(field, "")).strip() for row in r.json() if row.get(field)}
    return sorted(v for v in values if v)


def _cached(name: str, resource: str, field: str) -> tuple[list[str], str]:
    hit = _cache.get(name)
    if hit and (time.monotonic() - hit[0]) < TTL:
        return hit[1], hit[2]
    try:
        values, error = _fetch_distinct(resource, field), ""
    except (requests.RequestException, LookupError, ValueError) as exc:
        # Keep a previous good answer if we have one; a blip should not empty
        # the dropdown mid-session.
        if hit:
            return hit[1], f"cellDB lookup failed ({exc}); showing the last known list."
        values, error = [], f"cellDB lookup failed: {exc}"
    _cache[name] = (time.monotonic(), values, error)
    return values, error


def cb_ladders() -> tuple[list[str], str]:
    """Control barcode ladder names, as CONTROL_BARCODE_META accepts them.

    make_config_file.groovy: "If the CBs exist in cellDB, this can simply be
    the lowercase cb_ladder name (ie, h-a). Otherwise, this must be a csv file
    located in the build directory." So these are suggestions, not a closed
    set -- a build using a local CSV must still be able to name it.
    """
    return _cached("cb_ladders", "v_control_barcodes", "set")


# Params whose text input gets a dropdown of real values, and where from.
SUGGESTION_SOURCES = {"CONTROL_BARCODE_META": cb_ladders}


def suggestions() -> dict:
    out = {}
    for param, source in SUGGESTION_SOURCES.items():
        values, error = source()
        out[param] = {"values": values, "error": error}
    return out
