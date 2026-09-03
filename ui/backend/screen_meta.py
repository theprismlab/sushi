"""Per-screen values from the same view create_sample_meta reads.

scripts/create_sample_meta/create_sample_meta.py builds the sample metadata
from https://api.theprismlab.org/api/v_seq_metadata filtered on project_code,
and renames `control_barcodes` to `cb_ladder` on the way. Reading that view
here means the form can offer the ladder the screen actually uses rather than
whatever the previous build happened to be set to.

A screen has exactly one ladder. Checked against MTS032, EPS008, APS007 and
CPS016: one distinct value each, matching their config.json. So more than one
value is a genuine anomaly worth surfacing, not something to pick from.
"""

import os
import time

import requests

from celldb import api_key  # same user_key; one key file for both services

VIEW = os.environ.get("SEQ_META_URL", "https://api.theprismlab.org/api/v_seq_metadata")
TIMEOUT = float(os.environ.get("SEQ_META_TIMEOUT", "40"))

# The metadata for a registered screen does not change during a session.
TTL = float(os.environ.get("SUSHI_UI_SCREEN_META_TTL", "900"))

# The API field name, before create_sample_meta renames it to cb_ladder.
FIELD = "control_barcodes"

_cache: dict[str, tuple[float, list[str]]] = {}


def control_barcodes(screen: str) -> list[str]:
    """Distinct control-barcode ladders recorded for a screen. Raises on failure."""
    hit = _cache.get(screen)
    if hit and (time.monotonic() - hit[0]) < TTL:
        return hit[1]

    # `fields` keeps the payload to one column; a screen can be a few thousand
    # rows and we only ever want the distinct set.
    query = '{"where":{"project_code":"%s"},"fields":["%s"]}' % (screen, FIELD)
    r = requests.get(VIEW, params={"filter": query},
                     headers={"user_key": api_key(), "Accept": "application/json"},
                     timeout=TIMEOUT)
    if not r.ok:
        raise LookupError(f"{VIEW} returned {r.status_code}")

    rows = r.json()
    if not isinstance(rows, list):
        raise LookupError(f"unexpected response for screen {screen!r}")

    values = sorted({str(row.get(FIELD, "")).strip() for row in rows if row.get(FIELD)})
    _cache[screen] = (time.monotonic(), values)
    return values
