"""Real branches and commits for the pipeline-version fields.

GIT_BRANCH and COMMIT_ID decide which code the pipeline runs, and typing them
freehand means a typo is only discovered when Jenkins fails to check out -- or
worse, silently runs the wrong revision.

theprismlab/sushi is public, so this needs no credentials. GITHUB_TOKEN is
honoured if set, which is only needed if the repo goes private or the
unauthenticated 60 requests/hour becomes a problem. With the cache below a
whole day of form-filling is a handful of requests.
"""

import os
import time

import requests

API = os.environ.get("GITHUB_API", "https://api.github.com").rstrip("/")
REPO = os.environ.get("GITHUB_REPO", "theprismlab/sushi")
TOKEN = os.environ.get("GITHUB_TOKEN", "")
TIMEOUT = float(os.environ.get("GITHUB_TIMEOUT", "15"))

# Branches move rarely enough that five minutes of staleness is invisible, and
# the fields stay free text, so a missing entry is never a dead end.
TTL = float(os.environ.get("SUSHI_UI_GITHUB_TTL", "300"))

# Enough to cover "the commit I made this morning" without paging.
COMMIT_LIMIT = int(os.environ.get("SUSHI_UI_GITHUB_COMMITS", "30"))

_cache: dict[str, tuple[float, object]] = {}


class GitHubError(RuntimeError):
    pass


def _get(path: str, not_found: str | None = None, **params):
    headers = {"Accept": "application/vnd.github+json"}
    if TOKEN:
        headers["Authorization"] = f"Bearer {TOKEN}"
    # No trailing slash when path is empty: GitHub 404s /repos/<owner>/<repo>/
    # even though /repos/<owner>/<repo> is fine.
    url = f"{API}/repos/{REPO}" + (f"/{path.lstrip('/')}" if path else "")
    try:
        r = requests.get(url, headers=headers, params=params, timeout=TIMEOUT)
    except requests.RequestException as exc:
        raise GitHubError(f"could not reach {API}: {exc}") from exc

    if r.status_code == 403 and "rate limit" in r.text.lower():
        raise GitHubError("GitHub rate limit reached; set GITHUB_TOKEN to raise it")
    if r.status_code == 404:
        raise GitHubError(not_found or f"{REPO} not found (private repo without GITHUB_TOKEN?)")
    if not r.ok:
        raise GitHubError(f"GitHub returned {r.status_code} for {path}")
    try:
        return r.json()
    except ValueError as exc:
        raise GitHubError("GitHub did not return JSON") from exc


def _cached(key: str, produce):
    hit = _cache.get(key)
    if hit and time.time() - hit[0] < TTL:
        return hit[1]
    value = produce()
    _cache[key] = (time.time(), value)
    return value


def branches() -> list[str]:
    """Branch names, default branch first, then alphabetical."""
    def produce():
        names = [b["name"] for b in _get("branches", per_page=100)]
        default = _get("")["default_branch"]
        rest = sorted(n for n in names if n != default)
        return ([default] + rest) if default in names else rest
    return _cached("branches", produce)


def commits(branch: str) -> list[dict]:
    """Recent commits on a branch, newest first."""
    def produce():
        return [
            {
                "sha": c["sha"][:7],
                "full_sha": c["sha"],
                # First line only: the form shows these in a dropdown.
                "subject": (c["commit"]["message"] or "").splitlines()[0][:100],
                "author": (c["commit"].get("author") or {}).get("name", ""),
                "date": (c["commit"].get("author") or {}).get("date", ""),
            }
            for c in _get("commits", not_found=f"branch {branch!r} does not exist on GitHub",
                          sha=branch, per_page=COMMIT_LIMIT)
        ]
    return _cached(f"commits:{branch}", produce)


def refs(branch: str = "") -> dict:
    """Everything the pipeline-version section needs, in one call.

    Never raises: the fields fall back to free text, so a GitHub outage must
    not stop anyone launching a build.
    """
    result = {"repo": REPO, "branches": [], "commits": [], "branch": branch, "error": None}
    try:
        result["branches"] = branches()
        if branch:
            result["commits"] = commits(branch)
    except GitHubError as exc:
        result["error"] = str(exc)
    return result
