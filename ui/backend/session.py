"""Identity via Jenkins credentials.

There is no user table here. The Jenkins instance decides who someone is, and
the same credentials that prove it are the ones that authorize the build
trigger -- which this Jenkins requires, since anonymous has read access but not
Job/Build. That kills two birds: nobody has to type their name into a box (and
be trusted to type it honestly), and launching works at all.

Sessions are held in memory, so a restart signs everyone out. That is
deliberate: the alternative is writing Jenkins credentials to disk. It also
means the service must run as a single worker -- with several uvicorn workers a
request could land on a process that has never seen the session.
"""

import os
import secrets
import time

import requests

import jenkins

# Long enough to cover a working day without a second login.
TTL = float(os.environ.get("SUSHI_UI_SESSION_TTL", 12 * 3600))
COOKIE = "sushi_sid"

_sessions: dict[str, dict] = {}


class AuthError(RuntimeError):
    pass


def verify(user: str, password: str) -> dict:
    """Confirm credentials against Jenkins and return the display identity.

    /me/api/json answers 403 for anonymous and 401 for bad credentials, so a
    200 here means genuinely authenticated -- anonymous read access cannot be
    mistaken for a successful login.
    """
    if not user.strip() or not password.strip():
        raise AuthError("Both a Jenkins user and a password are required.")
    try:
        r = requests.get(f"{jenkins.BASE}/me/api/json", auth=(user, password),
                         timeout=jenkins.TIMEOUT)
    except requests.RequestException as exc:
        raise AuthError(f"Could not reach Jenkins at {jenkins.BASE}: {exc}") from exc

    if r.status_code in (401, 403):
        raise AuthError("Jenkins rejected those credentials. If your password is correct, this "
                        "Jenkins may only accept API tokens over the API -- generate one at "
                        "/user/<you>/configure and use that instead.")
    if not r.ok:
        raise AuthError(f"Jenkins returned {r.status_code} for the credential check.")

    try:
        me = r.json()
    except ValueError as exc:
        raise AuthError("Jenkins did not return JSON for the credential check.") from exc

    identity = me.get("id") or user
    if identity == "anonymous":
        raise AuthError("Those credentials resolve to the anonymous user.")
    return {"id": identity, "full_name": me.get("fullName") or identity}


def create(user: str, password: str) -> tuple[str, dict]:
    """Verify, then remember the credentials against an opaque session id."""
    identity = verify(user, password)
    _evict_expired()
    sid = secrets.token_urlsafe(32)
    _sessions[sid] = {**identity, "user": user, "password": password, "created": time.time()}
    return sid, identity


def get(sid: str | None) -> dict | None:
    entry = _sessions.get(sid or "")
    if not entry:
        return None
    if time.time() - entry["created"] > TTL:
        _sessions.pop(sid, None)
        return None
    return entry


def identity(sid: str | None) -> dict | None:
    entry = get(sid)
    return {"id": entry["id"], "full_name": entry["full_name"]} if entry else None


def auth_for(sid: str | None) -> tuple[str, str] | None:
    """The Jenkins credentials to act as, for a POST that needs permission."""
    entry = get(sid)
    return (entry["user"], entry["password"]) if entry else None


def destroy(sid: str | None) -> None:
    _sessions.pop(sid or "", None)


def _evict_expired() -> None:
    now = time.time()
    for sid in [s for s, e in _sessions.items() if now - e["created"] > TTL]:
        _sessions.pop(sid, None)
