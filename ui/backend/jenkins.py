"""Thin Jenkins REST client.

Jenkins stays the executor: it does the git checkout, writes config.json via
make_config_file.groovy, and runs the modules in podman. This module only
triggers builds and reads status/logs back.
"""

import os

import requests

BASE = os.environ.get("JENKINS_URL", "http://localhost:8080").rstrip("/")
# URL path of the job, so nested folders work: e.g. "job/prism/job/sushi".
JOB_PATH = os.environ.get("JENKINS_JOB_PATH", "job/sushi").strip("/")
TIMEOUT = float(os.environ.get("JENKINS_TIMEOUT", "30"))

_user = os.environ.get("JENKINS_USER")
_token = os.environ.get("JENKINS_TOKEN")
AUTH = (_user, _token) if _user and _token else None


class JenkinsError(RuntimeError):
    pass


_session = requests.Session()
_crumb: dict | None = None


def _url(*parts: str) -> str:
    return "/".join([BASE, JOB_PATH, *[p.strip("/") for p in parts if p]])


def _crumb_headers() -> dict:
    """Jenkins requires a CSRF crumb for POST unless protection is disabled."""
    global _crumb
    if _crumb is None:
        try:
            r = _session.get(f"{BASE}/crumbIssuer/api/json", auth=AUTH, timeout=TIMEOUT)
            _crumb = {r.json()["crumbRequestField"]: r.json()["crumb"]} if r.ok else {}
        except requests.RequestException:
            _crumb = {}
    return dict(_crumb)


def declared_params() -> list[str]:
    """Parameter names the job declares. Used to catch drift between the
    pipeline's parameters{} block and params.yml."""
    r = _session.get(_url("api/json"), auth=AUTH, timeout=TIMEOUT,
                     params={"tree": "property[parameterDefinitions[name]]"})
    if not r.ok:
        raise JenkinsError(f"GET job: {r.status_code} {r.text[:200]}")
    names = []
    for prop in r.json().get("property", []):
        for definition in prop.get("parameterDefinitions", []) or []:
            names.append(definition["name"])
    return names


def trigger(params: dict) -> str:
    """Queue a build. Returns the queue item URL.

    Jenkins does not hand back a build number here -- the build has not been
    assigned an executor yet -- so the queue URL is what we persist and resolve
    later via queue_build_number().
    """
    form = {k: ("true" if v else "false") if isinstance(v, bool) else str(v) for k, v in params.items()}
    r = _session.post(_url("buildWithParameters"), data=form, auth=AUTH,
                      headers=_crumb_headers(), timeout=TIMEOUT, allow_redirects=False)
    if r.status_code not in (200, 201, 302, 303):
        raise JenkinsError(f"trigger failed: {r.status_code} {r.text[:500]}")
    queue_url = r.headers.get("Location")
    if not queue_url:
        raise JenkinsError("trigger succeeded but Jenkins returned no queue Location header")
    return queue_url.rstrip("/")


def queue_build_number(queue_url: str) -> tuple[int | None, str | None]:
    """Resolve a queue item to a build number.

    Returns (number, cancel_reason). Both None means still waiting in the queue.
    """
    r = _session.get(f"{queue_url}/api/json", auth=AUTH, timeout=TIMEOUT)
    if r.status_code == 404:
        # Jenkins drops queue items ~5 min after they start; without a build
        # number by then we have permanently lost track of this build.
        return None, "queue item expired before a build number was recorded"
    if not r.ok:
        raise JenkinsError(f"queue lookup: {r.status_code} {r.text[:200]}")
    item = r.json()
    if item.get("cancelled"):
        return None, "cancelled in queue"
    executable = item.get("executable")
    return (executable.get("number") if executable else None), None


def build(number: int) -> dict:
    r = _session.get(_url(str(number), "api/json"), auth=AUTH, timeout=TIMEOUT,
                     params={"tree": "building,result,timestamp,duration,url,fullDisplayName"})
    if r.status_code == 404:
        raise JenkinsError(f"build {number} not found (rotated out of Jenkins history?)")
    if not r.ok:
        raise JenkinsError(f"build lookup: {r.status_code} {r.text[:200]}")
    return r.json()


def stop(number: int) -> None:
    r = _session.post(_url(str(number), "stop"), auth=AUTH, headers=_crumb_headers(), timeout=TIMEOUT)
    if not r.ok and r.status_code not in (302, 303):
        raise JenkinsError(f"stop failed: {r.status_code} {r.text[:200]}")


def console(number: int, start: int = 0) -> dict:
    """Incremental console text, so the UI can tail a running build."""
    r = _session.get(_url(str(number), "logText/progressiveText"), auth=AUTH,
                     timeout=TIMEOUT, params={"start": start})
    if not r.ok:
        raise JenkinsError(f"console fetch: {r.status_code} {r.text[:200]}")
    return {
        "text": r.text,
        "start": start,
        "next": int(r.headers.get("X-Text-Size", start + len(r.text))),
        "more": r.headers.get("X-More-Data", "false").lower() == "true",
    }


def build_url(number: int) -> str:
    return _url(str(number)) + "/"
