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


def _url(*parts: str) -> str:
    return "/".join([BASE, JOB_PATH, *[p.strip("/") for p in parts if p]])


def _posting_session(auth: tuple[str, str] | None) -> requests.Session:
    """A fresh session carrying a CSRF crumb issued for these credentials.

    Jenkins ties the crumb to the authenticated identity and to the session
    cookie handed out with it, so a crumb fetched anonymously (or as another
    user) will not validate the POST. Both must come from the same session as
    the request that uses them, which is why this is not cached.
    """
    s = requests.Session()
    try:
        r = s.get(f"{BASE}/crumbIssuer/api/json", auth=auth, timeout=TIMEOUT)
        if r.ok:
            body = r.json()
            s.headers[body["crumbRequestField"]] = body["crumb"]
    except (requests.RequestException, ValueError, KeyError):
        pass  # crumb protection may be off; let the POST itself report failure
    return s


# The parameter-separator plugin declares section headings as parameter
# definitions. They hold no value and can never cause a rejected build, so
# counting them as drift is a permanent false warning.
SEPARATOR_MARKER = "ParameterSeparator"


def real_parameters(properties: list[dict]) -> list[str]:
    """Names of parameters that actually take a value, separators excluded."""
    names = []
    for prop in properties or []:
        for definition in prop.get("parameterDefinitions") or []:
            kind = f"{definition.get('_class', '')} {definition.get('type', '')}"
            if SEPARATOR_MARKER in kind:
                continue
            names.append(definition["name"])
    return names


def declared_params() -> list[str]:
    """Parameter names the job declares. Used to catch drift between the
    pipeline's parameters{} block and params.yml."""
    r = _session.get(_url("api/json"), auth=AUTH, timeout=TIMEOUT,
                     params={"tree": "property[parameterDefinitions[name,type]]"})
    if not r.ok:
        raise JenkinsError(f"GET job: {r.status_code} {r.text[:200]}")
    return real_parameters(r.json().get("property", []))


def trigger(params: dict, auth: tuple[str, str] | None = None) -> str:
    """Queue a build as `auth`. Returns the queue item URL.

    Jenkins does not hand back a build number here -- the build has not been
    assigned an executor yet -- so the queue URL is what we persist and resolve
    later via queue_build_number().
    """
    auth = auth or AUTH
    form = {k: ("true" if v else "false") if isinstance(v, bool) else str(v) for k, v in params.items()}
    session = _posting_session(auth)
    r = session.post(_url("buildWithParameters"), data=form, auth=auth,
                     timeout=TIMEOUT, allow_redirects=False)
    if r.status_code in (401, 403):
        raise JenkinsError(
            f"Jenkins refused the trigger ({r.status_code}). The account needs Build permission "
            f"on {JOB_PATH}."
        )
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


def stop(number: int, auth: tuple[str, str] | None = None) -> None:
    auth = auth or AUTH
    r = _posting_session(auth).post(_url(str(number), "stop"), auth=auth, timeout=TIMEOUT)
    if r.status_code in (401, 403):
        raise JenkinsError(f"Jenkins refused the stop ({r.status_code}); the account needs "
                           f"Cancel permission on {JOB_PATH}.")
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
