"""Parameter catalog and screen-type presets, loaded from the yml files."""

import json
import os
import re
import time
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
from functools import lru_cache
from pathlib import Path, PurePosixPath

import yaml

HERE = Path(__file__).parent
PRISMSEQ_ROOT = os.environ.get("PRISMSEQ_ROOT", "/cmap/obelix/pod/prismSeq")


SHIPPED_PRESETS = HERE / "presets.yml"

# Screen-type defaults are editable from the UI, so the live copy cannot be the
# file in git. Once this store exists it is authoritative; until then the
# shipped file seeds it. Put it somewhere backed up, like the database.
PRESETS_STORE = Path(os.environ.get("SUSHI_UI_PRESETS", HERE / "presets.local.yml"))


@lru_cache(maxsize=1)
def _params_file() -> dict:
    """params.yml, which is not editable from the UI."""
    return yaml.safe_load((HERE / "params.yml").read_text())


def presets_source() -> str:
    return "store" if PRESETS_STORE.exists() else "shipped"


def _load_presets() -> list[dict]:
    path = PRESETS_STORE if PRESETS_STORE.exists() else SHIPPED_PRESETS
    data = yaml.safe_load(path.read_text()) or {}
    return data.get("presets") or []


@lru_cache(maxsize=1)
def catalog() -> dict:
    params = _params_file()
    presets = _load_presets()

    known = {p["name"] for p in params["params"]}
    for preset in presets:
        unknown = set(preset["params"]) - known
        if unknown:
            raise ValueError(f"preset {preset['id']} sets unknown params: {sorted(unknown)}")

    group_ids = {g["id"] for g in params["groups"]}
    orphans = {p["name"] for p in params["params"] if p["group"] not in group_ids}
    if orphans:
        raise ValueError(f"params in undeclared groups: {sorted(orphans)}")

    return {"groups": params["groups"], "params": params["params"], "presets": presets}


def validate_presets(presets: list[dict]) -> list[str]:
    """Problems with an incoming screen-type list. Empty means safe to save.

    Checked before writing: a bad store would otherwise make catalog() raise on
    every request, taking the whole UI down until someone edited a file by hand.
    """
    problems = []
    if not isinstance(presets, list) or not presets:
        return ["At least one screen type is required."]

    index = {p["name"]: p for p in _params_file()["params"]}
    seen = set()
    for i, preset in enumerate(presets):
        where = preset.get("id") or f"entry {i + 1}"
        for field in ("id", "label", "dir"):
            if not str(preset.get(field, "")).strip():
                problems.append(f"{where}: {field} is required.")

        pid = str(preset.get("id", "")).strip()
        if pid:
            if pid in seen:
                problems.append(f"{pid}: duplicate screen type id.")
            seen.add(pid)
            if not re.fullmatch(r"[A-Za-z0-9_-]+", pid):
                problems.append(f"{pid}: id may only contain letters, digits, _ and -.")

        # dir is joined onto PRISMSEQ_ROOT, so it must not climb out of it.
        directory = str(preset.get("dir", "")).strip()
        if directory and (directory.startswith("/") or ".." in Path(directory).parts):
            problems.append(f"{where}: dir must be a relative path inside the prismSeq root.")

        params = preset.get("params") or {}
        if not isinstance(params, dict):
            problems.append(f"{where}: params must be a mapping of parameter to value.")
            continue
        for name, value in params.items():
            spec = index.get(name)
            if spec is None:
                problems.append(f"{where}: {name} is not a parameter this pipeline declares.")
                continue
            if spec["type"] == "choice" and str(value) not in [str(c) for c in spec["choices"]]:
                problems.append(f"{where}: {name} must be one of {spec['choices']}, got {value!r}.")
            if spec["type"] == "bool" and not isinstance(value, bool):
                problems.append(f"{where}: {name} must be true or false, got {value!r}.")

    return problems


def save_presets(presets: list[dict]) -> None:
    problems = validate_presets(presets)
    if problems:
        raise ValueError("; ".join(problems))
    PRESETS_STORE.parent.mkdir(parents=True, exist_ok=True)
    # Write-then-rename so a failure cannot leave a half-written store that
    # breaks catalog() for everyone.
    tmp = PRESETS_STORE.with_suffix(".tmp")
    tmp.write_text(yaml.safe_dump({"presets": presets}, sort_keys=False, allow_unicode=True))
    tmp.replace(PRESETS_STORE)
    catalog.cache_clear()


def reset_presets() -> None:
    """Drop the store and fall back to the shipped presets.yml."""
    PRESETS_STORE.unlink(missing_ok=True)
    catalog.cache_clear()


def param_index() -> dict:
    return {p["name"]: p for p in catalog()["params"]}


def preset(preset_id: str) -> dict:
    for p in catalog()["presets"]:
        if p["id"] == preset_id:
            return p
    raise KeyError(preset_id)


def build_dir_for(preset_id: str, build_name: str) -> str:
    """Where a build of this type lives. Empty build name -> the type's parent dir."""
    return os.path.join(PRISMSEQ_ROOT, preset(preset_id)["dir"], build_name or "")


def defaults_for(preset_id: str, build_name: str = "", build_dir: str = "") -> dict:
    """Catalog defaults, overlaid with the preset, overlaid with name-derived values."""
    values = {p["name"]: p["default"] for p in catalog()["params"]}
    values.update(preset(preset_id)["params"])
    if build_name:
        values["BUILD_NAME"] = build_name
        values["SCREEN"] = build_name
        values["BUILD_DIR"] = build_dir or build_dir_for(preset_id, build_name)
    elif build_dir:
        values["BUILD_DIR"] = build_dir
    return values


# ------------------------------------------------------- build path discovery

# A directory is a build if it holds either of these. config.json is written by
# make_config_file.groovy; the counts file is what a collate run leaves behind,
# so a build part-way through its first run still registers.
BUILD_MARKERS = ("config.json", "raw_counts_uncollapsed.csv.gz")

# The scan is ~600 directory reads on an NFS mount, so it is cached. Typing a
# path that is not in the list still works, so a stale list costs suggestions,
# never the ability to launch.
SCAN_TTL = float(os.environ.get("SUSHI_UI_SCAN_TTL", "600"))
_scan: tuple[float, list[dict]] | None = None


def resolve_rel(rel_path: str) -> str:
    """Normalize a client-supplied build path; refuse anything outside the root.

    This value becomes BUILD_DIR, which drives file writes (archiving an
    existing config.json), so a path escaping PRISMSEQ_ROOT must not survive
    this function. normpath rather than resolve: '..' is collapsed textually,
    without following symlinks that may legitimately point off the mount.
    """
    rel = (rel_path or "").strip().strip("/")
    if not rel:
        return ""
    root = Path(PRISMSEQ_ROOT)
    target = Path(os.path.normpath(str(root / rel)))
    if target != root and root not in target.parents:
        raise ValueError(f"build path must stay under {PRISMSEQ_ROOT}: {rel_path!r}")
    return target.relative_to(root).as_posix()


def _scan_dir(path: Path) -> tuple[set[str], list[Path]]:
    """One pass over a directory: child names, and which of them are dirs.

    Unreadable directories exist on this mount; they are skipped quietly.
    """
    names: set[str] = set()
    subdirs: list[Path] = []
    try:
        for entry in os.scandir(path):
            names.add(entry.name)
            if not entry.name.startswith(".") and entry.is_dir():
                subdirs.append(Path(entry.path))
    except OSError:
        pass
    return names, subdirs


def _is_build(names: set[str]) -> bool:
    return any(m in names for m in BUILD_MARKERS)


def _iso_mtime(path: Path) -> str:
    try:
        mtime = path.stat().st_mtime
    except OSError:
        mtime = 0.0
    return datetime.fromtimestamp(mtime, timezone.utc).isoformat(timespec="seconds")


def _build_entry(root: Path, path: Path, names: set[str]) -> dict:
    parent = path.parent
    return {
        "path": path.relative_to(root).as_posix(),
        "name": path.name,
        "parent": parent.relative_to(root).as_posix() if parent != root else "",
        "has_config": "config.json" in names,
        "has_raw_counts": "raw_counts_uncollapsed.csv.gz" in names,
        "modified": _iso_mtime(path),
    }


def build_paths(refresh: bool = False) -> list[dict]:
    """Every directory under PRISMSEQ_ROOT that looks like a build directory.

    Builds sit either one level down (PDEV021/) or two (MTS_SEQ/MTS032/), so
    the walk stops at depth 2. A build never contains another build, which is
    why a depth-1 hit short-circuits its own subtree.
    """
    global _scan
    if not refresh and _scan and (time.monotonic() - _scan[0]) < SCAN_TTL:
        return _scan[1]

    root = Path(PRISMSEQ_ROOT)
    _, tops = _scan_dir(root)

    # NFS round-trip latency dominates here, not CPU, so the reads are fanned
    # out: serially this takes ~12s, which is too slow to block a page load.
    with ThreadPoolExecutor(max_workers=16) as pool:
        top_scans = dict(zip(tops, pool.map(_scan_dir, tops)))
        nested = [sub for top in tops
                  if not _is_build(top_scans[top][0])
                  for sub in top_scans[top][1]]
        sub_scans = dict(zip(nested, pool.map(_scan_dir, nested)))

        hits = [(top, top_scans[top][0]) for top in tops if _is_build(top_scans[top][0])]
        hits += [(sub, sub_scans[sub][0]) for sub in nested if _is_build(sub_scans[sub][0])]
        found = list(pool.map(lambda h: _build_entry(root, h[0], h[1]), hits))

    # Most recently touched first: the build someone is about to launch is
    # almost always the one they just created or last worked in.
    found.sort(key=lambda b: b["modified"], reverse=True)
    _scan = (time.monotonic(), found)
    return found


def parent_of(rel_path: str) -> str | None:
    """The level above a relative path; None when already at the root."""
    if not rel_path:
        return None
    parent = PurePosixPath(rel_path).parent
    return "" if str(parent) == "." else str(parent)


def list_dir(rel_path: str = "") -> dict:
    """One level of the build tree, for the directory browser.

    Each child is classified as a build directory or a plain folder to descend
    into, so the browser can offer 'select' or 'open' without a second call.
    """
    rel = resolve_rel(rel_path)
    root = Path(PRISMSEQ_ROOT)
    target = root / rel if rel else root
    names, subdirs = _scan_dir(target)

    # One directory read per child, fanned out: a screen-type folder holds ~80
    # builds and the root holds ~180 entries.
    with ThreadPoolExecutor(max_workers=32) as pool:
        scans = list(pool.map(_scan_dir, subdirs))
        mtimes = list(pool.map(_iso_mtime, subdirs))

    entries = [
        {
            "name": sub.name,
            "path": sub.relative_to(root).as_posix(),
            "is_build": _is_build(child_names),
            "has_config": "config.json" in child_names,
            "has_raw_counts": "raw_counts_uncollapsed.csv.gz" in child_names,
            "child_dirs": len(child_subdirs),
            "modified": mtime,
        }
        for sub, (child_names, child_subdirs), mtime in zip(subdirs, scans, mtimes)
    ]
    entries.sort(key=lambda e: e["name"].lower())

    return {
        "root": PRISMSEQ_ROOT,
        "path": rel,
        "parent": parent_of(rel),
        "exists": target.is_dir(),
        "is_build": _is_build(names),
        "entries": entries,
    }


def _declared_screen_type(build_dir: Path) -> str:
    try:
        return str(json.loads((build_dir / "config.json").read_text()).get("SCREEN_TYPE", ""))
    except (OSError, json.JSONDecodeError, AttributeError):
        return ""


def infer_preset(rel_path: str) -> tuple[str, str]:
    """Best guess at the screen type for a build path, plus why.

    Directory before config.json deliberately: presets.yml notes that two
    reference builds carry a SCREEN_TYPE that disagrees with the directory they
    live in, and the directory is what actually locates the data.
    """
    presets = catalog()["presets"]
    parts = PurePosixPath(rel_path).parts if rel_path else ()

    if len(parts) > 1:
        for p in presets:
            if p["dir"] == parts[0]:
                return p["id"], f"directory {parts[0]}/"

    if rel_path:
        declared = _declared_screen_type(Path(PRISMSEQ_ROOT) / rel_path)
        for p in presets:
            if declared and p["id"] == declared:
                return p["id"], f"SCREEN_TYPE={declared} in config.json"

    if parts:
        prefix = parts[-1][:3].upper()
        for p in presets:
            if prefix and p["id"].startswith(prefix):
                return p["id"], f"build name starts with {prefix}"

    return presets[0]["id"], "nothing in the path identifies a screen type"


def defaults_for_path(rel_path: str, preset_override: str = "") -> dict:
    """Resolve a build path to a screen type and a full set of form values."""
    rel = resolve_rel(rel_path)
    if preset_override:
        preset(preset_override)  # KeyError -> unknown screen type
        preset_id, reason = preset_override, "chosen manually"
    else:
        preset_id, reason = infer_preset(rel)

    build_dir = os.path.join(PRISMSEQ_ROOT, rel) if rel else PRISMSEQ_ROOT + "/"
    target = Path(build_dir)
    return {
        "path": rel,
        "preset": preset_id,
        "preset_reason": reason,
        "inferred": not preset_override,
        "values": defaults_for(preset_id, PurePosixPath(rel).name if rel else "", build_dir),
        "from_preset": sorted(preset(preset_id)["params"]),
        "build_dir": build_dir,
        "exists": target.is_dir(),
        "has_config": (target / "config.json").is_file(),
        "has_raw_counts": (target / "raw_counts_uncollapsed.csv.gz").is_file(),
    }


def coerce(values: dict) -> dict:
    """Drop unknown keys and normalize types so what we store is what Jenkins gets.

    Jenkins rejects buildWithParameters calls containing parameters the job does
    not declare, so silently dropping unknown keys here is the safe direction.
    """
    index = param_index()
    out = {}
    for name, value in values.items():
        spec = index.get(name)
        if spec is None:
            continue
        if spec["type"] == "bool":
            out[name] = value if isinstance(value, bool) else str(value).lower() in ("true", "1", "yes", "on")
        else:
            out[name] = "" if value is None else str(value)
    return out


def validate(values: dict) -> list[str]:
    """Return human-readable problems. Empty list means good to launch."""
    index = param_index()
    problems = []

    for name, spec in index.items():
        if spec.get("required") and not str(values.get(name, "")).strip():
            problems.append(f"{spec['label']} ({name}) is required.")
        if spec["type"] == "choice" and name in values:
            if str(values[name]) not in [str(c) for c in spec["choices"]]:
                problems.append(f"{name} must be one of {spec['choices']}, got {values[name]!r}.")

    for name in ("COUNT_THRESHOLD", "CHUNK_SIZE", "N_SAMPLES", "PSEUDOCOUNT", "VIABILITY_CAP"):
        raw = str(values.get(name, "")).strip()
        if raw:
            try:
                float(raw)
            except ValueError:
                problems.append(f"{name} must be a number, got {raw!r}.")

    # These are the dependencies the pipeline itself does not check; a run that
    # violates them fails deep inside a module with a confusing R error.
    if values.get("GENERATE_QC_TABLES_2") and not values.get("COMPUTE_LFC"):
        problems.append("QC tables (post-LFC) requires Compute log fold change.")
    if values.get("LFC_BIOMARKER") and not values.get("COMPUTE_LFC"):
        problems.append("Biomarker on LFC requires Compute log fold change.")
    if values.get("AUC_BIOMARKER") and not values.get("DRC"):
        problems.append("Biomarker on AUC requires Dose response curves.")
    if values.get("SYNERGY_BIOMARKER") and not values.get("SYNERGY"):
        problems.append("Biomarker on synergy requires Synergy.")
    if values.get("TEST_DATASET") and not values.get("CONVERT_SUSHI"):
        problems.append("Upload to test bucket has no effect unless Convert for MTS pipeline is on.")
    if not values.get("USE_LATEST") and not str(values.get("COMMIT_ID", "")).strip():
        problems.append("Pin to commit is required when 'Use latest commit on branch' is off.")

    return problems


MODULE_GROUPS = ("modules", "qc_modules", "analytics", "portal", "metadata")


def module_names() -> list[str]:
    """Params that decide which scripts run. Jenkins never writes these into
    config.json, so the run record is the only place they are preserved."""
    return [p["name"] for p in catalog()["params"] if p["group"] in MODULE_GROUPS]
