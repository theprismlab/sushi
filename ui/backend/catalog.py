"""Parameter catalog and screen-type presets, loaded from the yml files."""

import os
from functools import lru_cache
from pathlib import Path

import yaml

HERE = Path(__file__).parent
PRISMSEQ_ROOT = os.environ.get("PRISMSEQ_ROOT", "/cmap/obelix/pod/prismSeq")


@lru_cache(maxsize=1)
def catalog() -> dict:
    params = yaml.safe_load((HERE / "params.yml").read_text())
    presets = yaml.safe_load((HERE / "presets.yml").read_text())["presets"]

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


def defaults_for(preset_id: str, build_name: str = "") -> dict:
    """Catalog defaults, overlaid with the preset, overlaid with name-derived values."""
    values = {p["name"]: p["default"] for p in catalog()["params"]}
    values.update(preset(preset_id)["params"])
    if build_name:
        values["BUILD_NAME"] = build_name
        values["SCREEN"] = build_name
        values["BUILD_DIR"] = build_dir_for(preset_id, build_name)
    return values


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
