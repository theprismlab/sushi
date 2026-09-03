# Minimal sushi metadata: nori → `filtered_counts` (no control barcodes)

This template is the smallest set of metadata that lets a custom-PRISM collaborator run **only**
the two sushi steps needed to reach `filtered_counts.csv`:

1. **collate_counts** — nori output (`raw_counts_uncollapsed.csv.gz`) → `prism_barcode_counts.csv`
2. **filter_counts** — `prism_barcode_counts.csv` → `filtered_counts.csv`

No control barcode normalization, no QC, and no LFC computation — so **no control barcodes** are used. The config keys people
usually worry about (`sig_cols`, `control_cols`, `cell_line_cols`, `count_col_name`, `ctl_types`)
are consumed only in `compute_l2fc` and later; they are **irrelevant** here. Keep the defaults for
the only three that matter: `ID_COLS=pcr_plate,pcr_well`,
`SEQUENCING_INDEX_COLS=flowcell_names,flowcell_lanes,index_1,index_2`, `BARCODE_COL=forward_read_barcode`.

## Files in this template

| File | Purpose |
|---|---|
| `sample_meta.csv` | the 7 required sample-metadata columns (below) |
| `pool_ids.csv` | list of PRISM pools used in the experiment |

No `CB_meta.csv` is needed — the pipeline tolerates its absence (see below).

## Pipeline arguments to change from their defaults

These are the `make_config_file.groovy` parameters that must differ from their defaults for a
minimal, no-controls, `filtered_counts`-only run. Everything not listed here keeps its default.

### Why the column config keys are NOT in this list

`collate_counts.sh` and `filter_counts.sh` pass only three column keys to the R scripts —
`ID_COLS`, `SEQUENCING_INDEX_COLS`, `BARCODE_COL` (see their `args=(...)` arrays). The
"condition" column keys — `SIG_COLS`, `CONTROL_COLS`, `CELL_LINE_COLS`, `COUNT_COL_NAME`,
`CTL_TYPES` — are **never passed to either step**; they're read only by `compute_l2fc`/normalize/QC,
which are off. So they are inert here regardless of their value.

The three keys that *are* used keep their defaults **only because this template's
`sample_meta.csv` uses the canonical column names those defaults point at** (`pcr_plate`,
`pcr_well`, `index_1`, `index_2`, `flowcell_names`). The column requirement is real — it's just
expressed in the sample-meta column *names* (next section), not in a config change. If your data
uses different headers (e.g. `plate`/`well`), either rename to canonical **or** set the matching
key, e.g. `ID_COLS=plate,well`.

**Point the build at your data:**

| Argument | Default | Set to | Why |
|---|---|---|---|
| `BUILD_DIR` | `/cmap/obelix/pod/prismSeq/` | your build directory | Must contain the Nori `raw_counts_uncollapsed.csv.gz` **and** the three metadata files from this template |
| `BUILD_NAME` | `` (empty) | your build name | Names the output files; match the `BUILD_DIR` name |
| `GIT_BRANCH` | `main` | `develop` | The missing-`CB_meta` tolerance (`read_cb_meta`) lives on `develop`; use `main` only once it is merged there |

**Turn off every stage except collate + filter** (all default `true`, each pulls in control
barcodes, QC, or LFCs that this run does not do):

| Argument | Default | Set to | Why |
|---|---|---|---|
| `CREATE_CELLDB_METADATA` | `true` | `false` | You supply `cell_line_meta.csv` / `cell_set_and_pool_meta.csv`; don't fetch from cellDB |
| `CBNORMALIZE` | `true` | `false` | Normalization requires control barcodes + vehicle controls |
| `COMPUTE_LFC` | `true` | `false` | No fold changes |
| `BIAS_CORRECTION` | `true` | `false` | Operates on LFCs |
| `COLLAPSE` | `true` | `false` | Collapses LFC replicates |
| `GENERATE_WELL_METRICS` | `true` | `false` | Uses control barcodes |
| `GENERATE_QC_TABLES` | `true` | `false` | QC we don't run |
| `FILTER_SKIPPED_WELLS` | `true` | `false` | `skipped_wells.csv` comes from `create_sample_meta` (not run); a no-op without the file, but set `false` for clarity |

Keep on: `COLLATE_FASTQ_READS` and `FILTER_COUNTS` (both default `true`). `CREATE_SAMPLE_META`
is already `false` by default — leave it (you provide `sample_meta.csv`). `FILTER_QC_FLAGS` and
`FILTER_FAILED_LINES` only affect the QC/collapse stages you've turned off, so their values don't
matter. `CONTROL_BARCODE_META` is irrelevant here (only used by the disabled cellDB/QC stages).

## The 7 required `sample_meta.csv` columns

| Column | Why it is required |
|---|---|
| `pcr_plate` | `id_cols`: existence + uniqueness enforced in both steps; join key onto counts |
| `pcr_well` | other half of `id_cols` |
| `index_1` | sequencing index; maps reads → wells |
| `index_2` | sequencing index; maps reads → plates |
| `flowcell_names` | flowcell-filter key; also checked unconditionally at collate startup |
| `flowcell_lanes` | paired with `flowcell_names` to enumerate expected flowcell+lane combos |
| `cell_set` | filter_counts merge key to `cell_set_and_pool_meta`; empty/NA rows are dropped from the expected-reads template |

Value rules:
- `(flowcell_names, index_1, index_2)` must **uniquely identify** each well and map **1:1** to
  `(pcr_plate, pcr_well)`.
- every `cell_set` value must appear in `cell_set_and_pool_meta.csv`.

**Do not add `cb_ladder`.** Omitting it makes the control-barcode template block skip safely.
Any other sushi sample-meta column (`pert_name`, `pert_dose`, `day`, `bio_rep`, `x_project_id`, …)
is pure pass-through for these two steps and only matters if normalize/QC/LFC are run later.

## Companion files — required columns

- `pool_ids.csv` — `pool_id`

## Control barcodes: no `CB_meta.csv` required

Both steps read a **hardcoded** `$BUILD_DIR/CB_meta.csv` (collate_counts.sh:111,
filter_counts.sh:79), ignoring the `CONTROL_BARCODE_META` config value. Historically an absent
file errored at `fread`.

Both reads now go through `read_cb_meta()` (`scripts/utils/kitchen_utensils.R`), which returns an
empty CB_meta (with `forward_read_barcode` + `cb_log2_dose`) when the file is missing — so a
custom-PRISM run with no control barcodes just works, no `CB_meta.csv` needed:

- `collate_counts.R:49` — `read_cb_meta(...)`; empty barcode set is merely additive to `known_barcodes`
- `filter_counts.R:60` — `read_cb_meta(...)`

The empty CB_meta carries `cb_log2_dose`, so the dose-column check at
`filter_counts_functions.R:56-66` passes and the downstream CB joins become no-ops — no changes
were needed in `filter_counts_functions.R`.

## Verify a run

1. Put the three metadata files + a small `raw_counts_uncollapsed.csv.gz` (columns
   `flowcell_name,flowcell_lane,index_1,index_2,forward_read_barcode,n`) in a `BUILD_DIR`.
2. Run collate → expect `prism_barcode_counts.csv`, non-empty, no CB_meta error.
3. Run filter → expect `filtered_counts.csv`; the CB template block is skipped (no "control
   barcode" template rows in the log).
4. `filtered_counts.csv` has one row per (well × expected cell line) with `n` populated, and the
   cell-line-purity line prints without error.
