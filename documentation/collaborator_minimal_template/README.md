# Minimal sushi metadata: nori → `filtered_counts` (no control barcodes)

This template is the smallest set of metadata that lets a custom-PRISM collaborator run **only**
the two sushi steps needed to reach `filtered_counts.csv`:

1. **collate_counts** — nori output (`raw_counts_uncollapsed.csv.gz`) → `prism_barcode_counts.csv`
2. **filter_counts** — `prism_barcode_counts.csv` → `filtered_counts.csv`

No CB normalization, no QC, no LFCs — so **no control barcodes** are used. The config keys people
usually worry about (`sig_cols`, `control_cols`, `cell_line_cols`, `count_col_name`, `ctl_types`)
are consumed only in `compute_l2fc` and later; they are **irrelevant** here. Keep the defaults for
the only three that matter: `ID_COLS=pcr_plate,pcr_well`,
`SEQUENCING_INDEX_COLS=flowcell_names,index_1,index_2`, `BARCODE_COL=forward_read_barcode`.

## Files in this template

| File | Purpose |
|---|---|
| `sample_meta.csv` | the 7 required sample-metadata columns (below) |
| `cell_line_meta.csv` | barcode → cell-line lookup |
| `cell_set_and_pool_meta.csv` | which cell lines each `cell_set` expects |
| `CB_meta.csv` | **optional** header-only stub — the code now tolerates a missing CB_meta (see below) |

## The 7 required `sample_meta.csv` columns

| Column | Why it is required |
|---|---|
| `pcr_plate` | `id_cols`: existence + uniqueness enforced in both steps; join key onto counts |
| `pcr_well` | other half of `id_cols` |
| `index_1` | sequencing index; maps reads → wells; must also exist in the nori file |
| `index_2` | sequencing index; reverse-complemented if `reverse_index2=TRUE` |
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

- `cell_line_meta.csv`: `depmap_id`, `lua`, `forward_read_barcode`
- `cell_set_and_pool_meta.csv`: `cell_set`, `depmap_id`, `lua` (`pool_id` optional)
- `CB_meta.csv`: `forward_read_barcode` + `cb_log2_dose` (or `cb_log10_dose`) — **rows may be empty**

## Control barcodes: `CB_meta.csv` is now optional

Both steps read a **hardcoded** `$BUILD_DIR/CB_meta.csv` (collate_counts.sh:111,
filter_counts.sh:79), ignoring the `CONTROL_BARCODE_META` config value. Historically an absent
file errored at `fread`.

Both reads now go through `read_cb_meta()` (`scripts/utils/kitchen_utensils.R`), which returns an
empty CB_meta (with `forward_read_barcode` + `cb_log2_dose`) when the file is missing — so a
custom-PRISM run with no control barcodes just works:

- `collate_counts.R:49` — `read_cb_meta(...)`; empty barcode set is merely additive to `known_barcodes`
- `filter_counts.R:60` — `read_cb_meta(...)`

The empty CB_meta carries `cb_log2_dose`, so the dose-column check at
`filter_counts_functions.R:56-66` passes and the downstream CB joins become no-ops — no changes
were needed in `filter_counts_functions.R`. The `CB_meta.csv` in this template is kept only as an
explicit, self-documenting stub; you can delete it and the run still succeeds.

## Verify a run

1. Put the four files + a small `raw_counts_uncollapsed.csv.gz` (columns
   `flowcell_name,flowcell_lane,index_1,index_2,forward_read_barcode,n`) in a `BUILD_DIR`.
2. Run collate → expect `prism_barcode_counts.csv`, non-empty, no CB_meta error.
3. Run filter → expect `filtered_counts.csv`; the CB template block is skipped (no "control
   barcode" template rows in the log).
4. `filtered_counts.csv` has one row per (well × expected cell line) with `n` populated, and the
   cell-line-purity line prints without error.
