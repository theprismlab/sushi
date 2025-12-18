import polars as pl
from pathlib import Path
from typing import Callable
from sushilib.sushi_io.reader import read_polars


class SushiBuild:
    """
    Container for all of the CSV/JSON tables in a build directory.
    Usage:
        b = SushiBuild("/path/to/build")
        df = b.normalized_counts_filtered
        for name, df in b:
            pass
    """

    # manifest: (attribute_name, relative_path, needs_pert_dose_override)
    _FILES = [
        ("normalized_counts_filtered", "normalized_counts.csv"),
        ("normalized_counts_original", "normalized_counts_original.csv"),
        ("annotated_counts", "annotated_counts.csv"),
        ("filtered_counts_filtered", "filtered_counts.csv"),
        ("filtered_counts_original", "filtered_counts_original.csv"),
        ("l2fc_filtered", "l2fc.csv"),
        ("l2fc_original", "l2fc_original.csv"),
        ("collapsed_l2fc", "collapsed_l2fc.csv"),
        ("qc_table", "qc_tables/plate_cell_qc_table_internal.csv"),
        ("qc_flags", "qc_tables/plate_cell_qc_flags.csv"),
        ("id_cols_table", "qc_tables/id_cols_qc_table.csv"),
        ("id_cols_flags", "qc_tables/id_cols_qc_flags.csv"),
        ("pool_well_table", "qc_tables/pool_well_qc_table.csv"),
        ("sample_meta", "sample_meta.csv"),
        ("eps_qc_table", "qc_tables/eps_qc_table.csv"),
    ]


    def __init__(self, build_path):
        self.build_path = Path(build_path).expanduser().resolve()
        self._load_all()

    def _load_all(self):
        for attr, rel in self._FILES:
            p = self.build_path / rel
            if not p.exists():
                setattr(self, attr, None)
                continue
            df = read_polars(p)

            setattr(self, attr, df)

    def __repr__(self):
        # Show single build or combined builds
        if hasattr(self, "build_path"):
            path_repr = self.build_path
        elif hasattr(self, "build_paths"):
            path_repr = list(self.build_paths)
        else:
            path_repr = None
        loaded = [
            attr for attr, *_ in self._FILES if getattr(self, attr, None) is not None
        ]
        return f"<SushiBuild path={path_repr!r} tables={loaded}>"

    def __iter__(self):
        """Iterate over (attribute_name, DataFrame) pairs for all loaded tables."""
        for attr, rel in self._FILES:
            df = getattr(self, attr, None)
            if df is not None:
                yield attr, df

    def __getitem__(self, key):
        if isinstance(key, str):
            return getattr(self, key)
        if isinstance(key, int):
            attr = self._FILES[key][0]
            return getattr(self, attr)
        raise KeyError(f"Invalid key {key!r}; must be str or int.")

    def apply_to_tables(self, fn):
        """Apply fn(df) to every non-None table in this instance, replacing it."""
        for attr, rel in self._FILES:
            df = getattr(self, attr, None)
            if df is not None:
                setattr(self, attr, fn(df))
        return self

    @classmethod
    def concat_builds(cls, builds: list["SushiBuild"]) -> "SushiBuild":
        combined = object.__new__(cls)
        combined.build_paths = []
        for b in builds:
            if hasattr(b, "build_paths"):
                combined.build_paths.extend(b.build_paths)
            else:
                combined.build_paths.append(b.build_path)
        for attr, rel in cls._FILES:
            dfs = [
                getattr(b, attr, None)
                for b in builds
                if getattr(b, attr, None) is not None
            ]
            if len(dfs) >= 2:
                ref = dfs[0]
                ref_cols, ref_dtypes = ref.columns, ref.dtypes
                aligned = []
                for df in dfs:
                    for col, dtype in zip(ref_cols, ref_dtypes):
                        if col not in df.columns:
                            df = df.with_columns(pl.lit(None).cast(dtype).alias(col))
                    df = df.select(ref_cols)
                    df = df.with_columns(
                        [
                            pl.col(col).cast(dtype)
                            for col, dtype in zip(ref_cols, ref_dtypes)
                        ]
                    )
                    aligned.append(df)
                setattr(combined, attr, pl.concat(aligned, how="vertical"))
            else:
                setattr(combined, attr, None)
        return combined

    def __add__(self, other: object) -> "SushiBuild":
        if not isinstance(other, SushiBuild):
            return NotImplemented
        return self.concat_builds([self, other])

    def __radd__(self, other: object) -> "SushiBuild":
        if other == 0:
            return self
        if not isinstance(other, SushiBuild):
            return NotImplemented
        return self.concat_builds([other, self])

    def update_tables(self, fn: Callable[[pl.DataFrame], pl.DataFrame]) -> "SushiBuild":
        skipped = []
        for attr, rel in self._FILES:
            df = getattr(self, attr, None)
            if df is None:
                continue

            try:
                new_df = fn(df)
                n_rows_modified = abs(len(new_df) - len(df))
                print(f" ↳ {n_rows_modified} rows modified in {attr}")
                print(f" ↳ updating {attr} and saving as {rel}")

            except KeyError as e:
                skipped.append((attr, f"missing column {e}"))
                continue
            except Exception as e:
                skipped.append((attr, str(e)))
                continue

            # write it back
            out_path = (self.build_path / rel).resolve()
            out_path.parent.mkdir(parents=True, exist_ok=True)
            new_df.write_csv(out_path)
            setattr(self, attr, new_df)

        if skipped:
            for name, reason in skipped:
                print(f" ↳ skipped {name}: {reason}")
        return self
