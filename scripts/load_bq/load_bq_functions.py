# load_bq_functions.py
import logging
from google.cloud import bigquery
import os
import tempfile
import polars as pl
from typing import Dict, Any
import farmhash
import ctypes

import sys
import os

#TODO: Update how sushilib is imported once we have proper packaging
this_file_path = os.path.abspath(__file__)
project_root = os.path.dirname(os.path.dirname(os.path.dirname(this_file_path)))
if project_root not in sys.path:
    sys.path.append(project_root)
from sushilib.sushi_io.reader import read_polars

def get_polars_dtype_from_bq_field(field: bigquery.SchemaField) -> pl.DataType:
    """Map BigQuery field types to Polars data types."""
    type_mapping = {
        "STRING": pl.Utf8,
        "INTEGER": pl.Int64,
        "INT64": pl.Int64,
        "FLOAT": pl.Float64,
        "FLOAT64": pl.Float64,
        "BOOLEAN": pl.Boolean,
        "BOOL": pl.Boolean,
        "DATE": pl.Date,
        "DATETIME": pl.Datetime,
        "TIMESTAMP": pl.Datetime,
        "TIME": pl.Time,
        "NUMERIC": pl.Float64,
        "BIGNUMERIC": pl.Float64,
        "BYTES": pl.Binary,
        "JSON": pl.Utf8,  # JSON stored as string in Polars
    }

    bq_type = field.field_type.upper()
    polars_type = type_mapping.get(bq_type, pl.Utf8)  # Default to string if unknown

    if field.mode == "REPEATED":
        # For repeated fields, wrap in List
        return pl.List(polars_type)

    return polars_type


def coerce_dataframe_types(df: pl.DataFrame, table: bigquery.Table) -> pl.DataFrame:
    """Coerce DataFrame column types to match BigQuery table schema."""
    schema_map = {field.name: field for field in table.schema}

    cast_expressions = []

    for col in df.columns:
        if col in schema_map:
            field = schema_map[col]
            target_dtype = get_polars_dtype_from_bq_field(field)
            current_dtype = df[col].dtype

            # Skip if types already match
            if current_dtype == target_dtype:
                cast_expressions.append(pl.col(col))
                continue

            try:
                logging.info(
                    f"Coercing column '{col}' from {current_dtype} to {target_dtype}"
                )

                # Handle special cases for type coercion
                if target_dtype == pl.Boolean:
                    # Handle various boolean representations
                    cast_expr = (
                        pl.col(col)
                        .map_elements(
                            lambda x: _coerce_to_boolean(x), return_dtype=pl.Boolean
                        )
                        .alias(col)
                    )
                elif target_dtype in [pl.Int64, pl.Float64]:
                    # Handle numeric coercion with null handling
                    cast_expr = (
                        pl.col(col)
                        .map_elements(
                            lambda x: _coerce_to_numeric(x, target_dtype),
                            return_dtype=target_dtype,
                        )
                        .alias(col)
                    )
                elif target_dtype == pl.Date:
                    # Handle date parsing
                    cast_expr = (
                        pl.col(col)
                        .str.strptime(pl.Date, format=None, strict=False)
                        .alias(col)
                    )
                elif target_dtype == pl.Datetime:
                    # Handle datetime parsing
                    cast_expr = (
                        pl.col(col)
                        .str.strptime(pl.Datetime, format=None, strict=False)
                        .alias(col)
                    )
                else:
                    # Default casting
                    cast_expr = pl.col(col).cast(target_dtype, strict=False).alias(col)

                cast_expressions.append(cast_expr)

            except Exception as e:
                logging.warning(
                    f"Failed to coerce column '{col}' to {target_dtype}: {e}"
                )
                # Keep original column if coercion fails
                cast_expressions.append(pl.col(col))
        else:
            # Keep columns not in schema as-is
            cast_expressions.append(pl.col(col))

    # Apply all transformations at once
    df = df.with_columns(cast_expressions)
    return df


def _coerce_to_boolean(value) -> bool:
    """Helper function to coerce various values to boolean."""
    if value is None or (
        isinstance(value, str) and value.lower() in ["", "na", "null"]
    ):
        return None
    if isinstance(value, bool):
        return value
    if isinstance(value, str):
        return value.lower() in ["true", "1", "yes", "t", "y"]
    if isinstance(value, (int, float)):
        return bool(value)
    return None


def _coerce_to_numeric(value, target_dtype):
    """Helper function to coerce values to numeric types."""
    if value is None or (
        isinstance(value, str) and value.lower() in ["", "na", "null"]
    ):
        return None
    try:
        if target_dtype == pl.Int64:
            return int(float(value))  # Convert through float to handle "1.0" -> 1
        elif target_dtype == pl.Float64:
            return float(value)
    except (ValueError, TypeError):
        return None
    return None


def filter_csv_to_matching_columns(
    file_path: str, table: bigquery.Table, build_name: str, screen: str
) -> str:
    """Creates a temp CSV with columns correctly aligned to the BigQuery table schema."""
    schema_map = {field.name: field for field in table.schema}

    # Read CSV with more robust null handling
    null_values = ["", "NA", "NULL", "null", "N/A", "n/a", "NaN", "nan"]

    try:
        df = read_polars(
            file_path,
        )
        logging.info(
            f"Successfully read CSV with {len(df)} rows and columns: {df.columns}"
        )
    except Exception as e:
        logging.error(f"Failed to read CSV {file_path}: {e}")
        raise

    # 1. Get expected columns from BQ schema (excluding metadata columns we'll add)
    metadata_cols = {"sushi_build", "screen", "screen_hash"}
    bq_columns = [
        field.name for field in table.schema if field.name not in metadata_cols
    ]
    original_columns = set(df.columns)

    missing_cols = [col for col in bq_columns if col not in original_columns]
    extra_cols = original_columns - set(bq_columns) - metadata_cols

    if missing_cols:
        logging.info(f"Missing columns (will be added as null): {missing_cols}")
    if extra_cols:
        logging.info(f"Extra columns (will be dropped): {extra_cols}")

    # 2. Add missing columns with appropriate null values
    add_expressions = []
    for col in missing_cols:
        field = schema_map.get(col)
        if field:
            dtype = get_polars_dtype_from_bq_field(field)
            add_expressions.append(pl.lit(None).cast(dtype).alias(col))
        else:
            logging.warning(f"Could not find schema for missing column: {col}")
            add_expressions.append(pl.lit(None).alias(col))

    if add_expressions:
        df = df.with_columns(add_expressions)

    # 3. Add metadata columns
    metadata_expressions = []
    if build_name and "sushi_build" in schema_map:
        metadata_expressions.append(pl.lit(build_name).alias("sushi_build"))
    if screen and "screen" in schema_map:
        metadata_expressions.append(pl.lit(screen).alias("screen"))
    if screen and "screen_hash" in schema_map:
        unsigned_hash = farmhash.hash64(screen)
        signed_hash = ctypes.c_int64(unsigned_hash).value
        metadata_expressions.append(
            pl.lit(signed_hash, dtype=pl.Int64).alias("screen_hash")
        )

    if metadata_expressions:
        df = df.with_columns(metadata_expressions)

    # 4. Select and order columns exactly as in BQ schema
    final_columns = [field.name for field in table.schema]
    available_columns = [col for col in final_columns if col in df.columns]

    if len(available_columns) != len(final_columns):
        missing_final = set(final_columns) - set(available_columns)
        logging.warning(f"Still missing columns after processing: {missing_final}")

    df = df.select(available_columns)

    # 5. Coerce data types to match BQ schema
    df = coerce_dataframe_types(df, table)

    # 6. Validate the final DataFrame
    logging.info(f"Final DataFrame shape: {df.shape}")
    logging.info(f"Final columns: {df.columns}")
    logging.info(f"Column types: {dict(zip(df.columns, df.dtypes))}")

    # 7. Write to temporary file
    temp_fd, temp_path = tempfile.mkstemp(suffix=".csv")
    os.close(temp_fd)

    try:
        df.write_csv(temp_path, null_value="")  # Use empty string for nulls in CSV
        logging.info(f"Successfully wrote processed data to {temp_path}")
    except Exception as e:
        logging.error(f"Failed to write CSV to {temp_path}: {e}")
        raise

    return temp_path


def init_logger():
    logging.basicConfig(
        level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
    )


def get_bigquery_client():
    return bigquery.Client(project="prism-359612")


def table_exists(client, dataset_id, table_id):
    try:
        client.get_table(f"{dataset_id}.{table_id}")
        return True
    except Exception:
        return False


def delete_rows_for_build(client, dataset_id, table_id, build_name):
    query = f"""
    DELETE FROM `{dataset_id}.{table_id}`
    WHERE sushi_build = @build
    """
    job_config = bigquery.QueryJobConfig(
        query_parameters=[bigquery.ScalarQueryParameter("build", "STRING", build_name)]
    )
    logging.info(
        f"Deleting rows from {dataset_id}.{table_id} for sushi_build='{build_name}'"
    )
    client.query(query, job_config=job_config).result()


def add_screen_to_lookup_if_not_exists(client, dataset_id, screen_name):
    """Adds a screen to the lookup table if it doesn't already exist."""
    table_id = "sushi_screens_lookup"

    merge_query = f"""
    MERGE `{dataset_id}.{table_id}` T
    USING (SELECT @screen AS screen) S
    ON T.screen = S.screen
    WHEN NOT MATCHED THEN
      INSERT (screen) VALUES (screen)
    """

    job_config = bigquery.QueryJobConfig(
        query_parameters=[bigquery.ScalarQueryParameter("screen", "STRING", screen_name)]
    )

    logging.info(f"Ensuring screen '{screen_name}' exists in {dataset_id}.{table_id}")
    client.query(merge_query, job_config=job_config).result()  # Wait for completion
    logging.info(f"Screen '{screen_name}' is present in {dataset_id}.{table_id}")


def load_csv_to_bigquery(client, dataset_id, table_id, file_path, build_name, screen):
    """Load CSV to BigQuery with robust schema matching and error handling."""
    try:
        # Delete existing rows for this build if specified
        if build_name:
            delete_rows_for_build(client, dataset_id, table_id, build_name)

        table_ref = client.dataset(dataset_id).table(table_id)
        table = client.get_table(table_ref)  # Get full table object including schema

        logging.info(f"BQ table schema has {len(table.schema)} fields")
        for field in table.schema:
            logging.debug(f"  {field.name}: {field.field_type} ({field.mode})")

        # Filter and align CSV to match BQ schema
        filtered_csv = filter_csv_to_matching_columns(
            file_path, table, build_name, screen
        )

        # Configure load job with strict schema enforcement
        job_config = bigquery.LoadJobConfig(
            source_format=bigquery.SourceFormat.CSV,
            skip_leading_rows=1,
            autodetect=False,  # Use explicit schema
            schema=table.schema,
            write_disposition=bigquery.WriteDisposition.WRITE_APPEND,
            ignore_unknown_values=False,  # Be strict about unknown values
            allow_jagged_rows=False,  # Require consistent row structure
            allow_quoted_newlines=True,  # Handle quoted newlines in data
            field_delimiter=",",
            null_marker="",  # Empty string represents null
        )

        # Load the data
        with open(filtered_csv, "rb") as source_file:
            job = client.load_table_from_file(
                source_file, table_ref, job_config=job_config
            )

        # Wait for job completion and handle errors
        try:
            job.result()  # Wait for the job to complete
            logging.info(f"Successfully loaded {file_path} to {dataset_id}.{table_id}")
            logging.info(f"Loaded {job.output_rows} rows")
        except Exception as load_error:
            logging.error(f"BigQuery load job failed: {load_error}")
            if hasattr(job, "errors") and job.errors:
                for error in job.errors:
                    logging.error(f"  Error: {error}")
            raise

    except Exception as e:
        logging.error(f"Failed to load {file_path} to BigQuery: {e}")
        raise
    finally:
        # Clean up temporary file
        if "filtered_csv" in locals() and os.path.exists(filtered_csv):
            try:
                os.unlink(filtered_csv)
                logging.debug(f"Cleaned up temporary file: {filtered_csv}")
            except Exception as cleanup_error:
                logging.warning(f"Failed to clean up {filtered_csv}: {cleanup_error}")
