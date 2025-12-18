import polars as pl
import pandas as pd
import yaml
from pathlib import Path
from typing import Dict, Any, Optional

# --- Define the default path to the schema file ---
# This finds the directory of the current script and then navigates to the schema file.
# This makes the script work regardless of where it's run from.
try:
    # This will work when the script is run as a file
    DEFAULT_SCHEMA_PATH = Path(__file__).parent.parent / "sushi_io" / "schema.yml"
except NameError:
    # This is a fallback for interactive environments (like Jupyter) where __file__ is not defined.
    # It assumes a project structure relative to the current working directory.
    DEFAULT_SCHEMA_PATH = Path.cwd()  / "sushi_io" / "schema.yml"


def read_polars(csv_path: str, schema_path: Optional[str] = None) -> pl.DataFrame:
    """
    Reads a CSV into a Polars DataFrame, enforcing data types from a YAML schema.

    Args:
        csv_path: The path to the input CSV file.
        schema_path: Optional path to the master YAML schema file.
                     If None, defaults to the path defined in this module.

    Returns:
        A Polars DataFrame with the specified data types.
    """
    print("--- Running Polars CSV Reader ---")
    if schema_path is None:
        schema_path = DEFAULT_SCHEMA_PATH
        print(f"📘 Using default schema path: {schema_path}")

    try:
        # Define the mapping from your abstract YAML types to Polars types
        type_mapping = {
            'string': pl.Utf8,
            'integer': pl.Int64,
            'float': pl.Float64,
            'boolean': pl.Boolean,
            'datetime': pl.Datetime,  # For future use if needed
        }

        # Read the master schema file
        with open(schema_path, 'r') as f:
            schema = yaml.safe_load(f)

        # Create the dtypes dictionary for Polars by mapping schema types
        polars_dtypes = {col: type_mapping[dtype] for col, dtype in schema.items()}
        print(f"✅ Schema loaded and mapped for {len(polars_dtypes)} columns.")

        # Read the CSV file, enforcing the schema
        df = pl.read_csv(csv_path, dtypes=polars_dtypes, null_values="NA")
        print(f"✅ Success! Loaded {df.shape[0]} rows into a Polars DataFrame.\n")
        return df

    except FileNotFoundError as e:
        print(f"❌ Error: File not found - {e}")
        return pl.DataFrame()
    except Exception as e:
        print(f"❌ An unexpected error occurred: {e}")
        return pl.DataFrame()


def read_pandas(csv_path: str, schema_path: Optional[str] = None) -> pd.DataFrame:
    """
    Reads a CSV into a Pandas DataFrame, enforcing data types from a YAML schema.

    Args:
        csv_path: The path to the input CSV file.
        schema_path: Optional path to the master YAML schema file.
                     If None, defaults to the path defined in this module.

    Returns:
        A Pandas DataFrame with the specified data types.
    """
    print("--- Running Pandas CSV Reader ---")
    if schema_path is None:
        schema_path = DEFAULT_SCHEMA_PATH
        print(f"📘 Using default schema path: {schema_path}")

    try:
        # Define the mapping from your abstract YAML types to modern, nullable Pandas dtypes
        type_mapping = {
            'string': pd.StringDtype(),
            'integer': pd.Int64Dtype(),
            'float': pd.Float64Dtype(),
            'boolean': pd.BooleanDtype(),
        }

        # Read the master schema file
        with open(schema_path, 'r') as f:
            schema = yaml.safe_load(f)

        # Create the dtypes dictionary for Pandas
        # Note: We separate date columns as they are handled by a different parameter
        pandas_dtypes = {}
        date_columns = []
        for col, dtype in schema.items():
            if dtype == 'datetime':
                date_columns.append(col)
            else:
                pandas_dtypes[col] = type_mapping.get(dtype)

        print(f"✅ Schema loaded and mapped for {len(pandas_dtypes)} columns.")
        if date_columns:
            print(f"   (Will parse {len(date_columns)} date columns separately)")

        # Read the CSV file, enforcing the schema
        df = pd.read_csv(
            csv_path,
            dtype=pandas_dtypes,
            parse_dates=date_columns  # Use the dedicated parameter for dates
        )
        print(f"✅ Success! Loaded {len(df)} rows into a Pandas DataFrame.\n")
        return df

    except FileNotFoundError as e:
        print(f"❌ Error: File not found - {e}")
        return pd.DataFrame()
    except Exception as e:
        print(f"❌ An unexpected error occurred: {e}")
        return pd.DataFrame()


if __name__ == '__main__':
    # --- Create Dummy Files for Demonstration ---

    # In a real scenario, the schema would live at the default path.
    # For this demo, we'll create the directory structure and the dummy files.

    # 1. Define where the dummy schema should go
    demo_schema_dir = Path("./sushi/utils/io")
    demo_schema_path = demo_schema_dir / "schema.yaml"

    # Create the directories if they don't exist
    demo_schema_dir.mkdir(parents=True, exist_ok=True)

    # Set the default path for this demo script to point to our dummy file
    DEFAULT_SCHEMA_PATH = demo_schema_path

    # 2. Create a dummy schema.yaml
    dummy_schema_content = """
    pert_id: string
    pert_dose: float
    day: integer
    is_combination: boolean
    """
    demo_schema_path.write_text(dummy_schema_content.strip())

    # 3. Create a dummy data.csv file that matches the dummy schema
    dummy_csv_path = Path("data.csv")
    dummy_csv_content = """pert_id,pert_dose,day,is_combination
BRD-K12345,10.0,7,True
BRD-K67890,5.5,3,False
BRD-A11223,,5,true
"""
    dummy_csv_path.write_text(dummy_csv_content.strip())

    print("--- DEMONSTRATION ---")
    print(f"Created dummy '{dummy_csv_path}' and '{demo_schema_path}' for this example.\n")

    # --- Use the Polars Function (without specifying schema_path) ---
    polars_df = read_polars(csv_path="data.csv")
    if not polars_df.is_empty():
        print("Polars DataFrame Info:")
        print(polars_df)
        print("Polars Schema:")
        print(polars_df.schema, "\n")

    # --- Use the Pandas Function (without specifying schema_path) ---
    pandas_df = read_pandas(csv_path="data.csv")
    if not pandas_df.empty:
        print("Pandas DataFrame Info:")
        print(pandas_df)
        print("Pandas dtypes:")
        print(pandas_df.dtypes)

