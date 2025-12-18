import argparse
import requests
import pandas as pd
import os
import urllib.parse
import json
import logging


def fetch_metadata(filter_dict, base_url, api_key):
    # Convert the filter to a JSON string
    filter_json = json.dumps(filter_dict)
    # URL-encode the JSON string
    encoded_filter = urllib.parse.quote(filter_json)
    # Build the full URL with the encoded filter
    url = f"{base_url}?filter={encoded_filter}"

    headers = {
        "Accept": "application/json",
        "user_key": api_key,
        "prism_key": "prism_mts",
    }

    response = requests.get(url, headers=headers)

    if response.status_code == 200:
        data = response.json()
        df = pd.DataFrame(data)
        return df
    else:
        print(f"Error: {response.status_code}")
        print(response.text)  # Print the response content for more insight
        response.raise_for_status()


def save_dataframe(df, filename, build_dir):
    # Ensure the directory exists
    os.makedirs(build_dir, exist_ok=True)

    # Define the file path
    filepath = os.path.join(build_dir, f"{filename}.csv")

    # Write the DataFrame to a CSV file
    df.to_csv(filepath, index=False)
    print(f"Sample metadata written to {filepath}")


def rename_sample_meta(df):
    rename_map = {
        "index1": "index_1",
        "index2": "index_2",
        "trt_type": "pert_type",
        "treatment": "pert_name",
        "treatment_2": "pert2_name",
        "dose": "pert_dose",
        "dose_unit": "pert_dose_unit",
        "control_barcodes": "cb_ladder",
        "replicate": "replicate_plate",
    }
    return df.rename(columns=rename_map)


def remove_sample_meta_columns(df, columns_to_remove):
    # Remove specified columns from the DataFrame
    return df.drop(columns=columns_to_remove, errors="ignore")


def filter_nan_flowcells(df, build_dir):
    # Count rows with NaN in "flowcell_names"
    nan_count = df["flowcell_names"].isna().sum()

    # Print the count if any rows are dropped and append that output to a file
    if nan_count > 0:
        # If the build_dir/logs directory does not exist, create it
        os.makedirs(os.path.join(build_dir, "logs"), exist_ok=True)
        warning_message = (
            f"Warning: {nan_count} rows with NaN in 'flowcell_names' will be dropped."
        )
        print(warning_message)
        output_path = os.path.join(build_dir, "logs/critical_output.txt")
        with open(output_path, "a") as f:
            f.write(warning_message + "\n")

        # Write the filtered DataFrame to a CSV file
        df_filtered = df.dropna(subset=["flowcell_names"])
        filtered_path = os.path.join(build_dir, "logs/na_flowcells.csv")
        df_filtered.to_csv(filtered_path, index=False, mode="a", header=False)

    # Drop rows with NaN in "flowcell_names"
    return df.dropna(subset=["flowcell_names"])


def rename_val_projects(df, screen):
    # Where pert_type is 'trt_val', rename 'x_project_id' to {screen}_VALIDATION_COMPOUNDS
    df.loc[df["pert_type"] == "trt_val", "x_project_id"] = f"{screen}_VALIDATION_COMPOUNDS"
    return df


def update_project_code(df, screen):
    # Rename the 'screen' column to have consistent naming; used when screens have overlapping treatments
    df["project_code"] = screen
    return df


def main():
    parser = argparse.ArgumentParser(
        description="Fetch data from API and save to a CSV file."
    )
    parser.add_argument(
        "--screen", "-s", type=str, required=True, help="Screen code to query the API"
    )
    parser.add_argument(
        "--api_key", "-k", type=str, required=True, help="API key for authentication"
    )
    parser.add_argument(
        "--build_dir",
        "-b",
        type=str,
        required=True,
        help="Directory to save the CSV file",
    )
    parser.add_argument(
        "--pert_plates",
        "-p",
        type=str,
        required=False,
        help="Pert plates in the screen to use. If not provided, all pert plates will be used.",
    )
    parser.add_argument(
        "--screen_type",
        "-t",
        type=str,
        required=False,
        default="MTS_SEQ",
    )

    args = parser.parse_args()
    screen = args.screen
    screen_type = args.screen_type
    api_key = args.api_key
    build_dir = args.build_dir
    pert_plates = (
        args.pert_plates.replace(" ", "").split(",") if args.pert_plates else None
    )

    try:
        df = (
            fetch_metadata(
                filter_dict={"where": {"project_code": {"like": f"%{screen}%"}}},
                base_url="https://api.clue.io/api/v_seq_metadata",
                api_key=api_key,
            )
            .pipe(rename_sample_meta)
            .pipe(rename_val_projects, screen)
            .pipe(update_project_code, screen)
            .pipe(
                remove_sample_meta_columns,
                [
                    "id",
                    "pert_platemap_id",
                    "prism_pcr_plate_id",
                    "pcr_plate_well_id",
                    "index_set",
                ],
            )
            .pipe(filter_nan_flowcells, build_dir)
        )
        # Subset to provided pert plates if any
        if pert_plates:
            print(
                f"pert_plate provided, filtering sample_metadata to pert plates: {pert_plates}"
            )
            df = df[df["pert_plate"].isin(pert_plates)]

        print("Retrieved following sample_metadata from COMET: ")
        print(df.head())

        # Filter pert2* columns if screen type is not CPS_SEQ
        if screen_type != "CPS_SEQ":
            logging.info("Screen type is not CPS_SEQ, filtering pert2* columns.")
            logging.info(f"Columns before filtering: {df.columns.tolist()}")

            # Store original columns
            original_cols = set(df.columns)

            # Remove pert2* columns (only if they exist)
            pert2_cols = [col for col in df.columns if col.startswith("pert2")]
            if pert2_cols:
                df = df.drop(columns=pert2_cols)

            # Remove additional columns (only if they exist)
            additional_cols = ["treatment_2", "is_combination"]
            existing_additional_cols = [
                col for col in additional_cols if col in df.columns
            ]
            if existing_additional_cols:
                df = df.drop(columns=existing_additional_cols)

            # Calculate what was actually removed
            final_cols = set(df.columns)
            cols_removed = list(original_cols - final_cols)

            logging.info(f"Columns removed: {cols_removed}")
            logging.info(f"Columns after filtering: {df.columns.tolist()}")

        # Save the DataFrame to a CSV file
        save_dataframe(df=df, filename="sample_meta", build_dir=build_dir)
    except Exception as e:
        print(f"An error occurred: {e}")

    try:
        df = fetch_metadata(
            filter_dict={"where": {"screen": {"like": f"%{screen}%"}}},
            base_url="https://api.clue.io/api/v_c_assay_plate_skipped_well",
            api_key=api_key,
        )
        # If skipped wells is not empty, save it to a csv file
        if not df.empty:
            print("Retrieved following skipped_wells from COMET: ")
            print(df.head())
            save_dataframe(df=df, filename="skipped_wells", build_dir=build_dir)
        else:
            print(f"No skipped wells found for screen {screen}.")
    except Exception as e:
        print(f"An error occurred: {e}")


if __name__ == "__main__":
    logging.basicConfig(
        level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
    )
    main()
