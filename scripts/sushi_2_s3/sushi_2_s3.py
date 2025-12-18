import sys
import os

# This line gets the absolute path to the current script file
current_file_path = os.path.abspath(__file__)

# This line goes "up" two levels (from sushi_2_s3.py -> sushi_2_s3/ -> scripts/)
project_root = os.path.dirname(os.path.dirname(os.path.dirname(current_file_path)))

# Add the project root to the list of places Python looks for modules
sys.path.insert(0, project_root)

import pandas as pd
import argparse
import glob
import json
import boto3
from botocore.exceptions import NoCredentialsError
import logging
from sushilib.prism_tools import SushiBuild
import polars as pl

logger = logging.getLogger("sushi_2_s3")


class Config:
    def __init__(self, config_dict):
        self.MERGE_PATTERNS = None
        for key, value in config_dict.items():
            setattr(self, key, value)


def build_parser():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        "--s3_bucket", "-s", help="S3 bucket for output.", default="macchiato.clue.io"
    )
    parser.add_argument(
        "--build_path", "-b", help="Path to the build directory.", required=True
    )
    parser.add_argument(
        "--verbose",
        "-v",
        action="store_true",
        help="Enable verbose logging.",
    )
    return parser


def read_config(config_path):
    with open(config_path, "r") as f:
        config_dict = json.load(f)
    return Config(config_dict)


def read_build_file(search_pattern, args):
    fstr = os.path.join(args.build_path, search_pattern)
    fmatch = glob.glob(fstr)
    assert len(fmatch) == 1, "Too many files found: {}".format(fmatch)
    return pd.read_csv(fmatch[0])


def generate_compound_key(df, filt_cols: list):
    # Filter out the positive and vehicle controls
    df = df[~df["pert_type"].isin(["trt_poscon", "ctl_vehicle"])]

    # Select the columns we need
    df = df[filt_cols]

    # Drop duplicates
    distinct_df = df.drop_duplicates().reset_index(drop=True)
    return distinct_df


def generate_merge_key(df, merge_patterns):
    # Filter out the positive and vehicle controls
    df = df[~df["pert_type"].isin(["trt_poscon", "ctl_vehicle"])]

    # Get a list of x_project_ids
    project_ids = df["x_project_id"].unique()

    # Creat a dataframe with each row containing a project_id and each of the merge_patterns
    merge_keys = []
    for project_id in project_ids:
        for pattern in merge_patterns:
            merge_keys.append([project_id, pattern])

    # Create a dataframe from the merge_keys
    merge_df = pd.DataFrame(merge_keys, columns=["x_project_id", "merge_pattern"])

    # Ensure uniqueness
    distinct_df = merge_df.drop_duplicates().reset_index(drop=True)

    # Ensure that x_project_id is a string
    distinct_df["x_project_id"] = distinct_df["x_project_id"].astype(str)

    # Return the merge_df
    return distinct_df


def generate_pert_plate_project_list(df):
    # Remove "x_project_id" when it contains the string "CONTROLS"
    df = df[~df["x_project_id"].str.contains("CONTROLS")]

    # Select only the relevant columns and remove any duplicate rows.
    df = df[["x_project_id", "pert_plate"]].drop_duplicates()

    # Group by 'x_project_id' and aggregate the 'pert_plate' values into lists
    grouped = df.groupby("x_project_id")["pert_plate"].apply(list).reset_index()

    # Convert the DataFrame to a list of dictionaries (records)
    output_list = grouped.to_dict(orient="records")
    return output_list


# Turn the key_df into a json object
def key_to_json(
    df: pd.DataFrame, unique_by: list = None, output_path: str = None
) -> str:
    """
    Convert a pandas DataFrame to a pretty-printed JSON string (records orientation) and write it to a file.

    Parameters
    ----------
    df : pd.DataFrame
        The input DataFrame containing columns like "pert_name", "pert_id", "pert_plate",
        "x_project_id", "pert_dose", etc.
    unique_by : list, optional
        If provided, duplicate rows based on these columns will be dropped before conversion.
    output_path : str, optional
        If provided, the resulting JSON string will be written to this file. Necessary directories
        will be created if they do not exist.

    Returns
    -------
    str
        A pretty-printed JSON string representing the DataFrame.
    """
    if unique_by is not None:
        df = df.drop_duplicates(subset=unique_by)

    json_str = df.to_json(orient="records", indent=2)

    if output_path:
        # Create the directory if it doesn't exist
        directory = os.path.dirname(output_path)
        if directory and not os.path.exists(directory):
            os.makedirs(directory, exist_ok=True)
        with open(output_path, "w") as f:
            f.write(json_str)

    return json_str


# Sync all contents of current directory to S3
def sync_to_s3(local_dir, s3_bucket, s3_prefix, exclude_pattern=None):
    """Sync the local directory to S3."""
    s3 = boto3.client("s3")

    for root, dirs, files in os.walk(local_dir):
        for file in files:
            logger.info(f"Checking {file} for upload...")
            if exclude_pattern not in file and not file.startswith("."):
                logger.info(f"Uploading {file} to S3...")
                local_path = os.path.join(root, file)
                relative_path = os.path.relpath(local_path, local_dir)
                s3_path = os.path.join(s3_prefix, relative_path)
                try:
                    logger.info(f"Uploading {local_path} to s3://{s3_bucket}/{s3_path}")
                    s3.upload_file(local_path, s3_bucket, s3_path)
                    logger.info(f"Uploaded {local_path} to s3://{s3_bucket}/{s3_path}")
                except NoCredentialsError:
                    logger.info("AWS credentials not found.")
                    raise
            else:
                logger.info(f"Excluding upload of {file}")


def main(args):
    # Get the build path
    build_path = args.build_path

    # Read the config file
    config_path = args.build_path + "/config.json"
    config = read_config(config_path)

    # Get build name from config
    build_name = config.BUILD_NAME

    # Construct the s3 prefix
    s3_prefix = f"builds/{build_name}/"

    # Read the sample_meta file
    sample_meta_path = os.path.join(build_path, config.SAMPLE_META)
    sample_meta = pd.read_csv(sample_meta_path)

    # Generate the key_df
    logger.info("Generating compound key...")

    # From screen type, decide on which columns to use for compound key
    if config.SCREEN_TYPE == "CPS_SEQ":
        compound_key_cols = [
            "x_project_id",
            "pert_plate",
            "pert_name",
            "pert_id",
            "pert2_name",
            "pert2_id",
        ]
    else:
        compound_key_cols = ["x_project_id", "pert_plate", "pert_name", "pert_id"]

    key_df = generate_compound_key(sample_meta, filt_cols=compound_key_cols)

    # Convert the key_df to json and write to the build directory
    logger.info("Writing compound key to JSON...")
    key_to_json(key_df, output_path=os.path.join(build_path, "compound_key.json"))

    # Create the unique projects key
    unique_df = key_df[["x_project_id"]].drop_duplicates().reset_index(drop=True)
    key_to_json(unique_df, output_path=os.path.join(build_path, "unique_projects.json"))

    # Get the list of merge patterns
    merge_patterns = config.MERGE_PATTERNS
    merge_patterns = merge_patterns.split(",")

    # Generate the merge key
    logger.info("Generating merge key...")
    merge_key_df = generate_merge_key(sample_meta, merge_patterns)

    # Convert the merge_key_df to json and write to the build directory
    logger.info("Writing merge key to JSON...")
    key_to_json(merge_key_df, output_path=os.path.join(build_path, "merge_key.json"))

    # Generate the pert_plate_project key
    logger.info("Generating pert_plate_project key...")
    pert_plate_project_list = generate_pert_plate_project_list(sample_meta)
    with open(os.path.join(build_path, "pert_plate_project_key.json"), "w") as f:
        json.dump(pert_plate_project_list, f, indent=2)

    # Sync the build directory to S3
    logger.info(f"Syncing {build_path} to {s3_prefix}...")
    sync_to_s3(
        local_dir=build_path,
        s3_bucket=args.s3_bucket,
        s3_prefix=s3_prefix,
        exclude_pattern="raw_counts",
    )

    logger.info("Completed sync to S3.")


if __name__ == "__main__":
    args = build_parser().parse_args(sys.argv[1:])
    level = logging.DEBUG if args.verbose else logging.INFO
    logging.basicConfig(level=level)
    logger.info("args:  {}".format(args))

    main(args)
