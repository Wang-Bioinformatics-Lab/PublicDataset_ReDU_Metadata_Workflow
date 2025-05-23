#!/usr/bin/env python3
import argparse
import requests
import pandas as pd
from io import StringIO
from tqdm import tqdm
import re
from urllib.parse import unquote

def create_usi(row):
    return f"mzspec:NORMAN-{row['uuid']}:{row['file_paths']}"


def _get_existing_datasets(path_to_file):
    try:
        existing_datasets = pd.read_csv(path_to_file, sep="\t")
        existing_datasets = set(existing_datasets['datasets'])
        return existing_datasets
    except:
        return set()
    
def process_dataset_files(
    df_files,
    internal_id=None,
    uuid=None,
    title=None,
    filter_extensions=True,
):
    file_cols = ["data_independent", "data_dependent", "data_fullscan"]
    missing_cols = [c for c in file_cols if c not in df_files.columns]
    if missing_cols:
        print(f"Dataset {internal_id or '(unknown)'} missing columns: {missing_cols}. Skipping.")
        return pd.DataFrame()

    # reshape → one file per row
    df = (
        df_files[["sample_id"] + file_cols]
        .melt(id_vars=["sample_id"], value_name="file_urls")
        .dropna(subset=["file_urls"])
    )
    df = df[df["file_urls"].astype(str).str.strip() != ""]

    # filenames, paths
    df["file_names"] = df["file_urls"].apply(extract_file_name)
    df["file_paths"] = df["file_urls"].apply(extract_file_path)

    # optional metadata columns
    if internal_id is not None:
        df["internal_id"] = internal_id
    if uuid is not None:
        df["uuid"] = uuid
    if title is not None:
        df["title"] = title

    # extension filter
    if filter_extensions:
        exts = (".mzml", ".mzxml", ".cdf", ".raw", ".wiff", ".d")
        df = df[df["file_names"].str.lower().str.endswith(exts)]

    # USI
    df["USI"] = df.apply(create_usi, axis=1)

    return df

def extract_file_name(url):
    # Step 1: Remove everything before and including the first occurrence of "sample/{some_number}/"
    url_after_sample = re.sub(r'.*sample/\d+/','', url)
    
    # Step 2: Remove everything from "?VersionId" onwards if it exists
    url_after_version = re.sub(r'\?VersionId.*', '', url_after_sample)
    
    return unquote(url_after_version)

def extract_file_path(url):
    # Step 1: Remove everything before and including the first occurrence of "sample/{some_number}/"
    file_path = re.sub(r'https://files.dsfp.norman-data.eu/','', url)
        
    return unquote(file_path)

def main(output_filename, study_id, filter_extensions, existing_datasets):
    # Step 1: Fetch the list of datasets
    datasets_url = "https://dsfp.norman-data.eu/api/1/metastore/schemas/dataset/all"
    print(f"Fetching datasets from {datasets_url}")
    response = requests.get(datasets_url)
    response.raise_for_status()  # Raise an error for bad status codes
    datasets = response.json()
    print(f"Fetched {len(datasets)} datasets")

    # Filter datasets based on study_id
    if study_id != "ALL":
        datasets = [ds for ds in datasets if ds['uuid'] == study_id]
        print(f"Processing dataset with UUID: {study_id}")
    else:
        print("Processing all datasets")

    dfs = []  # List to hold DataFrames for each dataset
    errors = []  # List to hold errors

    # Process each dataset with a progress bar
    for ds in tqdm(datasets, desc="Processing datasets", unit="dataset"):
        try:
            uuid = ds['uuid']
            internal_id = ds['internal_id']
            title = ds['title']
            print(f"Processing dataset: {title} (UUID: {uuid}, Internal ID: {internal_id})")

            if internal_id in existing_datasets:
                print(f"Dataset {internal_id} already indexed. Skipping.")
                continue

            # Build URL to get the CSV with file info
            file_url = f"https://dsfp.norman-data.eu/data/{internal_id}/files.csv"
            print(f"Fetching file data from {file_url}")
            file_response = requests.get(file_url)

            if file_response.status_code == 200:
                # Read CSV data into DataFrame
                csv_data = StringIO(file_response.text)
                df_files = pd.read_csv(csv_data)
                print(f"Fetched file data for dataset {internal_id}")

                # Process the DataFrame
                df_melted = process_dataset_files(
                    df_files,
                    internal_id=internal_id,
                    uuid=uuid,
                    title=title,
                    filter_extensions=filter_extensions
                )

                # Append the processed DataFrame to our list
                dfs.append(df_melted)
                print(f"Processed dataset {internal_id} successfully")
            else:
                error_message = f"Failed to download files for dataset {internal_id}. Status code: {file_response.status_code}"
                print(error_message)
                errors.append(error_message)
        except Exception as e:
            error_message = f"Error processing dataset {internal_id}: {str(e)}"
            print(error_message)
            errors.append(error_message)

    # Step 3: Combine all per-dataset DataFrames (if any)
    if dfs:
        print("Combining all dataset DataFrames")
        combined_df = pd.concat(dfs, ignore_index=True)

        # rename columns
        combined_df = combined_df.rename(columns={'uuid': 'study_id', 'file_names': 'file_path', 'file_urls': 'file_url'})
        combined_df = combined_df[['study_id', 'internal_id', 'sample_id', 'file_path', 'file_url', 'USI']]

        # Step 4: Save the result as a TSV file
        combined_df.to_csv(output_filename, sep='\t', index=False)
        print(f"\nCombined data saved to {output_filename}")
    else:
        print("No file data found across datasets.")

    # Print errors 
    if errors:
        print("\nErrors encountered during processing:")
        for error in errors:
            print(error)



if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Download and combine dataset file information into a TSV file."
    )
    parser.add_argument(
        "-o", "--output",
        type=str,
        required=True,
        help="Output TSV filename (e.g., combined_files.tsv)"
    )
    parser.add_argument(
        "--study_id",
        type=str,
        default="ALL",
        help="UUID of the study to process, or 'ALL' to process all datasets."
    )
    parser.add_argument(
        "--filter_extensions", 
        type=str, 
        help='Filter extensions to ".mzml", ".mzxml", ".cdf", ".raw", ".wiff", and ".d".', 
        default='False'
        )
    parser.add_argument(
        "--existing_datasets", 
        type=str, 
        help="path to a file of datasets already indexed", 
        default="none"
        )
    
    args = parser.parse_args()

    existing_datasets = _get_existing_datasets(args.existing_datasets)

    main(args.output, args.study_id, args.filter_extensions, existing_datasets)
