import os
import pandas as pd
import argparse
import collections
import requests
from time import sleep
import glob
from io import StringIO
from subprocess import PIPE, run
import json
from pathlib import Path

ccms_peak_link = "https://datasetcache.gnps2.org/datasette/datasette/database/uniquemri.csv?_sort=usi&dataset__exact=" # MSV000081468&filepath__endswith=%25.mz%25ML&_size=max"
gnps_column_names_added = ['USI']

def _make_usi_from_filename(filename, dataset_id):
    # Replace initial "f." with "mzspec:"
    if filename.startswith("f.MSV"):
        usi = "mzspec:" + filename[2:]
        # Replace the first "/" with ":"
        first_slash_index = usi.find('/')
        if first_slash_index != -1:
            usi = filename[:first_slash_index] + ':' + filename[first_slash_index+1:]
    elif filename.startswith("f."):
        usi = "mzspec:" + dataset_id + ':' + filename[2:]
    else:
        print('Do filenames not start with f. anymore?')

    return usi

def longest_common_dir_suffix_score(meta_dirs, cand_dirs):
    score = 0
    for md, cd in zip(reversed(meta_dirs), reversed(cand_dirs)):
        if md == cd:
            score += 1
        else:
            break
    return score


def _match_filenames_and_add_usi(dataset_metadata_df):
         
    dataset = dataset_metadata_df['ATTRIBUTE_DatasetAccession'].iloc[0]
    print(f"Checking dataset {dataset} (value extracting from column.)")

    url_f = "{}{}&filepath__endswith=%25.mz%25ML&_size=max".format(ccms_peak_link, dataset) 
    print("Fetching {}".format(url_f))

    # Request URL and retry 3 times after waiting 5 seconds if ccms_df does not have a column named 'filepath'
    retry_attempts = 3
    for attempt in range(retry_attempts):
        print(f"Checking dataset {dataset} by going to the dataset cache")
        dataset_files_response = requests.get(url_f)

        csvStringIO = StringIO(dataset_files_response.text)
        ccms_df = pd.read_csv(csvStringIO)
        
        if 'filepath' in ccms_df.columns:
            break
        else:
            print(f"Attempt {attempt + 1} failed. Retrying in 5 seconds...")
            sleep(5)

    if len(ccms_df) > 0:
        print("Received a dataframe with {} rows.".format(len(ccms_df)))
        print(ccms_df)
        pass
    else:
        print("Error: the size of the input is too small, length:", len(ccms_df))
        return None

    ccms_df["query_path"] = ccms_df["filepath"].apply(lambda x: os.path.basename(x))
    ccms_df["_dir_parts"] = ccms_df["filepath"].apply(lambda p: tuple(Path(p).parts[:-1]))

    metadata_row_list = dataset_metadata_df.to_dict('records')
    output_row_list = []

    for metadata_row in metadata_row_list:

        print(f"Processing {metadata_row['filename']}")

        # Metadata path pieces
        meta_path = str(metadata_row["filename"])
        meta_basename = os.path.basename(meta_path)
        # extension-only case handling (keep as you wanted)
        meta_basename2 = meta_basename[:-3] + "ML" if len(meta_basename) >= 3 else meta_basename
        meta_dir_parts = tuple(Path(meta_path).parts[:-1])  # () if no dirs present

        # 1) Collect candidates by basename or basename2
        cand_mask = (ccms_df["query_path"] == meta_basename) | (ccms_df["query_path"] == meta_basename2)
        candidates = ccms_df.loc[cand_mask, ["filepath", "_dir_parts", "query_path"]].copy()

        if candidates.empty:
            found_file_paths = []
            print(f"Found file paths: {found_file_paths}")
        else:
            if len(candidates) == 1:
                # single candidate → just use it
                found_file_paths = candidates["filepath"].tolist()
                print(f"Found file paths: {found_file_paths}")
            else:
                if len(meta_dir_parts) == 0:
                    # multiple candidates but no directory info in metadata → keep all for preference filter
                    found_file_paths = candidates["filepath"].tolist()
                    print(f"Found file paths (no dir context): {found_file_paths}")
                else:
                    # 2) Score candidates by longest common directory suffix with the metadata path
                    candidates["_dir_score"] = candidates["_dir_parts"].apply(
                        lambda dp: longest_common_dir_suffix_score(meta_dir_parts, dp)
                    )
                    # 3) Keep only the best-scoring candidate(s)
                    best_score = candidates["_dir_score"].max()
                    found_file_paths = candidates.loc[candidates["_dir_score"] == best_score, "filepath"].tolist()
                    print(f"Found file paths (scored): {found_file_paths}")


        # 2) Score candidates by longest common directory suffix with the metadata path
        candidates["_dir_score"] = candidates["_dir_parts"].apply(
            lambda dp: longest_common_dir_suffix_score(meta_dir_parts, dp)
        )

        # 3) Keep only the best-scoring candidate(s)
        best_score = candidates["_dir_score"].max()
        best = candidates.loc[candidates["_dir_score"] == best_score, "filepath"].tolist()

        found_file_paths = best

        print(f"Found file paths: {found_file_paths}")

        if len(found_file_paths) > 0:
            preferred_directories = ['ccms_peak', 'peak', 'raw']
            selected_path = None
            for preferred_dir in preferred_directories:
                for path in found_file_paths:
                    if preferred_dir in path:
                        selected_path = path
                        break
                if selected_path:
                    break
            if not selected_path:
                selected_path = found_file_paths[0]

            print("Found match", selected_path)
            metadata_row["filename"] = "f." + selected_path
            metadata_row["USI"] = _make_usi_from_filename(
                metadata_row["filename"],
                metadata_row["ATTRIBUTE_DatasetAccession"]
            )
            output_row_list.append(metadata_row)
        else:
            # Didn't find or is ambiguous
            continue

    return pd.DataFrame(output_row_list)

def main():
    # parsing args
    parser = argparse.ArgumentParser(description='GNPS Name Matcher')
    parser.add_argument('passed_file_names', help='Input TSV file')
    parser.add_argument('metadata_folder')
    parser.add_argument('output_filename')
    parser.add_argument('path_allowed_terms_json')

    args = parser.parse_args()
     
    ccms_filenames = collections.defaultdict(set)
    
    if args.passed_file_names == 'all':
        passed_file_names = glob.glob(f"{args.metadata_folder}/*.tsv")
    elif args.passed_file_names == 'single':
        passed_file_names = []
    else:
        # Read the TSV file and specify the delimiter as a tab
        df = pd.read_csv(args.passed_file_names, delimiter='\t', header=None, names=['Name'])
        # Extract the names from a specific column (e.g., column 'Name')
        passed_file_names = df['Name'].tolist()

    
    with open(args.path_allowed_terms_json, 'r') as file:
        allowed_terms_json = json.load(file)

    gnps_column_names = ["ATTRIBUTE_DatasetAccession"] + list(set(allowed_terms_json.keys()) - {'USI', 'MassiveID'})

    print("echo Iterating though rows now")
    
    all_metadata_list = []

    for file in passed_file_names:
        # print("Length of visit is ",len(visit))
        print("Working on ", file)
        # print("echo Working on")
        # csv_path = os.path.join(current_dir, './data.csv' file)
        dataset_metadata_df = pd.read_csv( file , delimiter='\t')

        print(f"reading potential ReDU data with {len(dataset_metadata_df)} rows.")
        
        #Renaming the coloumn, Matching common columns and rearranging them in same order to final file
        dataset_metadata_df = dataset_metadata_df.rename(columns={'MassiveID': 'ATTRIBUTE_DatasetAccession'})
        common_cols = list(set(gnps_column_names).intersection(set(dataset_metadata_df.columns)))
        dataset_metadata_df = dataset_metadata_df.loc[:, common_cols]
        try:
            dataset_metadata_df = dataset_metadata_df[gnps_column_names]
        except KeyError:
            print(f"Skipping file {file} due to a TypeError.")
            continue

        # Matching the metadata
        enriched_metadata_df = _match_filenames_and_add_usi(dataset_metadata_df)
        if enriched_metadata_df is not None:

            print(f"Selected {len(enriched_metadata_df)} rows with columns: {list(enriched_metadata_df.columns)}")

            all_metadata_list.append(enriched_metadata_df)

    # Create a DataFrame from the list with headers
    merged_metadata_df = pd.concat(all_metadata_list)
    merged_metadata_df = merged_metadata_df[gnps_column_names + gnps_column_names_added]

    # Save the DataFrame to a TSV file without column names
    merged_metadata_df.to_csv(args.output_filename, sep='\t', index=False)


if __name__ == '__main__':
    main()