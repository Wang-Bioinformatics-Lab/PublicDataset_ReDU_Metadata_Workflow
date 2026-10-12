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
            usi = usi[:first_slash_index] + ':' + usi[first_slash_index+1:]
    elif filename.startswith("f."):
        usi = "mzspec:" + dataset_id + ':' + filename[2:]
    else:
        print('Do filenames not start with f. anymore?')
        usi = "mzspec:" + dataset_id + ':' + filename

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
    ccms_df = None
    for attempt in range(retry_attempts):
        print(f"Checking dataset {dataset} by going to the dataset cache")
        try:
            dataset_files_response = requests.get(url_f, timeout=300)
            dataset_files_response.raise_for_status()
            attempt_df = pd.read_csv(StringIO(dataset_files_response.text))
            if 'filepath' in attempt_df.columns:
                ccms_df = attempt_df
                break
        except (requests.RequestException, pd.errors.ParserError, pd.errors.EmptyDataError) as e:
            print(f"Attempt {attempt + 1} error: {e}")
        print(f"Attempt {attempt + 1} failed. Retrying in 5 seconds...")
        sleep(5)

    if ccms_df is None:
        print(f"Error: could not fetch the file list for dataset {dataset}, skipping it")
        return None

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

        # 2) Score candidates by longest common directory suffix with the metadata path
        #    (all score 0 when the metadata has no directory info)
        candidates["_dir_score"] = candidates["_dir_parts"].apply(
            lambda dp: longest_common_dir_suffix_score(meta_dir_parts, dp)
        )

        # 3) Keep only the best-scoring candidate(s)
        best_score = candidates["_dir_score"].max()
        found_file_paths = candidates.loc[candidates["_dir_score"] == best_score, "filepath"].tolist()

        print(f"Found file paths: {found_file_paths}")

        if len(found_file_paths) > 0:
            # The same file can be present in several top-level folders (e.g. the
            # converted ccms_peak copy and the uploaded peak/raw copy); prefer them in
            # this order. Folders are compared as whole path segments.
            preferred_directories = ['ccms_peak', 'peak', 'raw']
            selected_paths = found_file_paths
            for preferred_dir in preferred_directories:
                in_dir = [path for path in found_file_paths if preferred_dir in Path(path).parts[:-1]]
                if in_dir:
                    selected_paths = in_dir
                    break

            if len(selected_paths) > 1:
                # Several distinct files fit equally well (e.g. pos/s1.mzML and
                # neg/s1.mzML for metadata 's1.mzML'); do not guess
                print(f"Ambiguous match for {meta_path}, skipping: {selected_paths}")
                continue
            selected_path = selected_paths[0]

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
        except KeyError as e:
            print(f"Skipping file {file} due to missing columns: {e}")
            continue

        # Matching the metadata
        enriched_metadata_df = _match_filenames_and_add_usi(dataset_metadata_df)
        if enriched_metadata_df is not None:

            print(f"Selected {len(enriched_metadata_df)} rows with columns: {list(enriched_metadata_df.columns)}")

            all_metadata_list.append(enriched_metadata_df)

    # Create a DataFrame from the list with headers (header only if nothing matched,
    # e.g. no new MASST rows)
    if all_metadata_list:
        merged_metadata_df = pd.concat(all_metadata_list)
    else:
        print("No metadata rows matched to files")
        merged_metadata_df = pd.DataFrame(columns=gnps_column_names + gnps_column_names_added)
    merged_metadata_df = merged_metadata_df[gnps_column_names + gnps_column_names_added]

    # Save the DataFrame to a TSV file without column names
    merged_metadata_df.to_csv(args.output_filename, sep='\t', index=False)


if __name__ == '__main__':
    main()