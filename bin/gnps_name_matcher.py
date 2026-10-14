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
import re
from pathlib import Path

ccms_peak_link = "https://datasetcache.gnps2.org/datasette/datasette/database/uniquemri.csv?_sort=usi&dataset__exact=" # MSV000081468&filepath__endswith=%25.mz%25ML&_size=max"
# MetadataFileDataset: accession of the dataset the metadata file was found in (MassIVE)
# or named after (GitHub); merge_metadata.py prefers a dataset's own file on conflicts
gnps_column_names_added = ['USI', 'MetadataFileDataset']

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


# dataset accessions of the dataset cache: MSV/RMSV (MassIVE), MTBLS, ST (MASST rows)
_ACCESSION = re.compile(r'^[A-Z]+\d+$')
_MZ_EXTENSIONS = ('.mzml', '.mzxml')
PREFERRED_DIRECTORIES = ['ccms_peak', 'peak', 'raw']

# matching rules, in the order they are tried for a metadata row; a later rule is only
# used when the earlier ones find no file of that name at all
#   exact        same file name (or the x.mzXML -> x.mzML converted copy, as before)
#   stem_ext     same name, the other of .mzML/.mzXML (x.mzML listed, only x.mzXML uploaded)
#   case         same name ignoring upper/lower case (also with the other of .mzML/.mzXML)
MATCH_STATS = collections.Counter()
AMBIGUOUS = []
# one entry per metadata row (written with --match_report)
MATCH_REPORT = []


# The dataset cache's whole file table (https://datasetcache.gnps2.org/dataset/uniquemri, a
# ~1 GB export the cache regenerates about every 3 h), downloaded once per run by the
# workflow (--files_table): dataset -> its mz(X)ML file rows sorted by usi, i.e. what the
# per-dataset query below returns. One query per dataset took 20-90 min per run and a query
# that timed out dropped all rows of its dataset (MSV000087612, 192 rows, on 2026-10-14).
# Set by load_files_table().
FILES_TABLE = None
# the condition of the per-dataset query, filepath__endswith=%.mz%ML, i.e. SQLite
# "filepath LIKE '%.mz%ML'" (LIKE ignores upper/lower case of ASCII letters)
_LIKE_MZ_ML = r'\.[mM][zZ].*[mM][lL]\Z'


def load_files_table(path):
    """Reads the uniquemri export once and indexes its mz(X)ML rows by dataset."""
    global FILES_TABLE
    df = pd.read_csv(path, usecols=['usi', 'filepath', 'dataset'], dtype=str, keep_default_na=False, na_values=[])
    df = df[df['filepath'].str.contains(_LIKE_MZ_ML, regex=True, flags=re.DOTALL)]
    # the per-dataset query orders by usi (SQLite BINARY collation = code point order)
    df = df.sort_values('usi', kind='stable')
    FILES_TABLE = {dataset: group.reset_index(drop=True) for dataset, group in df.groupby('dataset', sort=False)}
    print(f"Dataset cache table {path}: {len(df)} mz(X)ML files in {len(FILES_TABLE)} datasets")


def _fetch_dataset_files(dataset):
    """mz(X)ML files of a dataset from the dataset cache, or None if they could not be fetched."""
    if FILES_TABLE is not None and dataset in FILES_TABLE:
        MATCH_STATS['datasets_from_files_table'] += 1
        return FILES_TABLE[dataset]
    # Not in the table (no mz(X)ML files when the export was made, or no table): ask the cache.
    # _stream=on returns all rows (without it at most max_returned_rows = 50000); the
    # timeout is short, as the cache answers within its 60 s SQL time limit or not at all
    url_f = "{}{}&filepath__endswith=%25.mz%25ML&_size=max&_stream=on".format(ccms_peak_link, dataset)
    print("Fetching {}".format(url_f))
    MATCH_STATS['datasets_queried'] += 1

    retry_attempts = 3
    for attempt in range(retry_attempts):
        print(f"Checking dataset {dataset} by going to the dataset cache")
        try:
            dataset_files_response = requests.get(url_f, timeout=(15, 120))
            dataset_files_response.raise_for_status()
            attempt_df = pd.read_csv(StringIO(dataset_files_response.text))
            if 'filepath' in attempt_df.columns:
                return attempt_df
        except (requests.RequestException, pd.errors.ParserError, pd.errors.EmptyDataError) as e:
            print(f"Attempt {attempt + 1} error: {e}")
        print(f"Attempt {attempt + 1} failed. Retrying in 5 seconds...")
        sleep(5)
    MATCH_STATS['datasets_query_failed'] += 1
    return None


class DatasetFiles:
    """File lookups of one dataset (by name, name without .mzML/.mzXML, lower-case name)."""

    def __init__(self, ccms_df):
        ccms_df = ccms_df.reset_index(drop=True)
        ccms_df["query_path"] = ccms_df["filepath"].apply(lambda x: os.path.basename(x))
        ccms_df["_dir_parts"] = ccms_df["filepath"].apply(lambda p: tuple(Path(p).parts[:-1]))
        self.df = ccms_df
        # plain lists (row i = position i) for the per-row candidate scoring in _select
        self.filepaths = ccms_df["filepath"].tolist()
        self.dir_parts = ccms_df["_dir_parts"].tolist()
        self.by_name = collections.defaultdict(list)
        self.by_stem = collections.defaultdict(list)
        self.by_lower = collections.defaultdict(list)
        self.by_lower_stem = collections.defaultdict(list)
        for i, name in zip(ccms_df.index, ccms_df["query_path"]):
            name = str(name)
            self.by_name[name].append(i)
            self.by_lower[name.lower()].append(i)
            stem, ext = os.path.splitext(name)
            if ext.lower() in _MZ_EXTENSIONS:
                self.by_stem[stem].append(i)
                self.by_lower_stem[stem.lower()].append(i)

    def candidates(self, meta_basename):
        """(rule, candidate rows) for the first rule that finds any file."""
        meta_basename2 = meta_basename[:-3] + "ML" if len(meta_basename) >= 3 else meta_basename
        exact = self.by_name.get(meta_basename, []) + \
            ([] if meta_basename2 == meta_basename else self.by_name.get(meta_basename2, []))
        if exact:
            return 'exact', exact
        stem, ext = os.path.splitext(meta_basename)
        if ext.lower() in _MZ_EXTENSIONS:
            same_stem = self.by_stem.get(stem, [])
            if same_stem:
                return 'stem_ext', same_stem
        lower = self.by_lower.get(meta_basename.lower(), [])
        if ext.lower() in _MZ_EXTENSIONS:
            lower = list(dict.fromkeys(lower + self.by_lower_stem.get(stem.lower(), [])))
        if lower:
            return 'case', lower
        return None, []


def _select(files, rows, meta_path):
    """The best candidate file path, or None (with a reason) if several fit equally well."""
    meta_dir_parts = tuple(Path(meta_path).parts[:-1])  # () if no dirs present

    # Score candidates by longest common directory suffix with the metadata path
    # (all score 0 when the metadata has no directory info); keep the best ones
    scores = [longest_common_dir_suffix_score(meta_dir_parts, files.dir_parts[i]) for i in rows]
    best_score = max(scores)
    found_file_paths = list(dict.fromkeys(files.filepaths[i] for i, score in zip(rows, scores) if score == best_score))

    # The same file can be present in several top-level folders (e.g. the
    # converted ccms_peak copy and the uploaded peak/raw copy); prefer them in
    # this order. Folders are compared as whole path segments.
    selected_paths = found_file_paths
    for preferred_dir in PREFERRED_DIRECTORIES:
        in_dir = [path for path in found_file_paths if preferred_dir in Path(path).parts[:-1]]
        if in_dir:
            selected_paths = in_dir
            break

    if len(selected_paths) > 1:
        # Several distinct files fit equally well (e.g. pos/s1.mzML and
        # neg/s1.mzML for metadata 's1.mzML'); do not guess
        return None, selected_paths
    return selected_paths[0], selected_paths


def _match_filenames_and_add_usi(dataset_metadata_df, file_accession=None, files_cache=None, source=''):
    """Matches each metadata row to a file of its own dataset (ATTRIBUTE_DatasetAccession;
    files can describe several datasets). Rows without a valid accession are matched in the
    dataset of the metadata file. Returns the matched rows with 'filename' set to the file
    path and a USI, or None if nothing could be matched."""
    files_cache = {} if files_cache is None else files_cache
    accessions = dataset_metadata_df['ATTRIBUTE_DatasetAccession'].astype(str).str.strip()
    valid = accessions.str.match(_ACCESSION)
    fallback = file_accession if file_accession and _ACCESSION.match(file_accession) else \
        (accessions[valid].iloc[0] if valid.any() else None)
    if (~valid).any():
        print(f"{(~valid).sum()} rows without a valid dataset accession; matching them in {fallback}")
        MATCH_STATS['rows_invalid_accession'] += int((~valid).sum())
        if fallback is None:
            dataset_metadata_df = dataset_metadata_df[valid]
            accessions = accessions[valid]
        else:
            accessions = accessions.where(valid, fallback)
            MATCH_STATS['rows_invalid_accession_set_to_file_dataset'] += int((~valid).sum())
    dataset_metadata_df = dataset_metadata_df.copy()
    dataset_metadata_df['ATTRIBUTE_DatasetAccession'] = accessions

    output_row_list = []
    for dataset, group in dataset_metadata_df.groupby('ATTRIBUTE_DatasetAccession', sort=False):
        print(f"Checking dataset {dataset} ({len(group)} rows)")
        if dataset not in files_cache:
            ccms_df = _fetch_dataset_files(dataset)
            if ccms_df is None:
                print(f"Error: could not fetch the file list for dataset {dataset}, skipping it")
                files_cache[dataset] = None
            elif len(ccms_df) == 0:
                print("Error: the size of the input is too small, length:", len(ccms_df))
                files_cache[dataset] = None
            else:
                print("Received a dataframe with {} rows.".format(len(ccms_df)))
                files_cache[dataset] = DatasetFiles(ccms_df)
        files = files_cache[dataset]
        if files is None:
            MATCH_STATS['rows_dataset_without_files'] += len(group)
            MATCH_REPORT.extend((source, dataset, f, 'dataset_without_files', '') for f in group['filename'].astype(str))
            continue

        matched = []
        for metadata_row in group.to_dict('records'):
            meta_path = str(metadata_row["filename"])
            meta_basename = os.path.basename(meta_path)
            rule, rows = files.candidates(meta_basename)
            if not rows:
                MATCH_STATS['rows_not_found'] += 1
                MATCH_REPORT.append((source, dataset, meta_path, 'not_found', ''))
                continue
            selected_path, options = _select(files, rows, meta_path)
            if selected_path is None:
                print(f"Ambiguous match for {meta_path}, skipping: {options}")
                MATCH_STATS[f'rows_ambiguous_{rule}'] += 1
                AMBIGUOUS.append((dataset, meta_path, rule, ' ; '.join(options)))
                MATCH_REPORT.append((source, dataset, meta_path, f'ambiguous_{rule}', ' ; '.join(options)))
                continue
            if rule == 'exact' and os.path.basename(selected_path) != meta_basename:
                rule = 'exact_mzxml_to_converted_mzml'
            matched.append((metadata_row, selected_path, rule, os.path.basename(selected_path) == meta_basename, meta_path))

        # A file named exactly by one row is not also given to rows that only reach it
        # through another extension (x.mzXML and x.mzML listed for the same sample: the
        # x.mzXML row would otherwise also get the converted ccms_peak/x.mzML)
        exact_targets = {path for _, path, _, is_exact, _ in matched if is_exact}
        for metadata_row, selected_path, rule, is_exact, meta_path in matched:
            if not is_exact and selected_path in exact_targets:
                MATCH_STATS['rows_twin_of_exactly_named_file'] += 1
                MATCH_REPORT.append((source, dataset, meta_path, 'twin_of_exactly_named_file', selected_path))
                continue
            MATCH_STATS[f'rows_matched_{rule}'] += 1
            MATCH_REPORT.append((source, dataset, meta_path, f'matched_{rule}', selected_path))
            metadata_row["filename"] = "f." + selected_path
            metadata_row["USI"] = _make_usi_from_filename(
                metadata_row["filename"],
                metadata_row["ATTRIBUTE_DatasetAccession"]
            )
            output_row_list.append(metadata_row)

    if not output_row_list:
        return None
    return pd.DataFrame(output_row_list)

def main():
    # parsing args
    parser = argparse.ArgumentParser(description='GNPS Name Matcher')
    parser.add_argument('passed_file_names', help='Input TSV file')
    parser.add_argument('metadata_folder')
    parser.add_argument('output_filename')
    parser.add_argument('path_allowed_terms_json')
    parser.add_argument('--match_report', default=None, help='optional TSV with the matching result of every metadata row')
    parser.add_argument('--files_table', default=None,
                        help='local copy of https://datasetcache.gnps2.org/dataset/uniquemri; without it (or if it '
                             'cannot be read, e.g. empty after a failed download) every dataset is queried from the dataset cache')

    args = parser.parse_args()
    if args.files_table:
        try:
            load_files_table(args.files_table)
        except Exception as e:
            print(f"WARNING: could not read the dataset cache table {args.files_table} ({type(e).__name__}: {e}); "
                  f"querying the dataset cache per dataset")
     
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
    # file lists of the datasets, fetched once (a dataset can be described by several files)
    files_cache = {}

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

        # Matching the metadata, each row within its own dataset
        file_accession = re.match(r'([A-Z]+\d+)', os.path.basename(str(file)))
        MATCH_STATS['rows_in'] += len(dataset_metadata_df)
        enriched_metadata_df = _match_filenames_and_add_usi(
            dataset_metadata_df, file_accession.group(1) if file_accession else None, files_cache,
            source=os.path.basename(str(file)))
        if enriched_metadata_df is not None:
            enriched_metadata_df['MetadataFileDataset'] = file_accession.group(1) if file_accession else ''

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

    print("Matching summary (rows):")
    for key, value in sorted(MATCH_STATS.items()):
        print(f"  {key}: {value}")
    for dataset, meta_path, rule, options in AMBIGUOUS:
        print(f"  ambiguous ({rule}) {dataset} {meta_path}: {options}")

    if args.match_report:
        pd.DataFrame(MATCH_REPORT, columns=['source_file', 'dataset', 'filename', 'result', 'file_path']).to_csv(
            args.match_report, sep='\t', index=False)

    # Save the DataFrame to a TSV file without column names
    merged_metadata_df.to_csv(args.output_filename, sep='\t', index=False)


if __name__ == '__main__':
    main()