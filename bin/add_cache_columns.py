import pandas as pd
import argparse
import requests
from time import sleep
import numpy as np
import os
import csv
import sys
import traceback


# only these columns of the export are used
CACHE_COLUMNS = ['usi', 'classification', 'spectra_ms2']
DOWNLOAD_PATH = 'uniquemri.part.csv'

# Parquet copy of the final table (read by the web app; the task fails without it); row groups of 131072
# rows let readers that filter (e.g. on a dataset) skip most of the file
PARQUET_INT_COLUMNS = ['MS2spectra_count']
PARQUET_ROW_GROUP_SIZE = 131072


def write_parquet_copy(tsv_path, parquet_path, expected_rows=None):
    """Parquet copy of the TSV just written, read back from that file so both hold the same values:
    every column is the TSV's text (nothing becomes null: 'missing value', '', 'NA' stay as they
    are), except MS2spectra_count (int64). Written under a temporary name and then renamed, so a
    failed write never leaves a partial file. The output only depends on the TSV's bytes."""
    import pyarrow as pa
    import pyarrow.csv as pa_csv
    import pyarrow.parquet as pq

    with open(tsv_path, encoding='utf-8', newline='') as handle:
        columns = next(csv.reader(handle, delimiter='\t'))
    column_types = {c: (pa.int64() if c in PARQUET_INT_COLUMNS else pa.string()) for c in columns}
    table = pa_csv.read_csv(
        tsv_path,
        parse_options=pa_csv.ParseOptions(delimiter='\t', newlines_in_values=True),
        convert_options=pa_csv.ConvertOptions(column_types=column_types, null_values=[],
                                              strings_can_be_null=False, quoted_strings_can_be_null=False))
    if table.column_names != columns or (expected_rows is not None and table.num_rows != expected_rows):
        raise ValueError("the TSV read back does not have the rows and columns written")
    # (lets a reader check that this file belongs to the TSV next to it)
    table = table.replace_schema_metadata({'redu_tsv_bytes': str(os.path.getsize(tsv_path))})

    temporary_path = parquet_path + '.part'
    try:
        pq.write_table(table, temporary_path, compression='zstd', version='2.6',
                       row_group_size=PARQUET_ROW_GROUP_SIZE)
        os.replace(temporary_path, parquet_path)
    finally:
        if os.path.exists(temporary_path):
            os.remove(temporary_path)


def fetch_cache_table(url, timeout=1800):
    """Downloads the dataset cache export to a file and reads its used columns. Same values
    as reading the whole response text (the undecodable bytes replaced, as requests' .text
    does), but with ~2 GB instead of ~10 GB of memory."""
    try:
        with requests.get(url, timeout=timeout, stream=True) as response:
            response.raise_for_status()
            # (with a Content-Encoding the length is that of the encoded body)
            expected = None if response.headers.get('Content-Encoding') else response.headers.get('Content-Length')
            with open(DOWNLOAD_PATH, 'wb') as handle:
                for chunk in response.iter_content(chunk_size=1 << 20):
                    handle.write(chunk)
        size = os.path.getsize(DOWNLOAD_PATH)
        if expected and int(expected) != size:
            raise IOError(f"incomplete download: {size} of {expected} bytes")
        return read_cache_table(DOWNLOAD_PATH)
    finally:
        if os.path.exists(DOWNLOAD_PATH):
            os.remove(DOWNLOAD_PATH)


def read_cache_table(path):
    """The used columns of a dataset cache export file."""
    return pd.read_csv(path, usecols=lambda c: c in CACHE_COLUMNS,
                       encoding='utf-8', encoding_errors='replace')


def usable_cache(cache_df):
    return set(CACHE_COLUMNS) <= set(cache_df.columns) and len(cache_df) > 0


if __name__ == "__main__":

    # parsing args
    parser = argparse.ArgumentParser(description='GNPS Name Matcher')
    parser.add_argument('--parquet_output', help='also write the table as Parquet to this path')
    parser.add_argument('merged_metadata_path', help='Input TSV file')
    parser.add_argument('output_filename')
    parser.add_argument('--files_table', default=None,
                        help='the dataset cache export fetched earlier in this run (fetchDatasetCache; empty if '
                             'that download failed), used when the dataset cache cannot be fetched now')

    args = parser.parse_args()


    url_cache_datasette = "https://datasetcache.gnps2.org/dataset/uniquemri"
    print("Fetching {}".format(url_cache_datasette))

    # Request URL and retry 3 times after waiting 5 seconds if ccms_df does not have a column named 'filepath'
    retry_attempts = 3
    cache_df = pd.DataFrame()
    for attempt in range(retry_attempts):
        print(f"Attempt {attempt + 1} to fetch the dataset cache")
        try:
            cache_df = fetch_cache_table(url_cache_datasette)
        except (requests.RequestException, IOError, ValueError, pd.errors.ParserError, pd.errors.EmptyDataError) as e:
            print(f"Attempt {attempt + 1} failed: {e}")
            cache_df = pd.DataFrame()

        if 'usi' in cache_df.columns:
            break
        else:
            print(f"Attempt {attempt + 1} failed. Retrying in 5 seconds...")
            sleep(5)


    # If the dataset cache cannot be fetched now, the export this run fetched a few hours earlier is used;
    # without either, rows keep the classification/MS2 counts they already carry (rows from earlier runs)
    # instead of blocking the whole update
    cache_available = usable_cache(cache_df)
    if not cache_available and args.files_table:
        try:
            cache_df = read_cache_table(args.files_table)
        except (IOError, ValueError) as e:  # (an empty file: pandas' EmptyDataError, a ValueError)
            print(f"The dataset cache copy of this run ({args.files_table}) is not usable: {e}")
            cache_df = pd.DataFrame()
        cache_available = usable_cache(cache_df)
        if cache_available:
            print(f"WARNING: could not fetch the dataset cache; using the copy fetched earlier in this run "
                  f"({args.files_table}, {len(cache_df)} files)")
    if not cache_available:
        cache_df = pd.DataFrame(columns=['usi', 'classification', 'spectra_ms2'])

    # subset the cache dataframe to only the columns we need
    cache_df = cache_df[['usi', 'classification', 'spectra_ms2']]

    #set name usi to USI
    cache_df.rename(columns={'usi':'USI'}, inplace=True)
    cache_df.rename(columns={'spectra_ms2':'MS2spectra_count'}, inplace=True)


    cache_df['MS2spectra_count'] = pd.to_numeric(cache_df['MS2spectra_count'], errors='coerce').fillna(-1).astype(int)

    # making nan or inf to -1 in the MS2spectra_count column
    cache_df['MS2spectra_count'] = cache_df['MS2spectra_count'].replace([np.inf, -np.inf], -1)
    # making nan to -1
    cache_df['MS2spectra_count'] = cache_df['MS2spectra_count'].fillna(-1)
    # casting to int
    cache_df['MS2spectra_count'] = cache_df['MS2spectra_count'].astype(int)
    


    # read the metadata file
    metadata_df = pd.read_csv(args.merged_metadata_path, sep='\t')


    if cache_available:
        if 'classification' in metadata_df.columns:
            metadata_df = metadata_df.drop(columns=['classification'])

        if 'MS2spectra_count' in metadata_df.columns:
            metadata_df = metadata_df.drop(columns=['MS2spectra_count'])

        # merge the metadata with the cache dataframe keeping all metadata rows
        metadata_df = metadata_df.merge(cache_df, on='USI', how='left')
    else:
        # A row without an earlier value carries 'missing value' (text): it gets what the merge
        # gives a row the cache does not know (-1 MS2 spectra, 'Unclassified'), and the two
        # columns go last, as after the merge
        for column in ['classification', 'MS2spectra_count']:
            if column in metadata_df.columns:
                values = metadata_df.pop(column).replace('missing value', np.nan)
            else:
                values = np.nan
            metadata_df[column] = values
        metadata_df['MS2spectra_count'] = pd.to_numeric(metadata_df['MS2spectra_count'], errors='coerce')
        no_count = int((~np.isfinite(metadata_df['MS2spectra_count'].astype(float))).sum())
        no_class = int(metadata_df['classification'].isna().sum())
        print(f"WARNING: could not fetch the dataset cache and no copy of it from this run is usable: of "
              f"{len(metadata_df)} rows, {no_count} get MS2spectra_count -1 and {no_class} classification "
              f"'Unclassified'; the others keep the values they carry from earlier runs")


    # fill NA value in classification with Unclassified
    metadata_df['classification'].fillna('Unclassified', inplace=True)

    # making nan or inf to -1 in the MS2spectra_count column
    metadata_df['MS2spectra_count'] = metadata_df['MS2spectra_count'].replace([np.inf, -np.inf], -1)
    # making nan to -1
    metadata_df['MS2spectra_count'] = metadata_df['MS2spectra_count'].fillna(-1)
    # casting to int
    metadata_df['MS2spectra_count'] = metadata_df['MS2spectra_count'].astype(int)


    # write the metadata file with the new columns
    metadata_df.to_csv(args.output_filename, sep='\t', index=False)

    if args.parquet_output:
        # (frees ~3 GB before the TSV is read back)
        rows = len(metadata_df)
        del metadata_df, cache_df
        try:
            write_parquet_copy(args.output_filename, args.parquet_output, expected_rows=rows)
        except Exception as e:
            # the web app reads the Parquet copy: without it the task fails (non-zero exit), so Nextflow
            # publishes neither file and the table of the previous run stays in place
            traceback.print_exc()
            print(f"ERROR: could not write the Parquet copy of {args.output_filename} to {args.parquet_output}: "
                  f"{e!r}", file=sys.stderr)
            sys.stdout.flush()
            sys.stderr.flush()
            # (os._exit, not sys.exit: after rejecting a CSV value pyarrow 10 can wait forever for its reader
            # threads when the interpreter exits, so the task would hang instead of failing)
            os._exit(1)
