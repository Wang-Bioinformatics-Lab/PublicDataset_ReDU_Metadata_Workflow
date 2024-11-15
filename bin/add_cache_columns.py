import pandas as pd
import argparse
import requests
from io import StringIO
from time import sleep



if __name__ == "__main__":

    # parsing args
    parser = argparse.ArgumentParser(description='GNPS Name Matcher')
    parser.add_argument('merged_metadata_path', help='Input TSV file')
    parser.add_argument('output_filename')

    args = parser.parse_args()


    url_cache_datasette = "https://datasetcache.gnps2.org/dataset/uniquemri"
    print("Fetching {}".format(url_cache_datasette))

    # Request URL and retry 3 times after waiting 5 seconds if ccms_df does not have a column named 'filepath'
    retry_attempts = 3
    for attempt in range(retry_attempts):
        print(f"Attempt {attempt + 1} to fetch the dataset cache")
        dataset_files_response = requests.get(url_cache_datasette)

        csvStringIO = StringIO(dataset_files_response.text)
        cache_df = pd.read_csv(csvStringIO)

        if 'usi' in cache_df.columns:
            break
        else:
            print(f"Attempt {attempt + 1} failed. Retrying in 5 seconds...")
            sleep(5)


    if len(cache_df) == 0:
        print("Cache dataframe is empty. Exiting...")
        exit(1)

    
    # subset the cache dataframe to only the columns we need
    cache_df = cache_df[['usi', 'classification', 'spectra_ms2']]

    #set name usi to USI
    cache_df.rename(columns={'usi':'USI'}, inplace=True)
    cache_df.rename(columns={'spectra_ms2':'MS2spectra_count'}, inplace=True)


    # read the metadata file
    metadata_df = pd.read_csv(args.merged_metadata_path, sep='\t')


    # merge the metadata with the cache dataframe keeping all metadata rows
    metadata_df = metadata_df.merge(cache_df, on='USI', how='left')


    # fill NA value in classification with Unclassified
    metadata_df['classification'].fillna('Unclassified', inplace=True)


    # write the metadata file with the new columns
    metadata_df.to_csv(args.output_filename, sep='\t', index=False)

