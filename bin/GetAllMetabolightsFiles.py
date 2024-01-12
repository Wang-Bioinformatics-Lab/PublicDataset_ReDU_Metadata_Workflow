import os
import requests
import pandas as pd
from bs4 import BeautifulSoup
import argparse
import json
import time
import numpy as np
from tqdm import tqdm



def safe_api_request(url, retries=3, expected_codes={200}):

  """
  Safely requests JSON data from an API and checks for errors.

  Args:
    url: The URL of the API endpoint.
    retries: The number of retries to attempt in case of errors.
    expected_codes: A set of expected HTTP status codes indicating success.

  Returns:
    A dictionary containing the JSON data if successful,
    or None if all retries fail.
  """
  for _ in range(retries):
    try:
      response = requests.get(url)
      print(f"Requested: {url}")
      if response.status_code in expected_codes:
        data = response.json()
        print("Request successful!")
        return data
      else:
        print(f"Error: Unexpected status code {response.status_code}")
        time.sleep(10)
    except Exception as e:
      print(f"Error requesting data: {e}")
  print(f"All retries failed for {url}.")
  return None


def GetMetabolightsFilePaths(study_id):
    """
    Converts Metabolights study data to a REDU table format.

    Args:
    study_id: The ID of the Metabolights study to convert.

    Returns:
    A DataFrame in the REDU table format with processed and aligned data from the Metabolights study,
    or an empty DataFrame if no applicable data is found.
    """
    study_url = "https://www.ebi.ac.uk:443/metabolights/ws/studies/public/study/" + study_id
    study_details = safe_api_request(study_url)

    with open('filename.json', 'w') as file:
        json.dump(study_details, file, indent=4)

    study_assays = study_details['content']['assays']

    df_assays = pd.DataFrame()
    ms_study_assays = []

    #get study assays if they are MS. There can be multiple assay tables in the same study
    for index, assay_json in enumerate(study_assays):
        if assay_json['technology'] == 'mass spectrometry':  

            # Extract headers in the correct order
            headers = [None] * len(assay_json['assayTable']['fields'])
            for key, value in assay_json['assayTable']['fields'].items():
                headers[value['index']] = value['header']        

            df = pd.DataFrame(assay_json['assayTable']['data'])
            df.columns = headers

            ms_study_assays.append(df)

    if len(ms_study_assays) > 0:

        all_columns = set()
        for df in ms_study_assays:
            all_columns.update(df.columns)

        aligned_dfs = []
        for df in ms_study_assays:
            # Add missing columns with NaN values
            for col in all_columns:
                if col not in df.columns:
                    df[col] = np.nan

        # Reorder columns and add to the aligned list
        aligned_dfs.append(df[list(all_columns)])

        df_assays = pd.concat(aligned_dfs, ignore_index=True)

        # Duplicate rows if we have mzml AND raw files
        df_study_raw = pd.DataFrame()
        if 'Raw Spectral Data File' in df_assays.columns:
            df_study_raw = df_assays[df_assays['Raw Spectral Data File'].str.contains('\.', regex=True, na=False)].copy()
            df_study_raw['filepath'] = df_study_raw['Raw Spectral Data File']
            df_study_raw.drop(columns=['Raw Spectral Data File', 'Derived Spectral Data File'], inplace=True)
        
        df_study_mzml = pd.DataFrame()
        if 'Derived Spectral Data File' in df_assays.columns:
            df_study_mzml = df_assays[df_assays['Derived Spectral Data File'].str.contains('\.', regex=True, na=False)].copy()
            df_study_mzml['filepath'] = df_study_mzml['Derived Spectral Data File']
            df_study_mzml.drop(columns=['Raw Spectral Data File', 'Derived Spectral Data File'], inplace=True)

        if len(df_study_raw) > 0 and len(df_study_mzml) > 0:
            df_study = pd.concat([df_study_mzml, df_study_raw], ignore_index=True)
        elif len(df_study_raw) > 0:
            df_study = df_study_raw
        elif len(df_study_mzml) > 0:
            df_study = df_study_mzml
        elif len(df_study_raw) == 0 and len(df_study_mzml) == 0:
            return pd.DataFrame()
        

        df_study['filename'] = df_study['filepath'].apply(lambda x: os.path.basename(x))


        # List of allowed extensions
        allowed_extensions = [".mzml", ".mzxml", ".cdf", ".raw", ".wiff", ".d"]

        # Use a tuple of allowed extensions for the endswith method
        df_study = df_study[df_study['filepath'].str.lower().str.endswith(tuple(allowed_extensions))]
        
        df_study['Metabolights_ID'] = study_id

        df_study = df_study[['Metabolights_ID', 'filepath']]

        return df_study


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Give an Metabolights study ID and a tsv with all file paths.')
    parser.add_argument("--study_id", type=str, help='An Metabolights study ID such as "MTBLS1015". If "ALL" all studys are requested.', required=True)
    
    
    args = parser.parse_args()

    if args.study_id == 'ALL':
        public_metabolights_studies = safe_api_request('https://www.ebi.ac.uk:443/metabolights/ws/studies', retries = 1)
        public_metabolights_studies = public_metabolights_studies['content']
    else:
        public_metabolights_studies = [args.study_id]


    REDU_dataframes = []
    redu_table_single = pd.DataFrame()
    for study_id in tqdm(public_metabolights_studies):
        try:
            redu_table_single = GetMetabolightsFilePaths(study_id)
        except Exception as e:
            print(f"An error occurred with study_id {study_id}: {e}")
            continue
        if redu_table_single is not None and len(redu_table_single) > 0:
            REDU_dataframes.append(redu_table_single)

    if len(REDU_dataframes) > 0:
        redu_tables_all = pd.concat(REDU_dataframes, ignore_index=True)
        redu_tables_all.to_csv('MetabolightsFilePaths_' + args.study_id + '.tsv', sep='\t', index=False, header=True)
        print(f'Output has been saved to MetabolightsFilePaths_{args.study_id}.tsv!')
    else:
        print('nothing to return!')




