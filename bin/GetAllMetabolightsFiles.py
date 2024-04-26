import os
import requests
import pandas as pd
from bs4 import BeautifulSoup
import argparse
import json
import time
import numpy as np
from tqdm import tqdm



def GetMetabolightsFilePaths(study_id):
    study_url = "https://www.ebi.ac.uk:443/metabolights/ws/studies/public/study/" + study_id

    response = requests.get(study_url)
    study_details = response.json()

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
        extensions = [".mzml", ".mzxml", ".cdf", ".raw", ".wiff", ".d"]
        raw_files = df_assays[df_assays['Raw Spectral Data File'].str.lower().str.endswith(tuple(extensions), na=False)]['Raw Spectral Data File'].tolist() if 'Raw Spectral Data File' in df_assays.columns else []
        mzml_files = df_assays[df_assays['Derived Spectral Data File'].str.lower().str.endswith(tuple(extensions), na=False)]['Derived Spectral Data File'].tolist() if 'Derived Spectral Data File' in df_assays.columns else []

        all_files = raw_files + mzml_files

        df_output = pd.DataFrame(all_files, columns=['filepath'])
        df_output['mtbls_id'] = study_id

        return df_output


def safe_api_request(url, retries=3, expected_codes={200}):

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




