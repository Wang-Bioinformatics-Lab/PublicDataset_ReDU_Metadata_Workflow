import requests
import pandas as pd
from tqdm import tqdm
import requests
import pandas as pd
import argparse 
import time
import os 

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


def get_all_files(study_id, headers):
    base_url = f'https://www.ebi.ac.uk:443/metabolights/ws/studies/{study_id}/data-files'
    search_patterns = [
        'FILES/**/*.*',
        'FILES/**/*.zip',
        'FILES/**/*.raw',
        'FILES/**/*.RAW',
        'FILES/**/*.d',
        'FILES/**/*.D',
        'FILES/**/*.lcd',
        'FILES/**/*.mzML',
        'FILES/**/*.mzXML',
        'FILES/**/*.cdf',
        'FILES/**/*.CDF',
        'FILES/**/*.wiff',
        'FILES/**/*.wiff.scan'
    ]
    
    all_files = []
    initial_try = True

    for pattern in search_patterns:
        url = f"{base_url}?search_pattern={pattern}&file_match=true&folder_match=true"
        try:
            response = requests.get(url, headers=headers)
            response.raise_for_status()  # Raises HTTPError for 4XX/5XX responses
            files_data = response.json()
            if files_data['files']:
                all_files.extend(files_data['files'])
                if initial_try:  # If initial broad pattern works, no need to try more specific patterns
                    break
        except requests.exceptions.HTTPError as http_err:
            if response.status_code == 504 and initial_try:  # Gateway Timeout only matters on first try
                print(f"Timeout encountered in {study_id} with pattern: {pattern}. Trying more specific patterns.")
                initial_try = False  # Subsequent attempts will not break after first success
            else:
                print(f"HTTP error occurred in {study_id} : {http_err} - Status code: {response.status_code}")
        except requests.RequestException as err:
            print(f"Request failed: {err}")
        except ValueError as json_err:
            print(f"JSON decoding failed in {study_id}: {json_err}")
            print("Response content:", response.text)  # Print the response text that was not JSON

    if not all_files:
        print(f"No files found for study {study_id}.")
    return {'files': all_files}

            
       
def create_usi(row):
    return f"mzspec:{row['study_id']}:{row['file_path']}"

def process_filename(filename):
    path_parts = filename.split('/')
    for part in path_parts:
        print(part)
        if part.endswith('.d'):
            return '/'.join(path_parts[:path_parts.index(part)+1])
    return filename
            
def _get_existing_datasets(path_to_file):
    try:
        existing_datasets = pd.read_csv(path_to_file, sep="\t")
        existing_datasets = set(existing_datasets['dataset'])
        return existing_datasets
    except:
        return set()

if __name__ == "__main__":
            
    parser = argparse.ArgumentParser(description='Get all Metabolights file paths')
    parser.add_argument("--output_filename", type=str, help="tsv file name for output", default="none")
    parser.add_argument("--user_token", type=str, help="user token you can get from metabolights account", default="none")
    parser.add_argument("--existing_datasets", type=str, help="path to a file of datasets already indexed", default="none")
    
    args = parser.parse_args()
        
    existing_datasets = _get_existing_datasets(args.existing_datasets)

    public_metabolights_studies = safe_api_request('https://www.ebi.ac.uk:443/metabolights/ws/studies/technology', retries = 1)


    # Initialize variables and headers with your API token
    headers = {'user_token': args.user_token}
    data = []

    # Loop over each study and gather files
    for study in tqdm(public_metabolights_studies, desc="Processing studies"):
        study_id = study['accession']

        if study_id in existing_datasets:
            print("Skipping", study_id, "Already indexed")
            continue

        files_list = get_all_files(study_id, headers)  
        
        # Extract each file and append to the list with the study_id
        for file in files_list['files']:
            data.append({'study_id': study_id, 'file_path': file['name']})

    # Convert the list to a DataFrame
    files_df = pd.DataFrame(data)

    #add  USI
    files_df['file_path'] = files_df['file_path'].apply(process_filename)
    files_df = files_df.drop_duplicates(keep='first')
    files_df['USI'] = files_df.apply(create_usi, axis=1)

    #export tsv file
    files_df.to_csv(args.output_filename, index = False, sep = "\t")






