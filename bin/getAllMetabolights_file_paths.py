import requests
import pandas as pd
from tqdm import tqdm
import requests
import pandas as pd
import argparse 
import time
import os 

# (connect, read) timeouts in seconds; the recursive listing of large studies is slow
REQUEST_TIMEOUT = (10, 300)
# Status codes worth retrying: throttling and transient server errors. 504 is the
# server-side timeout of a (too) broad search pattern and is handled separately.
RETRY_STATUS = {429, 500, 502, 503}

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
      response = requests.get(url, timeout=REQUEST_TIMEOUT)
      print(f"Requested: {url}")
      if response.status_code in expected_codes:
        data = response.json()
        print("Request successful!")
        return data
      else:
        print(f"Error: Unexpected status code {response.status_code}")
    except Exception as e:
      print(f"Error requesting data: {e}")
    time.sleep(10)
  print(f"All retries failed for {url}.")
  return None


def get_all_files(study_id, headers):
    base_url = f'https://www.ebi.ac.uk:443/metabolights/ws/studies/{study_id}/public-data-files'
    # Two complementary broad patterns are BOTH required to list a study completely:
    #   - 'FILES/*'     matches files that sit directly under FILES/
    #   - 'FILES/**/*.*' matches files nested inside subdirectories of FILES/
    # The recursive '**' glob does NOT match top-level files, so studies whose data
    # lives directly in FILES/ return nothing from '**' alone (this silently dropped
    # hundreds of studies). We union both and dedup by name.
    broad_patterns = ['FILES/*', 'FILES/**/*.*']
    # Extension-specific recursive fallbacks, only used when the recursive broad
    # pattern times out (504) on very large studies.
    fallback_patterns = [
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
    seen = set()

    def _add(files):
        for f in files:
            name = f['name']
            if name not in seen:
                seen.add(name)
                all_files.append(f)

    def _query(pattern, retries=4):
        url = f"{base_url}?search_pattern={pattern}&file_match=true&folder_match=true"
        for attempt in range(1, retries + 1):
            try:
                response = requests.get(url, headers=headers, timeout=REQUEST_TIMEOUT)
                if response.status_code not in RETRY_STATUS or attempt == retries:
                    response.raise_for_status()  # Raises HTTPError for 4XX/5XX responses
                    return response.json()['files']
                print(f"[{study_id}] {pattern}: status {response.status_code}, retrying (attempt {attempt})")
            except (requests.ConnectionError, requests.Timeout) as err:
                if attempt == retries:
                    raise
                print(f"[{study_id}] {pattern}: {err}, retrying (attempt {attempt})")
            time.sleep(min(5 * 2 ** attempt, 120))

    # A study is incomplete if any pattern could not be listed; it is reported so
    # the caller can tell an empty study from a failed one
    failed = False
    timed_out = []
    for pattern in broad_patterns:
        try:
            _add(_query(pattern))
        except requests.exceptions.HTTPError as http_err:
            status = getattr(http_err.response, 'status_code', None)
            if status == 504:  # the broad scans can time out on very large studies
                print(f"Timeout encountered in {study_id} with pattern: {pattern}. Trying more specific patterns.")
                timed_out.append(pattern)
            else:
                print(f"HTTP error occurred in {study_id} : {http_err} - Status code: {status}")
                failed = True
        except requests.RequestException as err:
            print(f"Request failed for {study_id}: {err}")
            failed = True
        except (ValueError, KeyError) as json_err:
            print(f"JSON decoding failed in {study_id}: {json_err}")
            failed = True

    # Recursive fallbacks do not match top-level files, so a timed out 'FILES/*'
    # falls back to the same extensions directly under FILES/
    retry_patterns = []
    if 'FILES/**/*.*' in timed_out:
        retry_patterns += fallback_patterns
    if 'FILES/*' in timed_out:
        retry_patterns += [p.replace('FILES/**/', 'FILES/') for p in fallback_patterns]

    if retry_patterns:
        for pattern in retry_patterns:
            try:
                _add(_query(pattern))
            except requests.exceptions.HTTPError as http_err:
                status = getattr(http_err.response, 'status_code', None)
                print(f"HTTP error occurred in {study_id} : {http_err} - Status code: {status}")
                failed = True
            except requests.RequestException as err:
                print(f"Request failed for {study_id}: {err}")
                failed = True
            except (ValueError, KeyError) as json_err:
                print(f"JSON decoding failed in {study_id}: {json_err}")
                failed = True

    if not all_files:
        print(f"No files found for study {study_id}.")
    return {'files': all_files, 'failed': failed}

            
       
def create_usi(row):
    return f"mzspec:{row['study_id']}:{row['file_path']}"

def process_filename(filename):
    path_parts = filename.split('/')
    for part in path_parts:
        if part.lower().endswith('.d'):
            return '/'.join(path_parts[:path_parts.index(part)+1])
    return filename
            
def _get_existing_datasets(path_to_file):
    try:
        existing_datasets = pd.read_csv(path_to_file, sep="\t")
        existing_datasets = set(existing_datasets['datasets'])
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

    public_metabolights_studies = safe_api_request('https://www.ebi.ac.uk:443/metabolights/ws/studies/technology', retries = 5)
    if not public_metabolights_studies:
        raise RuntimeError("Could not retrieve the list of MetaboLights studies.")


    # Initialize variables and headers with your API token (public files need none)
    headers = {'user_token': args.user_token} if args.user_token != "none" else {}
    data = []
    failed_studies = []

    # Loop over each study and gather files
    for study in tqdm(public_metabolights_studies, desc="Processing studies"):
        study_id = study['accession']

        if study_id in existing_datasets:
            print("Skipping", study_id, "Already indexed")
            continue

        files_list = get_all_files(study_id, headers)  
        if files_list['failed']:
            failed_studies.append(study_id)
        
        # Extract each file and append to the list with the study_id
        for file in files_list['files']:
            data.append({'study_id': study_id, 'file_path': file['name']})

    if failed_studies:
        print(f"File listing incomplete for {len(failed_studies)} of {len(public_metabolights_studies)} studies:")
        print(" ".join(failed_studies))

    if not data:
        raise RuntimeError("No MetaboLights files were listed for any study.")

    # Convert the list to a DataFrame
    files_df = pd.DataFrame(data)

    #add  USI
    files_df['file_path'] = files_df['file_path'].apply(process_filename)
    files_df = files_df.drop_duplicates(keep='first')
    files_df['USI'] = files_df.apply(create_usi, axis=1)

    #export tsv file
    files_df.to_csv(args.output_filename, index = False, sep = "\t")






