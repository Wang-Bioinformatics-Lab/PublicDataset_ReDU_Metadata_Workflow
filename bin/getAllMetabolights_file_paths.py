import requests
import pandas as pd
from tqdm import tqdm
import requests
import pandas as pd
import argparse 
import time

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


def get_files_recursive(study_id, directory, headers, files_list, depth=0, max_depth=4):
    if depth > max_depth or 'nmr' in directory.lower():
        return

    url = f'https://www.ebi.ac.uk/metabolights/ws/studies/{study_id}/fileslist?directory_name={directory}'
    response = requests.get(url, headers=headers)
    filelist_archive = response.json()
    
    for file_info in filelist_archive.get('files', []):
        if '.' in file_info['file']:  # Check if it's a file
            file_path = file_info.get('path', '')
            files_list.append({'study_id': study_id, 'file_path': file_path})
    
    for subdirectory_info in filelist_archive.get('directories', []):
        subdirectory_name = subdirectory_info['directory']
        if 'nmr' not in subdirectory_name.lower():
            new_directory = f"{directory}/{subdirectory_name}" if directory else subdirectory_name
            get_files_recursive(study_id, new_directory, headers, files_list, depth=depth+1, max_depth=max_depth)




def get_all_files(study_id, headers):
    try:
        url = f'https://www.ebi.ac.uk:443/metabolights/ws/studies/{study_id}/data-files?search_pattern=FILES%2F**%2F*.*&file_match=true&folder_match=true'
        response = requests.get(url, headers=headers)
        response.raise_for_status()  # Raises an HTTPError for bad responses (4XX, 5XX)
        return response.json()  # Attempt to return the JSON
    except requests.HTTPError as http_err:
        print(f"HTTP error occurred: {http_err} - Status code: {response.status_code}")
    except requests.RequestException as err:
        print(f"Request failed: {err}")
    except ValueError as json_err:  # Includes JSONDecodeError
        print(f"JSON decoding failed: {json_err}")
        print("Response content:", response.text)  # Print the response text that was not JSON
    return {'files': []}





    return filelist_archive

            
       
def create_usi(row):
    return f"mzspec:{row['study_id']}:{row['file_path']}"
     
            
if __name__ == "__main__":
            
    parser = argparse.ArgumentParser(description='Get all Metabolights file paths')
    parser.add_argument("--output_filename", type=str, help="tsv file name for output", default="none")
    parser.add_argument("--user_token", type=str, help="user token you can get from metabolights account", default="none")
    
    
    args = parser.parse_args()
        

    public_metabolights_studies = safe_api_request('https://www.ebi.ac.uk:443/metabolights/ws/studies/technology', retries = 1)


    # Initialize variables and headers with your API token
    headers = {'user_token': args.user_token}
    data = []

    # Loop over each study and gather files
    for study in tqdm(public_metabolights_studies, desc="Processing studies"):
        study_id = study['accession']
        files_list = get_all_files(study_id, headers)  
        
        # Extract each file and append to the list with the study_id
        for file in files_list['files']:
            data.append({'study_id': study_id, 'file_path': file['name']})

    # Convert the list to a DataFrame
    files_df = pd.DataFrame(data)

    #add  USI
    files_df['USI'] = files_df.apply(create_usi, axis=1)

    #export tsv file
    files_df.to_csv(args.output_filename, index = False, sep = "\t")






