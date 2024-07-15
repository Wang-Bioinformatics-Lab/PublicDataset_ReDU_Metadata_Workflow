import os
import pandas as pd
import requests
import argparse
import tqdm
from urllib.parse import urlparse, parse_qs


def clean_mac_path(path):
    if "__MACOSX/" in path and '/._' in path:
        # Remove the __MACOSX part and the ._ prefix from the filename
        return path.replace("__MACOSX/", "").replace("/._", "/")
    return path

# Define a function to extract the first folder with an extension in order to not get the individual files in e.g. agilent "files"
def extract_first_folder_with_extension(url):
    filepath = parse_qs(urlparse(url).query).get('F', [None])[0]
    parts = filepath.split("/")

    # Check for extensions in each part
    for i, part in enumerate(parts):
        if "." in part:
            # If an extension is found, return the path up to and including that part
            return "/".join(parts[:i + 1])

    # If no extension is found, return the entire filepath
    return filepath

def _get_metabolomicsworkbench_filepaths(study_id):

    try:
        dataset_list_url = "https://www.metabolomicsworkbench.org/data/show_archive_contents_json.php?STUDY_ID={}".format(
            study_id)
        mw_file_list = requests.get(dataset_list_url).json()
        workbench_df = pd.DataFrame(mw_file_list)

        workbench_df['raw_sample_name'] = workbench_df['URL'].apply(extract_first_folder_with_extension)


        workbench_df["USI_file"] = workbench_df["URL"].apply(
            lambda url: f"mzspec:{study_id}:{parse_qs(urlparse(url).query).get('F', [None])[0]}"
        )





        workbench_df["USI_sample"] = workbench_df.apply(
            lambda
                row: f"mzspec:{study_id}:{parse_qs(urlparse(row['URL']).query).get('A', [None])[0]}-{row['raw_sample_name']}",
            axis=1
        )
        
    except KeyboardInterrupt:
        raise
    except:
        workbench_df = pd.DataFrame()

    return workbench_df

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
        existing_datasets = set(existing_datasets['datasets'])
        return existing_datasets
    except:
        return set()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Give an MWB study ID and get a tsv with file paths present in the study.')
    parser.add_argument("--study_id", "-mwb_id", type=str, help='An MWB study ID such as "ST002050", "ALL" for every study', required=True)
    parser.add_argument("--output_path", type=str, help='Output file path to tsv file.')
    parser.add_argument("--filter_extensions", type=str, help='Filter extensions to ".mzml", ".mzxml", ".cdf", ".raw", ".wiff", and ".d".', default='False')
    parser.add_argument("--existing_datasets", type=str, help="path to a file of datasets already indexed", default="none")


    args = parser.parse_args()

    existing_datasets = _get_existing_datasets(args.existing_datasets)

    if args.study_id == "ALL":
        # Getting all files
        url = "https://www.metabolomicsworkbench.org/rest/study/study_id/ST/available"
        studies_dict = requests.get(url).json()

        study_list = []
        for key in studies_dict.keys():
            study_dict = studies_dict[key]
            study_list.append(study_dict['study_id'])

        study_list = list(set(study_list))

        all_results_list = []
        for study_id in tqdm.tqdm(study_list):
            if study_id in existing_datasets:
                print("Skipping", study_id, "Already indexed")
                continue

            try:
                temp_result_df = _get_metabolomicsworkbench_filepaths(study_id=study_id)
                all_results_list.append(temp_result_df)
            except KeyboardInterrupt:
                raise
            except:
                pass

        result_df = pd.concat(all_results_list, axis=0)

    else:
        result_df = _get_metabolomicsworkbench_filepaths(study_id=args.study_id)

    result_df['study_id'] = result_df['STUDY_ID']
    result_df['file_path'] = result_df['FILENAME'].apply(process_filename)
    result_df['file_path'] = result_df['file_path'].apply(clean_mac_path)
    
    result_df['USI'] = result_df['USI_file'].apply(process_filename)
    result_df['USI'] = result_df['USI'].apply(clean_mac_path)

    result_df = result_df.drop_duplicates(keep='first')

    if args.filter_extensions == 'True':
        extensions = [".mzml", ".mzxml", ".cdf", ".raw", ".wiff", ".d"]
        result_df = result_df[result_df['FILENAME'].str.lower().str.endswith(tuple(extensions))]


    result_df = result_df[['study_id', 'file_path', 'USI']]
     

    result_df.to_csv(args.output_path, sep='\t', index=False, header=True)

    print(f"Output written to {args.output_path}")