import os
import pandas as pd
import requests
import argparse
from urllib.parse import urlparse, parse_qs


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
            lambda url: f"mzspec:{study_id}:{parse_qs(urlparse(url).query).get('A', [None])[0]}-{parse_qs(urlparse(url).query).get('F', [None])[0]}"
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



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Give an MWB study ID and get a tsv with file paths present in the study.')
    parser.add_argument("--study_id", "-mwb_id", type=str, help='An MWB study ID such as "ST002050", "ALL" for every study', required=True)
    parser.add_argument("--output_path", type=str, help='Output file path to tsv file.')

    args = parser.parse_args()

    extensions = [".mzml", ".mzxml", ".cdf", ".raw", ".wiff", ".d"]

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
        for study_id in study_list:
            print("Downloading", study_id)

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

    result_df = result_df[result_df['FILENAME'].str.lower().str.endswith(tuple(extensions))]

    result_df.to_csv(args.output_path, sep='\t', index=False, header=True)

    print(f"Output written to {args.output_path}")