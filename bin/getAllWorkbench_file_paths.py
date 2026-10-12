import os
import re
import pandas as pd
import requests
import argparse
import tqdm
from urllib.parse import urlparse, parse_qs
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry


# Shared session: retries throttling/server errors with backoff instead of losing the study
session = requests.Session()
session.mount("https://", HTTPAdapter(max_retries=Retry(total=5, backoff_factor=2,
                                                         status_forcelist=[429, 500, 502, 503, 504],
                                                         allowed_methods=["GET"])))
REQUEST_TIMEOUT = (10, 300)


def _get_json(url):
    response = session.get(url, timeout=REQUEST_TIMEOUT)
    response.raise_for_status()
    return response.json()



def is_apple_double(path):
    # macOS resource-fork copies (__MACOSX/.../._x.mzML) are not data files
    return "__MACOSX/" in path or os.path.basename(path).startswith("._")

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
        mw_file_list = _get_json(dataset_list_url)
        workbench_df = pd.DataFrame(mw_file_list)
        if len(workbench_df) == 0:
            return workbench_df
        workbench_df = workbench_df[~workbench_df['FILENAME'].apply(is_apple_double)]

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
    except Exception as e:
        print(f"Failed to list files for {study_id}: {type(e).__name__}: {e}")
        workbench_df = pd.DataFrame()

    return workbench_df

# Some archive listings carry raw 7z listing columns in front of the path,
# e.g. '2014-10-05 10:45:02 ..H.A           80               MX183900 .../x.d'
_7Z_LISTING_PREFIX = re.compile(r'^\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2} \S{5} +\d+ +(?:\d+ +)?')


def clean_listing_path(path):
    return _7Z_LISTING_PREFIX.sub('', path)


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
        studies_dict = _get_json(url)

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

            temp_result_df = _get_metabolomicsworkbench_filepaths(study_id=study_id)
            all_results_list.append(temp_result_df)

        result_df = pd.concat(all_results_list, axis=0)

    else:
        result_df = _get_metabolomicsworkbench_filepaths(study_id=args.study_id)

    result_df['study_id'] = result_df['STUDY_ID']
    result_df['file_path'] = result_df['FILENAME'].apply(clean_listing_path).apply(process_filename)

    result_df['USI'] = result_df.apply(
        lambda row: f"mzspec:{row['study_id']}:{process_filename(clean_listing_path(parse_qs(urlparse(row['URL']).query).get('F', [''])[0]))}",
        axis=1)

    if args.filter_extensions == 'True':
        extensions = [".mzml", ".mzxml", ".cdf", ".raw", ".wiff", ".d"]
        # filter on the collapsed path so .d folders listed only by their inner files are kept
        result_df = result_df[result_df['file_path'].str.lower().str.endswith(tuple(extensions))]


    result_df = result_df[['study_id', 'file_path', 'USI']]

    # one row per file; .d folders collapse many listed files into one path
    result_df = result_df.drop_duplicates(subset=['study_id', 'file_path'], keep='first')
     

    result_df.to_csv(args.output_path, sep='\t', index=False, header=True)

    print(f"Output written to {args.output_path}")