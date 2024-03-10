import os
import argparse
import pandas as pd
import re
import numpy as np
import requests
from bs4 import BeautifulSoup
from urllib.parse import urlparse, parse_qs
import json
from tqdm import tqdm
import time
from extend_allowed_terms import adapt_allowed_terms 
from REDU_conversion_functions import get_taxonomy_id_from_name__allowedTerms
from read_and_validate_redu_from_github import complete_and_fill_REDU_table
from REDU_conversion_functions import update_unassigned_terms
from REDU_conversion_functions import get_taxonomy_id_from_name
from REDU_conversion_functions import get_uberon_table
from REDU_conversion_functions import get_ontology_table
from REDU_conversion_functions import age_category
from REDU_conversion_functions import get_taxonomy_info
from collections import Counter
import traceback


def clean_path(path):
    if path.startswith("__MACOSX/") and '/._' in path:
        # Remove the __MACOSX part and the ._ prefix from the filename
        return path.replace("__MACOSX/", "").replace("/._", "/")
    return path


def merge_repeated_fileobservations_across_mwatb(df, **kwargs):

    if 'polarity_table' in kwargs.keys():
        polarity_table = kwargs['polarity_table']
    else:
        polarity_table = pd.DataFrame(columns=['USI', 'negative_polarity_count', 'positive_polarity_count'])


    def process_group(group):
        # Initialize the list for columns to make the most common
        li_make_most_common = []
        
        # Check if all entries in 'NCBITaxonomy' for the group are the same
        if len(set(group['NCBITaxonomy'])) == 1:
            li_make_most_common = ['SampleType', 'SampleTypeSub1']
        
            # Find the most common value for specified columns and update the group
            for col in li_make_most_common:
                most_common_value = Counter(group[col]).most_common(1)[0][0]
                if most_common_value is None or most_common_value == 'None':
                    most_common_value = 'ML import: not available'
                group[col] = most_common_value

        # Handle discrepancies for columns not in 'li_make_most_common'
        li_rm_discrepancies = [col for col in group.columns if col not in li_make_most_common]
        for col in li_rm_discrepancies:
            if len(set(group[col])) > 1 or all(value is None for value in group[col]):
                group[col] = 'ML import: not available'
            else:
                # If there's no discrepancy, keep the original value
                group[col] = group[col].iloc[0]
        
        return group

    df = df.merge(polarity_table, left_on='USI', right_on='USI', how='left')

    df = df[~((df['IonizationSourceAndPolarity'].str.contains('negative')) & (df['positive_polarity_count'] > 0))]
    df = df[~((df['IonizationSourceAndPolarity'].str.contains('positive')) & (df['negative_polarity_count'] > 0))]
    df = df[~((df['IonizationSourceAndPolarity'].str.contains('alternating')) & ((df['negative_polarity_count'] == 0) | (df['positive_polarity_count'] == 0)))]


    # Group by 'USI' and apply the processing function to each group
    df = df.groupby('filename').apply(process_group)

    # Remove duplicates based on 'USI' column
    df = df.drop_duplicates(subset='filename')

    return df

def handle_duplicates(ordered_pairs):
    """
    Append numbers to duplicate keys. This is necessary because the RAW_FILE_NAME key often appears more than once
    in the mwTAb files
    """
    d = {}
    for k, v in ordered_pairs:
        if k in d:
            count = 1
            new_key = f"{k}_{count}"
            while new_key in d:
                count += 1
                new_key = f"{k}_{count}"
            d[new_key] = v
        else:
            d[k] = v
    return d


def _get_metabolomicsworkbench_files(dataset_accession):
    # Lets see if it is in massive
    try:
        msv_accession = _accession_to_msv_accession(dataset_accession)
        files_df = _get_massive_files(msv_accession)
    except:
        msv_accession = None
        files_df = pd.DataFrame()
    try:
        dataset_list_url = "https://www.metabolomicsworkbench.org/data/show_archive_contents_json.php?STUDY_ID={}".format(
            dataset_accession)
        mw_file_list = requests.get(dataset_list_url).json()

        acceptable_extensions = [".mzml", ".mzxml", ".cdf", ".raw", ".wiff", ".d"]

        mw_file_list = [file_obj for file_obj in mw_file_list if
                        file_obj["FILENAME"].lower().endswith(tuple(acceptable_extensions))]
        workbench_df = pd.DataFrame(mw_file_list)
        workbench_df["filename"] = workbench_df["FILENAME"]
        workbench_df['filename'] = workbench_df['filename'].apply(clean_path)
        workbench_df["size_mb"] = (workbench_df["FILESIZE"].astype(int) / 1024 / 1024).astype(int)
        workbench_df["URL"] = workbench_df["URL"]        

        # Generate USI column
        workbench_df["USI"] = workbench_df["URL"].apply(
            lambda url: "mzspec:{}:{}".format(
                dataset_accession,
                parse_qs(urlparse(url).query).get('F', [None])[0]
            )
        )

        workbench_df = workbench_df[["filename", "size_mb", "USI"]]
    except:
        workbench_df = pd.DataFrame()

    merged_df = pd.concat([files_df, workbench_df])
    #print(merged_df)
    return merged_df, msv_accession



def get_taxonomy_id_from_name(species_name):
    # Query NCBI's Entrez API to get the NCBI ID for the species based on its Latin name
    if species_name is not None and species_name != "NA" and species_name != "N/A":
        try:
            response = requests.get(
                f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=taxonomy&term={species_name}"
            )
            soup = BeautifulSoup(response.text, "xml")
            ncbi_id = soup.find("IdList").find("Id").text
            return ncbi_id
        except Exception:
            return None
    else:
        return None

def convert_to_numeric_or_original(x):
    try:
        numeric_x = pd.to_numeric(x)
        return 'numeric'
    except ValueError:
        return x


def print_dict_with_exclude(d, exclude=[]):
    for key, value in d.items():
        if key in exclude:
            continue
        if isinstance(value, list):
            print(f"{key}: {value[:3]}")
        else:
            print(f"{key}: {value}")


def get_enviromental_water(x):
    x = x.lower()

    if 'water' in x or 'sewerage' in x:
        if 'waste' in x or 'sewerage' in x:
            return ['environmental', 'water_waste']
        if 'surface' in x:
            return ['environmental', 'water_surface']
        if 'ground' in x:
            return ['environmental', 'water_ground']
        if 'storm' in x:
            return ['environmental', 'water_storm']
        if 'sea' in x or 'ocean' in x or 'coast':
            return ['environmental', 'water_seawater']
        else:
            return [None, None]
    else:
        return [None, None]


def get_blanks(x):
    x = x.lower()
    if 'blank' in x:
        if 'blank' == x:
            return ['blank_analysis', 'blank_analysis']
        if 'method' in x or 'extraction' in x:
            return ['blank_extraction', 'blank_extraction']
        if 'media' in x:
            return ['blank_culturemedia', 'blank_culturemedia']
        else:
            return [None, None]
    else:
        return [None, None]

def get_raw_file_names(x, y, raw_file_name_df):
    raw_file_names = list(raw_file_name_df['filename_base'])
    raw_file_nameswo_extension = [os.path.splitext(filename)[0] for filename in raw_file_names]

    if x in raw_file_names:
        filtered_df = raw_file_name_df[raw_file_name_df['filename_base'] == x]
        if not filtered_df.empty:
            return [x, filtered_df.iloc[0]['filename']]
    elif x in raw_file_nameswo_extension:
        filtered_df = raw_file_name_df[raw_file_name_df['filename_base'].apply(lambda f: os.path.splitext(f)[0]) == x]
        if not filtered_df.empty:
            return [x, filtered_df.iloc[0]['filename']]
    elif y in raw_file_names:
        filtered_df = raw_file_name_df[raw_file_name_df['filename_base'] == y]
        if not filtered_df.empty:
            return [y, filtered_df.iloc[0]['filename']]
    elif y in raw_file_nameswo_extension:
        filtered_df = raw_file_name_df[raw_file_name_df['filename_base'].apply(lambda f: os.path.splitext(f)[0]) == y]
        if not filtered_df.empty:
            return [y, filtered_df.iloc[0]['filename']]

    return [None, None]

def convert_sex(x):
    x = x.lower()
    if x in ['f', 'female']:
        return 'female'
    elif x in ['m', 'male']:
        return 'male'
    elif x in ['as', 'asexual']:
        return 'asexual'
    else:
        return 'ML import: not available'


def convert_to_years(age):
    num = int(re.findall(r'\b\d+\b', age)[0]) if re.findall(r'\b\d+\b', age) else None
    if num is None:
        return None
    if "week" in age:
        return num / 52.1429
    elif "day" in age:
        return num / 365
    elif "year" in age:
        return num
    else:
        return num


def extract_years(dates):
    try:
        pattern = re.compile(r'\d{1,2}[-/.]\d{1,2}[-/.]\d{4}')
        matches = pattern.findall(dates)
        if not matches:
            pattern = re.compile(r'\d{4}')
            matches = pattern.findall(dates)
        years = max([int(date.split('/')[-1]) if '/' in date else int(date.split('.')[-1]) for date in matches])
        if isinstance(years, int):
            if years >= 1998 and years <= 2030:
                return years
            else:
                return None
        else:
            return None
    except Exception:
        return None


def get_key_info_into_outer(df, key_vars, new_col):
    df[new_col] = df.loc[df['Key'].isin(key_vars), 'Value']
    df[new_col] = df.groupby('filename')[new_col].transform('first')
    if new_col == 'Latitude' or new_col == 'Longitude':
        df[new_col] = pd.to_numeric(df[new_col], errors='coerce')
    else:
        df[new_col] = df[new_col].fillna(value='NA')
    return df

def get_rawFile_names(df, key_vars, new_col):
    # Create a pattern to match key_vars and their numbered versions
    pattern = r'(' + '|'.join(key_vars) + r')(_\d+)?'

    # Extract the base key name without the number
    df['base_key'] = df['Key'].str.extract(pattern, expand=False)[0]

    filename_mapping = df[df['base_key'].notna()].set_index('script_id')['Value'].to_dict()

    # Apply the mapping to create the new filename column
    df[new_col] = df['script_id'].map(filename_mapping)

    # Handle numeric conversion for specific columns
    if new_col in ['Latitude', 'Longitude']:
        df[new_col] = pd.to_numeric(df[new_col], errors='coerce')
    else:
        df[new_col] = df[new_col].fillna(value='NA')

    # Drop the temporary 'base_key' column
    df.drop('base_key', axis=1, inplace=True)
    return df



def find_unique_match_index(row, target_df, target_column):
    # Find all potential matches based on the condition
    matches = [target_row[target_column] for _, target_row in target_df.iterrows()
               if target_row[target_column] in row['filename_raw_lower'] or row['filename_raw_lower'] in target_row[target_column]]
        
    # Deduplicate the matches to consider only unique values
    unique_matches = set(matches)
        
    # Check if there's exactly one unique match
    if len(unique_matches) == 1:
        match_value = unique_matches.pop()
        match_index = target_df[target_df[target_column] == match_value].index[0]
        return match_index, False  # No ambiguity
    else:
        return pd.NA, True  # Ambiguity exists



def map_potential_matches(df, target_df, target_column, source_columns):
    from tqdm import tqdm
    # Temporary storage for potential matches
    potential_matches = {}

    # Preprocess to get unique values from 'filename_base_wo_extension'
    unique_values = df['filename_base_wo_extension'].unique()

    # Identify potential matches for each unique value with a progress bar
    for value in tqdm(unique_values, desc="Mapping potential matches"):
        matches = [target_row[target_column] for _, target_row in target_df.iterrows()
                   if value in target_row[target_column] or target_row[target_column] in value]
        
        # Store unique matches only, mapping back to all rows with this value
        potential_matches[value] = set(matches)

    return potential_matches



def analyze_for_contradictions(potential_matches):
    # Reverse mapping to see which df values map to the same raw_file_name_df value
    reverse_mapping = {}
    for df_value, matches in potential_matches.items():
        for match in matches:
            if match in reverse_mapping:
                reverse_mapping[match].add(df_value)
            else:
                reverse_mapping[match] = {df_value}
    
    # Identify contradictions
    contradictions = {raw_value: df_values for raw_value, df_values in reverse_mapping.items() if len(df_values) > 1}
    
    return contradictions


def resolve_contradictions(contradictions, potential_matches):
    for raw_value, conflicting_df_values in contradictions.items():
        # Determine the "better" match for each contradiction
        better_match = None
        max_overlap = 0
        
        for df_value in conflicting_df_values:
            overlap = len(set(df_value) & set(raw_value))  # Simple way to estimate overlap; consider refining based on your actual data structure
            if overlap > max_overlap:
                max_overlap = overlap
                better_match = df_value
        
        # Update potential_matches to only include the better match for this raw_value
        if better_match:
            for df_value in conflicting_df_values:
                if df_value != better_match:
                    # Remove the raw_value from the non-better matches
                    potential_matches[df_value].discard(raw_value)
    
    return potential_matches

def filter_unresolved_contradictions(potential_matches, contradictions):
    # This function assumes 'contradictions' is updated to reflect any unresolved contradictions
    # Remove or mark potential matches that are part of unresolved contradictions
    for conflict_value in contradictions:
        for df_value in contradictions[conflict_value]:
            if df_value in potential_matches:
                potential_matches[df_value] = set()  # Remove match or mark as unresolved

    return potential_matches

def create_dataframe_from_SUBJECT_SAMPLE_FACTORS(data, raw_file_name_df=None, path_to_csvs='translation_sheets', attempt_taxa_extraction_from_inner = False, add_to_cols = [], **kwargs):
    data_list = []
    for id, item in enumerate(data):
        subject_id = item.get('Subject ID', 'NA')
        sample_id = item.get('Sample ID', 'NA')
        factors = item.get('Factors', {})
        additional_data = item.get('Additional sample data', {})

        for key, value in factors.items():
            data_list.append([id, subject_id, sample_id, key, value])
        for key, value in additional_data.items():
            data_list.append([id, subject_id, sample_id, key, value])
    df = pd.DataFrame(data_list, columns=['script_id', 'SubjectIdentifierAsRecorded', 'filename', 'Key', 'Value'])
    df['Key'] = df['Key'].str.lower()

    ontology_table = kwargs['ontology_table']

    #raw data name might be in the Sample ID, Factors
    expected_raw_file_keys = ["raw_file_name",  "rawfilename", "raw_file",
                              "datafile name", "raw files", "raw file name", 
                              "rawfile name", "rpm_filename", "zhp_filename",
                              "data file"]


    if not any(key in df['Key'].values for key in expected_raw_file_keys):
        print('no raw file key found.')

        if any("file" in key for key in df['Key'].values):
            print("Found data files with unknown name")
            df['Key'] = df['Key'].apply(lambda x: "raw_file_name" if "file" in x else x)


        if not any(key in df['Key'].values for key in expected_raw_file_keys):
            #if we have no raw file key try to use the SampleID as raw_file_name
            df_new_row = df.drop_duplicates(['script_id', 'SubjectIdentifierAsRecorded', 'filename']).copy()
            df_new_row['Key'] = 'raw_file_name'
            df_new_row['Value'] = df_new_row['filename']
            df = pd.concat([df, df_new_row], ignore_index=True)

    if any(key in df['Key'].values for key in expected_raw_file_keys):
        print('raw file key found.')
        #attempt to get raw files from associated keys
        raw_file_name_df = raw_file_name_df.rename(columns={'filename': 'filename_raw_path'})
        df = get_rawFile_names(df, key_vars=expected_raw_file_keys, new_col="filename_raw")
        #explode cases where we have multiple files by separators
        if raw_file_name_df['filename_base'].str.contains(" ").any():
            split_by_this = '[;,]+'
        else:
            split_by_this = '[; ,]+'

        #df = df.assign(filename_raw=df['filename_raw'].str.split(split_by_this)).explode('filename_raw')
        df = df.assign(filename_raw=df['filename_raw'].str.split(split_by_this).apply(lambda x: list(set(x)))).explode('filename_raw')


        df.reset_index(drop=True, inplace=True)

        raw_file_name_df['filename_base_lower'] = raw_file_name_df['filename_base'].str.lower()
        df['filename_raw_lower'] = df['filename_raw'].str.lower()
        df['filename_lower'] = df['filename'].str.lower()

        if bool(set(df['filename_raw_lower']) & set(raw_file_name_df['filename_base_lower'])):
            df = df.merge(raw_file_name_df, left_on='filename_raw_lower', right_on='filename_base_lower', how='left')
            df.drop(['filename_raw_lower', 'filename_lower', 'filename_base_lower'], axis=1, inplace=True)
        elif bool(set(df['filename_lower']) & set(raw_file_name_df['filename_base_lower'])):
            df = df.merge(raw_file_name_df, left_on='filename_lower', right_on='filename_base_lower', how='left')
            df.drop(['filename_raw_lower', 'filename_lower', 'filename_base_lower'], axis=1, inplace=True)
    
        else:
            df['filename_base_wo_extension'] = df['filename_raw_lower'].str.split('.').str[0]
            raw_file_name_df['filename_base_wo_extension'] = raw_file_name_df['filename_base'].str.split('.').str[0].str.lower()
            raw_file_name_df['filename_raw'] = raw_file_name_df['filename_raw_path'].str.lower()
            raw_file_name_df['filename_raw_lower'] = raw_file_name_df['filename_raw'].str.lower()

            
            match_df = 'filename_raw_lower'
            match_raw_file_name_df = 'filename_base_wo_extension'
            potential_matches = set(df[match_df]).intersection(set(raw_file_name_df[match_raw_file_name_df]))

            if not potential_matches:
                match_df = 'filename_raw_lower'
                match_raw_file_name_df = 'filename_raw_lower'
                potential_matches = set(df[match_df]).intersection(set(raw_file_name_df[match_raw_file_name_df]))

            if not potential_matches:
                match_df = 'filename_base_wo_extension'
                match_raw_file_name_df = 'filename_base_wo_extension'
                potential_matches = set(df[match_df]).intersection(set(raw_file_name_df[match_raw_file_name_df]))

            if potential_matches:
                df = df[['script_id', 'SubjectIdentifierAsRecorded', 'filename', 'Key', 'Value', match_df]]
                df = df.merge(raw_file_name_df, left_on=match_df, right_on=match_raw_file_name_df, how='left')
            else:

                print('getting substring matches:')
                potential_matches = map_potential_matches(df, raw_file_name_df, target_column='filename_base_wo_extension', source_columns = ['filename_base_wo_extension', ''])
                print('analyzing contradictions:')
                contradictions = analyze_for_contradictions(potential_matches)
                print('resolving contradictions:')
                potential_matches = resolve_contradictions(contradictions, potential_matches)
                print('filtering contradictions:')
                potential_matches = filter_unresolved_contradictions(potential_matches, contradictions)

                
                match_indices = {}
                for df_value, matches in potential_matches.items():
                    if matches:  
                        
                        match_value = matches.pop()  # Get the match value
                        match_index = raw_file_name_df.index[raw_file_name_df['filename_base_wo_extension'] == match_value].tolist()[0]  # Assuming this finds the index
                        match_indices[df_value] = match_index

                df['match_index'] = df['filename_base_wo_extension'].map(match_indices)
                
                df = pd.merge(df, raw_file_name_df, left_on='match_index', right_index=True, how='left', suffixes=('', '_matched'))

                df['filename_raw'] = df['filename_raw_path']
                df['filename'] = df['filename_raw_path']
                df.drop(columns=['match_index'], inplace=True)

            df['filename_raw_lower'] = df['filename_base']
            df = df.drop(columns=['filename_base'])

    df['Value'] = df['Value'].str.lower()
    df['SubjectIdentifierAsRecorded'] = df['SubjectIdentifierAsRecorded'].replace('-', '')
    df = get_key_info_into_outer(df, key_vars=["gender", "sex", "gender (f/m)", "biological_sex"], new_col="MWB_sex")
    df = get_key_info_into_outer(df, key_vars=["age", "age (years)"], new_col="MWB_age")
    df = get_key_info_into_outer(df, key_vars=["collection_country", "collection country", "country", "site", "location"],
                                 new_col="Country")
    df = get_key_info_into_outer(df, key_vars=["latitude"], new_col="Latitude")
    df = get_key_info_into_outer(df, key_vars=["longitude"], new_col="Longitude")
    df['LatitudeandLongitude'] = df.apply(
        lambda x: f'{x.Latitude}|{x.Longitude}' if x.Latitude is not None and x.Longitude is not None and not np.isnan(
            x.Latitude) and not np.isnan(x.Longitude) else None, axis=1)

    df[['SampleType_inner', 'SampleTypeSub1_inner']] = df.apply(lambda x: get_blanks(x['Value']), axis=1,
                                                                result_type='expand')

    df = translate_MWB_to_REDU_from_csv(df, case='inner', path_to_csvs=path_to_csvs, add_to_cols=add_to_cols, ontology_table=ontology_table, allowedTerm_dict=allowedTerm_dict)
    df = df.drop(columns=['Key', 'Value', 'Longitude', 'Latitude'])
    df = df.drop_duplicates().reset_index(drop=True)

    df['filename'] = df['filename_raw_path']

    return df


def create_dataframe_from_SUBJECT_SAMPLE_FACTORS_collapsed_factors(data):
    data_list = []
    for item in data:
        subject_id = item.get('Subject ID', 'NA')
        sample_id = item.get('Sample ID', 'NA')
        factors = " ".join([f"{key}: {value}" for key, value in item.get('Factors', {}).items()])
        additional_data = " ".join([f"{key}: {value}" for key, value in item.get('Additional sample data', {}).items()])
        data_list.append([subject_id, sample_id, factors + " " + additional_data])
    df = pd.DataFrame(data_list, columns=['SubjectIdentifierAsRecorded', 'filename', 'Factors'])
    return df


def create_dataframe_outer_dict(MWB_mwTAB_dict, raw_file_name_df=None, path_to_csvs = 'translation_sheets', **kwargs):
    entry = {
        'ANALYSIS_TYPE': MWB_mwTAB_dict.get('ANALYSIS', {}).get('ANALYSIS_TYPE', 'NA'),
        'PROJECT_TITLE': MWB_mwTAB_dict.get('PROJECT', {}).get('PROJECT_TITLE', 'NA'),
        'SUBJECT_TYPE': MWB_mwTAB_dict.get('SUBJECT', {}).get('SUBJECT_TYPE', 'NA'),
        'SPECIES_GROUP': MWB_mwTAB_dict.get('SUBJECT', {}).get('SPECIES_GROUP', 'NA'),
        'SUBJECT_SPECIES': MWB_mwTAB_dict.get('SUBJECT', {}).get('SUBJECT_SPECIES', 'NA'),
        'SAMPLE_TYPE': MWB_mwTAB_dict.get('COLLECTION', {}).get('SAMPLE_TYPE', 'NA'),
        'COLLECTION_LOCATION': MWB_mwTAB_dict.get('COLLECTION', {}).get('COLLECTION_LOCATION', 'NA'),
        'SampleCollectionMethod': MWB_mwTAB_dict.get('COLLECTION', {}).get('COLLECTION_METHOD', 'NA'),
        'COLLECTION_TIME': MWB_mwTAB_dict.get('COLLECTION', {}).get('COLLECTION_TIME', 'NA'),
        'MassSpectrometer': MWB_mwTAB_dict.get('MS', {}).get('INSTRUMENT_NAME', 'NA'),
        'MS_TYPE': MWB_mwTAB_dict.get('MS', {}).get('MS_TYPE', 'NA'),
        'ION_MODE': MWB_mwTAB_dict.get('MS', {}).get('ION_MODE', 'NA'),
        'InternalStandardsUsed': MWB_mwTAB_dict.get('CHROMATOGRAPHY', {}).get('INTERNAL_STANDARD', 'NA'),
        'COLUMN_NAME': MWB_mwTAB_dict.get('CHROMATOGRAPHY', {}).get('COLUMN_NAME', 'NA'),
        'CHROMATOGRAPHY_TYPE': MWB_mwTAB_dict.get('CHROMATOGRAPHY', {}).get('CHROMATOGRAPHY_TYPE', 'NA'),
        'SampleExtractionMethod': MWB_mwTAB_dict.get('SAMPLEPREP', {}).get('EXTRACTION_METHOD', 'NA'),
        'TAXONOMY_ID': MWB_mwTAB_dict.get('SUBJECT', {}).get('TAXONOMY_ID', 'NA'),
        'ACQUISITION_DATE': MWB_mwTAB_dict.get('ANALYSIS', {}).get('ACQUISITION_DATE', 'NA'),
        'CREATED_ON': MWB_mwTAB_dict.get('METABOLOMICS WORKBENCH', {}).get('CREATED_ON', 'NA'),
        'STUDY_ID': MWB_mwTAB_dict.get('METABOLOMICS WORKBENCH', {}).get('STUDY_ID', 'NA')

    }

    ontology_table = kwargs['ontology_table']

    if entry['ANALYSIS_TYPE'] != 'MS':
        raise ValueError("This is not an MS analysis!")

    df_outer = pd.DataFrame(entry, index=[0])

    #add NCBITaxonomy
    #######
    unique_species = set(df_outer['SUBJECT_SPECIES'])
    processed_species = {species: get_taxonomy_id_from_name__allowedTerms(species, allowedTerm_dict=allowedTerm_dict) for species in unique_species}
    df_outer['NCBITaxonomy'] = df_outer['SUBJECT_SPECIES'].map(processed_species)

    attempt_taxa_extraction_from_inner = False
    if all(value is None for value in processed_species.values()):
        attempt_taxa_extraction_from_inner = True


    #add LC method
    #######
    if 'CHROMATOGRAPHY_TYPE' in df_outer.columns:
        df_outer['CHROMATOGRAPHY_TYPE'] = df_outer['CHROMATOGRAPHY_TYPE'].fillna('')
        df_outer['CHROMATOGRAPHY_TYPE'] = df_outer['CHROMATOGRAPHY_TYPE'].str.lower()
        df_outer['ChromatographyAndPhase'] = ''
        df_outer.loc[df_outer['CHROMATOGRAPHY_TYPE'].str.contains('reverse'), 'ChromatographyAndPhase'] = 'reverse phase'
        df_outer.loc[df_outer['CHROMATOGRAPHY_TYPE'].str.contains('hilic'), 'ChromatographyAndPhase'] = 'normal phase (HILIC)'
        df_outer.loc[df_outer['CHROMATOGRAPHY_TYPE'].str.contains('normal phase'), 'ChromatographyAndPhase'] = 'normal phase (HILIC)'
        df_outer.loc[df_outer['CHROMATOGRAPHY_TYPE'].str.contains('gc'), 'ChromatographyAndPhase'] = 'gas chromatography (DB-5)'
        
        if 'COLUMN_NAME' in df_outer.columns:
            df_outer['COLUMN_NAME'] = df_outer['COLUMN_NAME'].str.lower()
            df_outer.loc[df_outer['COLUMN_NAME'].str.contains('Phenyl', case=False) & df_outer['COLUMN_NAME'].str.contains('Hexyl', case=False), 'ChromatographyAndPhase'] = ' (Phenyl-Hexyl)'
            df_outer.loc[df_outer['COLUMN_NAME'].str.contains('C18') & df_outer['COLUMN_NAME'].str.contains('polar', case=False), 'ChromatographyAndPhase'] = ' (polar-C18)'
            
            df_outer.loc[df_outer['COLUMN_NAME'].str.contains('HSS T3', case=False) & ~df_outer['ChromatographyAndPhase'].str.contains(r"\("), 'ChromatographyAndPhase'] += ' (C18)'
            df_outer.loc[df_outer['COLUMN_NAME'].str.contains('C18') & ~df_outer['ChromatographyAndPhase'].str.contains(r"\("), 'ChromatographyAndPhase'] += ' (C18)'
            df_outer.loc[df_outer['COLUMN_NAME'].str.contains('C30') & ~df_outer['ChromatographyAndPhase'].str.contains(r"\("), 'ChromatographyAndPhase'] += ' (C30)'
            df_outer.loc[df_outer['COLUMN_NAME'].str.contains('C8') & ~df_outer['ChromatographyAndPhase'].str.contains(r"\("), 'ChromatographyAndPhase'] += ' (C8)'
            df_outer.loc[df_outer['ChromatographyAndPhase'].str.contains('reverse phase') & ~df_outer['ChromatographyAndPhase'].str.contains(r"\("), 'ChromatographyAndPhase'] += ' (NOS)'
        else:
            df_outer.loc[df_outer['ChromatographyAndPhase'].str.contains('reverse phase') & ~df_outer['ChromatographyAndPhase'].str.contains(r"\("), 'ChromatographyAndPhase'] += ' (NOS)'


    df_outer['YearOfAnalysis'] = df_outer['ACQUISITION_DATE'].apply(lambda x: extract_years(x))
    df_outer['YearOfAnalysis'] = df_outer.apply(
        lambda x: extract_years(x['CREATED_ON']) if x['YearOfAnalysis'] is None else x['YearOfAnalysis'], axis=1)

    df_outer['IonizationSourceAndPolarity'] = df_outer['MS_TYPE'] + '|' + df_outer['ION_MODE']

    df_outer.loc[:, ['SampleType', 'SampleTypeSub1']] = 'ML import: not available'

    processed_taxonomy = {taxonomy.split('|')[0]: get_taxonomy_info(taxonomy.split('|')[0])
                        for taxonomy in df_outer['NCBITaxonomy'].unique()
                        if taxonomy is not None and '|' in taxonomy and 'None' not in taxonomy}


    df_outer.loc[:, 'SampleType'] = df_outer['NCBITaxonomy'].map(
        lambda x: processed_taxonomy.get(x.split('|')[0], [pd.NA, pd.NA])[0]
        if x is not None and '|' in x and 'None' not in x else pd.NA
    )

    df_outer.loc[:, 'SampleTypeSub1'] = df_outer['NCBITaxonomy'].map(
        lambda x: processed_taxonomy.get(x.split('|')[0], [pd.NA, pd.NA])[1]
        if x is not None and '|' in x and 'None' not in x else pd.NA
    )


    df_outer[['SampleType', 'SampleTypeSub1']] = df_outer.apply(
        lambda x: get_enviromental_water(x['SAMPLE_TYPE'] + x['SUBJECT_TYPE']) 
        if pd.isnull(x['SampleType']) and pd.isnull(x['SampleTypeSub1']) 
        else (x['SampleType'], x['SampleTypeSub1']), axis=1, result_type='expand'
    )

    if df_outer['NCBITaxonomy'].iloc[0] is None:
        add_to_cols = ['NCBITaxonomy']
    else:
        add_to_cols = []

    df_inner_SUBJECT_SAMPLE_FACTORS = create_dataframe_from_SUBJECT_SAMPLE_FACTORS(
        MWB_mwTAB_dict['SUBJECT_SAMPLE_FACTORS'],
        raw_file_name_df=raw_file_name_df,
        path_to_csvs=path_to_csvs,
        ontology_table=ontology_table,
        attempt_taxa_extraction_from_inner=attempt_taxa_extraction_from_inner,
        add_to_cols=add_to_cols
    )

    if 'NCBITaxonomy' in df_inner_SUBJECT_SAMPLE_FACTORS.columns:
        df_outer.drop('NCBITaxonomy', axis=1, inplace=True)


    df_outer = pd.concat([df_outer] * len(df_inner_SUBJECT_SAMPLE_FACTORS), ignore_index=True)
    df_outer_inner = pd.concat([df_outer, df_inner_SUBJECT_SAMPLE_FACTORS], axis=1)

    df_outer_inner.loc[df_outer_inner['SampleType_inner'].notnull(), 'SampleType'] = df_outer_inner['SampleType_inner']
    df_outer_inner.loc[df_outer_inner['SampleTypeSub1_inner'].notnull(), 'SampleTypeSub1'] = df_outer_inner[
        'SampleTypeSub1_inner']

    df_outer_inner = df_outer_inner.drop(columns=['SampleType_inner', 'SampleTypeSub1_inner'])
    
    return df_outer_inner


def translate_MWB_to_REDU_from_csv(MWB_table,
                                   column_and_csv_names_outer=['MassSpectrometer',
                                                               'ChromatographyAndPhase',
                                                               'InternalStandardsUsed',
                                                               'SampleExtractionMethod',
                                                               'NCBITaxonomy',
                                                               'SampleCollectionMethod',
                                                               'IonizationSourceAndPolarity',
                                                               'Country'],
                                   column_and_csv_names_inner=['UBERONBodyPartName',
                                                               'DOIDCommonName',
                                                               'HumanPopulationDensity',
                                                               'HealthStatus'],
                                   fill_col_from=[['UBERONBodyPartName', 'SAMPLE_TYPE']],
                                   case='outer',
                                   path_to_csvs='translation_sheets',
                                   add_to_cols=[],
                                   **kwargs):

    ontology_table = kwargs['ontology_table']
    allowedTerm_dict = kwargs['allowedTerm_dict']

    if case == 'outer':
        column_and_csv_names = column_and_csv_names_outer
    elif case == 'inner':
        column_and_csv_names = column_and_csv_names_inner + add_to_cols
        present_keys = MWB_table[['Key', 'Value']].copy()
        present_keys['Value'] = present_keys['Value'].apply(convert_to_numeric_or_original)
        present_keys.drop_duplicates(inplace=True)
    elif case == 'fill':
        column_and_csv_names = [x[0] for x in fill_col_from]
        origin_cols = [x[1] for x in fill_col_from]

    for index, col_csv_name in enumerate(column_and_csv_names):

        df_translations = pd.read_csv(path_to_csvs + "/{}.csv".format(str(col_csv_name)), encoding="ISO-8859-1",
                                      dtype=str)
        df_translations['MWB'] = df_translations['MWB'].str.lower()
        df_translations = df_translations.drop_duplicates()

        if case == 'outer':

            if col_csv_name == 'NCBITaxonomy' and 'NCBITaxonomy' in MWB_table.columns:
                continue

            MWB_table[col_csv_name] = MWB_table[col_csv_name].str.lower()
            MWB_table = pd.merge(MWB_table, df_translations, left_on=col_csv_name, right_on='MWB', how='left')
            MWB_table = MWB_table.drop(columns=['MWB', col_csv_name])
            MWB_table = MWB_table.rename(columns={'REDU': col_csv_name})
               
        if case == 'inner':

            if col_csv_name != 'NCBITaxonomy':
                MWB_table['match'] = MWB_table['Value'].isin(df_translations['MWB'].tolist())
                MWB_table = MWB_table.merge(df_translations, left_on='Value', right_on='MWB', how='left')
                MWB_table[col_csv_name] = MWB_table.groupby('filename')['REDU'].transform('first')
                if col_csv_name == 'UBERONBodyPartName':
                    MWB_table = pd.merge(MWB_table, ontology_table, left_on='UBERONBodyPartName', right_on='Label', how='left')
                    MWB_table = MWB_table.drop(columns=['REDU_UBERONOntologyIndex'])
                if col_csv_name == 'DOIDCommonName':
                    MWB_table['DOIDOntologyIndex'] = MWB_table.groupby('filename')['REDU_DOIDOntologyIndex'].transform('first')
                    MWB_table = MWB_table.drop(columns=['REDU_DOIDOntologyIndex'])

                MWB_table = MWB_table.drop(['MWB', 'match', 'REDU'], axis=1)

            elif col_csv_name == 'NCBITaxonomy':


                processed_species = {}
                none_count = 0
                unique_species = set(MWB_table['Value'])

                #processed_species = {species: get_taxonomy_id_from_name__allowedTerms(species, allowedTerm_dict=allowedTerm_dict) for species in unique_species} 
                for species in unique_species:
                    result = get_taxonomy_id_from_name__allowedTerms(species, allowedTerm_dict=allowedTerm_dict)
                    if result is None:
                        none_count += 1
                        if none_count >= 100:
                            processed_species = {}
                            break
                    else:
                        none_count = 0  # Reset count if a non-None result is found
                        processed_species[species] = result

                if len(processed_species) > 0:

                    MWB_table['REDU'] = MWB_table['Value'].map(processed_species)

                    def get_longest_value(series):
                        lengths = series.str.len()  # Calculate length of each string
                        if lengths.isnull().all():  # Check if all values are None
                            return None
                        return series.loc[lengths.idxmax()]
                    
                    #MWB_table['NCBITaxonomy'] = MWB_table.groupby('filename')['REDU'].transform('first')
                    MWB_table['NCBITaxonomy'] = MWB_table.groupby('filename')['REDU'].transform(get_longest_value)
                    MWB_table = MWB_table.drop(['REDU'], axis=1)

        if case == 'fill':
            MWB_table[origin_cols[index]] = MWB_table[origin_cols[index]].str.lower()
            MWB_table = pd.merge(MWB_table, df_translations, left_on=origin_cols[index], right_on='MWB', how='left')
            MWB_table['REDU'] = MWB_table['REDU'].fillna('ML import: not available')
            MWB_table[col_csv_name] = MWB_table[col_csv_name].fillna(MWB_table['REDU'])
            if col_csv_name == 'UBERONBodyPartName':
                MWB_table = MWB_table.drop(columns=['REDU_UBERONOntologyIndex'])
                MWB_table = MWB_table.drop(columns=['MWB', 'REDU', origin_cols[index]])


    return MWB_table


def translate_MWB_to_REDU_by_logic(MWB_table, path_to_csvs='translation_sheets'):
    if 'MWB_age' in MWB_table.columns:
        MWB_table['AgeInYears'] = MWB_table['MWB_age'].apply(convert_to_years)
    if 'MWB_sex' in MWB_table.columns:
        MWB_table['BiologicalSex'] = MWB_table['MWB_sex'].map(convert_sex)

    df_translations = pd.read_csv(path_to_csvs + "/biofluid_tissue_distinction.csv", encoding="ISO-8859-1", dtype=str)

    MWB_table = MWB_table.merge(df_translations, left_on='UBERONBodyPartName', right_on='sampletype', how='left')
    
    # Convert both columns to string type before filling NA values
    MWB_table['SampleTypeSub1'] = MWB_table['SampleTypeSub1'].astype(str)
    MWB_table['tissue_vs_biofluid'] = MWB_table['tissue_vs_biofluid'].astype(str)

    # Now you can fill NA values without worrying about datatype issues
    MWB_table['SampleTypeSub1'] = MWB_table['SampleTypeSub1'].fillna(MWB_table['tissue_vs_biofluid'])

    MWB_table = MWB_table.drop(columns=['sampletype', 'tissue_vs_biofluid'])


    # List of columns to check and add if not present
    columns_to_check = ['NCBITaxonomy', 'SampleCollectionMethod', 'SampleExtractionMethod',
                        'InternalStandardsUsed', 'SubjectIdentifierAsRecorded', 'AgeInYears',
                        'BiologicalSex', 'UBERONBodyPartName', 'HealthStatus', 'DOIDCommonName',
                        'ComorbidityListDOIDIndex', 'Country', 'HumanPopulationDensity',
                        'LatitudeandLongitude']

    # Iterate through each column in the list
    for column in columns_to_check:
        # Check if the column exists in the DataFrame
        if column not in MWB_table.columns:
            # If the column does not exist, create it with all values set to 'ML import: not available'
            MWB_table[column] = 'ML import: not available'
        else:
            # If the column exists, replace NaNs with 'ML import: not available'
            MWB_table[column] = MWB_table[column].replace(np.nan, 'ML import: not available')

    return MWB_table



def MWB_to_REDU_study_wrapper(study_id, path_to_csvs='translation_sheets',
                              duplicate_raw_file_handling='remove_duplicates', export_to_tsv = False, **kwargs):


    allowedTerm_dict = kwargs['allowedTerm_dict']
    ontology_table = kwargs['ontology_table']
    NCBIRankDivision_table = kwargs['NCBIRankDivision_table']

    raw_file_name_tupple = _get_metabolomicsworkbench_files(study_id)

    raw_file_name_df = pd.DataFrame(raw_file_name_tupple[0])
    if len(raw_file_name_df) == 0:
        print('No raw data detected. Ignoring {}'.format(str(study_id)))
        return None

    raw_file_name_df['filename_base'] = raw_file_name_df['filename'].apply(lambda x: os.path.basename(x))

    stdy_info_req = requests.get(
        'https://www.metabolomicsworkbench.org/rest/study/study_id/{}/analysis'.format(str(study_id)))
    stdy_info = stdy_info_req.json()

    if not isinstance(next(iter(stdy_info.values())), dict):
        stdy_info = {'1': stdy_info}

    redu_dfs = []
    for analysis_id, analysis_details in stdy_info.items():

        try:
            if analysis_details['analysis_type'] == 'MS':
                print('Initiating REDU-table generation for ' + str(analysis_details["analysis_id"]))
            else:
                print('This is not an MS-analysis. Ignoring ' + str(analysis_details["analysis_id"]))
                continue
        
        except Exception as e:
            traceback_info = traceback.format_exc()
            print(f"No analysis_type-key in {study_id}: {e}\nTraceback:\n{traceback_info}")
            continue

        
        redu_df = MWB_to_REDU_wrapper(MWB_analysis_ID=analysis_details["analysis_id"],
                                      raw_file_name_df=raw_file_name_df[['filename', 'filename_base', 'USI']],
                                      path_to_csvs=path_to_csvs,
                                      Massive_ID=study_id + '|' + str(analysis_details["analysis_id"]),
                                      allowedTerm_dict=allowedTerm_dict,
                                      ontology_table=ontology_table)
        
        if isinstance(redu_df, pd.DataFrame):
            redu_dfs.append(redu_df)

    
    if len(redu_dfs) > 0:
        redu_df_final = pd.concat(redu_dfs)
        
        
        if duplicate_raw_file_handling == 'remove_duplicates':
            duplicates = redu_df_final.filename.value_counts() > 1
            redu_df_final = redu_df_final[~redu_df_final.filename.isin(duplicates.index[duplicates])]

        if duplicate_raw_file_handling == 'keep_pols_dupl':
            mask = redu_df_final.duplicated(subset='filename', keep=False)

            remove_mask = mask.copy()
            for i, row in redu_df_final[mask].iterrows():
                if ((redu_df_final['filename'] == row['filename']) & mask).sum() == 2:

                    negative_mask = (redu_df_final['filename'] == row['filename']) & (
                        redu_df_final['IonizationSourceAndPolarity'].str.contains('negative'))
                    positive_mask = (redu_df_final['filename'] == row['filename']) & (
                        redu_df_final['IonizationSourceAndPolarity'].str.contains('positive'))

                    if negative_mask.sum() == 1 and positive_mask.sum() == 1:
                        remove_mask &= ~((redu_df_final['filename'] == row['filename']) & mask)

            # Remove all rows with entries occurring more than once in col1 (but not if one has 'negative' and the other has 'positive' in col2)
            redu_df_final = redu_df_final[~remove_mask]
        
    else:
        return None


    #ensure correct order of columns
    redu_df_final = redu_df_final[["MassiveID",
                                   "filename",
                                   "SampleType",
                                   "SampleTypeSub1",
                                   "NCBITaxonomy",
                                   "YearOfAnalysis",
                                   "SampleCollectionMethod",
                                   "SampleExtractionMethod",
                                   "InternalStandardsUsed",
                                   "MassSpectrometer",
                                   "IonizationSourceAndPolarity",
                                   "ChromatographyAndPhase",
                                   "SubjectIdentifierAsRecorded",
                                   "AgeInYears",
                                   "BiologicalSex",
                                   "UBERONBodyPartName",
                                   "TermsofPosition",
                                   "HealthStatus",
                                   "DOIDCommonName",
                                   "ComorbidityListDOIDIndex",
                                   "SampleCollectionDateandTime",
                                   "Country",
                                   "HumanPopulationDensity",
                                   "LatitudeandLongitude",
                                   "DepthorAltitudeMeters",
                                   "qiita_sample_name",
                                   "UniqueSubjectID",
                                   "LifeStage",
                                   "UBERONOntologyIndex",
                                   "DOIDOntologyIndex",
                                   "USI"
                                   ]]


    redu_df_final = merge_repeated_fileobservations_across_mwatb(redu_df_final, polarity_table=polarity_table)

    ontology_table = ontology_table.drop_duplicates(subset=['Label'])
    redu_df_final = complete_and_fill_REDU_table(redu_df_final, allowedTerm_dict, UBERONOntologyIndex_table=ontology_table, NCBIRankDivision_table=NCBIRankDivision_table,
                                                 add_usi = True, other_allowed_file_extensions = ['.raw', '.cdf', '.wiff', '.d'])


    if export_to_tsv == True:
        redu_df_final.to_csv('{}_REDU_from_MWB.tsv'.format(study_id), sep='\t', index=False, header=True)
        return None

    return redu_df_final


def MWB_to_REDU_wrapper(mwTab_json=None, MWB_analysis_ID=None, raw_file_name_df=None, Massive_ID='',
                        path_to_csvs='translation_sheets', **kwargs):
    if mwTab_json is None and MWB_analysis_ID is None:
        raise SystemExit("mwTab_json or MWB_analysis_ID has to be provided!")

    if mwTab_json is not None and MWB_analysis_ID is not None:
        raise SystemExit("Only mwTab_json OR MWB_analysis_ID should be provided!")

    if MWB_analysis_ID is not None:
        mwb_anlysis_req = requests.get(
            "https://www.metabolomicsworkbench.org/rest/study/analysis_id/{}/mwtab".format(str(MWB_analysis_ID)))

        try:
            mwTab_json = json.loads(mwb_anlysis_req.text, object_pairs_hook=handle_duplicates)
        except json.JSONDecodeError:
            print("Did not receive valid mwTab json for {}!".format(str(MWB_analysis_ID)))
            return None
        
    allowedTerm_dict = kwargs['allowedTerm_dict']
    ontology_table = kwargs['ontology_table']

    try:

        complete_df = translate_MWB_to_REDU_by_logic(
            translate_MWB_to_REDU_from_csv(
                translate_MWB_to_REDU_from_csv(
                    create_dataframe_outer_dict(mwTab_json, raw_file_name_df=raw_file_name_df, path_to_csvs=path_to_csvs, allowedTerm_dict=allowedTerm_dict, ontology_table=ontology_table),
                    path_to_csvs=path_to_csvs,
                ontology_table=ontology_table,
                allowedTerm_dict=allowedTerm_dict),
                case='fill',
                path_to_csvs=path_to_csvs,
                ontology_table=ontology_table,
                allowedTerm_dict=allowedTerm_dict),
            path_to_csvs=path_to_csvs)

    except Exception as e:
        traceback_info = traceback.format_exc()
        print(f"An error occurred in MWB_to_REDU_wrapper: {e}\nTraceback:\n{traceback_info}")

        return None
    
    

    complete_df['filename'] = complete_df['filename_raw_path']
    complete_df = complete_df[complete_df['filename'] != ''].dropna(subset=['filename'])


    complete_df['MassiveID'] = Massive_ID
    complete_df['MassiveID'] = complete_df['MassiveID'].apply(lambda x: x.split('|')[0])

    complete_df['UniqueSubjectID'] = complete_df.apply(
        lambda x: str(x['MassiveID']) + '_' + str(x['SubjectIdentifierAsRecorded']) 
        if x['SubjectIdentifierAsRecorded'] != '' else None,
        axis=1)
    
    

    complete_df['LifeStage'] = complete_df['AgeInYears'].apply(lambda x: age_category(x))

    complete_df[['SampleCollectionDateandTime',
                 'DepthorAltitudeMeters',
                 'qiita_sample_name']] = 'ML import: not available'

    missing_not_imported = ["SampleType",
                            "SampleTypeSub1",
                            "NCBITaxonomy",
                            "YearOfAnalysis",
                            "SampleCollectionMethod",
                            "SampleExtractionMethod",
                            "InternalStandardsUsed",
                            "MassSpectrometer",
                            "IonizationSourceAndPolarity",
                            "ChromatographyAndPhase",
                            "LifeStage",
                            "Country",
                            "UBERONOntologyIndex",
                            "DOIDOntologyIndex",
                            "SubjectIdentifierAsRecorded",
                            "BiologicalSex",
                            "UBERONBodyPartName",
                            "HealthStatus",
                            "DOIDCommonName",
                            "ComorbidityListDOIDIndex",
                            "SampleCollectionDateandTime",
                            "HumanPopulationDensity",
                            "LatitudeandLongitude",
                            "DepthorAltitudeMeters",
                            "qiita_sample_name"
                            ]
    
    complete_df[missing_not_imported] = complete_df[missing_not_imported].replace(
        {np.nan: 'ML import: not available', None: 'ML import: not available'})

    missing_not_collected = ["AgeInYears"]
    complete_df[missing_not_collected] = complete_df[missing_not_collected].replace(
        {np.nan: 'not collected', None: 'not collected'})

    complete_df[['TermsofPosition',
                 'SampleCollectionDateandTime',
                 'DepthorAltitudeMeters',
                 'qiita_sample_name']] = 'ML import: not available'

    REDU_df = complete_df[["MassiveID",
                           "filename",
                           "SampleType",
                           "SampleTypeSub1",
                           "NCBITaxonomy",
                           "YearOfAnalysis",
                           "SampleCollectionMethod",
                           "SampleExtractionMethod",
                           "InternalStandardsUsed",
                           "MassSpectrometer",
                           "IonizationSourceAndPolarity",
                           "ChromatographyAndPhase",
                           "SubjectIdentifierAsRecorded",
                           "AgeInYears",
                           "BiologicalSex",
                           "UBERONBodyPartName",
                           "TermsofPosition",
                           "HealthStatus",
                           "DOIDCommonName",
                           "ComorbidityListDOIDIndex",
                           "SampleCollectionDateandTime",
                           "Country",
                           "HumanPopulationDensity",
                           "LatitudeandLongitude",
                           "DepthorAltitudeMeters",
                           "qiita_sample_name",
                           "UniqueSubjectID",
                           "LifeStage",
                           "UBERONOntologyIndex",
                           "DOIDOntologyIndex",
                           "USI"
                           ]]
    
    return REDU_df


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Give an MWB study ID and get a REDU table tsv.')
    parser.add_argument("--study_id", "-mwb_id", type=str, help='An MWB study ID such as "ST002050". If "ALL" all study IDs are requested.', required=True)
    parser.add_argument("--path_to_csvs", "-csvs", type=str, help="Path to the translation csvs holding translations from MWB to REDU vocabulary (optional)", default="translation_sheets")
    parser.add_argument("--path_to_allowed_term_json", "-json", type=str, help="Allowed vocabulary json (optional)", default="allowed_terms")
    parser.add_argument("--path_ncbi_rank_division", type=str, help="Path to the path_ncbi_rank_division")
    parser.add_argument("--path_to_uberon_cl_po_csv", type=str, help="Path to the prepared uberon_cl_po ontology csv")
    parser.add_argument("--duplicate_raw_file_handling", "-duplStrat", type=str, help="What should be done with duplicate filenames across studies? Can be 'keep_pols_dupl' to keep cases where files can be distinguished by their polarity or 'remove_duplicates' to only keep cases where files can be assigned unambiguously (i.e. cases with only one analysis per study_id)(optional)", default='remove_duplicates')
    parser.add_argument("--path_to_polarity_info", type=str, help="Path to the polarity file.", default='none')

    #python3.8 workflows/PublicDataset_ReDU_Metadata_Workflow/bin/MWB_to_REDU.py --study_id  ST000107 --path_to_csvs /home/yasin/projects/ReDU-MS2-GNPS2/workflows/PublicDataset_ReDU_Metadata_Workflow/bin/translation_sheets --path_to_allowed_term_json /home/yasin/projects/ReDU-MS2-GNPS2/work/03/59b15db86cef31d8a337537d4cff04/allowed_terms.json --duplicate_raw_file_handling keep_all --path_to_uberon_cl_po_csv /home/yasin/projects/ReDU-MS2-GNPS2/work/ca/535593beabd1e7034883820a96e530/UBERON_CL_PO_ontology.csv --path_to_polarity_info /home/yasin/projects/ReDU-MS2-GNPS2/workflows/PublicDataset_ReDU_Metadata_Workflow/data/MWB_polarity_table.csv 
    print('Starting MWB2REDU script,..')

    args = parser.parse_args()

    study_id = args.study_id
    path_to_csvs = args.path_to_csvs
    duplicate_raw_file_handling = args.duplicate_raw_file_handling
    path_to_allowed_term_json = args.path_to_allowed_term_json

    # Read allowed terms json
    with open(path_to_allowed_term_json, 'r') as file:
        allowedTerm_dict = json.load(file)

    # Read ontology
    ontology_table = pd.read_csv(args.path_to_uberon_cl_po_csv)

    # Read polarity file
    polarity_table = pd.read_csv(args.path_to_polarity_info)
   
    # Read NCBI rank and division csv
    NCBIRankDivision_table = pd.read_csv(args.path_ncbi_rank_division, index_col = False)
    NCBIRankDivision_table = NCBIRankDivision_table.drop_duplicates(subset=['TaxonID'])

    # result
    if study_id == "ALL":
        # Getting all files
        url = "https://www.metabolomicsworkbench.org/rest/study/study_id/ST/available"
        studies_dict = requests.get(url).json()

        study_list = []
        for key in studies_dict.keys():
            study_dict = studies_dict[key]
            study_list.append(study_dict['study_id'])

        study_list = list(set(study_list))

        all_results_list = []
        for study_id in tqdm(study_list):
            print("Processing ", study_id)

            try:
                result = MWB_to_REDU_study_wrapper(study_id=study_id,
                                          path_to_csvs=path_to_csvs,
                                          duplicate_raw_file_handling=duplicate_raw_file_handling,
                                          export_to_tsv=False,
                                          allowedTerm_dict=allowedTerm_dict,
                                          ontology_table=ontology_table,
                                          polarity_table=polarity_table,
                                          NCBIRankDivision_table=NCBIRankDivision_table)
                print('Extracted information for {} samples.'.format(len(result)))
                if len(result) > 1:
                    all_results_list.append(result)
            except KeyboardInterrupt:
                raise
            except:
                pass
        merged_df = pd.concat(all_results_list, ignore_index=True)
        merged_df.to_csv('REDU_from_MWB_all.tsv', sep='\t', index=False, header=True)

    else:
        MWB_to_REDU_study_wrapper(study_id=study_id,
                                  path_to_csvs=path_to_csvs,
                                  duplicate_raw_file_handling=duplicate_raw_file_handling,
                                  export_to_tsv=True,
                                  allowedTerm_dict=allowedTerm_dict,
                                  ontology_table=ontology_table,
                                  polarity_table=polarity_table,
                                  NCBIRankDivision_table=NCBIRankDivision_table)

    print("Output files written to working directory")
