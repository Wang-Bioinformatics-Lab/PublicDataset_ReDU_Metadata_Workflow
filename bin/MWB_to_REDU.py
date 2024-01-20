import os
import argparse
import pandas as pd
import re
import numpy as np
import requests
from bs4 import BeautifulSoup
from urllib.parse import urlparse, parse_qs
import json
import time


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
        workbench_df["size_mb"] = (workbench_df["FILESIZE"].astype(int) / 1024 / 1024).astype(int)
        workbench_df["URL"] = workbench_df["URL"]        

        # Generate USI column
        workbench_df["USI"] = workbench_df["URL"].apply(
            lambda url: "mzspec:{}:{}-{}".format(
                dataset_accession,
                parse_qs(urlparse(url).query).get('A', [None])[0],
                parse_qs(urlparse(url).query).get('F', [None])[0]
            )
        )

        workbench_df = workbench_df[["filename", "size_mb", "USI"]]
    except:
        workbench_df = pd.DataFrame()

    merged_df = pd.concat([files_df, workbench_df])

    return merged_df, msv_accession


def age_category(age):
    if age is None:
        return ''
    try:
        age = float(age)
    except ValueError:
        return ''
    if age < 2:
        return 'Infancy (<2 yrs)'
    elif age <= 8:
        return 'Early Childhood (2 yrs < x <=8 yrs)'
    elif age <= 18:
        return 'Adolescence (8 yrs < x <= 18 yrs)'
    elif age <= 45:
        return 'Early Adulthood (18 yrs < x <= 45 yrs)'
    elif age <= 65:
        return 'Middle Adulthood (45 yrs < x <= 65 yrs)'
    elif age > 65:
        return 'Later Adulthood (>65 yrs)'
    else:
        return ''


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


def get_taxonomy_info(ncbi_id, cell_culture_key1, cell_culture_key2):
    if ncbi_id is not None and ncbi_id != "NA":
        cell_culture_key_words = ["cell", "media", "culture"]
        try:
            #try to get taxa via API
            for attempt in range(3):
                try:
                    response = requests.get(
                        f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=taxonomy&id={ncbi_id}")
                    soup = BeautifulSoup(response.text, "xml")
                    classification = [s.lower() for s in soup.find("Taxon").find("Lineage").text.split("; ")]
                    break

                except requests.exceptions.RequestException as e:
                    print(f"Taxonomy request failed. Retrying (attempt {attempt + 1}/3): {e}")
                    time.sleep(2)

            else:
                # If the loop completes without breaking, raise an exception
                raise Exception("All retries failed. Unable to fetch taxonomy data.")

            SampleType = None
            SampleTypeSub1 = None

            if "viridiplantae" in classification:
                SampleType = "plant"
                SampleTypeSub1 = "plant_NOS"
                if "Algae" in classification or 'rhodophyta' in classification or "phaeophyceae" in classification:
                    SampleType = "algae"
                if "chlorophyta" in classification or "microalgae" in classification or "microalga" in classification:
                    SampleType = "microalgae"
                if "streptophyta" in classification:
                    SampleTypeSub1 = "plant_angiospermae"
                if "cyanobacteria" in classification:
                    SampleTypeSub1 = "marine_cyanobacteria_insitu"
                if "bacillariophyta" in classification:
                    SampleTypeSub1 = "marine_diatom"
            elif "metazoa" in classification:
                SampleType = "animal"
                if "mammalia" in classification and any(
                        word in cell_culture_key1.lower() for word in cell_culture_key_words) or any(
                        word in cell_culture_key2.lower() for word in cell_culture_key_words):
                    SampleType = "culture_mammalian"
                    SampleTypeSub1 = "culture_mammalian"
                if "amphibia" in classification:
                    if "Caudata" in classification or "urodela" in classification or "echinodermata" in classification:
                        SampleTypeSub1 = "salamander"
                    else:
                        SampleTypeSub1 = "frog"
                if "insecta" in classification:
                    SampleTypeSub1 = "insect"
                if "porifera" in classification or "mollusca" in classification:
                    SampleTypeSub1 = "marine_invertebrates"
                if "cnidaria" in classification:
                    SampleTypeSub1 = "marine_coral"
            elif "fungi" in classification:
                SampleType = "culture_fungal"
                SampleTypeSub1 = "culture_fungal"
            elif "bacteria" in classification:
                SampleType = "culture_bacterial"
                SampleTypeSub1 = "culture_bacterial"
            return [SampleType, SampleTypeSub1]
        except Exception:
            return [None, None]
    else:
        return [None, None]


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
            if years >= 2000 and years <= 2030:
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
        df[new_col].fillna(value='NA', inplace=True)
    return df

def get_rawFile_names(df, key_vars, new_col):
    # Create a pattern to match key_vars and their numbered versions
    pattern = r'(' + '|'.join(key_vars) + r')(_\d+)?'

    # Extract the base key name without the number
    df['base_key'] = df['Key'].str.extract(pattern, expand=False)[0]

    # Filter the DataFrame to get rows where 'base_key' is not null (matches the pattern)
    key_df = df[df['base_key'].notnull()].copy()

    # Create a new DataFrame to hold the expanded rows
    expanded_rows = []

    # Iterate over each filename group
    for filename, group in key_df.groupby('filename'):
        # Iterate over each base_key within the group
        for base_key in group['base_key'].unique():
            # Create a row for each value associated with the base_key
            for value in group[group['base_key'] == base_key]['Value']:
                new_row = df[df['filename'] == filename].iloc[0].to_dict()
                new_row[new_col] = value  # Add the new column
                expanded_rows.append(new_row)

    # Create a new DataFrame from the expanded rows
    expanded_df = pd.DataFrame(expanded_rows)

    # Drop duplicate rows, if any
    expanded_df.drop_duplicates(inplace=True)

    # Handle numeric conversion for specific columns
    if new_col in ['Latitude', 'Longitude']:
        expanded_df[new_col] = pd.to_numeric(expanded_df[new_col], errors='coerce')
    else:
        expanded_df[new_col].fillna(value='NA', inplace=True)

    # Drop the temporary 'base_key' column
    expanded_df.drop('base_key', axis=1, inplace=True)

    return expanded_df

def create_dataframe_from_SUBJECT_SAMPLE_FACTORS(data, raw_file_name_df=None, path_to_csvs='translation_sheets'):
    data_list = []
    for item in data:
        subject_id = item.get('Subject ID', 'NA')
        sample_id = item.get('Sample ID', 'NA')
        factors = item.get('Factors', {})
        additional_data = item.get('Additional sample data', {})

        for key, value in factors.items():
            data_list.append([subject_id, sample_id, key, value])
        for key, value in additional_data.items():
            data_list.append([subject_id, sample_id, key, value])
    df = pd.DataFrame(data_list, columns=['SubjectIdentifierAsRecorded', 'filename', 'Key', 'Value'])
    df['Key'] = df['Key'].str.lower()


    #raw data name might be in the Sample ID, Factors
    expected_raw_file_keys = ["raw_file_name",  "rawfilename", "raw_file",
                              "datafile name", "raw files"]


    if not any(key in df['Key'].values for key in expected_raw_file_keys):
        df[['filename_raw', 'filename_raw_path']] = df.apply(lambda x: get_raw_file_names(x['Value'], x['filename'], raw_file_name_df),
                                                             axis=1, result_type='expand')
        df['filename_raw'] = df.groupby('filename')['filename_raw'].transform('first')
        df['filename_raw_path'] = df.groupby('filename')['filename_raw_path'].transform('first')
    else:
        raw_file_name_df = raw_file_name_df.rename(columns={'filename': 'filename_raw_path'})
        df = get_rawFile_names(df, key_vars=expected_raw_file_keys, new_col="filename_raw")

        if bool(set(list(df['filename_raw'])) & set(list(raw_file_name_df['filename_base']))):
            df = df.merge(raw_file_name_df, left_on='filename_raw', right_on='filename_base', how='left')
        else:
            raw_file_name_df['filename_base_wo_extension'] = raw_file_name_df['filename_base'].str.split('.').str[0]
            df = df.merge(raw_file_name_df, left_on='filename_raw', right_on='filename_base_wo_extension', how='left')
            df['filename_raw'] = df['filename_base']
            df = df.drop(columns=['filename_base'])

    df['Value'] = df['Value'].str.lower()

    df['SubjectIdentifierAsRecorded'] = df['SubjectIdentifierAsRecorded'].replace('-', '')
    df = get_key_info_into_outer(df, key_vars=["gender", "sex"], new_col="MWB_sex")
    df = get_key_info_into_outer(df, key_vars=["age", "age (years)"], new_col="MWB_age")
    df = get_key_info_into_outer(df, key_vars=["collection_country", "collection country", "country"],
                                 new_col="Country")
    df = get_key_info_into_outer(df, key_vars=["latitude"], new_col="Latitude")
    df = get_key_info_into_outer(df, key_vars=["longitude"], new_col="Longitude")
    df['LatitudeandLongitude'] = df.apply(
        lambda x: f'{x.Latitude}|{x.Longitude}' if x.Latitude is not None and x.Longitude is not None and not np.isnan(
            x.Latitude) and not np.isnan(x.Longitude) else None, axis=1)

    df[['SampleType_inner', 'SampleTypeSub1_inner']] = df.apply(lambda x: get_blanks(x['Value']), axis=1,
                                                                result_type='expand')

    df = translate_MWB_to_REDU_from_csv(df, case='inner', path_to_csvs=path_to_csvs)
    df = df.drop(columns=['Key', 'Value', 'Longitude', 'Latitude'])
    df = df.drop_duplicates().reset_index(drop=True)

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


def create_dataframe_outer_dict(MWB_mwTAB_dict, raw_file_name_df=None, path_to_csvs = 'translation_sheets'):
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

    if entry['ANALYSIS_TYPE'] != 'MS':
        raise ValueError("This is not an MS analysis!")

    df_outer = pd.DataFrame(entry, index=[0])

    if entry['TAXONOMY_ID'] == 'NA' and entry['SUBJECT_SPECIES'] != 'NA':
        df_outer['TAXONOMY_ID'] = df_outer['SUBJECT_SPECIES'].apply(lambda x: get_taxonomy_id_from_name(x))
    df_outer['NCBITaxonomy'] = df_outer['TAXONOMY_ID'] + '|' + df_outer['SUBJECT_SPECIES']
    df_outer['ChromatographyAndPhase'] = df_outer['CHROMATOGRAPHY_TYPE'] + '|' + df_outer['COLUMN_NAME']
    df_outer['YearOfAnalysis'] = df_outer['ACQUISITION_DATE'].apply(lambda x: extract_years(x))
    df_outer['YearOfAnalysis'] = df_outer.apply(
        lambda x: extract_years(x['CREATED_ON']) if x['YearOfAnalysis'] is None else x['YearOfAnalysis'], axis=1)

    if df_outer['YearOfAnalysis'][0] is None:
        raise ValueError("YearOfAnalysis was not collected!")

    df_outer['ChromatographyAndPhase'] = df_outer['CHROMATOGRAPHY_TYPE'] + '|' + df_outer['COLUMN_NAME']
    df_outer['IonizationSourceAndPolarity'] = df_outer['MS_TYPE'] + '|' + df_outer['ION_MODE']

    df_outer['SampleType'] = None
    df_outer['SampleTypeSub1'] = None
    df_outer[['SampleType', 'SampleTypeSub1']] = df_outer.apply(
        lambda x: get_taxonomy_info(x['TAXONOMY_ID'], x['SAMPLE_TYPE'], x['SUBJECT_TYPE']), axis=1,
        result_type='expand')
    df_outer[['SampleType', 'SampleTypeSub1']] = df_outer.apply(
        lambda x: get_enviromental_water(x['SAMPLE_TYPE'] + x['SUBJECT_TYPE']) if x['SampleType'] is None and x[
            'SampleTypeSub1'] is None else (x['SampleType'], x['SampleTypeSub1']), axis=1, result_type='expand')

    df_inner_SUBJECT_SAMPLE_FACTORS = create_dataframe_from_SUBJECT_SAMPLE_FACTORS(
        MWB_mwTAB_dict['SUBJECT_SAMPLE_FACTORS'],
        raw_file_name_df=raw_file_name_df,
        path_to_csvs=path_to_csvs
    )

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
                                                               'ChromatographyAndPhase',
                                                               'SampleCollectionMethod',
                                                               'IonizationSourceAndPolarity',
                                                               'Country'],
                                   column_and_csv_names_inner=['UBERONBodyPartName',
                                                               'DOIDCommonName',
                                                               'HumanPopulationDensity',
                                                               'HealthStatus'],
                                   fill_col_from=[['UBERONBodyPartName', 'SAMPLE_TYPE']],
                                   case='outer',
                                   path_to_csvs='translation_sheets'):


    if case == 'outer':
        column_and_csv_names = column_and_csv_names_outer
    elif case == 'inner':
        column_and_csv_names = column_and_csv_names_inner
    elif case == 'fill':
        column_and_csv_names = [x[0] for x in fill_col_from]
        origin_cols = [x[1] for x in fill_col_from]

    for index, col_csv_name in enumerate(column_and_csv_names):

        df_translations = pd.read_csv(path_to_csvs + "/{}.csv".format(str(col_csv_name)), encoding="ISO-8859-1",
                                      dtype=str)
        df_translations['MWB'] = df_translations['MWB'].str.lower()
        df_translations = df_translations.drop_duplicates()

        if case == 'outer':
            MWB_table[col_csv_name] = MWB_table[col_csv_name].str.lower()
            MWB_table = pd.merge(MWB_table, df_translations, left_on=col_csv_name, right_on='MWB', how='left')
            MWB_table = MWB_table.drop(columns=['MWB', col_csv_name])
            MWB_table = MWB_table.rename(columns={'REDU': col_csv_name})

        if case == 'inner':
            df_translations['MWB'] = df_translations['MWB'].str.lower()
            MWB_table['match'] = MWB_table['Value'].isin(df_translations['MWB'].tolist())
            MWB_table = MWB_table.merge(df_translations, left_on='Value', right_on='MWB', how='left')
            MWB_table[col_csv_name] = MWB_table.groupby('filename')['REDU'].transform('first')
            if col_csv_name == 'UBERONBodyPartName':
                MWB_table['UBERONOntologyIndex'] = MWB_table.groupby('filename')['REDU_UBERONOntologyIndex'].transform(
                    'first')
                MWB_table = MWB_table.drop(columns=['REDU_UBERONOntologyIndex'])
            if col_csv_name == 'DOIDCommonName':
                MWB_table['DOIDOntologyIndex'] = MWB_table.groupby('filename')['REDU_DOIDOntologyIndex'].transform(
                    'first')
                MWB_table = MWB_table.drop(columns=['REDU_DOIDOntologyIndex'])
            MWB_table = MWB_table.drop(['MWB', 'match', 'REDU'], axis=1)

        if case == 'fill':
            MWB_table[origin_cols[index]] = MWB_table[origin_cols[index]].str.lower()
            MWB_table = pd.merge(MWB_table, df_translations, left_on=origin_cols[index], right_on='MWB', how='left')
            MWB_table['REDU'] = MWB_table['REDU'].fillna('ML import: not available')
            MWB_table[col_csv_name] = MWB_table[col_csv_name].fillna(MWB_table['REDU'])
            if col_csv_name == 'UBERONBodyPartName':
                MWB_table['UBERONOntologyIndex'] = MWB_table['UBERONOntologyIndex'].fillna(
                    MWB_table['REDU_UBERONOntologyIndex'])
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
                              duplicate_raw_file_handling='remove_duplicates', export_to_tsv = False):

    raw_file_name_tupple = _get_metabolomicsworkbench_files(study_id)
    raw_file_name_df = pd.DataFrame(raw_file_name_tupple[0])

    if len(raw_file_name_df) == 0:
        print('No raw data detected. Ignoring {}'.format(str(study_id)))
        return None

    raw_file_name_df['filename_base'] = raw_file_name_df['filename'].apply(lambda x: os.path.basename(x))

    if raw_file_name_df['filename_base'].duplicated().any():
        print('Duplicate raw file names present in study {}. (Could attempt to read analysis_id from path)'.format(
            str(study_id)))
        return None

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
            print('No analysis_type-key in ' + study_id)
            continue
        
        redu_df = MWB_to_REDU_wrapper(MWB_analysis_ID=analysis_details["analysis_id"],
                                      raw_file_name_df=raw_file_name_df[['filename', 'filename_base']],
                                      path_to_csvs=path_to_csvs,
                                      Massive_ID = study_id + '|' + str(analysis_details["analysis_id"]))

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

    #add USI and correct filename (at this point filenames in redu table could be with our without extension)
    raw_file_name_df['filename_match'] = raw_file_name_df['filename'].apply(lambda x: os.path.splitext(os.path.basename(x))[0])
    raw_file_name_df['filename'] = raw_file_name_df['filename'].apply(lambda x: os.path.basename(x))

    redu_df_final['filename_match'] = redu_df_final['filename'].apply(lambda x: os.path.splitext(os.path.basename(x))[0])
    redu_df_final = redu_df_final.drop(columns=['filename'])

    merged_df = pd.merge(redu_df_final, raw_file_name_df[['filename_match', 'filename', 'USI']],
                         left_on='filename_match', right_on='filename_match', how='inner')

    redu_df_final = merged_df.drop(columns=['filename_match'])

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


    if export_to_tsv == True:
        redu_df_final.to_csv('{}_REDU_from_MWB.tsv'.format(study_id), sep='\t', index=False, header=True)
        return None

    return redu_df_final


def MWB_to_REDU_wrapper(mwTab_json=None, MWB_analysis_ID=None, raw_file_name_df=None, Massive_ID='',
                        path_to_csvs='translation_sheets'):
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

    try:

        complete_df = translate_MWB_to_REDU_by_logic(
            translate_MWB_to_REDU_from_csv(
                translate_MWB_to_REDU_from_csv(
                    create_dataframe_outer_dict(mwTab_json, raw_file_name_df=raw_file_name_df, path_to_csvs=path_to_csvs),
                    path_to_csvs=path_to_csvs),
                case='fill',
                path_to_csvs=path_to_csvs),
            path_to_csvs=path_to_csvs)

    except Exception as e:
        print(f"Error: {e}")
        return None

    # remove cases where filenames_raw appear more than once
    duplicates = complete_df.filename_raw.value_counts() > 1
    complete_df = complete_df[~complete_df.filename_raw.isin(duplicates.index[duplicates])]

    complete_df['filename'] = complete_df['filename_raw']

    complete_df = complete_df[complete_df['filename'] != ''].dropna(subset=['filename'])


    complete_df['MassiveID'] = Massive_ID
    complete_df['UniqueSubjectID'] = complete_df.apply(
        lambda x: str(x['MassiveID']) + '_' + str(x['SubjectIdentifierAsRecorded']) if x[
                                                                                           'SubjectIdentifierAsRecorded'] != '' else None,
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
                           #'filename_raw_path',
                           #'STUDY_ID',
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
                           "DOIDOntologyIndex"
                           ]]

    return REDU_df


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Give an MWB study ID and get a REDU table tsv.')
    parser.add_argument("--study_id", "-mwb_id", type=str, help='An MWB study ID such as "ST002050". If "ALL" all study IDs are requested.', required=True)
    parser.add_argument("--path_to_csvs", "-csvs", type=str, help="Path to the translation csvs holding translations from MWB to REDU vocabulary (optional)", default="translation_sheets")
    parser.add_argument("--duplicate_raw_file_handling", "-duplStrat", type=str, help="What should be done with duplicate filenames across studies? Can be 'keep_pols_dupl' to keep cases where files can be distinguished by their polarity or 'remove_duplicates' to only keep cases where files can be assigned unambiguously (i.e. cases with only one analysis per study_id)(optional)", default='remove_duplicates')

    args = parser.parse_args()

    study_id = args.study_id
    path_to_csvs = args.path_to_csvs
    duplicate_raw_file_handling = args.duplicate_raw_file_handling

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
        for study_id in study_list:
            print("Processing ", study_id)

            try:
                result = MWB_to_REDU_study_wrapper(study_id=study_id,
                                          path_to_csvs=path_to_csvs,
                                          duplicate_raw_file_handling=duplicate_raw_file_handling,
                                          export_to_tsv=False)
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
                                  export_to_tsv=True)

    print("Output files written to working directory")
