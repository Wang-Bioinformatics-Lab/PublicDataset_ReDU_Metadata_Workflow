import os
import requests
import pandas as pd
from bs4 import BeautifulSoup
import argparse
import json
import time
import numpy as np
from tqdm import tqdm
from extend_allowed_terms import adapt_allowed_terms 
from REDU_conversion_functions import get_taxonomy_id_from_name__allowedTerms
from REDU_conversion_functions import update_unassigned_terms
from REDU_conversion_functions import get_taxonomy_id_from_name
import traceback


def complete_and_fill_REDU_table(df, allowedTerm_dict):
    """
    Completes and fills a REDU table with values based on a dictionary of allowed terms and missing values.

    Args:
    df: A pandas DataFrame containing the initial data.
    allowedTerm_dict: A dictionary containing allowed terms and missing values for each column.

    Returns:
    A DataFrame that has been filled with default values for missing columns,
    with values replaced by the corresponding "missing" value from the dictionary 
    if they are not in the allowed terms or are missing/empty, except for specific columns.
    """
    # Convert year to string for comparison
    if 'YearOfAnalysis' in df.columns:
        df['YearOfAnalysis'] = df['YearOfAnalysis'].astype(str)

    # Add missing columns with their respective default missing value from the dictionary
    for key, value in allowedTerm_dict.items():
        if key not in df.columns:
            df[key] = value['missing']

    # Replace values with the respective "missing" value if they're not in the allowed terms or are missing/empty
    for key, value in allowedTerm_dict.items():
        allowed_terms = value['allowed_values']
        missing_value = value['missing']
        if key in df.columns:
            if key not in ["MassiveID", "filename", "AgeInYears"]:
                df[key] = df[key].apply(lambda x: x if x in allowed_terms else missing_value).fillna(missing_value).replace("", missing_value)
            else:
                df[key] = df[key].fillna(missing_value).replace("", missing_value)

    # Ensure the dataframe contains only the columns specified in the dictionary plus 'USI'
    return df[[*allowedTerm_dict.keys(), 'USI']]

def find_column_after_target_column(df, target_column='', search_column_prefix='Samples_Unit'):
    if target_column in df.columns:
        target_index = df.columns.get_loc(target_column)
        # Check the next 3 columns after the target column
        for i in range(1, 4):
            if target_index + i < len(df.columns):
                next_column = df.columns[target_index + i]
                # Check if the next column starts with the search_column_prefix
                if next_column.startswith(search_column_prefix):
                    return next_column
    return ''  # Return None if no matching column is found


def get_taxonomy_id_from_name__allowedTerms(organism_name, **kwargs):

    allowedTerm_dict = kwargs['allowedTerm_dict']
    unassigned_term_json = kwargs['unassigned_term_json']


    taxonomy_data = allowedTerm_dict["NCBITaxonomy"]["allowed_values"]
    for entry in taxonomy_data:
        parts = entry.split('|')
        if len(parts) == 2:
            ncbi_id, name = parts
            if name.lower() == organism_name.lower():
                return ncbi_id + '|' + str(name)
        else:
            continue

    req_ncbi_name = get_taxonomy_id_from_name(organism_name)
    autoupdate = False


    if req_ncbi_name is not None and unassigned_term_json != 'none':

        print([req_ncbi_name + '__AUTOUPDATE'])

        path, ext = os.path.splitext(allowedTerm_dict)

        # Add the suffix to the path and recombine with the extension
        allowedTerm_dict_AUTOUPDATE = f"{path}_AUTOUPDATE.json"

        allowedTerm_dict = adapt_allowed_terms(terms_dict = allowedTerm_dict, 
                                               redu_variable = 'NCBITaxonomy', 
                                               term_list = [req_ncbi_name + '__AUTOUPDATE'], 
                                               add_or_remove = 'add', 
                                               load_dict_from_path = 'no', 
                                               save_dict_to_path = allowedTerm_dict_AUTOUPDATE)

        autoupdate = True

    if unassigned_term_json != 'none':
        update_unassigned_terms(organism_name, autoupdated=autoupdate, unassigned_file=unassigned_term_json)
    return None

def update_unassigned_terms(organism_name, autoupdated=False, unassigned_file='unassigned_terms.json'):
    if os.path.exists(unassigned_file):
        with open(unassigned_file, 'r') as file:
            unassigned_data = json.load(file)
    else:
        unassigned_data = {"Samples_Organism": {}}
    
    if organism_name in unassigned_data["Samples_Organism"]:
        unassigned_data["Samples_Organism"][organism_name]["count"] += 1
    else:
        unassigned_data["Samples_Organism"][organism_name] = {"count": 1, "autoupdated": autoupdated}

    with open(unassigned_file, 'w') as file:
        json.dump(unassigned_data, file, indent=4)


def get_taxonomy_info(ncbi_id, cell_culture_key1 = '', cell_culture_key2 = ''):
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
                    time.sleep(10)

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

def get_enviromental_water(x):
    x = x.lower()

    if 'water' in x or 'sewerage' in x:
        if 'waste' in x or 'sewerage' in x or 'sewage' in x:
            return ['environmental', 'water_waste']
        if 'surface' in x:
            return ['environmental', 'water_surface']
        if 'ground' in x:
            return ['environmental', 'water_ground']
        if 'storm' in x:
            return ['environmental', 'water_storm']
        if 'sea' in x or 'ocean' in x or 'coast' in x:
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


def get_taxonomy_id_from_name(species_name, retries=3):
    if species_name is None or species_name in ["NA", "N/A"]:
        return None

    attempts = 0
    while attempts < retries:
        try:
            response = requests.get(
                f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=taxonomy&term={species_name}&retmode=xml"
            )

            if response.status_code == 200:
                soup = BeautifulSoup(response.text, "xml")
                id_list = soup.find("IdList")
                if id_list is not None and id_list.find("Id") is not None:
                    ncbi_id = id_list.find("Id").text
                    term = soup.find("Term").text.split("[")[0].strip()
                    if term.lower() == species_name.lower():
                        return str(ncbi_id) + '|' + species_name
                    else:
                        return None
                else:
                    return  None
            else:
                print(f"Server responded with status code {response.status_code} for {species_name}.")
                time.sleep(10)

        except Exception as e:
            print(f"Attempt {attempts + 1} failed for {species_name}: {e}")
            time.sleep(10)

        attempts += 1

    print(f'{species_name} returned no NCBI-ID after {retries} attempts')
    return  None

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

def rename_duplicated_column_names(df):
    # We need this to get rid of duplicated columns
    new_column_names = []

    # Dictionary to keep track of the count of duplicate column names
    dup_count = {}

    # Iterate through the columns using index
    for idx, col in enumerate(df.columns):
        if col in df.columns[:idx]:
            # If the column name is a duplicate, append a count to its name
            dup_count[col] = dup_count.get(col, 0) + 1
            new_name = f"{col}_{dup_count[col]}"
            new_column_names.append(new_name)
        else:
            # If not a duplicate, keep the original name
            new_column_names.append(col)

    # Set the new column names to the DataFrame
    df.columns = new_column_names

    return df


def Metabolights2REDU(study_id, **kwargs):
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

    study_assays = study_details['content']['assays']

    submissionYear = study_details['content']['derivedData']['submissionYear']

    df_assays = pd.DataFrame()
    ms_study_assays = []

    #get study assays if they are MS. There can be multiple assay tables in the same study
    for index, assay_json in enumerate(study_assays):
        if assay_json['technology'] == 'mass spectrometry':          

            # Extract headers in the correct order
            headers = [None] * len(assay_json['assayTable']['fields'])
            for key, value in assay_json['assayTable']['fields'].items():
                headers[value['index']] = 'Assay_' + value['header']

            df = pd.DataFrame(assay_json['assayTable']['data'])
            df = rename_duplicated_column_names(df)

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
            df = rename_duplicated_column_names(df)
            aligned_dfs.append(df[list(all_columns)])

        # Concatenate all dataframes into a single dataframe
        df_assays = pd.concat(aligned_dfs, ignore_index=True)
        
        #get info on samples (there is only one sample table per study)
        headers = [None] * len(study_details['content']['sampleTable']['fields'])
        for key, value in study_details['content']['sampleTable']['fields'].items():
            headers[value['index']] = 'Samples_' + value['header']

        df_samples = pd.DataFrame(study_details['content']['sampleTable']['data'])
        df_samples = rename_duplicated_column_names(df_samples)
        df_samples.columns = headers
    
        df_study = df_assays.merge(df_samples, left_on='Assay_Sample Name', right_on='Samples_Sample Name', how='inner')

        # Duplicate rows if we have mzml AND raw files
        df_study_raw = pd.DataFrame()
        
        if 'Assay_Raw Spectral Data File' in df_study.columns:
            df_study_raw = df_study[df_study['Assay_Raw Spectral Data File'].str.contains('\.', regex=True, na=False)].copy()
            df_study_raw['filepath'] = df_study_raw['Assay_Raw Spectral Data File']
            # Check if the column exists before dropping
            if 'Assay_Raw Spectral Data File' in df_study_raw.columns:
                df_study_raw.drop(columns=['Assay_Raw Spectral Data File'], inplace=True)
            if 'Assay_Derived Spectral Data File' in df_study_raw.columns:
                df_study_raw.drop(columns=['Assay_Derived Spectral Data File'], inplace=True)

        df_study_mzml = pd.DataFrame()
        if 'Assay_Derived Spectral Data File' in df_study.columns:
            df_study_mzml = df_study[df_study['Assay_Derived Spectral Data File'].str.contains('\.', regex=True, na=False)].copy()
            df_study_mzml['filepath'] = df_study_mzml['Assay_Derived Spectral Data File']
            # Check if the column exists before dropping
            if 'Assay_Raw Spectral Data File' in df_study_mzml.columns:
                df_study_mzml.drop(columns=['Assay_Raw Spectral Data File'], inplace=True)
            if 'Assay_Derived Spectral Data File' in df_study_mzml.columns:
                df_study_mzml.drop(columns=['Assay_Derived Spectral Data File'], inplace=True)

        if len(df_study_raw) > 0 and len(df_study_mzml) > 0:
            df_study = pd.concat([df_study_mzml, df_study_raw], ignore_index=True).copy(deep=True)
        elif len(df_study_raw) > 0:
            df_study = df_study_raw.copy(deep=True)
        elif len(df_study_mzml) > 0:
            df_study = df_study_mzml.copy(deep=True)
        elif len(df_study_raw) == 0 and len(df_study_mzml) == 0:
            return pd.DataFrame()
        

        df_study = rename_duplicated_column_names(df_study)
        df_study['filename'] = df_study['filepath'].apply(lambda x: os.path.basename(x))

        # List of allowed extensions
        allowed_extensions = [".mzml", ".mzxml", ".cdf", ".raw", ".wiff", ".d"]

        # Use a tuple of allowed extensions for the endswith method
        df_study = df_study[df_study['filepath'].str.lower().str.endswith(tuple(allowed_extensions))]
        

        if len(df_study) > 0:
            
            allowedTerm_dict = kwargs['allowedTerm_dict']
            unassigned_term_json = kwargs['unassigned_term_json']


            df_study.loc[:, 'YearOfAnalysis'] = submissionYear

            #add NCBITaxonomy and Sampletype & SampleTypeSub1
            #######
            if 'Samples_Organism' in df_study.columns:
                #processed_organisms = {org: str(get_taxonomy_id_from_name(org)) + '|' + str(org) for org in df_study['Samples_Organism'].unique()}
                processed_organisms = {org: str(get_taxonomy_id_from_name__allowedTerms(org, allowedTerm_dict = allowedTerm_dict, unassigned_term_json=unassigned_term_json)) for org in df_study['Samples_Organism'].unique()}

                df_study.loc[:, 'NCBITaxonomy'] = df_study['Samples_Organism'].map(processed_organisms)
                df_study.loc[:, 'NCBITaxonomy'] = df_study['NCBITaxonomy'].replace(to_replace=r'^.*None.*$', value='ML import: not available', regex=True)
                
                df_study.loc[:, ['SampleType', 'SampleTypeSub1']] = 'ML import: not available'

                processed_taxonomy = {taxonomy.split('|')[0]: get_taxonomy_info(taxonomy.split('|')[0])
                                    for taxonomy in df_study['NCBITaxonomy'].unique()
                                    if '|' in taxonomy and 'None' not in taxonomy}


                df_study.loc[:,'SampleType'] = df_study['NCBITaxonomy'].map(lambda x: processed_taxonomy.get(x.split('|')[0], [pd.NA, pd.NA])[0] 
                                                                    if '|' in x and 'None' not in x else pd.NA)
                df_study.loc[:,'SampleTypeSub1'] = df_study['NCBITaxonomy'].map(lambda x: processed_taxonomy.get(x.split('|')[0], [pd.NA, pd.NA])[1]
                                                                        if '|' in x and 'None' not in x else pd.NA)

                df_study[['SampleType', 'SampleTypeSub1']] = df_study.apply(
                    lambda row: get_blanks(row.Samples_Organism) 
                    if (pd.isna(row.SampleType) or row.SampleType == 'ML import: not available') 
                    else [row.SampleType, row.SampleTypeSub1], 
                    axis=1
                ).apply(pd.Series)

                df_study[['SampleType', 'SampleTypeSub1']] = df_study.apply(
                    lambda row: get_enviromental_water(row.Samples_Organism) 
                    if (pd.isna(row.SampleType) or row.SampleType == 'ML import: not available') 
                    else [row.SampleType, row.SampleTypeSub1], 
                    axis=1
                ).apply(pd.Series)


                df_biofluid_vs_tissue = pd.read_csv(os.path.join(transSheet_dir, 'biofluid_vs_tissue_Metabolights.csv'))

                for index, row in df_study.iterrows():
                    if pd.isna(row['SampleTypeSub1']):
                        org_part = row['Samples_Organism part'].lower()
                        matching_row = df_biofluid_vs_tissue[df_biofluid_vs_tissue['ML'] == org_part]
                        if not matching_row.empty:
                            df_study.at[index, 'SampleTypeSub1'] = matching_row.iloc[0]['tissueVSbiofluid']

                #add UBERON bodypart column
                #######
                df_bodypart = pd.read_csv(os.path.join(transSheet_dir, 'bodypart_Metabolights.csv'))
                df_study = df_study.merge(df_bodypart[['ML', 'bodypart']], left_on='Samples_Organism part', right_on='ML', how='left')
                df_study.rename(columns={'bodypart': 'UBERONBodyPartName'}, inplace=True)
                df_study.drop('ML', axis=1, inplace=True)

            #add MassSpectrometer column
            #######
            if 'Assay_Instrument' in df_study.columns:
                df_study['Assay_Instrument'].fillna('', inplace=True)
                df_ms = pd.read_csv(os.path.join(transSheet_dir, 'MassSpectrometer_Metabolights.csv'))
                df_study = df_study.merge(df_ms[['ML', 'MassSpectrometer']], left_on='Assay_Instrument', right_on='ML', how='left')
                df_study.drop('ML', axis=1, inplace=True)

            #add IonizationMethod/Polarity column
            #######
            if 'Assay_Ion source' in df_study.columns and 'Assay_Scan polarity' in df_study.columns:
                df_study['Assay_Ion source'].fillna('', inplace=True)
                df_study['Assay_Scan polarity'].fillna('', inplace=True)
                df_study['Assay_Ion source'] = df_study['Assay_Ion source'].str.lower()
                df_study['IonizationSourceAndPolarity'] = ''
                df_study.loc[df_study['Assay_Ion source'].str.contains('electrospray'), 'IonizationSourceAndPolarity'] = 'electrospray ionization'
                df_study.loc[df_study['Assay_Ion source'].str.contains('electron ionization'), 'IonizationSourceAndPolarity'] = 'electron ionization'      
                df_study.loc[df_study['Assay_Ion source'].str.contains('chemical ionization'), 'IonizationSourceAndPolarity'] = 'atmospheric pressure chemical ionization'
                df_study.loc[df_study['Assay_Scan polarity'].str.contains('negative'), 'IonizationSourceAndPolarity'] += ' (negative)'
                df_study.loc[df_study['Assay_Scan polarity'].str.contains('positive'), 'IonizationSourceAndPolarity'] += ' (positive)'
                df_study.loc[df_study['Assay_Scan polarity'].str.contains('alternating'), 'IonizationSourceAndPolarity'] += ' (alternating)'


            #add LC method
            #######
            if 'Assay_Column type' in df_study.columns:
                df_study['Assay_Column type'].fillna('', inplace=True)
                df_study['Assay_Column type'] = df_study['Assay_Column type'].str.lower()
                df_study['ChromatographyAndPhase'] = ''
                df_study.loc[df_study['Assay_Column type'].str.contains('reverse'), 'ChromatographyAndPhase'] = 'reverse phase'
                df_study.loc[df_study['Assay_Column type'].str.contains('hilic'), 'ChromatographyAndPhase'] = 'normal phase (HILIC)'
                df_study.loc[df_study['Assay_Column type'].str.contains('normal phase'), 'ChromatographyAndPhase'] = 'normal phase (HILIC)'
                
                if 'Assay_Column model' in df_study.columns:
                    df_study['Assay_Column model'].fillna('', inplace=True)
                    df_study.loc[df_study['Assay_Column model'].str.contains('Phenyl', case=False) & df_study['Assay_Column model'].str.contains('Hexyl', case=False), 'ChromatographyAndPhase'] = ' (Phenyl-Hexyl)'
                    df_study.loc[df_study['Assay_Column model'].str.contains('C18') & df_study['Assay_Column model'].str.contains('polar', case=False), 'ChromatographyAndPhase'] = ' (polar-C18)'
                    
                    df_study.loc[df_study['Assay_Column model'].str.contains('HSS T3') & ~df_study['ChromatographyAndPhase'].str.contains(r"\("), 'ChromatographyAndPhase'] += ' (C18)'
                    df_study.loc[df_study['Assay_Column model'].str.contains('C18') & ~df_study['ChromatographyAndPhase'].str.contains(r"\("), 'ChromatographyAndPhase'] += ' (C18)'
                    df_study.loc[df_study['Assay_Column model'].str.contains('C30') & ~df_study['ChromatographyAndPhase'].str.contains(r"\("), 'ChromatographyAndPhase'] += ' (C30)'
                    df_study.loc[df_study['Assay_Column model'].str.contains('C8') & ~df_study['ChromatographyAndPhase'].str.contains(r"\("), 'ChromatographyAndPhase'] += ' (C8)'
                    df_study.loc[df_study['ChromatographyAndPhase'].str.contains('reverse phase') & ~df_study['ChromatographyAndPhase'].str.contains(r"\("), 'ChromatographyAndPhase'] += ' (NOS)'
                else:
                    df_study.loc[df_study['ChromatographyAndPhase'].str.contains('reverse phase') & ~df_study['ChromatographyAndPhase'].str.contains(r"\("), 'ChromatographyAndPhase'] += ' (NOS)'

            #add AgeInYears 
            #######
            if 'samples_age' in df_study.columns.str.lower():
                if 'Samples_Age' in df_study.columns:
                    df_study.rename(columns={'Samples_Age': 'Samples_age'}, inplace=True)
                if 'Samples_AGE' in df_study.columns:
                    df_study.rename(columns={'Samples_AGE': 'Samples_age'}, inplace=True)

                unit_column = find_column_after_target_column(df_study, 'Samples_age', 'Samples_Unit')

                if unit_column != '':
                    df_study['Samples_age'] = pd.to_numeric(df_study['Samples_age'], errors='coerce')
                    df_study['AgeInYears'] = np.nan
                    df_study.loc[df_study[unit_column] == 'year', 'AgeInYears'] = df_study.loc[df_study[unit_column] == 'year', 'Samples_age']
                    df_study.loc[df_study[unit_column] == 'month', 'AgeInYears'] = df_study.loc[df_study[unit_column] == 'month', 'Samples_age'] / 12
                    df_study.loc[df_study[unit_column] == 'week', 'AgeInYears'] = df_study.loc[df_study[unit_column] == 'week', 'Samples_age'] / 52.1429
                    df_study.loc[df_study[unit_column] == 'day', 'AgeInYears'] = df_study.loc[df_study[unit_column] == 'day', 'Samples_age'] / 365
                    df_study.loc[df_study[unit_column] == 'hour', 'AgeInYears'] = df_study.loc[df_study[unit_column] == 'hour', 'Samples_age'] / 8760

                    df_study['AgeInYears'] = df_study['AgeInYears'].astype(str).replace('nan', 'ML import: not available')

            #add Sex
            #######
            if 'samples_sex' in df_study.columns.str.lower():
                if 'Samples_Sex' in df_study.columns:
                    df_study.rename(columns={'Samples_Sex': 'Samples_sex'}, inplace=True)
                if 'Samples_SEX' in df_study.columns:
                    df_study.rename(columns={'Samples_SEX': 'Samples_sex'}, inplace=True)
            
            if 'samples_gender' in df_study.columns.str.lower():
                if 'Samples_Gender' in df_study.columns:
                    df_study.rename(columns={'Samples_Gender': 'Samples_gender'}, inplace=True)
                if 'Samples_GENDER' in df_study.columns:
                    df_study.rename(columns={'Samples_GENDER': 'Samples_gender'}, inplace=True)
            
            if 'Samples_gender' in df_study.columns:

                df_study["BiologicalSex"] = 'ML import: not available'
                df_study.loc[df_study['Samples_gender'].str.lower().str.contains('female'), 'BiologicalSex'] = 'female'
                df_study.loc[(df_study['Samples_gender'].str.lower().str.contains('male')) & (df_study['BiologicalSex'] != 'female'), 'BiologicalSex'] = 'male'
            

            #add MassiveID and USIs
            #######
            df_study['MassiveID'] = study_id
            df_study['USI'] = 'mzspec:' + study_id + ':' + df_study['filepath']

            df_study = complete_and_fill_REDU_table(df_study, allowedTerm_dict)

            return df_study




if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Give an Metabolights study ID and get a REDU table tsv.')
    parser.add_argument("--study_id", type=str, help='An Metabolights study ID such as "MTBLS1015". If "ALL" all studys are requested.', required=True)
    parser.add_argument("--path_to_translation_sheet_csvs", "-csvs", type=str, help="Path to the translation csvs holding translations from MWB to REDU vocabulary (optional)", default="none")
    parser.add_argument("--path_to_allowed_term_json", type=str, help="Path to the json with allowed REDU terms (optional)", default="none")
    parser.add_argument("--path_to_unassigned_term_json", type=str, help="If you want to export a json with terms that have not been associated with anything (optional)", default="none")
    
    
    args = parser.parse_args()

    if args.study_id == 'ALL':
        public_metabolights_studies = safe_api_request('https://www.ebi.ac.uk:443/metabolights/ws/studies', retries = 1)
        public_metabolights_studies = public_metabolights_studies['content']
    else:
        public_metabolights_studies = [args.study_id]

    if args.path_to_translation_sheet_csvs == 'none':
        script_dir = os.path.dirname(os.path.realpath(__file__))
        transSheet_dir = os.path.join(script_dir, 'translation_sheets_metabolights')
    else:
        transSheet_dir = args.path_to_translation_sheet_csvs
        script_dir = os.path.dirname(transSheet_dir)

    if args.path_to_allowed_term_json == 'none':
        allowedTermSheet_json = os.path.join(script_dir, 'allowed_terms', 'allowed_terms.json')
    else:
        print('set')
        allowedTermSheet_json = args.path_to_allowed_term_json
    print('here')
    print(allowedTermSheet_json)
    print('there')
    # Read allowed terms json
    with open(allowedTermSheet_json, 'r') as file:
        allowedTerm_dict = json.load(file)


    REDU_dataframes = []
    redu_table_single = pd.DataFrame()
    for study_id in tqdm(public_metabolights_studies):
        try:
            redu_table_single = Metabolights2REDU(study_id, allowedTerm_dict = allowedTerm_dict, unassigned_term_json = args.path_to_unassigned_term_json)
        except Exception as e:
            traceback_info = traceback.format_exc()
            print(f"An error occurred with study_id {study_id}: {e}\nTraceback:\n{traceback_info}")
            continue
        if redu_table_single is not None and len(redu_table_single) > 0:
            REDU_dataframes.append(redu_table_single)

    if len(REDU_dataframes) > 0:
        redu_tables_all = pd.concat(REDU_dataframes, ignore_index=True)
        redu_tables_all.to_csv('Metabolights2REDU_' + args.study_id + '.tsv', sep='\t', index=False, header=True)
        print(f'Output has been saved to Metabolights2REDU_{args.study_id}.tsv!')
    else:
        print('nothing to return!')




