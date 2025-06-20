import os
import requests
import pandas as pd
from bs4 import BeautifulSoup
import argparse
import json
import time
import traceback
import numpy as np
from tqdm import tqdm
from REDU_conversion_functions import age_category
from REDU_conversion_functions import get_taxonomy_id_from_name__allowedTerms
from REDU_conversion_functions import get_taxonomy_info
from REDU_conversion_functions import merge_repeated_fileobservations
from read_and_validate_redu_from_github import complete_and_fill_REDU_table
from REDU_conversion_functions import find_column_after_target_column




def prefer_extension(group):
    extensions = group['extension'].values
    if 'mzml' in extensions or 'mzxml' in extensions:
        group['keep'] = (group['extension'] == 'mzml') | (group['extension'] == 'mzxml')
    else:
        # Mark only the first row to keep if no preferred extension found
        group['keep'] = [True] + [False] * (len(group) - 1)
    return group


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
        if assay_json['technology'] == 'mass spectrometry' or assay_json['technology'] == 'mass spectrometry assay':        

            # Extract headers in the correct order
            headers = [None] * len(assay_json['assayTable']['fields'])
            for key, value in assay_json['assayTable']['fields'].items():
                headers[value['index']] = 'Assay_' + value['header']

            df = pd.DataFrame(assay_json['assayTable']['data'])
            df = rename_duplicated_column_names(df)

            df.columns = headers

            print(f'Assay {index} has {len(df)} rows and {len(df.columns)} columns.')
            
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

        print(f'Assay table has {len(df_assays)} rows and {len(df_assays.columns)} columns.')

        #get info on samples (there is only one sample table per study)
        headers = [None] * len(study_details['content']['sampleTable']['fields'])
        for key, value in study_details['content']['sampleTable']['fields'].items():
            headers[value['index']] = 'Samples_' + value['header']

        df_samples = pd.DataFrame(study_details['content']['sampleTable']['data'])
        df_samples = rename_duplicated_column_names(df_samples)
        df_samples.columns = headers
    
        df_study = df_assays.merge(df_samples, left_on='Assay_Sample Name', right_on='Samples_Sample Name', how='inner')
        df_study['row_id'] = range(1, len(df_study) + 1)

        print(f'Sample table has {len(df_study)} rows and {len(df_study.columns)} columns.')

        # Duplicate rows if we have mzml AND raw files
        df_study_raw = pd.DataFrame()
        
        if 'Assay_Raw Spectral Data File' in df_study.columns:
            df_study_raw = df_study[df_study['Assay_Raw Spectral Data File'].str.contains(r'\.', regex=True, na=False)].copy()

            df_study_raw['filepath'] = df_study_raw['Assay_Raw Spectral Data File']
            # Check if the column exists before dropping
            if 'Assay_Raw Spectral Data File' in df_study_raw.columns:
                df_study_raw.drop(columns=['Assay_Raw Spectral Data File'], inplace=True)
            if 'Assay_Derived Spectral Data File' in df_study_raw.columns:
                df_study_raw.drop(columns=['Assay_Derived Spectral Data File'], inplace=True)
            print(f'Assay table has {len(df_study_raw)} rows and {len(df_study_raw.columns)} columns.')

        df_study_mzml = pd.DataFrame()
        if 'Assay_Derived Spectral Data File' in df_study.columns:
            df_study_mzml = df_study[df_study['Assay_Derived Spectral Data File'].str.contains(r'\.', regex=True, na=False)].copy()
            df_study_mzml['filepath'] = df_study_mzml['Assay_Derived Spectral Data File']
            # Check if the column exists before dropping
            if 'Assay_Raw Spectral Data File' in df_study_mzml.columns:
                df_study_mzml.drop(columns=['Assay_Raw Spectral Data File'], inplace=True)
            if 'Assay_Derived Spectral Data File' in df_study_mzml.columns:
                df_study_mzml.drop(columns=['Assay_Derived Spectral Data File'], inplace=True)
            print(f'Assay table has {len(df_study_mzml)} rows and {len(df_study_mzml.columns)} columns.')

        if len(df_study_raw) > 0 and len(df_study_mzml) > 0:
            df_study = pd.concat([df_study_mzml, df_study_raw], ignore_index=True).copy(deep=True)
        elif len(df_study_raw) > 0:
            df_study = df_study_raw.copy(deep=True)
        elif len(df_study_mzml) > 0:
            df_study = df_study_mzml.copy(deep=True)
        elif len(df_study_raw) == 0 and len(df_study_mzml) == 0:
            return pd.DataFrame()
        
        df_study = rename_duplicated_column_names(df_study)
        df_study['filename'] = df_study['filepath']#.apply(lambda x: os.path.basename(x))

        # Ensure the filename is a string and lowercase
        df_study['filename_lower'] = df_study['filename'].astype(str).str.lower()

        # Perform the split operation with expansion to get two columns
        split_values = df_study['filename_lower'].str.rsplit('.', n=1, expand=True)

        # List of allowed extensions
        allowed_extensions = [".mzml", ".mzxml", ".cdf", ".raw", ".wiff", ".d"]
        df_study = df_study[df_study['filename_lower'].apply(lambda x: any(x.endswith(ext) for ext in allowed_extensions))]

        print(f'Assay table has {len(df_study)} rows and {len(df_study.columns)} columns.')

        if len(df_study) == 0:
            return None

        # Assign the split columns to 'base' and 'extension'
        df_study['base'] = split_values[0]

        # Initially set 'extension' to None to handle rows without an extension
        df_study['extension'] = None

        # Only update 'extension' for rows where a split occurred
        df_study.loc[split_values[1].notna(), 'extension'] = split_values[1]

        df_study = df_study.groupby('row_id').apply(prefer_extension)
        df_study = df_study[df_study['keep']]

        if len(df_study) > 0:
            
            allowedTerm_dict = kwargs['allowedTerm_dict']
            ontology_table = kwargs['ontology_table']
            ENVOEnvironmentMaterialIndex_table = kwargs['ENVOEnvironmentMaterialIndex_table']
            ENVOEnvironmentBiomeIndex_table = kwargs['ENVOEnvironmentBiomeIndex_table']
            NCBIRankDivision_table  = kwargs['NCBIRankDivision_table']

            
            df_study.loc[:, 'YearOfAnalysis'] = submissionYear

            #add NCBITaxonomy and Sampletype & SampleTypeSub1
            #######
            if 'Samples_Organism' in df_study.columns:
                processed_organisms = {org: str(get_taxonomy_id_from_name__allowedTerms(org, allowedTerm_dict = allowedTerm_dict)) for org in df_study['Samples_Organism'].unique()}

                df_study.loc[:, 'NCBITaxonomy'] = df_study['Samples_Organism'].map(processed_organisms)
                df_study.loc[:, 'NCBITaxonomy'] = df_study['NCBITaxonomy'].replace(to_replace=r'^.*None.*$', value='missing value', regex=True)
                
                df_study.loc[:, ['SampleType', 'SampleTypeSub1']] = 'missing value'

                processed_taxonomy = {taxonomy.split('|')[0]: get_taxonomy_info(taxonomy.split('|')[0])
                                    for taxonomy in df_study['NCBITaxonomy'].unique()
                                    if '|' in taxonomy and 'None' not in taxonomy}
                
                print(f'Found {len(processed_taxonomy)} unique taxonomy IDs.')


                df_study.loc[:,'SampleType'] = df_study['NCBITaxonomy'].map(lambda x: processed_taxonomy.get(x.split('|')[0], [pd.NA, pd.NA])[0] 
                                                                    if '|' in x and 'None' not in x else pd.NA)
                df_study.loc[:,'SampleTypeSub1'] = df_study['NCBITaxonomy'].map(lambda x: processed_taxonomy.get(x.split('|')[0], [pd.NA, pd.NA])[1]
                                                                        if '|' in x and 'None' not in x else pd.NA)

                df_study[['SampleType', 'SampleTypeSub1']] = df_study.apply(
                    lambda row: get_blanks(row.Samples_Organism) 
                    if (pd.isna(row.SampleType) or row.SampleType == 'missing value') 
                    else [row.SampleType, row.SampleTypeSub1], 
                    axis=1
                ).apply(pd.Series)

                df_study[['SampleType', 'SampleTypeSub1']] = df_study.apply(
                    lambda row: get_enviromental_water(row.Samples_Organism) 
                    if (pd.isna(row.SampleType) or row.SampleType == 'missing value') 
                    else [row.SampleType, row.SampleTypeSub1], 
                    axis=1
                ).apply(pd.Series)


                # df_biofluid_vs_tissue = pd.read_csv(os.path.join(transSheet_dir, 'biofluid_vs_tissue_Metabolights.csv'))

                # df_study['Samples_Organism part'] = df_study['Samples_Organism part'].str.lower()

                # for index, row in df_study.iterrows():
                #     if pd.isna(row['SampleTypeSub1']):
                #         org_part = row['Samples_Organism part']
                #         matching_row = df_biofluid_vs_tissue[df_biofluid_vs_tissue['ML'] == org_part]
                #         if not matching_row.empty:
                #             df_study.at[index, 'SampleTypeSub1'] = matching_row.iloc[0]['tissueVSbiofluid']


                #add ENV material column
                #######
                updated_species_dict = {}

                for key, value in processed_organisms.items():
                    if value == 'None':
                        match = ENVOEnvironmentMaterialIndex_table[(ENVOEnvironmentMaterialIndex_table['Label'] == key) | (ENVOEnvironmentMaterialIndex_table['Synonym'] == key)]['Label'].values
                        if match.size > 0:
                            updated_species_dict[key] = match[0]
                        else:
                            updated_species_dict[key] = 'missing value'

                df_study.loc[:, 'ENVOEnvironmentMaterial'] = df_study['Samples_Organism'].map(updated_species_dict)


                #add ENV biome column
                #######
                updated_species_dict = {}

                for key, value in processed_organisms.items():
                    if value == 'None':
                        match = ENVOEnvironmentBiomeIndex_table[(ENVOEnvironmentBiomeIndex_table['Label'] == key) | (ENVOEnvironmentBiomeIndex_table['Synonym'] == key)]['Label'].values
                        if match.size > 0:
                            updated_species_dict[key] = match[0]
                        else:
                            updated_species_dict[key] = 'missing value'

                df_study.loc[:, 'ENVOEnvironmentBiome'] = df_study['Samples_Organism'].map(updated_species_dict)


                #add UBERON bodypart column
                #######
                df_study['UBERONBodyPartName'] = 'missing value'

                labels_in_ontology_table = set(ontology_table['Label'])
                synonyms_in_ontology_table = set(ontology_table['Synonym'])
                allowed_bodyparts = set(allowedTerm_dict['UBERONBodyPartName']['allowed_values'])
                allowed_bodypart_ids = set(allowedTerm_dict['UBERONOntologyIndex']['allowed_values'])

                for index, row in df_study.iterrows():
                    organism_part = row['Samples_Organism part']

                    if organism_part in allowed_bodyparts:
                        df_study.at[index, 'UBERONBodyPartName'] = organism_part
                        continue 

                    if organism_part in synonyms_in_ontology_table:

                        matching_rows = ontology_table[ontology_table['Synonym'] == organism_part]

                        plant = True if row['SampleType'] == 'plant' else False

                        if len(matching_rows) > 1 and plant == True:
                            matching_rows = matching_rows[matching_rows['UBERONOntologyIndex'].str.startswith('PO', na=False)]

                        if len(matching_rows) > 1 and plant == False:
                            matching_rows = matching_rows[(~matching_rows['UBERONOntologyIndex'].str.startswith('PO', na=False))]

                        if len(matching_rows) == 1:
                            
                            label = matching_rows['Label'].iloc[0] 
                            df_study.at[index, 'UBERONBodyPartName'] = label

                        
                ontology_table_unique = ontology_table.drop_duplicates(subset=['Label', 'UBERONOntologyIndex']).drop(columns=['Synonym'])

                df_study = pd.merge(df_study, ontology_table_unique, left_on='UBERONBodyPartName', right_on='Label', how='left')
                df_study.drop('Label', axis=1, inplace=True)

                df_study.loc[df_study['Is Fluid'] == True, 'SampleTypeSub1'] = 'biofluid'
                df_study.loc[df_study['Is Multicellular'] == True, 'SampleTypeSub1'] = 'tissue'
                df_study['UBERONOntologyIndex'] = df_study['UBERONOntologyIndex'].str.replace("_", ":", regex=False)


            #add MassSpectrometer column
            #######
            if 'Assay_Instrument' in df_study.columns:
                df_study['Assay_Instrument'] = df_study['Assay_Instrument'].fillna('')

                prefixes = [
                    "Thermo Scientific ",
                    "Waters ",
                    "Agilent ",
                    "AB SCIEX ",
                    "LECO ",
                    "Bruker ",
                    "Shimadzu "
                ]

                # Function to remove prefixes
                def remove_prefixes(value, prefixes):
                    for prefix in prefixes:
                        if value.startswith(prefix):
                            return value[len(prefix):]
                    return value

                # Apply the function to the column
                df_study['Assay_Instrument'] = df_study['Assay_Instrument'].apply(lambda x: remove_prefixes(x, prefixes))

                allowed_values_map = dict(zip(
                    allowedTerm_dict["MassSpectrometer"]["allowed_values_matching_0"],
                    allowedTerm_dict["MassSpectrometer"]["allowed_values"]
                ))

                df_study["MassSpectrometer"] = df_study["Assay_Instrument"].map(allowed_values_map)



            #add IonizationMethod/Polarity column
            #######
            if 'Assay_Ion source' in df_study.columns and 'Assay_Scan polarity' in df_study.columns:
                df_study['Assay_Ion source'] = df_study['Assay_Ion source'].fillna('')
                df_study['Assay_Scan polarity'] = df_study['Assay_Scan polarity'].fillna('')

                df_study['Assay_Ion source'] = df_study['Assay_Ion source'].str.lower()
                df_study['IonizationSourceAndPolarity'] = ''
                df_study.loc[df_study['Assay_Ion source'].str.contains('electrospray'), 'IonizationSourceAndPolarity'] = 'electrospray ionization'
                df_study.loc[df_study['Assay_Ion source'].str.contains('electron ionization'), 'IonizationSourceAndPolarity'] = 'electron ionization'      
                df_study.loc[df_study['Assay_Ion source'].str.contains('chemical ionization'), 'IonizationSourceAndPolarity'] = 'atmospheric pressure chemical ionization'
                df_study.loc[df_study['Assay_Scan polarity'].str.contains('negative') & ~df_study['IonizationSourceAndPolarity'].str.contains("electron ionization"), 'IonizationSourceAndPolarity'] += ' (negative)'
                df_study.loc[df_study['Assay_Scan polarity'].str.contains('positive') & ~df_study['IonizationSourceAndPolarity'].str.contains("electron ionization"), 'IonizationSourceAndPolarity'] += ' (positive)'
                df_study.loc[df_study['Assay_Scan polarity'].str.contains('alternating') & ~df_study['IonizationSourceAndPolarity'].str.contains("electron ionization"), 'IonizationSourceAndPolarity'] += ' (alternating)'


            #add LC method
            #######
            if 'Assay_Column type' in df_study.columns:
                df_study['Assay_Column type'] = df_study['Assay_Column type'].fillna('')
                df_study['Assay_Column type'] = df_study['Assay_Column type'].str.lower()
                df_study['ChromatographyAndPhase'] = ''
                df_study.loc[df_study['Assay_Column type'].str.contains('reverse'), 'ChromatographyAndPhase'] = 'reverse phase'
                df_study.loc[df_study['Assay_Column type'].str.contains('hilic'), 'ChromatographyAndPhase'] = 'normal phase (HILIC)'
                df_study.loc[df_study['Assay_Column type'].str.contains('normal phase'), 'ChromatographyAndPhase'] = 'normal phase (HILIC)'
                
                if 'Assay_Column model' in df_study.columns:
                    df_study['Assay_Column model'] = df_study['Assay_Column model'].fillna('')
                    df_study.loc[df_study['Assay_Column model'].str.contains('Phenyl', case=False) & df_study['Assay_Column model'].str.contains('Hexyl', case=False), 'ChromatographyAndPhase'] += ' (Phenyl-Hexyl)'
                    df_study.loc[df_study['Assay_Column model'].str.contains('C18') & df_study['Assay_Column model'].str.contains('polar', case=False), 'ChromatographyAndPhase'] += ' (polar-C18)'
                    
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

                    df_study['AgeInYears'] = df_study['AgeInYears'].astype(str).replace('nan', 'missing value')

            #add AgeInYears 
            #######
                    df_study['LifeStage'] = df_study['AgeInYears'].apply(lambda x: age_category(x))

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

                df_study["BiologicalSex"] = 'missing value'
                df_study.loc[df_study['Samples_gender'].str.lower().str.contains('female'), 'BiologicalSex'] = 'female'
                df_study.loc[(df_study['Samples_gender'].str.lower().str.contains('male')) & (df_study['BiologicalSex'] != 'female'), 'BiologicalSex'] = 'male'


            #add MassiveID and USIs
            #######
            df_study['MassiveID'] = study_id
            
            ontology_table = ontology_table.drop_duplicates(subset=['Label'])
            df_study = merge_repeated_fileobservations(df_study)
            df_study = complete_and_fill_REDU_table(df_study, allowedTerm_dict, UBERONOntologyIndex_table=ontology_table, ENVOEnvironmentBiomeIndex_table=ENVOEnvironmentBiomeIndex_table,
                                                    ENVOEnvironmentMaterialIndex_table=ENVOEnvironmentMaterialIndex_table,NCBIRankDivision_table=NCBIRankDivision_table, add_usi = True, 
                                                    other_allowed_file_extensions = ['.raw', '.cdf', '.wiff', '.d'])
            
            df_study = df_study.drop_duplicates() 

            #remove files if they are assigned multiple times as we cannot tell which sample they belong to (this is probably because people make mistakes when creating their study)
            df_study['count'] = df_study.groupby('USI')['USI'].transform('size')
            df_study = df_study[df_study['count'] == 1]
            df_study = df_study.drop(columns=['count'])

            print(f'We are adding {len(df_study)} samples to the REDU table for study {study_id}.')

            return df_study




if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Give an Metabolights study ID and get a REDU table tsv.')
    parser.add_argument("--study_id", type=str, help='An Metabolights study ID such as "MTBLS1015". If "ALL" all studys are requested.', required=True)
    parser.add_argument("--path_to_translation_sheet_csvs", "-csvs", type=str, help="Path to the translation csvs holding translations from MWB to REDU vocabulary", default="none")
    parser.add_argument("--path_to_allowed_term_json", type=str, help="Path to the json with allowed REDU terms")
    parser.add_argument("--path_to_uberon_cl_po_csv", type=str, help="Path to the prepared uberon_cl_po ontology csv")
    parser.add_argument("--path_to_envo_biome_csv", type=str, help="Path to the prepared uberon_cl_po ontology csv")
    parser.add_argument("--path_to_envo_material_csv", type=str, help="Path to the prepared uberon_cl_po ontology csv")
    parser.add_argument("--path_ncbi_rank_division", type=str, help="Path to the path_ncbi_rank_division")
            
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
        allowedTermSheet_json = args.path_to_allowed_term_json


    # Read allowed terms json
    with open(allowedTermSheet_json, 'r') as file:
        allowedTerm_dict = json.load(file)

    # Add terms for matching
    allowed_values = allowedTerm_dict["MassSpectrometer"]["allowed_values"]
    allowedTerm_dict["MassSpectrometer"]["allowed_values_matching_0"] = [value.split('|')[0] for value in allowed_values]


    # Read ontology tables
    ontology_table = pd.read_csv(args.path_to_uberon_cl_po_csv)
    ENVOEnvironmentBiomeIndex_table = pd.read_csv(args.path_to_envo_biome_csv)
    ENVOEnvironmentMaterialIndex_table = pd.read_csv(args.path_to_envo_material_csv)

    # Read NCBI rank and division csv
    NCBIRankDivision_table = pd.read_csv(args.path_ncbi_rank_division, index_col = False)
    NCBIRankDivision_table = NCBIRankDivision_table.drop_duplicates(subset=['TaxonID'])


    REDU_dataframes = []
    redu_table_single = pd.DataFrame()
    for study_id in tqdm(public_metabolights_studies):
        try:
            print(f'Processing study {study_id}...')
            redu_table_single = Metabolights2REDU(study_id, allowedTerm_dict = allowedTerm_dict, ontology_table = ontology_table, ENVOEnvironmentBiomeIndex_table=ENVOEnvironmentBiomeIndex_table,
                                                  ENVOEnvironmentMaterialIndex_table=ENVOEnvironmentMaterialIndex_table, NCBIRankDivision_table=NCBIRankDivision_table)
        except Exception as e:
            traceback_info = traceback.format_exc()
            print(f"An error occurred with study_id {study_id}: {e}\nTraceback:\n{traceback_info}")
            continue
        if redu_table_single is not None and len(redu_table_single) > 0:
            print(f'Added {len(redu_table_single)} samples.')
            REDU_dataframes.append(redu_table_single)
        else:
            print(f'Added {0} samples.')

    if len(REDU_dataframes) > 0:
        redu_tables_all = pd.concat(REDU_dataframes, ignore_index=True)
        redu_tables_all = redu_tables_all.drop_duplicates()
        redu_tables_all.to_csv('Metabolights2REDU_' + args.study_id + '.tsv', sep='\t', index=False, header=True)
        print(f'Output of {len(redu_tables_all)} samples has been saved to Metabolights2REDU_{args.study_id}.tsv!')
    else:
        print('nothing to return!')

