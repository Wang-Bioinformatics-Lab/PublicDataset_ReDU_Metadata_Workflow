import argparse
import pandas as pd
import glob
import re
import os
import json
from io import StringIO
from REDU_conversion_functions import age_category
from REDU_conversion_functions import get_uberon_table
from REDU_conversion_functions import get_ontology_table

def complete_and_fill_REDU_table(df, allowedTerm_dict, **kwargs):
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

    if 'UBERONOntologyIndex_table' in kwargs.keys():
        uberon_ontology_table = kwargs['UBERONOntologyIndex_table']
    if 'DOIDOntologyIndex_table' in kwargs.keys():
        doid_ontology_table = kwargs['DOIDOntologyIndex_table']
        if 'UBERONOntologyIndex' in doid_ontology_table.columns:
            doid_ontology_table.rename(columns={'UBERONOntologyIndex': 'DOIDOntologyIndex'}, inplace=True)

    # Convert year to string for comparison
    if 'YearOfAnalysis' in df.columns:
        df['YearOfAnalysis'] = df['YearOfAnalysis'].astype(str)
    if 'AgeInYears' in df.columns:
        df['AgeInYears'] = df['AgeInYears'].astype(str)



    #remove autogenerated columns if present
    columns_to_remove = ['LifeStage', 'UniqueSubjectID', 'UBERONOntologyIndex', 'DOIDOntologyIndex', 'USI']
    df.drop(columns=[col for col in columns_to_remove if col in df.columns], inplace=True)

    # Add missing columns with their respective default missing value from the dictionary or generate values if generate == True
    for key, value in allowedTerm_dict.items():
        if key not in df.columns:
            if value['generate'] == 'False':
                if value['missing'] == 'not allowed':
                    return pd.DataFrame()
                else:
                    df[key] = value['missing']
            elif value['generate'] == 'True':
                if key == 'LifeStage':
                    df[key] = df.apply(lambda x: age_category(x['AgeInYears']) if x['NCBITaxonomy'] == "9606|Homo sapiens" else value['missing'], axis=1)
                if key == 'UniqueSubjectID':
                    df[key] = df.apply(lambda x: str(x['MassiveID']) + '_' + str(x['SubjectIdentifierAsRecorded']) if x['SubjectIdentifierAsRecorded'] != allowedTerm_dict['SubjectIdentifierAsRecorded']['missing'] else value['missing'], axis=1)
                if key == 'UBERONOntologyIndex':
                    df = df.merge(uberon_ontology_table[['Label', 'UBERONOntologyIndex']], left_on='UBERONBodyPartName', right_on='Label', how='left')
                    df.drop(columns=['Label'], inplace=True)
                if key == 'DOIDOntologyIndex':
                    df = df.merge(doid_ontology_table[['Label', 'DOIDOntologyIndex']], left_on='DOIDCommonName', right_on='Label', how='left')
                    df.drop(columns=['Label'], inplace=True)
                if key == 'USI':
                    pass
                    # df['USI'] = 'mzspec:' + df['MassiveID'] + ':' + df['filename']


    # Replace values with the respective "missing" value if they're not in the allowed terms or are missing/empty
    for key, value in allowedTerm_dict.items():
        allowed_terms = value['allowed_values']
        missing_value = value['missing']
        if key in df.columns:
            if value['generate'] == 'False':
                if key != 'filename' and len(allowed_terms) > 1:
                    df[key] = df[key].apply(lambda x: x if x in allowed_terms else missing_value).fillna(missing_value).replace("", missing_value)
                elif allowed_terms[0] == '00':
                    df[key] = df[key].fillna(missing_value).replace("", missing_value)
                elif allowed_terms[0] == 'numeric':
                    df[key] = pd.to_numeric(df[key], errors='coerce').fillna(missing_value).replace("", missing_value)
                elif allowed_terms[0] == 'numeric|numeric':
                    df[key] = df[key].apply(
                        lambda x: x if all(part.replace('.', '', 1).isdigit() or part.lstrip('-').replace('.', '', 1).isdigit() for part in x.split('|')) else missing_value
                    ).replace("", missing_value)
                elif key == 'filename':
                    df['filename'] = df['filename'].apply(lambda x: x if any(x.endswith(ext) for ext in allowed_terms) else missing_value)


    # Ensure the dataframe contains only the columns specified in the dictionary
    keys_to_include = [key for key in allowedTerm_dict.keys() if key != 'USI']
    return df[keys_to_include]


def clean_read_tsv_quoted_lines(path):
    cleaned_data = ''
    with open(path, 'r') as file:
        for line in file:
            # Strip leading and trailing whitespace and then strip double quotes
            cleaned_line = line.strip().strip('"')
            cleaned_data += cleaned_line + '\n'

    # Use StringIO to simulate a file object for pandas
    cleaned_data_io = StringIO(cleaned_data)

    # Now read the cleaned data into a DataFrame
    df = pd.read_csv(cleaned_data_io, sep='\t')

    return df

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='GNPS Validator')
    parser.add_argument('path_to_github_metadata')    
    parser.add_argument('output_metadata_folder')
    parser.add_argument("--AllowedTermJson_path", type=str, help="Path to json with allowed terms")
    parser.add_argument('--path_to_uberon_cl_po_csv')
    parser.add_argument('--path_to_doid_csv')
    # parser.add_argument('--path_to_uberon_owl', default = 'none')
    # parser.add_argument('--path_to_po_owl', default = 'none')
    # parser.add_argument('--path_to_cl_owl', default = 'none')
    # parser.add_argument('--path_to_doid_owl', default = 'none')
    args = parser.parse_args()

    #python3.8 ../ReDU-MS2-GNPS2/workflows/PublicDataset_ReDU_Metadata_Workflow/bin/read_and_validate_redu_from_github.py /home/yasin/projects/ReDU_metadata/metadata output/ --AllowedTermJson_path /home/yasin/projects/ReDU-MS2-GNPS2/workflows/PublicDataset_ReDU_Metadata_Workflow/bin/allowed_terms/allowed_terms.json --path_to_uberon_owl /home/yasin/projects/ReDU-MS2-GNPS2/workflows/PublicDataset_ReDU_Metadata_Workflow/bin/allowed_terms/uberon.owl --path_to_po_owl /home/yasin/projects/ReDU-MS2-GNPS2/workflows/PublicDataset_ReDU_Metadata_Workflow/bin/allowed_terms/po.owl --path_to_cl_owl /home/yasin/projects/ReDU-MS2-GNPS2/workflows/PublicDataset_ReDU_Metadata_Workflow/bin/allowed_terms/cl.owl  --path_to_doid_owl /home/yasin/projects/ReDU-MS2-GNPS2/workflows/PublicDataset_ReDU_Metadata_Workflow/bin/allowed_terms/doid.owl

    # Read the JSON file
    with open(args.AllowedTermJson_path, 'r') as json_file:
        allowed_terms = json.load(json_file)

    #create ontology table to fill UBERONOntologyIndex from bodyparts
    # uberon_onto = get_uberon_table(args.path_to_uberon_owl)
    # cl_onto = get_ontology_table(args.path_to_cl_owl, ont_prefix = 'CL_')
    # po_onto = get_ontology_table(args.path_to_po_owl, ont_prefix = 'PO_', rm_synonym_info = True)

    # uberon_ontology_table = pd.concat([uberon_onto, cl_onto, po_onto], ignore_index=True, sort=False)
    # uberon_ontology_table['UBERONOntologyIndex'] = uberon_ontology_table['UBERONOntologyIndex'].str.replace('_', ':')

    # #create ontology table to fill DOIDOntologyIndex from bodyparts
    # doid_ontology_table = get_ontology_table(args.path_to_doid_owl, ont_prefix = 'DOID_', index_column_name = 'DOIDOntologyIndex')

    # doid_ontology_table['DOIDOntologyIndex'] = doid_ontology_table['DOIDOntologyIndex'].str.replace('_', ':')


    uberon_ontology_table = pd.read_csv(args.path_to_uberon_cl_po_csv, index_col=False)
    doid_ontology_table = pd.read_csv(args.path_to_doid_csv, index_col=False)


    # Loop through all .tsv files in the given directory
    #for file_path in glob.glob(f"{args.path_to_github_metadata}/redu_*.tsv"):
    for file_path in ["/home/yasin/projects/ReDU_metadata/metadata/redu_MSV000093329.tsv"]:
        
        if re.match(r'.*redu_MSV\d+\.tsv$', file_path):

            df = pd.read_csv(file_path, sep='\t')

            if 'ATTRIBUTE_MassiveID' in doid_ontology_table.columns:
                doid_ontology_table.rename(columns={'ATTRIBUTE_MassiveID': 'MassiveID'}, inplace=True)

            #generate extra columns, add missing columns and remove values which are not in allowed terms
            df = complete_and_fill_REDU_table(df, allowed_terms, UBERONOntologyIndex_table=uberon_ontology_table, DOIDOntologyIndex_table=doid_ontology_table)
            if len(df) > 0:
            
                file_name = os.path.join(args.output_metadata_folder, os.path.basename(file_path))
                df.to_csv(file_name, sep='\t', index=False)
