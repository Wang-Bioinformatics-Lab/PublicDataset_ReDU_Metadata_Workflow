import argparse
import pandas as pd
import json
import os
from read_and_validate_redu_from_github import complete_and_fill_REDU_table






if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='GNPS Validator')
    parser.add_argument('path_REDUtable')
    parser.add_argument('output_folder')
    parser.add_argument('--path_microbeMASST')
    parser.add_argument('--path_plantMASST')
    parser.add_argument("--AllowedTermJson_path", type=str, help="Path to json with allowed terms")
    args = parser.parse_args()



    df_redu = pd.read_csv(args.path_REDUtable, sep = '\t')
    df_redu['filename'] = df_redu['filename'].str.slice(2)

    df_mm = pd.read_csv(args.path_microbeMASST, dtype=str)
    print(df_mm)
    df_mm['filename'] = df_mm['Filepath'].apply(lambda x: '/'.join(x.split('/')[1:]))
    df_mm = df_mm[~df_mm['filename'].isin(df_redu['filename'])]

    df_pm = pd.read_csv(args.path_plantMASST, dtype=str)
    print(df_pm)
    df_pm['filename'] = df_pm['Filepath'].apply(lambda x: '/'.join(x.split('/')[1:]))
    df_pm = df_pm[~df_pm['filename'].isin(df_redu['filename'])]


    with open(args.AllowedTermJson_path, 'r') as json_file:
        allowed_terms = json.load(json_file)

    

    df_list = []

    if len(df_mm) > 0:
        df_mm['NCBITaxonomy'] = df_mm['Taxa_NCBI'] + '|' + df_mm['Taxaname_file']
        df_mm = df_mm.rename(columns={'MassIVE': 'MassiveID'})
        df_mm = df_mm[['MassiveID', 'filename', 'NCBITaxonomy']]
        df_list.append(df_mm)


    if len(df_pm) > 0:
        df_pm['NCBITaxonomy'] = df_pm['Taxa_NCBI'] + '|' + df_pm['Taxaname_file']
        df_pm = df_pm.rename(columns={'MassIVE': 'MassiveID'})
        df_pm['SampleType'] = 'plant'
        df_pm = df_pm[['MassiveID', 'filename', 'SampleType', 'NCBITaxonomy']]
        df_list.append(df_pm)


    if len(df_list) > 0:

        df_massts = pd.concat(df_list, ignore_index=True)
        df_massts_filled = complete_and_fill_REDU_table(df_massts, allowedTerm_dict=allowed_terms)

        #save output to csv
        for massive_id in df_massts_filled['MassiveID'].unique():
            # Filter the DataFrame for the current MassiveID
            df_filtered = df_massts_filled[df_massts_filled['MassiveID'] == massive_id]
            
            # Construct the file name
            file_name = f"{args.output_folder}/{massive_id}.tsv"
            
            # Save the filtered DataFrame to a TSV file
            df_filtered.to_csv(file_name, sep='\t', index=False)


