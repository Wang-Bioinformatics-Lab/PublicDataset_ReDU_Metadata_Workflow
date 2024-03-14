import pandas as pd
import argparse
import json
from REDU_conversion_functions import get_uberon_table
from REDU_conversion_functions import get_ontology_table

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='update allowed terms json')
    parser.add_argument('path_to_allowed_terms_json')
    parser.add_argument('--path_to_ncbi_dump', default = 'none')
    parser.add_argument('--path_to_uberon_owl', default = 'none')
    parser.add_argument('--path_to_po_owl', default = 'none')
    parser.add_argument('--path_to_cl_owl', default = 'none')
    parser.add_argument('--path_to_doid_owl', default = 'none')
    parser.add_argument('--path_to_ms_owl', default = 'none')
    parser.add_argument('--path_to_biome_envs_owl', default = 'none')
    parser.add_argument('--path_to_material_envs_owl', default = 'none')
    args = parser.parse_args()

    #ncbi_dump can be downloaded from https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdmp.zip 
    #after unzipping the file we need is named names.dmp

    #uberon_owl can be downloaded from http://aber-owl.net/ontology/UBERON/#/Download
    #always take the latest verion!

    #po_owl can be downloaded from http://aber-owl.net/ontology/PO/#/Download
    #always take the latest verion!

    #cl_owl can be downloaded from http://aber-owl.net/ontology/CL/#/Download
    #always take the latest verion!

    #doid_owl can be downloaded from http://aber-owl.net/ontology/DOID/#/Download
    #always take the latest verion!

    #ms_owl can be downloaded from http://aber-owl.net/ontology/MS/#/Download
    #always take the latest verion!

    #envo_bimoe can be downloaded from https://github.com/EnvironmentOntology/envo/blob/master/subsets/biome-hierarchy.owl
    #always take the latest verion!

    #envo_material can be downloaded from https://github.com/EnvironmentOntology/envo/blob/master/subsets/material-hierarchy.owl
    #always take the latest verion!


    # Load allowed Terms json
    print('Reading allowed terms json!')
    with open(args.path_to_allowed_terms_json, 'r') as file:
        allowedTerm_dict = json.load(file)


    #update from path_to_ncbi_dump
    ############
    if args.path_to_ncbi_dump != 'none':

        
        print('Updating NCBIs!')

        df = pd.read_csv(
            args.path_to_ncbi_dump,
            sep='\t\|\t',  
            engine='python',
            header=None,
            usecols=[0, 1, 2, 3],
            names=['ncbi_id', 'taxaname', 'unique_name', 'name_class'],
            dtype=str,
            lineterminator='\n'
    )
        
        df['name_class'] = df['name_class'].str.rstrip('\t|')

        # Now, filter for scientific names
        scientific_names_df = df[df['name_class'].str.strip() == 'scientific name'].copy()

        # Drop the 'unique_name' and 'name_class' columns as they are no longer needed
        scientific_names_df = scientific_names_df.drop(['unique_name', 'name_class'], axis=1)

        # Reset index
        scientific_names_df.reset_index(drop=True, inplace=True)

        # Make REDU formatted taxas 
        redu_formatted_ncbis = scientific_names_df.apply(lambda x: f'{x["ncbi_id"]}|{x["taxaname"]}', axis=1).to_frame(name='redu_ncbis')


        allowedTerm_dict['NCBITaxonomy']['allowed_values'] = redu_formatted_ncbis['redu_ncbis'].tolist()


    if args.path_to_uberon_owl != 'none' and args.path_to_po_owl != 'none' and args.path_to_cl_owl != 'none':
        
        print('Updating bodyparts!')

        uberon_df = get_uberon_table(args.path_to_uberon_owl)
        cl_df = get_ontology_table(args.path_to_cl_owl, ont_prefix = 'CL_')
        po_df = get_ontology_table(args.path_to_po_owl, ont_prefix = 'PO_', rm_synonym_info = True)


        combined_df = pd.concat([uberon_df, cl_df, po_df])

        # Get unique values for 'Label'
        unique_labels = combined_df['Label'].unique().tolist()

        # Get unique values for 'UBERONOntologyIndex'
        unique_ontology_indexes = combined_df['UBERONOntologyIndex'].str.replace("_", ":").unique().tolist()


        allowedTerm_dict['UBERONBodyPartName']['allowed_values'] = unique_labels
        allowedTerm_dict['UBERONOntologyIndex']['allowed_values'] = unique_ontology_indexes


    #environment ontology
    if args.path_to_biome_envs_owl != 'none' and args.path_to_material_envs_owl != 'none':

        print('Processing environmental biome ontology,..')
        envBiome_onto = get_ontology_table(args.path_to_biome_envs_owl, ont_prefix = 'ENVO_', index_column_name = 'ENVOEnvironmentBiomeIndex')
        unique_labels = envBiome_onto['Label'].unique().tolist()
        unique_indexes = envBiome_onto['ENVOEnvironmentBiomeIndex'].str.replace("_", ":").unique().tolist()

        allowedTerm_dict['ENVOEnvironmentBiome']['allowed_values'] = unique_labels
        allowedTerm_dict['ENVOEnvironmentBiomeIndex']['allowed_values'] = unique_indexes


        print('Processing environmental material ontology,..')
        envMaterial_onto = get_ontology_table(args.path_to_material_envs_owl, ont_prefix = 'ENVO_', index_column_name = 'ENVOEnvironmentMaterialIndex')
        unique_labels = envMaterial_onto['Label'].unique().tolist()
        unique_indexes = envMaterial_onto['ENVOEnvironmentMaterialIndex'].str.replace("_", ":").unique().tolist()

        allowedTerm_dict['ENVOEnvironmentMaterial']['allowed_values'] = unique_labels
        allowedTerm_dict['ENVOEnvironmentMaterialIndex']['allowed_values'] = unique_indexes


    if args.path_to_doid_owl != 'none':
        
        print('Updating diseases!')

        doid_df = get_ontology_table(args.path_to_doid_owl, ont_prefix = 'DOID_')

        # Get unique values for 'Label'
        unique_labels = doid_df['Label'].unique().tolist()

        # Get unique values for 'DOIDOntologyIndex' (ignore the wrong UBERONOntologyIndex name)
        unique_doid_indexes = doid_df['UBERONOntologyIndex'].str.replace("_", ":").unique().tolist()


        allowedTerm_dict['DOIDCommonName']['allowed_values'] = unique_labels
        allowedTerm_dict['DOIDOntologyIndex']['allowed_values'] = unique_doid_indexes



    if args.path_to_ms_owl != 'none':

        print('Updating mass spectrometers!')

        ms_df = get_ontology_table(args.path_to_ms_owl, ont_prefix = 'MS_', descendant_node='MS_1000031')
        ms_df.reset_index(drop=True, inplace=True)

        # Make REDU formatted taxas 
        ms_df = ms_df.apply(lambda x: f'{x["Label"]}|{x["UBERONOntologyIndex"].replace("_", ":")}', axis=1).to_frame(name='redu_ms')

        unique_labels = ms_df['redu_ms'].unique().tolist()

        allowedTerm_dict['MassSpectrometer']['allowed_values'] = unique_labels


    #print(f"Saving new allowed terms to {args.path_to_allowed_terms_json}")

    with open('allowed_terms.json' , 'w') as json_file:
        json.dump(allowedTerm_dict, json_file, indent=4)