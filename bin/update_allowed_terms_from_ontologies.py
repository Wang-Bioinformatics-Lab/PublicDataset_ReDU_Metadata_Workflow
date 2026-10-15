import pandas as pd
import argparse
import json
from REDU_conversion_functions import get_uberon_table
from REDU_conversion_functions import get_ontology_table
from obo_tables import hpo_table, exposure_table
from diet_columns import NOT_DIET_FOODS

def exact_synonym_labels(table):
    """{exact synonym (lower case): label} for the synonyms of a Label/Synonym table that
    belong to exactly one label and are not (in any case) the label of another class, and
    {label (lower case): label} for labels no other label equals in lower case.
    Synonyms of 3 characters or less are left out: they are abbreviations ('NA' = nasal
    artery) that collide with placeholders. Sorted by synonym, so the json does not depend
    on the order of the ontology files."""
    labels_of = {}
    lower_labels = {}
    for label in table['Label'].dropna().astype(str).unique():
        lower_labels.setdefault(label.strip().lower(), set()).add(label)
    for label, synonym in zip(table['Label'], table['Synonym']):
        if not isinstance(label, str) or not isinstance(synonym, str):
            continue
        key = synonym.strip().lower()
        if len(key) <= 3 or key == label.strip().lower():
            continue
        labels_of.setdefault(key, set()).add(label)
    out = {key: next(iter(labels)) for key, labels in labels_of.items()
           if len(labels) == 1 and lower_labels.get(key, labels) == labels}
    # the label in other case ('islet of langerhans' -> 'islet of Langerhans'), when only
    # one label is written so
    for key, labels in lower_labels.items():
        if len(labels) == 1:
            out[key] = next(iter(labels))
    return dict(sorted(out.items()))


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
    parser.add_argument('--path_to_hp_obo', default = 'none')
    parser.add_argument('--path_to_ecto_obo', default = 'none')
    parser.add_argument('--path_to_peco_obo', default = 'none')
    parser.add_argument('--path_to_cellosaurus_tsv', default = 'none')
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
        # exact synonyms -> label, so a value that is an exact synonym of one class is kept
        # under its current label (e.g. 'mononuclear cell', the former label of CL:0000842,
        # now its exact synonym under the label 'mononuclear leukocyte')
        allowedTerm_dict['UBERONBodyPartName']['synonyms'] = exact_synonym_labels(combined_df)
        print(f"UBERONBodyPartName: {len(allowedTerm_dict['UBERONBodyPartName']['synonyms'])} exact synonyms of one class")


    #environment ontology
    if args.path_to_biome_envs_owl != 'none' and args.path_to_material_envs_owl != 'none':

        print('Processing environmental biome ontology,..')
        envBiome_onto = get_ontology_table(args.path_to_biome_envs_owl, ont_prefix = 'ENVO_', index_column_name = 'ENVOEnvironmentBiomeIndex')
        unique_labels = envBiome_onto['Label'].unique().tolist()
        unique_indexes = envBiome_onto['ENVOEnvironmentBiomeIndex'].str.replace("_", ":").unique().tolist()

        allowedTerm_dict['ENVOEnvironmentBiome']['allowed_values'] = unique_labels
        allowedTerm_dict['ENVOEnvironmentBiomeIndex']['allowed_values'] = unique_indexes

        print('Processing environmental biome local scale ontology,..')
        allowedTerm_dict['ENVOLocalScale']['allowed_values'] = unique_labels
        allowedTerm_dict['ENVOLocalScaleIndex']['allowed_values'] = unique_indexes

        print('Processing environmental biome broad scale ontology,..')
        allowedTerm_dict['ENVOBroadScale']['allowed_values'] = unique_labels
        allowedTerm_dict['ENVOBroadScaleIndex']['allowed_values'] = unique_indexes


        print('Processing environmental material ontology,..')
        envMaterial_onto = get_ontology_table(args.path_to_material_envs_owl, ont_prefix = 'ENVO_', index_column_name = 'ENVOEnvironmentMaterialIndex')
        unique_labels = envMaterial_onto['Label'].unique().tolist()
        unique_indexes = envMaterial_onto['ENVOEnvironmentMaterialIndex'].str.replace("_", ":").unique().tolist()

        allowedTerm_dict['ENVOEnvironmentMaterial']['allowed_values'] = unique_labels
        allowedTerm_dict['ENVOEnvironmentMaterialIndex']['allowed_values'] = unique_indexes

        print('Processing environmental material medium scale ontology,..')
        allowedTerm_dict['ENVOMediumScale']['allowed_values'] = unique_labels
        allowedTerm_dict['ENVOMediumScaleIndex']['allowed_values'] = unique_indexes


    if args.path_to_doid_owl != 'none':
        
        print('Updating diseases!')

        doid_df = get_ontology_table(args.path_to_doid_owl, ont_prefix = 'DOID_')

        # Get unique values for 'Label'
        unique_labels = doid_df['Label'].unique().tolist()

        # Get unique values for 'DOIDOntologyIndex' (ignore the wrong UBERONOntologyIndex name)
        unique_doid_indexes = doid_df['UBERONOntologyIndex'].str.replace("_", ":").unique().tolist()


        allowedTerm_dict['DOIDCommonName']['allowed_values'] = unique_labels
        allowedTerm_dict['DOIDOntologyIndex']['allowed_values'] = unique_doid_indexes



    #hp_obo, ecto_obo and peco_obo: http://purl.obolibrary.org/obo/{hp,ecto,peco}.obo (fetched by the workflow)
    if args.path_to_hp_obo != 'none':
        print('Updating HPO health conditions!')
        hpo_df = hpo_table(args.path_to_hp_obo)
        if len(hpo_df) > 0:
            allowedTerm_dict['HPOCommonName']['allowed_values'] = hpo_df['Label'].unique().tolist()
            allowedTerm_dict['HPOOntologyIndex']['allowed_values'] = hpo_df['HPOOntologyIndex'].unique().tolist()
        else:
            print('WARNING: no HPO terms found, HPOCommonName keeps its previous allowed values')

    if args.path_to_ecto_obo != 'none' or args.path_to_peco_obo != 'none':
        print('Updating exposures!')
        exposure_df = exposure_table(args.path_to_ecto_obo, args.path_to_peco_obo)
        if len(exposure_df) > 0:
            allowedTerm_dict['Exposure']['allowed_values'] = exposure_df['Label'].unique().tolist()
            allowedTerm_dict['ExposureIndex']['allowed_values'] = exposure_df['ExposureIndex'].unique().tolist()
            agents = exposure_df.dropna(subset=['OrganismExposureAgentIndex']) if 'OrganismExposureAgentIndex' in exposure_df.columns else exposure_df.iloc[0:0]
            for key, column in (('OrganismExposureAgent', 'OrganismExposureAgent'), ('OrganismExposureAgentIndex', 'OrganismExposureAgentIndex')):
                if key in allowedTerm_dict and len(agents):
                    allowedTerm_dict[key]['allowed_values'] = agents[column].unique().tolist()
        else:
            print('WARNING: no ECTO/PECO terms found, Exposure keeps its previous allowed values')

    #cellosaurus_tsv: the pinned table shipped in data/ (built by build_cellosaurus_table.py)
    if args.path_to_cellosaurus_tsv != 'none':
        print('Updating cell lines!')
        from cellosaurus import load_cellosaurus_table
        cello_df = load_cellosaurus_table(args.path_to_cellosaurus_tsv)
        if cello_df is not None and len(cello_df) > 0:
            first = cello_df.drop_duplicates(subset=['CellLineIndex'])
            allowedTerm_dict['CellLine']['allowed_values'] = first['CellLine'].unique().tolist()
            allowedTerm_dict['CellLineIndex']['allowed_values'] = first['CellLineIndex'].tolist()
        else:
            print('WARNING: no cell lines found, CellLine keeps its previous allowed values')

    if args.path_to_ms_owl != 'none':

        print('Updating mass spectrometers!')

        ms_df = get_ontology_table(args.path_to_ms_owl, ont_prefix = 'MS_', descendant_node='MS_1000031')
        ms_df.reset_index(drop=True, inplace=True)

        # Make REDU formatted taxas 
        ms_df = ms_df.apply(lambda x: f'{x["Label"]}|{x["UBERONOntologyIndex"].replace("_", ":")}', axis=1).to_frame(name='redu_ms')

        unique_labels = ms_df['redu_ms'].unique().tolist()

        allowedTerm_dict['MassSpectrometer']['allowed_values'] = unique_labels


    # Diet: the regimen terms of the seed, then the food terms (the food a subject ate), which are
    # the FoodOntologyTerm (GFOP) terms, so both columns use the same food vocabulary; without the
    # GFOP classes too general for a diet and water (diet_columns.NOT_DIET_FOODS)
    if 'Diet' in allowedTerm_dict and 'FoodOntologyTerm' in allowedTerm_dict:
        gfop = allowedTerm_dict['FoodOntologyTerm']['allowed_values']
        regimens = [t for t in allowedTerm_dict['Diet']['allowed_values'] if t not in set(gfop)]
        foods = [t for t in gfop if t not in NOT_DIET_FOODS]
        allowedTerm_dict['Diet']['allowed_values'] = regimens + foods
        print(f"Diet: {len(regimens)} regimen terms and {len(foods)} food terms")

    #print(f"Saving new allowed terms to {args.path_to_allowed_terms_json}")

    with open('allowed_terms.json' , 'w') as json_file:
        json.dump(allowedTerm_dict, json_file, indent=4)