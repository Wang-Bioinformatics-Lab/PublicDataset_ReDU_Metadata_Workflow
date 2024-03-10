
from owlready2 import get_ontology
import owlready2
import pandas as pd
import tqdm
import os
import argparse
import json

def get_uberon_table(owl_path):
    onto = get_ontology(owl_path).load()

    multi_cellular = onto.search(iri="http://purl.obolibrary.org/obo/UBERON_0010000")[0]
    organ = onto.search(iri="http://purl.obolibrary.org/obo/UBERON_0000062")[0]
    bodily_fluid = onto.search(iri="http://purl.obolibrary.org/obo/UBERON_0006314")[0]

    downstream_classes_multi_cellular = set(multi_cellular.descendants())
    downstream_classes_organ = set(organ.descendants())
    downstream_classes_bodily_fluid = set(bodily_fluid.descendants())
            
    data = []
    for cls in tqdm.tqdm(onto.classes(), desc="Processing classes"):
        label = cls.label.first() if cls.label else None
        synonyms = [synonym for synonym in cls.hasExactSynonym] if hasattr(cls, 'hasExactSynonym') else []

        if not synonyms:
            synonyms = [label] if label else []
            
        is_multi_cellular = cls in downstream_classes_multi_cellular
        is_organ = cls in downstream_classes_organ
        is_bodily_fluid = cls in downstream_classes_bodily_fluid

        for synonym in synonyms:
            data.append({
                "UBERONOntologyIndex": cls.name,
                "Label": label,
                "Synonym": synonym,
                "Is Multicellular": is_multi_cellular,
                "Is Organ": is_organ,
                "Is Fluid": is_bodily_fluid
            })

    df = pd.DataFrame(data)

    df = df[(df["Label"].notna() | df["Synonym"].notna()) & ((df["Label"] != '') | (df["Synonym"] != ''))]
    df = df[df['UBERONOntologyIndex'].str.startswith('UBERON_', na=False)]

    return df



def get_ontology_table(owl_path, ont_prefix, rm_synonym_info=False, descendant_node=None, index_column_name = 'UBERONOntologyIndex'):
    onto = get_ontology(owl_path).load()

    filter_class = None
    descendants = set()
    if descendant_node:
        # Search for the specific class using its node ID
        filter_class = onto.search_one(iri="*" + descendant_node)
        if filter_class:
            # Get all its descendants
            descendants = set(filter_class.descendants())
            # Optionally remove the filter_class itself from the set of descendants to exclude it
            descendants.discard(filter_class)

    data = []
    for cls in tqdm.tqdm(onto.classes(), desc="Processing classes"):
        # Skip the class if it is the filter_class itself and descendant_node is provided
        if descendant_node and cls == filter_class:
            continue

        label = cls.label.first() if cls.label else None
        synonyms = [synonym for synonym in cls.hasExactSynonym] if hasattr(cls, 'hasExactSynonym') else []

        if not synonyms:
            synonyms = [label] if label else []

        # Determine if the current class is a descendant of the specified node
        is_descendant = cls in descendants

        for synonym in synonyms:
            data.append({
                index_column_name: cls.name,
                "Label": label,
                "Synonym": synonym,
                "Is Descendant": is_descendant  # Add this information to each entry
            })

    df = pd.DataFrame(data)
    df = df[(df["Label"].notna() | df["Synonym"].notna()) & ((df["Label"] != '') | (df["Synonym"] != ''))]
    df = df[df[index_column_name].str.startswith(ont_prefix, na=False)]

    if rm_synonym_info:
        df['Synonym'] = df['Synonym'].str.replace(r' \([^)]*\)', '', regex=True)

    # Optionally, you might want to filter the dataframe to include only the descendants
    if descendant_node:
        df = df[df["Is Descendant"] == True]

    return df


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Prepare ontologies')
    parser.add_argument('--path_to_uberon_owl')
    parser.add_argument('--path_to_po_owl')
    parser.add_argument('--path_to_cl_owl')
    parser.add_argument('--path_to_doid_owl')
    parser.add_argument('--path_to_biome_envs_owl')
    parser.add_argument('--path_to_material_envs_owl')
    parser.add_argument('--path_to_ncbi_nodes_dmp')
    parser.add_argument('--path_to_ncbi_division_dmp')

    args = parser.parse_args()

    #python3.8 ../ReDU-MS2-GNPS2/workflows/PublicDataset_ReDU_Metadata_Workflow/bin/read_and_validate_redu_from_github.py /home/yasin/projects/ReDU_metadata/metadata output/ --AllowedTermJson_path /home/yasin/projects/ReDU-MS2-GNPS2/workflows/PublicDataset_ReDU_Metadata_Workflow/bin/allowed_terms/allowed_terms.json --path_to_uberon_owl /home/yasin/projects/ReDU-MS2-GNPS2/workflows/PublicDataset_ReDU_Metadata_Workflow/bin/allowed_terms/uberon.owl --path_to_po_owl /home/yasin/projects/ReDU-MS2-GNPS2/workflows/PublicDataset_ReDU_Metadata_Workflow/bin/allowed_terms/po.owl --path_to_cl_owl /home/yasin/projects/ReDU-MS2-GNPS2/workflows/PublicDataset_ReDU_Metadata_Workflow/bin/allowed_terms/cl.owl  --path_to_doid_owl /home/yasin/projects/ReDU-MS2-GNPS2/workflows/PublicDataset_ReDU_Metadata_Workflow/bin/allowed_terms/doid.owl


    #create ontology table to fill UBERONOntologyIndex from bodyparts
    print('Processing body part ontology,..')
    uberon_onto = get_uberon_table(args.path_to_uberon_owl)

    print('Processing cell ontology,..')
    cl_onto = get_ontology_table(args.path_to_cl_owl, ont_prefix = 'CL_')

    print('Processing plant ontology,..')
    po_onto = get_ontology_table(args.path_to_po_owl, ont_prefix = 'PO_', rm_synonym_info = True)

    uberon_ontology_table = pd.concat([uberon_onto, cl_onto, po_onto], ignore_index=True, sort=False)
    uberon_ontology_table['UBERONOntologyIndex'] = uberon_ontology_table['UBERONOntologyIndex'].str.replace('_', ':')




    #environment ontology
    print('Processing environmental biome ontology,..')
    envBiome_onto = get_uberon_table(args.path_to_biome_envs_owl, ont_prefix = 'ENVO_', index_column_name = 'ENVOEnvironmentBiomeIndex')
    envBiome_onto['ENVOEnvironmentBiomeIndex'] = envBiome_onto['ENVOEnvironmentBiomeIndex'].str.replace('_', ':')



    print('Processing environmental material ontology,..')
    envMaterial_onto = get_ontology_table(args.path_to_material_envs_owl, ont_prefix = 'ENVO_', index_column_name = 'ENVOEnvironmentMaterialIndex')
    envMaterial_onto['ENVOEnvironmentMaterialIndex'] = envMaterial_onto['ENVOEnvironmentMaterialIndex'].str.replace('_', ':')


    #create ontology table to fill DOIDOntologyIndex from bodyparts
    print('Processing disease ontology,..')
    doid_ontology_table = get_ontology_table(args.path_to_doid_owl, ont_prefix = 'DOID_', index_column_name = 'DOIDOntologyIndex')
    doid_ontology_table['DOIDOntologyIndex'] = doid_ontology_table['DOIDOntologyIndex'].str.replace('_', ':')




    #create ncbi - rank - division sheet
    print('Loading and processing NCBI rank file,..')
    df_ncbi_rank = pd.read_csv(args.path_to_ncbi_nodes_dmp, sep='\t\|\t', engine='python', header=None, index_col=False, comment='#')
    df_ncbi_rank = df_ncbi_rank.dropna(axis=1, how='all')
    df_ncbi_rank.columns = ['TaxonID', 'ParentID', 'NCBIRank', 'Empties1', 'DivID', 'InheritedDivFlag', 'GeneticCodeID', 'InheritedGCFlag', 'MitochondrialGeneticCodeID', 'InheritedMitoGCFlag', 'GenBankHiddenFlag', 'HiddenSubtreeRootFlag', 'Comments']
    df_ncbi_rank = df_ncbi_rank[['TaxonID', 'NCBIRank', 'DivID']]

    print('Loading and processing NCBI divion file,..')
    df_ncbi_divisions = pd.read_csv(args.path_to_ncbi_division_dmp, sep='\t\|\t', engine='python', header=None, index_col=False)
    df_ncbi_divisions[3] = df_ncbi_divisions[3].str.replace('\t\|', '', regex=False).str.strip()
    df_ncbi_divisions.columns = ['DivID', 'Abbreviation', 'NCBIDivision', 'Notes']
    df_ncbi_divisions = df_ncbi_divisions[['DivID', 'NCBIDivision']]

    print('Merging NCBI files,..')
    df_ncbi_rank_divisions = pd.merge(df_ncbi_rank, df_ncbi_divisions, on='DivID', how='left')
    df_ncbi_rank_divisions = df_ncbi_rank_divisions[['TaxonID', 'NCBIRank', 'NCBIDivision']]


    uberon_ontology_table.to_csv('UBERON_CL_PO_ontology.csv', index=False)
    doid_ontology_table.to_csv('DOID_ontology.csv', index=False)
    envBiome_onto.to_csv('ENVO_biome_ontology.csv', index=False)
    envMaterial_onto.to_csv('ENVO_material_ontology.csv', index=False)
    df_ncbi_rank_divisions.to_csv('NCBI_Rank_Division.csv', index=False)