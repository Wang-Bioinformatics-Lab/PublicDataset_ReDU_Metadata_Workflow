
from owlready2 import get_ontology
import owlready2
import owlready2.namespace as _owl_ns
import pandas as pd
import tqdm
import os
import argparse
import json
import re
from obo_tables import hpo_table, exposure_table

# owlready2 raises TypeError when an OWL file (e.g. new cl.owl via STATO import)
# declares an IRI as both a property and a class/individual. Patch _load_properties
# to skip only those conflicting entries instead of crashing. (Catching the TypeError
# around the original method would also skip loading every property after the
# conflicting one.) Mirrors Ontology._load_properties of owlready2 0.46.
def _patched_load_properties(self):
    # Update props from other ontologies, if needed
    for prop in list(self.world._props.values()):
        if prop.namespace.world is _owl_ns.owl_world: continue
        prop._check_update(self)

    # Loads new props
    for prop_storid in _owl_ns.itertools.chain(self._get_obj_triples_po_s(_owl_ns.rdf_type, _owl_ns.owl_object_property),
                                               self._get_obj_triples_po_s(_owl_ns.rdf_type, _owl_ns.owl_data_property),
                                               self._get_obj_triples_po_s(_owl_ns.rdf_type, _owl_ns.owl_annotation_property)):
        Prop = self.world._get_by_storid(prop_storid)
        if not isinstance(Prop, _owl_ns.PropertyClass):
            print(f"Skipping '{Prop.iri}': it is both a property and a class/an individual")
            continue
        python_name_d = self.world._get_data_triple_sp_od(prop_storid, _owl_ns.owlready_python_name)
        if python_name_d is not None:
            with _owl_ns.LOADING: Prop.python_name = python_name_d[0]
_owl_ns.Ontology._load_properties = _patched_load_properties

# Body substances UBERON does not place under 'bodily fluid' (UBERON_0006314) but that ReDU
# treats as biofluid (GNPS submitters label urine and saliva biofluid 100% of the time, feces 93%)
EXTRA_FLUIDS = {
    "UBERON_0001088",  # urine
    "UBERON_0001988",  # feces
    "UBERON_0001836",  # saliva
    "UBERON_0000173",  # amniotic fluid
    "UBERON_0000038",  # follicular fluid
    "UBERON_0001968",  # semen
    "UBERON_0006530",  # seminal fluid
    "UBERON_0001971",  # gastric juice
}

# Common names that are not UBERON synonyms (e.g. MetaboLights and Workbench use 'plasma')
EXTRA_SYNONYMS = {
    "UBERON_0001969": ["plasma"],
    "UBERON_0000178": ["whole blood"],
    "UBERON_0001988": ["fecal", "faecal"],
}

def get_uberon_table(owl_path):
    # Each ontology gets its own owlready2 World: in the shared default world, classes an
    # earlier file already declared (e.g. UBERON declares ENVO's 'marine biome') belong to
    # that earlier ontology and silently disappear from this one's classes()
    onto = owlready2.World().get_ontology(owl_path).load()

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
        synonyms = synonyms + [s for s in EXTRA_SYNONYMS.get(cls.name, []) if s not in synonyms]

        is_multi_cellular = cls in downstream_classes_multi_cellular
        is_organ = cls in downstream_classes_organ
        is_bodily_fluid = cls in downstream_classes_bodily_fluid or cls.name in EXTRA_FLUIDS

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
    # obsolete classes keep their label as 'obsolete <name>'; they must not be assigned
    df = df[~df['Label'].astype(str).str.lower().str.startswith('obsolete')]

    return df


# Body parts only one sex has (BiologicalSex from the body part, bin/final_table_rules.py): the
# classes below the female / male reproductive system or organism (UBERON) or germ cell (CL), by
# is_a, part_of, produced_by (cervical mucus is produced by the uterine cervix) or has_part (semen
# has part sperm). Embryonic and extraembryonic structures belong to the conceptus, which can be of
# either sex (UBERON puts the placenta in the female reproductive system): left out, as are classes
# below both sexes (e.g. 'synchronous hermaphroditic organism').
SEX_ROOTS = {
    'female': {'UBERON_0000474', 'UBERON_0003100', 'CL_0000021'},  # reproductive system, organism, germ cell
    'male': {'UBERON_0000079', 'UBERON_0003101', 'CL_0000015'},
}
EITHER_SEX_ROOTS = {'UBERON_0000922', 'UBERON_0000478', 'UBERON_0016887'}  # embryo, extraembryonic structure and component
SEX_RELATIONS = {'BFO_0000050', 'RO_0003001', 'BFO_0000051'}  # part_of, produced_by, has_part
# milk and colostrum (mammary gland fluid/secretion) come from the mammary gland, which both sexes
# have, but only females lactate; 'Is Milk' lets the final table rule tell the subject's own milk
# from the milk an infant subject drank
MILK = 'UBERON_0006539'
# errors of the pinned UBERON, left out: both sexes have a cranial suspensory ligament (UBERON puts it
# in the male reproductive system), a mating plug is made by the male but sampled from the female,
# 'pars recta' is an amphibian oviduct part but in mammals the name of the kidney's proximal straight
# tubule. A label that names the other sex is left out too ('female deep inguinal ring' is part_of
# the male reproductive system)
SEX_EXCLUDED = {'UBERON_0008848', 'UBERON_0010148', 'UBERON_3010780'}
SEX_WORD = re.compile(r'\b(female|male)\b', re.I)


def get_sex_specific_parts(uberon_owl_path, cl_owl_path):
    """UBERON and CL classes only one sex has (SEX_ROOTS, EITHER_SEX_ROOTS, MILK, SEX_EXCLUDED):
    UBERONOntologyIndex, Label, BiologicalSex, Is Milk"""
    rows = []
    for owl_path, prefix in ((uberon_owl_path, 'UBERON_'), (cl_owl_path, 'CL_')):
        onto = owlready2.World().get_ontology(owl_path).load()
        children, labels = {}, {}
        for cls in onto.classes():
            labels[cls.name] = cls.label.first() if cls.label else None
            for parent in cls.is_a:
                if isinstance(parent, owlready2.Restriction):
                    if parent.type != owlready2.SOME or getattr(parent.property, 'name', None) not in SEX_RELATIONS:
                        continue
                    parent = parent.value
                if isinstance(parent, owlready2.ThingClass):
                    children.setdefault(parent.name, set()).add(cls.name)

        def below(roots):
            seen, todo = set(roots), list(roots)
            while todo:
                for child in children.get(todo.pop(), ()):
                    if child not in seen:
                        seen.add(child)
                        todo.append(child)
            return seen
        milk = below({MILK})
        female, male, either = below(SEX_ROOTS['female']) | milk, below(SEX_ROOTS['male']), below(EITHER_SEX_ROOTS)
        for name, label in labels.items():
            if not name.startswith(prefix) or not label or str(label).lower().startswith('obsolete') or name in either \
                    or name in SEX_EXCLUDED:
                continue
            if (name in female) != (name in male):
                sex = 'female' if name in female else 'male'
                if {w.lower() for w in SEX_WORD.findall(str(label))} - {sex}:
                    continue
                rows.append({'UBERONOntologyIndex': name.replace('_', ':'), 'Label': str(label),
                             'BiologicalSex': sex, 'Is Milk': name in milk})
    df = pd.DataFrame(rows, columns=['UBERONOntologyIndex', 'Label', 'BiologicalSex', 'Is Milk'])
    return df.sort_values('UBERONOntologyIndex', kind='stable').reset_index(drop=True)


def get_doid_parents(owl_path):
    """Direct is_a links between DOID classes (used to pick the most specific disease term):
    asserted subclass links plus the genus of equivalent-class definitions (e.g.
    adenocarcinoma = cancer and has_material_basis_in ..., so adenocarcinoma is_a cancer)."""
    onto = owlready2.World().get_ontology(owl_path).load()
    rows = []
    for cls in onto.classes():
        if not cls.name.startswith('DOID_') or True in getattr(cls, 'deprecated', []):
            continue
        parents = list(cls.is_a)
        for definition in cls.equivalent_to:
            parents += list(getattr(definition, 'Classes', []))
        for parent in parents:
            name = getattr(parent, 'name', None)
            if isinstance(name, str) and name.startswith('DOID_') and name != cls.name:
                rows.append({'DOIDOntologyIndex': cls.name.replace('_', ':'), 'ParentDOIDOntologyIndex': name.replace('_', ':')})
    return pd.DataFrame(rows, columns=['DOIDOntologyIndex', 'ParentDOIDOntologyIndex']).drop_duplicates()


def get_ontology_table(owl_path, ont_prefix, rm_synonym_info=False, descendant_node=None, index_column_name = 'UBERONOntologyIndex'):
    # Each ontology gets its own owlready2 World: in the shared default world, classes an
    # earlier file already declared (e.g. UBERON declares ENVO's 'marine biome') belong to
    # that earlier ontology and silently disappear from this one's classes()
    onto = owlready2.World().get_ontology(owl_path).load()
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
    # obsolete classes keep their label as 'obsolete <name>'; they must not be assigned
    df = df[~df['Label'].astype(str).str.lower().str.startswith('obsolete')]

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
    parser.add_argument('--path_to_hp_obo', default='none')
    parser.add_argument('--path_to_ecto_obo', default='none')
    parser.add_argument('--path_to_peco_obo', default='none')

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

    print('Processing body parts only one sex has,..')
    sex_specific_parts = get_sex_specific_parts(args.path_to_uberon_owl, args.path_to_cl_owl)

    #environment ontology
    print('Processing environmental biome ontology,..')
    envBiome_onto = get_ontology_table(args.path_to_biome_envs_owl, ont_prefix = 'ENVO_', index_column_name = 'ENVOEnvironmentBiomeIndex')
    envBiome_onto['ENVOEnvironmentBiomeIndex'] = envBiome_onto['ENVOEnvironmentBiomeIndex'].str.replace('_', ':')


    print('Processing environmental material ontology,..')
    envMaterial_onto = get_ontology_table(args.path_to_material_envs_owl, ont_prefix = 'ENVO_', index_column_name = 'ENVOEnvironmentMaterialIndex')
    envMaterial_onto['ENVOEnvironmentMaterialIndex'] = envMaterial_onto['ENVOEnvironmentMaterialIndex'].str.replace('_', ':')


    #create ontology table to fill DOIDOntologyIndex from bodyparts
    print('Processing disease ontology,..')
    doid_ontology_table = get_ontology_table(args.path_to_doid_owl, ont_prefix = 'DOID_', index_column_name = 'DOIDOntologyIndex')
    doid_ontology_table['DOIDOntologyIndex'] = doid_ontology_table['DOIDOntologyIndex'].str.replace('_', ':')
    doid_parents_table = get_doid_parents(args.path_to_doid_owl)




    # health conditions (HPO) and exposures (ECTO, PECO); empty tables if a file is missing
    print('Processing HPO, ECTO and PECO,..')
    hpo_ontology_table = hpo_table(args.path_to_hp_obo)
    exposure_ontology_table = exposure_table(args.path_to_ecto_obo, args.path_to_peco_obo)

    #create ncbi - rank - division sheet
    print('Loading and processing NCBI rank file,..')
    df_ncbi_rank = pd.read_csv(args.path_to_ncbi_nodes_dmp, sep='\t\|\t', engine='python', header=None, index_col=False, comment='#')
    df_ncbi_rank = df_ncbi_rank.dropna(axis=1, how='all')
    df_ncbi_rank.columns = ['TaxonID', 'ParentID', 'NCBIRank', 'Empties1', 'DivID', 'InheritedDivFlag', 'GeneticCodeID', 'InheritedGCFlag', 'MitochondrialGeneticCodeID', 'InheritedMitoGCFlag', 'GenBankHiddenFlag', 'HiddenSubtreeRootFlag', 'Comments']
    df_ncbi_rank = df_ncbi_rank[['TaxonID', 'NCBIRank', 'DivID']]
    # NCBI renamed the 'superkingdom' rank to 'domain' (2025); ReDU keeps the established term
    df_ncbi_rank['NCBIRank'] = df_ncbi_rank['NCBIRank'].replace('domain', 'superkingdom')

    print('Loading and processing NCBI divion file,..')
    df_ncbi_divisions = pd.read_csv(args.path_to_ncbi_division_dmp, sep='\t\|\t', engine='python', header=None, index_col=False)
    df_ncbi_divisions[3] = df_ncbi_divisions[3].str.replace('\t\|', '', regex=False).str.strip()
    df_ncbi_divisions.columns = ['DivID', 'Abbreviation', 'NCBIDivision', 'Notes']
    df_ncbi_divisions = df_ncbi_divisions[['DivID', 'NCBIDivision']]

    print('Merging NCBI files,..')
    df_ncbi_rank_divisions = pd.merge(df_ncbi_rank, df_ncbi_divisions, on='DivID', how='left')
    df_ncbi_rank_divisions = df_ncbi_rank_divisions[['TaxonID', 'NCBIRank', 'NCBIDivision']]


    uberon_ontology_table.to_csv('UBERON_CL_PO_ontology.csv', index=False)
    sex_specific_parts.to_csv('UBERON_sex_specific_parts.csv', index=False)
    doid_ontology_table.to_csv('DOID_ontology.csv', index=False)
    doid_parents_table.to_csv('DOID_parents.csv', index=False)
    hpo_ontology_table.to_csv('HPO_ontology.csv', index=False)
    exposure_ontology_table.to_csv('Exposure_ontology.csv', index=False)
    envBiome_onto.to_csv('ENVO_biome_ontology.csv', index=False)
    envMaterial_onto.to_csv('ENVO_material_ontology.csv', index=False)
    df_ncbi_rank_divisions.to_csv('NCBI_Rank_Division.csv', index=False)