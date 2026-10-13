"""MASST tables (https://github.com/robinschmid/microbe_masst/tree/master/data) -> ReDU rows
for files that no GNPS/ReDU metadata file describes.

Tables (each optional: a missing, empty or unreadable table is skipped with a warning):
  microbe_masst_table.csv               microbial cultures (NCBI taxonomy of the strain)
  plant_masst_table.csv                 plants
  food_masst_table.csv + _full.tsv      Global FoodOmics (foodMASST): food/beverage, organism,
                                        plant part, country, coordinates, date, altitude,
                                        GFOP food ontology term
  personalCareProduct_masst_table.csv   personal care / household products and their blanks
  microbiome_masst_table.tsv            microbiomeMASST, only rows not copied from ReDU
                                        (imported_from_redu != yes)
tissueMASST and plant_masst_table_full are ReDU exports and are not read.

Precedence for a file listed in several tables (the more specific table wins):
  food > personal care > microbe > microbiome > plant
e.g. the Global FoodOmics files (MSV000084900) that plantMASST also lists are food. Files
described by GNPS/ReDU metadata files are left out here, and GNPS rows win the merge
(merge_metadata.py) for the same USI.
"""
import argparse
import collections
import json
import os

import pandas as pd

from read_and_validate_redu_from_github import complete_and_fill_REDU_table, build_ncbi_name_resolver
from masst_tables import (read_table, food_rows, personal_care_rows, microbiome_rows, Vocab, load_gfop)

# NCBI taxids used to tell microbes apart (Bacteria division also holds archaea,
# 'Plants and Fungi' holds both)
TAXID_ARCHAEA = '2157'
TAXID_FUNGI = '4751'

# most specific first
TABLE_PRECEDENCE = ['food', 'personal_care', 'microbe', 'microbiome', 'plant']


def load_ncbi_parents(path_to_nodes_dmp):
    """taxid -> parent taxid from NCBI nodes.dmp"""
    parents = {}
    with open(path_to_nodes_dmp) as f:
        for line in f:
            taxid, parent = line.split('\t|\t', 2)[:2]
            parents[taxid] = parent
    return parents


def has_ancestor(taxid, ancestor, parents):
    seen = set()
    while taxid and taxid not in seen:
        if taxid == ancestor:
            return True
        seen.add(taxid)
        parent = parents.get(taxid)
        if parent == taxid:
            break
        taxid = parent
    return False


def microbe_sample_type(taxid, division_by_taxid, parents):
    """SampleType/SampleTypeSub1 of a microbeMASST culture, or ('', '') if unclear"""
    division = division_by_taxid.get(taxid)
    if division == 'Bacteria' and parents is not None:
        if has_ancestor(taxid, TAXID_ARCHAEA, parents):
            return 'archaea', ''
        return 'culture_bacterial', 'culture_bacterial'
    if division == 'Plants and Fungi' and parents is not None and has_ancestor(taxid, TAXID_FUNGI, parents):
        return 'culture_fungal', 'culture_fungal'
    # microbeMASST also holds mammalian cell cultures (e.g. human cell lines)
    if division in ('Primates', 'Rodents', 'Mammals'):
        return 'culture_mammalian', 'culture_mammalian'
    return '', ''


def apply_blank_qc_flags(df):
    """Blank/QC files carry no organism: label them as blank/QC and drop the taxonomy"""
    blank = df['Blank'].astype(str).str.strip().str.lower().eq('yes')
    qc = df['QC'].astype(str).str.strip().str.lower().eq('yes') & ~blank
    df.loc[blank, 'SampleType'] = 'blank_analysis'
    df.loc[blank, 'SampleTypeSub1'] = 'blank_analysis'
    df.loc[qc, 'SampleType'] = 'pool_QC'
    df.loc[qc, 'SampleTypeSub1'] = 'pool_QC'
    df.loc[blank | qc, 'NCBITaxonomy'] = ''
    print(f"Blank/QC files: {blank.sum()} blanks, {qc.sum()} QCs")
    return df


def _stem(filename):
    base = os.path.basename(str(filename))
    root, ext = os.path.splitext(base)
    return (root if ext.lower() in ('.mzml', '.mzxml') else base).lower()


def _collapse_repeats(df, name):
    """One row per file of the microbe/plant tables: a file listed several times (e.g.
    composite foods under each ingredient) keeps a column only where its rows agree."""
    key = ['MassiveID', 'filename']
    df = df.fillna('')
    repeated = df.duplicated(key, keep=False)
    if not repeated.any():
        return df
    def agree(values):
        vals = list(dict.fromkeys(v for v in values if str(v).strip()))
        return vals[0] if len(vals) == 1 else ''
    collapsed = df[repeated].groupby(key, sort=False).agg(agree).reset_index()
    print(f"{name}: {repeated.sum()} rows for {len(collapsed)} repeated files collapsed (conflicting fields left empty)")
    return pd.concat([df[~repeated], collapsed[df.columns]], ignore_index=True)


def main():
    parser = argparse.ArgumentParser(description='MASST tables to ReDU')
    parser.add_argument('path_REDUtable')
    parser.add_argument('output_folder')
    parser.add_argument('--path_microbeMASST')
    parser.add_argument('--path_plantMASST')
    parser.add_argument('--path_foodMASST', default=None, help='food_masst_table.csv (dataset of each file)')
    parser.add_argument('--path_foodMASST_full', default=None, help='food_masst_table_full.tsv (Global FoodOmics sample sheet)')
    parser.add_argument('--path_personalCareProductMASST', default=None)
    parser.add_argument('--path_microbiomeMASST', default=None)
    parser.add_argument('--path_ncbiRanksDivisions')
    parser.add_argument('--path_ncbi_nodes', default=None, help='NCBI nodes.dmp, to tell bacteria/archaea and fungi/plants apart')
    parser.add_argument('--path_ncbi_names', default=None, help='NCBI names.dmp: common names of food animals, synonyms')
    parser.add_argument('--path_to_uberon_cl_po_csv', default=None)
    parser.add_argument('--path_to_doid_csv', default=None)
    parser.add_argument('--path_to_exposure_csv', default=None)
    parser.add_argument("--AllowedTermJson_path", type=str, help="Path to json with allowed terms")
    args = parser.parse_args()

    df_redu = pd.read_csv(args.path_REDUtable, sep='\t', dtype=str)
    df_redu['filename'] = df_redu['filename'].str.slice(2)
    # Files already described by GNPS metadata, identified by dataset + path (the same
    # relative path, e.g. ccms_peak/Blank.mzML, exists in many datasets); tables that only
    # give the file name are compared by dataset + name
    redu_files = set(zip(df_redu['ATTRIBUTE_DatasetAccession'], df_redu['filename']))
    redu_stems = set(zip(df_redu['ATTRIBUTE_DatasetAccession'], df_redu['filename'].map(_stem)))

    with open(args.AllowedTermJson_path, 'r') as json_file:
        allowed_terms = json.load(json_file)

    NCBIRankDivision_table = pd.read_csv(args.path_ncbiRanksDivisions, index_col=False)
    NCBIRankDivision_table = NCBIRankDivision_table.drop_duplicates(subset=['TaxonID'])
    division_by_taxid = dict(zip(NCBIRankDivision_table['TaxonID'].astype(str), NCBIRankDivision_table['NCBIDivision']))
    parents = load_ncbi_parents(args.path_ncbi_nodes) if args.path_ncbi_nodes else None

    def _csv(path):
        return pd.read_csv(path) if path and path != 'none' and os.path.exists(path) else None
    uberon_table, doid_table, exposure_table = (_csv(args.path_to_uberon_cl_po_csv), _csv(args.path_to_doid_csv),
                                                _csv(args.path_to_exposure_csv))
    vocab = Vocab(allowed_terms, uberon_table, doid_table, exposure_table,
                  names_dmp=args.path_ncbi_names if args.path_ncbi_names and os.path.exists(args.path_ncbi_names) else None)
    stats = collections.Counter()

    tables = {}
    # microbeMASST
    df_mm = read_table(args.path_microbeMASST, ',', keep_default_na=True)
    if df_mm is not None:
        df_mm = df_mm.dropna(subset=['Filepath', 'MassIVE']).copy()
        df_mm['filename'] = df_mm['Filepath'].apply(lambda x: '/'.join(x.split('/')[1:]))
        df_mm['NCBITaxonomy'] = df_mm['Taxa_NCBI'] + '|' + df_mm['Taxaname_file']
        df_mm = df_mm.rename(columns={'MassIVE': 'MassiveID'})
        types = [microbe_sample_type(str(t).strip(), division_by_taxid, parents) for t in df_mm['Taxa_NCBI']]
        df_mm['SampleType'] = [t[0] for t in types]
        df_mm['SampleTypeSub1'] = [t[1] for t in types]
        df_mm = apply_blank_qc_flags(df_mm)
        tables['microbe'] = _collapse_repeats(df_mm[['MassiveID', 'filename', 'SampleType', 'SampleTypeSub1', 'NCBITaxonomy']], 'microbeMASST')
    # plantMASST
    df_pm = read_table(args.path_plantMASST, ',', keep_default_na=True)
    if df_pm is not None:
        df_pm = df_pm.dropna(subset=['Filepath', 'MassIVE']).copy()
        df_pm['filename'] = df_pm['Filepath'].apply(lambda x: '/'.join(x.split('/')[1:]))
        df_pm['NCBITaxonomy'] = df_pm['Taxa_NCBI'] + '|' + df_pm['Taxaname_file']
        df_pm = df_pm.rename(columns={'MassIVE': 'MassiveID'})
        df_pm['SampleType'] = 'plant'
        df_pm['SampleTypeSub1'] = ''
        df_pm = apply_blank_qc_flags(df_pm)
        tables['plant'] = _collapse_repeats(df_pm[['MassiveID', 'filename', 'SampleType', 'SampleTypeSub1', 'NCBITaxonomy']], 'plantMASST')
    # foodMASST
    df_food_csv = read_table(args.path_foodMASST, ',')
    df_food_full = read_table(args.path_foodMASST_full, '\t')
    if df_food_csv is not None and df_food_full is not None:
        tables['food'] = food_rows(df_food_csv, df_food_full, vocab, load_gfop(), division_by_taxid)
    elif (df_food_csv is None) != (df_food_full is None):
        print("WARNING: foodMASST needs both food_masst_table.csv and food_masst_table_full.tsv; skipping it")
    # personalCareProductMASST
    tables['personal_care'] = personal_care_rows(read_table(args.path_personalCareProductMASST, ','))
    # microbiomeMASST
    tables['microbiome'] = microbiome_rows(
        read_table(args.path_microbiomeMASST, '\t'), vocab, division_by_taxid,
        microbe_sample_type=lambda t: microbe_sample_type(t, division_by_taxid, parents), stats=stats)

    # foodMASST names the plant of a food by its botanical name; where that name is not an
    # NCBI name (typos such as 'rosmarinus'), the plantMASST taxon of the same file is used
    if tables.get('food') is not None and len(tables['food']) and tables.get('plant') is not None:
        plant_taxon = {(m, _stem(f)): t for m, f, t in zip(tables['plant']['MassiveID'], tables['plant']['filename'],
                                                           tables['plant']['NCBITaxonomy'].fillna(''))
                       if '|' in str(t) and str(t).split('|')[0].strip()}
        food = tables['food']
        keys = list(zip(food['MassiveID'], food['filename'].map(_stem)))
        fallback = [('|' not in str(t)) and k in plant_taxon and g in ('plant', 'fungi', 'algae')
                    for t, k, g in zip(food['NCBITaxonomy'], keys, food['_food_group'])]
        food.loc[fallback, 'NCBITaxonomy'] = [plant_taxon[k] for k, f in zip(keys, fallback) if f]
        print(f"foodMASST: {sum(fallback)} files without an NCBI name take the plantMASST taxon")

    # files of GNPS/ReDU metadata are left out; a file in several tables is taken from the
    # most specific one
    basename_only = {'food', 'personal_care'}
    taken_paths, taken_stems = set(), set()
    df_list = []
    for name in TABLE_PRECEDENCE:
        df = tables.get(name)
        if df is None or len(df) == 0:
            continue
        df = df.copy()
        df['filename'] = df['filename'].astype(str)
        stems = list(zip(df['MassiveID'], df['filename'].map(_stem)))
        paths = list(zip(df['MassiveID'], df['filename']))
        if name in basename_only:
            in_redu = [s in redu_stems for s in stems]
            shadowed = [s in taken_stems for s in stems]
        else:
            in_redu = [p in redu_files for p in paths]
            shadowed = [p in taken_paths or s in taken_stems for p, s in zip(paths, stems)]
        keep = [not a and not b for a, b in zip(in_redu, shadowed)]
        print(f"{name}MASST: {len(df)} files, {sum(in_redu)} described by GNPS metadata, "
              f"{sum(b and not a for a, b in zip(in_redu, shadowed))} taken from a more specific MASST table, "
              f"{sum(keep)} kept")
        # every listed file (also the GNPS-described ones) shadows the less specific tables
        if name in basename_only:
            taken_stems.update(stems)
        else:
            taken_paths.update(paths)
        df = df[keep]
        df['_masst_table'] = name
        df_list.append(df)

    for key, value in sorted(stats.items()):
        print(f"  {key}: {value}")

    if len(df_list) == 0:
        print("No MASST rows")
        return

    df_massts = pd.concat(df_list, ignore_index=True).fillna('').drop(columns=['_food_group'], errors='ignore')
    print('MASST rows by table:', df_massts['_masst_table'].value_counts().to_dict())
    df_massts = df_massts.drop(columns=['_masst_table'])

    ncbi_name_resolver = None
    if args.path_ncbi_names and os.path.exists(args.path_ncbi_names):
        names = {str(v).split('|', 1)[-1].strip() for v in df_massts['NCBITaxonomy'].unique() if str(v).strip()}
        ncbi_name_resolver = build_ncbi_name_resolver(args.path_ncbi_names, names,
                                                      allowed_terms['NCBITaxonomy']['allowed_values'], NCBIRankDivision_table)
    kwargs = {}
    if uberon_table is not None:
        kwargs['UBERONOntologyIndex_table'] = uberon_table
    if doid_table is not None:
        kwargs['DOIDOntologyIndex_table'] = doid_table
    if exposure_table is not None:
        kwargs['ExposureIndex_table'] = exposure_table
    df_massts_filled = complete_and_fill_REDU_table(df_massts, allowedTerm_dict=allowed_terms,
                                                    NCBIRankDivision_table=NCBIRankDivision_table,
                                                    NCBI_name_resolver=ncbi_name_resolver, **kwargs)

    # save output, one file per dataset
    for massive_id in df_massts_filled['MassiveID'].unique():
        df_filtered = df_massts_filled[df_massts_filled['MassiveID'] == massive_id]
        df_filtered.to_csv(f"{args.output_folder}/{massive_id}.tsv", sep='\t', index=False)


if __name__ == '__main__':
    main()
