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
described by GNPS/ReDU metadata files are left out here (also when the metadata describes
another copy of the file, e.g. ccms_peak/A/x.mzML for MASST's peak/A/x.mzXML: same_file_keys
in merge_metadata.py); their MASST rows only fill missing fields (--supplement_output), and
GNPS rows win the merge (merge_metadata.py) for the same file.
NCBI taxids that NCBI has merged into another taxon (merged.dmp, --path_ncbi_merged) give the
sample type by the current taxid, and replace the taxid where the name alone does not resolve.
"""
import argparse
import collections
import json
import os

import pandas as pd

from read_and_validate_redu_from_github import complete_and_fill_REDU_table, build_ncbi_name_resolver
from masst_tables import (read_table, food_rows, personal_care_rows, microbiome_rows, Vocab, load_gfop)
from merge_metadata import link_same_files, NameCopies, name_copy_datasets

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


def load_ncbi_merged(path_to_merged_dmp):
    """old taxid -> current taxid from NCBI merged.dmp ({} without the file)"""
    merged = {}
    if not path_to_merged_dmp or path_to_merged_dmp == 'none' or not os.path.exists(path_to_merged_dmp):
        return merged
    with open(path_to_merged_dmp) as f:
        for line in f:
            parts = line.split('\t|\t')
            if len(parts) >= 2:
                merged[parts[0].strip()] = parts[1].strip().rstrip('\t|').strip()
    return merged


def current_taxid(taxid, merged):
    """the current NCBI taxid of a taxid NCBI merged into another taxon (merged.dmp), else the taxid"""
    t = str(taxid).strip()
    return merged.get(t, t) if merged else t


def follow_merged_taxa(values, merged, name_resolves, stats=None):
    """'old|name' -> 'current|name' for NCBI taxids that NCBI merged into another taxon, where the
    name does not identify an allowed taxon on its own (the validator then finds the allowed term by
    the current id; a value whose name resolves is left to the validator as before)"""
    if not merged:
        return values
    def one(value):
        text = '' if value is None or (isinstance(value, float) and value != value) else str(value)
        tid, sep, name = text.partition('|')
        current = merged.get(tid.strip())
        if not current or not sep:
            return value
        if name_resolves(name):
            if stats is not None:
                stats['ncbi_merged_taxid_name_resolves_(left_as_is)'] += 1
            return value
        if stats is not None:
            stats['ncbi_merged_taxid_followed'] += 1
        return current + sep + name
    return [one(v) for v in values]


def ncbi_name_resolves(allowed_ncbi_terms, ncbi_name_resolver):
    """Whether the validator resolves a taxon name without its id (an allowed scientific name of one
    taxon, 'Genus sp.' of an allowed genus, or a names.dmp synonym / strain of an allowed taxon; not
    'Genus sp. <a lab's numbering>', which names the genus only without its id:
    REDU_conversion_functions.sp_generic_label)."""
    from REDU_conversion_functions import sp_generic_label
    by_name = collections.Counter(t.split('|', 1)[1].strip().lower() for t in allowed_ncbi_terms if '|' in t)
    def resolves(name):
        name = str(name).strip()
        key = name.lower()
        words = key.split()
        if sp_generic_label(name):
            return False
        if len(words) == 2 and words[1] in ('sp.', 'spp.', 'sp', 'spp') and by_name.get(words[0]) == 1:
            return True
        if by_name.get(key) == 1:
            return True
        return bool(ncbi_name_resolver is not None and ncbi_name_resolver(name))
    return resolves


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
    """Blank/QC files: label them as blank/QC. Blanks carry no organism (the taxonomy is dropped);
    a pooled QC keeps the organism of its row (a pool of the study's samples)"""
    blank = df['Blank'].astype(str).str.strip().str.lower().eq('yes')
    qc = df['QC'].astype(str).str.strip().str.lower().eq('yes') & ~blank
    df.loc[blank, 'SampleType'] = 'blank_analysis'
    df.loc[blank, 'SampleTypeSub1'] = 'blank_analysis'
    df.loc[qc, 'SampleType'] = 'pool_QC'
    df.loc[qc, 'SampleTypeSub1'] = 'pool_QC'
    df.loc[blank, 'NCBITaxonomy'] = ''
    print(f"Blank/QC files: {blank.sum()} blanks, {qc.sum()} QCs")
    return df


def write_refmat_evidence(df_massts, path):
    """The reference-material evidence (reference_materials.EVIDENCE_COLUMN) of the MASST rows of
    GNPS-described files (USI of the GNPS row), for final_table_rules.reference_materials_pass: the merge
    keeps the GNPS row of such a file, which may call a reference material a sample."""
    from reference_materials import EVIDENCE_COLUMN
    out = pd.DataFrame(columns=['USI', 'evidence'])
    if EVIDENCE_COLUMN in df_massts.columns and 'USI' in df_massts.columns:
        usi = df_massts['USI'].fillna('').astype(str)
        ev = df_massts[EVIDENCE_COLUMN].fillna('').astype(str)
        keep = usi.str.startswith('mzspec:') & (ev != '')
        out = pd.DataFrame({'USI': usi[keep], 'evidence': ev[keep]}).drop_duplicates('USI').sort_values('USI')
    out.to_csv(path, sep='\t', index=False)
    print(f"Reference-material evidence for GNPS-described files: {len(out)} USIs")


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
    parser.add_argument('--path_ncbi_merged', default=None,
                        help='NCBI merged.dmp: taxids NCBI merged into another taxon are replaced by the current one')
    parser.add_argument('--files_table', default=None,
                        help='dataset cache file table (uniquemri.csv): recognises copies of a GNPS-described file in '
                             'another folder layout (merge_metadata.NameCopies); without it only the same path or a copy '
                             'in another top folder counts')
    parser.add_argument('--path_to_uberon_cl_po_csv', default=None)
    parser.add_argument('--path_to_doid_csv', default=None)
    parser.add_argument('--path_to_exposure_csv', default=None)
    parser.add_argument("--AllowedTermJson_path", type=str, help="Path to json with allowed terms")
    parser.add_argument('--supplement_output', default=None,
                        help='TSV for the MASST rows of files that GNPS metadata already describes (with the '
                             'USI of the GNPS row); merge_metadata.py fills only their missing fields from it')
    parser.add_argument('--refmat_evidence_output', default=None,
                        help='TSV (USI, evidence) of the GNPS-described files whose MASST rows name a reference '
                             'material (e.g. microbiomeMASST "human reference material"); final_table_rules.py reads '
                             'it for the GNPS rows of these files')
    args = parser.parse_args()
    if args.supplement_output:
        open(args.supplement_output, 'w').close()
    if args.refmat_evidence_output:
        pd.DataFrame(columns=['USI', 'evidence']).to_csv(args.refmat_evidence_output, sep='\t', index=False)

    df_redu = pd.read_csv(args.path_REDUtable, sep='\t', dtype=str)
    df_redu['filename'] = df_redu['filename'].str.slice(2)
    # Files already described by GNPS metadata, identified by dataset + path (the same
    # relative path, e.g. ccms_peak/Blank.mzML, exists in many datasets) or as another copy
    # of the same file (same_file_keys); tables that only give the file name are compared by
    # dataset + name
    redu_stems = set(zip(df_redu['ATTRIBUTE_DatasetAccession'], df_redu['filename'].map(_stem)))
    redu_usis = ('mzspec:' + df_redu['ATTRIBUTE_DatasetAccession'].fillna('') + ':' + df_redu['filename'].fillna(''))
    # USI of each GNPS-described file, for the supplement rows (a name that fits several
    # files of the dataset, e.g. pos/ and neg/ copies, is not used)
    usi_by_redu_usi, usis_by_stem = {}, collections.defaultdict(set)
    if 'USI' in df_redu.columns:
        for m, f, u, ru in zip(df_redu['ATTRIBUTE_DatasetAccession'], df_redu['filename'], df_redu['USI'], redu_usis):
            usi_by_redu_usi[ru] = u
            usis_by_stem[(m, _stem(f))].add(u)

    with open(args.AllowedTermJson_path, 'r') as json_file:
        allowed_terms = json.load(json_file)

    NCBIRankDivision_table = pd.read_csv(args.path_ncbiRanksDivisions, index_col=False)
    NCBIRankDivision_table = NCBIRankDivision_table.drop_duplicates(subset=['TaxonID'])
    division_by_taxid = dict(zip(NCBIRankDivision_table['TaxonID'].astype(str), NCBIRankDivision_table['NCBIDivision']))
    parents = load_ncbi_parents(args.path_ncbi_nodes) if args.path_ncbi_nodes else None
    merged_taxa = load_ncbi_merged(args.path_ncbi_merged)

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
        # a taxid NCBI merged into another taxon is looked up by its current taxid
        types = [microbe_sample_type(current_taxid(t, merged_taxa), division_by_taxid, parents) for t in df_mm['Taxa_NCBI']]
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
        microbe_sample_type=lambda t: microbe_sample_type(current_taxid(t, merged_taxa), division_by_taxid, parents), stats=stats)

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
    name_copies = None
    if args.files_table:
        masst_usis = pd.Series(['mzspec:' + str(m) + ':' + str(f) for name in TABLE_PRECEDENCE if name not in basename_only
                                and tables.get(name) is not None and len(tables[name])
                                for m, f in zip(tables[name]['MassiveID'], tables[name]['filename'])], dtype=object)
        datasets = name_copy_datasets(masst_usis, redu_usis, link_same_files(masst_usis, redu_usis))
        name_copies = NameCopies(args.files_table, list(redu_usis) + list(masst_usis), datasets)
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
            # the GNPS-described copies of each file (the same path, the same path with another
            # extension, or a copy in another top folder)
            query = pd.Series(['mzspec:' + str(m) + ':' + str(f) for m, f in paths])
            same_files = link_same_files(query, redu_usis)
            if name_copies is not None:
                same_files = name_copies.link(query, redu_usis, same_files)
            in_redu = [bool(t) for t in same_files]
            shadowed = [p in taken_paths or s in taken_stems for p, s in zip(paths, stems)]
        keep = [not a and not b for a, b in zip(in_redu, shadowed)]
        # GNPS-described files of the most specific table listing them: supplement rows
        if name in basename_only:
            target = [next(iter(usis_by_stem[s])) if a and not b and len(usis_by_stem.get(s, ())) == 1 else ''
                      for s, a, b in zip(stems, in_redu, shadowed)]
        else:
            # only a file with one described copy fills that copy's missing fields
            target = [usi_by_redu_usi.get(t[0], '') if a and not b and len(t) == 1 else ''
                      for t, a, b in zip(same_files, in_redu, shadowed)]
        supplement = [bool(t) for t in target]
        print(f"{name}MASST: {len(df)} files, {sum(in_redu)} described by GNPS metadata "
              f"({sum(supplement)} as supplement rows), "
              f"{sum(b and not a for a, b in zip(in_redu, shadowed))} taken from a more specific MASST table, "
              f"{sum(keep)} kept")
        # every listed file (also the GNPS-described ones) shadows the less specific tables
        if name in basename_only:
            taken_stems.update(stems)
        else:
            taken_paths.update(paths)
        df['USI'] = target
        df = df[[k or u for k, u in zip(keep, supplement)]]
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
    # taxids NCBI merged into another taxon (e.g. 137779 -> 184924): the current taxid, where the
    # name alone does not give an allowed taxon (otherwise the row lost its taxonomy or was dropped)
    if merged_taxa:
        resolves = ncbi_name_resolves(allowed_terms['NCBITaxonomy']['allowed_values'], ncbi_name_resolver)
        df_massts['NCBITaxonomy'] = follow_merged_taxa(df_massts['NCBITaxonomy'], merged_taxa, resolves, stats)
        print(f"NCBI taxids merged into another taxon: {stats.get('ncbi_merged_taxid_followed', 0)} rows take the current "
              f"taxid, {stats.get('ncbi_merged_taxid_name_resolves_(left_as_is)', 0)} rows keep it (their name resolves)")
    if args.refmat_evidence_output:
        write_refmat_evidence(df_massts, args.refmat_evidence_output)
    kwargs = {}
    if uberon_table is not None:
        kwargs['UBERONOntologyIndex_table'] = uberon_table
    if doid_table is not None:
        kwargs['DOIDOntologyIndex_table'] = doid_table
    if exposure_table is not None:
        kwargs['ExposureIndex_table'] = exposure_table
    df_massts_filled = complete_and_fill_REDU_table(df_massts, allowedTerm_dict=allowed_terms,
                                                    NCBIRankDivision_table=NCBIRankDivision_table,
                                                    NCBI_name_resolver=ncbi_name_resolver, keep_usi=True, **kwargs)
    if len(df_massts_filled) == 0:
        print("No MASST rows left after validation")
        return
    is_supplement = df_massts_filled['USI'].fillna('').astype(str).str.startswith('mzspec:')
    if args.supplement_output:
        df_massts_filled[is_supplement].to_csv(args.supplement_output, sep='\t', index=False)
        print(f"Supplement rows for GNPS-described files: {int(is_supplement.sum())}")
    df_massts_filled = df_massts_filled[~is_supplement].drop(columns=['USI'])

    # save output, one file per dataset
    for massive_id in df_massts_filled['MassiveID'].unique():
        df_filtered = df_massts_filled[df_massts_filled['MassiveID'] == massive_id]
        df_filtered.to_csv(f"{args.output_folder}/{massive_id}.tsv", sep='\t', index=False)


if __name__ == '__main__':
    main()
