import argparse
import collections
import math
import numpy as np
import pandas as pd
import glob
import re
import os
import json
import sys
import traceback
import weakref
from io import StringIO
from REDU_conversion_functions import age_category, normalize_redu_datetime
from REDU_conversion_functions import pinned_genus_term, sp_generic_label, sp_other_strain
from REDU_conversion_functions import bare_sp_id_term, _GENUS_LEVEL_SP, NCBI_COLLECTING_TAXON
from blank_qc_names import apply_filename_blank_qc, blank_qc_type_from_filename
import reference_materials
# several organisms in one NCBITaxonomy value ('ID|Name||ID|Name')
from ncbi_database import MULTI_SEP, aligned, is_multi, join_parts, submitted_taxonomy_parts, taxon_ids
# the organism an isolate was isolated from (a value like NCBITaxonomy's), without the isolate's own organisms
from ncbi_database import ISOLATION_SOURCE, without_taxa_of_series

def _canonical_term(value):
    """Case-, space- and separator-insensitive form used to match submitted values to
    allowed terms ('Male' -> 'male', 'methanol:water (7:3)' -> 'methanol-water (7:3)').
    Ratios are reduced, so 'methanol-water (8:2)' matches 'methanol-water (4:1)'."""
    s = str(value).strip().lower()
    head, sep, tail = s.partition('(')
    head = re.sub(r'\s*[:/]\s*', '-', head)
    ratio = re.fullmatch(r'\s*(\d+(?::\d+)+)\s*\)(.*)', tail) if sep else None
    if ratio:
        parts = [int(x) for x in ratio.group(1).split(':')]
        g = 0
        for x in parts:
            g = math.gcd(g, x)
        if g > 1:
            parts = [x // g for x in parts]
        tail = ':'.join(str(x) for x in parts) + ')' + ratio.group(2)
    return re.sub(r'\s+', '', head + sep + tail)


# Lookups derived from the (large, unchanging) allowed-term lists, built once per run
# instead of once per metadata file. An entry keeps its source, so the id() it is keyed on
# cannot be taken by another object while the entry exists; at most _LOOKUP_CACHE_MAX
# entries are kept.
_LOOKUP_CACHE = {}
_LOOKUP_CACHE_MAX = 512


def _cached(name, source, build):
    key = (name, id(source), len(source))
    hit = _LOOKUP_CACHE.get(key)
    if hit is None or hit[0] is not source:
        # entries whose source only this cache still holds (one reference per entry of it)
        # can never be hit again: drop them, so a caller that passes a new allowed-term list
        # per call (e.g. the NCBI list, whose lookups take ~1.2 GB) does not pin a copy per call
        entries = collections.Counter(id(v[0]) for v in _LOOKUP_CACHE.values() if isinstance(v, tuple))
        for k in [k for k, v in _LOOKUP_CACHE.items()
                  if isinstance(v, tuple) and sys.getrefcount(v[0]) <= 1 + entries[id(v[0])]]:
            del _LOOKUP_CACHE[k]
        if len(_LOOKUP_CACHE) >= _LOOKUP_CACHE_MAX:
            _LOOKUP_CACHE.clear()
        hit = (source, build())
        _LOOKUP_CACHE[key] = hit
    return hit[1]


# Lookups derived from a table argument: one entry per name, for the table object last
# passed. The entry holds only a weak reference to the table, so tables that callers pass
# as a new copy per call (Workbench, MetaboLights) are rebuilt, are not kept alive and do
# not make the cache grow.
_TABLE_CACHE = {}


def _cached_for_table(name, table, build):
    hit = _TABLE_CACHE.get(name)
    if hit is not None and hit[0]() is table and hit[1] == len(table):
        return hit[2]
    value = build()
    try:
        _TABLE_CACHE[name] = (weakref.ref(table), len(table), value)
    except TypeError:
        _TABLE_CACHE.pop(name, None)  # not weak-referenceable: not cached
    return value


class PinnedNcbiRelations:
    """Lineage relations of the pinned NCBI taxonomy (REDU_conversion_functions.use_pinned_ncbi_taxonomy):
    is_ancestor(a, t) whether taxon ID a is a proper ancestor of t; microbe_group(t) 'Bacteria', 'Archaea',
    'Fungi' or 'Viruses' when t is that group or below it, else None (merged IDs followed; the same rules as
    Metabolights2REDU._ncbi_is_ancestor / _ncbi_microbe_group)."""
    MICROBE_GROUPS = {2: 'Bacteria', 2157: 'Archaea', 4751: 'Fungi', 10239: 'Viruses'}

    def __init__(self, pinned):
        self.taxonomy = pinned.taxonomy

    def _node(self, taxid):
        try:
            return self.taxonomy.resolve(int(str(taxid).strip()))
        except (TypeError, ValueError):
            return None

    def is_ancestor(self, ancestor, taxon):
        a, node = self._node(ancestor), self._node(taxon)
        if a is None or node is None or a == node:
            return False
        parent = self.taxonomy.parent
        for _ in range(1000):
            up = parent[node]
            if up == node or up < 0:
                return False
            if up == a:
                return True
            node = up
        return False

    def microbe_group(self, taxon):
        node = self._node(taxon)
        if node is None:
            return None
        parent = self.taxonomy.parent
        for _ in range(1000):
            if node in self.MICROBE_GROUPS:
                return self.MICROBE_GROUPS[node]
            up = parent[node]
            if up == node or up < 0:
                return None
            node = up
        return None


def pinned_ncbi_relations():
    """PinnedNcbiRelations of the pinned NCBI taxonomy, or None when the caller did not load one."""
    import REDU_conversion_functions
    pinned = REDU_conversion_functions._PINNED_LINEAGE
    return PinnedNcbiRelations(pinned) if pinned is not None else None


_NCBI_COLUMNS = ('TaxonID', 'NCBIRank', 'NCBIDivision')
# The NCBI rank/division table last prepared by _prepare_ncbi_table (one entry)
_NCBI_PREPARED = {}


def _column_buffers(table):
    # the arrays of the three columns and the address of their data
    arrays = [table[c].to_numpy() for c in _NCBI_COLUMNS]
    return arrays, [a.__array_interface__['data'][0] for a in arrays]


def _forget_ncbi_table(ref, prepared=_NCBI_PREPARED):
    # weakref callback: the table is gone, so is its entry (the dict is bound here, as
    # module globals may already be cleared when this runs at interpreter exit)
    hit = prepared.get('entry')
    if hit is not None and hit['table'] is ref:
        prepared.clear()


def _prepare_ncbi_table(table):
    """Converts TaxonID, NCBIRank and NCBIDivision of the caller's NCBI table to str in
    place, as complete_and_fill_REDU_table always did, and returns an Index of its TaxonIDs
    (None for an empty table) to take only the rows of the taxa in a metadata table from it.
    The ~3M-row table is converted and indexed once per table object: the entry holds a
    weak reference to the table (it is not kept alive) and the converted column arrays (so
    their memory cannot be reused), and is used again only for the same table object with
    the same length and the very same column arrays; a call with another table, or after a
    column of this one was replaced, converts and indexes again. Values edited in place
    inside the table after the first call (no caller does that) are not converted again
    (a NaN written into NCBIRank stays NaN and gives 'missing value', where the str
    conversion on every call gave 'nan'), and TaxonIDs edited in place are not noticed: a
    row given the TaxonID of another row is not found under it."""
    hit = _NCBI_PREPARED.get('entry')
    if hit is not None and hit['table']() is table and hit['len'] == len(table) \
            and _column_buffers(table)[1] == hit['addresses']:
        return hit['index']
    for c in _NCBI_COLUMNS:
        table[c] = table[c].astype(str)
    if not len(table):
        return None
    arrays, addresses = _column_buffers(table)
    # the Index owns a copy of the TaxonIDs: its hash table only borrows references to
    # the values, so it must not share the table's array (a value replaced in place there
    # would be freed under it)
    index = pd.Index(arrays[0], dtype=object, copy=True)
    try:
        ref = weakref.ref(table, _forget_ncbi_table)
    except TypeError:
        return index
    _NCBI_PREPARED['entry'] = {'table': ref, 'len': len(table), 'arrays': arrays,
                               'addresses': addresses, 'index': index}
    return index


def _ncbi_rows_of(table, taxon_index, taxon_ids):
    """The rows of the NCBI table whose TaxonID is one of taxon_ids, in table order: a
    left merge with them gives the same rows, values and columns as with the whole table."""
    needed = pd.unique(taxon_ids.to_numpy(dtype=object))
    if taxon_index.is_unique:
        rows = taxon_index.get_indexer(needed)
    else:
        rows, _ = taxon_index.get_indexer_non_unique(needed)
    return table.iloc[np.sort(rows[rows >= 0])]


# an organism name at genus level: 'Genus sp.', 'Genus spp.', also with NCBI's note on the genus name ('Bacillus sp.
# (in: firmicutes)' for a name of two genera, 'Saccharimonas sp. (SeqCode)') (complete_and_fill_REDU_table,
# NCBITaxonomy): REDU_conversion_functions._GENUS_LEVEL_SP, imported above.
# NCBI's taxa that collect the unnamed members of a taxon ('unclassified Sulfitobacter', 'Pseudomonas environmental
# samples') and its placeholders for organisms known by a sequence only ('uncultured Streptomyces sp.', 'unidentified
# Pseudomonas'): no one organism either (REDU_conversion_functions.NCBI_COLLECTING_TAXON)
_NCBI_COLLECTING_TAXON = NCBI_COLLECTING_TAXON

_DOID_ITEM = re.compile(r'(?i)^DOID[:_]?(\d+)$')


def normalize_doid_list(value, doid_ids, missing_value='missing value'):
    """'DOID:526|DOID:1596' style list of DOID IDs. Items may be separated by | , ; or
    whitespace and written DOID_526 / doid:526; each is written DOID:<digits>. IDs that are
    not in doid_ids (unknown or obsolete) are dropped, as are repeats (the submitter's order
    is kept). Returns missing_value when no valid ID remains."""
    text = str(value).strip()
    if not text:
        return missing_value
    text = re.sub(r'(?i)\bDOID\s*([:_])\s*', r'DOID\1', text)
    ids = []
    for item in re.split(r'[|,;\s]+', text):
        m = _DOID_ITEM.match(item.strip())
        if not m:
            continue
        doid = 'DOID:' + m.group(1)
        if doid in doid_ids and doid not in ids:
            ids.append(doid)
    return '|'.join(ids) if ids else missing_value


# Rows that are no samples of the study's subjects (blanks, pooled QCs, reference materials, standards,
# products) lose the subject fields (STRIP_FOR_NOT_SUBJECT). Pooled QCs (a pool of the study's samples)
# and reference materials keep their organism and body part (STRIP_ORGANISM): a pooled QC its own
# row's (not one whose file is named as a blank), a reference material the one of the material
# (reference_materials.py, rows it placed), every other reference material none. Blanks, standards and
# products keep none.
STRIP_FOR_NOT_SUBJECT = ['AgeInYears', 'UniqueSubjectID', 'SubjectIdentifierAsRecorded', 'DOIDCommonName', 'HPOCommonName',
                         'Exposure', 'TimePointHoursRelativeToEvent', 'FastingStatus', 'ENVOEnvironmentBiome',
                         'ENVOEnvironmentMaterial', 'BiologicalSex', 'HealthStatus', 'LatitudeandLongitude', 'SmokingStatus',
                         'SmokingPackYears', 'BodyMassIndex', 'Diet', 'CellLine', 'CellLineIndex', 'ColonizationStatus']
STRIP_ORGANISM = ['NCBITaxonomy', 'UBERONBodyPartName']
# the full list as it was (blanks): both
strip_for_blank = STRIP_ORGANISM + STRIP_FOR_NOT_SUBJECT
# where a cultured organism was isolated from: kept only with the organism of a pooled QC (a reference material has
# the material's organism)
STRIP_ISOLATION_SOURCE = [ISOLATION_SOURCE, 'IsolationSourceUBERONBodyPartName', 'IsolationSourceENVOEnvironmentMaterial']


def columns_with_terms_of(allowedTerm_dict, key):
    """key and the columns that have its allowed terms (allowed_values_from in the allowed terms json: the
    IsolationSource columns have those of NCBITaxonomy, UBERONBodyPartName and ENVOEnvironmentMaterial)"""
    return [key] + [k for k, v in allowedTerm_dict.items() if v.get('allowed_values_from') == key]


def _not_subject_rows(df):
    """(rows that are no samples of the study's subjects, those of them that keep NCBITaxonomy and
    UBERONBodyPartName: pooled QCs and the reference materials reference_materials.py placed)"""
    sample_type = df['SampleType'].fillna('').astype(str)
    mask = sample_type.str.contains('blank', case=False) | \
        sample_type.isin(['pool_QC', 'reference material', 'personal care product', 'molecular_standards'])
    placed = df[reference_materials.PLACED_COLUMN].astype(str).str.lower().eq('true') \
        if reference_materials.PLACED_COLUMN in df.columns else pd.Series(False, index=df.index)
    # (a pooled QC whose file is named as a blank, e.g. a blank file listed for a QC sample, keeps none)
    blank_file = df['filename'].map(lambda f: (blank_qc_type_from_filename(f) or '').startswith('blank')) \
        if 'filename' in df.columns and (sample_type == 'pool_QC').any() else pd.Series(False, index=df.index)
    keep = ((sample_type == 'pool_QC') & ~blank_file) | ((sample_type == 'reference material') & placed)
    return mask, keep & mask


def _reference_materials_and_pools(df, kwargs):
    """reference_materials.annotate on the rows (dataset column MassiveID): the rows' ReDU text columns
    plus the source's native sample fields (column reference_materials.EVIDENCE_COLUMN, if the converter
    adds it) and the catalog numbers its study-level text names (kwargs['refmat_study_text']: a list,
    or {dataset: list}). A problem here never costs the rows."""
    try:
        stats = collections.Counter()
        df, changed = reference_materials.annotate(
            df, dataset_column='MassiveID', evidence=reference_materials.EVIDENCE_COLUMN,
            study_text=kwargs.get('refmat_study_text'), stats=stats)
        df = reference_materials.place_labelled(df, stats=stats)
        if stats:
            # one line per call: what was decided, applied and kept (guards), by dataset
            datasets = sorted(set(df['MassiveID'].astype(str))) if 'MassiveID' in df.columns else []
            print(f"REFMAT {len(changed)} rows changed; datasets {','.join(datasets[:5])}{' ...' if len(datasets) > 5 else ''}; "
                  f"{json.dumps(dict(sorted(stats.items())))}")
    except Exception as e:
        print(f"WARNING: reference materials not checked: {type(e).__name__}: {e}")
        traceback.print_exc()
    return df


def complete_and_fill_REDU_table(df, allowedTerm_dict, add_usi = False, keep_usi = False, other_allowed_file_extensions = [], attempt_adding_file_extensions = False, case_insensitive_terms = False, **kwargs):
    """
    Completes and fills a REDU table with values based on a dictionary of allowed terms and missing values.

    Args:
    df: A pandas DataFrame containing the initial data.
    allowedTerm_dict: A dictionary containing allowed terms and missing values for each column.

    Returns:
    A DataFrame that has been filled with default values for missing columns,
    with values replaced by the corresponding "missing" value from the dictionary
    if they are not in the allowed terms or are missing/empty, except for specific columns.

    The ontology and NCBI tables passed in are changed in place: their 'Label' columns
    (and TaxonID, NCBIRank, NCBIDivision of the NCBI table) become str, so NaN becomes
    'nan', and a DOID table's 'UBERONOntologyIndex' column is renamed 'DOIDOntologyIndex'.
    Callers read these tables after the call (e.g. Metabolights2REDU matches ENVO labels
    of the converted table in the next study), so this is kept as it is.
    """

    if len(df) == 0:
        # On an empty (0-row) input, pandas' .map({}).fillna(...) chain below
        # produces float64-typed empty columns (no data to infer a dtype from),
        # which then fails when merged against the ontology tables' object-typed
        # 'Label' columns ("You are trying to merge on float64 and object
        # columns"). There's nothing to process for an empty table anyway, so
        # short-circuit before any of that logic runs.
        return df

    #prepare UberonOntologyIndex
    if 'UBERONOntologyIndex_table' in kwargs.keys():
        uberon_ontology_table = kwargs['UBERONOntologyIndex_table']
    else:
        uberon_ontology_table = pd.DataFrame(columns=['Label', 'UBERONOntologyIndex', 'Is Multicellular','Is Organ', 'Is Fluid'])

    uberon_ontology_table['Label'] = uberon_ontology_table['Label'].astype(str)

    #prepare DOIDOntologyTable
    if 'DOIDOntologyIndex_table' in kwargs.keys():
        doid_ontology_table = kwargs['DOIDOntologyIndex_table']
        if 'UBERONOntologyIndex' in doid_ontology_table.columns:
            doid_ontology_table.rename(columns={'UBERONOntologyIndex': 'DOIDOntologyIndex'}, inplace=True)
    else:
        doid_ontology_table = pd.DataFrame(columns=['Label', 'DOIDOntologyIndex'])
    
    doid_ontology_table['Label'] = doid_ontology_table['Label'].astype(str)

    #prepare HPO (health conditions) and exposure (ECTO/PECO) tables
    hpo_ontology_table = kwargs.get('HPOOntologyIndex_table')
    if hpo_ontology_table is None:
        hpo_ontology_table = pd.DataFrame(columns=['Label', 'HPOOntologyIndex'])
    hpo_ontology_table['Label'] = hpo_ontology_table['Label'].astype(str)
    exposure_ontology_table = kwargs.get('ExposureIndex_table')
    if exposure_ontology_table is None:
        exposure_ontology_table = pd.DataFrame(columns=['Label', 'ExposureIndex'])
    exposure_ontology_table['Label'] = exposure_ontology_table['Label'].astype(str)

    #prepare environment biome table
    if 'ENVOEnvironmentBiomeIndex_table' in kwargs.keys():
        envome_biome_ontology_table = kwargs['ENVOEnvironmentBiomeIndex_table']
    else:
        envome_biome_ontology_table = pd.DataFrame(columns=['Label', 'ENVOEnvironmentBiomeIndex'])

    envome_biome_ontology_table['Label'] = envome_biome_ontology_table['Label'].astype(str)

    #prepare environment material table
    if 'ENVOEnvironmentMaterialIndex_table' in kwargs.keys():
        envome_material_ontology_table = kwargs['ENVOEnvironmentMaterialIndex_table']
    else:
        envome_material_ontology_table = pd.DataFrame(columns=['Label', 'ENVOEnvironmentMaterialIndex'])

    envome_material_ontology_table['Label'] = envome_material_ontology_table['Label'].astype(str)

    #prepare environment medium scale table
    if 'ENVOEnvironmentMaterialIndex_table' in kwargs.keys():
        envome_medium_scale_ontology_table = kwargs['ENVOEnvironmentMaterialIndex_table']
        envome_medium_scale_ontology_table = envome_medium_scale_ontology_table.rename(
            columns={'ENVOEnvironmentMaterialIndex': 'ENVOMediumScaleIndex'}
            )
    else:
        envome_medium_scale_ontology_table = pd.DataFrame(columns=['Label', 'ENVOMediumScaleIndex'])

    envome_medium_scale_ontology_table['Label'] = envome_medium_scale_ontology_table['Label'].astype(str)


    #prepare environment local scale table
    if 'ENVOEnvironmentBiomeIndex_table' in kwargs.keys():
        envome_local_scale_ontology_table = kwargs['ENVOEnvironmentBiomeIndex_table']
        envome_local_scale_ontology_table = envome_local_scale_ontology_table.rename(
            columns={'ENVOEnvironmentBiomeIndex': 'ENVOLocalScaleIndex'}
            )
    else:
        envome_local_scale_ontology_table = pd.DataFrame(columns=['Label', 'ENVOLocalScaleIndex'])

    envome_local_scale_ontology_table['Label'] = envome_local_scale_ontology_table['Label'].astype(str)


    #prepare environment broad scale table
    if 'ENVOEnvironmentBiomeIndex_table' in kwargs.keys():
        envome_broad_scale_ontology_table = kwargs['ENVOEnvironmentBiomeIndex_table']
        envome_broad_scale_ontology_table = envome_broad_scale_ontology_table.rename(
            columns={'ENVOEnvironmentBiomeIndex': 'ENVOBroadScaleIndex'}
            )
    else:
        envome_broad_scale_ontology_table = pd.DataFrame(columns=['Label', 'ENVOBroadScaleIndex'])

    envome_broad_scale_ontology_table['Label'] = envome_broad_scale_ontology_table['Label'].astype(str)


    #prepare NCBI table
    if 'NCBIRankDivision_table' in kwargs.keys():
        NCBIRankDivision_table = kwargs['NCBIRankDivision_table']
    else:
        NCBIRankDivision_table = pd.DataFrame(columns=['TaxonID', 'NCBIRank', 'NCBIDivision'])

    # str in place, once per table object; index of its TaxonIDs for the NCBIRank merge
    ncbi_taxon_index = _prepare_ncbi_table(NCBIRankDivision_table)

    # resolves NCBI names that are not allowed scientific names (synonyms, merged taxa,
    # strains) to an allowed term, or returns None
    ncbi_name_resolver = kwargs.get('NCBI_name_resolver')
    # lineage relations of the pinned NCBI taxonomy for an 'ID|Name' whose two parts name different taxa
    # (name_over_id): only where the caller asks for it (the GNPS harmonize step, which loads the pinned taxonomy
    # from --path_to_ncbi_nodes_dmp); the converters' terms are allowed terms already
    ncbi_lineage = pinned_ncbi_relations() if kwargs.get('ncbi_name_over_id') else None
    # the pinned NCBI taxonomy for a bare 'Genus sp.' with the ID of a species or strain of that genus
    # (REDU_conversion_functions.bare_sp_id_term): the caller's (kwargs 'ncbi_taxonomy', a PinnedNcbiLineage:
    # MASST_to_REDU), else the one the caller loaded with use_pinned_ncbi_taxonomy (GNPS harmonize, the converters)
    import REDU_conversion_functions as _rcf_pinned
    bare_sp_taxonomy = kwargs.get('ncbi_taxonomy') or _rcf_pinned._PINNED_LINEAGE
    # names no lookup resolves: other white space, a strain designation at the end, one wrong character
    # (REDU_conversion_functions.ncbi_name_fallback); only where the caller asks for it (the GNPS harmonize
    # step, with the pinned NCBI taxonomy)
    use_ncbi_name_fallback = bool(kwargs.get('ncbi_name_fallback'))

    # Cellosaurus table (build_cellosaurus_table.py): checks CellLine, gives CellLineIndex
    cellosaurus_table = kwargs.get('CellLineIndex_table')
    cellosaurus_index = None
    if cellosaurus_table is not None and len(cellosaurus_table):
        from cellosaurus import CellosaurusIndex
        cellosaurus_index = _cached('cellosaurus', cellosaurus_table, lambda: CellosaurusIndex(cellosaurus_table))

    # Convert all columns to String (missing values become '' rather than the string 'nan')
    df = df.fillna('').astype(str)

    def _split_id_name(term: str):
        term = "" if term is None else str(term)
        if "|" in term:
            tid, tname = term.split("|", 1)
            return tid.strip(), tname.strip()
        return None, term.strip()

    def _is_sp_or_spp(name: str) -> bool:
        toks = name.strip().split()
        return len(toks) >= 2 and toks[1].rstrip(".").lower() in {"sp", "spp"}

    def _genus_from_name(name: str) -> str:
        return name.strip().split()[0] if name.strip() else ""

    def _genus_level_sp(name: str) -> bool:
        # 'Genus sp.' / 'Genus spp.' with nothing after it but NCBI's note on the genus name ('Bacillus sp.
        # (in: firmicutes)'): an answer at genus level. 'Genus sp. <strain, isolate, clone>' names one organism
        # ('Synechocystis sp. PCC 6803', 'Streptomyces sp. CNQ-509').
        return _GENUS_LEVEL_SP.fullmatch(" ".join(str(name).split())) is not None

    #remove autogenerated columns if present
    keys_to_drop = [key for key, value in allowedTerm_dict.items() if value.get("generate") == "True"]
    columns_to_remove = [col for col in keys_to_drop + ['Is Multicellular', 'Is Organ', 'Is Fluid'] if col in df.columns]
    if keep_usi:
        print("Keeping USI column")
        columns_to_remove = [c for c in columns_to_remove if c != "USI"]
    # LifeStage is regenerated from the age below; keep the submitter's value for when it can't be
    if 'LifeStage' in df.columns:
        df['_recorded_LifeStage'] = df['LifeStage']
    # the ages as given, to tell a LifeStage a converter derived from them (see LifeStage below)
    ages_as_given = df['AgeInYears'].copy() if 'AgeInYears' in df.columns else None
    df = df.drop(columns=columns_to_remove)


    # reference materials and in-house pools named in the rows' own fields (reference_materials.py)
    df = _reference_materials_and_pools(df, kwargs)

    #Resolve contradictory fields
    if 'SampleType' in df.columns:
        mask, keep_organism = _not_subject_rows(df)
        df.loc[mask, [c for c in STRIP_FOR_NOT_SUBJECT if c in df.columns]] = ''
        clear_organism = mask & ~keep_organism
        df.loc[clear_organism, [c for c in STRIP_ORGANISM if c in df.columns]] = ''
        # (the isolation source stays with the organism of a pooled QC only)
        df.loc[mask & ~(keep_organism & (df['SampleType'] == 'pool_QC')),
               [c for c in STRIP_ISOLATION_SOURCE if c in df.columns]] = ''
        if '_recorded_LifeStage' in df.columns:
            df.loc[mask, '_recorded_LifeStage'] = ''

    if 'filename' in df.columns and attempt_adding_file_extensions == True:
        df.loc[df['filename'].str.contains(r'^[^.]+$') & (df['filename'].str.len() >= 1), 'filename'] += '.mzML'
        

    input_columns = set(df.columns)

    # Add missing columns with their respective default missing value from the dictionary or generate values if generate == True
    for key, value in allowedTerm_dict.items():
        if key not in df.columns:
            if value['generate'] == 'False':
                if value['missing'] == 'not allowed':
                    print(f'{key}: MISSING. STOPPING INTERPRETATION!')
                    return pd.DataFrame()
                else:
                    df[key] = value['missing']
                    print(f'{key}: ADDED!')


    for key, value in allowedTerm_dict.items():
        if key in df.columns and value['generate'] == 'False':
            # a column with the allowed terms of another one (allowed_values_from) is checked as that one
            rule = value.get('allowed_values_from', key)
            terms_entry = allowedTerm_dict.get(rule, value)
            allowed_terms = terms_entry['allowed_values']
            missing_value = value['missing']

            # Extract unique values from the DataFrame column
            unique_values = df[key].unique()
            
            if key == 'ComorbidityListDOIDIndex':
                # a list of DOID IDs, not a fixed set of combinations: every ID must be a
                # current (non-obsolete) DOID class of the DOID table
                doid_ids = _cached_for_table('doid_ids', doid_ontology_table, lambda: set(
                    doid_ontology_table['DOIDOntologyIndex'].dropna().astype(str).str.strip())
                    if 'DOIDOntologyIndex' in doid_ontology_table.columns else set())
                if not doid_ids and any(str(v).strip() not in ('', missing_value) for v in unique_values):
                    print(f'{key}: no DOID table given; IDs cannot be checked and are set to {missing_value}')
                value_map = {v: normalize_doid_list(v, doid_ids, missing_value) for v in unique_values}
                df[key] = df[key].map(value_map).fillna(missing_value).replace("", missing_value)
                filtered_map = {k: v for k, v in value_map.items() if k != v}
                if filtered_map:
                    print(f"{key}:")
                    for k, v in filtered_map.items():
                        print(f"  {k}: {v}")
                continue

            if key == 'CellLine':
                # the Cellosaurus recommended name: the name itself, its accession, the name
                # in other case or a synonym of exactly one line; else missing. Without the
                # table only allowed names are kept.
                def _cell_line(v):
                    if str(v).strip() in ('', missing_value):
                        return missing_value
                    if cellosaurus_index is not None:
                        hit = cellosaurus_index.canonical(v)
                        return hit[0] if hit else missing_value
                    allowed_set = _cached('cellline_names', allowed_terms, lambda: set(allowed_terms))
                    return v if v in allowed_set else missing_value
                value_map = {v: _cell_line(v) for v in unique_values}
                df[key] = df[key].map(value_map).fillna(missing_value).replace("", missing_value)
                filtered_map = {k: v for k, v in value_map.items() if k != v and k != ''}
                if filtered_map:
                    print(f"{key}:")
                    for k, v in filtered_map.items():
                        print(f"  {k}: {v}")
                continue

            # Initialize value_map for each condition
            if rule not in ['NCBITaxonomy', 'MassSpectrometer', 'filename'] and len(allowed_terms) >= 1 \
                    and allowed_terms[0] not in ['00', 'datetime', 'numeric', 'numeric|numeric']:
                # General case for non-specific keys with a list of allowed terms: the first
                # allowed term that equals the value when spaces are ignored
                def _build_without_spaces(terms=allowed_terms):
                    if not all(isinstance(t, str) for t in terms):
                        return None
                    first = {}
                    for allowed in terms:
                        first.setdefault(allowed.replace(" ", ""), allowed)
                    return first
                without_spaces = _cached('without_spaces', allowed_terms, _build_without_spaces)
                if without_spaces is not None:
                    value_map = {
                        observed_value: without_spaces.get(str(observed_value).replace(" ", ""), missing_value)
                        for observed_value in unique_values
                    }
                else:
                    value_map = {
                        observed_value: next((allowed for allowed in allowed_terms if allowed.replace(" ", "") == str(observed_value).replace(" ", "")), missing_value)
                        for observed_value in unique_values
                    }
                if rule == 'UBERONBodyPartName' and isinstance(terms_entry.get('synonyms'), dict):
                    # an exact synonym of exactly one class (built from the pinned UBERON/CL/PO by
                    # update_allowed_terms_from_ontologies.py) -> its label, e.g. the former CL
                    # label 'mononuclear cell' -> 'mononuclear leukocyte' (CL:0000842)
                    synonyms = terms_entry['synonyms']
                    allowed_labels = _cached('allowed_set', allowed_terms, lambda: set(allowed_terms))
                    for observed_value, mapped in value_map.items():
                        if mapped == missing_value:
                            label = synonyms.get(str(observed_value).strip().lower())
                            if label is not None and label in allowed_labels:
                                value_map[observed_value] = label
                if case_insensitive_terms:
                    # Fall back to a case/separator-insensitive match when it is unambiguous
                    def _build_canonical(terms=allowed_terms):
                        canonical = {}
                        for allowed in terms:
                            canonical.setdefault(_canonical_term(allowed), set()).add(allowed)
                        return canonical
                    canonical = _cached('canonical', allowed_terms, _build_canonical)
                    for observed_value, mapped in value_map.items():
                        if mapped == missing_value:
                            hits = canonical.get(_canonical_term(observed_value), set())
                            if len(hits) == 1:
                                value_map[observed_value] = next(iter(hits))


            elif rule in ['NCBITaxonomy', 'MassSpectrometer']:
                if rule == 'NCBITaxonomy':
                    # Build fast lookups from allowed terms: 'ID|Name'
                    # Names are not unique in NCBI (e.g. Bacillus the bacterium and the stick
                    # insect), so a shared name only resolves when the provided ID disambiguates it
                    def _build_ncbi(terms=allowed_terms):
                        allowed_by_id = {}
                        allowed_by_name_ci = {}
                        shared_names_ci = {}
                        for t in terms:
                            if "|" in t:
                                tid, tname = t.split("|", 1)
                                tid = tid.strip()
                                tname = tname.strip().lower()
                                allowed_by_id[tid] = t
                                first = allowed_by_name_ci.setdefault(tname, t)
                                if first != t:
                                    shared_names_ci.setdefault(tname, {first}).add(t)
                        return set(terms), allowed_by_id, allowed_by_name_ci, shared_names_ci
                    allowed_set, allowed_by_id, allowed_by_name_ci, shared_names_ci = _cached('ncbi', allowed_terms, _build_ncbi)
                    division_by_id = {}

                    def division_of(taxid):
                        # NCBI division of a TaxonID (None if not in the table)
                        if ncbi_taxon_index is not None and ncbi_taxon_index.is_unique:
                            try:
                                position = ncbi_taxon_index.get_loc(taxid)
                            except KeyError:
                                position = None
                            if position is not None and NCBIRankDivision_table['TaxonID'].iat[position] == taxid:
                                return NCBIRankDivision_table['NCBIDivision'].iat[position]
                            if position is None:
                                return None
                        if not division_by_id:
                            division_by_id.update(zip(NCBIRankDivision_table['TaxonID'], NCBIRankDivision_table['NCBIDivision']))
                        return division_by_id.get(taxid)

                    def name_hit(name, oid):
                        terms = shared_names_ci.get(name.lower())
                        if terms is None:
                            return allowed_by_name_ci.get(name.lower())
                        # Shared name: use the NCBI division of the provided ID to pick one
                        # (e.g. '1409|Bacillus sp.' is in Bacteria -> 1386|Bacillus, not the insect)
                        if not oid:
                            return None
                        division = division_of(oid)
                        same_division = [t for t in terms if division and division_of(t.split("|", 1)[0].strip()) == division]
                        return same_division[0] if len(same_division) == 1 else None

                    def fallback_hit(name, oid):
                        # a name no lookup resolves, read with its white space normalised, without a strain designation
                        # at its end or with one wrong character fixed (REDU_conversion_functions.ncbi_name_fallback;
                        # caller's choice)
                        if not use_ncbi_name_fallback or not name:
                            return None
                        import REDU_conversion_functions
                        term, how = REDU_conversion_functions.ncbi_name_fallback(
                            name, allowed_by_id.get, lambda n: name_hit(n, oid))
                        if term:
                            print(f"[{key}] organism name {name!r} read as {term} ({how})")
                        return term

                    def name_over_id(oid, oname, a_term):
                        # An 'ID|Name' whose name names another taxon than its ID (a_term): the name wins, as
                        # Metabolights2REDU.organism_accession_terms decides organism text vs NCBITaxon accession
                        # (user decision: trust the name, for microbes the ID). The ID stays when it is below the
                        # named taxon (a strain or subspecies of the species), and for microbes (both in Bacteria,
                        # Archaea, Fungi or Viruses, the ID not above the name's taxon: strains of a collection
                        # reclassified since, 'Gardnerella vaginalis' with the ID of 'Gardnerella pickettii JCP7659';
                        # unlike MetaboLights, one species name with many strain IDs is no copy-paste in GNPS strain
                        # collections). A name that only names.dmp lists (a synonym) wins only over an ID above it
                        # or over the genus it names; the ID also stays without the pinned NCBI taxonomy
                        # (--path_to_ncbi_nodes_dmp). A name read by fallback_hit (white space, strain designation,
                        # one wrong character) counts as written that way.
                        if not oname or ncbi_lineage is None:
                            return a_term
                        direct = name_hit(oname, oid)
                        t_term = direct or (ncbi_name_resolver(oname) if ncbi_name_resolver is not None else None)
                        if not t_term:
                            t_term = direct = fallback_hit(oname, oid)
                        if not t_term or t_term == a_term:
                            return a_term
                        a_id, t_id = a_term.split("|", 1)[0].strip(), t_term.split("|", 1)[0].strip()
                        if ncbi_lineage.is_ancestor(t_id, a_id):
                            return a_term
                        group = ncbi_lineage.microbe_group(t_id)
                        if group is not None and ncbi_lineage.microbe_group(a_id) == group \
                                and not ncbi_lineage.is_ancestor(a_id, t_id):
                            return a_term
                        if direct:
                            return t_term
                        a_name = a_term.split("|", 1)[1].strip().lower()
                        if ncbi_lineage.is_ancestor(a_id, t_id) or a_name == oname.split()[0].lower():
                            return t_term
                        return a_term

                    def genus_of_genus_level(term, oid):
                        # the genus of an allowed term at genus level ('1931|Streptomyces sp.' -> '1883|Streptomyces')
                        name = term.split("|", 1)[1].strip() if "|" in term else ""
                        return name_hit(_genus_from_name(name), oid) if name and _genus_level_sp(name) else None

                    genus_relations = []

                    def in_genus(term, genus_term, genus_name):
                        # whether an allowed term is a taxon of the genus a text names: its name starts with the genus
                        # name, or it lies below the genus in the pinned NCBI taxonomy (where the caller loaded it)
                        name = term.split("|", 1)[-1].split()
                        if name and name[0].lower() == genus_name.lower():
                            return True
                        if not genus_relations:
                            genus_relations.append(ncbi_lineage if ncbi_lineage is not None else pinned_ncbi_relations())
                        relations = genus_relations[0]
                        return relations is not None and relations.is_ancestor(genus_term.split("|", 1)[0].strip(),
                                                                               term.split("|", 1)[0].strip())

                    def resolve_ncbi(observed_value: str) -> str:
                        # observed_value like '1931|Streptomyces sp.' or '1883|Streptomyces' or just 'Streptomyces sp.'
                        oid, oname = _split_id_name(observed_value)

                        # 1) 'Genus sp.' / 'Genus spp.', an answer at genus level: the genus, even if the exact value
                        #    is allowed ('1931|Streptomyces sp.' -> '1883|Streptomyces'). 'Genus sp. <strain, isolate,
                        #    clone>' names one organism: the allowed taxon the value names by steps 2-6 (the exact term
                        #    '1148|Synechocystis sp. PCC 6803', the name alone, a names.dmp name of it, the ID by
                        #    name_over_id; an ID of the genus-level 'Genus sp.' taxon stands for its genus there, the name
                        #    tells which organism of it) when the text names that taxon, or it is one of the genus and
                        #    no 'Genus sp.' of another designation; else the genus (a strain NCBI does not list, a taxon
                        #    of another genus, one at genus level, NCBI's 'unclassified Genus' or 'uncultured Genus sp.',
                        #    an ID of another strain). A designation that is a lab's numbering ('Bacillus sp. 1',
                        #    'Streptomyces sp. 12': REDU_conversion_functions.sp_generic_label) names the genus only,
                        #    unless the given ID is that taxon (NCBI lists such names for other submitters' isolates).
                        #    A 'Genus sp.' / 'Genus spp.' with the ID of a species or strain of that genus names that
                        #    organism, the ID being more specific ('502558|Eggerthella sp' is 'Eggerthella sp. YY7918';
                        #    REDU_conversion_functions.bare_sp_id_term, pinned NCBI taxonomy)
                        sp_genus = None
                        generic = False
                        if oname and _is_sp_or_spp(oname):
                            if oid and _genus_level_sp(oname):
                                own = bare_sp_id_term(oname, oid, allowed_by_id.get,
                                                      lambda n: allowed_by_name_ci.get(str(n).strip().lower()),
                                                      bare_sp_taxonomy)
                                if own:
                                    return own
                            sp_genus = name_hit(_genus_from_name(oname), oid)
                            if sp_genus and _genus_level_sp(oname):
                                return sp_genus  # e.g., '1883|Streptomyces'
                            generic = sp_generic_label(oname)
                        term = resolve_named(observed_value, oid, oname, sp_genus is not None)
                        t_id = term.split("|", 1)[0].strip()
                        if generic and term != missing_value and t_id != oid:
                            # (a genus name of two taxa without an ID to tell which: the NCBI genus of the taxon found,
                            # in the pinned taxonomy, when it has the text's genus name)
                            return sp_genus or pinned_genus_term(t_id, _genus_from_name(oname), allowed_by_id.get) \
                                or missing_value
                        if not sp_genus:
                            return term
                        t_name = term.split("|", 1)[-1].strip()
                        if term == missing_value or _genus_level_sp(t_name) or _NCBI_COLLECTING_TAXON.fullmatch(t_name):
                            return sp_genus
                        # a taxon the text itself names (an allowed name of it, a names.dmp name of it: 'Synechococcus
                        # sp. PCC 7502', now Pseudanabaenococcus habilitatus PCC 7502); else one of the genus the text
                        # names that is no strain of another designation
                        named = name_hit(oname, oid) or (ncbi_name_resolver(oname) if ncbi_name_resolver is not None else None)
                        if named != term and (not in_genus(term, sp_genus, _genus_from_name(oname))
                                              or other_designation(term, oid, oname)):
                            return sp_genus
                        return term

                    def other_designation(term, oid, oname):
                        # the taxon of the ID, whose NCBI name is 'Genus sp. Y', for a text 'Genus sp. X' of another
                        # designation (the ID of another strain: 'Pseudomonas sp. 2156' for 'Pseudomonas sp. 2 1 26'):
                        # the two name different organisms of the genus. The same designation written otherwise
                        # ('CNQ525' / 'CNQ-525', 'F0373' / 'oral taxon 186 str. F0373') is none; one designation
                        # inside another ('2' / '2156', 'PCC 680' / 'PCC 6803') is another strain
                        # (REDU_conversion_functions.sp_other_strain)
                        if term.split("|", 1)[0].strip() != oid:
                            return False
                        return sp_other_strain(oname, term.split("|", 1)[-1])

                    def resolve_named(observed_value, oid, oname, one_of_genus):
                        # 2) Exact allowed term
                        if observed_value in allowed_set:
                            return observed_value

                        # 3) Same ID → allowed term (unless the name names another taxon: name_over_id)
                        if oid and oid in allowed_by_id:
                            id_term = allowed_by_id[oid]
                            if one_of_genus:
                                id_term = genus_of_genus_level(id_term, oid) or id_term
                            return name_over_id(oid, oname, id_term)

                        # 4) Same Name (case-insensitive) → allowed term
                        if oname:
                            hit = name_hit(oname, oid)
                            if hit:
                                return hit

                        # 5) Name that NCBI lists for another (current) taxon, e.g. a
                        #    synonym or the name of a merged taxon (GNPS only; optional)
                        if ncbi_name_resolver is not None and oname:
                            hit = ncbi_name_resolver(oname)
                            if hit:
                                return hit

                        # 6) The name with its white space normalised, without a strain designation at its end, or
                        #    with one wrong character fixed (the GNPS harmonize step only)
                        hit = fallback_hit(oname, oid)
                        if hit:
                            return hit

                        # 7) No match
                        return missing_value

                    def resolve_ncbi_value(observed_value):
                        # several organisms ('ID|Name||ID|Name', also 'ID|Name, ID|Name'): each part as a value of
                        # its own; all parts must resolve (a sample of A and B is not recorded as A alone); each
                        # taxon once, in the submitter's order (ncbi_database.join_parts); placeholder parts
                        # ('ID|Name||missing value') are no organism
                        parts = submitted_taxonomy_parts(observed_value)
                        if len(parts) == 1 and parts[0] != str(observed_value).strip():
                            return resolve_ncbi(parts[0])
                        if not parts and MULTI_SEP in str(observed_value):
                            return missing_value
                        if len(parts) <= 1:
                            return resolve_ncbi(observed_value)
                        terms = [resolve_ncbi(p) for p in parts]
                        if any(t == missing_value for t in terms):
                            return missing_value
                        return join_parts(terms)

                    value_map = {ov: resolve_ncbi_value(ov) for ov in unique_values}

                else:
                    # MassSpectrometer: keep your original direct matching
                    value_map = {
                        observed_value: (
                            observed_value if observed_value in allowed_terms 
                            else next(
                                (y for y in allowed_terms if "|" in str(observed_value) and y.split("|", 1)[0] == str(observed_value).split("|", 1)[0]),
                                next(
                                    (y for y in allowed_terms if "|" in str(observed_value) and y.split("|", 1)[1] == str(observed_value).split("|", 1)[1]),
                                    missing_value
                                )
                            )
                        )
                        for observed_value in unique_values
                    }
                
            elif not allowed_terms:
                raise ValueError(f"No allowed values for {key}; was the allowed terms json updated from the ontologies?")
            elif allowed_terms[0] == '00':
                # Handling for terms where '00' indicates a placeholder for any value
                value_map = {x: x if pd.notna(x) and x != "" else missing_value for x in unique_values}
                if key in _IDENTIFIER_COLUMNS:
                    # placeholder words ('not applicable', 'ML import: not available', a Qiita
                    # study prefix with 'not applicable' as sample) identify nothing
                    value_map = {x: missing_value if _is_identifier_placeholder(x) else v
                                 for x, v in value_map.items()}
            elif allowed_terms[0] == 'datetime':
                # Flexible date/time normalisation to precision-preserving ISO-8601
                # (year / year-month / date / date+time / interval). Ambiguous slash
                # dates (both fields <= 12) are disambiguated per dataset: a value that
                # itself proves an order (a field > 12) always wins; otherwise the
                # dataset's other dates decide; NORMAN (European, ISO) defaults to
                # day-first, everything else (US MassIVE/Workbench) to month-first.
                _slash = re.compile(r'\s*(\d{1,2})[/-](\d{1,2})[/-]\d')

                def _infer_dayfirst(mid, values):
                    if str(mid).startswith('NORMAN'):
                        return True
                    dmy = mdy = False
                    for v in values:
                        m = _slash.match(str(v))
                        if not m:
                            continue
                        a, b = int(m.group(1)), int(m.group(2))
                        if a > 12 and b <= 12:
                            dmy = True
                        elif b > 12 and a <= 12:
                            mdy = True
                    if dmy and not mdy:
                        return True
                    return False  # month-first: US platform default (also for conflicts)

                out = df[key].copy()
                gid = (df['MassiveID'] if 'MassiveID' in df.columns
                       else pd.Series('', index=df.index)).fillna('')
                for mid, idx in out.groupby(gid).groups.items():
                    sub = df.loc[idx, key]
                    dayfirst = _infer_dayfirst(mid, sub.unique())
                    vmap = {
                        x: (normalize_redu_datetime(x, dayfirst=dayfirst) or missing_value)
                        if (pd.notna(x) and x != "") else missing_value
                        for x in sub.unique()
                    }
                    out.loc[idx] = sub.map(vmap)
                df[key] = out.fillna(missing_value).replace("", missing_value)
                continue
            elif allowed_terms[0] == 'numeric':
                # Numeric handling; optionally clamp to a plausible [min, max] range so
                # nonsensical values (e.g. an auto-miscalculated BMI) drop to missing.
                nums = pd.to_numeric(df[key], errors='coerce')
                lo, hi = value.get('min'), value.get('max')
                if lo is not None:
                    nums = nums.where(nums >= lo)
                if hi is not None:
                    nums = nums.where(nums <= hi)
                if key in _WHOLE_NUMBERS_WITHOUT_FRACTION:
                    # whole numbers are written '72', not '72.0' (other values as before)
                    nums = nums.astype(object).map(
                        lambda x: str(int(x)) if isinstance(x, (int, float, np.integer, np.floating))
                        and math.isfinite(x) and float(x).is_integer() else x)
                df[key] = nums.fillna(missing_value).replace("", missing_value)
                continue
            elif allowed_terms[0] == 'numeric|numeric':
                # Handling for numeric ranges split by '|'
                value_map = {
                    x: x if all(part.replace('.', '', 1).isdigit() or part.lstrip('-').replace('.', '', 1).isdigit() for part in str(x).split('|')) else missing_value
                    for x in unique_values
                }
            elif key == 'filename':
                # Specific handling for filenames with allowed extensions
                value_map = {
                    x: x if any(str(x).lower().endswith(ext.lower()) for ext in allowed_terms + other_allowed_file_extensions) 
                    else missing_value 
                    for x in unique_values
                }
            else:
                raise ValueError(f"No validation rule for {key} with allowed values {allowed_terms}")
            # Apply the mapping for columns, except 'numeric' which is handled separately
            if key != 'numeric':
                df[key] = df[key].map(value_map).fillna(missing_value).replace("", missing_value)
            
            # preparing valuemap for print into log
            filtered_map = {k: v for k, v in value_map.items() if k != v}
            if filtered_map: 
                print(f"{key}:")
                for k, v in filtered_map.items():
                    print(f"  {k}: {v}")
        elif key not in df.columns and value['generate'] == 'False':
            df[key] = value['missing']

    # where an isolate was isolated from is another organism than the isolate: an organism of
    # IsolationSourceNCBITaxonomy that is also in NCBITaxonomy (the same taxon ID) is left out, per organism
    if ISOLATION_SOURCE in df.columns and 'NCBITaxonomy' in df.columns:
        source = without_taxa_of_series(df[ISOLATION_SOURCE], df['NCBITaxonomy'])
        same = source != df[ISOLATION_SOURCE]
        if same.any():
            print(f"{ISOLATION_SOURCE}: the organisms that are also in NCBITaxonomy left out on {int(same.sum())} rows")
            df[ISOLATION_SOURCE] = source


    for key, value in allowedTerm_dict.items():
        if value['generate'] == 'True':
            missing_value = value['missing']

            if key == 'DataSource':
                df.loc[df['MassiveID'].str.startswith('NORMAN', na=False), 'DataSource'] = 'NORMAN'
                df.loc[df['MassiveID'].str.startswith('MSV', na=False), 'DataSource'] = 'GNPS'
                df.loc[df['MassiveID'].str.startswith('MTBLS', na=False), 'DataSource'] = 'MetaboLights'
                df.loc[df['MassiveID'].str.startswith('ST', na=False), 'DataSource'] = 'Workbench'


            if key == 'LifeStage':
                if df.columns.is_unique:
                    # the row-wise rule on the two columns it reads (df.apply(axis=1) builds a
                    # Series of all columns per row)
                    computed = pd.Series([age_category(age) if taxon == "9606|Homo sapiens" else value['missing']
                                          for age, taxon in zip(df['AgeInYears'], df['NCBITaxonomy'])],
                                         index=df.index, dtype=object)
                else:
                    computed = df.apply(lambda x: age_category(x['AgeInYears']) if x['NCBITaxonomy'] == "9606|Homo sapiens" else value['missing'], axis=1)
                computed = computed.replace('', value['missing'])  # age_category gives '' for non-numeric ages
                # a life stage recorded by the submitter is kept when the age gives none
                if '_recorded_LifeStage' in df.columns and len(value.get('allowed_values') or []) > 0:
                    recorded = df['_recorded_LifeStage'].where(df['_recorded_LifeStage'].isin(value['allowed_values']))
                    if not kwargs.get('life_stage_as_submitted', False) and ages_as_given is not None \
                            and ages_as_given.index.equals(df.index):
                        # The converters (Workbench, MetaboLights) fill LifeStage with the human
                        # age bracket of AgeInYears for every organism. LifeStage is defined for
                        # humans only, so on rows of other (or unknown) organisms a value that
                        # is just the bracket of the row's own age is not kept; a life stage
                        # stated otherwise (no age, or not the age's bracket) still is. GNPS
                        # sheets give the submitter's own value (life_stage_as_submitted=True).
                        brackets = {a: age_category(a) for a in ages_as_given.unique()}
                        derived = recorded.notna() & (df['NCBITaxonomy'] != "9606|Homo sapiens") \
                            & (recorded == ages_as_given.map(brackets))
                        recorded = recorded.mask(derived)
                    computed = computed.where(computed != value['missing'], recorded).fillna(value['missing'])
                df[key] = computed
            if key == 'UniqueSubjectID' and df.columns.is_unique:
                # the row-wise rule below on the two columns it reads
                # no dataset (MassiveID missing: 'not allowed'), no subject ID across datasets
                no_dataset = ('', allowedTerm_dict.get('MassiveID', {}).get('missing', 'not allowed'))
                df[key] = pd.Series([str(massive_id) + '_' + str(subject)
                                     if subject != allowedTerm_dict['SubjectIdentifierAsRecorded']['missing']
                                     and subject != ''
                                     and subject != 'not applicable'
                                     and not pd.isna(subject)
                                     and massive_id not in no_dataset
                                     else value['missing']
                                     for massive_id, subject in zip(df['MassiveID'], df['SubjectIdentifierAsRecorded'])],
                                    index=df.index, dtype=object)
            elif key == 'UniqueSubjectID':
                df[key] = df.apply(lambda x: str(x['MassiveID']) + '_' + str(x['SubjectIdentifierAsRecorded']) 
                                if x['SubjectIdentifierAsRecorded'] != allowedTerm_dict['SubjectIdentifierAsRecorded']['missing'] 
                                and x['SubjectIdentifierAsRecorded'] != ''
                                and x['SubjectIdentifierAsRecorded'] != 'not applicable' 
                                and not pd.isna(x['SubjectIdentifierAsRecorded']) 
                                else value['missing'], axis=1)
            if key == 'UBERONOntologyIndex':
                df = df.merge(uberon_ontology_table[['Label', 'UBERONOntologyIndex', 'Is Multicellular','Is Organ', 'Is Fluid']].drop_duplicates('Label'), left_on='UBERONBodyPartName', right_on='Label', how='left')
                # QC, blank, standard and product rows keep their own SampleTypeSub1 (e.g. pool_QC of plasma)
                sample_type = df['SampleType'] if 'SampleType' in df.columns else pd.Series('', index=df.index)
                not_qc = ~sample_type.astype(str).str.match(r'blank_|pool_QC|molecular_standards|reference material|personal care product')
                # only where Sub1 is still open: a given Sub1 (tissue, biofluid or a finer value such
                # as insect, plant_angiospermae) is kept. Both 'tissue' and 'biofluid' are valid
                # SampleTypeSub1 terms; the body part's flags (e.g. feces, which UBERON does not
                # place under bodily fluid but this table flags as fluid) only fill a missing Sub1
                sub1 = df['SampleTypeSub1'] if 'SampleTypeSub1' in df.columns else pd.Series(pd.NA, index=df.index, dtype=object)
                sub1_open = sub1.isna() | sub1.isin(['', value['missing'], allowedTerm_dict['SampleTypeSub1']['missing']])
                # a given 'tissue' on a fluid body part (blood plasma, urine) is still corrected;
                # feces keeps the given value (see above)
                tissue_on_fluid = (sub1 == 'tissue') & (df['Is Fluid'] == True) & (df['Label'] != 'feces')
                df.loc[(df['Is Multicellular'] == True) & not_qc & sub1_open, 'SampleTypeSub1'] = 'tissue'
                df.loc[(df['Is Fluid'] == True) & not_qc & (sub1_open | tissue_on_fluid), 'SampleTypeSub1'] = 'biofluid'
                df.drop(columns=['Label', 'Is Multicellular', 'Is Organ', 'Is Fluid'], inplace=True)
            if key == 'DOIDOntologyIndex':
                df = df.merge(doid_ontology_table[['Label', 'DOIDOntologyIndex']].drop_duplicates('Label'), left_on='DOIDCommonName', right_on='Label', how='left')
                df.drop(columns=['Label'], inplace=True)
            if key == 'CellLineIndex':
                if 'CellLine' not in df.columns:
                    df['CellLine'] = allowedTerm_dict.get('CellLine', {}).get('missing', missing_value)
                if cellosaurus_index is not None:
                    # looked up per distinct value: mapping with the ~160k-name dict itself
                    # turns it into a Series on every call
                    by_name = cellosaurus_index.by_name
                    accession = {v: by_name[v] for v in df['CellLine'].unique() if v in by_name}
                    df[key] = df['CellLine'].map(accession).fillna(missing_value)
                else:
                    df[key] = missing_value
            if key == 'HPOOntologyIndex':
                if 'HPOCommonName' not in df.columns:
                    df['HPOCommonName'] = allowedTerm_dict.get('HPOCommonName', {}).get('missing', missing_value)
                df = df.merge(hpo_ontology_table[['Label', 'HPOOntologyIndex']].drop_duplicates('Label'), left_on='HPOCommonName', right_on='Label', how='left')
                df.drop(columns=['Label'], inplace=True)
            if key == 'ExposureIndex':
                if 'Exposure' not in df.columns:
                    df['Exposure'] = allowedTerm_dict.get('Exposure', {}).get('missing', missing_value)
                # PECO is for plants: on other samples a PECO term becomes the ECTO term with
                # that name as label or exact synonym (e.g. PECO 'radiation exposure' ->
                # ECTO 'exposure to radiation'), or missing if ECTO has none
                if 'Synonym' in exposure_ontology_table.columns and 'SampleType' in df.columns:
                    ecto = exposure_ontology_table[exposure_ontology_table['ExposureIndex'].astype(str).str.startswith('ECTO:')]
                    ecto_by_name = {}
                    for label, synonym in zip(ecto['Label'].astype(str), ecto['Synonym'].astype(str)):
                        for name in (label, synonym):
                            ecto_by_name.setdefault(name.strip().lower(), set()).add(label)
                    peco_labels = set(exposure_ontology_table.loc[exposure_ontology_table['ExposureIndex'].astype(str).str.startswith('PECO:'), 'Label'].astype(str)) \
                        - set(ecto['Label'].astype(str))
                    not_plant = ~df['SampleType'].astype(str).str.lower().str.startswith('plant')
                    on_peco = not_plant & df['Exposure'].isin(peco_labels)
                    def _to_ecto(value):
                        labels = ecto_by_name.get(str(value).strip().lower(), set())
                        return next(iter(labels)) if len(labels) == 1 else allowedTerm_dict['Exposure']['missing']
                    df.loc[on_peco, 'Exposure'] = df.loc[on_peco, 'Exposure'].map(_to_ecto)
                # ECTO terms also give what the organism was exposed to (OrganismExposureAgent)
                agent_columns = [c for c in ('OrganismExposureAgentIndex', 'OrganismExposureAgent')
                                 if c in exposure_ontology_table.columns and c in allowedTerm_dict]
                df = df.drop(columns=[c for c in agent_columns if c in df.columns])
                df = df.merge(exposure_ontology_table[['Label', 'ExposureIndex'] + agent_columns].drop_duplicates('Label'), left_on='Exposure', right_on='Label', how='left')
                df.drop(columns=['Label'], inplace=True)
                for c in ('OrganismExposureAgentIndex', 'OrganismExposureAgent'):
                    if c in allowedTerm_dict and c not in df.columns:
                        df[c] = allowedTerm_dict[c]['missing']
            if key == 'ENVOEnvironmentBiomeIndex':
                df = df.merge(envome_biome_ontology_table[['Label', 'ENVOEnvironmentBiomeIndex']].drop_duplicates('Label'), left_on='ENVOEnvironmentBiome', right_on='Label', how='left')
                df.drop(columns=['Label'], inplace=True)
            if key == 'ENVOEnvironmentMaterialIndex':
                df = df.merge(envome_material_ontology_table[['Label', 'ENVOEnvironmentMaterialIndex']].drop_duplicates('Label'), left_on='ENVOEnvironmentMaterial', right_on='Label', how='left')
                df.drop(columns=['Label'], inplace=True)
            if key == 'IsolationSourceUBERONOntologyIndex':
                # (only the index: SampleTypeSub1 is that of the sample analyzed)
                df = df.merge(uberon_ontology_table[['Label', 'UBERONOntologyIndex']].drop_duplicates('Label').rename(columns={'UBERONOntologyIndex': key}),
                              left_on='IsolationSourceUBERONBodyPartName', right_on='Label', how='left')
                df.drop(columns=['Label'], inplace=True)
            if key == 'IsolationSourceENVOEnvironmentMaterialIndex':
                df = df.merge(envome_material_ontology_table[['Label', 'ENVOEnvironmentMaterialIndex']].drop_duplicates('Label').rename(columns={'ENVOEnvironmentMaterialIndex': key}),
                              left_on='IsolationSourceENVOEnvironmentMaterial', right_on='Label', how='left')
                df.drop(columns=['Label'], inplace=True)
            if key == 'ENVOMediumScaleIndex':
                df = df.merge(envome_medium_scale_ontology_table[['Label', 'ENVOMediumScaleIndex']].drop_duplicates('Label'), left_on='ENVOMediumScale', right_on='Label', how='left')
                df.drop(columns=['Label'], inplace=True)
            if key == 'ENVOLocalScaleIndex':
                df = df.merge(envome_local_scale_ontology_table[['Label', 'ENVOLocalScaleIndex']].drop_duplicates('Label'), left_on='ENVOLocalScale', right_on='Label', how='left')
                df.drop(columns=['Label'], inplace=True)
            if key == 'ENVOBroadScaleIndex':
                df = df.merge(envome_broad_scale_ontology_table[['Label', 'ENVOBroadScaleIndex']].drop_duplicates('Label'), left_on='ENVOBroadScale', right_on='Label', how='left')
                df.drop(columns=['Label'], inplace=True)
            if key == 'USI' and add_usi == True:
                df['filename'] = df['filename'].apply(process_filename)
                df['USI'] = 'mzspec:' + df['MassiveID'] + ':' + df['filename']
            if key == 'NCBIRank':
                df['TaxonID'] = df['NCBITaxonomy'].apply(lambda x: x.split('|')[0] if '|' in x else 'missing value')
                df['TaxonID'] = df['TaxonID'].astype(str)
                # several organisms: no single taxon (rank and division per part, below)
                multi_values = [v for v in df['NCBITaxonomy'].unique() if is_multi(v)]
                multi_rows = df['NCBITaxonomy'].isin(multi_values)
                df.loc[multi_rows, 'TaxonID'] = 'missing value'
                part_ids = pd.Series([t for v in multi_values for t in taxon_ids(v)], dtype=object)
                ncbi_rows = NCBIRankDivision_table
                if ncbi_taxon_index is not None:
                    # only the table rows of these taxa: the same left merge as with the
                    # whole table, without hashing its ~3M TaxonIDs on every call
                    ncbi_rows = _ncbi_rows_of(NCBIRankDivision_table, ncbi_taxon_index, df['TaxonID'])
                df = pd.merge(df, ncbi_rows, on='TaxonID', how='left')
                df.drop(columns=['TaxonID'], inplace=True)

                if 'NCBIRank' in df.columns:
                    df['NCBIRank'] = df['NCBIRank'].fillna(missing_value)
                else:
                    df['NCBIRank'] = missing_value

                if 'NCBIDivision' in df.columns:
                    df['NCBIDivision'] = df['NCBIDivision'].fillna(missing_value)
                else:
                    df['NCBIDivision'] = missing_value

                if multi_values:
                    # one rank / division per organism, in the order of the NCBITaxonomy parts
                    part_rows = (_ncbi_rows_of(NCBIRankDivision_table, ncbi_taxon_index, part_ids)
                                 if ncbi_taxon_index is not None else NCBIRankDivision_table[NCBIRankDivision_table['TaxonID'].isin(set(part_ids))])
                    rank_of = dict(zip(part_rows['TaxonID'], part_rows['NCBIRank']))
                    division_of_id = dict(zip(part_rows['TaxonID'], part_rows['NCBIDivision']))
                    multi_rows = df['NCBITaxonomy'].isin(multi_values)
                    for column, lookup in (('NCBIRank', rank_of), ('NCBIDivision', division_of_id)):
                        per_value = {v: aligned([lookup.get(t) for t in taxon_ids(v)]) for v in multi_values}
                        df.loc[multi_rows, column] = df.loc[multi_rows, 'NCBITaxonomy'].map(per_value)



    # generated columns without their own rule (e.g. OrganismExposureAgent without an
    # exposure table) still have to exist
    for key, value in allowedTerm_dict.items():
        if value.get('generate') == 'True' and key not in df.columns and key not in ('USI',):
            if key in ('OrganismExposureAgent', 'OrganismExposureAgentIndex'):
                df[key] = value['missing']
            elif key == 'FoodOntologyIndex':
                # WebProtege ID of the GFOP food ontology term (translation_sheets/GFOP_food_ontology.csv)
                terms = df['FoodOntologyTerm'] if 'FoodOntologyTerm' in df.columns else pd.Series('', index=df.index)
                df[key] = terms.map(_gfop_index()).fillna(value['missing'])

    # ontology indexes of unmatched labels come out of the merges as NaN
    for key, value in allowedTerm_dict.items():
        if key.endswith('OntologyIndex') or key in ('ExposureIndex', 'OrganismExposureAgent', 'OrganismExposureAgentIndex',
                                                    'IsolationSourceENVOEnvironmentMaterialIndex') \
                or (key.startswith('ENVO') and key.endswith('Index')):
            if key in df.columns:
                df[key] = df[key].fillna(value['missing']).replace('', value['missing'])

    # Ensure the dataframe contains only the columns specified in the dictionary
    keys_to_include = [key for key in allowedTerm_dict.keys() if key != 'USI' or add_usi or keep_usi]

    ignored_columns = input_columns - set(keys_to_include)
    
    for column in ignored_columns:
        print(f"IGNORED:  {column}")

    #remove rows if not enough metadata are present
    def count_non_missing_specific(row, cols):
        return (row[cols] != "missing value").sum()
    
    columns_to_consider = ["SampleType", "SampleTypeSub1", "NCBITaxonomy", "UBERONBodyPartName", "BiologicalSex", 
                           "AgeInYears", "LifeStage", "Country", "HealthStatus", "SampleExtractionMethod", 
                           "SampleCollectionMethod", "ComorbidityListDOIDIndex", "DOIDCommonName", "HPOCommonName", "Exposure", 
                           "DepthorAltitudeMeters", "HumanPopulationDensity", "LatitudeandLongitude",
                           "ENVOEnvironmentBiome", "ENVOEnvironmentMaterial"] + STRIP_ISOLATION_SOURCE
    
    original_row_count = df.shape[0]
    columns_to_check = [col for col in columns_to_consider if col in df.columns]
    if len(df) == 0 or isinstance(df.iloc[0].dtype, pd.api.extensions.ExtensionDtype):
        df = df[(df.apply(count_non_missing_specific, cols=columns_to_check, axis=1) >= 1)]
    else:
        # the same count for all rows at once. The rows df.apply passes hold the cells as
        # objects (the frame's columns are not all of one extension dtype), so the cells are
        # compared as objects here too: NaN, None and pd.NA count as not missing
        df = df[(df[columns_to_check].astype(object) != "missing value").sum(axis=1) >= 1]

    print(f"Number of rows removed due to not enough metadata: {original_row_count - df.shape[0]}")
    print(f"Returning {len(df)} rows!")

    return df[keys_to_include]

# Values submitters commonly use for a term the vocabulary spells differently
# (keys are lower case). Only mappings that do not change the meaning.
_VALUE_SYNONYMS = {
    'Country': {
        'usa': 'United States of America', 'us': 'United States of America', 'u.s.a.': 'United States of America',
        'u.s.': 'United States of America', 'united states': 'United States of America',
        'uk': 'United Kingdom', 'u.k.': 'United Kingdom', 'great britain': 'United Kingdom',
        'england': 'United Kingdom', 'scotland': 'United Kingdom', 'wales': 'United Kingdom',
        'northern ireland': 'United Kingdom',
    },
    'IonizationSourceAndPolarity': {
        'heated electrospray ionization (positive)': 'electrospray ionization (positive)',
        'heated electrospray ionization (negative)': 'electrospray ionization (negative)',
        'heated electrospray ionization (alternating)': 'electrospray ionization (alternating)',
    },
    'ChromatographyAndPhase': {
        'reverse phase (polar c18)': 'reverse phase (C18)',
    },
    # first-generation template values
    'SampleType': {
        'standard': 'molecular_standards', 'culture_mamalian': 'culture_mammalian',
        'deceased animal': 'animal', 'blank': 'blank_analysis', 'pool': 'pool_QC', 'qc': 'pool_QC',
        'animal_mammalia_nos': 'animal', 'human': 'animal',
    },
    'SmokingStatus': {
        'nonsmoker': 'non-smoker', 'non smoker': 'non-smoker',
        'ex-smoker': 'former smoker', 'exsmoker': 'former smoker', 'ex smoker': 'former smoker',
    },
}

# Prefixes of the first-generation ReDU template (e.g. 'ATTRIBUTE_Curated_SampleType',
# 'Subject_BiologicalSex', 'ATTRIBUTE_Analysis_MassSpectrometer') in front of current
# field names, and first-generation names that differ from the current ones (lower case,
# after the prefixes). Same rules as template_column in gnps_downloader.py.
_LEGACY_COLUMN_PREFIXES = ('ATTRIBUTE_Curated_', 'ATTRIBUTE_', 'Subject_', 'Analysis_', 'Study_', 'Curated_')
_LEGACY_COLUMN_ALIASES = {
    'sex': 'BiologicalSex', 'health': 'HealthStatus', 'bodypartontologyname': 'UBERONBodyPartName',
    'bmi': 'BodyMassIndex', 'diseasecommonname': 'DOIDCommonName', 'smoking': 'SmokingStatus',
    'sampletype_sub1': 'SampleTypeSub1', 'comorbidity_listdoids': 'ComorbidityListDOIDIndex',
    'massiveid': 'MassiveID', 'datasetaccession': 'MassiveID',
}


_GFOP_SHEET = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'translation_sheets', 'GFOP_food_ontology.csv')


def _gfop_index():
    """GFOP food ontology label -> FoodOntologyIndex (empty if the sheet is missing)."""
    if 'gfop_index' not in _LOOKUP_CACHE:
        try:
            sheet = pd.read_csv(_GFOP_SHEET, dtype=str, keep_default_na=False)
            _LOOKUP_CACHE['gfop_index'] = dict(zip(sheet['Label'], sheet['FoodOntologyIndex']))
        except (OSError, KeyError, pd.errors.EmptyDataError):
            print(f"WARNING: {_GFOP_SHEET} not readable; FoodOntologyIndex stays missing")
            _LOOKUP_CACHE['gfop_index'] = {}
    return _LOOKUP_CACHE['gfop_index']


def _template_column(column, allowedTerm_dict):
    """The current field a (first-generation or differently cased) column name stands for."""
    by_lower = _cached('keys_lower', allowedTerm_dict, lambda: {k.lower(): k for k in allowedTerm_dict})
    name = str(column).strip()
    while True:
        if name in allowedTerm_dict:
            return name
        if name.lower() in by_lower:
            return by_lower[name.lower()]
        if name.lower() in _LEGACY_COLUMN_ALIASES:
            return _LEGACY_COLUMN_ALIASES[name.lower()]
        for prefix in _LEGACY_COLUMN_PREFIXES:
            if name.lower().startswith(prefix.lower()) and len(name) > len(prefix):
                name = name[len(prefix):]
                break
        else:
            return None


_LATLON_NUMBER = r'([+-]?\d{1,3}(?:\.\d+)?)'
_LATLON = re.compile(r'^' + _LATLON_NUMBER + r'\s*°?\s*([NnSs])?\s*(AND|and|And|[;,|/]|\s)\s*'
                     + _LATLON_NUMBER + r'\s*°?\s*([EeWw])?$')


def _is_redu_latlon(value):
    """The format the ReDU validation accepts ('lat|lon', numbers without a '+')."""
    parts = str(value).split('|')
    return len(parts) == 2 and all(re.fullmatch(r'-?\d+(\.\d+)?', p) for p in parts)


def parse_latitude_longitude(value):
    """'lat|lon' from '32.842;-117.258', '48.2,-124.3', '32.842AND-117.258', '-33.99|+25.70',
    '41.8719 N, 12.5674 E' or '48.27 -124.45' (latitude first; N/S/E/W give the sign and
    exclude a minus sign). The numbers keep their digits. None if the value is not such a
    pair or out of range."""
    m = _LATLON.match(str(value).strip())
    if not m:
        return None
    lat, ns, separator, lon, ew = m.groups()
    if separator == ',' and '.' not in lat + lon:
        return None  # '41,87' may be one number with a decimal comma
    out = []
    for number, hemisphere, negative, limit in ((lat, ns, 'sS', 90), (lon, ew, 'wW', 180)):
        number = number.lstrip('+')
        if hemisphere:
            if number.startswith('-'):
                return None
            if hemisphere in negative:
                number = '-' + number
        if abs(float(number)) > limit:
            return None
        out.append(number)
    return '|'.join(out)


def _load_ms_aliases():
    path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'translation_sheets', 'GNPS_MassSpectrometer_aliases.csv')
    if not os.path.exists(path):
        return {}
    table = pd.read_csv(path, dtype=str).dropna()
    return {g.strip().lower(): r.strip() for g, r in zip(table['GNPS'], table['REDU'])}


def _diet_key(value):
    return re.sub(r'[\s_\-]+', '', str(value).strip().lower())


def _solvent_key(value):
    return str(value).strip().lower()

_MISSING_LIKE = {'', 'nan', 'na', 'n/a', 'none', 'missing value', 'missing_value', 'not specified', 'not collected',
                 'not applicable', 'not_applicable', 'ml import: not available', 'not available'}


def _is_missing_like(value):
    return str(value).strip().lower() in _MISSING_LIKE


# Numeric columns whose whole numbers are written without '.0' (a time point of 72 h is
# '72'); the other numeric columns keep their former representation
_WHOLE_NUMBERS_WITHOUT_FRACTION = ('TimePointHoursRelativeToEvent',)


# Identifier columns ('00' = any text) in which placeholder words stand for no identifier
_IDENTIFIER_COLUMNS = ('SubjectIdentifierAsRecorded', 'qiita_sample_name')
_IDENTIFIER_PLACEHOLDERS = _MISSING_LIKE | {
    'not provided', 'not reported', 'not recorded', 'not determined', 'not known', 'unknown',
    'missing', 'null', 'n.a.', 'n.a', 'restricted access', '-', '--', '?'}


def _is_identifier_placeholder(value):
    """True for a placeholder instead of an identifier: 'not applicable', 'Not_Applicable',
    'ML import: not available', or a Qiita study prefix followed by placeholders only
    ('14577.not applicable .not applicable', '14577.T.not reported' is not: 'T' is a value)."""
    text = re.sub(r'[\s_]+', ' ', str(value)).strip().lower()
    if text in _IDENTIFIER_PLACEHOLDERS:
        return True
    m = re.fullmatch(r'\d+\.(.+)', text)
    if m:
        parts = [p.strip() for p in m.group(1).split('.')]
        return all(p in _IDENTIFIER_PLACEHOLDERS for p in parts)
    return False


# First-generation sheets name the organism as SampleType (ATTRIBUTE_Curated_SampleType
# 'mouse', 'rat', 'human', 'Solanum lycopersicum'). Common names -> (SampleType, organism),
# as the curated (current-template) rows of the same datasets have them ('animal',
# '10088|Mus' / '10114|Rattus'); the organism fills an empty NCBITaxonomy.
_LEGACY_SAMPLE_TYPES = {
    'mouse': ('animal', '10088|Mus'), 'mice': ('animal', '10088|Mus'),
    'rat': ('animal', '10114|Rattus'), 'human': ('animal', '9606|Homo sapiens'),
    'animal_insecta_nos': ('animal', None), 'soil': ('environmental', None),
}
_SCIENTIFIC_NAME = re.compile(r'^[A-Z][a-z]+(?: [a-z][a-z-]+)?$')


def _legacy_sample_type_organisms(df, column, allowedTerm_dict):
    """Organisms given as SampleType in a first-generation column (see _LEGACY_SAMPLE_TYPES):
    a common name becomes its SampleType and fills an empty NCBITaxonomy with its taxon; a
    scientific name ('Zea mays', 'Penicillium') that is no SampleType term fills an empty
    NCBITaxonomy (the NCBI validation resolves it or drops it) and gives no SampleType."""
    values = df[column]
    text = values.astype(str).str.strip()
    lower = text.str.lower()
    common = lower.map({k: v for k, v in _LEGACY_SAMPLE_TYPES.items()})
    sample_types = {str(t).lower() for t in allowedTerm_dict['SampleType']['allowed_values']} \
        | {str(t).lower() for t in allowedTerm_dict['SampleTypeSub1']['allowed_values']} \
        | set(_VALUE_SYNONYMS['SampleType'])
    scientific = text.str.match(_SCIENTIFIC_NAME) & ~lower.isin(sample_types) & common.isna() \
        & ~values.map(_is_missing_like)
    if not common.notna().any() and not scientific.any():
        return df
    organism = common.map(lambda v: v[1] if isinstance(v, tuple) else None)
    organism = organism.where(organism.notna(), text.where(scientific))
    if organism.notna().any():
        if 'NCBITaxonomy' not in df.columns:
            df['NCBITaxonomy'] = ''
        fill = organism.notna() & df['NCBITaxonomy'].map(_is_missing_like)
        if fill.any():
            df.loc[fill, 'NCBITaxonomy'] = organism[fill]
            print(f"NCBITaxonomy of {int(fill.sum())} rows from the organism in legacy column {column}")
    is_common = common.notna()
    if is_common.any():
        df.loc[is_common, column] = common[is_common].map(lambda v: v[0])
    return df


def normalize_gnps_metadata(df, allowedTerm_dict, dataset_accession=None, uberon_synonym_table=None,
                            envo_tables=(), NCBIRankDivision_table=None):
    """Repairs common, unambiguous deviations in submitted GNPS/ReDU metadata files before
    validation, so the values are kept instead of being replaced by 'missing value'."""
    df = df.copy()

    # first-generation template column names (and case variants) -> current names; when
    # the current column exists too, its empty values are taken from the old one
    for column in list(df.columns):
        if column in allowedTerm_dict:
            continue
        name = _template_column(column, allowedTerm_dict)
        if name is None or name == column:
            continue
        if name == 'SampleType':
            df = _legacy_sample_type_organisms(df, column, allowedTerm_dict)
        if name not in df.columns:
            df = df.rename(columns={column: name})
            print(f"Renamed legacy column {column} -> {name}")
        else:
            empty = df[name].map(_is_missing_like) & ~df[column].map(_is_missing_like)
            if empty.any():
                df.loc[empty, name] = df.loc[empty, column]
                print(f"Filled {int(empty.sum())} empty {name} values from legacy column {column}")
            df = df.drop(columns=[column])

    # files without a dataset column: the accession is known from the file name
    if 'MassiveID' not in df.columns and dataset_accession:
        df['MassiveID'] = dataset_accession
        print(f"Added MassiveID {dataset_accession} from the file name")
    elif dataset_accession:
        # empty dataset cells are the file's dataset too (the matcher matches such rows in it;
        # the accession also builds UniqueSubjectID, which was 'not allowed_<subject>')
        empty = df['MassiveID'].map(_is_missing_like)
        if empty.any():
            df.loc[empty, 'MassiveID'] = dataset_accession
            print(f"MassiveID {dataset_accession} from the file name for {int(empty.sum())} rows without one")

    def _ci_lookup(mapping, value):
        return mapping.get(str(value).strip().lower())

    for key, synonyms in _VALUE_SYNONYMS.items():
        if key in df.columns:
            df[key] = df[key].map(lambda v: _ci_lookup(synonyms, v) or v)

    # SampleTypeSub1 values (e.g. biofluid, tissue, reference material_chemicalstandard)
    # entered as SampleType: move them, and infer SampleType only when unambiguous
    if 'SampleType' in df.columns:
        sample_types = set(allowedTerm_dict['SampleType']['allowed_values'])
        sub1_types = set(allowedTerm_dict['SampleTypeSub1']['allowed_values'])
        if 'SampleTypeSub1' not in df.columns:
            df['SampleTypeSub1'] = ''
        is_sub1 = df['SampleType'].isin(sub1_types - sample_types) & df['SampleTypeSub1'].map(_is_missing_like)
        if is_sub1.any():
            moved = df.loc[is_sub1, 'SampleType']
            df.loc[is_sub1, 'SampleTypeSub1'] = moved
            prefix = moved.str.rsplit('_', n=1).str[0]
            inferred = prefix.where(prefix.isin(sample_types) & (prefix != moved), '')
            # biofluid/tissue: animal when the organism is an animal
            if NCBIRankDivision_table is not None and 'NCBITaxonomy' in df.columns:
                division = _cached('division', NCBIRankDivision_table, lambda: dict(zip(
                    NCBIRankDivision_table['TaxonID'].astype(str), NCBIRankDivision_table['NCBIDivision'])))
                animal_divisions = {'Primates', 'Rodents', 'Mammals', 'Vertebrates', 'Invertebrates'}
                values = df.loc[is_sub1, 'NCBITaxonomy'].astype(str)
                taxid = values.str.split('|').str[0].str.strip()
                # the divisions of these taxa only: mapping with the ~3M-entry dict itself
                # turns it into a Series on every call
                division_here = {t: division[t] for t in taxid.unique() if t in division}
                is_animal = taxid.map(division_here).isin(animal_divisions)
                # several organisms ('ID|Name||ID|Name'): animal only when every organism is one
                several = {v: submitted_taxonomy_parts(v) for v in values.unique()}
                several = {v: all(division.get(p.split('|', 1)[0].strip()) in animal_divisions for p in parts)
                           for v, parts in several.items() if len(parts) > 1}
                if several:
                    is_animal = is_animal.mask(values.isin(several), values.map(several)).astype(bool)
                inferred = inferred.mask(moved.isin(['biofluid', 'tissue']) & is_animal, 'animal')
            df.loc[is_sub1, 'SampleType'] = inferred
            print(f"Moved {is_sub1.sum()} SampleTypeSub1 values out of SampleType")

    # 'healthy' entered as a disease
    if 'DOIDCommonName' in df.columns:
        healthy = df['DOIDCommonName'].astype(str).str.strip().str.lower() == 'healthy'
        if healthy.any():
            if 'HealthStatus' not in df.columns:
                df['HealthStatus'] = ''
            no_status = df['HealthStatus'].map(_is_missing_like)
            df.loc[healthy & no_status, 'HealthStatus'] = 'healthy'
            df.loc[healthy, 'DOIDCommonName'] = ''

    # 'disease NOS' says the subject is ill without naming the disease
    if 'DOIDCommonName' in df.columns:
        nos = df['DOIDCommonName'].astype(str).str.strip().str.lower() == 'disease nos'
        if nos.any():
            if 'HealthStatus' not in df.columns:
                df['HealthStatus'] = ''
            fill = nos & df['HealthStatus'].map(_is_missing_like)
            df.loc[fill, 'HealthStatus'] = 'unhealthy (NOS)'
            if fill.any():
                print(f"HealthStatus 'unhealthy (NOS)' for {int(fill.sum())} rows with DOIDCommonName 'disease NOS'")

    # 'lat | lon' with spaces, Excel float years
    if 'LatitudeandLongitude' in df.columns:
        df['LatitudeandLongitude'] = df['LatitudeandLongitude'].astype(str).str.replace(r'\s*\|\s*', '|', regex=True).str.strip()
        # other separators / hemisphere letters (values already in ReDU format are kept)
        def _latlon(v):
            if _is_redu_latlon(v) or _is_missing_like(v):
                return v
            return parse_latitude_longitude(v) or v
        df['LatitudeandLongitude'] = df['LatitudeandLongitude'].map(_latlon)

    # MassSpectrometer given as a bare MS accession, as the bare label of an allowed term
    # or as a first-generation name ('Maxis_ImpactHD'; translation_sheets/GNPS_MassSpectrometer_aliases.csv)
    if 'MassSpectrometer' in df.columns and 'MassSpectrometer' in allowedTerm_dict:
        ms_terms = allowedTerm_dict['MassSpectrometer']['allowed_values']

        def _build_ms():
            by_id, by_label = {}, {}
            for t in ms_terms:
                label, _, acc = str(t).partition('|')
                by_id[acc.strip().upper()] = t
                by_label.setdefault(label.strip().lower(), set()).add(t)
            return set(ms_terms), by_id, by_label, _load_ms_aliases()
        ms_set, ms_by_id, ms_by_label, ms_aliases = _cached('ms_lookup', ms_terms, _build_ms)

        def _ms(v):
            if v in ms_set or _is_missing_like(v):
                return v
            text = str(v).strip()
            if re.fullmatch(r'(?i)MS[:_]\d{7}', text):
                return ms_by_id.get('MS:' + text[3:], v)
            label = ms_by_label.get(text.lower())
            if label and len(label) == 1:
                return next(iter(label))
            alias = ms_aliases.get(text.lower())
            return alias if alias in ms_set else v
        df['MassSpectrometer'] = df['MassSpectrometer'].map(_ms)

    # SampleExtractionMethod: a bare solvent is the 100% term ('methanol' -> 'methanol (100%)'),
    # '<solvent> (80%)' the solvent-water ratio ('acetonitrile-water (4:1)'); only when that
    # term is allowed
    if 'SampleExtractionMethod' in df.columns and 'SampleExtractionMethod' in allowedTerm_dict:
        sem_terms = allowedTerm_dict['SampleExtractionMethod']['allowed_values']
        sem_set = set(sem_terms)

        def _sem(v):
            if v in sem_set or _is_missing_like(v):
                return v
            text = _solvent_key(v)
            if re.fullmatch(r'[a-z][a-z \-]*[a-z]', text) and f'{text} (100%)' in sem_set:
                return f'{text} (100%)'
            m = re.fullmatch(r'([a-z][a-z \-]*[a-z])\s*\((\d{1,2})\s*%\)', text)
            if m and m.group(1) != 'water':
                part = int(m.group(2))
                g = math.gcd(part, 100 - part)
                term = f'{m.group(1)}-water ({part // g}:{(100 - part) // g})'
                if term in sem_set:
                    return term
            return v
        df['SampleExtractionMethod'] = df['SampleExtractionMethod'].map(_sem)

    # Diet: spelling variants that differ only in case, spaces, '_' or '-' ('High_Fat_Diet')
    if 'Diet' in df.columns and 'Diet' in allowedTerm_dict:
        diet_terms = allowedTerm_dict['Diet']['allowed_values']

        def _build_diet():
            keys = {}
            for t in diet_terms:
                keys.setdefault(_diet_key(t), set()).add(t)
            return keys
        diet_keys = _cached('diet_keys', diet_terms, _build_diet)
        diet_set = set(diet_terms)

        def _diet(v):
            if v in diet_set or _is_missing_like(v):
                return v
            hits = diet_keys.get(_diet_key(v), set())
            return next(iter(hits)) if len(hits) == 1 else v
        df['Diet'] = df['Diet'].map(_diet)
    if 'YearOfAnalysis' in df.columns:
        df['YearOfAnalysis'] = df['YearOfAnalysis'].astype(str).str.replace(r'^\s*(\d{4})\.0+\s*$', r'\1', regex=True)

    # UBERON: exact synonyms and case variants of labels -> label (also in the columns with its terms)
    uberon_columns = [c for c in columns_with_terms_of(allowedTerm_dict, 'UBERONBodyPartName') if c in df.columns]
    if uberon_synonym_table is not None and uberon_columns:
        labels = _cached('uberon_labels', uberon_synonym_table,
                         lambda: set(allowedTerm_dict['UBERONBodyPartName']['allowed_values']))

        def _build_to_labels():
            to_labels = {}
            for label, synonym in zip(uberon_synonym_table['Label'].astype(str), uberon_synonym_table['Synonym'].astype(str)):
                if label in labels:
                    to_labels.setdefault(label.strip().lower(), set()).add(label)
                    # short synonyms are abbreviations ('NA' = nasal artery) that collide
                    # with placeholders, so only longer synonyms are used
                    if len(synonym.strip()) > 3:
                        to_labels.setdefault(synonym.strip().lower(), set()).add(label)
            return to_labels
        to_labels = _cached('uberon_synonyms', uberon_synonym_table, _build_to_labels)
        for column in uberon_columns:
            df[column] = df[column].map(
                lambda v: v if v in labels or _is_missing_like(v) else (
                    next(iter(to_labels[str(v).strip().lower()]))
                    if len(to_labels.get(str(v).strip().lower(), ())) == 1 else v))

    # ENVO fields given as ontology IDs (ENVO:00000446) -> label
    envo_by_id = {}
    for table in envo_tables:
        index_column = [c for c in table.columns if c.endswith('Index')][0]
        for index, label in zip(table[index_column].astype(str), table['Label'].astype(str)):
            envo_by_id[index.replace('_', ':')] = label
    for key in ['ENVOEnvironmentBiome', 'ENVOEnvironmentMaterial', 'ENVOBroadScale', 'ENVOLocalScale', 'ENVOMediumScale',
                'IsolationSourceENVOEnvironmentMaterial']:
        if key in df.columns:
            df[key] = df[key].map(lambda v: envo_by_id.get(str(v).strip(), v))

    return df


_NCBI_NAME_CLASSES = ('scientific name', 'synonym', 'equivalent name')


def build_ncbi_name_resolver(names_dmp, names, allowed_ncbi_terms, NCBIRankDivision_table=None):
    """Resolver for NCBI names that are not the scientific name of an allowed taxon:
    - a synonym / equivalent name / scientific name of exactly one taxon (in names.dmp)
      that is allowed, e.g. 'Clostridium coccoides' -> '33035|Blautia producta', including
      names of merged taxa ('37480|Streptomyces pluricolorescens');
    - a strain or other name below species ('Staphylococcus aureus JE3') that NCBI does not
      list at all -> its species, when the first two words are the name of exactly one
      allowed taxon of rank species. Hybrids ('A b x c') are not reduced.
    Only names written like a scientific name are considered (at least two words, the
    first a capitalised genus: not placeholders such as 'not specified', which NCBI lists
    as a synonym of 32644 'unidentified'). names: the names that occur in the metadata, so
    only those are kept from names.dmp."""
    def scientific_like(name):
        name = str(name).strip()
        return re.match(r'^[A-Z][a-z]+(\s+\S+)+$', name) is not None and not _is_missing_like(name)
    wanted = {str(n).strip().lower() for n in names if scientific_like(n)}
    if not wanted:
        return None
    resolving = collections.defaultdict(set)
    known = set()
    with open(names_dmp, encoding='utf-8', errors='replace') as f:
        for line in f:
            parts = line.split('\t|\t')
            if len(parts) < 4:
                continue
            name = parts[1].strip().lower()
            if name not in wanted:
                continue
            known.add(name)
            name_class = parts[3].rstrip('\t|\n').strip()
            if name_class in _NCBI_NAME_CLASSES:
                resolving[name].add(parts[0].strip())
    allowed_by_id = {}
    allowed_by_name = collections.defaultdict(set)
    for term in allowed_ncbi_terms:
        if '|' in term:
            tid, tname = term.split('|', 1)
            allowed_by_id[tid.strip()] = term
            allowed_by_name[tname.strip().lower()].add(term)
    rank = {}
    if NCBIRankDivision_table is not None:
        rank = dict(zip(NCBIRankDivision_table['TaxonID'].astype(str), NCBIRankDivision_table['NCBIRank'].astype(str)))
    print(f"NCBI names.dmp: {len(resolving)} of {len(wanted)} metadata names are NCBI names")

    def resolve(name):
        name = str(name).strip()
        key = name.lower()
        if key not in wanted or not scientific_like(name):
            return None
        ids = resolving.get(key)
        if ids:
            taxid = next(iter(ids))
            # never to the catch-all 'unidentified' taxon
            return allowed_by_id.get(taxid) if len(ids) == 1 and taxid != '32644' else None
        if key in known:
            return None  # listed by NCBI, but not as a name of one taxon
        words = name.split()
        if len(words) < 3 or re.search(r'(?i)\sx\s|×', name):
            return None
        species = allowed_by_name.get(' '.join(words[:2]).lower(), set())
        if len(species) != 1:
            return None
        term = next(iter(species))
        return term if rank.get(term.split('|', 1)[0].strip()) == 'species' else None
    return resolve


def _ncbi_names_in_files(paths, allowedTerm_dict):
    """Organism names in the NCBITaxonomy columns (any template name) of metadata files, and in the columns with its
    terms (IsolationSourceNCBITaxonomy)."""
    names = set()
    ncbi_columns = set(columns_with_terms_of(allowedTerm_dict, 'NCBITaxonomy'))
    for path in paths:
        try:
            header = pd.read_csv(path, sep='\t', dtype=str, nrows=0).columns
            columns = [c for c in header if _template_column(c, allowedTerm_dict) in ncbi_columns]
            if not columns:
                continue
            table = pd.read_csv(path, sep='\t', dtype=str, usecols=columns)
        except Exception:
            continue
        for column in columns:
            for value in table[column].dropna().unique():
                names.add(str(value).split('|', 1)[-1].strip())
                # several organisms: the name of each part
                names.update(p.split('|', 1)[-1].strip() for p in submitted_taxonomy_parts(value))
    return names


# sheet columns that are no sample-level evidence of a reference material: ReDU's sample type vocabulary
# (a GNPS 'reference material' is any standard or product), the file name (its directory may name a
# material; the base name is read from 'filename' itself) and where an isolate was isolated from (not the
# sample's matrix)
_GNPS_NOT_EVIDENCE = re.compile(r'(?i)^(?:attribute_|curated_)*(?:sample_?type(?:_?sub1)?|study_sampletypeasrecorded|'
                                r'file_?name|filepath|usi|isolationsource\w*)$')


def gnps_refmat_evidence(df):
    """The cells of a GNPS metadata file that can name a reference material, as one
    reference_materials.EVIDENCE_COLUMN cell per row (columns other than the sample type and file name
    columns; reference_materials skips protocol columns), or None on an error."""
    try:
        cols = [c for c in df.columns if not _GNPS_NOT_EVIDENCE.match(str(c).strip())]
        if not cols:
            return pd.Series('', index=df.index, dtype=object)
        values = df[cols].astype(object).where(df[cols].notna(), '').astype(str)
        out = pd.Series('', index=df.index, dtype=object)
        # only the columns and rows with a candidate word
        hits = {c: values[c].str.contains(reference_materials.CANDIDATE_RE, regex=True) for c in cols}
        # (also the rows whose file name has a candidate word: reference_materials reads the name from 'filename')
        name_hit = df['filename'].astype(str).str.contains(reference_materials.CANDIDATE_RE, regex=True) \
            if 'filename' in df.columns else pd.Series(False, index=df.index)
        if not any(h.any() for h in hits.values()) and not name_hit.any():
            return out
        rows = pd.concat([hits[c] for c in cols if hits[c].any()] + [name_hit], axis=1).any(axis=1)
        # (all columns of these rows: encode_evidence keeps the cells with a candidate word and those that only
        # name a matrix, e.g. UBERONBodyPartName 'blood plasma', the row's own matrix)
        out[rows] = [reference_materials.encode_evidence(zip(cols, row))
                     for row in values.loc[rows, cols].itertuples(index=False, name=None)]
        return out
    except Exception as e:
        print(f"WARNING: no reference-material evidence from the sheet: {type(e).__name__}: {e}")
        return None


def process_filename(filename):
    path_parts = filename.split('/')
    for part in path_parts:
        print(part)
        if part.endswith('.d'):
            return '/'.join(path_parts[:path_parts.index(part)+1])
    return filename


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
    parser.add_argument('--path_to_envo_biome_csv')
    parser.add_argument('--path_to_envo_material_csv')
    parser.add_argument('--path_ncbi_rank_division')
    parser.add_argument('--path_to_doid_csv')
    parser.add_argument('--path_to_hpo_csv', default='none')
    parser.add_argument('--path_to_exposure_csv', default='none')
    parser.add_argument('--path_to_ncbi_names_dmp', default='none',
                        help='NCBI names.dmp; resolves synonyms, merged taxa and strains of NCBITaxonomy names')
    parser.add_argument('--path_to_ncbi_nodes_dmp', default='none',
                        help='NCBI nodes.dmp of the pinned taxdump (with --path_to_ncbi_names_dmp and '
                             '--path_to_ncbi_merged_dmp): an ID|Name whose name names another taxon takes the name '
                             '(except strains of the named species and microbes); without it the ID wins')
    parser.add_argument('--path_to_ncbi_merged_dmp', default='none', help='NCBI merged.dmp of the pinned taxdump')
    parser.add_argument('--path_to_cellosaurus_tsv', default='none',
                        help='Cellosaurus table from build_cellosaurus_table.py; checks CellLine and gives CellLineIndex')
    args = parser.parse_args()

    with open(args.AllowedTermJson_path, 'r') as json_file:
        allowed_terms = json.load(json_file)

    uberon_ontology_table = pd.read_csv(args.path_to_uberon_cl_po_csv, index_col=False)
    uberon_ontology_table = uberon_ontology_table.drop_duplicates(subset=['Label'])

    doid_ontology_table = pd.read_csv(args.path_to_doid_csv, index_col=False)
    # HPO (HPOCommonName) and ECTO/PECO (Exposure) tables; optional
    hpo_ontology_table = pd.read_csv(args.path_to_hpo_csv, dtype=str) if args.path_to_hpo_csv != 'none' else None
    exposure_ontology_table = pd.read_csv(args.path_to_exposure_csv, dtype=str) if args.path_to_exposure_csv != 'none' else None
    from cellosaurus import load_cellosaurus_table
    cellosaurus_table = load_cellosaurus_table(args.path_to_cellosaurus_tsv)
    doid_ontology_table = doid_ontology_table.drop_duplicates(subset=['Label'])

    ENVOEnvironmentBiomeIndex_table = pd.read_csv(args.path_to_envo_biome_csv, index_col=False)
    ENVOEnvironmentBiomeIndex_table = ENVOEnvironmentBiomeIndex_table.drop_duplicates(subset=['Label'])

    ENVOEnvironmentMaterialIndex_table = pd.read_csv(args.path_to_envo_material_csv, index_col=False)
    ENVOEnvironmentMaterialIndex_table = ENVOEnvironmentMaterialIndex_table.drop_duplicates(subset=['Label'])

    NCBIRankDivision_table = pd.read_csv(args.path_ncbi_rank_division, index_col = False)
    NCBIRankDivision_table = NCBIRankDivision_table.drop_duplicates(subset=['TaxonID'])

    # all synonyms (the table above keeps one row per label)
    uberon_synonym_table = pd.read_csv(args.path_to_uberon_cl_po_csv, index_col=False, dtype=str)


    print('Starting tsv processing.')
    #for file_path in glob.glob(f"{args.path_to_github_metadata}/redu_*.tsv"):

    all_metadata_files = sorted(glob.glob(f"{args.path_to_github_metadata}/*.tsv"))

    print("Processing", len(all_metadata_files), "files")

    if 'none' not in (args.path_to_ncbi_nodes_dmp, args.path_to_ncbi_names_dmp, args.path_to_ncbi_merged_dmp):
        import REDU_conversion_functions
        REDU_conversion_functions.use_pinned_ncbi_taxonomy(args.path_to_ncbi_nodes_dmp, args.path_to_ncbi_names_dmp,
                                                           args.path_to_ncbi_merged_dmp)

    ncbi_name_resolver = None
    if args.path_to_ncbi_names_dmp != 'none':
        ncbi_name_resolver = build_ncbi_name_resolver(
            args.path_to_ncbi_names_dmp, _ncbi_names_in_files(all_metadata_files, allowed_terms),
            allowed_terms['NCBITaxonomy']['allowed_values'], NCBIRankDivision_table)
   

    # A file that cannot be read or harmonised is reported and skipped; the others are
    # still harmonised (one bad file used to stop the loop, and the step's '|| true' passed
    # the partial folder on)
    failed_files = []
    written_files = 0
    for file_path in all_metadata_files:
        print(f"Processing: {file_path}")
        partial_name = None
        try:
            # Read as strings so e.g. a YearOfAnalysis column with blanks is not parsed as
            # float ('2019.0' would not match the allowed '2019')
            df = pd.read_csv(file_path, sep='\t', dtype=str)

            if 'ATTRIBUTE_MassiveID' in df.columns:
                df.rename(columns={'ATTRIBUTE_MassiveID': 'MassiveID'}, inplace=True)
            if 'ATTRIBUTE_DatasetAccession' in df.columns:
                df.rename(columns={'ATTRIBUTE_DatasetAccession': 'MassiveID'}, inplace=True)

            # the sheet's own sample fields, for the reference materials (complete_and_fill_REDU_table)
            refmat_evidence = gnps_refmat_evidence(df)
            # GitHub files are named <MSV>.tsv, MassIVE downloads <MSV>_<n>_gnps_metadata.tsv
            accession = re.match(r'(MSV\d+)', os.path.basename(file_path))
            df = normalize_gnps_metadata(df, allowed_terms,
                                         dataset_accession=accession.group(1) if accession else None,
                                         uberon_synonym_table=uberon_synonym_table,
                                         envo_tables=(ENVOEnvironmentBiomeIndex_table, ENVOEnvironmentMaterialIndex_table),
                                         NCBIRankDivision_table=NCBIRankDivision_table)
            # files named as blanks / pooled QCs (BLANK_5_9_B, PoolQC_03) that are labelled as samples
            df, _ = apply_filename_blank_qc(df)
            if refmat_evidence is not None:
                df[reference_materials.EVIDENCE_COLUMN] = refmat_evidence.reindex(df.index)

            #generate extra columns, add missing columns and remove values which are not in allowed terms
            df = complete_and_fill_REDU_table(df, 
                                              allowed_terms, 
                                              UBERONOntologyIndex_table=uberon_ontology_table, 
                                              DOIDOntologyIndex_table=doid_ontology_table, 
                                              HPOOntologyIndex_table=hpo_ontology_table,
                                              ExposureIndex_table=exposure_ontology_table,
                                              NCBIRankDivision_table=NCBIRankDivision_table,
                                              ENVOEnvironmentBiomeIndex_table=ENVOEnvironmentBiomeIndex_table,
                                              ENVOEnvironmentMaterialIndex_table=ENVOEnvironmentMaterialIndex_table,
                                              attempt_adding_file_extensions=True,
                                              case_insensitive_terms=True,
                                              NCBI_name_resolver=ncbi_name_resolver,
                                              ncbi_name_over_id=True,
                                              ncbi_name_fallback=True,
                                              CellLineIndex_table=cellosaurus_table,
                                              life_stage_as_submitted=True)
            
            if len(df) > 0:
                # written under a name the next steps do not read (*.tsv), then renamed, so
                # a failed write leaves no partial file
                file_name = os.path.join(args.output_metadata_folder, os.path.basename(file_path))
                partial_name = os.path.join(args.output_metadata_folder, '.' + os.path.basename(file_path) + '.partial')
                df.to_csv(partial_name, sep='\t', index=False)
                os.replace(partial_name, file_name)
                partial_name = None
                written_files += 1
        except Exception as e:
            failed_files.append(file_path)
            sys.stdout.flush()
            print(f"WARNING: skipping {file_path}, it could not be harmonised: {type(e).__name__}: {e}",
                  file=sys.stderr)
            traceback.print_exc(file=sys.stderr)
            sys.stderr.flush()
            if partial_name is not None:
                try:
                    os.remove(partial_name)
                except OSError:
                    pass  # no partial file, or it cannot be removed; it is not read as *.tsv

    if failed_files:
        # one line on stdout (the step's log) that names the files missing from the output
        print(f"WARNING: harmonised only {len(all_metadata_files) - len(failed_files)} of {len(all_metadata_files)} files "
              f"({written_files} written, {len(failed_files)} failed, not in the output: "
              f"{', '.join(sorted(os.path.basename(p) for p in failed_files))})")
    else:
        print(f"Harmonised {len(all_metadata_files) - len(failed_files)} of {len(all_metadata_files)} files "
              f"({written_files} written, {len(failed_files)} failed)")
    if failed_files:
        print(f"WARNING: {len(failed_files)} of {len(all_metadata_files)} files could not be harmonised "
              f"and are not in the output: {', '.join(os.path.basename(p) for p in failed_files)}", file=sys.stderr)
        if len(failed_files) == len(all_metadata_files):
            sys.exit(1)
