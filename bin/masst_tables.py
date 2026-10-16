"""Readers for the foodMASST, personalCareProductMASST and microbiomeMASST tables
(https://github.com/robinschmid/microbe_masst/tree/master/data), used by MASST_to_REDU.py.

Each reader returns one row per file with MassiveID, filename (path within the dataset, or
the bare file name for tables that only give it) and ReDU columns in submitted form; the
values are validated afterwards by complete_and_fill_REDU_table like any other source.

tissueMASST and plant_masst_table_full are not read: they are ReDU exports.
"""
import io
import json
import os
import re
import collections

import pandas as pd

import diet_columns
from blank_qc_names import blank_qc_type_from_filename, name_tokens
from reference_materials import EVIDENCE_COLUMN as REFMAT_EVIDENCE_COLUMN, encode_evidence as refmat_encode_evidence

_MISSING = {'', 'nan', 'na', 'n/a', 'none', 'null', 'unknown', 'not applicable', 'not_applicable', 'not entered',
            'not collected', 'not provided', 'not specified', 'not_specified', 'not available', 'missing value',
            'missing', '#n/d', 'not_collected', 'other'}

TRANSLATION_SHEETS = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'translation_sheets')
GFOP_SHEET = os.path.join(TRANSLATION_SHEETS, 'GFOP_food_ontology.csv')
COUNTRY_ALIASES = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'translation_sheets_metabolights',
                               'MTBLS_Country_aliases.csv')


def clean(value):
    text = str(value).strip() if value is not None else ''
    return '' if text.lower() in _MISSING else text


def read_table(path, sep, keep_default_na=False):
    """The table as strings, or None (with a warning) when it is missing or empty."""
    if not path or path == 'none' or not os.path.exists(path) or os.path.getsize(path) == 0:
        print(f"WARNING: MASST table {path} is missing or empty; skipping it")
        return None
    try:
        if path.endswith('_full.tsv'):
            df = _read_food_full(path)
        else:
            df = pd.read_csv(path, sep=sep, dtype=str, keep_default_na=keep_default_na)
    except Exception as e:  # a broken download must not stop the other tables
        print(f"WARNING: could not read MASST table {path} ({e}); skipping it")
        return None
    if len(df) == 0:
        print(f"WARNING: MASST table {path} has no rows; skipping it")
        return None
    return df


def _read_food_full(path):
    """food_masst_table_full.tsv: the header line is broken by a newline after 'sample_type'
    (7 + 87 fields); the two header lines are joined when they only add up together to the
    number of fields of the data rows."""
    with open(path, encoding='utf-8', errors='replace') as f:
        text = f.read()
    lines = text.split('\n', 3)
    if len(lines) >= 3:
        n0, n1, n2 = (len(lines[i].split('\t')) for i in range(3))
        if n0 < n2 and n0 + n1 == n2:
            text = lines[0] + '\t' + text.split('\n', 1)[1]
    return pd.read_csv(io.StringIO(text), sep='\t', dtype=str, keep_default_na=False)


def _split_filepath(filepath):
    accession, _, rest = str(filepath).strip().partition('/')
    return accession, rest


# ---------------------------------------------------------------- shared lookups

class Vocab:
    """Exact (case-insensitive) label/synonym lookups against the ReDU vocabularies."""

    def __init__(self, allowed_terms, uberon_table=None, doid_table=None, exposure_table=None, names_dmp=None):
        self.allowed = allowed_terms
        self.countries = {t.lower(): t for t in allowed_terms['Country']['allowed_values']}
        self.country_aliases = {}
        if os.path.exists(COUNTRY_ALIASES):
            a = pd.read_csv(COUNTRY_ALIASES, dtype=str, keep_default_na=False)
            self.country_aliases = {k.strip().lower(): v for k, v in zip(a['alias'], a['Country'])
                                    if v in allowed_terms['Country']['allowed_values']}
        self.uberon, self.po = self._label_synonyms(uberon_table, 'UBERONOntologyIndex',
                                                    allowed_terms['UBERONBodyPartName']['allowed_values'])
        self.doid, _ = self._label_synonyms(doid_table, 'DOIDOntologyIndex', allowed_terms['DOIDCommonName']['allowed_values'])
        # ECTO 'exposure to X' terms whose agent is a chemical (ChEBI)
        self.ecto_chemical = {}
        if exposure_table is not None and len(exposure_table) and 'OrganismExposureAgentIndex' in exposure_table.columns:
            e = exposure_table[exposure_table['ExposureIndex'].astype(str).str.startswith('ECTO:')
                               & exposure_table['OrganismExposureAgentIndex'].astype(str).str.startswith('CHEBI:')]
            allowed_exposure = set(allowed_terms['Exposure']['allowed_values'])
            by = collections.defaultdict(set)
            for label, synonym in zip(e['Label'].astype(str), e['Synonym'].astype(str)):
                if label in allowed_exposure:
                    for name in (label, synonym):
                        by[name.strip().lower()].add(label)
            self.ecto_chemical = {k: next(iter(v)) for k, v in by.items() if len(v) == 1}
        self.common_names = self._common_names(names_dmp) if names_dmp else {}

    @staticmethod
    def _label_synonyms(table, index_column, allowed):
        """lower label/exact synonym -> label (unique only), and the PO subset of labels."""
        if table is None or len(table) == 0:
            return {t.lower(): t for t in allowed}, {}
        allowed = set(allowed)
        by = collections.defaultdict(set)
        po = {}
        for idx, label, synonym in zip(table[index_column].astype(str), table['Label'].astype(str),
                                       table['Synonym'].astype(str) if 'Synonym' in table.columns else table['Label'].astype(str)):
            if label not in allowed:
                continue
            by[label.strip().lower()].add(label)
            if len(synonym.strip()) > 3:
                by[synonym.strip().lower()].add(label)
            if idx.startswith('PO:'):
                po[label.strip().lower()] = label
        # a label always wins over another term's synonym
        labels = {t.lower(): t for t in allowed}
        out = {k: next(iter(v)) for k, v in by.items() if len(v) == 1}
        out.update({k: v for k, v in labels.items()})
        return out, po

    def _common_names(self, names_dmp):
        """NCBI common names ('cow' -> 9913|Bos taurus) that name exactly one taxon."""
        by_id = {}
        for t in self.allowed['NCBITaxonomy']['allowed_values']:
            tid, _, _ = t.partition('|')
            by_id[tid.strip()] = t
        hits = collections.defaultdict(set)
        with open(names_dmp, encoding='utf-8', errors='replace') as f:
            for line in f:
                parts = line.split('\t|\t')
                if len(parts) < 4:
                    continue
                cls = parts[3].rstrip('\t|\n')
                if cls in ('common name', 'genbank common name'):
                    hits[parts[1].strip().lower()].add(parts[0])
        return {k: by_id[next(iter(v))] for k, v in hits.items() if len(v) == 1 and next(iter(v)) in by_id}

    def plant_taxon(self, name, division_by_taxid):
        """'id|name' of the plant/fungus/alga NCBI name (case-insensitive), also for names
        NCBI uses for an animal too (Ficus, Olea, Pisum, Thymus): only the taxon in the
        'Plants and Fungi' division counts; the name itself when it is not unique there."""
        if not name:
            return ''
        if not hasattr(self, '_plant_names'):
            index = collections.defaultdict(set)
            for t in self.allowed['NCBITaxonomy']['allowed_values']:
                tid, _, tname = t.partition('|')
                if division_by_taxid.get(tid.strip()) == 'Plants and Fungi':
                    index[tname.strip().lower()].add(t)
            self._plant_names = index
        hits = self._plant_names.get(name.strip().lower(), set())
        return next(iter(hits)) if len(hits) == 1 else name

    def country(self, value):
        v = clean(value).lower()
        return self.countries.get(v) or self.country_aliases.get(v, '')

    def body_part(self, value):
        v = clean(value)
        if not v:
            return ''
        key = re.sub(r'[_\s]+', ' ', v).strip().lower()
        return self.uberon.get(key, '')

    def plant_part(self, value):
        return self.po.get(clean(value).lower(), '')

    def disease(self, value):
        v = clean(value)
        if not v:
            return ''
        key = re.sub(r'[_\s]+', ' ', v).strip().lower()
        return self.doid.get(key, '')

    def chemical_exposure(self, value):
        v = clean(value).lower()
        if not v or ',' in v or ';' in v:
            return ''
        return self.ecto_chemical.get('exposure to ' + v, '')


def _number(value):
    v = clean(value)
    try:
        x = float(v)
    except ValueError:
        return ''
    if x != x:
        return ''
    return str(int(x)) if x == int(x) else str(x)


# ---------------------------------------------------------------- GFOP (foodMASST tree)

def gfop_terms_from_tree(tree_json_path):
    """GFOP classes of the foodMASST tree: Label, FoodOntologyIndex (WebProtege local ID), IRI, Path."""
    tree = json.load(open(tree_json_path))
    rows = []

    def walk(node, path):
        name = node.get('name')
        iri = node.get('id', '')
        if iri:
            rows.append({'Label': name, 'FoodOntologyIndex': iri.rsplit('/', 1)[-1], 'IRI': iri,
                         'Path': ' > '.join(path + [name])})
        for child in node.get('children', []):
            walk(child, path + ([name] if iri else []))
    walk(tree, [])
    return pd.DataFrame(rows).drop_duplicates('Label')


def load_gfop(path=GFOP_SHEET):
    if not os.path.exists(path):
        return {}
    t = pd.read_csv(path, dtype=str, keep_default_na=False)
    return dict(zip(t['Label'], t['FoodOntologyIndex']))


# ---------------------------------------------------------------- foodMASST

def food_rows(df_csv, df_full, vocab, gfop_labels, division_by_taxid=None):
    """Global FoodOmics (foodMASST) sample sheet -> ReDU rows (one per file)."""
    if df_csv is None or df_full is None:
        return pd.DataFrame()
    t = df_full.merge(df_csv[['Filename', 'MassIVE']].drop_duplicates('Filename'), on='Filename', how='inner')
    st = t['sample_type'].str.strip().str.lower()
    g1 = t['sample_type_group1'].str.strip().str.lower()
    b1 = t['sample_type_groupB1'].str.strip().str.lower().str.replace('-', ' ')
    out = pd.DataFrame({'MassiveID': t['MassIVE'], 'filename': t['Filename']})
    out['SampleType'] = st.where(st.isin(['food', 'beverage']), '')
    food_sub = g1.map({'plant': 'food_source_plant', 'animal': 'food_source_animal', 'fungi': 'food_source_fungi',
                       'complex': 'food_source_complex', 'mineral': 'food_source_NOS', 'water': 'food_source_NOS',
                       'algae': 'food_source_NOS', 'supplement': 'food_source_NOS'}).fillna('')
    bev_sub = b1.map({'alcoholic': 'beverage_alcoholic', 'non alcoholic': 'beverage_nonalcoholic'}).fillna('')
    out['SampleTypeSub1'] = [b if s == 'beverage' else (f if s == 'food' else '') for s, f, b in zip(st, food_sub, bev_sub)]

    # organism: the plant/fungus/alga by its botanical name, an animal source by its NCBI common name
    def organism(group, species, genus, animal):
        if group in ('plant', 'fungi', 'algae'):
            name = re.sub(r'\s+spp?\.?$', '', clean(species) or clean(genus)).strip()
            return vocab.plant_taxon(name, division_by_taxid or {})
        if group == 'animal':
            return vocab.common_names.get(clean(animal).lower(), '')
        return ''
    out['NCBITaxonomy'] = [organism(g, s, ge, a) for g, s, ge, a in
                           zip(g1, t['botanical_genus_species'], t['botanical_genus'], t['animal_source'])]
    out['UBERONBodyPartName'] = [vocab.plant_part(v) if g in ('plant', 'fungi', 'algae') else ''
                                 for g, v in zip(g1, t['botanical_anatomy'])]
    out['Country'] = t['geo_loc_name'].map(lambda v: vocab.country(str(v).split(':')[0]))
    lat, lon = t['latitude'].map(_number), t['longitude'].map(_number)
    out['LatitudeandLongitude'] = [f'{a}|{b}' if a and b else '' for a, b in zip(lat, lon)]
    out['DepthorAltitudeMeters'] = t['elevation'].map(_number)
    out['SampleCollectionDateandTime'] = t['collection_timestamp'].map(clean)
    out['FoodOntologyTerm'] = t['ontology_terminal_leaf'].map(lambda v: v if v in gfop_labels else '')
    out['_food_group'] = g1.values
    return out


# ---------------------------------------------------------------- personalCareProductMASST

def personal_care_rows(df):
    """personalCareProduct_masst_table -> ReDU rows. Products (personal care and the few
    household cleaners) are SampleType personal care product (SampleTypeSub1
    reference material_personalcareproduct) with the last, most specific category of the
    product tree ('Root > Personal Care > Cosmetics > Skin Care > Sunscreen' -> 'Sunscreen');
    Root > Blank > Solvent/Method Blank are analysis/extraction blanks."""
    if df is None:
        return pd.DataFrame()
    out = pd.DataFrame({'MassiveID': df['MassIVE'], 'filename': df['Filename']})
    path = df['Category'].astype(str).str.strip()
    path = path.str.replace(r'^(Root|Health & Beauty)\s*>\s*', '', regex=True)
    path = path.map(lambda p: ' > '.join(x.strip() for x in p.split('>') if x.strip()))
    blank = path.str.match(r'(?i)^blank(\s|$)')
    method = blank & path.str.contains(r'(?i)method|extraction|procedur')
    out['SampleType'] = 'personal care product'
    out['SampleTypeSub1'] = 'reference material_personalcareproduct'
    out.loc[blank, ['SampleType', 'SampleTypeSub1']] = 'blank_analysis'
    out.loc[method, ['SampleType', 'SampleTypeSub1']] = 'blank_extraction'
    category = path.str.split(' > ').str[-1]
    # the same category name under several paths of the tree (none in the pinned table) is one value
    paths = path[~blank].groupby(category[~blank]).unique()
    for name, ps in paths[paths.map(len) > 1].items():
        print(f"personalCareProductMASST: category {name!r} is under {len(ps)} paths: {sorted(ps)}")
    out['ProductCategory'] = category.where(~blank, '')
    return out


# ---------------------------------------------------------------- microbiomeMASST

_MB_BLANK_TYPES = {'blank', 'qc_blank_monoculture', 'qc_blank_community', 'qc_blank_bsh', 'skin_blank',
                   'human cell lines qc blank'}
_MB_QC_TYPES = {'qc_pool_monoculture', 'qc_pool_community', 'qc_pool_bsh', 'skin_pool', 'colon_simulator_qc_pool',
                'human cell lines qc pool'}
_MB_MEDIA_TYPES = {'media_control_community', 'media_control_monoculture', 'skin_media_control', 'colon_simulator_medium'}
_MB_REFERENCE_TYPES = {'food reference material', 'human reference material'}
# kept on the files of pooled QCs and reference materials (the other subject fields are cleared)
MB_KEEP_ON_QC_AND_REFERENCE = ('NCBITaxonomy', 'UBERONBodyPartName')
# Taxaname_file values that are tree categories, not organisms
MB_CATEGORIES = {'bacteria', 'space', 'animal', 'food', 'human cell lines', 'protein screening', 'environment'}
_TAXON = re.compile(r'^\s*(\d+)\s*\|\s*(\S.*)$')
_ANIMAL_DIVISIONS = {'Primates', 'Rodents', 'Mammals', 'Vertebrates', 'Invertebrates'}
_CONTROL_TOKEN = re.compile(r'(?i)^(?:ctrl|control|controls|veh|vehicle|vehcontrol|sham|placebo|untreated|mock|naive)\d*$')
_COLONIZATION = {'germ_free': 'germ-free', 'germ-free': 'germ-free', 'germ free': 'germ-free', 'gf': 'germ-free',
                 'specific pathogen free': 'specific-pathogen-free', 'specific-pathogen-free': 'specific-pathogen-free',
                 'spf': 'specific-pathogen-free', 'conventional': 'conventional', 'conventionally raised': 'conventional'}
_HEALTH = {'healthy': 'healthy', 'unhealthy (nos)': 'unhealthy (NOS)', 'acute illness': 'acute illness',
           'chronic illness': 'chronic illness'}
_ECHO_SUBSPECIFIC = {'', 'not available', 'diet', 'human diet', 'gnotobiotic'}


def _row_kind(itype, blank_col, qc_col, community):
    t = itype.strip().lower()
    if t in _MB_MEDIA_TYPES or ('media' in t and 'control' in t):
        return 'media'
    if t in _MB_REFERENCE_TYPES:
        return 'reference'
    blank = t in _MB_BLANK_TYPES or 'blank' in t or blank_col.strip().lower() == 'yes' \
        or community.lower().startswith('qc_blank')
    qc = t in _MB_QC_TYPES or ('qc' in t and 'pool' in t) or qc_col.strip().lower() == 'yes' \
        or community.lower().startswith('qc_pool')
    if blank and qc:
        return 'blank_QC'
    if blank:
        return 'blank_analysis'
    if qc:
        return 'pool_QC'
    return ''


def microbiome_rows(df, vocab, division_by_taxid, microbe_sample_type=None, stats=None):
    """microbiomeMASST -> ReDU rows (one per file). Only rows not copied from ReDU
    (imported_from_redu != yes). The rows of one file (tree nodes / interventions) are
    merged: per-file fields that disagree are left empty; interventions are collected."""
    stats = stats if stats is not None else collections.Counter()
    if df is None:
        return pd.DataFrame()
    df = df[df['imported_from_redu'].str.strip().str.lower() != 'yes'].copy()
    stats['microbiome_rows_curated'] = len(df)
    df = df[df['Filepath'].str.contains('/', na=False)]
    acc_path = df['Filepath'].map(_split_filepath)
    df['MassiveID'] = [a for a, _ in acc_path]
    df['filename'] = [p for _, p in acc_path]
    df['_kind'] = [_row_kind(a, b, c, d) for a, b, c, d in
                   zip(df['Intervention_type'], df['Blank'], df['QC'], df['microbe_community_type'])]

    files = []
    for (mid, fname), g in df.groupby(['MassiveID', 'filename'], sort=False):
        rec = {'MassiveID': mid, 'filename': fname}
        conflicts = []

        def one(values, column):
            vals = list(dict.fromkeys(v for v in values if v))
            if len(vals) > 1:
                conflicts.append(column)
                return ''
            return vals[0] if vals else ''

        kinds = set(g['_kind']) - {''}
        name_type = blank_qc_type_from_filename(fname)
        key_type = None
        for key in g['SampleType_key']:
            key_type = key_type or blank_qc_type_from_filename(str(key).replace(' ', '_'))
        kind = ''
        if kinds & {'blank_analysis', 'blank_QC', 'pool_QC'} or name_type or key_type:
            table_kinds = kinds & {'blank_analysis', 'blank_QC', 'pool_QC'}
            if len(table_kinds) == 1:
                kind = next(iter(table_kinds))
            elif table_kinds:
                kind = 'blank_QC'
            else:
                kind = name_type or key_type
            # an extraction blank by name is more specific than the table's 'blank'
            if kind == 'blank_analysis' and name_type == 'blank_extraction':
                kind = 'blank_extraction'
            stats['microbiome_files_blank_qc'] += 1
            stats[f'microbiome_files_blank_qc_from_{"table" if table_kinds else ("name" if name_type else "sample_key")}'] += 1
        elif 'media' in kinds:
            kind = 'blank_culturemedia'
        elif 'reference' in kinds:
            kind = 'reference material'
        rec['_kind'] = kind

        # host / organism
        hosts = []
        category = ''
        for tax, tax_ncbi, itype in zip(g['Taxaname_file'], g['Taxa_NCBI'], g['Intervention_type']):
            tax = tax.strip()
            m = _TAXON.match(tax)
            if m:
                hosts.append(f'{m.group(1)}|{m.group(2).strip()}')
            elif tax.lower() in MB_CATEGORIES:
                category = category or tax.lower()
                if tax.lower() == 'bacteria' and _number(tax_ncbi):
                    hosts.append(_number(tax_ncbi) + '|')
                elif tax.lower() == 'animal' and _TAXON.match(itype):
                    m2 = _TAXON.match(itype)
                    hosts.append(f'{m2.group(1)}|{m2.group(2).strip()}')
                elif tax.lower() == 'human cell lines':
                    hosts.append('9606|Homo sapiens')
        host = one(hosts, 'NCBITaxonomy')
        rec['NCBITaxonomy'] = host
        taxid = host.split('|')[0] if host else ''

        if kind:
            sample_type = kind
        elif category == 'human cell lines':
            sample_type = 'culture_mammalian'
        elif category == 'bacteria':
            community = ' '.join(g['microbe_community_type']) + ' ' + ' '.join(g['Intervention_type']) + ' ' + fname
            if taxid and microbe_sample_type is not None:
                sample_type = microbe_sample_type(taxid)[0] or 'culture_bacterial'
            elif re.search(r'(?i)comm|syncom|consorti|mix', community):
                sample_type = 'culture_multiplespecies'
            else:
                sample_type = 'culture_bacterial'
        elif category == 'space':
            sample_type = 'built_environment'
        elif category == 'environment':
            sample_type = 'environmental'
        elif category == 'food':
            sample_type = 'food'
        elif category == 'animal' or division_by_taxid.get(taxid) in _ANIMAL_DIVISIONS:
            sample_type = 'animal'
        else:
            sample_type = ''
        rec['SampleType'] = sample_type
        rec['SampleTypeSub1'] = sample_type if sample_type in ('blank_analysis', 'blank_extraction', 'blank_QC', 'pool_QC',
                                                                'blank_culturemedia', 'culture_bacterial', 'culture_fungal',
                                                                'culture_mammalian', 'culture_multiplespecies') else ''

        rec['UBERONBodyPartName'] = one([vocab.body_part(v) for v in g['UBERONBodyPartName']], 'UBERONBodyPartName')
        doid, health = [], []
        for v in g['DOIDCommonName']:
            v = clean(v)
            if v.lower() in _HEALTH:
                health.append(_HEALTH[v.lower()])
            elif v and v.lower() != 'blank_qc':
                doid.append(vocab.disease(v) or v)
        rec['DOIDCommonName'] = one(doid, 'DOIDCommonName')
        rec['HealthStatus'] = one(health, 'HealthStatus')
        # a file listed as healthy (under another tree node) carries no disease
        if rec['HealthStatus'] == 'healthy' and rec['DOIDCommonName']:
            stats['microbiome_files_healthy_disease_dropped'] += 1
            rec['DOIDCommonName'] = ''
        rec['BiologicalSex'] = one([clean(v).lower() for v in g['BiologicalSex']], 'BiologicalSex')
        rec['AgeInYears'] = one([_number(v) for v in g['AgeInYears_MSV000092757']], 'AgeInYears')
        rec['TimePointHoursRelativeToEvent'] = one([_number(v) for v in g['Timepoint_h']], 'TimePointHoursRelativeToEvent')

        # interventions (graph layers): colonization status, chemical exposures, other groups
        interventions, colonization, raw = [], [], []
        for itype, spec, sub, rk in zip(g['Intervention_type'], g['Intervention_specific'],
                                        g['Intervention_subspecific'], g['_kind']):
            itype, spec, sub = itype.strip(), clean(spec), clean(sub)
            if rk or not clean(itype) or _TAXON.match(itype):
                continue
            if itype.lower() == 'gnotobiotic':
                status = _COLONIZATION.get(spec.lower()) or ('gnotobiotic' if spec else '')
                if status:
                    colonization.append(status)
                    continue
            raw.append((itype, spec, sub))
        # Intervention_subspecific is given per file (the same on all its rows, e.g. the
        # genotype WT/HET/KO): it goes to the genetic intervention if there is one
        has_genetic = any('genetic' in t.lower() for t, _, _ in raw)
        for itype, spec, sub in raw:
            value = spec
            attach = 'genetic' in itype.lower() or not has_genetic
            if attach and sub.lower() not in _ECHO_SUBSPECIFIC and sub.lower() != itype.lower() \
                    and sub.lower() not in spec.lower():
                value = f'{spec}, {sub}' if spec else sub
            interventions.append((itype, value, spec if spec else sub))
        rec['ColonizationStatus'] = one(colonization, 'ColonizationStatus')
        if kind:
            # blanks, QCs, media and reference material carry no subject metadata; pooled QCs and
            # reference materials keep the organism and body part (complete_and_fill_REDU_table: a
            # reference material then takes those of the material, reference_materials.py)
            for c in ('NCBITaxonomy', 'UBERONBodyPartName', 'DOIDCommonName', 'HealthStatus', 'BiologicalSex',
                      'AgeInYears', 'TimePointHoursRelativeToEvent', 'ColonizationStatus'):
                if c in MB_KEEP_ON_QC_AND_REFERENCE and kind in ('pool_QC', 'reference material'):
                    continue
                if rec.get(c):
                    stats[f'microbiome_blank_qc_stripped_{c}'] += 1
                rec[c] = ''
        # the intervention types (e.g. 'human reference material'), for reference_materials.py
        # (and, where they name one, the file's body part, which encode_evidence keeps if it only names a matrix: the
        # material's matrix)
        _rm = refmat_encode_evidence(
            ('Intervention_type', t) for t in dict.fromkeys(str(x).strip() for x in g['Intervention_type']) if t)
        if _rm and rec.get('UBERONBodyPartName'):
            _rm = refmat_encode_evidence([('Intervention_type', t) for t in dict.fromkeys(
                str(x).strip() for x in g['Intervention_type']) if t] + [('UBERONBodyPartName', rec['UBERONBodyPartName'])])
        rec[REFMAT_EVIDENCE_COLUMN] = _rm
        rec['_interventions'] = list(dict.fromkeys(interventions))
        rec['_control_name'] = any(_CONTROL_TOKEN.match(t) for t in name_tokens(fname))
        for c in conflicts:
            stats[f'microbiome_conflict_blanked_{c}'] += 1
        if len(g) > 1:
            stats['microbiome_files_with_several_rows'] += 1
        files.append(rec)
    out = pd.DataFrame(files)
    stats['microbiome_files'] = len(out)
    if len(out) == 0:
        return out
    out = _assign_interventions(out, vocab, stats)
    return out


def _intervention_value(k):
    """an intervention by what was given (its value, the type where it has none), any type"""
    return ' '.join((k[1] or k[0]).lower().split())


def _is_diet(k):
    """a diet intervention ('mice diet', 'human diet'; not the 'DrugsNDiet' communities)"""
    return diet_columns.masst_diet_type(k[0])


def _assign_interventions(out, vocab, stats):
    """Exposure (ECTO 'exposure to <chemical>' exact match), Diet (diet interventions, the
    mapping of diet_columns.py) or SampleGroupFreeText ('Intervention: <type>: <value>') from
    the interventions of each file.
    - an intervention given for every real (non blank/QC) file of its dataset is a
      property of the study, not of a group of files (e.g. the drug of a study whose
      vehicle controls are listed under the drug too): not used. Compared by value under
      any type (MSV000083777 lists its antibiotics under 'mice bile acids' and 'mice drug');
    - files named as controls (control, vehicle, sham, placebo, untreated, ...) keep only
      genetic and diet interventions (their genotype and diet: the vehicles of a diet study,
      MSV000089518);
    - a diet intervention also gives the file's Diet (animal files only); not one given for every
      real file, nor a Diet that every real file would get under different wordings ('high-fat,
      dark' / 'high-fat, light'; diet_columns.study_level_diet): microbiomeMASST lists a study's
      diet for its control files too (MSV000085879, MSV000097318); its text stays, the final
      table's pass (final_text_columns.py) removes it where Diet holds it on every file;
    - an intervention type that is an Exposure for some real files of the dataset and not
      for others keeps its text on all of them, as a Workbench / MetaboLights factor that the
      ReDU columns tell only partly apart (MSV000084783 'mice diet': iron / iron-deficient /
      normal chow diet);
    - SampleGroupFreeText only where the dataset's real files form >= 2 groups with on
      average >= 2 files per group."""
    gfop = diet_columns.gfop_index(vocab.allowed['FoodOntologyTerm']['allowed_values'])
    diet_terms = set(vocab.allowed['Diet']['allowed_values'])
    real = out['_kind'] == ''
    out['Exposure'] = ''
    out['Diet'] = ''
    out['SampleGroupFreeText'] = ''
    for mid, idx in out.groupby('MassiveID').groups.items():
        sub = out.loc[idx]
        real_idx = sub.index[real[sub.index]]
        n_real = len(real_idx)
        counts = collections.Counter(v for ivs in sub.loc[real_idx, '_interventions'] for v in {_intervention_value(k) for k in ivs})
        study_level = {k for k, n in counts.items() if n >= n_real}
        stats['microbiome_interventions_study_level'] += len(study_level)
        files, diets = {}, {}
        for i in real_idx:
            ivs = [k for k in out.at[i, '_interventions'] if _intervention_value(k) not in study_level]
            if out.at[i, '_control_name']:
                kept = [k for k in ivs if 'genetic' in k[0].lower() or _is_diet(k)]
                if len(kept) < len(ivs):
                    stats['microbiome_files_control_name_interventions_dropped'] += 1
                ivs = kept
            files[i] = (ivs, [vocab.chemical_exposure(k[2]) for k in ivs])
            diets[i], _ = diet_columns.row_diet([diet_columns.map_value('Intervention', k[1], gfop, masst_type=k[0])
                                                 for k in ivs if _is_diet(k)],
                                                human=str(out.at[i, 'NCBITaxonomy']).startswith('9606|'))
        if diet_columns.study_level_diet(diets.values()):
            stats['microbiome_datasets_study_level_Diet_not_used'] += 1
            diets = {}
        for i, diet in diets.items():
            if diet in diet_terms and out.at[i, 'SampleType'] == 'animal':
                out.at[i, 'Diet'] = diet
                stats['microbiome_files_Diet'] += 1
        is_exposure = collections.defaultdict(set)
        for ivs, exposures in files.values():
            for k, e in zip(ivs, exposures):
                is_exposure[k[0]].add(bool(e))
        partly = {t for t, s in is_exposure.items() if len(s) == 2}
        groups = {}
        for i, (ivs, exposures) in files.items():
            matched = [e for e in exposures if e]
            if len(set(matched)) == 1:
                out.at[i, 'Exposure'] = matched[0]
            elif len(set(matched)) > 1:
                stats['microbiome_conflict_blanked_Exposure'] += 1
            rest = [k for k, e in zip(ivs, exposures) if not e or k[0] in partly]
            text = '; '.join(sorted(f'Intervention: {t}: {v}' if v else f'Intervention: {t}' for t, v, _ in rest))
            if text:
                groups[i] = text
        n_groups = len(set(groups.values()))
        if n_groups >= 2 and len(groups) / n_groups >= 2:
            for i, text in groups.items():
                out.at[i, 'SampleGroupFreeText'] = text
            stats['microbiome_datasets_sample_groups_filled'] += 1
        elif groups:
            stats['microbiome_datasets_sample_groups_not_filled'] += 1
    return out.drop(columns=['_interventions', '_control_name'])
