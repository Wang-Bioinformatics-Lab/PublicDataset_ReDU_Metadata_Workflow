"""Rules applied to every row of the final table, after the NCBI lineage (enrich_ncbi_information) and
before the text columns (final_text_columns.py), so also to the rows saveOlderData keeps from earlier runs:
- collapse_cell_breaks: line breaks and tabs inside a cell (multi-line cells of a sheet) become one
  space, so every row of the TSV is one line (line-based readers: awk, cut, split by line)
- personal_care_products: personal care products are SampleType 'personal care product' with the most
  specific product category
- life_stage_humans_only: LifeStage (human age brackets) only on rows of humans
- sex_from_body_part: BiologicalSex of mammal samples with no recorded sex from a body part only one
  sex has (UBERON_sex_specific_parts.csv of prepare_ontologies.py)
- reference_materials_pass: reference materials and in-house pools named in the row's text
  (reference_materials.py; for the rows of earlier runs, and the GNPS rows of files whose MASST rows name
  one: --masst_refmat_evidence of MASST_to_REDU.py)
- pool_qc_study_values: a pooled QC without organism / body part takes the one all real samples of its
  dataset share
- exudate_samples: SampleTypeSub1 'exudate' on the samples whose body part is UBERON 'exudate' (any organism)
- plant_body_parts: Plant Ontology terms for the body parts of plants (the PO term of the same name, 'ear'
  of grasses), other UBERON terms cleared on plant rows, UBERON 'exudate' only on animal rows, PO terms
  cleared on animal rows (UBERON_CL_PO_ontology.csv of prepare_ontologies.py); runs after the pooled-QC
  passes, so a pooled QC of plants takes the same PO term as the samples (and loses an animal term)
- water_kind_from_material: SampleTypeSub1 water_* of environmental water samples whose kind of water is
  stated only in ENVOEnvironmentMaterial ('river water' -> water_surface)
- placeholder_factor_names: SampleGroupFreeText parts whose factor name is a template placeholder
  ('<enter factor name>: no drug') without the name
Values are only filled or renamed, never removed (but LifeStage on rows not of humans, the subject fields
of rows that become reference materials or pooled QCs, and the body parts plant_body_parts clears); a
second pass changes nothing.

    python final_table_rules.py merged_with_ncbi.tsv merged_with_rules.tsv \\
        --path_to_sex_specific_parts UBERON_sex_specific_parts.csv --masst_refmat_evidence masst_refmat_evidence.tsv \\
        --path_to_uberon_cl_po_csv UBERON_CL_PO_ontology.csv
"""
import argparse
import collections
import re

import pandas as pd

from ncbi_database import write_tsv_atomic
import reference_materials
# (the food track's marker of a digest of a food: in vitro digestion, colonic fermentation)
from REDU_conversion_functions import _FOOD_DIGEST_RE

MISSING = 'missing value'
# SampleTypeSub1 of personal care products (GNPS template, personalCareProductMASST)
PERSONAL_CARE_SUB1 = 'reference material_personalcareproduct'
# LifeStage (humans only) of a subject that does not lactate: the row of a milk sample then
# describes the infant who drank it
CHILD_LIFE_STAGES = {'Infancy (<2 yrs)', 'Early Childhood (2 yrs < x <=8 yrs)'}
# NCBI taxon of Homo sapiens (NCBITaxonomy is '<taxon id>|<name>')
HUMAN_TAXON_ID = '9606'


def _missing(series):
    return series.fillna('').astype(str).str.strip().str.lower().isin(('', 'nan', MISSING))


# a run of line breaks / tabs (with the spaces around it) inside a cell
_CELL_BREAK = re.compile(r' *[\r\n\t][\r\n\t ]*')
# columns that name a file: a USI or file name is kept as it is (it must match the file)
FILE_COLUMNS = ('USI', 'filename')


def collapse_cell_breaks(df):
    """\\r, \\n and \\t inside a text cell (sheets with multi-line cells, e.g. a subject identifier with
    a line feed at its end) become one space, and the cell is stripped; a cell left empty becomes
    'missing value'. Cells without these characters are not changed, so a second pass changes nothing.
    USI and filename are kept as they are (a file's name must match the file), their cells with such
    characters are only counted. Returns the frame and a Counter."""
    stats = collections.Counter()
    dataset = df['ATTRIBUTE_DatasetAccession'] if 'ATTRIBUTE_DatasetAccession' in df.columns else pd.Series('', index=df.index)
    changed = pd.Series(False, index=df.index)
    for column in df.columns:
        values = df[column]
        if values.dtype != object:
            continue
        # checked per distinct value (most columns have few)
        broken = [v for v in pd.unique(values) if isinstance(v, str) and _CELL_BREAK.search(v)]
        if not broken:
            continue
        hit = values.isin(broken)
        if column in FILE_COLUMNS:
            stats[f'line_breaks_in_file_column_kept {column}'] = int(hit.sum())
            continue
        fixed = {v: (_CELL_BREAK.sub(' ', v).strip() or MISSING) for v in broken}
        df.loc[hit, column] = values[hit].map(fixed)
        stats[f'line_breaks_collapsed {column}'] = int(hit.sum())
        changed |= hit
    for d, n in dataset[changed].value_counts().sort_index().items():
        stats[f'line_breaks_collapsed_rows {d}'] = int(n)
    stats['line_breaks_collapsed_rows'] = int(changed.sum())
    return df, stats


def personal_care_products(df):
    """A personal care product is not a reference material: SampleType 'reference material' with
    SampleTypeSub1 reference material_personalcareproduct (GNPS metadata files, MASST rows of earlier
    runs) becomes 'personal care product'; a ProductCategory path ('Personal Care > Cosmetics > Skin
    Care > Sunscreen', earlier runs) keeps its last, most specific category. Returns the frame and a
    Counter."""
    stats = collections.Counter()
    if 'SampleType' in df.columns and 'SampleTypeSub1' in df.columns:
        product = (df['SampleType'] == 'reference material') & (df['SampleTypeSub1'] == PERSONAL_CARE_SUB1)
        df.loc[product, 'SampleType'] = 'personal care product'
        stats['SampleType_personal_care_product'] = int(product.sum())
    if 'ProductCategory' in df.columns:
        path = df['ProductCategory'].fillna('').astype(str)
        is_path = path.str.contains('>', regex=False)
        df.loc[is_path, 'ProductCategory'] = path[is_path].map(
            lambda p: ([x.strip() for x in p.split('>') if x.strip()] or [MISSING])[-1])
        stats['ProductCategory_path_to_last_category'] = int(is_path.sum())
    return df, stats


def life_stage_humans_only(df):
    """LifeStage holds human age brackets ('Infancy (<2 yrs)', ...): 'missing value' on every row
    whose organism is not Homo sapiens, i.e. NCBITaxonomy not taxon 9606 ('9606|Homo sapiens') and
    NCBISpecies (the lineage of ncbi_database.py) not Homo sapiens (its subspecies, e.g.
    '63221|Homo sapiens neanderthalensis'); a row without an organism is not taken for human. GNPS
    sheets give the submitter's LifeStage for any organism (mice, cats, fish, bacteria, offices),
    Workbench a stated one ('newborn', 'adolescent') for any organism. Run before
    sex_from_body_part, which reads LifeStage (a second pass then changes nothing). Returns the
    frame and a Counter."""
    stats = collections.Counter()
    if 'LifeStage' not in df.columns or 'NCBITaxonomy' not in df.columns:
        return df, stats
    human = df['NCBITaxonomy'].fillna('').astype(str).str.split('|', n=1).str[0].str.strip() == HUMAN_TAXON_ID
    if 'NCBISpecies' in df.columns:
        human |= df['NCBISpecies'] == 'Homo sapiens'
    clear = ~human & (df['LifeStage'] != MISSING)
    source = df['DataSource'] if 'DataSource' in df.columns else pd.Series('', index=df.index)
    for s, n in source[clear].value_counts().sort_index().items():
        stats[f'LifeStage_not_human_cleared {s}'] = int(n)
    stats['LifeStage_not_human_cleared'] = int(clear.sum())
    df.loc[clear, 'LifeStage'] = MISSING
    return df, stats


def sex_from_body_part(df, parts):
    """Fills a missing BiologicalSex with the sex of the body part (parts: prepare_ontologies.py's
    UBERON_sex_specific_parts.csv; UBERONBodyPartName = its Label), only on rows of
    - mammals (NCBIClass Mammalia; other animals include hermaphrodites, e.g. oysters and many fish)
      sampled as SampleType 'animal' (not blanks, QCs, standards, foods, cultures), without a
      FoodOntologyTerm (milk analysed as food) or CellLine;
    - milk: not where the row's subject is a child (LifeStage), whose sex the row then gives;
    - not on digests (in vitro digestion, colonic fermentation: _FOOD_DIGEST_RE of
      REDU_conversion_functions.py in SampleGroupFreeText), nor on the rows of the same body part of a
      dataset that has such rows (MTBLS6768: in vitro digests of cow-milk infant formula, one of them a
      file listed for two samples, whose text the merge left empty): no animal was sampled;
    - not in a dataset whose rows of the body part record another sex (the column then describes
      someone else, e.g. the infants of a breast milk study).
    Recorded values are never changed: rows whose recorded sex differs from the body part's are only
    counted. Returns the frame and a Counter."""
    stats = collections.Counter()
    if len(df) == 0 or any(c not in df.columns for c in ('BiologicalSex', 'UBERONBodyPartName', 'NCBIClass', 'SampleType')):
        return df, stats
    by_label = parts.groupby('Label')['BiologicalSex'].agg(lambda s: s.iloc[0] if s.nunique() == 1 else '')
    milk_labels = set(parts.loc[parts['Is Milk'].astype(str).str.lower() == 'true', 'Label'])
    part = df['UBERONBodyPartName'].fillna('').astype(str)
    sex = part.map(by_label).fillna('')
    recorded = df['BiologicalSex'].fillna('').astype(str).str.strip().str.lower()
    unrecorded = _missing(df['BiologicalSex'])
    dataset = df['ATTRIBUTE_DatasetAccession'] if 'ATTRIBUTE_DatasetAccession' in df.columns else pd.Series('', index=df.index)
    other_sex = (sex != '') & ~unrecorded & (recorded != sex)
    for p, n in part[other_sex].value_counts().sort_index().items():
        stats[f'recorded_sex_differs_from_body_part_(kept) {p}'] = int(n)
    conflicting = set(zip(dataset[other_sex], part[other_sex]))
    digest = df['SampleGroupFreeText'].fillna('').astype(str).str.contains(_FOOD_DIGEST_RE) \
        if 'SampleGroupFreeText' in df.columns else pd.Series(False, index=df.index)
    digested = set(zip(dataset[digest], part[digest]))

    candidate = (sex != '') & unrecorded
    reasons = [
        ('not_a_mammal', df['NCBIClass'] != 'Mammalia'),
        ('not_SampleType_animal', df['SampleType'] != 'animal'),
        ('food', ~_missing(df['FoodOntologyTerm']) if 'FoodOntologyTerm' in df.columns else None),
        ('cell_line', ~_missing(df['CellLine']) if 'CellLine' in df.columns else None),
        ('milk_of_a_child_subject', part.isin(milk_labels) & df['LifeStage'].isin(CHILD_LIFE_STAGES)
         if 'LifeStage' in df.columns else None),
        ('dataset_records_another_sex', pd.Series([k in conflicting for k in zip(dataset, part)], index=df.index)
         if conflicting else None),
        ('digest', digest | pd.Series([k in digested for k in zip(dataset, part)], index=df.index) if digested else None),
    ]
    fill = candidate.copy()
    for reason, mask in reasons:
        if mask is None:
            continue
        stats[f'not_filled_{reason}'] = int((fill & mask).sum())
        fill &= ~mask
    df.loc[fill, 'BiologicalSex'] = sex[fill]
    for (p, s), n in pd.DataFrame({'p': part[fill], 's': sex[fill]}).value_counts().sort_index().items():
        stats[f'filled {s} {p}'] = int(n)
    stats['filled_rows'] = int(fill.sum())
    return df, stats


# ---- reference materials and pooled QCs (reference_materials.py) ----
# columns made from another column (ncbi_database.py, complete_and_fill_REDU_table): when the source column
# changes, they take the values other rows of the table with the same value have
DERIVED = {
    'NCBITaxonomy': ['NCBIRank', 'NCBIDivision', 'NCBISuperkingdom', 'NCBIKingdom', 'NCBIPhylum', 'NCBIClass', 'NCBIOrder',
                     'NCBIFamily', 'NCBIGenus', 'NCBISpecies', 'OpenTreeOfLifeTaxonomyID'],
    'UBERONBodyPartName': ['UBERONOntologyIndex'],
}
# subject fields of rows that are no samples of the study's subjects (complete_and_fill_REDU_table's
# STRIP_FOR_NOT_SUBJECT) and the columns made from them
NOT_SUBJECT_COLUMNS = ['AgeInYears', 'UniqueSubjectID', 'SubjectIdentifierAsRecorded', 'DOIDCommonName', 'DOIDOntologyIndex',
                       'HPOCommonName', 'HPOOntologyIndex', 'Exposure', 'ExposureIndex', 'OrganismExposureAgent',
                       'OrganismExposureAgentIndex', 'TimePointHoursRelativeToEvent', 'FastingStatus', 'ENVOEnvironmentBiome',
                       'ENVOEnvironmentBiomeIndex', 'ENVOEnvironmentMaterial', 'ENVOEnvironmentMaterialIndex', 'ENVOBroadScale',
                       'ENVOBroadScaleIndex', 'ENVOLocalScale', 'ENVOLocalScaleIndex', 'ENVOMediumScale', 'ENVOMediumScaleIndex',
                       'BiologicalSex', 'HealthStatus', 'LatitudeandLongitude', 'SmokingStatus', 'BodyMassIndex', 'Diet',
                       'CellLine', 'CellLineIndex', 'ColonizationStatus', 'LifeStage']


def _derived_lookup(df, column):
    """value of `column` -> the most common values of its DERIVED columns among the rows that have it
    (ties: the first in sorted order); only values whose rows all agree are used"""
    cols = [c for c in DERIVED[column] if c in df.columns]
    if not cols:
        return {}, cols
    sub = df[[column] + cols].astype(str)
    counts = sub.value_counts(sort=False).reset_index(name='_n')
    counts = counts.sort_values(['_n'] + [column] + cols, ascending=[False] + [True] * (len(cols) + 1), kind='stable')
    first = counts.drop_duplicates(column)
    return {r[0]: tuple(r[1:1 + len(cols)]) for r in first.itertuples(index=False, name=None)}, cols


def _set_with_derived(df, rows, column, values, lookups, stats, label):
    """sets `column` on `rows` to `values` (a Series by index) and its DERIVED columns as lookups[column]
    gives them; a value without a lookup entry (no other row has it) is not set"""
    lookup, cols = lookups[column]
    ok = [i for i in rows if values[i] == MISSING and not cols or values[i] in lookup]
    for i in ok:
        df.at[i, column] = values[i]
        for c, v in zip(cols, lookup.get(values[i], ())):
            df.at[i, c] = v
    if len(ok) < len(rows):
        stats[f'{label}_not_set_value_unknown_to_the_table'] += len(rows) - len(ok)
    return ok


# solvent injections by their file name: a solvent word and nothing else but run words ('20230426-MeOH-1.mzML',
# 'ACN_2.raw', 'H2O_pos_3.mzML'); not samples named with their extraction solvent ('005_Sa04_Water_POS.mzML',
# 'Sa12_Methanol_NEG.mzML': sample Sa04 extracted with water)
_SOLVENT_FILE_WORDS = {'meoh', 'methanol', 'acn', 'acetonitrile', 'ipa', 'isopropanol', 'solvent', 'solvents', 'water'}
_RUN_WORDS = {'pos', 'neg', 'positive', 'negative', 'ms', 'msms', 'dda', 'dia', 'full', 'scan', 'run', 'inj', 'injection',
              'rep', 'replicate', 'only', 'pure'}


def solvent_file_name(filename):
    """whether a file's base name says it is a solvent injection (see _SOLVENT_FILE_WORDS)"""
    base = str(filename).rstrip('/').rsplit('/', 1)[-1].lower()
    base = re.sub(r'(\.(mzml|mzxml|mzdata|mgf|raw|cdf|wiff2?|d|gz|zip|lzma|bz2|7z|xz|rar|scan))+$', '', base)
    base = re.sub(r'h2o(?!\d)', ' water ', base)
    runs = set(re.findall(r'[a-z]+', base))
    return bool(runs & _SOLVENT_FILE_WORDS) and runs <= _SOLVENT_FILE_WORDS | _RUN_WORDS


def reference_materials_pass(df, masst_evidence=None):
    """reference_materials.annotate on the final table (dataset ATTRIBUTE_DatasetAccession; the text of
    filename, SubjectIdentifierAsRecorded and SampleGroupFreeText, plus masst_evidence: a {USI: evidence}
    of MASST rows of the GNPS-described files). Rows already labelled as specifically keep their label, so
    this changes the rows of earlier runs and the GNPS rows of files whose MASST rows name a reference
    material. A row that becomes a reference material takes the material's organism and body part (and
    the lineage / index columns other rows with these values have); one that becomes a reference material
    or an in-house pool (pool_QC) loses its subject fields. Returns the frame and a Counter."""
    stats = collections.Counter()
    if len(df) == 0 or 'SampleType' not in df.columns:
        return df, stats
    evidence = None
    if masst_evidence and 'USI' in df.columns:
        evidence = df['USI'].map(masst_evidence)
        # a MASST table's dataset-level label does not make a solvent injection of the dataset a reference
        # material ('instrument_qc/20230426-MeOH-1.mzML' listed as 'human reference material')
        if 'filename' in df.columns and evidence.notna().any():
            solvent = pd.Series(False, index=df.index)
            has = evidence.notna()
            solvent[has] = [solvent_file_name(f) for f in df.loc[has, 'filename']]
            stats['rows_with_masst_evidence_solvent_file_not_used'] = int(solvent.sum())
            evidence = evidence.where(~solvent)
        stats['rows_with_masst_evidence'] = int(evidence.notna().sum())
    decisions = reference_materials.decide_frame(df, dataset_column='ATTRIBUTE_DatasetAccession', evidence=evidence)
    if not decisions:
        return df, stats
    before = df.loc[list(decisions), ['SampleType', 'SampleTypeSub1', 'NCBITaxonomy', 'UBERONBodyPartName']].copy()
    work = before.copy()
    work, changed = reference_materials.apply_decisions(work, decisions, stats=stats)
    if not changed:
        return df, stats
    lookups = {c: _derived_lookup(df, c) for c in DERIVED if c in df.columns}
    dataset = df['ATTRIBUTE_DatasetAccession'] if 'ATTRIBUTE_DatasetAccession' in df.columns else pd.Series('', index=df.index)
    rm = [i for i in changed if work.at[i, 'SampleType'] == reference_materials.REFERENCE_MATERIAL]
    df.loc[changed, 'SampleType'] = work.loc[changed, 'SampleType']
    df.loc[changed, 'SampleTypeSub1'] = work.loc[changed, 'SampleTypeSub1']
    for column in ('NCBITaxonomy', 'UBERONBodyPartName'):
        if column not in df.columns:
            continue
        new = work.loc[rm, column]
        diff = [i for i in rm if str(df.at[i, column]) != str(new[i])]
        set_ = _set_with_derived(df, diff, column, new, lookups, stats, f'reference_material_{column}')
        # a value the table does not know (no lineage to copy) is left out: the row keeps no organism
        unset = [i for i in diff if i not in set(set_)]
        if unset:
            _set_with_derived(df, unset, column, pd.Series(MISSING, index=unset), lookups, stats, f'reference_material_{column}')
        stats[f'reference_material_{column}_changed'] += len(set_)
    cols = [c for c in NOT_SUBJECT_COLUMNS if c in df.columns]
    cleared = df.loc[changed, cols].astype(str).ne(MISSING).any(axis=1)
    df.loc[changed, cols] = MISSING
    stats['rows_subject_fields_cleared'] = int(cleared.sum())
    for (d, st, s1), n in pd.DataFrame({'d': dataset[changed], 'st': df.loc[changed, 'SampleType'],
                                        's1': df.loc[changed, 'SampleTypeSub1']}).value_counts().sort_index().items():
        stats[f'reference_material_or_pool {d} {st} {s1}'] = int(n)
    stats['rows_changed'] = len(changed)
    return df, stats


# file-name words of QC injections that are no pool of the study's samples (blank_qc_names.NOT_POOL_WORDS: sources
# label them pool_QC too: 'calcurve_6mix', 'Std_Mix', 'QC_ESTD10000nM', 'PBQC_ISTDS', 'QCSolvent', 'acn1',
# 'MediaControl', 'SST', 'triazines_10ngml'): their study's organism is not filled in; nor that of QCs named as
# culture media ('pooled_med_1'), controls ('Ctrl1', 'Control_2') or vehicle injections ('DMSO_1')
from blank_qc_names import NOT_POOL_WORDS, not_a_pool_by_name as _not_a_pool_name
FILL_BLOCK_WORDS = NOT_POOL_WORDS | {'med', 'ctrl', 'ctl', 'control', 'controls', 'dmso', 'vehicle'}


def not_a_pool_by_name(filename):
    """whether a pooled QC's file name says it is no pool of the study's samples (blank, standard, calibration,
    solvent, medium, control or vehicle injection: FILL_BLOCK_WORDS as letter runs of the base name)"""
    return _not_a_pool_name(filename, FILL_BLOCK_WORDS)


# body parts a pooled QC's file name can name ('Pool_qc_liver_1'): a study's single body part is not filled into a
# pool named after another one (the words of the filled UBERON label, also as these synonyms)
_POOL_PART_WORDS = {'liver', 'kidney', 'brain', 'heart', 'lung', 'muscle', 'plasma', 'serum', 'urine', 'feces', 'fecal',
                    'faeces', 'faecal', 'stool', 'blood', 'milk', 'saliva', 'csf', 'cecum', 'caecum', 'cecal', 'caecal',
                    'colon', 'ileum', 'jejunum', 'duodenum', 'stomach', 'skin', 'bone', 'adipose', 'fat', 'spleen',
                    'pancreas', 'tumor', 'tumour', 'cortex', 'hippocampus', 'retina', 'testis', 'ovary', 'placenta',
                    'bile', 'sweat', 'tear', 'tears', 'breath', 'leaf', 'root', 'seed', 'fruit', 'stem', 'flower'}
_PART_SYNONYMS = {'fecal': 'feces', 'faeces': 'feces', 'faecal': 'feces', 'stool': 'feces', 'caecum': 'cecum',
                  'cecal': 'cecum', 'caecal': 'cecum', 'tumour': 'tumor', 'tears': 'tear', 'csf': 'cerebrospinal'}


def _part_word(w):
    return _PART_SYNONYMS.get(w, w)


def pool_name_names_other_part(filename, part):
    """whether a pooled QC's base name names a body part (_POOL_PART_WORDS, as letter runs) and none of the words
    of `part` (the body part to fill: 'Pool_qc_liver' and 'caecum')"""
    base = str(filename).rstrip('/').rsplit('/', 1)[-1].lower()
    named = {_part_word(w) for w in re.findall(r'[a-z]+', base) if w in _POOL_PART_WORDS}
    if not named:
        return False
    label = {_part_word(w) for w in re.findall(r'[a-z]+', str(part).lower())}
    label |= {'blood'} if label & {'plasma', 'serum'} else set()
    return not (named & label)


def pool_qc_study_values(df):
    """A pooled QC (pool of a study's samples) without NCBITaxonomy / UBERONBodyPartName takes the value
    every real sample of its dataset has (a SampleType that is no blank, QC, reference material, standard
    or product, or no SampleType but a value of the column; at least one such row, all with the same value),
    with its lineage / index columns as
    those rows have them; not a QC whose file name says it is a blank, standard, calibration, solvent,
    medium, control or vehicle injection (not_a_pool_by_name), and no body part into a QC whose file name
    names another body part (pool_name_names_other_part). Returns the frame and a Counter."""
    stats = collections.Counter()
    if len(df) == 0 or 'SampleType' not in df.columns or 'ATTRIBUTE_DatasetAccession' not in df.columns:
        return df, stats
    st = df['SampleType'].fillna('').astype(str)
    pool = st == 'pool_QC'
    if not pool.any():
        return df, stats
    names = df.loc[pool, 'filename'] if 'filename' in df.columns else pd.Series('', index=df.index[pool])
    not_pool = pd.Series(False, index=df.index)
    not_pool[pool] = [not_a_pool_by_name(f) for f in names]
    stats['pool_QC_rows_named_as_no_pool_(not_filled)'] = int(not_pool.sum())
    real = ~st.isin(reference_materials.NOT_SAMPLE_TYPES) & ~_missing(df['SampleType'])
    # (rows without a SampleType count where they have a value of the column: MTBLS15166's bacterial rows of no type
    # next to its mouse samples make the study's organism not one)
    untyped = _missing(df['SampleType'])
    dataset = df['ATTRIBUTE_DatasetAccession'].fillna('').astype(str)
    pool_ds = set(dataset[pool])
    real &= dataset.isin(pool_ds)
    untyped &= dataset.isin(pool_ds)
    for column in ('NCBITaxonomy', 'UBERONBodyPartName'):
        if column not in df.columns:
            continue
        cols = [c for c in DERIVED[column] if c in df.columns]
        sel = real | (untyped & ~_missing(df[column]))
        r = df.loc[sel, [column] + cols].astype(str)
        rd = dataset[sel]
        one = {}
        for d, g in r.groupby(rd, sort=True):
            vals = g[column]
            if (vals == MISSING).any() or vals.nunique() != 1:
                stats[f'{column}_study_not_one_value'] += 1
                continue
            derived = g[cols].drop_duplicates()
            if len(derived) != 1:
                stats[f'{column}_study_derived_columns_differ'] += 1
                continue
            one[d] = (vals.iloc[0], tuple(derived.iloc[0]))
        fill = pool & _missing(df[column]) & dataset.isin(one)
        stats[f'pool_QC_{column}_not_filled_name_says_no_pool'] = int((fill & not_pool).sum())
        fill &= ~not_pool
        if column == 'UBERONBodyPartName' and fill.any() and 'filename' in df.columns:
            other = pd.Series([pool_name_names_other_part(f, one[d][0]) for f, d in
                               zip(df.loc[fill, 'filename'], dataset[fill])], index=df.index[fill])
            other = other.reindex(df.index, fill_value=False)
            stats['pool_QC_UBERONBodyPartName_not_filled_name_names_another_part'] = int(other.sum())
            fill &= ~other
        for i in df.index[fill]:
            v, der = one[dataset.at[i]]
            df.at[i, column] = v
            for c, x in zip(cols, der):
                df.at[i, c] = x
        own = pool & ~_missing(df[column])
        stats[f'pool_QC_{column}_own_value'] = int(own.sum())
        stats[f'pool_QC_{column}_from_study'] = int(fill.sum())
        stats[f'pool_QC_{column}_none'] = int((pool & _missing(df[column])).sum())
        for d, n in dataset[fill].value_counts().sort_index().items():
            stats[f'pool_QC_{column}_from_study {d}'] = int(n)
    return df, stats


def read_masst_refmat_evidence(path):
    """{USI: evidence} from MASST_to_REDU.py's --refmat_evidence_output ({} for no or an empty file)"""
    try:
        t = pd.read_csv(path, sep='\t', dtype=str, keep_default_na=False)
    except (FileNotFoundError, pd.errors.EmptyDataError):
        return {}
    if 'USI' not in t.columns or 'evidence' not in t.columns:
        return {}
    return dict(zip(t['USI'], t['evidence']))
# ---- body parts of plants, exudates (round 2, user decisions of 2026-10-16) ----
EXUDATE_SUB1 = 'exudate'
EXUDATE_UBERON = ('exudate', 'UBERON:0007780')
# 'ear' (UBERON:0001690, the organ of hearing) written for the ear of a grass: PO inflorescence, on maize
# (Zea) PO ear inflorescence
EAR_UBERON = 'UBERON:0001690'
GRASS_FAMILY, MAIZE_GENUS = 'Poaceae', 'Zea'
EAR_PO, MAIZE_EAR_PO = 'PO:0009049', 'PO:0020136'
# SampleType of rows that are no sample of an organism (exudate never set on them)
_NOT_ORGANISM_SAMPLE = re.compile(r'(?i)blank|^pool_QC$|^reference material$|^personal care product$|^molecular_standards$')


def _organism_group(df):
    """'plant' (NCBI lineage Viridiplantae; no organism and SampleType 'plant'), 'animal' (Metazoa; no
    organism and SampleType 'animal') or '' per row"""
    kingdom = df['NCBIKingdom'] if 'NCBIKingdom' in df.columns else pd.Series('', index=df.index)
    sample_type = df['SampleType'] if 'SampleType' in df.columns else pd.Series('', index=df.index)
    no_organism = _missing(df['NCBITaxonomy']) if 'NCBITaxonomy' in df.columns else pd.Series(True, index=df.index)
    group = pd.Series('', index=df.index)
    group[(kingdom == 'Viridiplantae') | (no_organism & (sample_type == 'plant'))] = 'plant'
    group[(kingdom == 'Metazoa') | (no_organism & (sample_type == 'animal'))] = 'animal'
    return group


def plant_ontology_terms(ontology):
    """From prepare_ontologies.py's UBERON_CL_PO_ontology.csv: (PO id -> label, lower-case PO label or synonym
    -> set of PO ids, UBERON/CL label -> its index)"""
    import html
    index = ontology['UBERONOntologyIndex'].astype(str)
    po = ontology[index.str.startswith('PO:')]
    labels = dict(zip(po['UBERONOntologyIndex'], po['Label']))
    names = collections.defaultdict(set)
    for term, label, synonym in zip(po['UBERONOntologyIndex'], po['Label'], po['Synonym'].fillna('')):
        names[str(label).strip().lower()].add(term)
        if str(synonym).strip():
            names[html.unescape(str(synonym)).strip().lower()].add(term)
    other = ontology[~index.str.startswith('PO:')]
    animal = dict(zip(other['Label'], other['UBERONOntologyIndex']))
    return labels, names, animal


def exudate_samples(df):
    """Exudates exist across life (root exudates, coral exudates, culture exudates): SampleTypeSub1 'exudate'
    on every row whose body part is UBERON 'exudate', of any organism, not on blanks / QCs / standards /
    reference materials (SampleType). Run before plant_body_parts, which takes the UBERON term off the rows
    that are not of animals. Returns the frame and a Counter."""
    stats = collections.Counter()
    if 'UBERONBodyPartName' not in df.columns or 'SampleTypeSub1' not in df.columns:
        return df, stats
    sample_type = df['SampleType'].fillna('').astype(str) if 'SampleType' in df.columns else pd.Series('', index=df.index)
    exudate = (df['UBERONBodyPartName'] == EXUDATE_UBERON[0]) & ~sample_type.str.contains(_NOT_ORGANISM_SAMPLE)
    change = exudate & (df['SampleTypeSub1'] != EXUDATE_SUB1)
    dataset = df['ATTRIBUTE_DatasetAccession'] if 'ATTRIBUTE_DatasetAccession' in df.columns else pd.Series('', index=df.index)
    for (d, old), n in pd.DataFrame({'d': dataset[change], 'o': df.loc[change, 'SampleTypeSub1']}).value_counts().sort_index().items():
        stats[f'SampleTypeSub1_exudate {d} (was {old})'] = int(n)
    stats['SampleTypeSub1_exudate_rows'] = int(change.sum())
    df.loc[change, 'SampleTypeSub1'] = EXUDATE_SUB1
    return df, stats


def plant_body_parts(df, ontology):
    """Body parts of plants are Plant Ontology terms, UBERON (animal anatomy) terms are for animals:
    - on plant rows (_organism_group) a UBERON term whose label is also a PO term's label or synonym, or is
      the label of a PO 'plant <label>' term (PO's names for the plant structure: lamina, placenta,
      plant embryo, plant ovary) becomes that PO term when exactly one PO term fits; UBERON synonyms are
      not used (they give pairs such as 'anterior macula' -> 'apical meristem');
    - 'ear' on a grass (NCBIFamily Poaceae): PO 'inflorescence', on maize (NCBIGenus Zea) PO 'ear
      inflorescence';
    - other UBERON terms on plant rows are animal anatomy (olivary body for an olive, honey, milk for a
      plant milk, exudate: see exudate_samples) and are cleared;
    - UBERON 'exudate' is kept only on animal rows (other rows: cleared, their SampleTypeSub1 says exudate);
    - PO terms on animal rows are cleared, except where the row's SampleType is 'plant' (then the organism
      contradicts SampleType and body part, MSV000097043: a spider's taxon for a rhizome).
    CL terms are taxon-neutral cell types and are left as they are. A PO term is only used if the pinned
    ontology table (ontology: UBERON_CL_PO_ontology.csv) has it; UBERONOntologyIndex follows the name.
    Returns the frame and a Counter."""
    stats = collections.Counter()
    if len(df) == 0 or any(c not in df.columns for c in ('UBERONBodyPartName', 'UBERONOntologyIndex')):
        return df, stats
    labels, names, animal_terms = plant_ontology_terms(ontology)
    group = _organism_group(df)
    name = df['UBERONBodyPartName'].fillna('').astype(str)
    index = df['UBERONOntologyIndex'].fillna('').astype(str)
    # (a row without an index: the term of its name, a UBERON / CL label first: 'lamina' is both)
    by_label = {**{label: term for term, label in sorted(labels.items(), reverse=True)}, **animal_terms}
    no_index = _missing(index) & ~_missing(name)
    index[no_index] = name[no_index].map(lambda n: by_label.get(n, ''))
    is_uberon = index.str.startswith('UBERON:')
    family = df['NCBIFamily'] if 'NCBIFamily' in df.columns else pd.Series('', index=df.index)
    genus = df['NCBIGenus'] if 'NCBIGenus' in df.columns else pd.Series('', index=df.index)
    sample_type = df['SampleType'] if 'SampleType' in df.columns else pd.Series('', index=df.index)

    def po_of(label):
        key = label.strip().lower()
        found = names.get(key, set()) | {t for t in names.get('plant ' + key, set()) if labels[t].lower() == 'plant ' + key}
        return next(iter(found)) if len(found) == 1 else None

    # (the index only changes on the rows a rule below decides on)
    new_index, new_name = df['UBERONOntologyIndex'].fillna('').astype(str).copy(), name.copy()
    # plant rows with a UBERON term
    rows = (group == 'plant') & is_uberon
    for (label, term), idx in df[rows].groupby([name[rows], index[rows]], sort=True).groups.items():
        if term == EAR_UBERON:
            for i in idx:
                po = (MAIZE_EAR_PO if genus[i] == MAIZE_GENUS else EAR_PO) if family[i] == GRASS_FAMILY else None
                new_index[i], new_name[i] = (po, labels[po]) if po in labels else (MISSING, MISSING)
            continue
        po = po_of(label) if term != EXUDATE_UBERON[1] else None
        new_index[idx], new_name[idx] = (po, labels[po]) if po else (MISSING, MISSING)
    # UBERON exudate on rows that are not of animals
    other_exudate = (group != 'animal') & (group != 'plant') & (index == EXUDATE_UBERON[1])
    new_index[other_exudate], new_name[other_exudate] = MISSING, MISSING
    # PO terms on animal rows
    po_on_animal = (group == 'animal') & index.str.startswith('PO:')
    kept = po_on_animal & (sample_type == 'plant')
    stats['PO_on_animal_rows_kept_SampleType_plant'] = int(kept.sum())
    clear = po_on_animal & ~kept
    new_index[clear], new_name[clear] = MISSING, MISSING

    change = (new_name != name) | (new_index != df['UBERONOntologyIndex'].fillna('').astype(str))
    dataset = df['ATTRIBUTE_DatasetAccession'] if 'ATTRIBUTE_DatasetAccession' in df.columns else pd.Series('', index=df.index)
    log = pd.DataFrame({'d': dataset[change], 'g': group[change], 'f': name[change] + ' ' + df.loc[change, 'UBERONOntologyIndex'].astype(str),
                        't': new_name[change] + ' ' + new_index[change]})
    for (d, g, f, t), n in log.value_counts().sort_index().items():
        stats[f'body_part {d} {g or "other"}: {f} -> {t}'] = int(n)
    stats['body_part_rows_changed'] = int(change.sum())
    df.loc[change, 'UBERONBodyPartName'] = new_name[change]
    df.loc[change, 'UBERONOntologyIndex'] = new_index[change]
    return df, stats


# ENVO material of a water sample -> its SampleTypeSub1 (the ENVO labels of ENVO_material_ontology.csv;
# not 'fresh water', 'brackish water', 'estuarine water', 'spring water', drinking or tap water)
WATER_MATERIAL_SUB1 = {'surface water': 'water_surface', 'river water': 'water_surface', 'lake water': 'water_surface',
                       'pond water': 'water_surface', 'stream water': 'water_surface', 'ditch water': 'water_surface',
                       'sea water': 'water_seawater', 'coastal sea water': 'water_seawater',
                       'waste water': 'water_waste', 'treated wastewater': 'water_waste', 'sewage': 'water_waste',
                       'groundwater': 'water_ground', 'stormwater': 'water_storm'}


def water_kind_from_material(df):
    """SampleTypeSub1 of environmental water samples whose sheet states the kind of water only in
    ENVOEnvironmentMaterial (GNPS sheets: 'river water' with SampleTypeSub1 missing): the water_* value of the
    material (WATER_MATERIAL_SUB1), only where SampleType is 'environmental' and SampleTypeSub1 is missing.
    Returns the frame and a Counter."""
    stats = collections.Counter()
    if any(c not in df.columns for c in ('SampleType', 'SampleTypeSub1', 'ENVOEnvironmentMaterial')):
        return df, stats
    kind = df['ENVOEnvironmentMaterial'].map(WATER_MATERIAL_SUB1)
    fill = (df['SampleType'] == 'environmental') & _missing(df['SampleTypeSub1']) & kind.notna()
    dataset = df['ATTRIBUTE_DatasetAccession'] if 'ATTRIBUTE_DatasetAccession' in df.columns else pd.Series('', index=df.index)
    for (d, m), n in pd.DataFrame({'d': dataset[fill], 'm': df.loc[fill, 'ENVOEnvironmentMaterial']}).value_counts().sort_index().items():
        stats[f'SampleTypeSub1_from_water_material {d} {m}'] = int(n)
    stats['SampleTypeSub1_from_water_material_rows'] = int(fill.sum())
    df.loc[fill, 'SampleTypeSub1'] = kind[fill]
    return df, stats


def placeholder_factor_names(df):
    """SampleGroupFreeText parts whose factor name is a template placeholder ('<enter factor name>: no drug',
    'Factor name: ...', 'NA: Responder') lose the name (sample_groups.strip_placeholder_factor_names), on all
    rows (also the rows of earlier runs). Runs before final_text_columns.py, so its Diet / group rules see
    the texts without placeholders. Returns the frame and a Counter."""
    from sample_groups import strip_placeholder_factor_names
    stats = collections.Counter()
    if 'SampleGroupFreeText' not in df.columns:
        return df, stats
    texts = df['SampleGroupFreeText']
    new = strip_placeholder_factor_names(texts)
    change = new != texts
    dataset = df['ATTRIBUTE_DatasetAccession'] if 'ATTRIBUTE_DatasetAccession' in df.columns else pd.Series('', index=df.index)
    for d, n in dataset[change].value_counts().sort_index().items():
        stats[f'SampleGroupFreeText_placeholder_names {d}'] = int(n)
    stats['SampleGroupFreeText_placeholder_names_rows'] = int(change.sum())
    df.loc[change, 'SampleGroupFreeText'] = new[change]
    return df, stats


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Rules applied to every row of the final table')
    parser.add_argument('input_tsv')
    parser.add_argument('output_tsv')
    parser.add_argument('--path_to_sex_specific_parts', default=None,
                        help='UBERON_sex_specific_parts.csv of prepare_ontologies.py (BiologicalSex from the body part)')
    parser.add_argument('--masst_refmat_evidence', default=None,
                        help='masst_refmat_evidence.tsv of MASST_to_REDU.py (reference materials MASST names for '
                             'GNPS-described files)')
    parser.add_argument('--path_to_uberon_cl_po_csv', default=None,
                        help='UBERON_CL_PO_ontology.csv of prepare_ontologies.py (Plant Ontology terms of plant body parts)')
    args = parser.parse_args()

    # all cells as the previous step wrote them (no NA parsing: e.g. 'NA' stays 'NA')
    df = pd.read_csv(args.input_tsv, sep='\t', dtype=str, keep_default_na=False, na_filter=False, low_memory=False)
    df, stats = collapse_cell_breaks(df)
    df, care_stats = personal_care_products(df)
    stats.update(care_stats)
    df, rm_stats = reference_materials_pass(
        df, read_masst_refmat_evidence(args.masst_refmat_evidence) if args.masst_refmat_evidence else None)
    stats.update({f'refmat: {k}': v for k, v in rm_stats.items()})
    df, pool_stats = pool_qc_study_values(df)
    stats.update({f'pool_QC: {k}': v for k, v in pool_stats.items()})
    df, life_stage_stats = life_stage_humans_only(df)
    stats.update(life_stage_stats)
    if args.path_to_sex_specific_parts:
        parts = pd.read_csv(args.path_to_sex_specific_parts, dtype=str, keep_default_na=False)
        df, sex_stats = sex_from_body_part(df, parts)
        stats.update(sex_stats)
    df, exudate_stats = exudate_samples(df)
    stats.update(exudate_stats)
    if args.path_to_uberon_cl_po_csv:
        ontology = pd.read_csv(args.path_to_uberon_cl_po_csv, dtype=str, keep_default_na=False)
        df, plant_stats = plant_body_parts(df, ontology)
        stats.update(plant_stats)
    df, water_stats = water_kind_from_material(df)
    stats.update(water_stats)
    df, placeholder_stats = placeholder_factor_names(df)
    stats.update(placeholder_stats)
    for key, value in sorted(stats.items()):
        print(f"  {key}: {value}", flush=True)
    write_tsv_atomic(df, args.output_tsv)
