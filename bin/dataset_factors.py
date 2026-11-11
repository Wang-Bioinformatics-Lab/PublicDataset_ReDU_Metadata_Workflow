"""SampleGroupFreeText of MassIVE sheets from the datasets' own metadata (gnps_downloader.py).

Many MassIVE datasets describe their sample groups in columns ReDU has no field for: extra columns of the ReDU
sheet itself (the validator drops them) and GNPS-format metadata files next to it (a file-name column and
ATTRIBUTE_ / other columns, e.g. 'other/Turek Pilot Study metadata.txt' with ATTRIBUTE_diet_group). Their
factors are added to the sheet's SampleGroupFreeText as 'Factor: value' parts, joined by file name without
folders and extension (gnps_downloader._sample_key):
  - a factor is added only where it tells the sheet's samples apart: 2 to n/2 different values over the n
    samples it gives a value, numbers with more than MAX_NUMERIC_LEVELS values are measurements (no groups).
    The samples are the real samples (sample_groups.is_real_sample); a sheet without any has its reference
    material and personal care product rows as samples (sample_groups.OBJECT_TYPES: MSV000083605, tablets of
    6 drugs, all 'reference material'), not its standards (dilution series, adducts, collision energies);
  - columns that stand for a ReDU template column, file / sample / subject identifiers, cages, dates and
    technical, storage, sequencing-preparation and administrative columns, numbers as names (SKIP_NAME) and
    numberings with more than MAX_NUMERIC_LEVELS values (_NUMBERED: identifiers) are left out, and so are
    columns whose values are dates, also day-month dates Excel made of IDs ('1-Jan', _DATE; calendar years are
    groups: a harvest year),
    file names (_FILE_VALUE), database identifiers ('gbif:3790961', _DB_ID) and the ionization mode
    (_POLARITY_NAME and only polarity values: 'Metadata Table: posmode'; a value 'positive' alone is no mode,
    'RT PCR diagnostic: positive');
  - a column that stands for a ReDU field (_FIELD_PROXY: species, taxonomy, organism, sex, country, ..., or a
    template column under another name, e.g. ATTRIBUTE_Diet) is still added where the field does not hold
    what the column says (_proxy_adds): the column splits a value of the field (Species: 9 Salinispora species,
    NCBITaxonomy: the genus; each column value has one field value), or the field is empty there and the
    column's values are none of the field's allowed terms (ATTRIBUTE_Diet: Before / After). Not where the
    column and the field disagree, where the field is empty and the column could fill it, nor for SampleType /
    SampleTypeSub1, numeric fields, or first-generation ReDU columns (Curated_ / Subject_ / Study_ /
    Analysis_: ReDU fields under their old names); blank / media / not-applicable values of such a column are
    no value ('media blank', 'no_apply');
  - values that differ only in case and spaces are one value ('4wks', '4 wks': the most frequent spelling); a ';'
    in a value (it separates the parts of the text) is written as ',' ('sugar;citric acid'); Qiita / INSDC
    missing-value terms ('not entered', 'missing: not provided') are no value;
  - a factor whose values are a renaming of an already added factor's (the same rows, one to one, at least 3
    values) is a copy and left out (ATRRIBUTE_SampleType = Timepoint);
  - files with the same columns are one table (per-strain / per-polarity files of one design: MSV000085261, 56
    files of 36-60 of its 1277 files); a table is used when its file-name column names at least
    MIN_MATCHED_SHARE of the sheet's files;
  - a part whose factor the row's text already has is not added; the sheet's own columns come first, then the
    files by path.
Measured on the 126 MassIVE datasets of the run15b missing-metadata report (2026-11-09): the factors the report
found missing were added for 164 of 204 columns; the remaining ones are ReDU fields or not groups. The field,
reference material, same-columns and copy rules: run 18 inputs, 2026-11-10 (tmp_redu_run15/diag19).
"""
import collections
import re

import pandas as pd

from sample_groups import OBJECT_TYPES, is_real_sample

MISSING = 'missing value'
NA = {'', 'na', 'n/a', 'nan', 'none', 'null', 'not applicable', 'not available', 'unknown', 'not specified',
      'not collected', 'missing', 'missing value', '-', '--', '?', 'not provided', 'not reported', 'unspecified',
      'not determined', 'nd', 'n.d.',
      # Qiita / INSDC missing-value terms (MSV000086158's Global FoodOmics sheet: 'not entered')
      'not entered', 'restricted access', 'missing: not collected', 'missing: not provided',
      'missing: restricted access', 'missing: not reported', 'missing: control sample', 'missing: lab stock',
      'missing: third party data'}
MIN_MATCHED_SHARE = 0.3
MAX_NUMERIC_LEVELS = 6
# file / sample / subject identifiers, technical, storage and administrative columns, and what ReDU's own fields
# hold (MIxS env_* / geo_loc_name / lat-lon / taxon / age / sex / organism / country)
SKIP_NAME = re.compile(
    r'(?i)file|path|usi\b|^scan|^unnamed|^index$|^row|^#|^id$|sample[ _]?(id|name|code|number|no)|injection|'
    r'vial|well|plate|tube|barcode|position|order|^run|batch|datasetaccession|massiveid|weight|volume|'
    r'amount|timestamp|description|comment|note|qiita|dna|sequenc|specimen|^sample type$|taxon|'
    r'scientific name|latitude|longitude|elevation|altitude|depth|geo loc|^env |^empo|title|library|'
    r'instrument|operator|^extraction|ms2|ms1|polarity|ionization|column|chromatograph|^mass|'
    r'^age( |$)|^sex$|^gender$|^organism$|^species$|^country$|replicate|^rep( |$)|ontologyindex|room|rat number|'
    r'mouse number|^tank$|^public$|units?$|^year$|^global$|host body|body (habitat|product|site)|'
    r'sampletype ?as ?recorded|sampletype source|freezer|box|rack|shelf|storage|^date|^ms |anonymized|taxid|'
    r'common name|center|lab code|pf code|animal number|gps|coordinates|bodypart|\bid$|code$|catalog|#|ear tag|'
    r'experimenter|funding|^labels?$|line color|strain number|ibt number|'
    r'^(host )?subject( id)?$|^(host|patient|participant|individual|animal|mouse|donor|person) ?id$|'
    r'^(patient|participant|individual|donor)$|'
    # cages (housing, not a group) and dates (SampleCollectionDateandTime, not a group)
    r'(?<![a-z])(cages?|dates?)(?![a-z])|'
    # numbers as names (feature tables: 'm/z RT' columns), catalogue numbers, sequencing preparation (Qiita)
    r'^[-\d.\s]+$|(?<![a-z])lot( number)?$|primer|linker|robot|^platform$|^plating$|'
    r'target (gene|subfragment)|tool$|^project name$|emp status|samp size|assigned from geo|has extracted data|'
    r'required sample info')
# values that are dates ('3/1/19', '2019-03-07', '07-Mar-2019', '1/31/2025 8:45 AM'; not a calendar year '2019': a harvest year is a group)
_MONTH = r'(jan|feb|mar|apr|may|jun|jul|aug|sep|oct|nov|dec)[a-z]*'
_DATE = re.compile(r'(?i)(\d{1,4}[-/]\d{1,2}[-/]\d{1,4}|\d{1,2}\.\d{1,2}\.\d{4}|\d{4}\.\d{1,2}\.\d{1,2})([ t]\d{1,2}:\d{2}(:\d{2})?( ?[ap]m)?)?|'
                   r'\d{1,2}[- ]?' + _MONTH + r'[- ]?\d{2,4}|' + _MONTH + r' \d{1,2},? \d{4}|'
                   # day-month without a year: what Excel makes of IDs like 1-1, 10-11 (MSV000083605 '1-Jan')
                   r'\d{1,2}-' + _MONTH + r'|' + _MONTH + r'-\d{1,2}')
# a numbering ('Fish Number', 'Museum Number update', 'IBT-number') with more than MAX_NUMERIC_LEVELS values
# is an identifier ('cohort number' with 3 values is a group)
_NUMBERED = re.compile(r'(?i)number( update)?$')
_NUMBER = re.compile(r'[<>]?\s*-?\d+(?:\.\d+)?(?:e-?\d+)?', re.I)
# values of a measurement below / above the range ('LOD', '<LOQ', 'n.d.'): counted with its numbers
_CENSORED = re.compile(r'(?i)[<>]?\s*(lod|loq|lld|bld|blq|bdl|nd|n\.d\.|below .*|above .*|not detected)')
# file names as values (instrument methods, raw files)
_FILE_VALUE = re.compile(r'(?i)\S+\.(meth|method|raw|mzml|mzxml|d|wiff2?|txt|csv|tsv|xlsx?|m|qgd|lcd)')


# database identifiers as values ("['gbif:3790961']", 'ncbi:1008910', 'UBERON:0001988' are ReDU's index columns)
_DB_ID = re.compile(r"(?i)\W*[a-z][a-z_]*:\s*\d+\W*(,\s*\W*[a-z][a-z_]*:\s*\d+\W*)*")
# the ionization mode: a column with such a name and only polarity values
_POLARITY_NAME = re.compile(r'(?i)mode|polar|ioni[sz]ation|^esi$|metadata ?table')
_POLARITY_VALUE = re.compile(r'(?i)\s*(esi)?\s*[+-]?\s*(pos|neg|positive|negative)?\s*(_?mode|ion mode|ionization)?\s*[+-]?\s*')
# columns that stand for a ReDU field, by name ('Species', 'taxonomy', 'scientific name', 'Gender'; the other
# names SKIP_NAME leaves out for a field are measurements or places: no groups ReDU's field misses)
_FIELD_PROXY = [(re.compile(r'(?i)taxon|species|organism|scientific name|^genus$|^strain$'), 'NCBITaxonomy'),
                (re.compile(r'(?i)^sex$|^gender$'), 'BiologicalSex'),
                (re.compile(r'(?i)^country$'), 'Country'),
                (re.compile(r'(?i)host body|body (habitat|product|site)|bodypart'), 'UBERONBodyPartName')]
# fields a column is compared with (not SampleType / SampleTypeSub1: the template's types next to the
# submitter's own words are no groups; not the numeric fields: measurements)
PROXY_FIELDS = {'NCBITaxonomy', 'BiologicalSex', 'Diet', 'LifeStage', 'UBERONBodyPartName', 'HealthStatus',
                'DOIDCommonName', 'Country', 'Exposure', 'CellLine', 'FastingStatus', 'SmokingStatus', 'HPOCommonName',
                'ColonizationStatus', 'ENVOEnvironmentBiome', 'ENVOEnvironmentMaterial'}
# names of an identifier, never a field's value
_ID_NAME = re.compile(r'(?i)number|\bid$|code|file|date|index|sample[ _]?(id|name|no)|source')
# first-generation ReDU columns (the old names of ReDU fields; gnps_downloader.LEGACY_PREFIXES)
_LEGACY = re.compile(r'(?i)^(attribute_)?(curated|subject|analysis|study)_')
# values of a field column that are no value of the field: blanks, media, QCs, 'does not apply'
_NOT_A_VALUE = re.compile(r'(?i)blank|media\b|medium|^qcs?\b|_qc\b|solvent|no_?apply|not_?applicable|does_?not')
_MISSING_TERMS = {'missing value', 'not collected', 'not applicable', 'not specified', 'unknown', ''}


def factor_name(column):
    """'ATTRIBUTE_FungalStrain' -> 'Fungal Strain', 'diet_group' -> 'diet group' (also misspelt prefixes:
    'ATRRIBUTE_SampleType' -> 'Sample Type')"""
    name = re.sub(r'(?i)^at{1,2}r{1,2}[a-z]{3,5}_', '', str(column).strip())
    name = re.sub(r'(?<=[a-z])(?=[A-Z])', ' ', name)
    return ' '.join(name.replace('_', ' ').split())


def _norm(text):
    return re.sub(r'[^a-z0-9]+', '', str(text).lower())


def _skipped(name):
    return bool(SKIP_NAME.search(name) or SKIP_NAME.search(name.replace(' ', '')))


def _value(v):
    v = ' '.join(str(v).split()) if isinstance(v, str) else ('' if pd.isna(v) else str(v))
    return '' if v.lower() in NA else v


def _part_value(v):
    """a value as the value of a 'Factor: value' part: no ';' (it separates the parts; ingredient lists
    'sugar;citric acid' -> 'sugar, citric acid')"""
    return ' '.join(re.sub(r'\s*;\s*', ', ', v).split()).strip(', ') if ';' in v else v


def proxy_field(column, name, template_column):
    """The ReDU field (PROXY_FIELDS) a column stands for, or None: a template column under another name
    ('ATTRIBUTE_Diet', 'Sex'; not the field's own column), or a name SKIP_NAME leaves out for a field
    (_FIELD_PROXY; not identifiers, not first-generation ReDU columns)"""
    if _LEGACY.match(str(column)):
        return None
    field = template_column(column)
    if field is not None:
        return field if field in PROXY_FIELDS and field != str(column).strip() else None
    if _ID_NAME.search(name):
        return None
    for rx, field in _FIELD_PROXY:
        if rx.search(name) or rx.search(name.replace(' ', '')):
            return field
    return None


def factors_of(df, file_column, template_column, sample_key):
    """{factor name: {sample key: value}} of a table's factor columns, values conflicting within a key left out,
    and {factor name: ReDU field} of the columns that stand for a field (proxy_field)"""
    keys = df[file_column].fillna('').astype(str).map(sample_key)
    out, fields = {}, {}
    for c in df.columns:
        name = factor_name(c)
        if c == file_column or not name or name in out:
            continue
        field = None
        if template_column(c) is not None or _skipped(name):
            field = proxy_field(c, name, template_column)
            if field is None:
                continue
        values = collections.defaultdict(set)
        for k, v in zip(keys, df[c]):
            v = _part_value(_value(v))
            if field is not None and _NOT_A_VALUE.search(v):
                v = ''
            if k and v:
                values[k].add(v)
        out[name] = {k: next(iter(v)) for k, v in values.items() if len(v) == 1}
        if field is not None:
            fields[name] = field
    return out, fields


def file_column_of(df, sheet_keys, sample_key):
    """the column whose values name the most of the sheet's files, and how many"""
    best, n = None, 0
    for c in df.columns:
        hit = len(set(df[c].dropna().astype(str).map(sample_key)) & sheet_keys)
        if hit > n:
            best, n = c, hit
    return best, n


def same_column_tables(tables):
    """[(name, table)] with the tables of identical columns concatenated (in their order; the name of the
    first, '(+k files with the same columns)')"""
    by = collections.OrderedDict()
    for name, table in tables:
        by.setdefault(tuple(table.columns), []).append((name, table))
    return [(group[0][0] + (f' (+{len(group) - 1} files with the same columns)' if len(group) > 1 else ''),
             pd.concat([t for _, t in group], ignore_index=True) if len(group) > 1 else group[0][1])
            for group in by.values()]


def sample_rows(sheet):
    """the rows whose groups count: the real samples (is_real_sample), or in a sheet without any, its
    reference material and personal care product rows"""
    if 'SampleType' not in sheet.columns:
        return pd.Series(True, index=sheet.index)
    real = sheet['SampleType'].map(is_real_sample)
    if not real.any():
        real = sheet['SampleType'].fillna('').astype(str).str.strip().isin(OBJECT_TYPES)
    return real


def _field_value(field, value, allowed):
    v = '' if pd.isna(value) else ' '.join(str(value).split()).lower()
    if v in _MISSING_TERMS:
        return ''
    terms = allowed.get(field)
    return v if not terms or v in terms else ''


def _proxy_adds(values, real, field_values, closed_terms):
    """Whether a column standing for a ReDU field adds what the field does not hold: each of its values has
    one field value and some field value carries 2 or more of its values (the field is coarser: a genus for
    species), or the field is empty there and the column's values are not the field's allowed terms (they
    would be blanked). Not where they disagree (a column value with several field values) or where the
    field is empty and the column could fill it."""
    to_field, to_column = collections.defaultdict(set), collections.defaultdict(set)
    for v, r, f in zip(values, real, field_values):
        if v and r:
            to_field[v.lower()].add(f)
            to_column[f].add(v.lower())
    if not to_field or any(len(s) > 1 for s in to_field.values()) or not any(len(s) > 1 for s in to_column.values()):
        return False
    if set(to_column) != {''}:
        return '' not in to_column
    return bool(closed_terms) and sum(v not in closed_terms for v in to_field) >= 0.9 * len(to_field)


def _spelling(values):
    """values that differ only in case and spaces as their most frequent spelling"""
    counts = collections.Counter(v for v in values if v)
    best = {}
    for v, n in counts.most_common():
        best.setdefault(re.sub(r'\s+', '', v.lower()), v)
    return [best[re.sub(r'\s+', '', v.lower())] if v else v for v in values]


# values of yes / no flags (_is_copy)
_FLAG_VALUES = {'yes', 'no', 'y', 'n', 'true', 'false', '0', '1', 'present', 'absent', 'positive', 'negative', '+', '-',
                'x'}


def _is_copy(values, added, real):
    """the values (per row) a renaming of an added factor's (the same rows, one to one, at least 3 values), or the
    added factor's own values on the same rows (any number: 'susceptiblility' of a second metadata file =
    'susceptibility', High / Low, MSV000089633; not yes / no flags, which several factors share: the drugs of one
    patient)"""
    for other in added.values():
        pairs = {(a, b) for a, b, r in zip(values, other, real) if r and (a or b)}
        if not pairs or any(not a or not b for a, b in pairs):
            continue
        if all(a == b for a, b in pairs) and not {a.lower() for a, _ in pairs} <= _FLAG_VALUES:
            return True
        left, right = collections.defaultdict(set), collections.defaultdict(set)
        for a, b in pairs:
            left[a].add(b)
            right[b].add(a)
        if len(left) >= 3 and all(len(s) == 1 for s in left.values()) and all(len(s) == 1 for s in right.values()):
            return True
    return False


def add_factors(sheet, tables, template_column, sample_key, allowed_terms=None):
    """The sheet (a ReDU-format table with a 'filename' column) with SampleGroupFreeText completed from its own
    extra columns and the tables (other metadata files of its dataset, [(name, table)]), and [(source, factor,
    rows, levels)] of the factors added. allowed_terms (allowed_terms.json): the fields' allowed terms, for the
    columns that stand for a field (without it, only where the field is coarser)."""
    allowed = {k: {str(t).lower() for t in v.get('allowed_values') or []}
               for k, v in (allowed_terms or {}).items() if len(v.get('allowed_values') or []) >= 2}
    keys = sheet['filename'].fillna('').astype(str).map(sample_key)
    sheet_keys = set(keys) - {''}
    real = sample_rows(sheet)
    sources = [('own columns',) + factors_of(sheet, 'filename', template_column, sample_key)]
    for name, table in same_column_tables(tables):
        column, n = file_column_of(table, sheet_keys, sample_key)
        if column is not None and n >= max(2, MIN_MATCHED_SHARE * len(sheet_keys)):
            sources.append((name,) + factors_of(table, column, template_column, sample_key))
    added, log = collections.OrderedDict(), []
    for source, factors, fields in sources:
        for name, by_key in factors.items():
            if _norm(name) in {_norm(a) for a in added}:
                continue
            values = _spelling([by_key.get(k) for k in keys])
            given = [v for v, r in zip(values, real) if v and r]
            levels = len(set(given))
            numeric = sum(bool(_NUMBER.fullmatch(g) or _CENSORED.fullmatch(g)) for g in given) >= 0.9 * len(given) \
                if given else False
            if given and sum(bool(_DATE.fullmatch(g) or _FILE_VALUE.fullmatch(g) or _DB_ID.fullmatch(g))
                             for g in given) >= 0.9 * len(given):
                continue
            if given and _POLARITY_NAME.search(name) and all(_POLARITY_VALUE.fullmatch(g) for g in set(given)):
                continue
            if levels > MAX_NUMERIC_LEVELS and _NUMBERED.search(name):
                continue
            if not (len(given) >= 2 and 2 <= levels <= len(given) / 2 and not (numeric and levels > MAX_NUMERIC_LEVELS)):
                continue
            field = fields.get(name)
            if field is not None:
                if numeric:
                    continue
                field_values = (sheet[field].map(lambda v: _field_value(field, v, allowed)) if field in sheet.columns
                                else pd.Series('', index=sheet.index))
                if not _proxy_adds(values, real, field_values, allowed.get(field)):
                    continue
            if _is_copy(values, added, real):
                continue
            added[name] = values
            log.append((source, name, len(given), levels))
    if not added:
        return sheet, log
    sheet = sheet.copy()
    texts = sheet['SampleGroupFreeText'] if 'SampleGroupFreeText' in sheet.columns else pd.Series('', index=sheet.index)
    out = []
    for i, text in enumerate(texts):
        text = _value(text)
        have = {_norm(p.split(': ', 1)[0]) for p in text.split('; ') if p}
        parts = [f'{n}: {v[i]}' for n, v in added.items() if v[i] and _norm(n) not in have]
        out.append('; '.join([text] * bool(text) + parts) or MISSING)
    sheet['SampleGroupFreeText'] = out
    return sheet, log
