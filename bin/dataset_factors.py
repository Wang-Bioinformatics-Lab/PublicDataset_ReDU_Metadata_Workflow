"""SampleGroupFreeText of MassIVE sheets from the datasets' own metadata (gnps_downloader.py).

Many MassIVE datasets describe their sample groups in columns ReDU has no field for: extra columns of the ReDU
sheet itself (the validator drops them) and GNPS-format metadata files next to it (a file-name column and
ATTRIBUTE_ / other columns, e.g. 'other/Turek Pilot Study metadata.txt' with ATTRIBUTE_diet_group). Their
factors are added to the sheet's SampleGroupFreeText as 'Factor: value' parts, joined by file name without
folders and extension (gnps_downloader._sample_key):
  - a factor is added only where it tells the sheet's real samples apart (sample_groups.is_real_sample): 2 to n/2
    different values over the n real samples it gives a value, numbers with more than MAX_NUMERIC_LEVELS values
    are measurements (no groups);
  - columns that stand for a ReDU template column, file / sample / subject identifiers, cages, dates and
    technical, storage and administrative columns (SKIP_NAME) are left out, and so are columns whose values are
    dates or calendar years (_DATE);
  - a file is used when its file-name column names at least MIN_MATCHED_SHARE of the sheet's files;
  - a part whose factor the row's text already has is not added; the sheet's own columns come first, then the
    files by path.
Measured on the 126 MassIVE datasets of the run15b missing-metadata report (2026-11-09): the factors the report
found missing were added for 164 of 204 columns; the remaining ones are ReDU fields or not groups.
"""
import collections
import re

import pandas as pd

from sample_groups import is_real_sample

MISSING = 'missing value'
NA = {'', 'na', 'n/a', 'nan', 'none', 'null', 'not applicable', 'not available', 'unknown', 'not specified',
      'not collected', 'missing', 'missing value', '-', '--', '?', 'not provided', 'not reported', 'unspecified',
      'not determined', 'nd', 'n.d.'}
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
    r'(?<![a-z])(cages?|dates?)(?![a-z])')
# values that are dates or calendar years ('3/1/19', '2019-03-07', '07-Mar-2019', '2019')
_DATE = re.compile(r'(?i)(\d{1,4}[-/]\d{1,2}[-/]\d{1,4}|\d{1,2}\.\d{1,2}\.\d{4}|\d{4}\.\d{1,2}\.\d{1,2})([ t]\d{1,2}:\d{2}(:\d{2})?)?|'
                   r'\d{1,2}[- ]?(jan|feb|mar|apr|may|jun|jul|aug|sep|oct|nov|dec)[a-z]*[- ]?\d{2,4}|'
                   r'(jan|feb|mar|apr|may|jun|jul|aug|sep|oct|nov|dec)[a-z]* \d{1,2},? \d{4}|(19|20)\d{2}')
_NUMBER = re.compile(r'-?\d+(?:\.\d+)?(?:e-?\d+)?', re.I)


def factor_name(column):
    """'ATTRIBUTE_FungalStrain' -> 'Fungal Strain', 'diet_group' -> 'diet group'"""
    name = re.sub(r'(?i)^attribute_', '', str(column).strip())
    name = re.sub(r'(?<=[a-z])(?=[A-Z])', ' ', name)
    return ' '.join(name.replace('_', ' ').split())


def _norm(text):
    return re.sub(r'[^a-z0-9]+', '', str(text).lower())


def _skipped(name):
    return bool(SKIP_NAME.search(name) or SKIP_NAME.search(name.replace(' ', '')))


def _value(v):
    v = ' '.join(str(v).split()) if isinstance(v, str) else ('' if pd.isna(v) else str(v))
    return '' if v.lower() in NA else v


def factors_of(df, file_column, template_column, sample_key):
    """{factor name: {sample key: value}} of a table's factor columns, values conflicting within a key left out"""
    keys = df[file_column].fillna('').astype(str).map(sample_key)
    out = {}
    for c in df.columns:
        name = factor_name(c)
        if c == file_column or not name or template_column(c) is not None or _skipped(name) or name in out:
            continue
        values = collections.defaultdict(set)
        for k, v in zip(keys, df[c]):
            v = _value(v)
            if k and v:
                values[k].add(v)
        out[name] = {k: next(iter(v)) for k, v in values.items() if len(v) == 1}
    return out


def file_column_of(df, sheet_keys, sample_key):
    """the column whose values name the most of the sheet's files, and how many"""
    best, n = None, 0
    for c in df.columns:
        hit = len(set(df[c].dropna().astype(str).map(sample_key)) & sheet_keys)
        if hit > n:
            best, n = c, hit
    return best, n


def add_factors(sheet, tables, template_column, sample_key):
    """The sheet (a ReDU-format table with a 'filename' column) with SampleGroupFreeText completed from its own
    extra columns and the tables (other metadata files of its dataset, [(name, table)]), and [(source, factor,
    rows, levels)] of the factors added"""
    keys = sheet['filename'].fillna('').astype(str).map(sample_key)
    sheet_keys = set(keys) - {''}
    real = sheet['SampleType'].map(is_real_sample) if 'SampleType' in sheet.columns else pd.Series(True, index=sheet.index)
    sources = [('own columns', factors_of(sheet, 'filename', template_column, sample_key))]
    for name, table in tables:
        column, n = file_column_of(table, sheet_keys, sample_key)
        if column is not None and n >= max(2, MIN_MATCHED_SHARE * len(sheet_keys)):
            sources.append((name, factors_of(table, column, template_column, sample_key)))
    added, log = collections.OrderedDict(), []
    for source, factors in sources:
        for name, by_key in factors.items():
            if _norm(name) in {_norm(a) for a in added}:
                continue
            values = [by_key.get(k) for k in keys]
            given = [v for v, r in zip(values, real) if v and r]
            levels = len(set(given))
            numeric = sum(bool(_NUMBER.fullmatch(g)) for g in given) >= 0.9 * len(given) if given else False
            if given and sum(bool(_DATE.fullmatch(g)) for g in given) >= 0.9 * len(given):
                continue
            if len(given) >= 2 and 2 <= levels <= len(given) / 2 and not (numeric and levels > MAX_NUMERIC_LEVELS):
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
