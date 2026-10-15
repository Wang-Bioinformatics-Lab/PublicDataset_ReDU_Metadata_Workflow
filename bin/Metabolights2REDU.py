import os
import re
import requests
import pandas as pd
from bs4 import BeautifulSoup
import argparse
import json
import time
import html
import math
import itertools
import traceback
import csv
import io
from urllib.parse import quote
import numpy as np
from tqdm import tqdm
from REDU_conversion_functions import age_category
from REDU_conversion_functions import get_taxonomy_id_from_name__allowedTerms
from REDU_conversion_functions import get_taxonomy_info
from REDU_conversion_functions import merge_repeated_fileobservations
from read_and_validate_redu_from_github import complete_and_fill_REDU_table
from REDU_conversion_functions import find_column_after_target_column
from REDU_conversion_functions import map_instrument_to_allowed
from REDU_conversion_functions import convert_smoking, bmi_to_numeric, convert_diet
from REDU_conversion_functions import normalize_redu_datetime
from ms_instrument_matching import match_mass_spectrometer
# time point, fasting and sample group rules shared with the Workbench converter
from MWB_to_REDU import (_norm_key, _TP_UNIT_WORDS, timepoint_key_tier, timepoint_hours, _timepoint_key_unit,
                         sample_timepoint_hours, sample_fasting_status,
                         _GROUP_KEY_EXCLUDE_RE, _GROUP_ID_KEY_RE, _GROUP_NA_VALUES, _GROUP_SID, _GROUP_KV, _GROUP_PART,
                         finalize_sample_groups)




# MS file formats REDU indexes, and the compression suffixes MetaboLights uses for them
# (x.raw.zip, x.d.zip, x.cdf.gz, x.mzXML.lzma ...). A compressed file is listed under its
# real path; only single files of a supported format qualify (a bare x.zip is not indexed).
# one session (kept-alive connection) for all requests to EBI, made one at a time: EBI bans
# IPs that open many concurrent connections
EBI_SESSION = requests.Session()

MS_FILE_EXTENSIONS = ['.mzml', '.mzxml', '.cdf', '.raw', '.wiff', '.d']
COMPRESSION_SUFFIXES = ['.zip', '.gz', '.lzma', '.bz2', '.7z', '.xz', '.rar']
COMPRESSED_MS_FILE_EXTENSIONS = [e + c for e in MS_FILE_EXTENSIONS for c in COMPRESSION_SUFFIXES]
_MS_FILE_RE = re.compile(r'\.(%s)(%s)?$' % ('|'.join(e[1:] for e in MS_FILE_EXTENSIONS),
                                            '|'.join(re.escape(c) for c in COMPRESSION_SUFFIXES)), re.I)


def prefer_extension(group):
    # An uncompressed file is always preferred over a compressed one, so assay rows that name
    # an uncompressed file are resolved exactly as before; among compressed files only, the
    # same rule is applied to the inner format (x.mzML.gz over x.raw.zip)
    candidates = ~group['compressed']
    if not candidates.any():
        candidates = pd.Series(True, index=group.index)
    preferred = candidates & group['extension'].isin(['mzml', 'mzxml'])
    if preferred.any():
        group['keep'] = preferred
    else:
        # Mark only the first candidate to keep if no preferred extension found
        group['keep'] = candidates & (candidates.cumsum() == 1)
    return group


def get_enviromental_water(x):
    x = x.lower()

    if 'water' in x or 'sewerage' in x:
        if 'waste' in x or 'sewerage' in x or 'sewage' in x:
            return ['environmental', 'water_waste']
        if 'surface' in x:
            return ['environmental', 'water_surface']
        if 'ground' in x:
            return ['environmental', 'water_ground']
        if 'storm' in x:
            return ['environmental', 'water_storm']
        if 'sea' in x or 'ocean' in x or 'coast' in x:
            return ['environmental', 'water_seawater']
        else:
            return [None, None]
    else:
        return [None, None]


_AGE_UNIT_DIVISOR = {'year': 1, 'yr': 1, 'month': 12, 'mo': 12, 'week': 52.1429, 'wk': 52.1429, 'day': 365, 'd': 365, 'hour': 8760, 'h': 8760}


def _age_in_years_from_text(value):
    """Age in years from a value that carries its own unit ('73 year', '7 weeks'), else None."""
    m = re.fullmatch(r'\s*(\d+(?:\.\d+)?)\s*(years?|yrs?|months?|mo|weeks?|wks?|days?|d|hours?|h)\s*(?:old)?\s*', str(value), re.I)
    if not m:
        return None
    unit = m.group(2).lower()
    unit = unit if unit in _AGE_UNIT_DIVISOR else unit.rstrip('s')
    return float(m.group(1)) / _AGE_UNIT_DIVISOR[unit]


def _sample_type_qc(x):
    """REDU SampleType for a MetaboLights 'Sample type' value that marks a pooled QC or a
    blank, otherwise None."""
    x = re.sub(r'[^a-z]+', ' ', str(x).lower()).strip()
    if 'pool' in x and ('quality control' in x or re.search(r'\bqc\b', x)):
        return 'pool_QC'
    if 'blank' in x:
        if re.search(r'extraction|process|procedur|method', x):
            return 'blank_extraction'
        if re.search(r'media|medium', x):
            return 'blank_culturemedia'
        if re.search(r'solvent|instrument|analytic|injection|^blank$', x):
            return 'blank_analysis'
    return None


def get_blanks(x):
    x = x.lower()
    if 'blank' in x:
        if 'blank' == x:
            return ['blank_analysis', 'blank_analysis']
        if 'method' in x or 'extraction' in x:
            return ['blank_extraction', 'blank_extraction']
        if 'media' in x:
            return ['blank_culturemedia', 'blank_culturemedia']
        else:
            return [None, None]
    else:
        return [None, None]


def _label(value):
    """Lower case, punctuation as spaces: 'Pooled QC-sample' -> 'pooled qc sample'."""
    return re.sub(r'[^a-z0-9]+', ' ', str(value).lower()).strip()


# MetaboLights labels for QC samples, blanks and standards (exact, after _label). They are
# read from the 'Sample type' characteristic, from factor values (submitters write
# 'Treatment = QC') and from the organism column ('reference compound', 'solvent blank').
_QC_BLANK_LABELS = {}
for _term, _labels in {
        'pool_QC': ['pooled quality control sample', 'pooled quality control samples', 'pooled quality control',
                    'pooled qc', 'pooled qc sample', 'pooled qc samples', 'qc pool', 'pool qc', 'pooled sample',
                    'pooled samples', 'pool', 'pooled', 'study reference', 'long term reference',
                    'quality control', 'quality control sample', 'quality control samples', 'qc', 'qc sample',
                    'qc samples'],
        'blank_analysis': ['blank', 'blanks', 'blank sample', 'blank samples', 'solvent blank', 'solvent blanks',
                           'instrument blank', 'injection blank', 'analytical blank', 'water blank'],
        'blank_extraction': ['extraction blank', 'extraction blanks', 'blank extraction', 'process blank',
                             'process blanks', 'processing blank', 'procedural blank', 'procedure blank',
                             'method blank', 'sample preparation blank', 'preparation blank'],
        'blank_culturemedia': ['media blank', 'medium blank', 'culture medium blank', 'culture media blank',
                               'blank medium', 'blank media', 'growth medium blank'],
        'molecular_standards': ['reference compound', 'reference compounds', 'reference compound mix',
                                'reference compound mixture', 'reference compound single', 'standard',
                                'standards', 'standard mix', 'standard mixture', 'reference standard',
                                'reference standards', 'chemical standard', 'authentic standard',
                                'internal standard mixture', 'calibration curve', 'calibration point',
                                'calibration standard', 'calibration standards'],
        }.items():
    for _l in _labels:
        _QC_BLANK_LABELS[_l] = _term
_SUBTYPE = {'molecular_standards': 'defined_standards'}
# organism values that are not organisms but say what the sample is
_ORGANISM_BLANK_LABELS = {'experimental blank': 'blank_analysis', 'solvent': 'blank_analysis', 'acetonitrile': 'blank_analysis',
                          'methanol': 'blank_analysis', 'non sample control': 'blank_analysis',
                          'compound': 'molecular_standards', 'compounds': 'molecular_standards', 'chemical': 'molecular_standards'}

# Organism part values that name an UBERON term without being its label or exact synonym
_BODY_PART_ALIASES = {
    'plasma': 'blood plasma', 'plasma sample': 'blood plasma', 'plasma samples': 'blood plasma',
    'edta plasma': 'blood plasma', 'heparin plasma': 'blood plasma', 'human plasma': 'blood plasma',
    'serum': 'blood serum', 'serum sample': 'blood serum', 'serum samples': 'blood serum', 'human serum': 'blood serum',
    'whole blood': 'blood', 'blood sample': 'blood', 'blood samples': 'blood',
    'stool': 'feces', 'stool sample': 'feces', 'stool samples': 'feces', 'faeces': 'feces', 'fecal': 'feces',
    'faecal': 'feces', 'feces sample': 'feces', 'feces samples': 'feces', 'fecal sample': 'feces',
    'fecal samples': 'feces', 'faecal sample': 'feces', 'faeces sample': 'feces',
    'urine sample': 'urine', 'urine samples': 'urine', 'csf': 'cerebrospinal fluid',
}


def _qc_blank_term(value):
    """ReDU SampleType for an exact MetaboLights QC/blank/standard label, else None."""
    return _QC_BLANK_LABELS.get(_label(value))


# Columns holding the country where the sample was taken. 'Country of origin' and 'Country
# of collection' are left to the exact-name match above: in MetaboLights they describe the
# provenance of germplasm accessions (cultivar/landrace origin), not where the analysed
# sample was collected.
_COUNTRY_CODE_COLUMNS = re.compile(r'Samples_(?:country|country of sampling|sampling country|geographic(?:al)? location)(?:_\d+)?$', re.I)


def _country_term(value, vocab, aliases):
    """Country term for a country name/code ('USA', 'NZ', 'ITA', 'USA (Oregon)'), else None."""
    key = str(value).strip().lower()
    if not key or key in _MISSING_VALUES:
        return None
    for k in (key, re.sub(r'\s*\(.*?\)\s*$', '', key).strip()):
        if k in vocab:
            return vocab[k]
        if k in aliases:
            return aliases[k]
    return None


_UNIT_YEARS = {'year': 1, 'month': 12, 'week': 52.1429, 'day': 365, 'hour': 8760}


def _unit_name(unit):
    """'Month', 'weeks', 'yr' -> 'month', 'week', 'year'; None if not a time unit."""
    u = str(unit).strip().lower()
    u = {'yr': 'year', 'yrs': 'year', 'y': 'year', 'years old': 'year', 'mo': 'month', 'wk': 'week', 'wks': 'week',
         'd': 'day', 'h': 'hour', 'hr': 'hour', 'hrs': 'hour'}.get(u, u)
    u = u[:-1] if u.endswith('s') and u[:-1] in _UNIT_YEARS else u
    return u if u in _UNIT_YEARS else None


# Age of the sampled subject, other than 'Age' itself (which is read first)
_OTHER_AGE_COLUMNS = re.compile(r'Samples_(?:mouse age|rat age|animal age|subject age|patient age|age at sampling|age at sacrifice'
                                r'|age at collection|age at sample collection|exact age|current chronological age|chronological age'
                                r'|postnatal age|age in (?:years|months|weeks|days))(?:_\d+)?$', re.I)


def _age_years(value, unit='', header=''):
    """Age in years from a number with its unit (unit column, in the value, or in the header)."""
    s = str(value).strip()
    m = re.fullmatch(r'(\d+(?:\.\d+)?)\s*([A-Za-z]*)\s*(?:old)?', s)
    if not m:
        return None
    u = _unit_name(m.group(2)) if m.group(2) else None
    if m.group(2) and u is None:
        return None
    if u is None:
        u = _unit_name(unit) if str(unit).strip() else None
    if u is None:
        h = re.search(r'\bin (years|months|weeks|days)\b', str(header), re.I)
        u = _unit_name(h.group(1)) if h else None
    if u is None:
        return None
    return float(m.group(1)) / _UNIT_YEARS[u]


_LENGTH_M = {'meter': 1, 'meters': 1, 'metre': 1, 'metres': 1, 'm': 1, 'centimeter': 0.01, 'centimeters': 0.01,
             'centimetre': 0.01, 'cm': 0.01, 'millimeter': 0.001, 'millimeters': 0.001, 'mm': 0.001,
             'kilometer': 1000, 'kilometers': 1000, 'km': 1000}


def _meters(value, unit):
    """Number with a length unit in meters, else None."""
    try:
        x = float(str(value).strip())
    except ValueError:
        return None
    f = _LENGTH_M.get(str(unit).strip().lower())
    return None if f is None or not np.isfinite(x) else x * f


_SUBJECT_COLUMNS = re.compile(r'Samples_(?:subject|subject id|subject_id|subject identifier|subject number|patient|patient id'
                              r'|patient identifier|patient number|participant|participant id|participant identifier'
                              r'|study participant|animal id|animal identifier|animal number|animal name|mouse id'
                              r'|mouse identifier|rat id|donor|donor id|sample donor|stool donor|individual'
                              r'|individual id|person|volunteer|volunteer id|tree id|plant id|host id)(?:_\d+)?$', re.I)
# values in subject columns that mark QC/blank/pooled samples rather than a subject
_NOT_A_SUBJECT = re.compile(r'\b(?:qc|quality control|pool(?:ed)?|blank|standard|reference|mix(?:ture)?|not applicable|none|na'
                            r'|control|healthy|patient|disease|accession)\b', re.I)

_DATE_COLUMNS = re.compile(r'Samples_(?:collection date|sample collection date|sampling date|sample date|date of sampling'
                           r'|date of collection|date of sample collection|sampling datetime|collection datetime'
                           r'|sampling date and time|collection date and time|utc date and time|utc sampling time'
                           r'|harvest date|harvest year|harvesting date|year of harvest|year of harvesting'
                           r'|year of collection|collection year|sampling year|sampling year date)(?:_\d+)?$', re.I)


def _collection_dates(series):
    """ISO-8601 collection dates (precision kept) for one column, None where not parseable.
    Dotted dates are day-first (d.m.Y). Slash/dash day-month dates are only used when the
    column itself shows the order (a field > 12); otherwise they are left out."""
    values = series.astype(str).str.strip()
    slash = values.str.extract(r'^(\d{1,2})[/-](\d{1,2})[/-]\d{2,4}(?:\s|T|$)').dropna().astype(int)
    dmy = bool(((slash[0] > 12) & (slash[1] <= 12)).any())
    mdy = bool(((slash[1] > 12) & (slash[0] <= 12)).any())

    def one(v):
        if not v or v.lower() in _MISSING_VALUES:
            return None
        m = re.fullmatch(r'(\d{1,2})\.(\d{1,2})\.(\d{4})', v)
        if m:
            v = f'{m.group(3)}-{int(m.group(2)):02d}-{int(m.group(1)):02d}'
        elif re.match(r'\d{1,2}[/-]\d{1,2}[/-]\d', v):
            if dmy == mdy:
                return None
            return normalize_redu_datetime(v, dayfirst=dmy)
        return normalize_redu_datetime(v)
    cache = {}
    return values.map(lambda v: cache.setdefault(v, one(v)))


# ---------------------------------------------------------------------------------------
# SampleExtractionMethod from the 'Extraction' protocol text. Only a solvent system written
# with its ratio ('methanol:water (4:1, v/v)', '80% methanol', '4:4:2 ACN:MeOH:H2O') that is
# an allowed term, or a single named solvent with an exact '(100%)' term, is used. The text
# is cut where the extract is dried and from reconstitution/LC/derivatisation keywords on;
# protocols naming other solvents, several systems or per-sample-type procedures are skipped.
# ---------------------------------------------------------------------------------------
_EXT_SOLVENTS = [
    ('methyltertbutylether', r'methyl[\s-]*tert(?:iary)?[\s-]*butyl[\s-]*ether|tert[\s-]*butyl[\s-]*methyl[\s-]*ether|MTBE'),
    ('acetonitrile', r'acetonitrile|ACN|MeCN|CH3CN'),
    ('isopropanol', r'isopropanol|iso-propanol|2-propanol|isopropyl[\s-]+alcohol|IPA'),
    ('dichloromethane', r'dichloromethane|DCM|CH2Cl2|methylene[\s-]+chloride'),
    ('ethyl acetate', r'ethyl[\s-]+acetate|EtOAc'),
    ('chloroform', r'chloroform|CHCl3'),
    ('methanol', r'methanol|MeOH|CH3OH'),
    ('ethanol', r'ethanol|EtOH'),
    ('butanol', r'n-butanol|1-butanol|butanol|BuOH'),
    ('hexane', r'n-hexane|hexane'),
    ('acetone', r'acetone'),
    ('water', r'water|H2O|dd\s?H2O|Milli-?Q|ultrapure\s+water'),
    # other organic solvents: not in the vocabulary, only block a match
    ('other', r'petroleum\s+ether|diethyl\s+ether|heptane|pentane|toluene|cyclohexane|benzene|iso-?octane|n-?propanol|1-propanol|pyridine|tetrahydrofuran|THF|DMSO|dimethyl\s+sulfoxide|formamide'),
]
_EXT_SOLV_ALT = '|'.join('(?:%s)' % p for _, p in _EXT_SOLVENTS)
_EXT_SOLV_RE = re.compile(r'(?<![A-Za-z0-9])(?<!\d-)(?:%s)(?![A-Za-z0-9])' % _EXT_SOLV_ALT, re.I)
def _ext_solvent(token):
    for name, p in _EXT_SOLVENTS:
        if re.fullmatch(p, token.strip(), re.I):
            return name
    return None
_EXT_NUM = r'\d+(?:\.\d+)?'
_EXT_SEP = r'\s*(?:[:/]|-(?!\d)|\band\b|,)\s*'
_EXT_RSEP = r'\s*[:/]\s*'
_EXT_S = r'(%s)' % _EXT_SOLV_ALT
# 'methanol:water (4:1, v/v)', 'MeOH/H2O 80:20', 'acetonitrile:methanol:water (2:2:1)'
_EXT_NAMES_RATIO = re.compile(r'(?<![A-Za-z0-9])%s%s%s(?:%s%s)?(?![A-Za-z0-9])\s*[\(\[,]?\s*(%s)%s(%s)(?:%s(%s))?(?!\s*(?:%%|mm|µ|u[lLm]|m[lLgM]|°))' % (_EXT_S, _EXT_SEP, _EXT_S, _EXT_SEP, _EXT_S, _EXT_NUM, _EXT_RSEP, _EXT_NUM, _EXT_RSEP, _EXT_NUM), re.I)
# '4:1 methanol:water', '80:20 (v/v) MeOH/H2O'
_EXT_RATIO_NAMES = re.compile(r'(?<![\d.:/])(%s)%s(%s)(?:%s(%s))?\s*(?:\(?v/v(?:/v)?\)?,?\s*)?(?:of\s+)?%s%s%s(?:%s%s)?(?![A-Za-z0-9])' % (_EXT_NUM, _EXT_RSEP, _EXT_NUM, _EXT_RSEP, _EXT_NUM, _EXT_S, _EXT_SEP, _EXT_S, _EXT_SEP, _EXT_S), re.I)
# '80% methanol', '80 % (v/v) aqueous methanol', 'methanol (80%)'
_EXT_PCT = re.compile(r'(?<![\d.])(%s)\s*%%\s*(?:\(?\s*v\s*[:/]\s*v\s*\)?\s*)?(?:(?:aqueous|aq\.|cold|ice[\s-]cold|pre-?chilled|pre-?cooled|LC-?MS[\s-]grade|HPLC[\s-]grade)\s*)*%s(?![A-Za-z0-9])(?!\s*(?:/|:|-(?!\d)|and)\s*(?:%s))' % (_EXT_NUM, _EXT_S, _EXT_SOLV_ALT), re.I)
_EXT_PCT_AFTER = re.compile(r'(?<![A-Za-z0-9])%s\s*\(\s*(%s)\s*%%\s*(?:,?\s*v/v)?\s*\)' % (_EXT_S, _EXT_NUM), re.I)
# sentences about what happens after/around the extraction, not the extraction solvent
_EXT_MASK = re.compile(r'reconstitut|re-?dissol|dissolved\s+in|mobile\s+phase|eluent|elut|solvent\s+[AB]\b|buffer\s+[AB]\b|dilut|inject|derivati|stock|calibra|SPE\b|solid[\s-]phase|cartridge|column', re.I)
_EXT_NOT_SOLVENT_WATER = re.compile(r'water[\s-]*bath|ice[\s-]*water|ice/water|water[\s-]*cooled|tap\s+water|drinking\s+water|water\s+content', re.I)

_EXT_PER_TYPE = re.compile(r'(?:^|[.;:]\s+|\)\s+)For\s+(?:[\w-]+\s+){0,4}?[\w-]+\s*,|\b[\w-]+\s+(?:protocol|samples?)\s*:', re.I)
_EXT_DRY = re.compile(r'(?:evaporated|dried|dry)\s+(?:down\s+|completely\s+|it\s+)?(?:under|in\s+(?:a\s+)?(?:vacuo|vacuum|speed|centrifugal|nitrogen|N2)|using|with\s+(?:a\s+)?(?:nitrogen|N2|speed|vacuum|centrifugal)|by\s+(?:nitrogen|vacuum|speed|centrifugal|lyophili))|to\s+(?:complete\s+)?dryness|until\s+(?:complete\s+)?dry|speed-?vac|vacuum\s+concentrat|centrifugal\s+(?:evaporat|concentrat|vacuum)|were\s+(?:then\s+)?(?:dried|evaporated)\b|was\s+(?:then\s+)?(?:dried|evaporated)\b', re.I)
def _ext_clean(text):
    text = html.unescape(re.sub(r'<[^>]+>', ' ', str(text)))
    text = re.sub(r'\s+', ' ', text).strip()
    # literature references cited at the end of the protocol
    return re.split(r'\bRefs?\s*:|\bReferences?\s*:', text)[0]

def _ext_sentences(text):
    return [s for s in re.split(r'(?<=[.;])\s+(?=[A-Z0-9(])', text) if s.strip()]

def _ext_canon_ratio(names, ratio):
    try:
        vals = [float(x) for x in ratio]
    except ValueError:
        return None
    if any(v <= 0 for v in vals) or any(v != int(v) for v in vals):
        return None
    vals = [int(v) for v in vals]
    g = 0
    for v in vals: g = math.gcd(g, v)
    return names, [v // g for v in vals]

def _ext_lookup(names, ratio, allowed):
    """Allowed term for solvents+ratio in any order, else None."""
    if len(set(names)) != len(names):
        return None
    for perm in itertools.permutations(range(len(names))):
        term = '-'.join(names[i] for i in perm) + ' (' + ':'.join(str(ratio[i]) for i in perm) + ')'
        if term in allowed:
            return term
    return None

def _extraction_systems(text):
    """(solvent-system mentions [(names, ratio)], solvents mentioned) in the unmasked sentences."""
    systems, mentioned = [], set()
    text = _ext_clean(text)
    # the extraction solvent is named before the extract is dried; what follows is reconstitution
    dm = _EXT_DRY.search(text)
    if dm:
        text = text[:dm.start()]
    for sent in _ext_sentences(text):
        # drop the part of the sentence from a reconstitution/LC/derivatisation keyword on
        mm = _EXT_MASK.search(sent)
        if mm:
            sent = sent[:mm.start()]
        s = _EXT_NOT_SOLVENT_WATER.sub(' ', sent)
        for m in _EXT_SOLV_RE.finditer(s):
            mentioned.add(_ext_solvent(m.group(0)))
        spans = []
        for rx, kind in ((_EXT_NAMES_RATIO, 'nr'), (_EXT_RATIO_NAMES, 'rn')):
            for m in rx.finditer(s):
                g = m.groups()
                if kind == 'nr':
                    names = [x for x in g[:3] if x]; ratio = [x for x in g[3:] if x]
                else:
                    ratio = [x for x in g[:3] if x]; names = [x for x in g[3:] if x]
                if len(names) != len(ratio):
                    systems.append((None, m.group(0))); spans.append(m.span()); continue
                c = _ext_canon_ratio([_ext_solvent(n) for n in names], ratio)
                systems.append((c, m.group(0)) if c else (None, m.group(0))); spans.append(m.span())
        for rx in (_EXT_PCT, _EXT_PCT_AFTER):
            for m in rx.finditer(s):
                if any(a <= m.start() < b for a, b in spans):
                    continue
                pct, name = (m.group(1), m.group(2)) if rx is _EXT_PCT else (m.group(2), m.group(1))
                name = _ext_solvent(name); pct = float(pct)
                if name == 'water' or pct <= 0 or pct > 100 or pct != int(pct):
                    systems.append((None, m.group(0))); continue
                if pct == 100:
                    systems.append((([name], [100]), m.group(0)))
                else:
                    systems.append((_ext_canon_ratio([name, 'water'], [pct, 100 - pct]), m.group(0)))
    return systems, mentioned

def _extraction_method(text, allowed):
    """Allowed SampleExtractionMethod term if the extraction protocol names exactly one
    solvent system that is in the vocabulary and no other solvent; else None."""
    if not text or not str(text).strip():
        return None, 'empty'
    # several sample-type specific procedures ('For serum samples, ... For fecal samples, ...')
    if len(_EXT_PER_TYPE.findall(_ext_clean(text))) >= 2:
        return None, 'per_sample_type'
    systems, mentioned = _extraction_systems(text)
    mentioned.discard(None)
    terms = set()
    for c, raw in systems:
        if c is None:
            return None, 'unparsed:' + raw
        names, ratio = c
        if ratio == [100]:
            t = names[0] + ' (100%)'
            terms.add(t if t in allowed else 'NOVOCAB:' + t)
        else:
            t = _ext_lookup(names, ratio, allowed)
            terms.add(t or 'NOVOCAB:' + '-'.join(names) + ' (' + ':'.join(map(str, ratio)) + ')')
    if not systems:
        # a single solvent named without ratio: only pure solvents with an exact '(100%)' term
        if len(mentioned) == 1 and mentioned != {'water'} and mentioned != {'other'}:
            t = next(iter(mentioned)) + ' (100%)'
            if re.search(r'\b(?:RNA|DNA|TRIzol)\b', text) or re.search(r'\b(?:hot|boiling|heated)\s+(?:\w+\s+)?(?:ethanol|methanol|water)', _ext_clean(text), re.I):
                return None, 'single_context:' + t
            # a ratio or percentage next to the solvent that was not understood
            for sent in _ext_sentences(_ext_clean(text)):
                mm = _EXT_MASK.search(sent)
                sent = sent[:mm.start()] if mm else sent
                if _EXT_SOLV_RE.search(sent) and re.search(r'\d\s*%|\d\s*[:/]\s*\d', sent):
                    return None, 'single_with_numbers:' + t
            return (t, 'single') if t in allowed else (None, 'single_novocab:' + t)
        return None, 'none' if not mentioned else 'solvents_no_ratio:' + ','.join(sorted(mentioned))
    if len(terms) != 1:
        return None, 'multiple:' + ';'.join(sorted(terms))
    t = next(iter(terms))
    if t.startswith('NOVOCAB:'):
        return None, t
    term_solvents = set(re.split(r'-(?=[a-z])', t.split(' (')[0]))
    if mentioned != term_solvents:
        return None, 'other_solvents:' + t + ' vs ' + ','.join(sorted(mentioned))
    return t, 'ratio'


def _internal_standards(protocol_texts, allowed_terms):
    """InternalStandardsUsed term for the vocabulary compounds named in a sentence about
    internal standards ('ribitol was added as internal standard'), else None. Several
    compounds are only used when their combination is itself an allowed term."""
    singles = [t for t in allowed_terms if ';' not in t and t != 'none']
    found = set()
    for text in protocol_texts:
        for sent in _ext_sentences(_ext_clean(text)):
            if not re.search(r'internal\s+standard', sent, re.I):
                continue
            for t in singles:
                for m in re.finditer(r'(?<![A-Za-z0-9-])' + re.escape(t) + r'(?![A-Za-z0-9-])', sent, re.I):
                    # an isotope-labelled form ('deuterated biotin', '[13C5]ribitol') is another compound
                    after = sent[m.end():m.end() + 60]
                    if after.lstrip().startswith('('):
                        after = after[:after.index(')') + 1] if ')' in after else after
                    else:
                        after = after[:15]
                    around = sent[max(0, m.start() - 25):m.start()] + ' ' + after
                    if not re.search(r'deuter|label|isotop|13C|15N|\bD\d|-d\d|\(d\d|\bd\d+\b|\]\s*$', around, re.I):
                        found.add(t)
    if not found:
        return None
    if len(found) == 1:
        return next(iter(found))
    combos = [t for t in allowed_terms if ';' in t and set(t.split(';')) == found]
    return combos[0] if len(combos) == 1 else None


def load_metabolights_sheets(transSheet_dir, allowedTerm_dict):
    """Optional MetaboLights translation sheets: explicit instrument aliases and country
    codes/names. Rows whose target is not an allowed term are ignored."""
    sheets = {}
    path = os.path.join(transSheet_dir, 'MTBLS_MassSpectrometer_aliases.csv')
    if os.path.exists(path):
        t = pd.read_csv(path, dtype=str)
        t = t[t['MassSpectrometer'].isin(set(allowedTerm_dict['MassSpectrometer']['allowed_values']))]
        sheets['ms_aliases'] = dict(zip(t['ML'].str.strip(), t['MassSpectrometer']))
    path = os.path.join(transSheet_dir, 'MTBLS_Country_aliases.csv')
    if os.path.exists(path):
        t = pd.read_csv(path, dtype=str, keep_default_na=False)
        t = t[t['Country'].isin(set(allowedTerm_dict['Country']['allowed_values']))]
        sheets['country_aliases'] = dict(zip(t['alias'].str.strip().str.lower(), t['Country']))
    path = os.path.join(transSheet_dir, 'MTBLS_Organism_aliases.csv')
    if os.path.exists(path):
        t = pd.read_csv(path, dtype=str, keep_default_na=False)
        known = set(allowedTerm_dict['NCBITaxonomy']['allowed_values'])
        t = t[t['NCBITaxonomy'].isin(known)]
        sheets['organism_aliases'] = dict(zip(t['organism'].str.strip().str.lower(), t['NCBITaxonomy']))
    # study diseases: exact-meaning aliases of DOID terms (empty target = never map), DOID is_a edges
    path = os.path.join(transSheet_dir, 'MTBLS_disease_aliases.csv')
    if os.path.exists(path):
        t = pd.read_csv(path, dtype=str, keep_default_na=False)
        sheets['disease_aliases'] = dict(zip(t['alias'].str.strip().str.lower(), t['DOIDCommonName'].str.strip()))
    # health conditions and exposures: exact-meaning aliases of HPO / ECTO / PECO terms
    # (column HPOCommonName or Exposure; empty term = never map)
    # standard disease-model procedures (only where the study text names the procedure)
    path = os.path.join(transSheet_dir, 'MTBLS_procedure_aliases.csv')
    if os.path.exists(path):
        t = pd.read_csv(path, dtype=str, keep_default_na=False)
        t = t[t['column'].str.strip() == 'HPOCommonName']
        sheets['procedure_aliases'] = [(_dlabel(p), a.strip().lower(), x.strip()) for p, a, x in zip(t['procedure'], t['abbreviation'], t['term'])]
    path = os.path.join(transSheet_dir, 'MTBLS_condition_aliases.csv')
    if os.path.exists(path):
        t = pd.read_csv(path, dtype=str, keep_default_na=False)
        for column, key in (('HPOCommonName', 'hpo_aliases'), ('Exposure', 'exposure_aliases')):
            r = t[t['column'].str.strip() == column]
            sheets[key] = dict(zip(r['alias'].map(_dlabel), r['term'].str.strip()))
    return sheets


class StudyNotAvailable(Exception):
    """The study exists but is not public (MetaboLights returns no content for it)."""


def safe_api_request(url, retries=3, expected_codes={200}):

  for _ in range(retries):
    try:
      response = EBI_SESSION.get(url, timeout=(10, 300))
      print(f"Requested: {url}")
      if response.status_code in expected_codes:
        data = response.json()
        print("Request successful!")
        return data
      elif response.status_code in {401, 403, 404}:
        # private or unknown study: retrying will not help
        print(f"Error: status code {response.status_code}, not retrying")
        return {'content': None, 'message': f'HTTP {response.status_code}'}
      else:
        print(f"Error: Unexpected status code {response.status_code}")
    except Exception as e:
      print(f"Error requesting data: {e}")
    time.sleep(10)
  print(f"All retries failed for {url}.")
  return None

# ---------------------------------------------------------------------------------------------
# Study record when the per-study API is down
#
# The converter reads the study record of 'ws/studies/public/study/<MTBLS>' (assays with their
# technology and assay table, the sample table, protocols, submission year, title, description,
# design descriptors and factors). When that endpoint fails with server errors or timeouts
# (EBI maintenance: HTTP 503 for every study), the same structure is built from the study's
# ISA-Tab metadata, taken from
#   1. the API's working copy ('ws/studies/<MTBLS>/investigation' and
#      'ws/studies/<MTBLS>/<s_ or a_ file>'), which is what the study record is made from, else
#   2. the public FTP folder (i_Investigation.txt, s_ and a_ files over HTTPS), which can lag
#      behind the working copy (assay sheets naming files EBI has compressed since, ...).
# Headers follow the API: 'Characteristics[Organism part]' -> header 'Organism part' with
# fieldType 'Characteristics' (also 'Factor Value', 'Parameter Value', 'Comment'; a repeated
# one keeps its name), other columns keep their name with the pandas suffix for repeats
# ('Term Source REF.1', 'Unit.2'), and the pandas NA words ('NA', 'None', 'n/a') become ''.
# ---------------------------------------------------------------------------------------------
METABOLIGHTS_WS = 'https://www.ebi.ac.uk/metabolights/ws/studies/'
METABOLIGHTS_FTP = 'https://ftp.ebi.ac.uk/pub/databases/metabolights/studies/public/'
# consecutive failures (server errors / timeouts) after which a source is taken as down and
# skipped; while down, it is tried again once every REPROBE_EVERY studies
SOURCE_DOWN_AFTER = 2
REPROBE_EVERY = 100
_SOURCE_STATE = {name: {'consecutive_failures': 0, 'down': False, 'since_probe': 0} for name in ('api', 'workcopy')}
# 'Kind[name]' headers (anything after the first ']' is ignored, as in the API: 'Factor
# Value[Storage location]]', the pandas suffix of repeats '.1')
_ISA_BRACKET_HEADER = re.compile(r'([^\[\]]*?)\s*\[([^\]]*)\].*', re.S)
_ISA_CANONICAL_KINDS = ('Characteristics', 'Factor Value', 'Parameter Value', 'Comment')
_ISA_IN_VENDOR_FOLDER = re.compile(r'(.*\.(?:d|raw))/[^/]*', re.I)
_MS_TECHNOLOGIES = ('mass spectrometry', 'mass spectrometry assay')


class SourceUnavailable(Exception):
    """A metadata source answered with server errors or timed out."""


def _ebi_get(url, retries=3, pause=5):
    """GET with the shared session: the response for 200, None for 401/403/404; raises
    SourceUnavailable after repeated server errors or timeouts."""
    last = None
    for attempt in range(retries):
        try:
            response = EBI_SESSION.get(url, timeout=(10, 300))
            if response.status_code == 200:
                return response
            if response.status_code in (401, 403, 404):
                return None
            last = f'HTTP {response.status_code}'
        except requests.RequestException as e:
            last = f'{type(e).__name__}: {e}'
        if attempt + 1 < retries:
            time.sleep(pause * (attempt + 1))
    raise SourceUnavailable(f'{url}: {last}')


def _decode(content):
    try:
        return content.decode('utf-8-sig')
    except UnicodeDecodeError:
        return content.decode('latin-1')


def fetch_workcopy_investigation(study_id):
    """i_Investigation.txt of the API's working copy (text), None if not available."""
    response = _ebi_get(METABOLIGHTS_WS + study_id + '/investigation')
    return None if response is None else _decode(response.content)


def fetch_workcopy_table(study_id, file_name):
    """An s_ / a_ table of the API's working copy ({'header': {name: index}, 'data': {'rows':
    [...]}}), None if not available."""
    response = _ebi_get(METABOLIGHTS_WS + study_id + '/' + quote(file_name))
    return None if response is None else response.json()


def fetch_isa_file(study_id, file_name):
    """Text of one ISA-Tab file from the public FTP folder (over HTTPS), None if absent."""
    response = _ebi_get(METABOLIGHTS_FTP + study_id + '/' + quote(file_name))
    return None if response is None else _decode(response.content)


def _parse_investigation(text):
    """'Study ...' rows of an i_Investigation.txt (the STUDY part): key -> list of values.
    Quoted values may span lines; trailing empty values are dropped."""
    out, in_study = {}, False
    for row in csv.reader(io.StringIO(text), delimiter='\t', strict=False):
        if not row:
            continue
        key = row[0].strip()
        if key == 'STUDY':
            in_study = True
            continue
        if in_study and key and key not in out:
            values = [v.strip() for v in row[1:]]
            while values and values[-1] == '':
                values.pop()
            out[key] = values
    return out


def _table_record(raw_headers, rows):
    """{'fields', 'data'} as in the API study record from the raw column names (pandas
    suffixes '.1', ... on repeats) and the rows (lists of values in column order), with the
    API's reading of ISA-Tab: a repeated 'Characteristics[...]' etc. keeps its name; a
    bracketed name with 'comment' in it is a Comment ('Factor Value[Comments]'); other
    bracketed names count only with the exact kind ('Factor value[Sex]', 'Factor[Treatment]'
    are left out); columns without a name are left out; rows whose first value starts with
    '#' are comment lines."""
    fields, keep = {}, []
    for i, raw in enumerate(raw_headers):
        if not isinstance(raw, str) or not raw.strip():
            continue  # column without a header
        m = _ISA_BRACKET_HEADER.fullmatch(raw.strip())
        if m:
            kind, header = m.group(1).strip(), m.group(2)
            if 'comment' in raw.lower():
                kind = 'Comment'
            elif kind not in _ISA_CANONICAL_KINDS:
                continue
            raw = f'{kind}[{header}]'
        else:
            header, kind = raw, 'basic'
        fields[f'{i}~{raw.lower()}'] = {'index': i, 'header': header, 'fieldType': kind, 'description': '',
                                       'cleanHeader': header}
        keep.append(i)
    rows = [row for row in rows if not (row and str(row[0]).lstrip().startswith('#'))]
    return {'fields': fields, 'data': [[row[i] for i in keep] for row in rows]}


def _isa_table(text):
    """An ISA-Tab table file (FTP) as the API record gives it."""
    df = pd.read_csv(io.StringIO(text), sep='\t', dtype=str, header=None)
    if len(df) == 0:
        return {'fields': {}, 'data': []}
    raw, seen = [], {}
    for h in df.iloc[0].tolist():
        if isinstance(h, str) and h.strip():
            n = seen.get(h, 0)  # pandas-style suffix of repeated names, as in the working copy
            seen[h] = n + 1
            h = h if n == 0 else f'{h}.{n}'
        raw.append(h)
    return _table_record(raw, df.iloc[1:].fillna('').values.tolist())


_NA_WORDS = {'#N/A', '#N/A N/A', '#NA', '-1.#IND', '-1.#QNAN', '-NaN', '-nan', '1.#IND', '1.#QNAN', '<NA>', 'N/A',
             'NA', 'NULL', 'NaN', 'None', 'n/a', 'nan', 'null'}


def _ws_table(js):
    """An s_ / a_ table of the working copy (JSON) as the API record gives it."""
    header = (js or {}).get('header') or {}
    # sample tables give name -> index, assay tables name -> {'index': ..., ...}
    position = lambda v: v.get('index', 0) if isinstance(v, dict) else v
    names = [n for n, _ in sorted(header.items(), key=lambda kv: position(kv[1]))]
    rows = sorted(((js or {}).get('data') or {}).get('rows') or [], key=lambda r: r.get('index', 0))
    value = lambda v: '' if v is None or str(v) in _NA_WORDS else str(v)
    return _table_record(names, [[value(r.get(n)) for n in names] for r in rows])


def _record_from_isa(study_id, inv_text, get_table):
    """The study record ({'content': ...}) from the investigation text and a function giving
    the API-like table of an s_ / a_ file name (None if missing)."""
    inv = _parse_investigation(inv_text)
    first = lambda key: (inv.get(key) or [''])[0]
    at = lambda key, i: (inv.get(key) or [])[i] if i < len(inv.get(key) or []) else ''

    year = None
    for key in ('Study Submission Date', 'Study Public Release Date'):
        m = re.search(r'(19|20)\d{2}', first(key))
        if m:
            year = int(m.group(0))
            break
    names = inv.get('Study Protocol Name') or []
    protocols = [{'name': n, 'description': at('Study Protocol Description', i)} for i, n in enumerate(names)]
    types = inv.get('Study Design Type') or []
    descriptors = [{'description': (f"{at('Study Design Type Term Source REF', i)}:{t}"
                                    if at('Study Design Type Term Source REF', i) else t)}
                   for i, t in enumerate(types) if t]
    factors = [{'name': n} for n in inv.get('Study Factor Name') or [] if n]

    sample_file = first('Study File Name')
    sample_table = get_table(sample_file) if sample_file else None
    if sample_table is None:
        raise RuntimeError(f'{study_id}: sample file {sample_file!r} not available')
    assays = []
    for i, name in enumerate(inv.get('Study Assay File Name') or []):
        if not name:
            continue
        technology = at('Study Assay Technology Type', i)
        table = {'fields': {}, 'data': []}
        # only mass spectrometry assays are read by the converter
        if technology.strip().lower() in _MS_TECHNOLOGIES:
            table = get_table(name)
            if table is None:
                print(f'{study_id}: assay file {name} not available')
                continue
        assays.append({'measurement': at('Study Assay Measurement Type', i), 'technology': technology,
                       'platform': at('Study Assay Technology Platform', i), 'fileName': name,
                       'assayNumber': str(i + 1), 'assayTable': table})
    content = {'studyIdentifier': first('Study Identifier') or study_id, 'title': first('Study Title'),
               'description': first('Study Description'), 'studyDescription': None,
               'factors': factors, 'descriptors': descriptors, 'protocols': protocols,
               'derivedData': {'submissionYear': year}, 'assays': assays, 'sampleTable': sample_table}
    return {'content': content}


def study_details_from_workcopy(study_id):
    """The study record built from the API's working copy of the ISA-Tab files; None if the
    working copy has no investigation for the study."""
    inv_text = fetch_workcopy_investigation(study_id)
    if inv_text is None or not inv_text.strip():
        return None

    def table(name):
        js = fetch_workcopy_table(study_id, name)
        return None if js is None else _ws_table(js)
    details = _record_from_isa(study_id, inv_text, table)
    details['content']['source'] = 'API working copy (ISA-Tab)'
    return details


def study_details_from_isa(study_id):
    """The study record built from the public FTP ISA-Tab files; {'content': None, ...} if
    the study has no public investigation file."""
    inv_text = fetch_isa_file(study_id, 'i_Investigation.txt')
    if inv_text is None:
        return {'content': None, 'message': 'no public ISA-Tab investigation file'}

    def table(name):
        text = fetch_isa_file(study_id, name)
        if text is None:
            return None
        t = _isa_table(text)
        # the FTP copy of an assay sheet can still name a file inside a vendor folder
        # ('x.D/data.ms', 'x.raw/_FUNC001.DAT', 'x.D/') where the working copy (and the FTP
        # data) have the folder, often compressed since ('x.D.zip'): take the folder
        fields = sorted(t['fields'].values(), key=lambda f: f['index'])
        for j, f in enumerate(fields):
            if name.startswith('a_') and re.match(r'(?:Raw|Derived) Spectral Data File', f['header']):
                for row in t['data']:
                    m = _ISA_IN_VENDOR_FOLDER.fullmatch(str(row[j]))
                    if m:
                        row[j] = m.group(1)
        return t
    details = _record_from_isa(study_id, inv_text, table)
    details['content']['source'] = 'ISA-Tab files (EBI FTP)'
    return details


def _source_usable(name):
    """Whether to try a source for this study (always while up; one probe every
    REPROBE_EVERY studies while down), and the number of attempts."""
    s = _SOURCE_STATE[name]
    if not s['down']:
        return True, 3
    s['since_probe'] += 1
    if s['since_probe'] >= REPROBE_EVERY:
        s['since_probe'] = 0
        return True, 1
    return False, 0


def _source_result(name, ok):
    s = _SOURCE_STATE[name]
    if ok:
        if s['down']:
            print(f'MetaboLights {name} answers again; using it for the next studies')
        s.update(consecutive_failures=0, down=False, since_probe=0)
        return
    s['consecutive_failures'] += 1
    if not s['down'] and s['consecutive_failures'] >= SOURCE_DOWN_AFTER:
        s['down'] = True
        print(f'MetaboLights {name} failed for {SOURCE_DOWN_AFTER} studies in a row; skipping it '
              f'(probed again every {REPROBE_EVERY} studies)')


def get_study_details(study_id):
    """The study record from the per-study API; when it fails with server errors / timeouts
    (not for 401/403/404: private or unknown studies), from the API's working copy of the
    ISA-Tab files, else from the public FTP ISA-Tab files. A source that fails for
    SOURCE_DOWN_AFTER studies in a row is skipped (one probe every REPROBE_EVERY studies),
    so a down API does not cost 30 s per study. Requests are sequential, on one session."""
    use, retries = _source_usable('api')
    if use:
        details = safe_api_request(METABOLIGHTS_WS + 'public/study/' + study_id, retries=retries)
        _source_result('api', details is not None)
        if details is not None:
            return details
    use, retries = _source_usable('workcopy')
    if use:
        try:
            details = study_details_from_workcopy(study_id)
            _source_result('workcopy', True)
            if details is not None:
                print(f'{study_id}: study record from the API working copy of the ISA-Tab files')
                return details
        except SourceUnavailable as e:
            _source_result('workcopy', False)
            print(f'{study_id}: working copy not available ({e})')
    print(f'{study_id}: study record from the ISA-Tab files on the EBI FTP server')
    try:
        return study_details_from_isa(study_id)
    except SourceUnavailable as e:
        raise RuntimeError(f'Could not retrieve study details for {study_id}: {e}')


def rename_duplicated_column_names(df):
    # We need this to get rid of duplicated columns
    new_column_names = []

    # Dictionary to keep track of the count of duplicate column names
    dup_count = {}

    # Iterate through the columns using index
    for idx, col in enumerate(df.columns):
        if col in df.columns[:idx]:
            # If the column name is a duplicate, append a count to its name
            dup_count[col] = dup_count.get(col, 0) + 1
            new_name = f"{col}_{dup_count[col]}"
            new_column_names.append(new_name)
        else:
            # If not a duplicate, keep the original name
            new_column_names.append(col)

    # Set the new column names to the DataFrame
    df.columns = new_column_names

    return df


def _table_headers(table, prefix):
    """Column headers for a MetaboLights assay/sample table, in the order of its data rows.

    Field 'index' values are not always contiguous (MetaboLights leaves gaps, e.g. where
    columns without a header were dropped), while the rows hold one value per field in
    index order. So headers are assigned by rank of index, not by the index itself.
    """
    fields = sorted(table['fields'].values(), key=lambda field: field['index'])
    widths = {len(row) for row in table['data']}
    if widths and widths != {len(fields)}:
        if widths == {fields[-1]['index'] + 1}:
            # rows padded to the highest index: position == index
            headers = [None] * widths.pop()
            for field in fields:
                headers[field['index']] = prefix + field['header']
            return headers
        raise ValueError(f"Table rows have {sorted(widths)} values but {len(fields)} fields are declared")
    return [prefix + field['header'] for field in fields]


_ACCESSION_COLUMN = re.compile(r'Term Accession Number(?:[._]\d+)*$')
_ANNOTATION_COLUMN = re.compile(r'(?:Unit|Term Source REF|Term Accession Number)(?:[._]\d+)*$')
_MISSING_VALUES = ('', 'nan', 'none', 'missing value', 'not applicable', 'not available', 'n/a', 'na')


def _is_missing(series):
    return series.isna() | series.astype(str).str.strip().str.lower().isin(_MISSING_VALUES)


def _accession_column(df, value_column):
    """The 'Term Accession Number' column belonging to value_column. In ISA tables it
    follows the value directly, after an optional Unit / Term Source REF column."""
    if value_column not in df.columns:
        return None
    columns = list(df.columns)
    position = columns.index(value_column)
    prefix = value_column.split('_', 1)[0] + '_'
    for column in columns[position + 1:position + 4]:
        if not column.startswith(prefix) or not _ANNOTATION_COLUMN.match(column[len(prefix):]):
            return None
        if _ACCESSION_COLUMN.match(column[len(prefix):]):
            return column
    return None


def _accession_ids(series, ontology):
    """Ontology ids ('<ONTOLOGY>:<id>') from accession values such as
    'http://purl.obolibrary.org/obo/NCBITaxon_9606', '.../ontology/NCBITAXON/9606' (BioPortal)
    or 'UBERON:0001088'; '' if none."""
    pattern = r'(?i)(?:^|[/#:\s])' + ontology + r'[_:/]([A-Za-z0-9]+)\s*$'
    ids = series.astype(str).str.strip().str.extract(pattern, expand=False)
    return ids.fillna('').map(lambda x: f'{ontology}:{x}' if x else '')


_ncbi_by_id = {}


def _ncbi_term_by_id(allowedTerm_dict):
    """taxid -> allowed 'taxid|name' term, built once (the allowed list has ~3M entries)."""
    if not _ncbi_by_id:
        for term in allowedTerm_dict['NCBITaxonomy']['allowed_values']:
            taxid, _, _ = str(term).partition('|')
            _ncbi_by_id.setdefault(taxid, term)
    return _ncbi_by_id


_ncbi_by_name = {}


def _ncbi_term_by_name(allowedTerm_dict):
    """lower-case name -> first allowed 'taxid|name' term (as get_taxonomy_id_from_name__allowedTerms), built once."""
    if not _ncbi_by_name:
        for term in allowedTerm_dict['NCBITaxonomy']['allowed_values']:
            taxid, sep, name = str(term).partition('|')
            if sep and '|' not in name:
                _ncbi_by_name.setdefault(name.lower(), term)
    return _ncbi_by_name


def _investigation_submission_year(study_id):
    """Submission year from the study's public ISA investigation file, used when the
    API has no submissionYear. Returns None if it cannot be read."""
    url = f"https://ftp.ebi.ac.uk/pub/databases/metabolights/studies/public/{study_id}/i_Investigation.txt"
    try:
        response = EBI_SESSION.get(url, timeout=(10, 120))
        if response.status_code != 200:
            return None
    except requests.RequestException:
        return None
    for line in response.text.splitlines():
        for key in ('Study Submission Date', 'Study Public Release Date'):
            if line.startswith(key):
                match = re.search(r'(19|20)\d{2}', line)
                if match:
                    return int(match.group(0))
    return None



# ---------------------------------------------------------------------------------------------
# Study-level disease -> DOIDCommonName / HealthStatus of the samples
#
# MetaboLights records the disease of a study at study level (design descriptors such as
# 'EFO:ulcerative colitis', factor names such as 'Diabetes', disease names in the title) and the
# groups at sample level ('Disease = control', 'Cohort = Healthy Control / Mild / Severe').
# Cases get the study disease (or the narrower DOID term their label names) and
# 'unhealthy (NOS)'; individuals labelled healthy (or plain controls of a human case-control
# factor) get 'healthy'. Health is a property of the individual, so all samples of a subject
# with a case sample are patient samples. Unlabelled samples, other groups, treatment /
# genotype / exposure controls and tissue-only 'normal' labels never become healthy.
# ---------------------------------------------------------------------------------------------

# DOID terms too general to stand for the disease of a study: disease, its top-level classes,
# disease by infectious agent, syndrome, physical disorder, genetic disease
_DOID_GENERIC = {'DOID:4', 'DOID:7', 'DOID:14566', 'DOID:150', 'DOID:0014667', 'DOID:0050117', 'DOID:630',
                 'DOID:0080015', 'DOID:225'}
_DOID_NEOPLASM_ROOT = 'DOID:14566'     # disease of cellular proliferation
_DOID_INFECTION_ROOT = 'DOID:0050117'  # disease by infectious agent
# NCBI divisions whose samples can carry a (human) disease ontology term
_DISEASE_DIVISIONS = {'Primates', 'Rodents', 'Mammals'}


def _dlabel(value):
    return re.sub(r'\s+', ' ', str(value).replace('_', ' ')).strip().lower()


def _doid_abbreviation(synonym):
    """DOID synonyms such as 'MDS', 'DSS', 'FACE', 'CLN3' are ambiguous outside the ontology."""
    s = str(synonym).strip()
    return len(s) <= 4 or (s == s.upper() and len(s) <= 12)


class _DoidIndex:
    """Exact DOID label / synonym lookup, is_a ancestors and the MetaboLights disease aliases."""

    def __init__(self, doid_table, allowed_labels, parents, aliases):
        self.allowed = set(allowed_labels)
        self.by_label = {str(v).strip().lower(): str(v) for v in self.allowed if not str(v).lower().startswith('obsolete')}
        self.id_of, self.label_of = {}, {}
        syn = {}
        if doid_table is not None and len(doid_table):
            for i, l, s in zip(doid_table['DOIDOntologyIndex'].astype(str), doid_table['Label'].astype(str),
                               doid_table['Synonym'].astype(str)):
                if l in self.allowed and not l.lower().startswith('obsolete'):
                    self.id_of.setdefault(l, i)
                    self.label_of.setdefault(i, l)
                    if not _doid_abbreviation(s):
                        syn.setdefault(s.strip().lower(), set()).add(l)
        self.by_synonym = {k: next(iter(v)) for k, v in syn.items() if len(v) == 1}
        self.parents = parents or {}
        # alias -> DOID label; an empty label blocks a DOID synonym ('insulin resistance')
        self.aliases = {a: l for a, l in (aliases or {}).items() if not l or l in self.allowed}
        self._anc = {}

    def term(self, value):
        k = _dlabel(value)
        if len(k) <= 2 or k in _MISSING_VALUES:
            return None
        if k in self.aliases:
            return self.aliases[k] or None
        return self.by_label.get(k) or self.by_synonym.get(k)

    def ancestors(self, label):
        if label not in self._anc:
            seen, stack = set(), list(self.parents.get(self.id_of.get(label), []))
            while stack:
                p = stack.pop()
                if p not in seen:
                    seen.add(p)
                    stack.extend(self.parents.get(p, []))
            self._anc[label] = {self.label_of[p] for p in seen if p in self.label_of}
        return self._anc[label]

    def under(self, label, root_id):
        return self.id_of.get(label) == root_id or self.label_of.get(root_id) in self.ancestors(label)

    def generic(self, label):
        return self.id_of.get(label) in _DOID_GENERIC

    def foreign(self, label):
        """A name of another vocabulary (a disease for the HPO index): never a study term."""
        return False

    def related(self, a, b):
        """The more specific of a and b if one equals or descends from the other, else None."""
        if a == b or b in self.ancestors(a):
            return a
        if a in self.ancestors(b):
            return b
        return None

    def single(self, terms):
        """The most specific term if all terms lie on one is_a chain, else None."""
        terms = [t for t in terms if t]
        best = terms[0] if terms else None
        for t in terms[1:]:
            best = self.related(best, t)
            if best is None:
                return None
        return best


_mammal_taxids_cache = {}


def _mammal_taxids(NCBIRankDivision_table):
    """NCBI taxon ids (str) of the divisions whose samples can carry a disease term."""
    key = id(NCBIRankDivision_table)
    if key not in _mammal_taxids_cache:
        _mammal_taxids_cache.clear()
        t = NCBIRankDivision_table
        _mammal_taxids_cache[key] = set(t.loc[t['NCBIDivision'].isin(_DISEASE_DIVISIONS), 'TaxonID'].astype(str))
    return _mammal_taxids_cache[key]


def _doid_synonyms(idx):
    """lower-case synonym -> DOID label: exact synonyms naming one term, without abbreviations,
    with the MetaboLights aliases (an alias without a term removes the synonym)."""
    out = dict(idx.by_synonym)
    for alias, label in idx.aliases.items():
        if label:
            out[alias] = label
        else:
            out.pop(alias, None)
    return out


_doid_index_cache = {}


def _doid_index(doid_table, allowedTerm_dict, parents, aliases):
    key = (id(doid_table), id(parents), id(aliases))
    if key not in _doid_index_cache:
        _doid_index_cache.clear()
        _doid_index_cache[key] = _DoidIndex(doid_table, allowedTerm_dict['DOIDCommonName']['allowed_values'], parents, aliases)
    return _doid_index_cache[key]


# ---------------------------------------------------------------------------------------------
# Abbreviations a study defines for itself
#
# Ontology synonyms that are abbreviations ('AKI', 'LPS', 'DSS') are ambiguous across studies
# and never matched on their own. Within one study they stand for the full name if the study
# text (title, description, design descriptors, protocols, factor names) defines them:
# (a) '<full name> (<ABBR>)' or '<ABBR> (<full name>)', the full name being a term of the
# index and starting with the first letter of the abbreviation, or (b) the abbreviation is an
# ontology synonym of a term whose label or synonym the text names.
# ---------------------------------------------------------------------------------------------

def _abbreviation_like(token):
    t = str(token)
    return (2 <= len(t) <= 8 and re.fullmatch(r'[A-Za-z][A-Za-z0-9\-]*', t) is not None
            and (sum(c.isupper() for c in t) >= 2 or (any(c.isupper() for c in t) and any(c.isdigit() for c in t))))


def _abbreviation_synonyms(table, label_col, synonym_col, allowed):
    """lower-case abbreviation synonym -> (label, names of the term) for abbreviations of one term."""
    names, abbr = {}, {}
    for l, s in zip(table[label_col].astype(str), table[synonym_col].astype(str)):
        if l not in allowed:
            continue
        names.setdefault(l, {_dlabel(l)})
        if _abbreviation_like(s.strip()):
            abbr.setdefault(s.strip().lower(), set()).add(l)
        elif not _doid_abbreviation(s):
            names[l].add(_dlabel(s))
    return {a: (next(iter(ls)), sorted(names[next(iter(ls))], key=len, reverse=True)) for a, ls in abbr.items() if len(ls) == 1}


def _study_text(content, protocols):
    parts = [content.get('title'), content.get('description') or content.get('studyDescription')]
    parts += [re.sub(r'^[A-Za-z]+:', '', str(d.get('description') or '')) for d in content.get('descriptors') or []]
    parts += [f.get('name') for f in content.get('factors') or []]
    parts += [f'{n}: {d}' for n, d in protocols]
    return re.sub(r'<[^>]+>', ' ', '. '.join(str(p) for p in parts if p))


def _value_tokens(df, columns):
    """Abbreviation-like tokens of the sample labels (lower case)."""
    out = set()
    for c in columns:
        for v in df[c].dropna().astype(str).unique():
            out.update(w.lower() for w in re.split(r'[^A-Za-z0-9\-]+', v) if _abbreviation_like(w))
            out.update(w.lower() for w in re.split(r'[^A-Za-z0-9]+', v) if _abbreviation_like(w))
    return out


def _defined_abbreviations(text, resolve, abbr_syn, wanted):
    """lower-case abbreviation -> lower-case full name, for the abbreviations in wanted that the
    study text defines; resolve(full name) is the term it names (or None)."""
    found = {}
    text = str(text or '')
    if not wanted:
        return found
    for m in re.finditer(r"([A-Za-z][\w'\-, ]{2,120}?)\s*\(\s*([A-Za-z][A-Za-z0-9\-]{1,7}?)s?\s*\)", text):
        ab = m.group(2)
        if ab.lower() not in wanted or not _abbreviation_like(ab):
            continue
        words = m.group(1).split()
        for L in range(min(8, len(words)), 0, -1):
            full = _dlabel(' '.join(words[-L:]).strip(' ,'))
            if full[:1] == ab[0].lower() and resolve(full):
                found.setdefault(ab.lower(), set()).add(full)
                break
    for m in re.finditer(r"(?<![\w-])([A-Za-z][A-Za-z0-9\-]{1,7})\s*\(\s*([^()]{3,120}?)\s*\)", text):
        ab, full = m.group(1), _dlabel(m.group(2))
        if ab.lower() in wanted and _abbreviation_like(ab) and full[:1] == ab[0].lower() and resolve(full):
            found.setdefault(ab.lower(), set()).add(full)
    low = _dlabel(text)
    for ab in wanted:
        if ab in found or ab not in abbr_syn:
            continue
        label, names = abbr_syn[ab]
        if any(re.search(r'(?<![a-z0-9])' + re.escape(n) + r'(?![a-z0-9])', low) for n in names if len(n) >= 4):
            found[ab] = {_dlabel(label)}
    # an abbreviation the text defines in two ways is left out
    return {ab: next(iter(f)) for ab, f in found.items() if len(f) == 1 and resolve(next(iter(f)))}


def _expand_abbreviations(label, abbr):
    """The lower-case label with the study's abbreviations written out ('lps 1 mg' -> 'lipopolysaccharide 1 mg')."""
    for ab, full in (abbr or {}).items():
        if ab in label:
            label = re.sub(r'(?<![a-z0-9])' + re.escape(ab) + r'(?![a-z0-9])', full, label)
    return label


class _StudyTerms:
    """An index that also reads the abbreviations one study defines."""

    def __init__(self, idx, abbr):
        self._idx, self._abbr = idx, abbr

    def __getattr__(self, name):
        return getattr(self._idx, name)

    def term(self, value):
        t = self._idx.term(value)
        if t is None and self._abbr:
            k = _dlabel(value)
            if k in self._abbr:
                t = self._idx.term(self._abbr[k])
            else:
                k2 = _d_strip_number(k)
                if k2 in self._abbr:
                    t = self._idx.term(self._abbr[k2])
        return t


# ---------------------------------------------------------------------------------------------
# Study-level health condition -> HPOCommonName / HealthStatus of the samples
#
# A condition of the individual that is not a disease (sepsis, cachexia, insulin resistance,
# acute kidney injury): labels of the HPO 'Phenotypic abnormality' branch. The same study /
# sample-group / subject rules as for the study disease apply, with an HPO index in place of
# the DOID one. Diseases stay in DOIDCommonName: HPO terms named like a DOID term, neoplasms
# and the generic 'Abnormality of ...' classes are left out, and a DOID name in a label makes
# that group 'another disease'.
# ---------------------------------------------------------------------------------------------

_HPO_NOT_CONDITION = re.compile(r'^(?:abnormal|phenotypic abnormality)|neoplas|cancer|tumou?r|leuka?emia|metasta|carcinoid|'
                                r'(?:carcin|lymph|sarc|gli|blast|melan|myel|aden|cyt|endym|ingi|papill|fibr|lip|terat|neur|'
                                r'schwann|chord|hamart|thym|insulin|semin|theli|angi)omas?\b', re.I)
_FOREIGN_TERM = '\x00disease'


class _HpoIndex(_DoidIndex):
    """HPO terms that are not diseases; names of DOID terms are 'foreign'."""

    def __init__(self, hpo_table, allowed_labels, aliases, doid_idx):
        self.doid = doid_idx
        allowed = {str(l) for l in allowed_labels if not _HPO_NOT_CONDITION.search(str(l))
                   and not (doid_idx is not None and self._names_disease(doid_idx, l))}
        t = hpo_table.rename(columns={'HPOOntologyIndex': 'DOIDOntologyIndex'})
        # names of the left-out terms ('Acute colitis', 'Astrocytoma') mark another condition
        out = t[~t['Label'].astype(str).isin(allowed)]
        self.excluded = {_dlabel(x) for x in pd.concat([out['Label'], out['Synonym']]).astype(str)}
        t = t[t['Label'].astype(str).isin(allowed)]
        if doid_idx is not None:
            # a synonym that names a disease is left to DOIDCommonName
            t = t[[not doid_idx.term(x) for x in t['Synonym'].astype(str)]]
        super().__init__(t, allowed, None, aliases)
        self.abbr_syn = _abbreviation_synonyms(t, 'Label', 'Synonym', allowed)

    @staticmethod
    def _names_disease(doid_idx, label):
        """The label is a DOID term, or one with modifiers ('Acute colitis', 'Cerebral ischemia')."""
        words = _dlabel(label).split()
        return any(doid_idx.term(' '.join(words[i:])) for i in range(len(words)) if i == 0 or len(' '.join(words[i:])) >= 5)

    def term(self, value):
        t = super().term(value)
        if t is None and _dlabel(value) not in self.aliases and (
                _dlabel(value) in self.excluded or (self.doid is not None and self.doid.term(value))):
            return _FOREIGN_TERM
        return t

    def foreign(self, label):
        return label == _FOREIGN_TERM


_hpo_index_cache = {}


def _hpo_index(hpo_table, allowedTerm_dict, aliases, doid_idx):
    key = (id(hpo_table), id(aliases), id(doid_idx))
    if key not in _hpo_index_cache:
        _hpo_index_cache.clear()
        _hpo_index_cache[key] = _HpoIndex(hpo_table, allowedTerm_dict['HPOCommonName']['allowed_values'], aliases, doid_idx)
    return _hpo_index_cache[key]


# ---------------------------------------------------------------------------------------------
# Exposure (ECTO; PECO for plants) of the exposed groups
#
# Only a sample label that names the exposure gets a term: a treatment / exposure / stress
# column value that is an ECTO/PECO label or exact synonym, or that is one with the ontology's
# own wording ('cisplatin' -> 'exposure to cisplatin', plants: 'drought' -> 'drought
# exposure'), or a 'yes' in a column named after the exposure ('Cold stress = yes'). Controls,
# vehicles, sham, untreated, ambient and recovery groups name no exposure and get nothing.
# Animals (and humans) get ECTO terms, plants PECO terms (ECTO where PECO has none). Cultured
# cells, microbes and diets (own column) are left out.
# ---------------------------------------------------------------------------------------------

_E_COL = re.compile(r'treat|expos|drug|compound|chemical|toxi|pollut|pesticide|herbicide|fungicide|insecticide|stress|stimul|'
                    r'challenge|infect|inocul|administ|interven|therap|medicat|irradiat|radiation|\blight\b|temperature|'
                    r'irrigation|watering|salin|ozone|antibiotic|inhibitor|agent|condition|group|stimulus', re.I)
# columns about something else: diet (own column), extraction, media, analytes, the organism itself
_E_SKIP_COL = re.compile(r'diet|feed|food|nutri|supplement|solvent|extract|standard|medium|media|substrate|carbon source|'
                         r'source|analyte|organism|strain|genotype|species|genus|cell|tissue|time|duration|dose|concentration|'
                         r'amount|quantity|replicate|batch|sex|gender|age\b|cultivar|variety|location|site|disease|diagnos|'
                         r'phenotype|health|blood group|family|maternal|mother|father|donor|history|term source|term accession|'
                         r'^unit|age group|ethnic', re.I)
# vehicles, controls, unexposed and recovering groups
_E_NEGATIVE = re.compile(r'\b(?:control|controls|ctrl|vehicle|sham|mock|placebo|untreated|unexposed|unstressed|non-?treated|'
                         r'non-?exposed|without|no|not|none|normal|ambient|baseline|before|pre|recover(?:y|ed)?|rewater(?:ing|ed)?|'
                         r'well-?watered|wild type|wt|blank|qc)\b|^-|^non-|^un', re.I)
# ECTO/PECO terms too general to stand for an exposure, or vehicles and materials of the study
# rather than exposures
_E_GENERIC = {'exposure to water', 'exposure to solvent', 'exposure to dimethyl sulfoxide', 'dimethyl sulfoxide exposure',
              'exposure to drug', 'exposure to chemical', 'exposure to chemical entity', 'exposure to molecular entity',
              'exposure to lifestyle', 'exposure to lipid', 'exposure to peptide', 'exposure to amino acid',
              'exposure to polysaccharide', 'exposure to blood', 'exposure to feces', 'exposure to bile', 'exposure to soil',
              'exposure to rock', 'exposure to sand', 'exposure to tension', 'exposure to stress', 'exposure to environment',
              'exposure to treatment', 'exposure to material entity', 'exposure to organism', 'chemical exposure',
              'exposure to physical object', 'exposure to excreta', 'exposure to acids', 'exposure to acid',
              'exposure to carbohydrate', 'carbohydrate exposure', 'water environment exposure', 'watering exposure', 'temperature exposure'}
_E_BLOOD_CELLS = re.compile(r'blood|mononuclear|pbmc|erythrocyte|platelet|leukocyte|neutrophil|lymphocyte|monocyte|sperm|oocyte', re.I)
_E_CULTURE_COL = re.compile(r'cultivation|culture|incubation|passage', re.I)
_E_CULTURE_VALUE = re.compile(r'\b(?:medium|media|culture|cultured|in vitro)\b', re.I)
_E_DIET = re.compile(r'\bdiet\b|via ingestion|\bfood\b|\bmeal\b|feeding', re.I)
_E_POSITIVE = re.compile(r'(?:yes|y|true|positive|\+|exposed|treated|stressed|stress|treatment|with stress(?: treatment)?|with treatment)')
# wording around the exposure in a sample label ('drought stress', 'treated with cisplatin')
_E_AFFIX = re.compile(r'^(?:treated with|treatment with|exposed to|exposure to|\+)\s*|\s*(?:stress(?:ed)?|treatment|treated|exposure|exposed|'
                      r'regimen|group|administration|challenge|condition|conditions)$')


class _ExposureIndex:
    """Exact ECTO/PECO label / synonym lookup with the MetaboLights condition aliases."""

    def __init__(self, exposure_table, allowed_labels, aliases):
        allowed = {str(l) for l in allowed_labels}
        self.allowed = allowed
        self.ont = {}
        syn = {}
        for i, l, s in zip(exposure_table['ExposureIndex'].astype(str), exposure_table['Label'].astype(str),
                           exposure_table['Synonym'].astype(str)):
            if l not in allowed or l.lower() in _E_GENERIC or _E_DIET.search(l) or l.lower().startswith('obsolete'):
                continue
            self.ont.setdefault(l, i.split(':')[0])
            if not _doid_abbreviation(s):
                syn.setdefault(_dlabel(s), set()).add(l)
        self.by_label = {_dlabel(l): l for l in self.ont}
        self.by_synonym = {k: next(iter(v)) for k, v in syn.items() if len(v) == 1}
        self.aliases = {a: l for a, l in (aliases or {}).items() if not l or l in self.ont}
        self.abbr_syn = _abbreviation_synonyms(exposure_table, 'Label', 'Synonym', set(self.ont))

    def _exact(self, key, plant):
        """(term or None, whether the key decides): PECO terms only for plants."""
        t = self.aliases[key] if key in self.aliases else (self.by_label.get(key) or self.by_synonym.get(key))
        if not t:
            return None, key in self.aliases
        return (t if self.ont[t] == 'ECTO' or plant else None), True

    def term(self, value, plant, abbr=None):
        """Exposure term a sample label names (PECO only for plants), else None. abbr: the
        abbreviations the study defines (lower case -> full name), read as the full name."""
        k = _expand_abbreviations(_dlabel(value), abbr)
        if len(k) <= 3 or k in _MISSING_VALUES or _E_NEGATIVE.search(k):
            return None
        keys = [k]
        stripped = _E_AFFIX.sub('', k).strip()
        if stripped != k and len(stripped) > 3:
            keys.append(stripped)
        for key in keys:
            t, decided = self._exact(key, plant)
            if decided:
                return t
        # the ontology's own wording: PECO 'X exposure' for plants, ECTO 'exposure to X'
        for key in keys:
            for cand in ([key + ' exposure'] if plant else []) + ['exposure to ' + key]:
                t, decided = self._exact(cand, plant)
                if decided:
                    return t
        return None


_exposure_index_cache = {}


def _exposure_index(exposure_table, allowedTerm_dict, aliases):
    key = (id(exposure_table), id(aliases))
    if key not in _exposure_index_cache:
        _exposure_index_cache.clear()
        _exposure_index_cache[key] = _ExposureIndex(exposure_table, allowedTerm_dict['Exposure']['allowed_values'], aliases)
    return _exposure_index_cache[key]


def _sample_exposures(eidx, df, colkind, eligible, plant, abbr=None):
    """Exposure term per row from the treatment / exposure columns (None where unset).

    df: rows with the 'Samples_*' columns; colkind: column -> 'Characteristics' / 'Factor Value';
    eligible: samples of an individual animal or plant; plant: plant rows. A row whose
    columns name different exposures gets none."""
    found = []
    # tissue or cells kept in culture ('cultivation time' recorded): the treatment is in vitro
    for c in df.columns:
        if _E_CULTURE_COL.search(re.sub(r'_\d+$', '', str(c)[len('Samples_'):])) and colkind.get(c) in ('Characteristics', 'Factor Value'):
            eligible = eligible & _is_missing(df[c])
    for c, kind in colkind.items():
        if c not in df.columns or kind not in ('Characteristics', 'Factor Value'):
            continue
        name = re.sub(r'_\d+$', '', str(c)[len('Samples_'):])
        if _E_SKIP_COL.search(name) or not _E_COL.search(name):
            continue
        # a treatment given in a culture medium ('Blank medium' among the groups)
        if df.loc[eligible, c].astype(str).str.contains(_E_CULTURE_VALUE).any():
            continue
        # a column named after an exposure ('Cold stress', 'Drought stress'): 'yes' names it
        col_term = {p: eidx.term(name, p, abbr) for p in (False, True)}
        cache = {}

        def one(v, p):
            if (v, p) not in cache:
                lab = _dlabel(v)
                t = eidx.term(v, p, abbr)
                if t is None and col_term[p] and _E_POSITIVE.fullmatch(lab):
                    t = col_term[p]
                cache[(v, p)] = t
            return cache[(v, p)]
        vals = pd.Series([one(v, bool(p)) if e and isinstance(v, str) else None
                          for v, p, e in zip(df[c], plant, eligible)], index=df.index, dtype=object)
        if vals.notna().any():
            found.append(vals)
    out = pd.Series(None, index=df.index, dtype=object)
    for vals in found:
        clash = out.notna() & vals.notna() & (out != vals)
        out = out.where(out.notna(), vals)
        out[clash] = '\x00clash'
    return out.where(out != '\x00clash')


# ---- sample labels (matched on the lower-case label, '_' as space) ----
_D_PEOPLE = r'(?:subjects?|individuals?|volunteers?|donors?|participants?|people|persons?|adults?|children|women|men|group|cohort|controls?|population|samples?)'
_D_HEALTHY = re.compile(r'(?:healthy|health|healthy controls?|health controls?|hc|hcs|hv|hvs|non-?diseased?|disease-? ?free|'
                        r'healthy ' + _D_PEOPLE + r'(?: ' + _D_PEOPLE + r')?|normal healthy ' + _D_PEOPLE + r'|'
                        r'(?:control|controls) \(healthy\)|healthy \((?:control|controls|hc)\)|'
                        r'healthy (?:volunteers?|controls?|subjects?|donors?) \((?:hcs?|hvs?|hs|hd|ctrl|controls?)\)|'
                        r'normal (?:subjects?|volunteers?|individuals?|donors?|participants?|people|persons?|adults?|children)|'
                        r'control \(normal\))')
_D_CONTROL = re.compile(r'(?:controls?|control group|control subjects?|control individuals?|control participants?|control cohort|'
                        r'ctrl|ctr|cont|cntrl|con|ctl|non-?cases?|nc|control \(c\)|controls \(c\)|c)')
_D_NORMAL = re.compile(r'(?:normal|normals|normal group)')
# words that make a 'control' label a treatment / genotype / exposure / tissue / QC control
_D_NOT_PLAIN_CONTROL = re.compile(r'patient|case|treat|vehicle|sham|diet|wild|positive|negative|model|drug|knock|mock|placebo|dmso|pbs|'
                                  r'saline|infect|sirna|shrna|vector|empty|scrambl|\bnc\b|medium|media|chow|day|week|hour|time|dose|'
                                  r'cell|blank|qc|quality|pool|standard|reference|gene|mutant|wt|\bko\b|siblings?|littermates?|'
                                  r'tumou?r|lesion|tissue|exposure|stress|baseline|before|after|pre|post|\+|\bwith\b|disease|'
                                  r'unwell|sick|symptom|other|non-?responder')
_D_GENERIC_CASE = re.compile(r'(?:patients?|cases?|diseased|affected|disease|sick|patient group|case group|cases group|'
                             r'disease group|patients group|patient cohort|case subjects?|lesions?|lesional(?: skin| tissues?)?|'
                             r'lesion tissues?)')
_D_SEVERITY = re.compile(r'(?:mild|moderate|severe|critical|critically ill|very severe|mild-moderate|moderate-severe|'
                         r'stage (?:[0-4]|i{1,3}|iv)[abc]?|early stage|late stage|advanced(?: stage)?|early|advanced)'
                         r'(?: (?:disease|cases?|patients?|group|cohort|stage))?')
_D_POSITIVE = re.compile(r'(?:positive|pos|\+|infected|yes|present|affected|y|true)')
_D_MODEL = re.compile(r'(?:model|models|model group|disease model|mouse model|rat model|model mice|model rats|model control)')
_D_TUMOUR = re.compile(r'(?:tumou?rs?|tumou?r tissues?|tumou?r samples?|cancer|cancers|cancer tissues?|cancerous(?: tissues?)?|'
                       r'carcinoma tissues?|malignant(?: tissues?)?|neoplastic tissues?|primary tumou?rs?|tumou?r core)')
# normal tissue of a patient ('adjacent normal', 'para-carcinoma', 'non-lesional skin')
_D_TISSUE_NORMAL = re.compile(r'adjacent|para-?(?:carcinoma|cancer|tumou?r)|paracancer|peri-?tumou?r|non-?tumou?r|non-?cancer|'
                              r'normal (?:tissues?|mucosa)|tumou?r margin|\bmargins?\b|non-?lesion|peri-?lesion|'
                              r'unaffected (?:skin|tissues?|side)|distant normal|histologically normal')
_D_NEGATED = re.compile(r'\b(?:without|non|no|not|free|negative|control|controls|healthy|normal|pre|post|risk|history|relatives?|'
                        r'family|suspected|possible|probable|excluded|recovered|cured|remission)\b|-free|non-|\bun|'
                        r'vaccin|immuni[sz]|inactivat|de-?activ|attenuat|killed|mock|sham')
_D_FILLER = re.compile(r'\b(?:patients?|pts|cases?|subjects?|individuals?|participants?|group|cohort|positive|diagnosed|'
                       r'samples?|people|donors?|mice|rats|model|with|newly|confirmed|the)\b')

# ---- sample columns ----
_D_STATUS_COL = re.compile(r'disease|diagnos|condition|status|state|phenotype|health|case|patient|clinical|cohort|group|'
                           r'class|category|population|severity|stage|staging|subject type|participant type|sample group|'
                           r'type of subject|outcome', re.I)
_D_TREATMENT_COL = re.compile(r'treat|drug|dose|compound|diet|feed|food|genotype|genetype|mutant|mutation|knock|gene|strain|'
                              r'variant|exposure|infection|inocul|challenge|stimul|vaccin|medium|media|time|temperature|'
                              r'cell line|cell type|cultivar|replicate|batch|sex|gender|age|bmi|weight|location|site|species', re.I)
_D_TISSUE_COL = re.compile(r'(?:organism part|tissue|tissue type|sample type|biospecimen|sample material|material|specimen|'
                           r'tissue location|tissue origin|organ|sample origin|sample source|source tissue|tumou?r type|'
                           r'tissue status|lesion)', re.I)
# columns about someone other than the sampled individual, or not about health at all
_D_SKIP_COL = re.compile(r'^(?:organism|sample name|source name)$|family|maternal|mother|father|parent|donor (?:disease|status)|'
                         r'history|comorbid|medication|term source|term accession|^unit', re.I)

# cultured cells, organoids, media: the subject's disease does not apply; milk is often
# described by the infant's condition
_D_CULTURE = re.compile(r'\bcells?\b|cell lines?|cultur|organoid|spheroid|ipsc|pluripotent|immortal|xenograft|in vitro|'
                        r'supernatant|\bmedium\b|\bmedia\b|lysate|\bclones?\b|transfect', re.I)
_D_PRIMARY_CELLS = re.compile(r'blood|mononuclear|pbmc|erythrocyte|platelet|leukocyte|neutrophil|lymphocyte|sperm|semen|oocyte|'
                              r'bone marrow|exfoliat|urine|stool|feces|faeces|islet|sputum|lavage|swab', re.I)
_D_MILK = re.compile(r'\bmilk\b|colostrum', re.I)
_D_MICROBIOTA_RECIPIENT = re.compile(r'microbiota (?:from|of) (?:human|patient|donor)|(?:human|patient|donor)[- ]derived (?:fecal |faecal |gut )?microbiota|'
                                     r'fecal microbiota transplant|faecal microbiota transplant|\bFMT\b|humani[sz]ed (?:gut|microbiota|mice)|'
                                     r'colonized with (?:human|patient)|gnotobiotic', re.I)
_D_CELL_COLS = re.compile(r'Samples_(?:[^_]*cell lines?[^_]*|cell line name|cell type|cell)(?:_\d+)?', re.I)
_D_SPECIMEN_COLS = re.compile(r'Samples_(?:organism part|sample type|biospecimen|tissue|tissue type|variant|sample material|'
                              r'material|specimen|sample source|sample origin)(?:_\d+)?', re.I)


def _disease_not_subject(df):
    """Rows that are cultured cells (cell lines, organoids, culture media) or milk."""
    out = pd.Series(False, index=df.index)
    for c in df.columns:
        if _D_CELL_COLS.fullmatch(str(c)):
            out |= ~_is_missing(df[c])
        elif _D_SPECIMEN_COLS.fullmatch(str(c)):
            v = df[c].astype(str)
            out |= (v.str.contains(_D_CULTURE) & ~v.str.contains(_D_PRIMARY_CELLS)) | v.str.contains(_D_MILK)
    return out


def _d_strip_number(v):
    """Replicate numbering: 'control-1', 'con 3', 'HC_12' -> 'control', 'con', 'hc'."""
    w = re.sub(r'[\s\-#.]*\d{1,3}$', '', v).strip()
    return w if re.search(r'[a-z]', w) else v


def _study_disease_terms(idx, descriptors, factor_names, title):
    """DOID terms of the study (label -> source): design descriptors and factor names that
    name a DOID term exactly; disease names in the title only when neither has one."""
    out = {}
    for d in descriptors:
        s = str(d or '').strip()
        m = re.match(r'^[A-Za-z]+:(.*)$', s)
        t = idx.term(m.group(1) if m else s)
        if t and not idx.generic(t) and not idx.foreign(t):
            out.setdefault(t, 'descriptor')
    for f in factor_names:
        t = idx.term(f)
        if t and not idx.generic(t) and not idx.foreign(t):
            out.setdefault(t, 'factor')
    if not out and title:
        words = re.findall(r"[\w'\-]+", str(title).lower())
        i = 0
        while i < len(words):
            hit = None
            for L in range(min(8, len(words) - i), 0, -1):
                t = idx.term(' '.join(words[i:i + L]))
                if t and not idx.generic(t) and len(' '.join(words[i:i + L])) >= 4:
                    hit = (t, L)
                    break
            if hit:
                # a name of another vocabulary uses up its words ('irritable bowel syndrome')
                if not idx.foreign(hit[0]):
                    out.setdefault(hit[0], 'title')
                i += hit[1]
            else:
                i += 1
    return out


def _text_abbreviations(idx, terms, text):
    """Abbreviations the study text defines for its diseases: 'type 2 diabetes mellitus (T2DM)'."""
    out = {}
    for m in re.finditer(r"([A-Za-z][\w'\- ,]{2,80}?)\s*\(\s*([A-Za-z][A-Za-z0-9\-]{1,9}?)s?\s*\)", str(text or '')):
        ab = m.group(2)
        # an acronym: two capitals at least ('T2DM', 'AAU'), not '(Cases)' or '(n=20)'
        if sum(ch.isupper() for ch in ab) < 2:
            continue
        words = m.group(1).lower().split()
        for L in range(min(8, len(words)), 0, -1):
            t = idx.term(' '.join(words[-L:]))
            if t and not idx.generic(t) and any(idx.related(t, d) for d in terms):
                # its first letter starts the disease name ('HSP' for 'Henoch-Schonlein purpura')
                if words[-L][0] == ab[0].lower():
                    out.setdefault(ab.lower(), t)
                break
    return out


def _disease_phrases(idx, D, abbr):
    """(lower-case name, term): labels, exact synonyms and aliases of the study terms, text abbreviations."""
    out = set()
    for d in D:
        out.add((d.lower(), d))
        out.update((k, d) for k, l in idx.by_synonym.items() if l == d and len(k) >= 5)
        out.update((k, d) for k, l in idx.aliases.items() if l == d)
    out.update((k, t) for k, t in abbr.items() if len(k) >= 2)
    return sorted(out, key=lambda x: -len(x[0]))


def _disease_label_kind(idx, value, coltype, D, abbr, phrases, col_term, human):
    """Kind of one sample label and, for 'case', the DOID term it names.

    'healthy' (explicitly healthy), 'control' / 'normal' (plain control), 'case' (names a study
    disease or a narrower term), 'case_broad' (names a broader term), 'case_generic'
    (patient / case / severity word), 'case_tumour', 'case_model', 'patient' (normal tissue
    of a patient), 'other_disease' (another disease), 'other'."""
    v0 = _dlabel(value)
    if v0 in _MISSING_VALUES:
        return None, None
    v = _d_strip_number(v0)
    if _D_HEALTHY.fullmatch(v) or _D_HEALTHY.fullmatch(v0):
        return 'healthy', None
    if _D_CONTROL.fullmatch(v):
        return 'control', None
    if _D_NORMAL.fullmatch(v):
        return ('patient', None) if coltype == 'tissue' else ('normal', None)
    if re.search(r'\bcontrols?\b', v) and not _D_NOT_PLAIN_CONTROL.search(v):
        return 'control', None      # 'non-IBS (control)', 'Study Control'
    if _D_TISSUE_NORMAL.search(v) and (coltype == 'tissue' or re.search(r'tissue|mucosa|adjacent|para|peri|margin|lesion|skin', v)):
        return 'patient', None
    if _D_TUMOUR.fullmatch(v):
        return ('case_tumour', None) if any(idx.under(d, _DOID_NEOPLASM_ROOT) for d in D) else ('other', None)
    if col_term and _D_POSITIVE.fullmatch(v):
        return 'case', col_term     # 'Diabetes = yes'
    t = idx.term(v) or idx.term(v0)
    if t is None:
        stripped = re.sub(r'\s+', ' ', _D_FILLER.sub(' ', re.sub(r'\([^)]*\)', ' ', v))).strip(' -,')
        if stripped and stripped != v:
            t = idx.term(stripped) or abbr.get(stripped)
    if t is None:
        t = abbr.get(v) or abbr.get(re.sub(r'\s*\([^)]*\)', '', v).strip())
    if t is not None:
        if idx.foreign(t):
            return 'other_disease', None
        if idx.generic(t):
            return 'case_generic', None
        if any(idx.related(t, d) == t for d in D):
            return 'case', t
        if any(idx.related(t, d) == d for d in D):
            return 'case_broad', t
        return 'other_disease', t
    # the label contains a study disease ("Parkinson's Disease on Medication", 'IBS - diarrhea')
    if not _D_NEGATED.search(v):
        # short abbreviations only at the start ('CRC_Rectum', 'IBS - diarrhea'), not 'SLE without LN'
        hits = {term for name, term in phrases
                if re.search(('^' if len(name) <= 3 else r'(?<![\w-])') + re.escape(name) + r'(?![\w-])', v)}
        if hits:
            return 'case', idx.single(list(hits))
    if _D_GENERIC_CASE.fullmatch(v):
        return 'case_generic', None
    if re.search(r'\b(?:patients?|cases?)\b', v) and not re.search(r'\d|cohort|\bvs\b|versus', v) and not _D_NEGATED.search(v):
        return 'case_generic', None     # 'Patient (+) SARS-CoV-2', 'Rhinitis patient'
    if coltype in ('status', 'disease_factor') and _D_SEVERITY.fullmatch(v):
        return 'case_generic', None
    if v in ('infected', 'infection', 'infected patients') and any(idx.under(d, _DOID_INFECTION_ROOT) for d in D):
        return 'case_generic', None
    if not human and _D_MODEL.fullmatch(v):
        return 'case_model', None
    return 'other', None


_D_CASE_KINDS = ('case', 'case_generic', 'case_tumour', 'case_model')


def _study_disease_labels(idx, df, colkind, descriptors, factor_names, title, text, eligible, human, subject=None,
                          plain_control_healthy=True, extra_abbr=None):
    """DOIDCommonName and HealthStatus per row from the study disease and the sample groups.

    df: rows with the 'Samples_*' columns; colkind: column -> 'Characteristics' / 'Factor Value';
    eligible: rows that are samples of an individual; subject: subject identifier per row or
    None; extra_abbr: further abbreviation -> term the study defines (the index reads them too);
    plain_control_healthy: plain controls of a human case-control factor are healthy
    (otherwise only explicitly healthy labels). Returns (doid, health, info), doid and health
    being Series with None where unset."""
    info = {}
    none = pd.Series(None, index=df.index, dtype=object)
    D = _study_disease_terms(idx, descriptors, factor_names, title)
    # several unrelated study diseases ('HCC in a cirrhosis cohort'): a plain 'control' may be a patient
    several = len(D) > 1 and idx.single(list(D)) is None
    cols = []
    for c, kind in colkind.items():
        if c not in df.columns or kind not in ('Characteristics', 'Factor Value'):
            continue
        name = re.sub(r'_\d+$', '', str(c)[len('Samples_'):])
        if _D_SKIP_COL.search(name):
            continue
        nt = idx.term(name)
        nt = None if nt and idx.foreign(nt) else nt
        if _D_TISSUE_COL.fullmatch(name):
            ct = 'tissue'
        elif nt and not idx.generic(nt):
            ct = 'disease_factor'
        elif _D_STATUS_COL.search(name) and not _D_TREATMENT_COL.search(name):
            ct = 'status'
        elif _D_TREATMENT_COL.search(name):
            ct = 'treatment'
        else:
            ct = 'other'
        cols.append((c, name, ct))
    # diseases named in disease / status columns are study diseases too
    for c, name, ct in cols:
        if ct in ('status', 'disease_factor') or re.search(r'disease|diagnos|condition|phenotype', name, re.I):
            for v in df.loc[eligible, c].dropna().unique():
                t = idx.term(_d_strip_number(_dlabel(v))) or idx.term(v)
                if t and not idx.generic(t) and not idx.foreign(t):
                    D.setdefault(t, 'value:' + name)
    info['D'] = D
    if not D or not eligible.any():
        return none, none.copy(), info
    abbr = _text_abbreviations(idx, list(D), f'{title or ""}. {text or ""}')
    for k, t in (extra_abbr or {}).items():
        abbr.setdefault(k, t)
    phrases = _disease_phrases(idx, D, abbr)
    ctype = {c: ct for c, _, ct in cols}
    lab, has = {}, {}
    for c, name, ct in cols:
        col_term = idx.term(name) if ct == 'disease_factor' else None
        lab[c] = {v: _disease_label_kind(idx, v, ct, D, abbr, phrases, col_term, human)
                  for v in df.loc[eligible, c].dropna().astype(str).unique()}
        # a broader term ('liver cancer' in an HCC study) stands for the study disease, unless
        # another label of the column names that disease or a narrower one (a different group)
        named = {t for k, t in lab[c].values() if k == 'case'}
        for v, (k, t) in list(lab[c].items()):
            if k == 'case_broad':
                narrower = [d for d in D if idx.related(t, d) == d]
                if any(idx.related(n, d) == n for n in named for d in narrower):
                    lab[c][v] = ('other_disease', t)
                else:
                    lab[c][v] = ('case', idx.single(narrower) or t)
        # in animal studies, one disease label for all samples outside a disease column is the
        # study topic ('Sample type = non-alcoholic fatty liver disease' for every diet group)
        if not human and ct not in ('status', 'disease_factor') and len(lab[c]) == 1:
            v, (k, t) = next(iter(lab[c].items()))
            if k in _D_CASE_KINDS or k == 'case_broad':
                lab[c][v] = ('other', None)
        has[c] = {k for k, _ in lab[c].values() if k}
    group_cols = [c for c, _, ct in cols if has[c] & set(_D_CASE_KINDS + ('healthy', 'control', 'normal', 'patient'))
                  or (has[c] and ct in ('status', 'disease_factor'))]
    info['group_cols'] = group_cols
    if not group_cols:
        return none, none.copy(), info
    generic_term = idx.single(list(D))
    neo_term = idx.single([d for d in D if idx.under(d, _DOID_NEOPLASM_ROOT)])
    # a plain control is a healthy individual only in a human case-control factor with named cases
    cc_cols = {c for c in group_cols if ctype[c] in ('status', 'disease_factor') and not several
               and ('case' in has[c] or ('case_generic' in has[c] and generic_term))}
    # a bare 'normal' only in a disease / health column of a study without a tumour disease
    normal_ok = {c for c in cc_cols if re.search(r'disease|diagnos|health|condition', c, re.I)
                 and not any(idx.under(d, _DOID_NEOPLASM_ROOT) for d in D)}

    def decide(values):
        kinds = []
        for c, x in zip(group_cols, values):
            if isinstance(x, str):
                k, t = lab[c].get(x, (None, None))
                if k:
                    kinds.append((k, t, c))
        # disease / status columns (and tissue labels) decide; treatment and other columns only without them
        primary = [x for x in kinds if ctype[x[2]] in ('status', 'disease_factor') or (ctype[x[2]] == 'tissue' and x[0] != 'other')]
        kinds = primary or kinds
        ks = {k for k, _, _ in kinds}
        cases = [(k, t) for k, t, _ in kinds if k in _D_CASE_KINDS]
        st_ctrl = [(k, c) for k, _, c in kinds if k in ('control', 'normal') and ctype[c] in ('status', 'disease_factor')]
        st_other = [c for k, _, c in kinds if k in ('other', 'other_disease') and ctype[c] in ('status', 'disease_factor')]
        if cases or 'patient' in ks:
            if 'healthy' in ks or st_ctrl:
                return 'conflict', None
            if not cases and st_other:
                return 'other', None
            ts = [t if k == 'case' else neo_term if k == 'case_tumour' else generic_term for k, t in cases]
            if 'patient' in ks:
                ts.append(neo_term or generic_term)
            return ('case' if cases else 'patient'), idx.single([x for x in ts if x])
        if 'other_disease' in ks:
            return 'other', None
        if 'healthy' in ks:
            return 'healthy', None
        if st_ctrl:
            ok = plain_control_healthy and human and any(c in cc_cols and (k == 'control' or c in normal_ok) for k, c in st_ctrl)
            return ('healthy' if ok else 'control'), None
        return ('other' if kinds else None), None

    keys = pd.Series(list(zip(*[df[c].map(lambda x: None if pd.isna(x) else str(x)) for c in group_cols])), index=df.index)
    decided = {k: decide(k) for k in keys[eligible].unique()}
    # rows whose labels name the term itself (not only 'patients', a severity or 'model')
    def named(values):
        return any(lab[c].get(x, (None, None))[0] == 'case' for c, x in zip(group_cols, values) if isinstance(x, str))
    _named = {k: named(k) for k in decided}
    info['named'] = keys.map(lambda k: _named.get(k, False)) & eligible
    status = keys.map(lambda k: decided.get(k, (None, None))[0]).where(eligible)
    term = keys.map(lambda k: decided.get(k, (None, None))[1]).where(eligible)
    # health belongs to the individual: every sample of a patient is a patient sample, and a
    # subject labelled both control and case keeps neither label on its control samples
    if subject is not None:
        sub = subject.where(eligible & ~_is_missing(subject))
        grp = pd.DataFrame({'s': sub, 'st': status, 't': term})[sub.notna()]
        per = grp.groupby('s')['st'].agg(lambda x: set(x.dropna()))
        conflicted = set(per[per.map(lambda x: bool(x & {'case', 'patient'}) and bool(x & {'healthy', 'control', 'conflict'}))].index)
        info['subjects_conflicting'] = len(conflicted)
        demote = grp.index[grp['s'].isin(conflicted) & grp['st'].isin(['healthy', 'control'])]
        status[demote] = 'conflict'
        # subject ids of animals are often reused per group ('1'..'6'), so only human subjects
        if human and not conflicted:
            sick = grp[grp['st'].isin(['case', 'patient'])]
            subject_term = sick.groupby('s')['t'].agg(lambda x: idx.single(list(x)) if x.notna().all() else None)
            add = grp.index[grp['s'].isin(subject_term.index) & ~grp['st'].isin(['case', 'patient'])]
            status[add] = 'case_by_subject'
            term[add] = grp.loc[add, 's'].map(subject_term)
    info['status'] = status.value_counts().to_dict()
    sick = status.isin(['case', 'patient', 'case_by_subject'])
    doid = pd.Series(None, index=df.index, dtype=object)
    doid[sick & term.notna()] = term[sick & term.notna()]
    health = pd.Series(None, index=df.index, dtype=object)
    health[sick] = 'unhealthy (NOS)'
    health[status == 'healthy'] = 'healthy'
    return doid, health, info


# ---------------------------------------------------------------------------
# TimePointHoursRelativeToEvent, FastingStatus, SampleGroupFreeText
# ---------------------------------------------------------------------------
# The rules are those of the Workbench converter (MWB_to_REDU.py), applied to the
# Characteristics / Factor Value columns of the ISA sample table. ISA-Tab keeps the unit
# of a value in the 'Unit' column right after it ('24' + 'hour').
# time-like sample columns that are no sampling time relative to an event: light / dark
# periods and other growth conditions, symptom counts, plant / embryo ages (days after
# sowing, flowering, fertilisation), birth months, urine collection periods, ...
_ML_TP_KEY_EXCLUDE_RE = re.compile(
    r'\b(?:light|dark|photoperiod|dawn|dusk|temperature|humidity|salt|quantity|watering|container|medium|'
    r'sonication|partition|active|average|mean|total|fever|diarrh\w*|allergy|birth|sowing|sown|planting|'
    r'germination|anthesis|pollination|flowering|fertili[sz]ation|urine|sputter|coating|hospitali[sz]\w*|disease|illness|'
    r'smok\w*|month of)\b'
    # keys that hold a time themselves ('4 day post-exposure sampling', 'Salt_Quantity_after_08d'):
    # their values are no time ('0' / '1', concentrations)
    r'|\d\s*-?\s*(?:s|sec|min|mins|minutes?|h|hr|hrs|hours?|d|days?|w|wk|wks|weeks?|mo|months?|y|yrs?|years?)\b')
# keys that name the reference event ('Time post intervention', 'Timepoint postpartum')
_ML_TP_EVENT_KEY_RE = re.compile(r'\bpost\w*|\b(?:after|since|from|following|intervention|treatment|treated|therapy|'
                                 r'surgery|transplant\w*|infection|inoculation|injury|exposure|dosing|supplement\w*|'
                                 r'diet|baseline|follow)\b')
# an ISA unit that is one time unit ('hour', 'Days', 'hour post inoculation')
_ML_TIME_UNIT_RE = re.compile(r'(' + '|'.join(sorted(_TP_UNIT_WORDS, key=len, reverse=True)) + r')'
                              r'(?:\s+(?:post|after|following)\b.*)?')
_ML_BARE_NUMBER_RE = re.compile(r'[-+]?\d+(?:\.\d+)?')
_ML_GROUP_NA_VALUES = set(_GROUP_NA_VALUES) | {'not applicable', 'not available', 'missing value'}
# subject identifiers under other names than the Workbench ones ('Identification Number', 'Patient')
_ML_GROUP_ID_KEY_RE = re.compile(r'identification|identifier|^(?:patients?|animals?|individuals?|donors?|mice|mouse|rats?|'
                                 r'pigs?|participants?|persons?|volunteers?|subject (?:no|number|code))$')
# technical sample attributes that are no groups of the study design (identifiers, analysis
# and preparation settings)
_ML_GROUP_TECH_KEY_RE = re.compile(r'\b(?:barcodes?|urls?|uri|accessions?|operators?|dilutions?|aliquots?|volumes?|'
                                   r'evaporation|dda|dia|analysis|analyses|sample roles?|storage|freezer|tubes?|box)\b')


def _ml_unit_column(df, column):
    """The ISA 'Unit' column directly after a value column, else None."""
    i = df.columns.get_loc(column)
    if isinstance(i, (int, np.integer)) and i + 1 < len(df.columns) and str(df.columns[i + 1]).startswith('Samples_Unit'):
        return df.columns[i + 1]
    return None


def _ml_time_unit(unit):
    """'hour' for the ISA unit 'hour' / 'Hours' / 'hour post inoculation', '' for an
    empty unit, None for any other unit."""
    u = ' '.join(str(unit).strip().lower().split())
    if u in _MISSING_VALUES or u == '-':
        return ''
    m = _ML_TIME_UNIT_RE.fullmatch(u)
    return m.group(1) if m else None


def _ml_key(column):
    """Factor / characteristic name of a sample column ('Samples_Time point_1' -> 'Time point')."""
    return re.sub(r'_\d+$', '', str(column)[len('Samples_'):] if str(column).startswith('Samples_') else str(column))


def _ml_long_values(df, columns, unit_in_key=False):
    """The sample columns as a long table (row position as 'filename', Key, Value), as the
    Workbench functions read them. A bare number gets the time unit of its 'Unit' column
    ('24' + 'hour' -> '24 hour', or the key 'Fasting (hour)' with unit_in_key); a number
    with another unit ('10' + 'millimolar') is no time."""
    parts = []
    pos = np.arange(len(df))
    for c in columns:
        v = df[c].astype(str).str.strip().where(~_is_missing(df[c]), '')
        key = pd.Series(_ml_key(c), index=df.index)
        uc = _ml_unit_column(df, c)
        if uc is not None:
            unit = df[uc].map(_ml_time_unit)
            bare = v.str.fullmatch(_ML_BARE_NUMBER_RE)
            if unit_in_key:
                key = key.where(~(bare & unit.fillna('').ne('')), key + ' (' + unit.fillna('') + ')')
            else:
                v = v.where(~(bare & unit.fillna('').ne('')), v + ' ' + unit.fillna(''))
            v = v.where(~(bare & unit.isna()), '')
        parts.append(pd.DataFrame({'filename': pos, 'Key': key.values, 'Value': v.values}))
    if not parts:
        return pd.DataFrame(columns=['filename', 'Key', 'Value'])
    long = pd.concat(parts, ignore_index=True)
    return long[long['Value'] != ''].reset_index(drop=True)


def ml_sample_timepoints(df, columns, study_text=''):
    """Per-row TimePointHoursRelativeToEvent (Series on df's index) from the time-point
    columns (Workbench rules: explicit number and time unit, in the value or its 'Unit'
    column; ages, storage / incubation times, clock times, dates, ranges excluded)."""
    columns = [c for c in columns if timepoint_key_tier(_ml_key(c)) is not None
               and not _ML_TP_KEY_EXCLUDE_RE.search(_norm_key(_ml_key(c)))]
    # months / years without a named event and without a baseline 0 are mostly ages
    # ('Collection time' 6 / 24 / 36 / 48 months in a longitudinal animal study)
    keep = []
    for c in columns:
        part = _ml_long_values(df, [c])
        key = _ml_key(c)
        long_unit = _timepoint_key_unit(key) in ('mo', 'y') or \
            part['Value'].str.lower().str.contains(r'\b(?:mo|months?|y|yrs?|years?)\b', regex=True).any()
        hours = {timepoint_hours(v, key) for v in part['Value']} - {None}
        if long_unit and 0 not in hours and not _ML_TP_EVENT_KEY_RE.search(_norm_key(key)):
            print(f"[TimePointHoursRelativeToEvent] column {key!r} not used: months / years without an event or baseline")
            continue
        keep.append(c)
    columns = keep
    out = pd.Series(np.nan, index=df.index, dtype=object)
    long = _ml_long_values(df, columns)
    if len(long) == 0:
        return out
    hours = sample_timepoint_hours(long, study_text)
    per_row = hours.groupby(long['filename']).first().dropna()
    out.iloc[per_row.index.to_numpy()] = per_row.values
    return out


def ml_sample_fasting(df, columns):
    """Per-row FastingStatus (Series on df's index) from fasting / feeding columns and
    fasting-state values of design columns (Workbench rules)."""
    out = pd.Series(np.nan, index=df.index, dtype=object)
    long = _ml_long_values(df, columns, unit_in_key=True)
    if len(long) == 0:
        return out
    # 'Non-fasting' / 'not fasting' as the Workbench words 'non-fasted' / 'not fasted'
    long['Value'] = long['Value'].str.replace(r'(?i)\bnon[- ]?fasting\b', 'non-fasted', regex=True) \
        .str.replace(r'(?i)\bnot fasting\b', 'not fasted', regex=True)
    status = sample_fasting_status(long)
    per_row = status.groupby(long['filename']).first().dropna()
    out.iloc[per_row.index.to_numpy()] = per_row.values
    return out


# ---- study-level fasting statement in the protocols ----
# protocols that do not describe the subjects
_FAST_PROTOCOL_SKIP_RE = re.compile(r'chromatograph|mass spec|data transformation|metabolite identification|nmr|'
                                    r'flow injection|histolog|magnetic resonance|assay|normali[sz]ation', re.I)
# fasting words (a bare 'fast' only in 'an overnight fast', 'a 12-h fast', 'to fast')
_FAST_WORD_RE = re.compile(r'\b(?:non-?|un)?fast(?:ed|ing)\b|\bto fast\b|\bfast overnight\b|'
                           r'\b(?:an?|overnight|h|hr|hrs|hours?|\d+)[- ]fast(?![\w-])(?!\s*,)(?!\s+(?:and|gas|gc|lc|liquid|'
                           r'polarity|switching|filtration|scan|untargeted|lipid|screening|robust|sampling|isolation|data|'
                           r'fourier|spin|heating|way|mode|growing|growth|prep|method|protocol|freezing|cooling)\b)')
# fasting values as measurement names ('fasting blood glucose', 'impaired fasting glucose')
_FAST_MEASURE_RE = re.compile(r'(?:impaired )?fasting (?:blood |plasma |serum |venous )?(?:glucose|insulin|sugar|'
                              r'blood sugar|plasma glucose|lipids?|lipid profiles?|triglycerides?|tg|cholesterol|hdl|'
                              r'ldl|c-?peptide|levels?|values?)\b(?: \([a-z]+\))?|\b(?:fbg|fpg|fbs)\b')
# the statement is not about one fasting condition of all sampled subjects
_FAST_UNCLEAR_RE = re.compile(
    r'intermittent|ramadan|re-?fe(?:d|eding)|\b(?:fed|fasted|fasting|fast)\s+or\b|\bor\s+(?:fed|fasted|fasting|non)|'
    r'\b(?:fed|fasted|fasting)\s+and\s+(?:fed|fasted|fasting|non)|\bvs\b|versus|\beither\b|\bsome\b|subset|\bhalf\b|'
    r'majority|\bmost\b|\bgroups?\b|\bexcept|\bunless\b|\bif\b|possible|not specified|not required|restriction|'
    r'regardless|optional|respective|durations|different|various|varied|fasting experiment|challeng|\bdose|administ|'
    r'received|inject|infus|gavage|weigh|larva|starv|pre-?fasting|post-?fasting|before (?:the )?fasting|baseline|'
    r'infect|inocul|tolerance|\bgtt|\bogtt|\bipgtt|\bitt\b|clamp|meal|breakfast|lunch|postprandial|post-prandial|'
    r'ingest|\bdrink|caffein|\bwithout water\b')
_FAST_H = r'(?:h|hr|hrs|hours?)'
_FAST_NUM = r'\d+(?:\.\d+)?(?:\s*(?:-|–|to)\s*\d+(?:\.\d+)?)?'
_FAST_DUR = r'(?:at least |minimum (?:of )?|more than |over |≥\s*)?' + _FAST_NUM + r'\s*-?\s*' + _FAST_H
_FASTED_STATEMENT_RE = re.compile('|'.join([
    # 'after an overnight fast', 'following a 12-hour overnight fast', 'after a 20-h fasting period'
    r'\b(?:after|following|under|with)\s+(?:an?\s+|the\s+)?(?:(?:overnight|' + _FAST_DUR + r')\s*-?\s*)+'
    r'(?:\(\s*' + _FAST_DUR + r'\s*\)\s*)?(?:fast|fasting)\b',
    # 'after 12 hours of fasting', 'after fasting for at least 8 h', 'after overnight fasting (8 h)'
    r'\b(?:after|following|under|with)\s+' + _FAST_DUR + r'\s+of\s+fasting\b',
    r'\b(?:after|following|under)\s+(?:a\s+)?fasting\s+(?:period\s+of\s+)?(?:for\s+)?(?:overnight|' + _FAST_DUR + r')',
    r'\b(?:after|following|under|underwent)\s+(?:overnight\s+)?fasting\s+overnight',
    # 'mice were fasted overnight', 'fasted for 12 h', 'birds were fasted 12 h'
    r'(?<![-\w])fasted\s+(?:overnight|(?:for\s+)?' + _FAST_DUR + r'|for\s+(?:at least\s+)?\d+\s+days?)',
    # 'overnight fasted animals', '12 h-fasted', 'from fasted individuals', 'from fasting participants'
    r'\b(?:overnight|\d+\s*-?\s*' + _FAST_H + r')[- ]?fasted\b',
    r'\b(?:from|in|of)\s+(?:the\s+)?(?:fasted|fasting)\s+(?:individuals|patients|subjects|participants|animals|mice|'
    r'rats|donors|volunteers|pigs|dogs|cats)\b',
    # 'in the fasting state', 'under fasting conditions', 'while they were fasting'
    r'\b(?:in|at|under|during)\s+(?:a\s+|the\s+)?(?:early morning\s+)?(?:fasting|fasted)\s+(?:state|states|condition|'
    r'conditions|status)\b',
    r'\bwhile\s+(?:they\s+were\s+|subjects\s+were\s+|participants\s+were\s+)?fasting\b',
    # 'fasting blood samples were collected', 'fasting peripheral venous blood was drawn'
    r'(?<![-\w])fasting\s+(?:\(\s*' + _FAST_DUR + r'\s*\)\s+)?(?:[a-z0-9]+\s+){0,2}?(?:blood|plasma|serum|sera|urine)\s*'
    r'(?:samples?|specimens?)?\s+(?:were|was)\s+(?:\w+\s+)?(?:collected|drawn|obtained|taken|withdrawn)',
    # 'instructed to fast overnight', 'required to fast for more than 8 h', 'to ensure overnight fasting'
    r'\bto\s+fast\s+(?:overnight|for\s+' + _FAST_DUR + r')',
    r'\bto\s+ensure\s+(?:at least\s+\d+\s*' + _FAST_H + r'\s+of\s+|overnight\s+)fasting\b',
    r'\b(?:subjects|participants|patients|volunteers|individuals)\s+fasted\s+overnight\b',
]))
_NONFASTED_STATEMENT_RE = re.compile('|'.join([
    r'\bnon-?fasting\s+(?:\w+\s+)?(?:blood\s+|plasma\s+|serum\s+)?samples?\b',
    r'\b(?:were|was)\s+not\s+fasting\b',
    r'\bwithout\s+(?:prior\s+)?fasting\b',
    r'\bnon-?fasted\s+(?:male\s+|female\s+)?(?:\S+\s+)?(?:mice|rats|animals|subjects|participants|patients)\b',
    r'\b(?:during|under|in)\s+(?:a\s+|the\s+)?non-?fasting\s+(?:conditions?|state)\b',
    r'\bfrom\s+(?:\S+\s+){0,3}non-?fasting\s+(?:\S+\s+){0,2}(?:volunteers|participants|subjects|donors|individuals)\b',
]))
# specimens a statement may be limited to
_FAST_SPECIMENS = {'blood': re.compile(r'\b(?:blood|plasma|serum|sera|venipuncture|phlebotomy|venous)\b'),
                   'urine': re.compile(r'\burin(?:e|ary)\b'),
                   'feces': re.compile(r'\b(?:feces|faeces|fecal|faecal|stool|stools)\b'),
                   'saliva': re.compile(r'\bsaliva\b')}


def _fasting_hours_in(text):
    """Hours of the fasting durations in a statement (minimum of a range); [] if none."""
    hours = [float(a) for a in re.findall(r'(\d+(?:\.\d+)?)\s*(?:(?:-|–|to)\s*\d+(?:\.\d+)?\s*)?-?\s*' + _FAST_H + r'\b', text)]
    hours += [24 * float(a) for a in re.findall(r'(\d+(?:\.\d+)?)\s*days?\b', text)]
    words = {'six': 6, 'eight': 8, 'ten': 10, 'twelve': 12, 'sixteen': 16, 'four': 4, 'five': 5}
    hours += [words[w] for w in re.findall(r'\b(' + '|'.join(words) + r')[- ]?' + _FAST_H + r'\b', text)]
    return hours


def protocol_fasting_status(protocols):
    """One FastingStatus for all subjects of a study from its protocol texts, with the
    specimens it is limited to (None: all) and the statements, or None.

    Every sentence that speaks about the subjects' fasting must be an explicit statement
    of one condition ('after an overnight fast', 'fasted for 12 h', 'in the fasting state',
    'non-fasting blood samples'); a sentence with a fasting word that is no such statement,
    or that names groups, conditions, doses, tests or meals, leaves the study open. A
    fasting time shorter than 8 h is not 'fasted'. Fasting values as measurement names
    ('fasting blood glucose') are ignored."""
    found, specimens, statements = set(), set(), []
    all_specimens = False
    for name, description in protocols:
        if _FAST_PROTOCOL_SKIP_RE.search(str(name)):
            continue
        text = html.unescape(re.sub(r'<[^>]+>', ' ', str(description or ''))).replace('\xa0', ' ')
        # (no split after 'approx. 55%' or 'e.g. (...)')
        for sentence in re.split(r'(?<=[.!?])\s+(?=[A-Z])', text):
            s = ' '.join(sentence.lower().split())
            s2 = _FAST_MEASURE_RE.sub(' ', s)
            if not _FAST_WORD_RE.search(s2):
                continue
            if _FAST_UNCLEAR_RE.search(s2):
                return None
            fasted = list(_FASTED_STATEMENT_RE.finditer(s2))
            nonfasted = _NONFASTED_STATEMENT_RE.search(s2)
            if bool(fasted) == bool(nonfasted):
                return None  # no explicit statement, or both
            if fasted:
                # durations given with the statement (in it or right after it)
                hours = [h for m in fasted for h in _fasting_hours_in(s2[m.start():m.end() + 30])]
                if hours and min(hours) < 8:
                    return None
                found.add('fasted')
            else:
                found.add('non-fasted')
            # specimens named in the statement itself ('urine and fasting serum samples':
            # serum), else in the sentence
            named = {k for k, rx in _FAST_SPECIMENS.items() for m in fasted if rx.search(m.group(0))} \
                or {k for k, rx in _FAST_SPECIMENS.items() if rx.search(s2)}
            if named:
                specimens |= named
            else:
                all_specimens = True
            statements.append(' '.join(sentence.split()))
    if len(found) != 1:
        return None
    return found.pop(), (None if all_specimens else specimens), statements


def _specimen_family(body_part):
    b = str(body_part).lower()
    if re.search(r'blood|plasma|serum', b):
        return 'blood'
    if 'urine' in b:
        return 'urine'
    if re.search(r'feces|faeces|stool', b):
        return 'feces'
    if 'saliva' in b:
        return 'saliva'
    return None


# ---- SampleGroupFreeText ----
def ml_sample_group_texts(df, factor_columns, factor_order=(), sample_column='Samples_Sample Name', study_id=''):
    """Per-row sample group in the internal form of finalize_sample_groups (sample ID,
    then 'factor -> value' parts) from the Factor Value columns, in the study's factor
    order and wording (values with their unit). Factors that have a ReDU column of their
    own or are no groups (Workbench rules), and factors with a different value for nearly
    every sample, are left out. ';' inside a name or value becomes ','."""
    out = pd.Series(np.nan, index=df.index, dtype=object)
    if not factor_columns or sample_column not in df.columns:
        return out
    order = {str(n).strip().lower(): i for i, n in enumerate(factor_order)}
    columns = sorted(factor_columns, key=lambda c: (order.get(_ml_key(c).strip().lower(), len(order)),
                                                    list(df.columns).index(c)))
    clean = lambda x: ' '.join(str(x).replace(';', ',').replace(_GROUP_SID, ' ').replace(_GROUP_KV, ' ')
                               .replace(_GROUP_PART, ' ').split())
    sid = df[sample_column].astype(str).str.strip()
    n = sid.nunique()
    values, drop = {}, []
    for c in columns:
        key = clean(_ml_key(c))
        nk = _norm_key(key)
        v = df[c].map(clean)
        uc = _ml_unit_column(df, c)
        if uc is not None:
            unit = df[uc].map(clean).where(~_is_missing(df[uc]), '')
            v = v.where(unit == '', v + ' ' + unit)
        v = v.where(~v.str.lower().isin(_ML_GROUP_NA_VALUES) & ~_is_missing(df[c]))
        if not nk or _GROUP_KEY_EXCLUDE_RE.search(nk) or _GROUP_ID_KEY_RE.fullmatch(nk) or _ML_GROUP_ID_KEY_RE.search(nk) \
                or _ML_GROUP_TECH_KEY_RE.search(nk):
            drop.append(key)
            continue
        # the recorded subject identifier under another name
        if 'SubjectIdentifierAsRecorded' in df.columns:
            subject = df['SubjectIdentifierAsRecorded'].map(clean).where(~_is_missing(df['SubjectIdentifierAsRecorded']))
            both = v.notna() & subject.notna()
            if both.sum() >= 3 and (v[both] == subject[both]).mean() >= 0.9:
                drop.append(key)
                continue
        per_sample = pd.DataFrame({'sid': sid, 'v': v}).dropna().drop_duplicates()
        if n >= 8 and per_sample['v'].nunique() >= 0.9 * n:
            drop.append(key)
            continue
        values[key] = v
    if drop:
        print(f"[SampleGroupFreeText] {study_id}: factors not used (ReDU column or no groups): {drop}")
    if not values:
        print(f"[SampleGroupFreeText] {study_id}: no grouping factors left")
        return out
    text = pd.Series('', index=df.index, dtype=object)
    for key, v in values.items():
        part = (key + _GROUP_KV + v).where(v.notna(), '')
        text = text.where(part == '', text.where(text == '', text + _GROUP_PART) + part)
    # one sample listed with different factors gets none
    ok = pd.DataFrame({'sid': sid, 't': text}).groupby('sid')['t'].transform('nunique') == 1
    text = text.where(ok & (text != ''))
    return (sid.map(clean) + _GROUP_SID + text).where(text.notna())


def Metabolights2REDU(study_id, **kwargs):
    """
    Converts Metabolights study data to a REDU table format.

    Args:
    study_id: The ID of the Metabolights study to convert.

    Returns:
    A DataFrame in the REDU table format with processed and aligned data from the Metabolights study,
    or an empty DataFrame if no applicable data is found.
    """
    # the per-study API, or the study's ISA-Tab files when the API is down
    study_details = get_study_details(study_id)
    if study_details is None:
        raise RuntimeError(f"Could not retrieve study details for {study_id}")
    if study_details.get('content') is None:
        raise StudyNotAvailable(study_details.get('message', ''))

    study_assays = study_details['content']['assays']
    # protocol name -> description, for the study-level extraction/internal standard fields
    protocols = [(str(p.get('name') or ''), str(p.get('description') or ''))
                 for p in (study_details['content'].get('protocols') or [])]

    submissionYear = (study_details['content'].get('derivedData') or {}).get('submissionYear')
    if submissionYear in (None, '', 0):
        submissionYear = _investigation_submission_year(study_id)

    df_assays = pd.DataFrame()
    ms_study_assays = []

    #get study assays if they are MS. There can be multiple assay tables in the same study
    for index, assay_json in enumerate(study_assays):
        if assay_json['technology'] == 'mass spectrometry' or assay_json['technology'] == 'mass spectrometry assay':        

            # Extract headers in the correct order
            headers = _table_headers(assay_json['assayTable'], 'Assay_')

            df = pd.DataFrame(assay_json['assayTable']['data'])
            df = rename_duplicated_column_names(df)

            df.columns = headers

            print(f'Assay {index} has {len(df)} rows and {len(df.columns)} columns.')
            
            ms_study_assays.append(df)

    if len(ms_study_assays) > 0:
        # Union of the columns in first-seen order: the position of a column matters, as
        # an ontology accession sits in the 'Term Accession Number' column after its value.
        # Repeated headers are renamed before the union so none of them is dropped.
        ms_study_assays = [rename_duplicated_column_names(df) for df in ms_study_assays]
        all_columns = list(dict.fromkeys(col for df in ms_study_assays for col in df.columns))

        aligned_dfs = []
        for df in ms_study_assays:
            # Add missing columns with NaN values
            for col in all_columns:
                if col not in df.columns:
                    df[col] = np.nan
            # Reorder columns and add to the aligned list
            aligned_dfs.append(df[all_columns])

        # Concatenate all dataframes into a single dataframe
        df_assays = pd.concat(aligned_dfs, ignore_index=True)

        print(f'Assay table has {len(df_assays)} rows and {len(df_assays.columns)} columns.')

        #get info on samples (there is only one sample table per study)
        headers = _table_headers(study_details['content']['sampleTable'], 'Samples_')

        df_samples = pd.DataFrame(study_details['content']['sampleTable']['data'])
        df_samples = rename_duplicated_column_names(df_samples)
        df_samples.columns = headers
        # 'Characteristics' / 'Factor Value' per sample column (the headers drop the kind)
        _sample_kind = {'Samples_' + f['header']: f.get('fieldType')
                        for f in study_details['content']['sampleTable']['fields'].values()}
        factor_columns = [c for c, k in _sample_kind.items() if k == 'Factor Value']
        # the kind of every sample column by position (a characteristic and a factor can
        # have the same header, e.g. 'Sample type' and 'Sample type_1')
        _sample_table = study_details['content']['sampleTable']
        _sample_kind_pos = dict(zip(df_samples.columns, _table_headers(
            {'fields': {k: dict(f, header=str(f.get('fieldType'))) for k, f in _sample_table['fields'].items()},
             'data': _sample_table['data']}, '')))
    
        # Sample names are typed by hand in both sheets; ignore surrounding whitespace
        df_assays['Assay_Sample Name'] = df_assays['Assay_Sample Name'].astype(str).str.strip()
        df_samples['Samples_Sample Name'] = df_samples['Samples_Sample Name'].astype(str).str.strip()
        n_unmatched = (~df_assays['Assay_Sample Name'].isin(df_samples['Samples_Sample Name'])).sum()
        if n_unmatched > 0:
            print(f'{n_unmatched} of {len(df_assays)} assay rows have no matching sample in the sample table and are dropped.')

        df_study = df_assays.merge(df_samples, left_on='Assay_Sample Name', right_on='Samples_Sample Name', how='inner')
        df_study['row_id'] = range(1, len(df_study) + 1)

        print(f'Sample table has {len(df_study)} rows and {len(df_study.columns)} columns.')

        # Duplicate rows if we have mzml AND raw files. A study can have several raw/derived file
        # columns (e.g. 'Raw Spectral Data File' holding .wiff.scan and 'Raw Spectral Data File.1'
        # holding the .wiff), so every such column contributes rows; derived (mzML) files come first
        file_columns = [col for prefix in ['Assay_Derived Spectral Data File', 'Assay_Raw Spectral Data File']
                        for col in df_study.columns if re.fullmatch(re.escape(prefix) + r'(\.\d+)?', col)]
        file_dfs = []
        for col in file_columns:
            df_files = df_study[df_study[col].str.contains(r'\.', regex=True, na=False)].copy()
            df_files['filepath'] = df_files[col]
            df_files = df_files.drop(columns=file_columns)
            print(f'Assay table has {len(df_files)} rows and {len(df_files.columns)} columns.')
            if len(df_files) > 0:
                file_dfs.append(df_files)

        if len(file_dfs) == 0:
            return pd.DataFrame()
        df_study = pd.concat(file_dfs, ignore_index=True).copy(deep=True)
        
        df_study = rename_duplicated_column_names(df_study)
        df_study['filename'] = df_study['filepath']#.apply(lambda x: os.path.basename(x))

        # Ensure the filename is a string and lowercase
        df_study['filename_lower'] = df_study['filename'].astype(str).str.lower()

        # Keep supported MS files, uncompressed or compressed (x.raw, x.raw.zip, x.mzML.gz ...)
        _ms_match = df_study['filename_lower'].str.extract(_MS_FILE_RE)
        df_study = df_study[_ms_match[0].notna()].copy()

        print(f'Assay table has {len(df_study)} rows and {len(df_study.columns)} columns.')

        if len(df_study) == 0:
            # e.g. only bare .zip archives or .xls tables listed, nothing REDU can index
            return pd.DataFrame()

        # 'extension' is the MS format (inner format for a compressed file)
        df_study['extension'] = _ms_match.loc[df_study.index, 0]
        df_study['compressed'] = _ms_match.loc[df_study.index, 1].notna()

        df_study = df_study.groupby('row_id').apply(prefer_extension)
        df_study = df_study[df_study['keep']]

        if len(df_study) > 0:
            
            allowedTerm_dict = kwargs['allowedTerm_dict']
            ontology_table = kwargs['ontology_table']
            ENVOEnvironmentMaterialIndex_table = kwargs['ENVOEnvironmentMaterialIndex_table']
            ENVOEnvironmentBiomeIndex_table = kwargs['ENVOEnvironmentBiomeIndex_table']
            NCBIRankDivision_table  = kwargs['NCBIRankDivision_table']

            
            df_study.loc[:, 'YearOfAnalysis'] = submissionYear

            #add NCBITaxonomy and Sampletype & SampleTypeSub1
            #######
            # Coalesce the organism source: the block below is driven by
            # 'Samples_Organism', so studies that record the organism under a
            # different column (Species/Genus/Bacteria/Strain) would otherwise get no
            # NCBITaxonomy at all. Fill 'Samples_Organism' per-sample from those
            # aliases ONLY where it is blank, so an explicit Organism value always
            # wins (no existing value is changed).
            _tax_alias_cols = ['Samples_Species', 'Samples_species', 'Samples_Bacteria',
                               'Samples_Genus', 'Samples_genus', 'Samples_organism']
            # NCBITaxon accession recorded with the organism (taken from the same column
            # as the organism value)
            _tax_accession = pd.Series('', index=df_study.index, dtype=object)
            _acc_col = _accession_column(df_study, 'Samples_Organism')
            if _acc_col:
                _tax_accession = _accession_ids(df_study[_acc_col], 'NCBITaxon')
            if any(c in df_study.columns for c in _tax_alias_cols):
                if 'Samples_Organism' not in df_study.columns:
                    df_study['Samples_Organism'] = np.nan
                for _c in _tax_alias_cols:
                    if _c in df_study.columns:
                        _blank = df_study['Samples_Organism'].isna() | df_study['Samples_Organism'].astype(str).str.strip().str.lower().isin(['', 'nan', 'none', 'missing value'])
                        df_study.loc[_blank, 'Samples_Organism'] = df_study.loc[_blank, _c]
                        _alias_acc = _accession_column(df_study, _c)
                        if _alias_acc:
                            _tax_accession.loc[_blank] = _accession_ids(df_study.loc[_blank, _alias_acc], 'NCBITaxon')

            if 'Samples_Organism' in df_study.columns:
                processed_organisms = {org: str(get_taxonomy_id_from_name__allowedTerms(org, allowedTerm_dict = allowedTerm_dict)) for org in df_study['Samples_Organism'].unique()}
                # Organism names that NCBI lists as synonym or common name of exactly one taxon
                # ('Clostridium bolteae', 'House Mouse'): MetaboLights alias sheet built from names.dmp
                _org_alias = kwargs.get('organism_aliases') or {}
                for _org, _tax in processed_organisms.items():
                    if _tax == 'None' and str(_org).strip().lower() in _org_alias:
                        processed_organisms[_org] = _org_alias[str(_org).strip().lower()]

                df_study.loc[:, 'NCBITaxonomy'] = df_study['Samples_Organism'].map(processed_organisms)
                # 'None' is the name lookup without a match; taxa whose name contains 'None'
                # ('203745|Nonea', '1068978|Nonomuraea ...') are allowed terms
                df_study.loc[:, 'NCBITaxonomy'] = df_study['NCBITaxonomy'].replace('None', 'missing value')

                # The curator-chosen NCBITaxon accession is more reliable than resolving the
                # free-text organism name, so it is used whenever it is a valid taxon
                if (_tax_accession != '').any():
                    _by_id = _ncbi_term_by_id(allowedTerm_dict)
                    _from_acc = _tax_accession.str.replace('NCBITaxon:', '', regex=False).map(_by_id)
                    _use = _from_acc.notna()
                    df_study.loc[_use, 'NCBITaxonomy'] = _from_acc[_use]
                
                df_study.loc[:, ['SampleType', 'SampleTypeSub1']] = 'missing value'

                processed_taxonomy = {taxonomy.split('|')[0]: get_taxonomy_info(taxonomy.split('|')[0])
                                    for taxonomy in df_study['NCBITaxonomy'].unique()
                                    if '|' in taxonomy}
                
                print(f'Found {len(processed_taxonomy)} unique taxonomy IDs.')


                df_study.loc[:,'SampleType'] = df_study['NCBITaxonomy'].map(lambda x: processed_taxonomy.get(x.split('|')[0], [pd.NA, pd.NA])[0] 
                                                                    if '|' in x else pd.NA)
                df_study.loc[:,'SampleTypeSub1'] = df_study['NCBITaxonomy'].map(lambda x: processed_taxonomy.get(x.split('|')[0], [pd.NA, pd.NA])[1]
                                                                        if '|' in x else pd.NA)

                df_study[['SampleType', 'SampleTypeSub1']] = df_study.apply(
                    lambda row: get_blanks(row.Samples_Organism) 
                    if (pd.isna(row.SampleType) or row.SampleType == 'missing value') 
                    else [row.SampleType, row.SampleTypeSub1], 
                    axis=1
                ).apply(pd.Series)

                df_study[['SampleType', 'SampleTypeSub1']] = df_study.apply(
                    lambda row: get_enviromental_water(row.Samples_Organism) 
                    if (pd.isna(row.SampleType) or row.SampleType == 'missing value') 
                    else [row.SampleType, row.SampleTypeSub1], 
                    axis=1
                ).apply(pd.Series)

                # Pooled QCs, blanks and standards are labelled in the 'Sample type'
                # characteristic (MetaboLights vocabulary), while their organism is often that
                # of the study (e.g. a pooled human plasma QC), so this label takes precedence.
                # Labels the converter did not read before are matched exactly; for blanks and
                # standards they are only taken where the organism column does not name an
                # organism (empty, 'blank', 'reference compound', NCBI 'blank sample'): real
                # samples are labelled 'experimental blank' or 'sample preparation blank' in
                # many studies ('experimental blank' is the only value in 119 studies), and
                # 'reference compound mix' is used for pooled QCs of the study material.
                # Factor values ('Treatment = QC') are matched the same way.
                _org_lab = df_study['Samples_Organism'].map(_label)
                _not_organism = (_is_missing(df_study['Samples_Organism']) | _org_lab.isin(_QC_BLANK_LABELS) | _org_lab.isin(_ORGANISM_BLANK_LABELS)
                                 | df_study['NCBITaxonomy'].astype(str).str.startswith('2582415|'))
                # a 'QC' of solvent or of a compound mixture is not a pooled QC of the samples
                _part_lab = df_study['Samples_Organism part'].map(_label) if 'Samples_Organism part' in df_study.columns else pd.Series('', index=df_study.index)
                _not_pool = _part_lab.isin({'solvent', 'mixture', 'water', 'buffer'}) | _part_lab.map(lambda v: _QC_BLANK_LABELS.get(v, '').startswith('blank'))
                def _new_label(series):
                    t = series.map(lambda v: _qc_blank_term(v) or ('blank_analysis' if _label(v) == 'experimental blank' else None))
                    return t.where(((t == 'pool_QC') & ~_not_pool) | _not_organism)
                qc_label = pd.Series(None, index=df_study.index, dtype=object)
                for _c in [c for c in df_study.columns if re.fullmatch(r'Samples_Sample type(?:_\d+)?', str(c), re.I)]:
                    _qc = df_study[_c].map(_sample_type_qc)
                    qc_label = qc_label.where(qc_label.notna(), _qc.where(_qc.notna(), _new_label(df_study[_c])))
                for _c in [c for c in factor_columns if c in df_study.columns and not re.search(r'diet|feed|food', c, re.I)]:
                    qc_label = qc_label.where(qc_label.notna(), _new_label(df_study[_c]))
                _qc = df_study['Samples_Organism'].map(lambda v: _qc_blank_term(v) or _ORGANISM_BLANK_LABELS.get(_label(v)))
                qc_label = qc_label.where(qc_label.notna(), _qc.where(_not_organism))
                if 'Samples_Organism part' in df_study.columns:
                    qc_label = qc_label.where(qc_label.notna(), _new_label(df_study['Samples_Organism part']))
                # a generic blank label does not replace a more specific blank type ('Blank Extraction')
                qc_label = qc_label.where(~((qc_label == 'blank_analysis') & df_study['SampleType'].isin(['blank_extraction', 'blank_culturemedia'])),
                                          df_study['SampleType'])
                _set = qc_label.notna()
                if _set.any():
                    df_study.loc[_set, 'SampleType'] = qc_label[_set]
                    df_study.loc[_set, 'SampleTypeSub1'] = qc_label[_set].map(lambda t: _SUBTYPE.get(t, t))
                # kept as a column: the index is reset by the ontology merges below
                df_study['_qc_label'] = qc_label


                # df_biofluid_vs_tissue = pd.read_csv(os.path.join(transSheet_dir, 'biofluid_vs_tissue_Metabolights.csv'))

                # df_study['Samples_Organism part'] = df_study['Samples_Organism part'].str.lower()

                # for index, row in df_study.iterrows():
                #     if pd.isna(row['SampleTypeSub1']):
                #         org_part = row['Samples_Organism part']
                #         matching_row = df_biofluid_vs_tissue[df_biofluid_vs_tissue['ML'] == org_part]
                #         if not matching_row.empty:
                #             df_study.at[index, 'SampleTypeSub1'] = matching_row.iloc[0]['tissueVSbiofluid']


                #add ENV material column
                #######
                updated_species_dict = {}

                envo_mat_labels_norm = ENVOEnvironmentMaterialIndex_table['Label'].str.replace(' ', '').str.lower()
                envo_mat_synonyms_norm = ENVOEnvironmentMaterialIndex_table['Synonym'].str.replace(' ', '').str.lower()

                for key, value in processed_organisms.items():
                    if value == 'None':
                        match = ENVOEnvironmentMaterialIndex_table[(ENVOEnvironmentMaterialIndex_table['Label'] == key) | (ENVOEnvironmentMaterialIndex_table['Synonym'] == key)]['Label'].values
                        if match.size > 0:
                            updated_species_dict[key] = match[0]
                        else:
                            key_norm = key.replace(' ', '').lower()
                            match_norm = ENVOEnvironmentMaterialIndex_table[(envo_mat_labels_norm == key_norm) | (envo_mat_synonyms_norm == key_norm)]['Label'].values
                            updated_species_dict[key] = match_norm[0] if match_norm.size > 0 else 'missing value'

                df_study.loc[:, 'ENVOEnvironmentMaterial'] = df_study['Samples_Organism'].map(updated_species_dict)


                #add ENV biome column
                #######
                updated_species_dict = {}

                envo_biome_labels_norm = ENVOEnvironmentBiomeIndex_table['Label'].str.replace(' ', '').str.lower()
                envo_biome_synonyms_norm = ENVOEnvironmentBiomeIndex_table['Synonym'].str.replace(' ', '').str.lower()

                for key, value in processed_organisms.items():
                    if value == 'None':
                        match = ENVOEnvironmentBiomeIndex_table[(ENVOEnvironmentBiomeIndex_table['Label'] == key) | (ENVOEnvironmentBiomeIndex_table['Synonym'] == key)]['Label'].values
                        if match.size > 0:
                            updated_species_dict[key] = match[0]
                        else:
                            key_norm = key.replace(' ', '').lower()
                            match_norm = ENVOEnvironmentBiomeIndex_table[(envo_biome_labels_norm == key_norm) | (envo_biome_synonyms_norm == key_norm)]['Label'].values
                            updated_species_dict[key] = match_norm[0] if match_norm.size > 0 else 'missing value'

                df_study.loc[:, 'ENVOEnvironmentBiome'] = df_study['Samples_Organism'].map(updated_species_dict)


                #add UBERON bodypart column
                #######
                # Coalesce the body-part source: only 'Samples_Organism part' was read,
                # so organ/tissue/plant-part recorded under other columns was dropped.
                # Fill it per-sample from those aliases where blank (Organism part wins).
                _bp_alias_cols = ['Samples_Organ', 'Samples_organ', 'Samples_tissue type',
                                  'Samples_Tissue', 'Samples_tissue', 'Samples_Plant part',
                                  'Samples_Biospecimen Type', 'Samples_Cell type']
                if 'Samples_Organism part' not in df_study.columns:
                    df_study['Samples_Organism part'] = np.nan
                # UBERON/CL/PO accession recorded with the body part (from the same column
                # as the body part value)
                def _body_part_accession(column):
                    _acc = _accession_column(df_study, column)
                    if not _acc:
                        return pd.Series('', index=df_study.index, dtype=object)
                    ids = _accession_ids(df_study[_acc], 'UBERON')
                    for _ont in ('CL', 'PO'):
                        ids = ids.where(ids != '', _accession_ids(df_study[_acc], _ont))
                    return ids
                _bp_accession = _body_part_accession('Samples_Organism part')
                for _c in _bp_alias_cols:
                    if _c in df_study.columns:
                        _blank = df_study['Samples_Organism part'].isna() | df_study['Samples_Organism part'].astype(str).str.strip().str.lower().isin(['', 'nan', 'none', 'missing value'])
                        df_study.loc[_blank, 'Samples_Organism part'] = df_study.loc[_blank, _c]
                        _bp_accession.loc[_blank] = _body_part_accession(_c).loc[_blank]

                df_study['UBERONBodyPartName'] = 'missing value'

                labels_in_ontology_table = set(ontology_table['Label'])
                synonyms_in_ontology_table = set(ontology_table['Synonym'])
                allowed_bodyparts = set(allowedTerm_dict['UBERONBodyPartName']['allowed_values'])
                allowed_bodypart_ids = set(allowedTerm_dict['UBERONOntologyIndex']['allowed_values'])
                _label_by_id = dict(zip(ontology_table['UBERONOntologyIndex'].astype(str).str.replace('_', ':', regex=False),
                                        ontology_table['Label']))
                _synonyms_lower = ontology_table['Synonym'].astype(str).str.strip().str.lower()
                _labels_lower = ontology_table['Label'].astype(str).str.strip().str.lower()
                _case_insensitive_matches = {}

                for index, row in df_study.iterrows():
                    organism_part = row['Samples_Organism part']

                    _acc_id = _bp_accession.get(index, '')
                    if _acc_id and _acc_id in allowed_bodypart_ids and _label_by_id.get(_acc_id) in allowed_bodyparts:
                        df_study.at[index, 'UBERONBodyPartName'] = _label_by_id[_acc_id]
                        continue

                    if organism_part in allowed_bodyparts:
                        df_study.at[index, 'UBERONBodyPartName'] = organism_part
                        continue 

                    matching_rows = ontology_table.iloc[0:0]
                    if organism_part in synonyms_in_ontology_table:
                        matching_rows = ontology_table[ontology_table['Synonym'] == organism_part]
                    elif isinstance(organism_part, str) and organism_part.strip():
                        # same label/exact synonym, written with different capitalisation
                        _key = organism_part.strip().lower()
                        if _key not in _case_insensitive_matches:
                            _case_insensitive_matches[_key] = ontology_table[(_labels_lower == _key) | (_synonyms_lower == _key)].drop_duplicates(subset=['UBERONOntologyIndex'])
                        matching_rows = _case_insensitive_matches[_key]

                    if len(matching_rows) > 0:

                        plant = True if row['SampleType'] == 'plant' else False

                        if len(matching_rows) > 1 and plant == True:
                            matching_rows = matching_rows[matching_rows['UBERONOntologyIndex'].str.startswith('PO', na=False)]

                        if len(matching_rows) > 1 and plant == False:
                            matching_rows = matching_rows[(~matching_rows['UBERONOntologyIndex'].str.startswith('PO', na=False))]

                        if len(matching_rows) == 1:
                            
                            label = matching_rows['Label'].iloc[0] 
                            df_study.at[index, 'UBERONBodyPartName'] = label

                        
                # Common names that are neither UBERON label nor exact synonym ('plasma' is not a
                # synonym of 'blood plasma'), only for animal samples
                _alias = df_study['Samples_Organism part'].astype(str).str.strip().str.lower().map(_BODY_PART_ALIASES)
                _alias = _alias.where(_alias.isin(allowed_bodyparts))
                _use = (df_study['UBERONBodyPartName'] == 'missing value') & _alias.notna() & (df_study['SampleType'] == 'animal')
                df_study.loc[_use, 'UBERONBodyPartName'] = _alias[_use]

                ontology_table_unique = ontology_table.drop_duplicates(subset=['Label', 'UBERONOntologyIndex']).drop(columns=['Synonym'])

                df_study = pd.merge(df_study, ontology_table_unique, left_on='UBERONBodyPartName', right_on='Label', how='left')
                df_study.drop('Label', axis=1, inplace=True)

                df_study.loc[df_study['Is Fluid'] == True, 'SampleTypeSub1'] = 'biofluid'
                df_study.loc[df_study['Is Multicellular'] == True, 'SampleTypeSub1'] = 'tissue'
                # QC samples, blanks and standards keep their own SampleTypeSub1
                _qc_rows = df_study['_qc_label'].notna()
                df_study.loc[_qc_rows, 'SampleTypeSub1'] = df_study.loc[_qc_rows, '_qc_label'].map(lambda t: _SUBTYPE.get(t, t))
                # Chemical standards are not a sample of the study organism/body part/environment
                _std = df_study['_qc_label'] == 'molecular_standards'
                for _c in ['NCBITaxonomy', 'UBERONBodyPartName', 'ENVOEnvironmentMaterial', 'ENVOEnvironmentBiome']:
                    df_study.loc[_std, _c] = 'missing value'
                df_study.loc[_std, 'UBERONOntologyIndex'] = np.nan

                # SampleCollectionMethod terms that a blood or urine sample implies by definition
                # ('blood NOS', 'urine, NOS': collected, method not specified). Feces is left
                # out: 'feces' is one of several feces collection methods (swabs, frozen, ...)
                _scm = df_study['UBERONBodyPartName'].map({'blood': 'blood NOS', 'blood plasma': 'blood NOS',
                                                           'blood serum': 'blood NOS', 'urine': 'urine, NOS'})
                _scm = _scm.where(_scm.isin(allowedTerm_dict['SampleCollectionMethod']['allowed_values']) & df_study['_qc_label'].isna())
                if _scm.notna().any():
                    if 'SampleCollectionMethod' not in df_study.columns:
                        df_study['SampleCollectionMethod'] = 'missing value'
                    _need = _is_missing(df_study['SampleCollectionMethod']) & _scm.notna()
                    df_study.loc[_need, 'SampleCollectionMethod'] = _scm[_need]

                # UBERON does not flag every body fluid or excretion (urine, feces, follicular
                # fluid, ...) as such; the MetaboLights translation sheet covers those
                _bvt = kwargs.get('biofluid_vs_tissue')
                if _bvt is not None and len(_bvt):
                    _bvt_map = dict(zip(_bvt['ML'].astype(str).str.strip().str.lower(), _bvt['tissueVSbiofluid']))
                    _sheet = df_study['UBERONBodyPartName'].astype(str).str.strip().str.lower().map(_bvt_map)
                    _sheet = _sheet.where(_sheet.notna(), df_study['Samples_Organism part'].astype(str).str.strip().str.lower().map(_bvt_map))
                    _need = _is_missing(df_study['SampleTypeSub1']) & _sheet.notna() & (df_study['SampleType'] == 'animal')
                    df_study.loc[_need, 'SampleTypeSub1'] = _sheet[_need]
                df_study['UBERONOntologyIndex'] = df_study['UBERONOntologyIndex'].str.replace("_", ":", regex=False)


            if '_qc_label' not in df_study.columns:
                df_study['_qc_label'] = None

            #add MassSpectrometer column
            #######
            if 'Assay_Instrument' in df_study.columns:
                df_study['Assay_Instrument'] = df_study['Assay_Instrument'].fillna('')
                _instrument_name = df_study['Assay_Instrument'].astype(str)
                _instrument_acc_col = _accession_column(df_study, 'Assay_Instrument')
                _instrument_acc = df_study[_instrument_acc_col].astype(str) if _instrument_acc_col else pd.Series('', index=df_study.index)

                prefixes = [
                    "Thermo Scientific ",
                    "Waters ",
                    "Agilent ",
                    "AB SCIEX ",
                    "LECO ",
                    "Bruker ",
                    "Shimadzu "
                ]

                # Function to remove prefixes
                def remove_prefixes(value, prefixes):
                    for prefix in prefixes:
                        if value.startswith(prefix):
                            return value[len(prefix):]
                    return value

                # Apply the function to the column
                df_study['Assay_Instrument'] = df_study['Assay_Instrument'].apply(lambda x: remove_prefixes(x, prefixes))

                allowed_values_map = dict(zip(
                    allowedTerm_dict["MassSpectrometer"]["allowed_values_matching_0"],
                    allowedTerm_dict["MassSpectrometer"]["allowed_values"]
                ))

                df_study["MassSpectrometer"] = df_study["Assay_Instrument"].map(allowed_values_map)

                # Then the PSI-MS accession recorded with the instrument, then the shared
                # name matcher (vendor prefixes, abbreviations; only unambiguous matches)
                _unres = df_study["MassSpectrometer"].isna()
                if _unres.any():
                    _allowed_ms_terms = allowedTerm_dict["MassSpectrometer"]["allowed_values"]
                    _ms_matches = {}
                    def _match_ms(name, accession):
                        if (name, accession) not in _ms_matches:
                            _ms_matches[(name, accession)] = (match_mass_spectrometer(None, _allowed_ms_terms, accession=accession)
                                                              or match_mass_spectrometer(name, _allowed_ms_terms))
                        return _ms_matches[(name, accession)]
                    df_study.loc[_unres, "MassSpectrometer"] = [
                        _match_ms(n, a) for n, a in zip(_instrument_name[_unres], _instrument_acc[_unres])]

                # Fallback: the exact map misses ~40% of instruments (Q-Exactive vs
                # "Q Exactive", Agilent iFunnel Q-TOF, Synapt/maXis variants, ...). Map
                # the (already prefix-stripped) instrument string onto the allowed
                # vocabulary by token-subset match. Analysis-level, no over-assignment.
                _allowed_ms = allowedTerm_dict["MassSpectrometer"]["allowed_values"]
                _unres = df_study["MassSpectrometer"].isna()
                if _unres.any():
                    _ms_cache = {}
                    def _resolve_ms(raw):
                        if raw not in _ms_cache:
                            _ms_cache[raw] = map_instrument_to_allowed(raw, _allowed_ms)
                        return _ms_cache[raw]
                    df_study.loc[_unres, "MassSpectrometer"] = df_study.loc[_unres, "Assay_Instrument"].map(_resolve_ms)

                # Instrument names written in ways none of the matchers recognise
                # ('QE Plus', 'Agilent 6230 Accurate-Mass TOF'): explicit alias sheet
                _ms_alias = kwargs.get('ms_aliases') or {}
                _unres = df_study["MassSpectrometer"].isna()
                if _ms_alias and _unres.any():
                    df_study.loc[_unres, "MassSpectrometer"] = _instrument_name[_unres].str.strip().map(_ms_alias)



            #add IonizationMethod/Polarity column
            #######
            if 'Assay_Ion source' in df_study.columns and 'Assay_Scan polarity' in df_study.columns:
                df_study['Assay_Ion source'] = df_study['Assay_Ion source'].fillna('')
                df_study['Assay_Scan polarity'] = df_study['Assay_Scan polarity'].fillna('')
                _ion_acc_col = _accession_column(df_study, 'Assay_Ion source')
                _ion_acc = _accession_ids(df_study[_ion_acc_col], 'MS') if _ion_acc_col else pd.Series('', index=df_study.index)

                df_study['Assay_Ion source'] = df_study['Assay_Ion source'].str.lower()
                df_study['Assay_Scan polarity'] = df_study['Assay_Scan polarity'].astype(str).str.lower()
                df_study['IonizationSourceAndPolarity'] = ''
                df_study.loc[df_study['Assay_Ion source'].str.contains('electrospray'), 'IonizationSourceAndPolarity'] = 'electrospray ionization'
                df_study.loc[df_study['Assay_Ion source'].str.contains('electron ionization'), 'IonizationSourceAndPolarity'] = 'electron ionization'      
                df_study.loc[df_study['Assay_Ion source'].str.contains('chemical ionization'), 'IonizationSourceAndPolarity'] = 'atmospheric pressure chemical ionization'
                # Abbreviations and vendor source names (ESI, HESI, Turbo V/TurboIonSpray,
                # IonDrive, APCI, APPI, EI), then the PSI-MS accession of the ion source
                _unset = df_study['IonizationSourceAndPolarity'] == ''
                _src = df_study['Assay_Ion source']
                for _pattern, _term in [(r'\bn?h?esi\b|nano ?spray|turbo ?(?:v|ion ?spray)|ion ?spray|iondrive|jet stream', 'electrospray ionization'),
                                        (r'\bapci\b', 'atmospheric pressure chemical ionization'),
                                        (r'\bappi\b|photoionization|photo-ionization', 'atmospheric pressure photoionization'),
                                        (r'\bei\b|electron impact', 'electron ionization')]:
                    _hit = _unset & _src.str.contains(_pattern, regex=True)
                    df_study.loc[_hit, 'IonizationSourceAndPolarity'] = _term
                    _unset = _unset & ~_hit
                for _accession, _term in [('MS:1000073', 'electrospray ionization'), ('MS:1000398', 'electrospray ionization'),
                                          ('MS:1000070', 'atmospheric pressure chemical ionization'),
                                          ('MS:1000382', 'atmospheric pressure photoionization'),
                                          ('MS:1000389', 'electron ionization')]:
                    _hit = _unset & (_ion_acc == _accession)
                    df_study.loc[_hit, 'IonizationSourceAndPolarity'] = _term
                    _unset = _unset & ~_hit
                _negative = df_study['Assay_Scan polarity'].str.contains(r'negative|\bneg\b', regex=True)
                _positive = df_study['Assay_Scan polarity'].str.contains(r'posi?tive|\bpos\b', regex=True)
                _with_polarity = (df_study['IonizationSourceAndPolarity'] != '') & ~df_study['IonizationSourceAndPolarity'].str.contains("electron ionization")
                df_study.loc[_negative & _with_polarity, 'IonizationSourceAndPolarity'] += ' (negative)'
                df_study.loc[_positive & _with_polarity, 'IonizationSourceAndPolarity'] += ' (positive)'
                df_study.loc[df_study['Assay_Scan polarity'].str.contains('alternating') & _with_polarity, 'IonizationSourceAndPolarity'] += ' (alternating)'


            #add LC method
            #######
            if 'Assay_Column type' in df_study.columns:
                df_study['Assay_Column type'] = df_study['Assay_Column type'].fillna('')
                df_study['Assay_Column type'] = df_study['Assay_Column type'].str.lower()
                df_study['ChromatographyAndPhase'] = ''
                df_study.loc[df_study['Assay_Column type'].str.contains('reverse'), 'ChromatographyAndPhase'] = 'reverse phase'
                df_study.loc[df_study['Assay_Column type'].str.contains('hilic'), 'ChromatographyAndPhase'] = 'normal phase (HILIC)'
                df_study.loc[df_study['Assay_Column type'].str.contains('normal phase'), 'ChromatographyAndPhase'] = 'normal phase (HILIC)'
                
                if 'Assay_Column model' in df_study.columns:
                    df_study['Assay_Column model'] = df_study['Assay_Column model'].fillna('')
                    df_study.loc[df_study['Assay_Column model'].str.contains('Phenyl', case=False) & df_study['Assay_Column model'].str.contains('Hexyl', case=False), 'ChromatographyAndPhase'] += ' (Phenyl-Hexyl)'
                    df_study.loc[df_study['Assay_Column model'].str.contains('C18') & df_study['Assay_Column model'].str.contains('polar', case=False), 'ChromatographyAndPhase'] += ' (polar-C18)'
                    
                    df_study.loc[df_study['Assay_Column model'].str.contains('HSS T3') & ~df_study['ChromatographyAndPhase'].str.contains(r"\("), 'ChromatographyAndPhase'] += ' (C18)'
                    df_study.loc[df_study['Assay_Column model'].str.contains('C18') & ~df_study['ChromatographyAndPhase'].str.contains(r"\("), 'ChromatographyAndPhase'] += ' (C18)'
                    df_study.loc[df_study['Assay_Column model'].str.contains('C30') & ~df_study['ChromatographyAndPhase'].str.contains(r"\("), 'ChromatographyAndPhase'] += ' (C30)'
                    df_study.loc[df_study['Assay_Column model'].str.contains('C8') & ~df_study['ChromatographyAndPhase'].str.contains(r"\("), 'ChromatographyAndPhase'] += ' (C8)'
                    df_study.loc[df_study['ChromatographyAndPhase'].str.contains('reverse phase') & ~df_study['ChromatographyAndPhase'].str.contains(r"\("), 'ChromatographyAndPhase'] += ' (NOS)'
                else:
                    df_study.loc[df_study['ChromatographyAndPhase'].str.contains('reverse phase') & ~df_study['ChromatographyAndPhase'].str.contains(r"\("), 'ChromatographyAndPhase'] += ' (NOS)'

            # Gas chromatography: GC studies often leave 'Column type' empty, but name a GC
            # in the chromatography instrument or use GC column/type wording
            _gc_text = pd.Series('', index=df_study.index, dtype=object)
            for _c in ['Assay_Chromatography Instrument', 'Assay_Column type', 'Assay_Column model']:
                if _c in df_study.columns:
                    _gc_text = _gc_text + ' ' + df_study[_c].fillna('').astype(str)
            _is_gc = _gc_text.str.contains(r'\bgc\b|gas chromatograph|\bgc[- ]?ms\b|\bgcxgc\b', case=False, regex=True)
            # GC capillary columns named by model (5%-phenyl, wax, FFAP, 1/17-type phases), or the
            # MetaboLights GC column types 'low/medium/high polarity', unless an LC column is named
            _model = df_study['Assay_Column model'].fillna('').astype(str) if 'Assay_Column model' in df_study.columns else pd.Series('', index=df_study.index)
            _ctype = df_study['Assay_Column type'].fillna('').astype(str).map(_label) if 'Assay_Column type' in df_study.columns else pd.Series('', index=df_study.index)
            _lc_words = _model.str.contains(r'c18|\bc8\b|hilic|\bbeh\b|\bhss\b|amide|reverse|uplc|hplc|zorbax|kinetex|acquity|luna|atlantis|hypersil|accucore|synergi|cortecs|ion ?pac|\bt3\b', case=False, regex=True)
            _gc_model = _model.str.contains(r'(?<![a-z0-9])(?:db|hp|rxi|rtx|zb|vf|dp|tg|bpx?|sh-rxi|elite)-?(?:1|5|17|35|wax)(?:[a-z]{0,3}|-?ms|sil ?ms)?(?![0-9])|\bffap\b|\bdb-?wax|\bstabilwax', case=False, regex=True)
            _is_gc = _is_gc | ((_gc_model | _ctype.isin(['low polarity', 'medium polarity', 'high polarity'])) & ~_lc_words)
            if _is_gc.any():
                if 'ChromatographyAndPhase' not in df_study.columns:
                    df_study['ChromatographyAndPhase'] = ''
                _need = _is_missing(df_study['ChromatographyAndPhase']) & _is_gc
                _db5 = _gc_text.str.contains(r'\bdb-?5', case=False, regex=True)
                df_study.loc[_need & _db5, 'ChromatographyAndPhase'] = 'gas chromatography (DB-5)'
                df_study.loc[_need & ~_db5, 'ChromatographyAndPhase'] = 'gas chromatography (NOS)'

            #add AgeInYears 
            #######
            if 'samples_age' in df_study.columns.str.lower():
                if 'Samples_Age' in df_study.columns:
                    df_study.rename(columns={'Samples_Age': 'Samples_age'}, inplace=True)
                if 'Samples_AGE' in df_study.columns:
                    df_study.rename(columns={'Samples_AGE': 'Samples_age'}, inplace=True)

                _age_text = df_study['Samples_age'].copy() if 'Samples_age' in df_study.columns else pd.Series(np.nan, index=df_study.index)
                unit_column = find_column_after_target_column(df_study, 'Samples_age', 'Samples_Unit')

                if unit_column != '':
                    df_study['Samples_age'] = pd.to_numeric(df_study['Samples_age'], errors='coerce')
                    df_study['AgeInYears'] = np.nan
                    # units are written 'year', 'Year', 'Month', 'weeks', ...
                    _unit = df_study[unit_column].map(_unit_name)
                    for _u, _div in _UNIT_YEARS.items():
                        df_study.loc[_unit == _u, 'AgeInYears'] = df_study.loc[_unit == _u, 'Samples_age'] / _div

                    df_study['AgeInYears'] = df_study['AgeInYears'].astype(str).replace('nan', 'missing value')

            #add AgeInYears 
            #######
                    df_study['LifeStage'] = df_study['AgeInYears'].apply(lambda x: age_category(x))

                # The unit can also be written in the value itself ('73 year', '7 weeks')
                _age_from_value = _age_text.map(_age_in_years_from_text)
                if _age_from_value.notna().any():
                    if 'AgeInYears' not in df_study.columns:
                        df_study['AgeInYears'] = 'missing value'
                    _need = _is_missing(df_study['AgeInYears']) & _age_from_value.notna()
                    df_study.loc[_need, 'AgeInYears'] = _age_from_value[_need].astype(str)
                    if 'LifeStage' not in df_study.columns:
                        df_study['LifeStage'] = ''
                    df_study.loc[_need, 'LifeStage'] = _age_from_value[_need].apply(age_category)

            #add Sex
            #######
            # Consume any sex / gender / biological-sex column (the original only read
            # 'Samples_gender', silently dropping 'Samples_sex' and 'Samples_Biological
            # sex'). Exclude columns describing a DIFFERENT subject than the sample
            # (baby/fetus/maternal/paternal/donor/partner/contact).
            _sex_exclude = re.compile(r'baby|fetus|foetal|fetal|maternal|mother|paternal|father|donor|partner|contact|infant|neonat', re.I)
            def _tokens(col):
                return set(re.split(r'[^a-z0-9]+', str(col).lower()))
            _sex_cols = [c for c in df_study.columns
                         if str(c).startswith('Samples_')
                         and ({'sex', 'gender'} & _tokens(c))
                         and not _sex_exclude.search(str(c))]
            if _sex_cols:
                _sexsrc = pd.Series('', index=df_study.index, dtype=object)
                for c in _sex_cols:
                    _b = _sexsrc.astype(str).str.strip().str.lower().isin(['', 'nan', 'none'])
                    _sexsrc.loc[_b] = df_study[c].astype(str).loc[_b]
                _sl = _sexsrc.astype(str).str.lower()
                if 'BiologicalSex' not in df_study.columns:
                    df_study['BiologicalSex'] = 'missing value'
                _blank = df_study['BiologicalSex'].astype(str).str.lower().isin(['', 'nan', 'missing value', 'none'])
                df_study.loc[_blank & _sl.str.contains('female'), 'BiologicalSex'] = 'female'
                df_study.loc[_blank & _sl.str.contains(r'\bmale\b') & (df_study['BiologicalSex'] != 'female'), 'BiologicalSex'] = 'male'
                # abbreviations and other words for the two sexes ('F', 'M', 'Woman', 'Man')
                _abbr = _sl.str.strip().map({'f': 'female', 'm': 'male', 'woman': 'female', 'women': 'female',
                                             'man': 'male', 'men': 'male', 'girl': 'female', 'boy': 'male'})
                _need = _is_missing(df_study['BiologicalSex']) & _abbr.notna()
                df_study.loc[_need, 'BiologicalSex'] = _abbr[_need]


            #add LatitudeandLongitude
            #######
            # Match latitude/longitude by column NAME (paired), covering MetaboLights
            # variants like "Geographic Location Latitude" / "Geographic location
            # (latitude)" - never by value (a bare number could be a plate/time).
            _lat_cols = [c for c in df_study.columns if 'latitude' in str(c).lower()]
            _lon_cols = [c for c in df_study.columns if 'longitude' in str(c).lower()]
            if _lat_cols and _lon_cols:
                _latc, _lonc = _lat_cols[0], _lon_cols[0]
                def combine_lat_lon(row):
                    lat = str(row[_latc]).strip()
                    lon = str(row[_lonc]).strip()
                    if lat and lon and lat.lower() not in ('nan', '') and lon.lower() not in ('nan', ''):
                        try:
                            la = float(lat); lo = float(lon)
                            if -90 <= la <= 90 and -180 <= lo <= 180 and not (la == 0 and lo == 0):
                                return f"{lat}|{lon}"
                        except ValueError:
                            pass
                    return 'missing value'
                _ll = df_study.apply(combine_lat_lon, axis=1)
                if 'LatitudeandLongitude' not in df_study.columns:
                    df_study['LatitudeandLongitude'] = _ll
                else:
                    _blank = df_study['LatitudeandLongitude'].astype(str).str.lower().isin(['', 'nan', 'missing value', 'none'])
                    df_study.loc[_blank, 'LatitudeandLongitude'] = _ll[_blank]

            #add ENVOEnvironmentMaterial (value-matched from candidate columns)
            #######
            # Environmental studies record the material (soil/seawater/sediment/...) in
            # dedicated substrate/material columns. Match values against the ENVO
            # material vocabulary; fill only where currently missing. NOTE: we do NOT
            # scan the organism column here (organism -> ENVO material is handled
            # separately above via the ontology table) - matching organisms against the
            # material vocab risks false hits (e.g. a stray taxon name in the vocab).
            _envm_vocab = {str(v).strip().lower(): str(v).strip()
                           for v in allowedTerm_dict.get('ENVOEnvironmentMaterial', {}).get('allowed_values', []) if str(v).strip()}
            if _envm_vocab:
                _mat_cols = [c for c in ['Samples_Onsite Substrate', 'Samples_Overlaying substrate',
                                         'Samples_Material', 'Samples_material', 'Samples_Environmental material',
                                         'Samples_environment_material', 'Samples_env_material', 'Samples_Substrate',
                                         'Samples_Environment (material)'] if c in df_study.columns]
                if _mat_cols:
                    if 'ENVOEnvironmentMaterial' not in df_study.columns:
                        df_study['ENVOEnvironmentMaterial'] = 'missing value'
                    for _c in _mat_cols:
                        _m = df_study[_c].astype(str).str.strip().str.lower().map(_envm_vocab)
                        _need = df_study['ENVOEnvironmentMaterial'].astype(str).str.strip().str.lower().isin(['', 'nan', 'none', 'missing value']) & _m.notna()
                        df_study.loc[_need, 'ENVOEnvironmentMaterial'] = _m[_need]

            # Environmental samples without an organism name the material as 'Organism part'
            # ('soil', 'sediment', 'waste water'). Exact ENVO material labels only, and only
            # where the organism is not a taxon; names that are also anatomy (UBERON label or
            # synonym, e.g. 'plasma', which ENVO has as a state of matter) are never used.
            if _envm_vocab and 'Samples_Organism part' in df_study.columns and 'NCBITaxonomy' in df_study.columns:
                _anatomy = set(ontology_table['Label'].astype(str).str.strip().str.lower()) | set(ontology_table['Synonym'].astype(str).str.strip().str.lower())
                _part = df_study['Samples_Organism part'].astype(str).str.strip().str.lower()
                _m = _part.map(lambda v: None if v in _anatomy or v in ('plasma', 'cell culture') else _envm_vocab.get(v))
                _ok = (_is_missing(df_study['NCBITaxonomy']) | df_study['NCBITaxonomy'].astype(str).str.startswith('2582415|')) & df_study['_qc_label'].isna()
                if 'ENVOEnvironmentMaterial' not in df_study.columns:
                    df_study['ENVOEnvironmentMaterial'] = 'missing value'
                _need = _is_missing(df_study['ENVOEnvironmentMaterial']) & _m.notna() & _ok
                df_study.loc[_need, 'ENVOEnvironmentMaterial'] = _m[_need]

            #add HealthStatus (value-matched)
            #######
            _hs_map = {'healthy': 'healthy', 'health': 'healthy', 'unhealthy': 'unhealthy (NOS)',
                       'diseased': 'unhealthy (NOS)', 'chronic illness': 'chronic illness',
                       'acute illness': 'acute illness'}
            _hs_cols = [c for c in df_study.columns
                        if str(c).startswith('Samples_') and re.search(r'health\s*status|health state', str(c), re.I)]
            if _hs_cols:
                if 'HealthStatus' not in df_study.columns:
                    df_study['HealthStatus'] = 'missing value'
                for _c in _hs_cols:
                    _m = df_study[_c].astype(str).str.strip().str.lower().map(_hs_map)
                    _need = df_study['HealthStatus'].astype(str).str.strip().str.lower().isin(['', 'nan', 'none', 'missing value']) & _m.notna()
                    df_study.loc[_need, 'HealthStatus'] = _m[_need]

            #add DOIDCommonName / HealthStatus from disease columns
            #######
            # Values (or their DOID accessions) that exactly match a DOID label or exact
            # synonym; 'healthy control'-type values mark healthy subjects. Columns about
            # someone other than the subject are skipped.
            _disease_cols = [c for c in df_study.columns
                             if str(c).startswith('Samples_')
                             and re.search(r'disease|patholog|diagnos', str(c), re.I)
                             and not re.search(r'severity|stag|grade|score|duration|onset|family|maternal|mother|father|donor|Term Source|Term Accession|Unit', str(c), re.I)]
            if _disease_cols:
                _doid = kwargs.get('doid_table')
                _doid_allowed = set(allowedTerm_dict['DOIDCommonName']['allowed_values'])
                _doid_labels = {str(v).strip().lower(): str(v) for v in _doid_allowed}
                _doid_by_id, _doid_by_synonym = {}, {}
                if _doid is not None and len(_doid):
                    _doid_by_id = dict(zip(_doid['DOIDOntologyIndex'].astype(str), _doid['Label'].astype(str)))
                    # exact synonyms without the ambiguous abbreviations ('MDS', 'IBD'), plus the
                    # MetaboLights disease aliases (as for the study-level disease step)
                    _doid_by_synonym = _doid_synonyms(_doid_index(_doid, allowedTerm_dict, kwargs.get('doid_parents'), kwargs.get('disease_aliases')))
                def _to_doid(value):
                    key = str(value).strip().lower()
                    if not key or key in _MISSING_VALUES or key == 'disease':
                        return None
                    return _doid_labels.get(key) or (_doid_by_synonym.get(key) if _doid_by_synonym.get(key) in _doid_allowed else None)
                if 'DOIDCommonName' not in df_study.columns:
                    df_study['DOIDCommonName'] = 'missing value'
                if 'HealthStatus' not in df_study.columns:
                    df_study['HealthStatus'] = 'missing value'
                for _c in _disease_cols:
                    _m = df_study[_c].map(_to_doid)
                    _acc_col = _accession_column(df_study, _c)
                    if _acc_col and _doid_by_id:
                        _from_acc = _accession_ids(df_study[_acc_col], 'DOID').map(_doid_by_id)
                        _from_acc = _from_acc.where(_from_acc.isin(_doid_allowed))
                        _m = _from_acc.where(_from_acc.notna(), _m)
                    _need = _is_missing(df_study['DOIDCommonName']) & _m.notna()
                    df_study.loc[_need, 'DOIDCommonName'] = _m[_need]
                    _healthy = df_study[_c].astype(str).str.strip().str.lower().str.fullmatch(
                        r'(?:healthy|non-?diseased?)(?: (?:control|volunteer|subject|donor|individual|participant)s?)?|healthy controls?|control \(healthy\)')
                    _need = _is_missing(df_study['HealthStatus']) & _healthy
                    df_study.loc[_need, 'HealthStatus'] = 'healthy'

            #add DepthorAltitudeMeters
            #######
            if 'Samples_Nominal depth' in df_study.columns:
                def get_depth(val):
                    val = str(val).strip()
                    if val and val not in ('nan', ''):
                        try:
                            return str(-abs(float(val)))
                        except ValueError:
                            pass
                    return 'missing value'
                df_study['DepthorAltitudeMeters'] = df_study['Samples_Nominal depth'].apply(get_depth)

            #add Country
            #######
            # Country was never mapped from MetaboLights. Read geographic-origin
            # columns and match their values against the Country vocabulary
            # (case-insensitive). Sample-level; fills only where currently missing and
            # the value is a recognized country (e.g. "Spain", "Japan").
            _country_cols = ['Samples_Geographical origin', 'Samples_Geographic origin',
                             'Samples_Geographical location', 'Samples_Geographic location',
                             'Samples_Country', 'Samples_country', 'Samples_Country of origin',
                             'Samples_Rice orign', 'Samples_Nationality', 'Samples_Birth country']
            _cc_present = [c for c in _country_cols if c in df_study.columns]
            if _cc_present:
                _cvocab = {str(v).strip().lower(): str(v).strip()
                           for v in allowedTerm_dict.get('Country', {}).get('allowed_values', []) if str(v).strip()}
                if 'Country' not in df_study.columns:
                    df_study['Country'] = 'missing value'
                for _c in _cc_present:
                    _mapped = df_study[_c].astype(str).str.strip().str.lower().map(_cvocab)
                    _need = df_study['Country'].astype(str).str.strip().str.lower().isin(['', 'nan', 'none', 'missing value']) & _mapped.notna()
                    df_study.loc[_need, 'Country'] = _mapped[_need]

            #add SmokingStatus
            #######
            def _toks(col):
                return set(re.split(r'[^a-z0-9]+', str(col).lower()))
            _smoke_cols = [c for c in df_study.columns if str(c).startswith('Samples_') and _toks(c) & {'smoking', 'smoke', 'smoker', 'smokingstatus', 'tobacco'}]
            if _smoke_cols:
                if 'SmokingStatus' not in df_study.columns:
                    df_study['SmokingStatus'] = 'missing value'
                for _c in _smoke_cols:
                    _m = df_study[_c].apply(convert_smoking, key=_c)
                    _need = df_study['SmokingStatus'].astype(str).str.strip().str.lower().isin(['', 'nan', 'none', 'missing value']) & (_m != 'missing value')
                    df_study.loc[_need, 'SmokingStatus'] = _m[_need]

            #add BodyMassIndex
            #######
            _bmi_cols = [c for c in df_study.columns if str(c).startswith('Samples_') and ('bmi' in _toks(c) or {'body', 'mass', 'index'} <= _toks(c))]
            if _bmi_cols:
                if 'BodyMassIndex' not in df_study.columns:
                    df_study['BodyMassIndex'] = 'missing value'
                for _c in _bmi_cols:
                    _m = df_study[_c].apply(bmi_to_numeric).astype(str)
                    _need = df_study['BodyMassIndex'].astype(str).str.strip().str.lower().isin(['', 'nan', 'none', 'missing value']) & (_m != 'missing value')
                    df_study.loc[_need, 'BodyMassIndex'] = _m[_need]

            #add Diet
            #######
            _diet_cols = [c for c in df_study.columns if str(c).startswith('Samples_') and ('diet' in _toks(c) or {'feeding', 'regime'} <= _toks(c) or {'dietary', 'group'} <= _toks(c))]
            if _diet_cols:
                if 'Diet' not in df_study.columns:
                    df_study['Diet'] = 'missing value'
                for _c in _diet_cols:
                    _m = df_study[_c].apply(convert_diet)
                    _need = df_study['Diet'].astype(str).str.strip().str.lower().isin(['', 'nan', 'none', 'missing value']) & (_m != 'missing value')
                    df_study.loc[_need, 'Diet'] = _m[_need]

            if '_qc_label' not in df_study.columns:
                df_study['_qc_label'] = None
            _real_sample = df_study['_qc_label'].isna()

            def _fill(column, values):
                """Set column from values where it is missing and values is not None."""
                if column not in df_study.columns:
                    df_study[column] = 'missing value'
                _need = _is_missing(df_study[column]) & values.notna()
                df_study.loc[_need, column] = values[_need]

            #add SampleExtractionMethod / InternalStandardsUsed (study protocols)
            #######
            _extraction_text = ' '.join(d for n, d in protocols if re.match(r'\s*extraction', n, re.I))
            _ext_term = _extraction_method(_extraction_text, set(allowedTerm_dict['SampleExtractionMethod']['allowed_values']))[0]
            if _ext_term:
                _fill('SampleExtractionMethod', pd.Series(_ext_term, index=df_study.index))
            _is_term = _internal_standards([d for _, d in protocols], allowedTerm_dict['InternalStandardsUsed']['allowed_values'])
            if _is_term:
                _fill('InternalStandardsUsed', pd.Series(_is_term, index=df_study.index))

            #add SubjectIdentifierAsRecorded
            #######
            # Columns that name the sampled individual; a column with fewer than 3 distinct
            # values is a constant or a group, not an identifier. QC/blank/pool/group entries
            # ('Quality Control', 'Healthy Control') are not subjects.
            for _c in [c for c in df_study.columns if _SUBJECT_COLUMNS.fullmatch(str(c))]:
                _v = df_study[_c].astype(str).str.strip()
                _ok = ~_is_missing(df_study[_c]) & ~_v.str.contains(_NOT_A_SUBJECT) & (_v.str.len() <= 60) & _real_sample
                if _v[_ok].nunique() >= 3:
                    _fill('SubjectIdentifierAsRecorded', _v.where(_ok))
                    break

            #add SampleCollectionDateandTime
            #######
            for _c in [c for c in df_study.columns if _DATE_COLUMNS.fullmatch(str(c))]:
                _fill('SampleCollectionDateandTime', _collection_dates(df_study[_c]))

            #add Country from columns holding the collection country as a code or short name
            #######
            _country_vocab = {str(v).strip().lower(): str(v).strip() for v in allowedTerm_dict['Country']['allowed_values']}
            _country_aliases = kwargs.get('country_aliases') or {}
            for _c in [c for c in df_study.columns if _COUNTRY_CODE_COLUMNS.fullmatch(str(c))]:
                _fill('Country', df_study[_c].map(lambda v: _country_term(v, _country_vocab, _country_aliases)))

            #add AgeInYears from other columns holding the subject's age with a unit
            #######
            for _c in [c for c in df_study.columns if _OTHER_AGE_COLUMNS.fullmatch(str(c))]:
                _uc = find_column_after_target_column(df_study, _c, 'Samples_Unit')
                _units = df_study[_uc] if _uc else pd.Series('', index=df_study.index)
                _age = pd.Series([_age_years(v, u, _c) for v, u in zip(df_study[_c], _units)], index=df_study.index, dtype=object)
                _age = _age.where(_real_sample)
                if _age.notna().any():
                    if 'AgeInYears' not in df_study.columns:
                        df_study['AgeInYears'] = 'missing value'
                    if 'LifeStage' not in df_study.columns:
                        df_study['LifeStage'] = ''
                    _need = _is_missing(df_study['AgeInYears']) & _age.notna()
                    df_study.loc[_need, 'AgeInYears'] = _age[_need].astype(str)
                    df_study.loc[_need, 'LifeStage'] = _age[_need].apply(age_category)

            #add DepthorAltitudeMeters from depth/elevation columns with a length unit
            #######
            for _c in [c for c in df_study.columns if re.fullmatch(r'Samples_(?:depth|sediment depth|sampling depth|water depth|elevation|altitude)(?:_\d+)?', str(c), re.I)]:
                _uc = find_column_after_target_column(df_study, _c, 'Samples_Unit')
                if not _uc:
                    continue
                _m = pd.Series([_meters(v, u) for v, u in zip(df_study[_c], df_study[_uc])], index=df_study.index, dtype=object)
                _m = _m.map(lambda x: x if x is None else (abs(x) if re.search('elevation|altitude', _c, re.I) else -abs(x)) + 0.0)
                _fill('DepthorAltitudeMeters', _m.map(lambda x: x if x is None else str(round(x, 4))).where(_real_sample))

            #add BodyMassIndex from per-sample weight (kg) and height (cm/m) of human subjects
            #######
            _w_cols = [c for c in df_study.columns if re.fullmatch(r'Samples_(?:weight|body weight)(?:_\d+)?', str(c), re.I)]
            _h_cols = [c for c in df_study.columns if re.fullmatch(r'Samples_(?:height|body height)(?:_\d+)?', str(c), re.I)]
            if _w_cols and _h_cols and 'NCBITaxonomy' in df_study.columns:
                _wu = find_column_after_target_column(df_study, _w_cols[0], 'Samples_Unit')
                _hu = find_column_after_target_column(df_study, _h_cols[0], 'Samples_Unit')
                if _wu and _hu:
                    _kg = pd.to_numeric(df_study[_w_cols[0]], errors='coerce').where(df_study[_wu].astype(str).str.strip().str.lower().isin(['kilogram', 'kg']))
                    _hm = pd.to_numeric(df_study[_h_cols[0]], errors='coerce')
                    _hunit = df_study[_hu].astype(str).str.strip().str.lower()
                    _hm = (_hm / 100).where(_hunit.isin(['centimeter', 'cm']), _hm.where(_hunit.isin(['meter', 'm'])))
                    # weight and height entered as the same number are a data entry error
                    _same = pd.to_numeric(df_study[_w_cols[0]], errors='coerce') == pd.to_numeric(df_study[_h_cols[0]], errors='coerce')
                    _bmi = (_kg / _hm ** 2).where((df_study['NCBITaxonomy'] == '9606|Homo sapiens') & ~_same & _real_sample)
                    _fill('BodyMassIndex', _bmi.round(1).map(lambda x: None if pd.isna(x) else str(x)))

            #add DOIDCommonName from condition/phenotype/cohort columns (exact DOID label/synonym)
            #######
            _doid_cols = [c for c in df_study.columns if re.fullmatch(r'Samples_(?:condition|phenotype|cohort|patient group|infection|medical condition'
                                                                    r'|clinical condition|health condition|clinical phenotype|disease state)(?:_\d+)?', str(c), re.I)]
            if _doid_cols:
                _doid = kwargs.get('doid_table')
                _doid_allowed = set(allowedTerm_dict['DOIDCommonName']['allowed_values'])
                _doid_labels = {str(v).strip().lower(): str(v) for v in _doid_allowed}
                _doid_by_synonym = {}
                if _doid is not None and len(_doid):
                    _doid_by_synonym = _doid_synonyms(_doid_index(_doid, allowedTerm_dict, kwargs.get('doid_parents'), kwargs.get('disease_aliases')))
                def _to_doid_exact(value):
                    key = str(value).strip().lower()
                    if len(key) <= 3 or key in _MISSING_VALUES or key == 'disease':
                        return None
                    return _doid_labels.get(key) or _doid_by_synonym.get(key)
                for _c in _doid_cols:
                    _fill('DOIDCommonName', df_study[_c].map(_to_doid_exact).where(_real_sample))

            #add DOIDCommonName / HealthStatus from the study-level disease and the sample groups
            #######
            _doid_tab = kwargs.get('doid_table')
            if (_doid_tab is not None and len(_doid_tab) and allowedTerm_dict['DOIDCommonName']['allowed_values']
                    and 'NCBITaxonomy' in df_study.columns and 'SampleType' in df_study.columns):
                _idx = _doid_index(_doid_tab, allowedTerm_dict, kwargs.get('doid_parents'), kwargs.get('disease_aliases'))
                _taxid = df_study['NCBITaxonomy'].astype(str).str.split('|').str[0]
                # cell types (CL) other than blood and germ cells are cultured cells
                _cl = pd.Series(False, index=df_study.index)
                if 'UBERONOntologyIndex' in df_study.columns and 'UBERONBodyPartName' in df_study.columns:
                    _cl = (df_study['UBERONOntologyIndex'].astype(str).str.startswith('CL:')
                           & ~df_study['UBERONBodyPartName'].astype(str).str.contains(_D_PRIMARY_CELLS))
                _eligible = (_real_sample & (df_study['SampleType'] == 'animal') & _taxid.isin(_mammal_taxids(NCBIRankDivision_table))
                             & ~_disease_not_subject(df_study) & ~_cl)
                _human = bool((_taxid[_eligible] == '9606').mean() > 0.5) if _eligible.any() else False
                _content = study_details['content']
                # animals given a human microbiota: disease labels describe the donors
                if not _human and re.search(_D_MICROBIOTA_RECIPIENT, f"{_content.get('title') or ''} {_content.get('description') or ''}"):
                    _eligible = pd.Series(False, index=df_study.index)
                if _eligible.any():
                    _kind = {c: _sample_kind.get(re.sub(r'_\d+$', '', c), _sample_kind.get(c)) for c in df_study.columns
                             if str(c).startswith('Samples_')}
                    _sub = df_study['SubjectIdentifierAsRecorded'].reset_index(drop=True) if 'SubjectIdentifierAsRecorded' in df_study.columns else None
                    _dv, _hv, _ = _study_disease_labels(_idx, df_study[list(_kind)].reset_index(drop=True), _kind,
                                                        [d.get('description') for d in _content.get('descriptors') or []],
                                                        [f.get('name') for f in _content.get('factors') or []],
                                                        _content.get('title'), _content.get('description') or _content.get('studyDescription'),
                                                        _eligible.reset_index(drop=True), _human, _sub)
                    _dv.index = df_study.index
                    _hv.index = df_study.index
                    for _col in ('DOIDCommonName', 'HealthStatus'):
                        if _col not in df_study.columns:
                            df_study[_col] = 'missing value'
                    # never against a recorded 'healthy' or a recorded disease
                    _case = (_hv == 'unhealthy (NOS)') & (df_study['HealthStatus'].astype(str).str.strip() != 'healthy')
                    _healthy = (_hv == 'healthy') & _is_missing(df_study['DOIDCommonName'])
                    # a file listed for several samples gets the labels only if they agree for all of them
                    _final = pd.DataFrame({'f': df_study['filename'].astype(str),
                                           'd': df_study['DOIDCommonName'].where(~_is_missing(df_study['DOIDCommonName']), _dv.where(_case)).fillna('missing value').astype(str),
                                           'h': df_study['HealthStatus'].where(~_is_missing(df_study['HealthStatus']), _hv.where(_case | _healthy)).fillna('missing value').astype(str)})
                    _agree = (_final.groupby('f')['d'].transform('nunique') == 1) & (_final.groupby('f')['h'].transform('nunique') == 1)
                    _fill('DOIDCommonName', _dv.where(_case & _agree))
                    _fill('HealthStatus', _hv.where((_case | _healthy) & _agree))

            # rows that are samples of an individual: no QC, cultured cells or cell types (CL)
            # other than blood and germ cells
            _not_cells = _real_sample & ~_disease_not_subject(df_study)
            if 'UBERONOntologyIndex' in df_study.columns and 'UBERONBodyPartName' in df_study.columns:
                _not_cells &= ~(df_study['UBERONOntologyIndex'].astype(str).str.startswith('CL:')
                                & ~df_study['UBERONBodyPartName'].astype(str).str.contains(_D_PRIMARY_CELLS))
            _kind = {c: _sample_kind.get(re.sub(r'_\d+$', '', c), _sample_kind.get(c)) for c in df_study.columns
                     if str(c).startswith('Samples_')}
            _content = study_details['content']
            _text = _study_text(_content, protocols)
            _tokens = _value_tokens(df_study, [c for c, k in _kind.items() if k in ('Characteristics', 'Factor Value')])

            #add HPOCommonName / HealthStatus from the study-level health condition and the sample groups
            #######
            _hpo_tab = kwargs.get('hpo_table')
            if (_hpo_tab is not None and len(_hpo_tab) and allowedTerm_dict.get('HPOCommonName', {}).get('allowed_values')
                    and 'NCBITaxonomy' in df_study.columns and 'SampleType' in df_study.columns):
                _doid_tab = kwargs.get('doid_table')
                _didx = (_doid_index(_doid_tab, allowedTerm_dict, kwargs.get('doid_parents'), kwargs.get('disease_aliases'))
                         if _doid_tab is not None and len(_doid_tab) else None)
                _hidx = _hpo_index(_hpo_tab, allowedTerm_dict, kwargs.get('hpo_aliases'), _didx)
                # abbreviations the study defines ('acute kidney injury (AKI)'), and standard
                # disease-model procedures the study text names ('cecal ligation and puncture')
                _habbr = _defined_abbreviations(_text, _hidx.term, _hidx.abbr_syn, _tokens)
                _hdesc = [d.get('description') for d in _content.get('descriptors') or []]
                _tl = _dlabel(_text)
                for _proc, _pab, _pterm in kwargs.get('procedure_aliases') or []:
                    if _hidx.term(_pterm) == _pterm and re.search(r'(?<![a-z0-9])' + re.escape(_proc) + r'(?![a-z0-9])', _tl):
                        _habbr[_proc] = _dlabel(_pterm)
                        if _pab:
                            _habbr.setdefault(_pab, _dlabel(_pterm))
                        _hdesc.append(_pterm)
                _hextra = {k: _hidx.term(f) for k, f in _habbr.items()}
                _hextra = {k: t for k, t in _hextra.items() if t and not _hidx.foreign(t)}
                _hidx = _StudyTerms(_hidx, _habbr) if _habbr else _hidx
                _taxid = df_study['NCBITaxonomy'].astype(str).str.split('|').str[0]
                _eligible = _not_cells & (df_study['SampleType'] == 'animal') & _taxid.isin(_mammal_taxids(NCBIRankDivision_table))
                _human = bool((_taxid[_eligible] == '9606').mean() > 0.5) if _eligible.any() else False
                # animals given a human microbiota: condition labels describe the donors
                if not _human and re.search(_D_MICROBIOTA_RECIPIENT, f"{_content.get('title') or ''} {_content.get('description') or ''}"):
                    _eligible = pd.Series(False, index=df_study.index)
                if _eligible.any():
                    _sub = df_study['SubjectIdentifierAsRecorded'].reset_index(drop=True) if 'SubjectIdentifierAsRecorded' in df_study.columns else None
                    _pv, _hv, _hinfo = _study_disease_labels(_hidx, df_study[list(_kind)].reset_index(drop=True), _kind,
                                                        _hdesc,
                                                        [f.get('name') for f in _content.get('factors') or []],
                                                        _content.get('title'), _content.get('description') or _content.get('studyDescription'),
                                                        _eligible.reset_index(drop=True), _human, _sub, plain_control_healthy=False,
                                                        extra_abbr=_hextra)
                    _pv.index = df_study.index
                    _hv.index = df_study.index
                    for _col in ('HPOCommonName', 'HealthStatus'):
                        if _col not in df_study.columns:
                            df_study[_col] = 'missing value'
                    _case = (_hv == 'unhealthy (NOS)') & (df_study['HealthStatus'].astype(str).str.strip() != 'healthy') & _pv.notna()
                    # a patient / severity / model label of a study with a disease stands for the disease
                    _named = _hinfo['named'].set_axis(df_study.index) if 'named' in _hinfo else pd.Series(False, index=df_study.index)
                    _case &= _named | _is_missing(df_study['DOIDCommonName'] if 'DOIDCommonName' in df_study.columns else pd.Series('', index=df_study.index))
                    _healthy = (_hv == 'healthy') & _is_missing(df_study['HPOCommonName'])
                    # a file listed for several samples gets the labels only if they agree for all of them
                    _final = pd.DataFrame({'f': df_study['filename'].astype(str),
                                           'p': df_study['HPOCommonName'].where(~_is_missing(df_study['HPOCommonName']), _pv.where(_case)).fillna('missing value').astype(str),
                                           'h': df_study['HealthStatus'].where(~_is_missing(df_study['HealthStatus']), _hv.where(_case | _healthy)).fillna('missing value').astype(str)})
                    _agree = (_final.groupby('f')['p'].transform('nunique') == 1) & (_final.groupby('f')['h'].transform('nunique') == 1)
                    _fill('HPOCommonName', _pv.where(_case & _agree))
                    _fill('HealthStatus', _hv.where((_case | _healthy) & _agree))

            #add Exposure from the treatment / exposure / stress columns of the exposed groups
            #######
            _exp_tab = kwargs.get('exposure_table')
            if (_exp_tab is not None and len(_exp_tab) and allowedTerm_dict.get('Exposure', {}).get('allowed_values')
                    and 'SampleType' in df_study.columns):
                _eidx = _exposure_index(_exp_tab, allowedTerm_dict, kwargs.get('exposure_aliases'))
                _plant = df_study['SampleType'] == 'plant'
                # with plant samples (or a PECO design descriptor) in the study, the treatments are given
                # to the plants (insects feeding on inoculated plants are not inoculated themselves)
                _plant_study = bool((_plant & _real_sample).any()) or any(
                    str(d.get('description') or '').strip().upper().startswith('PECO:') for d in _content.get('descriptors') or [])
                _e_eligible = _not_cells & (_plant | ((df_study['SampleType'] == 'animal') & (not _plant_study)))
                # cells grown from a tissue ('bone marrow-derived macrophage', 'iPSC derived') are cultured
                for _c in ('Samples_Organism part', 'UBERONBodyPartName'):
                    if _c in df_study.columns:
                        _e_eligible &= ~df_study[_c].astype(str).str.contains(r'derived|\bcultured?\b', case=False)
                # isolated cell types (CL) other than blood cells are treated in culture ('bone marrow stromal cell')
                if 'UBERONOntologyIndex' in df_study.columns and 'UBERONBodyPartName' in df_study.columns:
                    _e_eligible &= ~(df_study['UBERONOntologyIndex'].astype(str).str.startswith('CL:')
                                     & ~df_study['UBERONBodyPartName'].astype(str).str.contains(_E_BLOOD_CELLS))
                if _e_eligible.any():
                    _eabbr = _defined_abbreviations(_text, lambda f: _eidx.term(f, True), _eidx.abbr_syn, _tokens)
                    _ev = _sample_exposures(_eidx, df_study[list(_kind)].reset_index(drop=True), _kind,
                                            _e_eligible.reset_index(drop=True), _plant.reset_index(drop=True), _eabbr)
                    _ev.index = df_study.index
                    if 'Exposure' not in df_study.columns:
                        df_study['Exposure'] = 'missing value'
                    _final = pd.DataFrame({'f': df_study['filename'].astype(str),
                                           'e': df_study['Exposure'].where(~_is_missing(df_study['Exposure']), _ev).fillna('missing value').astype(str)})
                    _fill('Exposure', _ev.where(_final.groupby('f')['e'].transform('nunique') == 1))

            #add TimePointHoursRelativeToEvent / FastingStatus / SampleGroupFreeText
            #######
            # Workbench rules on the Characteristics / Factor Value columns. Time point and
            # fasting only for samples of the study's subjects (no QC / blanks), fasting only
            # for whole animals and humans (no cultured or isolated cells, plants, microbes).
            _value_cols = [c for c in df_study.columns if _sample_kind_pos.get(c) in ('Characteristics', 'Factor Value')]
            try:
                _tp_text = re.sub(r'<[^>]+>', ' ', ' \n '.join(str(p) for p in (_content.get('title'), _content.get('description') or _content.get('studyDescription')) if p))
                _tp = ml_sample_timepoints(df_study, _value_cols, _tp_text)
                if _tp.notna().any():
                    _fill('TimePointHoursRelativeToEvent', _tp.where(_real_sample))
            except Exception as e:
                print(f"Could not derive time points: {type(e).__name__}: {e}")
            try:
                _animal = _not_cells & (df_study['SampleType'] == 'animal') if 'SampleType' in df_study.columns \
                    else pd.Series(False, index=df_study.index)
                for _c in ('Samples_Organism part', 'UBERONBodyPartName'):
                    if _c in df_study.columns:
                        _animal &= ~df_study[_c].astype(str).str.contains(r'derived|\bcultured?\b', case=False)
                _fs = ml_sample_fasting(df_study, _value_cols)
                if _fs.notna().any():
                    _fill('FastingStatus', _fs.where(_animal))
                elif _animal.any():
                    # one fasting statement in the protocols, only for studies without
                    # per-sample values and without fasting / feeding design columns
                    _design = any(re.search(r'fast|feed|\bfed\b|prandial', _ml_key(c), re.I)
                                  or df_study[c].astype(str).str.contains(r'fast|\bfed\b|prandial', case=False).any()
                                  for c in _value_cols)
                    _st = None if _design else protocol_fasting_status(protocols)
                    if _st:
                        _status, _spec, _stmts = _st
                        _use = _animal.copy()
                        if _spec is not None:
                            # a statement about some specimens ('fasting blood samples') only for these
                            _fam = df_study['UBERONBodyPartName'].map(_specimen_family) if 'UBERONBodyPartName' in df_study.columns \
                                else pd.Series(None, index=df_study.index, dtype=object)
                            _use &= _fam.isin(_spec)
                        if _use.any():
                            print(f"[FastingStatus] {study_id}: study level '{_status}' for {int(_use.sum())} rows from: {' | '.join(_stmts)[:500]}")
                            _fill('FastingStatus', pd.Series(_status, index=df_study.index).where(_use))
            except Exception as e:
                print(f"Could not derive fasting states: {type(e).__name__}: {e}")
            try:
                _groups = ml_sample_group_texts(df_study, [c for c in df_study.columns if _sample_kind_pos.get(c) == 'Factor Value'],
                                                [f.get('name') for f in _content.get('factors') or []], study_id=study_id)
                if _groups.notna().any():
                    df_study['SampleGroupFreeText'] = _groups.fillna('missing value')
            except Exception as e:
                print(f"Could not derive sample groups: {type(e).__name__}: {e}")

            #refine NCBITaxonomy from the Variant/Strain column
            #######
            # Only an NCBI name that starts with the recorded species name ('Saccharomyces
            # cerevisiae BY4741' for 'Saccharomyces cerevisiae'), i.e. a lower rank of it
            _var_cols = [c for c in df_study.columns if re.fullmatch(r'Samples_(?:variant|strain)(?:_\d+)?', str(c), re.I)]
            if _var_cols and 'NCBITaxonomy' in df_study.columns:
                for _c in _var_cols:
                    def _refine(term, variant):
                        if not isinstance(term, str) or '|' not in term:
                            return None
                        species = term.split('|', 1)[1].strip().lower()
                        v = re.sub(r'\s+', ' ', str(variant)).strip()
                        if not v.lower().startswith(species + ' ') or len(v) <= len(species) + 1:
                            return None
                        hit = _ncbi_term_by_name(allowedTerm_dict).get(v.lower())
                        return hit if hit and '|' in hit and hit != term and hit.split('|', 1)[1].strip().lower() == v.lower() else None
                    _cache = {}
                    _ref = pd.Series([_cache.setdefault((t, v), _refine(t, v)) for t, v in zip(df_study['NCBITaxonomy'], df_study[_c])],
                                     index=df_study.index, dtype=object)
                    _use = _ref.notna() & _real_sample
                    df_study.loc[_use, 'NCBITaxonomy'] = _ref[_use]

            #add MassiveID and USIs
            #######
            df_study['MassiveID'] = study_id
            
            ontology_table = ontology_table.drop_duplicates(subset=['Label'])
            df_study = merge_repeated_fileobservations(df_study)
            # HPO and exposure tables (one row per label) so complete_and_fill_REDU_table can generate their indexes
            _index_tables = {}
            for _arg, _tab, _id in (('DOIDOntologyIndex_table', kwargs.get('doid_table'), 'DOIDOntologyIndex'),
                                    ('HPOOntologyIndex_table', kwargs.get('hpo_table'), 'HPOOntologyIndex'),
                                    ('ExposureIndex_table', kwargs.get('exposure_table'), 'ExposureIndex')):
                if _tab is not None and len(_tab):
                    # the exposure table also carries the exposure agent (OrganismExposureAgent / Index)
                    _index_tables[_arg] = _tab.drop_duplicates(subset=['Label'])[['Label', _id] + [c for c in ('OrganismExposureAgentIndex', 'OrganismExposureAgent') if c in _tab.columns]].copy()
            df_study = complete_and_fill_REDU_table(df_study, allowedTerm_dict, UBERONOntologyIndex_table=ontology_table, ENVOEnvironmentBiomeIndex_table=ENVOEnvironmentBiomeIndex_table,
                                                    ENVOEnvironmentMaterialIndex_table=ENVOEnvironmentMaterialIndex_table,NCBIRankDivision_table=NCBIRankDivision_table, add_usi = True, 
                                                    other_allowed_file_extensions = ['.raw', '.cdf', '.wiff', '.d'] + COMPRESSED_MS_FILE_EXTENSIONS, **_index_tables)
            
            df_study = df_study.drop_duplicates() 

            #files assigned to several samples were already collapsed into one row by merge_repeated_fileobservations,
            #with conflicting fields set to 'missing value'. Rows can still share a USI here when different filenames
            #resolve to the same file (e.g. several paths inside one .d folder); drop those as we cannot tell them apart
            df_study['count'] = df_study.groupby('USI')['USI'].transform('size')
            df_study = df_study[df_study['count'] == 1]
            df_study = df_study.drop(columns=['count'])

            # sample groups: factors whose groups the final ReDU columns tell apart are dropped;
            # filled only for at least 2 groups and on average at least 2 samples per group
            if 'SampleGroupFreeText' in df_study.columns:
                try:
                    df_study = finalize_sample_groups(df_study, study_id=study_id)
                except Exception as e:
                    print(f"Could not finalize sample groups: {type(e).__name__}: {e}")
                    df_study['SampleGroupFreeText'] = 'missing value'

            print(f'We are adding {len(df_study)} samples to the REDU table for study {study_id}.')

            return df_study




if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Give an Metabolights study ID and get a REDU table tsv.')
    parser.add_argument("--study_id", type=str, help='An Metabolights study ID such as "MTBLS1015". If "ALL" all studys are requested.', required=True)
    parser.add_argument("--path_to_translation_sheet_csvs", "-csvs", type=str, help="Path to the translation csvs holding translations from MWB to REDU vocabulary", default="none")
    parser.add_argument("--path_to_allowed_term_json", type=str, help="Path to the json with allowed REDU terms")
    parser.add_argument("--path_to_uberon_cl_po_csv", type=str, help="Path to the prepared uberon_cl_po ontology csv")
    parser.add_argument("--path_to_envo_biome_csv", type=str, help="Path to the prepared uberon_cl_po ontology csv")
    parser.add_argument("--path_to_envo_material_csv", type=str, help="Path to the prepared uberon_cl_po ontology csv")
    parser.add_argument("--path_ncbi_rank_division", type=str, help="Path to the path_ncbi_rank_division")
    parser.add_argument("--path_to_doid_parents_csv", type=str, default="none", help="Optional DOID is_a links from prepare_ontologies.py (most specific disease term)")
    parser.add_argument("--path_to_doid_csv", type=str, default="none", help="Optional path to the prepared DOID ontology csv (exact synonyms and DOID accessions for DOIDCommonName)")
    parser.add_argument("--path_to_hpo_csv", type=str, default="none", help="Optional path to the prepared HPO ontology csv (HPOCommonName)")
    parser.add_argument("--path_to_exposure_csv", type=str, default="none", help="Optional path to the prepared ECTO/PECO exposure csv (Exposure)")
            
    args = parser.parse_args()

    if args.study_id == 'ALL':
        public_metabolights_studies = safe_api_request('https://www.ebi.ac.uk:443/metabolights/ws/studies', retries = 5)
        if not public_metabolights_studies or not public_metabolights_studies.get('content'):
            raise RuntimeError("Could not retrieve the list of MetaboLights studies.")
        public_metabolights_studies = public_metabolights_studies['content']
    else:
        public_metabolights_studies = [args.study_id]

    if args.path_to_translation_sheet_csvs == 'none':
        script_dir = os.path.dirname(os.path.realpath(__file__))
        transSheet_dir = os.path.join(script_dir, 'translation_sheets_metabolights')
    else:
        transSheet_dir = args.path_to_translation_sheet_csvs
        script_dir = os.path.dirname(transSheet_dir)

    if args.path_to_allowed_term_json == 'none':
        allowedTermSheet_json = os.path.join(script_dir, 'allowed_terms', 'allowed_terms.json')
    else:
        allowedTermSheet_json = args.path_to_allowed_term_json


    # Read allowed terms json
    with open(allowedTermSheet_json, 'r') as file:
        allowedTerm_dict = json.load(file)

    # Add terms for matching
    allowed_values = allowedTerm_dict["MassSpectrometer"]["allowed_values"]
    allowedTerm_dict["MassSpectrometer"]["allowed_values_matching_0"] = [value.split('|')[0] for value in allowed_values]


    # Read ontology tables
    ontology_table = pd.read_csv(args.path_to_uberon_cl_po_csv)
    ENVOEnvironmentBiomeIndex_table = pd.read_csv(args.path_to_envo_biome_csv)
    ENVOEnvironmentMaterialIndex_table = pd.read_csv(args.path_to_envo_material_csv)

    # Read NCBI rank and division csv
    NCBIRankDivision_table = pd.read_csv(args.path_ncbi_rank_division, index_col = False)
    NCBIRankDivision_table = NCBIRankDivision_table.drop_duplicates(subset=['TaxonID'])

    # Optional tables: DOID ontology, and the body part -> biofluid/tissue translation sheet
    doid_table = pd.read_csv(args.path_to_doid_csv) if args.path_to_doid_csv != 'none' else None
    # Optional HPO (health conditions) and ECTO/PECO (exposures) tables
    hpo_table = pd.read_csv(args.path_to_hpo_csv, dtype=str, keep_default_na=False) if args.path_to_hpo_csv != 'none' else None
    exposure_table = pd.read_csv(args.path_to_exposure_csv, dtype=str, keep_default_na=False) if args.path_to_exposure_csv != 'none' else None
    biofluid_vs_tissue_path = os.path.join(transSheet_dir, 'biofluid_vs_tissue_Metabolights.csv')
    biofluid_vs_tissue = pd.read_csv(biofluid_vs_tissue_path) if os.path.exists(biofluid_vs_tissue_path) else None
    sheets = load_metabolights_sheets(transSheet_dir, allowedTerm_dict)
    # DOID is_a links from prepare_ontologies.py, to pick the most specific disease term
    if args.path_to_doid_parents_csv != 'none' and os.path.exists(args.path_to_doid_parents_csv):
        t = pd.read_csv(args.path_to_doid_parents_csv, dtype=str, keep_default_na=False)
        sheets['doid_parents'] = t.groupby('DOIDOntologyIndex')['ParentDOIDOntologyIndex'].agg(list).to_dict()


    REDU_dataframes = []
    redu_table_single = pd.DataFrame()
    not_public_studies = []
    no_data_studies = []
    failed_studies = []
    for study_id in tqdm(public_metabolights_studies):
        try:
            print(f'Processing study {study_id}...')
            redu_table_single = Metabolights2REDU(study_id, allowedTerm_dict = allowedTerm_dict, ontology_table = ontology_table, ENVOEnvironmentBiomeIndex_table=ENVOEnvironmentBiomeIndex_table,
                                                  ENVOEnvironmentMaterialIndex_table=ENVOEnvironmentMaterialIndex_table, NCBIRankDivision_table=NCBIRankDivision_table,
                                                  doid_table=doid_table, hpo_table=hpo_table, exposure_table=exposure_table,
                                                  biofluid_vs_tissue=biofluid_vs_tissue, **sheets)
        except StudyNotAvailable as e:
            print(f"Study {study_id} is not available: {e}")
            not_public_studies.append(study_id)
            continue
        except Exception as e:
            traceback_info = traceback.format_exc()
            print(f"An error occurred with study_id {study_id}: {e}\nTraceback:\n{traceback_info}")
            failed_studies.append(study_id)
            continue
        if redu_table_single is not None and len(redu_table_single) > 0:
            print(f'Added {len(redu_table_single)} samples.')
            REDU_dataframes.append(redu_table_single)
        else:
            print(f'Added {0} samples.')
            no_data_studies.append(study_id)

    print(f'\nStudies with samples: {len(REDU_dataframes)}, without MS samples: {len(no_data_studies)}, '
          f'not public: {len(not_public_studies)}, failed: {len(failed_studies)} (of {len(public_metabolights_studies)})')
    if failed_studies:
        print('Failed studies: ' + ' '.join(failed_studies))

    if len(REDU_dataframes) > 0:
        redu_tables_all = pd.concat(REDU_dataframes, ignore_index=True)
        redu_tables_all = redu_tables_all.drop_duplicates()
        redu_tables_all.to_csv('Metabolights2REDU_' + args.study_id + '.tsv', sep='\t', index=False, header=True)
        print(f'Output of {len(redu_tables_all)} samples has been saved to Metabolights2REDU_{args.study_id}.tsv!')
    else:
        print('nothing to return!')

