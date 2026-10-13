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




def prefer_extension(group):
    extensions = group['extension'].values
    if 'mzml' in extensions or 'mzxml' in extensions:
        group['keep'] = (group['extension'] == 'mzml') | (group['extension'] == 'mzxml')
    else:
        # Mark only the first row to keep if no preferred extension found
        group['keep'] = [True] + [False] * (len(group) - 1)
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
    return sheets


class StudyNotAvailable(Exception):
    """The study exists but is not public (MetaboLights returns no content for it)."""


def safe_api_request(url, retries=3, expected_codes={200}):

  for _ in range(retries):
    try:
      response = requests.get(url, timeout=(10, 300))
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
        response = requests.get(url, timeout=(10, 120))
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


def Metabolights2REDU(study_id, **kwargs):
    """
    Converts Metabolights study data to a REDU table format.

    Args:
    study_id: The ID of the Metabolights study to convert.

    Returns:
    A DataFrame in the REDU table format with processed and aligned data from the Metabolights study,
    or an empty DataFrame if no applicable data is found.
    """
    study_url = "https://www.ebi.ac.uk:443/metabolights/ws/studies/public/study/" + study_id
    study_details = safe_api_request(study_url)
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

        # Perform the split operation with expansion to get two columns
        split_values = df_study['filename_lower'].str.rsplit('.', n=1, expand=True)

        # List of allowed extensions
        allowed_extensions = [".mzml", ".mzxml", ".cdf", ".raw", ".wiff", ".d"]
        df_study = df_study[df_study['filename_lower'].apply(lambda x: any(x.endswith(ext) for ext in allowed_extensions))]

        print(f'Assay table has {len(df_study)} rows and {len(df_study.columns)} columns.')

        if len(df_study) == 0:
            # e.g. only .zip archives or .xls tables listed, nothing REDU can index
            return pd.DataFrame()

        # Assign the split columns to 'base' and 'extension'
        df_study['base'] = split_values[0]

        # Initially set 'extension' to None to handle rows without an extension
        df_study['extension'] = None

        # Only update 'extension' for rows where a split occurred
        df_study.loc[split_values[1].notna(), 'extension'] = split_values[1]

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
                df_study.loc[:, 'NCBITaxonomy'] = df_study['NCBITaxonomy'].replace(to_replace=r'^.*None.*$', value='missing value', regex=True)

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
                                    if '|' in taxonomy and 'None' not in taxonomy}
                
                print(f'Found {len(processed_taxonomy)} unique taxonomy IDs.')


                df_study.loc[:,'SampleType'] = df_study['NCBITaxonomy'].map(lambda x: processed_taxonomy.get(x.split('|')[0], [pd.NA, pd.NA])[0] 
                                                                    if '|' in x and 'None' not in x else pd.NA)
                df_study.loc[:,'SampleTypeSub1'] = df_study['NCBITaxonomy'].map(lambda x: processed_taxonomy.get(x.split('|')[0], [pd.NA, pd.NA])[1]
                                                                        if '|' in x and 'None' not in x else pd.NA)

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
                    _syn = _doid.assign(_s=_doid['Synonym'].astype(str).str.strip().str.lower())
                    _syn = _syn.groupby('_s')['Label'].agg(lambda x: set(x))
                    _doid_by_synonym = {k: next(iter(v)) for k, v in _syn.items() if len(v) == 1}
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
                    _syn = _doid.assign(_s=_doid['Synonym'].astype(str).str.strip().str.lower()).groupby('_s')['Label'].agg(lambda x: set(x))
                    _doid_by_synonym = {k: next(iter(v)) for k, v in _syn.items() if len(v) == 1 and next(iter(v)) in _doid_allowed}
                def _to_doid_exact(value):
                    key = str(value).strip().lower()
                    if len(key) <= 3 or key in _MISSING_VALUES or key == 'disease':
                        return None
                    return _doid_labels.get(key) or _doid_by_synonym.get(key)
                for _c in _doid_cols:
                    _fill('DOIDCommonName', df_study[_c].map(_to_doid_exact).where(_real_sample))

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
            df_study = complete_and_fill_REDU_table(df_study, allowedTerm_dict, UBERONOntologyIndex_table=ontology_table, ENVOEnvironmentBiomeIndex_table=ENVOEnvironmentBiomeIndex_table,
                                                    ENVOEnvironmentMaterialIndex_table=ENVOEnvironmentMaterialIndex_table,NCBIRankDivision_table=NCBIRankDivision_table, add_usi = True, 
                                                    other_allowed_file_extensions = ['.raw', '.cdf', '.wiff', '.d'])
            
            df_study = df_study.drop_duplicates() 

            #files assigned to several samples were already collapsed into one row by merge_repeated_fileobservations,
            #with conflicting fields set to 'missing value'. Rows can still share a USI here when different filenames
            #resolve to the same file (e.g. several paths inside one .d folder); drop those as we cannot tell them apart
            df_study['count'] = df_study.groupby('USI')['USI'].transform('size')
            df_study = df_study[df_study['count'] == 1]
            df_study = df_study.drop(columns=['count'])

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
    parser.add_argument("--path_to_doid_csv", type=str, default="none", help="Optional path to the prepared DOID ontology csv (exact synonyms and DOID accessions for DOIDCommonName)")
            
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
    biofluid_vs_tissue_path = os.path.join(transSheet_dir, 'biofluid_vs_tissue_Metabolights.csv')
    biofluid_vs_tissue = pd.read_csv(biofluid_vs_tissue_path) if os.path.exists(biofluid_vs_tissue_path) else None
    sheets = load_metabolights_sheets(transSheet_dir, allowedTerm_dict)


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
                                                  doid_table=doid_table, biofluid_vs_tissue=biofluid_vs_tissue, **sheets)
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

