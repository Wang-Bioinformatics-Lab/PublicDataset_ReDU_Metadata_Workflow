"""Blanks and pooled QCs recognised from the file name.

A file is a blank or a pooled QC when its name (basename without extension) contains one
of these words as a whole token. Tokens are the parts of the name between '_', '-', '.',
spaces and other non-alphanumeric characters; a token may carry a number ('Blank11',
'QC03', 'pool2'), but a word inside a longer word does not count ('QCT', 'blanket',
'spool', 'Poole').

    blank:  blank, blk (also written together with solvent/method/process(ing)/
            procedural/extraction/instrument/injection, e.g. 'SolventBlank')
    QC:     qc, pool, pooled, poolqc, pooledqc, qcpool

Type:
    blank token + extraction/extract/ext/procedural/procedure/method/process(ing) token
        (or 'MethodBlank', 'ExtractionBlank', ...)      -> blank_extraction
    blank token + QC token                              -> blank_QC
    other blank tokens                                  -> blank_analysis
    QC token                                            -> pool_QC

'pool' is not a pooled QC where the name says it is a water body ('tide pool', 'rock
pool', 'swimming pool', 'pool water'): see _POOL_PLACE_TOKENS.

A token that mixes letters and digits is also read as its letter and digit runs, so a word
followed by a number and more letters counts ('qc1list2' -> qc 1 list 2, as 'qc1-list2';
'S3QC' -> s 3 qc); a word inside a longer word still does not ('QCT', 'blanket').
"""
import os
import re

import pandas as pd

_SPLIT = re.compile(r'[^A-Za-z0-9]+')
_NUM = r'\d*'
_BLANK_PREFIXES = r'(?:solvent|method|process|processing|procedural|procedure|extraction|extract|instrument|injection|inj|lc|system)'
_EXTRACTION_PREFIXES = r'(?:method|process|processing|procedural|procedure|extraction|extract)'
_BLANK_TOKEN = re.compile(rf'(?i)^(?:{_NUM}(?:blank|blk|{_BLANK_PREFIXES}blank){_NUM}|blank{_NUM}[a-z]?)$')
_EXTRACTION_BLANK_TOKEN = re.compile(rf'(?i)^{_NUM}{_EXTRACTION_PREFIXES}blank{_NUM}$')
_EXTRACTION_WORD = re.compile(r'(?i)^(?:extraction|extract|ext|extr|procedural|procedure|method|process|processing)$')
_QC_TOKEN = re.compile(rf'(?i)^{_NUM}(?:qc|pool|pooled|poolqc|pooledqc|qcpool){_NUM}$')
_POOL_TOKEN = re.compile(rf'(?i)^{_NUM}(?:pool|pooled){_NUM}$')
# 'pool' as a body of water / place, not a pooled sample
_POOL_PLACE_TOKENS = {'tide', 'tidal', 'rock', 'rockpool', 'swimming', 'swim', 'water', 'pond', 'spa', 'beach',
                      'river', 'lake', 'stream', 'puddle', 'marsh', 'seawater'}
# MS file extensions that are removed before the name is split (also double ones: x.mzML.gz)
_EXTENSIONS = re.compile(r'(?i)(\.(mzml|mzxml|mzdata|mgf|raw|cdf|wiff2?|d|gz|zip|lzma|bz2|7z|xz|rar))+$')

BLANK_QC_TYPES = ('blank_analysis', 'blank_extraction', 'blank_QC', 'pool_QC')


def name_tokens(filename):
    base = os.path.basename(str(filename).rstrip('/'))
    base = _EXTENSIONS.sub('', base)
    return [t for t in _SPLIT.split(base) if t]


_RUNS = re.compile(r'[A-Za-z]+|[0-9]+')


def _tokens_and_runs(tokens):
    """the tokens and, for tokens of letters and digits, their letter and digit runs ('qc1list2' ->
    'qc1list2', 'qc', '1', 'list', '2')"""
    out = list(tokens)
    for t in tokens:
        runs = _RUNS.findall(t)
        if len(runs) > 1:
            out.extend(runs)
    return out


def blank_qc_type_from_filename(filename):
    """'blank_analysis', 'blank_extraction', 'blank_QC', 'pool_QC' or None."""
    tokens = _tokens_and_runs(name_tokens(filename))
    if not tokens:
        return None
    lower = {t.lower() for t in tokens}
    blank = any(_BLANK_TOKEN.match(t) for t in tokens)
    qc = any(_QC_TOKEN.match(t) for t in tokens)
    if qc and not blank and lower & _POOL_PLACE_TOKENS and \
            not any(_QC_TOKEN.match(t) and not _POOL_TOKEN.match(t) for t in tokens):
        qc = False
    if blank:
        if any(_EXTRACTION_BLANK_TOKEN.match(t) for t in tokens) or any(_EXTRACTION_WORD.match(t) for t in tokens):
            return 'blank_extraction'
        return 'blank_QC' if qc else 'blank_analysis'
    if qc:
        return 'pool_QC'
    return None


def _is_set(value):
    return str(value).strip().lower() not in ('', 'nan', 'none', 'missing value', 'not applicable', 'not specified')


def apply_filename_blank_qc(df, filename_column='filename', verbose=True):
    """Labels rows whose file name says blank / pooled QC (see module docstring) and whose
    SampleType is not already a blank, QC, culture-medium blank, reference material, personal
    care product or standard. SampleType and SampleTypeSub1 are set to the blank/QC type; the
    subject fields are then cleared by complete_and_fill_REDU_table (a blank's organism and body
    part too; a pooled QC keeps its own). Returns the frame and the number of changed rows."""
    if filename_column not in df.columns or len(df) == 0:
        return df, 0
    df = df.copy()
    if 'SampleType' not in df.columns:
        df['SampleType'] = ''
    if 'SampleTypeSub1' not in df.columns:
        df['SampleTypeSub1'] = ''
    current = df['SampleType'].fillna('').astype(str).str.strip()
    already = current.str.contains('blank', case=False) | current.isin(
        ['pool_QC', 'reference material', 'personal care product', 'molecular_standards'])
    types = df[filename_column].map(lambda f: blank_qc_type_from_filename(f) if pd.notna(f) else None)
    change = types.notna() & ~already
    if change.any():
        df.loc[change, 'SampleType'] = types[change]
        df.loc[change, 'SampleTypeSub1'] = types[change]
        if verbose:
            counts = types[change].value_counts().to_dict()
            print(f"File-name blank/QC rule: {int(change.sum())} rows relabelled {counts}")
    return df, int(change.sum())


# ---- QC samples that are no pool of the study's samples (by their name or their values) ----
# File-name / sample-name words of QC injections that are no pool of the study's samples (sources label them
# pool_QC or QC too: 'calcurve_6mix', 'Std_Mix', 'QC_ESTD10000nM', 'PBQC_ISTDS', 'QCSolvent', 'acn1',
# 'MediaControl', 'SST', 'triazines_10ngml', 'ISTD_QAQC', 'p1Std'). As letter runs of the base name.
NOT_POOL_WORDS = {'blank', 'blk', 'std', 'stds', 'standard', 'standards', 'stdmix', 'estd', 'istd', 'istds', 'is', 'cal',
                  'calib', 'calibration', 'calcurve', 'curve', 'mix', 'mixture', 'solvent', 'acn', 'meoh', 'methanol',
                  'water', 'h', 'o', 'buffer', 'media', 'medium', 'sst', 'crude', 'reaction', 'spike', 'spiked', 'nm', 'um',
                  'mm', 'ngml', 'ugml', 'mgml', 'ppb', 'ppm', 'ngl', 'ugl'}
# the same for the values of a sample ('Sample source: Standard', 'Solvent Mixture', 'QC(Internal Standard)',
# 'Quality standards'): not the short words and units, which are also parts of ordinary values
NOT_POOL_VALUE_WORDS = {'blank', 'blk', 'std', 'stds', 'standard', 'standards', 'stdmix', 'estd', 'istd', 'istds',
                        'calib', 'calibration', 'calcurve', 'solvent', 'solvents', 'acn', 'meoh', 'methanol', 'sst'}
# the same for the names of QC samples a source calls QC (Workbench '_QC_' / '_POOLED_' sample sources): only
# the words that say blank, standard, calibration or solvent; not 'mix' ('QC mix', 'Tech-Mix': a mix of the
# study's samples), media ('Pool_QC_Media': a pool of the study's media samples), 'IS', units ('UM053')
NOT_POOL_QC_NAME_WORDS = {'blank', 'blk', 'std', 'stds', 'standard', 'standards', 'stdmix', 'estd', 'istd', 'istds', 'cal',
                          'calib', 'calibration', 'calcurve', 'solvent', 'solvents', 'acn', 'meoh', 'methanol', 'sst'}
# words found inside longer runs too ('C18negblank1', 'blankH2O')
NOT_POOL_SUBSTRINGS = ('blank', 'solvent', 'standard', 'calibrat')
# water as 'H2O' (not hydrogen peroxide 'H2O2')
_H2O = re.compile(r'h2o(?!\d)')
_LETTER_RUNS = re.compile(r'[a-z]+')
_NAME_EXTENSIONS = re.compile(r'(?i)(\.(mzml|mzxml|mzdata|mgf|raw|cdf|wiff2?|d|gz|zip|lzma|bz2|7z|xz|rar|scan))+$')


def not_a_pool_by_name(filename, words=None):
    """whether a file's (or sample's) base name says it is a blank, standard, calibration, solvent or medium
    injection (its letter runs, camel case split: 'QCSolvent', 'MediaControl'); words: the words to use
    (default NOT_POOL_WORDS)"""
    words = NOT_POOL_WORDS if words is None else words
    base = str(filename).rstrip('/').rsplit('/', 1)[-1]
    base = _NAME_EXTENSIONS.sub('', base)
    split = re.sub(r'(?<=[a-z])(?=[A-Z])', ' ', base).lower()
    base = base.lower()
    if any(w in base for w in NOT_POOL_SUBSTRINGS):
        return True
    # letter runs with and without the camel-case split ('MeOH': 'meoh', 'StdMix': 'std', 'mix', 'stdmix')
    runs = set(_LETTER_RUNS.findall(split)) | set(_LETTER_RUNS.findall(base))
    # 'h2o' -> 'h', 'o': only together
    if {'h', 'o'} <= words and {'h', 'o'} <= runs and _H2O.search(base):
        return True
    return bool(runs & (words - {'h', 'o'}))


def not_a_pool_by_value(value):
    """whether a sample's value names a standard, calibration, solvent or blank ('Standard', 'Solvent Mixture',
    'QC(Internal Standard)'), NOT_POOL_VALUE_WORDS as letter runs"""
    text = str(value)
    low = text.lower()
    if any(w in low for w in NOT_POOL_SUBSTRINGS):
        return True
    runs = set(_LETTER_RUNS.findall(re.sub(r'(?<=[a-z])(?=[A-Z])', ' ', text).lower())) | set(_LETTER_RUNS.findall(low))
    return bool(runs & NOT_POOL_VALUE_WORDS)
