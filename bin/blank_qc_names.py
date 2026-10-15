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


def blank_qc_type_from_filename(filename):
    """'blank_analysis', 'blank_extraction', 'blank_QC', 'pool_QC' or None."""
    tokens = name_tokens(filename)
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
    organism, body part and other subject fields are then cleared by complete_and_fill_REDU_table,
    as for any blank. Returns the frame and the number of changed rows."""
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
