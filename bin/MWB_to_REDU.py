import os
import argparse
import pandas as pd
import re
import numpy as np
import requests
from bs4 import BeautifulSoup
from urllib.parse import urlparse, parse_qs
import json
from tqdm import tqdm
from collections import Counter
import traceback
import time
from REDU_conversion_functions import get_taxonomy_id_from_name__allowedTerms
from REDU_conversion_functions import map_instrument_to_allowed
from REDU_conversion_functions import convert_smoking, bmi_to_numeric, convert_diet
from read_and_validate_redu_from_github import complete_and_fill_REDU_table
from REDU_conversion_functions import age_category
from REDU_conversion_functions import get_taxonomy_info
from REDU_conversion_functions import normalize_redu_datetime
from ms_instrument_matching import match_mass_spectrometer
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry


# Shared session for Metabolomics Workbench: retries throttling/server errors with backoff
mwb_session = requests.Session()
mwb_session.mount("https://", HTTPAdapter(max_retries=Retry(total=5, backoff_factor=2,
                                                             status_forcelist=[429, 500, 502, 503, 504],
                                                             allowed_methods=["GET"])))
MWB_TIMEOUT = (10, 300)


def _mwb_get(url):
    response = mwb_session.get(url, timeout=MWB_TIMEOUT)
    response.raise_for_status()
    return response


def is_apple_double(path):
    # macOS resource-fork copies (__MACOSX/.../._x.mzML) are not data files
    return "__MACOSX/" in path or os.path.basename(path).startswith("._")


_POS_RE = re.compile(r'(?:^|[^a-z]|hilic|rp|c18|esi|lipid|lipids)(?:pos|positive)(?:[^a-z]|$)')
_NEG_RE = re.compile(r'(?:^|[^a-z]|hilic|rp|c18|esi|lipid|lipids)(?:neg|negative)(?:[^a-z]|$)')
_CHROM_TOKENS = {
    'hilic': re.compile(r'(?:^|[^a-z])hilic|amide|zic'),
    'reverse': re.compile(r'(?:^|[^a-z])(?:c18|rp|rplc|reverse|reversed|c8|c30|t3)(?:[^a-z]|pos|neg|$)'),
}


def _path_indicates(path):
    """Polarity ('positive'/'negative'/None) and chromatography ('hilic'/'reverse'/None)
    stated by tokens in a raw file path, e.g. 'C18-NEG/0001.raw' or 'x_HILICpos_1.mzML'."""
    p = str(path).lower()
    pos = bool(_POS_RE.search(p) or re.search(r'(?<![a-z])(?:esi|apci|hesi)\s*\+', p))
    neg = bool(_NEG_RE.search(p) or re.search(r'(?<![a-z])(?:esi|apci|hesi)\s*-(?![a-z])', p))
    polarity = 'positive' if pos and not neg else 'negative' if neg and not pos else None
    chrom_hits = [k for k, rx in _CHROM_TOKENS.items() if rx.search(p)]
    chrom = chrom_hits[0] if len(chrom_hits) == 1 else None
    return polarity, chrom


# --- ChromatographyAndPhase from CHROMATOGRAPHY_TYPE + COLUMN_NAME ---------------------
# Column names are compared after lower-casing and replacing every non-alphanumeric
# character (hyphens of any kind, '(', '™', '_') by a space.
# 5%-phenyl methylpolysiloxane GC phases (the DB-5 family under its vendor names)
_GC_DB5_RE = re.compile(r' (?:db|hp|rtx|rxi|zb|vf|tg|elite|slb|equity|bpx|optima|inertcap|sh i|ov|velocity) ?5(?![0-9])'
                        r'|(?<![0-9])5 ?% ?(?:di)?phenyl')
# other GC phases: a column list naming one of these is not (only) a DB-5 type column
_GC_OTHER_RE = re.compile(r'wax|ffap|fame|aaa|porabond|semivolatiles|624|1701|cp sil'
                          r'| (?:db|hp|rtx|rxi|zb|vf|tg) ?(?:1|17|23|35|50|200|225)(?![0-9])')
# LC column words in a GC analysis: the analysis is mis-typed or combined, leave it empty
_LC_COLUMN_RE = re.compile(r' c ?(?:18|8|30)(?![0-9])|hilic|amide|hss|beh|kinetex|acquity')
_RP_FAMILIES = {
    'polar-C18': re.compile(r' polar ?c ?18(?![0-9])'),
    'Phenyl-Hexyl': re.compile(r'phenyl ?hexyl'),
    'C18': re.compile(r' [a-z]*c ?18(?![0-9])| t3 | ods | octadecyl'),
    'C8': re.compile(r' [a-z]*c ?8(?![0-9])'),
    'C30': re.compile(r' [a-z]*c ?30(?![0-9])'),
}
_HILIC_COLUMN_RE = re.compile(r'hilic|amide|zic')
# chromatography types that name no phase; only here does the column decide the phase
_NO_PHASE_TYPES = {'', '-', 'na', 'lc', 'hplc', 'uplc', 'uhplc', 'lc ms', 'other', 'unspecified', 'csh', 'rp', 'rplc'}


def _norm_column(text):
    return ' ' + re.sub(r'[^a-z0-9%]+', ' ', str(text).lower()).strip() + ' '


def _reverse_phase_subtype(col):
    if 'pfp' in col:  # pentafluorophenyl phases, also when sold as 'C18-PFP'
        return 'NOS'
    hits = [name for name, rx in _RP_FAMILIES.items() if rx.search(col)]
    if 'polar-C18' in hits:
        hits.remove('C18')
    return hits[0] if len(hits) == 1 else 'NOS'


def chromatography_and_phase(chrom_type, column_name):
    """ReDU ChromatographyAndPhase from the mwTab CHROMATOGRAPHY_TYPE and COLUMN_NAME.
    The type gives the separation (reverse phase / HILIC / GC); the column gives the
    phase. GC analyses only get the DB-5 term when the column is a DB-5 type column,
    otherwise the generic GC term. Returns '' where nothing can be said."""
    t = '' if chrom_type is None or (isinstance(chrom_type, float) and np.isnan(chrom_type)) else str(chrom_type).lower()
    col = _norm_column('' if column_name is None or (isinstance(column_name, float) and np.isnan(column_name)) else column_name)
    if 'gc' in t:
        if _LC_COLUMN_RE.search(col):
            return ''
        if _GC_DB5_RE.search(col) and not _GC_OTHER_RE.search(col):
            return 'gas chromatography (DB-5)'
        return 'gas chromatography (NOS)'
    if 'scherzo sm c18' in col:
        return 'mixed mode (Scherzo SM-C18)'
    if 'scherzo ss c18' in col:
        return 'mixed mode (Scherzo SS-C18)'
    is_rp = 'reverse' in t or t.strip() in ('rp', 'rplc')
    is_hilic = 'hilic' in t or 'normal phase' in t
    if is_rp and is_hilic:
        return ''  # e.g. 'Reversed phase/HILIC': cannot tell which
    if is_hilic:
        return 'normal phase (HILIC)'
    if not is_rp and t.strip() in _NO_PHASE_TYPES:
        # no phase in the type: a column that is clearly one family decides
        rp_hit = any(rx.search(col) for rx in _RP_FAMILIES.values()) and 'pfp' not in col
        hilic_hit = bool(_HILIC_COLUMN_RE.search(col))
        if hilic_hit and not rp_hit:
            return 'normal phase (HILIC)'
        if rp_hit and not hilic_hit:
            is_rp = True
    if is_rp:
        if 'hypercarb' in col or 'porous graphitic' in col:
            return 'porous graphitic carbon (PGC)'
        if 'exchange' in t or 'exchange' in col:
            return 'reverse phase (NOS)'  # combined with an ion-exchange separation
        return 'reverse phase ({})'.format(_reverse_phase_subtype(col))
    return ''


_MS_ALIAS_CACHE = {}


def mwb_ms_alias(name, instrument_type, allowed_ms, path_to_csvs='translation_sheets'):
    """Allowed MassSpectrometer term for a Workbench INSTRUMENT_NAME listed in
    MWB_MassSpectrometer_aliases.csv (exact name after lower-casing / whitespace
    normalisation; a row with an INSTRUMENT_TYPE only applies to that type). Only
    targets that exist in the allowed terms are returned."""
    if path_to_csvs not in _MS_ALIAS_CACHE:
        table = {}
        path = os.path.join(path_to_csvs, 'MWB_MassSpectrometer_aliases.csv')
        if os.path.exists(path):
            for row in pd.read_csv(path, dtype=str, keep_default_na=False).itertuples(index=False):
                table[(' '.join(row.MWB.lower().split()), row.INSTRUMENT_TYPE.strip().lower())] = row.REDU
        _MS_ALIAS_CACHE[path_to_csvs] = table
    table = _MS_ALIAS_CACHE[path_to_csvs]
    key = ' '.join(str(name).lower().split())
    itype = ' '.join(str(instrument_type).lower().split()) if isinstance(instrument_type, str) else ''
    term = table.get((key, itype)) or table.get((key, ''))
    return term if term in set(allowed_ms) else None


def study_level_sex(gender):
    """BiologicalSex from the study-level SUBJECT GENDER when it names one sex."""
    if not isinstance(gender, str) or not gender.strip():
        return None
    sex = convert_sex(gender)
    return sex if sex in ('male', 'female') else None


_AGE_UNIT_YEARS = {'year': 1, 'yr': 1, 'month': 1 / 12, 'mo': 1 / 12, 'week': 1 / 52.1429, 'wk': 1 / 52.1429, 'day': 1 / 365}
_STUDY_AGE_RE = re.compile(r'(\d+(?:\.\d+)?)[ -]?(years?|yrs?|months?|mo|weeks?|wks?|wk|days?)(?:[ -]old)?')


def study_level_age_years(age):
    """AgeInYears from the study-level SUBJECT AGE_OR_AGE_RANGE when it is exactly one
    age with a unit ('8 weeks', '6-weeks old', '4-month-old', '25 wk'). Ranges, bare
    numbers, approximate and described ages are not used."""
    if not isinstance(age, str):
        return None
    m = _STUDY_AGE_RE.fullmatch(' '.join(age.lower().split()))
    if not m:
        return None
    unit = m.group(2).rstrip('s')
    factor = _AGE_UNIT_YEARS.get(unit) or _AGE_UNIT_YEARS.get(unit[:2])
    return round(float(m.group(1)) * factor, 4) if factor else None


def _phase_family(values):
    v = values.astype(str)
    v = v.where(~v.str.startswith('reverse phase ('), 'reverse phase (NOS)')
    return v.where(~v.str.startswith('gas chromatography ('), 'gas chromatography (NOS)')


def resolve_analysis_conflicts_from_path(df):
    """A raw file listed under several analyses of a study (e.g. a positive and a
    negative mode analysis) gets conflicting analysis-level values, which are later
    blanked to 'missing value'. When the file path itself states the polarity and/or
    the chromatography (C18-NEG/..., ..._HILICpos...), keep only the analysis rows that
    agree with it - but only if that leaves rows that agree with each other on the
    analysis-level columns. Otherwise the group is left unchanged."""
    analysis_cols = [c for c in ['IonizationSourceAndPolarity', 'ChromatographyAndPhase', 'MassSpectrometer'] if c in df.columns]
    dup = df['filename'].duplicated(keep=False)
    if not dup.any() or not analysis_cols:
        return df
    drop_idx = []
    for filename, group in df[dup].groupby('filename'):
        if all(group[c].astype(str).nunique() == 1 for c in analysis_cols):
            continue
        polarity, chrom = _path_indicates(filename)
        if polarity is None and chrom is None:
            continue
        keep = pd.Series(True, index=group.index)
        if polarity is not None:
            keep &= group['IonizationSourceAndPolarity'].astype(str).str.lower().str.contains(polarity)
        if chrom is not None and 'ChromatographyAndPhase' in group.columns:
            keep &= group['ChromatographyAndPhase'].astype(str).str.lower().str.contains(chrom)
        kept = group[keep]
        if len(kept) == 0 or len(kept) == len(group):
            continue
        # analyses that differ only in the reverse-phase (or GC) column type agree on
        # the polarity; merge_repeated_fileobservations_across_mwatb then keeps the shared phase
        if all(_phase_family(kept[c]) .nunique() == 1 if c == 'ChromatographyAndPhase' else kept[c].astype(str).nunique() == 1
               for c in analysis_cols):
            drop_idx.extend(group.index[~keep])
    if drop_idx:
        print(f"Resolved analysis conflicts from file paths: {len(drop_idx)} conflicting analysis rows dropped")
    return df.drop(index=drop_idx)


def propagate_to_same_name_files(redu_df, raw_file_name_df):
    """Studies often hold one raw file per sample under several method folders
    (Method_1/Sample_1.raw, Method_2/Sample_1.raw, ...), but the sample-to-file match
    links each sample to one of them. Give the other same-name files (same basename and
    extension) the sample's metadata, when the basename belongs to exactly one sample.
    Analysis-level fields are copied only if the study has a single analysis setup;
    otherwise they come from the file path tokens, or stay missing."""
    if redu_df is None or len(redu_df) == 0 or raw_file_name_df is None or len(raw_file_name_df) == 0:
        return redu_df
    out = redu_df.copy()
    out['_base'] = out['filename'].astype(str).map(lambda x: os.path.basename(x).lower())
    base_counts = out['_base'].value_counts()
    unique_rows = out[out['_base'].map(base_counts) == 1].set_index('_base')
    if len(unique_rows) == 0:
        return redu_df
    files = raw_file_name_df[['filename', 'USI']].drop_duplicates('filename').copy()
    files['_base'] = files['filename'].astype(str).map(lambda x: os.path.basename(x).lower())
    files = files[~files['filename'].isin(set(out['filename'])) & files['_base'].isin(unique_rows.index)]
    if len(files) == 0:
        return redu_df
    analysis_cols = [c for c in ['IonizationSourceAndPolarity', 'ChromatographyAndPhase', 'MassSpectrometer'] if c in out.columns]
    single_setup = all(out[c].astype(str).nunique() == 1 for c in analysis_cols)
    copies = unique_rows.loc[files['_base']].reset_index(drop=True)
    copies['filename'] = files['filename'].values
    copies['USI'] = files['USI'].values
    if not single_setup:
        # a field with one value across the study's files (e.g. the instrument) is still
        # known for the copies
        for c in analysis_cols:
            if out[c].astype(str).nunique() != 1:
                copies[c] = 'missing value'
        for i, path in enumerate(copies['filename']):
            polarity, chrom = _path_indicates(path)
            match = out
            if polarity is not None:
                match = match[match['IonizationSourceAndPolarity'].astype(str).str.lower().str.contains(polarity)]
            if chrom is not None and 'ChromatographyAndPhase' in out.columns:
                match = match[match['ChromatographyAndPhase'].astype(str).str.lower().str.contains(chrom)]
            if (polarity is not None or chrom is not None) and len(match) and all(match[c].astype(str).nunique() == 1 for c in analysis_cols):
                for c in analysis_cols:
                    copies.at[i, c] = match[c].iloc[0]
    print(f"Copied sample metadata to {len(copies)} same-name raw files in other folders")
    return pd.concat([redu_df, copies[redu_df.columns]], ignore_index=True)


def merge_repeated_fileobservations_across_mwatb(df, **kwargs):

    if 'polarity_table' in kwargs.keys():
        polarity_table = kwargs['polarity_table']
    else:
        polarity_table = pd.DataFrame(columns=['USI', 'negative_polarity_count', 'positive_polarity_count'])


    def process_group(group):
        # Initialize the list for columns to make the most common
        li_make_most_common = []
        
        # Check if all entries in 'NCBITaxonomy' for the group are the same
        if len(set(group['NCBITaxonomy'])) == 1:
            li_make_most_common = ['SampleType', 'SampleTypeSub1']
        
            # Find the most common value for specified columns and update the group
            for col in li_make_most_common:
                most_common_value = Counter(group[col]).most_common(1)[0][0]
                if most_common_value is None or most_common_value == 'None':
                    most_common_value = 'missing value'
                group[col] = most_common_value

        # one file under several analyses with different reverse-phase (or GC) columns:
        # what they share is still known
        if 'ChromatographyAndPhase' in group.columns:
            phases = set(group['ChromatographyAndPhase'].dropna().astype(str)) - {'missing value'}
            if len(phases) > 1 and len(phases) == group['ChromatographyAndPhase'].astype(str).nunique():
                if all(p.startswith('reverse phase (') for p in phases):
                    group['ChromatographyAndPhase'] = 'reverse phase (NOS)'
                elif all(p.startswith('gas chromatography (') for p in phases):
                    group['ChromatographyAndPhase'] = 'gas chromatography (NOS)'

        # Handle discrepancies for columns not in 'li_make_most_common'
        li_rm_discrepancies = [col for col in group.columns if col not in li_make_most_common]
        for col in li_rm_discrepancies:
            if len(set(group[col])) > 1 or all(value is None for value in group[col]):
                group[col] = 'missing value'
            else:
                # If there's no discrepancy, keep the original value
                group[col] = group[col].iloc[0]
        
        return group

    df = df.merge(polarity_table, left_on='USI', right_on='USI', how='left')

    # A file with many scans of both polarities is a polarity-switching file even when
    # its analysis is declared positive or negative: call it alternating instead of
    # dropping it below (whole studies were lost this way)
    try:
        pos_n = pd.to_numeric(df['positive_polarity_count'], errors='coerce').fillna(0)
        neg_n = pd.to_numeric(df['negative_polarity_count'], errors='coerce').fillna(0)
        switching = (pos_n > 0) & (neg_n > 0) & (np.minimum(pos_n, neg_n) >= 0.05 * (pos_n + neg_n))
        allowed_ion = set(((globals().get('allowedTerm_dict') or {}).get('IonizationSourceAndPolarity') or {}).get('allowed_values') or [])
        ion = df['IonizationSourceAndPolarity'].astype(str)
        alt = ion.str.replace(r'\((?:positive|negative)\)$', '(alternating)', regex=True)
        relabel = switching & (alt != ion) & alt.isin(allowed_ion)
        if relabel.any():
            print(f"{int(relabel.sum())} polarity-switching files relabelled as alternating")
            df.loc[relabel, 'IonizationSourceAndPolarity'] = alt[relabel]
    except Exception as e:
        print(f"Could not check polarity-switching files: {type(e).__name__}: {e}")

    df = df[~((df['IonizationSourceAndPolarity'].str.contains('negative')) & (df['positive_polarity_count'] > 0))]
    df = df[~((df['IonizationSourceAndPolarity'].str.contains('positive')) & (df['negative_polarity_count'] > 0))]
    df = df[~((df['IonizationSourceAndPolarity'].str.contains('alternating')) & ((df['negative_polarity_count'] == 0) | (df['positive_polarity_count'] == 0)))]

    try:
        df = resolve_analysis_conflicts_from_path(df.reset_index(drop=True))
    except Exception as e:
        print(f"Could not resolve analysis conflicts from file paths: {type(e).__name__}: {e}")


    # Group by 'USI' and apply the processing function to each group
    df = df.groupby('filename').apply(process_group)

    # Remove duplicates based on 'USI' column
    df = df.drop_duplicates(subset='filename')

    return df

def handle_duplicates(ordered_pairs):
    """
    Append numbers to duplicate keys. This is necessary because the RAW_FILE_NAME key often appears more than once
    in the mwTAb files
    """
    d = {}
    for k, v in ordered_pairs:
        if k in d:
            count = 1
            new_key = f"{k}_{count}"
            while new_key in d:
                count += 1
                new_key = f"{k}_{count}"
            d[new_key] = v
        else:
            d[k] = v
    return d


# Some archive listings carry raw 7z listing columns in front of the path,
# e.g. '2014-10-05 10:45:02 ..H.A           80               MX183900 .../x.d'
_7Z_LISTING_PREFIX = re.compile(r'^\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2} \S{5} +\d+ +(?:\d+ +)?')


def clean_listing_path(path):
    return _7Z_LISTING_PREFIX.sub('', path)


def collapse_d_folder(path):
    """Agilent/Bruker .d raw data are folders, listed only as the files inside them;
    represent them by the .d folder (as getAllWorkbench_file_paths.py does)."""
    parts = path.split('/')
    for i, part in enumerate(parts):
        if part.lower().endswith('.d'):
            return '/'.join(parts[:i + 1])
    return path


def _get_metabolomicsworkbench_files(dataset_accession):
    try:
        dataset_list_url = "https://www.metabolomicsworkbench.org/data/show_archive_contents_json.php?STUDY_ID={}".format(
            dataset_accession)
        mw_file_list = _mwb_get(dataset_list_url).json()

        acceptable_extensions = [".mzml", ".mzxml", ".cdf", ".raw", ".wiff", ".d"]

        mw_file_list = [file_obj for file_obj in mw_file_list if not is_apple_double(file_obj["FILENAME"])]
        workbench_df = pd.DataFrame(mw_file_list, columns=["FILENAME", "FILESIZE", "URL"])
        workbench_df["filename"] = workbench_df["FILENAME"].apply(clean_listing_path).apply(collapse_d_folder)
        workbench_df = workbench_df[workbench_df["filename"].str.lower().str.endswith(tuple(acceptable_extensions))]
        workbench_df["FILESIZE"] = pd.to_numeric(workbench_df["FILESIZE"], errors='coerce').fillna(0)

        # Generate USI column
        workbench_df["USI"] = workbench_df["URL"].apply(
            lambda url: "mzspec:{}:{}".format(
                dataset_accession,
                collapse_d_folder(clean_listing_path(parse_qs(urlparse(url).query).get('F', [''])[0]))
            )
        )

        # one row per file / .d folder (summing the sizes of the files inside a folder)
        workbench_df = workbench_df.groupby("filename", as_index=False, sort=False).agg(
            FILESIZE=("FILESIZE", "sum"), USI=("USI", "first"))
        workbench_df["size_mb"] = (workbench_df["FILESIZE"] / 1024 / 1024).astype(int)

        workbench_df = workbench_df[["filename", "size_mb", "USI"]]
    except Exception as e:
        print(f"Failed to list files for {dataset_accession}: {type(e).__name__}: {e}")
        workbench_df = pd.DataFrame()

    return workbench_df



def convert_to_numeric_or_original(x):
    try:
        numeric_x = pd.to_numeric(x)
        return 'numeric'
    except ValueError:
        return x



def get_enviromental_water(x):
    x = x.lower()

    if 'water' in x or 'sewerage' in x:
        if 'waste' in x or 'sewerage' in x:
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


# --- QC / blank / reference material samples ------------------------------------------
# Exact values (after lower-casing and turning separators into spaces). Plain 'blank'
# and 'pool' only count under keys that state what kind of sample it is.
_QC_TYPE_KEYS = {'sample source', 'sample type', 'sampletype', 'samle type', 'type', 'sample class', 'class',
              'group', 'sample group', 'treatment', 'treatment group', 'factor', 'factor group', 'source',
              'origin', 'condition', 'category', 'sample category', 'stuy type', 'study type', 'material'}
_QC_VALUES = {'qc', 'qcs', 'quality control', 'quality control sample', 'qc sample', 'qc samples', 'test qc',
              'instrument qc', 'qc drift correction', 'qc driftcorrection', 'qc prerun', 'qc reference', 'qc recas',
              'qc externalstandard', 'qc external standard'}
_POOL_QC_VALUES = {'pooled qc', 'qc pooled', 'pool qc', 'qc pool', 'pooled qc sample', 'pooled qc samples', 'study pool',
                   'pooled quality control', 'pooled quality control sample', 'pooled biological quality control',
                   'sample qc pooled study samples', 'pooled plasma quality control sample', 'study total pool',
                   'pooled study plasma', 'pooled study sample', 'pooled study samples', 'qc pooled sample'}
_POOL_VALUES_QC_TYPE_KEYS = {'pool', 'sample pool'}
_BLANK_ANALYSIS_VALUES = {'solvent blank', 'instrument blank', 'pre run blank', 'post run blank', 'analytical blank',
                          'water blank', 'blank h2o', 'qc blank', 'blank solvent', 'solvent blanks', 'reagent blank'}
_BLANK_EXTRACTION_VALUES = {'process blank', 'process blanks', 'extraction blank', 'extraction blanks', 'method blank',
                            'method blanks', 'procedural blank', 'procedure blank', 'blank extraction', 'extracted blank'}
_BLANK_MEDIA_VALUES = {'media blank', 'blank media', 'medium blank', 'culture media blank'}
_REFERENCE_RE = re.compile(r'(?:pooled )?nist(?: standard)?(?: reference)?(?: material)?(?: srm)?(?: ?(?:1950|1957|1958|3672|3667|1953))?'
                           r'(?: (?:plasma|serum|urine|pool|urine pool))?(?: (?:srm )?(?:1950|1957|1958|3672|3667|1953))?(?: srm 1950)?'
                           r'|srm ?(?:1950|1957|1958|3672|3667|1953)(?: (?:plasma|serum|urine))?'
                           r'|nist plasma srm 1950|nist standard reference material srm 1950|pooled nist reference plasma')


def _norm_value(v):
    return re.sub(r'\s+', ' ', re.sub(r'[^a-z0-9+]+', ' ', str(v).lower())).strip()


def qc_sample_type(key, value):
    """ReDU SampleType for a QC / blank / reference-material sample, or None."""
    k, v = _norm_value(key), _norm_value(value)
    if not v:
        return None
    if v in _QC_VALUES or re.fullmatch(r'qc ?\d{1,2}', v):
        return 'blank_QC'
    if v in _POOL_QC_VALUES or (k in _QC_TYPE_KEYS and v in _POOL_VALUES_QC_TYPE_KEYS):
        return 'pool_QC'
    if v in _BLANK_ANALYSIS_VALUES or (k in _QC_TYPE_KEYS and v in ('blank', 'blanks')):
        return 'blank_analysis'
    if v in _BLANK_EXTRACTION_VALUES:
        return 'blank_extraction'
    if v in _BLANK_MEDIA_VALUES:
        return 'blank_culturemedia'
    if _REFERENCE_RE.fullmatch(v):
        return 'reference material'
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


_SEX_WORDS = {
    'female': {'f', 'female', 'females', 'woman', 'women', 'girl', 'girls', 'spayed', 'sow', 'hen', 'mare', 'doe', 'ewe'},
    'male': {'m', 'male', 'males', 'man', 'men', 'boy', 'boys', 'gelding', 'castrated', 'barrow', 'boar', 'stallion', 'rooster', 'buck', 'ram'},
}


def convert_sex(x):
    x = str(x).lower().strip()
    if x in ['f', 'female']:
        return 'female'
    elif x in ['m', 'male']:
        return 'male'
    elif x in ['as', 'asexual']:
        return 'asexual'
    # spelled-out / plural / animal-husbandry forms ('Males', 'BOY', 'Male (gelding)',
    # 'neutered male'): exactly one sex must be named
    tokens = set(re.split(r'[^a-z]+', x)) - {''}
    hits = [sex for sex, words in _SEX_WORDS.items() if tokens & words]
    if len(hits) == 1 and not tokens & {'and', 'mixed', 'both', 'pooled', 'pool'}:
        return hits[0]
    return 'missing value'


# Keys describing the sex of someone other than the sampled subject
_SEX_KEY_EXCLUDE_RE = r'pup|fetus|fetal|foetal|offspring|mother|maternal|father|paternal|donor|partner|infant|baby|child of|recipient'


def convert_sex_coded(value, key):
    """Numeric sex codes are only decoded when the key defines them, e.g.
    'sex (1 boy, 2 girl)' or 'gender: 1=male, 2=female'."""
    v = str(value).strip()
    if not re.fullmatch(r'\d+', v):
        return None
    codes = {}
    for code, word in re.findall(r'(\d+)\s*[=:\-]?\s*([a-z]+)', str(key).lower()):
        sex = convert_sex(word)
        if sex in ('male', 'female'):
            codes[code] = sex
    return codes.get(v)


# Depth / altitude in metres; depths are negative from sea level, altitudes positive
_DEPTH_KEY_RE = r'(?<![a-z])(depth|altitude|elevation)(?![a-z])'
_DEPTH_KEY_EXCLUDE_RE = r'sequenc|read|coverage|seq\b|snow|soil depth class'


def depth_or_altitude_meters(value, key):
    nums = re.findall(r'-?\d+(?:\.\d+)?', str(value))
    if len(nums) != 1:
        return None  # missing, or a range such as '0-10 cm'
    num = abs(float(nums[0]))
    unit_src = (str(key) + ' ' + str(value)).lower()
    if re.search(r'\bcm\b|centimet', unit_src):
        num /= 100
    elif re.search(r'\bmm\b|millimet', unit_src):
        num /= 1000
    elif re.search(r'\bkm\b|kilomet', unit_src):
        num *= 1000
    elif re.search(r'\bft\b|feet|foot', unit_src):
        num *= 0.3048
    if re.search(r'(?<![a-z])depth(?![a-z])', str(key).lower()) and num:
        num = -num
    return round(num, 3)


_DATE_KEYS = {'collection date', 'collection_date', 'sampling date', 'sampling_date', 'sample collection date',
              'date of collection', 'date collected', 'collection_datetime', 'sampling datetime', 'sample date',
              'date of sampling', 'collection date and time', 'date_collected', 'sampledate', 'date'}

# per-sample health keys; values are matched exactly (case-insensitive)
_HEALTH_KEY_RE = r'^(condition|phenotype|disease|diagnosis|disease state|disease status|disease_state|clinical diagnosis|health status|health_status|health state|pathology|group|subject group|cohort)$'
_HEALTHY_VALUES = {'healthy', 'healthy control', 'healthy controls', 'healthy volunteer', 'healthy volunteers',
                   'healthy subject', 'healthy subjects', 'healthy individual', 'healthy donor', 'healthy donors',
                   'normal healthy', 'healthy adult', 'healthy adults'}
# plain 'control' only means healthy under keys that are about disease
_DISEASE_KEY_RE = r'^(disease|diagnosis|disease state|disease status|disease_state|clinical diagnosis|health status|health_status|health state)$'

# lower-cased DOID label/exact synonym -> 'Label|DOID:x'; set from --path_to_doid_csv
DOID_LOOKUP = {}
# prepared DOID ontology table; set from --path_to_doid_csv
DOID_TABLE = None


def build_doid_lookup(doid_table):
    """Exact (case-insensitive) DOID label/synonym lookup. Names that point to more
    than one disease, the root term 'disease' and very short names are left out."""
    names = {}
    for idx, label, syn in doid_table[['DOIDOntologyIndex', 'Label', 'Synonym']].itertuples(index=False):
        if not isinstance(label, str) or idx == 'DOID:4':
            continue
        for name in {label, syn}:
            if isinstance(name, str) and len(name.strip()) >= 4:
                names.setdefault(name.strip().lower(), set()).add(f'{label}|{idx}')
    return {name: next(iter(v)) for name, v in names.items() if len(v) == 1}


_COUNTRY_ALIASES = {
    'usa': 'United States of America', 'u.s.a.': 'United States of America', 'u.s.a': 'United States of America',
    'united states': 'United States of America', 'u.s.': 'United States of America',
    'uk': 'United Kingdom', 'u.k.': 'United Kingdom', 'england': 'United Kingdom', 'scotland': 'United Kingdom',
    'wales': 'United Kingdom', 'great britain': 'United Kingdom', 'northern ireland': 'United Kingdom',
    'korea': 'South Korea', 'republic of korea': 'South Korea', 'the netherlands': 'Netherlands', 'holland': 'Netherlands',
    'p.r. china': 'China', "people's republic of china": 'China', 'pr china': 'China', 'viet nam': 'Vietnam',
    'russian federation': 'Russia', 'new mexico': 'United States of America',
}
# country names that are also common US states / names / words in location text
_COUNTRY_SKIP = {'georgia', 'jordan', 'chad', 'guinea', 'niger', 'dominica', 'turkey', 'mali', 'oman', 'peru'}


def country_from_location(text, allowed_countries):
    """Country named in a free-text location ('Children's Hospital, Boston, USA'):
    the longest country name/alias found as a whole phrase, if exactly one country."""
    if not isinstance(text, str) or not text.strip():
        return None
    t = ' ' + re.sub(r'[^a-z.\' ]+', ' ', text.lower()) + ' '
    names = {c.lower(): c for c in allowed_countries if ',' not in c and c.lower() not in _COUNTRY_SKIP}
    names.update({k: v for k, v in _COUNTRY_ALIASES.items() if v in allowed_countries})
    found = [(n, c) for n, c in names.items() if re.search(r'(?<![a-z])' + re.escape(n) + r'(?![a-z])', t)]
    # drop names contained in a longer match ('sudan' in 'south sudan')
    found = [(n, c) for n, c in found if not any(n != m and n in m for m, _ in found)]
    countries = {c for _, c in found}
    return next(iter(countries)) if len(countries) == 1 else None


_LIFESTAGE_VALUES = {
    'newborn': 'Infancy (<2 yrs)', 'neonate': 'Infancy (<2 yrs)', 'neonatal': 'Infancy (<2 yrs)',
    'infant': 'Infancy (<2 yrs)', 'infancy': 'Infancy (<2 yrs)',
    'adolescent': 'Adolescence (8 yrs < x <= 18 yrs)', 'adolescence': 'Adolescence (8 yrs < x <= 18 yrs)',
    'teenager': 'Adolescence (8 yrs < x <= 18 yrs)',
}
_LIFESTAGE_KEY_RE = r'^(life stage|life_stage|lifestage|developmental stage|development stage|age group|age_group)$'


def convert_to_years(age):
    num = int(re.findall(r'\b\d+\b', age)[0]) if re.findall(r'\b\d+\b', age) else None
    if num is None:
        return None
    if "week" in age:
        return num / 52.1429
    elif "day" in age:
        return num / 365
    elif "year" in age:
        return num
    else:
        return num


# Age Factor keys that describe a DIFFERENT subject than the sampled one (mother,
# fetus, household contact, index case, donor...) - never map these to the
# sample's AgeInYears.
AGE_EXCLUDE_RE = r'gestational|maternal|postnatal|postconceptional|foetal|\bfetal\b|\bmother\b|paternal|\bfather\b|\bdonor\b|contact|\bindex\b|\bbaby\b|neonat'
_AGE_CANONICAL = {'age', 'age (years)', 'age(years)', 'age (year)'}


def mwb_age_to_years(value, key):
    """Convert an age Factor value to years. Unit is taken from the KEY name first
    (e.g. 'age (months)', 'age in days'), falling back to any unit in the VALUE
    string, else years. Number extraction matches convert_to_years exactly, so for
    the previously-handled keys ('age', 'age (years)') the result is identical."""
    value = '' if value is None else str(value)
    nums = re.findall(r'\b\d+\b', value)
    if not nums:
        return None
    num = int(nums[0])
    src = str(key).lower()
    if 'day' not in src and 'week' not in src and 'month' not in src and 'year' not in src:
        src = value.lower()  # no unit in key -> use value (preserves old behavior)
    if 'day' in src:
        return num / 365
    if 'week' in src:
        return num / 52.1429
    if 'month' in src:
        return num / 12
    return num  # years (default)


def get_age_years_into_outer(df, new_col='MWB_age'):
    """Per-sample AgeInYears (in years) from any age Factor whose key isn't an
    excluded-subject one. Canonical 'age'/'age (years)' keys win over unit-named
    variants so studies that already have them produce identical output; studies
    that only have e.g. 'age (months)' now get a correctly-converted value."""
    keyl = df['Key'].astype(str).str.lower()
    is_age = keyl.str.contains(r'\bage\b', regex=True, na=False) & ~keyl.str.contains(AGE_EXCLUDE_RE, regex=True, na=False)
    tmp = df.loc[is_age, ['filename', 'Key', 'Value']].copy()
    if len(tmp) == 0:
        df[new_col] = 'NA'
        return df
    tmp['_yr'] = tmp.apply(lambda r: mwb_age_to_years(r['Value'], r['Key']), axis=1)
    tmp = tmp.dropna(subset=['_yr'])
    if len(tmp) == 0:
        df[new_col] = 'NA'
        return df
    tmp['_rank'] = tmp['Key'].astype(str).str.lower().apply(lambda k: 0 if k in _AGE_CANONICAL else 1)
    tmp = tmp.sort_values(['filename', '_rank'])
    best = tmp.drop_duplicates('filename', keep='first').set_index('filename')['_yr']
    df[new_col] = df['filename'].map(best)
    df[new_col] = df[new_col].fillna('NA')
    return df


def extract_years(dates):
    try:
        # the year of MM-DD-YYYY / MM/DD/YYYY / DD.MM.YYYY dates (the '-' form used to
        # raise and return no year), else any 4-digit number (YYYY-MM-DD, '2016')
        matches = re.findall(r'\d{1,2}[-/.]\d{1,2}[-/.](\d{4})', dates)
        if not matches:
            matches = re.findall(r'\d{4}', dates)
        years = max([int(year) for year in matches])
        if isinstance(years, int):
            if years >= 1998 and years <= 2030:
                return years
            else:
                return None
        else:
            return None
    except Exception:
        return None


def _norm_key(key):
    """Key with every run of non-alphanumeric characters ('_', '.', '-') as one space."""
    return re.sub(r'[^a-z0-9]+', ' ', str(key).lower()).strip()


def _fill_from_rows(df, values, col):
    """Fill df[col] where it is 'NA'/missing with the per-sample value of `values`
    (a Series on a subset of df's rows); only samples with exactly one distinct
    non-missing value are filled."""
    values = values.dropna()
    values = values[~values.astype(str).str.strip().str.lower().isin(['', 'na', 'nan', 'none', 'missing value'])]
    if len(values) == 0:
        return df
    per_sample = pd.DataFrame({'filename': df.loc[values.index, 'filename'], 'v': values.astype(object)})
    agreed = per_sample.groupby('filename')['v'].agg(lambda s: s.iloc[0] if s.astype(str).nunique() == 1 else None).dropna()
    current = df[col]
    missing = current.isna() | current.astype(str).str.strip().str.lower().isin(['', 'na', 'nan', 'none', 'missing value'])
    fill = missing & df['filename'].isin(agreed.index)
    df.loc[fill, col] = df.loc[fill, 'filename'].map(agreed)
    return df


# fallback age keys that do not give the age at sampling or are not one age
_AGE_FB_EXCLUDE_RE = r'diagnos|\bdx\b|onset|\bstart\b|\bunits?\b|\brange\b|\bgroup\b|\bcat\b|categor|\bdich\b|relative|worker'
_SINGLE_AGE_VALUE_RE = re.compile(r'\d+(?:\.\d+)?(?: ?(?:years?|months?|weeks?|days?))?', re.I)

_EXTRA_SMOKING_KEYS = {'smoke status', 'smoker participant', 'ever smoked?', 'ever smoker (cigarette, cigar, etc)',
                       'smoking condition', 'smoke now?', 'smoking (never 0, current 1, past 2)'}
_SMOKING_CODE_WORDS = {'never': 'never smoker', 'current': 'current smoker', 'past': 'former smoker', 'former': 'former smoker'}


def convert_smoking_coded(value, key):
    """Numeric smoking codes, only where the key defines them,
    e.g. 'smoking (never 0, current 1, past 2)'."""
    v = str(value).strip()
    if not re.fullmatch(r'\d+', v):
        return None
    codes = {}
    for word, code in re.findall(r'([a-z]+)\s*[=:]?\s*(\d+)', str(key).lower()):
        if word in _SMOKING_CODE_WORDS:
            codes[code] = _SMOKING_CODE_WORDS[word]
    return codes.get(v)


# Subject ID fields that hold a sample category instead of an identifier
_LABEL_SUBJECT_ID_RE = re.compile(r'(?<![a-z0-9])(?:samples?|study|qc|quality control|pool|pooled|blank|blanks|standards?|'
                                  r'reference|nist|unknown|na|n/a|none|not applicable)(?![a-z0-9])')


def _drop_label_subject_ids(df):
    """Subject IDs that are labels rather than identifiers become '': values naming a
    sample category ('Plasma Study Sample', 'Quality Control'), values equal to one of
    the sample's own metadata values (a group such as 'wild-type' or the sample source),
    and a study whose samples all carry one and the same Subject ID."""
    sid = df['SubjectIdentifierAsRecorded'].fillna('').astype(str).str.strip()
    low = sid.str.lower()
    bad = low.str.contains(_LABEL_SUBJECT_ID_RE, regex=True) & (low != '')
    key_l = df['Key'].astype(str).str.lower()
    is_meta = ~key_l.map(lambda k: k in _EXACT_RAW_KEYS or is_variant_raw_file_key(k)) & (key_l != 'sample id')
    own_values = (df.loc[is_meta, ['filename', 'Value']].assign(v=lambda t: t['Value'].astype(str).str.strip().str.lower())
                  .groupby('filename')['v'].agg(set))
    bad |= pd.Series([v != '' and v in own_values.get(f, set()) for f, v in zip(df['filename'], low)], index=df.index)
    per_sample = pd.DataFrame({'filename': df['filename'], 'sid': low}).query("sid != ''").drop_duplicates()
    if per_sample['filename'].nunique() > 1 and per_sample['sid'].nunique() == 1:
        bad |= low != ''
    return sid.where(~bad, '')


_NCBI_NAME_CACHE = {}


def _ncbi_name_lookup(allowedTerm_dict):
    """lower-cased NCBI name -> first allowed 'ID|Name' term with that name."""
    terms = allowedTerm_dict['NCBITaxonomy']['allowed_values']
    key = (id(terms), len(terms))
    if key not in _NCBI_NAME_CACHE:
        lookup = {}
        for entry in terms:
            parts = entry.split('|')
            if len(parts) == 2:
                lookup.setdefault(parts[1].lower(), parts[0] + '|' + parts[1])
        _NCBI_NAME_CACHE.clear()
        _NCBI_NAME_CACHE[key] = lookup
    return _NCBI_NAME_CACHE[key]


def get_key_info_into_outer(df, key_vars, new_col):
    # Match keys case-insensitively: MWB depositors use 'Sex'/'Gender'/'Age' etc.
    # with varied capitalization, and the fixed key_vars lists are lowercase, so a
    # case-sensitive match silently drops capitalized keys. Lower-casing both sides
    # only ever matches MORE of the same semantic keys (never fewer), so it cannot
    # remove a value we already capture.
    key_set = {str(k).lower() for k in key_vars}
    df[new_col] = df.loc[df['Key'].astype(str).str.lower().isin(key_set), 'Value']
    df[new_col] = df.groupby('filename')[new_col].transform('first')
    if new_col == 'Latitude' or new_col == 'Longitude':
        df[new_col] = pd.to_numeric(df[new_col], errors='coerce')
    else:
        df[new_col] = df[new_col].fillna(value='NA')
    print(f"Found {df[new_col].nunique()} unique values for {new_col}")
    return df

_RAW_KEY_RE = re.compile(r'^(?:raw_file_name|raw file name|rawfilename|raw_file|raw file|raw files|datafile name|'
                         r'rawfile name|rpm_filename|zhp_filename|data file|raw_data)'
                         r'(?:\s*\(([^)]*)\))?\s*(?:[_ ]?\d+)?$')
# the bracketed part of a variant key must be about the file (not e.g. 'RAW_FILE_NAME(sex)')
_RAW_KEY_QUALIFIER_RE = re.compile(r'raw|file|data|name|mz|pos|neg|hilic|c18|c8|c30|rp|ion|mode|source|lc|gc|ms|\d|^\s*$')
_EXACT_RAW_KEYS = {"raw_file_name", "rawfilename", "raw_file", "datafile name", "raw files", "raw file name",
                   "rawfile name", "rpm_filename", "zhp_filename", "data file", "raw_data"}


def is_variant_raw_file_key(key):
    key = str(key).lower().strip()
    if key in _EXACT_RAW_KEYS:
        return False
    m = _RAW_KEY_RE.match(key)
    if not m:
        return False
    return m.group(1) is None or bool(_RAW_KEY_QUALIFIER_RE.search(m.group(1)))


def _hint_text_indicates(text):
    """(polarity, chromatography) stated by a raw-file key qualifier or a polarity /
    chromatography factor value ('pos', 'HILIC-neg', 'ESI (+)', 'cation', ...)."""
    t = str(text).lower().replace('(', ' ').replace(')', ' ')
    t = re.sub(r'\bcation\b', 'positive', re.sub(r'\banion\b', 'negative', t))
    return _path_indicates(t)


def _raw_key_contradicts_analysis(key, hints):
    m = _RAW_KEY_RE.match(str(key).lower().strip())
    if not m or not m.group(1):
        return False
    polarity, chrom = _hint_text_indicates(m.group(1))
    return ((polarity is not None and hints.get('polarity') is not None and polarity != hints['polarity'])
            or (chrom is not None and hints.get('chrom') is not None and chrom != hints['chrom']))


_POLARITY_FACTOR_KEYS = {'polarity', 'ion mode', 'ionization mode', 'ionization', 'ion polarity', 'ms mode', 'esi mode',
                         'ms polarity', 'ionisation mode', 'mode', 'ms method'}
_CHROM_FACTOR_KEYS = {'chromatography', 'column', 'column type', 'lc method', 'separation', 'platform', 'mode',
                      'ms method', 'chromatography type'}


def _samples_of_other_analyses(df_additional, df_rest, hints):
    """Sample IDs whose polarity / chromatography factor (one unambiguous value per
    sample) contradicts this analysis' ION_MODE / chromatography."""
    rows = [df_additional[['filename', 'Key', 'Value']]]
    if len(df_rest):
        rows.append(df_rest[['Local Sample ID', 'Factor Name', 'Factor Value']].set_axis(['filename', 'Key', 'Value'], axis=1))
    kv = pd.concat(rows, ignore_index=True)
    kv['_k'] = kv['Key'].map(_norm_key)
    kv = kv[kv['_k'].isin(_POLARITY_FACTOR_KEYS | _CHROM_FACTOR_KEYS) & (kv['Value'].astype(str).str.len() <= 40)]
    out = set()
    for sample, g in kv.groupby('filename'):
        pols, chroms = set(), set()
        for k, v in zip(g['_k'], g['Value']):
            polarity, chrom = _hint_text_indicates(v)
            if k in _POLARITY_FACTOR_KEYS and polarity:
                pols.add(polarity)
            if k in _CHROM_FACTOR_KEYS and chrom:
                chroms.add(chrom)
        if (len(pols) == 1 and hints.get('polarity') and next(iter(pols)) != hints['polarity']) or \
                (len(chroms) == 1 and hints.get('chrom') and next(iter(chroms)) != hints['chrom']):
            out.add(sample)
    return out


def get_rawFile_names(df, key_vars, new_col):
    # Create a pattern to match key_vars and their numbered versions
    pattern = r'(' + '|'.join(key_vars) + r')(_\d+)?'

    # Extract the base key name without the number
    df['base_key'] = df['Key'].str.extract(pattern, expand=False)[0]

    # all raw file names of a sample (a sample may list several, e.g. one per
    # polarity), not only the last one
    filename_mapping = (df[df['base_key'].notna() & df['Value'].notna()]
                        .groupby('script_id')['Value']
                        .agg(lambda s: ';'.join(dict.fromkeys(str(v) for v in s))).to_dict())

    # Apply the mapping to create the new filename column
    df[new_col] = df['script_id'].map(filename_mapping)

    # Handle numeric conversion for specific columns
    if new_col in ['Latitude', 'Longitude']:
        df[new_col] = pd.to_numeric(df[new_col], errors='coerce')
    else:
        df[new_col] = df[new_col].fillna(value='NA')

    # Drop the temporary 'base_key' column
    df.drop('base_key', axis=1, inplace=True)
    return df



def map_potential_matches(df, target_df, target_column, source_columns):
    
    # Temporary storage for potential matches
    potential_matches = {}

    # Preprocess to get unique values from 'filename_base_wo_extension'
    unique_values = df['filename_base_wo_extension'].unique()

    # Identify potential matches for each unique value with a progress bar
    # plain list instead of DataFrame.iterrows (same result; iterrows made this take
    # hours for studies with thousands of files and samples)
    targets = target_df[target_column].tolist()
    for value in tqdm(unique_values, desc="Mapping potential matches", mininterval=30):
        matches = [t for t in targets if value in t or t in value]
        
        # Store unique matches only, mapping back to all rows with this value
        potential_matches[value] = set(matches)

    return potential_matches



def analyze_for_contradictions(potential_matches):
    # Reverse mapping to see which df values map to the same raw_file_name_df value
    reverse_mapping = {}
    for df_value, matches in potential_matches.items():
        for match in matches:
            if match in reverse_mapping:
                reverse_mapping[match].add(df_value)
            else:
                reverse_mapping[match] = {df_value}
    
    # Identify contradictions
    contradictions = {raw_value: df_values for raw_value, df_values in reverse_mapping.items() if len(df_values) > 1}
    
    return contradictions




def filter_unresolved_contradictions(potential_matches):
    
    # Step 1: Flatten all values into a single list and count occurrences
    all_values = [value for values in potential_matches.values() for value in values]
    
    value_counts = Counter(all_values)
    
    # Step 2: Identify values that occur more than once
    values_to_filter = {value for value, count in value_counts.items() if count > 1}
    
    # Step 3: Set the value of keys to an empty set if they contain any value to be filtered
    for key in potential_matches:
        if any(value in values_to_filter for value in potential_matches[key]):
            potential_matches[key] = set()
    
    return potential_matches

def resolve_contradictions(contradictions, potential_matches):
    for raw_value, conflicting_values in contradictions.items():
        better_match = None
        highest_percentage = 0
        
        raw_value_length = len(raw_value)
        
        for df_value in conflicting_values:
            # Find the percentage of the raw_value that matches the df_value
            # Assuming df_value represents the key in potential_matches that relates to raw_value
            percentage = len(df_value) / raw_value_length
            
            if percentage > highest_percentage:
                highest_percentage = percentage
                better_match = df_value
        
        # Now, keep only the better match in potential_matches for the conflicting raw_value
        if better_match:
            for key in list(potential_matches):  # Convert to list to avoid runtime error during iteration
                match_set = potential_matches[key]
                if raw_value in match_set and key != better_match:
                    match_set.remove(raw_value)

    return potential_matches

def create_dataframe_from_SUBJECT_SAMPLE_FACTORS(data, rest_response, raw_file_name_df=None, path_to_csvs='translation_sheets', attempt_taxa_extraction_from_inner = False, add_to_cols = [], **kwargs):
    
    ontology_table = kwargs['ontology_table']
    
    data_list = []
    for id, item in enumerate(data):
        subject_id = item.get('Subject ID', 'NA')
        sample_id = item.get('Sample ID', 'NA')
        factors = item.get('Factors', {})
        additional_data = item.get('Additional sample data', {})

        for key, value in factors.items():
            data_list.append([id, subject_id, sample_id, key, value])
        for key, value in additional_data.items():
            data_list.append([id, subject_id, sample_id, key, value])
    df_additional = pd.DataFrame(data_list, columns=['script_id', 'SubjectIdentifierAsRecorded', 'filename', 'Key', 'Value'])
    df_additional['Key'] = df_additional['Key'].str.lower()

    print(f"Rest length is {len(rest_response)}")
    df_rest = pd.DataFrame(columns=['ID', 'Study ID', 'MB Sample ID', 'Local Sample ID', 'Factor Name', 'Factor Value', 'raw_data_source'])
    if len(rest_response) > 0:
        data_list = []
        for id, item in rest_response.items():
            study_id = item.get('study_id', 'NA')
            local_sample_id = item.get('local_sample_id', 'NA')
            mb_sample_id = item.get('mb_sample_id', 'NA')
            sample_source = item.get('sample_source', 'NA')  # Corrected to sample_source
            factors = item.get('factors', "")
            raw_data = item.get('raw_data', "")

            # Adding sample source
            data_list.append([id, study_id, mb_sample_id, local_sample_id, "Sample Source", sample_source])

            # Splitting factors string into key-value pairs
            for factor in factors.split(' | '):
                key_value = factor.split(':')
                if len(key_value) == 2:
                    data_list.append([id, study_id, mb_sample_id, local_sample_id,  key_value[0], key_value[1]])
            
            # Adding raw_data if present
            if raw_data:
                data_list.append([id, study_id, mb_sample_id, local_sample_id, "raw_data", raw_data])
            
        df_rest = pd.DataFrame(data_list, columns=['ID', 'Study ID', 'MB Sample ID', 'Local Sample ID', 'Factor Name', 'Factor Value'])
        df_rest['raw_data_source'] = "Rest"


    if len(df_rest) == 0 and len(df_additional):
        # no REST factors at all: list the mwTab samples so that their Sample IDs can
        # still be matched to raw file names (this used to fail with KeyError 'USI')
        _ids = df_additional['filename'].drop_duplicates()
        df_rest = pd.DataFrame({'ID': None, 'Study ID': None, 'MB Sample ID': None, 'Local Sample ID': _ids.values,
                                'Factor Name': 'Sample ID', 'Factor Value': _ids.values, 'raw_data_source': 'mwTab'})

    #raw data name might be in the Sample ID, Factors
    expected_raw_file_keys = ["raw_file_name",  "rawfilename", "raw_file",
                              "datafile name", "raw files", "raw file name", 
                              "rawfile name", "rpm_filename", "zhp_filename",
                              "data file", "raw_data"]

    # Samples whose own polarity / chromatography factor says they belong to another
    # analysis of a multi-analysis study are left out of this analysis
    hints = kwargs.get('analysis_hints') or {}
    other_analysis_samples = set()
    if hints.get('multi'):
        other_analysis_samples = _samples_of_other_analyses(df_additional, df_rest, hints)
        if other_analysis_samples:
            print(f"{len(other_analysis_samples)} samples belong to another analysis by their polarity/chromatography factor")
            df_additional = df_additional[~df_additional['filename'].isin(other_analysis_samples)]
            df_rest = df_rest[~df_rest['Local Sample ID'].isin(other_analysis_samples)].reset_index(drop=True)

    _is_exact_raw_key = df_additional['Key'].isin(expected_raw_file_keys)
    # variant raw-file keys: RAW_FILE_NAME(Raw file pos), RAW_FILE_NAME_1 (a repeated
    # RAW_FILE_NAME key), RAW_FILE_NAME(mzXML file name), ...
    _is_variant_raw_key = ~_is_exact_raw_key & df_additional['Key'].map(is_variant_raw_file_key)
    df_filtered = df_additional[_is_exact_raw_key | _is_variant_raw_key]   # after this we only have raw data

    # Assigning unique script_id based on Local Sample ID
    df_rest['script_id'] = df_rest.groupby('Local Sample ID').ngroup()

    # Appending rows from `df` to `df_rest` where raw_data is missing
    samples_with_rest_raw = set(df_rest.loc[(df_rest['Factor Name'] == 'raw_data') & df_rest['Factor Value'].map(bool), 'Local Sample ID'])
    # samples with names under variant keys get all their mwTab raw file names (also
    # the plain RAW_FILE_NAME), so that none of a sample's files is left out
    samples_with_variant = set(df_filtered.loc[_is_variant_raw_key.reindex(df_filtered.index, fill_value=False), 'filename'])
    new_rows = []
    n_key_hint_skipped = 0
    for _, row in df_filtered.iterrows():
        is_variant = row['Key'] not in expected_raw_file_keys
        if is_variant and hints.get('multi') and _raw_key_contradicts_analysis(row['Key'], hints):
            n_key_hint_skipped += 1
            continue
        if is_variant or row['filename'] not in samples_with_rest_raw or row['filename'] in samples_with_variant:
            if len(str(row['Value'])) > 1:
                new_row = {
                    'ID': None,  # Assuming ID needs to be set or left as None
                    'Study ID': None,  # Assuming Study ID needs to be set or left as None
                    'MB Sample ID': None,  # Assuming MB Sample ID needs to be set or left as None
                    'Local Sample ID': row['filename'],
                    'Factor Name': row['Key'],
                    'Factor Value': row['Value'],
                    'raw_data_source': 'mwTab'
                }
                new_rows.append(new_row)
    if n_key_hint_skipped:
        print(f"{n_key_hint_skipped} raw file names left out: their key names another polarity/chromatography")
    print(f"Adding {len(new_rows)} raw file names from the mwTab")
    if new_rows:
        df_rest = pd.concat([df_rest, pd.DataFrame(new_rows)], ignore_index=True)
        # Re-number by Local Sample ID: mwTab rows carry their list position as id, which
        # does not match the REST numbering and would attach their raw file to another sample
        df_rest['script_id'] = df_rest.groupby('Local Sample ID').ngroup()

    df_rest['SubjectIdentifierAsRecorded'] = ""

    # Renaming columns to match the df DataFrame's column names
    df = df_rest.rename(columns={
        'Local Sample ID': 'filename',
        'Factor Name': 'Key',
        'Factor Value': 'Value'
    })
    # 'filename' holds the MWB Sample ID here but is overwritten with raw file paths by
    # some of the matching branches below; keep the Sample ID so the sample-level
    # metadata can still be joined back per sample afterwards
    df['_sample_id'] = df['filename']


    # Sample IDs are tried as raw file names when there is no raw file key; names
    # under variant raw-file keys do not replace that fallback
    only_variant_raw = (not any(key in df['Key'].values for key in expected_raw_file_keys)
                        and df['Key'].map(is_variant_raw_file_key).any())
    if only_variant_raw and not any("file" in key for key in df.loc[~df['Key'].map(is_variant_raw_file_key), 'Key'].values):
        df_new_row = df.drop_duplicates(['script_id', 'SubjectIdentifierAsRecorded', 'filename']).copy()
        df_new_row['Key'] = 'raw_file_name'
        df_new_row['Value'] = df_new_row['filename']
        df = pd.concat([df, df_new_row], ignore_index=True)

    if not any(key in df['Key'].values for key in expected_raw_file_keys):
        print('no raw file key found.')

        if any("file" in key for key in df['Key'].values):
            print("Found data files with unknown name")
            df['Key'] = df['Key'].apply(lambda x: "raw_file_name" if "file" in x else x)


        if not any(key in df['Key'].values for key in expected_raw_file_keys):
            #if we have no raw file key try to use the SampleID as raw_file_name
            df_new_row = df.drop_duplicates(['script_id', 'SubjectIdentifierAsRecorded', 'filename']).copy()
            df_new_row['Key'] = 'raw_file_name'
            df_new_row['Value'] = df_new_row['filename']
            df = pd.concat([df, df_new_row], ignore_index=True)

    if any(key in df['Key'].values for key in expected_raw_file_keys):
        print('raw file key found.')
        #attempt to get raw files from associated keys
        raw_file_name_df = raw_file_name_df.rename(columns={'filename': 'filename_raw_path'})


        # ---------- BEGIN: keep only preferred uploader file per stem (mzML > mzXML > others) ----------
        # This modifies raw_file_name_df so downstream matching will prefer mzML/mzXML automatically.

        RAWFILE_PREF_DEBUG_STEM = ""  # e.g. "p10_1_a6_qing_hao__dist" ; empty disables debug prints

        def _basename_lower(x):
            s = "" if pd.isna(x) else str(x)
            s = s.replace("\\", "/").split("/")[-1].lower()
            return s

        def _strip_gz(name_lower: str) -> str:
            return re.sub(r"\.gz$", "", name_lower)

        tmp = raw_file_name_df.copy()
        tmp["_base_lower"] = tmp["filename_base"].map(_basename_lower)
        tmp["_no_gz"] = tmp["_base_lower"].map(_strip_gz)

        # stem = filename without final extension (after optional .gz removal)
        tmp["_stem"] = tmp["_no_gz"].str.rsplit(".", n=1).str[0]

        # preference rank
        tmp["_rank"] = 2
        tmp.loc[tmp["_no_gz"].str.endswith(".mzml"), "_rank"] = 0
        tmp.loc[tmp["_no_gz"].str.endswith(".mzxml"), "_rank"] = 1

        before_n = len(tmp)

        # keep one row per stem: best rank, tie-break by name for determinism
        tmp = tmp.sort_values(["_stem", "_rank", "_base_lower"]).drop_duplicates("_stem", keep="first")

        after_n = len(tmp)

        # optional targeted debug (won't spam)
        if RAWFILE_PREF_DEBUG_STEM:
            cand = raw_file_name_df.copy()
            cand["_base_lower"] = cand["filename_base"].map(_basename_lower)
            cand["_no_gz"] = cand["_base_lower"].map(_strip_gz)
            cand["_stem"] = cand["_no_gz"].str.rsplit(".", n=1).str[0]
            cand["_rank"] = 2
            cand.loc[cand["_no_gz"].str.endswith(".mzml"), "_rank"] = 0
            cand.loc[cand["_no_gz"].str.endswith(".mzxml"), "_rank"] = 1

            cand = cand.loc[cand["_stem"] == RAWFILE_PREF_DEBUG_STEM, ["filename_base", "filename_raw_path", "_rank"]] \
                    .sort_values(["_rank", "filename_base"])

            print(f"\n[rawfile pref DEBUG] stem={RAWFILE_PREF_DEBUG_STEM}")
            if len(cand) == 0:
                print("No uploader candidates found for this stem (stem mismatch).")
            else:
                print(cand.to_string(index=False))
                chosen = tmp.loc[tmp["_stem"] == RAWFILE_PREF_DEBUG_STEM, ["filename_base", "filename_raw_path", "_rank"]]
                print("\n[rawfile pref DEBUG] chosen:")
                print(chosen.to_string(index=False))

        # write back cleaned uploader df
        raw_file_name_df = tmp.drop(columns=["_base_lower", "_no_gz", "_stem", "_rank"])

        # tiny summary (one line only)
        print(f"[rawfile pref] uploader list filtered {before_n} -> {after_n} using mzML > mzXML > others")
        # ---------- END: keep only preferred uploader file per stem ----------


        df = get_rawFile_names(df, key_vars=expected_raw_file_keys, new_col="filename_raw")

        #explode cases where we have multiple files by separators
        if raw_file_name_df['filename_base'].str.contains(" ").any():
            split_by_this = r'; ?|,'
        else:
            split_by_this = r'; ?|,| '

        df = df.assign(filename_raw=df['filename_raw'].str.split(split_by_this).apply(lambda x: list(set(x)))).explode('filename_raw')

        df.reset_index(drop=True, inplace=True)


        # --- Prefer mzML > mzXML > others BEFORE matching/merging ---

        # Set to {} to silence stem-specific debug
        DEBUG_STEMS = set()  # e.g. {"p10_1_a6_qing_hao__dist"}, lowercase stem WITHOUT extension

        def _basename_lower(s: str) -> str:
            # remove folder parts (handles / and \), lowercase
            return str(s).replace("\\", "/").split("/")[-1].lower()

        # build preferred filename_base per stem from uploader table
        raw_file_name_df['_base_lower'] = raw_file_name_df['filename_base'].map(_basename_lower)

        raw_file_name_df['_stem'] = raw_file_name_df['_base_lower'].str.rsplit('.', n=1).str[0]
        raw_file_name_df['_ext']  = raw_file_name_df['_base_lower'].str.rsplit('.', n=1).str[-1].fillna('')

        raw_file_name_df['_rank'] = (
            raw_file_name_df['_ext']
            .map({'mzml': 0, 'mzxml': 1})
            .fillna(2)
            .astype(int)
        )

        preferred_raw_by_stem = (
            raw_file_name_df
            .sort_values(['_stem', '_rank', '_base_lower'])
            .drop_duplicates('_stem', keep='first')
            [['_stem', 'filename_base']]
            .rename(columns={'filename_base': 'preferred_filename_base'})
        )

        # compute stem from metadata-provided filename_raw
        df['filename_raw_orig'] = df['filename_raw']

        df['_raw_base_lower'] = df['filename_raw'].map(_basename_lower)
        df['_raw_stem'] = df['_raw_base_lower'].str.rsplit('.', n=1).str[0]
        df['_raw_ext']  = df['_raw_base_lower'].str.rsplit('.', n=1).str[-1].fillna('')

        df = df.merge(preferred_raw_by_stem, left_on='_raw_stem', right_on='_stem', how='left')

        # replace filename_raw with preferred uploader filename_base when available
        has_pref = df['preferred_filename_base'].notna()
        df.loc[has_pref, 'filename_raw'] = df.loc[has_pref, 'preferred_filename_base']

        # compact debug
        changed = (df['filename_raw'] != df['filename_raw_orig'])
        upgraded = changed & (df['_raw_ext'] == 'raw')
        print(f"[rawfile preference] upgraded {int(upgraded.sum())} '.raw' -> preferred (mzML/mzXML/other); "
            f"total changed = {int(changed.sum())}")

        # stem-specific debug only
        if DEBUG_STEMS:
            dbg_mask = df['_raw_stem'].isin(DEBUG_STEMS)
            if dbg_mask.any():
                print("\n[rawfile preference DEBUG] affected rows (orig -> new):")
                print(df.loc[dbg_mask, ['filename_raw_orig', 'filename_raw', '_raw_stem', '_raw_ext', 'preferred_filename_base']]
                        .drop_duplicates()
                        .to_string(index=False))

                # show uploader candidates for these stems (small, ranked)
                cand = raw_file_name_df.loc[raw_file_name_df['_stem'].isin(DEBUG_STEMS),
                                            ['filename_raw_path', 'filename_base', '_ext', '_rank']] \
                                    .sort_values(['_rank', 'filename_base'])
                print("\n[rawfile preference DEBUG] uploader candidates for stem(s):")
                print(cand.to_string(index=False))

        # cleanup temp cols (keep filename_raw_orig if you want)
        df.drop(columns=['_raw_base_lower', '_raw_stem', '_raw_ext', '_stem', 'preferred_filename_base'], inplace=True, errors='ignore')
        raw_file_name_df.drop(columns=['_base_lower', '_stem', '_ext', '_rank'], inplace=True, errors='ignore')

        raw_file_name_df['filename_base_lower'] = raw_file_name_df['filename_base'].str.lower()


        df['filename_raw_lower'] = df['filename_raw'].str.lower()
        df['filename_lower'] = df['filename'].str.lower()

        if bool(set(df['filename_raw']) & set(raw_file_name_df['filename_base'])):
            print('found filename_raw to filename_base match')
            df = df.merge(raw_file_name_df, left_on='filename_raw', right_on='filename_base', how='left')
            df.drop(['filename_raw_lower', 'filename_lower', 'filename_base_lower'], axis=1, inplace=True)
        elif bool(set(df['filename_raw_lower']) & set(raw_file_name_df['filename_base_lower'])):
            df = df.merge(raw_file_name_df, left_on='filename_raw_lower', right_on='filename_base_lower', how='left')
            df.drop(['filename_raw_lower', 'filename_lower', 'filename_base_lower'], axis=1, inplace=True)
        elif bool(set(df['filename_lower']) & set(raw_file_name_df['filename_base_lower'])):
            df = df.merge(raw_file_name_df, left_on='filename_lower', right_on='filename_base_lower', how='left')
            df.drop(['filename_raw_lower', 'filename_lower', 'filename_base_lower'], axis=1, inplace=True)
    
        else:
            df['filename_base_wo_extension'] = df['filename_raw_lower'].str.split('.').str[0]
            raw_file_name_df['filename_base_wo_extension'] = raw_file_name_df['filename_base'].str.split('.').str[0].str.lower()
            raw_file_name_df['filename_raw'] = raw_file_name_df['filename_raw_path'].str.lower()
            raw_file_name_df['filename_raw_lower'] = raw_file_name_df['filename_raw'].str.lower()

            match_df = 'filename_raw_lower'
            match_raw_file_name_df = 'filename_base_wo_extension'
            potential_matches = set(df[match_df]).intersection(set(raw_file_name_df[match_raw_file_name_df]))

            if not potential_matches:
                match_df = 'filename_raw_lower'
                match_raw_file_name_df = 'filename_raw_lower'
                potential_matches = set(df[match_df]).intersection(set(raw_file_name_df[match_raw_file_name_df]))

            if not potential_matches:
                match_df = 'filename_base_wo_extension'
                match_raw_file_name_df = 'filename_base_wo_extension'
                potential_matches = set(df[match_df]).intersection(set(raw_file_name_df[match_raw_file_name_df]))

            if potential_matches:
                df = df[['script_id', 'SubjectIdentifierAsRecorded', 'filename', '_sample_id', 'Key', 'Value', match_df]]
                df = df.merge(raw_file_name_df, left_on=match_df, right_on=match_raw_file_name_df, how='left')
            else:
                # here we try to find raw file names which are represented as substrings oif actual file names in the uploader provided raw file names

                print(f"use from match_df: {match_df}")
                print(f"use from raw_file_name_df: {match_raw_file_name_df}")

                print('getting substring matches,..')
                potential_matches = map_potential_matches(df, raw_file_name_df, target_column='filename_base_wo_extension', source_columns = ['filename_base_wo_extension', ''])
                print('analyzing contradictions,..')
                contradictions = analyze_for_contradictions(potential_matches)
                print('resolving contradictions,..')
                potential_matches = resolve_contradictions(contradictions, potential_matches)
                print('filtering contradictions,..')
                potential_matches = filter_unresolved_contradictions(potential_matches)
                
                #create dict to assign matches
                match_indices = {}

                for df_value, matches in potential_matches.items():
                    if matches:
                        for match_value in matches:  
                            match_index = raw_file_name_df.index[raw_file_name_df['filename_base_wo_extension'] == match_value].tolist()[0] 
                            if df_value not in match_indices:
                                match_indices[df_value] = [match_index]
                            else:
                                match_indices[df_value].append(match_index)

                
                df['match_index'] = df['filename_base_wo_extension'].map(match_indices)
                df = df.explode('match_index')
                
                df = pd.merge(df, raw_file_name_df, left_on='match_index', right_index=True, how='left', suffixes=('', '_matched'))

                df['filename_raw'] = df['filename_raw_path']
                df['filename'] = df['filename_raw_path']
                df.drop(columns=['match_index'], inplace=True)

            if 'filename_base_wo_extension' in df.columns:
                df['_stem'] = df['filename_base_wo_extension'].astype(str).str.lower()
                missing_stem = df['_stem'].isin(['', 'nan', 'none']) | df['_stem'].isna()
                if missing_stem.any():
                    df.loc[missing_stem, '_stem'] = (
                        df.loc[missing_stem, 'filename_raw']
                        .astype(str).str.lower().str.split('.').str[0]
                    )
            else:
                df['_stem'] = df['filename_raw'].astype(str).str.lower().str.split('.').str[0]

            df = df.merge(preferred_raw_by_stem, on='_stem', how='left', suffixes=('', '_preferred'))

            # Guard the '_preferred' overrides: the preferred_raw_by_stem merged here
            # only carries 'preferred_filename_base', so the 'filename_raw_path_preferred'
            # / 'filename_base_preferred' columns the original code referenced never
            # exist -> an unconditional reference raised KeyError and dropped the whole
            # study (this is what silently lost ST001000 and other studies routed
            # through the fuzzy-match branch). This early override is in any case
            # superseded by the full mzML>mzXML>vendor preference block right below,
            # so no-op'ing it here when the columns are absent is safe.
            if 'filename_raw_path_preferred' in df.columns:
                if 'filename_raw_path' in df.columns:
                    df['filename_raw_path'] = df['filename_raw_path_preferred'].combine_first(df['filename_raw_path'])
                else:
                    df['filename_raw_path'] = df['filename_raw_path_preferred']

            if 'filename_base_preferred' in df.columns:
                if 'filename_base' in df.columns:
                    df['filename_base'] = df['filename_base_preferred'].combine_first(df['filename_base'])
                else:
                    df['filename_base'] = df['filename_base_preferred']

            # Drop every helper column this early-override block introduced (both the
            # '_preferred'-suffixed names AND the un-suffixed 'preferred_*' columns that
            # the merged preferred_raw_by_stem carries). Leaving 'preferred_filename_base'
            # behind here made the *next* preferred_raw_by_stem merge below collide and
            # rename to _x/_y, producing a 'preferred_filename_base' KeyError that again
            # dropped the whole study.
            df.drop(columns=['_stem', 'filename_base_preferred', 'filename_raw_path_preferred',
                             'preferred_filename_base', 'preferred_filename_raw_path'],
                    inplace=True, errors='ignore')


            # --- PREFER mzML > mzXML > others (and print debug for why it might not upgrade) ---
            print("Applying raw-file extension preference: mzML > mzXML > others")

            # helper: strip any folder from a path-like string (handles / and \)
            raw_file_name_df['_base_lower'] = (
                raw_file_name_df['filename_base']
                .astype(str)
                .str.replace(r'^.*[\\/]', '', regex=True)
                .str.lower()
            )

            raw_file_name_df['_stem'] = raw_file_name_df['_base_lower'].str.rsplit('.', n=1).str[0]
            raw_file_name_df['_ext']  = raw_file_name_df['_base_lower'].str.rsplit('.', n=1).str[-1].fillna('')

            raw_file_name_df['_rank'] = (
                raw_file_name_df['_ext']
                .map({'mzml': 0, 'mzxml': 1})
                .fillna(2)
                .astype(int)
            )

            preferred_raw_by_stem = (
                raw_file_name_df
                .sort_values(['_stem', '_rank', '_base_lower'])
                .drop_duplicates('_stem', keep='first')
                [['_stem', 'filename_base', 'filename_raw_path']]
                .rename(columns={
                    'filename_base': 'preferred_filename_base',
                    'filename_raw_path': 'preferred_filename_raw_path'
                })
            )

            # compute df stem from whatever you currently have
            if 'filename_base' in df.columns and df['filename_base'].notna().any():
                df['_df_base_lower'] = (
                    df['filename_base'].astype(str)
                    .str.replace(r'^.*[\\/]', '', regex=True)
                    .str.lower()
                )
            elif 'filename_raw' in df.columns and df['filename_raw'].notna().any():
                df['_df_base_lower'] = (
                    df['filename_raw'].astype(str)
                    .str.replace(r'^.*[\\/]', '', regex=True)
                    .str.lower()
                )
            else:
                df['_df_base_lower'] = (
                    df['filename'].astype(str)
                    .str.replace(r'^.*[\\/]', '', regex=True)
                    .str.lower()
                )

            df['_stem'] = df['_df_base_lower'].str.rsplit('.', n=1).str[0]

            # merge preferred candidate
            df = df.merge(preferred_raw_by_stem, on='_stem', how='left')

            # DEBUG: your specific problematic stem
            debug_stem = ""  # e.g. "p10_1_a6_qing_hao__dist"
            if debug_stem and (df['_stem'] == debug_stem).any():
                print("\n[DEBUG] stem =", debug_stem)

                cand_rows = raw_file_name_df.loc[
                    raw_file_name_df['_stem'] == debug_stem,
                    ['filename_raw_path', 'filename_base', '_ext', '_rank']
                ].sort_values(['_rank', 'filename_base'])

                if len(cand_rows) == 0:
                    print("[DEBUG] No uploader candidates found for this stem. That means your stem extraction didn't match.")
                    print("[DEBUG] Example df bases for this stem:")
                    print(df.loc[df['_stem'] == debug_stem, ['_df_base_lower']].head(10).to_string(index=False))
                else:
                    print("[DEBUG] Uploader candidates for this stem (sorted by preference):")
                    print(cand_rows.to_string(index=False))

                cols = [c for c in [
                    'script_id', 'SubjectIdentifierAsRecorded',
                    'filename', 'filename_raw', 'filename_base', 'filename_raw_path',
                    'preferred_filename_base', 'preferred_filename_raw_path'
                ] if c in df.columns]

                print("[DEBUG] Row(s) BEFORE override:")
                print(df.loc[df['_stem'] == debug_stem, cols].to_string(index=False))

            # apply override where we have a preferred file
            has_pref = df['preferred_filename_base'].notna()

            # these two are your merged uploader columns (if present)
            if 'filename_base' in df.columns:
                df.loc[has_pref, 'filename_base'] = df.loc[has_pref, 'preferred_filename_base']
            if 'filename_raw_path' in df.columns:
                df.loc[has_pref, 'filename_raw_path'] = df.loc[has_pref, 'preferred_filename_raw_path']

            # IMPORTANT: also overwrite the columns that often “carry through” to outputs/deduping
            # If in your pipeline `filename` is intended to represent the associated raw file, set it to the preferred *path*.
            if 'filename_raw' in df.columns:
                df.loc[has_pref, 'filename_raw'] = df.loc[has_pref, 'preferred_filename_base']

            if 'filename' in df.columns:
                # If `filename` is *metadata filename*, do NOT overwrite it.
                # In your shown code, you later treat filename as raw-associated in several places and dedupe on it,
                # so overwriting it is likely what you want.
                df.loc[has_pref, 'filename'] = df.loc[has_pref, 'preferred_filename_raw_path']

            if (df['_stem'] == debug_stem).any():
                cols = [c for c in [
                    'script_id', 'SubjectIdentifierAsRecorded',
                    'filename', 'filename_raw', 'filename_base', 'filename_raw_path'
                ] if c in df.columns]
                print("[DEBUG] Row(s) AFTER override:")
                print(df.loc[df['_stem'] == debug_stem, cols].to_string(index=False))
                print()

            # cleanup
            df.drop(columns=[
                'preferred_filename_base', 'preferred_filename_raw_path',
                '_df_base_lower', '_stem'
            ], inplace=True, errors='ignore')

            # also cleanup temp cols from raw_file_name_df to avoid polluting later logic
            raw_file_name_df.drop(columns=['_base_lower', '_stem', '_ext', '_rank'], inplace=True, errors='ignore')



            df['filename_raw_lower'] = df['filename_base']
            df = df.drop(columns=['filename_base'])



    if 'USI' not in df.columns:
        raise ValueError("No raw data could be associated with metadata. Stopping!")
    df = df[pd.notna(df['USI'])]
    df = df.drop_duplicates(subset=['filename', 'Key' , 'Value', 'USI'], keep='first').reset_index(drop=True)
    if len(df) == 0:
        raise ValueError("No raw data could be associated with metadata. Stopping!")
    

    # Group the sample-level metadata by Sample ID again (the final filename is taken
    # from filename_raw_path below, so this only affects the grouping in between)
    if '_sample_id' in df.columns:
        df['filename'] = df['_sample_id'].where(df['_sample_id'].notna(), df['filename'])
        df = df.drop(columns=['_sample_id'])

    # appending the additional data again to make sure we have all metadata
    try:
        df_update = pd.concat([df, df_additional], ignore_index=True)

        # This consolidates *sample-level metadata* per sample: group by 'filename'
        # (which at this stage is the MWB Sample ID) and, for each column, keep the
        # single agreed value or blank it to NaN if the sample's rows disagree.
        #
        # File-association columns are EXEMPT: a single sample legitimately maps to
        # several raw files (e.g. one measurement per chromatography/ionization
        # method, all listed under the same Sample ID). Those distinct file paths
        # are not a "discrepancy" to blank - blanking them destroys the file
        # association and silently drops the sample entirely (this was losing whole
        # multi-method studies like ST002470/ST002471). Exempting these columns is a
        # no-op for the common one-file-per-sample case (a single value is kept
        # either way), so it cannot change output for those studies.
        file_assoc_cols = {c for c in df_update.columns if 'filename' in c} | {'USI', 'size_mb'}

        def _agreed_value(s):
            # Empty strings are not a value: REST rows carry SubjectIdentifierAsRecorded
            # "" while the mwTab rows of the same sample carry the real Subject ID, and
            # counting "" as a second value blanked every Subject ID
            s = s.dropna()
            filled = s[s.astype(str) != '']
            if filled.nunique() == 1:
                return filled.iloc[0]
            if len(filled) == 0 and s.nunique() == 1:
                return s.iloc[0]
            return np.nan

        for col in df_update.columns:
            if col in ['filename', 'Key', 'Value'] or col in file_assoc_cols:
                continue
            # a sample with no value at all in this column (e.g. one listed only in the
            # mwTab additional data) stays NaN instead of raising on an empty group
            df_update[col] = df_update.groupby('filename')[col].transform(_agreed_value)
        df = df_update.copy()
    except Exception as e:
        print("Error during re-adding additional data:", e)


    df['Value'] = df['Value'].str.lower()
    df['SubjectIdentifierAsRecorded'] = df['SubjectIdentifierAsRecorded'].replace('-', '')
    try:
        df['SubjectIdentifierAsRecorded'] = _drop_label_subject_ids(df)
    except Exception as e:
        print(f"Could not check Subject IDs: {type(e).__name__}: {e}")
        df['SubjectIdentifierAsRecorded'] = ''
    df = get_key_info_into_outer(df, key_vars=["Sample Source"], new_col="MWB_sampleSource")
    # sex keys: the known ones plus other 'sex'/'gender' keys about the sampled subject;
    # numeric codes are decoded only where the key defines them
    _keyl = df['Key'].astype(str).str.lower()
    _sex_keys = ["gender", "sex", "gender (f/m)", "biological_sex"]
    _extra_sex = sorted(set(_keyl[_keyl.str.contains(r'\b(?:sex|gender)\b', regex=True)
                                  & ~_keyl.str.contains(_SEX_KEY_EXCLUDE_RE, regex=True)]) - set(_sex_keys))
    _is_sex = _keyl.isin(_sex_keys + _extra_sex)
    _decoded = [convert_sex_coded(v, k) for k, v in zip(df.loc[_is_sex, 'Key'], df.loc[_is_sex, 'Value'])]
    df.loc[_is_sex, 'Value'] = [d if d is not None else v for d, v in zip(_decoded, df.loc[_is_sex, 'Value'])]
    df = get_key_info_into_outer(df, key_vars=_sex_keys + _extra_sex, new_col="MWB_sex")
    # keys the patterns above miss because '_' / '.' are word characters for '\b'
    # (sex_participant, user_gender, mouse_sex): only fill samples without a value
    _keyn = _keyl.map(_norm_key)
    _sex_fb = (~_is_sex & _keyn.str.contains(r'\b(?:sex|gender)\b', regex=True)
               & ~_keyn.str.contains(_SEX_KEY_EXCLUDE_RE + r'|\bindex\b|\bcontact\b|\bchild\b', regex=True))
    if _sex_fb.any():
        _decoded = [convert_sex_coded(v, k) for k, v in zip(df.loc[_sex_fb, 'Key'], df.loc[_sex_fb, 'Value'])]
        _vals = [d if d is not None else convert_sex(v) for d, v in zip(_decoded, df.loc[_sex_fb, 'Value'])]
        df = _fill_from_rows(df, pd.Series(_vals, index=df.index[_sex_fb]), 'MWB_sex')
    df = get_age_years_into_outer(df, new_col="MWB_age")
    # age keys such as host_age, age_at_sample_collection, age_months, age_visit: a single
    # number (optionally with a unit) under a key about the sampled subject's age
    _age_fb = (~_keyl.str.contains(r'\bage\b', regex=True) & _keyn.str.contains(r'\bage\b', regex=True)
               & ~_keyn.str.contains(AGE_EXCLUDE_RE + r'|' + _AGE_FB_EXCLUDE_RE, regex=True))
    if _age_fb.any():
        _vals = [mwb_age_to_years(v, k.replace('_', ' ')) if _SINGLE_AGE_VALUE_RE.fullmatch(str(v).strip()) else None
                 for k, v in zip(df.loc[_age_fb, 'Key'], df.loc[_age_fb, 'Value'])]
        df = _fill_from_rows(df, pd.Series(_vals, index=df.index[_age_fb], dtype=object), 'MWB_age')
    smoking_keys = ["smoking_status", "smoking status", "smokstatus",
                    "smoking", "smoke", "tobacco use", "smoking history",
                    "cigarette smoking", "smoker", "smokingstatus",
                    "smoker status", "smoking_history", "subject_smoker",
                    "ever smoked", "ever_smoked", "current smoker", "current_smoker"]
    # yes/no answers only have a meaning together with the key ('Smoker: No' vs
    # 'Ever smoked: No'), so convert them here while the key is still known
    _is_smoking = df['Key'].astype(str).str.lower().isin(smoking_keys)
    _converted = [convert_smoking(v, key=k) for k, v in zip(df.loc[_is_smoking, 'Key'], df.loc[_is_smoking, 'Value'])]
    df.loc[_is_smoking, 'Value'] = [c if c != 'missing value' else v for c, v in zip(_converted, df.loc[_is_smoking, 'Value'])]
    df = get_key_info_into_outer(df, key_vars=smoking_keys, new_col="MWB_smoking")
    # further smoking keys (and the listed ones written with other separators); only
    # fill samples without a value
    _smoke_fb = ~_is_smoking & (_keyn.isin({_norm_key(k) for k in smoking_keys}) | _keyl.isin(_EXTRA_SMOKING_KEYS))
    if _smoke_fb.any():
        _vals = [convert_smoking_coded(v, k) or convert_smoking(v, key=k)
                 for k, v in zip(df.loc[_smoke_fb, 'Key'], df.loc[_smoke_fb, 'Value'])]
        df = _fill_from_rows(df, pd.Series(_vals, index=df.index[_smoke_fb]).replace('missing value', None), 'MWB_smoking')
    df = get_key_info_into_outer(df, key_vars=["bmi", "body mass index", "body_mass_index",
                                               "bmi (kg/m2)", "bmi_kg_m2"], new_col="MWB_bmi")
    df = get_key_info_into_outer(df, key_vars=["diet", "diet_type", "diet type", "diet group",
                                               "dietary group", "dietary pattern", "feeding regime",
                                               "diet_group"], new_col="MWB_diet")
    df = get_key_info_into_outer(df, key_vars=["collection_country", "collection country", "country", "site", "location", "country_of_origin"],
                                 new_col="Country")
    df = get_key_info_into_outer(df, key_vars=["latitude", "lat", "lat_dd", "decimal latitude", "latitude (decimal degrees)"], new_col="Latitude")
    df = get_key_info_into_outer(df, key_vars=["longitude", "lon", "long", "lon_dd", "decimal longitude", "longitude (decimal degrees)"], new_col="Longitude")

    # extra per-sample fields (depth, collection date, health, life stage); a problem
    # here must never cost the study its rows
    try:
        # depth / altitude
        _keyl = df['Key'].astype(str).str.lower()
        _is_depth = _keyl.str.contains(_DEPTH_KEY_RE, regex=True) & ~_keyl.str.contains(_DEPTH_KEY_EXCLUDE_RE, regex=True)
        _depth = pd.Series([depth_or_altitude_meters(v, k) for k, v in zip(df.loc[_is_depth, 'Key'], df.loc[_is_depth, 'Value'])],
                           index=df.index[_is_depth], dtype=object)
        df['MWB_depth'] = _depth
        df['MWB_depth'] = df.groupby('filename')['MWB_depth'].transform(lambda s: s.dropna().iloc[0] if s.dropna().nunique() == 1 else np.nan)

        # collection date: raw value, normalised (and rejected if not a real date) during validation
        df = get_key_info_into_outer(df, key_vars=sorted(_DATE_KEYS), new_col="MWB_collection_date")

        # health: 'healthy' keywords and exact DOID names under disease-like keys
        _is_health = _keyl.str.contains(_HEALTH_KEY_RE, regex=True)
        _vals = df['Value'].astype(str).str.strip().str.lower()
        _healthy = _is_health & (_vals.isin(_HEALTHY_VALUES) | (_keyl.str.contains(_DISEASE_KEY_RE, regex=True) & _vals.isin({'control', 'controls', 'normal', 'none'})))
        df['MWB_healthy'] = np.where(_healthy, 'healthy', None)
        df['MWB_healthy'] = df.groupby('filename')['MWB_healthy'].transform(lambda s: s.dropna().iloc[0] if s.dropna().nunique() == 1 else None)
        if DOID_LOOKUP:
            df['MWB_doid'] = _vals.where(_is_health).map(DOID_LOOKUP)
            df['MWB_doid'] = df.groupby('filename')['MWB_doid'].transform(lambda s: s.dropna().iloc[0] if s.dropna().nunique() == 1 else None)
        else:
            df['MWB_doid'] = None

        # life stage stated directly (newborn, adolescent, ...)
        _is_ls = _keyl.str.contains(_LIFESTAGE_KEY_RE, regex=True)
        df['MWB_lifestage'] = _vals.where(_is_ls).map(_LIFESTAGE_VALUES)
        df['MWB_lifestage'] = df.groupby('filename')['MWB_lifestage'].transform(lambda s: s.dropna().iloc[0] if s.dropna().nunique() == 1 else None)
    except Exception as e:
        print(f"Could not extract additional sample fields: {type(e).__name__}: {e}")
        for col in ['MWB_depth', 'MWB_collection_date', 'MWB_healthy', 'MWB_doid', 'MWB_lifestage']:
            if col in df.columns:
                df = df.drop(columns=[col])
    df['LatitudeandLongitude'] = df.apply(
        lambda x: f'{x.Latitude}|{x.Longitude}' if x.Latitude is not None and x.Longitude is not None and not np.isnan(
            x.Latitude) and not np.isnan(x.Longitude) else None, axis=1)
    
    df[['SampleType_inner', 'SampleTypeSub1_inner']] = df.apply(lambda x: get_blanks(x['Value']), axis=1,
                                                                result_type='expand')

    df = translate_MWB_to_REDU_from_csv(df, case='inner', path_to_csvs=path_to_csvs, add_to_cols=add_to_cols, ontology_table=ontology_table, allowedTerm_dict=allowedTerm_dict)

    
    df.loc[df['MWB_sampleSource'] == '_blank_', ['SampleType_inner', 'SampleTypeSub1_inner']] = ['blank_analysis', 'blank_analysis']
    df.loc[df['MWB_sampleSource'] == '_qc_', ['SampleType_inner', 'SampleTypeSub1_inner']] = ['blank_QC', 'blank_QC']
    df.loc[df['MWB_sampleSource'] == '_pooled_', ['SampleType_inner', 'SampleTypeSub1_inner']] = ['blank_QC', 'blank_QC']
    df.loc[(df['MWB_sampleSource'] == 'culture media') | 
           (df['MWB_sampleSource'] == 'bacterial media') | 
           (df['MWB_sampleSource'] == 'adipocyte media') |
           (df['MWB_sampleSource'] == 'media'), ['SampleType_inner', 'SampleTypeSub1_inner']] = ['blank_culturemedia', 'blank_culturemedia']
    df.loc[(df['MWB_sampleSource'] == 'bacterial cells'), ['SampleType_inner', 'SampleTypeSub1_inner']] = ['culture_bacterial', 'culture_bacterial']

    # QC pools, blanks and reference materials named by a sample's values (sample type
    # 'Study Pool', group 'QC', 'process blank', 'NIST SRM 1950', ...); only samples
    # that got no sample type above, and only one class per sample
    try:
        _qc = pd.Series([qc_sample_type(k, v) for k, v in zip(df['Key'], df['Value'])], index=df.index, dtype=object)
        _qc_per_sample = (pd.DataFrame({'filename': df['filename'], 'qc': _qc}).dropna()
                          .groupby('filename')['qc'].agg(lambda s: s.iloc[0] if s.nunique() == 1 else None).dropna())
        _typed = df.groupby('filename')['SampleType_inner'].transform(lambda s: s.notna().any())
        _set = ~_typed & df['filename'].isin(_qc_per_sample.index)
        if _set.any():
            _t = df.loc[_set, 'filename'].map(_qc_per_sample)
            df.loc[_set, 'SampleType_inner'] = _t
            df.loc[_set, 'SampleTypeSub1_inner'] = _t.map(lambda t: 'missing value' if t == 'reference material' else t)
            print(f"QC/blank/reference samples from their values: {df.loc[_set, 'filename'].nunique()}")
    except Exception as e:
        print(f"Could not detect QC/blank samples: {type(e).__name__}: {e}")
   
   #add environemntal columns
    
    df['MWB_envo_exact'] = envo_material_exact_name(df['MWB_sampleSource'], allowedTerm_dict).values
    df = df.drop(columns=['Key', 'Value', 'Longitude', 'Latitude', 'MWB_sampleSource'])
    df = df.drop_duplicates().reset_index(drop=True)

    df['filename'] = df['filename_raw_path']
    return df


def _section(mwtab, name):
    """A mwTab section as dict; sections without content come as [] (e.g. an empty
    CHROMATOGRAPHY section of a direct-infusion analysis), which used to crash."""
    section = mwtab.get(name, {})
    return section if isinstance(section, dict) else {}


_HUMAN_SUBJECT_TYPES = {'human', 'human clinical study', 'human subjects'}
# SUBJECT_TYPE -> (SampleType, SampleTypeSub1), used only when the species gave none
_SUBJECT_TYPE_SAMPLE_TYPE = {
    'synthetic sample': ('molecular_standards', None), 'synthetic': ('molecular_standards', None),
    'chemical': ('molecular_standards', None), 'metabolite standards': ('molecular_standards', None),
    'standards': ('molecular_standards', None), 'chemical standards': ('molecular_standards', None),
    'soil sample': ('environmental', None), 'water sample': ('environmental', None),
    'food item': ('food', None),
    'plant': ('plant', 'plant_NOS'), 'plants': ('plant', 'plant_NOS'),
    'bacteria': ('bacteria', None), 'bacterial cells': ('bacteria', None),
    'yeast': ('fungi', None), 'fungi': ('fungi', None), 'fungal': ('fungi', None),
    'human': ('animal', None), 'mammal': ('animal', None), 'animal': ('animal', None), 'mouse': ('animal', None),
    'rat': ('animal', None), 'fish': ('animal', None), 'insect': ('animal', None), 'invertebrate': ('animal', None),
    'bird': ('animal', None), 'amphibian': ('animal', None), 'human;animal': ('animal', None), 'sheep': ('animal', None),
    'human clinical study': ('animal', None), 'algae': ('algae', None),
}


_CULTURED_CELL_SUBJECT_TYPES = {'cultured cells', 'cells', 'cell culture', 'human cells', 'cell line', 'cell lines'}


def multi_species_sample_type(species, subject_type, allowedTerm_dict):
    """(SampleType, SampleTypeSub1) shared by every species of a multi-species
    SUBJECT_SPECIES, when each part is exactly an allowed NCBI name; else None.
    Mammals in a cultured-cell study are 'culture_mammalian'."""
    parts = [p.strip() for p in re.split(r'[;,/]| and ', str(species)) if p.strip()]
    if len(parts) < 2:
        return None
    types = set()
    for part in parts:
        term = _ncbi_name_lookup(allowedTerm_dict).get(part.lower())
        if term is None:
            return None
        info = get_taxonomy_info(term.split('|')[0])
        if not info or info[0] is None:
            return None
        lineage_mammal = info[0] == 'animal' and info[1] is None
        types.add((info[0], info[1], lineage_mammal))
    if len({(t, s1) for t, s1, _ in types}) != 1:
        return None
    sample_type, sub1, _ = next(iter(types))
    if sample_type == 'animal' and str(subject_type).strip().lower() in _CULTURED_CELL_SUBJECT_TYPES:
        if all(_is_mammal(p, allowedTerm_dict) for p in parts):
            return ('culture_mammalian', 'culture_mammalian')
        return None
    return (sample_type, sub1)


def _is_mammal(name, allowedTerm_dict):
    term = _ncbi_name_lookup(allowedTerm_dict).get(name.lower())
    if term is None:
        return False
    try:
        from REDU_conversion_functions import _get_ncbi_lineage
        lineage = _get_ncbi_lineage(term.split('|')[0])
    except Exception:
        return False
    return lineage is not None and 'mammalia' in lineage


def sample_type_from_subject_type(subject_type, sample_type=''):
    if not isinstance(subject_type, str):
        return None
    hit = _SUBJECT_TYPE_SAMPLE_TYPE.get(' '.join(subject_type.strip().lower().split()))
    # 'Synthetic sample' also covers in-vitro assays and mislabelled studies: standards
    # only when the sample type says so ('Synthetic Mixture', 'Standard phosphatidylcholines')
    if hit and hit[0] == 'molecular_standards' and not re.search(r'standard|synthetic', str(sample_type), re.I):
        return None
    return hit


_TXID_RE = re.compile(r'^\s*(?:ncbi\s*:\s*)?(?:txid|taxid|taxon\s*:)\s*(\d+)\s*$', re.I)
# trailing botanical author abbreviation: 'Lactuca sativa L.', 'Elaeis guineensis Jacq.'
_AUTHOR_RE = re.compile(r'^([A-Z][a-z]+ [a-z][a-z-]+) (?:[A-Z][a-z]{0,6}\.)$')


def study_species_variant_term(species, taxonomy_id, allowedTerm_dict):
    """Allowed NCBITaxonomy term for a study species/ID the plain lookup misses only
    because of its format; one exact name or ID, else None."""
    terms = allowedTerm_dict['NCBITaxonomy']['allowed_values']
    m = _TXID_RE.match(str(taxonomy_id))
    if m:
        hit = get_taxonomy_id_from_name__allowedTerms('', ncbi_id={m.group(1)}, allowedTerm_dict=allowedTerm_dict)
        if hit:
            return hit
    sp = str(species).strip()
    candidates = []
    paren = re.fullmatch(r'[^()]*\(([^()]+)\)\s*', sp)
    if paren:
        candidates.append(paren.group(1).strip())
    candidates.append(sp)
    for cand in candidates:
        a = _AUTHOR_RE.match(cand)
        for name in ([cand, a.group(1)] if a else [cand]):
            if cand is sp and name == sp:
                continue  # already tried by the plain lookup
            if re.fullmatch(r'[A-Z][a-z]+ [a-z][a-z-]+(?: [a-z][a-z-]+)?', name):
                hit = _ncbi_name_lookup(allowedTerm_dict).get(name.lower())
                if hit:
                    return hit
    return None


# number of MS analyses of the study being converted (set by MWB_to_REDU_study_wrapper)
_STUDY_N_MS_ANALYSES = 1


def _analysis_hints(df_outer):
    """Polarity ('positive'/'negative') and chromatography ('reverse'/'hilic') of the
    analysis, used to tell which raw files of a multi-analysis study belong to it."""
    ion = str(df_outer['ION_MODE'].iloc[0]).lower() if 'ION_MODE' in df_outer.columns else ''
    pos, neg = 'positive' in ion, 'negative' in ion
    chrom = str(df_outer['ChromatographyAndPhase'].iloc[0]).lower() if 'ChromatographyAndPhase' in df_outer.columns else ''
    return {'multi': _STUDY_N_MS_ANALYSES > 1,
            'polarity': 'positive' if pos and not neg else 'negative' if neg and not pos else None,
            'chrom': 'reverse' if 'reverse' in chrom else 'hilic' if 'hilic' in chrom else None}


def create_dataframe_outer_dict(MWB_mwTAB_dict, rest_response, raw_file_name_df=None, path_to_csvs = 'translation_sheets', **kwargs):
    entry = {
        'ANALYSIS_TYPE': _section(MWB_mwTAB_dict, 'ANALYSIS').get('ANALYSIS_TYPE', 'NA'),
        'PROJECT_TITLE': _section(MWB_mwTAB_dict, 'PROJECT').get('PROJECT_TITLE', 'NA'),
        'SUBJECT_TYPE': _section(MWB_mwTAB_dict, 'SUBJECT').get('SUBJECT_TYPE', 'NA'),
        'SPECIES_GROUP': _section(MWB_mwTAB_dict, 'SUBJECT').get('SPECIES_GROUP', 'NA'),
        'SUBJECT_GENDER': _section(MWB_mwTAB_dict, 'SUBJECT').get('GENDER', ''),
        'SUBJECT_AGE': _section(MWB_mwTAB_dict, 'SUBJECT').get('AGE_OR_AGE_RANGE', ''),
        'SUBJECT_SPECIES': _section(MWB_mwTAB_dict, 'SUBJECT').get('SUBJECT_SPECIES', 'NA'),
        'SAMPLE_TYPE': _section(MWB_mwTAB_dict, 'COLLECTION').get('SAMPLE_TYPE', 'NA'),
        'COLLECTION_LOCATION': _section(MWB_mwTAB_dict, 'COLLECTION').get('COLLECTION_LOCATION', 'NA'),
        'SampleCollectionMethod': _section(MWB_mwTAB_dict, 'COLLECTION').get('COLLECTION_METHOD', 'NA'),
        'COLLECTION_TIME': _section(MWB_mwTAB_dict, 'COLLECTION').get('COLLECTION_TIME', 'NA'),
        'MassSpectrometer': _section(MWB_mwTAB_dict, 'MS').get('INSTRUMENT_NAME', 'NA'),
        'MS_INSTRUMENT_TYPE': _section(MWB_mwTAB_dict, 'MS').get('INSTRUMENT_TYPE', ''),
        'MS_TYPE': _section(MWB_mwTAB_dict, 'MS').get('MS_TYPE', 'NA'),
        'ION_MODE': _section(MWB_mwTAB_dict, 'MS').get('ION_MODE', 'NA'),
        'InternalStandardsUsed': _section(MWB_mwTAB_dict, 'CHROMATOGRAPHY').get('INTERNAL_STANDARD', 'NA'),
        'COLUMN_NAME': _section(MWB_mwTAB_dict, 'CHROMATOGRAPHY').get('COLUMN_NAME', 'NA'),
        'CHROMATOGRAPHY_TYPE': _section(MWB_mwTAB_dict, 'CHROMATOGRAPHY').get('CHROMATOGRAPHY_TYPE', 'NA'),
        'SampleExtractionMethod': _section(MWB_mwTAB_dict, 'SAMPLEPREP').get('EXTRACTION_METHOD', 'NA'),
        'TAXONOMY_ID': _section(MWB_mwTAB_dict, 'SUBJECT').get('TAXONOMY_ID', 'NA'),
        'ACQUISITION_DATE': _section(MWB_mwTAB_dict, 'ANALYSIS').get('ACQUISITION_DATE', 'NA'),
        'CREATED_ON': _section(MWB_mwTAB_dict, 'METABOLOMICS WORKBENCH').get('CREATED_ON', 'NA'),
        'STUDY_ID': _section(MWB_mwTAB_dict, 'METABOLOMICS WORKBENCH').get('STUDY_ID', 'NA')

    }

    ontology_table = kwargs['ontology_table']

    if entry['ANALYSIS_TYPE'] != 'MS':
        raise ValueError("This is not an MS analysis!")

    df_outer = pd.DataFrame(entry, index=[0])
    df_outer['STUDY_uberon_exact'] = uberon_exact_name(df_outer['SAMPLE_TYPE'], ontology_table, allowedTerm_dict).values
    df_outer['STUDY_envo_exact'] = envo_material_exact_name(df_outer['SAMPLE_TYPE'], allowedTerm_dict).values

    #add NCBITaxonomy on dataset level
    #######
    unique_species = set(df_outer['SUBJECT_SPECIES'])
    unique_species_ncbiIDs = set(df_outer['TAXONOMY_ID'])
    processed_species = {species: get_taxonomy_id_from_name__allowedTerms(species, ncbi_id=unique_species_ncbiIDs, allowedTerm_dict=allowedTerm_dict) for species in unique_species}
    # names the lookup above misses for formatting reasons only: 'NCBI:txid51953',
    # 'Oat (Avena sativa L.)', 'Lactuca sativa L.'
    processed_species = {sp: (term if term is not None else study_species_variant_term(sp, df_outer['TAXONOMY_ID'].iloc[0], allowedTerm_dict))
                         for sp, term in processed_species.items()}
    df_outer['NCBITaxonomy'] = df_outer['SUBJECT_SPECIES'].map(processed_species)
    attempt_taxa_extraction_from_inner = False
    if all(value is None for value in processed_species.values()):
        attempt_taxa_extraction_from_inner = True


    #add LC method on dataset level
    #######
    if 'CHROMATOGRAPHY_TYPE' in df_outer.columns:
        df_outer['ChromatographyAndPhase'] = [
            chromatography_and_phase(t, c) for t, c in
            zip(df_outer['CHROMATOGRAPHY_TYPE'], df_outer['COLUMN_NAME'] if 'COLUMN_NAME' in df_outer.columns else [''] * len(df_outer))]
        df_outer['CHROMATOGRAPHY_TYPE'] = df_outer['CHROMATOGRAPHY_TYPE'].fillna('').astype(str).str.lower()


    #add YearOfAnalysis on dataset level
    df_outer['YearOfAnalysis'] = df_outer['ACQUISITION_DATE'].apply(lambda x: extract_years(x))
    df_outer['YearOfAnalysis'] = df_outer.apply(
        lambda x: extract_years(x['CREATED_ON']) if x['YearOfAnalysis'] is None else x['YearOfAnalysis'], axis=1)

    df_outer['IonizationSourceAndPolarity'] = df_outer['MS_TYPE'] + '|' + df_outer['ION_MODE']

    df_outer.loc[:, ['SampleType', 'SampleTypeSub1']] = 'missing value'

    processed_taxonomy = {taxonomy.split('|')[0]: get_taxonomy_info(taxonomy.split('|')[0])
                        for taxonomy in df_outer['NCBITaxonomy'].unique()
                        if taxonomy is not None and '|' in taxonomy and 'None' not in taxonomy}


    df_outer.loc[:, 'SampleType'] = df_outer['NCBITaxonomy'].map(
        lambda x: processed_taxonomy.get(x.split('|')[0], [pd.NA, pd.NA])[0]
        if x is not None and '|' in x and 'None' not in x else pd.NA
    )

    df_outer.loc[:, 'SampleTypeSub1'] = df_outer['NCBITaxonomy'].map(
        lambda x: processed_taxonomy.get(x.split('|')[0], [pd.NA, pd.NA])[1]
        if x is not None and '|' in x and 'None' not in x else pd.NA
    )


    df_outer[['SampleType', 'SampleTypeSub1']] = df_outer.apply(
        lambda x: get_enviromental_water(str(x['SAMPLE_TYPE'] or '') + str(x['SUBJECT_TYPE'] or '')) 
        if pd.isnull(x['SampleType']) and pd.isnull(x['SampleTypeSub1']) 
        else (x['SampleType'], x['SampleTypeSub1']), axis=1, result_type='expand'
    )

    # several species named ('Homo sapiens; Mus musculus', 'Quercus ilex, Quercus suber'):
    # no NCBITaxonomy, but the sample type they all share
    no_type = df_outer['SampleType'].isna() | df_outer['SampleType'].isin([None, 'missing value'])
    if no_type.any() and df_outer['NCBITaxonomy'].iloc[0] is None:
        shared = multi_species_sample_type(df_outer['SUBJECT_SPECIES'].iloc[0], df_outer['SUBJECT_TYPE'].iloc[0], allowedTerm_dict)
        if shared is not None:
            df_outer.loc[no_type, 'SampleType'] = shared[0]
            df_outer.loc[no_type, 'SampleTypeSub1'] = shared[1]

    # study-level SUBJECT_TYPE where the species gave no sample type ('Synthetic sample',
    # 'Soil sample', 'Plant', 'Bacteria', 'Human', ...)
    no_type = df_outer['SampleType'].isna() | df_outer['SampleType'].isin([None, 'missing value'])
    if no_type.any():
        from_subject = pd.Series([sample_type_from_subject_type(t, st) for t, st in zip(df_outer['SUBJECT_TYPE'], df_outer['SAMPLE_TYPE'])],
                                 index=df_outer.index, dtype=object)
        fill = no_type & from_subject.notna()
        if fill.any():
            df_outer.loc[fill, 'SampleType'] = from_subject[fill].str[0]
            df_outer.loc[fill, 'SampleTypeSub1'] = from_subject[fill].str[1]
    # human subjects without a species
    if df_outer['NCBITaxonomy'].iloc[0] is None and str(df_outer['SUBJECT_TYPE'].iloc[0]).strip().lower() in _HUMAN_SUBJECT_TYPES \
            and str(df_outer['SUBJECT_SPECIES'].iloc[0]).strip().lower() in ('', '-', 'na', 'n/a', 'none'):
        df_outer['NCBITaxonomy'] = '9606|Homo sapiens'

    if df_outer['NCBITaxonomy'].iloc[0] is None:
        add_to_cols = ['NCBITaxonomy']
    else:
        add_to_cols = []


    df_inner_SUBJECT_SAMPLE_FACTORS = create_dataframe_from_SUBJECT_SAMPLE_FACTORS(
        MWB_mwTAB_dict['SUBJECT_SAMPLE_FACTORS'],
        rest_response=rest_response,
        raw_file_name_df=raw_file_name_df,
        path_to_csvs=path_to_csvs,
        ontology_table=ontology_table,
        attempt_taxa_extraction_from_inner=attempt_taxa_extraction_from_inner,
        add_to_cols=add_to_cols,
        analysis_hints=_analysis_hints(df_outer)
    )

    if 'NCBITaxonomy' in df_inner_SUBJECT_SAMPLE_FACTORS.columns:
        df_outer.drop('NCBITaxonomy', axis=1, inplace=True)


    df_outer = pd.concat([df_outer] * len(df_inner_SUBJECT_SAMPLE_FACTORS), ignore_index=True)
    df_outer_inner = pd.concat([df_outer, df_inner_SUBJECT_SAMPLE_FACTORS], axis=1)

    df_outer_inner.loc[df_outer_inner['SampleType_inner'].notnull(), 'SampleType'] = df_outer_inner['SampleType_inner']
    df_outer_inner.loc[df_outer_inner['SampleTypeSub1_inner'].notnull(), 'SampleTypeSub1'] = df_outer_inner[
        'SampleTypeSub1_inner']

    df_outer_inner = df_outer_inner.drop(columns=['SampleType_inner', 'SampleTypeSub1_inner'])
    
    return df_outer_inner


# Per-sample Factor keys (besides "Sample Source") that commonly carry a body part /
# tissue / biofluid. Their values are matched against UBERONBodyPartName.csv to fill a
# body part when Sample Source didn't provide one. Environmental values (soil, water,
# ...) simply don't match the sheet and are left untouched.
UBERON_BODYPART_KEYS = {
    'sample.matrix', 'sample matrix', 'client matrix', 'matrix', 'material',
    'sample material', 'biological material', 'specimen_type', 'specimen type',
    'specimen', 'tissue', 'tissues', 'tissue type', 'tissue location', 'organ',
    'organism_part', 'organism part', 'body site', 'body_site', 'anatomical site',
}


def _fill_uberon_from_other_keys(MWB_table, df_translations):
    """Fill-only-if-missing body part from UBERON_BODYPART_KEYS. Never changes a value
    Sample Source already produced (so no loss); only adds where it was missing."""
    if 'Key' not in MWB_table.columns or 'Value' not in MWB_table.columns:
        return MWB_table
    key_l = MWB_table['Key'].astype(str).str.lower()
    cand = MWB_table.loc[key_l.isin(UBERON_BODYPART_KEYS), ['filename', 'Value']].copy()
    if cand.empty:
        return MWB_table
    cand['_mkey'] = cand['Value'].astype(str).str.lower()
    cand = cand.merge(df_translations[['MWB', 'REDU']], left_on='_mkey', right_on='MWB', how='left')
    fill = cand.dropna(subset=['REDU']).groupby('filename')['REDU'].first()
    if fill.empty:
        return MWB_table
    cur = MWB_table['UBERONBodyPartName']
    missing = cur.isna() | cur.astype(str).str.strip().str.lower().isin(['', 'nan', 'none', 'missing value'])
    filled = MWB_table['filename'].map(fill)
    MWB_table.loc[missing, 'UBERONBodyPartName'] = filled[missing].fillna(cur[missing])
    return MWB_table


# Sample Source values that are no body part, and labels too generic to be one
_UBERON_PLACEHOLDERS = {'', 'na', 'n/a', 'nan', 'none', 'other', 'others', 'unknown', 'not applicable', 'not available',
                        'missing', 'missing value', '-', '--', 'mixed', 'pool', 'pooled', 'qc', 'blank', 'sample', 'samples',
                        'control', 'standard', 'standards'}
_UBERON_TOO_GENERIC = {'cell', 'tissue', 'organ', 'anatomical entity', 'material anatomical entity', 'organism',
                       'multicellular organism', 'whole organism', 'plant', 'animal', 'body', 'organism substance',
                       'anatomical structure', 'portion of tissue', 'plant structure', 'plant anatomical entity'}
_UBERON_LOOKUP_CACHE = {}


def _uberon_exact_lookup(ontology_table, allowed):
    """lower-cased label / synonym (> 3 characters) -> label, for names that point to
    exactly one allowed, non-obsolete label. Labels win over synonyms."""
    key = (id(ontology_table), len(ontology_table))
    if key in _UBERON_LOOKUP_CACHE:
        return _UBERON_LOOKUP_CACHE[key]
    labels, synonyms = {}, {}
    for label, syn in ontology_table[['Label', 'Synonym']].itertuples(index=False):
        if not isinstance(label, str) or label not in allowed or label.lower().startswith('obsolete'):
            continue
        labels.setdefault(label.strip().lower(), set()).add(label)
        if isinstance(syn, str) and len(syn.strip()) > 3:
            synonyms.setdefault(syn.strip().lower(), set()).add(label)
    lookup = {n: next(iter(v)) for n, v in synonyms.items() if len(v) == 1}
    # a name that is any label's synonym but also an obsolete/disallowed label's name is
    # left to the label table below
    lookup.update({n: next(iter(v)) for n, v in labels.items() if len(v) == 1})
    for n, v in labels.items():
        if len(v) > 1:
            lookup.pop(n, None)
    lookup = {n: l for n, l in lookup.items()
              if n not in _UBERON_PLACEHOLDERS and l.lower() not in _UBERON_TOO_GENERIC}
    _UBERON_LOOKUP_CACHE.clear()
    _UBERON_LOOKUP_CACHE[key] = lookup
    return lookup


_ENVO_MATERIAL_TABLE = None
_ENV_SUBJECT_TYPES = {'water sample', 'soil sample', 'other abiotic sample', 'environmental sample', 'sediment sample'}
_ENVO_TOO_GENERIC = {'environmental material', 'anthropogenic environmental material', 'material', 'liquid water', 'water'}


def envo_material_exact_name(values, allowedTerm_dict):
    """Allowed ENVOEnvironmentMaterial label for each value that is exactly an ENVO
    material label or synonym (> 3 characters, one label only), else NaN. 'water' is
    left out: in lab studies it usually names a blank or solvent."""
    values = pd.Series(values)
    table = _ENVO_MATERIAL_TABLE
    if table is None or 'Synonym' not in table.columns:
        return pd.Series(np.nan, index=values.index, dtype=object)
    allowed = set((allowedTerm_dict.get('ENVOEnvironmentMaterial') or {}).get('allowed_values') or [])
    names = {}
    for label, syn in table[['Label', 'Synonym']].itertuples(index=False):
        if not isinstance(label, str) or label not in allowed:
            continue
        for n in {label, syn}:
            if isinstance(n, str) and len(n.strip()) > 3:
                names.setdefault(n.strip().lower(), set()).add(label)
    lookup = {n: next(iter(v)) for n, v in names.items()
              if len(v) == 1 and n not in _ENVO_TOO_GENERIC and next(iter(v)).lower() not in _ENVO_TOO_GENERIC}
    return values.astype(str).str.strip().str.lower().map(lookup)


def uberon_exact_name(values, ontology_table, allowedTerm_dict):
    """Allowed UBERON/CL/PO label for each value with an exact (case-insensitive)
    label or synonym match, else NaN."""
    values = pd.Series(values)
    if ontology_table is None or 'Synonym' not in ontology_table.columns:
        return pd.Series(np.nan, index=values.index, dtype=object)
    allowed = set((allowedTerm_dict.get('UBERONBodyPartName') or {}).get('allowed_values') or [])
    lookup = _uberon_exact_lookup(ontology_table, allowed)
    return values.astype(str).str.strip().str.lower().map(lookup)


def translate_MWB_to_REDU_from_csv(MWB_table,
                                   column_and_csv_names_outer=['MassSpectrometer',
                                                               'ChromatographyAndPhase',
                                                               'InternalStandardsUsed',
                                                               'SampleExtractionMethod',
                                                               'NCBITaxonomy',
                                                               'SampleCollectionMethod',
                                                               'IonizationSourceAndPolarity',
                                                               'Country'],
                                   column_and_csv_names_inner=['UBERONBodyPartName',
                                                               'DOIDCommonName',
                                                               'HumanPopulationDensity',
                                                               'HealthStatus',
                                                               'Country'],
                                   fill_col_from=[['UBERONBodyPartName', 'SAMPLE_TYPE']],
                                   case='outer',
                                   path_to_csvs='translation_sheets',
                                   add_to_cols=[],
                                   **kwargs):

    ontology_table = kwargs['ontology_table']
    allowedTerm_dict = kwargs['allowedTerm_dict']

    if case == 'outer':
        column_and_csv_names = column_and_csv_names_outer
        print(f"Translating MWB to REDU for outer columns: {column_and_csv_names}")
    elif case == 'inner':
        column_and_csv_names = column_and_csv_names_inner + add_to_cols
        print(f"Translating MWB to REDU for inner columns: {column_and_csv_names}")
        present_keys = MWB_table[['Key', 'Value']].copy()
        present_keys['Value'] = present_keys['Value'].apply(convert_to_numeric_or_original)
        present_keys.drop_duplicates(inplace=True)
    elif case == 'fill':
        column_and_csv_names = [x[0] for x in fill_col_from]
        origin_cols = [x[1] for x in fill_col_from]

    for index, col_csv_name in enumerate(column_and_csv_names):

        print(f"Translating column {col_csv_name}")

        df_translations = pd.read_csv(path_to_csvs + "/{}.csv".format(str(col_csv_name)), encoding="ISO-8859-1",
                                      dtype=str)
        df_translations['MWB'] = df_translations['MWB'].str.lower()
        df_translations = df_translations.drop_duplicates()

        if case == 'outer':

            if col_csv_name == 'NCBITaxonomy' and 'NCBITaxonomy' in MWB_table.columns:
                continue

            if col_csv_name == 'MassSpectrometer':
                MWB_table['_ms_raw'] = MWB_table[col_csv_name].astype(str)
            if col_csv_name == 'ChromatographyAndPhase':
                MWB_table['_chrom_derived'] = MWB_table[col_csv_name]
            MWB_table['_raw_value'] = MWB_table[col_csv_name].astype(str).str.strip().str.lower()
            MWB_table[col_csv_name] = MWB_table[col_csv_name].str.lower()
            MWB_table = pd.merge(MWB_table, df_translations, left_on=col_csv_name, right_on='MWB', how='left')
            MWB_table = MWB_table.drop(columns=['MWB', col_csv_name])
            MWB_table = MWB_table.rename(columns={'REDU': col_csv_name})

            # Workbench-only additions to a shared sheet (exact values), used only for
            # values the shared sheet does not translate
            _extra_sheet = os.path.join(path_to_csvs, 'MWB_{}.csv'.format(col_csv_name))
            if col_csv_name != 'MassSpectrometer' and os.path.exists(_extra_sheet):
                _extra = pd.read_csv(_extra_sheet, dtype=str, keep_default_na=False)
                _extra = dict(zip(_extra['MWB'].str.strip().str.lower(), _extra['REDU']))
                _allowed = set((allowedTerm_dict.get(col_csv_name) or {}).get('allowed_values') or [])
                _found = MWB_table['_raw_value'].map(_extra).where(lambda v: v.isin(_allowed))
                _fill = MWB_table[col_csv_name].isna() & _found.notna()
                MWB_table.loc[_fill, col_csv_name] = _found[_fill]
            MWB_table = MWB_table.drop(columns=['_raw_value'])

            if col_csv_name == 'ChromatographyAndPhase':
                # the value derived in create_dataframe_outer_dict already is a ReDU term;
                # keep it where the sheet has no row for it (e.g. 'gas chromatography (NOS)')
                allowed_chrom = set((allowedTerm_dict.get('ChromatographyAndPhase') or {}).get('allowed_values') or [])
                keep = MWB_table[col_csv_name].isna() & MWB_table['_chrom_derived'].isin(allowed_chrom)
                MWB_table.loc[keep, col_csv_name] = MWB_table.loc[keep, '_chrom_derived']
                MWB_table = MWB_table.drop(columns=['_chrom_derived'])

            if col_csv_name == 'MassSpectrometer':
                # Fallback for the ~58% of instrument strings the curated sheet misses
                # (Q Exactive HF-X, Orbitrap Exploris, Sciex QTrap, timsTOF, ...): map the
                # raw instrument name onto the allowed vocabulary by token-subset match.
                allowed_ms = (allowedTerm_dict.get('MassSpectrometer') or {}).get('allowed_values') or []
                unresolved = MWB_table[col_csv_name].isna()
                # explicit Workbench instrument-name aliases (exact names only)
                if allowed_ms and unresolved.any():
                    itype = MWB_table['MS_INSTRUMENT_TYPE'] if 'MS_INSTRUMENT_TYPE' in MWB_table.columns else pd.Series('', index=MWB_table.index)
                    aliased = [mwb_ms_alias(n, t, allowed_ms, path_to_csvs) for n, t in zip(MWB_table['_ms_raw'], itype)]
                    MWB_table.loc[unresolved, col_csv_name] = pd.Series(aliased, index=MWB_table.index)[unresolved]
                    unresolved = MWB_table[col_csv_name].isna()
                if allowed_ms and unresolved.any():
                    _cache = {}
                    def _resolve(raw):
                        # unambiguous vendor-normalised match first; the token-subset
                        # match only for names it cannot decide
                        if raw not in _cache:
                            _cache[raw] = match_mass_spectrometer(raw, allowed_ms) or map_instrument_to_allowed(raw, allowed_ms)
                        return _cache[raw]
                    MWB_table.loc[unresolved, col_csv_name] = MWB_table.loc[unresolved, '_ms_raw'].map(_resolve)
                MWB_table = MWB_table.drop(columns=['_ms_raw'])

            print(f"Unique values in {col_csv_name}: {MWB_table[col_csv_name].nunique()}")
               
        if case == 'inner':

            if col_csv_name != 'NCBITaxonomy':
                if col_csv_name == 'UBERONBodyPartName':
                    left = 'MWB_sampleSource'
                else:
                    left = 'Value'
                # df_translations['MWB'] is lowercased at load; match the factor value
                # case-insensitively so mixed-case values (e.g. "Crohn's disease",
                # "Asthma") also hit the sheet instead of being silently dropped.
                MWB_table['_mkey'] = MWB_table[left].astype(str).str.lower()
                MWB_table = MWB_table.merge(df_translations, left_on='_mkey', right_on='MWB', how='left')
                MWB_table = MWB_table.drop(columns=['_mkey'])
                MWB_table[col_csv_name] = MWB_table.groupby('filename')['REDU'].transform('first')
                if col_csv_name == 'UBERONBodyPartName':
                    # Sample Source is not the only key carrying a body part. Where it
                    # yielded nothing, fill from other per-sample body-part keys
                    # (matrix/material/tissue/organ/specimen/organism_part/body site).
                    MWB_table = _fill_uberon_from_other_keys(MWB_table, df_translations)
                    # Sample Source values the sheet does not list: exact UBERON/CL/PO label
                    # or synonym (used in translate_MWB_to_REDU_by_logic, after the
                    # study-level SAMPLE_TYPE fill, so it only fills what is left empty)
                    MWB_table['MWB_uberon_exact'] = uberon_exact_name(MWB_table['MWB_sampleSource'], ontology_table, allowedTerm_dict)
                    # ontology_table has one row per synonym; merge on unique labels so rows are not multiplied
                    MWB_table = pd.merge(MWB_table, ontology_table.drop_duplicates(subset=['Label']), left_on='UBERONBodyPartName', right_on='Label', how='left')
                    MWB_table = MWB_table.drop(columns=['REDU_UBERONOntologyIndex'])
                if col_csv_name == 'DOIDCommonName':
                    MWB_table['DOIDOntologyIndex'] = MWB_table.groupby('filename')['REDU_DOIDOntologyIndex'].transform('first')
                    MWB_table = MWB_table.drop(columns=['REDU_DOIDOntologyIndex'])

                MWB_table = MWB_table.drop(['MWB', 'REDU'], axis=1)

                print(f"Unique values in {col_csv_name}: {MWB_table[col_csv_name].nunique()}")

            elif col_csv_name == 'NCBITaxonomy':


                processed_species = {}
                none_count = 0
                # raw file names under variant keys and the Sample-ID placeholder rows are
                # no species names
                _not_species = (MWB_table['Key'] == 'raw_file_name') | (MWB_table['Key'] == 'Sample ID') \
                    | MWB_table['Key'].map(is_variant_raw_file_key)
                unique_species = set(MWB_table[~_not_species]['Value'])
                _names = _ncbi_name_lookup(allowedTerm_dict)

                for species in unique_species:
                    # same result as get_taxonomy_id_from_name__allowedTerms(species) (first
                    # allowed term with that name, case-insensitive), as a dictionary lookup
                    result = _names.get(str(species).lower())
                    if result is None:
                        none_count += 1
                        if none_count >= 100:
                            processed_species = {}
                            break
                    else:
                        none_count = 0  # Reset count if a non-None result is found
                        processed_species[species] = result

                if len(processed_species) > 0:

                    MWB_table['REDU'] = MWB_table['Value'].map(processed_species)

                    def get_longest_value(series):
                        lengths = series.str.len()  # Calculate length of each string
                        if lengths.isnull().all():  # Check if all values are None
                            return None
                        return series.loc[lengths.idxmax()]
                    
                    MWB_table['NCBITaxonomy'] = MWB_table.groupby('filename')['REDU'].transform(get_longest_value)
                    MWB_table = MWB_table.drop(['REDU'], axis=1)

        if case == 'fill':
            MWB_table[origin_cols[index]] = MWB_table[origin_cols[index]].str.lower()
            MWB_table = pd.merge(MWB_table, df_translations, left_on=origin_cols[index], right_on='MWB', how='left')
            MWB_table['REDU'] = MWB_table['REDU'].fillna('missing value')
            MWB_table[col_csv_name] = MWB_table[col_csv_name].fillna(MWB_table['REDU'])
            if col_csv_name == 'UBERONBodyPartName':
                MWB_table = MWB_table.drop(columns=['REDU_UBERONOntologyIndex'])
                MWB_table = MWB_table.drop(columns=['MWB', 'REDU', origin_cols[index]])


    return MWB_table


def translate_MWB_to_REDU_by_logic(MWB_table, path_to_csvs='translation_sheets'):
    if 'MWB_age' in MWB_table.columns:
        # MWB_age already holds AgeInYears (converted per-sample, unit-from-key) via
        # get_age_years_into_outer; 'NA' where no usable age. Route the 'NA' through
        # convert_to_years (-> None) but keep numeric years as-is.
        MWB_table['AgeInYears'] = MWB_table['MWB_age'].apply(
            lambda v: None if (isinstance(v, str) and not re.search(r'\d', v)) else v)
    if 'MWB_sex' in MWB_table.columns:
        MWB_table['BiologicalSex'] = MWB_table['MWB_sex'].map(convert_sex)
    if 'MWB_smoking' in MWB_table.columns:
        MWB_table['SmokingStatus'] = MWB_table['MWB_smoking'].map(convert_smoking)
    if 'MWB_bmi' in MWB_table.columns:
        MWB_table['BodyMassIndex'] = MWB_table['MWB_bmi'].apply(bmi_to_numeric)
    if 'MWB_diet' in MWB_table.columns:
        MWB_table['Diet'] = MWB_table['MWB_diet'].map(convert_diet)

    def _is_missing(col):
        if col not in MWB_table.columns:
            MWB_table[col] = np.nan
        return MWB_table[col].isna() | MWB_table[col].astype(str).str.strip().str.lower().isin(['', 'nan', 'none', 'na', 'missing value'])

    # values derived from extra factor keys only fill gaps; a problem here must never
    # cost the study its rows
    try:
        if 'MWB_healthy' in MWB_table.columns:
            fill = _is_missing('HealthStatus') & MWB_table['MWB_healthy'].notna()
            MWB_table.loc[fill, 'HealthStatus'] = MWB_table.loc[fill, 'MWB_healthy']
        if 'MWB_doid' in MWB_table.columns:
            fill = _is_missing('DOIDCommonName') & MWB_table['MWB_doid'].notna()
            if fill.any():
                doid = MWB_table.loc[fill, 'MWB_doid'].astype(str).str.split('|')
                MWB_table.loc[fill, 'DOIDCommonName'] = doid.str[0]
                MWB_table.loc[fill, 'DOIDOntologyIndex'] = doid.str[1]
        if 'MWB_depth' in MWB_table.columns:
            MWB_table['DepthorAltitudeMeters'] = MWB_table['MWB_depth']
        # study-level SUBJECT GENDER / AGE_OR_AGE_RANGE, only for studies without any
        # per-sample value (so a study-wide statement never mixes with per-sample ones)
        if 'SUBJECT_GENDER' in MWB_table.columns and _is_missing('BiologicalSex').all():
            study_sex = MWB_table['SUBJECT_GENDER'].map(study_level_sex)
            fill = _is_missing('BiologicalSex') & study_sex.notna()
            MWB_table.loc[fill, 'BiologicalSex'] = study_sex[fill]
        if 'SUBJECT_AGE' in MWB_table.columns and _is_missing('AgeInYears').all():
            study_age = MWB_table['SUBJECT_AGE'].map(study_level_age_years)
            fill = _is_missing('AgeInYears') & study_age.notna()
            MWB_table.loc[fill, 'AgeInYears'] = study_age[fill]
        if 'COLLECTION_LOCATION' in MWB_table.columns:
            allowed_countries = (allowedTerm_dict.get('Country') or {}).get('allowed_values') or []
            fill = _is_missing('Country')
            if fill.any() and allowed_countries:
                MWB_table.loc[fill, 'Country'] = MWB_table.loc[fill, 'COLLECTION_LOCATION'].map(
                    lambda v: country_from_location(v, allowed_countries))
        if 'MWB_collection_date' in MWB_table.columns:
            dates = MWB_table['MWB_collection_date'].where(MWB_table['MWB_collection_date'].astype(str) != 'NA')
            if 'COLLECTION_TIME' in MWB_table.columns:
                # study-level collection time only when it is a real calendar date
                study_date = MWB_table['COLLECTION_TIME'].map(lambda v: normalize_redu_datetime(v) if isinstance(v, str) and re.search(r'(19|20)\d\d', v) else None)
                dates = dates.fillna(study_date)
            MWB_table['SampleCollectionDateandTime'] = dates
    except Exception as e:
        print(f"Could not derive additional sample fields: {type(e).__name__}: {e}")

    # ENVO environmental material, exact label/synonym of the Sample Source (else the
    # study SAMPLE_TYPE), only in studies of environmental samples
    if 'SUBJECT_TYPE' in MWB_table.columns:
        env_study = (MWB_table['SUBJECT_TYPE'].astype(str).str.strip().str.lower().isin(_ENV_SUBJECT_TYPES)
                     | MWB_table['SampleType'].astype(str).eq('environmental'))
        envo = pd.Series(np.nan, index=MWB_table.index, dtype=object)
        for col in ['MWB_envo_exact', 'STUDY_envo_exact']:
            if col in MWB_table.columns:
                envo = envo.fillna(MWB_table[col])
        fill = env_study & envo.notna() & ~MWB_table['SampleType'].astype(str).str.contains('blank|QC', case=False)
        if fill.any():
            MWB_table.loc[fill, 'ENVOEnvironmentMaterial'] = envo[fill]

    # body part from an exact UBERON/CL/PO name: the sample's Sample Source, else the
    # study-level SAMPLE_TYPE (only where the translation sheets gave none)
    # (not where the sample type already has a finer SampleTypeSub1, e.g. 'insect', which
    # a body part would turn into 'tissue')
    _sub1_open = MWB_table['SampleTypeSub1'].isna() | MWB_table['SampleTypeSub1'].astype(str).isin(['', 'nan', 'None', 'missing value', 'tissue', 'biofluid']) \
        if 'SampleTypeSub1' in MWB_table.columns else pd.Series(True, index=MWB_table.index)
    for col in ['MWB_uberon_exact', 'STUDY_uberon_exact']:
        if col in MWB_table.columns:
            fill = _is_missing('UBERONBodyPartName') & MWB_table[col].notna() & _sub1_open
            MWB_table.loc[fill, 'UBERONBodyPartName'] = MWB_table.loc[fill, col]

    df_translations = pd.read_csv(path_to_csvs + "/biofluid_tissue_distinction.csv", encoding="ISO-8859-1", dtype=str)

    MWB_table = MWB_table.merge(df_translations, left_on='UBERONBodyPartName', right_on='sampletype', how='left')
    
    # Fill NA values before any string conversion (astype(str) turns NaN into 'nan', which fillna cannot see)
    MWB_table['SampleTypeSub1'] = MWB_table['SampleTypeSub1'].fillna(MWB_table['tissue_vs_biofluid'])

    MWB_table = MWB_table.drop(columns=['sampletype', 'tissue_vs_biofluid'])


    # List of columns to check and add if not present
    columns_to_check = ['NCBITaxonomy', 'SampleCollectionMethod', 'SampleExtractionMethod',
                        'InternalStandardsUsed', 'SubjectIdentifierAsRecorded', 'AgeInYears',
                        'BiologicalSex', 'UBERONBodyPartName', 'HealthStatus', 'DOIDCommonName',
                        'ComorbidityListDOIDIndex', 'Country', 'HumanPopulationDensity',
                        'LatitudeandLongitude', 'SmokingStatus', 'BodyMassIndex', 'Diet']

    # Iterate through each column in the list
    for column in columns_to_check:
        # Check if the column exists in the DataFrame
        if column not in MWB_table.columns:
            # If the column does not exist, create it with all values set to 'missing value'
            MWB_table[column] = 'missing value'
        else:
            # If the column exists, replace NaNs with 'missing value'
            MWB_table[column] = MWB_table[column].replace(np.nan, 'missing value')

    return MWB_table



# Tolerance for the sex-distribution guard: max |frac_female_output - frac_female_source|
# before the per-file sex assignment is considered untrustworthy for a study.
SEX_GUARD_TOLERANCE = 0.12


def _study_source_sex_fracF(rest_response):
    """Authoritative per-sample sex distribution from the study /factors response,
    using the pipeline's own convert_sex (so 'F'/'M', German 'w'/'m' via the same
    mapping the pipeline uses). Dedupe by sample id. Returns (n_with_sex, frac_female)
    or (0, None) when the study carries no recognizable sex."""
    recs = list(rest_response.values()) if isinstance(rest_response, dict) else rest_response
    if not isinstance(recs, list):
        return (0, None)
    seen = {}
    for rec in recs:
        if not isinstance(rec, dict):
            continue
        sid = rec.get('local_sample_id') or rec.get('mb_sample_id')
        fac = rec.get('factors', '')
        if sid is None or not isinstance(fac, str):
            continue
        for part in fac.split('|'):
            if ':' not in part:
                continue
            k, v = part.split(':', 1)
            if k.strip().lower() in ('sex', 'gender', 'biological_sex', 'biological sex'):
                s = convert_sex(v.strip())
                if s in ('male', 'female'):
                    seen[sid] = s
                break
    if not seen:
        return (0, None)
    n = len(seen)
    return (n, sum(1 for s in seen.values() if s == 'female') / n)


def apply_sex_distribution_guard(redu_df_final, rest_response, study_id=''):
    """Safety net against raw-file mis-attribution. The per-sample sex is attached to
    raw files via (sometimes fuzzy) filename matching; when that matching collapses
    many files onto few samples the per-file BiologicalSex comes out wrong (e.g. an
    all-female study where the source is 50/50). If the output BiologicalSex
    distribution disagrees with the authoritative per-sample distribution from the
    study /factors response beyond SEX_GUARD_TOLERANCE, the file<->sample mapping can't
    be trusted for sex, so BiologicalSex is suppressed for the whole study rather than
    emit wrong values. Never removes rows and never touches any other field."""
    if 'BiologicalSex' not in redu_df_final.columns or len(redu_df_final) == 0:
        return redu_df_final
    out = redu_df_final['BiologicalSex'].astype(str).str.lower()
    out = out[out.isin(['male', 'female'])]
    if len(out) == 0:
        return redu_df_final
    out_fracF = (out == 'female').mean()
    src_n, src_fracF = _study_source_sex_fracF(rest_response)
    if src_n == 0 or src_fracF is None:
        return redu_df_final
    if abs(src_fracF - out_fracF) > SEX_GUARD_TOLERANCE:
        print(f"[sex-guard] {study_id}: output BiologicalSex (fracF={out_fracF:.2f}, "
              f"n={len(out)}) disagrees with source (fracF={src_fracF:.2f}, n={src_n}); "
              f"suppressing BiologicalSex for this study to avoid mis-attributed values.")
        redu_df_final = redu_df_final.copy()
        redu_df_final['BiologicalSex'] = 'missing value'
    return redu_df_final


def MWB_to_REDU_study_wrapper(study_id, path_to_csvs='translation_sheets',
                              duplicate_raw_file_handling='remove_duplicates', export_to_tsv = False, **kwargs):


    allowedTerm_dict = kwargs['allowedTerm_dict']
    ontology_table = kwargs['ontology_table']
    ENVOEnvironmentBiomeIndex_table = kwargs['ENVOEnvironmentBiomeIndex_table']
    ENVOEnvironmentMaterialIndex_table = kwargs['ENVOEnvironmentMaterialIndex_table']
    global _ENVO_MATERIAL_TABLE
    _ENVO_MATERIAL_TABLE = ENVOEnvironmentMaterialIndex_table
    NCBIRankDivision_table = kwargs['NCBIRankDivision_table']

    raw_file_name_df = _get_metabolomicsworkbench_files(study_id)
    if len(raw_file_name_df) == 0:
        print('No raw data detected. Ignoring {}'.format(str(study_id)))
        return None

    raw_file_name_df['filename_base'] = raw_file_name_df['filename'].apply(lambda x: os.path.basename(x))

    stdy_info_req = _mwb_get(
        'https://www.metabolomicsworkbench.org/rest/study/study_id/{}/analysis'.format(str(study_id)))
    stdy_info = stdy_info_req.json()

    if not isinstance(next(iter(stdy_info.values())), dict):
        stdy_info = {'1': stdy_info}

    response = _mwb_get(f'https://www.metabolomicsworkbench.org/rest/study/study_id/{str(study_id)}/factors')
    rest_response = response.json()
    # a study with a single sample gets the record itself, not {'1': record}
    if isinstance(rest_response, dict) and 'local_sample_id' in rest_response:
        rest_response = {'1': rest_response}
    elif not isinstance(rest_response, dict):
        rest_response = {}

    global _STUDY_N_MS_ANALYSES
    _STUDY_N_MS_ANALYSES = sum(1 for a in stdy_info.values() if isinstance(a, dict) and a.get('analysis_type') == 'MS')

    redu_dfs = []
    for analysis_id, analysis_details in stdy_info.items():

        try:
            if analysis_details['analysis_type'] == 'MS':
                print('Initiating REDU-table generation for ' + str(analysis_details["analysis_id"]))
            else:
                print('This is not an MS-analysis. Ignoring ' + str(analysis_details["analysis_id"]))
                continue
        
        except Exception as e:
            traceback_info = traceback.format_exc()
            print(f"No analysis_type-key in {study_id}: {e}\nTraceback:\n{traceback_info}")
            continue

        
        redu_df = MWB_to_REDU_wrapper(MWB_analysis_ID=analysis_details["analysis_id"],
                                      rest_response=rest_response,
                                      raw_file_name_df=raw_file_name_df[['filename', 'filename_base', 'USI']],
                                      path_to_csvs=path_to_csvs,
                                      Massive_ID=study_id + '|' + str(analysis_details["analysis_id"]),
                                      allowedTerm_dict=allowedTerm_dict,
                                      ontology_table=ontology_table,
                                      ENVOEnvironmentBiomeIndex_table=ENVOEnvironmentBiomeIndex_table,
                                      ENVOEnvironmentMaterialIndex_table=ENVOEnvironmentMaterialIndex_table)
        
        if isinstance(redu_df, pd.DataFrame):
            redu_dfs.append(redu_df)

    
    if len(redu_dfs) > 0:
        redu_df_final = pd.concat(redu_dfs)
        
        
        if duplicate_raw_file_handling == 'remove_duplicates':
            duplicates = redu_df_final.filename.value_counts() > 1
            redu_df_final = redu_df_final[~redu_df_final.filename.isin(duplicates.index[duplicates])]

        if duplicate_raw_file_handling == 'keep_pols_dupl':
            mask = redu_df_final.duplicated(subset='filename', keep=False)

            remove_mask = mask.copy()
            for i, row in redu_df_final[mask].iterrows():
                if ((redu_df_final['filename'] == row['filename']) & mask).sum() == 2:

                    negative_mask = (redu_df_final['filename'] == row['filename']) & (
                        redu_df_final['IonizationSourceAndPolarity'].str.contains('negative'))
                    positive_mask = (redu_df_final['filename'] == row['filename']) & (
                        redu_df_final['IonizationSourceAndPolarity'].str.contains('positive'))

                    if negative_mask.sum() == 1 and positive_mask.sum() == 1:
                        remove_mask &= ~((redu_df_final['filename'] == row['filename']) & mask)

            # Remove all rows with entries occurring more than once in col1 (but not if one has 'negative' and the other has 'positive' in col2)
            redu_df_final = redu_df_final[~remove_mask]
        
    else:
        return None


    redu_df_final = merge_repeated_fileobservations_across_mwatb(redu_df_final, polarity_table=polarity_table)
    try:
        redu_df_final = propagate_to_same_name_files(redu_df_final, raw_file_name_df)
    except Exception as e:
        print(f"Could not copy metadata to same-name files: {type(e).__name__}: {e}")

    # QC pools and reference materials are not samples of the study's subjects: like
    # blanks (cleared in complete_and_fill_REDU_table) they get no organism, body part
    # or subject fields
    _not_subject = (redu_df_final['SampleType'].isin(['pool_QC', 'reference material'])
                    | redu_df_final['SampleType'].astype(str).str.contains('blank', case=False))
    if _not_subject.any():
        _subject_cols = ['NCBITaxonomy', 'UBERONBodyPartName', 'UBERONOntologyIndex', 'AgeInYears', 'UniqueSubjectID',
                         'SubjectIdentifierAsRecorded', 'DOIDCommonName', 'DOIDOntologyIndex', 'BiologicalSex',
                         'HealthStatus', 'LatitudeandLongitude', 'SmokingStatus', 'BodyMassIndex', 'Diet', 'LifeStage']
        redu_df_final.loc[_not_subject, [c for c in _subject_cols if c in redu_df_final.columns]] = 'missing value'

    ontology_table = ontology_table.drop_duplicates(subset=['Label'])
    # DOID table (one row per label) so complete_and_fill_REDU_table can generate DOIDOntologyIndex
    doid_table = kwargs.get('DOIDOntologyIndex_table')
    if doid_table is None:
        doid_table = DOID_TABLE
    doid_kwargs = {}
    if doid_table is not None and len(doid_table):
        doid_kwargs['DOIDOntologyIndex_table'] = doid_table.drop_duplicates(subset=['Label'])[['Label', 'DOIDOntologyIndex']].copy()
    redu_df_final = complete_and_fill_REDU_table(redu_df_final, allowedTerm_dict, UBERONOntologyIndex_table=ontology_table, **doid_kwargs,
                                                 ENVOEnvironmentBiomeIndex_table=ENVOEnvironmentBiomeIndex_table.drop_duplicates(subset=["Label"]),
                                                 ENVOEnvironmentMaterialIndex_table=ENVOEnvironmentMaterialIndex_table.drop_duplicates(subset=["Label"]),
                                                 NCBIRankDivision_table=NCBIRankDivision_table,
                                                 add_usi = True, other_allowed_file_extensions = ['.raw', '.cdf', '.wiff', '.d'])

    redu_df_final = apply_sex_distribution_guard(redu_df_final, rest_response, study_id=str(study_id))

    if export_to_tsv == True:
        redu_df_final.to_csv('{}_REDU_from_MWB.tsv'.format(study_id), sep='\t', index=False, header=True)
        return None

    return redu_df_final


def MWB_to_REDU_wrapper(mwTab_json=None, rest_response=None, MWB_analysis_ID=None, raw_file_name_df=None, Massive_ID='',
                        path_to_csvs='translation_sheets', **kwargs):
    if mwTab_json is None and MWB_analysis_ID is None:
        raise SystemExit("mwTab_json or MWB_analysis_ID has to be provided!")

    if mwTab_json is not None and MWB_analysis_ID is not None:
        raise SystemExit("Only mwTab_json OR MWB_analysis_ID should be provided!")

    if MWB_analysis_ID is not None:
        try:
            mwb_anlysis_req = _mwb_get(
                "https://www.metabolomicsworkbench.org/rest/study/analysis_id/{}/mwtab".format(str(MWB_analysis_ID)))
        except Exception as e:
            # one unavailable mwTab must not cost the study its other analyses
            print(f"Could not get the mwTab of {MWB_analysis_ID}: {type(e).__name__}: {e}")
            return None

        try:
            mwTab_json = json.loads(mwb_anlysis_req.text, object_pairs_hook=handle_duplicates)

        except json.JSONDecodeError:
            print("Did not receive valid mwTab json for {}!".format(str(MWB_analysis_ID)))
            return None
        
    allowedTerm_dict = kwargs['allowedTerm_dict']
    ontology_table = kwargs['ontology_table']

    try:

        #create dataframe from mwTab file only considering per study variables
        df_outer_dict = create_dataframe_outer_dict(mwTab_json, rest_response, raw_file_name_df=raw_file_name_df, path_to_csvs=path_to_csvs, allowedTerm_dict=allowedTerm_dict, ontology_table=ontology_table)

        #translate terms from per-study-variables to ReDU ontology 
        df_outer_dict_REDUfied = translate_MWB_to_REDU_from_csv(df_outer_dict,case='outer',path_to_csvs=path_to_csvs,ontology_table=ontology_table,allowedTerm_dict=allowedTerm_dict)

        #
        df_outer_dict_REDUfied_filled = translate_MWB_to_REDU_from_csv(df_outer_dict_REDUfied,case='fill',path_to_csvs=path_to_csvs,ontology_table=ontology_table,allowedTerm_dict=allowedTerm_dict)
        
        #utilize substring matches and other logic to get to harmonized ontology
        complete_df = translate_MWB_to_REDU_by_logic(df_outer_dict_REDUfied_filled,path_to_csvs=path_to_csvs)


    except Exception as e:
        traceback_info = traceback.format_exc()
        print(f"An error occurred in MWB_to_REDU_wrapper: {e}\nTraceback:\n{traceback_info}")

        return None
    
    

    complete_df['filename'] = complete_df['filename_raw_path']
    complete_df = complete_df[complete_df['filename'] != ''].dropna(subset=['filename'])


    complete_df['MassiveID'] = Massive_ID
    complete_df['MassiveID'] = complete_df['MassiveID'].apply(lambda x: x.split('|')[0])

    complete_df['UniqueSubjectID'] = complete_df.apply(
        lambda x: str(x['MassiveID']) + '_' + str(x['SubjectIdentifierAsRecorded']) 
        if x['SubjectIdentifierAsRecorded'] != '' else None,
        axis=1)
    
    

    complete_df['LifeStage'] = complete_df['AgeInYears'].apply(lambda x: age_category(x))
    if 'MWB_lifestage' in complete_df.columns:
        # LifeStage is defined for humans only
        is_human = complete_df['NCBITaxonomy'].astype(str).str.startswith('9606|')
        fill = (complete_df['LifeStage'].isin(['', None]) | complete_df['LifeStage'].isna()) & complete_df['MWB_lifestage'].notna() & is_human
        complete_df.loc[fill, 'LifeStage'] = complete_df.loc[fill, 'MWB_lifestage']

    for col in ['SampleCollectionDateandTime', 'DepthorAltitudeMeters', 'ENVOEnvironmentMaterial']:
        if col not in complete_df.columns:
            complete_df[col] = 'missing value'
    complete_df['qiita_sample_name'] = 'missing value'

    missing_not_imported = ["SampleType",
                            "SampleTypeSub1",
                            "NCBITaxonomy",
                            "YearOfAnalysis",
                            "SampleCollectionMethod",
                            "SampleExtractionMethod",
                            "InternalStandardsUsed",
                            "MassSpectrometer",
                            "IonizationSourceAndPolarity",
                            "ChromatographyAndPhase",
                            "LifeStage",
                            "Country",
                            "UBERONOntologyIndex",
                            "DOIDOntologyIndex",
                            "SubjectIdentifierAsRecorded",
                            "BiologicalSex",
                            "UBERONBodyPartName",
                            "HealthStatus",
                            "DOIDCommonName",
                            "ComorbidityListDOIDIndex",
                            "SampleCollectionDateandTime",
                            "HumanPopulationDensity",
                            "LatitudeandLongitude",
                            "DepthorAltitudeMeters",
                            "qiita_sample_name",
                            "SmokingStatus",
                            "BodyMassIndex",
                            "Diet"
                            ]
    
    complete_df['ENVOEnvironmentMaterial'] = complete_df['ENVOEnvironmentMaterial'].replace({np.nan: 'missing value', None: 'missing value'})
    complete_df[missing_not_imported] = complete_df[missing_not_imported].replace(
        {np.nan: 'missing value', None: 'missing value'})

    missing_not_collected = ["AgeInYears"]
    complete_df[missing_not_collected] = complete_df[missing_not_collected].replace(
        {np.nan: 'not collected', None: 'not collected'})

    complete_df[['TermsofPosition',
                 'qiita_sample_name']] = 'missing value'

    REDU_df = complete_df[["MassiveID",
                           "filename",
                           "SampleType",
                           "SampleTypeSub1",
                           "NCBITaxonomy",
                           "YearOfAnalysis",
                           "SampleCollectionMethod",
                           "SampleExtractionMethod",
                           "InternalStandardsUsed",
                           "MassSpectrometer",
                           "IonizationSourceAndPolarity",
                           "ChromatographyAndPhase",
                           "SubjectIdentifierAsRecorded",
                           "AgeInYears",
                           "BiologicalSex",
                           "UBERONBodyPartName",
                           "TermsofPosition",
                           "HealthStatus",
                           "DOIDCommonName",
                           "ComorbidityListDOIDIndex",
                           "SampleCollectionDateandTime",
                           "Country",
                           "HumanPopulationDensity",
                           "LatitudeandLongitude",
                           "DepthorAltitudeMeters",
                           "qiita_sample_name",
                           "SmokingStatus",
                           "BodyMassIndex",
                           "Diet",
                           "UniqueSubjectID",
                           "LifeStage",
                           "UBERONOntologyIndex",
                           "DOIDOntologyIndex",
                           "ENVOEnvironmentMaterial",
                           "USI"
                           ]]
    
    return REDU_df


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Give an MWB study ID and get a REDU table tsv.')
    parser.add_argument("--study_id", "-mwb_id", type=str, help='An MWB study ID such as "ST002050". If "ALL" all study IDs are requested.', required=True)
    parser.add_argument("--path_to_csvs", "-csvs", type=str, help="Path to the translation csvs holding translations from MWB to REDU vocabulary (optional)", default="translation_sheets")
    parser.add_argument("--path_to_allowed_term_json", "-json", type=str, help="Allowed vocabulary json (optional)", default="allowed_terms")
    parser.add_argument("--path_ncbi_rank_division", type=str, help="Path to the path_ncbi_rank_division")
    parser.add_argument("--path_to_uberon_cl_po_csv", type=str, help="Path to the prepared uberon_cl_po ontology csv")
    parser.add_argument("--path_to_envo_biome_csv", type=str, help="Path to the prepared uberon_cl_po ontology csv")
    parser.add_argument("--path_to_envo_material_csv", type=str, help="Path to the prepared uberon_cl_po ontology csv")
    parser.add_argument("--duplicate_raw_file_handling", "-duplStrat", type=str, help="What should be done with duplicate filenames across studies? Can be 'keep_pols_dupl' to keep cases where files can be distinguished by their polarity or 'remove_duplicates' to only keep cases where files can be assigned unambiguously (i.e. cases with only one analysis per study_id)(optional)", default='remove_duplicates')
    parser.add_argument("--path_to_polarity_info", type=str, help="Path to the polarity file.", default='none')
    parser.add_argument("--path_to_doid_csv", type=str, help="Prepared DOID ontology csv (optional); enables exact disease-name matching", default='none')

    print('Starting MWB2REDU script,..')

    args = parser.parse_args()

    study_id = args.study_id
    path_to_csvs = args.path_to_csvs
    duplicate_raw_file_handling = args.duplicate_raw_file_handling
    path_to_allowed_term_json = args.path_to_allowed_term_json

    # Read allowed terms json
    with open(path_to_allowed_term_json, 'r') as file:
        allowedTerm_dict = json.load(file)

    # Read ontology
    ontology_table = pd.read_csv(args.path_to_uberon_cl_po_csv)
    ENVOEnvironmentBiomeIndex_table = pd.read_csv(args.path_to_envo_biome_csv)
    ENVOEnvironmentMaterialIndex_table = pd.read_csv(args.path_to_envo_material_csv)

    # Read polarity file
    polarity_table = pd.read_csv(args.path_to_polarity_info)

    if args.path_to_doid_csv != 'none':
        DOID_TABLE = pd.read_csv(args.path_to_doid_csv, dtype=str)
        DOID_LOOKUP = build_doid_lookup(DOID_TABLE)
   
    # Read NCBI rank and division csv
    NCBIRankDivision_table = pd.read_csv(args.path_ncbi_rank_division, index_col = False)
    NCBIRankDivision_table = NCBIRankDivision_table.drop_duplicates(subset=['TaxonID'])

    # result
    if study_id == "ALL":
        # Getting all files
        url = "https://www.metabolomicsworkbench.org/rest/study/study_id/ST/available"
        studies_dict = _mwb_get(url).json()

        study_list = []
        for key in studies_dict.keys():
            study_dict = studies_dict[key]
            study_list.append(study_dict['study_id'])

        study_list = list(set(study_list))

        all_results_list = []
        failed_studies = []
        for study_id in tqdm(study_list):
            print("Processing ", study_id)

            try:
                result = MWB_to_REDU_study_wrapper(study_id=study_id,
                                          path_to_csvs=path_to_csvs,
                                          duplicate_raw_file_handling=duplicate_raw_file_handling,
                                          export_to_tsv=False,
                                          allowedTerm_dict=allowedTerm_dict,
                                          ontology_table=ontology_table,
                                          ENVOEnvironmentBiomeIndex_table=ENVOEnvironmentBiomeIndex_table,
                                          ENVOEnvironmentMaterialIndex_table=ENVOEnvironmentMaterialIndex_table,
                                          polarity_table=polarity_table,
                                          NCBIRankDivision_table=NCBIRankDivision_table)
                if result is not None and len(result) > 0:
                    print('Extracted information for {} samples.'.format(len(result)))
                    all_results_list.append(result)
            except KeyboardInterrupt:
                raise
            except Exception as e:
                print(f"Failed to process {study_id}: {type(e).__name__}: {e}")
                failed_studies.append(study_id)
        print(f"{len(failed_studies)} of {len(study_list)} studies failed to process")
        merged_df = pd.concat(all_results_list, ignore_index=True)
        merged_df.to_csv('REDU_from_MWB_all.tsv', sep='\t', index=False, header=True)

    else:
        MWB_to_REDU_study_wrapper(study_id=study_id,
                                  path_to_csvs=path_to_csvs,
                                  duplicate_raw_file_handling=duplicate_raw_file_handling,
                                  export_to_tsv=True,
                                  allowedTerm_dict=allowedTerm_dict,
                                  ontology_table=ontology_table,
                                  ENVOEnvironmentBiomeIndex_table=ENVOEnvironmentBiomeIndex_table,
                                  ENVOEnvironmentMaterialIndex_table=ENVOEnvironmentMaterialIndex_table,
                                  polarity_table=polarity_table,
                                  NCBIRankDivision_table=NCBIRankDivision_table)

    print("Output files written to working directory")
