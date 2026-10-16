import os
import sys
import traceback
import pandas as pd
import argparse
import collections
import requests
from time import sleep
import glob
from io import StringIO
from subprocess import PIPE, run
import json
import re
from pathlib import Path

from merge_metadata import FINGERPRINT_COLUMNS, fingerprint_column, file_identity, identical_files, without_method_fixed

ccms_peak_link = "https://datasetcache.gnps2.org/datasette/datasette/database/uniquemri.csv?_sort=usi&dataset__exact=" # MSV000081468&filepath__endswith=%25.mz%25ML&_size=max"
# MetadataFileDataset: accession of the dataset the metadata file was found in (MassIVE)
# or named after (GitHub); merge_metadata.py prefers a dataset's own file on conflicts
gnps_column_names_added = ['USI', 'MetadataFileDataset']
# metadata files whose rows (for other datasets, from older versions of a dataset's sheet; see
# FILL_ONLY_MARK in gnps_downloader.py) may only fill fields of rows from the dataset's own
# metadata: a row of such a file is kept only if a row from its dataset's own metadata names the
# same file (merge_metadata.py then fills only the fields those rows leave missing)
FILL_ONLY_MARK = '_fill_only_gnps_metadata.tsv'


# a fill-only row whose value in one of these columns differs from the value of the own rows
# of the same file describes another sample (e.g. a copied sheet saying Mus musculus for
# MSV000085862's Nepenthes files) and fills nothing
FILL_ONLY_IDENTITY_COLUMNS = ['SampleType', 'SampleTypeSub1', 'NCBITaxonomy', 'UBERONBodyPartName']


def _present(value):
    return not pd.isna(value) and str(value).strip().lower() not in ('', 'nan', 'missing value')


def keep_fill_only_rows_of_described_files(frames):
    """frames: [(is_fill_only, matched rows)] -> list of tables; a fill-only row is kept only if
    a row of its dataset's own metadata (MetadataFileDataset == dataset) has the same USI and
    it does not contradict the own rows in FILL_ONLY_IDENTITY_COLUMNS."""
    own_values = collections.defaultdict(lambda: collections.defaultdict(set))
    for fill_only, df in frames:
        if not fill_only:
            own = df[df['MetadataFileDataset'] == df['ATTRIBUTE_DatasetAccession']]
            for row in own.reindex(columns=['USI'] + FILL_ONLY_IDENTITY_COLUMNS).itertuples(index=False):
                values = own_values[row[0]]
                for column, value in zip(FILL_ONLY_IDENTITY_COLUMNS, row[1:]):
                    if _present(value):
                        values[column].add(str(value).strip())
                    else:
                        values[column]
    kept = []
    for fill_only, df in frames:
        if fill_only:
            keep = []
            for row in df.reindex(columns=['USI', 'MetadataFileDataset', 'ATTRIBUTE_DatasetAccession'] + FILL_ONLY_IDENTITY_COLUMNS).itertuples(index=False):
                if row[0] not in own_values or row[1] == row[2]:
                    MATCH_STATS['rows_fill_only_dropped_no_own_row'] += 1
                    keep.append(False)
                    continue
                values = own_values[row[0]]
                conflict = any(_present(v) and values[c] and str(v).strip() not in values[c]
                               for c, v in zip(FILL_ONLY_IDENTITY_COLUMNS, row[3:]))
                MATCH_STATS['rows_fill_only_dropped_contradicting_own_rows' if conflict else 'rows_fill_only_kept_to_fill_own_rows'] += 1
                keep.append(not conflict)
            df = df[keep]
        kept.append(df)
    return kept


def _make_usi_from_filename(filename, dataset_id):
    # Replace initial "f." with "mzspec:"
    if filename.startswith("f.MSV"):
        usi = "mzspec:" + filename[2:]
        # Replace the first "/" with ":"
        first_slash_index = usi.find('/')
        if first_slash_index != -1:
            usi = usi[:first_slash_index] + ':' + usi[first_slash_index+1:]
    elif filename.startswith("f."):
        usi = "mzspec:" + dataset_id + ':' + filename[2:]
    else:
        print('Do filenames not start with f. anymore?')
        usi = "mzspec:" + dataset_id + ':' + filename

    return usi

def longest_common_dir_suffix_score(meta_dirs, cand_dirs):
    score = 0
    for md, cd in zip(reversed(meta_dirs), reversed(cand_dirs)):
        if md == cd:
            score += 1
        else:
            break
    return score


# dataset accessions of the dataset cache: MSV/RMSV (MassIVE), MTBLS, ST (MASST rows)
_ACCESSION = re.compile(r'^[A-Z]+\d+$')
_MZ_EXTENSIONS = ('.mzml', '.mzxml')
PREFERRED_DIRECTORIES = ['ccms_peak', 'peak', 'raw']

# matching rules, in the order they are tried for a metadata row; a later rule is only
# used when the earlier ones find no file of that name at all
#   exact        same file name (or the x.mzXML -> x.mzML converted copy, as before)
#   stem_ext     same name, the other of .mzML/.mzXML (x.mzML listed, only x.mzXML uploaded)
#   case         same name ignoring upper/lower case (also with the other of .mzML/.mzXML)
MATCH_STATS = collections.Counter()
AMBIGUOUS = []
# one entry per metadata row (written with --match_report)
MATCH_REPORT = []


# The dataset cache's whole file table (https://datasetcache.gnps2.org/dataset/uniquemri, a
# ~1 GB export the cache regenerates about every 3 h), downloaded once per run by the
# workflow (--files_table): dataset -> its mz(X)ML file rows sorted by usi, i.e. what the
# per-dataset query below returns. One query per dataset took 20-90 min per run and a query
# that timed out dropped all rows of its dataset (MSV000087612, 192 rows, on 2026-10-14).
# Set by load_files_table().
FILES_TABLE = None
# datasets whose file list could not be fetched this run (not: fetched and empty)
FETCH_FAILED = set()
# the condition of the per-dataset query, filepath__endswith=%.mz%ML, i.e. SQLite
# "filepath LIKE '%.mz%ML'" (LIKE ignores upper/lower case of ASCII letters)
_LIKE_MZ_ML = r'\.[mM][zZ].*[mM][lL]\Z'


def load_files_table(path):
    """Reads the uniquemri export once and indexes its mz(X)ML rows by dataset (with the size and
    the spectra fingerprint of each file, merge_metadata.FINGERPRINT_COLUMNS, where the export has
    them: they tell identical copies apart, see _break_tie)."""
    global FILES_TABLE
    wanted = {'usi', 'filepath', 'dataset', 'size', 'classification'} | set(FINGERPRINT_COLUMNS)
    parts = []
    for chunk in pd.read_csv(path, usecols=lambda c: c in wanted, dtype=str, keep_default_na=False, na_values=[],
                             chunksize=1_000_000):
        chunk = chunk[chunk['filepath'].str.contains(_LIKE_MZ_ML, regex=True, flags=re.DOTALL)]
        fingerprints = fingerprint_column(chunk)
        out = chunk[['usi', 'filepath', 'dataset']].copy()
        out['_fingerprint'] = fingerprints if fingerprints is not None else ''
        out['_size'] = chunk['size'] if 'size' in chunk.columns else ''
        out['_classification'] = chunk['classification'] if 'classification' in chunk.columns else ''
        parts.append(out)
    df = pd.concat(parts, ignore_index=True) if parts else pd.DataFrame(columns=['usi', 'filepath', 'dataset'])
    # the per-dataset query orders by usi (SQLite BINARY collation = code point order)
    df = df.sort_values('usi', kind='stable')
    FILES_TABLE = {dataset: group.reset_index(drop=True) for dataset, group in df.groupby('dataset', sort=False)}
    print(f"Dataset cache table {path}: {len(df)} mz(X)ML files in {len(FILES_TABLE)} datasets")


def _incomplete_answer(text, df, dataset):
    """Why a per-dataset answer of the dataset cache is incomplete, or None if it is complete.
    datasette ends a streamed CSV (_stream=on) whose later page fails (SQL time limit) or that
    exceeds max_csv_mb with HTTP 200 and its error message as the last line; a connection cut
    off mid-answer leaves a partial last line. Either would silently drop files of the dataset
    (and a NaN file path crashed the matching of every dataset)."""
    if 'filepath' not in df.columns or 'usi' not in df.columns:
        return f"no usi/filepath columns ({text[:200]!r})"
    if text and not text.endswith('\n'):
        return f"answer does not end with a complete line (last line {text.splitlines()[-1][:200]!r})"
    bad = df['filepath'].isna() | ~df['usi'].astype(str).str.startswith('mzspec:')
    if 'dataset' in df.columns:
        bad |= df['dataset'].astype(str) != dataset
    if bad.any():
        return f"{int(bad.sum())} malformed rows (first: {str(df.loc[bad].iloc[0].tolist())[:200]})"
    return None


def _fetch_dataset_files(dataset, retry_attempts=3):
    """mz(X)ML files of a dataset from the dataset cache, or None if they could not be fetched
    completely (then the dataset is in FETCH_FAILED and its rows are skipped, loudly)."""
    if FILES_TABLE is not None and dataset in FILES_TABLE:
        MATCH_STATS['datasets_from_files_table'] += 1
        return FILES_TABLE[dataset]
    # Not in the table (no mz(X)ML files when the export was made, or no table): ask the cache.
    # _stream=on returns all rows (without it at most max_returned_rows = 50000); the
    # timeout is short, as the cache answers within its 60 s SQL time limit or not at all
    url_f = "{}{}&filepath__endswith=%25.mz%25ML&_size=max&_stream=on".format(ccms_peak_link, dataset)
    print("Fetching {}".format(url_f))
    MATCH_STATS['datasets_queried'] += 1

    for attempt in range(retry_attempts):
        print(f"Checking dataset {dataset} by going to the dataset cache")
        try:
            dataset_files_response = requests.get(url_f, timeout=(15, 120))
            dataset_files_response.raise_for_status()
            text = dataset_files_response.text
            attempt_df = pd.read_csv(StringIO(text))
            problem = _incomplete_answer(text, attempt_df, dataset)
            if problem is None:
                return attempt_df
            MATCH_STATS['dataset_answers_incomplete'] += 1
            print(f"WARNING: attempt {attempt + 1}: incomplete dataset cache answer for {dataset}: {problem}")
        except (requests.RequestException, pd.errors.ParserError, pd.errors.EmptyDataError) as e:
            print(f"Attempt {attempt + 1} error: {e}")
        if attempt + 1 < retry_attempts:
            wait = 5 * 3 ** attempt
            print(f"Attempt {attempt + 1} failed. Retrying in {wait} seconds...")
            sleep(wait)
    MATCH_STATS['datasets_query_failed'] += 1
    FETCH_FAILED.add(dataset)
    print(f"ERROR: could not fetch the complete file list of {dataset} after {retry_attempts} attempts; "
          f"its metadata rows are not matched this run")
    return None


def _ws_key(name):
    """a file name with runs of whitespace as one space and none at either end"""
    return re.sub(r'\s+', ' ', str(name)).strip()


def _strip_segments(path):
    """a metadata file path without whitespace at either end of each folder / file name: sheets write
    ' 05_Mastigocladus_laminosus_SPE.mzML' (MSV000088697) for the listed '05_Mastigocladus_laminosus_SPE.mzML'"""
    return '/'.join(part.strip() for part in str(path).split('/'))


class DatasetFiles:
    """File lookups of one dataset (by name, name without .mzML/.mzXML, name with runs of whitespace
    as one space, lower-case name)."""

    def __init__(self, ccms_df):
        # (a row without a file path cannot be matched; answers with such rows are retried
        # by _fetch_dataset_files, this only guards against a crash)
        ccms_df = ccms_df[ccms_df["filepath"].notna()].reset_index(drop=True)
        ccms_df["query_path"] = ccms_df["filepath"].apply(lambda x: os.path.basename(x))
        # (folder names compared with runs of whitespace as one space, none at either end: _ws_key)
        ccms_df["_dir_parts"] = ccms_df["filepath"].apply(lambda p: tuple(_ws_key(d) for d in Path(p).parts[:-1]))
        self.df = ccms_df
        # plain lists (row i = position i) for the per-row candidate scoring in _select
        self.filepaths = ccms_df["filepath"].tolist()
        self.dir_parts = ccms_df["_dir_parts"].tolist()
        # identical-copy identity of each file (merge_metadata.file_identity; None without a fingerprint):
        # from the export (load_files_table) or the columns of a per-dataset answer
        none = pd.Series('', index=ccms_df.index)
        if '_fingerprint' in ccms_df.columns:
            fingerprints = ccms_df['_fingerprint']
            sizes = ccms_df['_size'] if '_size' in ccms_df.columns else none
            classes = ccms_df['_classification'] if '_classification' in ccms_df.columns else none
        else:
            fingerprints = fingerprint_column(ccms_df)
            sizes = ccms_df['size'] if 'size' in ccms_df.columns else none
            classes = ccms_df['classification'] if 'classification' in ccms_df.columns else none
        # (fingerprints that files of other names share at other sizes are the method's: none)
        self.identity = {} if fingerprints is None else \
            without_method_fixed({p: file_identity(p, f, z, k) for p, f, z, k in zip(self.filepaths, fingerprints, sizes, classes)})
        self.by_name = collections.defaultdict(list)
        self.by_stem = collections.defaultdict(list)
        self.by_ws = collections.defaultdict(list)
        self.by_ws_stem = collections.defaultdict(list)
        self.by_lower = collections.defaultdict(list)
        self.by_lower_stem = collections.defaultdict(list)
        # (folders without polarity words, file name) -> paths, built on first use (polarity_twin)
        self.twin_index = None
        for i, name in zip(ccms_df.index, ccms_df["query_path"]):
            name = str(name)
            self.by_name[name].append(i)
            self.by_ws[_ws_key(name)].append(i)
            self.by_lower[_ws_key(name).lower()].append(i)
            stem, ext = os.path.splitext(name)
            if ext.lower() in _MZ_EXTENSIONS:
                self.by_stem[stem].append(i)
                self.by_ws_stem[_ws_key(stem)].append(i)
                self.by_lower_stem[_ws_key(stem).lower()].append(i)

    def candidates(self, meta_basename):
        """(rule, candidate rows) for the first rule that finds any file. Whitespace at either end of
        the name is ignored (also by 'exact'); 'whitespace': the same name once runs of whitespace are
        one space, also with the other of .mzML/.mzXML and without whitespace before the extension
        (MSV000097243 'Babayani_DOC_33 .mzML' for the listed 'Babayani_DOC_33.mzML'); 'case' also
        ignores such whitespace. A name listed exactly as written (whitespace included) is 'exact'."""
        meta_basename = str(meta_basename)
        if not self.by_name.get(meta_basename):
            meta_basename = meta_basename.strip()
        meta_basename2 = meta_basename[:-3] + "ML" if len(meta_basename) >= 3 else meta_basename
        exact = self.by_name.get(meta_basename, []) + \
            ([] if meta_basename2 == meta_basename else self.by_name.get(meta_basename2, []))
        if exact:
            return 'exact', exact
        stem, ext = os.path.splitext(meta_basename)
        if ext.lower() in _MZ_EXTENSIONS:
            same_stem = self.by_stem.get(stem, [])
            if same_stem:
                return 'stem_ext', same_stem
        ws = self.by_ws.get(_ws_key(meta_basename), [])
        if ext.lower() in _MZ_EXTENSIONS:
            ws = list(dict.fromkeys(ws + self.by_ws_stem.get(_ws_key(stem), [])))
        if ws:
            return 'whitespace', ws
        lower = self.by_lower.get(_ws_key(meta_basename).lower(), [])
        if ext.lower() in _MZ_EXTENSIONS:
            lower = list(dict.fromkeys(lower + self.by_lower_stem.get(_ws_key(stem).lower(), [])))
        if lower:
            return 'case', lower
        return None, []


# Ties: several files fit a metadata row equally well (same name, same folder score, same
# preferred top folder). A tie is decided only when that is safe, in this order:
#   exact_path   the row names a path (with folders) that is one of the candidates, also when
#                runs of whitespace differ (MASST tables write 'ATCC13124 C.perfringens' for the
#                folder 'ATCC13124  C.perfringens')
#   polarity     every candidate sits in a positive- or negative-mode folder (a 'pos'/'neg'/
#                'positive'/'negative' word in its folders) and exactly one has the row's
#                polarity (IonizationSourceAndPolarity, else a pos/neg word in the row's file name).
#                If the row has no polarity and the candidates are one positive and one negative
#                file at otherwise the same path (the two acquisitions of the sample the row
#                describes), the row describes both: each gets a row ('each_mode'; its
#                IonizationSourceAndPolarity stays missing)
#   copies       all candidates are copies of one file: the same path once an 'updates/<id>/'
#                prefix is removed (a re-upload), format folders (raw, mzML, mzXML, ...) are
#                ignored and x.mzML / x.mzXML count as one (a conversion). Preference, in order:
#                the file named exactly as in the row, the original upload before an update's
#                copy, .mzML before .mzXML, then the path in code-point order
#   identical_copies  all candidates are the same data uploaded twice (or more) under other folder
#                names: the dataset cache describes their spectra identically and their sizes agree
#                (merge_metadata.identical_files; e.g. peak/Plate01_mzXMLs/x.mzXML and
#                peak/Plate02_mzXMLs/x.mzXML). Preference as for copies
#   each_acquisition  the candidates are other acquisitions of the sample the row describes: their
#                paths differ in one folder only, and only by acquisition words (_ACQUISITION_WORDS:
#                aqueous / organic extract, chromatography, MS1 / MS/MS run; format words such as raw
#                or mzML are ignored), e.g. Aq/x.mzML and Org/x.mzML, or Pos_MSMS_raw/x.mzML and
#                pos_raw/x.mzML. As for each_mode each file gets the row; the fields that describe
#                what differs (_ACQUISITION_FIELDS, e.g. SampleExtractionMethod for Aq/Org) are left
#                missing, as one row cannot tell the acquisitions apart
# Otherwise (e.g. plate01/x.mzML and plate02/x.mzML, or blanks of several batches) the row is
# not matched (logged as ambiguous, 'different_folders').
_UPDATE_PREFIX = re.compile(r'^updates/[^/]+/')
_FORMAT_FOLDER = re.compile(r'(?i)^(raws?|raw[ _-]?(data|files?)|mzx?mls?|mzx?ml[ _-]?(data|files?))$')
# (not 'positive_controls' / 'negative ctrl': control samples, not acquisition modes)
_POLARITY_WORD = re.compile(r'(?i)(?<![a-z])(pos|positive|neg|negative)(?![a-z])(?![ _-]*(?:controls?|ctrls?|ctls?)(?![a-z]))')


# acquisition words in folder names (whole tokens between '_', '-', ' ', ...) and what they describe
_ACQUISITION_WORDS = {'aq': 'extract', 'aqueous': 'extract', 'org': 'extract', 'organic': 'extract',
                      'hilic': 'chromatography', 'c18': 'chromatography', 'c8': 'chromatography', 'rp': 'chromatography',
                      'ms1': 'ms_level', 'ms2': 'ms_level', 'msms': 'ms_level', 'dda': 'ms_level', 'dia': 'ms_level',
                      'fullscan': 'ms_level'}
# the fields a sheet row gives for both acquisitions but can only be true for one of them
_ACQUISITION_FIELDS = {'extract': ['SampleExtractionMethod'], 'chromatography': ['ChromatographyAndPhase'], 'ms_level': []}
# file format words, not part of what a folder says about the acquisition (Pos_MSMS_raw vs pos_raw)
_FORMAT_WORDS = {'raw', 'raws', 'mzml', 'mzmls', 'mzxml', 'mzxmls', 'data', 'file', 'files'}
_TOKEN_SPLIT = re.compile(r'[^A-Za-z0-9]+')


def _norm_ws(path):
    return '/'.join(re.sub(r'\s+', ' ', part).strip() for part in str(path).split('/'))


def _acquisition_of_folder(folder):
    """(the folder's words without acquisition and format words, its acquisition words)"""
    words = [w.lower() for w in _TOKEN_SPLIT.split(folder) if w]
    return (tuple(w for w in words if w not in _ACQUISITION_WORDS and w not in _FORMAT_WORDS),
            frozenset(w for w in words if w in _ACQUISITION_WORDS))


def _acquisition_difference(paths):
    """What the paths differ in (a set of _ACQUISITION_FIELDS keys) if they are other acquisitions
    of one sample (see each_acquisition), else None"""
    split = [p.split('/') for p in paths]
    if len({len(parts) for parts in split}) != 1:
        return None
    if len({os.path.splitext(parts[-1])[0].lower() for parts in split}) != 1:
        return None
    differing = [i for i in range(len(split[0]) - 1) if len({parts[i] for parts in split}) > 1]
    if len(differing) != 1:
        return None
    folders = [_acquisition_of_folder(parts[differing[0]]) for parts in split]
    if len({rest for rest, _ in folders}) != 1 or len({acq for _, acq in folders}) != len(folders):
        return None
    words = frozenset().union(*(acq for _, acq in folders)) - frozenset.intersection(*(acq for _, acq in folders))
    return {_ACQUISITION_WORDS[w] for w in words} or None


def _identical_copies(paths, identity):
    """True if every pair of paths is an identical copy (merge_metadata.identical_files)"""
    records = [identity.get(p) for p in paths]
    return all(r is not None for r in records) and \
        all(identical_files(a, b) for i, a in enumerate(records) for b in records[i + 1:])


def _copy_key(path):
    parts = _UPDATE_PREFIX.sub('', path).split('/')
    dirs = tuple('<format>' if _FORMAT_FOLDER.match(d) else d for d in parts[:-1])
    stem, ext = os.path.splitext(parts[-1])
    return dirs, (stem if ext.lower() in _MZ_EXTENSIONS else parts[-1])


def _copy_preference(path, meta_basename):
    name = path.rsplit('/', 1)[-1]
    return (name != meta_basename, bool(_UPDATE_PREFIX.match(path)), not name.lower().endswith('.mzml'), path)


def _polarity_of(words):
    found = {w.lower()[:3] for w in words}
    return found.pop() if len(found) == 1 else None


def _path_polarity(path):
    """'pos' or 'neg' from the folders of a file path, or None."""
    return _polarity_of(w for d in path.split('/')[:-1] for w in _POLARITY_WORD.findall(d))


def _row_polarity(metadata_row, meta_path):
    """'pos' or 'neg' of a metadata row, or None."""
    value = str(metadata_row.get('IonizationSourceAndPolarity', '')).lower()
    if ('positive' in value) != ('negative' in value):
        return 'pos' if 'positive' in value else 'neg'
    return _polarity_of(_POLARITY_WORD.findall(os.path.basename(meta_path)))


# Polarity twins: a sheet row that states one polarity (IonizationSourceAndPolarity positive or
# negative) and is matched to a file in a positive- or negative-mode folder of that polarity, where
# the dataset also has a file of exactly the same name in the folder of the other polarity, the
# folders differing only by the polarity word (peak/Positive mzml/x.mzML and peak/Negative mzml/x.mzML,
# ccms_peak/raw_pos/Blank_1.mzML and ccms_peak/raw_neg/Blank_1.mzML: the other acquisition of the
# sample the row describes; MSV000093514, MSV000085779, checked by the scans' polarity cvParams).
# The other file gets a copy of the row with the polarity of its folder (add_polarity_twins), unless
# the sheet also describes it (a row matched to it, or a row naming its file name with its polarity),
# several different rows would give it a row, or the row comes from a fill-only file. Fields that this
# sheet gives differently for its positive and its negative rows (_TWIN_GUARD_FIELDS) are acquisition-
# specific and stay missing on the copy.
_TWIN_GUARD_FIELDS = ['ChromatographyAndPhase', 'MassSpectrometer', 'InternalStandardsUsed', 'SampleExtractionMethod']
# hidden columns of a matched row: the twin's path and the fields left missing on its copy
_TWIN_PATH, _TWIN_BLANK = '_polarity_twin_path', '_polarity_twin_blank'
_POLARITY_TERM = re.compile(r'\((positive|negative)\)$')


def _stated_polarity(value):
    """'pos' / 'neg' of an IonizationSourceAndPolarity term that states exactly one polarity, else None"""
    value = str(value).strip().lower()
    if ('positive' in value) == ('negative' in value):
        return None
    return 'pos' if 'positive' in value else 'neg'


def _other_polarity_term(value, allowed=None):
    """'electrospray ionization (positive)' -> 'electrospray ionization (negative)' (and back), None if the
    term states no single polarity or the other term is not an allowed term (allowed: a set, or None)"""
    value = str(value).strip()
    m = _POLARITY_TERM.search(value)
    if not m:
        return None
    other = value[:m.start()] + ('(negative)' if m.group(1) == 'positive' else '(positive)')
    return other if allowed is None or other in allowed else None


def _unpolar_dirs(path):
    """the folders of a path with every polarity word replaced by one placeholder"""
    return tuple(_POLARITY_WORD.sub('<polarity>', d) for d in path.split('/')[:-1])


def polarity_twin(files, path):
    """The one file of the dataset (DatasetFiles) with exactly the same file name as path in the folder
    of the other polarity (folders equal but for the polarity word), or None"""
    polarity = _path_polarity(path)
    if polarity is None:
        return None
    if files.twin_index is None:
        files.twin_index = collections.defaultdict(list)
        for p in files.filepaths:
            if _path_polarity(p) is not None:
                files.twin_index[(_unpolar_dirs(p), p.rsplit('/', 1)[-1])].append(p)
    others = [p for p in files.twin_index.get((_unpolar_dirs(path), path.rsplit('/', 1)[-1]), [])
              if p != path and _path_polarity(p) not in (None, polarity)]
    return others[0] if len(others) == 1 else None


def _polarity_specific_fields(group):
    """the _TWIN_GUARD_FIELDS a sheet's rows (one dataset of one metadata file) give with other values on
    its positive than on its negative rows"""
    polarity = group['IonizationSourceAndPolarity'].map(_stated_polarity) if 'IonizationSourceAndPolarity' in group.columns \
        else pd.Series(None, index=group.index)
    out = []
    for field in _TWIN_GUARD_FIELDS:
        if field not in group.columns:
            continue
        values = {pol: {str(v).strip() for v in group.loc[polarity == pol, field] if _present(v)} for pol in ('pos', 'neg')}
        if values['pos'] and values['neg'] and values['pos'] != values['neg']:
            out.append(field)
    return out


def _sheet_names_by_polarity(group):
    """(file name, stated polarity) of every row of a sheet (one dataset of one metadata file)"""
    if 'IonizationSourceAndPolarity' not in group.columns:
        return set()
    return {(os.path.basename(_strip_segments(f)), _stated_polarity(p))
            for f, p in zip(group['filename'].astype(str), group['IonizationSourceAndPolarity']) if _stated_polarity(p)}


def add_polarity_twins(frames, allowed_polarities=None):
    """frames: the matched tables after keep_fill_only_rows_of_described_files (rows of metadata files that
    are not fill-only carry the hidden _TWIN_PATH / _TWIN_BLANK columns) -> the same list with one more
    table of the polarity-twin rows (see _TWIN_GUARD_FIELDS). A twin is not added where a matched row
    already has its USI, or where rows that differ (other than in their file) would give it a row.
    allowed_polarities: the allowed IonizationSourceAndPolarity terms (a set) or None."""
    described = set()
    for df in frames:
        if 'USI' in df.columns:
            described.update(df['USI'].astype(str))
    proposals = collections.defaultdict(list)
    for df in frames:
        if _TWIN_PATH not in df.columns:
            continue
        for row in df[df[_TWIN_PATH].fillna('') != ''].to_dict('records'):
            twin = dict(row)
            path = twin.pop(_TWIN_PATH)
            blank = [f for f in str(twin.pop(_TWIN_BLANK, '') or '').split('|') if f]
            other = _other_polarity_term(twin.get('IonizationSourceAndPolarity', ''), allowed_polarities)
            if other is None:
                MATCH_STATS['rows_polarity_twin_no_term_of_other_polarity'] += 1
                continue
            twin['IonizationSourceAndPolarity'] = other
            for field in blank:
                if _present(twin.get(field)):
                    twin[field] = 'missing value'
                    MATCH_STATS[f'fields_left_missing_polarity_twin_{field}'] += 1
            twin['filename'] = 'f.' + path
            twin['USI'] = _make_usi_from_filename(twin['filename'], twin['ATTRIBUTE_DatasetAccession'])
            proposals[twin['USI']].append((twin, row.get('USI'), row.get('MetadataFileDataset', '')))
    rows = []
    for usi in sorted(proposals):
        twins = proposals[usi]
        if usi in described:
            MATCH_STATS['rows_polarity_twin_file_already_described'] += len(twins)
            MATCH_REPORT.append(('polarity_twin', usi.split(':')[1], ' ; '.join(sorted(str(s) for _, s, _ in twins)),
                                 'polarity_twin_file_already_described', usi.split(':', 2)[2]))
            continue
        if len({tuple(sorted((k, str(v)) for k, v in t.items())) for t, _, _ in twins}) > 1:
            MATCH_STATS['rows_polarity_twin_conflicting_rows'] += len(twins)
            MATCH_REPORT.append(('polarity_twin', usi.split(':')[1], usi, 'polarity_twin_conflicting_rows',
                                 ' ; '.join(sorted(str(s) for _, s, _ in twins))))
            continue
        twin, source_usi, _ = twins[0]
        MATCH_STATS['rows_matched_polarity_twin'] += 1
        MATCH_REPORT.append(('polarity_twin', twin['ATTRIBUTE_DatasetAccession'], str(source_usi), 'matched_polarity_twin',
                             twin['filename'][2:]))
        rows.append(twin)
    if rows:
        frames = frames + [pd.DataFrame(rows)]
    return frames


def _break_tie(paths, meta_path, polarity, identity=None):
    """(path, rule) when a tie between paths can be decided safely, else (None, kind); a tuple of
    paths (each gets the row) for each_mode and each_acquisition. identity: path -> file_identity
    record (DatasetFiles.identity), for identical_copies."""
    meta = str(meta_path).strip()
    if '/' in meta:
        same = [p for p in paths if p == meta]
        if len(same) == 1:
            return same[0], 'exact_path'
        same = [p for p in paths if _norm_ws(p) == _norm_ws(meta)]
        if len(same) == 1:
            return same[0], 'exact_path_whitespace'
    pool, rules = paths, []
    pols = [_path_polarity(p) for p in paths]
    if None not in pols and len(set(pols)) > 1:
        if polarity is None:
            unpolar = {_copy_key('/'.join(_POLARITY_WORD.sub('<polarity>', d) for d in p.split('/')[:-1]) + '/' + p.rsplit('/', 1)[-1])
                       for p in paths}
            if len(paths) == 2 and len(unpolar) == 1:
                return tuple(sorted(paths)), 'each_mode'
            return None, 'polarity_unknown'
        pool = [p for p, pol in zip(paths, pols) if pol == polarity]
        if not pool:
            return None, 'polarity_no_candidate'
        if len(pool) == 1:
            return pool[0], 'polarity'
        rules.append('polarity')
    if len({_copy_key(p) for p in pool}) == 1:
        return min(pool, key=lambda p: _copy_preference(p, os.path.basename(meta))), '_'.join(rules + ['copies'])
    if identity and _identical_copies(pool, identity):
        return min(pool, key=lambda p: _copy_preference(p, os.path.basename(meta))), '_'.join(rules + ['identical_copies'])
    if _acquisition_difference(pool):
        return tuple(sorted(pool)), '_'.join(rules + ['each_acquisition'])
    return None, 'different_folders'


def _select(files, rows, meta_path, polarity=None):
    """(best candidate file path or None, the equally good candidates, tie rule or None)."""
    meta_dir_parts = tuple(_ws_key(d) for d in Path(meta_path).parts[:-1])  # () if no dirs present

    # Score candidates by longest common directory suffix with the metadata path
    # (all score 0 when the metadata has no directory info); keep the best ones
    scores = [longest_common_dir_suffix_score(meta_dir_parts, files.dir_parts[i]) for i in rows]
    best_score = max(scores)
    found_file_paths = list(dict.fromkeys(files.filepaths[i] for i, score in zip(rows, scores) if score == best_score))

    # The same file can be present in several top-level folders (e.g. the
    # converted ccms_peak copy and the uploaded peak/raw copy); prefer them in
    # this order. Folders are compared as whole path segments.
    selected_paths = found_file_paths
    for preferred_dir in PREFERRED_DIRECTORIES:
        in_dir = [path for path in found_file_paths if preferred_dir in Path(path).parts[:-1]]
        if in_dir:
            selected_paths = in_dir
            break

    if len(selected_paths) > 1:
        # Several distinct files fit equally well (e.g. pos/s1.mzML and neg/s1.mzML for
        # metadata 's1.mzML'): decided only by the safe rules of _break_tie
        selected, how = _break_tie(selected_paths, meta_path, polarity, files.identity)
        return selected, selected_paths, how
    return selected_paths[0], selected_paths, None


def _dataset_files(dataset, files_cache):
    """DatasetFiles of a dataset (fetched once), or None if it has no mz(X)ML files or they
    could not be fetched (then the dataset is in FETCH_FAILED)."""
    if dataset not in files_cache:
        ccms_df = _fetch_dataset_files(dataset)
        if ccms_df is None:
            print(f"Error: could not fetch the file list for dataset {dataset}, skipping it")
            files_cache[dataset] = None
        elif len(ccms_df) == 0:
            print("Error: the size of the input is too small, length:", len(ccms_df))
            files_cache[dataset] = None
        else:
            print("Received a dataframe with {} rows.".format(len(ccms_df)))
            files_cache[dataset] = DatasetFiles(ccms_df)
    return files_cache[dataset]


def _lookup_path(files, filename):
    """The metadata path a row is matched with: as written when its file name is listed exactly so
    (whitespace at either end included), else without whitespace at either end of each folder /
    file name (_strip_segments)"""
    raw = str(filename)
    stripped = _strip_segments(raw)
    if raw != stripped and files.by_name.get(os.path.basename(raw)):
        return raw
    return stripped


def _match_row(files, metadata_row):
    """(rule, selected path or None, candidate paths, tie rule) for one metadata row."""
    meta_path = _lookup_path(files, metadata_row["filename"])
    rule, rows = files.candidates(os.path.basename(meta_path))
    if not rows:
        return None, None, [], None
    selected_path, options, how = _select(files, rows, meta_path, _row_polarity(metadata_row, meta_path))
    return rule, selected_path, options, how


def _mark_polarity_twin(files, metadata_row, selected_path, rule, sheet_names, twin_blank, sheet_dirs=frozenset()):
    """Sets the hidden _TWIN_PATH / _TWIN_BLANK of a matched row whose file has a polarity twin (see
    _TWIN_GUARD_FIELDS): the row states one polarity, the one of its file's folder, the sheet has no
    row naming the twin's file name with the twin's polarity, and no row of the sheet is matched to a file
    in the twin's folder (sheet_dirs: the folders of the sheet's matched files): a sheet that describes
    files of that folder leaves the others out on purpose (MSV000097243: the positive re-runs DOC_23r /
    DOC_27r are described, the first positive runs DOC_23 / DOC_27 are not)."""
    polarity = _stated_polarity(metadata_row.get('IonizationSourceAndPolarity', ''))
    if polarity is None or 'each_' in rule or _path_polarity(selected_path) != polarity:
        return
    twin = polarity_twin(files, selected_path)
    if twin is None:
        return
    if (twin.rsplit('/', 1)[-1], _path_polarity(twin)) in sheet_names:
        MATCH_STATS['rows_polarity_twin_sheet_names_the_other_file'] += 1
        MATCH_STATS[f'rows_polarity_twin_sheet_names_the_other_file {metadata_row.get("ATTRIBUTE_DatasetAccession")}'] += 1
        return
    if twin.rsplit('/', 1)[0] in sheet_dirs:
        MATCH_STATS['rows_polarity_twin_folder_described_by_the_sheet'] += 1
        MATCH_STATS[f'rows_polarity_twin_folder_described_by_the_sheet {metadata_row.get("ATTRIBUTE_DatasetAccession")}'] += 1
        return
    MATCH_STATS['rows_polarity_twin_proposed'] += 1
    metadata_row[_TWIN_PATH] = twin
    metadata_row[_TWIN_BLANK] = twin_blank


def _match_filenames_and_add_usi(dataset_metadata_df, file_accession=None, files_cache=None, source='', rehome=True,
                                 polarity_twins=False):
    """Matches each metadata row to a file of its own dataset (ATTRIBUTE_DatasetAccession;
    files can describe several datasets). Rows without a valid accession are matched in the
    dataset of the metadata file. A row whose stated dataset has no file of that name (a wrong
    MassiveID, e.g. MSV000090180..239 written by dragging MSV000090179 down a spreadsheet
    column) is matched in the metadata file's own dataset when that has the file and no other
    row of this metadata file is matched to it. Returns the matched rows with 'filename' set to
    the file path and a USI, or None if nothing could be matched. polarity_twins: mark the rows whose
    file has a polarity twin (_mark_polarity_twin; the rows are added by add_polarity_twins)."""
    files_cache = {} if files_cache is None else files_cache
    accessions = dataset_metadata_df['ATTRIBUTE_DatasetAccession'].astype(str).str.strip()
    valid = accessions.str.match(_ACCESSION)
    fallback = file_accession if file_accession and _ACCESSION.match(file_accession) else \
        (accessions[valid].iloc[0] if valid.any() else None)
    if (~valid).any():
        print(f"{(~valid).sum()} rows without a valid dataset accession; matching them in {fallback}")
        MATCH_STATS['rows_invalid_accession'] += int((~valid).sum())
        if fallback is None:
            dataset_metadata_df = dataset_metadata_df[valid]
            accessions = accessions[valid]
        else:
            accessions = accessions.where(valid, fallback)
            MATCH_STATS['rows_invalid_accession_set_to_file_dataset'] += int((~valid).sum())
    dataset_metadata_df = dataset_metadata_df.copy()
    dataset_metadata_df['ATTRIBUTE_DatasetAccession'] = accessions
    # the dataset of the metadata file itself (from its name), where rows with a wrong
    # accession are looked for
    own_dataset = file_accession if file_accession and _ACCESSION.match(file_accession) else None

    output_row_list = []
    # rows not found in their stated dataset, to look for in own_dataset: (row, result)
    elsewhere = []
    # file paths of own_dataset that rows of this metadata file are matched to
    own_paths = set()
    for dataset, group in dataset_metadata_df.groupby('ATTRIBUTE_DatasetAccession', sort=False):
        print(f"Checking dataset {dataset} ({len(group)} rows)")
        files = _dataset_files(dataset, files_cache)
        retry_elsewhere = rehome and own_dataset is not None and dataset != own_dataset and dataset not in FETCH_FAILED
        if files is None:
            if retry_elsewhere:
                elsewhere.extend((row, 'dataset_without_files') for row in group.to_dict('records'))
                continue
            MATCH_STATS['rows_dataset_without_files'] += len(group)
            MATCH_REPORT.extend((source, dataset, f, 'dataset_without_files', '') for f in group['filename'].astype(str))
            continue

        # polarity twins (add_polarity_twins): what this sheet says per polarity
        twin_blank = '|'.join(_polarity_specific_fields(group)) if polarity_twins else ''
        sheet_names = _sheet_names_by_polarity(group) if polarity_twins else set()
        matched = []
        for metadata_row in group.to_dict('records'):
            meta_path = str(metadata_row["filename"])
            meta_basename = os.path.basename(_lookup_path(files, meta_path))
            rule, selected_path, options, how = _match_row(files, metadata_row)
            if rule is None:
                if retry_elsewhere:
                    elsewhere.append((metadata_row, 'not_found'))
                    continue
                MATCH_STATS['rows_not_found'] += 1
                MATCH_REPORT.append((source, dataset, meta_path, 'not_found', ''))
                continue
            if selected_path is None:
                print(f"Ambiguous match for {meta_path} ({how}), skipping: {options}")
                MATCH_STATS[f'rows_ambiguous_{rule}'] += 1
                MATCH_STATS[f'rows_ambiguous_{how}'] += 1
                AMBIGUOUS.append((dataset, meta_path, rule, ' ; '.join(options)))
                MATCH_REPORT.append((source, dataset, meta_path, f'ambiguous_{rule}', ' ; '.join(options)))
                continue
            if meta_path != _lookup_path(files, meta_path):
                MATCH_STATS['rows_filename_whitespace_stripped'] += 1
            base_rule = rule
            # each_acquisition: the fields that tell the acquisitions apart stay missing
            blank = [f for kind in sorted(_acquisition_difference(selected_path) or ())
                     for f in _ACQUISITION_FIELDS[kind]] if how and how.endswith('each_acquisition') else []
            for path in (selected_path if isinstance(selected_path, tuple) else (selected_path,)):
                rule = base_rule
                if rule == 'exact' and os.path.basename(path) != meta_basename:
                    rule = 'exact_mzxml_to_converted_mzml'
                if how:
                    rule = f'{rule}_tie_{how}'
                row = dict(metadata_row)
                for field in blank:
                    if field in row and str(row[field]).strip().lower() not in ('', 'nan', 'missing value'):
                        row[field] = 'missing value'
                        MATCH_STATS[f'fields_left_missing_each_acquisition_{field}'] += 1
                matched.append((row, path, rule, os.path.basename(path) == meta_basename, meta_path))

        # A file named exactly by one row is not also given to rows that only reach it
        # through another extension (x.mzXML and x.mzML listed for the same sample: the
        # x.mzXML row would otherwise also get the converted ccms_peak/x.mzML)
        exact_targets = {path for _, path, _, is_exact, _ in matched if is_exact}
        # the folders of the files this sheet describes (_mark_polarity_twin)
        sheet_dirs = {path.rsplit('/', 1)[0] for _, path, _, is_exact, _ in matched
                      if '/' in path and (is_exact or path not in exact_targets)} if polarity_twins else set()
        for metadata_row, selected_path, rule, is_exact, meta_path in matched:
            if not is_exact and selected_path in exact_targets:
                MATCH_STATS['rows_twin_of_exactly_named_file'] += 1
                MATCH_REPORT.append((source, dataset, meta_path, 'twin_of_exactly_named_file', selected_path))
                continue
            MATCH_STATS[f'rows_matched_{rule}'] += 1
            MATCH_REPORT.append((source, dataset, meta_path, f'matched_{rule}', selected_path))
            if polarity_twins:
                _mark_polarity_twin(files, metadata_row, selected_path, rule, sheet_names, twin_blank, sheet_dirs)
            metadata_row["filename"] = "f." + selected_path
            metadata_row["USI"] = _make_usi_from_filename(
                metadata_row["filename"],
                metadata_row["ATTRIBUTE_DatasetAccession"]
            )
            output_row_list.append(metadata_row)
            if dataset == own_dataset:
                own_paths.add(selected_path)

    # rows whose stated dataset has no file of that name: the metadata file's own dataset
    files = _dataset_files(own_dataset, files_cache) if elsewhere else None
    for metadata_row, result in elsewhere:
        dataset = metadata_row['ATTRIBUTE_DatasetAccession']
        meta_path = str(metadata_row["filename"])
        rule, selected_path, options, how = _match_row(files, metadata_row) if files is not None else (None, None, [], None)
        if isinstance(selected_path, tuple):
            selected_path = None  # (one row for the files of both modes: not for a row with a wrong accession)
        if selected_path is None or selected_path in own_paths:
            MATCH_STATS[f'rows_{result}'] += 1
            if selected_path is not None:
                MATCH_STATS['rows_wrong_accession_file_already_described'] += 1
            MATCH_REPORT.append((source, dataset, meta_path, result, ''))
            continue
        own_paths.add(selected_path)
        MATCH_STATS['rows_matched_in_metadata_file_dataset'] += 1
        MATCH_REPORT.append((source, own_dataset, meta_path, f'matched_in_metadata_file_dataset_not_{dataset}', selected_path))
        metadata_row['ATTRIBUTE_DatasetAccession'] = own_dataset
        # the harmonizer built UniqueSubjectID as <stated accession>_<subject>
        subject = metadata_row.get('UniqueSubjectID')
        if isinstance(subject, str) and subject.startswith(str(dataset) + '_'):
            metadata_row['UniqueSubjectID'] = own_dataset + subject[len(str(dataset)):]
            MATCH_STATS['rows_matched_in_metadata_file_dataset_subject_id_fixed'] += 1
        metadata_row["filename"] = "f." + selected_path
        metadata_row["USI"] = _make_usi_from_filename(metadata_row["filename"], own_dataset)
        output_row_list.append(metadata_row)

    if not output_row_list:
        return None
    return pd.DataFrame(output_row_list)

def _match_metadata_file(file, gnps_column_names, files_cache):
    """The matched rows of one metadata file, or None."""
    print("Working on ", file)
    dataset_metadata_df = pd.read_csv( file , delimiter='\t')

    print(f"reading potential ReDU data with {len(dataset_metadata_df)} rows.")

    #Renaming the coloumn, Matching common columns and rearranging them in same order to final file
    dataset_metadata_df = dataset_metadata_df.rename(columns={'MassiveID': 'ATTRIBUTE_DatasetAccession'})
    common_cols = list(set(gnps_column_names).intersection(set(dataset_metadata_df.columns)))
    dataset_metadata_df = dataset_metadata_df.loc[:, common_cols]
    try:
        dataset_metadata_df = dataset_metadata_df[gnps_column_names]
    except KeyError as e:
        print(f"Skipping file {file} due to missing columns: {e}")
        return None

    # Matching the metadata, each row within its own dataset
    file_accession = re.match(r'([A-Z]+\d+)', os.path.basename(str(file)))
    MATCH_STATS['rows_in'] += len(dataset_metadata_df)
    enriched_metadata_df = _match_filenames_and_add_usi(
        dataset_metadata_df, file_accession.group(1) if file_accession else None, files_cache,
        source=os.path.basename(str(file)), rehome=not os.path.basename(str(file)).endswith(FILL_ONLY_MARK),
        polarity_twins=not os.path.basename(str(file)).endswith(FILL_ONLY_MARK))
    if enriched_metadata_df is not None:
        enriched_metadata_df['MetadataFileDataset'] = file_accession.group(1) if file_accession else ''

        print(f"Selected {len(enriched_metadata_df)} rows with columns: {list(enriched_metadata_df.columns)}")
    return enriched_metadata_df


def main():
    # parsing args
    parser = argparse.ArgumentParser(description='GNPS Name Matcher')
    parser.add_argument('passed_file_names', help='Input TSV file')
    parser.add_argument('metadata_folder')
    parser.add_argument('output_filename')
    parser.add_argument('path_allowed_terms_json')
    parser.add_argument('--match_report', default=None, help='optional TSV with the matching result of every metadata row')
    parser.add_argument('--files_table', default=None,
                        help='local copy of https://datasetcache.gnps2.org/dataset/uniquemri; without it (or if it '
                             'cannot be read, e.g. empty after a failed download) every dataset is queried from the dataset cache')

    args = parser.parse_args()
    if args.files_table:
        try:
            load_files_table(args.files_table)
        except Exception as e:
            print(f"WARNING: could not read the dataset cache table {args.files_table} ({type(e).__name__}: {e}); "
                  f"querying the dataset cache per dataset")
     
    ccms_filenames = collections.defaultdict(set)
    
    if args.passed_file_names == 'all':
        # sorted: the output must not depend on the order the file system lists the folder in
        passed_file_names = sorted(glob.glob(f"{args.metadata_folder}/*.tsv"))
    elif args.passed_file_names == 'single':
        passed_file_names = []
    else:
        # Read the TSV file and specify the delimiter as a tab
        df = pd.read_csv(args.passed_file_names, delimiter='\t', header=None, names=['Name'])
        # Extract the names from a specific column (e.g., column 'Name')
        passed_file_names = df['Name'].tolist()

    
    with open(args.path_allowed_terms_json, 'r') as file:
        allowed_terms_json = json.load(file)

    # (in the order of the allowed terms, not of a set: that would depend on PYTHONHASHSEED)
    gnps_column_names = ["ATTRIBUTE_DatasetAccession"] + [k for k in allowed_terms_json if k not in ('USI', 'MassiveID')]

    print("echo Iterating though rows now")
    
    all_metadata_list = []
    # file lists of the datasets, fetched once (a dataset can be described by several files)
    files_cache = {}

    failed_files = []
    for file in passed_file_names:
        # one metadata file that cannot be read or matched (e.g. an unexpected answer of the
        # dataset cache) only loses its own rows, not those of every other file
        try:
            enriched_metadata_df = _match_metadata_file(file, gnps_column_names, files_cache)
        except Exception as e:
            failed_files.append(os.path.basename(str(file)))
            MATCH_STATS['metadata_files_failed'] += 1
            print(f"ERROR: matching {file} failed ({type(e).__name__}: {e}); its rows are not in the output")
            traceback.print_exc(file=sys.stdout)
            continue
        if enriched_metadata_df is not None:
            all_metadata_list.append((os.path.basename(str(file)).endswith(FILL_ONLY_MARK), enriched_metadata_df))
    try:
        all_metadata_list = keep_fill_only_rows_of_described_files(all_metadata_list)
    except Exception as e:
        print(f"ERROR: the fill-only step failed ({type(e).__name__}: {e}); fill-only rows are left out")
        traceback.print_exc(file=sys.stdout)
        all_metadata_list = [df for fill_only, df in all_metadata_list if not fill_only]
    # polarity twins: the file of the other polarity of a described file gets the row (add_polarity_twins);
    # on an error none are added
    try:
        allowed_polarities = set((allowed_terms_json.get('IonizationSourceAndPolarity') or {}).get('allowed_values') or []) or None
        all_metadata_list = add_polarity_twins(all_metadata_list, allowed_polarities)
    except Exception as e:
        print(f"ERROR: the polarity-twin step failed ({type(e).__name__}: {e}); no polarity-twin rows are added")
        traceback.print_exc(file=sys.stdout)
    all_metadata_list = [df for df in all_metadata_list if len(df)]

    # Create a DataFrame from the list with headers (header only if nothing matched,
    # e.g. no new MASST rows)
    if all_metadata_list:
        merged_metadata_df = pd.concat(all_metadata_list)
    else:
        print("No metadata rows matched to files")
        merged_metadata_df = pd.DataFrame(columns=gnps_column_names + gnps_column_names_added)
    merged_metadata_df = merged_metadata_df[gnps_column_names + gnps_column_names_added]

    print("Matching summary (rows):")
    for key, value in sorted(MATCH_STATS.items()):
        print(f"  {key}: {value}")
    for dataset, meta_path, rule, options in AMBIGUOUS:
        print(f"  ambiguous ({rule}) {dataset} {meta_path}: {options}")

    if FETCH_FAILED:
        print(f"ERROR: the file lists of {len(FETCH_FAILED)} datasets could not be fetched; their rows are not matched: "
              f"{', '.join(sorted(FETCH_FAILED))}")
    if failed_files:
        print(f"ERROR: {len(failed_files)} metadata files failed and are not in the output: {', '.join(failed_files)}")

    if args.match_report:
        pd.DataFrame(MATCH_REPORT, columns=['source_file', 'dataset', 'filename', 'result', 'file_path']).to_csv(
            args.match_report, sep='\t', index=False)

    # Save the DataFrame to a TSV file without column names
    merged_metadata_df.to_csv(args.output_filename, sep='\t', index=False)


if __name__ == '__main__':
    main()