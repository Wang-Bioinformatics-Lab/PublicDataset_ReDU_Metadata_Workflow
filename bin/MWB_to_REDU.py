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
from blank_qc_names import blank_qc_type_from_filename, name_tokens, not_a_pool_by_name, not_a_pool_by_value, NOT_POOL_QC_NAME_WORDS
import reference_materials
import traceback
import time
from REDU_conversion_functions import get_taxonomy_id_from_name__allowedTerms
from REDU_conversion_functions import map_instrument_to_allowed
from REDU_conversion_functions import convert_smoking, bmi_to_numeric, convert_diet
from REDU_conversion_functions import food_ontology_term, food_sample_marker, food_term_of_organism
from read_and_validate_redu_from_github import complete_and_fill_REDU_table
from REDU_conversion_functions import age_category
from REDU_conversion_functions import get_taxonomy_info
from REDU_conversion_functions import normalize_redu_datetime
from REDU_conversion_functions import multi_sample_type, ORGANISM_LIST
# several organisms in one NCBITaxonomy value ('ID|Name||ID|Name')
from ncbi_database import all_parts, is_multi, join_parts, single_taxid_series, taxonomy_parts
from ms_instrument_matching import match_mass_spectrometer
from water_kind import water_kind_of_unstated_water
import REDU_conversion_functions as _rcf
from getAllWorkbench_file_paths import fetch_archive_listing, listing_query, ArchiveListingError, collapse_vendor_folder
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
    # word starts only: 'ceramides', 'pantothenamide', 'nicotinamide' are no HILIC columns
    'hilic': re.compile(r'(?:^|[^a-z])(?:hilic|amide|zic)'),
    'reverse': re.compile(r'(?:^|[^a-z])(?:c18|rp|rplc|reverse|reversed|c8|c30|t3)(?:[^a-z]|pos|neg|$)'),
}


def _path_tokens(path):
    """The polarities and chromatographies tokens of a path name (sets; both polarities in
    'pos_neg/x.raw')."""
    p = str(path).lower()
    pols = set()
    if _POS_RE.search(p) or re.search(r'(?<![a-z])(?:esi|apci|hesi)\s*\+', p):
        pols.add('positive')
    if _NEG_RE.search(p) or re.search(r'(?<![a-z])(?:esi|apci|hesi)\s*-(?![a-z])', p):
        pols.add('negative')
    return pols, [k for k, rx in _CHROM_TOKENS.items() if rx.search(p)]


def _path_indicates(path):
    """Polarity ('positive'/'negative'/None) and chromatography ('hilic'/'reverse'/None)
    stated by tokens in a raw file path, e.g. 'C18-NEG/0001.raw' or 'x_HILICpos_1.mzML'."""
    pols, chrom_hits = _path_tokens(path)
    polarity = next(iter(pols)) if len(pols) == 1 else None
    chrom = chrom_hits[0] if len(chrom_hits) == 1 else None
    return polarity, chrom


# archive names of one listed path in the 'archive' column of the file listing
ARCHIVE_SEP = '\n'


def _archive_names(archives):
    if archives is None or (isinstance(archives, float) and np.isnan(archives)):
        return []
    if isinstance(archives, str):
        archives = archives.split(ARCHIVE_SEP)
    return sorted({str(a) for a in archives if str(a)})


def archive_context(archives):
    """(polarity, chromatography) the Workbench archive name(s) of a listed path state
    ('ST004470_HILIC_Rawfiles.zip', 'ST003799_HILIC-NEG_mzml.zip'): per dimension the value
    every archive the path lies in states, else None (a path in a pos and a neg archive is
    one USI for both files)."""
    return _archive_context(tuple(_archive_names(archives)))


_ARCHIVE_CTX_CACHE = {}
# in archive names: a range of sample names ('C1-C8': the control mice C1 to C8, 'QC1-QC8') and a
# column word followed by digits ('C3015', 'C3002': project codes, not a C30 column) name no column
_ARCHIVE_RANGE = re.compile(r'(?<![a-z0-9])([a-z]+)\d+\s*-\s*\1\d+(?![0-9])')
_ARCHIVE_CODE = re.compile(r'(?<![a-z])(?:c18|c8|c30|t3)(?=\d)')


def _archive_text(name):
    """An archive name (lower case) without sample-name ranges and codes that the path tokens would
    read as a column (_ARCHIVE_RANGE, _ARCHIVE_CODE)"""
    text = _ARCHIVE_RANGE.sub(' ', str(name).lower())
    return _ARCHIVE_CODE.sub(' ', text)


def _archive_context(names):
    if names not in _ARCHIVE_CTX_CACHE:
        if not names:
            _ARCHIVE_CTX_CACHE[names] = (None, None)
        else:
            ctx = [_path_indicates(_archive_text(a)) for a in names]
            pols, chroms = {c[0] for c in ctx}, {c[1] for c in ctx}
            _ARCHIVE_CTX_CACHE[names] = ((next(iter(pols)) if len(pols) == 1 else None),
                                         (next(iter(chroms)) if len(chroms) == 1 else None))
    return _ARCHIVE_CTX_CACHE[names]


def archive_dimensions(archives_of_files):
    """The dimensions (0 polarity, 1 chromatography) in which the archive contexts of a
    study's listed files differ (HILIC and RP archives; pos and neg archives; a HILIC archive
    next to one stating nothing): only there do archive names tell its files apart. Archive
    names that all state the same (one 'x_HILIC.zip', 'x_C3015_Serum.zip' and
    'x_C3015_Serum-QC.zip') say nothing about which file is which acquisition."""
    values = {archive_context(a) for a in archives_of_files}
    return frozenset(d for d in (0, 1) if len({v[d] for v in values}) > 1)


def file_context(path, archives=None, dims=(0, 1)):
    """(polarity, chromatography) of a listed file: what its path states (_path_indicates);
    a dimension the path has no token of at all is taken from the archive name(s) the file
    lies in (archive_context), for the dimensions dims (archive_dimensions of the study). The
    path keeps precedence: the archive is not used for a file whose path states a dimension
    otherwise than the archive (Blank_POS1.mzML in x_Negative_Rawdata.zip; x_1C8.mzML, a well
    C8, in x_HILIC_pos.zip) or names both values of it, nor where the path states a dimension
    only in its file name and the archive does not state it too: a column or polarity word
    in a file name is often a sample or well name ('C8_194.d', 'X0141-M001A-103_1C8.mzML'
    in x_ionexchange-neg.zip) and must not decide the analysis together with the archive."""
    pols, chrom_hits = _path_tokens(path)
    polarity = next(iter(pols)) if len(pols) == 1 else None
    chrom = chrom_hits[0] if len(chrom_hits) == 1 else None
    if not ((not pols and 0 in dims) or (not chrom_hits and 1 in dims)) or not _archive_names(archives):
        return polarity, chrom
    a_ctx = archive_context(archives)
    folder = None
    for d, tokens, value in ((0, pols, polarity), (1, chrom_hits, chrom)):
        if not tokens:
            continue
        if value is None or (a_ctx[d] is not None and a_ctx[d] != value):
            return polarity, chrom
        if a_ctx[d] is None:
            folder = folder if folder is not None else _path_indicates(str(path).replace('\\', '/').rsplit('/', 1)[0]
                                                                       if '/' in str(path).replace('\\', '/') else '')
            if folder[d] != value:
                return polarity, chrom
    polarity = a_ctx[0] if not pols and 0 in dims else polarity
    chrom = a_ctx[1] if not chrom_hits and 1 in dims else chrom
    return polarity, chrom


def listing_contexts(raw_file_name_df, path_col='filename', archive_col='archive'):
    """{listed path: file_context} of a study's file listing (archive names in archive_col,
    used in the study's archive_dimensions; path contexts only without the column)."""
    paths = raw_file_name_df[path_col].astype(str).tolist()
    archives = raw_file_name_df[archive_col].tolist() if archive_col in raw_file_name_df.columns else [None] * len(paths)
    dims = archive_dimensions(archives)
    return {p: file_context(p, a, dims) for p, a in zip(paths, archives)}


def same_name_copy_context(path, archives=None, dims=(0, 1)):
    """(polarity, chromatography) that tells a same-name copy of a matched file apart from it
    (propagate_to_same_name_files: Method_2/x.raw, a copy of the sample matched as Method_1/x.raw).
    Chromatography: what the copy's folders state, with the dimension taken from its archive name(s)
    where they leave it open, as in file_context; a column word in the file name says nothing about
    which acquisition the copy is ('08C8.mzML', the well C8 of a HILIC study, in a POS and a NEG
    folder). Polarity: the whole path, as for any listed file (file_context): a polarity word in the
    copy's own file name is evidence ('lcqtof_t3_neg_1.mzML' is no positive acquisition)."""
    p = str(path).replace('\\', '/')
    folder = p.rsplit('/', 1)[0] + '/' if '/' in p else ''
    return file_context(p, archives, dims)[0], file_context(folder, archives, dims)[1]


def listing_copy_contexts(raw_file_name_df, path_col='filename', archive_col='archive'):
    """{listed path: same_name_copy_context} of a study's file listing (archive names in
    archive_col, used in the study's archive_dimensions as in listing_contexts)."""
    paths = raw_file_name_df[path_col].astype(str).tolist()
    archives = raw_file_name_df[archive_col].tolist() if archive_col in raw_file_name_df.columns else [None] * len(paths)
    dims = archive_dimensions(archives)
    return {p: same_name_copy_context(p, a, dims) for p, a in zip(paths, archives)}


def _context_contradicts(context, hints):
    """True if a file path's (polarity, chromatography) names another polarity or
    chromatography than the analysis hints ({'polarity': .., 'chrom': ..})."""
    polarity, chrom = context
    return bool((polarity is not None and hints.get('polarity') is not None and polarity != hints['polarity'])
                or (chrom is not None and hints.get('chrom') is not None and chrom != hints['chrom']))


def acquisition_contexts(stems, paths, hints=None, ranks=None, archives=None):
    """Acquisition context of each file (lists aligned with stems / paths / ranks / archives), for
    keeping one file per sample and acquisition: files of one stem (sample name) are raw and
    converted copies of one acquisition unless their paths state different polarities or
    chromatographies (temp_pos/x.mzML vs temp_neg/x.mzXML, rp_pos/x.mzML vs hilic_pos/x.mzML),
    which are separate acquisitions. Returns per file a group number within its stem (0 for
    a stem whose files agree, i.e. one file per stem as before) or -1 for a file left out:
    a file whose path contradicts the analysis' own polarity / chromatography while another
    file of the stem does not; in a stem whose files lie in several stated contexts, a file
    that states less of the analysis' context than another of at least as good a format
    (GT_H_Negative/x vs GT_RP_Negative/x for a reverse-phase negative analysis); and -
    where a stem is split by context - a file whose path names no context in the split
    dimension (an unattributable copy). A file whose path states no context at all, in a stem whose other files lie in
    stated contexts, is kept as a group of its own only when it is a better format (rank,
    mzML 0 < mzXML 1 < vendor 2) than every file kept for the contexts (a converted x.mzML
    next to hilic_p/x.raw, rplc_n/x.raw, ...). A file's context is that of its path, with
    a dimension the path leaves open taken from its Workbench archive name(s) where the
    study's archives differ in it (file_context, archive_dimensions; paths / archives: the
    study's whole listing): mzML_workbench/601.mzML of ST004470_HILIC_Rawfiles.zip is HILIC,
    MzML_workbench/601.mzML of ST004470_RP_Rawfiles.zip reverse phase."""
    hints = hints or {}
    ranks = list(ranks) if ranks is not None else [0] * len(stems)
    archives = list(archives) if archives is not None else [None] * len(stems)
    dims = archive_dimensions(archives)
    ctx = [file_context(p, a, dims) for p, a in zip(paths, archives)]
    by_stem = {}
    for i, s in enumerate(stems):
        by_stem.setdefault(s, []).append(i)
    group = [0] * len(stems)

    def _agreement(c):
        return int(c[0] is not None and c[0] == hints.get('polarity')) + int(c[1] is not None and c[1] == hints.get('chrom'))

    for s, idx_all in by_stem.items():
        if len(idx_all) < 2:
            continue
        idx = idx_all
        stated = {ctx[i] for i in idx} - {(None, None)}
        keep = [i for i in idx if not _context_contradicts(ctx[i], hints)]
        # files left that state nothing in the dimension the left-out ones contradict and
        # are another format of them (x.wiff / x.raw / x.d next to positive/x.mzXML for a
        # negative analysis) are a raw / converted copy of that other acquisition, not this
        # analysis' file: the stem is left to the one-file-per-stem rule and the path
        # conflict resolution, as before. Same-format files (postive/x.raw next to
        # negative/x.raw) stay: they are the other acquisition in an unrecognised folder.
        _dims = {d for i in idx if i not in keep for d, h in ((0, 'polarity'), (1, 'chrom'))
                 if ctx[i][d] is not None and hints.get(h) is not None and ctx[i][d] != hints[h]}
        if keep and _dims and all(ctx[i][d] is None for i in keep for d in _dims):
            _out_ext = {os.path.splitext(paths[i])[1].lower() for i in idx if i not in keep}
            if all(os.path.splitext(paths[i])[1].lower() not in _out_ext for i in keep):
                keep = []
        if keep and len(stated) > 1:
            best = max(_agreement(ctx[i]) for i in keep)
            top = [i for i in keep if _agreement(ctx[i]) == best]
            # only where that does not give up a better format: x.mzML (HILIC) next to
            # HILIC_Neg/x.wiff is its conversion, not another acquisition
            if best > 0 and min(ranks[i] for i in top) <= min(ranks[i] for i in keep):
                keep = top
        if keep and len(keep) < len(idx):
            for i in set(idx) - set(keep):
                group[i] = -1
            idx = keep
        pols = {ctx[i][0] for i in idx} - {None}
        chroms = {ctx[i][1] for i in idx} - {None}
        if len(pols) > 1 or len(chroms) > 1:
            split_pol, split_chrom = len(pols) > 1, len(chroms) > 1
            keys = {}
            for i in idx:
                pol, chrom = ctx[i]
                if (split_pol and pol is None) or (split_chrom and chrom is None):
                    group[i] = -1
                    continue
                keys[i] = (pol if split_pol else '', chrom if split_chrom else '')
            order = {k: n for n, k in enumerate(sorted(set(keys.values())))}
            for i, k in keys.items():
                group[i] = order[k]
        if len(stated) > 1:
            # the unattributed copies (no context in the path) left out above: the best of
            # them stays as its own group if it is a better format than all kept files
            kept = [i for i in idx_all if group[i] >= 0 and ctx[i] != (None, None)]
            loose = sorted((ranks[i], i) for i in idx_all if ctx[i] == (None, None))
            if kept:
                for i in (j for _, j in loose):
                    group[i] = -1
                if loose and loose[0][0] < min(ranks[i] for i in kept):
                    group[loose[0][1]] = max(group[i] for i in kept) + 1
    return group


def contexts_told_apart_by_listed_folders(stems, paths, groups, listed_names):
    """Acquisition groups (acquisition_contexts) with the files of one stem and group that the samples' listed raw
    file names with folders tell apart ('Seed_check/49.D/' and 'Leaf_check/49.D/' of two samples) in groups of their
    own: they are files of different samples, not copies of one file of which one is kept. A file of such a stem
    whose folders none or several of those listed names give is left out (-1: a copy of one of them, which one is
    not known). listed_names: the raw file names the samples list. Other files keep their group."""
    by_last = {}
    for n in listed_names:
        parts = [p for p in re.split(r'[\\/]+', str(n).strip()) if p]
        if len(parts) < 2:
            continue
        last, folder = _match_tokens(parts[-1].split('.')[0]), _match_tokens(' '.join(parts[:-1]))
        if last and folder:
            by_last.setdefault(last, set()).add(folder)
    groups = list(groups)
    if not by_last:
        return groups
    by_key, top = {}, {}
    for i, (st, g) in enumerate(zip(stems, groups)):
        if g >= 0:
            by_key.setdefault((st, g), []).append(i)
            top[st] = max(top.get(st, 0), g)
    for (st, g), idx in sorted(by_key.items(), key=lambda x: (str(x[0][0]), x[0][1])):
        folders = by_last.get(_match_tokens(st))
        if len(idx) < 2 or not folders or len(folders) < 2:
            continue
        named, loose = {}, []
        for i in idx:
            f = _path_folder_tokens(paths[i])[0]
            hit = [x for x in sorted(folders) if _token_run_in(x, f)]
            if len(hit) == 1:
                named.setdefault(hit[0], []).append(i)
            else:
                loose.append(i)
        if len(named) < 2:
            continue
        for n, x in enumerate(sorted(named)):
            if n:
                top[st] += 1
            for i in named[x]:
                groups[i] = g if n == 0 else top[st]
        for i in loose:
            groups[i] = -1
    return groups


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


_REST_BLANK = {'', 'na', 'n/a', 'nan', 'none', 'null', '-', 'not applicable'}


def _rest_blank(value):
    return value is None or (isinstance(value, float) and np.isnan(value)) or str(value).strip().lower() in _REST_BLANK


def chromatography_from_rest_analysis(mwtab, rest_analysis):
    """The mwTab of an analysis whose CHROMATOGRAPHY section states neither a CHROMATOGRAPHY_TYPE
    nor a COLUMN_NAME (an empty section: ST004691 AN007923), with the CHROMATOGRAPHY_TYPE of the
    study's REST analysis record of the same analysis (/rest/study/study_id/<ST>/analysis:
    'chromatography_type'), where that record has one ('Reversed phase' -> reverse phase (NOS)).
    Not the record's 'column_name': it can be that of another setup of the study's series
    (AN007923 of the nano-flow upload names the column of the high-flow sister studies). Only
    fills: a mwTab stating either field is returned as it is, also where the REST record says
    otherwise. The mwTab passed in is not modified."""
    if not isinstance(mwtab, dict) or not isinstance(rest_analysis, dict):
        return mwtab
    # the REST record must be the one of this mwTab's analysis
    an = _section(mwtab, 'METABOLOMICS WORKBENCH').get('ANALYSIS_ID')
    if _rest_blank(rest_analysis.get('analysis_id')) or (not _rest_blank(an) and str(an) != str(rest_analysis.get('analysis_id'))):
        return mwtab
    ch = _section(mwtab, 'CHROMATOGRAPHY')
    if not (_rest_blank(ch.get('CHROMATOGRAPHY_TYPE')) and _rest_blank(ch.get('COLUMN_NAME'))):
        return mwtab
    fill = {k: str(rest_analysis.get(r)).strip() for k, r in (('CHROMATOGRAPHY_TYPE', 'chromatography_type'),)
            if not _rest_blank(rest_analysis.get(r))}
    if not fill:
        return mwtab
    out = dict(mwtab)
    out['CHROMATOGRAPHY'] = {**ch, **fill}
    print(f"Chromatography of {rest_analysis.get('analysis_id')} from the REST analysis record "
          f"(the mwTab states none): {fill}")
    return out


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
    if is_hilic and 'hilic' not in t and any(rx.search(col) for rx in _RP_FAMILIES.values()) \
            and not _HILIC_COLUMN_RE.search(col) and 'pfp' not in col:
        # 'Normal phase' with a reverse-phase column (ST003453 'BEH C8', 'BEH C18', 'HSS T3'):
        # the column decides
        is_hilic, is_rp = False, True
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


def resolve_analysis_conflicts_from_path(df, contexts=None):
    """A raw file listed under several analyses of a study (e.g. a positive and a
    negative mode analysis) gets conflicting analysis-level values, which are later
    blanked to 'missing value'. When the file path itself states the polarity and/or
    the chromatography (C18-NEG/..., ..._HILICpos...), keep only the analysis rows that
    agree with it - but only if that leaves rows that agree with each other on the
    analysis-level columns. Otherwise the group is left unchanged. contexts: {path:
    (polarity, chromatography)} of the listing (file_context, with the archive names);
    paths not in it get their path's context. A path without any polarity / chromatography
    token whose context comes from its archive name(s) alone (100.raw, only in an RPLC pos and
    an RPLC neg archive, named by RP, HILIC and direct-infusion analyses) also loses the rows of
    the analyses the archives contradict when the rows left still disagree (on the polarity):
    those analyses name a file of other archives; the merge keeps what the rows left agree on.
    Tokens of a path do not do this (a file name's 'C8' may be a well)."""
    analysis_cols = [c for c in ['IonizationSourceAndPolarity', 'ChromatographyAndPhase', 'MassSpectrometer'] if c in df.columns]
    dup = df['filename'].duplicated(keep=False)
    if not dup.any() or not analysis_cols:
        return df
    contexts = contexts or {}
    drop_idx = []
    n_partial = 0
    for filename, group in df[dup].groupby('filename'):
        if all(group[c].astype(str).nunique() == 1 for c in analysis_cols):
            continue
        polarity, chrom = contexts.get(filename) or _path_indicates(filename)
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
        elif filename in contexts and _path_tokens(filename) == (set(), []):
            drop_idx.extend(group.index[~keep])
            n_partial += 1
    if drop_idx:
        print(f"Resolved analysis conflicts from file paths: {len(drop_idx)} conflicting analysis rows dropped"
              + (f" ({n_partial} files by their archive names only, other analysis fields left to the merge)" if n_partial else ""))
    return df.drop(index=drop_idx)


def _value_contradicts_path(column, value, path_context):
    """True if an IonizationSourceAndPolarity / ChromatographyAndPhase value names another
    polarity / chromatography than a file path states ((polarity, chrom) of _path_indicates)."""
    v = str(value).lower()
    if v in ('', 'nan', 'none', 'missing value'):
        return False
    polarity, chrom = path_context
    if column == 'IonizationSourceAndPolarity' and polarity is not None:
        other = 'negative' if polarity == 'positive' else 'positive'
        return other in v and polarity not in v
    if column == 'ChromatographyAndPhase' and chrom is not None:
        return v.startswith('reverse phase') if chrom == 'hilic' else 'hilic' in v
    return False


_POL_WORD_RE = re.compile(r'positive|negative|pos|neg')
_POL_SWAP = {'positive': 'negative', 'negative': 'positive', 'pos': 'neg', 'neg': 'pos'}


def _path_without_extension(path):
    head, _, base = str(path).lower().rpartition('/')
    base = base.rsplit('.', 1)[0] if '.' in base else base
    return f'{head}/{base}' if head else base


def polarity_siblings(path):
    """The lower-cased path (without extension) with its polarity words swapped (pos<->neg,
    positive<->negative): each word alone, and all of them together."""
    p = _path_without_extension(path)
    words = list(_POL_WORD_RE.finditer(p))
    out = {p[:m.start()] + _POL_SWAP[m.group(0)] + p[m.end():] for m in words}
    if len(words) > 1:
        out.add(_POL_WORD_RE.sub(lambda m: _POL_SWAP[m.group(0)], p))
    return out


def clear_polarity_contradicted_by_sibling(df, raw_file_name_df, polarity_table=None):
    """IonizationSourceAndPolarity that the file's own path contradicts ('x_neg.mzML' with the
    positive mode of the study's only analysis, a sample listing its _pos and _neg files) is
    cleared to 'missing value' when the study lists the file's counterpart of the other
    polarity ('x_pos.mzML'): the polarity word in the path then tells two acquisitions
    apart. Without such a counterpart the word may belong to a sample or folder name
    ('Negative control1', 'triple negative breast cancer cells/') and the value stays; so
    does a value the file's scans (polarity_table) agree with. A path without a polarity
    word takes the polarity of its archive name(s) (file_context); the counterpart is then
    an archive of the other polarity ('x_NEG.zip' next to 'x_POS.zip')."""
    if df is None or len(df) == 0 or 'IonizationSourceAndPolarity' not in df.columns:
        return df
    listing = {_path_without_extension(p) for p in raw_file_name_df['filename'].astype(str)}
    contexts = listing_contexts(raw_file_name_df)
    archives_of = {}
    if 'archive' in raw_file_name_df.columns:
        archives_of = {p: _archive_names(a) for p, a in zip(raw_file_name_df['filename'].astype(str), raw_file_name_df['archive'])}
    study_archives = {_path_without_extension(a) for names in archives_of.values() for a in names}
    scans = {}
    if polarity_table is not None and len(polarity_table):
        _pt = polarity_table[polarity_table['USI'].isin(set(df['USI'].astype(str)))]
        for u, n_neg, n_pos in zip(_pt['USI'], pd.to_numeric(_pt['negative_polarity_count'], errors='coerce'),
                                   pd.to_numeric(_pt['positive_polarity_count'], errors='coerce')):
            if n_neg > 0 and not n_pos > 0:
                scans[u] = 'negative'
            elif n_pos > 0 and not n_neg > 0:
                scans[u] = 'positive'
    values = df['IonizationSourceAndPolarity'].copy()
    n = 0
    for i, (path, usi, v) in enumerate(zip(df['filename'].astype(str), df['USI'].astype(str), values.astype(str))):
        polarity, _ = _path_indicates(path)
        from_archive = polarity is None and not _path_tokens(path)[0]
        if from_archive:
            polarity = (contexts.get(path) or (None, None))[0]
        if polarity is None:
            continue
        other = 'negative' if polarity == 'positive' else 'positive'
        vl = v.lower()
        if other not in vl or polarity in vl or scans.get(usi) == other:
            continue
        if (polarity_siblings(path) & listing) if not from_archive else \
                any(polarity_siblings(a) & study_archives for a in archives_of.get(path, ())):
            values.iat[i] = 'missing value'
            n += 1
    if n:
        print(f"{n} polarity values cleared: the file's path names the other polarity and the study has its counterpart")
        df = df.copy()
        df['IonizationSourceAndPolarity'] = values.values
    return df


def _file_name_polarities(path):
    """The polarities the file name (last path part) of a listed path states (_path_tokens)."""
    return _path_tokens(str(path).replace('\\', '/').rsplit('/', 1)[-1])[0]


def _scan_polarities(df, polarity_table):
    """{USI: 'positive' / 'negative'} of the rows' files whose scans (polarity_table) are of one polarity"""
    scans = {}
    if polarity_table is not None and len(polarity_table) and 'USI' in df.columns:
        _pt = polarity_table[polarity_table['USI'].isin(set(df['USI'].astype(str)))]
        for u, n_neg, n_pos in zip(_pt['USI'], pd.to_numeric(_pt['negative_polarity_count'], errors='coerce'),
                                   pd.to_numeric(_pt['positive_polarity_count'], errors='coerce')):
            if n_neg > 0 and not n_pos > 0:
                scans[u] = 'negative'
            elif n_pos > 0 and not n_neg > 0:
                scans[u] = 'positive'
    return scans


def clear_polarity_contradicted_by_all_file_names(df, min_rows=3, polarity_table=None):
    """IonizationSourceAndPolarity of an analysis that the file name of every one of its files
    contradicts is cleared to 'missing value': when all (at least min_rows) rows of a study with
    the positive (negative) value have a file name stating only the other polarity, the mwTab's
    ION_MODE does not describe these files (ST004626: ION_MODE POSITIVE, title '... TOF MS
    (Negative)', every file 'lcqtof_t3_neg_*'; its positive twin study names its files
    'lcqtof_t3_pos_*'). Only the file name counts: a folder word can describe the samples of a
    whole study ('triple negative breast cancer cells/' in a positive-mode study). Values are
    only cleared, never set from the file name. Polarity-switching ('alternating') values and
    rows without a polarity are left alone. A file name word can also be a sample word ('Gram-negative_1',
    'HIV_neg_2', 'ER-negative_3'): the group is left alone when the scans of any of its files
    (polarity_table) agree with the value, and when its files are all blanks / QCs by their name."""
    if df is None or len(df) == 0 or 'IonizationSourceAndPolarity' not in df.columns or 'filename' not in df.columns:
        return df
    values = df['IonizationSourceAndPolarity'].astype(str).str.lower()
    names = df['filename'].astype(str).map(_file_name_polarities)
    scans = _scan_polarities(df, polarity_table)
    usis = df['USI'].astype(str) if 'USI' in df.columns else pd.Series('', index=df.index)
    clear = pd.Series(False, index=df.index)
    for pol, other in (('positive', 'negative'), ('negative', 'positive')):
        rows = values.str.contains(f'({pol})', regex=False) & ~values.str.contains(f'({other})', regex=False)
        if rows.sum() >= min_rows and all(n == {other} for n in names[rows]) \
                and not any(scans.get(u) == pol for u in usis[rows]) \
                and not all(blank_qc_type_from_filename(f) for f in df.loc[rows, 'filename'].astype(str)):
            clear |= rows
    if clear.any():
        print(f"{int(clear.sum())} polarity values cleared: the file name of every file with that polarity "
              f"names the other one")
        df = df.copy()
        df.loc[clear, 'IonizationSourceAndPolarity'] = 'missing value'
    return df


def clear_polarity_contradicted_by_folder_and_file_name(df, polarity_table=None):
    """IonizationSourceAndPolarity that both the folders and the file name of the file's path
    contradict (each states only the other polarity: 'Experiment 1 FA Positive/x_RW1_Pos.mzXML'
    with the negative mode of a study's only analysis, whose mwTab names the _Pos sample too) is
    cleared to 'missing value', also where the study lists no counterpart of the other polarity
    (clear_polarity_contradicted_by_sibling needs one: x_RW1_Neg is listed under another date).
    One word alone may be a sample or folder name ('Negative control1', 'triple negative breast
    cancer cells/'); two independent ones are not. A value the file's scans (polarity_table)
    agree with stays."""
    if df is None or len(df) == 0 or 'IonizationSourceAndPolarity' not in df.columns or 'filename' not in df.columns:
        return df
    scans = {}
    if polarity_table is not None and len(polarity_table) and 'USI' in df.columns:
        _pt = polarity_table[polarity_table['USI'].isin(set(df['USI'].astype(str)))]
        for u, n_neg, n_pos in zip(_pt['USI'], pd.to_numeric(_pt['negative_polarity_count'], errors='coerce'),
                                   pd.to_numeric(_pt['positive_polarity_count'], errors='coerce')):
            if n_neg > 0 and not n_pos > 0:
                scans[u] = 'negative'
            elif n_pos > 0 and not n_neg > 0:
                scans[u] = 'positive'
    usis = df['USI'].astype(str) if 'USI' in df.columns else pd.Series('', index=df.index)
    values = df['IonizationSourceAndPolarity'].copy()
    n = 0
    for i, (path, usi, v) in enumerate(zip(df['filename'].astype(str), usis, values.astype(str).str.lower())):
        p = path.replace('\\', '/')
        if '/' not in p:
            continue
        folder, name = p.rsplit('/', 1)
        f_pols, n_pols = _path_tokens(folder)[0], _path_tokens(name)[0]
        if len(f_pols) != 1 or f_pols != n_pols:
            continue
        stated = next(iter(f_pols))
        other = 'negative' if stated == 'positive' else 'positive'
        if f'({other})' in v and f'({stated})' not in v and scans.get(usi) != other:
            values.iat[i] = 'missing value'
            n += 1
    if n:
        print(f"{n} polarity values cleared: the file's folders and file name both name the other polarity")
        df = df.copy()
        df['IonizationSourceAndPolarity'] = values.values
    return df


def clear_chromatography_contradicted_by_path(df, raw_file_name_df):
    """ChromatographyAndPhase that the file's own path contradicts (a reverse-phase value on
    'C18/x.mzML', HILIC on 'HILIC_pos/x.raw' given a reverse phase) is cleared to 'missing
    value' in a study whose file paths name both separations (HILIC and C18 / RP folders):
    the path then tells the acquisitions apart. Elsewhere a token may belong to a sample
    name ('T3' time points, 'C8' wells) and the value stays. A path without a chromatography
    token takes that of its archive name(s) (file_context), in a study whose files' contexts
    (paths and archive names) name both separations."""
    if df is None or len(df) == 0 or 'ChromatographyAndPhase' not in df.columns:
        return df
    paths = raw_file_name_df['filename'].astype(str)
    path_seps = len({_path_indicates(p)[1] for p in paths} - {None}) >= 2
    contexts = listing_contexts(raw_file_name_df)
    archive_seps = len({c[1] for c in contexts.values()} - {None}) >= 2
    if not path_seps and not archive_seps:
        return df
    values = df['ChromatographyAndPhase'].copy()
    n = 0
    for i, (path, v) in enumerate(zip(df['filename'].astype(str), values.astype(str))):
        if _path_tokens(path)[1]:
            chrom = _path_indicates(path)[1] if path_seps else None
        else:
            chrom = (contexts.get(path) or (None, None))[1] if archive_seps else None
        if chrom is not None and _value_contradicts_path('ChromatographyAndPhase', v, (None, chrom)):
            values.iat[i] = 'missing value'
            n += 1
    if n:
        print(f"{n} chromatography values cleared: the file's path names the other separation")
        df = df.copy()
        df['ChromatographyAndPhase'] = values.values
    return df


def propagate_to_same_name_files(redu_df, raw_file_name_df, polarity_table=None):
    """Studies often hold one raw file per sample under several method folders
    (Method_1/Sample_1.raw, Method_2/Sample_1.raw, ...), but the sample-to-file match
    links each sample to one of them. Give the other same-name files (same basename and
    extension) the sample's metadata, when the basename belongs to exactly one sample.
    Analysis-level fields are copied only if the study has a single analysis setup;
    otherwise they come from the file path tokens, or stay missing. A polarity /
    chromatography that the copy's own folders (Negative/..., rp_neg/...) or its scans
    (polarity_table) contradict is never copied: such a file is another acquisition of
    the sample. A copy's context takes a dimension its folders leave open from the
    copy's archive name(s) (file_context); its file name, the same as the matched file's,
    is not read (same_name_copy_context)."""
    if redu_df is None or len(redu_df) == 0 or raw_file_name_df is None or len(raw_file_name_df) == 0:
        return redu_df
    out = redu_df.copy()
    out['_base'] = out['filename'].astype(str).map(lambda x: os.path.basename(x).lower())
    base_counts = out['_base'].value_counts()
    unique_rows = out[out['_base'].map(base_counts) == 1].set_index('_base')
    # a basename on several rows of one sample (one file per acquisition context, e.g.
    # Positive/x.raw and Negative/x.raw, each matched by its own analysis): the values its
    # rows agree on. Rows that differ in a sample field are different samples sharing a
    # file name; their basename is left alone.
    multi = out[out['_base'].map(base_counts) > 1]
    if len(multi):
        per_file = {'filename', 'USI', '_base', 'IonizationSourceAndPolarity', 'ChromatographyAndPhase', 'MassSpectrometer',
                    'YearOfAnalysis', 'SampleCollectionMethod', 'SampleExtractionMethod', 'InternalStandardsUsed',
                    'negative_polarity_count', 'positive_polarity_count'}
        sample_cols = [c for c in out.columns if c not in per_file]
        merged = []
        for _b, g in multi.groupby('_base', sort=True):
            if any(g[c].astype(str).nunique() > 1 for c in sample_cols):
                continue
            row = g.iloc[0].copy()
            for c in per_file - {'filename', 'USI', '_base'}:
                if c in g.columns and g[c].astype(str).nunique() > 1:
                    row[c] = np.nan if c.endswith('_polarity_count') else 'missing value'
            merged.append(row)
        if merged:
            unique_rows = pd.concat([unique_rows, pd.DataFrame(merged).set_index('_base')])
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
    contexts = listing_contexts(raw_file_name_df)
    if not single_setup:
        # a field with one value across the study's files (e.g. the instrument) is still
        # known for the copies
        for c in analysis_cols:
            if out[c].astype(str).nunique() != 1:
                copies[c] = 'missing value'
        for i, path in enumerate(copies['filename']):
            polarity, chrom = contexts.get(path) or _path_indicates(path)
            match = out
            if polarity is not None:
                match = match[match['IonizationSourceAndPolarity'].astype(str).str.lower().str.contains(polarity)]
            if chrom is not None and 'ChromatographyAndPhase' in out.columns:
                match = match[match['ChromatographyAndPhase'].astype(str).str.lower().str.contains(chrom)]
            if (polarity is not None or chrom is not None) and len(match) and all(match[c].astype(str).nunique() == 1 for c in analysis_cols):
                for c in analysis_cols:
                    copies.at[i, c] = match[c].iloc[0]
    # never a polarity / chromatography the copy's own path or scans contradict
    n_cleared = 0
    scans = {}
    if polarity_table is not None and len(polarity_table) and 'IonizationSourceAndPolarity' in copies.columns:
        _pt = polarity_table[polarity_table['USI'].isin(set(copies['USI']))]
        for u, n_neg, n_pos in zip(_pt['USI'], pd.to_numeric(_pt['negative_polarity_count'], errors='coerce'),
                                   pd.to_numeric(_pt['positive_polarity_count'], errors='coerce')):
            if n_neg > 0 and not n_pos > 0:
                scans[u] = 'negative'
            elif n_pos > 0 and not n_neg > 0:
                scans[u] = 'positive'
    # a copy is told apart from the matched file by its folders / archive only: both have the
    # same file name (same_name_copy_context)
    copy_contexts = listing_copy_contexts(raw_file_name_df)
    for i, path in enumerate(copies['filename']):
        context = copy_contexts.get(path) or same_name_copy_context(path)
        for c in ('IonizationSourceAndPolarity', 'ChromatographyAndPhase'):
            if c not in copies.columns:
                continue
            ctx = context
            if c == 'IonizationSourceAndPolarity' and copies.at[i, 'USI'] in scans:
                ctx = (scans[copies.at[i, 'USI']], context[1])
            if _value_contradicts_path(c, copies.at[i, c], ctx):
                copies.at[i, c] = 'missing value'
                n_cleared += 1
    print(f"Copied sample metadata to {len(copies)} same-name raw files in other folders"
          + (f" ({n_cleared} polarity/chromatography values their path or scans contradict not copied)" if n_cleared else ""))
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
        # file contexts of the listing (path, else archive names: file_context)
        df = resolve_analysis_conflicts_from_path(df.reset_index(drop=True), contexts=kwargs.get('file_contexts'))
    except Exception as e:
        print(f"Could not resolve analysis conflicts from file paths: {type(e).__name__}: {e}")


    # groupby drops rows without a filename (dropna=True); so does the fast path below
    df = df[df['filename'].notna()]
    dup = df['filename'].duplicated(keep=False)
    if len(df) == 0 or not (~dup).any():
        # Group by 'USI' and apply the processing function to each group
        df = df.groupby('filename').apply(process_group)
    else:
        # Files listed once (most of them) need no per-group work: for a one-row group
        # process_group only turns None into 'missing value' (and 'None' in SampleType /
        # SampleTypeSub1). Only files under several analyses go through process_group.
        # Rows come out in the order groupby.apply gives (sorted by filename), with its index.
        group_no = df.groupby('filename').ngroup()
        single = df[~dup].copy()
        for col in single.columns:
            if single[col].dtype != object:
                continue
            if col in ('SampleType', 'SampleTypeSub1'):
                is_none = np.fromiter((v is None or (isinstance(v, str) and v == 'None') for v in single[col].values), bool, len(single))
            else:
                is_none = np.fromiter((v is None for v in single[col].values), bool, len(single))
            if is_none.any():
                single.loc[is_none, col] = 'missing value'
        parts = [single]
        if dup.any():
            multi = df[dup].groupby('filename').apply(process_group)
            multi.index = multi.index.droplevel(0)
            parts.append(multi)
        df = pd.concat(parts)
        df = df.iloc[np.argsort(group_no.loc[df.index].values, kind='stable')]
        df.index = pd.MultiIndex.from_arrays([df['filename'].values, df.index], names=['filename', None])
        # groupby.apply gives a numeric column when no group holds a string (e.g. ints and
        # floats -> float64): same here
        for col in df.columns:
            if df[col].dtype == object and not any(isinstance(v, str) for v in df[col].values):
                df[col] = df[col].infer_objects()

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


# --- mwTab responses that are not valid JSON -------------------------------------------
# About 50 MS analyses get an mwTab JSON that json.loads rejects: control characters in
# strings; '=' or ':' inside a sample's factor / additional-data values, which the JSON
# writer turns into "key":"a":"b"; and two documents in one response (a copy of the same
# one, an empty stub plus the real one, or one per study for an analysis shared by two
# studies). The document of the study is taken; a broken SUBJECT_SAMPLE_FACTORS list is
# replaced by the one of the mwTab text format (/mwtab/txt). Everything else still comes
# from the JSON.
_MWTAB_FIELD_RE = re.compile(r'^[A-Z]{2}:([A-Za-z0-9_]+)[ ]*\t(.*)$')
_MWTAB_SECTIONS = {'PROJECT', 'STUDY', 'SUBJECT', 'COLLECTION', 'TREATMENT', 'SAMPLEPREP', 'CHROMATOGRAPHY',
                   'ANALYSIS', 'MS', 'NMR'}
_MWTAB_JSON_DOC_SPLIT_RE = re.compile(r'\n\}\n\n(?=\{\n"METABOLOMICS WORKBENCH")')
_MWTAB_SSF_END_RE = re.compile(r'\n\],\n"[A-Z][A-Z_ ]*":')


def _mwtab_value(v):
    # the JSON format has no double quotes and no surrounding spaces
    return v.replace('"', '').strip(' ')


def _mwtab_pairs(text, sep, kv):
    """[(key, value)] of 'k=v; k=v' (sep '; ', kv '='); a piece without kv continues the
    previous value ('-&gt; 64 ul' in 'Comment=a-&gt; 64 ul'). In additional sample data the
    text format sometimes has two spaces where '; ' separated two items
    ('RAW_FILE_NAME=a.raw  RAW_FILE_NAME=b.raw'); the JSON of such analyses is never valid
    (its writer cannot have '=' in a value), so this only applies to the fallback."""
    pieces = text.split(sep)
    if kv == '=':
        split = []
        for piece in pieces:
            k, found, v = piece.partition(kv)
            # only inside the value: keys may have two spaces ('Age at  Consent=45')
            more = re.split(r'\s{2,}(?=[^=;\s][^=;]*=)', v) if found else [v]
            split.append(k + found + more[0])
            split.extend(more[1:])
        pieces = split
    out = []
    for piece in pieces:
        k, found, v = piece.partition(kv)
        if found:
            out.append([k, v])
        elif out:
            out[-1][1] += sep + piece
    return [(_mwtab_value(k), _mwtab_value(v)) for k, v in out]


def mwtab_text_to_dict(text):
    """One mwTab text-format document in the layout of the JSON format: the
    'METABOLOMICS WORKBENCH' header, the key-value sections and SUBJECT_SAMPLE_FACTORS
    ([{'Subject ID', 'Sample ID', 'Factors', 'Additional sample data'}])."""
    sections, order, header, factors = {}, [], {}, []
    section, last_key = None, None
    for line in text.split('\n'):
        line = line.rstrip('\r')
        if line.startswith('#'):
            last_key = None
            if line.startswith('#METABOLOMICS WORKBENCH'):
                section = 'METABOLOMICS WORKBENCH'
                for tok in line.split():
                    k, sep, v = tok.partition(':')
                    if sep and k in ('STUDY_ID', 'ANALYSIS_ID') and k not in header:
                        header[k] = v
                continue
            name = line[1:].split(':')[0].split('\t')[0].strip()
            section = name if name in _MWTAB_SECTIONS or name == 'SUBJECT_SAMPLE_FACTORS' else None
            if section in _MWTAB_SECTIONS and section not in sections:
                sections[section] = {}
                order.append(section)
            continue
        if section == 'METABOLOMICS WORKBENCH':
            for key in ('VERSION', 'CREATED_ON'):
                if line.startswith(key):
                    header[key] = line[len(key):].replace('\t', ' ')
            continue
        if section == 'SUBJECT_SAMPLE_FACTORS':
            parts = line.split('\t')
            if not line.startswith('SUBJECT_SAMPLE_FACTORS') or len(parts) < 4:
                continue
            item = [('Subject ID', _mwtab_value(parts[1])), ('Sample ID', _mwtab_value(parts[2])),
                    ('Factors', handle_duplicates(_mwtab_pairs(parts[3], ' | ', ':')))]
            add = _mwtab_pairs('\t'.join(parts[4:]), '; ', '=') if len(parts) > 4 else []
            if add:
                item.append(('Additional sample data', handle_duplicates(add)))
            factors.append(handle_duplicates(item))
            continue
        m = _MWTAB_FIELD_RE.match(line) if section in _MWTAB_SECTIONS else None
        if not m:
            continue
        key, value = m.group(1), m.group(2).replace('\t', ' ')
        d = sections[section]
        # continuation lines of a field are joined; a field given again later replaces it
        d[key] = d[key] + ' ' + value if key == last_key else value
        last_key = key
    out = {'METABOLOMICS WORKBENCH': {k: _mwtab_value(header[k]) for k in ('STUDY_ID', 'ANALYSIS_ID', 'VERSION', 'CREATED_ON')
                                      if k in header and _mwtab_value(header[k]) != ''}}
    for sec in order:
        out[sec] = {k: _mwtab_value(v) for k, v in sections[sec].items() if _mwtab_value(v) != ''}
        if sec == 'SUBJECT':
            out['SUBJECT_SAMPLE_FACTORS'] = factors
    out.setdefault('SUBJECT_SAMPLE_FACTORS', factors)
    return out


def _mwtab_text_sample_factors(text, study_id, analysis_id):
    """SUBJECT_SAMPLE_FACTORS of the study's document of an mwTab text response, or None."""
    docs = [d for d in re.split(r'(?m)^(?=#METABOLOMICS WORKBENCH)', text) if d.startswith('#METABOLOMICS WORKBENCH')]
    found = []
    for doc in docs:
        head = doc.split('\n', 1)[0]
        # with several documents, the study's is the one whose header names it (a single
        # document is taken as it is, as the JSON path does)
        if len(docs) > 1 and not (study_id in re.findall(r'STUDY_ID:(\S+)', head)
                                  and analysis_id in re.findall(r'ANALYSIS_ID:(\S+)', head)):
            continue
        d = mwtab_text_to_dict(doc)
        if d['SUBJECT_SAMPLE_FACTORS']:
            found.append(d['SUBJECT_SAMPLE_FACTORS'])
    if found and all(f == found[0] for f in found):
        return found[0]
    return None


# analyses left out because the Workbench gave no usable mwTab, with the reason
DROPPED_ANALYSES = {}


def _drop_analysis(study_id, analysis_id, reason):
    print(f"DROPPED analysis {analysis_id} of {study_id}: {reason}")
    DROPPED_ANALYSES[f'{study_id}/{analysis_id}'] = reason.split(' (')[0]


def mwtab_from_invalid_json(text, study_id, analysis_id, error, get_text=None):
    """The mwTab of an analysis whose JSON response json.loads rejected, or None (with the
    reason in the log). get_text(url) returns the text of a Workbench response."""
    if not text.strip():
        _drop_analysis(study_id, analysis_id, 'the Workbench answered with an empty mwTab')
        return None
    docs = _MWTAB_JSON_DOC_SPLIT_RE.split(text)
    docs = [d + '\n}' for d in docs[:-1]] + docs[-1:]
    parsed, problems = [], []
    for doc in docs:
        try:
            parsed.append((json.loads(doc, strict=False, object_pairs_hook=handle_duplicates), False))
            continue
        except ValueError as e:
            doc_error = e
        start = doc.find('\n"SUBJECT_SAMPLE_FACTORS":[')
        end = _MWTAB_SSF_END_RE.search(doc, start) if start >= 0 else None
        if end is None:
            problems.append(f"{type(doc_error).__name__}: {doc_error}")
            continue
        try:
            rest = doc[:start] + '\n"SUBJECT_SAMPLE_FACTORS":[]' + doc[end.start() + 2:]
            parsed.append((json.loads(rest, strict=False, object_pairs_hook=handle_duplicates), True))
        except ValueError as e:
            problems.append(f"{type(e).__name__}: {e}")
    def _study_of(d):
        head = d.get('METABOLOMICS WORKBENCH') if isinstance(d, dict) else None
        return head.get('STUDY_ID') if isinstance(head, dict) else None
    # with several documents, the study's are those whose header names it; a single document
    # is taken as it is, as the JSON path does (some headers name another study or none)
    mine = [(d, broken) for d, broken in parsed if isinstance(d, dict) and (len(docs) == 1 or _study_of(d) == study_id)]
    with_samples = [(d, broken) for d, broken in mine if broken or d.get('SUBJECT_SAMPLE_FACTORS')]
    if with_samples:
        mine = with_samples
    distinct = []
    for d, broken in mine:
        if not any(d == e and broken == b for e, b in distinct):
            distinct.append((d, broken))
    if len(distinct) != 1:
        _drop_analysis(study_id, analysis_id, f"the mwTab JSON is not valid JSON and has {len(distinct)} usable "
                                              f"documents of the study ({error}; of {len(docs)} documents: {problems[:3]})")
        return None
    mwtab, broken = distinct[0]
    how = 'control characters in strings accepted' if len(docs) == 1 else f'the document of {study_id} of {len(docs)} documents'
    if broken:
        factors = None
        try:
            if get_text is not None:
                factors = _mwtab_text_sample_factors(
                    get_text(f"https://www.metabolomicsworkbench.org/rest/study/analysis_id/{analysis_id}/mwtab/txt"),
                    study_id, analysis_id)
        except Exception as e:
            print(f"Could not get the mwTab text of {analysis_id}: {type(e).__name__}: {e}")
        if not factors:
            _drop_analysis(study_id, analysis_id, f"the samples of the mwTab JSON are not valid JSON and the mwTab "
                                                  f"text gives none ({error})")
            return None
        mwtab['SUBJECT_SAMPLE_FACTORS'] = factors
        how = ('' if len(docs) == 1 else how + '; ') + f'its {len(factors)} samples (SUBJECT_SAMPLE_FACTORS) from the mwTab text'
    print(f"mwTab JSON of {analysis_id} is not valid JSON ({error}): read with {how}")
    return mwtab


# Some archive listings carry raw 7z listing columns in front of the path,
# e.g. '2014-10-05 10:45:02 ..H.A           80               MX183900 .../x.d'
_7Z_LISTING_PREFIX = re.compile(r'^\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2} \S{5} +\d+ +(?:\d+ +)?')


def clean_listing_path(path):
    return _7Z_LISTING_PREFIX.sub('', path)


def collapse_d_folder(path):
    """Agilent/Bruker .d and Waters .raw raw data are folders, listed only as the files
    inside them; represent them by the folder (the same function as
    getAllWorkbench_file_paths.py, whose paths the rows are joined with)."""
    return collapse_vendor_folder(path)


# studies left out because the Workbench gave no file list, with the reason
DROPPED_STUDIES = {}

# ---- files the converter withholds on purpose (side table for saveOlderData --withheld_files) -----------------
# Files whose sample the name matching does not pair (substring_file_matches: named after another or several
# samples, several names fit equally, no whole-token match, a blank / QC file, only a placeholder name; same-name
# copies of the files of several samples: contexts_told_apart_by_listed_folders) get no row; their rows of earlier
# runs (paired by the earlier rules, also copies of them in other folders) must not come back through saveOlderData.
# Per study, the withheld files and the study's other files of their names without a row of the same file (USI
# without extension) this run, for studies read without an error. Columns USI, dataset, reason.
WITHHELD_COLUMNS = ['USI', 'dataset', 'reason']
WITHHELD_FILES = []
_STUDY_WITHHELD = {}
_WITHHELD_REASONS = {'named after another sample', 'named after several samples', 'several names fit equally',
                     'no whole-token match', 'blank / QC file, study sample name', 'only a placeholder name'}
_STUDY_LISTING = []
_MS_EXT_RE = re.compile(r'\.(?:mzml|mzxml|mzdata|mgf|raw|cdf|wiff2?|scan|d)$', re.I)


def note_withheld(usi, reason='sample_pairing_withheld'):
    if isinstance(usi, str) and usi.strip():
        _STUDY_WITHHELD.setdefault(usi.strip(), reason)


def _file_stem(usi):
    """file name without folders and extensions, lower case ('mzspec:ST1:a/X.mzML' -> 'x')"""
    return re.split(r'[\\/:]', str(usi))[-1].split('.')[0].lower()


def flush_study_withheld(study_id, result):
    """The study's withheld files, and the other listed files of their names (copies in other folders or formats,
    which the matching leaves out as copies: _STUDY_LISTING), without a row of the same file in result (USI without
    extension), to WITHHELD_FILES"""
    out = set()
    if result is not None and len(result) and 'USI' in result.columns:
        out = {_MS_EXT_RE.sub('', str(u)).lower() for u in result['USI'].dropna()}
    withheld = dict(_STUDY_WITHHELD)
    stems = {_file_stem(u): r for u, r in withheld.items()}
    for u in _STUDY_LISTING:
        if isinstance(u, str) and u.strip() and _file_stem(u) in stems:
            withheld.setdefault(u.strip(), stems[_file_stem(u)])
    for usi, reason in sorted(withheld.items()):
        if _MS_EXT_RE.sub('', usi).lower() not in out:
            WITHHELD_FILES.append({'USI': usi, 'dataset': str(study_id), 'reason': reason})
    _STUDY_WITHHELD.clear()
    del _STUDY_LISTING[:]


def write_withheld_files(path, rows):
    """The side table of withheld files (header only when there are none), sorted by USI."""
    table = pd.DataFrame(rows, columns=WITHHELD_COLUMNS).drop_duplicates(subset=['USI'])
    table.sort_values('USI', kind='mergesort').to_csv(path, sep='\t', index=False)
    print(f'Withheld files (sample pairing withheld): {len(table)} of {table["dataset"].nunique()} studies -> {path}')


def _get_metabolomicsworkbench_files(dataset_accession):
    try:
        mw_file_list, _ = fetch_archive_listing(dataset_accession, _mwb_get)

        acceptable_extensions = [".mzml", ".mzxml", ".cdf", ".raw", ".wiff", ".d"]

        mw_file_list = [file_obj for file_obj in mw_file_list if not is_apple_double(file_obj["FILENAME"])]
        workbench_df = pd.DataFrame(mw_file_list, columns=["FILENAME", "FILESIZE", "URL", "ARCHIVE"])
        workbench_df["filename"] = workbench_df["FILENAME"].apply(clean_listing_path).apply(collapse_d_folder)
        workbench_df = workbench_df[workbench_df["filename"].str.lower().str.endswith(tuple(acceptable_extensions))]
        workbench_df["FILESIZE"] = pd.to_numeric(workbench_df["FILESIZE"], errors='coerce').fillna(0)
        # the archive each file lies in (ST004470_HILIC_Rawfiles.zip / ST004470_RP_Rawfiles.zip):
        # where the path states no polarity / chromatography, the archive name may
        # (file_context); a row without ARCHIVE has it as the URL's A= parameter
        workbench_df["ARCHIVE"] = [a if isinstance(a, str) and a else listing_query(u).get('A', [''])[0]
                                   for a, u in zip(workbench_df["ARCHIVE"], workbench_df["URL"].astype(str))]

        # Generate USI column
        workbench_df["USI"] = workbench_df["URL"].apply(
            lambda url: "mzspec:{}:{}".format(
                dataset_accession,
                collapse_d_folder(clean_listing_path(listing_query(url).get('F', [''])[0]))
            )
        )

        # one row per file / .d folder (summing the sizes of the files inside a folder), with
        # the names of the archives holding the path (several for one path in a pos and a neg
        # archive: one USI for both), sorted and joined by ARCHIVE_SEP
        archives = workbench_df.drop_duplicates(["filename", "ARCHIVE"]).groupby("filename", sort=False)["ARCHIVE"].agg(
            lambda s: ARCHIVE_SEP.join(_archive_names(s.tolist())))
        workbench_df = workbench_df.groupby("filename", as_index=False, sort=False).agg(
            FILESIZE=("FILESIZE", "sum"), USI=("USI", "first"))
        workbench_df["size_mb"] = (workbench_df["FILESIZE"] / 1024 / 1024).astype(int)
        workbench_df["archive"] = workbench_df["filename"].map(archives).fillna('')

        workbench_df = workbench_df[["filename", "size_mb", "USI", "archive"]]
    except ArchiveListingError as e:
        print(f"DROPPED {dataset_accession}: no archive listing ({e})")
        DROPPED_STUDIES[dataset_accession] = 'no archive listing'
        workbench_df = pd.DataFrame()
    except Exception as e:
        print(f"Failed to list files for {dataset_accession}: {type(e).__name__}: {e}")
        DROPPED_STUDIES[dataset_accession] = f'file list unreadable ({type(e).__name__})'
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
        if 'sea' in x or 'ocean' in x or 'coast' in x or 'marine' in x:
            return ['environmental', 'water_seawater']
        # a water sample of no stated kind: sea water when the study's title / summary
        # is about the sea and names no fresh water
        if _MARINE_TEXT_RE.search(_STUDY_TEXT_REST or '') and not _FRESHWATER_TEXT_RE.search(_STUDY_TEXT_REST or ''):
            return ['environmental', 'water_seawater']
        # surface water when the sample's text or the study's title / summary names lakes, rivers,
        # ponds ... and no other kind of water (water_kind.py)
        kind = water_kind_of_unstated_water(x, _STUDY_TEXT_REST)
        if kind:
            return ['environmental', kind]
        else:
            return [None, None]
    else:
        return [None, None]


_MARINE_TEXT_RE = re.compile(r'(?i)\b(?:marine|oceans?|oceanic|seas?|seawater|sea water|coastal|atlantic|antarctic|arctic|'
                             r'(?:north|south|eastern|western|central|subtropical|tropical|equatorial) pacific|pacific ocean|'
                             r'research cruise|cruise)\b')
_FRESHWATER_TEXT_RE = re.compile(r'(?i)\b(?:lakes?|rivers?|freshwater|fresh water|ponds?|reservoirs?|streams?|creeks?|'
                                 r'groundwater|wastewater|estuar\w*|brackish)\b')


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


# Workbench sample sources of QC samples ('_QC_', '_POOLED_', lower case here): pooled QCs
QC_SAMPLE_SOURCES = ('_qc_', '_pooled_')
# reference materials (reference_materials.py): the sample-level items of the files of the study being
# converted (raw file path -> evidence cell) and the catalog numbers its mwTab sections name; set per
# analysis (create_dataframe_from_SUBJECT_SAMPLE_FACTORS, MWB_to_REDU_wrapper), reset per study
_REFMAT_EVIDENCE = {}
_REFMAT_STUDY_NUMBERS = set()
# mwTab sections that hold no study text (sample factors, data)
_REFMAT_NOT_TEXT_SECTIONS = re.compile(r'(?i)^(?:SUBJECT_SAMPLE_FACTORS|MS_METABOLITE_DATA|NMR_.*DATA|METABOLITES|'
                                       r'.*_DATA|.*METABOLITE.*)$')
# raw-file keys (their values are file names: the file itself is read by reference_materials.py)
_REFMAT_RAW_KEY = re.compile(r'(?i)raw[ _]?(?:file|data)|data ?file|filename|^rpm_|^zhp_')


def mwb_mixture_design_samples(df):
    """rows of samples (column filename: the Sample ID here) whose factor / additional-data items are amount
    factors of a reference material ('SRM 1950: 500', 'SRM 8231: 0'): a mixture design"""
    try:
        if not {'filename', 'Key', 'Value'} <= set(df.columns):
            return pd.Series(False, index=df.index)
        item = df['Key'].fillna('').astype(str) + ': ' + df['Value'].fillna('').astype(str)
        hit = item.map(reference_materials.is_amount_factor)
        samples = set(df.loc[hit, 'filename'])
        return df['filename'].isin(samples)
    except Exception as e:
        print(f"WARNING: mixture designs not checked: {type(e).__name__}: {e}")
        return pd.Series(False, index=df.index)


def mwb_qc_source_not_pool(df):
    """rows of samples (column filename: the Sample ID here) of the QC sample sources '_QC_' / '_POOLED_' that are no
    pool of the study's samples: a value of the sample says blank (qc_sample_type: 'Sample Type: Blank',
    'Class: Water_blank', 'Group: QC-Blank'), its Sample ID or raw file name says blank, standard, calibration or
    solvent (blank_qc_names.not_a_pool_by_name with NOT_POOL_QC_NAME_WORDS: 'ISTD_QAQC', 'p1Std'; not 'QC mix',
    'Tech-Mix', 'Pool_QC_Media': pools of the study's samples), or a value of a sample-type key
    names a standard, calibration, solvent or blank (blank_qc_names.not_a_pool_by_value: 'Sample source: Standard',
    'Solvent Mixture', 'QC(Internal Standard)'). These stay blank_QC (no organism)."""
    try:
        if not {'filename', 'Key', 'Value'} <= set(df.columns):
            return pd.Series(False, index=df.index)
        sid = df['filename'].fillna('').astype(str)
        blank = pd.Series([qc_sample_type(k, v) in ('blank_analysis', 'blank_extraction', 'blank_culturemedia')
                           for k, v in zip(df['Key'], df['Value'])], index=df.index)
        type_key = df['Key'].map(lambda k: _norm_value(k) in _QC_TYPE_KEYS)
        by_value = type_key & df['Value'].map(lambda v: not_a_pool_by_value(v) if isinstance(v, str) else False)
        samples = set(sid[blank | by_value])
        words = NOT_POOL_QC_NAME_WORDS
        samples |= {s_ for s_ in set(sid) if s_ and not_a_pool_by_name(s_, words)}
        if 'filename_raw_path' in df.columns:
            raw = df['filename_raw_path'].fillna('').astype(str)
            samples |= set(sid[raw.map(lambda r: bool(r) and not_a_pool_by_name(r, words))])
        return sid.isin(samples)
    except Exception as e:
        print(f"WARNING: QC samples not checked for blanks / standards: {type(e).__name__}: {e}")
        return pd.Series(False, index=df.index)


def mwb_reset_refmat():
    _REFMAT_EVIDENCE.clear()
    _REFMAT_STUDY_NUMBERS.clear()


def mwb_note_refmat_study_text(mwtab):
    """catalog numbers of reference materials named in an analysis' mwTab sections (all but the
    sample factors and data), added to the study's (_REFMAT_STUDY_NUMBERS)"""
    try:
        parts = []
        if isinstance(mwtab, dict):
            for sec, body in mwtab.items():
                if _REFMAT_NOT_TEXT_SECTIONS.match(str(sec)):
                    continue
                if isinstance(body, dict):
                    parts += [str(v) for v in body.values() if isinstance(v, (str, int, float))]
                elif isinstance(body, str):
                    parts.append(body)
        _REFMAT_STUDY_NUMBERS.update(reference_materials.study_text_materials(*parts))
    except Exception as e:
        print(f"WARNING: mwTab text not read for reference materials: {type(e).__name__}: {e}")


def mwb_note_refmat_evidence(df):
    """Per sample (column filename: the Sample ID here) its factor and additional-data items
    ('key: value', not the raw-file keys; those with a candidate word or only a matrix) and its Sample ID, as one
    reference_materials.EVIDENCE_COLUMN
    cell, stored for each of its raw files (filename_raw_path) in _REFMAT_EVIDENCE; a file of samples
    with different items gets none."""
    try:
        if not {'filename', 'Key', 'Value'} <= set(df.columns):
            return
        raw = df['filename_raw_path'] if 'filename_raw_path' in df.columns else df['filename']
        key = df['Key'].fillna('').astype(str)
        item = key + ': ' + df['Value'].fillna('').astype(str)
        # items with a candidate word, items that only name a matrix ('Sample type: plasma': the sample's own matrix;
        # not the REST record's 'Sample Source', the Workbench's sample source of the study, which a QC sample of
        # another matrix shares: ST001335 'blood (serum)' on NIST plasma) and species items saying human ('Species: Human')
        use = ~key.str.contains(_REFMAT_RAW_KEY, regex=True) & (
            item.str.contains(reference_materials.CANDIDATE_RE, regex=True)
            | (item.map(reference_materials.is_matrix_value) & (key != 'Sample Source'))
            | pd.Series([reference_materials.is_human_species_item(k, v) for k, v in zip(key, item)], index=df.index))
        sid = df['filename'].fillna('').astype(str)
        per_sample = {}
        for s_, k, it in zip(sid[use], key[use], item[use]):
            per_sample.setdefault(s_, []).append((k, it))
        for s_ in set(sid):
            if reference_materials.CANDIDATE_RE.search(s_):
                per_sample.setdefault(s_, []).append(('sample_id', s_))
        by_file = {}
        for s_, r in zip(sid, raw.fillna('').astype(str)):
            if not r:
                continue
            by_file.setdefault(r, set()).add(reference_materials.encode_evidence(per_sample.get(s_, [])))
        for r, evs in by_file.items():
            prev = _REFMAT_EVIDENCE.get(r)
            evs = evs | ({prev} if prev is not None else set())
            _REFMAT_EVIDENCE[r] = next(iter(evs)) if len(evs) == 1 else ''
    except Exception as e:
        print(f"WARNING: no reference-material evidence: {type(e).__name__}: {e}")


def mwb_refmat_evidence_column(filenames):
    """reference_materials.EVIDENCE_COLUMN cells for the study's rows: by raw file path, else by a base
    name that only one recorded file has (the copies propagate_to_same_name_files adds)"""
    by_base = {}
    for r, e in _REFMAT_EVIDENCE.items():
        by_base.setdefault(os.path.basename(r), set()).add(e)
    def one(f):
        f = str(f)
        if f in _REFMAT_EVIDENCE:
            return _REFMAT_EVIDENCE[f]
        evs = by_base.get(os.path.basename(f), set())
        return next(iter(evs)) if len(evs) == 1 else ''
    return [one(f) for f in filenames]


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
# DOID is_a links (DOIDOntologyIndex, ParentDOIDOntologyIndex); set from --path_to_doid_parents_csv
DOID_PARENTS_CSV = 'none'


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


# ---------------------------------------------------------------------------------
# Study-level disease (Workbench REST 'disease', e.g. "Cancer", "Diabetes") given to
# the samples of a study's case groups, and HealthStatus 'healthy' to samples of
# explicitly healthy control groups. Health is decided per individual (Subject ID):
# one case sample makes all samples of that subject cases. Only fills gaps.
# ---------------------------------------------------------------------------------
_DOID_CANCER = 'DOID:162'
_DOID_INFECTIOUS = 'DOID:0050117'
# prepared study-level disease of the study being converted (MWB_to_REDU_study_wrapper)
_STUDY_DISEASE = None
_DOID_INDEX_CACHE = {}
_DISEASE_RE_CACHE = {}

# adjectives of single-word disease names ('diabetic' for diabetes)
_DISEASE_ADJECTIVES = {
    'diabetes': ['diabetic', 'diabetics'], 'diabetes mellitus': ['diabetic', 'diabetics'], 'obesity': ['obese'],
    'asthma': ['asthmatic', 'asthmatics'], 'anemia': ['anemic', 'anaemic'], 'autism': ['autistic'],
    'epilepsy': ['epileptic'], 'schizophrenia': ['schizophrenic'], 'tuberculosis': ['tuberculous'],
    'psoriasis': ['psoriatic'], 'malaria': ['malarial'], 'hypertension': ['hypertensive'],
    'arthritis': ['arthritic'], 'cirrhosis': ['cirrhotic'], 'atherosclerosis': ['atherosclerotic'],
    'ischemia': ['ischemic', 'ischaemic'], 'colitis': ['colitic'], 'neuropathy': ['neuropathic'],
    'cardiomyopathy': ['cardiomyopathic'], 'glaucoma': ['glaucomatous'], 'cancer': ['cancerous'],
    'dementia': ['demented'], 'hyperglycemia': ['hyperglycemic'], 'hypoglycemia': ['hypoglycemic'],
    'steatosis': ['steatotic'], 'preeclampsia': ['preeclamptic'], 'pre eclampsia': ['preeclamptic'],
}
# DOID names that are ordinary words in factor values
_DISEASE_NAME_SKIP = {'all', 'aids', 'cold', 'pain', 'fever', 'shock', 'stress', 'injury', 'disease', 'diseases',
                      'syndrome', 'disorder', 'infection', 'tumor', 'tumour', 'neoplasm', 'lesion', 'mass', 'growth',
                      'obese', 'lean', 'aging', 'ageing', 'sick', 'case', 'cases', 'control', 'normal', 'positive',
                      'negative', 'mild', 'severe', 'acute', 'chronic', 'primary', 'secondary', 'early', 'late',
                      'stage', 'grade', 'high', 'low', 'type', 'wild type', 'mutant', 'carrier', 'carriers'}

_DIS_NA_VALUES = {'', 'na', 'n a', 'nan', 'none', 'null', 'not applicable', 'not available', 'unknown', 'not specified',
                  'not reported', 'missing', 'missing value', 'other', 'others', 'unspecified'}
_DIS_QC_RE = re.compile(r'\b(?:qcs?|pool|pools|pooled|blanks?|standards?|nist|srm|solvent|water|buffer|calibra\w*|'
                        r'quality controls?|system suitability|process controls?|extraction controls?|method controls?|'
                        r'reference material|reference sample|internal standards?|tech(?:nical)? rep\w*|matrix blank)\b')
_DIS_KEY_TISSUE_RE = re.compile(r'\b(?:sample source|sample type|sampletype|sample types|tissue|tissues|tissue type|organ|'
                                r'organs|matrix|specimen|specimens|histology|biospecimen|material|body site|sample matrix|'
                                r'biofluid|fluid)\b')
_DIS_KEY_IGNORE_RE = re.compile(r'\b(?:time|times|timepoint|timepoints|day|days|week|weeks|month|months|hour|hours|hr|'
                                r'hrs|h|min|mins|minutes|age|ages|sex|gender|batch|order|replicate|replicates|rep|run|'
                                r'plate|well|position|weight|bmi|race|ethnicity|ethnic|smoking|smoker|smoke|diet|diets|'
                                r'concentration|conc|dose|dosage|volume|tracer|labeling|labelling|label|isotope|isotopic|'
                                r'medium|media|extraction|polarity|mode|instrument|date|year|years|site|location|region|'
                                r'country|center|centre|visit|fraction|injection|id|ids|raw|file|filename|height|length|'
                                r'ph|temperature|vial|tube|storage|cell|cells|line|lines|clone|organoid|organoids|culture|'
                                r'passage|zt|zeitgeber|compartment|ionization|column|sequence|well|bed|room|cage|tank|'
                                r'pen|farm|family|sibling|twin|mother|father|maternal|paternal|offspring|parity)\b')
_DIS_KEY_STRICT_RE = re.compile(r'\b(?:disease|diseases|diagnosis|diagnoses|diagnosed|health|healthy|case|cases|'
                                r'casestatus|patient|patients|clinical|pathology|pathological|disorder|infection|'
                                r'infected|illness)\b')
_DIS_KEY_TREAT_RE = re.compile(r'\b(?:treatment|treatments|treated|treat|drug|drugs|exposure|exposed|intervention|'
                               r'interventions|compound|compounds|stimulus|stimulation|genotype|genotypes|strain|strains|'
                               r'knockout|transgene|transgenic|mutation|mutations|mutant|operation|surgery|inhibitor|'
                               r'agent|irradiation|radiation|protocol|vaccine|vaccination|infusion|supplement|'
                               r'supplementation|challenge|transfection|overexpression|knockdown|variant|allele)\b')

# groups that are neither cases of the study disease nor healthy (other patients, risk
# groups, benign lesions, recovered individuals, test-negative individuals)
_DIS_OTHER_RE = re.compile(r'\b(?:disease controls?|patient controls?|hospital controls?|symptomatic controls?|'
                           r'other diseases?|non \w+ patients?|non\w+ patients?|pre ?diabet\w*|pre ?cancer\w*|'
                           r'pre ?malignan\w*|pre ?neoplas\w*|pre ?symptomatic|borderline|benign|polyps?|'
                           r'adenomas?|dysplas\w*|at risk|high risk|risk|suspected|suspect|possible|probable|remission|'
                           r'recovered|recovering|recovery|convalescen\w*|survivors?|history|former|resolved|cured|'
                           r'unhealthy|non ?healthy|carriers?|exposed|contacts?|relatives?|siblings?|'
                           r'unaffected relatives?|latent|prone|susceptible|discharged|survived|deceased|died|exitus)\b')
# 'negative' that is about the measurement, not a test result
_DIS_NOT_STATUS_RE = re.compile(r'\b(?:mode|ion|ions|esi|hilic|rp|c18|lipidomics|polarity|scaffold|control|controls)\b')
_DIS_MOCK_RE = re.compile(r'\b(?:mock|sham|placebo)\b')
_DIS_PLAIN_CONTROL_RE = re.compile(r'(?:the )?(?:controls?|ctrls?|ctl|ctr|ctrl group|con|cont|hc|hcs|nc|normal|normals|'
                                   r'normal controls?|nondiseased|non diseased)'
                                   r'(?: (?:group|groups|subjects?|individuals?|samples?|participants?|cohort|population|'
                                   r'people|persons?|volunteers?|donors?|mice|mouse|animals?|rats?|dogs?|women|men|'
                                   r'children|adults?))?(?: ?\d{1,3})?')
_DIS_MODEL_RE = re.compile(r'\b(?:wild ?type|wt|vehicle|veh|untreated|non ?treated|baseline|ambient|normoxi\w*|chow|'
                           r'naive|dmso|pbs|saline|unstimulated|diet|scrambled|scr|empty vector|knock ?out|ko|transgenic|'
                           r'treated|treatment|drug|dose|vector|sham|mock|placebo|pre|post|before|after|day|week|'
                           r'month|time|timepoint|surgery|operated|injected|exposure|exposed|air|water|food)\b')
_DIS_HEALTHY_RE = re.compile(r'(?<!un)(?<!non )(?<!non)(?<!not )\b(?:healthy|healty|heatlhy|helathy|healthly|heathy|'
                             r'disease ?free|free of disease|non ?diseased|'
                             r'normal (?:subjects?|volunteers?|individuals?|donors?|participants?|people|persons?|'
                             r'adults?|children|child|women|men|cohort|population))\b')
_DIS_HEALTHY_TISSUE_RE = re.compile(r'(?<!un)(?<!non )(?<!non)\b(?:healthy|healty|heatlhy) (?:mouse|mice|rats?|donors?|'
                                    r'subjects?|volunteers?|individuals?|controls?|people|persons?|participants?|'
                                    r'dogs?|animals?|humans?|adults?|children|women|men|cohort)\b')
# values describing a treatment, not the individual's health
_DIS_TREATED_RE = re.compile(r'\b(?:treated|treatment|injected|injection|dose|dosed|mg|kg|ug|um|nm|mm|exposed|exposure|'
                             r'infused|fed|diet|drug|vehicle|dmso|stimulated|transfected|knock ?out|ko)\b')
_DIS_CASE_WORD_RE = re.compile(r'\b(?:patients?|cases?|diseased?|affected|sick|ill)\b')
_DIS_HUMAN_CASE_WORD_RE = re.compile(r'\b(?:patients?|cases?)\b')
_DIS_CANCER_WORD_RE = re.compile(r'\b(?:tumou?rs?|tumou?ral|intratumou?ral|tumou?r ?bearing|cancers?|cancerous|carcinomas?|'
                                 r'adenocarcinomas?|malignan\w*|metasta\w*|neoplas\w*|sarcomas?|lymphomas?|leuka?emi\w*|'
                                 r'melanomas?|gliomas?|glioblastomas?|\w+blastomas?|myelomas?|mesotheliomas?|xenografts?)\b')
# samples of cultured / isolated cells, organoids or media: not an individual's sample
_DIS_CULTURE_RE = re.compile(r'\b(?:cells?|cell lines?|cultures?|cultured|organoids?|spheroids?|lysates?|supernatants?|'
                             r'media|medium|ipscs?|hipscs?|fibroblasts?|macrophages?|monocytes?|lymphocytes?|pbmcs?|'
                             r'neurons?|hepatocytes?|astrocytes?|microglia|myotubes?|explants?|in vitro|exosomes?|'
                             r'extracellular vesicles?|evs?|mitochondria|organelles?|tdcls?)\b')
# patient tissue next to a tumour ('adjacent normal', 'Normal area of the pancreas with tumor')
_DIS_ADJACENT_RE = re.compile(r'\b(?:adjacent|paired normal|matched normal|para ?(?:carcinoma|cancer|tumou?r)\w*|'
                              r'paracancerous|peri ?tumou?r\w*|distant normal|normal\b.*\bwith\b)\b')
_DIS_INFECTION_WORD_RE = re.compile(r'\b(?:infected|infection|infections|infested|infestation)\b')
_DIS_POSITIVE_RE = re.compile(r'\bpositive\b|\bpos\b|\+$')


def _norm_doid_name(name):
    return _norm_value(name).replace('+', ' ').strip()


def _doid_index(doid_table, path_to_csvs):
    """Name, label and hierarchy indexes of the DOID table (cached)."""
    key = (id(doid_table), len(doid_table), str(path_to_csvs))
    if key in _DOID_INDEX_CACHE:
        return _DOID_INDEX_CACHE[key]
    names, norm_names, labels, id_names = {}, {}, {}, {}
    for idx, label, syn in doid_table[['DOIDOntologyIndex', 'Label', 'Synonym']].itertuples(index=False):
        if not isinstance(label, str) or not isinstance(idx, str) or label.lower().startswith('obsolete'):
            continue
        labels[idx] = label
        for name in {label, syn}:
            if isinstance(name, str) and name.strip():
                names.setdefault(name.strip().lower(), set()).add(idx)
                n = _norm_doid_name(name)
                norm_names.setdefault(n, set()).add(idx)
                id_names.setdefault(idx, set()).add(n)
    children = {}
    try:
        if DOID_PARENTS_CSV == 'none':
            raise FileNotFoundError('no --path_to_doid_parents_csv given')
        parents = pd.read_csv(DOID_PARENTS_CSV, dtype=str)
        for child, parent in parents[['DOIDOntologyIndex', 'ParentDOIDOntologyIndex']].itertuples(index=False):
            children.setdefault(parent, set()).add(child)
    except Exception as e:
        print(f"Could not read the DOID hierarchy: {type(e).__name__}: {e}")
    aliases = {}
    try:
        al = pd.read_csv(os.path.join(path_to_csvs, 'MWB_disease_aliases.csv'), dtype=str, keep_default_na=False)
        for mwb, redu in al[['MWB', 'REDU']].itertuples(index=False):
            aliases[mwb.strip().lower()] = redu.strip()
    except Exception as e:
        print(f"Could not read the disease aliases: {type(e).__name__}: {e}")
    label_ids = {}
    for idx, label in labels.items():
        label_ids.setdefault(label.lower(), set()).add(idx)
    index = {'names': names, 'norm_names': norm_names, 'labels': labels, 'id_names': id_names,
             'children': children, 'aliases': aliases, 'label_ids': label_ids, 'desc': {}}
    _DOID_INDEX_CACHE.clear()
    _DOID_INDEX_CACHE[key] = index
    return index


def _doid_descendants(index, doid):
    """doid and all its DOID subclasses."""
    if doid not in index['desc']:
        seen, todo = {doid}, [doid]
        while todo:
            for c in index['children'].get(todo.pop(), ()):
                if c not in seen:
                    seen.add(c)
                    todo.append(c)
        index['desc'][doid] = seen
    return index['desc'][doid]


def mwb_disease_categories(disease_json):
    """Disease categories of a Workbench REST 'disease' response; 'A; B' gives two."""
    if isinstance(disease_json, dict):
        disease_json = [disease_json] if 'Disease' in disease_json else list(disease_json.values())
    if not isinstance(disease_json, list):
        return []
    out = []
    for item in disease_json:
        if isinstance(item, dict) and isinstance(item.get('Disease'), str):
            for name in item['Disease'].split(';'):
                name = re.sub(r'\s+', ' ', name).strip()
                if name and name not in out:
                    out.append(name)
    return out


def _disease_name_regex(index, doid):
    """Regex for the names (labels, exact synonyms, adjectives) of doid and its subclasses."""
    if doid in _DISEASE_RE_CACHE:
        return _DISEASE_RE_CACHE[doid]
    names = set()
    for d in _doid_descendants(index, doid):
        for n in index['id_names'].get(d, ()):
            names.add(n)
            names.update(_DISEASE_ADJECTIVES.get(n, []))
    names = {n for n in names if len(n) >= 3 and n not in _DISEASE_NAME_SKIP and not n.isdigit()}
    rx = None
    if names:
        rx = re.compile(r'(?<![a-z0-9])(?:' + '|'.join(re.escape(n) for n in sorted(names, key=len, reverse=True))
                        + r')(?![a-z0-9])')
    _DISEASE_RE_CACHE[doid] = rx
    return rx


def prepare_study_disease(raw_names, doid_table, path_to_csvs='translation_sheets'):
    """DOID terms of a study's disease categories (exact label/exact synonym, or an
    exact-meaning alias from MWB_disease_aliases.csv). Categories without a DOID
    term are skipped. Returns None if no category has one."""
    if doid_table is None or len(doid_table) == 0 or not raw_names:
        return None
    index = _doid_index(doid_table, path_to_csvs)
    cats, skipped = [], []
    for raw in raw_names:
        low = raw.strip().lower()
        if low in index['aliases']:
            target = index['aliases'][low]
            ids = index['label_ids'].get(target.lower(), set()) if target else set()
        else:
            ids = index['names'].get(low, set())
        if len(ids) != 1:
            skipped.append(raw)
            continue
        doid = next(iter(ids))
        if doid == 'DOID:4' or any(c['id'] == doid for c in cats):
            continue
        desc = _doid_descendants(index, doid)
        cats.append({'raw': raw, 'id': doid, 'label': index['labels'][doid], 'desc': desc,
                     'cancer': doid in _doid_descendants(index, _DOID_CANCER),
                     'infectious': doid in _doid_descendants(index, _DOID_INFECTIOUS),
                     'name_re': _disease_name_regex(index, doid)})
    if skipped:
        print(f"Study disease categories without a DOID term: {skipped}")
    if not cats:
        return None
    return {'cats': cats, 'index': index}


def _pick_disease(study, cat_ids, specific):
    """'label|DOID' for a case: the specific term, else the single (deepest) category."""
    index = study['index']
    if len(specific) == 1:
        d = next(iter(specific))
        return f"{index['labels'][d]}|{d}"
    cats = [c for c in study['cats'] if c['id'] in cat_ids]
    # keep the deepest of nested categories ('Cancer' and 'Lung cancer' -> lung cancer)
    deepest = [c for c in cats if all(c['id'] in o['desc'] for o in cats)]
    if len(deepest) == 1:
        return f"{deepest[0]['label']}|{deepest[0]['id']}"
    return None


def _disease_key_class(key):
    k = _norm_key(key)
    if not k:
        return 'ignore'
    if _DIS_KEY_TISSUE_RE.search(k):
        return 'tissue'
    if _DIS_KEY_IGNORE_RE.search(k):
        return 'ignore'
    if _DIS_KEY_STRICT_RE.search(k):
        return 'strict'
    if _DIS_KEY_TREAT_RE.search(k):
        return 'treat'
    return 'group'


def classify_disease_value(value, key_class, study, key=''):
    """Class of one factor value: ('case', category ids, specific DOID ids), ('healthy',),
    ('control',) (plain 'Control', healthy only in context), ('otherdis',) (other patients,
    risk groups, ...), ('other',), ('qc',) or ('na',)."""
    v = _norm_value(value)
    if v in _DIS_NA_VALUES or v.replace(' ', '').isdigit():
        return ('na',)
    if qc_sample_type(key, value) or _DIS_QC_RE.search(v):
        return ('qc',)
    cats, index = study['cats'], study['index']
    all_ids = {c['id'] for c in cats}
    # generic case words ('patients', 'cases') only give the study-level categories, not
    # terms found in the sample values (health conditions)
    word_ids = {c['id'] for c in cats if not c.get('from_values')}
    if _DIS_MOCK_RE.search(v):
        return ('other',)
    if key_class == 'tissue':
        cancer_ids = {c['id'] for c in cats if c['cancer']}
        if cancer_ids and _DIS_ADJACENT_RE.search(v):
            return ('case', cancer_ids, set())
    # 'non-diabetic', 'tumor-free', 'without asthma', 'uninfected': explicitly without the disease
    vn = v.replace('+', ' ').strip()
    neg = re.match(r'^(?:non|un|not|no|without|free of)\s?(.*)$', vn)
    free = re.match(r'^(.*?)\s?free(?: (?:subjects?|individuals?|controls?|mice|group|volunteers?|participants?))?$', vn)
    def _names_disease(text, full):
        m = (lambda rx: rx.fullmatch(text)) if full else (lambda rx: rx.match(text))
        return any(c['name_re'] is not None and m(c['name_re']) for c in cats) \
            or (any(c['cancer'] for c in cats) and m(_DIS_CANCER_WORD_RE)) \
            or (any(c['infectious'] for c in cats) and m(re.compile(r'(?:infected|infection)'))) \
            or m(re.compile(r'(?:diseased?|affected)'))
    for rest in [m.group(1) for m in (neg, free) if m]:
        rest = re.sub(r' (?:controls?|subjects?|individuals?|mice|group|volunteers?|participants?|donors?|women|men|'
                      r'children|adults?)$', '', rest).strip()
        if not rest or _DIS_OTHER_RE.search(rest) or _DIS_CASE_WORD_RE.search(rest):
            continue
        if key_class == 'tissue' or rest.startswith('malignan'):
            # non-tumour tissue / non-malignant: not a statement about a healthy individual
            if _names_disease(rest, full=False):
                return ('other',) if key_class == 'tissue' else ('otherdis',)
            continue
        if _names_disease(rest, full=True):
            # health conditions: 'non-septic' / 'angina-free' people are often other patients
            return ('other',) if study.get('no_negation_healthy') else ('healthy',)
    # 'patients without CKD', 'participant without the virus infection': not cases
    for m in re.finditer(r'\b(?:(?:without|w o|no|not|free of)\s+|non\s?)(?:the |any |a )?(.+)$', vn):
        if _names_disease(m.group(1), full=False):
            return ('otherdis',)
    if _DIS_OTHER_RE.search(v) or (re.search(r'\bnegative\b', v) and not _DIS_NOT_STATUS_RE.search(v)):
        return ('otherdis',)
    # exact DOID name, possibly with 'patients'/'group' around it
    forms = {vn, re.sub(r'\s*(?:patients?|cases?|group|subjects?|individuals?|samples?|cohort|mice|tissues?)$', '', vn),
             re.sub(r'^(?:patients? with|subjects? with|individuals? with)\s*', '', vn)}
    for f in forms:
        ids = index['norm_names'].get(f, set()) - {'DOID:4'}
        if len(ids) == 1 and f not in _DISEASE_NAME_SKIP:
            d = next(iter(ids))
            under = {c['id'] for c in cats if d in c['desc']}
            if under:
                return ('case', under, {d} if d not in all_ids else set())
            # another disease; short abbreviations ('con', 'dss') are too ambiguous to tell
            if len(f) >= 5:
                return ('otherdis',)
    hit = {c['id'] for c in cats if c['name_re'] is not None and c['name_re'].search(v)}
    if any(c['cancer'] for c in cats) and _DIS_CANCER_WORD_RE.search(v):
        hit |= {c['id'] for c in cats if c['cancer']}
    if any(c['infectious'] for c in cats) and _DIS_INFECTION_WORD_RE.search(v):
        hit |= {c['id'] for c in cats if c['infectious']}
    if _DIS_CASE_WORD_RE.search(v):
        hit |= word_ids
    if not hit and _DIS_POSITIVE_RE.search(v):
        k = _norm_value(key)
        if any(c['name_re'] is not None and c['name_re'].search(k) for c in cats) or \
                (any(c['infectious'] for c in cats) and _DIS_INFECTION_WORD_RE.search(k)):
            hit = word_ids
    if hit:
        return ('case', hit, set())
    if len(v) <= 60 and not _DIS_TREATED_RE.search(v):
        if key_class == 'tissue':
            if _DIS_HEALTHY_TISSUE_RE.search(v):
                return ('healthy',)
        elif _DIS_HEALTHY_RE.search(v):
            return ('healthy',)
    if _DIS_PLAIN_CONTROL_RE.fullmatch(v):
        return ('control',)
    return ('other',)


def study_disease_assignment(df, study):
    """Per sample ('filename'): (HealthStatus, 'DOID label|DOID id' or None) from the
    study-level disease. df has one row per sample factor (filename, Key, Value,
    SubjectIdentifierAsRecorded). Returns {filename: (health, doid)}; health is
    'unhealthy (NOS)', 'healthy' or 'healthy if human'."""
    if study is None or len(df) == 0:
        return {}
    rows = df[['filename', 'Key', 'Value']].dropna(subset=['filename']).copy()
    rows['Key'] = rows['Key'].astype(str).str.strip().str.lower()
    rows['Value'] = rows['Value'].astype(str).str.strip().str.lower()
    rows = rows.drop_duplicates()
    rows['kc'] = rows['Key'].map(_disease_key_class)
    rows = rows[rows['kc'] != 'ignore']
    if len(rows) == 0:
        return {}
    cls = {}
    for key, kc, val in rows[['Key', 'kc', 'Value']].drop_duplicates().itertuples(index=False):
        cls[(key, val)] = classify_disease_value(val, kc, study, key=key)
    # plain 'Control' is a healthy control only under a disease-status key, or when every
    # other group of that factor is a case group (in treatment/tissue factors: a human
    # study or a 'patient'/'case' group)
    control_ok = {}
    for (key, kc), g in rows.groupby(['Key', 'kc']):
        vals = set(g['Value'])
        kinds = [cls[(key, v)][0] for v in vals]
        has_case = 'case' in kinds
        pure = all(k in ('case', 'control', 'healthy', 'qc', 'na') for k in kinds)
        human_case = any(cls[(key, v)][0] == 'case' and _DIS_HUMAN_CASE_WORD_RE.search(v) for v in vals)
        if study.get('explicit_healthy_only'):
            # health conditions: plain 'Control' groups are often other patients
            control_ok[key] = None
        elif kc == 'strict':
            control_ok[key] = 'healthy'
        elif kc == 'group':
            control_ok[key] = 'healthy' if has_case and pure else None
        else:
            if has_case and pure and human_case:
                control_ok[key] = 'healthy'
            elif has_case and pure and kc == 'treat':
                control_ok[key] = 'healthy if human'
            else:
                control_ok[key] = None
    # per sample
    per_sample = {}
    cultured = set(rows.loc[(rows['kc'] == 'tissue') & rows['Value'].map(lambda x: bool(_DIS_CULTURE_RE.search(_norm_value(x)))),
                            'filename'])
    for fn, g in rows.groupby('filename'):
        if fn in cultured:
            continue
        case_ids, specific, healthy, otherdis, case = set(), set(), None, False, False
        for key, val in g[['Key', 'Value']].itertuples(index=False):
            c = cls[(key, val)]
            if c[0] == 'case':
                case = True
                case_ids |= c[1]
                specific |= c[2]
            elif c[0] == 'healthy':
                healthy = 'healthy'
            elif c[0] == 'control' and control_ok.get(key):
                healthy = healthy or control_ok[key]
            elif c[0] == 'otherdis':
                otherdis = True
        if case and healthy:
            per_sample[fn] = ('conflict',)
        elif case:
            per_sample[fn] = ('case', case_ids, specific)
        elif healthy and not otherdis:
            per_sample[fn] = ('healthy', healthy)
        elif otherdis:
            per_sample[fn] = ('otherdis',)
    # per individual
    subj = df[['filename', 'SubjectIdentifierAsRecorded']].dropna(subset=['filename']).copy() \
        if 'SubjectIdentifierAsRecorded' in df.columns else pd.DataFrame(columns=['filename', 'SubjectIdentifierAsRecorded'])
    subj['sid'] = subj['SubjectIdentifierAsRecorded'].fillna('').astype(str).str.strip().str.lower()
    subj = subj[~subj['sid'].isin(['', '-', 'na', 'n/a', 'nan', 'none'])].drop_duplicates(['filename', 'sid'])
    # a sample listed under several Subject IDs is judged alone
    subj = subj[~subj['filename'].duplicated(keep=False)]
    out = {}
    done = set()
    for sid, g in subj.groupby('sid'):
        files = list(g['filename'])
        if len(files) < 2:
            continue
        states = [per_sample.get(f, ('none',)) for f in files]
        kinds = {s[0] for s in states}
        done.update(files)
        if 'case' in kinds:
            ids = set().union(*[s[1] for s in states if s[0] == 'case'])
            spec = set().union(*[s[2] for s in states if s[0] == 'case'])
            doid = _pick_disease(study, ids, spec)
            if 'healthy' in kinds or 'conflict' in kinds:
                # contradicting group labels for one Subject ID: only the case samples
                for f, s in zip(files, states):
                    if s[0] == 'case':
                        out[f] = ('unhealthy (NOS)', doid)
            else:
                for f in files:
                    out[f] = ('unhealthy (NOS)', doid)
        elif 'healthy' in kinds and not kinds & {'otherdis', 'conflict'}:
            h = 'healthy' if any(s[0] == 'healthy' and s[1] == 'healthy' for s in states) else 'healthy if human'
            for f in files:
                out[f] = (h, None)
    for f, s in per_sample.items():
        if f in done:
            continue
        if s[0] == 'case':
            out[f] = ('unhealthy (NOS)', _pick_disease(study, s[1], s[2]))
        elif s[0] == 'healthy':
            out[f] = (s[1], None)
    return {f: r for f, r in out.items() if f not in cultured}


# ---------------------------------------------------------------------------------
# Health conditions (HPOCommonName: HPO 'Phenotypic abnormality' terms that are not
# diseases, e.g. Sepsis, Cachexia) and exposures (Exposure: ECTO, for plants PECO, e.g.
# 'exposure to ozone'). Sources: Workbench REST disease categories without a DOID term,
# and sample values naming a term exactly. Health conditions use the per-individual
# rules of the study-level disease above; an exposure only goes to the samples of the
# exposed group (never to controls, vehicle, sham, ambient air, pre-exposure samples).
# Exact labels / exact synonyms only (case and separators ignored); only fills gaps.
# ---------------------------------------------------------------------------------
# prepared HPO table; set from --path_to_hpo_csv
HPO_TABLE = None
# prepared ECTO/PECO table; set from --path_to_exposure_csv
EXPOSURE_TABLE = None
# health conditions / exposure of the study being converted (MWB_to_REDU_study_wrapper)
_STUDY_HPO = None
_STUDY_EXPOSURE = None
_COND_INDEX_CACHE = {}

_HPO_ADJECTIVES = {'sepsis': ['septic'], 'cachexia': ['cachectic'], 'malnutrition': ['malnourished'],
                   'anorexia': ['anorexic']}
# key of a sample value that names a health condition of the individual
_HPO_KEY_SKIP_RE = re.compile(r'\b(?:death|died|cause|history|family|maternal|paternal|mother|father|offspring|'
                              r'parent|parents|donor|recipient|comorbidit\w*|medications?|complications?|follow|'
                              r'outcomes?|response|toxicit\w*|adverse|events?)\b')
# HPO names that sample values also use for a sample material
_HPO_VALUE_SKIP = {'ascites', 'edema', 'oedema', 'pleural effusion', 'effusion', 'abscess', 'scar', 'wound',
                   # states of hibernating animals, no condition (HPO has 'torpor' as Lethargy)
                   'torpor', 'hibernation', 'hibernating'}
# keys whose values can name an exposure of the individual
_EXP_KEY_STRONG_RE = re.compile(r'\b(?:treat\w*|expos\w*|drugs?|compounds?|chemicals?|agents?|interventions?|stimul\w*|'
                                r'challeng\w*|toxi\w*|pollut\w*|irradiat\w*|radiation|inhal\w*|administ\w*|medicat\w*|'
                                r'inocul\w*|infect\w*|gavage|injections?|injected|dosing|therap\w*|regimens?)\b')
_EXP_KEY_GENERIC_RE = re.compile(r'^(?:groups?|factors?|conditions?|class|classes|categor(?:y|ies)|arms?|study arms?|'
                                 r'study groups?|sample groups?|experimental (?:factors?|groups?|conditions?)|'
                                 r'experiment groups?|cohort groups?)$')
_EXP_KEY_SKIP_RE = re.compile(r'\b(?:diets?|dietary|feed|feeding|food|chow|tracers?|label\w*|isotop\w*|medium|media|'
                              r'cells?|culture\w*|genotypes?|strains?|species|organisms?|batch|vehicle|solvent|'
                              r'extraction|controls?)\b')
# time-like keys: their numbers are no exposure dose
_EXP_TIME_KEY_RE = re.compile(r'\b(?:time|times|timepoints?|day|days|hours?|hr|hrs|h|weeks?|months?|post|after|since|'
                              r'duration|min|mins|minutes|age|visit)\b')
_EXP_STATE_KEY_RE = re.compile(r'\b(?:status|stage|phase|period|collection|sampling|state|visit|timing)\b')
_EXP_DOSE_RE = re.compile(r'(?<![a-z0-9.])(\d+(?:\.\d+)?)\s*(?:mg/kg(?:/d(?:ay)?)?|mg/kg/day|g/kg|mg|µg|ug|ng|g/l|mg/l|'
                          r'mg/ml|ug/ml|µg/ml|µg/m3|ug/m3|mg/m3|µm|um|nm|mm|m|mmol/l|µmol/l|umol/l|ppm|ppb|ppt|c?gy|'
                          r'%|mpk|x|h|hr|hrs|hours?|d|days?|weeks?|wks?|min)(?![a-z0-9])')
_EXP_GY_RE = re.compile(r'(?<![a-z0-9.])(\d+(?:\.\d+)?)\s*c?gy(?![a-z])')
# values that are a control, a vehicle or too generic to be an exposure
_EXP_VALUE_SKIP = {'water', 'air', 'oil', 'feces', 'faeces', 'blood', 'protein', 'proteins', 'solvent', 'drug', 'drugs',
                   'toxin', 'toxins', 'peptide', 'peptides', 'patient', 'patients', 'chemical', 'chemicals', 'food',
                   'saline', 'vehicle', 'control', 'placebo', 'sham', 'mock', 'media', 'medium', 'diet', 'ambient',
                   'filtered air', 'normal air', 'clean air', 'room air', 'ambient air', 'pbs', 'dmso', 'buffer',
                   'serum', 'plasma', 'urine', 'tissue', 'cells', 'culture', 'light', 'dark', 'sugar', 'salt',
                   'mixture', 'mix', 'extract', 'compound', 'agent', 'treatment', 'surgery', 'stress', 'exposure',
                   'plaster'}
# exposure terms too generic to say what the individual was exposed to
_EXP_LABEL_SKIP = {'control exposure', 'plant exposure', 'exposure to organism', 'exposure to soil',
                   'exposure to solvent', 'exposure to water', 'exposure to air', 'watering exposure',
                   'temperature exposure', 'light intensity exposure', 'light quantity exposure',
                   'water environment exposure', 'exposure to drug', 'exposure to chemical entity'}
# explicitly not exposed
_EXP_NEG_VALUES = {'control', 'controls', 'ctrl', 'ctl', 'con', 'cont', 'vehicle', 'veh', 'untreated', 'non treated',
                   'not treated', 'none', 'no', 'n', 'air', 'ambient', 'ambient air', 'filtered air', 'normal air',
                   'room air', 'clean air', 'fa', 'non irradiated', 'nonirradiated', 'unirradiated', 'not irradiated',
                   'unexposed', 'non exposed', 'not exposed', 'no exposure', 'naive', 'sham', 'mock', 'placebo', 'normoxia',
                   'normoxic', 'baseline'}
_EXP_PRE_RE = re.compile(r'(?:pre|pre ?(?:exposure|exposed|treatment|treated|irradiation|irradiated|dose|dosing|op|'
                         r'operative|intervention|challenge|infection|immune|injection|ir|rt)|before(?: \w+)?|baseline|'
                         r'day 0|d ?0|0 ?h|0 ?hr?s?|0 ?d|t ?0|time 0|week 0|w ?0)')
_EXP_FIRST_VISIT_RE = re.compile(r'(?:v|t|tp|visit|time|timepoint|time point|sample|s)? ?0?1')
# cell lines named without the word 'cells' (in vitro treatments are no exposure of an individual)
_EXP_CELL_LINE_RE = re.compile(r'\b(?:hek ?293\w*|hela|hep ?g2|huh ?7|mcf ?7|mcf ?10a|mda mb \d+|a549|caco ?2|ht ?29|'
                               r'hct ?116|sw ?480|sw ?620|jurkat|thp ?1|u937|raw ?264\w*|3t3\w*|cho|cos ?7|u87\w*|'
                               r'u251|sh ?sy5y|pc ?12|pc ?3|du ?145|lncap|22rv1|huvecs?|hl ?60|k562|c2c12|l6|ins ?1|'
                               r'min6|bv ?2|ipsc\w*|hepa ?1 ?6|panc ?1|mia paca ?2|a375|b16\w*|ct26|4t1|llc|e ?coli|'
                               r'saccharomyces|yeast)\b')
# studies of cultured / isolated material: their 'treatments' are no exposure of an individual
_EXP_IN_VITRO_RE = re.compile(r'\b(?:cells|cell lines?|cell cultures?|in vitro|cultured|cultures|organoids?|spheroids?|'
                              r'ipscs?|hipscs?|induced neurons|fibroblasts|macrophages|hepatocytes|myotubes|'
                              r'cardiomyocytes|adipocytes|neurons|microglia|astrocytes|keratinocytes|explants?)\b')
_EXP_IN_VIVO_RE = re.compile(r'\b(?:mice|mouse|murine|rats?|rodents?|patients|participants|subjects|volunteers|children|'
                             r'adults|women|men|humans|individuals|people|cohorts?|population|zebrafish|embryos|larvae|'
                             r'fish|pigs?|piglets|cows|cattle|dogs|sheep|monkeys|primates|animals|birds|chickens|broilers|'
                             r'bees|flies|drosophila|worms|elegans|plants?|seedlings|leaves|in vivo|serum|plasma|urine|'
                             r'feces|fecal|faecal|stool|blood|tissues?|biopsy|biopsies|milk|xenografts?)\b')
# set per analysis from the mwTab (create_dataframe_outer_dict)
_ANALYSIS_IN_VITRO = False


def analysis_in_vitro(title, sample_type='', subject_type=''):
    """True if the study title (or the subject type) describes cultured or isolated
    material only ('... Amlodipine in Lung Cancer Cells'). (The sample type is not used:
    'Media' / 'Liver cells' there also describe embryo media and tissue.)"""
    if re.search(r'\b(?:cells?|cell lines?|cultured|cultures?|organoids?|ipscs?)\b', _norm_value(subject_type)):
        return True
    t = re.sub(r'[^a-z0-9]+', ' ', str(title or '').lower())
    return bool(_EXP_IN_VITRO_RE.search(t)) and not _EXP_IN_VIVO_RE.search(t)


_EXP_POS_WORDS = {'yes', 'y', 'exposed', 'treated', 'irradiated', 'positive', 'true'}
# wording of the study-level exposure in sample values
_EXP_CATEGORY_FORMS = {'radiation': {'radiation', 'irradiated', 'irradiation', 'radiated', 'ionizing radiation',
                                     'ionising radiation', 'total body irradiation', 'whole body irradiation'},
                       'ozone': {'ozone'}}
_EXP_RADIATION_SPECIFIC = {'gamma': 'gamma radiation', 'gamma ray': 'gamma radiation', 'gamma rays': 'gamma radiation',
                           'gamma radiation': 'gamma radiation', 'gamma irradiation': 'gamma radiation',
                           'gamma irradiated': 'gamma radiation', 'x ray': 'x-ray radiation',
                           'x rays': 'x-ray radiation', 'xray': 'x-ray radiation', 'x ray radiation': 'x-ray radiation',
                           'x ray irradiation': 'x-ray radiation', 'x ray irradiated': 'x-ray radiation'}


def _allowed_values(column):
    return set(((globals().get('allowedTerm_dict') or {}).get(column) or {}).get('allowed_values') or [])


def _ontology_name_index(table, id_col, allowed):
    """normalised label/exact synonym -> {(label, id): is_label}, names longer than 3
    characters of allowed labels only; and id -> normalised names."""
    names, id_names = {}, {}
    if table is None or len(table) == 0:
        return names, id_names
    for idx, label, syn in table[[id_col, 'Label', 'Synonym']].itertuples(index=False):
        if not isinstance(label, str) or not isinstance(idx, str) or (allowed and label not in allowed):
            continue
        for name in {label, syn}:
            if not isinstance(name, str):
                continue
            n = _norm_doid_name(name)
            if len(n) > 3:
                names.setdefault(n, {})
                names[n][(label, idx)] = names[n].get((label, idx), False) or name == label
                id_names.setdefault(idx, set()).add(n)
    return names, id_names


def _condition_index(doid_table, path_to_csvs):
    """HPO / exposure name indexes and the condition aliases (cached)."""
    hpo_table = HPO_TABLE
    exp_table = EXPOSURE_TABLE
    key = (id(hpo_table), id(exp_table), id(doid_table), str(path_to_csvs))
    if key in _COND_INDEX_CACHE:
        return _COND_INDEX_CACHE[key]
    dindex = _doid_index(doid_table, path_to_csvs) if doid_table is not None and len(doid_table) else \
        {'names': {}, 'norm_names': {}, 'labels': {}, 'id_names': {}, 'children': {}, 'aliases': {}, 'label_ids': {},
         'desc': {}}
    hpo_names, hpo_id_names = _ontology_name_index(hpo_table, 'HPOOntologyIndex', _allowed_values('HPOCommonName'))
    # names that are also DOID names are diseases (DOIDCommonName); tumours likewise.
    # Sample values only give terms none of whose names is a DOID name ('heart failure'
    # is HPO 'Congestive heart failure', a disease in DOID)
    doid_names = set(dindex['norm_names'])
    disease_like = {i for i, ns in hpo_id_names.items() if ns & doid_names}
    for n in list(hpo_names):
        if n in doid_names or n in _DISEASE_NAME_SKIP or _DIS_CANCER_WORD_RE.search(n) or n.isdigit():
            del hpo_names[n]
    exp_names, _ = _ontology_name_index(exp_table, 'ExposureIndex', _allowed_values('Exposure'))
    for n in list(exp_names):
        exp_names[n] = {t: is_label for t, is_label in exp_names[n].items()
                        if not re.search(r'\b(?:diet|diets|feeding|food)\b', t[0].lower())}
        if not exp_names[n]:
            del exp_names[n]
    hpo_labels = {}
    for n, terms in hpo_names.items():
        for (label, idx) in terms:
            hpo_labels[idx] = label
    aliases = {}
    try:
        al = pd.read_csv(os.path.join(path_to_csvs, 'MWB_condition_aliases.csv'), dtype=str, keep_default_na=False)
        for mwb, col, redu in al[['MWB', 'Column', 'REDU']].itertuples(index=False):
            aliases[mwb.strip().lower()] = (col.strip(), redu.strip())
    except Exception as e:
        print(f"Could not read the condition aliases: {type(e).__name__}: {e}")
    index = {'doid': dindex, 'hpo_names': hpo_names, 'hpo_id_names': hpo_id_names, 'hpo_labels': hpo_labels,
             'hpo_disease_like': disease_like,
             'exp_names': exp_names, 'aliases': aliases}
    _COND_INDEX_CACHE.clear()
    _COND_INDEX_CACHE[key] = index
    return index


def _pick_exposure(cands, plant):
    """One exposure term of candidates {(label, id): is_label}: PECO for plants, else
    ECTO; then labels before synonyms. None if still ambiguous."""
    if not cands:
        return None
    # PECO (plant exposures) only for plants; generic terms ('control exposure') never
    pool = [t for t in cands if (plant or not t[1].startswith('PECO:')) and t[0].lower() not in _EXP_LABEL_SKIP]
    if not pool:
        return None
    for keep in (lambda t: t[1].startswith('PECO:') == plant, lambda t: cands[t]):
        sub = [t for t in pool if keep(t)]
        if sub:
            pool = sub
    if len({label.lower() for label, _ in pool}) == 1:
        label, idx = sorted(pool)[0]
        return f'{label}|{idx}'
    return None


def _exposure_term(name, cindex, plant):
    cands = {}
    for f in (name, 'exposure to ' + name, name + ' exposure'):
        cands.update(cindex['exp_names'].get(f, {}))
    return _pick_exposure(cands, plant)


def _doid_covers(raw, dindex):
    """True if the disease step gives this study category a DOID term."""
    low = raw.strip().lower()
    if low in dindex['aliases']:
        target = dindex['aliases'][low]
        return bool(target) and len(dindex['label_ids'].get(target.lower(), set())) == 1
    ids = dindex['names'].get(low, set())
    return len(ids) == 1 and ids != {'DOID:4'}


def _hpo_category(hp, cindex, from_values=False, extra_names=()):
    names = {n for n in cindex['hpo_id_names'].get(hp, set()) if n in cindex['hpo_names']} | set(extra_names)
    for n in list(names):
        names.update(_HPO_ADJECTIVES.get(n, []))
    names = {n for n in names if len(n) > 3 and n not in _DISEASE_NAME_SKIP}
    rx = re.compile(r'(?<![a-z0-9])(?:' + '|'.join(re.escape(n) for n in sorted(names, key=len, reverse=True))
                    + r')(?![a-z0-9])') if names else None
    return {'raw': hp, 'id': hp, 'label': cindex['hpo_labels'][hp], 'desc': {hp}, 'cancer': False,
            'infectious': False, 'name_re': rx, 'names': names, 'from_values': from_values}


def prepare_study_conditions(raw_names, doid_table, path_to_csvs='translation_sheets'):
    """HPO health conditions and the exposure of a study's REST disease categories that
    have no DOID term: exact HPO / ECTO / PECO label or exact synonym, or an exact-meaning
    alias from MWB_condition_aliases.csv (an empty alias target means no term).
    Returns ([(HPO id, alias names)], exposure dict or None)."""
    if not raw_names or ((HPO_TABLE is None or len(HPO_TABLE) == 0) and
                         (EXPOSURE_TABLE is None or len(EXPOSURE_TABLE) == 0)):
        return [], None
    cindex = _condition_index(doid_table, path_to_csvs)
    hpo_ids, exposures, alias_names = [], [], {}
    for raw in raw_names:
        if _doid_covers(raw, cindex['doid']):
            continue
        low = raw.strip().lower()
        n = _norm_doid_name(raw)
        hp, exp = set(), None
        if low in cindex['aliases']:
            col, target = cindex['aliases'][low]
            if not target:
                continue
            t = _norm_doid_name(target)
            if col == 'HPOCommonName':
                hp = {i for (label, i) in cindex['hpo_names'].get(t, {}) if label.lower() == target.lower()}
            elif col == 'Exposure':
                exp = t
        else:
            hp = {i for (_, i) in cindex['hpo_names'].get(n, {})}
            if not hp and _pick_exposure(cindex['exp_names'].get(n, {}), False):
                exp = n
        if len(hp) == 1:
            if next(iter(hp)) not in hpo_ids:
                hpo_ids.append(next(iter(hp)))
            if low in cindex['aliases'] and len(n) > 3 and n not in cindex['doid']['norm_names']:
                # the category's own wording names the condition in sample values ('Angina')
                alias_names.setdefault(next(iter(hp)), set()).add(n)
        elif exp is not None:
            exposures.append(exp)
    exposure = None
    if len(set(exposures)) == 1:
        name = exposures[0]
        term = _pick_exposure(cindex['exp_names'].get(name, {}), False)
        plant_term = _pick_exposure(cindex['exp_names'].get(name, {}), True)
        if term:
            words = set()
            for label in {term.split('|')[0], (plant_term or term).split('|')[0]}:
                w = re.sub(r'^exposure to |\s*exposure$', '', _norm_doid_name(label)).strip()
                words |= _EXP_CATEGORY_FORMS.get(w, {w})
            exposure = {'term': term, 'plant_term': plant_term or term, 'words': words,
                        'radiation': 'radiation' in words}
    if hpo_ids or exposure:
        print(f"Study health conditions: {[cindex['hpo_labels'][h] for h in hpo_ids]}; "
              f"exposure: {exposure['term'] if exposure else None}")
    return [(h, alias_names.get(h, set())) for h in hpo_ids], exposure


def _hpo_value_terms(rows, cindex):
    """HPO terms named exactly by a sample value under a disease/group-like key."""
    found = set()
    for key, kc, val in rows[['Key', 'kc', 'Value']].drop_duplicates().itertuples(index=False):
        if kc not in ('strict', 'group') or _HPO_KEY_SKIP_RE.search(_norm_key(key)):
            continue
        vn = _norm_doid_name(val)
        forms = {vn, re.sub(r'\s*(?:patients?|cases?|group|subjects?|individuals?|cohort|mice|rats?)$', '', vn),
                 re.sub(r'^(?:patients? with|subjects? with|individuals? with)\s*', '', vn)}
        for f in forms:
            if f in _HPO_VALUE_SKIP:
                continue
            ids = {i for (_, i) in cindex['hpo_names'].get(f, {})}
            if len(ids) == 1 and not ids & cindex['hpo_disease_like']:
                found |= ids
    return found


def study_condition_assignment(df, hpo_ids, doid_table, path_to_csvs='translation_sheets', study_has_doid=False):
    """Per sample ('filename'): (HealthStatus, 'HPO label|HPO id') from the study's HPO
    categories (hpo_ids: [(HPO id, extra names)]) and HPO terms named by sample values, with the per-individual rules of
    study_disease_assignment. 'healthy' only in studies without a DOID category."""
    if HPO_TABLE is None or len(HPO_TABLE) == 0 or len(df) == 0:
        return {}
    cindex = _condition_index(doid_table, path_to_csvs)
    # keys about other people, the past or outcomes say nothing about the condition
    df = df[~df['Key'].astype(str).str.lower().map(lambda k: bool(_HPO_KEY_SKIP_RE.search(_norm_key(k))))]
    rows = df[['filename', 'Key', 'Value']].dropna(subset=['filename']).copy()
    rows['Key'] = rows['Key'].astype(str).str.strip().str.lower()
    rows['Value'] = rows['Value'].astype(str).str.strip().str.lower()
    rows = rows.drop_duplicates()
    rows['kc'] = rows['Key'].map(_disease_key_class)
    rows = rows[rows['kc'] != 'ignore']
    if len(rows) == 0:
        return {}
    cats = [_hpo_category(h, cindex, extra_names=names) for h, names in hpo_ids]
    for h in sorted(_hpo_value_terms(rows, cindex) - {h for h, _ in hpo_ids}):
        cats.append(_hpo_category(h, cindex, from_values=True))
    if not cats:
        return {}
    dindex = cindex['doid']
    norm_names = dict(dindex['norm_names'])
    for c in cats:
        for n in c['names']:
            if n in cindex['hpo_names']:
                norm_names[n] = {c['id']}
    labels = dict(dindex['labels'])
    labels.update({c['id']: c['label'] for c in cats})
    index = dict(dindex, norm_names=norm_names, labels=labels)
    study = {'cats': cats, 'index': index, 'no_negation_healthy': True, 'explicit_healthy_only': True}
    out = study_disease_assignment(df, study)
    if study_has_doid:
        out = {f: r for f, r in out.items() if r[0] == 'unhealthy (NOS)'}
    return out


def _exposure_value(value):
    """(normalised value without doses and 'exposed to' wording, dose numbers)."""
    x = re.sub(r'\s+', ' ', str(value).strip().lower())
    nums = [float(m.group(1)) for m in _EXP_DOSE_RE.finditer(x)]
    y = _EXP_DOSE_RE.sub(' ', x)
    y = re.sub(r'\s+', ' ', y).strip(' -_,;:()')
    y = re.sub(r'^(?:exposed to|exposure to|treated with|treatment with|treated|exposed)\s+', '', y)
    y = re.sub(r'[\s_-]+(?:exposed|exposure|treated|treatment|group|administration|injection|injected)$', '', y)
    return _norm_doid_name(y), nums


def study_exposure_assignment(df, study_exposure, doid_table, path_to_csvs='translation_sheets'):
    """Per sample ('filename'): ('label|id' for non-plants, 'label|id' for plants) of what
    the sample's individual was exposed to: a sample value under a treatment/exposure key
    that names an ECTO/PECO term exactly ('ozone' -> 'exposure to ozone', 'Cadmium' ->
    'exposure to cadmium'), or the exposed group of the study-level exposure ('5 Gy
    Irradiated' in a 'Radiation exposure' study). Controls, vehicle, sham, ambient air,
    zero doses and pre-exposure samples get nothing; so do QC and cultured samples and
    samples with two different exposures."""
    if EXPOSURE_TABLE is None or len(EXPOSURE_TABLE) == 0 or len(df) == 0 or _ANALYSIS_IN_VITRO:
        return {}
    cindex = _condition_index(doid_table, path_to_csvs)
    rows = df[['filename', 'Key', 'Value']].dropna(subset=['filename']).copy()
    rows['Key'] = rows['Key'].astype(str).str.strip().str.lower()
    rows['Value'] = rows['Value'].astype(str).str.strip()
    rows = rows.drop_duplicates()
    if len(rows) == 0:
        return {}
    words = study_exposure['words'] if study_exposure else set()
    cache = {}

    def _classify(key, value):
        # ('pos', terms) / ('neg',) / ('pre',) / ('qc',) / ('culture',) / None
        kn = _norm_key(key)
        vn = _norm_value(value)
        if not vn or vn in _DIS_NA_VALUES:
            return None
        # ('water' is a vehicle here, not a blank)
        if qc_sample_type(key, value):
            return ('qc',)
        if _DIS_KEY_TISSUE_RE.search(kn):
            return ('culture',) if _DIS_CULTURE_RE.search(vn) or _EXP_CELL_LINE_RE.search(vn) else None
        strong = bool(_EXP_KEY_STRONG_RE.search(kn))
        key_word = any(re.search(r'(?<![a-z])' + re.escape(w), kn) for w in words)
        if _EXP_PRE_RE.fullmatch(vn) and (strong or key_word or _EXP_TIME_KEY_RE.search(kn) or
                                          _EXP_KEY_GENERIC_RE.fullmatch(kn) or _EXP_STATE_KEY_RE.search(kn)):
            return ('pre',)
        # the first visit / time point of a longitudinal design may be before the exposure
        if _EXP_FIRST_VISIT_RE.fullmatch(vn) and len(kn) <= 25 and \
                (_EXP_TIME_KEY_RE.search(kn) or _EXP_STATE_KEY_RE.search(kn)) and \
                not re.search(r'\b(?:censor\w*|survival|os|pfs|death|age|score|grade|stage)\b', kn):
            return ('pre',)
        if not (strong or key_word or _EXP_KEY_GENERIC_RE.fullmatch(kn)) or _HPO_KEY_SKIP_RE.search(kn) or \
                (_EXP_KEY_SKIP_RE.search(kn) and not key_word and not re.search(r'expos|irradiat|treat', kn)):
            return None
        name, nums = _exposure_value(value)
        if _DIS_MOCK_RE.search(vn) or vn in _EXP_NEG_VALUES or name in _EXP_NEG_VALUES or \
                any(n == 0 for n in nums) or re.match(r'(?:non|un|no|not|without|anti)\b', vn) or \
                (key_word and not _EXP_TIME_KEY_RE.search(kn) and re.fullmatch(r'0+(?:\.0+)?', vn.replace(' ', ''))):
            return ('neg', key_word) if (strong or key_word) else None
        terms = set()
        if name and name not in _EXP_VALUE_SKIP and not re.search(r'[+&/]| and | with ', str(value).lower()):
            t, tp = _exposure_term(name, cindex, False), _exposure_term(name, cindex, True)
            if t and (strong or not t.split('|')[1].startswith('ECTO:3')) and \
                    not (t.split('|')[1].startswith('ECTO:3') and not strong):
                terms.add((t, tp or t))
        if study_exposure:
            pos = name in words
            if study_exposure['radiation']:
                gy = [float(m.group(1)) for m in _EXP_GY_RE.finditer(str(value).lower())]
                pos = pos or (bool(gy) and all(g > 0 for g in gy))
                spec = _EXP_RADIATION_SPECIFIC.get(name)
                if spec and (key_word or strong):
                    st = _exposure_term(_norm_doid_name(spec), cindex, False)
                    if st:
                        terms.add((st, _exposure_term(_norm_doid_name(spec), cindex, True) or st))
                    pos = True
            if key_word and not _EXP_TIME_KEY_RE.search(kn):
                pos = pos or vn in _EXP_POS_WORDS or \
                    bool(re.fullmatch(r'\d+(?:\.\d+)?', vn.replace(' ', '.', 1)) and float(vn.replace(' ', '.', 1)) > 0)
            if pos:
                terms.add((study_exposure['term'], study_exposure['plant_term']))
        return ('pos', terms) if terms else None

    per_sample = {}
    for fn, g in rows.groupby('filename'):
        pos, pos_keys, neg_keys, skip = set(), set(), set(), False
        for key, val in g[['Key', 'Value']].itertuples(index=False):
            if (key, val) not in cache:
                cache[(key, val)] = _classify(key, val)
            c = cache[(key, val)]
            if c is None:
                continue
            if c[0] in ('qc', 'culture', 'pre'):
                skip = True
            elif c[0] == 'neg':
                # a control value cancels the exposure of the same factor, or any exposure
                # where it is about the study-level exposure ('radiation source: sham')
                neg_keys.add('*' if c[1] else key)
            else:
                pos |= c[1]
                pos_keys.add(key)
        if skip or not pos or '*' in neg_keys or neg_keys & pos_keys:
            continue
        if len(pos) > 1 and study_exposure:
            # the study-level term next to a narrower one ('exposure to gamma radiation')
            general = (study_exposure['term'], study_exposure['plant_term'])
            narrower = {t for t in pos if t != general and
                        any(w in _norm_doid_name(t[0].split('|')[0]) for w in words)}
            if general in pos and len(narrower) == 1:
                pos = narrower
        if len(pos) == 1:
            per_sample[fn] = next(iter(pos))
    return per_sample


# ---------------------------------------------------------------------------------
# Study-defined abbreviations ('acute kidney injury (AKI)'): within the study whose own
# text defines them, sample values naming the abbreviation ('AKI', 'AKI_d3', 'sham vs
# AKI') are read as the full name by the disease, health-condition and exposure steps,
# under all their rules, and only where the unexpanded values gave nothing. An
# abbreviation is defined if the study text has '<full name> (<ABBR>)' or '<ABBR> (<full
# name>)' with the full name an exact DOID / HPO / ECTO (PECO for plants) name and the
# abbreviation built from its initials (or an exact synonym of the term), or if it is an
# exact synonym of a term whose name appears in the study text. Standard disease-model
# procedures ('cecal ligation and puncture', CLP) mean their condition when the study
# text names the procedure (MWB_procedure_aliases.csv). Outside the defining study
# nothing changes.
# ---------------------------------------------------------------------------------
# REST title/summary of the study being converted (MWB_to_REDU_study_wrapper)
_STUDY_TEXT_REST = ''
# per analysis (create_dataframe_outer_dict): [(regex, full name)]
_ANALYSIS_ABBR = []
_ABBR_INDEX_CACHE = {}
# lab / method abbreviations and vehicles that must never become a condition or exposure
_ABBR_STOP = {'QC', 'QCS', 'WT', 'KO', 'KI', 'MS', 'LC', 'GC', 'NMR', 'HPLC', 'UPLC', 'UHPLC', 'ESI', 'ID', 'IDS', 'PBS',
              'DMSO', 'ETOH', 'MEOH', 'ACN', 'H2O', 'BSA', 'FBS', 'HFD', 'ND', 'CD', 'NC', 'SD', 'NA', 'NS', 'IP', 'IV',
              'PO', 'SC', 'IM', 'BW', 'BMI', 'RT', 'PCR', 'QPCR', 'RNA', 'DNA', 'CON', 'CTRL', 'CTL', 'VEH', 'ALL',
              'AND', 'FOR', 'THE', 'NOT', 'MIX', 'PRE', 'POST', 'LOW', 'HIGH', 'MAX', 'MIN', 'AGE', 'SEX', 'DAY',
              'USA', 'UK', 'EU', 'NIH', 'NIST', 'SRM', 'CSF', 'BAL', 'GI', 'IS', 'ISTD', 'OR', 'US', 'TIC', 'MRM',
              'SRM', 'PRM', 'DDA', 'DIA', 'API', 'HILIC', 'RP', 'NP', 'TOF', 'QTOF', 'FT', 'IT'}
# English words that are also abbreviations: factor values are compared in lower case
_ABBR_WORDS = {'ART', 'PER', 'GEN', 'HEX', 'MARCH', 'CANDLE', 'CAN', 'MAN', 'MEN', 'FAT', 'DOG', 'CAT', 'RAT', 'BAT', 'PIG',
               'COW', 'HEN', 'EGG', 'SUN', 'END', 'TOP', 'RED', 'HOT', 'ICE', 'AIR', 'SEA', 'OIL', 'GAS', 'TEA', 'JAM',
               'BAD', 'MAD', 'SAD', 'MID', 'LAB', 'MAP', 'PET', 'CUT', 'CAP', 'LID', 'TIP', 'SET', 'RUN', 'LAG', 'GAP',
               'BED', 'AID', 'ACE', 'APE', 'ARM', 'EAR', 'EYE', 'LEG', 'LIP', 'HIP', 'RIB', 'TOE', 'JAW', 'GUM', 'FUR',
               'HAY', 'LAW', 'WAR', 'ERA', 'ZOO', 'HAM', 'OAT', 'RYE', 'PEA', 'BEAN', 'CORN', 'RICE', 'MILK', 'MEAT',
               'FISH', 'BEEF', 'PORK', 'SALT', 'ACID', 'BASE', 'CORE', 'MASS', 'CELL', 'HOST', 'LINE', 'SITE', 'TYPE',
               'MODE', 'DOSE', 'ZERO', 'NONE', 'MALE', 'POOL', 'MOCK', 'SHAM', 'NAIVE', 'BLANK', 'MEDIA', 'STAGE',
               'PHASE', 'STOP', 'ARE', 'WAS', 'HAS', 'HAD', 'BUT', 'YES', 'ONE', 'TWO', 'SIX', 'TEN', 'OFF', 'OUT',
               'NEW', 'OLD', 'BIG', 'WET', 'DRY', 'RAW', 'MIXED', 'LEAN', 'OBESE', 'CASE', 'CASES', 'FED', 'FAST',
               'WEEK', 'YEAR', 'HOUR', 'LATE', 'EARLY', 'BOTH', 'LEFT', 'RIGHT', 'FREE', 'HALF', 'FULL', 'NEAR', 'FAR'}
_ABBR_CONTACT_KEY_RE = re.compile(r'(?:INSTITUTE|NAME|ADDRESS|EMAIL|PHONE|DEPARTMENT|LABORATORY|LAB_NAME|FUNDING|'
                                  r'SUBMIT_DATE|DATE|NUM_GROUPS|TOTAL_SUBJECTS|PUBLICATIONS|DOI|ID$|_ID$)')


def _abbr_like(s):
    """2-8 characters, mostly upper case or digits ('AKI', 'T2D', 'PM2.5')."""
    if not isinstance(s, str) or not (2 <= len(s) <= 8) or \
            not re.fullmatch(r'[A-Za-z0-9](?:[A-Za-z0-9.\-+]*[A-Za-z0-9+])?', s):
        return False
    core = re.sub(r'[^A-Za-z0-9]', '', s)
    if not re.search(r'[A-Za-z]', core) or s.upper() in _ABBR_STOP or s.upper() in _ABBR_WORDS:
        return False
    return sum(c.isupper() or c.isdigit() for c in core) >= 0.6 * len(core)


def _abbr_fits(abbr, full):
    """Schwartz-Hearst: the abbreviation's characters appear in order in the full name,
    its first character starting the full name's first word."""
    a = re.sub(r'[^a-z0-9]', '', abbr.lower())
    f = full.lower()
    if not a or not f or a[0] != f[0]:
        return False
    ai, fi = len(a) - 1, len(f) - 1
    while ai >= 0:
        c = a[ai]
        while fi >= 0 and (f[fi] != c or (ai == 0 and fi > 0 and f[fi - 1].isalnum())):
            fi -= 1
        if fi < 0:
            return False
        fi -= 1
        ai -= 1
    return True


def _abbr_index(doid_table, path_to_csvs, plant):
    """Long names -> terms and terms -> abbreviation synonyms of DOID, HPO and exposures."""
    key = (id(doid_table), id(HPO_TABLE), id(EXPOSURE_TABLE), str(path_to_csvs), plant)
    if key in _ABBR_INDEX_CACHE:
        return _ABBR_INDEX_CACHE[key]
    cindex = _condition_index(doid_table, path_to_csvs)
    dindex = cindex['doid']
    names, subst, term_names, term_abbrs = {}, {}, {}, {}
    for n, ids in dindex['norm_names'].items():
        ids = ids - {'DOID:4'}
        if len(n) > 3 and n not in _DISEASE_NAME_SKIP and len(ids) == 1:
            t = ('DOID', next(iter(ids)))
            names.setdefault(n, set()).add(t)
            term_names.setdefault(t, set()).add(n)
    for t in list(term_names):
        if t[1] in dindex['labels']:
            subst[t] = _norm_doid_name(dindex['labels'][t[1]])
    for n, terms in cindex['hpo_names'].items():
        ids = {i for _, i in terms}
        if len(ids) == 1:
            t = ('HPO', next(iter(ids)))
            names.setdefault(n, set()).add(t)
            term_names.setdefault(t, set()).add(n)
            subst[t] = _norm_doid_name(cindex['hpo_labels'][t[1]])
    for n in cindex['exp_names']:
        bare = re.sub(r'^exposure to |\s*exposure$', '', n).strip()
        if len(bare) <= 3 or bare in _EXP_VALUE_SKIP:
            continue
        term = _exposure_term(bare, cindex, plant)
        if term:
            t = ('EXP', term)
            names.setdefault(bare, set()).add(t)
            term_names.setdefault(t, set()).add(bare)
            subst.setdefault(t, bare)
    # exact synonyms that are abbreviations, in their own spelling
    for table, id_col, onto in [(doid_table, 'DOIDOntologyIndex', 'DOID'), (HPO_TABLE, 'HPOOntologyIndex', 'HPO')]:
        if table is None or len(table) == 0:
            continue
        for idx, syn in table[[id_col, 'Synonym']].itertuples(index=False):
            if _abbr_like(syn) and (onto, idx) in subst:
                term_abbrs.setdefault((onto, idx), set()).add(syn)
    procedures = []
    try:
        al = pd.read_csv(os.path.join(path_to_csvs, 'MWB_procedure_aliases.csv'), dtype=str, keep_default_na=False)
        for proc, abbr, redu in al[['Procedure', 'Abbreviation', 'REDU']].itertuples(index=False):
            procedures.append((_norm_doid_name(proc), abbr.strip(), _norm_doid_name(redu)))
    except Exception as e:
        print(f"Could not read the procedure aliases: {type(e).__name__}: {e}")
    # an abbreviation is no evidence of its own meaning: 'MED1' (a gene) names DOID's
    # 'multiple epiphyseal dysplasia 1' only through that very synonym
    abbr_norms = {_norm_doid_name(a) for abbrs in term_abbrs.values() for a in abbrs}
    index = {'names': names, 'subst': subst, 'term_names': term_names, 'term_abbrs': term_abbrs,
             'procedures': procedures, 'abbr_norms': abbr_norms}
    _ABBR_INDEX_CACHE[key] = index
    return index


def study_text(mwtab):
    """Text of the study's REST title/summary and its mwTab PROJECT, STUDY, SUBJECT,
    TREATMENT and COLLECTION sections (no contact fields)."""
    parts = [_STUDY_TEXT_REST]
    for sec in ['PROJECT', 'STUDY', 'SUBJECT', 'TREATMENT', 'COLLECTION']:
        s = mwtab.get(sec) if isinstance(mwtab, dict) else None
        if isinstance(s, dict):
            parts += [str(v) for k, v in s.items() if isinstance(v, (str, int, float)) and not _ABBR_CONTACT_KEY_RE.search(str(k))]
    return ' \n '.join(p for p in parts if p)


def study_abbreviations(text, doid_table, path_to_csvs='translation_sheets', plant=False):
    """[(regex for the abbreviation, full name)] defined by the study text."""
    if not text or not text.strip():
        return []
    index = _abbr_index(doid_table, path_to_csvs, plant)
    names, found = index['names'], {}

    def _add(abbr, full):
        found.setdefault(abbr, set()).add(full)

    # (a) 'full name (ABBR)' and 'ABBR (full name)'
    for m in re.finditer(r'\(\s*([A-Za-z0-9][A-Za-z0-9.\-+]{0,7})\s*[),;]', text):
        abbr = m.group(1).rstrip('.')
        if not _abbr_like(abbr):
            continue
        words = re.findall(r'[A-Za-z0-9\'\-]+', text[max(0, m.start() - 200):m.start()])[-(len(abbr) + 5):]
        for i in range(len(words)):
            full = ' '.join(words[i:])
            n = _norm_doid_name(full)
            if n in names and n not in index['abbr_norms'] and (_abbr_fits(abbr, n) or
                               any(abbr.lower() == _norm_doid_name(x) for t in names[n] for x in index['term_abbrs'].get(t, ()))):
                _add(abbr, n)
                break
    for m in re.finditer(r'(?<![A-Za-z0-9])([A-Za-z0-9][A-Za-z0-9.\-+]{1,7})\s*\(\s*([^()]{4,120}?)\s*\)', text):
        abbr, n = m.group(1), _norm_doid_name(m.group(2))
        if _abbr_like(abbr) and n in names and n not in index['abbr_norms'] and _abbr_fits(abbr, n):
            _add(abbr, n)
    # (b) abbreviation synonyms of terms named in the text
    words = _norm_doid_name(text).split()
    present = set()
    for i in range(len(words)):
        for j in range(i + 1, min(i + 9, len(words) + 1)):
            p = ' '.join(words[i:j])
            if p in names and p not in index['abbr_norms']:
                present |= names[p]
    for t in present:
        for abbr in index['term_abbrs'].get(t, ()):
            _add(abbr, index['subst'][t])
    # disease-model procedures named in the text
    ntext = ' ' + ' '.join(words) + ' '
    for proc, abbr, target in index['procedures']:
        if f' {proc} ' in ntext:
            _add(abbr, target)
            _add(proc, target)
    out = []
    for abbr, fulls in found.items():
        if len(fulls) != 1:
            print(f"Abbreviation {abbr} has several meanings in this study: {sorted(fulls)}")
            continue
        full = next(iter(fulls))
        if re.search(r'(?<![A-Za-z0-9])' + re.escape(abbr) + r'(?![A-Za-z0-9])(?:[\s\-/]+\S+){0,2}?[\s\-/]*'
                     r'(?:genes?|knock ?outs?|KO|-/-|\+/-|deficien\w*|mutants?|mutations?|mutated|variants?|alleles?|'
                     r'expression|expressing|overexpress\w*|mRNA|proteins?|levels?|knock ?down|siRNA|shRNA|'
                     r'transgenic|null|floxed|fl/fl|cohort|study|trial|consortium)\b', text):
            # a gene ('FADS2 levels', 'TSC mutation') or a study name ('MARCH cohort')
            print(f"Abbreviation {abbr} looks like a gene or study name here; not used")
            continue
        if abbr.lower() == abbr and ' ' in abbr:
            # a procedure's full name: any case, whole words
            rx = re.compile(r'(?<![a-z0-9])' + r'[\W_]+'.join(re.escape(w) for w in abbr.split()) + r'(?![a-z0-9])', re.I)
        else:
            # (values arrive in lower case) a separate token: 'AKI', 'AKI_d3', 'sham vs AKI';
            # not part of 'SARS-CoV-2'. Two-character abbreviations ('AD', 'PD') only as the
            # whole value: as tokens they are words ('ad libitum', 'as')
            if len(re.sub(r'[^A-Za-z0-9]', '', abbr)) <= 2:
                rx = re.compile(r'\s*' + re.escape(abbr) + r'\s*', re.I)
                rx = re.compile(r'^' + rx.pattern + r'$', re.I)
            else:
                rx = re.compile(r'(?<![A-Za-z0-9])(?<![A-Za-z]-)' + re.escape(abbr) + r'(?![A-Za-z0-9]|-[A-Za-z])', re.I)
        out.append((rx, full, abbr))
    if out:
        print(f"Study abbreviations: {sorted((a, f) for _, f, a in out)}")
    return [(rx, full) for rx, full, _ in sorted(out, key=lambda t: -len(t[2]))]


def expand_abbreviations(value, abbreviations):
    """The value with the study's abbreviations written out."""
    if not abbreviations or not isinstance(value, str):
        return value
    out = value
    for rx, full in abbreviations:
        out = rx.sub(' ' + full + ' ', out)
    return re.sub(r'\s+', ' ', out).strip() if out != value else value


_ABBR_ISOLATED_RE = re.compile(r'\b(?:\w+cytes?|\w+blasts?|cells?|cell lines?|chondrocytes?)\b')


def isolated_material_samples(df):
    """Samples whose tissue / sample-source value names isolated or cultured cells
    ('Human Chondrocytes', 'blood Tfh cells'): abbreviations give them nothing."""
    keys = df['Key'].astype(str).map(lambda k: bool(_DIS_KEY_TISSUE_RE.search(_norm_key(k))))
    vals = df['Value'].astype(str).map(_norm_value)
    iso = keys & vals.map(lambda v: bool(_DIS_CULTURE_RE.search(v) or _EXP_CELL_LINE_RE.search(v) or
                                         _ABBR_ISOLATED_RE.search(v)))
    return set(df.loc[iso, 'filename'])


def _merge_health(orig, alt):
    """Per-sample (health, term): the unexpanded result, completed by the expanded one."""
    out = dict(orig)
    for f, r in alt.items():
        if f not in out:
            out[f] = r
        elif out[f][1] is None and r[1] is not None and r[0] == out[f][0]:
            out[f] = r
    return out


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


# age units: units per year (the divisors used before for days, weeks, months)
_AGE_UNITS_PER_YEAR = {'mo': 12, 'w': 52.1429, 'd': 365, 'h': 8760}
_AGE_UNIT_WORDS = {'year': 'y', 'years': 'y', 'yr': 'y', 'yrs': 'y', 'y': 'y', 'yo': 'y', 'y/o': 'y',
                   'month': 'mo', 'months': 'mo', 'mo': 'mo', 'mos': 'mo', 'mth': 'mo', 'mths': 'mo', 'm': 'mo',
                   'week': 'w', 'weeks': 'w', 'wk': 'w', 'wks': 'w', 'w': 'w',
                   'day': 'd', 'days': 'd', 'd': 'd',
                   'hour': 'h', 'hours': 'h', 'hr': 'h', 'hrs': 'h', 'h': 'h'}
# unit words anywhere in a key or value (single letters only right after a number, see below)
_AGE_UNIT_WORD_RE = re.compile(r'(?<![a-z])(years?|yrs?|months?|mos|mths?|weeks?|wks?|days?|hours?|hrs?)(?![a-z])')
# a unit right after the number: '8w', '20M', '21 d', '7.5-mo', '9 h'
_AGE_NUM_UNIT_RE = re.compile(r'\d\s*-?\s*(years?|yrs?|y/o|yo|y|months?|mos?|mths?|m|weeks?|wks?|w|days?|d|hours?|hrs?|h)(?![a-z])')
# a key's unit in brackets: 'age (d)', 'age [w]'
_AGE_KEY_LETTER_RE = re.compile(r'[(\[]\s*(y|m|mo|w|wk|d|h)\s*[)\]]')
# the number of an age: not part of a word or code ('P21', 'E18.5', 'C57BL/6', 'SU001')
_AGE_NUMBER_RE = re.compile(r'(?<![a-z0-9_./])(\d+(?:\.\d+)?)(?![0-9]|[./][0-9]|[a-z]*[0-9])')
# 'less than 9 h', '>65', '~3 months': not one age
_AGE_BOUND_RE = re.compile(r'[<>≤≥~]|(?<![a-z])(?:over|under|above|below|less than|more than|greater than|older than|'
                           r'younger than|up to|at least|approx\w*|about|around|ca)(?![a-z])')

def _age_units(text, after_number):
    """Distinct units a text names (after_number: also single letters right after a number)."""
    t = str(text).lower()
    units = {_AGE_UNIT_WORDS[u] for u in _AGE_UNIT_WORD_RE.findall(t)}
    if after_number:
        units |= {_AGE_UNIT_WORDS.get(u, _AGE_UNIT_WORDS.get(u.rstrip('s'))) for u in _AGE_NUM_UNIT_RE.findall(t)}
    else:
        units |= {_AGE_UNIT_WORDS[u] for u in _AGE_KEY_LETTER_RE.findall(t)}
    return units - {None}


def study_age_unit(age_range):
    """The one unit the study-level SUBJECT AGE_OR_AGE_RANGE gives its ages in ('8-12 week
    old' -> 'w'), else None."""
    units = _age_units(age_range or '', after_number=True)
    return next(iter(units)) if len(units) == 1 else None


def mwb_age_to_years(value, key, default_unit=None):
    """Convert an age Factor value to years. The unit comes from the KEY ('age (months)',
    'age in days', 'age (d)'), else from the VALUE ('8 weeks', '20M', '9 h'), else
    default_unit (the unit of the study's AGE_OR_AGE_RANGE, or years for human studies);
    without one the age is unclear and None. Bounds ('< 9 h', '>65'), values with other
    words ('C57BL/6', 'P21') and keys or values naming two units are not an age. Decimals
    are kept ('16.8'); whole numbers of years stay integers."""
    value = '' if value is None else str(value)
    v = value.strip().lower()
    nums = _AGE_NUMBER_RE.findall(v)
    if not nums or _AGE_BOUND_RE.search(v):
        return None
    if re.search(r'[a-z][a-z0-9]*\d|\d+[/.]\d+[/.]', v):
        # a code with digits ('C57BL/6', 'P21', 'D15') or a date
        return None
    key_units = _age_units(str(key).replace('_', ' '), after_number=False)
    value_units = _age_units(v, after_number=True)
    units = key_units or value_units
    if len(units) > 1:
        return None
    unit = next(iter(units)) if units else default_unit
    if unit is None:
        return None
    num = float(nums[0])
    if unit == 'y':
        return int(num) if num.is_integer() else num
    if num.is_integer():
        num = int(num)
    return num / _AGE_UNITS_PER_YEAR[unit]


# unit for ages without one, set per analysis in create_dataframe_outer_dict: the unit of
# the study's AGE_OR_AGE_RANGE, else years for human studies, else None (unclear)
_ANALYSIS_AGE_UNIT = 'y'


def _is_human_study(entry, taxonomy):
    """The analysis' subjects are humans: taxonomy 9606, TAXONOMY_ID 9606, or a human
    SUBJECT_TYPE / SUBJECT_SPECIES."""
    if str(taxonomy or '').startswith('9606|') or str(entry.get('TAXONOMY_ID') or '').strip() == '9606':
        return True
    species = str(entry.get('SUBJECT_SPECIES') or '').strip().lower()
    return species in ('homo sapiens', 'human', 'humans') or \
        (str(entry.get('SUBJECT_TYPE') or '').strip().lower() in _HUMAN_SUBJECT_TYPES and species in ('', '-', 'na', 'n/a', 'none'))


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
    # a unit given in a key of its own ('age_unit: month'), else the analysis' default
    is_unit = keyl.str.fullmatch(r'(?:age[ _]?units?|units? of age|age[ _]unit[ _]of[ _]measure)', na=False)
    file_unit = {}
    for f, v in zip(df.loc[is_unit, 'filename'], df.loc[is_unit, 'Value']):
        units = _age_units(' ' + str(v), after_number=False)
        file_unit.setdefault(f, set()).update(units if len(units) == 1 else {None})
    file_unit = {f: next(iter(u)) for f, u in file_unit.items() if len(u) == 1 and None not in u}
    tmp['_yr'] = [mwb_age_to_years(v, k, file_unit.get(f, _ANALYSIS_AGE_UNIT))
                  for v, k, f in zip(tmp['Value'], tmp['Key'], tmp['filename'])]
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


# Excel stores dates as days since 1899-12-30; '41732' in an ACQUISITION_DATE is 2014-04-03
_EXCEL_SERIAL_RE = re.compile(r'\s*(\d{5})(?:\.\d+)?\s*')
# M/D/YY dates with a two-digit year ('1/14/21', '4/15/19-4/19/19')
_SHORT_YEAR_DATE_RE = re.compile(r'(?<![\d/.-])\d{1,2}[-/.]\d{1,2}[-/.](\d{2})(?![\d/.-])')


def analysis_year(dates, latest=None):
    """YearOfAnalysis from an acquisition date as written in the mwTab: 4-digit years in
    any date format, Excel serial dates ('41732') and M/D/YY dates. A year after
    `latest` (the study's submission year; data are acquired before they are submitted)
    is no year."""
    if not isinstance(dates, str) or not dates.strip():
        return None
    year = None
    m = _EXCEL_SERIAL_RE.fullmatch(dates)
    if m and 32874 <= int(m.group(1)) <= 47848:  # 1990-01-01 .. 2030-12-31
        import datetime as _dt
        year = (_dt.date(1899, 12, 30) + _dt.timedelta(days=int(m.group(1)))).year
    if year is None:
        year = extract_years(dates)
    if year is None:
        # also ranges of such dates ('4/15/19-4/19/19')
        parts = re.split(r'(?<=\d)\s*-\s*(?=\d{1,2}/)', dates)
        short = [2000 + int(y) for part in parts for y in _SHORT_YEAR_DATE_RE.findall(part)]
        short = [y for y in short if 1998 <= y <= 2030]
        if short:
            year = max(short)
    if year is not None and latest is not None and year > latest:
        return None
    return year


# submission year of the study being converted (MWB_to_REDU_study_wrapper): the latest
# YearOfAnalysis allowed
_STUDY_SUBMIT_YEAR = None


def submission_year(*dates):
    """Year of the first date that has one ('2014-09-15'), else None."""
    for d in dates:
        m = re.search(r'(?<!\d)((?:19|20)\d\d)(?!\d)', str(d or ''))
        if m:
            return int(m.group(1))
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


# keys that name the subject itself (their values may be the Subject IDs) and keys whose
# value cannot differ between the samples of one subject
_SUBJECT_KEY_RE = re.compile(r'subject|patient|participant|individual|donor|animal|mouse|mice|\brat\b|person|volunteer|'
                             r'\bid\b|identifier|\bname\b|\bno\b|number|\bcode\b|replicate|\brep\b|\bpair\b|family|'
                             r'litter|cage|\bpig\b|\bdog\b|\bcow\b|horse|monkey|\bplant\b|\btree\b|\bfish\b|\bhost\b')
_SUBJECT_LEVEL_KEY_RE = re.compile(r'genotype|strain|\bsex\b|gender|species|organism|breed|cell line|mutant|mutation|'
                                   r'knock ?out|genetic|transgen|\brace\b|ethnic')
# Subject IDs that only name a sample matrix or an analysis method ('serum', 'liver',
# 'HILIC NEG', 'RP POS')
_MATRIX_OR_METHOD_SID_RE = re.compile(
    r'(?:(?:whole |human |mouse |rat )?(?:serum|plasma|blood|urine|feces|faeces|stool|liver|kidney|brain|heart|lungs?|muscle|'
    r'tissue|cells?|media|medium|culture|supernatant|saliva|csf|milk|breast milk|bile|sweat|cecum|caecum|colon|'
    r'intestine|skin|adipose|fat|spleen|pancreas|retina|tumou?r)s?)'
    r'|(?:(?:hilic|rp|rplc|c18|c8|lc|gc|pos|neg|positive|negative|esi|apci|ms|ion|mode|polarity|reverse|reversed|'
    r'phase|normal|plus|minus)(?: |$))+')
# study-level SUBJECT / COLLECTION texts of the analysis (create_dataframe_outer_dict):
# a Subject ID they name ('Il10-/-' in GENOTYPE_STRAIN, 'liver' as SAMPLE_TYPE) is a label
_ANALYSIS_SUBJECT_TEXTS = []


def _norm_label(text):
    return ' '.join(re.sub(r'[\s_\-.,;:()\[\]{}]+', ' ', str(text).lower()).split())


def _drop_label_subject_ids(df):
    """Subject IDs that are labels rather than identifiers become '': values naming a
    sample category ('Plasma Study Sample', 'Quality Control'), values equal to one of
    the sample's own metadata values (a group such as 'wild-type' or the sample source),
    and a study whose samples all carry one and the same Subject ID. Also group labels
    (_group_label_subject_ids)."""
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
    try:
        bad |= _group_label_subject_ids(df, low, bad, is_meta, key_l)
    except Exception as e:
        print(f"Could not check Subject IDs for group labels: {type(e).__name__}: {e}")
    return sid.where(~bad, '')


def _group_label_subject_ids(df, low, bad, is_meta, key_l):
    """Rows whose Subject ID is a group label, not an individual:
    - a Subject ID (with a letter) equal to a metadata value of any sample of the study,
      made only of such values ('KO-ND' with Genotype KO and Diet ND), or named in the
      study's subject / collection description ('Il10-/-' in GENOTYPE_STRAIN, 'liver');
    - a Subject ID whose samples differ in a value one subject cannot change (genotype,
      strain, sex, species, cell line);
    - Subject IDs that are the levels of one (non-numeric) factor: each level one
      Subject ID, each Subject ID one level, on average >= 2 samples per Subject ID."""
    sample = (df['_sample_id'] if '_sample_id' in df.columns else df['filename']).astype(str)
    ok = (low != '') & ~bad
    if not ok.any():
        return pd.Series(False, index=df.index)
    norm_sid = low.map(_norm_label)
    has_letter = low.str.contains(r'[a-z]', regex=True)
    meta = is_meta & ~key_l.isin(['raw_data']) & ~key_l.map(lambda k: bool(_SUBJECT_KEY_RE.search(k)))
    study_values = {_norm_label(v) for v in df.loc[meta, 'Value'].astype(str).unique()} - {''}
    texts = [' ' + _norm_label(t) + ' ' for t in _ANALYSIS_SUBJECT_TEXTS if t]
    labels, why = set(), {}
    for s_raw, s_norm in sorted(set(zip(low[ok & has_letter], norm_sid[ok & has_letter]))):
        tokens = s_norm.split()
        if _MATRIX_OR_METHOD_SID_RE.fullmatch(s_norm):
            labels.add(s_raw); why.setdefault('sample matrix / method', []).append(s_raw)
        elif s_norm in study_values:
            labels.add(s_raw); why.setdefault('value', []).append(s_raw)
        elif len(tokens) >= 2 and all(t in study_values and re.search(r'[a-z]', t) for t in tokens):
            labels.add(s_raw); why.setdefault('values', []).append(s_raw)
        elif len(s_norm) >= 3 and any(' ' + s_norm + ' ' in t for t in texts):
            labels.add(s_raw); why.setdefault('text', []).append(s_raw)
    rows = pd.DataFrame({'sample': sample, 'sid': low, 'key': key_l, 'value': df['Value'].astype(str).str.strip().str.lower(),
                         'meta': meta})[ok]
    sid_of = rows.drop_duplicates(['sample', 'sid']).groupby('sample')['sid'].agg(lambda s: s.iloc[0] if s.nunique() == 1 else None).dropna()
    kv = rows[rows['meta'] & ~rows['value'].isin(['', 'na', 'n/a', '-', 'none', 'nan', 'not applicable', 'unknown'])]
    kv = kv.drop_duplicates(['sample', 'key', 'value'])
    # one subject cannot have two genotypes, sexes, species...
    subj_level = kv[kv['key'].map(lambda k: bool(_SUBJECT_LEVEL_KEY_RE.search(k)))].copy()
    # samples with one value of the key only (a sample listed as 'M' and 'Male' is no conflict)
    subj_level = subj_level[subj_level.groupby(['sample', 'key'])['value'].transform('nunique') == 1]
    if len(subj_level):
        subj_level['sid'] = subj_level['sample'].map(sid_of)
        n_vals = subj_level.dropna(subset=['sid']).groupby(['sid', 'key'])['value'].nunique()
        conflict = set(n_vals[n_vals > 1].index.get_level_values(0)) - labels
        labels |= conflict
        if conflict:
            why['subject-level ' + '/'.join(sorted(set(n_vals[n_vals > 1].index.get_level_values(1))))] = sorted(conflict)
    # the Subject IDs are the levels of one factor
    sids = set(sid_of.values)
    if 2 <= len(sids) <= 50 and len(sid_of) >= 2 * len(sids) and not sids <= labels:
        one_val = kv.groupby(['sample', 'key'])['value'].agg(lambda s: s.iloc[0] if len(s) == 1 else None).dropna()
        per_key = {}
        for (smp, key), val in one_val.items():
            if smp in sid_of.index:
                per_key.setdefault(key, {})[smp] = val
        for key in sorted(per_key):
            level_of = per_key[key]
            levels = set(level_of.values())
            if len(level_of) != len(sid_of) or len(levels) != len(sids) or not all(re.search(r'[a-z]', v) for v in levels):
                continue
            sid_levels, level_sids = {}, {}
            for smp, lev in level_of.items():
                sid_levels.setdefault(sid_of[smp], set()).add(lev)
                level_sids.setdefault(lev, set()).add(sid_of[smp])
            if all(len(v) == 1 for v in sid_levels.values()) and all(len(v) == 1 for v in level_sids.values()):
                why['levels of ' + key] = sorted(sids - labels)
                labels |= sids
                break
        else:
            # or of two factors combined ('KO-ND', 'WT-WD' for Genotype x Diet)
            keys = [k for k in sorted(per_key) if len(per_key[k]) == len(sid_of)
                    and all(re.search(r'[a-z]', v) for v in set(per_key[k].values()))
                    and 2 <= len(set(per_key[k].values())) < len(sids)]
            for i, k1 in enumerate(keys[:40]):
                for k2 in keys[i + 1:40]:
                    combo = {smp: (per_key[k1][smp], per_key[k2][smp]) for smp in sid_of.index}
                    if len(set(combo.values())) != len(sids):
                        continue
                    sid_levels, level_sids = {}, {}
                    for smp, lev in combo.items():
                        sid_levels.setdefault(sid_of[smp], set()).add(lev)
                        level_sids.setdefault(lev, set()).add(sid_of[smp])
                    if all(len(v) == 1 for v in sid_levels.values()) and all(len(v) == 1 for v in level_sids.values()):
                        why[f'levels of {k1} x {k2}'] = sorted(sids - labels)
                        labels |= sids
                        break
                else:
                    continue
                break
    # what is left is one Subject ID for all samples ('Breast_Milk' next to the NIST and
    # water samples' labels)
    left = set(sid_of.values) - labels
    if len(left) == 1 and (sid_of == next(iter(left))).sum() > 1 and len(labels) > 0:
        labels |= left
        why['only one left'] = sorted(left)
    if why:
        print('Subject IDs that are group labels: ' + '; '.join(f"{k}: {v[:5]}{' ...' if len(v) > 5 else ''} ({len(v)})"
                                                               for k, v in sorted(why.items())))
    return low.isin(labels) & ok


# ---------------------------------------------------------------------------
# TimePointHoursRelativeToEvent, SampleGroupFreeText, FastingStatus
# ---------------------------------------------------------------------------
# hours per unit; a month is 730 h (8760 h / 12), a year 8760 h (365 days)
_TP_UNIT_HOURS = {'s': 1 / 3600, 'min': 1 / 60, 'h': 1, 'd': 24, 'w': 168, 'mo': 730, 'y': 8760}
_TP_UNIT_WORDS = {
    's': 's', 'sec': 's', 'secs': 's', 'second': 's', 'seconds': 's',
    'min': 'min', 'mins': 'min', 'minute': 'min', 'minutes': 'min',
    'h': 'h', 'hr': 'h', 'hrs': 'h', 'hour': 'h', 'hours': 'h', 'hpi': 'h',
    'd': 'd', 'day': 'd', 'days': 'd', 'dpi': 'd',
    'w': 'w', 'wk': 'w', 'wks': 'w', 'week': 'w', 'weeks': 'w',
    'mo': 'mo', 'month': 'mo', 'months': 'mo',
    'y': 'y', 'yr': 'y', 'year': 'y', 'years': 'y'}
_TP_WORD_NUMBERS = {'zero': 0, 'one': 1, 'two': 2, 'three': 3, 'four': 4, 'five': 5, 'six': 6, 'seven': 7, 'eight': 8,
                    'nine': 9, 'ten': 10, 'eleven': 11, 'twelve': 12, 'fourteen': 14, 'twenty': 20, 'thirty': 30}
# keys about the time of sampling (also words glued to 'time': 'timedraw', 'samplingtimepoint')
_TP_KEY_RE = re.compile(r'time|\b(?:tp|timept|day|days|week|weeks|wk|wks|hour|hours|hr|hrs|minute|minutes|min|mins|'
                        r'month|months|hpi|dpi|seconds|post|after|since)\b')
# keys giving how long a treatment / exposure / diet lasted until sampling (used only
# where no time-point key gives a value)
_TP_DURATION_KEY_RE = re.compile(r'\b(?:duration|period|on diet|of exposure|exposure time)\b')
# time keys that are not the sampling time relative to an event: ages and development,
# storage, extraction / incubation / LC run, clock times and dates, clinical history and
# outcome, frequencies, fasting durations, ...
_TP_KEY_EXCLUDE_RE = re.compile(
    r'\b(?:age|aged|ages|ga|ega|gestation\w*|postnatal|postconceptional|postmenstrual|embryonic|developmental|'
    r'development|adulthood|hpf|dpf|in vitro|div|storage|store|stored|fridge|freezer|frozen|freez\w*|thaw\w*|ice|'
    r'incubation|incubated|extraction|extracted|retention|gradient|run|runs|acquisition|injection order|doubling|'
    r'ischemia|ischemic|pmi|postmortem|mortem|death|died|survival|os|censor\w*|mortality|outcome|rebound|follow|'
    r'followed|until|onset|sick|illness|dialysis|art|ati|diagnosis|btwn|between|date|dates|clock|zt|zeitgeber|'
    r'zietgeber|circadian|year|years|season|per|past|frequency|drinks|exercise|stools?|well being|interval|cycle|'
    r'caseindex|mgit|positive|difference|stab|subculture|preculture|sleep|fasting|fast|fasted|starvation|starved|'
    r'lag|delay|processing|centrifug\w*|transport|shipping|stability|elapsed|half life|ventilation|icu|copies|tb|'
    r'bmi|weight|fev1|fvc|fev1fvc|glucose|calprotectin|defensin|of day|hour of|draw hour|caffeine|use|kinetic|raw|'
    r'file|files|fraction|intensity|vo2\w*|mg|ml|kg|g|ug|um|nm|mm|mol|number)\b|ageat|ageof|µ'
    r'|\b(?:time|times|days?|weeks?|months?|hours?) to\b|\bfrom\b.*\bto\b')
_TP_NUM = r'(\d+(?:\.\d+)?|' + '|'.join(_TP_WORD_NUMBERS) + r')'
_TP_UNIT = r'(' + '|'.join(sorted(_TP_UNIT_WORDS, key=len, reverse=True)) + r')'
# '24 h', '2hr', '30_min', '0-h', '24hrpost', '1 wk after FI', '2 h before surgery'
_TP_NUM_UNIT_RE = re.compile(r'([-+]?)\s*' + _TP_NUM + r'\s*[-_]?\s*' + _TP_UNIT + r'\.?'
                             r'(?:\s*[-_]?\s*(post|after|following|later|pi|pre|before)\b.{0,40})?')
# 'day 7', 'D7', 'week 4', 'month 3', 'day -1', 'day one'
_TP_UNIT_NUM_RE = re.compile(r'(day|days|d|week|wk|month|hour|hr|h)\s*[-_]?\s*([-+]?)\s*' + _TP_NUM)
_TP_BARE_NUM_RE = re.compile(r'([-+]?)\s*(\d+(?:\.\d+)?)')


def _tp_number(text):
    return float(_TP_WORD_NUMBERS[text]) if text in _TP_WORD_NUMBERS else float(text)


def timepoint_key_tier(key):
    """1 for keys about the time point of sampling, 2 for treatment / exposure durations,
    None for other keys."""
    nk = _norm_key(key)
    if not nk or _TP_KEY_EXCLUDE_RE.search(nk):
        return None
    if _TP_DURATION_KEY_RE.search(nk):
        return 2
    if _TP_KEY_RE.search(nk):
        return 1
    return None


def _timepoint_key_unit(key):
    """The one time unit a key states ('time (h)', 'days post infection'), else None.
    A key that only is a unit word that may also mean a clock hour ('hour') gives none."""
    nk = _norm_key(key)
    if nk in {'hour', 'month', 'months', 'year', 'years'}:  # also clock hours / calendar months
        return None
    units = {_TP_UNIT_WORDS[t] for t in nk.split() if t in _TP_UNIT_WORDS}
    return units.pop() if len(units) == 1 else None


def timepoint_hours(value, key):
    """Hours of the sampling time relative to the study's reference event, from one
    factor value with an explicit number and time unit (in the value, or a number
    under a key that states the unit). None for anything else (baseline, pre, T1,
    ranges, several values, clock times)."""
    if timepoint_key_tier(key) is None:
        return None
    v = ' '.join(str(value).strip().lower().replace('−', '-').split())
    v = ' '.join(re.sub(r'\(\s*([a-z]+)\s*\)$', r' \1', v).split())  # '144(hrs)'
    v = re.sub(r'^(?:post|after)\s+(?=\d)', '', v)  # 'post 4 hr'
    if not v or ':' in v or '/' in v or '>' in v or '<' in v or '~' in v:
        return None
    m = _TP_NUM_UNIT_RE.fullmatch(v)
    if m:
        sign, num, unit, qual = m.group(1), m.group(2), m.group(3), m.group(4)
        hours = _tp_number(num) * _TP_UNIT_HOURS[_TP_UNIT_WORDS[unit]]
        if sign == '-' or qual in ('pre', 'before'):
            if sign == '-' and qual in ('pre', 'before'):
                return None
            hours = -hours
        return round(hours, 4)
    m = _TP_UNIT_NUM_RE.fullmatch(v)
    if m:
        unit = _TP_UNIT_WORDS[m.group(1)]
        hours = _tp_number(m.group(3)) * _TP_UNIT_HOURS[unit]
        return round(-hours if m.group(2) == '-' else hours, 4)
    m = _TP_BARE_NUM_RE.fullmatch(v)
    if m:
        unit = _timepoint_key_unit(key)
        if unit is None:
            return None
        hours = float(m.group(2)) * _TP_UNIT_HOURS[unit]
        return round(-hours if m.group(1) == '-' else hours, 4)
    return None


_TP_AGE_VALUE_RE = re.compile(r'postnatal|\bage\b|\bold\b|\bembryonic\b|^p ?\d+d?$|^e ?\d+(?:\.\d+)?$')
_TP_AGE_TEXT_RE = re.compile(r'\bage[- ]?matched|\b(?:months?|years?|yrs?)[- ]old\b|\b(?:months?|years?) of age\b', re.I)


def sample_timepoint_hours(df, study_text=''):
    """Per-sample TimePointHoursRelativeToEvent (Series on df's rows): time-point keys
    first, durations only for samples without one; the sample's values must agree.
    A key is not used where its values give ages ('postnatal day 7', 'P21', 'E12.5') or,
    for months / years, where the study describes its subjects by age ('6-month-old',
    'age matched')."""
    keys = df['Key'].astype(str)
    tiers = keys.map(timepoint_key_tier)
    age_keys = {k for k, v in zip(keys[tiers.notna()], df.loc[tiers.notna(), 'Value'].astype(str))
                if _TP_AGE_VALUE_RE.search(' '.join(v.strip().lower().split()))}
    tiers = tiers.where(~keys.isin(age_keys))
    rows = tiers.notna()
    if not rows.any():
        return pd.Series(np.nan, index=df.index, dtype=object)
    hours = pd.Series([timepoint_hours(v, k) for k, v in zip(keys[rows], df.loc[rows, 'Value'])],
                      index=df.index[rows], dtype=object)
    if study_text and _TP_AGE_TEXT_RE.search(study_text):
        # months / years (multiples of 730 h) in a study about subjects of some age
        long_keys = set(keys[rows][hours.map(lambda h: h is not None and abs(h) >= 730 and abs(h) % 730 == 0)])
        hours = hours.where(~keys[rows].isin(long_keys))
    t = pd.DataFrame({'filename': df.loc[rows, 'filename'], 'tier': tiers[rows], 'h': hours}).dropna()
    if len(t) == 0:
        return pd.Series(np.nan, index=df.index, dtype=object)
    t = t.drop_duplicates()
    t = t[t['tier'] == t.groupby('filename')['tier'].transform('min')]
    agreed = t.groupby('filename')['h'].agg(lambda s: s.iloc[0] if s.nunique() == 1 else None).dropna()
    return df['filename'].map(agreed)


# ---- FastingStatus ----
_FAST_KEY_RE = re.compile(r'\b(?:fast|fasted|fasting|feeding|feeding state|feeding status|nutritional state|'
                          r'nutritional status|prandial state|prandial status|fed or fasted|fed vs fasted)\b')
_FAST_KEY_EXCLUDE_RE = re.compile(r'glucose|insulin|glyc|lipid|cholesterol|triglycer|\btg\b|hba1c|level|concentration|'
                                  r'breast|period|\bdha\b|\bmg\b|\bml\b|after|re ?feed|refed|\btime\b|\bduration\b|'
                                  r'\bweeks?\b|\bwk\b|\bdays?\b')
# keys that name the groups of the study design (fasting values are taken there only
# when the value is a fasting state and nothing else)
_FAST_GENERIC_KEY_RE = re.compile(r'^(?:groups?|factors?|conditions?|class|classes|categor(?:y|ies)|arms?|study arms?|'
                                  r'treatments?|treatment groups?|treamtent group|experimental factors?|experimental groups?|'
                                  r'experimental conditions?|interventions?|diets?|diet groups?|sampling|sample groups?|'
                                  r'time ?points?|time|state|status|nutrition|nutritional status)$')
_FASTED_PART_RE = re.compile(r'(?:(?:\d+(?:\.\d+)? ?(?:h|hr|hrs|hours?) ?)?(?:fasted|fasting|fast)'
                             r'(?: ?(?:for )?\d+(?:\.\d+)? ?(?:h|hr|hrs|hours?))?(?: overnight)?|'
                             r'overnight (?:fasted|fasting|fast)|morning fasting|fasting overnight)')
_NONFASTED_PART_RE = re.compile(r'(?:fed|non ?fasted|nonfasted|unfasted|not fasted|postprandial|post prandial|sated|satiated)')


def _fasting_text(value):
    v = str(value).strip().lower().replace('-', ' ').replace('_', ' ')
    v = re.sub(r'(\d)([a-z])', r'\1 \2', v)
    v = re.sub(r'([a-z])(\d)', r'\1 \2', v)
    return ' '.join(v.split())


def fasting_status(value, key):
    """'fasted' / 'non-fasted' from one per-sample value, else None."""
    nk = _norm_key(key)
    if not nk:
        return None
    v = _fasting_text(value)
    if not v:
        return None
    if _is_fasting_key(key):
        if re.search(r'\bintermittent\b|\bpre\b|\bpost fast|\bre ?fe[de]', v):
            return None
        if _FASTED_PART_RE.fullmatch(v) or v in {'overnight', 'over night', 'overnight fast'}:
            return 'fasted'
        if _NONFASTED_PART_RE.fullmatch(v) or v in {'ad lib', 'ad libitum', 'ad lib fed', 'fed ad libitum'}:
            return 'non-fasted'
        # yes / no only where the key asks whether the subject fasted
        if re.fullmatch(r'(?:subject )?(?:fast|fasted|fasting)(?: status| state)?', nk):
            if v in {'yes', 'y', 'true'}:
                return 'fasted'
            if v in {'no', 'n', 'false'}:
                return 'non-fasted'
        return None
    if _FAST_GENERIC_KEY_RE.match(nk):
        # the value, or one part of it ('Fasted/Saline', 'WT - 16h fast', 'LKO_Fast'), is
        # only a fasting state; parts with both states leave it open
        parts = [p.strip() for p in re.split(r'\s*(?:/|;|,|\||\+| and )\s*', str(value).strip().lower()) if p.strip()]
        parts = [q for p in parts for q in re.split(r'\s+-\s+|_', p) if q.strip()]
        found = set()
        for p in parts:
            t = _fasting_text(p)
            if _FASTED_PART_RE.fullmatch(t):
                found.add('fasted')
            elif _NONFASTED_PART_RE.fullmatch(t):
                found.add('non-fasted')
            elif re.search(r'fast|\bfed\b|feed|prandial', t):
                return None
        if len(found) == 1:
            whole = _fasting_text(value)
            if re.search(r'\bintermittent\b|\bpre\b|\bpost fast|\bre ?fe[de]', whole):
                return None
            return found.pop()
    return None


def fasting_hours_status(value, key):
    """Fasting duration in hours under a key that gives it ('fasting prior to blood
    draw (hours)'): 0 -> non-fasted, >= 8 h -> fasted, else None."""
    nk = _norm_key(key)
    if not re.search(r'\b(?:fast|fasted|fasting)\b', nk) or re.search(r'glucose|insulin|glyc|lipid', nk):
        return None
    if _timepoint_key_unit(key) != 'h':
        return None
    v = str(value).strip()
    if not re.fullmatch(r'\d+(?:\.\d+)?', v):
        return None
    h = float(v)
    return 'non-fasted' if h == 0 else 'fasted' if h >= 8 else None


def _is_fasting_key(key):
    nk = _norm_key(key)
    return bool(_FAST_KEY_RE.search(nk)) and not _FAST_KEY_EXCLUDE_RE.search(nk)


def sample_fasting_status(df):
    """Per-sample FastingStatus (Series on df's rows); the sample's values must agree.
    Under a key that is not about fasting ('Treatment', 'Condition'), 'fed' / 'sated'
    count only where the same key also names fasted samples ('Fed' alone can be a
    feeding treatment)."""
    keys = df['Key'].astype(str)
    vals = [fasting_status(v, k) or fasting_hours_status(v, k) for k, v in zip(keys, df['Value'])]
    s = pd.Series(vals, index=df.index, dtype=object)
    nkeys = keys.map(_norm_key)
    generic = ~keys.map(_is_fasting_key)
    fasted_keys = set(nkeys[s.eq('fasted')])
    s = s.where(~(generic & s.eq('non-fasted') & ~nkeys.isin(fasted_keys)))
    t = pd.DataFrame({'filename': df['filename'], 'f': s}).dropna()
    if len(t) == 0:
        return pd.Series(np.nan, index=df.index, dtype=object)
    agreed = t.groupby('filename')['f'].agg(lambda x: x.iloc[0] if x.nunique() == 1 else None).dropna()
    return df['filename'].map(agreed)


_STUDY_FAST_UNCLEAR_RE = re.compile(r'intermittent|ramadan|likely|maybe|possibl|probabl|some|subset|\bfor (?:the )?\w+ '
                                    r'(?:study|test|experiment|group)|first sample|inoculation|challenge|larva|do not feed|'
                                    r'\bif\b|unless|\bor\b|varied|various|depend|not required|no fasting required|'
                                    r'not controlled|not specified|unknown|breakfast|\bbut\b|incubator|dietary|diet\b')
_STUDY_FAST_EVENT_RE = re.compile(r'\b(?:before|prior to)\s+(?:the\s+)?(?:\w+\s+)?(\w+)')
_STUDY_FAST_EVENT_OK = re.compile(r'(?:sampl|collect|sacrific|euthan|end|blood|serum|plasma|necrops|dissect|harvest|'
                                  r'experiment|kill|surgery|tissue|draw|study|being studied|test|measure)')


def study_fasting_status(text):
    """FastingStatus from the study-level ANIMAL_FASTING / HUMAN_FASTING text when it
    states one fasting condition for all subjects ('Overnight', '12 h', 'non-fasted',
    ...); None for mixed or unclear statements. A bare 'no'/'none' is None too: in this
    field it often means 'no information' rather than 'not fasted'."""
    if not isinstance(text, str):
        return None
    v = ' '.join(text.strip().lower().replace(' ', ' ').split()).rstrip('.')
    if not v or v in {'na', 'n/a', 'n a', 'not applicable', 'not available', '-', 'nan', 'unknown', 'not specified'}:
        return None
    if '|' in v or _STUDY_FAST_UNCLEAR_RE.search(v):
        return None
    if v in {'no', 'none'}:
        return None
    if re.fullmatch(r'(?:no fasting|non[- ]?fasted|not fasted|unfasted|fed|ad lib(?:itum)?|ad libitum fed|'
                    r'fed ad lib(?:itum)?|free feeding)', v):
        return 'non-fasted'
    if re.search(r'\b(?:no|not|non|without|un)[- ]?fast', v):
        return None
    for m in _STUDY_FAST_EVENT_RE.finditer(v):
        if not _STUDY_FAST_EVENT_OK.match(m.group(1)):
            return None
    nums = [float(x) for x in re.findall(r'(\d+(?:\.\d+)?)\s*(?:-\s*\d+(?:\.\d+)?\s*)?(?:h|hr|hrs|hours?)\b', v)]
    words = {'six': 6, 'eight': 8, 'twelve': 12, 'four': 4, 'five': 5, 'ten': 10, 'sixteen': 16}
    nums += [words[w] for w in re.findall(r'\b(' + '|'.join(words) + r')\s*(?:h|hr|hrs|hours?)\b', v)]
    if re.search(r'\b(?:0|zero)\s*(?:-|to)', v):
        return None
    if nums:
        return 'fasted' if min(nums) >= 4 else None
    if re.search(r'\b(?:yes|fasted|fasting|fast|overnight|over night|one night)\b', v) or v in {'y'}:
        return 'fasted'
    return None


# ---- SampleGroupFreeText ----
# factor keys whose information has a ReDU column of its own (sex, age, life stage,
# organism, body part / sample type, disease / health, smoking, BMI, diet, location,
# dates and times, fasting, instrument / LC / polarity, methods, QC / blank class) or
# that are not groups (batch, replicate, run order, injection, plate / well, raw files,
# subject IDs)
_GROUP_KEY_EXCLUDE_RE = re.compile(
    r'time|visit|\b(?:day|days|week|weeks|hour|hours|hr|hrs|month|months|tp|timept)\b|'
    r'\b(?:batch|batches|replicate|replicates|rep|reps|replication|run|runs|order|injection|inj|plate|well|position|'
    r'vial|tray|rack|sequence|seq|block|polarity|ion mode|ionization|ionisation|column|chromatography|lc|ms|'
    r'instrument|platform|method|methods|extraction|protocol|date|raw|file|files|mzml|id|ids|identifier|identifiers|'
    r'sex|gender|age|aged|ages|life stage|lifestage|developmental stage|species|organism|organisms|taxon|taxonomy|'
    r'tissue|tissues|organ|organs|body part|bodypart|body site|sample type|sample types|sampletype|sample source|'
    r'source|matrix|specimen|specimens|biofluid|fluid|material|biospecimen|disease|diseases|diagnosis|diagnoses|'
    r'health|healthy|smok\w*|tobacco|bmi|body mass|diet|diets|dietary|chow|feed|feeding|fed|fasting|fasted|fast|'
    r'country|location|site|region|latitude|longitude|lat|lon|depth|altitude|elevation|collection|season|year|'
    r'qc|blank|blanks|quality|sample class|samle type|subject|subjects|cage|cages|moment|duration|ion model)\b|'
    r'^age(?!nt)|\bage(?:months?|years?|weeks?|days?|at\w*)|collection\w*|ionmode|tissue\w*|\w*duration\w*')
_GROUP_ID_KEY_RE = re.compile(r'(?:subject|participant|individual|donor|person)s?|'
                              r'(?:subject|participant|patient|animal|mouse|mice|rat|individual|donor|sample|pig|cow|'
                              r'fish|plant|tree|person|host)s? (?:no|number|num|nr|name|code|label)')
# design keys that are dropped where the study's ReDU columns already tell their groups
# apart (treatment -> Exposure, group -> DOIDCommonName, phenotype -> HPOCommonName, ...)
# (every factor except genotype / strain / cell line / cultivar-like ones, which the ReDU
# columns never describe)
_GROUP_INDEPENDENT_KEY_RE = re.compile(
    r'genotyp\w*|\b(?:strains?|lines?|cell lines?|mutants?|mutations?|knock ?outs?|ko|clones?|cultivars?|variet(?:y|ies)|'
    r'accessions?|alleles?|transgen\w*|construct|plasmid|sirna|shrna|crispr|wt)\b')
_GROUP_NA_VALUES = {'', '-', '--', 'na', 'n/a', 'nan', 'null',
                     # no group information ('Factor: unknown')
                     'unknown', 'n.a.', 'n.a', 'not applicable', 'not available', 'not specified', 'unspecified',
                     'not provided', 'not reported', 'missing', 'missing value'}
# internal form until the study-level check: sample ID \x1d key \x1e value \x1f key \x1e value
_GROUP_SID, _GROUP_KV, _GROUP_PART = '\x1d', '\x1e', '\x1f'
# ReDU columns that may already tell a design factor's groups apart
_GROUP_REDU_COLUMNS = ['SampleType', 'SampleTypeSub1', 'NCBITaxonomy', 'UBERONBodyPartName', 'BiologicalSex',
                       'HealthStatus', 'DOIDCommonName', 'HPOCommonName', 'Exposure', 'Diet', 'SmokingStatus',
                       'LifeStage', 'Country', 'ENVOEnvironmentMaterial', 'FastingStatus',
                       'TimePointHoursRelativeToEvent']


def sample_group_texts(data):
    """Sample ID -> the sample's grouping factors from the mwTab SUBJECT_SAMPLE_FACTORS
    'Factors', in the study's order and wording (internal form, see
    finalize_sample_groups). Factors that have a ReDU column of their own, that are no
    groups, and factors with a different value for nearly every sample are left out.
    ';' inside a key or value becomes ','."""
    items = [(str(it.get('Sample ID', '')).strip(), it.get('Factors') or {})
             for it in (data or []) if isinstance(it, dict)]
    items = [(sid, f) for sid, f in items if sid and isinstance(f, dict)]
    if not items:
        return {}
    values = {}
    for _, f in items:
        for k, v in f.items():
            values.setdefault(str(k).strip(), []).append(str(v).strip())
    n = len({sid for sid, _ in items})
    drop = set()
    for k, vs in values.items():
        nk = _norm_key(k)
        if not nk or _GROUP_KEY_EXCLUDE_RE.search(nk) or _GROUP_ID_KEY_RE.fullmatch(nk) \
                or is_variant_raw_file_key(k.lower()) or k.lower() in _EXACT_RAW_KEYS:
            drop.add(k)
            continue
        filled = [x for x in vs if x.lower() not in _GROUP_NA_VALUES]
        if n >= 8 and len(set(filled)) >= 0.9 * n:
            drop.add(k)
    clean = lambda x: ' '.join(str(x).replace(';', ',').replace(_GROUP_SID, ' ').replace(_GROUP_KV, ' ')
                               .replace(_GROUP_PART, ' ').split())
    out = {}
    for sid, f in items:
        parts = []
        for k, v in f.items():
            k, v = str(k).strip(), clean(v)
            if k in drop or v.lower() in _GROUP_NA_VALUES:
                continue
            parts.append(clean(k) + _GROUP_KV + v)
        text = _GROUP_PART.join(parts)
        if sid in out and out[sid] != text:
            out[sid] = None  # one Sample ID listed with different factors
        elif sid not in out:
            out[sid] = text or None
    return {sid: clean(sid) + _GROUP_SID + t for sid, t in out.items() if t}


def finalize_sample_groups(df, column='SampleGroupFreeText', missing='missing value', study_id=''):
    """Turns the internal sample groups into 'Key: value; Key: value' for one study.
    Factors (other than genotype / strain / cell line-like ones) whose groups the
    study's ReDU columns already tell apart are dropped; the column is filled only where the real
    samples (no QC / blanks) form at least 2 groups and at most half as many groups as
    samples. Never removes rows."""
    if column not in df.columns or len(df) == 0:
        return df
    df = df.copy()
    raw = df[column].astype(str)
    has = raw.str.contains(_GROUP_SID, regex=False)
    if not has.any():
        df[column] = missing
        return df
    sid = raw.where(has).str.split(_GROUP_SID, n=1).str[0]
    parts = raw.where(has).str.split(_GROUP_SID, n=1).str[1].fillna('').map(
        lambda t: [p.split(_GROUP_KV, 1) for p in t.split(_GROUP_PART) if _GROUP_KV in p])
    sample_type = df['SampleType'].astype(str) if 'SampleType' in df.columns else pd.Series('', index=df.index)
    real = has & ~(sample_type.isin(['pool_QC', 'reference material']) | sample_type.str.contains('blank', case=False))
    # design keys whose value is a function of the ReDU columns over the real samples
    redu_cols = [c for c in _GROUP_REDU_COLUMNS if c in df.columns]
    redu_tuple = df[redu_cols].astype(str).agg('\x1f'.join, axis=1) if redu_cols else pd.Series('', index=df.index)
    keys = list(dict.fromkeys(k for ps in parts[real] for k, _ in ps))
    dropped = set()
    for k in keys:
        if _GROUP_INDEPENDENT_KEY_RE.search(_norm_key(k)):
            continue
        kv = parts[real].map(lambda ps: next((v for kk, v in ps if kk == k), None))
        t = pd.DataFrame({'redu': redu_tuple[real], 'v': kv, 'sid': sid[real]}).dropna().drop_duplicates(['sid'])
        if len(t) and t['v'].nunique() > 1 and (t.groupby('redu')['v'].nunique() <= 1).all():
            dropped.add(k)
    text = parts.map(lambda ps: '; '.join(f'{k}: {v}' for k, v in ps if k not in dropped))
    text = text.where(has & (text != ''))
    groups = pd.DataFrame({'sid': sid[real], 'g': text[real]}).dropna().drop_duplicates('sid')
    n_samples, n_groups = len(groups), groups['g'].nunique()
    reason = None
    if n_groups < 2:
        reason = 'fewer than 2 groups'
    elif n_groups > n_samples / 2:
        reason = 'groups > half of the samples'
    if dropped:
        print(f"[SampleGroupFreeText] {study_id}: dropped as covered by ReDU columns: {sorted(dropped)}")
    if reason:
        print(f"[SampleGroupFreeText] {study_id}: not filled ({reason}; {n_groups} groups, {n_samples} samples)")
        df[column] = missing
        return df
    print(f"[SampleGroupFreeText] {study_id}: filled ({n_groups} groups, {n_samples} samples)")
    df[column] = text.fillna(missing)
    return df


def study_fasting_text(mwtab):
    """The study-level ANIMAL_FASTING / HUMAN_FASTING statements (TREATMENT, else
    SUBJECT section), joined by ' || '."""
    texts = []
    for section in ('TREATMENT', 'SUBJECT'):
        for key in ('ANIMAL_FASTING', 'HUMAN_FASTING'):
            value = _section(mwtab, section).get(key)
            if isinstance(value, str) and value.strip():
                texts.append(' '.join(value.split()))
    return ' || '.join(dict.fromkeys(texts))


def study_fasting_from_texts(joined):
    """One FastingStatus from all study-level fasting statements, or None if any of
    them is unclear or they disagree."""
    if not isinstance(joined, str) or not joined.strip():
        return None
    found = set()
    for text in joined.split(' || '):
        t = text.strip().lower()
        if t in {'na', 'n/a', 'n a', 'not applicable', 'not available', '-', 'nan'}:
            continue
        status = study_fasting_status(text)
        if status is None:
            return None
        found.add(status)
    return found.pop() if len(found) == 1 else None


# ---- CellLine (Cellosaurus) and cell types (Cell Ontology) ----
# Cell-line names under sample factor keys ('Cell line', 'cell_type', ...) and in the
# study-level CELL_STRAIN_DETAILS; matched to the shipped Cellosaurus table (exact
# name / exact synonym / normalised name, checked against the study organism). Values
# that are cell types ('macrophages', 'PBMC') go to UBERONBodyPartName as Cell Ontology
# terms instead.
# Cellosaurus table (DataFrame); set from --path_to_cellosaurus_tsv
CELLOSAURUS_TABLE = None
_CELLOSAURUS_INDEX_CACHE = {}
_CELL_KEY_RE = re.compile(r'cell[ _-]?(?:lines?|types?|strains?|names?|models?|source)\b|^cells?$|cellosaurus|^cell[ _]?line')
_CELL_KEY_EXCLUDE_RE = re.compile(r'count|number|passage|density|viab|conc|volume|seed|size|cycle|confluen|amount|#|'
                                  r'harvest|date|time|fraction|stage|media|medium')
_CELL_VALUE_SEP = '\x1f'
# Cell Ontology labels too generic for a body part
_CL_TOO_GENERIC = {'cell', 'native cell', 'eukaryotic cell', 'animal cell', 'somatic cell', 'cultured cell',
                   'cell in vitro', 'primary cultured cell', 'experimentally modified cell in vitro', 'cell line cell',
                   'nucleate cell', 'mononuclear cell', 'single nucleate cell', 'motile cell', 'non-terminally differentiated cell'}
_CL_PREFIX_RE = re.compile(r'^(?:(?:primary|human|mouse|murine|rat|cultured|isolated)\s+)+')
_CL_LOOKUP_CACHE = {}
# allowed Cell Ontology labels (set by _cl_lookup)
_CL_LABELS = set()


def _cellosaurus_index():
    from cellosaurus import CellosaurusIndex
    table = CELLOSAURUS_TABLE
    if table is None or len(table) == 0:
        return None
    key = (id(table), len(table))
    if key not in _CELLOSAURUS_INDEX_CACHE:
        _CELLOSAURUS_INDEX_CACHE.clear()
        _CELLOSAURUS_INDEX_CACHE[key] = CellosaurusIndex(table)
    return _CELLOSAURUS_INDEX_CACHE[key]


def is_cell_line_key(key):
    k = str(key).strip().lower()
    return bool(_CELL_KEY_RE.search(k)) and not _CELL_KEY_EXCLUDE_RE.search(k)


def is_cultured_cell_study(subject_type):
    """SUBJECT_TYPE of cultured (non-microbial) cells: 'Cultured cells', 'Cells', 'Human cells', ..."""
    s = _norm_value(subject_type)
    return bool(re.search(r'\b(?:cells?|cell lines?|cultured|cultures?)\b', s)) and \
        not re.search(r'bacteri|yeast|fung|microb|plant|alga|archae', s)


def _cl_lookup(ontology_table, allowedTerm_dict):
    """lower-case Cell Ontology label / synonym (> 3 characters) -> label, for names of
    exactly one allowed, not too generic CL label."""
    key = (id(ontology_table), len(ontology_table))
    if key in _CL_LOOKUP_CACHE:
        return _CL_LOOKUP_CACHE[key]
    allowed = set((allowedTerm_dict.get('UBERONBodyPartName') or {}).get('allowed_values') or [])
    names = {}
    labels = set()
    if ontology_table is not None and 'Synonym' in ontology_table.columns:
        cl = ontology_table[ontology_table['UBERONOntologyIndex'].astype(str).str.startswith('CL:')]
        for label, syn in cl[['Label', 'Synonym']].itertuples(index=False):
            if not isinstance(label, str) or label not in allowed or label.lower().startswith('obsolete') \
                    or label.lower() in _CL_TOO_GENERIC:
                continue
            labels.add(label.strip().lower())
            for n in (label, syn):
                if isinstance(n, str) and len(n.strip()) > 3:
                    names.setdefault(n.strip().lower(), set()).add(label)
    # a name that is a UBERON/PO label too is left alone
    other = set()
    if ontology_table is not None:
        other = set(ontology_table.loc[~ontology_table['UBERONOntologyIndex'].astype(str).str.startswith('CL:'), 'Label']
                    .dropna().astype(str).str.strip().str.lower())
    lookup = {n: next(iter(v)) for n, v in names.items() if len(v) == 1 and (n in labels or n not in other)}
    _CL_LABELS.clear()
    _CL_LABELS.update(lookup.values())
    _CL_LOOKUP_CACHE.clear()
    _CL_LOOKUP_CACHE[key] = lookup
    return lookup


def cell_type_term(value, ontology_table, allowedTerm_dict):
    """Allowed Cell Ontology label a value names exactly (case-insensitive; also the
    plural and with a leading 'primary' / 'human' / 'mouse'), else None."""
    v = ' '.join(str(value).strip().lower().split())
    if len(v) <= 3:
        return None
    lookup = _cl_lookup(ontology_table, allowedTerm_dict)
    if not lookup:
        return None
    base = _CL_PREFIX_RE.sub('', v)
    for s in dict.fromkeys([v, base]):
        if len(s) <= 3:
            continue
        for cand in (s, re.sub(r'cells$', 'cell', s), re.sub(r'(?<=[a-z])s$', '', s)):
            if cand in lookup:
                return lookup[cand]
    return None


def sample_cell_values(df, ontology_table, allowedTerm_dict):
    """(cell-line values, cell type) per row of df: the sample's distinct values under
    cell-line keys that are no Cell Ontology cell type ('\\x1f'-joined), and the one CL
    label its values name (None if none or several)."""
    keys = df['Key'].astype(str)
    is_cell = keys.map(is_cell_line_key)
    if not is_cell.any():
        empty = pd.Series(None, index=df.index, dtype=object)
        return empty, empty.copy()
    sub = pd.DataFrame({'filename': df.loc[is_cell, 'filename'],
                        'v': df.loc[is_cell, 'Value'].astype(str).str.strip()})
    sub = sub[sub['v'] != '']
    sub['ct'] = sub['v'].map(lambda v: cell_type_term(v, ontology_table, allowedTerm_dict))
    lines = (sub[sub['ct'].isna()].groupby('filename')['v']
             .agg(lambda s: _CELL_VALUE_SEP.join(sorted(set(s)))))
    types = (sub.dropna(subset=['ct']).groupby('filename')['ct']
             .agg(lambda s: s.iloc[0] if s.nunique() == 1 else None).dropna())
    return df['filename'].map(lines), df['filename'].map(types)


def _lineage(taxid):
    import REDU_conversion_functions as _rcf
    return _rcf._get_ncbi_lineage(str(taxid))


def assign_cell_lines(MWB_table):
    """CellLine (Cellosaurus recommended name), Cell Ontology cell types for an empty
    UBERONBodyPartName and SampleType 'culture_mammalian' for rows given a mammalian
    cell line. Sample-level values first; the study-level CELL_STRAIN_DETAILS only when
    it names exactly one line and no sample-level value names another. Blanks/QC get
    nothing; a line has to come from the study organism (see CellosaurusIndex.match)."""
    t = MWB_table
    t['CellLine'] = np.nan
    n = len(t)
    if n == 0:
        return t
    st = t['SampleType'].astype(str) if 'SampleType' in t.columns else pd.Series('', index=t.index)
    not_subject = st.str.contains(r'blank|qc|pool', case=False) | st.isin(['reference material', 'molecular_standards'])
    subject_type = t['SUBJECT_TYPE'].iloc[0] if 'SUBJECT_TYPE' in t.columns else ''
    cultured = is_cultured_cell_study(subject_type)

    def _missing(col):
        if col not in t.columns:
            t[col] = np.nan
        return t[col].isna() | t[col].astype(str).str.strip().str.lower().isin(['', 'nan', 'none', 'na', 'missing value'])

    sample_type_cl = t['MWB_celltype'] if 'MWB_celltype' in t.columns else pd.Series(None, index=t.index, dtype=object)
    celltype = sample_type_cl
    if 'STUDY_celltype' in t.columns and celltype.isna().all():
        celltype = t['STUDY_celltype']
    values = t['MWB_cellline_values'] if 'MWB_cellline_values' in t.columns else pd.Series(None, index=t.index, dtype=object)

    def _fill_cell_types(line_ac):
        # cell types -> UBERONBodyPartName (Cell Ontology), only where empty. Sample-level
        # cell types only when every sample with a cell value got a cell type or a line
        # (no half-labelled designs, e.g. 'Reticulocyte' vs an unmatched 'Mature erythrocyte')
        from cellosaurus import is_rejected_value
        if sample_type_cl.notna().any():
            open_value = values.map(lambda v: isinstance(v, str) and any(
                not is_rejected_value(x) for x in v.split(_CELL_VALUE_SEP)))
            unresolved = open_value & sample_type_cl.isna() & line_ac.isna() & ~not_subject
            if unresolved.any():
                print(f"[CellLine] cell types not used: {int(unresolved.sum())} rows have other unmatched cell values")
                return
        fill = _missing('UBERONBodyPartName') & celltype.notna() & ~not_subject
        if fill.any():
            t.loc[fill, 'UBERONBodyPartName'] = celltype[fill]
            print(f"[CellLine] cell type as UBERONBodyPartName for {int(fill.sum())} rows: {sorted(set(celltype[fill]))}")

    index = _cellosaurus_index()
    if index is None:
        _fill_cell_types(pd.Series(None, index=t.index, dtype=object))
        return t
    taxon = t['NCBITaxonomy'] if 'NCBITaxonomy' in t.columns else pd.Series(None, index=t.index, dtype=object)
    # (a sample of several organisms gives no study organism to match cell lines with)
    taxon = taxon.where(taxon.astype(str).str.contains(r'^\d+\|', regex=True) & ~taxon.map(is_multi).astype(bool), '')
    cache = {}
    log = {}

    def _resolve(value, tax):
        k = (value, tax)
        if k not in cache:
            tid, _, tname = str(tax).partition('|')
            ac, reason = index.match(value, tid, tname, cultured=cultured, lineage_fn=_lineage)
            cache[k] = ac
            log[value] = (index.info[ac]['name'] if ac else None, reason)
        return cache[k]

    sample_ac = pd.Series(None, index=t.index, dtype=object)
    if 'MWB_cellline_values' in t.columns:
        for i, vals, tax in zip(t.index, t['MWB_cellline_values'], taxon):
            if not isinstance(vals, str) or not vals or not_subject[i]:
                continue
            found = {_resolve(v, tax) for v in vals.split(_CELL_VALUE_SEP)} - {None}
            if len(found) == 1:
                sample_ac[i] = next(iter(found))
            elif len(found) > 1:
                log[vals.replace(_CELL_VALUE_SEP, ' / ')] = (None, 'several cell lines in one sample')
    ac = sample_ac.copy()
    study_value = str(t['CELL_STRAIN_DETAILS'].iloc[0]).strip() if 'CELL_STRAIN_DETAILS' in t.columns else ''
    study_is_type = 'STUDY_celltype' in t.columns and t['STUDY_celltype'].notna().any()
    if study_value and study_value.lower() not in ('na', 'nan', 'none', '-', '') and not study_is_type:
        study_ac = pd.Series([_resolve(study_value, tax) for tax in taxon], index=t.index, dtype=object)
        lines = set(study_ac[~not_subject]) - {None}
        other = set(sample_ac.dropna()) - lines
        if len(lines) == 1 and study_ac[~not_subject].notna().all() and not other:
            fill = ac.isna() & ~not_subject
            ac[fill] = next(iter(lines))
        elif len(lines) == 1 and other:
            log[study_value] = (None, 'study-level line disagrees with sample-level lines')
    for value, (name, reason) in sorted(log.items()):
        print(f"[CellLine] {value!r}: {name if name else 'not assigned'} ({reason})")
    _fill_cell_types(ac)
    has = ac.notna()
    if not has.any():
        return t
    t.loc[has, 'CellLine'] = ac[has].map(lambda a: index.info[a]['name'])
    print(f"[CellLine] {int(has.sum())} rows: {dict(Counter(t.loc[has, 'CellLine']))}")

    # SampleType culture_mammalian for mammalian lines, only over 'animal' / no type;
    # outside cultured-cell studies only for rows that are cells (no body part, or a cell type)
    mammal = ac.map(lambda a: bool(a) and index.info[a]['mammal'] if isinstance(a, str) else False).astype(bool)
    open_type = t['SampleType'].isna() | t['SampleType'].astype(str).str.strip().str.lower().isin(
        ['animal', '', 'nan', 'none', 'missing value'])
    is_cells = pd.Series(cultured, index=t.index) | _missing('UBERONBodyPartName') | t['UBERONBodyPartName'].isin(
        _CL_LABELS | set(celltype.dropna()))
    change = has & mammal & open_type & ~not_subject & is_cells
    if change.any():
        t.loc[change, 'SampleType'] = 'culture_mammalian'
        t.loc[change, 'SampleTypeSub1'] = 'culture_mammalian'
        print(f"[CellLine] SampleType culture_mammalian for {int(change.sum())} rows")
    return t


# COLLECTION SAMPLE_TYPE and SUBJECT_TYPE of the analysis being converted (create_dataframe_outer_dict)
_ANALYSIS_SAMPLE_TYPE = ''
_ANALYSIS_SUBJECT_TYPE = ''
_ANALYSIS_TAXON = None   # the analysis' NCBITaxonomy term (create_dataframe_outer_dict), for sample_body_parts
_NOT_CULTURE_SAMPLE_RE = re.compile(r'(?i)fraction|microsom|\bs9\b|homogen|tissue|\borgans?\b|biops|plasma|serum|'
                                    r'\bblood\b|urine|feces|faeces|stool|\bmilk\b|saliva')


def cultured_study_sample_type(t):
    """SampleType culture_mammalian for the 'animal' rows of a cultured-cell study
    (SUBJECT_TYPE 'Cultured cells', ...) of a mammal that got no Cellosaurus line
    (primary cells, 'HEK cells', 'Fibroblasts'): rows that are cells (no body part, or a
    cell type), not blanks / QC."""
    if len(t) == 0 or 'SUBJECT_TYPE' not in t.columns or not is_cultured_cell_study(t['SUBJECT_TYPE'].iloc[0]):
        return t
    # subcellular fractions, tissues and body fluids of a 'Cultured cells' study are no cell cultures
    sample_type_text = str(t['SAMPLE_TYPE'].iloc[0]) if 'SAMPLE_TYPE' in t.columns else str(_ANALYSIS_SAMPLE_TYPE or '')
    if _NOT_CULTURE_SAMPLE_RE.search(sample_type_text):
        return t
    st = t['SampleType'].astype(str).str.strip().str.lower()
    body = t['UBERONBodyPartName'] if 'UBERONBodyPartName' in t.columns else pd.Series(np.nan, index=t.index)
    is_cells = body.isna() | body.astype(str).str.strip().str.lower().isin(['', 'nan', 'none', 'na', 'missing value']) \
        | body.isin(_CL_LABELS)
    tax = t['NCBITaxonomy'].astype(str) if 'NCBITaxonomy' in t.columns else pd.Series('', index=t.index)
    mammal_of = {}
    def _mammal(x):
        m = re.match(r'^(\d+)\|', x)
        lineage = _lineage(m.group(1)) if m else None
        return lineage is not None and 'mammalia' in lineage
    for x in sorted(set(tax[st.eq('animal') & is_cells])):
        # (several organisms: all of them mammals)
        mammal_of[x] = all_parts(x, _mammal) if is_multi(x) else _mammal(x)
    change = st.eq('animal') & is_cells & tax.map(lambda x: mammal_of.get(x, False)).astype(bool)
    if change.any():
        t.loc[change, 'SampleType'] = 'culture_mammalian'
        t.loc[change, 'SampleTypeSub1'] = 'culture_mammalian'
        print(f"[CellLine] SampleType culture_mammalian for {int(change.sum())} rows of a cultured-cell study")
    return t


def _ncbi_name_lookup(allowedTerm_dict):
    """lower-cased NCBI name -> first allowed 'ID|Name' term with that name (the index
    get_taxonomy_id_from_name__allowedTerms uses; one copy of the ~3M names)."""
    return _rcf._taxonomy_term_index(allowedTerm_dict['NCBITaxonomy']['allowed_values'])[1]


def _sp_name_term(text, term, given_id, allowedTerm_dict):
    """_rcf.sp_name_term of an organism text and the term the converter's lookups found for it (by the given ID
    given_id, '' if none, or by name): for 'Genus sp. <designation>', the genus where the designation is a lab's
    numbering ('Bacillus sp. 1') and the term is not the given ID's taxon, or where the given ID is a strain of another
    designation that the text's own name (allowed name, names.dmp synonym) does not name; else term."""
    if not term or not _rcf.sp_designation(text):
        return term
    by_id, by_name = _rcf._taxonomy_term_index(allowedTerm_dict['NCBITaxonomy']['allowed_values'])
    return _rcf.sp_name_term(text, term, given_id, by_id.get, lambda n: by_name.get(str(n).lower()),
                             _sp_text_term(text, allowedTerm_dict))


def _sp_text_term(text, allowedTerm_dict):
    """The term an organism text names by its allowed name or a names.dmp synonym (white space normalised), else None."""
    _, by_name = _rcf._taxonomy_term_index(allowedTerm_dict['NCBITaxonomy']['allowed_values'])
    return by_name.get(' '.join(str(text).lower().split())) or _rcf.pinned_ncbi_synonym_term(str(text), allowedTerm_dict)


def _sp_given_id(text, term, given_ids, allowedTerm_dict):
    """The ID of `term` where it is the term of the analysis' TAXONOMY_ID (its ID one of given_ids) and, for a text
    'Genus sp. <designation>', that ID is the text's own taxon (_sp_text_term) or lies in the text's genus (the pinned
    NCBI genus of its taxon has the text's genus name: _rcf.pinned_genus_term), as the validator and MetaboLights take an
    ID over such a text; else '' (_sp_name_term: an ID of an unrelated taxon, 'Streptomyces sp. 12' with 9606, is no
    evidence that the text names the isolate NCBI lists under that name)."""
    taxid = term.split('|', 1)[0].strip()
    if taxid not in given_ids:
        return ''
    if not _rcf.sp_designation(text):
        return taxid
    by_id, _ = _rcf._taxonomy_term_index(allowedTerm_dict['NCBITaxonomy']['allowed_values'])
    if _rcf.pinned_genus_term(taxid, str(text).split()[0], by_id.get) or _sp_text_term(text, allowedTerm_dict) == term:
        return taxid
    return ''


def _bare_sp_term(species, given_ids, allowedTerm_dict):
    """The term of the analysis' TAXONOMY_ID (one of given_ids) where SUBJECT_SPECIES is a bare 'Genus sp.' and the ID
    that of a species or strain of that genus (_rcf.bare_sp_id_term; one such term), else None."""
    if not _rcf.genus_level_sp(species):
        return None
    by_id, by_name = _rcf._taxonomy_term_index(allowedTerm_dict['NCBITaxonomy']['allowed_values'])
    terms = {_rcf.bare_sp_id_term(species, i, by_id.get, lambda n: by_name.get(str(n).strip().lower()))
             for i in sorted(given_ids) if i.isdigit()} - {None}
    return terms.pop() if len(terms) == 1 else None


_organism_fallback_memo = {}


def mwb_organism_fallback(text, allowedTerm_dict, where=''):
    """Allowed term of an organism text that neither the allowed names nor names.dmp synonyms resolve: other white
    space (a no-break space), a strain designation at the end ('... ATCC 31916'), one wrong character
    ('Arbidopsis thaliana'); see
    REDU_conversion_functions.ncbi_name_fallback (pinned taxonomy only). Else None. Logged where it is used."""
    if not isinstance(text, str) or _rcf._PINNED_LINEAGE is None:
        return None
    if text not in _organism_fallback_memo:
        by_id, by_name = _rcf._taxonomy_term_index(allowedTerm_dict['NCBITaxonomy']['allowed_values'])

        def exact(name):
            return by_name.get(name.lower()) or _rcf.pinned_ncbi_synonym_term(name, allowedTerm_dict)
        _organism_fallback_memo[text] = _rcf.ncbi_name_fallback(text, by_id.get, exact)
    term, how = _organism_fallback_memo[text]
    if term:
        print(f"[NCBITaxonomy] {where}: organism text {text!r} read as {term} ({how})")
    return term


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



# a name that ends in an MS file extension, then '.' and a space: 'a.raw, b.raw. c.raw'
_NAME_DOT_SEP_RE = re.compile(r'(\.(?:raw|mzml|mzxml|wiff|cdf|d))\.\s+', re.I)
# placeholders between listed names ('a.mzML b.mzML - -')
_RAW_NAME_PLACEHOLDERS = {'-', '--', 'na', 'n/a', 'none'}
# one raw file name with an MS file extension (also the .wiff.scan companion)
_NAMED_FILE_RE = re.compile(r'\S+\.(?:raw|mzml|mzxml|wiff|wiff\.scan|cdf|d)', re.I)


def split_raw_file_names(value, split_at_spaces):
    """The raw file names of one RAW_FILE_NAME value (sorted, distinct): split at ';' and
    ',' (and at spaces when no listed file name has a space), and also at '. ' after a
    file extension ('a.raw, b.raw. c.raw'), at '/' between names that all have an MS
    file extension ('x.RAW/x_diluted.RAW'), and at spaces between names that all have one
    ('ESI+1.wiff ESI+1.wiff.scan ESI-.wiff', also with '-' placeholders between them) when a
    listed file name has a space."""
    v = _NAME_DOT_SEP_RE.sub(r'\1;', str(value))
    out = []
    for p in re.split(r'; ?|,| ' if split_at_spaces else r'; ?|,', v):
        words = p.split()
        named = sum(1 for w in words if _NAMED_FILE_RE.fullmatch(w))
        if not split_at_spaces and named > 1 and \
                all(_NAMED_FILE_RE.fullmatch(w) or w.lower() in _RAW_NAME_PLACEHOLDERS for w in words):
            out.extend(words)
            continue
        pieces = p.split('/')
        if len(pieces) > 1 and all(_NAMED_FILE_RE.fullmatch(x) for x in pieces):
            out.extend(pieces)
            continue
        out.append(p)
    return sorted(set(out))


def _name_stem(name):
    return os.path.basename(str(name).replace('\\', '/')).lower().rsplit('.', 1)[0]


def _own_file_key(name):
    """A listed file name or a raw file name for matching own files: stem, lower case, '-' = '_'."""
    return _name_stem(name).replace('-', '_')


def _own_id_key(sample_id):
    """A Sample ID for matching own files: lower case, '-' = '_' (no extension removed)."""
    return str(sample_id).strip().lower().replace('-', '_')


# Sample ID tokens of standards (blanks and QCs: blank_qc_names)
_STANDARD_TOKEN_RE = re.compile(r'(?i)^(?:std|stds|standards?|istds?|calibrants?|calibration|calibrators?|nist|srm)\d*$')
# own files found by the Sample ID alone need an ID at least this long, with a letter: '1', 'A1',
# '0012' name unrelated files too easily
_OWN_ID_MIN_LEN = 5

OWN_FILE_REASONS = ('shared name not a listed file', 'blank / QC / standard', 'Sample ID and file name differ in case',
                    'Sample ID only: several MS analyses',
                    'Sample ID only: short ID or no letter', 'Sample ID only: sample names several shared files',
                    'several files fit', 'file named by a sample', 'other polarity / chromatography',
                    'file wanted by several samples')


def not_a_study_sample(sample_id, values=()):
    """Blank, QC or standard: by a token of the Sample ID (blank_qc_names; 'STD', 'NIST') or by a
    QC / blank / reference-material value of the sample (qc_sample_type)."""
    if blank_qc_type_from_filename(sample_id):
        return True
    if any(_STANDARD_TOKEN_RE.match(t) for t in name_tokens(sample_id)):
        return True
    return any(qc_sample_type(k, v) for k, v in values)


def own_file_matches(pairs, listing, hints=None, is_not_study_sample=None, contexts=None):
    """The own file of each sample that names a raw file several samples of the analysis name
    (a multi-sample SCIEX .wiff, or one file named for several samples).

    pairs: (raw file name, Sample ID) of every sample of the analysis; listing: (filename_base,
    path) of the files kept for the analysis (one per stem and acquisition context).
    Returns ({(raw name lower case, Sample ID): (filename_base, form)}, [candidate details]).

    Forms, first that fits:
      '<file>-<Sample ID>'    '<file stem>-<Sample ID>' or '_' (case-insensitive; the 6f46b12
                              rule, unchanged and without the guards below)
      'normalised'            '<file stem>_<Sample ID>' with '-' and '_' alike
      '<Sample ID> tied'      the Sample ID itself ('-' = '_'), where it starts with the file stem
                              ('<file stem>_x')
      '<Sample ID>'           the Sample ID itself, not tied to the shared file's name
    In these three forms the Sample ID part of the file name must have the Sample ID's case
    (the shared file's stem may differ in case). A candidate of them is left out (details give
    the first reason and every guard it fails) when the shared name is not a listed file
    (placeholders like '-', '.mzML'); for blanks, QCs and standards; when the Sample ID differs
    in case; for the untied form, in studies with several MS analyses (the file may be another
    analysis' acquisition), for Sample IDs shorter than 5 characters or without a letter, and
    for samples naming several shared files; when several listed files fit; when a sample of the
    analysis names the file (or another file of its stem, e.g. x.raw for x.mzML) as its raw
    file; when its path states another polarity / chromatography than the analysis or the
    shared file; and when the file would be the own file of several samples (or shared files of
    different stems). contexts: {path: (polarity, chromatography)} of the listing (file_context,
    with the archive names); a path not in it has its path's context."""
    hints = hints or {}
    contexts = contexts or {}
    pairs = sorted({(str(r).lower(), str(s)) for r, s in pairs})
    naming = {}
    for r, s in pairs:
        naming.setdefault(r, set()).add(s)
    shared = {r for r, ss in naming.items() if len(ss) > 1} - {'na', 'nan', ''}
    matches, details = {}, []
    if not shared:
        return matches, details
    listing = sorted({(str(b), str(p)) for b, p in listing})
    by_stem, by_key, paths_of_stem = {}, {}, {}
    for b, p in listing:
        by_stem.setdefault(_name_stem(b), b)
        by_key.setdefault(_own_file_key(b), []).append((b, p))
        paths_of_stem.setdefault(_name_stem(b), []).append(p)
    named = {}
    for r, s in pairs:
        named.setdefault(_own_file_key(r), set()).add(s)
    shared_of = {}   # listed shared files (stems) a sample names
    for r, s in pairs:
        if r in shared and _name_stem(r) in paths_of_stem:
            shared_of.setdefault(s, set()).add(_name_stem(r))
    proposals = {}
    for r, s in pairs:
        if r not in shared:
            continue
        stem = _name_stem(r)
        hit = next((by_stem[k] for k in (f'{stem}-{s.lower()}', f'{stem}_{s.lower()}') if k in by_stem), None)
        if hit is not None:
            matches[(r, s)] = (hit, '<file>-<Sample ID>')
            continue
        bkey, skey = _own_file_key(r), _own_id_key(s)
        if not skey or skey == bkey:
            continue
        if f'{bkey}_{skey}' in by_key:
            form, hits = 'normalised', by_key[f'{bkey}_{skey}']
        elif skey in by_key:
            tied = skey.startswith(bkey.rstrip('_') + '_')
            form, hits = '<Sample ID> tied' if tied else '<Sample ID>', by_key[skey]
        else:
            continue
        batch_paths = paths_of_stem.get(stem, [])
        untied = form == '<Sample ID>'
        # the Sample ID's part of the file name in the Sample ID's case ('-' = '_'): no own file
        # of the Workbench needed a case-insensitive match
        _id = str(s).strip().replace('-', '_')
        exact = [(b, p) for b, p in hits if os.path.basename(b.replace('\\', '/')).rsplit('.', 1)[0].replace('-', '_').endswith(_id)]
        hits = exact or hits
        ctx = contexts.get(hits[0][1]) or _path_indicates(hits[0][1])
        batch_ctx = {contexts.get(p) or _path_indicates(p) for p in batch_paths}
        # every guard the candidate fails (the first is the reason given)
        fails = [why for why, bad in (
            ('shared name not a listed file', not batch_paths),
            ('blank / QC / standard', (is_not_study_sample or not_a_study_sample)(s)),
            ('Sample ID and file name differ in case', not exact),
            ('Sample ID only: several MS analyses', untied and hints.get('multi')),
            ('Sample ID only: short ID or no letter', untied and (len(s.strip()) < _OWN_ID_MIN_LEN or not re.search(r'[a-z]', s, re.I))),
            ('Sample ID only: sample names several shared files', untied and len(shared_of.get(s, ())) > 1),
            ('several files fit', len(hits) > 1),
            ('file named by a sample', any(named.get(_own_file_key(b)) for b, _ in hits)),
            ('other polarity / chromatography', _context_contradicts(ctx, hints) or (len(batch_ctx) == 1 and _context_contradicts(
                ctx, dict(zip(('polarity', 'chrom'), next(iter(batch_ctx))))))),
        ) if bad]
        reason = fails[0] if fails else None
        details.append([r, s, form, hits[0][0] if len(hits) == 1 else ';'.join(b for b, _ in hits), reason, fails])
        if reason is None:
            proposals[(r, s)] = (hits[0][0], form, len(details) - 1)
    # one file per sample and shared file: a sample naming the shared file twice ('x.wiff' and
    # 'x.wiff - - -', one stem) wants the file once
    wanted = {}
    for (r, s), v in list(matches.items()) + list(proposals.items()):
        wanted.setdefault(v[0], set()).add((_name_stem(r), s))
    for k, (b, form, d) in proposals.items():
        if len(wanted[b]) > 1:
            details[d][4] = 'file wanted by several samples'
            details[d][5] = ['file wanted by several samples']
        else:
            matches[k] = (b, form)
    return matches, details


def per_sample_files_of_shared_names(df, raw_file_name_df, hints=None, sample_values=None, contexts=None):
    """filename_raw per df row, where several samples name one raw file (a multi-sample
    SCIEX .wiff) whose per-sample files are listed too: the sample's own file name instead of
    the shared one (own_file_matches). Other rows unchanged. sample_values: function giving the
    (Sample ID, key, value) rows of the samples, to leave out blanks / QCs / standards.
    contexts: {path: (polarity, chromatography)} of the study's listing (listing_contexts)."""
    raw = df['filename_raw'].astype(str)
    sid = df['_sample_id'].astype(str) if '_sample_id' in df.columns else df['filename'].astype(str)
    pairs = pd.DataFrame({'raw': raw.str.lower(), 'sid': sid}).drop_duplicates()
    n_samples = pairs.groupby('raw')['sid'].nunique()
    shared = set(n_samples.index[n_samples > 1]) - {'na', 'nan', ''}
    if not shared:
        return df['filename_raw']
    _values = {}

    def _not_study_sample(s):
        if sample_values is not None and not _values:
            _values[None] = []   # built once, at the first candidate
            for v_sid, k, v in sample_values():
                _values.setdefault(str(v_sid), []).append((k, v))
        return not_a_study_sample(s, _values.get(s, ()))

    matches, details = own_file_matches(zip(raw, sid), zip(raw_file_name_df['filename_base'], raw_file_name_df['filename_raw_path']),
                                        hints=hints, is_not_study_sample=_not_study_sample, contexts=contexts)
    out = df['filename_raw'].copy()
    n = 0
    for i, (r, s) in enumerate(zip(raw, sid)):
        m = matches.get((r.lower(), s))
        if m is not None:
            out.iat[i] = m[0]
            n += 1
    if n:
        forms = Counter(f for _, f in matches.values())
        print(f"[rawfile] {n} rows ({len(matches)} samples x shared files): a raw file named by several samples replaced by the "
              f"sample's own file (" + ', '.join(f"{f}: {c}" for f, c in sorted(forms.items())) + ")")
    left = Counter(d[4] for d in details if d[4])
    if left:
        print("[rawfile] own-file candidates left out: " + ', '.join(f"{r}: {left[r]}" for r in OWN_FILE_REASONS if left[r]))
    return out


def match_unmatched_names_by_prefix(df, raw_file_name_df):
    """After an exact file-name match of an analysis' raw file names (some names matched): a
    listed name that matched nothing (samples listing one raw file per platform,
    'x-N;x-P;20160915_S1;20160913_S1', whose lipidomics names lack the '_N' / '_P' suffix)
    gets the one file nothing else matched whose name starts with the listed name and a
    separator ('20160915_S1' -> NEG/BATCH1/20160915_S1_N.wiff). Names with a file extension,
    names that fit no such file or several files, names fitting a file another name also
    fits, and names fitting a file named after another sample of the analysis (its Sample ID or
    listed name is the file name, compared as in substring_file_matches) stay unmatched."""
    if 'USI' not in df.columns:
        return df
    todo = df['USI'].isna()
    if not todo.any():
        return df
    # names without a file extension only: a full file name that is not listed is another file
    names = sorted({str(n) for n in df.loc[todo, 'filename_raw']
                    if len(str(n)) >= 4 and str(n).lower() not in ('na', 'nan') and not _NAMED_FILE_RE.fullmatch(str(n))})
    claimed = set(df['USI'].dropna())
    free = raw_file_name_df[~raw_file_name_df['USI'].isin(claimed)]
    by_prefix = {}
    for lbl, base in zip(free.index, free['filename_base'].astype(str)):
        stem = _name_stem(base)
        for m in re.finditer(r'[-_.]', stem):
            by_prefix.setdefault(stem[:m.start()], []).append(lbl)
    hits = {}
    for n in names:
        lbls = by_prefix.get(os.path.basename(n.replace('\\', '/')).lower(), [])
        if len(lbls) == 1:
            hits[n] = lbls[0]
    used = Counter(hits.values())
    hits = {n: l for n, l in hits.items() if used[l] == 1}
    # never a file named after another sample of the analysis (its Sample ID or listed name is the file name)
    if hits and '_sample_id' in df.columns:
        named = {}
        for r, s in zip(df['filename_raw'].astype(str), df['_sample_id'].astype(str)):
            for t in (_match_tokens(_name_stem(r)), _match_tokens(s)):
                if t:
                    named.setdefault(t, set()).add(s)
        own = {}
        for r, s in zip(df.loc[todo, 'filename_raw'].astype(str), df.loc[todo, '_sample_id'].astype(str)):
            own.setdefault(r, set()).add(s)
        hits = {n: l for n, l in hits.items()
                if not (named.get(_match_tokens(_name_stem(raw_file_name_df.loc[l, 'filename_base']))) or set()) - own.get(n, set())}
    if not hits:
        return df
    df = df.copy()
    rows = todo & df['filename_raw'].astype(str).isin(hits)
    cols = [c for c in raw_file_name_df.columns if c in df.columns]
    fill = raw_file_name_df.loc[df.loc[rows, 'filename_raw'].astype(str).map(hits), cols]
    fill.index = df.index[rows]
    for c in cols:
        df.loc[rows, c] = fill[c]
    print(f"[rawfile] {len(hits)} listed raw file names without extension matched by name prefix ('<name>_N.wiff')")
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

# a trailing format word of a listed raw file name: 'L32650_MzML' is the file L32650.mzML
_LISTED_FORMAT_WORD_RE = re.compile(r'[\s._-]+(?:mzml|mzxml|mzdata|raw|wiff|cdf|mgf)$', re.I)
# listed raw file names that name no file (fillna('NA'), placeholders)
_NO_LISTED_NAME = {'', 'na', 'nan', 'none', 'null', 'n/a', '-', '--'}
_NAME_TOKEN_RE = re.compile(r'[a-z]+|[0-9]+')
# word boundaries inside a name written in camel case: 'fooBar', 'CSHneg'
_CAMEL_RE = re.compile(r'(?<=[a-z])(?=[A-Z])|(?<=[A-Z]{2})(?=[a-z])')


def _match_tokens(text):
    """letter and digit runs of a name, lower case ('QC-1' and 'qc1' -> ('qc', '1'); 'KPF_1_1' is not 'KPF_11')"""
    return tuple(_NAME_TOKEN_RE.findall(str(text).lower()))


def _listed_name_tokens(name):
    """tokens of a listed raw file name (before the first '.') without a trailing format word"""
    s = str(name)
    return _match_tokens(_LISTED_FORMAT_WORD_RE.sub('', s) or s)


def _token_run_in(short, long):
    n = len(short)
    return n > 0 and any(long[i:i + n] == short for i in range(len(long) - n + 1))


def _run_positions(short, long):
    n = len(short)
    return [i for i in range(len(long) - n + 1) if n and long[i:i + n] == short]


def _camel_tokens(text):
    """letter and digit runs of a name as written, also split at 'fooBar' and 'CSHneg' (lower case)"""
    s = _LISTED_FORMAT_WORD_RE.sub('', str(text)) or str(text)
    return _match_tokens(_CAMEL_RE.sub(' ', s))


_POLARITY_WORDS = {'pos', 'neg', 'positive', 'negative'}


def _camel_words(text):
    """(tokens, word start flags, word end flags) of a name as written: its letter and digit runs (lower case) and
    whether each starts / ends a word; words are split at non-alphanumerics, at 'fooBar' / 'CSHneg', and around a
    polarity word ('1NEG' is 1 + neg)"""
    s = _LISTED_FORMAT_WORD_RE.sub('', str(text)) or str(text)
    toks, starts, ends = [], [], []
    for w in re.findall(r'[A-Za-z0-9]+', _CAMEL_RE.sub(' ', s)):
        t = _NAME_TOKEN_RE.findall(w.lower())
        for i, x in enumerate(t):
            toks.append(x)
            starts.append(i == 0 or x in _POLARITY_WORDS or t[i - 1] in _POLARITY_WORDS)
            ends.append(i == len(t) - 1 or x in _POLARITY_WORDS or t[i + 1] in _POLARITY_WORDS)
    return tuple(toks), starts, ends


def _weak_name(tokens):
    """a name that pairs a file by containment only as whole words ('T1' in 'T1_Neg', not '1' in '1ppmCoc'):
    shorter than 3 characters, or fewer than 4 digits only"""
    j = ''.join(tokens)
    return len(j) < 3 or (j.isdigit() and len(j) < 4)


SUBSTRING_MATCH_REASONS = ('named after the sample', 'folder named after the sample',
                           'named after the sample without the common part', 'listed name is the folders and the file name',
                           'closest name', 'named after another sample', 'named after several samples',
                           'several names fit equally', 'no whole-token match', 'blank / QC file, study sample name',
                           'only a placeholder name')


def _path_folder_tokens(path):
    """tokens of the folders of a path ('a/B_1/x.mzML' -> a b 1) and of the parent folder alone (b 1)"""
    parts = [p for p in re.split(r'[\\/]+', str(path)) if p][:-1]
    return _match_tokens(' '.join(parts)), (_match_tokens(parts[-1]) if parts else ())


def _common_affix(token_lists):
    """the leading and the trailing tokens every name has (() when one name is all of it)"""
    lists = [t for t in token_lists if t]
    if len(lists) < 2:
        return (), ()
    pre = lists[0]
    for t in lists[1:]:
        n = 0
        while n < min(len(pre), len(t)) and pre[n] == t[n]:
            n += 1
        pre = pre[:n]
    suf = lists[0][len(pre):]
    for t in lists[1:]:
        t = t[len(pre):]
        n = 0
        while n < min(len(suf), len(t)) and suf[len(suf) - 1 - n] == t[len(t) - 1 - n]:
            n += 1
        suf = suf[len(suf) - n:] if n else ()
    if any(len(t) <= len(pre) + len(suf) for t in lists):
        return (), ()
    return pre, suf


def _folder_key(path):
    """the folders of a path, lower case ('Raw Data/Pos/S1.d' -> 'raw data/pos')"""
    parts = [p for p in re.split(r'[\\/]+', str(path).lower()) if p]
    return '/'.join(parts[:-1])


def substring_file_matches(potential_matches, sample_ids_of, files, listed_names=None, reasons=None):
    """The samples of each listing row in the substring matching of listed raw file names (map_potential_matches: a
    listed name and a file name, lower case without extension, of which one contains the other), one sample per file.

    potential_matches: {listed name: set of file names}; sample_ids_of: {listed name: Sample IDs whose rows list it};
    files: (row label, file name, path) of the listing rows the matching uses; listed_names: {listed name: the name as
    listed (letter case)}; reasons: a dict that gets {row label: reason}. Names are equal when their runs of letters
    and digits are ('QC-1' = 'qc1'; 'Sample1' is not 'Sample10'), a listed name without a trailing format word
    ('L32650_MzML' = 'L32650'), a listed name with folders ('Seed/49.D/') also by its last part, its folders telling
    apart samples of one name. Per file, the first rule that names samples decides:
      1. the file is named after samples (a listed name or a Sample ID is the file name): it goes to that sample
         ('32650.mzML' is sample '32650', not 'L32650');
      2. its folder is named after samples ('FB0001/FB000101.CDF');
      3. it is named after samples once the leading / trailing parts all file names (of the analysis, or of the
         file's folder) share are left out ('23M034_M002' is 'M002' where every file starts with '23M034_');
      4. else, of the listed names that are the file's folders and its name ('Supernatant_neg13_G_t2_1' for
         Supernatant/neg/13_G_t2_1.raw) or that contain the file name or are contained in it as a run of whole tokens
         ('Batch1_S1_pos': 'S1', not 'S11'; not 'S1' in 'MAS1'; also in camel-case words: 'BruinCSH' in
         'BruinCSHneg'), the one all of whose tokens the file's path has, else the longest name in the file name,
         else the one listed name the file name is in. A name shorter than 3 characters or of fewer than 4 digits
         only pairs as whole words ('T1' in 'T1_Neg', '1' in '1NEG', not '1' in 'Inj01' or '1ppmCoc'); never by a
         name inside the part all files of the analysis or of the folder share where those files are named after
         other samples ('M034' in 'GC-23M034-M168' where every file of the folder starts with 'GC-23M034-M'); never a
         file named as a blank or QC (blank_qc_names) to samples that are none (not_a_study_sample, nor a blank / QC
         listed name).
    Rules 1-3 give the file to the named sample alone (also where other samples list the same name); rule 4 to every
    sample listing the name. The file goes to nobody when the samples a rule names are several, or none of them has a
    listed name that matched the file, or (rule 4) several names fit equally or none in whole tokens. Placeholders
    ('NA', '', '-') match nothing (a file only they contain: 'only a placeholder name'). Returns ({(listed name, Sample ID): [row labels]}, Counter of reasons: how files
    were paired, and why the others were not)."""
    listed_names = listed_names or {}
    names, last, folders, samples_of = {}, {}, {}, {}
    for v in potential_matches:
        s = str(v)
        parts = [p for p in re.split(r'[\\/]+', s) if p]
        names[v] = _listed_name_tokens(s)
        last[v] = _listed_name_tokens(parts[-1]) if parts else ()
        folders[v] = _match_tokens(' '.join(parts[:-1]))
        samples_of[v] = frozenset(str(x) for x in (sample_ids_of or {}).get(v, ()) if str(x).strip()) or frozenset([str(v)])

    def placeholder(v):
        return str(v).strip().lower() in _NO_LISTED_NAME or not names[v]

    # samples by the names that name them: listed names (whole, and their last part) and Sample IDs
    named_by = {}
    for v in sorted(potential_matches, key=str):
        if not placeholder(v):
            for t in {names[v], last[v]}:
                if t:
                    named_by.setdefault(t, set()).update(samples_of[v])
        for x in samples_of[v] if (sample_ids_of or {}).get(v) else ():
            t = _match_tokens(x)
            if t:
                named_by.setdefault(t, set()).add(x)
    claims, placeholder_claims = {}, set()
    for v in sorted(potential_matches, key=str):
        for t in potential_matches[v]:
            if placeholder(v):
                placeholder_claims.add(t)
            else:
                claims.setdefault(t, []).append(v)
    pre, suf = _common_affix([_match_tokens(t) for _, t, _ in files])
    by_folder = {}
    for _, t, path in files:
        by_folder.setdefault(_folder_key(path), []).append(_match_tokens(t))
    folder_affix = {k: _common_affix(v) for k, v in by_folder.items()}
    words = {label: _camel_words(re.split(r'[\\/]+', str(path).rstrip('/\\'))[-1].split('.')[0]) for label, _, path in files}
    camel = {label: w[0] for label, w in words.items()}
    by_folder_c = {}
    for label, _, path in files:
        by_folder_c.setdefault(_folder_key(path), []).append(camel[label])
    pre_c, suf_c = _common_affix(list(camel.values()))
    folder_affix_c = {k: _common_affix(v) for k, v in by_folder_c.items()}
    listed_w = {v: _camel_words(str(listed_names.get(v, v)).replace('\\', '/').split('/')[-1]) for v in potential_matches}
    listed_c = {v: w[0] for v, w in listed_w.items()}

    def cores(pre_, suf_):
        return [(i, j) for i in range(len(pre_) + 1) for j in range(len(suf_) + 1) if i or j]
    out, why = {}, Counter()

    def note(label, reason):
        why[reason] += 1
        if reasons is not None:
            reasons[label] = reason

    assigned = {}

    def give(pairs, label, reason):
        for p in sorted(pairs):
            out.setdefault(p, []).append(label)
            assigned.setdefault(label, set()).add(p[1])
        note(label, reason)

    def decide(owners, cands, label, reason, path_folders=None):
        """True when the samples a name names decide the file (paired or withheld)"""
        if not owners:
            return False
        fit = [v for v in cands if samples_of[v] & owners]
        if len(owners) > 1 and path_folders is not None:
            # samples of one name told apart by the folders their listed names give
            told = [v for v in fit if folders[v] and _token_run_in(folders[v], path_folders)]
            if told:
                owners = owners & set().union(*(samples_of[v] for v in told))
                fit = told
        if len(owners) == 1 and fit:
            sid = next(iter(owners))
            give({(v, sid) for v in fit}, label, reason)
        else:
            note(label, 'named after several samples' if len(owners) > 1 else 'named after another sample')
        return True

    # rules 1-3 first: which samples the files of each folder (and of the analysis) are named after
    later, named_in = [], {}
    for label, t, path in files:
        cands = claims.get(t, [])
        if not cands:
            if t in placeholder_claims:
                note(label, 'only a placeholder name')
            continue
        tt = _match_tokens(t)
        path_folders, parent = _path_folder_tokens(path)
        if decide(named_by.get(tt), cands, label, 'named after the sample', path_folders) or \
                (parent and decide(named_by.get(parent), cands, label, 'folder named after the sample')):
            pass
        else:
            # the leading / trailing tokens all files (of the analysis, or of the folder) share left out, fewest
            # first (the shared part may end in a token the sample names start with: '23M034_M002', '23M034_M040'
            # share '23 m 034 m')
            strips = sorted(set(cores(pre, suf)) | set(cores(*folder_affix.get(_folder_key(path), ((), ())))),
                            key=lambda x: (x[0] + x[1], x))
            core = next((tt[i:len(tt) - j] for i, j in strips if tt[i:len(tt) - j] and named_by.get(tt[i:len(tt) - j])), None)
            if not (core and decide(named_by.get(core), cands, label, 'named after the sample without the common part')):
                later.append((label, t, path, cands, tt, path_folders))
                continue
        for sid in assigned.get(label, ()):
            named_in.setdefault(_folder_key(path), set()).add(sid)
            named_in.setdefault(None, set()).add(sid)
    for label, t, path, cands, tt, path_folders in later:
        path_tokens = set(path_folders) | set(tt)
        in_folders = [v for v in cands if len(names[v]) > len(tt) and names[v][len(names[v]) - len(tt):] == tt
                      and set(names[v][:len(names[v]) - len(tt)]) <= set(path_folders)]
        if in_folders:
            reason, top = 'listed name is the folders and the file name', in_folders
        else:
            key = _folder_key(path)
            fpre, fsuf = folder_affix.get(key, ((), ()))
            tc = camel[label]
            fpre_c, fsuf_c = folder_affix_c.get(key, ((), ()))

            def run_inside(n, t_, affixes, v, flags=None):
                """n a run of whole tokens of t_ (of whole words, flags: word start / end flags of t_), not inside
                the part all files (of the analysis / of the folder) share where those files are named after other
                samples ('M034' in 'GC-23M034-M168')"""
                if not n or len(n) > len(t_):
                    return False
                sh = [(len(p_), len(t_) - len(s_)) for (p_, s_), k in affixes if named_in.get(k, set()) - samples_of[v]]
                return any(not any(i + len(n) <= a_ or i >= b_ for a_, b_ in sh)
                           and (flags is None or (flags[0][i] and flags[1][i + len(n) - 1])) for i in _run_positions(n, t_))

            def inside(v):
                """the listed name in the file name, as a run of whole tokens (also as written in camel case:
                'BruinCSH' in 'BruinCSHneg'); a short name ('T1', '706') as whole words"""
                c_aff = [((pre_c, suf_c), None), ((fpre_c, fsuf_c), key)]
                if _weak_name(names[v]):
                    return run_inside(listed_c[v], tc, c_aff, v, words[label][1:])
                return run_inside(names[v], tt, [((pre, suf), None), ((fpre, fsuf), key)], v) or \
                    run_inside(listed_c[v], tc, c_aff, v)

            def around(v):
                """the file name in the listed name, as a run of whole tokens (a short file name as whole words)"""
                if len(names[v]) <= len(tt):
                    return False
                if not _weak_name(tt):
                    return _token_run_in(tt, names[v])
                lc, ls, le = listed_w[v]
                return any(ls[i] and le[i + len(tc) - 1] for i in _run_positions(tc, lc))

            whole = [v for v in cands if inside(v) or around(v)]
            if not whole:
                note(label, 'no whole-token match')
                continue
            explained = [v for v in whole if set(names[v]) <= path_tokens]
            if explained and len(explained) < len(whole):
                whole = explained
            if any(inside(v) for v in whole):
                # several listed names in the file name: the longest ('S11', not 'S1', in 'Batch1_S11_pos')
                whole = [v for v in whole if inside(v)]
                best = max(len(_LISTED_FORMAT_WORD_RE.sub('', str(v)) or str(v)) for v in whole)
                top = [v for v in whole if len(_LISTED_FORMAT_WORD_RE.sub('', str(v)) or str(v)) == best]
            else:
                # the file name in several listed names (with other words each): only one sample's
                top = whole
            reason = 'closest name'
            # a file named as a blank or QC is not given to study samples by a part of its name
            if blank_qc_type_from_filename(path) and not any(
                    blank_qc_type_from_filename(str(listed_names.get(v, v))) or not_a_study_sample(x) for v in top for x in samples_of[v]):
                note(label, 'blank / QC file, study sample name')
                continue
        if len({samples_of[v] for v in top}) == 1:
            give({(v, s) for v in top for s in samples_of[v]}, label, reason)
        else:
            note(label, 'several names fit equally')
    # each name's files in the order of the listing rows given
    order = {label: i for i, (label, _, _) in enumerate(files)}
    return {k: sorted(v, key=order.get) for k, v in out.items()}, why


def resolve_contradictions(contradictions, potential_matches):
    for raw_value, conflicting_values in contradictions.items():
        better_match = None
        highest_percentage = 0
        
        raw_value_length = len(raw_value)
        
        # sorted: of equally long names the first wins, not the first in set order (which
        # changes with PYTHONHASHSEED)
        for df_value in sorted(conflicting_values):
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


# factor keys that say what a sample is (not what its subject was given or ate)
_FOOD_SAMPLE_KEY_RE = re.compile(r'sample[ _]?(?:type|matrix|material)|matrix|material|specimen(?:[ _]type)?|'
                                 r'food[ _](?:type|item|product)|product(?:[ _]type)?|commodity')
_NO_SAMPLE_SOURCE = {'', 'na', 'n/a', '-', 'none', 'nan', 'not applicable'}


def sample_food_terms(df):
    """GFOP label (FoodOntologyTerm) of the sample of each row (food_ontology_term), or None:
    from the sample source and the factors that say what the sample is ('Sample type',
    'Matrix', 'Food type'), then the analysis' SAMPLE_TYPE where the sample source says
    nothing. Samples of a 'Food item' study (SUBJECT_TYPE) and samples with 'Sample type:
    Food' are food samples: plain GFOP names count ('Food type: Beef'), and product names in
    their other factors ('Treatment: Ground Beef'). The other factors of other samples can
    only narrow the product (grape cultivar) or rule it out (digests). The organism is
    checked later (food_term_of_organism)."""
    gfop = allowedTerm_dict['FoodOntologyTerm']['allowed_values'] if 'FoodOntologyTerm' in allowedTerm_dict else []
    keyl = df['Key'].astype(str).str.strip().str.lower()
    is_source = keyl == 'sample source'
    described = keyl.str.fullmatch(_FOOD_SAMPLE_KEY_RE) & ~is_source
    other = ~is_source & ~described & ~keyl.str.contains(r'file|raw', regex=True)
    sample_type = _ANALYSIS_SAMPLE_TYPE if isinstance(_ANALYSIS_SAMPLE_TYPE, str) else ''
    food_study = str(_ANALYSIS_SUBJECT_TYPE).strip().lower() == 'food item'
    # most studies: no value names a food, no sample is declared food
    candidates = set(df.loc[is_source | described, 'Value'].dropna()) | {sample_type}
    if not food_study and not any(food_sample_marker(v) or food_ontology_term([v], gfop) for v in candidates):
        return pd.Series(None, index=df.index, dtype=object)
    terms = {}
    for f, g in pd.DataFrame({'f': df['filename'], 'v': df['Value'], 'source': is_source, 'described': described,
                              'other': other}).groupby('f', sort=False):
        source = list(g.loc[g['source'], 'v'])
        values = source + list(g.loc[g['described'], 'v'])
        if all(str(s).strip().lower() in _NO_SAMPLE_SOURCE or food_sample_marker(s) for s in source):
            values.append(sample_type)
        food = food_study or any(food_sample_marker(v) for v in values)
        others = list(g.loc[g['other'], 'v'])
        terms[f] = food_ontology_term(values, gfop, food_sample=food, extra_values=others if food else (),
                                      detail_values=() if food else others)
    return df['filename'].map(terms)


# ---- standard injections ----------------------------------------------------------------------------
# A sample whose own value says it is an injection of standards (a standard mixture, calibration standards, an
# external standard; the Workbench sample source '_STD_') is SampleType 'molecular_standards' with SampleTypeSub1
# 'defined_standards' (the terms Metabolights2REDU uses); complete_and_fill_REDU_table then clears its organism,
# body part and subject fields. Values as qc_sample_type reads them (_norm_value), text in brackets dropped
# ('standards (STD)'). Phrases that only name standards count under the keys that say what a sample is or
# which group it belongs to, and as the Subject ID ('Synthetic Standard'); bare words ('Standard', 'STD 3',
# 'Amino acid standard', 'Internal Std') only under the keys that say what the sample is: 'GROUP: Standard' (a
# standard-care group of patients), 'Genotype: Standard', 'timing of metabolite extract addition: standard' are
# no standards. Not after words that make it a
# sample of the study or a QC ('plasma based standard', 'pooled reference standard', 'quality control standard'), also
# in brackets ('Standard (plasma)'); 'reference standard' only under the sample-kind keys.
# Only for samples no QC / blank value typed (qc_sample_type; its 'QC externalstandard' stays blank_QC). Kept
# apart from qc_sample_type, which also decides which samples the file matching treats as no study samples.
STANDARD_SAMPLE_TYPE = ('molecular_standards', 'defined_standards')
STANDARD_SAMPLE_SOURCE = '_std_'
_STANDARD_KIND_KEYS = {'sample source', 'sample type', 'sampletype', 'samle type', 'type', 'sample class', 'class',
                       'matrix', 'sample matrix', 'material', 'sample material', 'specimen', 'specimen type',
                       'category', 'sample category'}
_STANDARD_PHRASE_KEYS = _STANDARD_KIND_KEYS | {'group', 'sample group', 'treatment', 'treatment group', 'factor',
                                               'factor group', 'condition', 'source', 'origin', 'subject id'}
_STD_WORD = r'(?:standards?|stds?)'
_STD_NUM = r'(?: [0-9]+[a-z]?)?'
_STANDARD_PHRASE_RE = re.compile(
    rf'{_STD_WORD} (?:mix|mixes|mixture|mixtures|cocktail|solution|solutions){_STD_NUM}'
    rf'|(?:mix|mixture) of {_STD_WORD}'
    rf'|(?:calibration|calibrant|external|authentic|analytical|chemical|synthetic|reference|retention time|alkane) '
    rf'{_STD_WORD}(?: (?:mix|mixture|solution))?{_STD_NUM}'
    rf'|(?:calibration|standard) curve(?: (?:point|sample|level))?{_STD_NUM}'
    rf'|calibration (?:point|level|sample|samples|solution|solutions){_STD_NUM}')
_STANDARD_WORD_RE = re.compile(rf'(?:[a-z0-9+]+ ){{0,3}}{_STD_WORD}{_STD_NUM}|{_STD_WORD} (?=[a-z]*[0-9])[0-9a-z]{{1,5}}'
                               rf'|istds?|internal {_STD_WORD}')
_NOT_STANDARD_WORDS = {'pool', 'pooled', 'plasma', 'serum', 'urine', 'blood', 'tissue', 'based', 'gold', 'qc', 'quality',
                       'diet', 'care', 'chow', 'growth', 'protocol', 'nist', 'srm', 'study', 'patient', 'patients',
                       'control', 'controls', 'treatment', 'sample', 'samples', 'group'}


def standard_sample_type(key, value):
    """STANDARD_SAMPLE_TYPE[0] if a sample's item (key, value) says it is a standard injection, else None."""
    k = _norm_value(key)
    if k not in _STANDARD_PHRASE_KEYS:
        return None
    # the words of the whole value, text in brackets too ('Standard (plasma)', 'Reference standard (NIST plasma)')
    if set(_norm_value(value).split()) & _NOT_STANDARD_WORDS:
        return None
    v = _norm_value(re.sub(r'\([^()]*\)', ' ', str(value)))
    if not v:
        return None
    if _STANDARD_PHRASE_RE.fullmatch(v):
        # 'reference standard' only says what a sample is under the sample-kind keys ('Group: Reference standard' is a
        # positive-control drug group)
        if v.startswith('reference ') and k not in _STANDARD_KIND_KEYS:
            return None
        return STANDARD_SAMPLE_TYPE[0]
    if k in _STANDARD_KIND_KEYS and _STANDARD_WORD_RE.fullmatch(v):
        return STANDARD_SAMPLE_TYPE[0]
    return None


# ---- the sample's own organism --------------------------------------------------------------------------
# The organism comes from the analysis' SUBJECT_SPECIES / TAXONOMY_ID; a sample's own items can name another
# one (a multi-species study: ST003799 cultures of several bacteria under one SUBJECT_SPECIES, ST005015 human
# and rat samples in a mouse study). Items: mwTab factors and additional sample data whose key is about the
# organism of the sample itself (_OWN_ORGANISM_KEY_RE: Species, Organism, Strain, Taxon, Genus, Bacteria, ...),
# not about another organism (_NOT_OWN_ORGANISM_KEY_RE: 'Infection', 'Colonized bacteria', 'Inoculum', 'Host',
# 'Treatment', 'Diet', co-culture partners, cultivars, 'mouse strain', ...). A value names a taxon when it is an
# allowed NCBI scientific name (case-insensitive) or a names.dmp synonym of one taxon (pinned taxonomy:
# 'Clostridium bolteae' -> Enterocloster bolteae), else (two words at least) the value with its white space
# normalised, without a strain designation at its end or with one wrong character fixed (mwb_organism_fallback:
# 'Bacteroides finegoldii DSMZ 17565', 'Muribaculun intestinale'); a one-word value only when written as NCBI writes
# the name ('Rattus', not the epithet 'melanogaster', which is a fungus genus), never under strain keys ('Strain:
# Valencia'); no virus, no placeholder taxon. Decided per analysis (apply_own_organisms), only when the
# analysis' species resolves (else the older per-sample lookup of translate_MWB_to_REDU_from_csv runs):
#   * a genus named by a one-word value that is the first word of another value of the same sample, whose taxon NCBI
#     has put in another genus since ('genus: Bacteroides' with 'taxonomy: Bacteroides plebeius DSMZ 17135', now
#     Phocaeicola plebeius: ST001688), names the same organism and is not counted as a taxon of its own;
#   * own taxa of another kind than the analysis' taxon (animal / plant / virus / microbe) are not used: a host
#     whose items name a microbe (an infection, its microbiota) keeps the host;
#   * a sample's own taxon above or equal to the analysis' taxon keeps the analysis' (more specific) one, one
#     below it (a strain of the species, a species of the genus) is used;
#   * the analysis is multi-species when its samples name two taxa of that kind that are not on one lineage;
#     then every sample with its own taxon gets it, and a sample without one (no value, a value that names no
#     taxon or several) gets no organism (the analysis' species is not the one of every sample: media controls,
#     co-cultures 'S. aureus and C. striatum', unresolved names), unless all the own taxa lie below the analysis'
#     taxon (ST001688 'Bacteria' over the genera of a strain collection; 'Enterobacter sp.' read as its genus):
#     then it holds for every sample and the sample keeps it;
#   * otherwise a sample's own taxon on another lineage is used, except for microbes (bacteria, archaea, fungi,
#     protists): the analysis' TAXONOMY_ID wins over a name (the user's decision: microbiology curators pick IDs
#     with care).
#   * for animals and plants the own taxon must lie in the analysis taxon's NCBI class ('Species: Schistosoma
#     mansoni' in a mouse study is a parasite; 'Genus: Major' an insect genus);
#   * pooled QCs, blanks and reference materials keep the analysis' organism (pooled QCs keep it: user decision).
#   * a value that lists organisms ('S. aureus and C. striatum': REDU_conversion_functions.ORGANISM_LIST) and names no
#     taxon as a whole names a sample of several organisms when each part names a taxon as above and the parts name two
#     taxa at least (own_taxa_term: 'ID|Name||ID|Name', ncbi_database.MULTI_SEP). It is used only in a multi-species
#     analysis, when every organism is of the analysis' kind (and class) and one of them is another sample's own taxon
#     (the co-culture of two of the analysis' monocultures); otherwise the sample counts as one without an own taxon.
#     A sample of several organisms is a value of its own (never read as its first organism).
# SampleType / SampleTypeSub1 follow the new taxon where they came from the analysis' taxon (several organisms:
# multi_sample_type).
_OWN_ORGANISM_KEY_RE = re.compile(
    r'(?:^|[^a-z])(?:species|organisms?|strains?|taxon|taxonomy|genus|isolates?|serovar|serotype|bacteria|bacterium|'
    r'microorganisms?|microbes?)(?:[^a-z]|$)')
_OWN_STRAIN_KEY_RE = re.compile(r'(?:^|[^a-z])(?:strains?|isolates?|serovar|serotype)(?:[^a-z]|$)')
_NOT_OWN_ORGANISM_KEY_RE = re.compile(
    r'treat|infect|coloni[sz]|inocul|diet|feed|\bfed\b|supplement|co ?cultur|partner|host|donor|recipient|gavage|'
    r'probiotic|prebiotic|expos|challeng|pathogen|parasit|virus|viral|phage|symbion|transplant|\bwith\b|cultivar|'
    r'variet|breed|\bpart\b|antibod|vector|plasmid|knock|mutant|genotype|background|reference|control|media\b|medium|'
    r'source|origin|prey|food|target|competitor|mouse|mice|\brats?\b|human|animal|plant|\bcells?\b|number|count|'
    r'\bid\b|code|abundance|percent|ratio|\bod\b|density|load')
_NO_ORGANISM_VALUES = {'', '-', '--', 'na', 'n/a', 'nan', 'none', 'null', 'not applicable', 'not available', 'unknown',
                       'not specified', 'unspecified', 'missing', 'missing value', 'nd', 'n.a.', 'blank', 'control',
                       'qc', 'pool', 'pooled', 'mixed', 'mix', 'all'}
# NCBI taxa whose descendants form a kind for apply_own_organisms ('microbe': every other cellular organism)
_ORGANISM_KIND_TAXA = {33208: 'animal', 33090: 'plant', 10239: 'virus'}
OWN_TAXA_SEP = ';'


def _taxon_node(taxid):
    pinned = _rcf._PINNED_LINEAGE
    if pinned is None or is_multi(taxid):
        # (a value of several organisms is no one taxon: never its first ID)
        return None
    try:
        return pinned.taxonomy.resolve(int(str(taxid).split('|')[0]))
    except (TypeError, ValueError):
        return None


def taxon_kind(taxid):
    """'animal' (Metazoa), 'plant' (Viridiplantae), 'virus' or 'microbe' (any other cellular organism:
    bacteria, archaea, fungi, protists) of an NCBI taxon ID or 'ID|Name' term (pinned taxonomy), else None."""
    node = _taxon_node(taxid)
    if node is None:
        return None
    parent = _rcf._PINNED_LINEAGE.taxonomy.parent
    for _ in range(1000):
        if node in _ORGANISM_KIND_TAXA:
            return _ORGANISM_KIND_TAXA[node]
        if node == 131567:   # cellular organisms
            return 'microbe'
        up = parent[node]
        if up == node or up < 0:
            return None
        node = up
    return None


def taxon_is_ancestor(ancestor, taxon):
    """Whether NCBI taxon `ancestor` is a proper ancestor of `taxon` (IDs or 'ID|Name' terms; pinned taxonomy)."""
    a, t = _taxon_node(ancestor), _taxon_node(taxon)
    if a is None or t is None or a == t:
        return False
    parent = _rcf._PINNED_LINEAGE.taxonomy.parent
    node = t
    for _ in range(1000):
        up = parent[node]
        if up == node or up < 0:
            return False
        if up == a:
            return True
        node = up
    return False


_CLASS_COLUMN = 3   # ncbi_database.DESIRED_RANKS: superkingdom, kingdom, phylum, class, ...


def taxon_class(taxid):
    """The NCBI class (taxon ID) a taxon ID / 'ID|Name' term lies in or is, else None (pinned taxonomy)."""
    node = _taxon_node(taxid)
    if node is None:
        return None
    tax = _rcf._PINNED_LINEAGE.taxonomy
    for _ in range(1000):
        code = tax.rank[node]
        if code >= 0 and tax.column[code] == _CLASS_COLUMN:
            return node
        up = tax.parent[node]
        if up == node or up < 0:
            return None
        node = up
    return None


def _qc_like_rows(df):
    """rows of pooled QCs, blanks and reference materials (the converter's own sample type, or the file name)"""
    st = df['SampleType_inner'].astype(str) if 'SampleType_inner' in df.columns else pd.Series('', index=df.index)
    qc = st.isin(['pool_QC', 'blank_QC', 'reference material']) | st.str.contains('blank', case=False)
    if 'filename' in df.columns:
        qc |= df['filename'].map(lambda f: bool(f) and isinstance(f, str) and blank_qc_type_from_filename(f) is not None)
    return qc.fillna(False).astype(bool)


def _same_taxon(a, b):
    return _taxon_node(a) is not None and _taxon_node(a) == _taxon_node(b)


def is_own_organism_key(key):
    k = _norm_key(key)
    return bool(k) and bool(_OWN_ORGANISM_KEY_RE.search(k)) and not _NOT_OWN_ORGANISM_KEY_RE.search(k) \
        and not is_variant_raw_file_key(str(key).lower()) and str(key).lower() not in _EXACT_RAW_KEYS


def own_taxon_term(key, value, allowedTerm_dict):
    """The allowed 'ID|Name' NCBITaxonomy term the value of an own-organism item names (see the block comment
    above), else None."""
    v = ' '.join(str(value).split())
    if v.lower() in _NO_ORGANISM_VALUES or not re.search(r'[A-Za-z]{3}', v):
        return None
    term = _ncbi_name_lookup(allowedTerm_dict).get(v.lower())
    if term is None and len(v.split()) > 1:
        term = _rcf.pinned_ncbi_synonym_term(v, allowedTerm_dict)
    if term is None and len(v.split()) > 1:
        # other white space, a strain designation at the end, one wrong character (two words at least)
        term = mwb_organism_fallback(v, allowedTerm_dict, f'sample item {key!r}')
    # 'Genus sp. <a lab's numbering>' ('Bacillus sp. 1') names the genus only (_rcf.sp_name_term)
    term = _sp_name_term(v, term, '', allowedTerm_dict)
    if term is None:
        return None
    if len(v.split()) == 1 and (_OWN_STRAIN_KEY_RE.search(_norm_key(key)) or term.split('|', 1)[1] != v):
        return None
    node = _taxon_node(term)
    if node is None or node in _rcf._NCBI_PLACEHOLDER_TAXA or taxon_kind(term) in (None, 'virus'):
        return None
    return term


def own_taxa_term(key, value, allowedTerm_dict):
    """own_taxon_term of an own-organism item's value; for a value that lists organisms (ORGANISM_LIST) and names no
    taxon itself, the 'ID|Name||ID|Name' term of its parts when each names a taxon (own_taxon_term, the same key), each
    taxon once in the value's order (one taxon left: that term); else None."""
    term = own_taxon_term(key, value, allowedTerm_dict)
    if term is not None:
        return term
    v = ' '.join(str(value).split())
    if not ORGANISM_LIST.search(v):
        return None
    parts = [p.strip() for p in ORGANISM_LIST.split(v) if p.strip()]
    if len(parts) < 2:
        return None
    terms = [own_taxon_term(key, p, allowedTerm_dict) for p in parts]
    if any(t is None for t in terms):
        return None
    # each taxon once: by its pinned node (merged IDs followed), by its ID where it has none
    return join_parts(terms, key=lambda t: _taxon_node(t) or t.split('|', 1)[0].strip())


def sample_own_taxa(df_items, allowedTerm_dict):
    """Sample ID -> the distinct taxa its own-organism items name (own_taxa_term), joined by OWN_TAXA_SEP in
    sorted order; df_items: the mwTab items (columns filename = Sample ID, Key, Value as given)."""
    out = {}
    if df_items is None or not len(df_items) or _rcf._PINNED_LINEAGE is None:
        return out
    keys = df_items['Key'].astype(str)
    use = keys.map(is_own_organism_key)
    memo = {}
    named = {}
    for sid, k, v in zip(df_items.loc[use, 'filename'], keys[use], df_items.loc[use, 'Value']):
        if (k, v) not in memo:
            memo[(k, v)] = own_taxa_term(k, v, allowedTerm_dict)
        t = memo[(k, v)]
        if t is not None:
            out.setdefault(str(sid), set()).add(t)
            named.setdefault(str(sid), []).append((' '.join(str(v).lower().split()), t))
        else:
            out.setdefault(str(sid), set())
    # a genus written as the first word of another value of the sample that names a taxon NCBI has put in another
    # genus since ('genus: Bacteroides' with 'taxonomy: Bacteroides plebeius DSMZ 17135', now Phocaeicola
    # plebeius) names the same organism: the genus is not counted as a taxon of its own
    for sid, pairs in named.items():
        for v1, t1 in pairs:
            if ' ' in v1:
                continue
            if any(v2.split(' ', 1)[0] == v1 and ' ' in v2 and not (_same_taxon(t1, t2) or taxon_is_ancestor(t1, t2)
                                                                      or taxon_is_ancestor(t2, t1))
                   for v2, t2 in pairs):
                out[sid].discard(t1)
    return {s: OWN_TAXA_SEP.join(sorted(ts)) for s, ts in out.items()}


def _one_lineage(terms):
    """The deepest of terms that all lie on one lineage, else None."""
    terms = sorted(set(terms))
    if not terms:
        return None
    deepest = terms[0]
    for t in terms[1:]:
        if _same_taxon(t, deepest) or taxon_is_ancestor(t, deepest):
            continue
        if taxon_is_ancestor(deepest, t):
            deepest = t
            continue
        return None
    return deepest


def apply_own_organisms(df, study_term, column='MWB_own_taxa'):
    """NCBITaxonomy (and SampleType / SampleTypeSub1 where they came from the analysis' taxon) of the rows of one
    analysis from the samples' own taxa (column `column`, from sample_own_taxa), by the rules of the block
    comment above. Returns the frame and the number of rows changed."""
    if column not in df.columns or not isinstance(study_term, str) or '|' not in study_term \
            or _rcf._PINNED_LINEAGE is None:
        return df, 0
    kind = taxon_kind(study_term)
    if kind is None:
        return df, 0
    own = df[column].fillna('').astype(str)
    # pooled QCs, blanks and reference materials keep the analysis' organism (pooled QCs keep it: the user's
    # decision; blanks lose it later): their own values are neither used nor counted
    exempt = _qc_like_rows(df)
    own = own.where(~exempt, '')
    study_class = taxon_class(study_term) if kind in ('animal', 'plant') else None

    def usable(t):
        # an own taxon of the analysis' kind; for animals and plants also of its NCBI class (a parasite named in a
        # mouse study, 'Species: Schistosoma mansoni', does not replace the host; nor 'Genus: Major', an insect)
        if is_multi(t):
            # several organisms: every one of them
            return all(usable(p) for p in taxonomy_parts(t))
        if taxon_kind(t) != kind:
            return False
        if kind not in ('animal', 'plant'):
            return True
        if study_class is None:
            return taxon_is_ancestor(study_term, t) or _same_taxon(t, study_term)
        return taxon_class(t) == study_class
    per_value = {}
    for value in own.unique():
        terms = [t for t in value.split(OWN_TAXA_SEP) if t and usable(t)]
        per_value[value] = (_one_lineage(terms) if terms else None, bool(terms))
    # (the samples of several organisms do not make the analysis multi-species)
    taxa = {t for t, _ in per_value.values() if t and not is_multi(t)}
    multi = any(not (_same_taxon(a, b) or taxon_is_ancestor(a, b) or taxon_is_ancestor(b, a))
                for a in taxa for b in taxa if a < b)
    own_nodes = {_taxon_node(t) for t in taxa} - {None}
    # the analysis' taxon holds for every sample when all own taxa lie below it ('Bacteria' over the genera of a
    # strain collection; an 'X sp.' taxon stands for its genus X here); then a sample without an own taxon keeps it
    anchor = _taxon_node(study_term)
    sp = re.fullmatch(r'([A-Z][a-z]+) spp?\.?', study_term.split('|', 1)[-1].strip())
    if anchor is not None and sp:
        node, parent, names = anchor, _rcf._PINNED_LINEAGE.taxonomy.parent, _rcf._PINNED_LINEAGE.names
        for _ in range(50):
            node = int(parent[node])
            if node < 0 or names.get(node) == sp.group(1):
                break
        anchor = node if node >= 0 and names.get(node) == sp.group(1) else anchor
    all_below = anchor is not None and all(_taxon_node(t) == anchor or taxon_is_ancestor(anchor, t) for t in taxa)
    study_type = get_taxonomy_info(study_term.split('|')[0])
    new = {}
    for value, (t, named) in per_value.items():
        if t is not None and is_multi(t):
            # several organisms: in a multi-species analysis, one of them another sample's own taxon
            if multi and any(_taxon_node(p) in own_nodes for p in taxonomy_parts(t)):
                new[value] = t
                continue
            t = None
        if t is not None:
            if _same_taxon(t, study_term) or taxon_is_ancestor(t, study_term):
                new[value] = study_term
            elif taxon_is_ancestor(study_term, t) or multi or kind != 'microbe':
                new[value] = t
            else:
                new[value] = study_term
        elif multi and not all_below:
            new[value] = None
        else:
            new[value] = study_term
    target = own.map(new).where(~exempt, study_term)
    changed = target.ne(study_term)
    if not changed.any():
        return df, 0
    df = df.copy()
    df.loc[changed, 'NCBITaxonomy'] = target[changed]
    if study_type and study_type[0] is not None and {'SampleType', 'SampleTypeSub1'} <= set(df.columns):
        sub1 = df['SampleTypeSub1']
        sub1_missing = sub1.isna() | sub1.astype(str).isin(['', 'nan', 'None', 'missing value'])
        from_taxon = df['SampleType'].astype(str).eq(str(study_type[0])) & (
            sub1.astype(str).eq(str(study_type[1])) if study_type[1] is not None else sub1_missing)
        for t in sorted(set(target[changed].dropna())):
            info = multi_sample_type(taxonomy_parts(t)) if is_multi(t) else get_taxonomy_info(t.split('|')[0])
            if not info or info[0] is None:
                continue
            rows = changed & from_taxon & target.eq(t)
            df.loc[rows, 'SampleType'] = info[0]
            df.loc[rows, 'SampleTypeSub1'] = info[1] if info[1] is not None else np.nan
    return df, int(changed.sum())


# ---- the sample's own body part ---------------------------------------------------------------------------
# UBERONBodyPartName came from the REST record's sample_source, which in many studies is one value for every
# sample (the study's sample type: ST002254 'Liver' on liver, plasma and muscle samples; ST001808 'Liver' on
# serum samples), and where that did not translate from the analysis' SAMPLE_TYPE (ST005015 'BAT' / 'eWAT' samples
# labelled pancreas). sample_body_parts() decides per sample from its own values, each translated by
# body_part_of_value (UBERONBodyPartName.csv, else an exact UBERON/CL/PO label or synonym (no taxon-specific class),
# else _BODY_PART_ABBREVIATIONS as written; also with '_' read as a space, and the sheet tissue left without species
# words; only allowed labels, generic ones such as tissue, organ, cell are no body part here):
#   own: the first of P1 the mwTab 'Sample source' item of the sample, P3 its items under body-part keys
#   (_OWN_BODY_PART_KEYS: Organ, Tissue, Sample type, Matrix, Specimen, ...; all of them agreeing), P4 its Subject ID
#   if that is exactly a body part ('serum'; only while the REST sample_source has one value for all of the
#   analysis' samples: a study-level value) that names one body part; REST: the REST sample_source's body part.
#   The own one wins over the REST one where its level tells the analysis' samples apart (two or more values over
#   the samples: per-sample information; ST001467 'Sample Source: muscle / plasma' under one REST 'Blood (plasma)')
#   or where it is the REST one made more specific (its words include the REST label's: 'blood plasma' over
#   'blood'); otherwise (a study-level own value against another study-level REST value: ST003112 'Brain' against
#   REST 'Brain cortex') the REST one stays, as before. An own value never replaces a more specific REST value, nor,
#   where the REST record gives none (pooled QCs 'Sample source: blood' under '_QC_'), the study-level SAMPLE_TYPE's
#   (broader / more specific by _broader_body_part's refinement table: the table of the run has no is_a / part_of).
# A REST value that differs between the analysis' samples is the sample's own record: an own value only refines it.
# The study-level SAMPLE_TYPE no longer fills a sample of a multi-tissue analysis (a sample's own tissue is neither
# the SAMPLE_TYPE's nor a refinement of it, and other samples name the SAMPLE_TYPE's tissue themselves) whose own body
# part is unknown (its source does not translate, or it has none): left missing rather than given another tissue. QC / blank / standard sample sources (_QC_, _POOLED_, _BLANK_, _STD_) keep the
# earlier handling (pooled QCs keep the study's body part).
_OWN_BODY_PART_KEYS = {'sample matrix', 'client matrix', 'matrix', 'material', 'sample material', 'biological material',
                       'specimen type', 'specimen', 'tissue', 'tissues', 'tissue type', 'tissuetype', 'tissue location',
                       'tissue region', 'organ', 'organs', 'organism part', 'body site', 'anatomical site', 'body part',
                       'sample type', 'sampletype', 'sample type matrix', 'tissue fluid', 'tissue serum', 'biofluid',
                       'biospecimen', 'biospecimen type', 'plant organ', 'plant organs', 'plant part', 'plant tissue',
                       'sample tissue', 'tissue sample', 'organ tissue', 'tissue organ'}
_OWN_SAMPLE_SOURCE_KEYS = {'sample source', 'samplesource'}
_QC_SOURCE_MARKERS = {'_qc_', '_pooled_', '_blank_', STANDARD_SAMPLE_SOURCE}
# common tissue abbreviations, written as they are written (case-sensitive: 'Bat' / 'bat' is the animal), that name
# one body part (eWAT: epididymal or epigonadal white fat, both a gonadal fat pad; gWAT gonadal; iWAT inguinal; sWAT /
# scWAT subcutaneous; BAT / iBAT brown, interscapular brown), and short names (any case)
_BODY_PART_ABBREVIATIONS = {'BAT': 'brown adipose tissue', 'iBAT': 'brown adipose tissue', 'WAT': 'white adipose tissue',
                            'eWAT': 'gonadal fat pad', 'gWAT': 'gonadal fat pad', 'iWAT': 'inguinal fat pad',
                            'sWAT': 'subcutaneous adipose tissue', 'scWAT': 'subcutaneous adipose tissue'}
_BODY_PART_SHORT_NAMES = {'adipose': 'adipose tissue'}
# species words around a sheet tissue ('human_liver', 'mouse spinal cord') and the NCBI taxon each names (the taxon or
# one below it): the rest is looked up in the sheet only, and only where every such word names the row's organism
# (the analysis' or the sample's own; 'Pooled_human_plasma' in a mouse study is not translated)
_BODY_PART_SPECIES_WORDS = {'human': 9606, 'humans': 9606, 'mouse': 10088, 'mice': 10088, 'murine': 10088, 'rat': 10114,
                            'rats': 10114, 'porcine': 9823, 'pig': 9823, 'bovine': 9903, 'cow': 9913, 'canine': 9611,
                            'dog': 9615, 'feline': 9681, 'cat': 9685, 'ovine': 9935, 'sheep': 9940, 'caprine': 9922,
                            'goat': 9925, 'equine': 9789, 'horse': 9796, 'rabbit': 9986, 'hamster': 10026,
                            'macaque': 9539, 'baboon': 9554, 'chicken': 9031, 'zebrafish': 7955}
# (with _UBERON_TOO_GENERIC) labels that name no particular body part
_GENERIC_BODY_PARTS = {'tissue', 'organ', 'cell', 'multicellular organism', 'organism substance', 'anatomical entity',
                       'material anatomical entity', 'adult organism', 'whole organism', 'post-embryonic organism',
                       'cell cluster', 'not applicable', 'not specified', 'ml import: not available'}
# body part -> its refinements (subtypes and parts a sample of it is often given as), for _broader_body_part: the
# UBERON/CL/PO table of the run has no is_a / part_of links, and label words alone mistake 'blood vessel' for a kind of
# blood. Dataset-agnostic; only allowed labels count.
_ADIPOSE_DEPOTS = {'subcutaneous adipose tissue', 'subcutaneous abdominal adipose tissue', 'visceral abdominal adipose tissue',
                   'gonadal fat pad', 'epididymal fat pad', 'parametrial fat pad', 'periovarian fat pad', 'inguinal fat pad',
                   'retroperitoneal fat pad', 'mesenteric fat pad', 'abdominal fat pad', 'perirenal fat', 'omental fat pad'}
_BODY_PART_REFINEMENTS = {
    'blood': {'blood plasma', 'blood serum', 'umbilical cord blood', 'venous blood', 'arterial blood'},
    'adipose tissue': {'brown adipose tissue', 'white adipose tissue', 'interscapular fat pad', 'mammary gland fat',
                       'epicardial fat'} | _ADIPOSE_DEPOTS,
    'white adipose tissue': _ADIPOSE_DEPOTS,
    'gonadal fat pad': {'epididymal fat pad', 'parametrial fat pad', 'periovarian fat pad'},
    'bone tissue': {'compact bone tissue', 'spongy bone tissue'},
    'lung': {'left lung', 'right lung'},
    'skeletal muscle organ, vertebrate': {'gastrocnemius', 'soleus muscle', 'quadriceps femoris', 'tibialis anterior',
                                          'vastus lateralis', 'extensor digitorum longus', 'plantaris'},
    'zone of skin': {'skin of back', 'skin of body', 'skin of leg', 'skin of trunk', 'arm skin', 'axilla skin', 'scalp',
                     'skin of manus', 'skin of pes', 'head or neck skin', 'skin of prepuce of penis'},
}


def _is_generic_body_part(label):
    low = str(label).strip().lower()
    return low in _GENERIC_BODY_PARTS or low in _UBERON_TOO_GENERIC


def _body_part_lookups(path_to_csvs, ontology_table, allowedTerm_dict):
    """(sheet, exact names, abbreviations, short names): value -> allowed specific UBERONBodyPartName label (keys
    lower-cased, abbreviations as written). A sheet value whose label is no allowed one or a generic one gets none
    here (False: the older sheet path keeps it); a sheet value with an empty label (None: 'BAT') is looked up by its
    exact name and abbreviation."""
    allowed = set((allowedTerm_dict.get('UBERONBodyPartName') or {}).get('allowed_values') or [])
    exact = _uberon_exact_lookup(ontology_table, allowed) if ontology_table is not None and 'Synonym' in ontology_table.columns else {}
    keep = lambda l: l if isinstance(l, str) and l in allowed and not _is_generic_body_part(l) else False
    sheet = _read_sheet(os.path.join(path_to_csvs, 'UBERONBodyPartName.csv'), encoding="ISO-8859-1", dtype=str)
    sheet_map = {}
    for k, v in zip(sheet['MWB'].astype(str).str.lower(), sheet['REDU']):
        label = v.strip() if isinstance(v, str) and v.strip() else None
        sheet_map.setdefault(' '.join(k.split()), None if label is None else keep(label))
    exact = {k: keep(v) for k, v in exact.items()}
    abbreviations = {k: keep(v) for k, v in _BODY_PART_ABBREVIATIONS.items()}
    short = {k: keep(v) for k, v in _BODY_PART_SHORT_NAMES.items()}
    return sheet_map, exact, abbreviations, short


def _species_word_fits(word, organisms):
    """Whether a species word names one of the organisms (NCBI terms / IDs: the taxon or one below it)"""
    taxon = _BODY_PART_SPECIES_WORDS[word]
    return any(_same_taxon(o, taxon) or taxon_is_ancestor(taxon, o) for o in organisms if o)


def body_part_of_value(value, lookups, organisms=()):
    """The allowed specific UBERONBodyPartName label a sample's value names (sheet, exact UBERON/CL/PO name,
    abbreviation as written, short name; also with '_' read as a space; else the sheet tissue left once species
    words that name the row's organism (organisms: NCBI terms) are taken out: 'human_liver' -> liver in a human
    sample), else None."""
    sheet_map, exact, abbreviations, short = lookups[:4]
    raw = ' '.join(str(value).split())
    low = raw.lower()
    if not low or low in _UBERON_PLACEHOLDERS or low in _QC_SOURCE_MARKERS:
        return None
    for cand, cand_raw in dict.fromkeys([(low, raw), (' '.join(low.replace('_', ' ').split()),
                                                      ' '.join(raw.replace('_', ' ').split()))]):
        hit = sheet_map.get(cand)
        if hit is False:   # the sheet's label for it is generic or no allowed label
            return None
        hit = hit or exact.get(cand) or abbreviations.get(cand_raw) or short.get(cand)
        if hit:
            return hit
    words = [w for w in re.split(r'[\s_\-./]+', low) if w]
    species = [w for w in words if w in _BODY_PART_SPECIES_WORDS]
    rest = [w for w in words if w not in _BODY_PART_SPECIES_WORDS]
    if rest and species and all(_species_word_fits(w, organisms) for w in species):
        hit = sheet_map.get(' '.join(rest))
        if hit:
            return hit
    return None


def _sheet_lists(value, lookups):
    """Whether the translation sheet gives the value a label (as the older sheet path reads it: lower case, '_' as a
    space), allowed or not"""
    low = ' '.join(str(value).lower().split())
    return any(lookups[0].get(c) is not None for c in (low, ' '.join(low.replace('_', ' ').split())))


def _broader_body_part(a, b):
    """Whether body part label b is a refinement of label a (_BODY_PART_REFINEMENTS: 'blood' / 'blood serum',
    'adipose tissue' / 'brown adipose tissue')."""
    return isinstance(a, str) and isinstance(b, str) and b in _BODY_PART_REFINEMENTS.get(a, ())


def sample_body_parts(df_items, rest_sources, study_sample_type, lookups, study_taxon=None, own_taxa=None):
    """Sample ID -> (UBERONBodyPartName label of the sample's own values or None, whether the study-level
    SAMPLE_TYPE may not fill it), by the rules of the block comment above. df_items: the analysis' mwTab items
    (filename = Sample ID, SubjectIdentifierAsRecorded, Key, Value); rest_sources: Sample ID -> REST
    sample_source (the analysis' samples); study_taxon: the analysis' NCBI term; own_taxa: Sample ID -> the taxa
    the sample's own items name (sample_own_taxa), for the species words."""
    if df_items is None:
        df_items = pd.DataFrame(columns=['filename', 'SubjectIdentifierAsRecorded', 'Key', 'Value'])
    keyn = df_items['Key'].astype(str).map(_norm_key)
    own_taxa = own_taxa or {}
    memo = {}

    def organisms_of(sid):
        found = [t for t in str(own_taxa.get(sid, '') or '').split(OWN_TAXA_SEP) if t]
        return tuple(sorted(set(found) | ({study_taxon} if isinstance(study_taxon, str) and '|' in study_taxon else set())))

    def bp(v, sid):
        key = (v, organisms_of(sid))
        if key not in memo:
            memo[key] = body_part_of_value(v, lookups, key[1])
        return memo[key]
    levels = {}
    for sid, k, v in zip(df_items['filename'].astype(str), keyn, df_items['Value']):
        if k in _OWN_SAMPLE_SOURCE_KEYS:
            levels.setdefault(sid, {}).setdefault('P1', set()).add(str(v))
        elif k in _OWN_BODY_PART_KEYS:
            levels.setdefault(sid, {}).setdefault('P3', set()).add(str(v))
    if 'SubjectIdentifierAsRecorded' in df_items.columns:
        for sid, subj in zip(df_items['filename'].astype(str), df_items['SubjectIdentifierAsRecorded']):
            if isinstance(subj, str) and subj.strip() and subj.strip() != '-':
                levels.setdefault(sid, {}).setdefault('P4', set()).add(subj.strip())
    rest_sources = {str(s): v for s, v in rest_sources.items()}
    marked = {sid for sid, src in rest_sources.items() if str(src).strip().lower() in _QC_SOURCE_MARKERS}
    real_rest = {str(src).strip().lower() for sid, src in rest_sources.items()
                 if sid not in marked and str(src).strip().lower() not in _UBERON_PLACEHOLDERS}
    rest_study_level = len(real_rest) <= 1
    # levels whose values tell the analysis' samples apart: P1 / P3 by their values, P4 by its body parts
    told_apart = {}
    for level in ('P1', 'P3', 'P4'):
        vals = set()
        for sid, lv in levels.items():
            if sid in marked:
                continue
            for v in lv.get(level, ()):
                low = ' '.join(str(v).lower().split())
                if level == 'P4':
                    low = bp(v, sid)
                if low and low not in _UBERON_PLACEHOLDERS:
                    vals.add(low)
        told_apart[level] = len(vals) >= 2
    out = {}
    for sid in set(levels) | set(rest_sources):
        lv = levels.get(sid, {})
        rest = bp(rest_sources[sid], sid) if isinstance(rest_sources.get(sid), str) else None
        own = own_level = None
        for level in ('P1', 'P3', 'P4'):
            if level == 'P4' and not rest_study_level:
                break
            labels = {bp(v, sid) for v in lv.get(level, ())} - {None}
            if len(labels) == 1:
                own, own_level = next(iter(labels)), level
                break
        if own is None or own == rest:
            out[sid] = rest
        elif rest is None or _broader_body_part(rest, own):
            out[sid] = own    # no REST body part, or the own value refines it ('blood serum' / REST 'blood')
        elif _broader_body_part(own, rest) or not rest_study_level:
            # the own value is a broader form of the REST one, or the REST value differs between the samples (it is
            # the sample's own record: ST003077 pools 'Heart' / 'Kidney' with a copied 'Sample type: Plasma')
            out[sid] = rest
        elif told_apart[own_level]:
            out[sid] = own
        else:
            out[sid] = rest
    study_bp = body_part_of_value(study_sample_type, lookups,
                                  (study_taxon,) if isinstance(study_taxon, str) else ()) \
        if isinstance(study_sample_type, str) else None
    # a value that is a broader form of the study-level SAMPLE_TYPE's tissue where the translation sheet gives the REST
    # record none (pooled QCs 'Sample source: blood' under '_QC_'; REST 'Mouse blood', read by the species-word rule, in
    # a 'Blood (serum)' study): the study-level fill stays, as before
    broader_than_study = set()
    if study_bp is not None:
        for sid, l in out.items():
            rest = rest_sources.get(sid)
            if l is not None and _broader_body_part(l, study_bp) \
                    and (not isinstance(rest, str) or not _sheet_lists(rest, lookups)):
                broader_than_study.add(sid)
    for sid in broader_than_study:
        out[sid] = None
    # the SAMPLE_TYPE does not fill: a multi-tissue analysis (a sample's own tissue is neither the study tissue nor a
    # refinement of it) where other samples name the study tissue (or a refinement) themselves; where none does, the
    # SAMPLE_TYPE describes the samples without a translated value (ST003472 'Ventricle' in a 'Heart' study with
    # plasma samples) and fills them
    named = {l for s, l in out.items() if l and s not in marked}
    multi = study_bp is not None and any(l != study_bp and not _broader_body_part(study_bp, l) for l in named) \
        and any(l == study_bp or _broader_body_part(study_bp, l) for l in named)
    return {sid: (l, bool(multi and l is None and sid not in marked and sid not in broader_than_study))
            for sid, l in out.items()}


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
        # the archive names ride along as a '_ctx' helper column (dropped with the others below)
        raw_file_name_df = raw_file_name_df.rename(columns={'filename': 'filename_raw_path', 'archive': '_ctx_archive'})


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

        # keep one row per stem and acquisition context (a stem's files in folders of
        # another polarity / chromatography are other acquisitions, not copies; files
        # whose path contradicts this analysis are left out where the stem has others):
        # best rank, tie-break by name for determinism
        tmp["_ctx"] = acquisition_contexts(tmp["_stem"].tolist(), tmp["filename_raw_path"].astype(str).tolist(), hints,
                                           tmp["_rank"].tolist(),
                                           tmp["_ctx_archive"].tolist() if "_ctx_archive" in tmp.columns else None)
        # same-name files of different samples, told apart by the folders of the listed raw file names
        _raw_key = df['Key'].isin(expected_raw_file_keys) | df['Key'].map(is_variant_raw_file_key)
        _listed_paths = sorted({n for v in df.loc[_raw_key, 'Value'].dropna().astype(str)
                                for n in split_raw_file_names(v, False) if '/' in n or '\\' in n})
        if _listed_paths:
            _before = tmp["_ctx"].tolist()
            tmp["_ctx"] = contexts_told_apart_by_listed_folders(tmp["_stem"].tolist(), tmp["filename_raw_path"].astype(str).tolist(),
                                                                _before, _listed_paths)
            for _u, _a, _b in zip(tmp['USI'], _before, tmp["_ctx"]):
                if _a >= 0 and _b < 0:
                    note_withheld(_u)
            _n_apart = sum(1 for a, b in zip(_before, tmp["_ctx"]) if a != b)
            if _n_apart:
                print(f"[rawfile pref] {_n_apart} same-name files told apart by the folders of listed raw file names")
        # (polarity, chromatography) of every listed file, with the study's archive names
        _file_contexts = listing_contexts(tmp, 'filename_raw_path', '_ctx_archive')
        n_ctx_out = int((tmp["_ctx"] < 0).sum())
        n_ctx_split = int(tmp.loc[tmp["_ctx"] > 0, "_stem"].nunique())
        if n_ctx_out or n_ctx_split:
            print(f"[rawfile pref] acquisition contexts: {n_ctx_split} stems kept per polarity/chromatography folder, "
                  f"{n_ctx_out} files of another polarity/chromatography left out")
        tmp = tmp[tmp["_ctx"] >= 0]
        tmp = tmp.sort_values(["_stem", "_ctx", "_rank", "_base_lower"]).drop_duplicates(["_stem", "_ctx"], keep="first")

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
        _split_at_spaces = not raw_file_name_df['filename_base'].str.contains(" ").any()
        df = df.assign(filename_raw=df['filename_raw'].map(lambda v: split_raw_file_names(v, _split_at_spaces))).explode('filename_raw')

        df.reset_index(drop=True, inplace=True)

        # a raw file several samples name, uploaded also as one file per sample
        df['filename_raw'] = per_sample_files_of_shared_names(
            df, raw_file_name_df, hints=hints,
            sample_values=lambda: zip(pd.concat([df['_sample_id'], df_additional['filename']], ignore_index=True),
                                      pd.concat([df['Key'], df_additional['Key']], ignore_index=True),
                                      pd.concat([df['Value'], df_additional['Value']], ignore_index=True)),
            contexts=_file_contexts)


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

        # a stem kept in several acquisition contexts (pos/x.mzML, neg/x.mzXML) gives each
        # of its file names; any other stem its one preferred name
        _pref_sorted = raw_file_name_df.sort_values(['_stem', '_rank', '_base_lower'])
        _multi_ctx = _pref_sorted.groupby('_stem')['_ctx'].transform('nunique') > 1
        preferred_raw_by_stem = (
            pd.concat([_pref_sorted[~_multi_ctx].drop_duplicates('_stem', keep='first'),
                       _pref_sorted[_multi_ctx].drop_duplicates(['_stem', 'filename_base'], keep='first')])
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
            df = match_unmatched_names_by_prefix(df, raw_file_name_df)
        elif bool(set(df['filename_raw_lower']) & set(raw_file_name_df['filename_base_lower'])):
            df = df.merge(raw_file_name_df, left_on='filename_raw_lower', right_on='filename_base_lower', how='left')
            df.drop(['filename_raw_lower', 'filename_lower', 'filename_base_lower'], axis=1, inplace=True)
            df = match_unmatched_names_by_prefix(df, raw_file_name_df)
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
                # one sample per file: a file named after a sample goes to that sample only, else to the one closest
                # whole-token match; ambiguous files to nobody (substring_file_matches; the longest listed name used
                # to win, which gave '32650.mzML' to sample 'L32650' and '1.mzML' to sample '11641')
                print('resolving contradictions,..')
                _ids_of = {}
                for _v, _s in zip(df['filename_base_wo_extension'], df['_sample_id']):
                    if pd.notna(_s):
                        _ids_of.setdefault(_v, set()).add(str(_s))
                # the listing rows used: the first row of each file stem per acquisition context (one stem kept in
                # pos and neg folders: both)
                _first_label = {}
                for _lbl, _stem_value, _c in zip(raw_file_name_df.index, raw_file_name_df['filename_base_wo_extension'],
                                                 raw_file_name_df['_ctx']):
                    _first_label.setdefault(_stem_value, {}).setdefault(_c, _lbl)
                _files = [(_lbl, _stem_value, raw_file_name_df.at[_lbl, 'filename_raw_path'])
                          for _stem_value in sorted(_first_label, key=str) for _lbl in _first_label[_stem_value].values()]
                # the listed names as listed (letter case: word boundaries 'BruinCSHneg')
                _listed = {}
                for _v, _r in zip(df['filename_base_wo_extension'], df['filename_raw'].astype(str)):
                    _listed.setdefault(_v, set()).add(_r.split('.')[0])
                _listed = {_v: sorted(_rs)[0] for _v, _rs in _listed.items()}
                _reasons = {}
                _pairs, _why = substring_file_matches(potential_matches, _ids_of, _files, _listed, reasons=_reasons)
                for _lbl, _r in _reasons.items():
                    if _r in _WITHHELD_REASONS:
                        note_withheld(raw_file_name_df.at[_lbl, 'USI'])
                print("[rawfile] substring matches: " + ', '.join(f"{r}: {_why[r]}" for r in SUBSTRING_MATCH_REASONS if _why[r]))

                df['match_index'] = [_pairs.get((_v, str(_s) if pd.notna(_s) else str(_v)))
                                     for _v, _s in zip(df['filename_base_wo_extension'], df['_sample_id'])]
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

            # stems kept in several acquisition contexts keep the file each row matched
            _multi_ctx = raw_file_name_df.groupby('_stem')['_ctx'].transform('nunique') > 1
            preferred_raw_by_stem = (
                raw_file_name_df[~_multi_ctx]
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



    df = df.drop(columns=[c for c in df.columns if str(c).startswith('_ctx')])
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
    # Subject ID and Sample ID as given, to tell pooled QC samples ('WB_SERUM_PLASMA_...')
    _id_text = (df['SubjectIdentifierAsRecorded'].fillna('').astype(str) + ' ' + df['filename'].astype(str)) \
        .groupby(df['filename'].astype(str)).agg(lambda s: ' '.join(sorted(set(s)))).str.lower()
    try:
        df['SubjectIdentifierAsRecorded'] = _drop_label_subject_ids(df)
    except Exception as e:
        print(f"Could not check Subject IDs: {type(e).__name__}: {e}")
        df['SubjectIdentifierAsRecorded'] = ''
    df = get_key_info_into_outer(df, key_vars=["Sample Source"], new_col="MWB_sampleSource")
    # the sample's own body part and whether the study-level SAMPLE_TYPE may fill it (sample_body_parts); a
    # problem here must never cost the study its rows
    try:
        _rest_sources = {}
        for _rec in (rest_response.values() if isinstance(rest_response, dict) else []):
            if isinstance(_rec, dict) and _rec.get('local_sample_id') is not None:
                _rest_sources.setdefault(str(_rec.get('local_sample_id')), _rec.get('sample_source'))
        _sids = set(df['filename'].astype(str))
        _rest_sources = {s: v for s, v in _rest_sources.items() if s in _sids}
        _bp = sample_body_parts(df_additional[df_additional['filename'].astype(str).isin(_sids)], _rest_sources,
                                _ANALYSIS_SAMPLE_TYPE, _body_part_lookups(path_to_csvs, ontology_table, allowedTerm_dict),
                                study_taxon=_ANALYSIS_TAXON, own_taxa=sample_own_taxa(df_additional, allowedTerm_dict))
        df['MWB_own_bodypart'] = df['filename'].astype(str).map(lambda s: _bp.get(s, (None, False))[0])
        df['MWB_no_study_bodypart'] = df['filename'].astype(str).map(lambda s: _bp.get(s, (None, False))[1])
        print(f"[own body part] {df.loc[df['MWB_own_bodypart'].notna(), 'filename'].nunique()} samples with their own "
              f"body part, {df.loc[df['MWB_no_study_bodypart'], 'filename'].nunique()} not filled from the study level")
    except Exception as e:
        print(f"Could not read the samples' own body parts: {type(e).__name__}: {e}")
        df['MWB_own_bodypart'] = None
        df['MWB_no_study_bodypart'] = False
    # the food or beverage product a sample is (FoodOntologyTerm); a problem here must never
    # cost the study its rows
    try:
        df['MWB_food'] = sample_food_terms(df)
    except Exception as e:
        print(f"Could not check for food samples: {type(e).__name__}: {e}")
        df['MWB_food'] = None
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
        _vals = [mwb_age_to_years(v, k.replace('_', ' '), _ANALYSIS_AGE_UNIT) if _SINGLE_AGE_VALUE_RE.fullmatch(str(v).strip()) else None
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
            if _ANALYSIS_ABBR:
                # study-defined abbreviations ('Diagnosis: COPD'), only where nothing was found
                _xvals = df['Value'].map(lambda v: expand_abbreviations(v, _ANALYSIS_ABBR)).astype(str).str.strip().str.lower()
                _xdoid = _xvals.where(_is_health & (_xvals != _vals) &
                                  ~df['filename'].isin(isolated_material_samples(df))).map(DOID_LOOKUP)
                _xdoid = _xdoid.groupby(df['filename']).transform(lambda s: s.dropna().iloc[0] if s.dropna().nunique() == 1 else None)
                df['MWB_doid'] = df['MWB_doid'].where(df['MWB_doid'].notna(), _xdoid)
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
    # time point, fasting state and the study's own sample group; a problem here must
    # never cost the study its rows
    for _col in ['MWB_timepoint', 'MWB_fasting', 'MWB_samplegroup']:
        df[_col] = None
    try:
        df['MWB_timepoint'] = sample_timepoint_hours(df, _STUDY_TEXT_REST)
    except Exception as e:
        print(f"Could not extract time points: {type(e).__name__}: {e}")
        df['MWB_timepoint'] = None
    try:
        df['MWB_fasting'] = sample_fasting_status(df)
    except Exception as e:
        print(f"Could not extract fasting states: {type(e).__name__}: {e}")
        df['MWB_fasting'] = None
    # cell lines and cell types named per sample (resolved in translate_MWB_to_REDU_by_logic)
    try:
        df['MWB_cellline_values'], df['MWB_celltype'] = sample_cell_values(df, ontology_table, allowedTerm_dict)
    except Exception as e:
        print(f"Could not extract cell lines: {type(e).__name__}: {e}")
        df['MWB_cellline_values'] = None
        df['MWB_celltype'] = None
    try:
        _groups = sample_group_texts(data)
        if _groups:
            df['MWB_samplegroup'] = df['filename'].map(_groups)
    except Exception as e:
        print(f"Could not build sample groups: {type(e).__name__}: {e}")
        df['MWB_samplegroup'] = None
    # study-level disease for the samples of case groups and healthy-control groups
    # (filled in translate_MWB_to_REDU_by_logic where nothing per-sample was found)
    df['MWB_study_health'] = None
    df['MWB_study_doid'] = None
    # the factor values with the study's own abbreviations written out (None if none apply)
    _df_x = None
    if _ANALYSIS_ABBR:
        try:
            _df_x = df.copy()
            _df_x['Value'] = df['Value'].map(lambda v: expand_abbreviations(v, _ANALYSIS_ABBR))
            if (_df_x['Value'].astype(str) == df['Value'].astype(str)).all():
                _df_x = None
            else:
                _iso = isolated_material_samples(df)
                _df_x = _df_x[~_df_x['filename'].isin(_iso)]
                if len(_df_x) == 0:
                    _df_x = None
        except Exception as e:
            print(f"Could not expand abbreviations: {type(e).__name__}: {e}")
            _df_x = None
    if _STUDY_DISEASE is not None:
        try:
            _assigned = study_disease_assignment(df, _STUDY_DISEASE)
            if _df_x is not None:
                _assigned = _merge_health(_assigned, study_disease_assignment(_df_x, _STUDY_DISEASE))
            if _assigned:
                df['MWB_study_health'] = df['filename'].map(lambda f: _assigned.get(f, (None, None))[0])
                df['MWB_study_doid'] = df['filename'].map(lambda f: _assigned.get(f, (None, None))[1])
        except Exception as e:
            print(f"Could not assign the study-level disease: {type(e).__name__}: {e}")
            df['MWB_study_health'] = None
            df['MWB_study_doid'] = None
    # health conditions (HPO) and exposures (ECTO/PECO), filled in translate_MWB_to_REDU_by_logic
    for _col in ['MWB_hpo_health', 'MWB_hpo', 'MWB_exposure', 'MWB_exposure_plant']:
        df[_col] = None
    if _STUDY_HPO is not None:
        try:
            _assigned = study_condition_assignment(df, _STUDY_HPO['hpo_ids'], _STUDY_HPO['doid_table'],
                                                   _STUDY_HPO['path_to_csvs'], study_has_doid=_STUDY_DISEASE is not None)
            if _df_x is not None:
                _assigned = _merge_health(_assigned, study_condition_assignment(
                    _df_x, _STUDY_HPO['hpo_ids'], _STUDY_HPO['doid_table'], _STUDY_HPO['path_to_csvs'],
                    study_has_doid=_STUDY_DISEASE is not None))
            if _assigned:
                df['MWB_hpo_health'] = df['filename'].map(lambda f: _assigned.get(f, (None, None))[0])
                df['MWB_hpo'] = df['filename'].map(lambda f: _assigned.get(f, (None, None))[1])
        except Exception as e:
            print(f"Could not assign health conditions: {type(e).__name__}: {e}")
            df['MWB_hpo_health'] = None
            df['MWB_hpo'] = None
        try:
            _assigned = study_exposure_assignment(df, _STUDY_HPO['exposure'], _STUDY_HPO['doid_table'],
                                                  _STUDY_HPO['path_to_csvs'])
            if _df_x is not None:
                _x = study_exposure_assignment(_df_x, _STUDY_HPO['exposure'], _STUDY_HPO['doid_table'],
                                               _STUDY_HPO['path_to_csvs'])
                _assigned = {**_x, **_assigned}
            if _assigned:
                df['MWB_exposure'] = df['filename'].map(lambda f: _assigned.get(f, (None, None))[0])
                df['MWB_exposure_plant'] = df['filename'].map(lambda f: _assigned.get(f, (None, None))[1])
        except Exception as e:
            print(f"Could not assign exposures: {type(e).__name__}: {e}")
            df['MWB_exposure'] = None
            df['MWB_exposure_plant'] = None
    if len(df) == 0:
        df['LatitudeandLongitude'] = df.apply(
            lambda x: f'{x.Latitude}|{x.Longitude}' if x.Latitude is not None and x.Longitude is not None and not np.isnan(
                x.Latitude) and not np.isnan(x.Longitude) else None, axis=1)

        df[['SampleType_inner', 'SampleTypeSub1_inner']] = df.apply(lambda x: get_blanks(x['Value']), axis=1,
                                                                    result_type='expand')
    else:
        # column-wise instead of DataFrame.apply(axis=1) (same values)
        df['LatitudeandLongitude'] = pd.Series(
            [f'{lat}|{lon}' if lat is not None and lon is not None and not np.isnan(lat) and not np.isnan(lon) else None
             for lat, lon in zip(df['Latitude'], df['Longitude'])], index=df.index, dtype=object)
        df[['SampleType_inner', 'SampleTypeSub1_inner']] = pd.DataFrame([get_blanks(v) for v in df['Value']], index=df.index)

    df = translate_MWB_to_REDU_from_csv(df, case='inner', path_to_csvs=path_to_csvs, add_to_cols=add_to_cols, ontology_table=ontology_table, allowedTerm_dict=allowedTerm_dict)

    
    df.loc[df['MWB_sampleSource'] == '_blank_', ['SampleType_inner', 'SampleTypeSub1_inner']] = ['blank_analysis', 'blank_analysis']
    # the Workbench's QC sample sources '_QC_' and '_POOLED_' are pooled QCs (pools of the study's samples, which
    # keep the study's organism), not blanks; samples of a mixture design (amount factors of reference materials,
    # 'SRM 1950: 500') stay as they were (blank_QC)
    # (nor are QC samples whose own values or names say blank, standard, solvent: mwb_qc_source_not_pool)
    _qc_source = df['MWB_sampleSource'].isin(QC_SAMPLE_SOURCES)
    _not_pool = mwb_mixture_design_samples(df)
    if _qc_source.any():
        _not_pool = _not_pool | mwb_qc_source_not_pool(df)
    df.loc[_qc_source & ~_not_pool, ['SampleType_inner', 'SampleTypeSub1_inner']] = ['pool_QC', 'pool_QC']
    df.loc[_qc_source & _not_pool, ['SampleType_inner', 'SampleTypeSub1_inner']] = ['blank_QC', 'blank_QC']
    df.loc[(df['MWB_sampleSource'] == 'culture media') | 
           (df['MWB_sampleSource'] == 'bacterial media') | 
           (df['MWB_sampleSource'] == 'adipocyte media') |
           (df['MWB_sampleSource'] == 'media'), ['SampleType_inner', 'SampleTypeSub1_inner']] = ['blank_culturemedia', 'blank_culturemedia']
    df.loc[(df['MWB_sampleSource'] == 'bacterial cells'), ['SampleType_inner', 'SampleTypeSub1_inner']] = ['culture_bacterial', 'culture_bacterial']
    # the Workbench's sample source of standard injections
    df.loc[df['MWB_sampleSource'] == STANDARD_SAMPLE_SOURCE, ['SampleType_inner', 'SampleTypeSub1_inner']] = list(STANDARD_SAMPLE_TYPE)

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
            # a 'QC' sample whose Subject / Sample ID names a pool or a body fluid is a pooled QC
            _pooled = _t.eq('blank_QC') & df.loc[_set, 'filename'].astype(str).map(_id_text).fillna('').str.contains(
                r'pool|(?<![a-z])(?:serum|plasma|blood|urine|saliva|csf|milk)(?![a-z])', regex=True)
            _t = _t.mask(_pooled, 'pool_QC')
            df.loc[_set, 'SampleType_inner'] = _t
            df.loc[_set, 'SampleTypeSub1_inner'] = _t.map(lambda t: 'missing value' if t == 'reference material' else t)
            print(f"QC/blank/reference samples from their values: {df.loc[_set, 'filename'].nunique()}")
    except Exception as e:
        print(f"Could not detect QC/blank samples: {type(e).__name__}: {e}")
    # standard injections named by a sample's values (standard_sample_type): only samples that got no sample type
    # above and have no QC / blank / reference-material value
    try:
        _std = pd.Series([standard_sample_type(k, v) for k, v in zip(df['Key'], df['Value'])], index=df.index, dtype=object)
        # the Subject ID as given ('Synthetic Standard'): the phrases only
        _std_samples = set(df.loc[_std.notna(), 'filename'])
        _std_samples |= {str(s_) for s_, subj in zip(df_additional['filename'], df_additional['SubjectIdentifierAsRecorded'])
                         if isinstance(subj, str) and standard_sample_type('subject id', subj)}
        if _std_samples:
            _qc_named = set(df.loc[[qc_sample_type(k, v) is not None for k, v in zip(df['Key'], df['Value'])], 'filename'])
            _typed = df.groupby('filename')['SampleType_inner'].transform(lambda s: s.notna().any())
            _set = ~_typed & df['filename'].astype(str).isin(_std_samples - _qc_named)
            if _set.any():
                df.loc[_set, ['SampleType_inner', 'SampleTypeSub1_inner']] = list(STANDARD_SAMPLE_TYPE)
                print(f"standard injections from their values: {df.loc[_set, 'filename'].nunique()}")
    except Exception as e:
        print(f"Could not detect standard injections: {type(e).__name__}: {e}")
    # the sample-level items that can name a reference material, by raw file (MWB_to_REDU_study_wrapper)
    mwb_note_refmat_evidence(df)
   
   #add environemntal columns
    
    df['MWB_envo_exact'] = envo_material_exact_name(df['MWB_sampleSource'], allowedTerm_dict).values
    # the taxa the samples' own organism items name (apply_own_organisms, in create_dataframe_outer_dict)
    try:
        _own_taxa = sample_own_taxa(df_additional, allowedTerm_dict)
        df['MWB_own_taxa'] = df['filename'].astype(str).map(_own_taxa).fillna('')
    except Exception as e:
        print(f"Could not read the samples' own organisms: {type(e).__name__}: {e}")
        df['MWB_own_taxa'] = ''
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
    lineage = _rcf._get_ncbi_lineage(term.split('|')[0])
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
        'STUDY_ID': _section(MWB_mwTAB_dict, 'METABOLOMICS WORKBENCH').get('STUDY_ID', 'NA'),
        'STUDY_FASTING': study_fasting_text(MWB_mwTAB_dict),
        'CELL_STRAIN_DETAILS': str(_section(MWB_mwTAB_dict, 'SUBJECT').get('CELL_STRAIN_DETAILS', '') or '').strip()

    }

    ontology_table = kwargs['ontology_table']
    # the Cell Ontology labels (_CL_LABELS, read by assign_cell_lines and
    # cultured_study_sample_type) are set by the first cell-type lookup: make it here, so
    # that a study's rows never depend on which studies were converted before it
    try:
        _cl_lookup(ontology_table, allowedTerm_dict)
    except Exception as e:
        print(f"Could not load the Cell Ontology labels: {type(e).__name__}: {e}")

    if entry['ANALYSIS_TYPE'] != 'MS':
        raise ValueError("This is not an MS analysis!")

    df_outer = pd.DataFrame(entry, index=[0])
    df_outer['STUDY_uberon_exact'] = uberon_exact_name(df_outer['SAMPLE_TYPE'], ontology_table, allowedTerm_dict).values
    df_outer['STUDY_envo_exact'] = envo_material_exact_name(df_outer['SAMPLE_TYPE'], allowedTerm_dict).values
    try:
        df_outer['STUDY_celltype'] = [cell_type_term(v, ontology_table, allowedTerm_dict) if v else None
                                      for v in df_outer['CELL_STRAIN_DETAILS']]
    except Exception as e:
        print(f"Could not check CELL_STRAIN_DETAILS for a cell type: {type(e).__name__}: {e}")
        df_outer['STUDY_celltype'] = None

    #add NCBITaxonomy on dataset level
    #######
    unique_species = set(df_outer['SUBJECT_SPECIES'])
    unique_species_ncbiIDs = set(df_outer['TAXONOMY_ID'])
    # (a copy of the IDs per lookup: get_taxonomy_id_from_name__allowedTerms pops the ID from the set it is given, and the
    # rules below need the analysis' TAXONOMY_ID)
    processed_species = {species: get_taxonomy_id_from_name__allowedTerms(species, ncbi_id=set(unique_species_ncbiIDs), allowedTerm_dict=allowedTerm_dict) for species in unique_species}
    # names the lookup above misses for formatting reasons only: 'NCBI:txid51953',
    # 'Oat (Avena sativa L.)', 'Lactuca sativa L.'
    processed_species = {sp: (term if term is not None else study_species_variant_term(sp, df_outer['TAXONOMY_ID'].iloc[0], allowedTerm_dict))
                         for sp, term in processed_species.items()}
    # renamed and merged taxa: a TAXONOMY_ID merged into another one, a synonym in the
    # pinned names.dmp ('Blautia coccoides')
    processed_species = {sp: (term if term is not None else
                              _rcf.pinned_ncbi_current_term(df_outer['TAXONOMY_ID'].iloc[0], allowedTerm_dict)
                              or _rcf.pinned_ncbi_synonym_term(sp, allowedTerm_dict))
                         for sp, term in processed_species.items()}
    # other white space, a strain designation at the end of the species, one wrong character (mwb_organism_fallback)
    processed_species = {sp: (term if term is not None else
                              mwb_organism_fallback(sp, allowedTerm_dict, f"{entry['STUDY_ID']} SUBJECT_SPECIES"))
                         for sp, term in processed_species.items()}
    # 'Genus sp. <designation>': the genus for a lab's numbering ('Bacillus sp. 1') without its own TAXONOMY_ID, and
    # for a TAXONOMY_ID of another strain (_sp_name_term); a TAXONOMY_ID of an unrelated taxon is no ID of the text's
    # isolate (_sp_given_id)
    _given_ids = {str(i).strip() for i in unique_species_ncbiIDs}
    _given_ids |= {m.group(1) for m in map(_TXID_RE.match, _given_ids) if m}   # ('NCBI:txid51953')
    _given_ids |= {str(_taxon_node(i)) for i in _given_ids if _taxon_node(i) is not None}   # (merged IDs: the current one)
    processed_species = {sp: (_sp_name_term(sp, term, _sp_given_id(sp, term, _given_ids, allowedTerm_dict), allowedTerm_dict)
                              if term is not None else None)
                         for sp, term in processed_species.items()}
    # a bare 'Genus sp.' ('Streptomyces sp.', 'Bacillus spp', 'Bacillus sp. (in: firmicutes)') with the TAXONOMY_ID of a
    # species or strain of that genus: that taxon, also where the ID is written 'NCBI:txid...' or was merged
    # (_bare_sp_term; user decision 2026-10-18)
    for sp, term in list(processed_species.items()):
        own = _bare_sp_term(sp, _given_ids, allowedTerm_dict)
        if own and own != term:
            print(f"[NCBITaxonomy] {entry['STUDY_ID']}: SUBJECT_SPECIES {sp!r} with TAXONOMY_ID "
                  f"{df_outer['TAXONOMY_ID'].iloc[0]!r}: {own} (was {term})")
        processed_species[sp] = own or term
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


    #add YearOfAnalysis on dataset level: the acquisition date of the analysis (ANALYSIS,
    # else MS section), never later than the study's submission. CREATED_ON is the date
    # the mwTab file was (re)generated, not when the samples were analysed: no fallback.
    _latest = _STUDY_SUBMIT_YEAR
    if _latest is None:
        _latest = submission_year(_section(MWB_mwTAB_dict, 'STUDY').get('SUBMIT_DATE', ''))
    _acq_ms = _section(MWB_mwTAB_dict, 'MS').get('ACQUISITION_DATE', '')
    # without an acquisition date: the year the mwTab was created (CREATED_ON), only when
    # it is not after the submission year (the original deposition; a regenerated mwTab
    # says 2025/2026 and is no year of analysis)
    _created = analysis_year(_section(MWB_mwTAB_dict, 'METABOLOMICS WORKBENCH').get('CREATED_ON', ''), _latest) \
        if _latest is not None else None
    df_outer['YearOfAnalysis'] = [analysis_year(x, _latest) if analysis_year(x, _latest) is not None
                                  else analysis_year(_acq_ms, _latest) if analysis_year(_acq_ms, _latest) is not None
                                  else _created for x in df_outer['ACQUISITION_DATE']]

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
        lambda x: get_enviromental_water(str(x['SAMPLE_TYPE'] or '') + ' ' + str(x['SUBJECT_TYPE'] or ''))
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


    global _ANALYSIS_IN_VITRO, _ANALYSIS_ABBR, _ANALYSIS_AGE_UNIT
    _ANALYSIS_ABBR = []
    # unit of ages given without one: the study's AGE_OR_AGE_RANGE unit, else years in
    # human studies, else unclear (no age)
    _ANALYSIS_AGE_UNIT = study_age_unit(entry['SUBJECT_AGE'])
    global _ANALYSIS_SUBJECT_TEXTS, _ANALYSIS_SAMPLE_TYPE, _ANALYSIS_SUBJECT_TYPE, _ANALYSIS_TAXON
    _ANALYSIS_SAMPLE_TYPE = entry['SAMPLE_TYPE'] if isinstance(entry['SAMPLE_TYPE'], str) else ''
    # the analysis' organism, for the species words of body-part values (sample_body_parts)
    _ANALYSIS_TAXON = df_outer['NCBITaxonomy'].iloc[0] if not add_to_cols else None
    _ANALYSIS_SUBJECT_TYPE = entry['SUBJECT_TYPE'] if isinstance(entry['SUBJECT_TYPE'], str) else ''
    _ANALYSIS_SUBJECT_TEXTS = [str(v) for v in (entry['SAMPLE_TYPE'], entry['SUBJECT_SPECIES'], entry['SUBJECT_TYPE'],
                                                _section(MWB_mwTAB_dict, 'SUBJECT').get('GENOTYPE_STRAIN', ''))
                               if isinstance(v, str) and v.strip() and v.strip().upper() != 'NA']
    if _ANALYSIS_AGE_UNIT is None and _is_human_study(entry, df_outer['NCBITaxonomy'].iloc[0]):
        _ANALYSIS_AGE_UNIT = 'y'
    try:
        _doid_t = _STUDY_HPO['doid_table'] if _STUDY_HPO is not None else DOID_TABLE
        if _doid_t is not None and len(_doid_t):
            _plant = bool(re.search(r'plant', f"{entry['SUBJECT_TYPE']} {entry['SPECIES_GROUP']}", re.I))
            _ANALYSIS_ABBR = study_abbreviations(study_text(MWB_mwTAB_dict), _doid_t, path_to_csvs, plant=_plant)
    except Exception as e:
        print(f"Could not read the study's abbreviations: {type(e).__name__}: {e}")
        _ANALYSIS_ABBR = []
    try:
        _ANALYSIS_IN_VITRO = analysis_in_vitro(_section(MWB_mwTAB_dict, 'STUDY').get('STUDY_TITLE', ''),
                                               entry['SAMPLE_TYPE'], entry['SUBJECT_TYPE'])
    except Exception as e:
        print(f"Could not tell whether the study is in vitro: {type(e).__name__}: {e}")
        _ANALYSIS_IN_VITRO = True
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

    # the samples' own organisms where the analysis' species resolved (apply_own_organisms); a problem here must
    # never cost the study its rows
    if not add_to_cols and 'MWB_own_taxa' in df_outer_inner.columns and 'NCBITaxonomy' in df_outer_inner.columns \
            and len(df_outer_inner):
        try:
            df_outer_inner, _n_own = apply_own_organisms(df_outer_inner, df_outer_inner['NCBITaxonomy'].iloc[0])
            if _n_own:
                print(f"[own organism] {_n_own} rows take their samples' own organism (analysis species "
                      f"{df_outer['NCBITaxonomy'].iloc[0]}): "
                      + ', '.join(f'{t}: {n}' for t, n in df_outer_inner['NCBITaxonomy'].fillna('none').value_counts().sort_index().items()))
        except Exception as e:
            print(f"Could not check the samples' own organisms: {type(e).__name__}: {e}")
    df_outer_inner = df_outer_inner.drop(columns=['MWB_own_taxa'], errors='ignore')

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


# taxon-specific UBERON ranges (teleost 2xxxxxx, amphibian 3xxxxxx, insect 6xxxxxx; 4xxxxxx are vertebrate skeletal
# terms and stay): 'corpus' there is an amphibian skull bone, not the corpus of the stomach. Left out of the
# exact-name lookup (the translation sheet may still name them)
_TAXON_SPECIFIC_UBERON = re.compile(r'^UBERON:[236]\d{6}$')


def _uberon_exact_lookup(ontology_table, allowed):
    """lower-cased label / synonym (> 3 characters) -> label, for names that point to
    exactly one allowed, non-obsolete label (no taxon-specific UBERON class). Labels win over synonyms."""
    key = (id(ontology_table), len(ontology_table))
    if key in _UBERON_LOOKUP_CACHE:
        return _UBERON_LOOKUP_CACHE[key]
    labels, synonyms = {}, {}
    index = ontology_table['UBERONOntologyIndex'] if 'UBERONOntologyIndex' in ontology_table.columns \
        else pd.Series('', index=ontology_table.index)
    for label, syn, idx in zip(ontology_table['Label'], ontology_table['Synonym'], index.astype(str)):
        if not isinstance(label, str) or label not in allowed or label.lower().startswith('obsolete') \
                or _TAXON_SPECIFIC_UBERON.match(idx.strip()):
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


_SHEET_CACHE = {}


def _read_sheet(path, **kwargs):
    """pd.read_csv of a translation sheet, read once per run (a copy is returned, so
    callers may change it)."""
    key = (path, tuple(sorted(kwargs.items())))
    if key not in _SHEET_CACHE:
        _SHEET_CACHE[key] = pd.read_csv(path, **kwargs)
    return _SHEET_CACHE[key].copy()


# keys whose values name the sampled organism ('Species', 'Organism', 'Bacterial strain')
_SPECIES_KEY_RE = r'(?:^|[^a-z])(?:species|organisms?|strains?|taxon|taxonomy|genus|isolates?|serovar|serotype)(?:[^a-z]|$)'


def _miss_run_in_most_orders(n_hits, n_misses, run=100):
    """True if, of all orders of n_hits found and n_misses not found names, more than half
    have run misses in a row somewhere (exact count: the misses form n_hits + 1 groups
    between the hits; orders with every group shorter than run, by inclusion-exclusion)."""
    if n_misses < run:
        return False
    if n_hits == 0:
        return True
    from math import comb
    groups = n_hits + 1
    short = sum((-1) ** j * comb(groups, j) * comb(n_misses - run * j + n_hits, n_hits)
                for j in range(0, min(groups, n_misses // run) + 1))
    return 2 * short < comb(n_misses + n_hits, n_hits)


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

        df_translations = _read_sheet(path_to_csvs + "/{}.csv".format(str(col_csv_name)), encoding="ISO-8859-1",
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
                _extra = _read_sheet(_extra_sheet, dtype=str, keep_default_na=False)
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
                if col_csv_name == 'UBERONBodyPartName' and 'MWB_own_bodypart' in MWB_table.columns:
                    # the sample's own body part (sample_body_parts) wins over the REST sample_source
                    _own = MWB_table['MWB_own_bodypart']
                    MWB_table.loc[_own.notna(), col_csv_name] = _own[_own.notna()]
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
                # raw file names under variant keys and the Sample-ID placeholder rows are
                # no species names
                _not_species = (MWB_table['Key'] == 'raw_file_name') | (MWB_table['Key'] == 'Sample ID') \
                    | MWB_table['Key'].map(is_variant_raw_file_key)
                _names = _ncbi_name_lookup(allowedTerm_dict)
                # same result as get_taxonomy_id_from_name__allowedTerms(species) (first
                # allowed term with that name, case-insensitive), as a dictionary lookup;
                # values tried in sorted order (the set order changes with PYTHONHASHSEED)
                _species_key = MWB_table['Key'].astype(str).str.lower().str.contains(_SPECIES_KEY_RE, regex=True)
                key_values = set(MWB_table[~_not_species & _species_key]['Value'])
                other_values = set(MWB_table[~_not_species & ~_species_key]['Value']) - key_values
                # values of a species / organism / strain key are species names wherever
                # they are one (a 'Species' factor with two bacteria among hundreds of other
                # values kept them only in some runs)
                for species in sorted(key_values, key=str):
                    result = _names.get(str(species).lower())
                    if result is None:
                        # surrounding blanks, renamed taxa (synonyms in the pinned names.dmp)
                        result = _names.get(' '.join(str(species).lower().split())) \
                            or _rcf.pinned_ncbi_synonym_term(str(species), allowedTerm_dict) \
                            or mwb_organism_fallback(str(species), allowedTerm_dict, 'sample species value')
                    # ('Genus sp. <a lab's numbering>': the genus, _sp_name_term)
                    result = _sp_name_term(species, result, '', allowedTerm_dict)
                    if result is not None:
                        processed_species[species] = result
                # values of other keys ('Acquisition method: IDA', 'group: Salmo') are no
                # species names if 100 of all the values in a row are none - decided as most
                # of the orders of trying them decide it (the result used to change with the
                # set order); this rule no longer drops the values of species keys
                other_hits = {}
                for species in sorted(other_values, key=str):
                    result = _sp_name_term(species, _names.get(str(species).lower()), '', allowedTerm_dict)
                    if result is not None:
                        other_hits[species] = result
                n_hits = len(processed_species) + len(other_hits)
                if not _miss_run_in_most_orders(n_hits, len(key_values) + len(other_values) - n_hits):
                    processed_species.update(other_hits)

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
            if col_csv_name == 'UBERONBodyPartName' and 'MWB_no_study_bodypart' in MWB_table.columns:
                # not for samples of a multi-tissue analysis whose own body part is unknown (sample_body_parts)
                MWB_table.loc[MWB_table['MWB_no_study_bodypart'].fillna(False).astype(bool), 'REDU'] = 'missing value'
            MWB_table[col_csv_name] = MWB_table[col_csv_name].fillna(MWB_table['REDU'])
            if col_csv_name == 'UBERONBodyPartName':
                MWB_table = MWB_table.drop(columns=['REDU_UBERONOntologyIndex'])
                MWB_table = MWB_table.drop(columns=['MWB', 'REDU', origin_cols[index]])


    return MWB_table


# UBERONBodyPartName.csv translates muscle values ('Muscle', 'muscles') to the vertebrate term; the muscle of any
# other organism (a shrimp) is a 'muscle organ' (UBERON:0001630; missing if that is no allowed label)
VERTEBRATE_MUSCLE = 'skeletal muscle organ, vertebrate'
OTHER_MUSCLE = 'muscle organ'
VERTEBRATA = 7742


def non_vertebrate_muscle(t):
    """rows with VERTEBRATE_MUSCLE whose NCBITaxonomy is no vertebrate (or none) get OTHER_MUSCLE"""
    if 'UBERONBodyPartName' not in t.columns or _rcf._PINNED_LINEAGE is None:
        return t
    rows = t['UBERONBodyPartName'].astype(str).eq(VERTEBRATE_MUSCLE)
    if not rows.any():
        return t
    tax = t['NCBITaxonomy'] if 'NCBITaxonomy' in t.columns else pd.Series(None, index=t.index, dtype=object)
    # (several organisms: all of them vertebrates)
    vertebrate = tax.map(lambda x: isinstance(x, str) and '|' in x and all_parts(x, lambda p: taxon_is_ancestor(VERTEBRATA, p)))
    change = rows & ~vertebrate.astype(bool)
    if change.any():
        allowed = set((allowedTerm_dict.get('UBERONBodyPartName') or {}).get('allowed_values') or [])
        t = t.copy()
        t.loc[change, 'UBERONBodyPartName'] = OTHER_MUSCLE if OTHER_MUSCLE in allowed else np.nan
        print(f"[muscle] {int(change.sum())} rows of no vertebrate: {VERTEBRATE_MUSCLE} -> "
              f"{OTHER_MUSCLE if OTHER_MUSCLE in allowed else 'missing'}")
    return t


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
        if 'MWB_study_health' in MWB_table.columns:
            # study-level disease: whole organisms only (no cultures, blanks, QC), and
            # never against a per-sample HealthStatus / DOIDCommonName
            organism = MWB_table['SampleType'].astype(str).isin(['animal', 'plant']) \
                if 'SampleType' in MWB_table.columns else pd.Series(False, index=MWB_table.index)
            if 'SUBJECT_TYPE' in MWB_table.columns:
                organism &= ~MWB_table['SUBJECT_TYPE'].astype(str).str.contains(r'cell|cultur|synthetic|water', case=False)
            study_h = MWB_table['MWB_study_health']
            human = single_taxid_series(MWB_table['NCBITaxonomy'].astype(str)).eq('9606') \
                if 'NCBITaxonomy' in MWB_table.columns else pd.Series(False, index=MWB_table.index)
            study_h = study_h.where(study_h != 'healthy if human', np.where(human, 'healthy', None))
            h_missing, d_missing = _is_missing('HealthStatus'), _is_missing('DOIDCommonName')
            current_h = MWB_table['HealthStatus'].astype(str)
            case = organism & (study_h == 'unhealthy (NOS)') & (h_missing | (current_h != 'healthy'))
            MWB_table.loc[case & h_missing, 'HealthStatus'] = 'unhealthy (NOS)'
            fill = case & d_missing & MWB_table['MWB_study_doid'].notna()
            if fill.any():
                doid = MWB_table.loc[fill, 'MWB_study_doid'].astype(str).str.split('|')
                MWB_table.loc[fill, 'DOIDCommonName'] = doid.str[0]
                MWB_table.loc[fill, 'DOIDOntologyIndex'] = doid.str[1]
            healthy = organism & (study_h == 'healthy') & h_missing & d_missing
            MWB_table.loc[healthy, 'HealthStatus'] = 'healthy'
        if 'MWB_hpo_health' in MWB_table.columns or 'MWB_exposure' in MWB_table.columns:
            # health conditions and exposures: whole organisms only, never against a
            # per-sample HealthStatus 'healthy', and never over existing values
            organism = MWB_table['SampleType'].astype(str).isin(['animal', 'plant']) \
                if 'SampleType' in MWB_table.columns else pd.Series(False, index=MWB_table.index)
            if 'SUBJECT_TYPE' in MWB_table.columns:
                organism &= ~MWB_table['SUBJECT_TYPE'].astype(str).str.contains(r'cell|cultur|synthetic|water', case=False)
            if 'MWB_hpo_health' in MWB_table.columns:
                human = single_taxid_series(MWB_table['NCBITaxonomy'].astype(str)).eq('9606') \
                    if 'NCBITaxonomy' in MWB_table.columns else pd.Series(False, index=MWB_table.index)
                hpo_h = MWB_table['MWB_hpo_health']
                hpo_h = hpo_h.where(hpo_h != 'healthy if human', np.where(human, 'healthy', None))
                h_missing, d_missing, p_missing = _is_missing('HealthStatus'), _is_missing('DOIDCommonName'), _is_missing('HPOCommonName')
                current_h = MWB_table['HealthStatus'].astype(str)
                # HPO describes humans and other animals, not plants
                animal = organism & MWB_table['SampleType'].astype(str).eq('animal')
                case = animal & (hpo_h == 'unhealthy (NOS)') & (h_missing | (current_h != 'healthy'))
                fill = case & p_missing & MWB_table['MWB_hpo'].notna()
                if fill.any():
                    MWB_table.loc[fill, 'HPOCommonName'] = MWB_table.loc[fill, 'MWB_hpo'].astype(str).str.split('|').str[0]
                MWB_table.loc[case & h_missing & (fill | ~p_missing), 'HealthStatus'] = 'unhealthy (NOS)'
                healthy = animal & (hpo_h == 'healthy') & h_missing & d_missing & p_missing
                MWB_table.loc[healthy, 'HealthStatus'] = 'healthy'
            if 'MWB_exposure' in MWB_table.columns:
                plant = MWB_table['SampleType'].astype(str).eq('plant') \
                    if 'SampleType' in MWB_table.columns else pd.Series(False, index=MWB_table.index)
                term = MWB_table['MWB_exposure'].where(~plant, MWB_table['MWB_exposure_plant'])
                fill = organism & _is_missing('Exposure') & term.notna()
                if fill.any():
                    MWB_table.loc[fill, 'Exposure'] = term[fill].astype(str).str.split('|').str[0]
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
                # the study-level location is the same on every row: look each distinct
                # text up once
                _country_memo = {}

                def _country(v):
                    if isinstance(v, str):
                        if v not in _country_memo:
                            _country_memo[v] = country_from_location(v, allowed_countries)
                        return _country_memo[v]
                    return country_from_location(v, allowed_countries)
                MWB_table.loc[fill, 'Country'] = MWB_table.loc[fill, 'COLLECTION_LOCATION'].map(_country)
        if 'MWB_collection_date' in MWB_table.columns:
            dates = MWB_table['MWB_collection_date'].where(MWB_table['MWB_collection_date'].astype(str) != 'NA')
            if 'COLLECTION_TIME' in MWB_table.columns:
                # study-level collection time only when it is a real calendar date
                study_date = MWB_table['COLLECTION_TIME'].map(lambda v: normalize_redu_datetime(v) if isinstance(v, str) and re.search(r'(19|20)\d\d', v) else None)
                dates = dates.fillna(study_date)
            MWB_table['SampleCollectionDateandTime'] = dates
    except Exception as e:
        print(f"Could not derive additional sample fields: {type(e).__name__}: {e}")

    # time point, sample group and fasting state
    try:
        if 'MWB_timepoint' in MWB_table.columns:
            MWB_table['TimePointHoursRelativeToEvent'] = MWB_table['MWB_timepoint']
        if 'MWB_samplegroup' in MWB_table.columns:
            MWB_table['SampleGroupFreeText'] = MWB_table['MWB_samplegroup']
        # fasting: whole animals (humans included) only, no cultured cells, plants or
        # microbes; the study-level statement only for studies without per-sample values
        animal = MWB_table['SampleType'].astype(str).eq('animal') \
            if 'SampleType' in MWB_table.columns else pd.Series(False, index=MWB_table.index)
        if 'SUBJECT_TYPE' in MWB_table.columns:
            animal &= ~MWB_table['SUBJECT_TYPE'].astype(str).str.contains(r'cell|cultur|synthetic|water', case=False)
        fasting = MWB_table['MWB_fasting'] if 'MWB_fasting' in MWB_table.columns \
            else pd.Series(None, index=MWB_table.index, dtype=object)
        if fasting.isna().all() and 'STUDY_FASTING' in MWB_table.columns:
            fasting = MWB_table['STUDY_FASTING'].map(study_fasting_from_texts)
        MWB_table['FastingStatus'] = fasting.where(animal)
    except Exception as e:
        print(f"Could not derive time points / sample groups / fasting states: {type(e).__name__}: {e}")
        for _col in ['TimePointHoursRelativeToEvent', 'SampleGroupFreeText', 'FastingStatus']:
            MWB_table[_col] = np.nan

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
            if col == 'STUDY_uberon_exact' and 'MWB_no_study_bodypart' in MWB_table.columns:
                fill &= ~MWB_table['MWB_no_study_bodypart'].fillna(False).astype(bool)
            MWB_table.loc[fill, 'UBERONBodyPartName'] = MWB_table.loc[fill, col]

    # the vertebrate skeletal-muscle term only on vertebrates (non_vertebrate_muscle)
    try:
        MWB_table = non_vertebrate_muscle(MWB_table)
    except Exception as e:
        print(f"Could not check muscle labels against the organism: {type(e).__name__}: {e}")

    df_translations = _read_sheet(path_to_csvs + "/biofluid_tissue_distinction.csv", encoding="ISO-8859-1", dtype=str)

    MWB_table = MWB_table.merge(df_translations, left_on='UBERONBodyPartName', right_on='sampletype', how='left')
    
    # Fill NA values before any string conversion (astype(str) turns NaN into 'nan', which fillna cannot see)
    MWB_table['SampleTypeSub1'] = MWB_table['SampleTypeSub1'].fillna(MWB_table['tissue_vs_biofluid'])

    MWB_table = MWB_table.drop(columns=['sampletype', 'tissue_vs_biofluid'])

    # cell lines (Cellosaurus) and cell types; a problem here must never cost the study its rows
    try:
        MWB_table = assign_cell_lines(MWB_table)
    except Exception as e:
        print(f"Could not assign cell lines: {type(e).__name__}: {e}")
        MWB_table['CellLine'] = np.nan
    try:
        MWB_table = cultured_study_sample_type(MWB_table)
    except Exception as e:
        print(f"Could not check the sample type of a cultured-cell study: {type(e).__name__}: {e}")


    # List of columns to check and add if not present
    columns_to_check = ['NCBITaxonomy', 'SampleCollectionMethod', 'SampleExtractionMethod',
                        'InternalStandardsUsed', 'SubjectIdentifierAsRecorded', 'AgeInYears',
                        'BiologicalSex', 'UBERONBodyPartName', 'HealthStatus', 'DOIDCommonName',
                        'ComorbidityListDOIDIndex', 'Country', 'HumanPopulationDensity',
                        'LatitudeandLongitude', 'SmokingStatus', 'BodyMassIndex', 'Diet', 'HPOCommonName', 'Exposure',
                        'TimePointHoursRelativeToEvent', 'SampleGroupFreeText', 'FastingStatus', 'CellLine']

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


# subject fields of the rows that are no samples of the study's subjects (blanks, pooled QCs, reference
# materials, standards); blanks and standards also lose the organism and body part, pooled QCs keep the study's species and their
# body part (a pool of the study's samples), reference materials get those of the material
# (reference_materials.py, in complete_and_fill_REDU_table)
_MWB_SUBJECT_COLUMNS = ['AgeInYears', 'UniqueSubjectID', 'SubjectIdentifierAsRecorded', 'DOIDCommonName', 'DOIDOntologyIndex',
                        'BiologicalSex', 'HealthStatus', 'LatitudeandLongitude', 'SmokingStatus', 'BodyMassIndex', 'Diet',
                        'LifeStage', 'HPOCommonName', 'Exposure', 'TimePointHoursRelativeToEvent', 'FastingStatus', 'CellLine']
_MWB_ORGANISM_COLUMNS = ['NCBITaxonomy', 'UBERONBodyPartName', 'UBERONOntologyIndex']


def mwb_clear_not_subject_rows(redu_df_final):
    sample_type = redu_df_final['SampleType'].astype(str)
    blank = sample_type.str.contains('blank', case=False)
    # standard injections (standard_sample_type) are no sample of an organism either
    no_organism = blank | sample_type.eq(STANDARD_SAMPLE_TYPE[0])
    not_subject = sample_type.isin(['pool_QC', 'reference material']) | no_organism
    if not_subject.any():
        redu_df_final.loc[not_subject, [c for c in _MWB_SUBJECT_COLUMNS if c in redu_df_final.columns]] = 'missing value'
        redu_df_final.loc[no_organism, [c for c in _MWB_ORGANISM_COLUMNS if c in redu_df_final.columns]] = 'missing value'
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

    mwb_reset_refmat()
    _STUDY_WITHHELD.clear()
    del _STUDY_LISTING[:]
    raw_file_name_df = _get_metabolomicsworkbench_files(study_id)
    if len(raw_file_name_df) and 'USI' in raw_file_name_df.columns:
        _STUDY_LISTING.extend(raw_file_name_df['USI'].astype(str))
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

    global _STUDY_TEXT_REST, _STUDY_SUBMIT_YEAR
    _STUDY_TEXT_REST = ''
    _STUDY_SUBMIT_YEAR = None
    try:
        _summ = _mwb_get(f'https://www.metabolomicsworkbench.org/rest/study/study_id/{str(study_id)}/summary').json()
        if isinstance(_summ, dict):
            _STUDY_TEXT_REST = ' \n '.join(str(_summ.get(k) or '') for k in ('study_title', 'study_summary'))
            _STUDY_SUBMIT_YEAR = submission_year(_summ.get('submission_date'))
    except Exception as e:
        print(f"No study summary: {type(e).__name__}: {e}")

    global _STUDY_N_MS_ANALYSES
    _STUDY_N_MS_ANALYSES = sum(1 for a in stdy_info.values() if isinstance(a, dict) and a.get('analysis_type') == 'MS')

    # curated study-level disease category ('Cancer', 'Diabetes', ...); a missing or
    # unreadable answer only means no study-level disease
    global _STUDY_DISEASE, _STUDY_HPO, HPO_TABLE, EXPOSURE_TABLE, CELLOSAURUS_TABLE
    if kwargs.get('CellLineIndex_table') is not None:
        CELLOSAURUS_TABLE = kwargs['CellLineIndex_table']
    _STUDY_DISEASE = None
    _STUDY_HPO = None
    if kwargs.get('HPOOntologyIndex_table') is not None:
        HPO_TABLE = kwargs['HPOOntologyIndex_table']
    if kwargs.get('ExposureIndex_table') is not None:
        EXPOSURE_TABLE = kwargs['ExposureIndex_table']
    _doid_full = kwargs.get('DOIDOntologyIndex_table')
    if _doid_full is None:
        _doid_full = DOID_TABLE
    _has_conditions = (HPO_TABLE is not None and len(HPO_TABLE) > 0) or \
        (EXPOSURE_TABLE is not None and len(EXPOSURE_TABLE) > 0)
    _categories = []
    try:
        if (_doid_full is not None and len(_doid_full)) or _has_conditions:
            _dis = _mwb_get(f'https://www.metabolomicsworkbench.org/rest/study/study_id/{str(study_id)}/disease').json()
            _categories = mwb_disease_categories(_dis)
    except Exception as e:
        print(f"No study-level disease: {type(e).__name__}: {e}")
    try:
        if _doid_full is not None and len(_doid_full) and _categories:
            _STUDY_DISEASE = prepare_study_disease(_categories, _doid_full, path_to_csvs)
            if _STUDY_DISEASE is not None:
                print(f"Study disease: {[c['label'] for c in _STUDY_DISEASE['cats']]}")
    except Exception as e:
        print(f"No study-level disease: {type(e).__name__}: {e}")
    # health conditions and exposures: from the categories without a DOID term and from
    # the sample values (also in studies without a category)
    try:
        if _has_conditions:
            _hpo_ids, _exposure = prepare_study_conditions(_categories, _doid_full, path_to_csvs)
            _STUDY_HPO = {'hpo_ids': _hpo_ids, 'exposure': _exposure, 'doid_table': _doid_full,
                          'path_to_csvs': path_to_csvs}
    except Exception as e:
        print(f"No study-level health conditions / exposures: {type(e).__name__}: {e}")
        _STUDY_HPO = None

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
                                      raw_file_name_df=raw_file_name_df[['filename', 'filename_base', 'USI', 'archive']],
                                      path_to_csvs=path_to_csvs,
                                      Massive_ID=study_id + '|' + str(analysis_details["analysis_id"]),
                                      rest_analysis=analysis_details,
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


    redu_df_final = merge_repeated_fileobservations_across_mwatb(redu_df_final, polarity_table=polarity_table,
                                                                 file_contexts=listing_contexts(raw_file_name_df))
    try:
        redu_df_final = clear_polarity_contradicted_by_all_file_names(redu_df_final, polarity_table=polarity_table)
    except Exception as e:
        print(f"Could not check polarities against the file names: {type(e).__name__}: {e}")
    try:
        redu_df_final = propagate_to_same_name_files(redu_df_final, raw_file_name_df, polarity_table=polarity_table)
    except Exception as e:
        print(f"Could not copy metadata to same-name files: {type(e).__name__}: {e}")
    try:
        redu_df_final = clear_polarity_contradicted_by_sibling(redu_df_final, raw_file_name_df, polarity_table=polarity_table)
        redu_df_final = clear_chromatography_contradicted_by_path(redu_df_final, raw_file_name_df)
    except Exception as e:
        print(f"Could not check polarities / chromatographies against file paths: {type(e).__name__}: {e}")
    try:
        redu_df_final = clear_polarity_contradicted_by_folder_and_file_name(redu_df_final, polarity_table=polarity_table)
    except Exception as e:
        print(f"Could not check polarities against folders and file names: {type(e).__name__}: {e}")

    redu_df_final = mwb_clear_not_subject_rows(redu_df_final)
    # reference materials named in the samples' items and the study's text (complete_and_fill_REDU_table)
    redu_df_final[reference_materials.EVIDENCE_COLUMN] = mwb_refmat_evidence_column(redu_df_final['filename'])
    _refmat_study_text = sorted(_REFMAT_STUDY_NUMBERS | set(reference_materials.study_text_materials(_STUDY_TEXT_REST)))

    ontology_table = ontology_table.drop_duplicates(subset=['Label'])
    # DOID table (one row per label) so complete_and_fill_REDU_table can generate DOIDOntologyIndex
    doid_table = kwargs.get('DOIDOntologyIndex_table')
    if doid_table is None:
        doid_table = DOID_TABLE
    doid_kwargs = {}
    if doid_table is not None and len(doid_table):
        doid_kwargs['DOIDOntologyIndex_table'] = doid_table.drop_duplicates(subset=['Label'])[['Label', 'DOIDOntologyIndex']].copy()
    # HPO and ECTO/PECO tables (one row per label) for HPOOntologyIndex / ExposureIndex
    for _key, _table, _col in [('HPOOntologyIndex_table', HPO_TABLE, 'HPOOntologyIndex'),
                               ('ExposureIndex_table', EXPOSURE_TABLE, 'ExposureIndex')]:
        if _table is not None and len(_table):
            # the exposure table also carries the exposure agent (OrganismExposureAgent / Index)
            doid_kwargs[_key] = _table.drop_duplicates(subset=['Label'])[['Label', _col] + [c for c in ('OrganismExposureAgentIndex', 'OrganismExposureAgent') if c in _table.columns]].copy()
    # Cellosaurus table, so complete_and_fill_REDU_table can check CellLine and generate CellLineIndex
    if CELLOSAURUS_TABLE is not None and len(CELLOSAURUS_TABLE):
        doid_kwargs['CellLineIndex_table'] = CELLOSAURUS_TABLE
    redu_df_final = complete_and_fill_REDU_table(redu_df_final, allowedTerm_dict, UBERONOntologyIndex_table=ontology_table, **doid_kwargs,
                                                 ENVOEnvironmentBiomeIndex_table=ENVOEnvironmentBiomeIndex_table.drop_duplicates(subset=["Label"]),
                                                 ENVOEnvironmentMaterialIndex_table=ENVOEnvironmentMaterialIndex_table.drop_duplicates(subset=["Label"]),
                                                 NCBIRankDivision_table=NCBIRankDivision_table,
                                                 add_usi = True, other_allowed_file_extensions = ['.raw', '.cdf', '.wiff', '.d'],
                                                 refmat_study_text=_refmat_study_text)

    redu_df_final = apply_sex_distribution_guard(redu_df_final, rest_response, study_id=str(study_id))
    try:
        redu_df_final = finalize_sample_groups(redu_df_final, study_id=str(study_id))
    except Exception as e:
        print(f"Could not finalize sample groups: {type(e).__name__}: {e}")
        redu_df_final['SampleGroupFreeText'] = 'missing value'

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
            _drop_analysis(str(Massive_ID).split('|')[0], str(MWB_analysis_ID),
                           f"could not get the mwTab ({type(e).__name__}: {e})")
            return None

        try:
            mwTab_json = json.loads(mwb_anlysis_req.text, object_pairs_hook=handle_duplicates)

        except json.JSONDecodeError as e:
            mwTab_json = mwtab_from_invalid_json(mwb_anlysis_req.text, str(Massive_ID).split('|')[0], str(MWB_analysis_ID),
                                                 e, get_text=lambda url: _mwb_get(url).text)
            if mwTab_json is None:
                return None
        
    # reference materials the analysis' mwTab sections name by number (reference_materials.py)
    mwb_note_refmat_study_text(mwTab_json)
    allowedTerm_dict = kwargs['allowedTerm_dict']
    ontology_table = kwargs['ontology_table']

    try:
        mwTab_json = chromatography_from_rest_analysis(mwTab_json, kwargs.get('rest_analysis'))
    except Exception as e:
        print(f"Could not check the REST analysis record for the chromatography: {type(e).__name__}: {e}")

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
    
    

    # LifeStage brackets are human ages: from the age only for humans; other organisms only
    # get a life stage the study states (newborn, adolescent, ...)
    # (a sample of several organisms is not one of humans)
    is_human = single_taxid_series(complete_df['NCBITaxonomy'].astype(str)).eq('9606')
    complete_df['LifeStage'] = [age_category(a) if h else '' for a, h in zip(complete_df['AgeInYears'], is_human)]
    if 'MWB_lifestage' in complete_df.columns:
        fill = (complete_df['LifeStage'].isin(['', None]) | complete_df['LifeStage'].isna()) & complete_df['MWB_lifestage'].notna()
        complete_df.loc[fill, 'LifeStage'] = complete_df.loc[fill, 'MWB_lifestage']

    # food and beverage samples (sample_food_terms), not blanks / QCs / standards, as fits their
    # organism (food_term_of_organism); SampleType and NCBITaxonomy stay as they are (honey has
    # SampleType 'food', no organism)
    _food = complete_df['MWB_food'] if 'MWB_food' in complete_df.columns else pd.Series(None, index=complete_df.index)
    _not_sample = complete_df['SampleType'].astype(str).str.contains('blank', case=False) | \
        complete_df['SampleType'].isin(['pool_QC', 'reference material', 'personal care product', 'molecular_standards'])
    complete_df['FoodOntologyTerm'] = [food_term_of_organism(t, o, s) if isinstance(t, str) and not q else None
                                       for t, q, o, s in zip(_food, _not_sample, complete_df['NCBITaxonomy'], complete_df['SampleType'])]

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
                            "Diet",
                            "HPOCommonName",
                            "Exposure",
                            "TimePointHoursRelativeToEvent",
                            "SampleGroupFreeText",
                            "FastingStatus",
                            "CellLine",
                            "FoodOntologyTerm"
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
                           "HPOCommonName",
                           "Exposure",
                           "TimePointHoursRelativeToEvent",
                           "SampleGroupFreeText",
                           "FastingStatus",
                           "CellLine",
                           "FoodOntologyTerm",
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
    parser.add_argument("--path_to_doid_parents_csv", type=str, help="DOID is_a links from prepare_ontologies.py (optional); enables the most specific disease term", default='none')
    parser.add_argument("--path_to_hpo_csv", type=str, help="Prepared HPO ontology csv (optional); enables health conditions (HPOCommonName)", default='none')
    parser.add_argument("--path_to_exposure_csv", type=str, help="Prepared ECTO/PECO ontology csv (optional); enables exposures (Exposure)", default='none')
    parser.add_argument("--path_to_cellosaurus_tsv", type=str, help="Cellosaurus table from build_cellosaurus_table.py (optional); enables CellLine / CellLineIndex", default='none')
    parser.add_argument("--path_to_ncbi_nodes_dmp", type=str, required=True, help="nodes.dmp of the pinned NCBI taxdump (fetchReferenceData); NCBI lineages for SampleType, cell lines")
    parser.add_argument("--path_to_ncbi_names_dmp", type=str, required=True, help="names.dmp of the pinned NCBI taxdump")
    parser.add_argument("--path_to_ncbi_merged_dmp", type=str, required=True, help="merged.dmp of the pinned NCBI taxdump")
    parser.add_argument("--withheld_output", type=str, default="none",
                        help="Side table (USI, dataset, reason) of the files whose sample the name matching withheld "
                             "(saveOlderData --withheld_files; written with --study_id ALL, also when empty)")

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

    DOID_PARENTS_CSV = args.path_to_doid_parents_csv
    if args.path_to_doid_csv != 'none':
        DOID_TABLE = pd.read_csv(args.path_to_doid_csv, dtype=str)
        DOID_LOOKUP = build_doid_lookup(DOID_TABLE)
    if args.path_to_hpo_csv != 'none':
        HPO_TABLE = pd.read_csv(args.path_to_hpo_csv, dtype=str)
    if args.path_to_exposure_csv != 'none':
        EXPOSURE_TABLE = pd.read_csv(args.path_to_exposure_csv, dtype=str)
    if args.path_to_cellosaurus_tsv != 'none':
        from cellosaurus import load_cellosaurus_table
        CELLOSAURUS_TABLE = load_cellosaurus_table(args.path_to_cellosaurus_tsv)
   
    # NCBI lineages from the pinned taxonomy (no E-utilities requests)
    _rcf.use_pinned_ncbi_taxonomy(args.path_to_ncbi_nodes_dmp, args.path_to_ncbi_names_dmp, args.path_to_ncbi_merged_dmp)

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

        study_list = sorted(set(study_list))

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
                flush_study_withheld(study_id, result)
                if result is not None and len(result) > 0:
                    print('Extracted information for {} samples.'.format(len(result)))
                    all_results_list.append(result)
            except KeyboardInterrupt:
                raise
            except Exception as e:
                print(f"Failed to process {study_id}: {type(e).__name__}: {e}")
                failed_studies.append(study_id)
        print(f"{len(failed_studies)} of {len(study_list)} studies failed to process")
        if DROPPED_STUDIES:
            print(f"{len(DROPPED_STUDIES)} studies dropped because the Workbench gave no file list: "
                  + ', '.join(f'{s} ({r})' for s, r in sorted(DROPPED_STUDIES.items())))
        if DROPPED_ANALYSES:
            print(f"{len(DROPPED_ANALYSES)} MS analyses dropped because the Workbench gave no usable mwTab: "
                  + ', '.join(f'{a} ({r})' for a, r in sorted(DROPPED_ANALYSES.items())))
        merged_df = pd.concat(all_results_list, ignore_index=True)
        merged_df.to_csv('REDU_from_MWB_all.tsv', sep='\t', index=False, header=True)
        if args.withheld_output != 'none':
            write_withheld_files(args.withheld_output, WITHHELD_FILES)

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
