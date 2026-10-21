import argparse
import collections
import os
import pandas as pd
import re

from merge_metadata import usi_noext, link_same_files, NameCopies, name_copy_datasets


def _usi_noext(usi):
    """USI without known MS file / compression extensions (x.raw and x.mzML are one file)."""
    return usi_noext(usi)


SOURCE_PREFIXES = {'Workbench': 'mzspec:ST', 'MetaboLights': 'mzspec:MTBLS', 'GNPS': 'mzspec:MSV', 'NORMAN': 'mzspec:NORMAN'}
# fields only submitter / converter metadata gives, never a MASST table: a row from an earlier run
# without any of them may come from MASST (the published table does not record it)
SUBMITTER_ONLY_FIELDS = ['YearOfAnalysis', 'SampleCollectionMethod', 'SampleExtractionMethod', 'InternalStandardsUsed',
                         'MassSpectrometer', 'IonizationSourceAndPolarity', 'ChromatographyAndPhase',
                         'SubjectIdentifierAsRecorded', 'UniqueSubjectID', 'TermsofPosition', 'ComorbidityListDOIDIndex',
                         'HumanPopulationDensity', 'qiita_sample_name', 'IsolationSourceNCBITaxonomy',
                         'IsolationSourceUBERONBodyPartName', 'IsolationSourceENVOEnvironmentMaterial']
# a rebuilt dataset whose rows of its own shrink below this share of its earlier rows keeps its
# earlier rows (as before per-dataset replacement): e.g. a sheet that was only partly read
SHRINK_GUARD = 0.5
# acquisition tokens in a path (polarity, chromatography, extract): an earlier file whose path differs from
# every file of this run with the same name in these tokens is another acquisition of the sample (e.g.
# temp_neg/x.mzXML next to temp_pos/x.mzML), not an older copy of a file of this run
_ACQ_TOKEN = re.compile(r'(?i)(?<![a-z0-9])(pos|neg|positive|negative|hilic|c18|c8|rp|aq|org|aqueous|organic)(?![a-z0-9])')
_ACQ_NORM = {'positive': 'pos', 'negative': 'neg', 'aqueous': 'aq', 'organic': 'org'}


def _acq_tokens(path):
    return frozenset(_ACQ_NORM.get(t.lower(), t.lower()) for t in _ACQ_TOKEN.findall(re.sub(r'[_\\/.+-]', ' ', path)))


def _missing_value(value):
    return pd.isna(value) or str(value).strip().lower() in ('', 'nan', 'missing value')


_POLARITY_TOKENS = frozenset({'pos', 'neg'})
_CHROMATOGRAPHY_TOKENS = frozenset({'hilic', 'c18', 'c8', 'rp'})


def _polarity_value_tokens(value):
    """IonizationSourceAndPolarity -> path tokens it agrees with ('alternating': both, no polarity: none)"""
    v = str(value).lower()
    tokens = set()
    if 'positive' in v or 'alternating' in v:
        tokens.add('pos')
    if 'negative' in v or 'alternating' in v:
        tokens.add('neg')
    return tokens


def _chromatography_value_tokens(value):
    """ChromatographyAndPhase -> path tokens it agrees with: HILIC, C18 / C8 (and rp), reverse phase
    (NOS) with any reverse phase token; other phases (GC, PGC, ...) with none"""
    v = str(value).lower()
    tokens = set()
    if 'hilic' in v:
        tokens.add('hilic')
    if 'c18' in v:
        tokens |= {'c18', 'rp'}
    if re.search(r'(?<![0-9])c8(?![0-9])', v):
        tokens |= {'c8', 'rp'}
    if 'reverse phase' in v:
        tokens.add('rp')
        if '(nos)' in v:
            tokens |= {'c18', 'c8'}
    return tokens


def _contradicts_path(path_tokens, polarity, chromatography):
    """True when a row's own IonizationSourceAndPolarity or ChromatographyAndPhase (if filled)
    contradicts the acquisition tokens of its path: a polarity none of the path's polarity tokens
    (neg/x.mzML labelled positive), a chromatography none of its chromatography tokens (hilic/x.mzML
    labelled reverse phase (C18)). A missing value, or a path without such tokens, contradicts nothing."""
    path_polarity = path_tokens & _POLARITY_TOKENS
    if path_polarity and not _missing_value(polarity):
        mine = _polarity_value_tokens(polarity)
        if mine and not (mine & path_polarity):
            return True
    path_chromatography = path_tokens & _CHROMATOGRAPHY_TOKENS
    if path_chromatography and not _missing_value(chromatography):
        if not (_chromatography_value_tokens(chromatography) & path_chromatography):
            return True
    return False


def _name_keys(usi):
    """(dataset, lower-case file name without MS file extension) per USI of a Series ('\\' and '/' separate folders)"""
    dataset = usi.str.extract(r'^mzspec:([^:]+):', expand=False).fillna('')
    name = usi.str.replace(r'^mzspec:[^:]+:', '', regex=True).str.replace(r'^.*[\\/]', '', regex=True)
    return pd.Series(list(zip(dataset, usi_noext(name).str.lower())), index=usi.index)


def _source_of(usi):
    out = pd.Series('', index=usi.index)
    for source, prefix in SOURCE_PREFIXES.items():
        out[usi.str.startswith(prefix, na=False)] = source
    return out


def _dataset_of(usi):
    return usi.str.extract(r'^mzspec:([^:]+):', expand=False).fillna('')


def _missing(series):
    return series.fillna('').astype(str).str.strip().str.lower().isin(('', 'nan', 'missing value'))


def _cache_files_noext(files_table, datasets):
    """USIs without MS file extension of the files of these datasets in the dataset cache table, or
    None when the table is not given, empty or cannot be read"""
    if not files_table or files_table == 'none' or not os.path.exists(files_table) or os.path.getsize(files_table) == 0:
        return None
    files = set()
    try:
        for chunk in pd.read_csv(files_table, usecols=['usi', 'dataset'], dtype=str, keep_default_na=False,
                                 chunksize=1_000_000):
            chunk = chunk[chunk['dataset'].isin(datasets)]
            files.update(usi_noext(chunk['usi']))
    except Exception as e:
        print(f"WARNING: could not read the dataset cache table {files_table} ({type(e).__name__}: {e})")
        return None
    return files


def _cache_datasets_of_paths(files_table, paths):
    """{path below the accession, without MS file extension: datasets that have a file of that path in the
    dataset cache table} for these paths, or None when the table is not given, empty or cannot be read"""
    if not files_table or files_table == 'none' or not os.path.exists(files_table) or os.path.getsize(files_table) == 0:
        return None
    out = collections.defaultdict(set)
    paths = set(paths)
    try:
        for chunk in pd.read_csv(files_table, usecols=['usi', 'dataset'], dtype=str, keep_default_na=False,
                                 chunksize=1_000_000):
            path = usi_noext(chunk['usi'].str.replace(r'^mzspec:[^:]+:', '', regex=True))
            hit = path.isin(paths)
            for p, d in zip(path[hit], chunk.loc[hit, 'dataset']):
                out[p].add(d)
    except Exception as e:
        print(f"WARNING: could not read the dataset cache table {files_table} ({type(e).__name__}: {e})")
        return None
    return out


def _specific_paths(old_ds, old_path, same_biology_datasets, cache_datasets):
    """Per earlier row (dataset old_ds, path old_path): whether its path is a file of exactly one other
    dataset in the dataset cache table (cache_datasets: _cache_datasets_of_paths) and that dataset is one
    whose row of the path this run describes the same biology (same_biology_datasets: per row, a set).
    A generic path (peak/Blank_1.mzML, ccms_peak/Raw/3.mzML) that several other datasets have is not
    specific: a dataset made private or removed upstream may share it with another dataset of the same
    biology and still be another file."""
    out = []
    for d, p, same in zip(old_ds, old_path, same_biology_datasets):
        others = (cache_datasets.get(p) or set()) - {d}
        out.append(len(others) == 1 and next(iter(others)) in same)
    return pd.Series(out, index=old_ds.index, dtype=bool)


# the fields an earlier row must share with this run's row of the same path in another dataset to be
# taken for that file under a wrong accession (_file_of_another_dataset)
BIOLOGY_FIELDS = ('NCBITaxonomy', 'UBERONBodyPartName', 'SampleType')


def _biology(df, index):
    """(NCBITaxonomy, UBERONBodyPartName, SampleType) per row of df (missing values as ''), for the rows of index"""
    columns = []
    for c in BIOLOGY_FIELDS:
        if c in df.columns:
            v = df.loc[index, c].fillna('').astype(str).str.strip()
            columns.append(v.where(~v.str.lower().isin(('', 'nan', 'missing value')), ''))
        else:
            columns.append(pd.Series('', index=index))
    return list(zip(*columns)) if columns else [()] * len(index)


def _file_of_another_dataset(old_usi, old_source, cur_usi, cur_source, candidates, files_table, df_older=None, df_current=None):
    """Earlier MassIVE rows (of `candidates`) whose USI is not a file of their dataset in the dataset
    cache table (files_table) and whose path below the accession (without MS file extension) is a
    file of another dataset among this run's rows that describes the same biology: a row of that path
    in the other dataset has the same NCBITaxonomy, UBERONBodyPartName and SampleType (BIOLOGY_FIELDS;
    missing on both counts as the same, not all three missing). Generic paths (ccms_peak/Raw/3.mzML,
    peak/mzML/Blank_17.mzML) are files of many datasets; an earlier row that describes another
    organism, body part or sample type there is another file and is kept (also when its dataset is no
    longer in the cache table: made private or removed upstream). The path must also be specific: a file of
    exactly one other dataset in the cache table, the one with the same biology (_specific_paths), so an
    earlier row of a generic path (peak/Blank_1.mzML of a removed dataset, blanks of two datasets of one
    lab) is kept. Without the cache table: none."""
    out = pd.Series(False, index=old_usi.index)
    massive = candidates & (old_source == 'GNPS') & (old_usi != '')
    cur_massive = cur_source == 'GNPS'
    if not massive.any() or not cur_massive.any():
        return out
    def dataset_and_path(usi):
        return _dataset_of(usi), usi_noext(usi.str.replace(r'^mzspec:[^:]+:', '', regex=True))
    old_ds, old_path = dataset_and_path(old_usi[massive])
    cur_ds, cur_path = dataset_and_path(cur_usi[cur_massive])
    wanted = set(old_path)
    in_wanted = cur_path.isin(wanted)
    cur_bio = _biology(df_current, cur_path.index[in_wanted]) if df_current is not None else [()] * int(in_wanted.sum())
    old_bio = _biology(df_older, old_path.index) if df_older is not None else [()] * len(old_path)
    # path -> other dataset -> the biology of its rows of that path this run
    datasets_of_path = collections.defaultdict(lambda: collections.defaultdict(set))
    for d, p, b in zip(cur_ds[in_wanted], cur_path[in_wanted], cur_bio):
        datasets_of_path[p][d].add(b)
    shared = pd.Series([bool(set(datasets_of_path.get(p, {})) - {d}) for d, p in zip(old_ds, old_path)], index=old_ds.index)
    elsewhere = pd.Series([any(o != d and b in bios and any(b) for o, bios in datasets_of_path.get(p, {}).items())
                           for d, p, b in zip(old_ds, old_path, old_bio)], index=old_ds.index)
    if shared.any():
        print(f"Earlier MassIVE rows with the path of a file of another dataset this run: {int(shared.sum())}, "
              f"of them with other biology than every row of that path there (kept): {int((shared & ~elsewhere).sum())}")
    if not elsewhere.any():
        return out
    cache = _cache_files_noext(files_table, set(old_ds[elsewhere]))
    if cache is None:
        print(f"No dataset cache table: {int(elsewhere.sum())} earlier MassIVE rows with the path of a file of another "
              f"dataset this run are kept (their own dataset's files cannot be checked)")
        return out
    not_listed = ~usi_noext(old_usi[massive][elsewhere]).isin(cache)
    # the path must be a file of exactly one other dataset in the cache table (_specific_paths)
    candidates = not_listed[not_listed].index
    cache_datasets = _cache_datasets_of_paths(files_table, set(old_path[candidates])) if len(candidates) else {}
    bio_of = dict(zip(old_path.index, old_bio))
    same_biology = [{o for o, bios in datasets_of_path.get(p, {}).items() if o != d and bio_of[i] in bios}
                    for i, d, p in zip(candidates, old_ds[candidates], old_path[candidates])]
    specific = (_specific_paths(old_ds[candidates], old_path[candidates], same_biology, cache_datasets)
                if cache_datasets is not None else pd.Series(False, index=candidates, dtype=bool))
    if (~specific).any():
        print(f"Earlier MassIVE rows with the path and biology of a file of another dataset this run, not a file of "
              f"their own dataset, but the path is not specific (not a file of exactly that one other dataset in the "
              f"dataset cache table; kept): {int((~specific).sum())} {old_ds[candidates][~specific].value_counts().sort_index().to_dict()}")
    hit = specific[specific].index
    out[hit] = True
    if len(hit):
        # per dataset this run that has the files: rows, and the earlier datasets they were under
        by_other = collections.defaultdict(list)
        for d, p in zip(old_ds[hit], old_path[hit]):
            by_other[','.join(sorted(set(datasets_of_path[p]) - {d}))].append(d)
        print(f"Earlier MassIVE rows dropped, not a file of their dataset (dataset cache) but a file of another "
              f"dataset this run: {len(hit)} rows of {old_ds[hit].nunique()} datasets; "
              + '; '.join(f"files of {o}: {len(ds)} rows under {len(set(ds))} datasets ({', '.join(sorted(set(ds))[:3])}"
                          f"{', ...' if len(set(ds)) > 3 else ''})" for o, ds in sorted(by_other.items())))
    print(f"Earlier MassIVE rows with the path and biology of a file of another dataset this run: {int(elsewhere.sum())}, "
          f"of them files of their own dataset (kept): {int(elsewhere.sum()) - len(candidates)}")
    return out


# ---- files a converter withheld on purpose this run (side table, --withheld_files) ----
# A converter can read a file and deliberately give it no row (Metabolights2REDU.py: data files whose
# sample pairing it withheld, MTBLS718 / MTBLS13295). The side table (USI, dataset, reason) tells
# replace_per_dataset that this run has seen these files: earlier rows of them (the same file as
# merge_metadata.link_same_files: same USI, same path without extension, or a MetaboLights file name
# unique on both sides) are dropped, unless this run has a row of the file itself, and the files count
# as files of their dataset for the shrink guard and make the dataset rebuilt even without rows of its
# own. Any source can add to it; a source without rows of its own this run (outage) is not touched.
WITHHELD_COLUMNS = ['USI', 'dataset', 'reason']
# reasons of files a converter did not withhold but joined to another sample than the earlier runs did
# (Metabolights2REDU.py 'sample_pairing_changed'). Their earlier rows are dropped only in datasets with
# non-MASST rows of their own this run, and they do not count as files of their dataset: once the earlier
# table comes from this code their earlier rows are paired as this run pairs them, so a dataset whose
# files this run delivered no rows (e.g. its file listing failed) keeps its earlier rows.
CHANGED_REASONS = {'sample_pairing_changed'}


def read_withheld_files(paths):
    """The side tables of withheld files (WITHHELD_COLUMNS, tab-separated with header) as one table, each USI
    once (first reason); missing or empty files, and tables without a USI column, are skipped (with a log
    line). paths: a path or a list of paths (None: none)."""
    if paths is None:
        paths = []
    elif isinstance(paths, str):
        paths = [paths]
    tables = []
    for path in paths:
        if not path or path == 'none' or not os.path.exists(path) or os.path.getsize(path) == 0:
            print(f"No withheld files in {path!r}")
            continue
        try:
            t = pd.read_csv(path, sep='\t', dtype=str, keep_default_na=False)
        except pd.errors.EmptyDataError:
            continue
        if 'USI' not in t.columns:
            print(f"WARNING: {path} has no USI column; not used")
            continue
        for c in WITHHELD_COLUMNS:
            if c not in t.columns:
                t[c] = ''
        t['USI'] = t['USI'].str.strip()
        t = t[t['USI'] != '']
        t.loc[t['dataset'] == '', 'dataset'] = _dataset_of(t.loc[t['dataset'] == '', 'USI'])
        t.loc[t['reason'] == '', 'reason'] = 'withheld'
        tables.append(t[WITHHELD_COLUMNS])
    if not tables:
        return pd.DataFrame(columns=WITHHELD_COLUMNS)
    return pd.concat(tables, ignore_index=True).drop_duplicates(subset=['USI']).reset_index(drop=True)


def _withheld_links(old_usi, withheld, silent):
    """(withheld USIs used, per earlier row the reason of the withheld file it is a copy of ('' if none)).
    Withheld files of sources without rows of their own this run (silent) are not used."""
    reason = pd.Series('', index=old_usi.index)
    if withheld is None or not len(withheld):
        return pd.Series([], dtype=str), reason
    w_usi = withheld['USI'].astype(str)
    w_source = _source_of(w_usi)
    unused = w_source.isin(silent) | (w_source == '')
    if unused.any():
        print(f"Withheld files of sources without rows of their own this run (or of no known source), not used: "
              f"{int(unused.sum())} {w_source[unused].value_counts().to_dict()}")
    w_usi, w_reason = w_usi[~unused], withheld['reason'][~unused]
    if not len(w_usi):
        return w_usi, reason
    reason_of = dict(zip(w_usi, w_reason))
    # (a copy of a file is always in the same dataset: only the earlier rows of these datasets are compared)
    query = old_usi[_dataset_of(old_usi).isin(set(_dataset_of(w_usi)))]
    links = link_same_files(query, w_usi)
    reason[query.index] = [reason_of[t[0]] if t else '' for t in links]
    return w_usi, reason


def replace_per_dataset(df_current, masst_only, df_older, files_table=None, withheld=None):
    """Rows of this run and of earlier runs, per file (same_file_keys in merge_metadata.py:
    the same path without extension or a copy in another top folder), by precedence
      this run's non-MASST row > earlier non-MASST row > this run's MASST row > earlier MASST row.
    A dataset with at least one non-MASST row this run is rebuilt: its earlier rows for files
    this run does not produce are dropped, unless its non-MASST rows this run are fewer than
    SHRINK_GUARD of its earlier non-MASST rows (then its earlier rows are kept as before, i.e. only
    the same USI without extension replaces them). A dataset without rows this run (embargoed,
    removed upstream, source outage) keeps all its earlier rows; so does a source that delivered
    no rows of its own this run, also against MASST rows of the same files. Earlier rows without
    any SUBMITTER_ONLY_FIELDS count as MASST rows.
    With the dataset cache table (files_table) copies of a MassIVE file in another folder layout
    count as the same file too (merge_metadata.NameCopies).
    An earlier row of a rebuilt dataset whose path differs from every file of this run with its name in
    the acquisition tokens (_ACQ_TOKEN) is kept as another acquisition of the sample, unless its own
    polarity / chromatography contradicts its path (_contradicts_path). An earlier MassIVE row that is
    not a file of its dataset in the dataset cache table but whose path is a file of another dataset
    this run is dropped (a wrong accession; _file_of_another_dataset; not without the cache table).
    withheld: the side table of files the converters withheld this run (read_withheld_files): earlier rows
    of these files are dropped ('withheld_this_run_<reason>') where this run has no non-MASST row of the
    file, and the files count as files of their dataset this run (shrink guard, rebuilt datasets);
    pairing-changed files (CHANGED_REASONS) only drop earlier rows in datasets with non-MASST rows of
    their own this run, and do not count as files.
    Returns (df_current, df_older) with the replaced rows removed, and a table of the dropped rows."""
    dropped = []
    cur_usi = df_current['USI'].fillna('').astype(str)
    old_usi = df_older['USI'].fillna('').astype(str)
    cur_source, old_source = _source_of(cur_usi), _source_of(old_usi)
    cur_ds, old_ds = _dataset_of(cur_usi), _dataset_of(old_usi)
    fields = [c for c in SUBMITTER_ONLY_FIELDS if c in df_older.columns]
    old_masst_like = (pd.concat([_missing(df_older[c]) for c in fields], axis=1).all(axis=1) if fields
                      else pd.Series(False, index=df_older.index)) & old_source.isin(['GNPS', 'MetaboLights', 'Workbench'])
    silent = [source for source in SOURCE_PREFIXES if not ((cur_source == source) & ~masst_only).any()]
    # earlier rows of files a converter withheld this run (side table), and the withheld files used
    withheld_usi, withheld_reason = _withheld_links(old_usi, withheld, silent)
    # (pairing-changed files only in datasets with non-MASST rows of their own this run: CHANGED_REASONS)
    fresh_datasets = set(cur_ds[~masst_only])
    withheld_reason[withheld_reason.isin(CHANGED_REASONS) & ~old_ds.isin(fresh_datasets)] = ''
    is_withheld = withheld_reason != ''

    # the files of this run each earlier row is a copy of
    fresh = ~masst_only
    link_fresh = link_same_files(old_usi, cur_usi[fresh])
    link_masst = link_same_files(old_usi, cur_usi[masst_only])
    if files_table:
        datasets = name_copy_datasets(old_usi, cur_usi, [f or m for f, m in zip(link_fresh, link_masst)])
        name_copies = NameCopies(files_table, list(old_usi[old_ds.isin(datasets)]) + list(cur_usi[cur_ds.isin(datasets)]), datasets)
        link_fresh = name_copies.link(old_usi, cur_usi[fresh], link_fresh)
        link_masst = name_copies.link(old_usi, cur_usi[masst_only], link_masst)
    link_fresh = pd.Series(link_fresh, index=df_older.index)
    link_masst = pd.Series(link_masst, index=df_older.index)
    has_fresh = link_fresh.map(bool)
    has_masst = link_masst.map(bool) & ~has_fresh

    # rebuilt datasets and the shrink guard (counted as files: USI without extension); the files a converter
    # withheld this run are files of their dataset this run too (not the pairing-changed ones, CHANGED_REASONS)
    # withheld files count once per file name (without extensions) of their dataset, and not where this run
    # has a row of a file of that name: a raw file and its derived twin in another folder are one file, as the
    # earlier runs counted them (MTBLS718: FILES/RAW_FILES/x.raw and FILES/DERIVED_FILES/x.mzML)
    n_fresh = usi_noext(cur_usi[fresh]).groupby(cur_ds[fresh]).nunique()
    n_withheld = pd.Series(dtype=int)
    reason_of = dict(zip(withheld['USI'].astype(str), withheld['reason'])) if withheld is not None and len(withheld) else {}
    counted_usi = [u for u in withheld_usi if reason_of.get(u) not in CHANGED_REASONS]
    if len(counted_usi):
        w_usi = pd.Series(counted_usi, dtype=str)
        w_ds = _dataset_of(w_usi)
        w_stem = usi_noext(w_usi).str.rsplit('/', n=1).str[-1].str.replace(r'^mzspec:[^:]+:', '', regex=True).str.lower()
        c_stem = usi_noext(cur_usi[fresh]).str.rsplit('/', n=1).str[-1].str.replace(r'^mzspec:[^:]+:', '', regex=True).str.lower()
        own = set(zip(cur_ds[fresh], c_stem))
        new = pd.Series([(d, s) not in own for d, s in zip(w_ds, w_stem)], index=w_usi.index)
        n_withheld = w_stem[new].groupby(w_ds[new]).nunique()
        n_fresh = n_fresh.add(n_withheld, fill_value=0).astype(int)
    old_nonmasst = ~old_masst_like & (old_usi != '')
    n_old = usi_noext(old_usi[old_nonmasst]).groupby(old_ds[old_nonmasst]).nunique()
    rebuilt = set(n_fresh.index) - {''}
    guarded = sorted(d for d in rebuilt if n_fresh.get(d, 0) < SHRINK_GUARD * n_old.get(d, 0))
    for d in guarded:
        print(f"shrink guard: {d} has {n_fresh.get(d, 0)} files of its own this run"
              f"{f' ({n_withheld[d]} of them withheld)' if d in n_withheld.index else ''}, {n_old.get(d, 0)} in earlier runs; "
              f"keeping its earlier rows")
    if len(n_withheld):
        print(f"Datasets with files withheld this run: {len(n_withheld)} ({int(n_withheld.sum())} files); "
              f"rebuilt without rows of their own: {sorted(set(n_withheld.index) - set(cur_ds[fresh]))}")
    rebuilt -= set(guarded)
    in_rebuilt = old_ds.isin(rebuilt)
    in_guarded = old_ds.isin(set(guarded))

    # earlier rows of a silent source win over this run's MASST rows of the same files (da65f05)
    old_silent = old_source.isin(silent)
    drop_masst_usis = set()
    for i in df_older.index[has_masst & old_silent]:
        drop_masst_usis.update(link_masst[i])
    # an earlier non-MASST row wins over this run's MASST rows of the same file; this run's MASST
    # row wins over an earlier MASST row (guarded datasets: as before)
    keep_old_over_masst = has_masst & ~old_silent & ~old_masst_like & ~in_guarded & ~is_withheld
    for i in df_older.index[keep_old_over_masst]:
        drop_masst_usis.update(link_masst[i])
    same_noext = usi_noext(old_usi).isin(set(usi_noext(cur_usi)))
    drop_old = pd.Series('', index=df_older.index)
    drop_old[has_fresh & ~in_guarded] = 'other_copy_of_a_file_of_this_run'
    drop_old[has_masst & ~old_silent & old_masst_like & ~in_guarded] = 'masst_row_of_this_run_for_the_file'
    not_produced = ~has_fresh & ~has_masst & in_rebuilt
    # this run without any MASST rows (MASST_to_REDU or gnpsmatchName_masst failed): earlier MASST rows are
    # kept, as for a source that delivered no rows (nf_workflow.nf: "MASST rows from earlier runs are kept")
    if not masst_only.any():
        old_masst_rows = old_masst_like & old_source.isin(['GNPS', 'MetaboLights'])
        n_keep = int((not_produced & old_masst_rows).sum())
        if n_keep:
            print(f"No MASST rows this run: keeping {n_keep} earlier MASST rows of rebuilt datasets")
        not_produced &= ~old_masst_rows
    # an earlier file that is another acquisition (polarity / chromatography / extract) of a file of this
    # run with the same name was missed by this run (converter / matcher), it is not a stale row
    if not_produced.any():
        old_keys = _name_keys(old_usi[not_produced])
        fresh_keys = _name_keys(cur_usi[fresh])
        wanted = set(old_keys)
        fresh_tokens = collections.defaultdict(list)
        for u, k in zip(cur_usi[fresh][fresh_keys.isin(wanted)], fresh_keys[fresh_keys.isin(wanted)]):
            fresh_tokens[k].append(_acq_tokens(u.split(':', 2)[-1]))
        other_acq = pd.Series(False, index=df_older.index)
        # (kept only if the row's own polarity / chromatography does not contradict its path: the matcher
        # of this run may have rejected the file, e.g. a positive-mode sheet row on neg/x.mzML)
        contradicted = pd.Series(False, index=df_older.index)
        no_value = pd.Series('', index=df_older.index)
        polarity = df_older['IonizationSourceAndPolarity'] if 'IonizationSourceAndPolarity' in df_older.columns else no_value
        chromatography = df_older['ChromatographyAndPhase'] if 'ChromatographyAndPhase' in df_older.columns else no_value
        for i, k in old_keys.items():
            tokens = fresh_tokens.get(k)
            mine = _acq_tokens(old_usi[i].split(':', 2)[-1])
            other_acq[i] = bool(tokens) and bool(mine) and all(t and t != mine for t in tokens)
            if other_acq[i] and _contradicts_path(mine, polarity[i], chromatography[i]):
                other_acq[i], contradicted[i] = False, True
        if other_acq.any():
            print(f"Earlier rows kept as another acquisition of a file of this run: {int(other_acq.sum())} "
                  f"{old_ds[other_acq].value_counts().sort_index().to_dict()}")
        if contradicted.any():
            print(f"Earlier rows of another acquisition not kept, their polarity / chromatography contradicts their path: "
                  f"{int(contradicted.sum())} {old_ds[contradicted].value_counts().sort_index().to_dict()}")
        not_produced &= ~other_acq
    else:
        contradicted = pd.Series(False, index=df_older.index)
    drop_old[not_produced] = 'file_not_produced_by_rebuilt_dataset'
    drop_old[not_produced & contradicted] = 'other_acquisition_contradicts_its_path'
    # earlier rows of a file a converter withheld this run (no non-MASST row of the file this run; guarded
    # datasets too: the earlier row is the one the converter does not trust)
    # (not the earlier MASST-like rows: their metadata does not come from the converter's sample pairing)
    _w = is_withheld & ~has_fresh & ~old_masst_like
    drop_old[_w] = 'withheld_this_run_' + withheld_reason[_w]
    # an earlier MassIVE row whose file is not in its dataset (dataset cache table) but whose path is a
    # file of another dataset this run: a wrong accession (e.g. a sheet's MassiveID column filled down,
    # MSV000090179 -> MSV000090180, 181, ...), the file is a row of its own dataset this run
    drop_old[_file_of_another_dataset(old_usi, old_source, cur_usi, cur_source, drop_old == '', files_table,
                                      df_older=df_older, df_current=df_current)] = \
        'file_of_another_dataset_this_run'
    # (the same USI without extension as a row of this run: replaced as before)
    drop_old[(drop_old != '') & same_noext] = 'same_usi_noext_as_this_run'

    is_dropped_masst = cur_usi.isin(drop_masst_usis) & masst_only
    dropped.append(pd.DataFrame({'USI': cur_usi[is_dropped_masst], 'origin': 'this_run_masst',
                                 'reason': ['source_silent_this_run_earlier_row_kept' if s in silent
                                            else 'earlier_non_masst_row_kept' for s in cur_source[is_dropped_masst]],
                                 'source': cur_source[is_dropped_masst], 'dataset': cur_ds[is_dropped_masst]}))
    is_dropped_old = drop_old != ''
    dropped.append(pd.DataFrame({'USI': old_usi[is_dropped_old], 'origin': 'earlier_runs', 'reason': drop_old[is_dropped_old],
                                 'source': old_source[is_dropped_old], 'dataset': old_ds[is_dropped_old]}))
    dropped = pd.concat(dropped, ignore_index=True)

    if silent:
        print(f"Sources without rows of their own this run (their earlier rows are kept): {silent}")
    print(f"Per-dataset replacement: {len(rebuilt)} datasets rebuilt, {len(guarded)} kept by the shrink guard")
    if len(dropped):
        summary = dropped.groupby(['source', 'origin', 'reason']).agg(rows=('USI', 'size'), datasets=('dataset', 'nunique'))
        print("Rows dropped (per source):")
        print(summary.to_string())
    guard_sources = collections.Counter(_source_of(pd.Series(['mzspec:' + d + ':' for d in guarded])).tolist())
    print(f"Shrink guard tripped per source: {dict(sorted(guard_sources.items())) or 0}")
    return (df_current[~is_dropped_masst.to_numpy()], masst_only[~is_dropped_masst.to_numpy()],
            df_older[~is_dropped_old.to_numpy()], dropped)


#    python $TOOL_FOLDER/save_older_data.py ${merged_ch} ${older_redu_data} merged_with_old.tsv
if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Save older data')
    parser.add_argument('input_new_data', help='Input file path')
    parser.add_argument('older_data', help='Older data file path')
    parser.add_argument('output_merged_data', help='Output file path for merged data')
    parser.add_argument('--metabolights_files', default=None,
                        help='Current MetaboLights file listing (study_id, file_path, USI). Used to drop '
                             'stale MetaboLights paths carried in from the older data when MetaboLights has '
                             'reorganised its folder layout (e.g. FILES/FILES/x -> FILES/DERIVED_FILES/x).')
    parser.add_argument('--summary_output', default=None,
                        help='Writes per-source counts of rows from this run vs. carried from earlier runs, '
                             'so a source that failed (0 rows this run) is visible, and own_this_run: the rows of this '
                             'run not only MASST describes (source_status.py).')
    parser.add_argument('--files_table', default=None,
                        help='dataset cache file table (uniquemri.csv): copies of a MassIVE file in another folder layout '
                             'count as the same file (merge_metadata.NameCopies)')
    parser.add_argument('--dropped_output', default=None,
                        help='Writes the rows dropped by the per-dataset replacement (USI, origin, reason, source, dataset).')
    parser.add_argument('--withheld_files', nargs='*', default=None,
                        help='Side tables of files the converters withheld this run (USI, dataset, reason): their earlier '
                             'rows are dropped and they count as files of their dataset (replace_per_dataset)')
    args = parser.parse_args()

    df_current = pd.read_csv(args.input_new_data, dtype=str, sep = '\t')
    # Only a missing/empty older file means "no history"; any other read error must fail
    # the run rather than silently discarding all older data.
    try:
        df_older = pd.read_csv(args.older_data, dtype=str, sep = '\t')
    except (FileNotFoundError, pd.errors.EmptyDataError):
        print(f"No older data found at {args.older_data}, continuing without it")
        df_older = pd.DataFrame()
    if 'DataSource' in df_older.columns:
        #strip ?VersionId... from USI values if DataSource is NORMAN
        is_norman = df_older['DataSource'] == 'NORMAN'
        df_older.loc[is_norman, 'USI'] = df_older.loc[is_norman, 'USI'].str.replace(r'\?VersionId.*', '', regex=True)

    # '_masst_only' (merge_metadata.py) marks the rows only MASST describes
    masst_only = (df_current.pop('_masst_only').fillna('') == 'yes') if '_masst_only' in df_current.columns \
        else pd.Series(False, index=df_current.index)
    dropped = pd.DataFrame(columns=['USI', 'origin', 'reason', 'source', 'dataset'])
    withheld = read_withheld_files(args.withheld_files)
    if len(withheld):
        print(f"Withheld files this run: {len(withheld)} of {withheld['dataset'].nunique()} datasets "
              f"{withheld['reason'].value_counts().to_dict()}")
    if len(df_older) and 'USI' in df_older.columns:
        df_current, masst_only, df_older, dropped = replace_per_dataset(df_current, masst_only, df_older, args.files_table,
                                                                        withheld=withheld)
    if args.dropped_output:
        dropped.to_csv(args.dropped_output, sep='\t', index=False)

    # rows of this run per source that not only MASST describes; none: the source delivered nothing this run and
    # keeps its earlier rows (replace_per_dataset, 'silent'). The summary's own_this_run, for source_status.py
    own_this_run = _source_of(df_current['USI'].fillna('').astype(str))[~masst_only.to_numpy()].value_counts()

    df_merged = pd.concat([df_current.assign(_origin='this_run'), df_older.assign(_origin='earlier_runs')],
                          ignore_index=True)

    # --- Resolve stale MetaboLights paths against the current file listing ---
    # MetaboLights occasionally reorganises where a study's files live (e.g. the same
    # file moves from FILES/FILES/<n>/x.mzML to FILES/DERIVED_FILES/<n>/x.mzML). The
    # older data then carries the stale path while the fresh run carries the current
    # one; because the USI paths differ, the plain USI dedup below keeps BOTH and the
    # same physical file appears twice. The freshly fetched listing is the source of
    # truth for which path exists right now, so for each MetaboLights file identity
    # (dataset + filename) we keep only the row(s) whose USI is in the current listing
    # and drop the stale alternative -- but only when a currently-valid path exists, so
    # studies that were not re-listed this run are never dropped.
    if args.metabolights_files:
        try:
            files_df = pd.read_csv(args.metabolights_files, sep='\t', dtype=str)
            valid_usi = set(files_df['USI'].dropna())
        except Exception:
            valid_usi = set()
        if valid_usi:
            is_mtbls = df_merged['USI'].str.startswith('mzspec:MTBLS', na=False)
            mt = df_merged.loc[is_mtbls, ['USI']].copy()
            # identity independent of the differing folder: dataset + basename (incl extension).
            # A compression suffix is not part of the identity: MetaboLights re-packed files
            # (FILES/RAW_FILES/x.d -> FILES/x.d.zip), so the older uncompressed path is stale
            # once the compressed file is listed
            id_key = mt['USI'].str.replace(r'^mzspec:([^:]+):.*/', r'\1/', regex=True)
            id_key = id_key.str.replace(r'(\.(mzml|mzxml|cdf|raw|wiff|d))\.(zip|gz|lzma|bz2|7z|xz|rar)$', r'\1',
                                        regex=True, case=False)
            valid = mt['USI'].isin(valid_usi)
            group_has_valid = valid.groupby(id_key).transform('any')
            stale = mt.index[group_has_valid & (~valid)]
            df_merged = df_merged.drop(index=stale).reset_index(drop=True)
            print(f"MetaboLights stale-path rows dropped: {len(stale)}")

    # rows from earlier runs may have empty ontology indexes; use the missing-value term
    # (and none for columns added later, e.g. HPOCommonName / Exposure)
    for col in [c for c in df_merged.columns if c.endswith('OntologyIndex') or (c.startswith('ENVO') and c.endswith('Index'))
                or c in ('HPOCommonName', 'Exposure', 'ExposureIndex', 'OrganismExposureAgent', 'OrganismExposureAgentIndex',
                         'TimePointHoursRelativeToEvent', 'SampleGroupFreeText', 'FastingStatus',
                         'CellLine', 'CellLineIndex', 'ProductCategory', 'FoodOntologyTerm', 'ColonizationStatus',
                         'IsolationSourceNCBITaxonomy', 'IsolationSourceUBERONBodyPartName',
                         'IsolationSourceENVOEnvironmentMaterial', 'IsolationSourceENVOEnvironmentMaterialIndex')]:
        df_merged[col] = df_merged[col].fillna('missing value').replace('', 'missing value')

    #unique by USI column, ignoring the file extension (e.g. x.raw vs x.mzML). Only known
    #MS file extensions are stripped: dots elsewhere in folder or file names (e.g.
    #'S.epidermidis/A.mzML', 'QC_1.5uL.mzML') belong to the identity of the file.
    df_merged["USI_noext"] = _usi_noext(df_merged["USI"])
    df_merged = df_merged.drop_duplicates(subset=["USI_noext"]).drop(columns=["USI_noext"])


    # make sure we have the dataset labels
    df_merged.loc[df_merged['USI'].str.startswith('mzspec:NORMAN', na=False), 'DataSource'] = 'NORMAN'
    df_merged.loc[df_merged['USI'].str.startswith('mzspec:MSV', na=False), 'DataSource'] = 'GNPS'
    df_merged.loc[df_merged['USI'].str.startswith('mzspec:MTBLS', na=False), 'DataSource'] = 'MetaboLights'
    df_merged.loc[df_merged['USI'].str.startswith('mzspec:ST', na=False), 'DataSource'] = 'Workbench'



    if args.summary_output:
        summary = (df_merged.groupby(['DataSource', '_origin'], dropna=False).size()
                   .unstack(fill_value=0).reindex(columns=['this_run', 'earlier_runs'], fill_value=0))
        summary['total'] = summary.sum(axis=1)
        summary['own_this_run'] = own_this_run.reindex(summary.index, fill_value=0)
        summary.to_csv(args.summary_output, sep='\t')
        print(summary.to_string())
    df_merged = df_merged.drop(columns=['_origin'])

    #reindex
    df_merged.reset_index(drop=True, inplace=True)


    df_merged.to_csv(args.output_merged_data, index=False, sep = '\t')


    