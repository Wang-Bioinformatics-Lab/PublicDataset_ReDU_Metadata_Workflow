import argparse
import collections
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
                         'HumanPopulationDensity', 'qiita_sample_name']
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


def replace_per_dataset(df_current, masst_only, df_older, files_table=None):
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

    # rebuilt datasets and the shrink guard (counted as files: USI without extension)
    n_fresh = usi_noext(cur_usi[fresh]).groupby(cur_ds[fresh]).nunique()
    old_nonmasst = ~old_masst_like & (old_usi != '')
    n_old = usi_noext(old_usi[old_nonmasst]).groupby(old_ds[old_nonmasst]).nunique()
    rebuilt = set(n_fresh.index) - {''}
    guarded = sorted(d for d in rebuilt if n_fresh.get(d, 0) < SHRINK_GUARD * n_old.get(d, 0))
    for d in guarded:
        print(f"shrink guard: {d} has {n_fresh.get(d, 0)} files of its own this run, {n_old.get(d, 0)} in earlier runs; "
              f"keeping its earlier rows")
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
    keep_old_over_masst = has_masst & ~old_silent & ~old_masst_like & ~in_guarded
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
        for i, k in old_keys.items():
            tokens = fresh_tokens.get(k)
            mine = _acq_tokens(old_usi[i].split(':', 2)[-1])
            other_acq[i] = bool(tokens) and bool(mine) and all(t and t != mine for t in tokens)
        if other_acq.any():
            print(f"Earlier rows kept as another acquisition of a file of this run: {int(other_acq.sum())} "
                  f"{old_ds[other_acq].value_counts().sort_index().to_dict()}")
        not_produced &= ~other_acq
    drop_old[not_produced] = 'file_not_produced_by_rebuilt_dataset'
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
                             'so a source that failed (0 rows this run) is visible.')
    parser.add_argument('--files_table', default=None,
                        help='dataset cache file table (uniquemri.csv): copies of a MassIVE file in another folder layout '
                             'count as the same file (merge_metadata.NameCopies)')
    parser.add_argument('--dropped_output', default=None,
                        help='Writes the rows dropped by the per-dataset replacement (USI, origin, reason, source, dataset).')
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
    if len(df_older) and 'USI' in df_older.columns:
        df_current, masst_only, df_older, dropped = replace_per_dataset(df_current, masst_only, df_older, args.files_table)
    if args.dropped_output:
        dropped.to_csv(args.dropped_output, sep='\t', index=False)

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
                         'CellLine', 'CellLineIndex', 'ProductCategory', 'FoodOntologyTerm', 'ColonizationStatus')]:
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
        summary.to_csv(args.summary_output, sep='\t')
        print(summary.to_string())
    df_merged = df_merged.drop(columns=['_origin'])

    #reindex
    df_merged.reset_index(drop=True, inplace=True)


    df_merged.to_csv(args.output_merged_data, index=False, sep = '\t')


    