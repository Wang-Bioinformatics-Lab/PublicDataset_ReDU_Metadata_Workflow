import argparse
import pandas as pd
import re

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
                         'CellLine', 'CellLineIndex')]:
        df_merged[col] = df_merged[col].fillna('missing value').replace('', 'missing value')

    #unique by USI column, ignoring the file extension (e.g. x.raw vs x.mzML). Only known
    #MS file extensions are stripped: dots elsewhere in folder or file names (e.g.
    #'S.epidermidis/A.mzML', 'QC_1.5uL.mzML') belong to the identity of the file.
    df_merged["USI_noext"] = df_merged["USI"].str.replace(
        r"(\.(mzml|mzxml|mzdata|mgf|raw|cdf|wiff2?|d|gz|zip|lzma|bz2|7z|xz|rar))+$", "", regex=True, case=False)
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


    