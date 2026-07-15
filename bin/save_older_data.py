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
    args = parser.parse_args()

    df_current = pd.read_csv(args.input_new_data, dtype=str, sep = '\t')
    try:
        df_older = pd.read_csv(args.older_data, dtype=str, sep = '\t')
        #use re.sub(r'\?VersionId.*', '', value) to alter USI values if DataSource is NORMAN
        df_older.loc[df_older['DataSource'] == 'NORMAN', 'USI'] = df_older.loc[df_older['DataSource'] == 'NORMAN', 'USI'].apply(lambda x: re.sub(r'\?VersionId.*', '', x))
    except:
        df_older = pd.DataFrame()

    df_merged = pd.concat([df_current, df_older], ignore_index=True)

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
            # identity independent of the differing folder: dataset + basename (incl extension)
            id_key = mt['USI'].str.replace(r'^mzspec:([^:]+):.*/', r'\1/', regex=True)
            valid = mt['USI'].isin(valid_usi)
            group_has_valid = valid.groupby(id_key).transform('any')
            stale = mt.index[group_has_valid & (~valid)]
            df_merged = df_merged.drop(index=stale).reset_index(drop=True)
            print(f"MetaboLights stale-path rows dropped: {len(stale)}")

    #unique by USI column
    df_merged["USI_noext"] = df_merged["USI"].str.replace(r"(\.[^.]+)+$", "", regex=True)
    df_merged = df_merged.drop_duplicates(subset=["USI_noext"]).drop(columns=["USI_noext"])


    # make sure we have the dataset labels
    df_merged.loc[df_merged['USI'].str.startswith('mzspec:NORMAN', na=False), 'DataSource'] = 'NORMAN'
    df_merged.loc[df_merged['USI'].str.startswith('mzspec:MSV', na=False), 'DataSource'] = 'GNPS'
    df_merged.loc[df_merged['USI'].str.startswith('mzspec:MTBLS', na=False), 'DataSource'] = 'MetaboLights'
    df_merged.loc[df_merged['USI'].str.startswith('mzspec:ST', na=False), 'DataSource'] = 'Workbench'



    #reindex
    df_merged.reset_index(drop=True, inplace=True)


    df_merged.to_csv(args.output_merged_data, index=False, sep = '\t')


    