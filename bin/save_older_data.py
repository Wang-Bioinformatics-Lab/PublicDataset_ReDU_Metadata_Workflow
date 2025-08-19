import argparse
import pandas as pd


#    python $TOOL_FOLDER/save_older_data.py ${merged_ch} ${older_redu_data} merged_with_old.tsv
if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Save older data')
    parser.add_argument('input_new_data', help='Input file path')
    parser.add_argument('older_data', help='Older data file path')
    parser.add_argument('output_merged_data', help='Output file path for merged data')
    args = parser.parse_args()

    df_current = pd.read_csv(args.input_new_data, dtype=str, sep = '\t')
    try:
        df_older = pd.read_csv(args.older_data, dtype=str, sep = '\t')
    except:
        df_older = pd.DataFrame()

    df_merged = pd.concat([df_current, df_older])

    #unique by USI column
    df_merged = df_merged.drop_duplicates(subset=['USI'])


    # make sure we have the dataset labels
    df_merged.loc[df_merged['USI'].str.startswith('mzspec:NORMAN', na=False), 'DataSource'] = 'NORMAN'
    df_merged.loc[df_merged['USI'].str.startswith('mzspec:MSV', na=False), 'DataSource'] = 'GNPS'
    df_merged.loc[df_merged['USI'].str.startswith('mzspec:MTBLS', na=False), 'DataSource'] = 'MetaboLights'
    df_merged.loc[df_merged['USI'].str.startswith('mzspec:ST', na=False), 'DataSource'] = 'Workbench'



    #reindex
    df_merged.reset_index(drop=True, inplace=True)


    df_merged.to_csv(args.output_merged_data, index=False, sep = '\t')


    