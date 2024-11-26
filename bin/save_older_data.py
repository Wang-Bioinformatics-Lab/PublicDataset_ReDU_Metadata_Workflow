import argparse
import pandas as pd


#    python $TOOL_FOLDER/save_older_data.py ${merged_ch} ${older_redu_data} merged_with_old.tsv
if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Save older data')
    parser.add_argument('input', help='Input file path')
    parser.add_argument('older_data', help='Older data file path')
    parser.add_argument('output', help='Output file path')
    args = parser.parse_args()

    

    df = pd.read_csv(args.input, dtype=str, sep = '\t')
    df_older = pd.read_csv(args.older_data, dtype=str, sep = '\t')

    df = pd.concat([df, df_older])

    #unique by USI column
    df = df.drop_duplicates(subset=['USI'])

    #reindex
    df.reset_index(drop=True, inplace=True)


    df.to_csv(args.output, index=False, sep = '\t')


    