import os
import argparse
import pandas as pd

def main():
    # parsing args
    parser = argparse.ArgumentParser()
    parser.add_argument('input_metadata')
    parser.add_argument('input_files')
    parser.add_argument('output_merged_file')

    args = parser.parse_args()

    print(args)

    # reading metadata
    metadata_df = pd.read_csv(args.input_metadata, sep='\t')
    # # reading files
    files_df = pd.read_csv(args.input_files, sep='\t')
    files_df = files_df[["USI"]]
    # #metadata_df["MassiveID"] = metadata_df["MassiveID"].apply(lambda x: x.split('|')[0])
    # print('meta:')
    # print(metadata_df)

    # print('files:')
    # print(files_df)

    # files_df["short_filename"] = files_df["filename"].apply(lambda x: os.path.basename(x))
    # files_df["key"] = files_df["STUDY_ID"] + ":" + files_df["short_filename"]

    #metadata_df["key"] = metadata_df["MassiveID"] + ":" + metadata_df["filename"]

    # # merging from both dataframes
    merged_df = pd.merge(metadata_df, files_df, on="USI", how="inner")
    # print('merged: ')
    # print(merged_df)
    # # Filtering columsn to original
    # merged_df = merged_df[metadata_df.columns]

    metadata_df.rename(columns={'MassiveID': 'ATTRIBUTE_DatasetAccession'}, inplace=True)
    #metadata_df.drop('key', axis=1, inplace=True)

    metadata_df['filename'] = 'f.' + metadata_df['filename']

    # Saving file
    metadata_df.to_csv(args.output_merged_file, sep='\t', index=False)


if __name__ == '__main__':
    main()

