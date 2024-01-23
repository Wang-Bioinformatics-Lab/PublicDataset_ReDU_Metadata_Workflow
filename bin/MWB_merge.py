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

    # reading files
    files_df = pd.read_csv(args.input_files, sep='\t')

    metadata_df["MassiveID"] = metadata_df["MassiveID"].apply(lambda x: x.split('|')[0])

    print(metadata_df)

    print(files_df)

    files_df["short_filename"] = files_df["filename"].apply(lambda x: os.path.basename(x))
    files_df["key"] = files_df["STUDY_ID"] + ":" + files_df["short_filename"]

    metadata_df["key"] = metadata_df["MassiveID"] + ":" + metadata_df["filename"]

    # merging from both dataframes
    merged_df = pd.merge(metadata_df, files_df, on="key", how="inner")

    # Filtering columsn to original
    merged_df = merged_df[metadata_df.columns]
    print(merged_df)

    # Saving file
    merged_df.to_csv(args.output_merged_file, sep='\t', index=False)


if __name__ == '__main__':
    main()

