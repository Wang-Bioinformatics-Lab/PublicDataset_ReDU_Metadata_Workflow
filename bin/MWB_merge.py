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
    metadata_df.drop('USI', axis=1, inplace=True)

    # reading files
    files_df = pd.read_csv(args.input_files, sep='\t') #'study_id', 'file_path', 'USI' are headers

    metadata_df["MassiveID"] = metadata_df["MassiveID"].apply(lambda x: x.split('|')[0])

    files_df["key"] = files_df["study_id"] + ":" + files_df["file_path"]

    metadata_df["key"] = metadata_df["MassiveID"] + ":" + metadata_df["filename"]



    # merging from both dataframes
    merged_df = pd.merge(metadata_df, files_df, on="key", how="inner")

    # Filtering columsn to original
    merged_df = merged_df[list(metadata_df.columns) + ["USI"]]

    
    merged_df.rename(columns={'MassiveID': 'ATTRIBUTE_DatasetAccession'}, inplace=True)
    merged_df.drop('key', axis=1, inplace=True)
   

    #adding f. to each filename
    merged_df['filename'] = 'f.' + merged_df['filename']

    # Saving file
    merged_df.to_csv(args.output_merged_file, sep='\t', index=False)


if __name__ == '__main__':
    main()

