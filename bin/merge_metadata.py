import pandas as pd
import argparse

def main():
    # parsing arguments
    parser = argparse.ArgumentParser(description='Merge GNPS and ReDU metadata')
    parser.add_argument('gnps_metadata')
    parser.add_argument('redu_metadata')
    parser.add_argument('metabolights_metadata')
    parser.add_argument('output_metadata')
    args = parser.parse_args()

    # read GNPS metadata
    gnps_df = pd.read_csv(args.gnps_metadata, sep='\t')
    gnps_df["DataSource"] = "GNPS"

    # read ReDU metadata
    redu_df = pd.read_csv(args.redu_metadata, sep='\t')
    redu_df["DataSource"] = "Workbench"

    # read MetaboLights metadata
    metabo_df = pd.read_csv(args.metabolights_metadata, sep='\t')
    metabo_df["DataSource"] = "MetaboLights"

    # merge GNPS and ReDU metadata
    merged_df = pd.concat([gnps_df, redu_df, metabo_df], ignore_index=True)

    # include only columns from gnps
    merged_df = merged_df[gnps_df.columns]

    # write merged metadata to file
    merged_df.to_csv(args.output_metadata, sep='\t', index=False)

if __name__ == '__main__':
    main()