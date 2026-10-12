import pandas as pd
import argparse
import json


def read_source(path, data_source, name=None):
    """Reads one source's table; an empty file (the source failed or timed out this run)
    gives an empty frame so the other sources are still merged."""
    name = name or data_source
    try:
        df = pd.read_csv(path, sep='\t', dtype=str)
    except pd.errors.EmptyDataError:
        print(f"WARNING: no new {name} metadata this run ({path} is empty)")
        df = pd.DataFrame()
    print(f"{name}: {len(df)} rows from {path}")
    df["DataSource"] = data_source
    return df


def main():
    # parsing arguments
    parser = argparse.ArgumentParser(description='Merge GNPS and ReDU metadata')
    parser.add_argument('--gnps_metadata')
    parser.add_argument('--mwb_metadata')
    parser.add_argument('--metabolights_metadata')
    parser.add_argument('--norman_metadata')
    parser.add_argument('--masst_metadata')
    parser.add_argument('--path_to_allowed_term_json')
    parser.add_argument('--output_metadata')
    args = parser.parse_args()

    # read allowed terms
    with open(args.path_to_allowed_term_json, 'r') as f:
        allowed_terms = json.load(f)

    columns_to_use = list(allowed_terms.keys())

    # replace value "MassiveID" with "ATTRIBUTE_DatasetAccession"
    columns_to_use = [col.replace("MassiveID", "ATTRIBUTE_DatasetAccession") for col in columns_to_use]

    # read GNPS metadata
    gnps_df = read_source(args.gnps_metadata, "GNPS")

    #drop duplicated files
    if 'USI' in gnps_df.columns:
        duplicates = gnps_df.duplicated(subset='USI', keep=False)
        gnps_df = gnps_df[~duplicates]

    # read Workbench, MetaboLights, NORMAN and MASST metadata
    mwb_df = read_source(args.mwb_metadata, "Workbench")
    metabo_df = read_source(args.metabolights_metadata, "MetaboLights")
    norman_df = read_source(args.norman_metadata, "NORMAN")
    masst_df = read_source(args.masst_metadata, "GNPS", name="MASST")

    # merge GNPS and ReDU metadata
    merged_df = pd.concat([gnps_df, mwb_df, metabo_df, masst_df, norman_df], ignore_index=True)

    # include only columns from gnps; a column no source provides is added empty
    # instead of failing the whole merge
    missing_columns = [col for col in columns_to_use if col not in merged_df.columns]
    if missing_columns:
        print(f"Columns missing from all sources, adding them empty: {missing_columns}")
    merged_df = merged_df.reindex(columns=columns_to_use)

    # drop eventual duplicates (could happen with MASST/GNPS and also for a different reason with MWB) MWB has duplicates across all rows GNPS/MASST we should keep the first entry (from GNPS)
    merged_df = merged_df.drop_duplicates(subset='USI')

    # write merged metadata to file
    merged_df.to_csv(args.output_metadata, sep='\t', index=False)

if __name__ == '__main__':
    main()