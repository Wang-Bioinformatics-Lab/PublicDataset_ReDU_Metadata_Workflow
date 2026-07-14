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


    # reading metadata
    metadata_df = pd.read_csv(args.input_metadata, sep='\t')
    metadata_df.drop('USI', axis=1, inplace=True)

    # reading files
    files_df = pd.read_csv(args.input_files, sep='\t') #'study_id', 'file_path', 'USI' are headers

    metadata_df["MassiveID"] = metadata_df["MassiveID"].apply(lambda x: x.split('|')[0])

    files_df["key"] = files_df["study_id"] + ":" + files_df["file_path"]

    metadata_df["key"] = metadata_df["MassiveID"] + ":" + metadata_df["filename"]

    meta_cols = list(metadata_df.columns)

    # Stage 1: exact match on study:full_path. This is the confident path where the
    # filename recorded in the metadata already equals the real file path.
    exact_df = pd.merge(metadata_df, files_df[["key", "USI"]], on="key", how="inner")

    # Stage 2: basename fallback. Some studies record a bare filename in the assay
    # sheet (e.g. 'FILES/sample.raw') while the real file lives in a subdirectory
    # (e.g. 'FILES/<batch>/sample.raw'), so the exact key never matches and the whole
    # study is dropped by the inner join. Recover these by matching on the basename
    # within the study, but ONLY when that basename is unambiguous (appears once) on
    # both sides, so we never mis-link two distinct files that share a name.
    matched_keys = set(exact_df["key"])
    unmatched_df = metadata_df[~metadata_df["key"].isin(matched_keys)].copy()

    recovered_df = None
    if len(unmatched_df):
        def _basename(p):
            return str(p).split('/')[-1]

        files_bn = files_df.copy()
        files_bn["basename"] = files_bn["file_path"].apply(_basename)
        # unambiguous basenames within a study on the file-listing side
        f_counts = files_bn.groupby(["study_id", "basename"]).size()
        f_unique = files_bn.set_index(["study_id", "basename"]).loc[f_counts[f_counts == 1].index]
        bn_to_usi = dict(zip(f_unique.index, f_unique["USI"]))
        bn_to_path = dict(zip(f_unique.index, f_unique["file_path"]))

        unmatched_df["basename"] = unmatched_df["filename"].apply(_basename)
        # unambiguous on the metadata side too
        m_counts = unmatched_df.groupby(["MassiveID", "basename"]).size()
        m_unique_idx = set(m_counts[m_counts == 1].index)

        lookup_key = list(zip(unmatched_df["MassiveID"], unmatched_df["basename"]))
        unmatched_df["USI"] = [bn_to_usi.get(k) for k in lookup_key]
        real_path = [bn_to_path.get(k) for k in lookup_key]
        is_unique_meta = [k in m_unique_idx for k in lookup_key]

        import numpy as np
        keep = unmatched_df["USI"].notna().to_numpy() & np.array(is_unique_meta, dtype=bool)
        recovered_df = unmatched_df[keep].copy()
        # align filename to the real (subdirectory) path so it is consistent with USI
        recovered_df["filename"] = [p for p, k in zip(real_path, keep) if k]
        recovered_df = recovered_df[meta_cols + ["USI"]]

    exact_df = exact_df[meta_cols + ["USI"]]
    merged_df = exact_df if recovered_df is None else pd.concat([exact_df, recovered_df], ignore_index=True)

    # Filtering columsn to original
    merged_df = merged_df[meta_cols + ["USI"]]


    merged_df.rename(columns={'MassiveID': 'ATTRIBUTE_DatasetAccession'}, inplace=True)
    merged_df.drop('key', axis=1, inplace=True)
   

    #adding f. to each filename
    merged_df['filename'] = 'f.' + merged_df['filename']

    # Saving file
    merged_df.to_csv(args.output_merged_file, sep='\t', index=False)

if __name__ == '__main__':
    main()

