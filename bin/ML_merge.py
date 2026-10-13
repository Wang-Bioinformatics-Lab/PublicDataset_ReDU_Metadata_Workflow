import os
import re
from collections import Counter
import argparse
import numpy as np
import pandas as pd

# A supported MS file with a compression suffix (x.raw.zip, x.d.zip, x.CDF.gz, x.mzXML.lzma).
# Group 1 is the uncompressed path ('core'). Keep in sync with Metabolights2REDU.py.
_COMPRESSED_MS_FILE = re.compile(r'^(.*\.(?:mzml|mzxml|cdf|raw|wiff|d))\.(?:zip|gz|lzma|bz2|7z|xz|rar)$', re.I)


def _core(path):
    """Path without a compression suffix, lowercased; the path itself (lowercased) otherwise."""
    p = str(path)
    m = _COMPRESSED_MS_FILE.match(p)
    return (m.group(1) if m else p).lower()


def _is_compressed(path):
    return _COMPRESSED_MS_FILE.match(str(path)) is not None

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

    # When a file exists both uncompressed and compressed (x.raw and x.raw.zip, compared
    # case-insensitively within the study), only the uncompressed one is listed
    files_df["_core"] = files_df["file_path"].apply(_core)
    files_df["_compressed"] = files_df["file_path"].apply(_is_compressed)
    uncompressed = set(zip(files_df.loc[~files_df["_compressed"], "study_id"], files_df.loc[~files_df["_compressed"], "_core"]))
    shadowed = files_df["_compressed"] & pd.Series([k in uncompressed for k in zip(files_df["study_id"], files_df["_core"])],
                                                   index=files_df.index, dtype=bool)
    if shadowed.any():
        print(f"Compressed files ignored because the uncompressed file exists: {int(shadowed.sum())}")
    files_df = files_df[~shadowed].copy()

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
    # both sides, so we never mis-link two distinct files that share a name. Basenames
    # are compared case-insensitively (e.g. 'x.CDF' in the assay sheet vs 'x.cdf' on FTP);
    # names that differ only in case count as the same name and are thus ambiguous.
    matched_keys = set(exact_df["key"])
    unmatched_df = metadata_df[~metadata_df["key"].isin(matched_keys)].copy()

    recovered_df = None
    if len(unmatched_df):
        def _basename(p):
            return str(p).split('/')[-1].lower()

        files_bn = files_df.copy()
        files_bn["basename"] = files_bn["file_path"].apply(_basename)
        # unambiguous basenames within a study on the file-listing side
        f_counts = files_bn.groupby(["study_id", "basename"]).size()
        f_unique = files_bn.set_index(["study_id", "basename"]).loc[f_counts[f_counts == 1].index]
        bn_to_usi = dict(zip(f_unique.index, f_unique["USI"]))
        bn_to_path = dict(zip(f_unique.index, f_unique["file_path"]))

        unmatched_df["basename"] = unmatched_df["filename"].apply(_basename)
        # unambiguous on the metadata side too, counted over ALL metadata rows: a basename
        # that also appears on an exactly matched row is ambiguous
        all_basenames = metadata_df["filename"].apply(_basename)
        m_counts = all_basenames.groupby([metadata_df["MassiveID"], all_basenames]).size()
        m_unique_idx = set(m_counts[m_counts == 1].index)

        lookup_key = list(zip(unmatched_df["MassiveID"], unmatched_df["basename"]))
        unmatched_df["USI"] = [bn_to_usi.get(k) for k in lookup_key]
        real_path = [bn_to_path.get(k) for k in lookup_key]
        is_unique_meta = [k in m_unique_idx for k in lookup_key]

        # never recover onto a file that an exact match already claimed
        not_claimed = ~unmatched_df["USI"].isin(set(exact_df["USI"])).to_numpy()
        keep = unmatched_df["USI"].notna().to_numpy() & np.array(is_unique_meta, dtype=bool) & not_claimed
        recovered_df = unmatched_df[keep].copy()
        # align filename to the real (subdirectory) path so it is consistent with USI
        recovered_df["filename"] = [p for p, k in zip(real_path, keep) if k]
        recovered_df = recovered_df[meta_cols + ["USI"]]

    exact_df = exact_df[meta_cols + ["USI"]]
    merged_df = exact_df if recovered_df is None else pd.concat([exact_df, recovered_df], ignore_index=True)

    # Stage 3: compressed files. The assay sheet often names the uncompressed file (x.raw)
    # while the FTP holds x.raw.zip, or names one compression and the FTP has another
    # (x.CDF.gz vs x.CDF.zip); a compressed name is also resolved onto the uncompressed file
    # when only that exists. Rows still unmatched are compared on the path without the
    # compression suffix (case-insensitive): (a) the same path, else (b) the same basename
    # within the study. As in stage 2, only names that are unique on both sides are matched,
    # and never onto a file already claimed.
    still = metadata_df[~metadata_df["key"].isin(set(merged_df["key"]))].copy()
    n_compressed = 0
    if len(still):
        f_core_path = list(zip(files_df["study_id"], files_df["_core"]))
        f_core_bn = list(zip(files_df["study_id"], files_df["_core"].str.split("/").str[-1]))
        m_core = metadata_df["filename"].apply(_core)
        m_compressed = metadata_df["filename"].apply(_is_compressed)
        m_core_path = list(zip(metadata_df["MassiveID"], m_core))
        m_core_bn = list(zip(metadata_df["MassiveID"], m_core.str.split("/").str[-1]))
        claimed = set(merged_df["USI"])
        recovered = []
        for f_keys, m_keys in ((f_core_path, m_core_path), (f_core_bn, m_core_bn)):
            f_counts = Counter(f_keys)
            m_counts = Counter(m_keys)
            lookup = {k: (usi, path, comp) for k, usi, path, comp in
                      zip(f_keys, files_df["USI"], files_df["file_path"], files_df["_compressed"]) if f_counts[k] == 1}
            m_key_of = dict(zip(metadata_df.index, m_keys))
            keep_idx, usis, paths = [], [], []
            for idx in still.index:
                k = m_key_of[idx]
                hit = lookup.get(k)
                # only pairs where at least one side is compressed (other pairs are stages 1 and 2)
                if hit is None or m_counts[k] != 1 or hit[0] in claimed or not (hit[2] or m_compressed[idx]):
                    continue
                keep_idx.append(idx); usis.append(hit[0]); paths.append(hit[1])
                claimed.add(hit[0])
            if keep_idx:
                rec = still.loc[keep_idx].copy()
                rec["USI"] = usis
                rec["filename"] = paths
                recovered.append(rec[meta_cols + ["USI"]])
                still = still.drop(index=keep_idx)
        if recovered:
            rec = pd.concat(recovered, ignore_index=True)
            n_compressed = len(rec)
            merged_df = pd.concat([merged_df, rec], ignore_index=True)
    print(f"Rows matched via compression-suffix rule: {n_compressed}")

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

