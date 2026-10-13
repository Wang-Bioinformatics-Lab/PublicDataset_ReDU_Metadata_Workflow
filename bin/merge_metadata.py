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


def _is_missing(value):
    return pd.isna(value) or str(value).strip().lower() in ('', 'nan', 'missing value')


def resolve_duplicate_usis(df, source_column='MetadataFileDataset'):
    """One row per USI. Exact duplicate rows are dropped first. For a USI that still has
    several rows (e.g. the same file described in two metadata files), the rows are
    collapsed as in merge_repeated_fileobservations: a column keeps its value when the
    rows agree (missing entries are not a conflict) and becomes 'missing value' when they
    give different values. Rows from a metadata file of the USI's own dataset (source_column
    equals the dataset accession) take precedence over rows from another dataset's file:
    a column takes the own rows' value when they have one, the other rows' only otherwise.
    Returns the frame and a table of the blanked conflicts."""
    empty_conflicts = pd.DataFrame(columns=['USI', 'column', 'values', 'kept_from_own_dataset'])
    if 'USI' not in df.columns or len(df) == 0:
        return df.drop(columns=[source_column], errors='ignore'), empty_conflicts
    if source_column in df.columns:
        own = (df[source_column].fillna('') == df['ATTRIBUTE_DatasetAccession'].fillna('')).to_numpy()
    else:
        own = pd.Series(True, index=df.index).to_numpy()
    df = df.drop(columns=[source_column], errors='ignore').assign(_own=own)
    n_rows = len(df)
    # an exact copy is dropped; if one copy comes from the dataset's own file, that one stays
    df = df.sort_values('_own', ascending=False, kind='stable')
    df = df[~df.drop(columns=['_own']).duplicated()].sort_index()
    n_exact = n_rows - len(df)
    repeated = df['USI'].duplicated(keep=False)
    if not repeated.any():
        print(f"GNPS: dropped {n_exact} exact duplicate rows; no conflicting USIs")
        return df.drop(columns=['_own']), empty_conflicts

    def collapse(values):
        present = [v for v in values if not _is_missing(v)]
        distinct = list(dict.fromkeys(str(v) for v in present))
        return (present[0] if len(distinct) == 1 else 'missing value'), distinct

    conflicts = []
    collapsed = []
    columns = [c for c in df.columns if c not in ('USI', '_own')]
    for usi, group in df[repeated].groupby('USI', sort=False):
        row = group.iloc[0].copy()
        own_rows, other_rows = group[group['_own']], group[~group['_own']]
        for col in columns:
            value, distinct = collapse(own_rows[col]) if len(own_rows) else ('missing value', [])
            from_own = bool(distinct)
            if not distinct:
                value, distinct = collapse(other_rows[col])
            if distinct:
                row[col] = value
            other = {str(v) for v in other_rows[col] if not _is_missing(v)} if from_own else set()
            if len(distinct) > 1 or (other - set(distinct)):
                conflicts.append((usi, col, ' || '.join(list(dict.fromkeys(distinct + sorted(other - set(distinct))))),
                                  from_own and len(distinct) == 1))
        collapsed.append(row)
    collapsed = pd.DataFrame(collapsed, columns=df.columns)
    # keep the position of the first row of each USI
    first = df[repeated].drop_duplicates('USI', keep='first')
    collapsed.index = first.index
    out = pd.concat([df[~repeated], collapsed]).sort_index().drop(columns=['_own'])
    conflicts = pd.DataFrame(conflicts, columns=['USI', 'column', 'values', 'kept_from_own_dataset'])
    blanked = conflicts[~conflicts['kept_from_own_dataset']]
    print(f"GNPS: dropped {n_exact} exact duplicate rows; {repeated.sum()} rows shared "
          f"{df.loc[repeated, 'USI'].nunique()} USIs and were collapsed to one row each "
          f"({len(conflicts) - len(blanked)} conflicting fields kept from the dataset's own metadata file, "
          f"{len(blanked)} fields on {blanked['USI'].nunique()} USIs set to 'missing value')")
    for column, n_usis in blanked.groupby('column').size().items():
        print(f"  blanked conflict in {column}: {n_usis} USIs")
    return out, conflicts


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

    # one row per file: exact duplicates are dropped, conflicting fields of a file
    # described more than once become 'missing value' (instead of dropping every copy)
    gnps_df, conflicts = resolve_duplicate_usis(gnps_df)
    if len(conflicts) and args.output_metadata:
        conflicts.to_csv(args.output_metadata + '.gnps_usi_conflicts.tsv', sep='\t', index=False)

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