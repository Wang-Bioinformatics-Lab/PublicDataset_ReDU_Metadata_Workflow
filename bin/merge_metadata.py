import pandas as pd
import argparse
import collections
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


# ---- MASST supplement: fill missing fields of files that a submitter/converter described ----
# fields filled together, only when all of them are missing on the ReDU row (name + index etc.)
SUPPLEMENT_GROUPS = [
    ['SampleType', 'SampleTypeSub1'],
    ['NCBITaxonomy', 'NCBIRank', 'NCBIDivision'],
    ['UBERONBodyPartName', 'UBERONOntologyIndex'],
    ['DOIDCommonName', 'DOIDOntologyIndex'],
    ['HPOCommonName', 'HPOOntologyIndex'],
    ['Exposure', 'ExposureIndex', 'OrganismExposureAgent', 'OrganismExposureAgentIndex'],
    ['FoodOntologyTerm', 'FoodOntologyIndex'],
    ['CellLine', 'CellLineIndex'],
    ['LatitudeandLongitude', 'DepthorAltitudeMeters'],
    ['AgeInYears', 'LifeStage'],
    ['SubjectIdentifierAsRecorded', 'UniqueSubjectID'],
    ['ENVOEnvironmentBiome', 'ENVOEnvironmentBiomeIndex'],
    ['ENVOEnvironmentMaterial', 'ENVOEnvironmentMaterialIndex'],
    ['ENVOBroadScale', 'ENVOBroadScaleIndex'],
    ['ENVOLocalScale', 'ENVOLocalScaleIndex'],
    ['ENVOMediumScale', 'ENVOMediumScaleIndex'],
]
# never filled from MASST
SUPPLEMENT_NEVER = {'USI', 'filename', 'ATTRIBUTE_DatasetAccession', 'MassiveID', 'DataSource'}
# subject/biology fields: not filled on blanks, QCs, standards or reference material
SUPPLEMENT_SUBJECT = {'NCBITaxonomy', 'UBERONBodyPartName', 'DOIDCommonName', 'HPOCommonName', 'Exposure',
                      'CellLine', 'LatitudeandLongitude', 'AgeInYears', 'SubjectIdentifierAsRecorded',
                      'ENVOEnvironmentBiome', 'ENVOEnvironmentMaterial', 'BiologicalSex', 'HealthStatus',
                      'SmokingStatus', 'BodyMassIndex', 'Diet', 'TimePointHoursRelativeToEvent', 'FastingStatus',
                      'ColonizationStatus', 'ComorbidityListDOIDIndex', 'SampleGroupFreeText', 'TermsofPosition'}
# fields that do not depend on what the sample is: filled even when the sample types disagree
SUPPLEMENT_TYPE_INDEPENDENT = {'ProductCategory', 'Country', 'SampleCollectionDateandTime'}


def _sample_class(value):
    v = '' if _is_missing(value) else str(value).strip()
    if not v:
        return 'missing'
    if 'blank' in v.lower() or v == 'pool_QC':
        return 'blank/QC'
    if v in ('reference material', 'molecular_standards'):
        return 'standard/reference'
    return 'sample'


def supplement_from_masst(primary, supplement, columns):
    """Fills missing fields of `primary` (one row per USI: the submitter/converter rows)
    from MASST rows of the same USI. Only missing fields are filled; linked fields
    (SUPPLEMENT_GROUPS) only together and only when all of them are missing. Rows where
    ReDU and MASST disagree on blank/QC vs sample are left as they are; where the sample
    types differ, only SUPPLEMENT_TYPE_INDEPENDENT fields are filled; subject fields are
    never filled on blanks, QCs, standards or reference material. Returns the frame, a
    table of the filled fields and a Counter of fills and conflicts."""
    stats = collections.Counter()
    fills = []
    if supplement is None or len(supplement) == 0 or len(primary) == 0:
        return primary, pd.DataFrame(columns=['USI', 'column', 'value']), stats
    groups = [[c for c in g if c in columns] for g in SUPPLEMENT_GROUPS]
    groups = [g for g in groups if g]
    grouped = {c for g in groups for c in g}
    units = groups + [[c] for c in columns if c not in grouped and c not in SUPPLEMENT_NEVER]
    position = {u: i for i, u in enumerate(primary['USI'])}
    supplement = supplement[supplement['USI'].isin(position)].drop_duplicates('USI')
    stats['supplement_rows_for_described_files'] = len(supplement)
    values = primary[columns].to_numpy(dtype=object)
    col_index = {c: i for i, c in enumerate(columns)}
    for row in supplement.to_dict('records'):
        i = position[row['USI']]
        redu = {c: values[i, col_index[c]] for c in columns}
        cr, cm = _sample_class(redu.get('SampleType')), _sample_class(row.get('SampleType'))
        redu_subject = any(not _is_missing(redu.get(c)) for c in SUPPLEMENT_SUBJECT if c in redu)
        if cr == 'blank/QC' and cm in ('sample', 'standard/reference'):
            stats['conflict_redu_blankqc_masst_sample'] += 1
            continue
        if cm == 'blank/QC' and (cr in ('sample', 'standard/reference') or (cr == 'missing' and redu_subject)):
            stats['conflict_redu_sample_masst_blankqc'] += 1
            continue
        same_type = cr == 'missing' or cm == 'missing' or \
            str(redu.get('SampleType')).strip() == str(row.get('SampleType')).strip()
        if not same_type:
            stats['conflict_sample_type_differs'] += 1
        subject_ok = cm == 'sample' and cr in ('sample', 'missing')
        if 'NCBITaxonomy' in redu and not _is_missing(redu['NCBITaxonomy']) and not _is_missing(row.get('NCBITaxonomy')) \
                and str(redu['NCBITaxonomy']).split('|')[0] != str(row.get('NCBITaxonomy')).split('|')[0]:
            stats['note_taxonomy_differs_(ReDU_kept)'] += 1
        filled_any = False
        for unit in units:
            head = unit[0]
            if _is_missing(row.get(head)):
                continue
            if not same_type and head not in SUPPLEMENT_TYPE_INDEPENDENT:
                continue
            if head in SUPPLEMENT_SUBJECT and not subject_ok:
                continue
            if not all(_is_missing(redu.get(c)) for c in unit):
                continue
            for c in unit:
                v = row.get(c)
                if not _is_missing(v):
                    values[i, col_index[c]] = v
                    fills.append((row['USI'], c, v))
                    stats[f'filled_{c}'] += 1
            filled_any = True
        # SampleTypeSub1 (grouped with SampleType) alone, below the same SampleType
        if cr != 'missing' and cr == cm and same_type and 'SampleTypeSub1' in col_index \
                and _is_missing(redu.get('SampleTypeSub1')) and not _is_missing(row.get('SampleTypeSub1')):
            values[i, col_index['SampleTypeSub1']] = row['SampleTypeSub1']
            fills.append((row['USI'], 'SampleTypeSub1', row['SampleTypeSub1']))
            stats['filled_SampleTypeSub1'] += 1
            filled_any = True
        if filled_any:
            stats['files_with_fills'] += 1
    out = pd.DataFrame(values, columns=columns, index=primary.index)
    return out[primary.columns], pd.DataFrame(fills, columns=['USI', 'column', 'value']), stats


def main():
    # parsing arguments
    parser = argparse.ArgumentParser(description='Merge GNPS and ReDU metadata')
    parser.add_argument('--gnps_metadata')
    parser.add_argument('--mwb_metadata')
    parser.add_argument('--metabolights_metadata')
    parser.add_argument('--norman_metadata')
    parser.add_argument('--masst_metadata')
    parser.add_argument('--masst_supplement', default=None,
                        help='MASST rows of files that GNPS metadata describes (MASST_to_REDU.py --supplement_output); '
                             'they only fill missing fields')
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
    from_masst = merged_df.index.isin(range(len(gnps_df) + len(mwb_df) + len(metabo_df),
                                            len(gnps_df) + len(mwb_df) + len(metabo_df) + len(masst_df)))

    # include only columns from gnps; a column no source provides is added empty
    # instead of failing the whole merge
    missing_columns = [col for col in columns_to_use if col not in merged_df.columns]
    if missing_columns:
        print(f"Columns missing from all sources, adding them empty: {missing_columns}")
    merged_df = merged_df.reindex(columns=columns_to_use)

    # MASST rows of files another source describes (the supplement rows of MASST_to_REDU and
    # MASST rows whose USI a GNPS/Workbench/MetaboLights row has) only fill missing fields
    supplement_parts = []
    if args.masst_supplement:
        supplement_parts.append(read_source(args.masst_supplement, "GNPS", name="MASST supplement")
                                .reindex(columns=columns_to_use))
    described = set(merged_df.loc[~from_masst, 'USI'].dropna())
    masst_only = set(merged_df.loc[from_masst, 'USI'].dropna()) - described
    supplement_parts.append(merged_df[from_masst & merged_df['USI'].isin(described).to_numpy()])

    # drop eventual duplicates (could happen with MASST/GNPS and also for a different reason with MWB) MWB has duplicates across all rows GNPS/MASST we should keep the first entry (from GNPS)
    merged_df = merged_df.drop_duplicates(subset='USI')

    supplement = pd.concat(supplement_parts, ignore_index=True)
    if len(supplement):
        n_before = len(merged_df)
        primary = ~merged_df['USI'].isin(masst_only)
        filled, fills, stats = supplement_from_masst(merged_df[primary], supplement, columns_to_use)
        merged_df = pd.concat([filled, merged_df[~primary]]).sort_index()
        assert len(merged_df) == n_before
        print(f"MASST supplement: {stats.get('files_with_fills', 0)} files gained fields from "
              f"{stats.get('supplement_rows_for_described_files', 0)} MASST rows of described files")
        for key, value in sorted(stats.items()):
            print(f"  {key}: {value}")
        if args.output_metadata:
            fills.to_csv(args.output_metadata + '.masst_supplement_fills.tsv', sep='\t', index=False)

    # rows only MASST describes; saveOlderData uses this (and drops the column) so that a source
    # that delivered nothing this run keeps its own earlier rows instead of MASST's
    merged_df['_masst_only'] = merged_df['USI'].isin(masst_only).map({True: 'yes', False: ''})

    # write merged metadata to file
    merged_df.to_csv(args.output_metadata, sep='\t', index=False)

if __name__ == '__main__':
    main()