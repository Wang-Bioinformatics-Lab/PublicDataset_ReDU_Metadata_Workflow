import pandas as pd
import argparse
import collections
import json
import os
import re


# ---- same file: one identity for the copies of a file under other paths ----
# Used by MASST_to_REDU.py (MASST files that GNPS metadata describes), the merge below (MASST
# rows of files another source describes) and save_older_data.py (rows from earlier runs).
# Known MS file / compression extensions (x.raw, x.mzML and x.mzML.gz are one file); dots
# elsewhere in folder or file names belong to the identity of the file
MS_FILE_EXTENSIONS = r"(\.(mzml|mzxml|mzdata|mgf|raw|cdf|wiff2?|d|gz|zip|lzma|bz2|7z|xz|rar))+$"
# MassIVE keeps copies of one uploaded file in these top folders: peak/ and raw/ (uploads),
# ccms_peak/ (GNPS's mzML conversion of them) and updates/<date>_<user>/<folder>/ (re-uploads)
_MASSIVE_UPDATES = r"^updates/[^/]+/"
_MASSIVE_COPY_FOLDERS = r"^(ccms_peak|peak|raw)/"


def usi_noext(usi):
    """USI without known MS file / compression extensions."""
    return usi.str.replace(MS_FILE_EXTENSIONS, "", regex=True, case=False)


def same_file_keys(usi):
    """Two identities per USI (a Series of USIs):
    - exact: dataset + path without MS file / compression extension;
    - copy: the identity across the copies of a file in other top folders. MassIVE (MSV,
      RMSV): the path below updates/<date>/ and below peak/, ccms_peak/ or raw/ (peak/A/x.mzXML,
      ccms_peak/A/x.mzML, raw/A/x.raw and updates/<date>/peak/A/x.mzML are one file).
      MetaboLights (MTBLS): dataset + file name (the files move between FILES/ subfolders), only
      valid where the name is unique (link_same_files). Other repositories: the exact identity.
    Returns (exact, copy, copy_needs_unique) as Series aligned with usi."""
    usi = usi.fillna('').astype(str)
    dataset = usi.str.extract(r'^mzspec:([^:]+):', expand=False).fillna('')
    path = usi.str.replace(r'^mzspec:[^:]+:', '', regex=True)
    exact = dataset + ':' + usi_noext(path)
    massive = dataset.str.match(r'^R?MSV\d')
    mtbls = dataset.str.startswith('MTBLS')
    copy = exact.copy()
    rel = path[massive].str.replace(_MASSIVE_UPDATES, '', regex=True).str.replace(_MASSIVE_COPY_FOLDERS, '', regex=True)
    copy[massive] = dataset[massive] + ':' + usi_noext(rel)
    copy[mtbls] = dataset[mtbls] + ':' + usi_noext(path[mtbls].str.rsplit('/', n=1).str[-1])
    return exact, copy, mtbls


def link_same_files(query_usi, target_usi):
    """For each query USI, the target USIs that are the same file (same_file_keys): the same
    USI if present, else the targets with the same exact identity, else those with the same
    copy identity. A MetaboLights file name links only where it names one file among the
    targets and one among the queries (FILES/pos/x.mzML and FILES/neg/x.mzML are two files).
    Returns a list (aligned with query_usi) of tuples of target USIs (empty: no copy)."""
    target_usi = pd.Series(pd.unique(target_usi.dropna().astype(str)))
    t_exact, t_copy, t_unique = same_file_keys(target_usi)
    by_usi = set(target_usi)
    by_exact = collections.defaultdict(list)
    by_copy = collections.defaultdict(list)
    for u, e, c in zip(target_usi, t_exact, t_copy):
        by_exact[e].append(u)
        by_copy[c].append(u)
    t_files_per_copy = pd.Series(t_exact.values).groupby(t_copy.values).nunique().to_dict()
    query_usi = query_usi.fillna('').astype(str)
    q_exact, q_copy, q_unique = same_file_keys(query_usi)
    q_files_per_copy = pd.Series(q_exact.values).groupby(q_copy.values).nunique().to_dict()
    out = []
    for u, e, c, needs_unique in zip(query_usi, q_exact, q_copy, q_unique):
        if not u:
            out.append(())
        elif u in by_usi:
            out.append((u,))
        elif e in by_exact:
            out.append(tuple(sorted(by_exact[e])))
        elif c in by_copy and not (needs_unique and (t_files_per_copy[c] > 1 or q_files_per_copy[c] > 1)):
            out.append(tuple(sorted(by_copy[c])))
        else:
            out.append(())
    return out


def _massive_top_and_stem(path):
    """(top folder, e.g. peak or updates/<date>/peak; lower-case file name without extension) of a MassIVE path"""
    m = re.match(_MASSIVE_UPDATES, path)
    prefix, rel = (m.group(0), path[m.end():]) if m else ('', path)
    top = prefix + rel.split('/', 1)[0] if '/' in rel else ''
    name = rel.rsplit('/', 1)[-1]
    return top, re.sub(MS_FILE_EXTENSIONS, '', name, flags=re.I).lower()


def name_copy_datasets(query_usi, target_usi, links):
    """MassIVE datasets where a query without a link (link_same_files) has a target of the same
    file name: the only datasets NameCopies needs from the dataset cache table"""
    def keys(usis):
        out = set()
        for u in usis:
            parts = str(u).split(':', 2)
            if len(parts) == 3 and re.match(r'^R?MSV\d', parts[1]):
                out.add((parts[1], _massive_top_and_stem(parts[2])[1]))
        return out
    unlinked = keys(u for u, t in zip(query_usi, links) if not t)
    return {d for d, _ in unlinked & keys(target_usi)}


class NameCopies:
    """Copies of one MassIVE file that sit in other folder layouts, e.g. the GNPS conversion
    ccms_peak/RAW/P1/x.mzML of raw/RAW/P1/x.raw and the submitter's own peak/mzML/x.mzML. Two
    groups of copies (same_file_keys 'copy' identity) are one file when the dataset's file listing
    (the dataset cache table, plus the given USIs) has exactly these two groups with that file
    name, they share no top folder (two files of one name in one top folder are two files, e.g.
    peak/Plate01/x and peak/Plate02/x), and where the cache counted spectra in both, the
    MS1/MS2 counts agree (pos/ and neg/ runs of one sample do not)."""

    def __init__(self, files_table, extra_usis, datasets):
        self.partner = {}
        datasets = {d for d in datasets if re.match(r'^R?MSV\d', d)}
        if not datasets:
            return
        usis, counts = [], {}
        if files_table and files_table != 'none' and os.path.exists(files_table) and os.path.getsize(files_table) > 0:
            try:
                for chunk in pd.read_csv(files_table, usecols=['usi', 'dataset', 'spectra_ms1', 'spectra_ms2'], dtype=str,
                                         keep_default_na=False, chunksize=1_000_000):
                    chunk = chunk[chunk['dataset'].isin(datasets)]
                    usis.extend(chunk['usi'])
                    counts.update({u: (a, b) for u, a, b in zip(chunk['usi'], chunk['spectra_ms1'], chunk['spectra_ms2'])
                                   if a not in ('', '0', 'nan') or b not in ('', '0', 'nan')})
            except Exception as e:  # without the listing only the other identities are used
                print(f"WARNING: could not read the dataset cache table {files_table} ({type(e).__name__}: {e}); "
                      f"copies in other folder layouts are not recognised")
                return
        else:
            print("No dataset cache table: copies of a file in other folder layouts are not recognised")
            return
        listed = set(usis)
        usis = pd.Series(list(dict.fromkeys(usis + [u for u in extra_usis
                                                     if u.split(':')[1:2] and u.split(':')[1] in datasets])))
        if len(usis) == 0:
            return
        _, copy, _ = same_file_keys(usis)
        groups = collections.defaultdict(lambda: collections.defaultdict(set))
        spectra = collections.defaultdict(set)
        for u, c in zip(usis, copy):
            d = u.split(':')[1]
            top, stem = _massive_top_and_stem(u.split(':', 2)[2])
            groups[(d, stem)][c].add(top)
            if u in counts:
                spectra[c].add(counts[u])
        for (d, stem), by_copy in groups.items():
            if len(by_copy) != 2:
                continue
            (c1, tops1), (c2, tops2) = sorted(by_copy.items())
            if tops1 & tops2 or '' in tops1 | tops2:
                continue
            if spectra[c1] and spectra[c2] and not (spectra[c1] & spectra[c2]):
                continue
            self.partner[c1], self.partner[c2] = c2, c1
        print(f"Dataset cache table: {len(listed)} files of {len(datasets)} datasets; "
              f"{len(self.partner) // 2} files with a copy in another folder layout")

    def link(self, query_usi, target_usi, links):
        """links (from link_same_files) completed with the name copies: a query without a link gets the
        targets that are the other copy group of its file"""
        if not self.partner:
            return links
        target_usi = pd.Series(pd.unique(target_usi.dropna().astype(str)))
        _, t_copy, _ = same_file_keys(target_usi)
        by_copy = collections.defaultdict(list)
        for u, c in zip(target_usi, t_copy):
            by_copy[c].append(u)
        _, q_copy, _ = same_file_keys(query_usi.fillna('').astype(str))
        return [t if t else tuple(sorted(by_copy.get(self.partner.get(c), ()))) for t, c in zip(links, q_copy)]


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
# disease fields: MASST tables often give one disease for all files of a study, so a disease is
# never filled into a row its own source marks as healthy or as a control / no-disease group
# (harmonisation turns 'no disease reported' into a missing DOIDCommonName)
SUPPLEMENT_DISEASE = {'DOIDCommonName', 'HPOCommonName'}
_NO_DISEASE_GROUP = re.compile(r'(?i)(?<![a-z0-9])(healthy|no[ _-]disease|disease[ _-]free|non[ _-]?diseased?|'
                               r'controls?|ctrl|normal[ _-]chow|ncd)(?![a-z0-9])')


def study_level_diseases(masst_rows):
    """(dataset, column) of the diseases (SUPPLEMENT_DISEASE) a MASST table gives to every real
    sample file of a dataset: a label of the study (whose controls, e.g. the wild-type mice of
    an Alzheimer's model study, are listed under it too), not of the file; these are not filled
    into the rows of other sources."""
    out = set()
    if masst_rows is None or len(masst_rows) == 0 or 'USI' not in masst_rows.columns:
        return out
    real = masst_rows[[_sample_class(v) == 'sample' for v in masst_rows.get('SampleType', pd.Series('', index=masst_rows.index))]]
    dataset = real['USI'].fillna('').astype(str).str.extract(r'^mzspec:([^:]+):', expand=False).fillna('')
    for column in sorted(SUPPLEMENT_DISEASE):
        if column not in real.columns:
            continue
        values = real[column].map(lambda v: '' if _is_missing(v) else str(v).strip())
        for d, vals in values.groupby(dataset):
            distinct = set(vals)
            if d and len(distinct) == 1 and '' not in distinct:
                out.add((d, column))
    return out


def _says_no_disease(redu):
    if str(redu.get('HealthStatus', '')).strip().lower() == 'healthy':
        return True
    group = redu.get('SampleGroupFreeText')
    return not _is_missing(group) and _NO_DISEASE_GROUP.search(str(group)) is not None


def _sample_class(value):
    v = '' if _is_missing(value) else str(value).strip()
    if not v:
        return 'missing'
    if 'blank' in v.lower() or v == 'pool_QC':
        return 'blank/QC'
    if v in ('reference material', 'molecular_standards'):
        return 'standard/reference'
    return 'sample'


def supplement_from_masst(primary, supplement, columns, study_level=frozenset()):
    """Fills missing fields of `primary` (one row per USI: the submitter/converter rows)
    from MASST rows of the same USI. Only missing fields are filled; linked fields
    (SUPPLEMENT_GROUPS) only together and only when all of them are missing. Rows where
    ReDU and MASST disagree on blank/QC vs sample are left as they are; where the sample
    types differ, only SUPPLEMENT_TYPE_INDEPENDENT fields are filled; subject fields are
    never filled on blanks, QCs, standards or reference material, a disease never on a row
    that is healthy or a control group (_says_no_disease) nor where it is a label of the whole
    study (study_level: study_level_diseases). Returns the frame, a table of the filled fields
    and a Counter of fills and conflicts."""
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
            if head in SUPPLEMENT_DISEASE and _says_no_disease(redu):
                stats[f'not_filled_{head}_row_healthy_or_control'] += 1
                continue
            if head in SUPPLEMENT_DISEASE and (row['USI'].split(':')[1] if row['USI'].count(':') >= 2 else '', head) in study_level:
                stats[f'not_filled_{head}_study_level_label'] += 1
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

    # merge GNPS and ReDU metadata. A MASST row never wins a USI tie: MASST rows of a USI (or
    # file) another source has are taken out below, before the duplicates are dropped
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
    # MASST rows of the same file as a GNPS/Workbench/MetaboLights/NORMAN row: the same USI, the
    # same path with another extension or a copy in another top folder, see same_file_keys) only
    # fill missing fields of that row, and only of a row that is the one such file
    supplement_parts = []
    if args.masst_supplement:
        supplement_parts.append(read_source(args.masst_supplement, "GNPS", name="MASST supplement")
                                .reindex(columns=columns_to_use))
    masst_rows = merged_df[from_masst]
    links = link_same_files(masst_rows['USI'], merged_df.loc[~from_masst, 'USI'])
    copy_of_described = pd.Series([bool(t) for t in links], index=masst_rows.index)
    fill_target = [t[0] if len(t) == 1 else '' for t in links]
    linked = masst_rows.assign(USI=fill_target)[[bool(t) for t in fill_target]]
    supplement_parts.append(linked)
    n_same_usi = sum(1 for u, t in zip(masst_rows['USI'], links) if t == (u,))
    print(f"MASST rows of files another source describes: {int(copy_of_described.sum())} "
          f"({n_same_usi} same USI, {int(copy_of_described.sum()) - n_same_usi} another copy of the file; "
          f"{len(linked)} fill missing fields of that file's row, "
          f"{int(copy_of_described.sum()) - len(linked)} match several rows and are only dropped)")
    if copy_of_described.any():
        by_source = (masst_rows.loc[copy_of_described.values, 'USI'].str.extract(r'^mzspec:([A-Z]+)', expand=False)
                     .value_counts().sort_index())
        print(f"  by dataset prefix: {by_source.to_dict()}")
    merged_df = merged_df[~merged_df.index.isin(masst_rows.index[copy_of_described.values])]
    masst_only = set(merged_df.loc[from_masst[merged_df.index], 'USI'].dropna())

    # drop eventual duplicates (could happen with MASST/GNPS and also for a different reason with MWB) MWB has duplicates across all rows GNPS/MASST we should keep the first entry (from GNPS)
    merged_df = merged_df.drop_duplicates(subset='USI')

    supplement = pd.concat(supplement_parts, ignore_index=True)
    if len(supplement):
        n_before = len(merged_df)
        primary = ~merged_df['USI'].isin(masst_only)
        # the MASST rows of all files (MASST-only rows and the rows of described files)
        study_level = study_level_diseases(pd.concat([masst_rows, supplement], ignore_index=True))
        if study_level:
            print(f"MASST diseases given to every sample file of a dataset (study labels, not filled into other "
                  f"sources' rows): {len(study_level)} datasets x columns")
        filled, fills, stats = supplement_from_masst(merged_df[primary], supplement, columns_to_use, study_level)
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