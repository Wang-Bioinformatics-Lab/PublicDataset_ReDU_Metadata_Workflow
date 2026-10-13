import argparse
import collections
import math
import pandas as pd
import glob
import re
import os
import json
from io import StringIO
from REDU_conversion_functions import age_category, normalize_redu_datetime

def _canonical_term(value):
    """Case-, space- and separator-insensitive form used to match submitted values to
    allowed terms ('Male' -> 'male', 'methanol:water (7:3)' -> 'methanol-water (7:3)').
    Ratios are reduced, so 'methanol-water (8:2)' matches 'methanol-water (4:1)'."""
    s = str(value).strip().lower()
    head, sep, tail = s.partition('(')
    head = re.sub(r'\s*[:/]\s*', '-', head)
    ratio = re.fullmatch(r'\s*(\d+(?::\d+)+)\s*\)(.*)', tail) if sep else None
    if ratio:
        parts = [int(x) for x in ratio.group(1).split(':')]
        g = 0
        for x in parts:
            g = math.gcd(g, x)
        if g > 1:
            parts = [x // g for x in parts]
        tail = ':'.join(str(x) for x in parts) + ')' + ratio.group(2)
    return re.sub(r'\s+', '', head + sep + tail)


# Lookups derived from the (large, unchanging) allowed-term lists, built once per run
# instead of once per metadata file
_LOOKUP_CACHE = {}


def _cached(name, source, build):
    key = (name, id(source), len(source))
    if key not in _LOOKUP_CACHE:
        _LOOKUP_CACHE[key] = build()
    return _LOOKUP_CACHE[key]


_DOID_ITEM = re.compile(r'(?i)^DOID[:_]?(\d+)$')


def normalize_doid_list(value, doid_ids, missing_value='missing value'):
    """'DOID:526|DOID:1596' style list of DOID IDs. Items may be separated by | , ; or
    whitespace and written DOID_526 / doid:526; each is written DOID:<digits>. IDs that are
    not in doid_ids (unknown or obsolete) are dropped, as are repeats (the submitter's order
    is kept). Returns missing_value when no valid ID remains."""
    text = str(value).strip()
    if not text:
        return missing_value
    text = re.sub(r'(?i)\bDOID\s*([:_])\s*', r'DOID\1', text)
    ids = []
    for item in re.split(r'[|,;\s]+', text):
        m = _DOID_ITEM.match(item.strip())
        if not m:
            continue
        doid = 'DOID:' + m.group(1)
        if doid in doid_ids and doid not in ids:
            ids.append(doid)
    return '|'.join(ids) if ids else missing_value


def complete_and_fill_REDU_table(df, allowedTerm_dict, add_usi = False, keep_usi = False, other_allowed_file_extensions = [], attempt_adding_file_extensions = False, case_insensitive_terms = False, **kwargs):
    """
    Completes and fills a REDU table with values based on a dictionary of allowed terms and missing values.

    Args:
    df: A pandas DataFrame containing the initial data.
    allowedTerm_dict: A dictionary containing allowed terms and missing values for each column.

    Returns:
    A DataFrame that has been filled with default values for missing columns,
    with values replaced by the corresponding "missing" value from the dictionary
    if they are not in the allowed terms or are missing/empty, except for specific columns.
    """

    if len(df) == 0:
        # On an empty (0-row) input, pandas' .map({}).fillna(...) chain below
        # produces float64-typed empty columns (no data to infer a dtype from),
        # which then fails when merged against the ontology tables' object-typed
        # 'Label' columns ("You are trying to merge on float64 and object
        # columns"). There's nothing to process for an empty table anyway, so
        # short-circuit before any of that logic runs.
        return df

    #prepare UberonOntologyIndex
    if 'UBERONOntologyIndex_table' in kwargs.keys():
        uberon_ontology_table = kwargs['UBERONOntologyIndex_table']
    else:
        uberon_ontology_table = pd.DataFrame(columns=['Label', 'UBERONOntologyIndex', 'Is Multicellular','Is Organ', 'Is Fluid'])

    uberon_ontology_table['Label'] = uberon_ontology_table['Label'].astype(str)

    #prepare DOIDOntologyTable
    if 'DOIDOntologyIndex_table' in kwargs.keys():
        doid_ontology_table = kwargs['DOIDOntologyIndex_table']
        if 'UBERONOntologyIndex' in doid_ontology_table.columns:
            doid_ontology_table.rename(columns={'UBERONOntologyIndex': 'DOIDOntologyIndex'}, inplace=True)
    else:
        doid_ontology_table = pd.DataFrame(columns=['Label', 'DOIDOntologyIndex'])
    
    doid_ontology_table['Label'] = doid_ontology_table['Label'].astype(str)

    #prepare HPO (health conditions) and exposure (ECTO/PECO) tables
    hpo_ontology_table = kwargs.get('HPOOntologyIndex_table')
    if hpo_ontology_table is None:
        hpo_ontology_table = pd.DataFrame(columns=['Label', 'HPOOntologyIndex'])
    hpo_ontology_table['Label'] = hpo_ontology_table['Label'].astype(str)
    exposure_ontology_table = kwargs.get('ExposureIndex_table')
    if exposure_ontology_table is None:
        exposure_ontology_table = pd.DataFrame(columns=['Label', 'ExposureIndex'])
    exposure_ontology_table['Label'] = exposure_ontology_table['Label'].astype(str)

    #prepare environment biome table
    if 'ENVOEnvironmentBiomeIndex_table' in kwargs.keys():
        envome_biome_ontology_table = kwargs['ENVOEnvironmentBiomeIndex_table']
    else:
        envome_biome_ontology_table = pd.DataFrame(columns=['Label', 'ENVOEnvironmentBiomeIndex'])

    envome_biome_ontology_table['Label'] = envome_biome_ontology_table['Label'].astype(str)

    #prepare environment material table
    if 'ENVOEnvironmentMaterialIndex_table' in kwargs.keys():
        envome_material_ontology_table = kwargs['ENVOEnvironmentMaterialIndex_table']
    else:
        envome_material_ontology_table = pd.DataFrame(columns=['Label', 'ENVOEnvironmentMaterialIndex'])

    envome_material_ontology_table['Label'] = envome_material_ontology_table['Label'].astype(str)

    #prepare environment medium scale table
    if 'ENVOEnvironmentMaterialIndex_table' in kwargs.keys():
        envome_medium_scale_ontology_table = kwargs['ENVOEnvironmentMaterialIndex_table']
        envome_medium_scale_ontology_table = envome_medium_scale_ontology_table.rename(
            columns={'ENVOEnvironmentMaterialIndex': 'ENVOMediumScaleIndex'}
            )
    else:
        envome_medium_scale_ontology_table = pd.DataFrame(columns=['Label', 'ENVOMediumScaleIndex'])

    envome_medium_scale_ontology_table['Label'] = envome_medium_scale_ontology_table['Label'].astype(str)


    #prepare environment local scale table
    if 'ENVOEnvironmentBiomeIndex_table' in kwargs.keys():
        envome_local_scale_ontology_table = kwargs['ENVOEnvironmentBiomeIndex_table']
        envome_local_scale_ontology_table = envome_local_scale_ontology_table.rename(
            columns={'ENVOEnvironmentBiomeIndex': 'ENVOLocalScaleIndex'}
            )
    else:
        envome_local_scale_ontology_table = pd.DataFrame(columns=['Label', 'ENVOLocalScaleIndex'])

    envome_local_scale_ontology_table['Label'] = envome_local_scale_ontology_table['Label'].astype(str)


    #prepare environment broad scale table
    if 'ENVOEnvironmentBiomeIndex_table' in kwargs.keys():
        envome_broad_scale_ontology_table = kwargs['ENVOEnvironmentBiomeIndex_table']
        envome_broad_scale_ontology_table = envome_broad_scale_ontology_table.rename(
            columns={'ENVOEnvironmentBiomeIndex': 'ENVOBroadScaleIndex'}
            )
    else:
        envome_broad_scale_ontology_table = pd.DataFrame(columns=['Label', 'ENVOBroadScaleIndex'])

    envome_broad_scale_ontology_table['Label'] = envome_broad_scale_ontology_table['Label'].astype(str)


    #prepare NCBI table
    if 'NCBIRankDivision_table' in kwargs.keys():
        NCBIRankDivision_table = kwargs['NCBIRankDivision_table']
    else:
        NCBIRankDivision_table = pd.DataFrame(columns=['TaxonID', 'NCBIRank', 'NCBIDivision'])

    NCBIRankDivision_table['TaxonID'] = NCBIRankDivision_table['TaxonID'].astype(str)
    NCBIRankDivision_table['NCBIRank'] = NCBIRankDivision_table['NCBIRank'].astype(str)
    NCBIRankDivision_table['NCBIDivision'] = NCBIRankDivision_table['NCBIDivision'].astype(str)

    # resolves NCBI names that are not allowed scientific names (synonyms, merged taxa,
    # strains) to an allowed term, or returns None
    ncbi_name_resolver = kwargs.get('NCBI_name_resolver')

    # Convert all columns to String (missing values become '' rather than the string 'nan')
    df = df.fillna('').astype(str)

    def _split_id_name(term: str):
        term = "" if term is None else str(term)
        if "|" in term:
            tid, tname = term.split("|", 1)
            return tid.strip(), tname.strip()
        return None, term.strip()

    def _is_sp_or_spp(name: str) -> bool:
        toks = name.strip().split()
        return len(toks) >= 2 and toks[1].rstrip(".").lower() in {"sp", "spp"}

    def _genus_from_name(name: str) -> str:
        return name.strip().split()[0] if name.strip() else ""

    #remove autogenerated columns if present
    keys_to_drop = [key for key, value in allowedTerm_dict.items() if value.get("generate") == "True"]
    columns_to_remove = [col for col in keys_to_drop + ['Is Multicellular', 'Is Organ', 'Is Fluid'] if col in df.columns]
    if keep_usi:
        print("Keeping USI column")
        columns_to_remove = [c for c in columns_to_remove if c != "USI"]
    # LifeStage is regenerated from the age below; keep the submitter's value for when it can't be
    if 'LifeStage' in df.columns:
        df['_recorded_LifeStage'] = df['LifeStage']
    df = df.drop(columns=columns_to_remove)


    #Resolve contradictory fields
    if 'SampleType' in df.columns:
        # blanks, pooled QCs and standards are not samples of the study's subjects/organisms
        mask = df['SampleType'].str.contains('blank', case=False) | \
            df['SampleType'].isin(['pool_QC', 'reference material', 'molecular_standards'])
        strip_for_blank = ['NCBITaxonomy', 'UBERONBodyPartName', 'AgeInYears', 'UniqueSubjectID', 'SubjectIdentifierAsRecorded', 'DOIDCommonName', 'HPOCommonName', 'Exposure', 'TimePointHoursRelativeToEvent', 'FastingStatus', 'ENVOEnvironmentBiome', 'ENVOEnvironmentMaterial', 'BiologicalSex', 'HealthStatus', 'LatitudeandLongitude', 'SmokingStatus', 'BodyMassIndex', 'Diet']
        df.loc[mask, [c for c in strip_for_blank if c in df.columns]] = ''
        if '_recorded_LifeStage' in df.columns:
            df.loc[mask, '_recorded_LifeStage'] = ''

    if 'filename' in df.columns and attempt_adding_file_extensions == True:
        df.loc[df['filename'].str.contains(r'^[^.]+$') & (df['filename'].str.len() >= 1), 'filename'] += '.mzML'
        

    input_columns = set(df.columns)

    # Add missing columns with their respective default missing value from the dictionary or generate values if generate == True
    for key, value in allowedTerm_dict.items():
        if key not in df.columns:
            if value['generate'] == 'False':
                if value['missing'] == 'not allowed':
                    print(f'{key}: MISSING. STOPPING INTERPRETATION!')
                    return pd.DataFrame()
                else:
                    df[key] = value['missing']
                    print(f'{key}: ADDED!')


    for key, value in allowedTerm_dict.items():
        if key in df.columns and value['generate'] == 'False':
            allowed_terms = value['allowed_values']
            missing_value = value['missing']

            # Extract unique values from the DataFrame column
            unique_values = df[key].unique()
            
            if key == 'ComorbidityListDOIDIndex':
                # a list of DOID IDs, not a fixed set of combinations: every ID must be a
                # current (non-obsolete) DOID class of the DOID table
                doid_ids = _cached('doid_ids', doid_ontology_table, lambda: set(
                    doid_ontology_table['DOIDOntologyIndex'].dropna().astype(str).str.strip())
                    if 'DOIDOntologyIndex' in doid_ontology_table.columns else set())
                if not doid_ids and any(str(v).strip() not in ('', missing_value) for v in unique_values):
                    print(f'{key}: no DOID table given; IDs cannot be checked and are set to {missing_value}')
                value_map = {v: normalize_doid_list(v, doid_ids, missing_value) for v in unique_values}
                df[key] = df[key].map(value_map).fillna(missing_value).replace("", missing_value)
                filtered_map = {k: v for k, v in value_map.items() if k != v}
                if filtered_map:
                    print(f"{key}:")
                    for k, v in filtered_map.items():
                        print(f"  {k}: {v}")
                continue

            # Initialize value_map for each condition
            if key not in ['NCBITaxonomy', 'MassSpectrometer', 'filename'] and len(allowed_terms) >= 1 \
                    and allowed_terms[0] not in ['00', 'datetime', 'numeric', 'numeric|numeric']:
                # General case for non-specific keys with a list of allowed terms
                value_map = {
                    observed_value: next((allowed for allowed in allowed_terms if allowed.replace(" ", "") == str(observed_value).replace(" ", "")), missing_value)
                    for observed_value in unique_values
                }
                if case_insensitive_terms:
                    # Fall back to a case/separator-insensitive match when it is unambiguous
                    def _build_canonical(terms=allowed_terms):
                        canonical = {}
                        for allowed in terms:
                            canonical.setdefault(_canonical_term(allowed), set()).add(allowed)
                        return canonical
                    canonical = _cached('canonical', allowed_terms, _build_canonical)
                    for observed_value, mapped in value_map.items():
                        if mapped == missing_value:
                            hits = canonical.get(_canonical_term(observed_value), set())
                            if len(hits) == 1:
                                value_map[observed_value] = next(iter(hits))


            elif key in ['NCBITaxonomy', 'MassSpectrometer']:
                if key == 'NCBITaxonomy':
                    # Build fast lookups from allowed terms: 'ID|Name'
                    # Names are not unique in NCBI (e.g. Bacillus the bacterium and the stick
                    # insect), so a shared name only resolves when the provided ID disambiguates it
                    def _build_ncbi(terms=allowed_terms):
                        allowed_by_id = {}
                        allowed_by_name_ci = {}
                        shared_names_ci = {}
                        for t in terms:
                            if "|" in t:
                                tid, tname = t.split("|", 1)
                                tid = tid.strip()
                                tname = tname.strip().lower()
                                allowed_by_id[tid] = t
                                first = allowed_by_name_ci.setdefault(tname, t)
                                if first != t:
                                    shared_names_ci.setdefault(tname, {first}).add(t)
                        return set(terms), allowed_by_id, allowed_by_name_ci, shared_names_ci
                    allowed_set, allowed_by_id, allowed_by_name_ci, shared_names_ci = _cached('ncbi', allowed_terms, _build_ncbi)
                    division_by_id = {}

                    def name_hit(name, oid):
                        terms = shared_names_ci.get(name.lower())
                        if terms is None:
                            return allowed_by_name_ci.get(name.lower())
                        # Shared name: use the NCBI division of the provided ID to pick one
                        # (e.g. '1409|Bacillus sp.' is in Bacteria -> 1386|Bacillus, not the insect)
                        if not oid:
                            return None
                        if not division_by_id:
                            division_by_id.update(zip(NCBIRankDivision_table['TaxonID'], NCBIRankDivision_table['NCBIDivision']))
                        division = division_by_id.get(oid)
                        same_division = [t for t in terms if division and division_by_id.get(t.split("|", 1)[0].strip()) == division]
                        return same_division[0] if len(same_division) == 1 else None

                    def resolve_ncbi(observed_value: str) -> str:
                        # observed_value like '1931|Streptomyces sp.' or '1883|Streptomyces' or just 'Streptomyces sp.'
                        oid, oname = _split_id_name(observed_value)

                        # 1) Prefer genus mapping if looks like 'Genus sp.' / 'Genus spp.' (even if exact value is allowed)
                        if oname and _is_sp_or_spp(oname):
                            genus = _genus_from_name(oname)
                            hit = name_hit(genus, oid)
                            if hit:
                                return hit  # e.g., '1883|Streptomyces'

                        # 2) Exact allowed term
                        if observed_value in allowed_set:
                            return observed_value

                        # 3) Same ID → allowed term
                        if oid and oid in allowed_by_id:
                            return allowed_by_id[oid]

                        # 4) Same Name (case-insensitive) → allowed term
                        if oname:
                            hit = name_hit(oname, oid)
                            if hit:
                                return hit

                        # 5) Name that NCBI lists for another (current) taxon, e.g. a
                        #    synonym or the name of a merged taxon (GNPS only; optional)
                        if ncbi_name_resolver is not None and oname:
                            hit = ncbi_name_resolver(oname)
                            if hit:
                                return hit

                        # 6) No match
                        return missing_value

                    value_map = {ov: resolve_ncbi(ov) for ov in unique_values}

                else:
                    # MassSpectrometer: keep your original direct matching
                    value_map = {
                        observed_value: (
                            observed_value if observed_value in allowed_terms 
                            else next(
                                (y for y in allowed_terms if "|" in str(observed_value) and y.split("|", 1)[0] == str(observed_value).split("|", 1)[0]),
                                next(
                                    (y for y in allowed_terms if "|" in str(observed_value) and y.split("|", 1)[1] == str(observed_value).split("|", 1)[1]),
                                    missing_value
                                )
                            )
                        )
                        for observed_value in unique_values
                    }
                
            elif not allowed_terms:
                raise ValueError(f"No allowed values for {key}; was the allowed terms json updated from the ontologies?")
            elif allowed_terms[0] == '00':
                # Handling for terms where '00' indicates a placeholder for any value
                value_map = {x: x if pd.notna(x) and x != "" else missing_value for x in unique_values}
            elif allowed_terms[0] == 'datetime':
                # Flexible date/time normalisation to precision-preserving ISO-8601
                # (year / year-month / date / date+time / interval). Ambiguous slash
                # dates (both fields <= 12) are disambiguated per dataset: a value that
                # itself proves an order (a field > 12) always wins; otherwise the
                # dataset's other dates decide; NORMAN (European, ISO) defaults to
                # day-first, everything else (US MassIVE/Workbench) to month-first.
                _slash = re.compile(r'\s*(\d{1,2})[/-](\d{1,2})[/-]\d')

                def _infer_dayfirst(mid, values):
                    if str(mid).startswith('NORMAN'):
                        return True
                    dmy = mdy = False
                    for v in values:
                        m = _slash.match(str(v))
                        if not m:
                            continue
                        a, b = int(m.group(1)), int(m.group(2))
                        if a > 12 and b <= 12:
                            dmy = True
                        elif b > 12 and a <= 12:
                            mdy = True
                    if dmy and not mdy:
                        return True
                    return False  # month-first: US platform default (also for conflicts)

                out = df[key].copy()
                gid = (df['MassiveID'] if 'MassiveID' in df.columns
                       else pd.Series('', index=df.index)).fillna('')
                for mid, idx in out.groupby(gid).groups.items():
                    sub = df.loc[idx, key]
                    dayfirst = _infer_dayfirst(mid, sub.unique())
                    vmap = {
                        x: (normalize_redu_datetime(x, dayfirst=dayfirst) or missing_value)
                        if (pd.notna(x) and x != "") else missing_value
                        for x in sub.unique()
                    }
                    out.loc[idx] = sub.map(vmap)
                df[key] = out.fillna(missing_value).replace("", missing_value)
                continue
            elif allowed_terms[0] == 'numeric':
                # Numeric handling; optionally clamp to a plausible [min, max] range so
                # nonsensical values (e.g. an auto-miscalculated BMI) drop to missing.
                nums = pd.to_numeric(df[key], errors='coerce')
                lo, hi = value.get('min'), value.get('max')
                if lo is not None:
                    nums = nums.where(nums >= lo)
                if hi is not None:
                    nums = nums.where(nums <= hi)
                df[key] = nums.fillna(missing_value).replace("", missing_value)
                continue
            elif allowed_terms[0] == 'numeric|numeric':
                # Handling for numeric ranges split by '|'
                value_map = {
                    x: x if all(part.replace('.', '', 1).isdigit() or part.lstrip('-').replace('.', '', 1).isdigit() for part in str(x).split('|')) else missing_value
                    for x in unique_values
                }
            elif key == 'filename':
                # Specific handling for filenames with allowed extensions
                value_map = {
                    x: x if any(str(x).lower().endswith(ext.lower()) for ext in allowed_terms + other_allowed_file_extensions) 
                    else missing_value 
                    for x in unique_values
                }
            else:
                raise ValueError(f"No validation rule for {key} with allowed values {allowed_terms}")
            # Apply the mapping for columns, except 'numeric' which is handled separately
            if key != 'numeric':
                df[key] = df[key].map(value_map).fillna(missing_value).replace("", missing_value)
            
            # preparing valuemap for print into log
            filtered_map = {k: v for k, v in value_map.items() if k != v}
            if filtered_map: 
                print(f"{key}:")
                for k, v in filtered_map.items():
                    print(f"  {k}: {v}")
        elif key not in df.columns and value['generate'] == 'False':
            df[key] = value['missing']



    for key, value in allowedTerm_dict.items():
        if value['generate'] == 'True':
            missing_value = value['missing']

            if key == 'DataSource':
                df.loc[df['MassiveID'].str.startswith('NORMAN', na=False), 'DataSource'] = 'NORMAN'
                df.loc[df['MassiveID'].str.startswith('MSV', na=False), 'DataSource'] = 'GNPS'
                df.loc[df['MassiveID'].str.startswith('MTBLS', na=False), 'DataSource'] = 'MetaboLights'
                df.loc[df['MassiveID'].str.startswith('ST', na=False), 'DataSource'] = 'Workbench'


            if key == 'LifeStage':
                computed = df.apply(lambda x: age_category(x['AgeInYears']) if x['NCBITaxonomy'] == "9606|Homo sapiens" else value['missing'], axis=1)
                computed = computed.replace('', value['missing'])  # age_category gives '' for non-numeric ages
                # a life stage recorded by the submitter is kept when the age gives none
                if '_recorded_LifeStage' in df.columns and len(value.get('allowed_values') or []) > 0:
                    recorded = df['_recorded_LifeStage'].where(df['_recorded_LifeStage'].isin(value['allowed_values']))
                    computed = computed.where(computed != value['missing'], recorded).fillna(value['missing'])
                df[key] = computed
            if key == 'UniqueSubjectID':
                df[key] = df.apply(lambda x: str(x['MassiveID']) + '_' + str(x['SubjectIdentifierAsRecorded']) 
                                if x['SubjectIdentifierAsRecorded'] != allowedTerm_dict['SubjectIdentifierAsRecorded']['missing'] 
                                and x['SubjectIdentifierAsRecorded'] != ''
                                and x['SubjectIdentifierAsRecorded'] != 'not applicable' 
                                and not pd.isna(x['SubjectIdentifierAsRecorded']) 
                                else value['missing'], axis=1)
            if key == 'UBERONOntologyIndex':
                df = df.merge(uberon_ontology_table[['Label', 'UBERONOntologyIndex', 'Is Multicellular','Is Organ', 'Is Fluid']].drop_duplicates('Label'), left_on='UBERONBodyPartName', right_on='Label', how='left')
                # QC, blank and standard rows keep their own SampleTypeSub1 (e.g. pool_QC of plasma)
                sample_type = df['SampleType'] if 'SampleType' in df.columns else pd.Series('', index=df.index)
                not_qc = ~sample_type.astype(str).str.match(r'blank_|pool_QC|molecular_standards|reference material')
                # only where Sub1 is still open; finer Sub1 values (insect, plant_angiospermae, ...) are kept
                sub1 = df['SampleTypeSub1'] if 'SampleTypeSub1' in df.columns else pd.Series(pd.NA, index=df.index, dtype=object)
                sub1_open = sub1.isna() | sub1.isin(['', value['missing'], allowedTerm_dict['SampleTypeSub1']['missing'], 'tissue', 'biofluid'])
                df.loc[(df['Is Multicellular'] == True) & not_qc & sub1_open, 'SampleTypeSub1'] = 'tissue'
                df.loc[(df['Is Fluid'] == True) & not_qc & sub1_open, 'SampleTypeSub1'] = 'biofluid'
                df.drop(columns=['Label', 'Is Multicellular', 'Is Organ', 'Is Fluid'], inplace=True)
            if key == 'DOIDOntologyIndex':
                df = df.merge(doid_ontology_table[['Label', 'DOIDOntologyIndex']].drop_duplicates('Label'), left_on='DOIDCommonName', right_on='Label', how='left')
                df.drop(columns=['Label'], inplace=True)
            if key == 'HPOOntologyIndex':
                if 'HPOCommonName' not in df.columns:
                    df['HPOCommonName'] = allowedTerm_dict.get('HPOCommonName', {}).get('missing', missing_value)
                df = df.merge(hpo_ontology_table[['Label', 'HPOOntologyIndex']].drop_duplicates('Label'), left_on='HPOCommonName', right_on='Label', how='left')
                df.drop(columns=['Label'], inplace=True)
            if key == 'ExposureIndex':
                if 'Exposure' not in df.columns:
                    df['Exposure'] = allowedTerm_dict.get('Exposure', {}).get('missing', missing_value)
                # PECO is for plants: on other samples a PECO term becomes the ECTO term with
                # that name as label or exact synonym (e.g. PECO 'radiation exposure' ->
                # ECTO 'exposure to radiation'), or missing if ECTO has none
                if 'Synonym' in exposure_ontology_table.columns and 'SampleType' in df.columns:
                    ecto = exposure_ontology_table[exposure_ontology_table['ExposureIndex'].astype(str).str.startswith('ECTO:')]
                    ecto_by_name = {}
                    for label, synonym in zip(ecto['Label'].astype(str), ecto['Synonym'].astype(str)):
                        for name in (label, synonym):
                            ecto_by_name.setdefault(name.strip().lower(), set()).add(label)
                    peco_labels = set(exposure_ontology_table.loc[exposure_ontology_table['ExposureIndex'].astype(str).str.startswith('PECO:'), 'Label'].astype(str)) \
                        - set(ecto['Label'].astype(str))
                    not_plant = ~df['SampleType'].astype(str).str.lower().str.startswith('plant')
                    on_peco = not_plant & df['Exposure'].isin(peco_labels)
                    def _to_ecto(value):
                        labels = ecto_by_name.get(str(value).strip().lower(), set())
                        return next(iter(labels)) if len(labels) == 1 else allowedTerm_dict['Exposure']['missing']
                    df.loc[on_peco, 'Exposure'] = df.loc[on_peco, 'Exposure'].map(_to_ecto)
                # ECTO terms also give what the organism was exposed to (OrganismExposureAgent)
                agent_columns = [c for c in ('OrganismExposureAgentIndex', 'OrganismExposureAgent')
                                 if c in exposure_ontology_table.columns and c in allowedTerm_dict]
                df = df.drop(columns=[c for c in agent_columns if c in df.columns])
                df = df.merge(exposure_ontology_table[['Label', 'ExposureIndex'] + agent_columns].drop_duplicates('Label'), left_on='Exposure', right_on='Label', how='left')
                df.drop(columns=['Label'], inplace=True)
                for c in ('OrganismExposureAgentIndex', 'OrganismExposureAgent'):
                    if c in allowedTerm_dict and c not in df.columns:
                        df[c] = allowedTerm_dict[c]['missing']
            if key == 'ENVOEnvironmentBiomeIndex':
                df = df.merge(envome_biome_ontology_table[['Label', 'ENVOEnvironmentBiomeIndex']].drop_duplicates('Label'), left_on='ENVOEnvironmentBiome', right_on='Label', how='left')
                df.drop(columns=['Label'], inplace=True)
            if key == 'ENVOEnvironmentMaterialIndex':
                df = df.merge(envome_material_ontology_table[['Label', 'ENVOEnvironmentMaterialIndex']].drop_duplicates('Label'), left_on='ENVOEnvironmentMaterial', right_on='Label', how='left')
                df.drop(columns=['Label'], inplace=True)
            if key == 'ENVOMediumScaleIndex':
                df = df.merge(envome_medium_scale_ontology_table[['Label', 'ENVOMediumScaleIndex']].drop_duplicates('Label'), left_on='ENVOMediumScale', right_on='Label', how='left')
                df.drop(columns=['Label'], inplace=True)
            if key == 'ENVOLocalScaleIndex':
                df = df.merge(envome_local_scale_ontology_table[['Label', 'ENVOLocalScaleIndex']].drop_duplicates('Label'), left_on='ENVOLocalScale', right_on='Label', how='left')
                df.drop(columns=['Label'], inplace=True)
            if key == 'ENVOBroadScaleIndex':
                df = df.merge(envome_broad_scale_ontology_table[['Label', 'ENVOBroadScaleIndex']].drop_duplicates('Label'), left_on='ENVOBroadScale', right_on='Label', how='left')
                df.drop(columns=['Label'], inplace=True)
            if key == 'USI' and add_usi == True:
                df['filename'] = df['filename'].apply(process_filename)
                df['USI'] = 'mzspec:' + df['MassiveID'] + ':' + df['filename']
            if key == 'NCBIRank':
                df['TaxonID'] = df['NCBITaxonomy'].apply(lambda x: x.split('|')[0] if '|' in x else 'missing value')
                df['TaxonID'] = df['TaxonID'].astype(str)
                df = pd.merge(df, NCBIRankDivision_table, on='TaxonID', how='left')
                df.drop(columns=['TaxonID'], inplace=True)

                if 'NCBIRank' in df.columns:
                    df['NCBIRank'] = df['NCBIRank'].fillna(missing_value)
                else:
                    df['NCBIRank'] = missing_value

                if 'NCBIDivision' in df.columns:
                    df['NCBIDivision'] = df['NCBIDivision'].fillna(missing_value)
                else:
                    df['NCBIDivision'] = missing_value



    # generated columns without their own rule (e.g. OrganismExposureAgent without an
    # exposure table) still have to exist
    for key, value in allowedTerm_dict.items():
        if value.get('generate') == 'True' and key not in df.columns and key not in ('USI',):
            if key in ('OrganismExposureAgent', 'OrganismExposureAgentIndex'):
                df[key] = value['missing']

    # ontology indexes of unmatched labels come out of the merges as NaN
    for key, value in allowedTerm_dict.items():
        if key.endswith('OntologyIndex') or key in ('ExposureIndex', 'OrganismExposureAgent', 'OrganismExposureAgentIndex') \
                or (key.startswith('ENVO') and key.endswith('Index')):
            if key in df.columns:
                df[key] = df[key].fillna(value['missing']).replace('', value['missing'])

    # Ensure the dataframe contains only the columns specified in the dictionary
    keys_to_include = [key for key in allowedTerm_dict.keys() if key != 'USI' or add_usi or keep_usi]

    ignored_columns = input_columns - set(keys_to_include)
    
    for column in ignored_columns:
        print(f"IGNORED:  {column}")

    #remove rows if not enough metadata are present
    def count_non_missing_specific(row, cols):
        return (row[cols] != "missing value").sum()
    
    columns_to_consider = ["SampleType", "SampleTypeSub1", "NCBITaxonomy", "UBERONBodyPartName", "BiologicalSex", 
                           "AgeInYears", "LifeStage", "Country", "HealthStatus", "SampleExtractionMethod", 
                           "SampleCollectionMethod", "ComorbidityListDOIDIndex", "DOIDCommonName", "HPOCommonName", "Exposure", 
                           "DepthorAltitudeMeters", "HumanPopulationDensity", "LatitudeandLongitude",
                           "ENVOEnvironmentBiome", "ENVOEnvironmentMaterial"]
    
    original_row_count = df.shape[0]
    columns_to_check = [col for col in columns_to_consider if col in df.columns]
    df = df[(df.apply(count_non_missing_specific, cols=columns_to_check, axis=1) >= 1)]

    print(f"Number of rows removed due to not enough metadata: {original_row_count - df.shape[0]}")
    print(f"Returning {len(df)} rows!")

    return df[keys_to_include]

# Values submitters commonly use for a term the vocabulary spells differently
# (keys are lower case). Only mappings that do not change the meaning.
_VALUE_SYNONYMS = {
    'Country': {
        'usa': 'United States of America', 'us': 'United States of America', 'u.s.a.': 'United States of America',
        'u.s.': 'United States of America', 'united states': 'United States of America',
        'uk': 'United Kingdom', 'u.k.': 'United Kingdom', 'great britain': 'United Kingdom',
        'england': 'United Kingdom', 'scotland': 'United Kingdom', 'wales': 'United Kingdom',
        'northern ireland': 'United Kingdom',
    },
    'IonizationSourceAndPolarity': {
        'heated electrospray ionization (positive)': 'electrospray ionization (positive)',
        'heated electrospray ionization (negative)': 'electrospray ionization (negative)',
        'heated electrospray ionization (alternating)': 'electrospray ionization (alternating)',
    },
    'ChromatographyAndPhase': {
        'reverse phase (polar c18)': 'reverse phase (C18)',
    },
    # first-generation template values
    'SampleType': {
        'standard': 'molecular_standards', 'culture_mamalian': 'culture_mammalian',
        'deceased animal': 'animal', 'blank': 'blank_analysis', 'pool': 'pool_QC', 'qc': 'pool_QC',
        'animal_mammalia_nos': 'animal', 'human': 'animal',
    },
    'SmokingStatus': {
        'nonsmoker': 'non-smoker', 'non smoker': 'non-smoker',
        'ex-smoker': 'former smoker', 'exsmoker': 'former smoker', 'ex smoker': 'former smoker',
    },
}

# Prefixes of the first-generation ReDU template (e.g. 'ATTRIBUTE_Curated_SampleType',
# 'Subject_BiologicalSex', 'ATTRIBUTE_Analysis_MassSpectrometer') in front of current
# field names, and first-generation names that differ from the current ones (lower case,
# after the prefixes). Same rules as template_column in gnps_downloader.py.
_LEGACY_COLUMN_PREFIXES = ('ATTRIBUTE_Curated_', 'ATTRIBUTE_', 'Subject_', 'Analysis_', 'Study_', 'Curated_')
_LEGACY_COLUMN_ALIASES = {
    'sex': 'BiologicalSex', 'health': 'HealthStatus', 'bodypartontologyname': 'UBERONBodyPartName',
    'bmi': 'BodyMassIndex', 'diseasecommonname': 'DOIDCommonName', 'smoking': 'SmokingStatus',
    'sampletype_sub1': 'SampleTypeSub1', 'comorbidity_listdoids': 'ComorbidityListDOIDIndex',
    'massiveid': 'MassiveID', 'datasetaccession': 'MassiveID',
}


def _template_column(column, allowedTerm_dict):
    """The current field a (first-generation or differently cased) column name stands for."""
    by_lower = _cached('keys_lower', allowedTerm_dict, lambda: {k.lower(): k for k in allowedTerm_dict})
    name = str(column).strip()
    while True:
        if name in allowedTerm_dict:
            return name
        if name.lower() in by_lower:
            return by_lower[name.lower()]
        if name.lower() in _LEGACY_COLUMN_ALIASES:
            return _LEGACY_COLUMN_ALIASES[name.lower()]
        for prefix in _LEGACY_COLUMN_PREFIXES:
            if name.lower().startswith(prefix.lower()) and len(name) > len(prefix):
                name = name[len(prefix):]
                break
        else:
            return None


_LATLON_NUMBER = r'([+-]?\d{1,3}(?:\.\d+)?)'
_LATLON = re.compile(r'^' + _LATLON_NUMBER + r'\s*°?\s*([NnSs])?\s*(AND|and|And|[;,|/]|\s)\s*'
                     + _LATLON_NUMBER + r'\s*°?\s*([EeWw])?$')


def _is_redu_latlon(value):
    """The format the ReDU validation accepts ('lat|lon', numbers without a '+')."""
    parts = str(value).split('|')
    return len(parts) == 2 and all(re.fullmatch(r'-?\d+(\.\d+)?', p) for p in parts)


def parse_latitude_longitude(value):
    """'lat|lon' from '32.842;-117.258', '48.2,-124.3', '32.842AND-117.258', '-33.99|+25.70',
    '41.8719 N, 12.5674 E' or '48.27 -124.45' (latitude first; N/S/E/W give the sign and
    exclude a minus sign). The numbers keep their digits. None if the value is not such a
    pair or out of range."""
    m = _LATLON.match(str(value).strip())
    if not m:
        return None
    lat, ns, separator, lon, ew = m.groups()
    if separator == ',' and '.' not in lat + lon:
        return None  # '41,87' may be one number with a decimal comma
    out = []
    for number, hemisphere, negative, limit in ((lat, ns, 'sS', 90), (lon, ew, 'wW', 180)):
        number = number.lstrip('+')
        if hemisphere:
            if number.startswith('-'):
                return None
            if hemisphere in negative:
                number = '-' + number
        if abs(float(number)) > limit:
            return None
        out.append(number)
    return '|'.join(out)


def _load_ms_aliases():
    path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'translation_sheets', 'GNPS_MassSpectrometer_aliases.csv')
    if not os.path.exists(path):
        return {}
    table = pd.read_csv(path, dtype=str).dropna()
    return {g.strip().lower(): r.strip() for g, r in zip(table['GNPS'], table['REDU'])}


def _diet_key(value):
    return re.sub(r'[\s_\-]+', '', str(value).strip().lower())


def _solvent_key(value):
    return str(value).strip().lower()

_MISSING_LIKE = {'', 'nan', 'na', 'n/a', 'none', 'missing value', 'missing_value', 'not specified', 'not collected',
                 'not applicable', 'not_applicable', 'ml import: not available', 'not available'}


def _is_missing_like(value):
    return str(value).strip().lower() in _MISSING_LIKE


def normalize_gnps_metadata(df, allowedTerm_dict, dataset_accession=None, uberon_synonym_table=None,
                            envo_tables=(), NCBIRankDivision_table=None):
    """Repairs common, unambiguous deviations in submitted GNPS/ReDU metadata files before
    validation, so the values are kept instead of being replaced by 'missing value'."""
    df = df.copy()

    # first-generation template column names (and case variants) -> current names; when
    # the current column exists too, its empty values are taken from the old one
    for column in list(df.columns):
        if column in allowedTerm_dict:
            continue
        name = _template_column(column, allowedTerm_dict)
        if name is None or name == column:
            continue
        if name not in df.columns:
            df = df.rename(columns={column: name})
            print(f"Renamed legacy column {column} -> {name}")
        else:
            empty = df[name].map(_is_missing_like) & ~df[column].map(_is_missing_like)
            if empty.any():
                df.loc[empty, name] = df.loc[empty, column]
                print(f"Filled {int(empty.sum())} empty {name} values from legacy column {column}")
            df = df.drop(columns=[column])

    # files without a dataset column: the accession is known from the file name
    if 'MassiveID' not in df.columns and dataset_accession:
        df['MassiveID'] = dataset_accession
        print(f"Added MassiveID {dataset_accession} from the file name")

    def _ci_lookup(mapping, value):
        return mapping.get(str(value).strip().lower())

    for key, synonyms in _VALUE_SYNONYMS.items():
        if key in df.columns:
            df[key] = df[key].map(lambda v: _ci_lookup(synonyms, v) or v)

    # SampleTypeSub1 values (e.g. biofluid, tissue, reference material_chemicalstandard)
    # entered as SampleType: move them, and infer SampleType only when unambiguous
    if 'SampleType' in df.columns:
        sample_types = set(allowedTerm_dict['SampleType']['allowed_values'])
        sub1_types = set(allowedTerm_dict['SampleTypeSub1']['allowed_values'])
        if 'SampleTypeSub1' not in df.columns:
            df['SampleTypeSub1'] = ''
        is_sub1 = df['SampleType'].isin(sub1_types - sample_types) & df['SampleTypeSub1'].map(_is_missing_like)
        if is_sub1.any():
            moved = df.loc[is_sub1, 'SampleType']
            df.loc[is_sub1, 'SampleTypeSub1'] = moved
            prefix = moved.str.rsplit('_', n=1).str[0]
            inferred = prefix.where(prefix.isin(sample_types) & (prefix != moved), '')
            # biofluid/tissue: animal when the organism is an animal
            if NCBIRankDivision_table is not None and 'NCBITaxonomy' in df.columns:
                division = _cached('division', NCBIRankDivision_table, lambda: dict(zip(
                    NCBIRankDivision_table['TaxonID'].astype(str), NCBIRankDivision_table['NCBIDivision'])))
                animal_divisions = {'Primates', 'Rodents', 'Mammals', 'Vertebrates', 'Invertebrates'}
                taxid = df.loc[is_sub1, 'NCBITaxonomy'].astype(str).str.split('|').str[0].str.strip()
                is_animal = taxid.map(division).isin(animal_divisions)
                inferred = inferred.mask(moved.isin(['biofluid', 'tissue']) & is_animal, 'animal')
            df.loc[is_sub1, 'SampleType'] = inferred
            print(f"Moved {is_sub1.sum()} SampleTypeSub1 values out of SampleType")

    # 'healthy' entered as a disease
    if 'DOIDCommonName' in df.columns:
        healthy = df['DOIDCommonName'].astype(str).str.strip().str.lower() == 'healthy'
        if healthy.any():
            if 'HealthStatus' not in df.columns:
                df['HealthStatus'] = ''
            no_status = df['HealthStatus'].map(_is_missing_like)
            df.loc[healthy & no_status, 'HealthStatus'] = 'healthy'
            df.loc[healthy, 'DOIDCommonName'] = ''

    # 'disease NOS' says the subject is ill without naming the disease
    if 'DOIDCommonName' in df.columns:
        nos = df['DOIDCommonName'].astype(str).str.strip().str.lower() == 'disease nos'
        if nos.any():
            if 'HealthStatus' not in df.columns:
                df['HealthStatus'] = ''
            fill = nos & df['HealthStatus'].map(_is_missing_like)
            df.loc[fill, 'HealthStatus'] = 'unhealthy (NOS)'
            if fill.any():
                print(f"HealthStatus 'unhealthy (NOS)' for {int(fill.sum())} rows with DOIDCommonName 'disease NOS'")

    # 'lat | lon' with spaces, Excel float years
    if 'LatitudeandLongitude' in df.columns:
        df['LatitudeandLongitude'] = df['LatitudeandLongitude'].astype(str).str.replace(r'\s*\|\s*', '|', regex=True).str.strip()
        # other separators / hemisphere letters (values already in ReDU format are kept)
        def _latlon(v):
            if _is_redu_latlon(v) or _is_missing_like(v):
                return v
            return parse_latitude_longitude(v) or v
        df['LatitudeandLongitude'] = df['LatitudeandLongitude'].map(_latlon)

    # MassSpectrometer given as a bare MS accession, as the bare label of an allowed term
    # or as a first-generation name ('Maxis_ImpactHD'; translation_sheets/GNPS_MassSpectrometer_aliases.csv)
    if 'MassSpectrometer' in df.columns and 'MassSpectrometer' in allowedTerm_dict:
        ms_terms = allowedTerm_dict['MassSpectrometer']['allowed_values']

        def _build_ms():
            by_id, by_label = {}, {}
            for t in ms_terms:
                label, _, acc = str(t).partition('|')
                by_id[acc.strip().upper()] = t
                by_label.setdefault(label.strip().lower(), set()).add(t)
            return set(ms_terms), by_id, by_label, _load_ms_aliases()
        ms_set, ms_by_id, ms_by_label, ms_aliases = _cached('ms_lookup', ms_terms, _build_ms)

        def _ms(v):
            if v in ms_set or _is_missing_like(v):
                return v
            text = str(v).strip()
            if re.fullmatch(r'(?i)MS[:_]\d{7}', text):
                return ms_by_id.get('MS:' + text[3:], v)
            label = ms_by_label.get(text.lower())
            if label and len(label) == 1:
                return next(iter(label))
            alias = ms_aliases.get(text.lower())
            return alias if alias in ms_set else v
        df['MassSpectrometer'] = df['MassSpectrometer'].map(_ms)

    # SampleExtractionMethod: a bare solvent is the 100% term ('methanol' -> 'methanol (100%)'),
    # '<solvent> (80%)' the solvent-water ratio ('acetonitrile-water (4:1)'); only when that
    # term is allowed
    if 'SampleExtractionMethod' in df.columns and 'SampleExtractionMethod' in allowedTerm_dict:
        sem_terms = allowedTerm_dict['SampleExtractionMethod']['allowed_values']
        sem_set = set(sem_terms)

        def _sem(v):
            if v in sem_set or _is_missing_like(v):
                return v
            text = _solvent_key(v)
            if re.fullmatch(r'[a-z][a-z \-]*[a-z]', text) and f'{text} (100%)' in sem_set:
                return f'{text} (100%)'
            m = re.fullmatch(r'([a-z][a-z \-]*[a-z])\s*\((\d{1,2})\s*%\)', text)
            if m and m.group(1) != 'water':
                part = int(m.group(2))
                g = math.gcd(part, 100 - part)
                term = f'{m.group(1)}-water ({part // g}:{(100 - part) // g})'
                if term in sem_set:
                    return term
            return v
        df['SampleExtractionMethod'] = df['SampleExtractionMethod'].map(_sem)

    # Diet: spelling variants that differ only in case, spaces, '_' or '-' ('High_Fat_Diet')
    if 'Diet' in df.columns and 'Diet' in allowedTerm_dict:
        diet_terms = allowedTerm_dict['Diet']['allowed_values']

        def _build_diet():
            keys = {}
            for t in diet_terms:
                keys.setdefault(_diet_key(t), set()).add(t)
            return keys
        diet_keys = _cached('diet_keys', diet_terms, _build_diet)
        diet_set = set(diet_terms)

        def _diet(v):
            if v in diet_set or _is_missing_like(v):
                return v
            hits = diet_keys.get(_diet_key(v), set())
            return next(iter(hits)) if len(hits) == 1 else v
        df['Diet'] = df['Diet'].map(_diet)
    if 'YearOfAnalysis' in df.columns:
        df['YearOfAnalysis'] = df['YearOfAnalysis'].astype(str).str.replace(r'^\s*(\d{4})\.0+\s*$', r'\1', regex=True)

    # UBERON: exact synonyms and case variants of labels -> label
    if uberon_synonym_table is not None and 'UBERONBodyPartName' in df.columns:
        labels = _cached('uberon_labels', uberon_synonym_table,
                         lambda: set(allowedTerm_dict['UBERONBodyPartName']['allowed_values']))

        def _build_to_labels():
            to_labels = {}
            for label, synonym in zip(uberon_synonym_table['Label'].astype(str), uberon_synonym_table['Synonym'].astype(str)):
                if label in labels:
                    to_labels.setdefault(label.strip().lower(), set()).add(label)
                    # short synonyms are abbreviations ('NA' = nasal artery) that collide
                    # with placeholders, so only longer synonyms are used
                    if len(synonym.strip()) > 3:
                        to_labels.setdefault(synonym.strip().lower(), set()).add(label)
            return to_labels
        to_labels = _cached('uberon_synonyms', uberon_synonym_table, _build_to_labels)
        df['UBERONBodyPartName'] = df['UBERONBodyPartName'].map(
            lambda v: v if v in labels or _is_missing_like(v) else (
                next(iter(to_labels[str(v).strip().lower()]))
                if len(to_labels.get(str(v).strip().lower(), ())) == 1 else v))

    # ENVO fields given as ontology IDs (ENVO:00000446) -> label
    envo_by_id = {}
    for table in envo_tables:
        index_column = [c for c in table.columns if c.endswith('Index')][0]
        for index, label in zip(table[index_column].astype(str), table['Label'].astype(str)):
            envo_by_id[index.replace('_', ':')] = label
    for key in ['ENVOEnvironmentBiome', 'ENVOEnvironmentMaterial', 'ENVOBroadScale', 'ENVOLocalScale', 'ENVOMediumScale']:
        if key in df.columns:
            df[key] = df[key].map(lambda v: envo_by_id.get(str(v).strip(), v))

    return df


_NCBI_NAME_CLASSES = ('scientific name', 'synonym', 'equivalent name')


def build_ncbi_name_resolver(names_dmp, names, allowed_ncbi_terms, NCBIRankDivision_table=None):
    """Resolver for NCBI names that are not the scientific name of an allowed taxon:
    - a synonym / equivalent name / scientific name of exactly one taxon (in names.dmp)
      that is allowed, e.g. 'Clostridium coccoides' -> '33035|Blautia producta', including
      names of merged taxa ('37480|Streptomyces pluricolorescens');
    - a strain or other name below species ('Staphylococcus aureus JE3') that NCBI does not
      list at all -> its species, when the first two words are the name of exactly one
      allowed taxon of rank species. Hybrids ('A b x c') are not reduced.
    Only names written like a scientific name are considered (at least two words, the
    first a capitalised genus: not placeholders such as 'not specified', which NCBI lists
    as a synonym of 32644 'unidentified'). names: the names that occur in the metadata, so
    only those are kept from names.dmp."""
    def scientific_like(name):
        name = str(name).strip()
        return re.match(r'^[A-Z][a-z]+(\s+\S+)+$', name) is not None and not _is_missing_like(name)
    wanted = {str(n).strip().lower() for n in names if scientific_like(n)}
    if not wanted:
        return None
    resolving = collections.defaultdict(set)
    known = set()
    with open(names_dmp, encoding='utf-8', errors='replace') as f:
        for line in f:
            parts = line.split('\t|\t')
            if len(parts) < 4:
                continue
            name = parts[1].strip().lower()
            if name not in wanted:
                continue
            known.add(name)
            name_class = parts[3].rstrip('\t|\n').strip()
            if name_class in _NCBI_NAME_CLASSES:
                resolving[name].add(parts[0].strip())
    allowed_by_id = {}
    allowed_by_name = collections.defaultdict(set)
    for term in allowed_ncbi_terms:
        if '|' in term:
            tid, tname = term.split('|', 1)
            allowed_by_id[tid.strip()] = term
            allowed_by_name[tname.strip().lower()].add(term)
    rank = {}
    if NCBIRankDivision_table is not None:
        rank = dict(zip(NCBIRankDivision_table['TaxonID'].astype(str), NCBIRankDivision_table['NCBIRank'].astype(str)))
    print(f"NCBI names.dmp: {len(resolving)} of {len(wanted)} metadata names are NCBI names")

    def resolve(name):
        name = str(name).strip()
        key = name.lower()
        if key not in wanted or not scientific_like(name):
            return None
        ids = resolving.get(key)
        if ids:
            taxid = next(iter(ids))
            # never to the catch-all 'unidentified' taxon
            return allowed_by_id.get(taxid) if len(ids) == 1 and taxid != '32644' else None
        if key in known:
            return None  # listed by NCBI, but not as a name of one taxon
        words = name.split()
        if len(words) < 3 or re.search(r'(?i)\sx\s|×', name):
            return None
        species = allowed_by_name.get(' '.join(words[:2]).lower(), set())
        if len(species) != 1:
            return None
        term = next(iter(species))
        return term if rank.get(term.split('|', 1)[0].strip()) == 'species' else None
    return resolve


def _ncbi_names_in_files(paths, allowedTerm_dict):
    """Organism names in the NCBITaxonomy columns (any template name) of metadata files."""
    names = set()
    for path in paths:
        try:
            header = pd.read_csv(path, sep='\t', dtype=str, nrows=0).columns
            columns = [c for c in header if _template_column(c, allowedTerm_dict) == 'NCBITaxonomy']
            if not columns:
                continue
            table = pd.read_csv(path, sep='\t', dtype=str, usecols=columns)
        except Exception:
            continue
        for column in columns:
            for value in table[column].dropna().unique():
                names.add(str(value).split('|', 1)[-1].strip())
    return names


def process_filename(filename):
    path_parts = filename.split('/')
    for part in path_parts:
        print(part)
        if part.endswith('.d'):
            return '/'.join(path_parts[:path_parts.index(part)+1])
    return filename


def clean_read_tsv_quoted_lines(path):
    cleaned_data = ''
    with open(path, 'r') as file:
        for line in file:
            # Strip leading and trailing whitespace and then strip double quotes
            cleaned_line = line.strip().strip('"')
            cleaned_data += cleaned_line + '\n'

    # Use StringIO to simulate a file object for pandas
    cleaned_data_io = StringIO(cleaned_data)

    # Now read the cleaned data into a DataFrame
    df = pd.read_csv(cleaned_data_io, sep='\t')

    return df

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='GNPS Validator')
    parser.add_argument('path_to_github_metadata')    
    parser.add_argument('output_metadata_folder')
    parser.add_argument("--AllowedTermJson_path", type=str, help="Path to json with allowed terms")
    parser.add_argument('--path_to_uberon_cl_po_csv')
    parser.add_argument('--path_to_envo_biome_csv')
    parser.add_argument('--path_to_envo_material_csv')
    parser.add_argument('--path_ncbi_rank_division')
    parser.add_argument('--path_to_doid_csv')
    parser.add_argument('--path_to_hpo_csv', default='none')
    parser.add_argument('--path_to_exposure_csv', default='none')
    parser.add_argument('--path_to_ncbi_names_dmp', default='none',
                        help='NCBI names.dmp; resolves synonyms, merged taxa and strains of NCBITaxonomy names')
    args = parser.parse_args()

    with open(args.AllowedTermJson_path, 'r') as json_file:
        allowed_terms = json.load(json_file)

    uberon_ontology_table = pd.read_csv(args.path_to_uberon_cl_po_csv, index_col=False)
    uberon_ontology_table = uberon_ontology_table.drop_duplicates(subset=['Label'])

    doid_ontology_table = pd.read_csv(args.path_to_doid_csv, index_col=False)
    # HPO (HPOCommonName) and ECTO/PECO (Exposure) tables; optional
    hpo_ontology_table = pd.read_csv(args.path_to_hpo_csv, dtype=str) if args.path_to_hpo_csv != 'none' else None
    exposure_ontology_table = pd.read_csv(args.path_to_exposure_csv, dtype=str) if args.path_to_exposure_csv != 'none' else None
    doid_ontology_table = doid_ontology_table.drop_duplicates(subset=['Label'])

    ENVOEnvironmentBiomeIndex_table = pd.read_csv(args.path_to_envo_biome_csv, index_col=False)
    ENVOEnvironmentBiomeIndex_table = ENVOEnvironmentBiomeIndex_table.drop_duplicates(subset=['Label'])

    ENVOEnvironmentMaterialIndex_table = pd.read_csv(args.path_to_envo_material_csv, index_col=False)
    ENVOEnvironmentMaterialIndex_table = ENVOEnvironmentMaterialIndex_table.drop_duplicates(subset=['Label'])

    NCBIRankDivision_table = pd.read_csv(args.path_ncbi_rank_division, index_col = False)
    NCBIRankDivision_table = NCBIRankDivision_table.drop_duplicates(subset=['TaxonID'])

    # all synonyms (the table above keeps one row per label)
    uberon_synonym_table = pd.read_csv(args.path_to_uberon_cl_po_csv, index_col=False, dtype=str)


    print('Starting tsv processing.')
    #for file_path in glob.glob(f"{args.path_to_github_metadata}/redu_*.tsv"):

    all_metadata_files = glob.glob(f"{args.path_to_github_metadata}/*.tsv")

    print("Processing", len(all_metadata_files), "files")

    ncbi_name_resolver = None
    if args.path_to_ncbi_names_dmp != 'none':
        ncbi_name_resolver = build_ncbi_name_resolver(
            args.path_to_ncbi_names_dmp, _ncbi_names_in_files(all_metadata_files, allowed_terms),
            allowed_terms['NCBITaxonomy']['allowed_values'], NCBIRankDivision_table)
   

    for file_path in all_metadata_files:
    #for file_path in ["/home/yasin/projects/ReDU_metadata/metadata/redu_MSV000093329.tsv"]:
        
        #if re.match(r'.*redu_MSV\d+\.tsv$', file_path):
        if True:

            print(f"Processing: {file_path}")

            # Read as strings so e.g. a YearOfAnalysis column with blanks is not parsed as
            # float ('2019.0' would not match the allowed '2019')
            df = pd.read_csv(file_path, sep='\t', dtype=str)

            if 'ATTRIBUTE_MassiveID' in df.columns:
                df.rename(columns={'ATTRIBUTE_MassiveID': 'MassiveID'}, inplace=True)
            if 'ATTRIBUTE_DatasetAccession' in df.columns:
                df.rename(columns={'ATTRIBUTE_DatasetAccession': 'MassiveID'}, inplace=True)

            # GitHub files are named <MSV>.tsv, MassIVE downloads <MSV>_<n>_gnps_metadata.tsv
            accession = re.match(r'(MSV\d+)', os.path.basename(file_path))
            df = normalize_gnps_metadata(df, allowed_terms,
                                         dataset_accession=accession.group(1) if accession else None,
                                         uberon_synonym_table=uberon_synonym_table,
                                         envo_tables=(ENVOEnvironmentBiomeIndex_table, ENVOEnvironmentMaterialIndex_table),
                                         NCBIRankDivision_table=NCBIRankDivision_table)

            #generate extra columns, add missing columns and remove values which are not in allowed terms
            df = complete_and_fill_REDU_table(df, 
                                              allowed_terms, 
                                              UBERONOntologyIndex_table=uberon_ontology_table, 
                                              DOIDOntologyIndex_table=doid_ontology_table, 
                                              HPOOntologyIndex_table=hpo_ontology_table,
                                              ExposureIndex_table=exposure_ontology_table,
                                              NCBIRankDivision_table=NCBIRankDivision_table,
                                              ENVOEnvironmentBiomeIndex_table=ENVOEnvironmentBiomeIndex_table,
                                              ENVOEnvironmentMaterialIndex_table=ENVOEnvironmentMaterialIndex_table,
                                              attempt_adding_file_extensions=True,
                                              case_insensitive_terms=True,
                                              NCBI_name_resolver=ncbi_name_resolver)
            
            if len(df) > 0:
            
                file_name = os.path.join(args.output_metadata_folder, os.path.basename(file_path))
                df.to_csv(file_name, sep='\t', index=False)
