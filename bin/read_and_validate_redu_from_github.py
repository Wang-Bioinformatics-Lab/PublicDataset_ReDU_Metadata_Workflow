import argparse
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
    df = df.drop(columns=columns_to_remove)


    #Resolve contradictory fields
    if 'SampleType' in df.columns:
        mask = df['SampleType'].str.contains('blank', case=False)
        strip_for_blank = ['NCBITaxonomy', 'UBERONBodyPartName', 'AgeInYears', 'UniqueSubjectID', 'DOIDCommonName', 'ENVOEnvironmentBiome', 'ENVOEnvironmentMaterial', 'BiologicalSex', 'HealthStatus', 'LatitudeandLongitude', 'SmokingStatus', 'BodyMassIndex', 'Diet']
        df.loc[mask, [c for c in strip_for_blank if c in df.columns]] = ''

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

                        # 5) No match
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
                df[key] = df.apply(lambda x: age_category(x['AgeInYears']) if x['NCBITaxonomy'] == "9606|Homo sapiens" else value['missing'], axis=1)
            if key == 'UniqueSubjectID':
                df[key] = df.apply(lambda x: str(x['MassiveID']) + '_' + str(x['SubjectIdentifierAsRecorded']) 
                                if x['SubjectIdentifierAsRecorded'] != allowedTerm_dict['SubjectIdentifierAsRecorded']['missing'] 
                                and x['SubjectIdentifierAsRecorded'] != ''
                                and x['SubjectIdentifierAsRecorded'] != 'not applicable' 
                                and not pd.isna(x['SubjectIdentifierAsRecorded']) 
                                else value['missing'], axis=1)
            if key == 'UBERONOntologyIndex':
                df = df.merge(uberon_ontology_table[['Label', 'UBERONOntologyIndex', 'Is Multicellular','Is Organ', 'Is Fluid']], left_on='UBERONBodyPartName', right_on='Label', how='left')
                df.loc[df['Is Multicellular'] == True, 'SampleTypeSub1'] = 'tissue'
                df.loc[df['Is Fluid'] == True, 'SampleTypeSub1'] = 'biofluid'
                df.drop(columns=['Label', 'Is Multicellular', 'Is Organ', 'Is Fluid'], inplace=True)
            if key == 'DOIDOntologyIndex':
                df = df.merge(doid_ontology_table[['Label', 'DOIDOntologyIndex']], left_on='DOIDCommonName', right_on='Label', how='left')
                df.drop(columns=['Label'], inplace=True)
            if key == 'ENVOEnvironmentBiomeIndex':
                df = df.merge(envome_biome_ontology_table[['Label', 'ENVOEnvironmentBiomeIndex']], left_on='ENVOEnvironmentBiome', right_on='Label', how='left')
                df.drop(columns=['Label'], inplace=True)
            if key == 'ENVOEnvironmentMaterialIndex':
                df = df.merge(envome_material_ontology_table[['Label', 'ENVOEnvironmentMaterialIndex']], left_on='ENVOEnvironmentMaterial', right_on='Label', how='left')
                df.drop(columns=['Label'], inplace=True)
            if key == 'ENVOMediumScaleIndex':
                df = df.merge(envome_medium_scale_ontology_table[['Label', 'ENVOMediumScaleIndex']], left_on='ENVOMediumScale', right_on='Label', how='left')
                df.drop(columns=['Label'], inplace=True)
            if key == 'ENVOLocalScaleIndex':
                df = df.merge(envome_local_scale_ontology_table[['Label', 'ENVOLocalScaleIndex']], left_on='ENVOLocalScale', right_on='Label', how='left')
                df.drop(columns=['Label'], inplace=True)
            if key == 'ENVOBroadScaleIndex':
                df = df.merge(envome_broad_scale_ontology_table[['Label', 'ENVOBroadScaleIndex']], left_on='ENVOBroadScale', right_on='Label', how='left')
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
                           "SampleCollectionMethod", "ComorbidityListDOIDIndex", "DOIDCommonName", 
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
    'SampleType': {
        'standard': 'molecular_standards', 'culture_mamalian': 'culture_mammalian',
        'deceased animal': 'animal', 'blank': 'blank_analysis', 'pool': 'pool_QC', 'qc': 'pool_QC',
    },
    'IonizationSourceAndPolarity': {
        'heated electrospray ionization (positive)': 'electrospray ionization (positive)',
        'heated electrospray ionization (negative)': 'electrospray ionization (negative)',
        'heated electrospray ionization (alternating)': 'electrospray ionization (alternating)',
    },
    'ChromatographyAndPhase': {
        'reverse phase (polar c18)': 'reverse phase (C18)',
    },
}

# Prefixes of the first-generation ReDU template (e.g. 'ATTRIBUTE_Curated_SampleType',
# 'Subject_BiologicalSex', 'Analysis_MassSpectrometer') in front of current field names
_LEGACY_COLUMN_PREFIXES = ('ATTRIBUTE_Curated_', 'ATTRIBUTE_', 'Subject_', 'Analysis_', 'Study_', 'Curated_')

_MISSING_LIKE = {'', 'nan', 'na', 'n/a', 'none', 'missing value', 'missing_value', 'not specified', 'not collected',
                 'not applicable', 'not_applicable', 'ml import: not available', 'not available'}


def _is_missing_like(value):
    return str(value).strip().lower() in _MISSING_LIKE


def normalize_gnps_metadata(df, allowedTerm_dict, dataset_accession=None, uberon_synonym_table=None,
                            envo_tables=(), NCBIRankDivision_table=None):
    """Repairs common, unambiguous deviations in submitted GNPS/ReDU metadata files before
    validation, so the values are kept instead of being replaced by 'missing value'."""
    df = df.copy()

    # first-generation template column names -> current names
    for column in list(df.columns):
        if column in allowedTerm_dict:
            continue
        for prefix in _LEGACY_COLUMN_PREFIXES:
            name = column[len(prefix):]
            if column.startswith(prefix) and name in allowedTerm_dict and name not in df.columns:
                df = df.rename(columns={column: name})
                print(f"Renamed legacy column {column} -> {name}")
                break

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

    # 'lat | lon' with spaces, Excel float years
    if 'LatitudeandLongitude' in df.columns:
        df['LatitudeandLongitude'] = df['LatitudeandLongitude'].astype(str).str.replace(r'\s*\|\s*', '|', regex=True).str.strip()
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
    args = parser.parse_args()

    with open(args.AllowedTermJson_path, 'r') as json_file:
        allowed_terms = json.load(json_file)

    uberon_ontology_table = pd.read_csv(args.path_to_uberon_cl_po_csv, index_col=False)
    uberon_ontology_table = uberon_ontology_table.drop_duplicates(subset=['Label'])

    doid_ontology_table = pd.read_csv(args.path_to_doid_csv, index_col=False)
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
                                              NCBIRankDivision_table=NCBIRankDivision_table,
                                              ENVOEnvironmentBiomeIndex_table=ENVOEnvironmentBiomeIndex_table,
                                              ENVOEnvironmentMaterialIndex_table=ENVOEnvironmentMaterialIndex_table,
                                              attempt_adding_file_extensions=True,
                                              case_insensitive_terms=True)
            
            if len(df) > 0:
            
                file_name = os.path.join(args.output_metadata_folder, os.path.basename(file_path))
                df.to_csv(file_name, sep='\t', index=False)
