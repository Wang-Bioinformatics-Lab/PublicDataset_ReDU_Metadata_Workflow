#!/usr/bin/env python3
import argparse
import requests
import pandas as pd
from io import StringIO
from tqdm import tqdm
import re
from urllib.parse import unquote
import numpy as np
from REDU_conversion_functions import get_taxonomy_id_from_name__allowedTerms
from REDU_conversion_functions import map_instrument_to_allowed
from ms_instrument_matching import match_mass_spectrometer
import json
import os
import sys
import traceback
from datetime import datetime
from getAllNORMAN_file_paths import process_dataset_files, SESSION, HTTP_TIMEOUT
from read_and_validate_redu_from_github import complete_and_fill_REDU_table


def norm(s):
    """Lower-case and keep only letters/digits, used to match instrument names."""
    return re.sub(r'[^a-z0-9]', '', str(s).lower())

def create_usi(row):
    return f"mzspec:DSFP_{row['uuid']}:{row['file_names']}"

def get_metadata_sheet(url):
    print(f"Fetching file data from {url}")
    file_response = SESSION.get(url, timeout=HTTP_TIMEOUT)

    if file_response.status_code == 200:
        
        csv_data = StringIO(file_response.text)
        # read as str so ids (sample_id, instrument_id) join consistently across sheets
        try:
            df_files = pd.read_csv(csv_data, dtype=str)
        except pd.errors.EmptyDataError:
            print(f"Empty sheet at {url}")
            return pd.DataFrame()
        return df_files
    else:
        error_message = f"Failed to download files for {url}. Status code: {file_response.status_code}"
        print(error_message)
        return pd.DataFrame()  


def str_contains(series, pattern):
    """Case-insensitive regex search that never crashes on NaNs or numerics."""
    return (
        series.fillna('')          
              .astype(str)         
              .str.replace(' ', '', regex=True)
              .str.contains(pattern, case=False, regex=True)
    )

def collapse_conflicting_rows(df, key, missing_value='missing value'):
    """One row per key; columns whose values disagree within a key are set to missing_value.

    Used for instrument setups: NORMAN does not link a raw file to one setup, so when a
    sample has several setups we cannot tell which file used which polarity/column.
    """
    if not df[key].duplicated().any():
        return df

    def _collapse(series):
        values = series.dropna().unique()
        if len(values) == 0:
            return np.nan
        if len(values) == 1:
            return values[0]
        return missing_value

    n_before = len(df)
    df = df.groupby(key, as_index=False, sort=False).agg(_collapse)
    print(f"Collapsed {n_before} instrument setup rows into {len(df)} (one per {key})")
    return df


def clean(series):
    """Lower-case, drop spaces, tolerate NaNs/numbers → always returns a string Series."""
    return (
        series.fillna('')        # NaN → ''
              .astype(str)       # numeric → '123'
              .str.replace(r'\s+', '', regex=True)
              .str.lower()
    )

# NORMAN (Drupal) dates look like 'Mon, 07/01/2019 - 22:00' (month first)
NORMAN_DATE_FORMAT = '%a, %m/%d/%Y - %H:%M'


def parse_norman_date(value):
    """Parsed date, or None if not in NORMAN format or implausible (typos like year 0003)."""
    try:
        parsed = datetime.strptime(str(value).strip(), NORMAN_DATE_FORMAT)
    except ValueError:
        return None
    if parsed.year < 1900 or parsed > datetime.now():
        return None
    return parsed


def norman_collection_date(value):
    """Sampling date as an ISO date. The time is dropped (the stored hours are UTC
    shifts such as 22:00/03:00, not sampling times), and 01/01/YYYY, which submitters
    use as a placeholder for 'some time in YYYY', keeps only the year. Values in other
    formats pass through for the shared datetime validation."""
    if pd.isna(value):
        return value
    parsed = parse_norman_date(value)
    if parsed is None:
        # implausible NORMAN-format dates are dropped; other formats go to the validation
        return np.nan if re.match(r'^[A-Za-z]{3}, \d', str(value).strip()) else str(value)
    if parsed.month == 1 and parsed.day == 1:
        return parsed.strftime('%Y')
    return parsed.strftime('%Y-%m-%d')


def norman_year(value):
    parsed = parse_norman_date(value)
    if parsed is not None:
        return parsed.strftime('%Y')
    parsed = pd.to_datetime(value, errors='coerce')
    return parsed.strftime('%Y') if pd.notna(parsed) else np.nan


# Sheets use both 'Real Sample'/'Blank Sample' and the short 'Sample'/'Blank'/'QC'
SAMPLE_TYPE_NAMES = {'Sample': 'Real Sample', 'Blank': 'Blank Sample', 'QC': 'QC Sample'}

# Files CSV matrix -> SampleTypeSub1 (water types only; sediment, air, biota have no
# matching SampleTypeSub1 and are described by ENVOEnvironmentMaterial / SampleType)
MATRIX_TO_SAMPLETYPESUB1 = {'Surface water': 'water_surface', 'Wastewater': 'water_waste',
                            'Groundwater': 'water_ground', 'Ground water': 'water_ground'}

# Water body / habitat descriptions -> ENVO biome. 'Coastal region' is left out: it is
# used for marine and for terrestrial coastal samples alike.
HABITAT_TO_BIOME = {'marine environment': 'marine biome', 'marine water': 'marine biome',
                    'coastal water': 'marine biome', 'river': 'freshwater river biome',
                    'river water': 'freshwater river biome', 'lake': 'freshwater lake biome',
                    'lake water': 'freshwater lake biome', 'transitional water': 'estuarine biome',
                    'terrestrial environment': 'terrestrial biome'}

# Species group (the number is NORMAN's term id) -> SampleType
SPECIES_GROUP_TO_SAMPLETYPE = {'fish': 'animal', 'bird': 'animal', 'other': 'animal', 'molluscs': 'animal',
                               'mussel': 'animal', 'adult mussel': 'animal', 'juvenile mussel': 'animal',
                               'rat': 'animal', 'deer': 'animal', 'earthworm': 'animal', 'eggs': 'animal',
                               'crab': 'animal', 'barnacle': 'animal', 'dab': 'animal',
                               'vegetation': 'plant', 'macrophytes': 'plant', 'algae': 'algae'}
VERTEBRATE_GROUPS = {'fish', 'bird', 'rat', 'deer', 'dab'}
VERTEBRATE_DIVISIONS = {'Vertebrates', 'Mammals', 'Primates', 'Rodents'}
PLANT_GROUPS = {'vegetation', 'macrophytes', 'algae'}

COUNTRY_NAMES = {'UK': 'United Kingdom', 'England': 'United Kingdom', 'Scotland': 'United Kingdom',
                 'Wales': 'United Kingdom', 'Northern Ireland': 'United Kingdom',
                 'USA': 'United States of America', 'US': 'United States of America'}

# Files CSV matrix of the human biomonitoring sheets -> SampleTypeSub1 (SampleType 'animal',
# as for human samples in the other sources). Urine is set here because the ontology
# table does not flag it as a fluid.
HUMAN_MATRIX_TO_SAMPLETYPESUB1 = {'Blood': 'biofluid', 'Urine': 'biofluid', 'Placenta': 'tissue',
                                  'Human Biomonitoring': ''}

# Instrument names the vocabulary matching cannot resolve on its own (model field, or
# the instrument title when the model is empty) -> MassSpectrometer label
NORMAN_INSTRUMENT_NAMES = {'exploris 120': 'Orbitrap Exploris 120',
                           'orbitrap velos': 'LTQ Orbitrap Velos',
                           'maxis hd': 'Bruker Daltonics maXis series',
                           'uhplc-orbitrap-id-x-tribrid': 'Orbitrap ID-X'}

# GC columns with a 5%-phenyl methylpolysiloxane phase (DB-5 equivalents)
GC_DB5_COLUMN = r'(?<![a-z0-9])(?:db|hp|tg|rxi|rtx|zb|vf|bpx|inertcap|intercap)[- ]?5(?![0-9])'

# Sample sheet columns naming the water body of a surface-water sample
WATER_BODY_COLUMNS = ['Name of river, estuary, lake, reservoir or sea', 'Name of river, estuary, lake, reservoir, sea',
                      'Name of surface water (river/estuary/lake/reservoir/sea)']

# Whole words in a water-body name that identify its biome. German 'See' is a lake and
# English 'sea' the sea; Dutch/German 'meer' (lake or sea) and bays are left out.
WATER_BODY_WORDS = {'river': 'freshwater river biome', 'rivers': 'freshwater river biome', 'rio': 'freshwater river biome',
                    'río': 'freshwater river biome', 'rivière': 'freshwater river biome', 'riviere': 'freshwater river biome',
                    'fluss': 'freshwater river biome', 'fiume': 'freshwater river biome',
                    'lake': 'freshwater lake biome', 'lago': 'freshwater lake biome', 'lac': 'freshwater lake biome',
                    'see': 'freshwater lake biome',
                    'sea': 'marine biome', 'ocean': 'marine biome', 'mer': 'marine biome',
                    'estuary': 'estuarine biome', 'estuaire': 'estuarine biome'}
# names of catchments or artificial waters say nothing certain about the sampled water body
WATER_BODY_EXCLUDE = {'basin', 'catchment', 'watershed', 'district', 'canal', 'channel', 'kanal', 'gracht', 'reservoir'}

# 'Sediment Type' -> ENVO material (only exact ENVO labels; river sediment has none)
SEDIMENT_TYPE_TO_MATERIAL = {'marine sediment': 'marine sediment', 'sediment-lake': 'lake sediment',
                             'sedimenet-lake': 'lake sediment'}

# Country names that are also common names of rivers, lakes or regions elsewhere
AMBIGUOUS_COUNTRY_NAMES = {'Georgia', 'Jordan', 'Chad', 'Niger'}


def water_body_biome(name):
    """'Danube river' -> 'freshwater river biome', 'Black Sea' -> 'marine biome'; names
    without exactly one biome word (bare names, basins, canals) give ''."""
    if pd.isna(name):
        return ''
    words = set(re.findall(r'[^\W\d_]+', str(name).lower()))
    if words & WATER_BODY_EXCLUDE:
        return ''
    biomes = {WATER_BODY_WORDS[w] for w in words if w in WATER_BODY_WORDS}
    return biomes.pop() if len(biomes) == 1 else ''


def dataset_country(dataset_record, allowed_countries):
    """The one country named in the dataset keywords / spatial description, or None
    when none or several are named."""
    texts = list(dataset_record.get('keyword') or [])
    if isinstance(dataset_record.get('spatial'), str):
        texts.append(dataset_record['spatial'])
    names = {c: c for c in allowed_countries if c not in AMBIGUOUS_COUNTRY_NAMES}
    names.update({k: v for k, v in COUNTRY_NAMES.items() if k != 'US'})
    found = set()
    for text in texts:
        for name, country in names.items():
            if re.search(r'(?<![A-Za-z])' + re.escape(name) + r'(?![A-Za-z])', str(text)):
                found.add(country)
    # samples from Antarctica (e.g. a national research station) are not in that country
    if any(re.search(r'antarctic', str(text), re.IGNORECASE) for text in texts):
        return None
    return found.pop() if len(found) == 1 else None


def norman_chromatography(column_model, instrument_type=None):
    """instrument-setups column_model (+ instrument_type, e.g. 'GC-EI-QFT') ->
    ChromatographyAndPhase, '' when unknown."""
    c = '' if pd.isna(column_model) else str(column_model).lower()
    is_hilic = 'hilic' in c
    if is_hilic and 'c18' in c:
        # serially coupled C18 and HILIC columns ('Poroshell 120 EC C18 + ZIC-HILIC')
        return '2D:HILIC and reverse phase (C18)'
    if is_hilic and re.search('c8|c30|reverse', c):
        return 'conflict (HILIC cannot coexist with C18/C8/C30)'
    if is_hilic:
        return 'normal phase (HILIC)'
    for token, term in [('c18', 'reverse phase (C18)'), ('c8', 'reverse phase (C8)'), ('c30', 'reverse phase (C30)')]:
        if token in c:
            return term
    # Waters T3 (trifunctional C18) columns: HSS T3, Atlantis T3
    if re.search(r'(?:hss|atlantis)\s*t3', c):
        return 'reverse phase (C18)'
    # reverse phase columns without a C18/C8 label (Accucore RP-MS, ZORBAX SB-Aq)
    if 'reverse phase' in c or re.search(r'(?<![a-z])rp(?![a-z])', c) or 'sb-aq' in c:
        return 'reverse phase (NOS)'
    is_gc = bool(re.search('gc|gas', c)) or str(instrument_type).strip().lower().startswith('gc')
    if is_gc:
        return 'gas chromatography (DB-5)' if re.search(GC_DB5_COLUMN, c) else 'gas chromatography (NOS)'
    return ''


_ENVO_SYNONYMS = {}


def envo_material_label(value, material_table, allowedTerm_dict):
    """An allowed ENVO material label for a value: the value itself, or the one label
    whose synonym it is. Other values are returned unchanged."""
    if pd.isna(value) or value == '':
        return value
    allowed = allowedTerm_dict['ENVOEnvironmentMaterial']['allowed_values']
    if not _ENVO_SYNONYMS:
        labels = {}
        for label, synonym in zip(material_table['Label'].astype(str), material_table['Synonym'].astype(str)):
            labels.setdefault(synonym.strip().lower(), set()).add(label)
        _ENVO_SYNONYMS.update(labels)
        _ENVO_SYNONYMS['__allowed__'] = set(allowed)
    if value in _ENVO_SYNONYMS['__allowed__']:
        return value
    labels = {l for l in _ENVO_SYNONYMS.get(str(value).strip().lower(), set()) if l in _ENVO_SYNONYMS['__allowed__']}
    return labels.pop() if len(labels) == 1 else value


def norm_plus(s):
    """norm() that keeps '+', so 'TripleTOF 5600' and 'TripleTOF 5600+' stay apart."""
    return re.sub(r'[^a-z0-9+]', '', str(s).lower())


def norman_mass_spectrometer(model, title, allowed_ms):
    """MassSpectrometer term for an instruments.csv row: the model, or the instrument title
    when the model is empty. Curated names first, then the shared vocabulary matcher,
    then the former NORMAN normalisation (leading 'Orbitrap '/'maXis ' dropped,
    '6550 iFunnel Q-TOF' -> 'iFunnel Q-TOF 6550') with an exact or token match."""
    name = model if isinstance(model, str) and model.strip() else title
    if not isinstance(name, str) or not name.strip():
        return np.nan
    name = ' '.join(name.split())
    labels = {v.split('|')[0].lower(): v for v in allowed_ms}
    curated = NORMAN_INSTRUMENT_NAMES.get(name.lower())
    if curated and curated.lower() in labels:
        return labels[curated.lower()]
    hit = match_mass_spectrometer(name, allowed_ms)
    if hit:
        return hit
    stripped = re.sub(r'^(orbitrap|maxis)\s+', '', name, flags=re.IGNORECASE)
    stripped = re.sub(r'^\s*(\d+)\s+(.*)$', r'\2 \1', stripped).strip()
    exact = {}
    for v in allowed_ms:
        exact.setdefault(norm_plus(v.split('|')[0]), set()).add(v)
    hits = exact.get(norm_plus(stripped), set())
    if len(hits) == 1:
        return next(iter(hits))
    return map_instrument_to_allowed(stripped, allowed_ms) or np.nan


def species_group(value):
    """'Fish (22486)' -> 'fish'"""
    if pd.isna(value):
        return ''
    return re.sub(r'\s*\(\d+\)\s*$', '', str(value)).strip().lower()


def tissue_for_group(tissue, group, division=None):
    """NORMAN Tissue -> a body part label ReDU resolves; None keeps the original value.
    'Muscle' alone is ambiguous in the ontology synonyms, so it is resolved by organism:
    the NCBI division of the species if known, else the species group ('Other' mixes
    mammals and invertebrates, so it falls back to the generic muscle tissue)."""
    if pd.isna(tissue):
        return None
    t = str(tissue).strip().lower()
    if t in ('muscle', 'muscel'):
        if group in PLANT_GROUPS:
            return ''
        vertebrate = division in VERTEBRATE_DIVISIONS if division else group in VERTEBRATE_GROUPS
        return 'skeletal muscle tissue' if vertebrate else 'muscle tissue'
    if t in ('whole body', 'soft body'):
        return 'multicellular organism'
    if t == 'egg':
        return 'shelled egg'
    if t == 'leaves':
        return 'leaf'
    return None


def load_species_names():
    path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'translation_sheets', 'NORMAN_species_names.csv')
    table = pd.read_csv(path, dtype=str)
    return {k.strip().lower(): v.strip() for k, v in zip(table['name_in_source'], table['scientific_name'])}


def scientific_name_candidates(value, species_names):
    """Names to try for a 'Biota species name (in Latin)' value, most specific first:
    curated common name / typo fix, 'Common, Latin', 'Common (Latin)', the raw value."""
    s = str(value).strip()
    candidates = []
    if s.lower() in species_names:
        candidates.append(species_names[s.lower()])
    m = re.fullmatch(r'(.+?)\s*\(([^)]+)\)', s)
    if m:
        candidates += [species_names.get(m.group(2).strip().lower(), m.group(2).strip()),
                       species_names.get(m.group(1).strip().lower(), '')]
    if ',' in s:
        latin = s.rsplit(',', 1)[1].strip()
        candidates.append(species_names.get(latin.lower(), latin))
    candidates.append(s)
    return [c for c in candidates if c]


def main(output_filename, study_id, **kwargs):


    # Load allowed terms and ontology tables
    allowedTerm_dict = kwargs['allowedTerm_dict']
    ontology_table = kwargs['ontology_table']
    ENVOEnvironmentBiomeIndex_table = kwargs['ENVOEnvironmentBiomeIndex_table']
    ENVOEnvironmentMaterialIndex_table = kwargs['ENVOEnvironmentMaterialIndex_table']
    NCBIRankDivision_table = kwargs['NCBIRankDivision_table']
    species_names = load_species_names()
    division_by_taxid = dict(zip(NCBIRankDivision_table['TaxonID'].astype(str), NCBIRankDivision_table['NCBIDivision']))


    # Fetch the list of datasets
    datasets_url = "https://dsfp.norman-data.eu/api/1/metastore/schemas/dataset/all"
    print(f"Fetching datasets from {datasets_url}", flush = True)
    response = SESSION.get(datasets_url, timeout=HTTP_TIMEOUT)
    response.raise_for_status() 
    datasets = response.json()
    print(f"Fetched {len(datasets)} datasets", flush = True)

    # Filter datasets based on study_id
    if study_id != "ALL":
        datasets = [ds for ds in datasets if ds['uuid'] == study_id]
        print(f"Processing dataset with UUID: {study_id}", flush = True)
    else:
        print("Processing all datasets", flush = True)

    dfs = [] 
    errors = [] 
    skipped = []


    for ds in tqdm(datasets, desc="Processing datasets", unit="dataset"):
        internal_id = ds.get('internal_id')
        try:
            uuid = ds['uuid']
            internal_id = ds['internal_id']
            title = ds['title']
            print(f"Processing dataset: {title} (UUID: {uuid}, Internal ID: {internal_id})", flush = True)

            # Build URL to get the CSV with file info
            metadata_collection_url = f"https://dsfp.norman-data.eu/api/1/metastore/schemas/dataset/items/{uuid}"
            print(f"Fetching file data from {metadata_collection_url}", flush = True)
            metadata_collection_response = SESSION.get(metadata_collection_url, timeout=HTTP_TIMEOUT)

            if metadata_collection_response.status_code == 200:
                
                metadata_collection_dict = metadata_collection_response.json()
                
                print(f"Fetched metadata collection dict for dataset {uuid}", flush = True)


                # Records without sheets (no 'distribution') hold no files to index
                if not metadata_collection_dict.get('distribution'):
                    message = f"Skipping dataset {internal_id}: no sheets listed in the dataset record"
                    print(message)
                    skipped.append(message)
                    continue
                metadata_collection_links = metadata_collection_dict['distribution']

                nomran_sheets_for_single_dataset_sheets = {'Sample info': [],
                                                           'File info': [],
                                                           'Spectral raw files': [],
                                                           'Instrument info': [],
                                                           'Instrument setup': []}


                # Get all Samples CSV data
                for metadata_sheet in metadata_collection_links:

                    download_url = metadata_sheet['downloadURL']
                    df_metadata_sheet = pd.DataFrame()

                    if metadata_sheet['title'].endswith(' - Samples CSV'):
                        df_metadata_sheet = get_metadata_sheet(download_url)

                        required_columns = ['ID', 'Sample type', 'Instrument setup used']

                        if df_metadata_sheet.empty or len(df_metadata_sheet) == 0 or not all(col in df_metadata_sheet.columns for col in required_columns): 
                            print(f"Skipping {metadata_sheet['title']} as it is empty or can not be downloaded or required columns not present.")
                            print(f"URL: {download_url}")
                            continue
                                            

                        df_metadata_sheet = df_metadata_sheet.rename(columns={'ID': 'sample_id', 'Sample type': 'sample_type'})
                        df_metadata_sheet['sample_type'] = df_metadata_sheet['sample_type'].replace(SAMPLE_TYPE_NAMES)
                        is_real = df_metadata_sheet['sample_type'] == 'Real Sample'

                        # 'Mon, 07/01/2019 - 22:00' -> '2019-07-01' (date only, 01/01 placeholders
                        # keep only the year); the shared validation then checks the ISO value
                        if 'Sampling date' in df_metadata_sheet.columns:
                            df_metadata_sheet['SampleCollectionDateandTime'] = df_metadata_sheet['Sampling date'].map(norman_collection_date)

                        # in YearOfAnalysis reformat 'Sat, 11/27/1999 - 22:00' to 1999
                        if 'Analysis date' in df_metadata_sheet.columns:
                            df_metadata_sheet['YearOfAnalysis'] = df_metadata_sheet['Analysis date'].map(norman_year)

                        # Make LatitudeandLongitude column in format latitude and longitude separated by a |
                        if 'Latitude' in df_metadata_sheet.columns and 'Longitude' in df_metadata_sheet.columns:
                            df_metadata_sheet['LatitudeandLongitude'] = df_metadata_sheet['Latitude'].astype(str) + '|' + df_metadata_sheet['Longitude'].astype(str)

                        # Make DepthorAltitudeMeters first to numeric to make sign negative and then to string
                        if 'Depth (m)' in df_metadata_sheet.columns:
                            df_metadata_sheet['DepthorAltitudeMeters'] = pd.to_numeric(df_metadata_sheet['Depth (m)'], errors='coerce').apply(lambda x: str(-x) if pd.notna(x) else "").astype(str)

                        # Waste water overwrite if Type of wastewater sample column exists
                        if 'Type of wastewater sample' in df_metadata_sheet.columns:
                            df_metadata_sheet['ENVOEnvironmentMaterial_overwrite'] = np.where(
                                df_metadata_sheet['sample_type'] == 'Real Sample',  
                                df_metadata_sheet['Type of wastewater sample']
                                .map({'Effluent wastewater': 'treated wastewater', 'Influent wastewater':  'waste water',
                                    'Primary sludge': 'raw primary sludge', 'Sewage sludge': 'sewage', 'Bio sludge': 'biological waste material'}) 
                                .fillna(''),
                                ''
                                )

                        # Surface water overwrite if Type of Surface water column exists
                        if 'Type of Surface water' in df_metadata_sheet.columns:
                            df_metadata_sheet['ENVOEnvironmentMaterial_overwrite'] = np.where(
                                df_metadata_sheet['sample_type'] == 'Real Sample',  
                                df_metadata_sheet['Type of Surface water']
                                .map({'River water': 'river water', 'Marine water':  'sea water', 'Lake water': 'lake water',
                                    'Transitional water': 'brackish water', 'Coastal water': 'coastal sea water'})
                                .fillna(''),
                                ''
                                )

                        # Soil overwrite if Soil Types column exists
                        if 'Soil Types' in df_metadata_sheet.columns:
                            df_metadata_sheet['ENVOEnvironmentMaterial_overwrite'] = np.where(
                                df_metadata_sheet['sample_type'] == 'Real Sample',  
                                df_metadata_sheet['Soil Types']
                                .map({'Acrisols': 'acrisol', 'Alisols': 'alisol', 'Andosols': 'andosol', 'Anthrosols': 'anthrosol',
                                    'Arenosols': 'arenosol', 'Calcisols': 'calcisol', 'Cambisols': 'cambisol', 'Chernozems': 'chernozem',
                                    'Cryosols': 'cryosol', 'Durisols': 'durisol', 'Ferralsols': 'ferralsol', 'Fluvisols': 'fluvisol',
                                    'Gleysols': 'gleysol', 'Gypsisols': 'gypsisol', 'Histosols': 'histosol', 'Kastanozems': 'kastanozem',
                                    'Leptosols': 'leptosol', 'Lixisols': 'lixisol', 'Luvisols': 'luvisol', 'Phaeozems': 'phaeozem',
                                    'Planosols': 'planosol', 'Plinthosols': 'plinthosol', 'Podzols': 'podzol', 'Regosols': 'regosol',
                                    'Solonchaks': 'solonchak', 'Solonetz': 'solonetz', 'Umbrisols': 'umbrisol', 'Vertisols': 'vertisol'})
                                .fillna(''),
                                ''
                                )                                               

                        # Biome from the water body / habitat of the sample
                        for habitat_column in ['Type of Surface water', 'Origin of biota', 'Origin of SPM']:
                            if habitat_column in df_metadata_sheet.columns:
                                biome = df_metadata_sheet[habitat_column].astype('string').str.strip().str.lower().map(HABITAT_TO_BIOME)
                                df_metadata_sheet['ENVOEnvironmentBiome_overwrite'] = np.where(is_real & biome.notna(), biome, '')

                        # Biome of surface water from the water body name ('Danube river', 'Black Sea');
                        # applied after the merge to Surface water files only
                        water_body = None
                        for water_body_column in WATER_BODY_COLUMNS:
                            if water_body_column in df_metadata_sheet.columns:
                                water_body = (df_metadata_sheet[water_body_column] if water_body is None
                                              else water_body.fillna(df_metadata_sheet[water_body_column]))
                        if water_body is not None:
                            df_metadata_sheet['water_body_biome'] = np.where(is_real, water_body.map(water_body_biome), '')

                        # Sediment overwrite from 'Sediment Type' ('Sediment-Lake', 'Marine sediment'),
                        # else 'Origin of Sediment'
                        if 'Sediment Type' in df_metadata_sheet.columns or 'Origin of Sediment' in df_metadata_sheet.columns:
                            sediment = pd.Series(np.nan, index=df_metadata_sheet.index, dtype=object)
                            if 'Sediment Type' in df_metadata_sheet.columns:
                                sediment = df_metadata_sheet['Sediment Type'].astype('string').str.strip().str.lower().map(SEDIMENT_TYPE_TO_MATERIAL)
                            if 'Origin of Sediment' in df_metadata_sheet.columns:
                                sediment = sediment.fillna(df_metadata_sheet['Origin of Sediment']
                                    .map({'Lake': 'lake sediment', 'Marine environment':  'marine sediment', 'Coastal region': 'marine sediment',
                                        'River': 'sediment permeated by freshwater'}))
                            df_metadata_sheet['ENVOEnvironmentMaterial_overwrite'] = np.where(is_real, sediment.fillna(''), '')
                        
                        # Biota: SampleType from the species group (fish, birds, mussels -> animal, ...)
                        groups = (df_metadata_sheet['Species group'].map(species_group) if 'Species group' in df_metadata_sheet.columns
                                  else pd.Series('', index=df_metadata_sheet.index))
                        if 'Species group' in df_metadata_sheet.columns:
                            df_metadata_sheet['SampleType_overwrite'] = np.where(
                                is_real, groups.map(SPECIES_GROUP_TO_SAMPLETYPE).fillna(''), '')

                        # NCBI Taxonomy overwrite if Biota species name (in Latin) column exists.
                        # The column also holds 'Common, Latin', 'Common (Latin)', common names
                        # only and typos; candidates are tried in order against the allowed terms.
                        if 'Biota species name (in Latin)' in df_metadata_sheet.columns:
                            def _taxonomy(org):
                                if pd.isna(org):
                                    return ''
                                for candidate in scientific_name_candidates(org, species_names):
                                    hit = get_taxonomy_id_from_name__allowedTerms(candidate, allowedTerm_dict = allowedTerm_dict)
                                    if hit:
                                        return hit
                                return ''

                            df_metadata_sheet['NCBITaxonomy'] = np.where(
                                is_real,
                                df_metadata_sheet['Biota species name (in Latin)']
                                .map({org: _taxonomy(org) for org in df_metadata_sheet['Biota species name (in Latin)'].unique()})
                                .fillna(''),
                                ''
                                )


                        # if the column names include 'Occupation' and	'Family annual income bracket' assign NCBITaxonomy to "9606|Homo sapiens" and overwrite the column
                        if 'Occupation' in df_metadata_sheet.columns and 'Family annual income bracket' in df_metadata_sheet.columns:
                            df_metadata_sheet['NCBITaxonomy'] = np.where(
                                df_metadata_sheet['sample_type'] == 'Real Sample',          
                                '9606|Homo sapiens',       
                                ''                                                            
                            )



                        # Blood subcategory overwrite if Blood subcategory column exists
                        if 'Blood subcategory' in df_metadata_sheet.columns:
                            df_metadata_sheet['Tissue'] = df_metadata_sheet['Blood subcategory'].str.lower().replace({'serum': 'blood serum', 'plasma': 'blood plasma'})
                            df_metadata_sheet['Tissue'] = df_metadata_sheet['Tissue'].where(df_metadata_sheet['Tissue'].isin(['blood serum', 'blood plasma']), df_metadata_sheet['Blood subcategory'].str.lower())


                        

                        # if the tissue column does not exist split metadata_sheet['title']  by " - " and take the first part. then assign that value to the tissue column for all real samples
                        if 'Tissue' not in df_metadata_sheet.columns:
                            tissue = metadata_sheet['title'].split(' - ')[0]
                            df_metadata_sheet['Tissue'] = np.where(
                                df_metadata_sheet['sample_type'] == 'Real Sample',  
                                tissue,
                                ''
                                )


                        # Resolve NORMAN tissue names ('Muscle', 'Muscel', 'Soft body', ...) by organism
                        if 'Tissue' in df_metadata_sheet.columns:
                            divisions = (df_metadata_sheet['NCBITaxonomy'].astype(str).str.split('|').str[0].map(division_by_taxid)
                                         if 'NCBITaxonomy' in df_metadata_sheet.columns else pd.Series(np.nan, index=df_metadata_sheet.index))
                            resolved = [tissue_for_group(t, g, d if isinstance(d, str) else None)
                                        for t, g, d in zip(df_metadata_sheet['Tissue'], groups, divisions)]
                            df_metadata_sheet['Tissue'] = [t if r is None else r for t, r in zip(df_metadata_sheet['Tissue'], resolved)]

                        # If Tissue exists renamethe column to UBERONBodyPartName
                        if 'Tissue' in df_metadata_sheet.columns:
                            synonyms_in_ontology_table = set(ontology_table['Synonym'].str.lower())  
                            allowed_bodyparts = set(allowedTerm_dict['UBERONBodyPartName']['allowed_values'])
                            df_metadata_sheet['UBERONBodyPartName'] = ''

                            
                            for index, row in df_metadata_sheet.iterrows():
                                tissue = row['Tissue']  
                                tissue = '' if pd.isna(tissue) else str(tissue).lower()
                                
                                if tissue in allowed_bodyparts:
                                    matching_label = next((term for term in allowed_bodyparts if term.lower() == tissue), tissue)
                                    df_metadata_sheet.at[index, 'UBERONBodyPartName'] = matching_label
                                    continue 

                                # If the tissue is a synonym, get the matching rows from the ontology
                                if tissue in synonyms_in_ontology_table:
                                    matching_rows = ontology_table[ontology_table['Synonym'].str.lower() == tissue]

                                    # Filter matching rows to ensure a unique match
                                    matching_rows = matching_rows.drop_duplicates(subset=['Label', 'UBERONOntologyIndex'])

                                    # If there's exactly one match, assign it
                                    if len(matching_rows) == 1:
                                        label = matching_rows['Label'].iloc[0]  
                                        df_metadata_sheet.at[index, 'UBERONBodyPartName'] = label

                            # Now merge the harmonized terms with df_metadata_sheet
                            ontology_table_unique = ontology_table.drop_duplicates(subset=['Label', 'UBERONOntologyIndex']).drop(columns=['Synonym'])
                            df_metadata_sheet = pd.merge(df_metadata_sheet, ontology_table_unique, left_on='UBERONBodyPartName', right_on='Label', how='left')
                            df_metadata_sheet.drop('Label', axis=1, inplace=True)

                            # Assign sample types based on additional conditions
                            df_metadata_sheet.loc[df_metadata_sheet['Is Fluid'] == True, 'SampleTypeSub1'] = 'biofluid'
                            df_metadata_sheet.loc[df_metadata_sheet['Is Multicellular'] == True, 'SampleTypeSub1'] = 'tissue'


                        # Country from 'Monitored country' (or 'field_country'), with UK/USA spelled
                        # as in the vocabulary
                        if 'Monitored country' in df_metadata_sheet.columns or 'field_country' in df_metadata_sheet.columns:
                            country = df_metadata_sheet.get('Monitored country', pd.Series(np.nan, index=df_metadata_sheet.index))
                            if 'field_country' in df_metadata_sheet.columns:
                                country = country.fillna(df_metadata_sheet['field_country'])
                            df_metadata_sheet = df_metadata_sheet.drop(columns=[c for c in ['Monitored country', 'field_country'] if c in df_metadata_sheet.columns])
                            df_metadata_sheet['Country'] = country.astype('string').str.strip().replace(COUNTRY_NAMES)


                        nomran_sheets_for_single_dataset_sheets['Sample info'].append(df_metadata_sheet)


                    elif download_url.endswith('files.csv'):
                        fileinfo_redu_sheet = get_metadata_sheet(download_url)

                        if fileinfo_redu_sheet.empty or len(fileinfo_redu_sheet) == 0: 
                            print(f"Skipping {metadata_sheet['title']} as it is empty or can not be downloaded.")
                            continue               

                        required_columns = ['sample_id', 'sample_type', 'matrix']

                        if not all(col in fileinfo_redu_sheet.columns for col in required_columns):
                            print(f"Skipping {metadata_sheet['title']} as required columns are not present.")
                            print(f"URL: {download_url}")
                            continue
                        
                        fileinfo_redu_sheet['sample_type'] = fileinfo_redu_sheet['sample_type'].replace(SAMPLE_TYPE_NAMES)
                        spectral_raw_files_sheet  = fileinfo_redu_sheet.copy()

                        spectral_raw_files_sheet = process_dataset_files(spectral_raw_files_sheet, internal_id=internal_id, uuid=uuid, filter_extensions=True)
                        spectral_raw_files_sheet = spectral_raw_files_sheet[['sample_id', 'file_names', 'USI']]


                        matrix_type = (fileinfo_redu_sheet['matrix_type'] if 'matrix_type' in fileinfo_redu_sheet.columns
                                       else pd.Series(np.nan, index=fileinfo_redu_sheet.index))
                        fileinfo_redu_sheet = fileinfo_redu_sheet[required_columns].copy()
                        is_real_file = fileinfo_redu_sheet['sample_type'] == 'Real Sample'

                        fileinfo_redu_sheet['SampleType'] = fileinfo_redu_sheet['sample_type'].replace({'Blank Sample': 'blank_extraction'})
                        fileinfo_redu_sheet['SampleTypeSub1'] = np.where(
                            is_real_file,
                            fileinfo_redu_sheet['matrix'].map(MATRIX_TO_SAMPLETYPESUB1).fillna('Real Sample'),
                            fileinfo_redu_sheet['sample_type'].replace({'Blank Sample': 'blank_extraction'}))

                        # Human biomonitoring (blood, urine, placenta): human samples are 'animal'
                        is_human = is_real_file & fileinfo_redu_sheet['matrix'].isin(HUMAN_MATRIX_TO_SAMPLETYPESUB1.keys())
                        fileinfo_redu_sheet.loc[is_human, 'SampleType'] = 'animal'
                        human_sub1 = fileinfo_redu_sheet['matrix'].map(HUMAN_MATRIX_TO_SAMPLETYPESUB1).fillna('')
                        fileinfo_redu_sheet.loc[is_human & (human_sub1 != ''), 'SampleTypeSub1'] = human_sub1

                        # Food samples are 'food', not environmental material
                        fileinfo_redu_sheet.loc[is_real_file & (fileinfo_redu_sheet['matrix'] == 'Foodomics'), 'SampleType'] = 'food'

                        # Biome from the matrix type ('River', 'Marine environment', ...)
                        biome = matrix_type.astype('string').str.strip().str.lower().map(HABITAT_TO_BIOME)
                        fileinfo_redu_sheet['ENVOEnvironmentBiome'] = np.where(is_real_file & biome.notna(), biome, '')

                        fileinfo_redu_sheet['ENVOEnvironmentMaterial'] = np.where(
                            fileinfo_redu_sheet['sample_type'] == 'Real Sample',  
                            fileinfo_redu_sheet['matrix']
                            .map({'Wastewater': 'waste water', 'Groundwater': 'groundwater', 'Ground water': 'groundwater', 'Surface water': 'surface water', 'Sludge': 'sludge',
                                  'Soil': 'soil', 'Sediment': 'sediment', 'Indoor Air': 'air', 'Air-Indoor': 'air', 'SPM': 'waterborne particulate matter',
                                  'Ambient Air': 'air', 'Workplace Air': 'air', 'Emission Air': 'air'})
                            .fillna(''),
                            ''
                            )                                                           

                        fileinfo_redu_sheet['NCBITaxonomy_overwrite'] = np.where(
                            fileinfo_redu_sheet['sample_type'] == 'Real Sample',          
                            fileinfo_redu_sheet['matrix'].map({'Human Biomonitoring': '9606|Homo sapiens'}).fillna(''),       
                            ''                                                            
                        )

                    
                        nomran_sheets_for_single_dataset_sheets['Spectral raw files'].append(spectral_raw_files_sheet)
                        nomran_sheets_for_single_dataset_sheets['File info'].append(fileinfo_redu_sheet)

                    elif download_url.endswith('instruments.csv'):

                        instrumentinfo_sheet = get_metadata_sheet(download_url)

                        required_columns = ['instrument_id', 'instrument_model']

                        if instrumentinfo_sheet.empty or len(instrumentinfo_sheet) == 0 or not all(col in instrumentinfo_sheet.columns for col in required_columns):
                            print(f"Skipping {metadata_sheet['title']} as it is empty or can not be downloaded or required columns not present.")
                            print(f"URL: {download_url}")
                            continue

                        # Model name (or the instrument title when the model is empty) -> vocabulary
                        _allowed_ms = allowedTerm_dict["MassSpectrometer"]["allowed_values"]
                        _titles = (instrumentinfo_sheet['instrument_title'] if 'instrument_title' in instrumentinfo_sheet.columns
                                   else pd.Series(np.nan, index=instrumentinfo_sheet.index))
                        instrumentinfo_sheet["MassSpectrometer"] = [
                            norman_mass_spectrometer(model, title, _allowed_ms)
                            for model, title in zip(instrumentinfo_sheet['instrument_model'], _titles)]


                        nomran_sheets_for_single_dataset_sheets['Instrument info'].append(instrumentinfo_sheet)

                    elif download_url.endswith('instrument-setups.csv'):

                        instrumentsetup_sheet = get_metadata_sheet(download_url)

                        some_columns = ['sample_id', 'setup_id', 'instrument', 'ionization_type', 'ionization_esi_apci_appi', 'column_model']

                        required_columns = ['sample_id']

                        if instrumentsetup_sheet.empty or len(instrumentsetup_sheet) == 0 or not all(col in instrumentsetup_sheet.columns for col in required_columns):
                            print(f"Skipping {metadata_sheet['title']} as it is empty or can not be downloaded or required columns not present.")
                            print(f"URL: {download_url}")
                            continue

                        
                        if 'instrument' in instrumentsetup_sheet.columns:
                            instrumentsetup_sheet = instrumentsetup_sheet.rename(columns={'instrument': 'instrument_id'})

                        
                        if 'ionization_type' in instrumentsetup_sheet.columns and 'ionization_esi_apci_appi' in instrumentsetup_sheet.columns:

                            # first if "esi" or "electrospray" in ionization_esi_apci_appi (ignoring case and spaces) put "electrospray ionization"
                            # and if "apci" or "chemical ionization" or "chemical ionization" in ionization_esi_apci_appi (ignoring case and spaces) put "atmospheric pressure chemical ionization"
                            # and if "ei" or "electron impact" in ionization_esi_apci_appi (ignoring case and spaces) put "electron ionization"
                            # and if "pi" or "photoionization" in ionization_esi_apci_appi (ignoring case and spaces) put "atmospheric pressure photoionization"
                            s1 = str_contains(instrumentsetup_sheet['ionization_type'], r'esi|electrospray')
                            s2 = str_contains(instrumentsetup_sheet['ionization_esi_apci_appi'], r'esi|electrospray')

                            s3 = str_contains(instrumentsetup_sheet['ionization_type'], r'apci|atmosphericpressurechemicalionization')
                            s4 = str_contains(instrumentsetup_sheet['ionization_esi_apci_appi'], r'apci|atmosphericpressurechemicalionization')

                            s5 = str_contains(instrumentsetup_sheet['ionization_type'], r'ei|electronimpact')
                            s6 = str_contains(instrumentsetup_sheet['ionization_esi_apci_appi'], r'ei|electronimpact')

                            s7 = str_contains(instrumentsetup_sheet['ionization_type'], r'pi|photoionization')
                            s8 = str_contains(instrumentsetup_sheet['ionization_esi_apci_appi'], r'pi|photoionization')

                            instrumentsetup_sheet['IonizationSourceAndPolarity'] = np.select(
                                [
                                    s1 | s2,
                                    s3 | s4,
                                    s5 | s6,
                                    s7 | s8,
                                ],
                                [
                                    'electrospray ionization',
                                    'atmospheric pressure chemical ionization',
                                    'electron ionization',
                                    'atmospheric pressure photoionization',
                                ],
                                default=''
                            )


                            # if ionization_type is "positive" or "negative" append that value to IonizationSourceAndPolarity as " (positive)" or " (negative)" except for electron ionization
                            ion_type  = clean(instrumentsetup_sheet['ionization_type'])
                            ion_appi  = clean(instrumentsetup_sheet['ionization_esi_apci_appi'])
                            sourcepol = clean(instrumentsetup_sheet['IonizationSourceAndPolarity'])

                            # polarity masks (excluding EI)
                            pos_mask = ion_type.str.contains('positive') & ~sourcepol.str.contains('electronionization')
                            neg_mask = ion_type.str.contains('negative') & ~sourcepol.str.contains('electronionization')

                            # start from a string-safe column, then append in place
                            instrumentsetup_sheet['IonizationSourceAndPolarity'] = (
                                instrumentsetup_sheet['IonizationSourceAndPolarity'].fillna('').astype(str)
                            )

                            instrumentsetup_sheet.loc[pos_mask, 'IonizationSourceAndPolarity'] += ' (positive)'
                            instrumentsetup_sheet.loc[neg_mask, 'IonizationSourceAndPolarity'] += ' (negative)'

                        # Column model (+ instrument type for GC) -> ChromatographyAndPhase
                        if 'column_model' in instrumentsetup_sheet.columns:
                            _types = (instrumentsetup_sheet['instrument_type'] if 'instrument_type' in instrumentsetup_sheet.columns
                                      else pd.Series(np.nan, index=instrumentsetup_sheet.index))
                            instrumentsetup_sheet['ChromatographyAndPhase'] = [
                                norman_chromatography(column, instrument_type)
                                for column, instrument_type in zip(instrumentsetup_sheet['column_model'], _types)]

                        
                        nomran_sheets_for_single_dataset_sheets['Instrument setup'].append(instrumentsetup_sheet)

                        

                
                # Without file info and raw files there is nothing to index for this dataset
                missing_sheets = [name for name in ['File info', 'Spectral raw files']
                                  if len(nomran_sheets_for_single_dataset_sheets[name]) == 0]
                if missing_sheets:
                    message = f"Skipping dataset {internal_id}: no usable {', '.join(missing_sheets)}"
                    print(message)
                    skipped.append(message)
                    continue

                # Combine individual sheet types (instrument sheets are optional)
                combined_file_info = pd.concat(nomran_sheets_for_single_dataset_sheets['File info'], ignore_index=True)
                if nomran_sheets_for_single_dataset_sheets['Sample info']:
                    combined_sample_info = pd.concat(nomran_sheets_for_single_dataset_sheets['Sample info'], ignore_index=True)
                else:
                    # No usable Samples CSV: the Files CSV still gives sample type and matrix
                    print(f"Dataset {internal_id}: no usable Samples CSV, using the Files CSV only")
                    combined_sample_info = pd.DataFrame({'sample_id': combined_file_info['sample_id'].dropna().unique()})
                combined_instrument_info = (pd.concat(nomran_sheets_for_single_dataset_sheets['Instrument info'], ignore_index=True)
                                            if nomran_sheets_for_single_dataset_sheets['Instrument info'] else None)
                combined_instrument_setup_info = (pd.concat(nomran_sheets_for_single_dataset_sheets['Instrument setup'], ignore_index=True)
                                                  if nomran_sheets_for_single_dataset_sheets['Instrument setup'] else None)


                # Combine sample info and file info
                ########
                columns_to_drop = [col for col in combined_sample_info.columns if col in combined_file_info.columns and col != 'sample_id']
                combined_sample_info = combined_sample_info.drop(columns=columns_to_drop)
                # Keep every file of the Files CSV, also when its sample is missing from the
                # Samples CSVs (sheet not downloadable, blanks/QCs listed only in the Files CSV)
                merged_df = pd.merge(combined_file_info, combined_sample_info, on='sample_id', how='left')
                print(f"Dataset {internal_id}: {len(combined_sample_info)} samples, {len(combined_file_info)} file info rows, {len(merged_df)} matched")

                # Gain most specific values; empty or missing overwrites (NaN when only some
                # sample sheets have the column) keep the original value
                for column in merged_df.columns:
                    if column.endswith('_overwrite'):
                        original_column = column.replace('_overwrite', '')
                        if original_column in merged_df.columns:
                            merged_df[original_column] = merged_df[column].replace('', np.nan).fillna(merged_df[original_column])
                
                # Drop the overwrite columns
                merged_df = merged_df.drop(columns=[col for col in merged_df.columns if col.endswith('_overwrite')])

                # Surface water without a biome: biome from the water body name
                if 'water_body_biome' in merged_df.columns:
                    no_biome = merged_df['ENVOEnvironmentBiome'].fillna('') == ''
                    is_surface_water = (merged_df['matrix'] == 'Surface water') & (merged_df['sample_type'] == 'Real Sample')
                    use = no_biome & is_surface_water & (merged_df['water_body_biome'].fillna('') != '')
                    merged_df.loc[use, 'ENVOEnvironmentBiome'] = merged_df.loc[use, 'water_body_biome']
                    merged_df = merged_df.drop(columns=['water_body_biome'])

                # Real samples without a country: the one country named in the dataset record,
                # if no sample of the dataset names another country
                country = dataset_country(metadata_collection_dict, allowedTerm_dict['Country']['allowed_values'])
                if country:
                    if 'Country' not in merged_df.columns:
                        merged_df['Country'] = np.nan
                    known = merged_df['Country'].dropna().astype(str).str.strip()
                    known = set(known[known != ''])
                    if known <= {country}:
                        no_country = merged_df['Country'].fillna('').astype(str).str.strip() == ''
                        merged_df.loc[no_country & (merged_df['sample_type'] == 'Real Sample'), 'Country'] = country

                # Material names that are ENVO synonyms ('ocean water') -> the ENVO label
                merged_df['ENVOEnvironmentMaterial'] = merged_df['ENVOEnvironmentMaterial'].map(
                    lambda v: envo_material_label(v, ENVOEnvironmentMaterialIndex_table, allowedTerm_dict))

                # The biome is the broad-scale environmental context (MIxS env_broad_scale)
                merged_df['ENVOBroadScale'] = merged_df['ENVOEnvironmentBiome']

                # If ENVOEnvironmentMaterial has values other than "" or NA set SampleType to environmental if sample_type is "Real Sample"
                merged_df['SampleType'] = np.where(
                    merged_df['ENVOEnvironmentMaterial'].notna() & (merged_df['ENVOEnvironmentMaterial'] != ''),
                    np.where(merged_df['sample_type'] == 'Real Sample', 'environmental', merged_df['SampleType']),
                    merged_df['SampleType']
                )

                # Combine main DataFrame with spectral raw files
                ########

                # Merge with the spectral raw files DataFrame
                spectral_raw_files_df = pd.concat(nomran_sheets_for_single_dataset_sheets['Spectral raw files'], ignore_index=True)
                spectral_raw_files_df = spectral_raw_files_df.rename(columns={'file_names': 'filename'})

                # Merge with the main DataFrame
                merged_df = pd.merge(merged_df, spectral_raw_files_df, on='sample_id', how='inner')
                print(f"Dataset {internal_id}: {len(spectral_raw_files_df)} raw files, {len(merged_df)} rows after matching raw files")

                # Combine instrument info and instrument setup info
                ########
                instrument_and_setup_df = None
                if combined_instrument_setup_info is not None:
                    if combined_instrument_info is not None and 'instrument_id' in combined_instrument_setup_info.columns:
                        columns_to_drop = [col for col in combined_instrument_info.columns if col in combined_instrument_setup_info.columns and col != 'instrument_id']
                        combined_instrument_info = combined_instrument_info.drop(columns=columns_to_drop)
                        instrument_and_setup_df = pd.merge(combined_instrument_setup_info, combined_instrument_info, on='instrument_id', how='left')
                    else:
                        print(f"Dataset {internal_id}: no instrument info to join, using instrument setups only")
                        instrument_and_setup_df = combined_instrument_setup_info

                # Merge with the main DataFrame (one setup row per sample, so files are not duplicated)
                if instrument_and_setup_df is not None:
                    instrument_and_setup_df = collapse_conflicting_rows(instrument_and_setup_df, 'sample_id')
                    merged_df = pd.merge(merged_df, instrument_and_setup_df, on='sample_id', how='left')

                # Add dataset identifiers
                merged_df['MassiveID'] = 'NORMAN-' + uuid

                dfs.append(merged_df)

            else:
                error_message = f"Failed to download files for dataset {internal_id}. Status code: {metadata_collection_response.status_code}"
                print(error_message)
                errors.append(error_message)
        except Exception as e:
            error_message = f"Error processing dataset {internal_id}: {str(e)}"
            print(error_message)
            traceback.print_exc() 
            errors.append(error_message)


    if dfs:
        
        combined_df = pd.concat(dfs, ignore_index=True)


        combined_df = complete_and_fill_REDU_table(combined_df, allowedTerm_dict, UBERONOntologyIndex_table=ontology_table, ENVOEnvironmentBiomeIndex_table=ENVOEnvironmentBiomeIndex_table,
                                                    ENVOEnvironmentMaterialIndex_table=ENVOEnvironmentMaterialIndex_table,NCBIRankDivision_table=NCBIRankDivision_table, add_usi = False, 
                                                    other_allowed_file_extensions = ['.raw', '.cdf', '.wiff', '.d'], keep_usi = True)
        
        # Make unique by USI
        combined_df = combined_df.drop_duplicates(subset=['USI'], keep='first')
        
        # Make REDU compliant
        combined_df.rename(columns={'MassiveID': 'ATTRIBUTE_DatasetAccession'}, inplace=True)
        combined_df['filename'] = 'f.' + combined_df['filename']


        # Step 4: Save the result as a TSV file
        combined_df.to_csv(output_filename, sep='\t', index=False)
        print(f"\nCombined data saved to {output_filename}")
    else:
        # Still write an (empty) output so the workflow can merge the other sources and
        # carry the previous NORMAN rows forward
        print("No file data found across datasets.")
        open(output_filename, 'w').close()

    # Print skipped datasets and errors
    if skipped:
        print(f"\nSkipped {len(skipped)} datasets without usable sheets:")
        for message in skipped:
            print(message)

    if errors:
        print("\nErrors encountered during processing:")
        for error in errors:
            print(error)

    print(f"\n{len(dfs)} datasets processed, {len(skipped)} skipped, {len(errors)} failed out of {len(datasets)}")

    # Signal a total failure (the output file has already been written)
    if datasets and len(errors) == len(datasets):
        sys.exit(1)



if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Download and combine dataset file information into a TSV file."
        )
    parser.add_argument(
        "-o", "--output",
        type=str,
        required=True,
        help="Output TSV filename (e.g., combined_files.tsv)"
        )
    parser.add_argument(
        "--study_id",
        type=str,
        default="ALL",
        help="UUID of the study to process, or 'ALL' to process all datasets."
        )
    parser.add_argument(
        "--path_to_allowed_term_json", 
        type=str, 
        help="Path to the json with allowed REDU terms"
        )
    parser.add_argument(
        "--path_to_uberon_cl_po_csv", 
        type=str, 
        help="Path to the prepared uberon_cl_po ontology csv"
        )
    parser.add_argument(
        "--path_to_envo_biome_csv", 
        type=str, 
        help="Path to the prepared uberon_cl_po ontology csv"
        )
    parser.add_argument(
        "--path_to_envo_material_csv", 
        type=str, 
        help="Path to the prepared uberon_cl_po ontology csv"
        )
    parser.add_argument(
        "--path_ncbi_rank_division", 
        type=str, 
        help="Path to the path_ncbi_rank_division"
        )
    
    args = parser.parse_args()
    
    path_to_allowed_term_json = args.path_to_allowed_term_json

    # Read allowed terms json
    with open(path_to_allowed_term_json, 'r') as file:
        allowedTerm_dict = json.load(file)

    allowed_values = allowedTerm_dict["MassSpectrometer"]["allowed_values"]
    allowedTerm_dict["MassSpectrometer"]["allowed_values_matching_0"] = [
    norm(v.split('|')[0]) for v in allowed_values]

    # Read ontology
    ontology_table = pd.read_csv(args.path_to_uberon_cl_po_csv)
    ENVOEnvironmentBiomeIndex_table = pd.read_csv(args.path_to_envo_biome_csv)
    ENVOEnvironmentMaterialIndex_table = pd.read_csv(args.path_to_envo_material_csv)

   
    # Read NCBI rank and division csv
    NCBIRankDivision_table = pd.read_csv(args.path_ncbi_rank_division, index_col = False)
    NCBIRankDivision_table = NCBIRankDivision_table.drop_duplicates(subset=['TaxonID'])

    main(args.output, args.study_id,
         allowedTerm_dict = allowedTerm_dict,
         ontology_table = ontology_table,
         ENVOEnvironmentBiomeIndex_table = ENVOEnvironmentBiomeIndex_table,
         ENVOEnvironmentMaterialIndex_table = ENVOEnvironmentMaterialIndex_table,
         NCBIRankDivision_table = NCBIRankDivision_table)
