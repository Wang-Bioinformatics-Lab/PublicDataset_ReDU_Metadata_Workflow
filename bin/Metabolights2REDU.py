import os
import re
import requests
import pandas as pd
from bs4 import BeautifulSoup
import argparse
import json
import time
import traceback
import numpy as np
from tqdm import tqdm
from REDU_conversion_functions import age_category
from REDU_conversion_functions import get_taxonomy_id_from_name__allowedTerms
from REDU_conversion_functions import get_taxonomy_info
from REDU_conversion_functions import merge_repeated_fileobservations
from read_and_validate_redu_from_github import complete_and_fill_REDU_table
from REDU_conversion_functions import find_column_after_target_column
from REDU_conversion_functions import map_instrument_to_allowed
from REDU_conversion_functions import convert_smoking, bmi_to_numeric, convert_diet
from ms_instrument_matching import match_mass_spectrometer




def prefer_extension(group):
    extensions = group['extension'].values
    if 'mzml' in extensions or 'mzxml' in extensions:
        group['keep'] = (group['extension'] == 'mzml') | (group['extension'] == 'mzxml')
    else:
        # Mark only the first row to keep if no preferred extension found
        group['keep'] = [True] + [False] * (len(group) - 1)
    return group


def get_enviromental_water(x):
    x = x.lower()

    if 'water' in x or 'sewerage' in x:
        if 'waste' in x or 'sewerage' in x or 'sewage' in x:
            return ['environmental', 'water_waste']
        if 'surface' in x:
            return ['environmental', 'water_surface']
        if 'ground' in x:
            return ['environmental', 'water_ground']
        if 'storm' in x:
            return ['environmental', 'water_storm']
        if 'sea' in x or 'ocean' in x or 'coast' in x:
            return ['environmental', 'water_seawater']
        else:
            return [None, None]
    else:
        return [None, None]


_AGE_UNIT_DIVISOR = {'year': 1, 'yr': 1, 'month': 12, 'mo': 12, 'week': 52.1429, 'wk': 52.1429, 'day': 365, 'd': 365, 'hour': 8760, 'h': 8760}


def _age_in_years_from_text(value):
    """Age in years from a value that carries its own unit ('73 year', '7 weeks'), else None."""
    m = re.fullmatch(r'\s*(\d+(?:\.\d+)?)\s*(years?|yrs?|months?|mo|weeks?|wks?|days?|d|hours?|h)\s*(?:old)?\s*', str(value), re.I)
    if not m:
        return None
    unit = m.group(2).lower()
    unit = unit if unit in _AGE_UNIT_DIVISOR else unit.rstrip('s')
    return float(m.group(1)) / _AGE_UNIT_DIVISOR[unit]


def _sample_type_qc(x):
    """REDU SampleType for a MetaboLights 'Sample type' value that marks a pooled QC or a
    blank, otherwise None."""
    x = re.sub(r'[^a-z]+', ' ', str(x).lower()).strip()
    if 'pool' in x and ('quality control' in x or re.search(r'\bqc\b', x)):
        return 'pool_QC'
    if 'blank' in x:
        if re.search(r'extraction|process|procedur|method', x):
            return 'blank_extraction'
        if re.search(r'media|medium', x):
            return 'blank_culturemedia'
        if re.search(r'solvent|instrument|analytic|injection|^blank$', x):
            return 'blank_analysis'
    return None


def get_blanks(x):
    x = x.lower()
    if 'blank' in x:
        if 'blank' == x:
            return ['blank_analysis', 'blank_analysis']
        if 'method' in x or 'extraction' in x:
            return ['blank_extraction', 'blank_extraction']
        if 'media' in x:
            return ['blank_culturemedia', 'blank_culturemedia']
        else:
            return [None, None]
    else:
        return [None, None]


class StudyNotAvailable(Exception):
    """The study exists but is not public (MetaboLights returns no content for it)."""


def safe_api_request(url, retries=3, expected_codes={200}):

  for _ in range(retries):
    try:
      response = requests.get(url, timeout=(10, 300))
      print(f"Requested: {url}")
      if response.status_code in expected_codes:
        data = response.json()
        print("Request successful!")
        return data
      elif response.status_code in {401, 403, 404}:
        # private or unknown study: retrying will not help
        print(f"Error: status code {response.status_code}, not retrying")
        return {'content': None, 'message': f'HTTP {response.status_code}'}
      else:
        print(f"Error: Unexpected status code {response.status_code}")
    except Exception as e:
      print(f"Error requesting data: {e}")
    time.sleep(10)
  print(f"All retries failed for {url}.")
  return None

def rename_duplicated_column_names(df):
    # We need this to get rid of duplicated columns
    new_column_names = []

    # Dictionary to keep track of the count of duplicate column names
    dup_count = {}

    # Iterate through the columns using index
    for idx, col in enumerate(df.columns):
        if col in df.columns[:idx]:
            # If the column name is a duplicate, append a count to its name
            dup_count[col] = dup_count.get(col, 0) + 1
            new_name = f"{col}_{dup_count[col]}"
            new_column_names.append(new_name)
        else:
            # If not a duplicate, keep the original name
            new_column_names.append(col)

    # Set the new column names to the DataFrame
    df.columns = new_column_names

    return df


def _table_headers(table, prefix):
    """Column headers for a MetaboLights assay/sample table, in the order of its data rows.

    Field 'index' values are not always contiguous (MetaboLights leaves gaps, e.g. where
    columns without a header were dropped), while the rows hold one value per field in
    index order. So headers are assigned by rank of index, not by the index itself.
    """
    fields = sorted(table['fields'].values(), key=lambda field: field['index'])
    widths = {len(row) for row in table['data']}
    if widths and widths != {len(fields)}:
        if widths == {fields[-1]['index'] + 1}:
            # rows padded to the highest index: position == index
            headers = [None] * widths.pop()
            for field in fields:
                headers[field['index']] = prefix + field['header']
            return headers
        raise ValueError(f"Table rows have {sorted(widths)} values but {len(fields)} fields are declared")
    return [prefix + field['header'] for field in fields]


_ACCESSION_COLUMN = re.compile(r'Term Accession Number(?:[._]\d+)*$')
_ANNOTATION_COLUMN = re.compile(r'(?:Unit|Term Source REF|Term Accession Number)(?:[._]\d+)*$')
_MISSING_VALUES = ('', 'nan', 'none', 'missing value', 'not applicable', 'not available', 'n/a', 'na')


def _is_missing(series):
    return series.isna() | series.astype(str).str.strip().str.lower().isin(_MISSING_VALUES)


def _accession_column(df, value_column):
    """The 'Term Accession Number' column belonging to value_column. In ISA tables it
    follows the value directly, after an optional Unit / Term Source REF column."""
    if value_column not in df.columns:
        return None
    columns = list(df.columns)
    position = columns.index(value_column)
    prefix = value_column.split('_', 1)[0] + '_'
    for column in columns[position + 1:position + 4]:
        if not column.startswith(prefix) or not _ANNOTATION_COLUMN.match(column[len(prefix):]):
            return None
        if _ACCESSION_COLUMN.match(column[len(prefix):]):
            return column
    return None


def _accession_ids(series, ontology):
    """Ontology ids ('<ONTOLOGY>:<id>') from accession values such as
    'http://purl.obolibrary.org/obo/NCBITaxon_9606', '.../ontology/NCBITAXON/9606' (BioPortal)
    or 'UBERON:0001088'; '' if none."""
    pattern = r'(?i)(?:^|[/#:\s])' + ontology + r'[_:/]([A-Za-z0-9]+)\s*$'
    ids = series.astype(str).str.strip().str.extract(pattern, expand=False)
    return ids.fillna('').map(lambda x: f'{ontology}:{x}' if x else '')


_ncbi_by_id = {}


def _ncbi_term_by_id(allowedTerm_dict):
    """taxid -> allowed 'taxid|name' term, built once (the allowed list has ~3M entries)."""
    if not _ncbi_by_id:
        for term in allowedTerm_dict['NCBITaxonomy']['allowed_values']:
            taxid, _, _ = str(term).partition('|')
            _ncbi_by_id.setdefault(taxid, term)
    return _ncbi_by_id


def _investigation_submission_year(study_id):
    """Submission year from the study's public ISA investigation file, used when the
    API has no submissionYear. Returns None if it cannot be read."""
    url = f"https://ftp.ebi.ac.uk/pub/databases/metabolights/studies/public/{study_id}/i_Investigation.txt"
    try:
        response = requests.get(url, timeout=(10, 120))
        if response.status_code != 200:
            return None
    except requests.RequestException:
        return None
    for line in response.text.splitlines():
        for key in ('Study Submission Date', 'Study Public Release Date'):
            if line.startswith(key):
                match = re.search(r'(19|20)\d{2}', line)
                if match:
                    return int(match.group(0))
    return None


def Metabolights2REDU(study_id, **kwargs):
    """
    Converts Metabolights study data to a REDU table format.

    Args:
    study_id: The ID of the Metabolights study to convert.

    Returns:
    A DataFrame in the REDU table format with processed and aligned data from the Metabolights study,
    or an empty DataFrame if no applicable data is found.
    """
    study_url = "https://www.ebi.ac.uk:443/metabolights/ws/studies/public/study/" + study_id
    study_details = safe_api_request(study_url)
    if study_details is None:
        raise RuntimeError(f"Could not retrieve study details for {study_id}")
    if study_details.get('content') is None:
        raise StudyNotAvailable(study_details.get('message', ''))

    study_assays = study_details['content']['assays']

    submissionYear = (study_details['content'].get('derivedData') or {}).get('submissionYear')
    if submissionYear in (None, '', 0):
        submissionYear = _investigation_submission_year(study_id)

    df_assays = pd.DataFrame()
    ms_study_assays = []

    #get study assays if they are MS. There can be multiple assay tables in the same study
    for index, assay_json in enumerate(study_assays):
        if assay_json['technology'] == 'mass spectrometry' or assay_json['technology'] == 'mass spectrometry assay':        

            # Extract headers in the correct order
            headers = _table_headers(assay_json['assayTable'], 'Assay_')

            df = pd.DataFrame(assay_json['assayTable']['data'])
            df = rename_duplicated_column_names(df)

            df.columns = headers

            print(f'Assay {index} has {len(df)} rows and {len(df.columns)} columns.')
            
            ms_study_assays.append(df)

    if len(ms_study_assays) > 0:
        # Union of the columns in first-seen order: the position of a column matters, as
        # an ontology accession sits in the 'Term Accession Number' column after its value.
        # Repeated headers are renamed before the union so none of them is dropped.
        ms_study_assays = [rename_duplicated_column_names(df) for df in ms_study_assays]
        all_columns = list(dict.fromkeys(col for df in ms_study_assays for col in df.columns))

        aligned_dfs = []
        for df in ms_study_assays:
            # Add missing columns with NaN values
            for col in all_columns:
                if col not in df.columns:
                    df[col] = np.nan
            # Reorder columns and add to the aligned list
            aligned_dfs.append(df[all_columns])

        # Concatenate all dataframes into a single dataframe
        df_assays = pd.concat(aligned_dfs, ignore_index=True)

        print(f'Assay table has {len(df_assays)} rows and {len(df_assays.columns)} columns.')

        #get info on samples (there is only one sample table per study)
        headers = _table_headers(study_details['content']['sampleTable'], 'Samples_')

        df_samples = pd.DataFrame(study_details['content']['sampleTable']['data'])
        df_samples = rename_duplicated_column_names(df_samples)
        df_samples.columns = headers
    
        # Sample names are typed by hand in both sheets; ignore surrounding whitespace
        df_assays['Assay_Sample Name'] = df_assays['Assay_Sample Name'].astype(str).str.strip()
        df_samples['Samples_Sample Name'] = df_samples['Samples_Sample Name'].astype(str).str.strip()
        n_unmatched = (~df_assays['Assay_Sample Name'].isin(df_samples['Samples_Sample Name'])).sum()
        if n_unmatched > 0:
            print(f'{n_unmatched} of {len(df_assays)} assay rows have no matching sample in the sample table and are dropped.')

        df_study = df_assays.merge(df_samples, left_on='Assay_Sample Name', right_on='Samples_Sample Name', how='inner')
        df_study['row_id'] = range(1, len(df_study) + 1)

        print(f'Sample table has {len(df_study)} rows and {len(df_study.columns)} columns.')

        # Duplicate rows if we have mzml AND raw files. A study can have several raw/derived file
        # columns (e.g. 'Raw Spectral Data File' holding .wiff.scan and 'Raw Spectral Data File.1'
        # holding the .wiff), so every such column contributes rows; derived (mzML) files come first
        file_columns = [col for prefix in ['Assay_Derived Spectral Data File', 'Assay_Raw Spectral Data File']
                        for col in df_study.columns if re.fullmatch(re.escape(prefix) + r'(\.\d+)?', col)]
        file_dfs = []
        for col in file_columns:
            df_files = df_study[df_study[col].str.contains(r'\.', regex=True, na=False)].copy()
            df_files['filepath'] = df_files[col]
            df_files = df_files.drop(columns=file_columns)
            print(f'Assay table has {len(df_files)} rows and {len(df_files.columns)} columns.')
            if len(df_files) > 0:
                file_dfs.append(df_files)

        if len(file_dfs) == 0:
            return pd.DataFrame()
        df_study = pd.concat(file_dfs, ignore_index=True).copy(deep=True)
        
        df_study = rename_duplicated_column_names(df_study)
        df_study['filename'] = df_study['filepath']#.apply(lambda x: os.path.basename(x))

        # Ensure the filename is a string and lowercase
        df_study['filename_lower'] = df_study['filename'].astype(str).str.lower()

        # Perform the split operation with expansion to get two columns
        split_values = df_study['filename_lower'].str.rsplit('.', n=1, expand=True)

        # List of allowed extensions
        allowed_extensions = [".mzml", ".mzxml", ".cdf", ".raw", ".wiff", ".d"]
        df_study = df_study[df_study['filename_lower'].apply(lambda x: any(x.endswith(ext) for ext in allowed_extensions))]

        print(f'Assay table has {len(df_study)} rows and {len(df_study.columns)} columns.')

        if len(df_study) == 0:
            # e.g. only .zip archives or .xls tables listed, nothing REDU can index
            return pd.DataFrame()

        # Assign the split columns to 'base' and 'extension'
        df_study['base'] = split_values[0]

        # Initially set 'extension' to None to handle rows without an extension
        df_study['extension'] = None

        # Only update 'extension' for rows where a split occurred
        df_study.loc[split_values[1].notna(), 'extension'] = split_values[1]

        df_study = df_study.groupby('row_id').apply(prefer_extension)
        df_study = df_study[df_study['keep']]

        if len(df_study) > 0:
            
            allowedTerm_dict = kwargs['allowedTerm_dict']
            ontology_table = kwargs['ontology_table']
            ENVOEnvironmentMaterialIndex_table = kwargs['ENVOEnvironmentMaterialIndex_table']
            ENVOEnvironmentBiomeIndex_table = kwargs['ENVOEnvironmentBiomeIndex_table']
            NCBIRankDivision_table  = kwargs['NCBIRankDivision_table']

            
            df_study.loc[:, 'YearOfAnalysis'] = submissionYear

            #add NCBITaxonomy and Sampletype & SampleTypeSub1
            #######
            # Coalesce the organism source: the block below is driven by
            # 'Samples_Organism', so studies that record the organism under a
            # different column (Species/Genus/Bacteria/Strain) would otherwise get no
            # NCBITaxonomy at all. Fill 'Samples_Organism' per-sample from those
            # aliases ONLY where it is blank, so an explicit Organism value always
            # wins (no existing value is changed).
            _tax_alias_cols = ['Samples_Species', 'Samples_species', 'Samples_Bacteria',
                               'Samples_Genus', 'Samples_genus', 'Samples_organism']
            # NCBITaxon accession recorded with the organism (taken from the same column
            # as the organism value)
            _tax_accession = pd.Series('', index=df_study.index, dtype=object)
            _acc_col = _accession_column(df_study, 'Samples_Organism')
            if _acc_col:
                _tax_accession = _accession_ids(df_study[_acc_col], 'NCBITaxon')
            if any(c in df_study.columns for c in _tax_alias_cols):
                if 'Samples_Organism' not in df_study.columns:
                    df_study['Samples_Organism'] = np.nan
                for _c in _tax_alias_cols:
                    if _c in df_study.columns:
                        _blank = df_study['Samples_Organism'].isna() | df_study['Samples_Organism'].astype(str).str.strip().str.lower().isin(['', 'nan', 'none', 'missing value'])
                        df_study.loc[_blank, 'Samples_Organism'] = df_study.loc[_blank, _c]
                        _alias_acc = _accession_column(df_study, _c)
                        if _alias_acc:
                            _tax_accession.loc[_blank] = _accession_ids(df_study.loc[_blank, _alias_acc], 'NCBITaxon')

            if 'Samples_Organism' in df_study.columns:
                processed_organisms = {org: str(get_taxonomy_id_from_name__allowedTerms(org, allowedTerm_dict = allowedTerm_dict)) for org in df_study['Samples_Organism'].unique()}

                df_study.loc[:, 'NCBITaxonomy'] = df_study['Samples_Organism'].map(processed_organisms)
                df_study.loc[:, 'NCBITaxonomy'] = df_study['NCBITaxonomy'].replace(to_replace=r'^.*None.*$', value='missing value', regex=True)

                # The curator-chosen NCBITaxon accession is more reliable than resolving the
                # free-text organism name, so it is used whenever it is a valid taxon
                if (_tax_accession != '').any():
                    _by_id = _ncbi_term_by_id(allowedTerm_dict)
                    _from_acc = _tax_accession.str.replace('NCBITaxon:', '', regex=False).map(_by_id)
                    _use = _from_acc.notna()
                    df_study.loc[_use, 'NCBITaxonomy'] = _from_acc[_use]
                
                df_study.loc[:, ['SampleType', 'SampleTypeSub1']] = 'missing value'

                processed_taxonomy = {taxonomy.split('|')[0]: get_taxonomy_info(taxonomy.split('|')[0])
                                    for taxonomy in df_study['NCBITaxonomy'].unique()
                                    if '|' in taxonomy and 'None' not in taxonomy}
                
                print(f'Found {len(processed_taxonomy)} unique taxonomy IDs.')


                df_study.loc[:,'SampleType'] = df_study['NCBITaxonomy'].map(lambda x: processed_taxonomy.get(x.split('|')[0], [pd.NA, pd.NA])[0] 
                                                                    if '|' in x and 'None' not in x else pd.NA)
                df_study.loc[:,'SampleTypeSub1'] = df_study['NCBITaxonomy'].map(lambda x: processed_taxonomy.get(x.split('|')[0], [pd.NA, pd.NA])[1]
                                                                        if '|' in x and 'None' not in x else pd.NA)

                df_study[['SampleType', 'SampleTypeSub1']] = df_study.apply(
                    lambda row: get_blanks(row.Samples_Organism) 
                    if (pd.isna(row.SampleType) or row.SampleType == 'missing value') 
                    else [row.SampleType, row.SampleTypeSub1], 
                    axis=1
                ).apply(pd.Series)

                df_study[['SampleType', 'SampleTypeSub1']] = df_study.apply(
                    lambda row: get_enviromental_water(row.Samples_Organism) 
                    if (pd.isna(row.SampleType) or row.SampleType == 'missing value') 
                    else [row.SampleType, row.SampleTypeSub1], 
                    axis=1
                ).apply(pd.Series)

                # Pooled QCs and blanks are labelled in the 'Sample type' characteristic
                # (MetaboLights vocabulary), while their organism is often that of the study
                # (e.g. a pooled human plasma QC), so this label takes precedence
                for _c in [c for c in df_study.columns if re.fullmatch(r'Samples_Sample type(?:_\d+)?', str(c), re.I)]:
                    _qc = df_study[_c].map(_sample_type_qc)
                    _set = _qc.notna()
                    if _set.any():
                        df_study.loc[_set, 'SampleType'] = _qc[_set]
                        df_study.loc[_set, 'SampleTypeSub1'] = _qc[_set]


                # df_biofluid_vs_tissue = pd.read_csv(os.path.join(transSheet_dir, 'biofluid_vs_tissue_Metabolights.csv'))

                # df_study['Samples_Organism part'] = df_study['Samples_Organism part'].str.lower()

                # for index, row in df_study.iterrows():
                #     if pd.isna(row['SampleTypeSub1']):
                #         org_part = row['Samples_Organism part']
                #         matching_row = df_biofluid_vs_tissue[df_biofluid_vs_tissue['ML'] == org_part]
                #         if not matching_row.empty:
                #             df_study.at[index, 'SampleTypeSub1'] = matching_row.iloc[0]['tissueVSbiofluid']


                #add ENV material column
                #######
                updated_species_dict = {}

                envo_mat_labels_norm = ENVOEnvironmentMaterialIndex_table['Label'].str.replace(' ', '').str.lower()
                envo_mat_synonyms_norm = ENVOEnvironmentMaterialIndex_table['Synonym'].str.replace(' ', '').str.lower()

                for key, value in processed_organisms.items():
                    if value == 'None':
                        match = ENVOEnvironmentMaterialIndex_table[(ENVOEnvironmentMaterialIndex_table['Label'] == key) | (ENVOEnvironmentMaterialIndex_table['Synonym'] == key)]['Label'].values
                        if match.size > 0:
                            updated_species_dict[key] = match[0]
                        else:
                            key_norm = key.replace(' ', '').lower()
                            match_norm = ENVOEnvironmentMaterialIndex_table[(envo_mat_labels_norm == key_norm) | (envo_mat_synonyms_norm == key_norm)]['Label'].values
                            updated_species_dict[key] = match_norm[0] if match_norm.size > 0 else 'missing value'

                df_study.loc[:, 'ENVOEnvironmentMaterial'] = df_study['Samples_Organism'].map(updated_species_dict)


                #add ENV biome column
                #######
                updated_species_dict = {}

                envo_biome_labels_norm = ENVOEnvironmentBiomeIndex_table['Label'].str.replace(' ', '').str.lower()
                envo_biome_synonyms_norm = ENVOEnvironmentBiomeIndex_table['Synonym'].str.replace(' ', '').str.lower()

                for key, value in processed_organisms.items():
                    if value == 'None':
                        match = ENVOEnvironmentBiomeIndex_table[(ENVOEnvironmentBiomeIndex_table['Label'] == key) | (ENVOEnvironmentBiomeIndex_table['Synonym'] == key)]['Label'].values
                        if match.size > 0:
                            updated_species_dict[key] = match[0]
                        else:
                            key_norm = key.replace(' ', '').lower()
                            match_norm = ENVOEnvironmentBiomeIndex_table[(envo_biome_labels_norm == key_norm) | (envo_biome_synonyms_norm == key_norm)]['Label'].values
                            updated_species_dict[key] = match_norm[0] if match_norm.size > 0 else 'missing value'

                df_study.loc[:, 'ENVOEnvironmentBiome'] = df_study['Samples_Organism'].map(updated_species_dict)


                #add UBERON bodypart column
                #######
                # Coalesce the body-part source: only 'Samples_Organism part' was read,
                # so organ/tissue/plant-part recorded under other columns was dropped.
                # Fill it per-sample from those aliases where blank (Organism part wins).
                _bp_alias_cols = ['Samples_Organ', 'Samples_organ', 'Samples_tissue type',
                                  'Samples_Tissue', 'Samples_tissue', 'Samples_Plant part',
                                  'Samples_Biospecimen Type', 'Samples_Cell type']
                if 'Samples_Organism part' not in df_study.columns:
                    df_study['Samples_Organism part'] = np.nan
                # UBERON/CL/PO accession recorded with the body part (from the same column
                # as the body part value)
                def _body_part_accession(column):
                    _acc = _accession_column(df_study, column)
                    if not _acc:
                        return pd.Series('', index=df_study.index, dtype=object)
                    ids = _accession_ids(df_study[_acc], 'UBERON')
                    for _ont in ('CL', 'PO'):
                        ids = ids.where(ids != '', _accession_ids(df_study[_acc], _ont))
                    return ids
                _bp_accession = _body_part_accession('Samples_Organism part')
                for _c in _bp_alias_cols:
                    if _c in df_study.columns:
                        _blank = df_study['Samples_Organism part'].isna() | df_study['Samples_Organism part'].astype(str).str.strip().str.lower().isin(['', 'nan', 'none', 'missing value'])
                        df_study.loc[_blank, 'Samples_Organism part'] = df_study.loc[_blank, _c]
                        _bp_accession.loc[_blank] = _body_part_accession(_c).loc[_blank]

                df_study['UBERONBodyPartName'] = 'missing value'

                labels_in_ontology_table = set(ontology_table['Label'])
                synonyms_in_ontology_table = set(ontology_table['Synonym'])
                allowed_bodyparts = set(allowedTerm_dict['UBERONBodyPartName']['allowed_values'])
                allowed_bodypart_ids = set(allowedTerm_dict['UBERONOntologyIndex']['allowed_values'])
                _label_by_id = dict(zip(ontology_table['UBERONOntologyIndex'].astype(str).str.replace('_', ':', regex=False),
                                        ontology_table['Label']))
                _synonyms_lower = ontology_table['Synonym'].astype(str).str.strip().str.lower()
                _labels_lower = ontology_table['Label'].astype(str).str.strip().str.lower()
                _case_insensitive_matches = {}

                for index, row in df_study.iterrows():
                    organism_part = row['Samples_Organism part']

                    _acc_id = _bp_accession.get(index, '')
                    if _acc_id and _acc_id in allowed_bodypart_ids and _label_by_id.get(_acc_id) in allowed_bodyparts:
                        df_study.at[index, 'UBERONBodyPartName'] = _label_by_id[_acc_id]
                        continue

                    if organism_part in allowed_bodyparts:
                        df_study.at[index, 'UBERONBodyPartName'] = organism_part
                        continue 

                    matching_rows = ontology_table.iloc[0:0]
                    if organism_part in synonyms_in_ontology_table:
                        matching_rows = ontology_table[ontology_table['Synonym'] == organism_part]
                    elif isinstance(organism_part, str) and organism_part.strip():
                        # same label/exact synonym, written with different capitalisation
                        _key = organism_part.strip().lower()
                        if _key not in _case_insensitive_matches:
                            _case_insensitive_matches[_key] = ontology_table[(_labels_lower == _key) | (_synonyms_lower == _key)].drop_duplicates(subset=['UBERONOntologyIndex'])
                        matching_rows = _case_insensitive_matches[_key]

                    if len(matching_rows) > 0:

                        plant = True if row['SampleType'] == 'plant' else False

                        if len(matching_rows) > 1 and plant == True:
                            matching_rows = matching_rows[matching_rows['UBERONOntologyIndex'].str.startswith('PO', na=False)]

                        if len(matching_rows) > 1 and plant == False:
                            matching_rows = matching_rows[(~matching_rows['UBERONOntologyIndex'].str.startswith('PO', na=False))]

                        if len(matching_rows) == 1:
                            
                            label = matching_rows['Label'].iloc[0] 
                            df_study.at[index, 'UBERONBodyPartName'] = label

                        
                ontology_table_unique = ontology_table.drop_duplicates(subset=['Label', 'UBERONOntologyIndex']).drop(columns=['Synonym'])

                df_study = pd.merge(df_study, ontology_table_unique, left_on='UBERONBodyPartName', right_on='Label', how='left')
                df_study.drop('Label', axis=1, inplace=True)

                df_study.loc[df_study['Is Fluid'] == True, 'SampleTypeSub1'] = 'biofluid'
                df_study.loc[df_study['Is Multicellular'] == True, 'SampleTypeSub1'] = 'tissue'

                # UBERON does not flag every body fluid or excretion (urine, feces, follicular
                # fluid, ...) as such; the MetaboLights translation sheet covers those
                _bvt = kwargs.get('biofluid_vs_tissue')
                if _bvt is not None and len(_bvt):
                    _bvt_map = dict(zip(_bvt['ML'].astype(str).str.strip().str.lower(), _bvt['tissueVSbiofluid']))
                    _sheet = df_study['UBERONBodyPartName'].astype(str).str.strip().str.lower().map(_bvt_map)
                    _sheet = _sheet.where(_sheet.notna(), df_study['Samples_Organism part'].astype(str).str.strip().str.lower().map(_bvt_map))
                    _need = _is_missing(df_study['SampleTypeSub1']) & _sheet.notna() & (df_study['SampleType'] == 'animal')
                    df_study.loc[_need, 'SampleTypeSub1'] = _sheet[_need]
                df_study['UBERONOntologyIndex'] = df_study['UBERONOntologyIndex'].str.replace("_", ":", regex=False)


            #add MassSpectrometer column
            #######
            if 'Assay_Instrument' in df_study.columns:
                df_study['Assay_Instrument'] = df_study['Assay_Instrument'].fillna('')
                _instrument_name = df_study['Assay_Instrument'].astype(str)
                _instrument_acc_col = _accession_column(df_study, 'Assay_Instrument')
                _instrument_acc = df_study[_instrument_acc_col].astype(str) if _instrument_acc_col else pd.Series('', index=df_study.index)

                prefixes = [
                    "Thermo Scientific ",
                    "Waters ",
                    "Agilent ",
                    "AB SCIEX ",
                    "LECO ",
                    "Bruker ",
                    "Shimadzu "
                ]

                # Function to remove prefixes
                def remove_prefixes(value, prefixes):
                    for prefix in prefixes:
                        if value.startswith(prefix):
                            return value[len(prefix):]
                    return value

                # Apply the function to the column
                df_study['Assay_Instrument'] = df_study['Assay_Instrument'].apply(lambda x: remove_prefixes(x, prefixes))

                allowed_values_map = dict(zip(
                    allowedTerm_dict["MassSpectrometer"]["allowed_values_matching_0"],
                    allowedTerm_dict["MassSpectrometer"]["allowed_values"]
                ))

                df_study["MassSpectrometer"] = df_study["Assay_Instrument"].map(allowed_values_map)

                # Then the PSI-MS accession recorded with the instrument, then the shared
                # name matcher (vendor prefixes, abbreviations; only unambiguous matches)
                _unres = df_study["MassSpectrometer"].isna()
                if _unres.any():
                    _allowed_ms_terms = allowedTerm_dict["MassSpectrometer"]["allowed_values"]
                    _ms_matches = {}
                    def _match_ms(name, accession):
                        if (name, accession) not in _ms_matches:
                            _ms_matches[(name, accession)] = (match_mass_spectrometer(None, _allowed_ms_terms, accession=accession)
                                                              or match_mass_spectrometer(name, _allowed_ms_terms))
                        return _ms_matches[(name, accession)]
                    df_study.loc[_unres, "MassSpectrometer"] = [
                        _match_ms(n, a) for n, a in zip(_instrument_name[_unres], _instrument_acc[_unres])]

                # Fallback: the exact map misses ~40% of instruments (Q-Exactive vs
                # "Q Exactive", Agilent iFunnel Q-TOF, Synapt/maXis variants, ...). Map
                # the (already prefix-stripped) instrument string onto the allowed
                # vocabulary by token-subset match. Analysis-level, no over-assignment.
                _allowed_ms = allowedTerm_dict["MassSpectrometer"]["allowed_values"]
                _unres = df_study["MassSpectrometer"].isna()
                if _unres.any():
                    _ms_cache = {}
                    def _resolve_ms(raw):
                        if raw not in _ms_cache:
                            _ms_cache[raw] = map_instrument_to_allowed(raw, _allowed_ms)
                        return _ms_cache[raw]
                    df_study.loc[_unres, "MassSpectrometer"] = df_study.loc[_unres, "Assay_Instrument"].map(_resolve_ms)



            #add IonizationMethod/Polarity column
            #######
            if 'Assay_Ion source' in df_study.columns and 'Assay_Scan polarity' in df_study.columns:
                df_study['Assay_Ion source'] = df_study['Assay_Ion source'].fillna('')
                df_study['Assay_Scan polarity'] = df_study['Assay_Scan polarity'].fillna('')
                _ion_acc_col = _accession_column(df_study, 'Assay_Ion source')
                _ion_acc = _accession_ids(df_study[_ion_acc_col], 'MS') if _ion_acc_col else pd.Series('', index=df_study.index)

                df_study['Assay_Ion source'] = df_study['Assay_Ion source'].str.lower()
                df_study['Assay_Scan polarity'] = df_study['Assay_Scan polarity'].astype(str).str.lower()
                df_study['IonizationSourceAndPolarity'] = ''
                df_study.loc[df_study['Assay_Ion source'].str.contains('electrospray'), 'IonizationSourceAndPolarity'] = 'electrospray ionization'
                df_study.loc[df_study['Assay_Ion source'].str.contains('electron ionization'), 'IonizationSourceAndPolarity'] = 'electron ionization'      
                df_study.loc[df_study['Assay_Ion source'].str.contains('chemical ionization'), 'IonizationSourceAndPolarity'] = 'atmospheric pressure chemical ionization'
                # Abbreviations and vendor source names (ESI, HESI, Turbo V/TurboIonSpray,
                # IonDrive, APCI, APPI, EI), then the PSI-MS accession of the ion source
                _unset = df_study['IonizationSourceAndPolarity'] == ''
                _src = df_study['Assay_Ion source']
                for _pattern, _term in [(r'\bn?h?esi\b|nano ?spray|turbo ?(?:v|ion ?spray)|ion ?spray|iondrive|jet stream', 'electrospray ionization'),
                                        (r'\bapci\b', 'atmospheric pressure chemical ionization'),
                                        (r'\bappi\b|photoionization|photo-ionization', 'atmospheric pressure photoionization'),
                                        (r'\bei\b|electron impact', 'electron ionization')]:
                    _hit = _unset & _src.str.contains(_pattern, regex=True)
                    df_study.loc[_hit, 'IonizationSourceAndPolarity'] = _term
                    _unset = _unset & ~_hit
                for _accession, _term in [('MS:1000073', 'electrospray ionization'), ('MS:1000398', 'electrospray ionization'),
                                          ('MS:1000070', 'atmospheric pressure chemical ionization'),
                                          ('MS:1000382', 'atmospheric pressure photoionization'),
                                          ('MS:1000389', 'electron ionization')]:
                    _hit = _unset & (_ion_acc == _accession)
                    df_study.loc[_hit, 'IonizationSourceAndPolarity'] = _term
                    _unset = _unset & ~_hit
                _negative = df_study['Assay_Scan polarity'].str.contains(r'negative|\bneg\b', regex=True)
                _positive = df_study['Assay_Scan polarity'].str.contains(r'posi?tive|\bpos\b', regex=True)
                _with_polarity = (df_study['IonizationSourceAndPolarity'] != '') & ~df_study['IonizationSourceAndPolarity'].str.contains("electron ionization")
                df_study.loc[_negative & _with_polarity, 'IonizationSourceAndPolarity'] += ' (negative)'
                df_study.loc[_positive & _with_polarity, 'IonizationSourceAndPolarity'] += ' (positive)'
                df_study.loc[df_study['Assay_Scan polarity'].str.contains('alternating') & _with_polarity, 'IonizationSourceAndPolarity'] += ' (alternating)'


            #add LC method
            #######
            if 'Assay_Column type' in df_study.columns:
                df_study['Assay_Column type'] = df_study['Assay_Column type'].fillna('')
                df_study['Assay_Column type'] = df_study['Assay_Column type'].str.lower()
                df_study['ChromatographyAndPhase'] = ''
                df_study.loc[df_study['Assay_Column type'].str.contains('reverse'), 'ChromatographyAndPhase'] = 'reverse phase'
                df_study.loc[df_study['Assay_Column type'].str.contains('hilic'), 'ChromatographyAndPhase'] = 'normal phase (HILIC)'
                df_study.loc[df_study['Assay_Column type'].str.contains('normal phase'), 'ChromatographyAndPhase'] = 'normal phase (HILIC)'
                
                if 'Assay_Column model' in df_study.columns:
                    df_study['Assay_Column model'] = df_study['Assay_Column model'].fillna('')
                    df_study.loc[df_study['Assay_Column model'].str.contains('Phenyl', case=False) & df_study['Assay_Column model'].str.contains('Hexyl', case=False), 'ChromatographyAndPhase'] += ' (Phenyl-Hexyl)'
                    df_study.loc[df_study['Assay_Column model'].str.contains('C18') & df_study['Assay_Column model'].str.contains('polar', case=False), 'ChromatographyAndPhase'] += ' (polar-C18)'
                    
                    df_study.loc[df_study['Assay_Column model'].str.contains('HSS T3') & ~df_study['ChromatographyAndPhase'].str.contains(r"\("), 'ChromatographyAndPhase'] += ' (C18)'
                    df_study.loc[df_study['Assay_Column model'].str.contains('C18') & ~df_study['ChromatographyAndPhase'].str.contains(r"\("), 'ChromatographyAndPhase'] += ' (C18)'
                    df_study.loc[df_study['Assay_Column model'].str.contains('C30') & ~df_study['ChromatographyAndPhase'].str.contains(r"\("), 'ChromatographyAndPhase'] += ' (C30)'
                    df_study.loc[df_study['Assay_Column model'].str.contains('C8') & ~df_study['ChromatographyAndPhase'].str.contains(r"\("), 'ChromatographyAndPhase'] += ' (C8)'
                    df_study.loc[df_study['ChromatographyAndPhase'].str.contains('reverse phase') & ~df_study['ChromatographyAndPhase'].str.contains(r"\("), 'ChromatographyAndPhase'] += ' (NOS)'
                else:
                    df_study.loc[df_study['ChromatographyAndPhase'].str.contains('reverse phase') & ~df_study['ChromatographyAndPhase'].str.contains(r"\("), 'ChromatographyAndPhase'] += ' (NOS)'

            # Gas chromatography: GC studies often leave 'Column type' empty, but name a GC
            # in the chromatography instrument or use GC column/type wording
            _gc_text = pd.Series('', index=df_study.index, dtype=object)
            for _c in ['Assay_Chromatography Instrument', 'Assay_Column type', 'Assay_Column model']:
                if _c in df_study.columns:
                    _gc_text = _gc_text + ' ' + df_study[_c].fillna('').astype(str)
            _is_gc = _gc_text.str.contains(r'\bgc\b|gas chromatograph|\bgc[- ]?ms\b|\bgcxgc\b', case=False, regex=True)
            if _is_gc.any():
                if 'ChromatographyAndPhase' not in df_study.columns:
                    df_study['ChromatographyAndPhase'] = ''
                _need = _is_missing(df_study['ChromatographyAndPhase']) & _is_gc
                _db5 = _gc_text.str.contains(r'\bdb-?5', case=False, regex=True)
                df_study.loc[_need & _db5, 'ChromatographyAndPhase'] = 'gas chromatography (DB-5)'
                df_study.loc[_need & ~_db5, 'ChromatographyAndPhase'] = 'gas chromatography (NOS)'

            #add AgeInYears 
            #######
            if 'samples_age' in df_study.columns.str.lower():
                if 'Samples_Age' in df_study.columns:
                    df_study.rename(columns={'Samples_Age': 'Samples_age'}, inplace=True)
                if 'Samples_AGE' in df_study.columns:
                    df_study.rename(columns={'Samples_AGE': 'Samples_age'}, inplace=True)

                _age_text = df_study['Samples_age'].copy() if 'Samples_age' in df_study.columns else pd.Series(np.nan, index=df_study.index)
                unit_column = find_column_after_target_column(df_study, 'Samples_age', 'Samples_Unit')

                if unit_column != '':
                    df_study['Samples_age'] = pd.to_numeric(df_study['Samples_age'], errors='coerce')
                    df_study['AgeInYears'] = np.nan
                    df_study.loc[df_study[unit_column] == 'year', 'AgeInYears'] = df_study.loc[df_study[unit_column] == 'year', 'Samples_age']
                    df_study.loc[df_study[unit_column] == 'month', 'AgeInYears'] = df_study.loc[df_study[unit_column] == 'month', 'Samples_age'] / 12
                    df_study.loc[df_study[unit_column] == 'week', 'AgeInYears'] = df_study.loc[df_study[unit_column] == 'week', 'Samples_age'] / 52.1429
                    df_study.loc[df_study[unit_column] == 'day', 'AgeInYears'] = df_study.loc[df_study[unit_column] == 'day', 'Samples_age'] / 365
                    df_study.loc[df_study[unit_column] == 'hour', 'AgeInYears'] = df_study.loc[df_study[unit_column] == 'hour', 'Samples_age'] / 8760

                    df_study['AgeInYears'] = df_study['AgeInYears'].astype(str).replace('nan', 'missing value')

            #add AgeInYears 
            #######
                    df_study['LifeStage'] = df_study['AgeInYears'].apply(lambda x: age_category(x))

                # The unit can also be written in the value itself ('73 year', '7 weeks')
                _age_from_value = _age_text.map(_age_in_years_from_text)
                if _age_from_value.notna().any():
                    if 'AgeInYears' not in df_study.columns:
                        df_study['AgeInYears'] = 'missing value'
                    _need = _is_missing(df_study['AgeInYears']) & _age_from_value.notna()
                    df_study.loc[_need, 'AgeInYears'] = _age_from_value[_need].astype(str)
                    if 'LifeStage' not in df_study.columns:
                        df_study['LifeStage'] = ''
                    df_study.loc[_need, 'LifeStage'] = _age_from_value[_need].apply(age_category)

            #add Sex
            #######
            # Consume any sex / gender / biological-sex column (the original only read
            # 'Samples_gender', silently dropping 'Samples_sex' and 'Samples_Biological
            # sex'). Exclude columns describing a DIFFERENT subject than the sample
            # (baby/fetus/maternal/paternal/donor/partner/contact).
            _sex_exclude = re.compile(r'baby|fetus|foetal|fetal|maternal|mother|paternal|father|donor|partner|contact|infant|neonat', re.I)
            def _tokens(col):
                return set(re.split(r'[^a-z0-9]+', str(col).lower()))
            _sex_cols = [c for c in df_study.columns
                         if str(c).startswith('Samples_')
                         and ({'sex', 'gender'} & _tokens(c))
                         and not _sex_exclude.search(str(c))]
            if _sex_cols:
                _sexsrc = pd.Series('', index=df_study.index, dtype=object)
                for c in _sex_cols:
                    _b = _sexsrc.astype(str).str.strip().str.lower().isin(['', 'nan', 'none'])
                    _sexsrc.loc[_b] = df_study[c].astype(str).loc[_b]
                _sl = _sexsrc.astype(str).str.lower()
                if 'BiologicalSex' not in df_study.columns:
                    df_study['BiologicalSex'] = 'missing value'
                _blank = df_study['BiologicalSex'].astype(str).str.lower().isin(['', 'nan', 'missing value', 'none'])
                df_study.loc[_blank & _sl.str.contains('female'), 'BiologicalSex'] = 'female'
                df_study.loc[_blank & _sl.str.contains(r'\bmale\b') & (df_study['BiologicalSex'] != 'female'), 'BiologicalSex'] = 'male'


            #add LatitudeandLongitude
            #######
            # Match latitude/longitude by column NAME (paired), covering MetaboLights
            # variants like "Geographic Location Latitude" / "Geographic location
            # (latitude)" - never by value (a bare number could be a plate/time).
            _lat_cols = [c for c in df_study.columns if 'latitude' in str(c).lower()]
            _lon_cols = [c for c in df_study.columns if 'longitude' in str(c).lower()]
            if _lat_cols and _lon_cols:
                _latc, _lonc = _lat_cols[0], _lon_cols[0]
                def combine_lat_lon(row):
                    lat = str(row[_latc]).strip()
                    lon = str(row[_lonc]).strip()
                    if lat and lon and lat.lower() not in ('nan', '') and lon.lower() not in ('nan', ''):
                        try:
                            la = float(lat); lo = float(lon)
                            if -90 <= la <= 90 and -180 <= lo <= 180 and not (la == 0 and lo == 0):
                                return f"{lat}|{lon}"
                        except ValueError:
                            pass
                    return 'missing value'
                _ll = df_study.apply(combine_lat_lon, axis=1)
                if 'LatitudeandLongitude' not in df_study.columns:
                    df_study['LatitudeandLongitude'] = _ll
                else:
                    _blank = df_study['LatitudeandLongitude'].astype(str).str.lower().isin(['', 'nan', 'missing value', 'none'])
                    df_study.loc[_blank, 'LatitudeandLongitude'] = _ll[_blank]

            #add ENVOEnvironmentMaterial (value-matched from candidate columns)
            #######
            # Environmental studies record the material (soil/seawater/sediment/...) in
            # dedicated substrate/material columns. Match values against the ENVO
            # material vocabulary; fill only where currently missing. NOTE: we do NOT
            # scan the organism column here (organism -> ENVO material is handled
            # separately above via the ontology table) - matching organisms against the
            # material vocab risks false hits (e.g. a stray taxon name in the vocab).
            _envm_vocab = {str(v).strip().lower(): str(v).strip()
                           for v in allowedTerm_dict.get('ENVOEnvironmentMaterial', {}).get('allowed_values', []) if str(v).strip()}
            if _envm_vocab:
                _mat_cols = [c for c in ['Samples_Onsite Substrate', 'Samples_Overlaying substrate',
                                         'Samples_Material', 'Samples_material', 'Samples_Environmental material',
                                         'Samples_environment_material', 'Samples_env_material', 'Samples_Substrate',
                                         'Samples_Environment (material)'] if c in df_study.columns]
                if _mat_cols:
                    if 'ENVOEnvironmentMaterial' not in df_study.columns:
                        df_study['ENVOEnvironmentMaterial'] = 'missing value'
                    for _c in _mat_cols:
                        _m = df_study[_c].astype(str).str.strip().str.lower().map(_envm_vocab)
                        _need = df_study['ENVOEnvironmentMaterial'].astype(str).str.strip().str.lower().isin(['', 'nan', 'none', 'missing value']) & _m.notna()
                        df_study.loc[_need, 'ENVOEnvironmentMaterial'] = _m[_need]

            #add HealthStatus (value-matched)
            #######
            _hs_map = {'healthy': 'healthy', 'health': 'healthy', 'unhealthy': 'unhealthy (NOS)',
                       'diseased': 'unhealthy (NOS)', 'chronic illness': 'chronic illness',
                       'acute illness': 'acute illness'}
            _hs_cols = [c for c in df_study.columns
                        if str(c).startswith('Samples_') and re.search(r'health\s*status|health state', str(c), re.I)]
            if _hs_cols:
                if 'HealthStatus' not in df_study.columns:
                    df_study['HealthStatus'] = 'missing value'
                for _c in _hs_cols:
                    _m = df_study[_c].astype(str).str.strip().str.lower().map(_hs_map)
                    _need = df_study['HealthStatus'].astype(str).str.strip().str.lower().isin(['', 'nan', 'none', 'missing value']) & _m.notna()
                    df_study.loc[_need, 'HealthStatus'] = _m[_need]

            #add DOIDCommonName / HealthStatus from disease columns
            #######
            # Values (or their DOID accessions) that exactly match a DOID label or exact
            # synonym; 'healthy control'-type values mark healthy subjects. Columns about
            # someone other than the subject are skipped.
            _disease_cols = [c for c in df_study.columns
                             if str(c).startswith('Samples_')
                             and re.search(r'disease|patholog|diagnos', str(c), re.I)
                             and not re.search(r'severity|stag|grade|score|duration|onset|family|maternal|mother|father|donor|Term Source|Term Accession|Unit', str(c), re.I)]
            if _disease_cols:
                _doid = kwargs.get('doid_table')
                _doid_allowed = set(allowedTerm_dict['DOIDCommonName']['allowed_values'])
                _doid_labels = {str(v).strip().lower(): str(v) for v in _doid_allowed}
                _doid_by_id, _doid_by_synonym = {}, {}
                if _doid is not None and len(_doid):
                    _doid_by_id = dict(zip(_doid['DOIDOntologyIndex'].astype(str), _doid['Label'].astype(str)))
                    _syn = _doid.assign(_s=_doid['Synonym'].astype(str).str.strip().str.lower())
                    _syn = _syn.groupby('_s')['Label'].agg(lambda x: set(x))
                    _doid_by_synonym = {k: next(iter(v)) for k, v in _syn.items() if len(v) == 1}
                def _to_doid(value):
                    key = str(value).strip().lower()
                    if not key or key in _MISSING_VALUES or key == 'disease':
                        return None
                    return _doid_labels.get(key) or (_doid_by_synonym.get(key) if _doid_by_synonym.get(key) in _doid_allowed else None)
                if 'DOIDCommonName' not in df_study.columns:
                    df_study['DOIDCommonName'] = 'missing value'
                if 'HealthStatus' not in df_study.columns:
                    df_study['HealthStatus'] = 'missing value'
                for _c in _disease_cols:
                    _m = df_study[_c].map(_to_doid)
                    _acc_col = _accession_column(df_study, _c)
                    if _acc_col and _doid_by_id:
                        _from_acc = _accession_ids(df_study[_acc_col], 'DOID').map(_doid_by_id)
                        _from_acc = _from_acc.where(_from_acc.isin(_doid_allowed))
                        _m = _from_acc.where(_from_acc.notna(), _m)
                    _need = _is_missing(df_study['DOIDCommonName']) & _m.notna()
                    df_study.loc[_need, 'DOIDCommonName'] = _m[_need]
                    _healthy = df_study[_c].astype(str).str.strip().str.lower().str.fullmatch(
                        r'(?:healthy|non-?diseased?)(?: (?:control|volunteer|subject|donor|individual|participant)s?)?|healthy controls?|control \(healthy\)')
                    _need = _is_missing(df_study['HealthStatus']) & _healthy
                    df_study.loc[_need, 'HealthStatus'] = 'healthy'

            #add DepthorAltitudeMeters
            #######
            if 'Samples_Nominal depth' in df_study.columns:
                def get_depth(val):
                    val = str(val).strip()
                    if val and val not in ('nan', ''):
                        try:
                            return str(-abs(float(val)))
                        except ValueError:
                            pass
                    return 'missing value'
                df_study['DepthorAltitudeMeters'] = df_study['Samples_Nominal depth'].apply(get_depth)

            #add Country
            #######
            # Country was never mapped from MetaboLights. Read geographic-origin
            # columns and match their values against the Country vocabulary
            # (case-insensitive). Sample-level; fills only where currently missing and
            # the value is a recognized country (e.g. "Spain", "Japan").
            _country_cols = ['Samples_Geographical origin', 'Samples_Geographic origin',
                             'Samples_Geographical location', 'Samples_Geographic location',
                             'Samples_Country', 'Samples_country', 'Samples_Country of origin',
                             'Samples_Rice orign', 'Samples_Nationality', 'Samples_Birth country']
            _cc_present = [c for c in _country_cols if c in df_study.columns]
            if _cc_present:
                _cvocab = {str(v).strip().lower(): str(v).strip()
                           for v in allowedTerm_dict.get('Country', {}).get('allowed_values', []) if str(v).strip()}
                if 'Country' not in df_study.columns:
                    df_study['Country'] = 'missing value'
                for _c in _cc_present:
                    _mapped = df_study[_c].astype(str).str.strip().str.lower().map(_cvocab)
                    _need = df_study['Country'].astype(str).str.strip().str.lower().isin(['', 'nan', 'none', 'missing value']) & _mapped.notna()
                    df_study.loc[_need, 'Country'] = _mapped[_need]

            #add SmokingStatus
            #######
            def _toks(col):
                return set(re.split(r'[^a-z0-9]+', str(col).lower()))
            _smoke_cols = [c for c in df_study.columns if str(c).startswith('Samples_') and _toks(c) & {'smoking', 'smoke', 'smoker', 'smokingstatus', 'tobacco'}]
            if _smoke_cols:
                if 'SmokingStatus' not in df_study.columns:
                    df_study['SmokingStatus'] = 'missing value'
                for _c in _smoke_cols:
                    _m = df_study[_c].apply(convert_smoking, key=_c)
                    _need = df_study['SmokingStatus'].astype(str).str.strip().str.lower().isin(['', 'nan', 'none', 'missing value']) & (_m != 'missing value')
                    df_study.loc[_need, 'SmokingStatus'] = _m[_need]

            #add BodyMassIndex
            #######
            _bmi_cols = [c for c in df_study.columns if str(c).startswith('Samples_') and ('bmi' in _toks(c) or {'body', 'mass', 'index'} <= _toks(c))]
            if _bmi_cols:
                if 'BodyMassIndex' not in df_study.columns:
                    df_study['BodyMassIndex'] = 'missing value'
                for _c in _bmi_cols:
                    _m = df_study[_c].apply(bmi_to_numeric).astype(str)
                    _need = df_study['BodyMassIndex'].astype(str).str.strip().str.lower().isin(['', 'nan', 'none', 'missing value']) & (_m != 'missing value')
                    df_study.loc[_need, 'BodyMassIndex'] = _m[_need]

            #add Diet
            #######
            _diet_cols = [c for c in df_study.columns if str(c).startswith('Samples_') and ('diet' in _toks(c) or {'feeding', 'regime'} <= _toks(c) or {'dietary', 'group'} <= _toks(c))]
            if _diet_cols:
                if 'Diet' not in df_study.columns:
                    df_study['Diet'] = 'missing value'
                for _c in _diet_cols:
                    _m = df_study[_c].apply(convert_diet)
                    _need = df_study['Diet'].astype(str).str.strip().str.lower().isin(['', 'nan', 'none', 'missing value']) & (_m != 'missing value')
                    df_study.loc[_need, 'Diet'] = _m[_need]

            #add MassiveID and USIs
            #######
            df_study['MassiveID'] = study_id
            
            ontology_table = ontology_table.drop_duplicates(subset=['Label'])
            df_study = merge_repeated_fileobservations(df_study)
            df_study = complete_and_fill_REDU_table(df_study, allowedTerm_dict, UBERONOntologyIndex_table=ontology_table, ENVOEnvironmentBiomeIndex_table=ENVOEnvironmentBiomeIndex_table,
                                                    ENVOEnvironmentMaterialIndex_table=ENVOEnvironmentMaterialIndex_table,NCBIRankDivision_table=NCBIRankDivision_table, add_usi = True, 
                                                    other_allowed_file_extensions = ['.raw', '.cdf', '.wiff', '.d'])
            
            df_study = df_study.drop_duplicates() 

            #files assigned to several samples were already collapsed into one row by merge_repeated_fileobservations,
            #with conflicting fields set to 'missing value'. Rows can still share a USI here when different filenames
            #resolve to the same file (e.g. several paths inside one .d folder); drop those as we cannot tell them apart
            df_study['count'] = df_study.groupby('USI')['USI'].transform('size')
            df_study = df_study[df_study['count'] == 1]
            df_study = df_study.drop(columns=['count'])

            print(f'We are adding {len(df_study)} samples to the REDU table for study {study_id}.')

            return df_study




if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Give an Metabolights study ID and get a REDU table tsv.')
    parser.add_argument("--study_id", type=str, help='An Metabolights study ID such as "MTBLS1015". If "ALL" all studys are requested.', required=True)
    parser.add_argument("--path_to_translation_sheet_csvs", "-csvs", type=str, help="Path to the translation csvs holding translations from MWB to REDU vocabulary", default="none")
    parser.add_argument("--path_to_allowed_term_json", type=str, help="Path to the json with allowed REDU terms")
    parser.add_argument("--path_to_uberon_cl_po_csv", type=str, help="Path to the prepared uberon_cl_po ontology csv")
    parser.add_argument("--path_to_envo_biome_csv", type=str, help="Path to the prepared uberon_cl_po ontology csv")
    parser.add_argument("--path_to_envo_material_csv", type=str, help="Path to the prepared uberon_cl_po ontology csv")
    parser.add_argument("--path_ncbi_rank_division", type=str, help="Path to the path_ncbi_rank_division")
    parser.add_argument("--path_to_doid_csv", type=str, default="none", help="Optional path to the prepared DOID ontology csv (exact synonyms and DOID accessions for DOIDCommonName)")
            
    args = parser.parse_args()

    if args.study_id == 'ALL':
        public_metabolights_studies = safe_api_request('https://www.ebi.ac.uk:443/metabolights/ws/studies', retries = 5)
        if not public_metabolights_studies or not public_metabolights_studies.get('content'):
            raise RuntimeError("Could not retrieve the list of MetaboLights studies.")
        public_metabolights_studies = public_metabolights_studies['content']
    else:
        public_metabolights_studies = [args.study_id]

    if args.path_to_translation_sheet_csvs == 'none':
        script_dir = os.path.dirname(os.path.realpath(__file__))
        transSheet_dir = os.path.join(script_dir, 'translation_sheets_metabolights')
    else:
        transSheet_dir = args.path_to_translation_sheet_csvs
        script_dir = os.path.dirname(transSheet_dir)

    if args.path_to_allowed_term_json == 'none':
        allowedTermSheet_json = os.path.join(script_dir, 'allowed_terms', 'allowed_terms.json')
    else:
        allowedTermSheet_json = args.path_to_allowed_term_json


    # Read allowed terms json
    with open(allowedTermSheet_json, 'r') as file:
        allowedTerm_dict = json.load(file)

    # Add terms for matching
    allowed_values = allowedTerm_dict["MassSpectrometer"]["allowed_values"]
    allowedTerm_dict["MassSpectrometer"]["allowed_values_matching_0"] = [value.split('|')[0] for value in allowed_values]


    # Read ontology tables
    ontology_table = pd.read_csv(args.path_to_uberon_cl_po_csv)
    ENVOEnvironmentBiomeIndex_table = pd.read_csv(args.path_to_envo_biome_csv)
    ENVOEnvironmentMaterialIndex_table = pd.read_csv(args.path_to_envo_material_csv)

    # Read NCBI rank and division csv
    NCBIRankDivision_table = pd.read_csv(args.path_ncbi_rank_division, index_col = False)
    NCBIRankDivision_table = NCBIRankDivision_table.drop_duplicates(subset=['TaxonID'])

    # Optional tables: DOID ontology, and the body part -> biofluid/tissue translation sheet
    doid_table = pd.read_csv(args.path_to_doid_csv) if args.path_to_doid_csv != 'none' else None
    biofluid_vs_tissue_path = os.path.join(transSheet_dir, 'biofluid_vs_tissue_Metabolights.csv')
    biofluid_vs_tissue = pd.read_csv(biofluid_vs_tissue_path) if os.path.exists(biofluid_vs_tissue_path) else None


    REDU_dataframes = []
    redu_table_single = pd.DataFrame()
    not_public_studies = []
    no_data_studies = []
    failed_studies = []
    for study_id in tqdm(public_metabolights_studies):
        try:
            print(f'Processing study {study_id}...')
            redu_table_single = Metabolights2REDU(study_id, allowedTerm_dict = allowedTerm_dict, ontology_table = ontology_table, ENVOEnvironmentBiomeIndex_table=ENVOEnvironmentBiomeIndex_table,
                                                  ENVOEnvironmentMaterialIndex_table=ENVOEnvironmentMaterialIndex_table, NCBIRankDivision_table=NCBIRankDivision_table,
                                                  doid_table=doid_table, biofluid_vs_tissue=biofluid_vs_tissue)
        except StudyNotAvailable as e:
            print(f"Study {study_id} is not available: {e}")
            not_public_studies.append(study_id)
            continue
        except Exception as e:
            traceback_info = traceback.format_exc()
            print(f"An error occurred with study_id {study_id}: {e}\nTraceback:\n{traceback_info}")
            failed_studies.append(study_id)
            continue
        if redu_table_single is not None and len(redu_table_single) > 0:
            print(f'Added {len(redu_table_single)} samples.')
            REDU_dataframes.append(redu_table_single)
        else:
            print(f'Added {0} samples.')
            no_data_studies.append(study_id)

    print(f'\nStudies with samples: {len(REDU_dataframes)}, without MS samples: {len(no_data_studies)}, '
          f'not public: {len(not_public_studies)}, failed: {len(failed_studies)} (of {len(public_metabolights_studies)})')
    if failed_studies:
        print('Failed studies: ' + ' '.join(failed_studies))

    if len(REDU_dataframes) > 0:
        redu_tables_all = pd.concat(REDU_dataframes, ignore_index=True)
        redu_tables_all = redu_tables_all.drop_duplicates()
        redu_tables_all.to_csv('Metabolights2REDU_' + args.study_id + '.tsv', sep='\t', index=False, header=True)
        print(f'Output of {len(redu_tables_all)} samples has been saved to Metabolights2REDU_{args.study_id}.tsv!')
    else:
        print('nothing to return!')

