import json
import os
import re
import requests
import time
import functools
from datetime import datetime
from bs4 import BeautifulSoup
from owlready2 import get_ontology
import owlready2
import owlready2.namespace as _owl_ns
import pandas as pd
import tqdm

# owlready2 raises TypeError when an OWL file (e.g. new cl.owl via STATO import)
# declares an IRI as both a property and a class/individual. Patch _load_properties
# to skip only those conflicting entries instead of crashing. (Catching the TypeError
# around the original method would also skip loading every property after the
# conflicting one.) Mirrors Ontology._load_properties of owlready2 0.46.
def _patched_load_properties(self):
    # Update props from other ontologies, if needed
    for prop in list(self.world._props.values()):
        if prop.namespace.world is _owl_ns.owl_world: continue
        prop._check_update(self)

    # Loads new props
    for prop_storid in _owl_ns.itertools.chain(self._get_obj_triples_po_s(_owl_ns.rdf_type, _owl_ns.owl_object_property),
                                               self._get_obj_triples_po_s(_owl_ns.rdf_type, _owl_ns.owl_data_property),
                                               self._get_obj_triples_po_s(_owl_ns.rdf_type, _owl_ns.owl_annotation_property)):
        Prop = self.world._get_by_storid(prop_storid)
        if not isinstance(Prop, _owl_ns.PropertyClass):
            print(f"Skipping '{Prop.iri}': it is both a property and a class/an individual")
            continue
        python_name_d = self.world._get_data_triple_sp_od(prop_storid, _owl_ns.owlready_python_name)
        if python_name_d is not None:
            with _owl_ns.LOADING: Prop.python_name = python_name_d[0]
_owl_ns.Ontology._load_properties = _patched_load_properties



def _is_missing_observation(value):
    return value is None or (isinstance(value, float) and pd.isna(value)) or str(value).strip().lower() in ('', 'nan', 'missing value')


def merge_repeated_fileobservations(df):
    """Collapses rows describing the same file into one.

    A file listed for several samples or assays cannot be attributed to one of them, so
    columns with conflicting values become 'missing value'. Blank/missing entries are not
    a conflict: if exactly one real value is present it is kept (e.g. a sample listed twice,
    once with 'France' and once with an empty country).
    """

    def process_group(group):

        for col in group.columns:
            present = [value for value in group[col] if not _is_missing_observation(value)]
            if len({str(value) for value in present}) == 1:
                group[col] = present[0]
            elif present or len(set(group[col])) > 1 or all(value is None for value in group[col]):
                group[col] = 'missing value'
            else:
                # If there's no discrepancy, keep the original value
                group[col] = group[col].iloc[0]
        
        return group

    # Group by 'USI' and apply the processing function to each group
    df = df.groupby('filename').apply(process_group)

    # Remove duplicates based on 'filename' column
    df = df.drop_duplicates(subset='filename')

    return df

def find_column_after_target_column(df, target_column='', search_column_prefix='Samples_Unit'):
    if target_column in df.columns:
        target_index = df.columns.get_loc(target_column)
        # Check the next 3 columns after the target column
        for i in range(1, 4):
            if target_index + i < len(df.columns):
                next_column = df.columns[target_index + i]
                # Check if the next column starts with the search_column_prefix
                if next_column.startswith(search_column_prefix):
                    return next_column
    return ''  # Return None if no matching column is found

def age_category(age):
    if age is None:
        return ''
    try:
        age = float(age)
    except ValueError:
        return ''
    if age < 2:
        return 'Infancy (<2 yrs)'
    elif age <= 8:
        return 'Early Childhood (2 yrs < x <=8 yrs)'
    elif age <= 18:
        return 'Adolescence (8 yrs < x <= 18 yrs)'
    elif age <= 45:
        return 'Early Adulthood (18 yrs < x <= 45 yrs)'
    elif age <= 65:
        return 'Middle Adulthood (45 yrs < x <= 65 yrs)'
    elif age > 65:
        return 'Later Adulthood (>65 yrs)'
    else:
        return ''
    

# NCBI allows 3 requests/s without an API key; going faster returns HTTP 429
_NCBI_MIN_INTERVAL = 0.4
_ncbi_last_request = [0.0]
NCBI_TIMEOUT = (10, 60)


def _ncbi_throttle():
    wait = _NCBI_MIN_INTERVAL - (time.time() - _ncbi_last_request[0])
    if wait > 0:
        time.sleep(wait)
    _ncbi_last_request[0] = time.time()


def _ncbi_get(url, retries=4):
    """GET an NCBI URL, throttled, retrying on throttling, server and connection errors.

    Returns the response if it has status 200, raises otherwise.
    """
    for attempt in range(1, retries + 1):
        _ncbi_throttle()
        try:
            response = requests.get(url, timeout=NCBI_TIMEOUT)
        except (requests.ConnectionError, requests.Timeout) as e:
            if attempt == retries:
                raise
            print(f"NCBI request failed (attempt {attempt}/{retries}): {e}")
        else:
            if response.status_code == 200:
                return response
            if attempt == retries or (response.status_code != 429 and response.status_code < 500):
                response.raise_for_status()
            print(f"NCBI responded with status {response.status_code} (attempt {attempt}/{retries})")
        time.sleep(min(2 ** attempt, 30))


@functools.lru_cache(maxsize=None)
def _get_ncbi_lineage(ncbi_id, retries=3):
    """Returns the lowercased lineage of an NCBI taxon id, or None if NCBI has no lineage for it.

    Cached, as the same ids recur across studies. Network errors propagate (and are not cached).
    """
    for attempt in range(1, retries + 1):
        response = _ncbi_get(f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=taxonomy&id={ncbi_id}")
        soup = BeautifulSoup(response.text, "xml")
        taxon = soup.find("Taxon")
        lineage = taxon.find("Lineage") if taxon is not None else None
        if lineage is not None:
            return tuple(s.lower() for s in lineage.text.split("; "))
        # NCBI occasionally answers 200 with an error document; retry before giving up
        print(f"No lineage in NCBI response for {ncbi_id} (attempt {attempt}/{retries})")
        time.sleep(2)
    return None


@functools.lru_cache(maxsize=None)
def _get_ncbi_scientific_name(ncbi_id):
    """Cached esummary lookup; network errors propagate and are not cached."""
    response = _ncbi_get(f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=taxonomy&id={ncbi_id}&retmode=json")
    summary = response.json().get("result", {}).get(ncbi_id, {})
    return summary.get("scientificname") or None


def get_taxonomic_name_from_id(ncbi_id):
    """Returns the scientific name of an NCBI taxon id, or None if NCBI does not know it.

    Uses the E-utilities esummary endpoint; the Taxonomy Browser page this used to scrape
    no longer carries the name in its title.
    """
    try:
        return _get_ncbi_scientific_name(str(ncbi_id).strip())
    except (requests.RequestException, ValueError) as e:
        print(f"Could not get the taxonomic name for {ncbi_id}: {e}")
        return None


def get_taxonomy_info(ncbi_id, cell_culture_key1 = '', cell_culture_key2 = ''):
    if ncbi_id is not None and ncbi_id != "NA":
        cell_culture_key_words = ["cell", "media", "culture"]
        try:
            #try to get taxa via API
            classification = _get_ncbi_lineage(str(ncbi_id))
            if classification is None:
                raise Exception(f"No lineage found for NCBI id {ncbi_id}.")

            SampleType = None
            SampleTypeSub1 = None

            if "viridiplantae" in classification:
                SampleType = "plant"
                SampleTypeSub1 = "plant_NOS"
                if "Algae" in classification or 'rhodophyta' in classification or "phaeophyceae" in classification:
                    SampleType = "algae"
                if "chlorophyta" in classification or "microalgae" in classification or "microalga" in classification:
                    SampleType = "microalgae"
                if "streptophyta" in classification:
                    SampleTypeSub1 = "plant_angiospermae"
                if "cyanobacteria" in classification:
                    SampleTypeSub1 = "marine_cyanobacteria_insitu"
                if "bacillariophyta" in classification:
                    SampleTypeSub1 = "marine_diatom"
            elif "metazoa" in classification:
                SampleType = "animal"
                if "mammalia" in classification and any(
                        word in cell_culture_key1.lower() for word in cell_culture_key_words) or any(
                        word in cell_culture_key2.lower() for word in cell_culture_key_words):
                    SampleType = "culture_mammalian"
                    SampleTypeSub1 = "culture_mammalian"
                if "amphibia" in classification:
                    if "Caudata" in classification or "urodela" in classification or "echinodermata" in classification:
                        SampleTypeSub1 = "salamander"
                    else:
                        SampleTypeSub1 = "frog"
                if "insecta" in classification:
                    SampleTypeSub1 = "insect"
                if "porifera" in classification or "mollusca" in classification:
                    SampleTypeSub1 = "marine_invertebrates"
                if "cnidaria" in classification:
                    SampleTypeSub1 = "marine_coral"
            elif "fungi" in classification:
                SampleType = "culture_fungal"
                SampleTypeSub1 = "culture_fungal"
            elif "bacteria" in classification:
                SampleType = "culture_bacterial"
                SampleTypeSub1 = "culture_bacterial"
            return [SampleType, SampleTypeSub1]
        except Exception as e:
            print(f"Could not get taxonomy info for {ncbi_id}: {e}")
            return [None, None]
    else:
        return [None, None]

# Vendor / marketing words that carry no model information; dropped before matching an
# instrument string against the allowed MassSpectrometer vocabulary. Shared by the MWB
# and MetaboLights converters.
_MS_STOPWORDS = {
    'thermo', 'abi', 'agilent', 'bruker', 'waters', 'sciex', 'leco', 'shimadzu', 'ab',
    'daltonics', 'hybrid', 'series', 'system', 'lc', 'ms', 'lcms', 'msms', 'scientific',
    'fisher', 'the', 'and', 'with', 'gc',
}


def _ms_tokens(s):
    s = str(s).lower().replace('q-exactive', 'q exactive').replace('qexactive', 'q exactive')
    s = re.sub(r'[^a-z0-9+]+', ' ', s)
    return {t for t in s.split() if t and t not in _MS_STOPWORDS}


def convert_smoking(x):
    """Normalize a free-text smoking value to the SmokingStatus vocabulary
    (current/former/never smoker). Only unambiguous text is mapped; coded values
    (0/1/yes/no) stay missing to avoid guessing polarity. Separators are normalized
    to spaces first so word boundaries match (e.g. 'never_used')."""
    x = re.sub(r'[^a-z0-9]+', ' ', str(x).strip().lower()).strip()
    if re.search(r'\bnever\b', x) or re.search(r'\bnon smoker\b|\bnonsmoker\b', x):
        return 'never smoker'
    if re.search(r'\bformer\b|\bex smoker\b|\bprevious\b|\bpast\b|\bquit\b', x):
        return 'former smoker'
    if re.search(r'\bcurrent\b', x) or x == 'smoker':
        return 'current smoker'
    return 'missing value'


def convert_diet(x):
    """Normalize a free-text diet value (from a diet-named field) to the Diet
    vocabulary. Human patterns and experimental/composition diet types; the input
    is assumed to come from a diet field, so short codes (cd/hf/hsd) decode safely.
    Bespoke study-specific formulas (e.g. 'ser/gly free') stay missing. Order is
    significant: most specific compositions first."""
    x = re.sub(r'[^a-z0-9]+', ' ', str(x).strip().lower()).strip()
    if not x:
        return 'missing value'
    # human dietary patterns
    if re.search(r'\bvegan\b', x):
        return 'vegan'
    if re.search(r'\bpesc[ae]tarian\b', x):
        return 'pescatarian'
    if re.search(r'\bvegetarian\b|lacto ovo|ovo vegetarian', x):
        return 'vegetarian'
    if re.search(r'\bomnivor', x):
        return 'omnivore'
    # experimental / composition diets (specific first)
    if re.search(r'\bhfhs\b|\bhchf\b|high fat high (sugar|sucrose|carb)', x):
        return 'high-fat high-sugar diet'
    if re.search(r'\bhfd\b|\bhf\b|high[ _]?fat|highfat', x):
        return 'high-fat diet'
    if re.search(r'\bhsd\b|high[ _]?sugar|high[ _]?sucrose|\bh sucrose\b|high carb', x):
        return 'high-sugar diet'
    if re.search(r'\blfd\b|low[ _]?fat', x):
        return 'low-fat diet'
    if re.search(r'high[ _]?protein|\bh protein\b', x):
        return 'high-protein diet'
    if re.search(r'ketogenic|\bketo\b', x):
        return 'ketogenic diet'
    if re.search(r'mediterranean', x):
        return 'Mediterranean diet'
    if re.search(r'western|regular american', x):
        return 'Western diet'
    if re.search(r'prebiotic', x):
        return 'prebiotic diet'
    if re.search(r'calorie restrict|caloric restrict|calorie restricted|\brestricted\b', x):
        return 'calorie-restricted diet'
    # control / standard chow
    if re.search(r'\bchow\b|\bnormal\b|\bcontrol\b|\bstandard\b|\bncd\b|\bcd\b|\bcon\b|\bctrl\b|ain ?76|regular chow|standard chow|\bsc\b|normal diet', x):
        return 'control diet'
    return 'missing value'


def bmi_to_numeric(x):
    """Extract a plausible BMI number (kg/m^2). Rejects z-scores, codes and other
    non-BMI values by requiring a human-plausible range (10-60); the same range is
    enforced centrally in complete_and_fill_REDU_table for all sources."""
    nums = re.findall(r'-?\d+\.?\d*', str(x))
    if not nums:
        return 'missing value'
    try:
        v = float(nums[0])
    except ValueError:
        return 'missing value'
    return v if 10 <= v <= 60 else 'missing value'


def map_instrument_to_allowed(instrument, allowed_values):
    """Map a free-text instrument name (e.g. "Thermo Scientific Q-Exactive") onto an
    allowed MassSpectrometer vocabulary entry ("Q Exactive|MS:1001911") by token-subset
    match: every (non-stopword) token of a vocab label must be present in the instrument
    string. The most specific label wins (most tokens), so "HF" never steals "HF-X".
    Returns the allowed value or None if nothing matches confidently."""
    it = _ms_tokens(instrument)
    if not it:
        return None
    best, best_n = None, 0
    for val in allowed_values:
        lab_tokens = _ms_tokens(str(val).split('|')[0])
        if lab_tokens and lab_tokens.issubset(it) and len(lab_tokens) > best_n:
            best, best_n = val, len(lab_tokens)
    return best


def get_taxonomy_id_from_name__allowedTerms(organism_name, **kwargs):

    allowedTerm_dict = kwargs['allowedTerm_dict']
    taxonomy_data = allowedTerm_dict["NCBITaxonomy"]["allowed_values"]
    
    ncbi_id_input = ''
    if 'ncbi_id' in kwargs.keys():
        ncbi_id_input = str(kwargs['ncbi_id'].pop())

    try:
        ncbi_id_numeric = int(ncbi_id_input)
    except ValueError:
        # ncbi_id_input cannot be made numeric, return None or an appropriate error message
        ncbi_id_input = ''

    
    if ncbi_id_input != '':
        ncbi_id_input = str(ncbi_id_input)
        for entry in taxonomy_data:
            parts = entry.split('|')
            if len(parts) == 2:
                ncbi_id, name = parts
                if str(ncbi_id) == ncbi_id_input:
                    return ncbi_id + '|' + str(name)
            else:
                continue


    organism_name = str(organism_name)
    for entry in taxonomy_data:
        parts = entry.split('|')
        if len(parts) == 2:
            ncbi_id, name = parts
            if name.lower() == organism_name.lower():
                return ncbi_id + '|' + str(name)
        else:
            continue

    return None
    req_ncbi_name = get_taxonomy_id_from_name(organism_name)

    if req_ncbi_name is not None:

        req_name = get_taxonomic_name_from_id(req_ncbi_name)
        return ncbi_id + '|' + str(req_name)

    return None


def get_taxonomy_id_from_name(species_name, retries=3):
    if species_name is None or species_name in ["NA", "N/A"]:
        return None

    attempts = 0
    while attempts < retries:
        try:
            _ncbi_throttle()
            response = requests.get(
                f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=taxonomy&term={species_name}&retmode=xml",
                timeout=NCBI_TIMEOUT
            )

            if response.status_code == 200:
                soup = BeautifulSoup(response.text, "xml")
                id_list = soup.find("IdList")
                if id_list is not None and id_list.find("Id") is not None:
                    ncbi_id = id_list.find("Id").text
                    term = soup.find("Term").text.split("[")[0].strip()
                    if term.lower() == species_name.lower():
                        return str(ncbi_id) + '|' + species_name
                    else:
                        return None
                else:
                    return  None
            else:
                print(f"Server responded with status code {response.status_code} for {species_name}.")
                time.sleep(10)

        except Exception as e:
            print(f"Attempt {attempts + 1} failed for {species_name}: {e}")
            time.sleep(10)

        attempts += 1

    print(f'{species_name} returned no NCBI-ID after {retries} attempts')
    return  None



def get_uberon_table(owl_path):
    onto = get_ontology(owl_path).load()

    multi_cellular = onto.search(iri="http://purl.obolibrary.org/obo/UBERON_0010000")[0]
    organ = onto.search(iri="http://purl.obolibrary.org/obo/UBERON_0000062")[0]
    bodily_fluid = onto.search(iri="http://purl.obolibrary.org/obo/UBERON_0006314")[0]

    downstream_classes_multi_cellular = set(multi_cellular.descendants())
    downstream_classes_organ = set(organ.descendants())
    downstream_classes_bodily_fluid = set(bodily_fluid.descendants())
            
    data = []
    for cls in onto.classes():
        label = cls.label.first() if cls.label else None
        synonyms = [synonym for synonym in cls.hasExactSynonym] if hasattr(cls, 'hasExactSynonym') else []

        if not synonyms:
            synonyms = [label] if label else []
            
        is_multi_cellular = cls in downstream_classes_multi_cellular
        is_organ = cls in downstream_classes_organ
        is_bodily_fluid = cls in downstream_classes_bodily_fluid

        for synonym in synonyms:
            data.append({
                "UBERONOntologyIndex": cls.name,
                "Label": label,
                "Synonym": synonym,
                "Is Multicellular": is_multi_cellular,
                "Is Organ": is_organ,
                "Is Fluid": is_bodily_fluid
            })

    df = pd.DataFrame(data)

    df = df[(df["Label"].notna() | df["Synonym"].notna()) & ((df["Label"] != '') | (df["Synonym"] != ''))]
    df = df[df['UBERONOntologyIndex'].str.startswith('UBERON_', na=False)]

    return df



def get_ontology_table(owl_path, ont_prefix, rm_synonym_info=False, descendant_node=None, index_column_name = 'UBERONOntologyIndex'):
    onto = get_ontology(owl_path).load()

    filter_class = None
    descendants = set()
    if descendant_node:
        # Search for the specific class using its node ID
        filter_class = onto.search_one(iri="*" + descendant_node)
        if filter_class:
            # Get all its descendants
            descendants = set(filter_class.descendants())
            # Optionally remove the filter_class itself from the set of descendants to exclude it
            descendants.discard(filter_class)

    data = []
    for cls in tqdm.tqdm(onto.classes(), desc="Processing classes"):
        # Skip the class if it is the filter_class itself and descendant_node is provided
        if descendant_node and cls == filter_class:
            continue

        label = cls.label.first() if cls.label else None
        synonyms = [synonym for synonym in cls.hasExactSynonym] if hasattr(cls, 'hasExactSynonym') else []

        if not synonyms:
            synonyms = [label] if label else []

        # Determine if the current class is a descendant of the specified node
        is_descendant = cls in descendants

        for synonym in synonyms:
            data.append({
                index_column_name: cls.name,
                "Label": label,
                "Synonym": synonym,
                "Is Descendant": is_descendant  # Add this information to each entry
            })

    df = pd.DataFrame(data)
    df = df[(df["Label"].notna() | df["Synonym"].notna()) & ((df["Label"] != '') | (df["Synonym"] != ''))]
    df = df[df[index_column_name].str.startswith(ont_prefix, na=False)]

    if rm_synonym_info:
        df['Synonym'] = df['Synonym'].str.replace(r' \([^)]*\)', '', regex=True)

    # Optionally, you might want to filter the dataframe to include only the descendants
    if descendant_node:
        df = df[df["Is Descendant"] == True]

    return df


# ---------------------------------------------------------------------------
# SampleCollectionDateandTime normalisation
# ---------------------------------------------------------------------------
# Submitters use a wide diversity of date formats and precisions: year only,
# year-month, full date, date+time, textual months, US M/D/Y, ISO Y-M-D, 2-digit
# years, ISO with timezone, and ranges. normalize_redu_datetime() accepts all of
# these and emits ISO-8601 at the SAME precision that was submitted (it never
# fabricates a day or a time that was not given):
#   2018 | 2018-07 | 2018-07-12 | 2018-07-12T22:00 | 2018-07-12T22:00:00
# and intervals as "<start>/<end>" (e.g. 2022-06-14/2022-06-30). Junk / missing
# markers / genuinely unparseable values return None so the caller can substitute
# the column's missing value.

_DATE_MONTHS = {}
for _i, _names in enumerate(
    [('jan', 'january'), ('feb', 'february'), ('mar', 'march'), ('apr', 'april'),
     ('may',), ('jun', 'june'), ('jul', 'july'), ('aug', 'august'),
     ('sep', 'sept', 'september'), ('oct', 'october'), ('nov', 'november'),
     ('dec', 'december')], start=1):
    for _n in _names:
        _DATE_MONTHS[_n] = _i

_DATE_MISSING_TOKENS = {
    '', 'nan', 'nat', 'none', 'null', 'na', 'n/a', 'missing value', 'not applicable',
    'not collected', 'not specified', 'not_applicable', 'ml import: not available',
    'unknown', 'unspecified', 'not available', 'not known', 'not provided',
    'not recorded', 'tbd',
}


def _date_yy(y):
    """Expand a 2-digit year: 00-68 -> 2000-2068, 69-99 -> 1969-1999."""
    y = int(y)
    if y < 100:
        y = 2000 + y if y <= 68 else 1900 + y
    return y


def _date_valid(y, m, d, H, M, S):
    try:
        datetime(y, m or 1, d or 1, H or 0, M or 0, S or 0)
        return 1900 <= y <= 2099
    except (ValueError, TypeError):
        return False


def _date_fmt(y, m, d, H, M, S, prec):
    if prec == 'year':
        return "%04d" % y
    if prec == 'month':
        return "%04d-%02d" % (y, m)
    if prec == 'day':
        return "%04d-%02d-%02d" % (y, m, d)
    if prec == 'minute':
        return "%04d-%02d-%02dT%02d:%02d" % (y, m, d, H, M)
    if prec == 'second':
        return "%04d-%02d-%02dT%02d:%02d:%02d" % (y, m, d, H, M, S)


def _parse_single_date(s, dayfirst=False):
    s = str(s).strip()
    low = s.lower()
    if low in _DATE_MISSING_TOKENS or not s or set(s) <= set('#?-/ .:'):
        return None

    # optional trailing time (':' or '.' separator), preceded by space/T/_, with an
    # optional ISO timezone designator. The time prefix is required so a bare
    # "04-04-2023" is never mistaken for a tz offset.
    H = M = S = None
    tm = re.search(r'[ T_](\d{1,2})[:.](\d{2})(?:[:.](\d{2}))?(?:\s*(?:Z|[+-]\d{2}:?\d{2}))?\s*$', s)
    if tm:
        H = int(tm.group(1)); M = int(tm.group(2))
        S = int(tm.group(3)) if tm.group(3) else None
        if H > 23 or M > 59 or (S is not None and S > 59):
            H = M = S = None
        else:
            s = s[:tm.start()].strip()

    y = m = d = None
    prec = None
    mo = None
    if re.fullmatch(r'\d{4}', s):
        y = int(s); prec = 'year'
    elif (mo := re.fullmatch(r'(\d{4})[-/.](\d{1,2})[-/.](\d{1,2})', s)):
        y, m, d = int(mo.group(1)), int(mo.group(2)), int(mo.group(3)); prec = 'day'
    elif (mo := re.fullmatch(r'(\d{4})[-/.](\d{1,2})', s)):
        y, m = int(mo.group(1)), int(mo.group(2)); prec = 'month'
    elif (mo := re.fullmatch(r'(\d{1,2})[/-](\d{1,2})[/-](\d{2,4})', s)):
        a, b, c = int(mo.group(1)), int(mo.group(2)), _date_yy(mo.group(3))
        if a > 12 and b <= 12:            # unambiguously D/M/Y
            d, m = a, b
        elif b > 12 and a <= 12:          # unambiguously M/D/Y
            m, d = a, b
        elif dayfirst:                     # ambiguous -> caller's derived preference
            d, m = a, b
        else:                              # ambiguous -> platform default (US M/D/Y)
            m, d = a, b
        y = c; prec = 'day'
    elif (mo := re.fullmatch(r'(\d{1,2})[/-](\d{4})', s)):
        m, y = int(mo.group(1)), int(mo.group(2)); prec = 'month'
    elif (mo := re.fullmatch(r'(\d{1,2})[-\s]([A-Za-z]{3,})[-\s,]+(\d{2,4})', s)):
        mon = _DATE_MONTHS.get(mo.group(2).lower()) or _DATE_MONTHS.get(mo.group(2).lower()[:3])
        if mon:
            d, m, y = int(mo.group(1)), mon, _date_yy(mo.group(3)); prec = 'day'
    elif (mo := re.fullmatch(r'([A-Za-z]{3,})[-\s,]+(\d{1,2})[-\s,]+(\d{2,4})', s)):
        mon = _DATE_MONTHS.get(mo.group(1).lower()) or _DATE_MONTHS.get(mo.group(1).lower()[:3])
        if mon:
            m, d, y = mon, int(mo.group(2)), _date_yy(mo.group(3)); prec = 'day'
    elif (mo := re.fullmatch(r'([A-Za-z]{3,})[-\s,]+(\d{2,4})', s)):
        mon = _DATE_MONTHS.get(mo.group(1).lower()) or _DATE_MONTHS.get(mo.group(1).lower()[:3])
        if mon:
            m, y = mon, _date_yy(mo.group(2)); prec = 'month'

    if y is None:
        return None
    if H is not None and prec == 'day':
        prec = 'second' if S is not None else 'minute'
    if not _date_valid(y, m, d, H or 0, M or 0, S or 0):
        return None
    return _date_fmt(y, m, d, H or 0, M or 0, S or 0, prec)


def normalize_redu_datetime(raw, dayfirst=False):
    """Normalise a submitted collection date/time to precision-preserving ISO-8601,
    or return None if it is missing / junk / unparseable.

    dayfirst disambiguates ONLY genuinely ambiguous slash dates (both fields <= 12);
    a value that itself proves an order (a field > 12) always wins over dayfirst.
    """
    if raw is None:
        return None
    s = str(raw).strip()
    if not s:
        return None
    parts = re.split(r'\s+-\s+|\s+to\s+', s)
    if len(parts) == 1:
        mo = re.fullmatch(r'(\d{1,2}/\d{1,2}/\d{2,4})\s*-\s*(\d{1,2}/\d{1,2}/\d{2,4})', s)
        if mo:
            parts = [mo.group(1), mo.group(2)]
    if len(parts) == 2:
        a, b = _parse_single_date(parts[0], dayfirst), _parse_single_date(parts[1], dayfirst)
        if a and b:
            return a + '/' + b
        return a or b
    return _parse_single_date(s, dayfirst)


if __name__ == '__main__':
    pass