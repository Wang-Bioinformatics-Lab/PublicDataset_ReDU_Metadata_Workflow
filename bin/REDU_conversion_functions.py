import json
import os
import re
import requests
import time
from datetime import datetime
from bs4 import BeautifulSoup
from owlready2 import get_ontology
import owlready2
import owlready2.namespace as _owl_ns
import pandas as pd
import tqdm

# owlready2 raises TypeError when an OWL file (e.g. new cl.owl via STATO import)
# declares an IRI as both a property and a class/individual. Patch _load_properties
# to skip those conflicting entries instead of crashing.
_original_load_properties = _owl_ns.Ontology._load_properties
def _patched_load_properties(self):
    try:
        _original_load_properties(self)
    except TypeError:
        pass
_owl_ns.Ontology._load_properties = _patched_load_properties



def merge_repeated_fileobservations(df):


    def process_group(group):

        # Handle discrepancies for columns not in 'li_make_most_common'
        li_rm_discrepancies = [col for col in group.columns]
        for col in li_rm_discrepancies:
            if len(set(group[col])) > 1 or all(value is None for value in group[col]):
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
    

def get_taxonomic_name_from_id(ncbi_id):
    url = f"https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id={ncbi_id}"
    attempts = 0
    max_attempts = 3

    while attempts < max_attempts:
        try:
            response = requests.get(url)
            response.raise_for_status()  # Check for HTTP errors.
            soup = BeautifulSoup(response.text, 'html.parser')

            # Attempt to find the taxonomic name more reliably.
            title_content = soup.title.string
            # Adjusting strategy to account for potential differences in page structure.
            if "Taxonomy browser" in title_content:
                taxonomic_name = title_content.split("(")[-1].split(")")[0]
            else:
                # Fallback if the expected pattern is not found
                header = soup.find('h2')
                if header and "Taxonomy browser" in header.text:
                    taxonomic_name = header.text.split("(")[-1].split(")")[0]
                else:
                    raise ValueError("Taxonomic name pattern not recognized.")

            if taxonomic_name:  # Check if a name was found
                return taxonomic_name
            else:
                raise ValueError("Taxonomic name not found.")
        except Exception as e:
            print(f"Attempt {attempts + 1}: An error occurred - {e}")
            time.sleep(6)  # Wait before retrying
            attempts += 1

    return None
                      

def get_taxonomy_info(ncbi_id, cell_culture_key1 = '', cell_culture_key2 = ''):
    if ncbi_id is not None and ncbi_id != "NA":
        cell_culture_key_words = ["cell", "media", "culture"]
        try:
            #try to get taxa via API
            for attempt in range(3):
                try:
                    response = requests.get(
                        f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=taxonomy&id={ncbi_id}")
                    soup = BeautifulSoup(response.text, "xml")
                    classification = [s.lower() for s in soup.find("Taxon").find("Lineage").text.split("; ")]
                    break

                except requests.exceptions.RequestException as e:
                    print(f"Taxonomy request failed. Retrying (attempt {attempt + 1}/3): {e}")
                    time.sleep(10)

            else:
                # If the loop completes without breaking, raise an exception
                raise Exception("All retries failed. Unable to fetch taxonomy data.")

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
        except Exception:
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
    # named clinical / nutrition dietary patterns (MIND before DASH so the
    # MIND-DASH hybrid resolves to MIND, its standard name)
    if re.search(r'low fodmap|\bfodmap\b', x):
        return 'low FODMAP diet'
    if re.search(r'gluten free|\bgluten\b free', x):
        return 'gluten-free diet'
    if re.search(r'\bmind\b', x):
        return 'MIND diet'
    if re.search(r'\bdash\b', x):
        return 'DASH diet'
    if re.search(r'\bbrat\b', x):
        return 'BRAT diet'
    if re.search(r'ultra ?processed|\bupf\b', x):
        return 'ultra-processed diet'
    if re.search(r'minimally processed|minimal(ly)? process', x):
        return 'minimally-processed diet'
    # infant feeding regimens (breast milk / formula / solids); check before the
    # experimental diets so "formula" is not misread as a composition diet
    _bm = re.search(r'breast ?milk|breast ?fed|breast ?feed', x)
    _formula = re.search(r'formula', x)
    if re.search(r'\bno breast|without breast|never breast', x):
        return 'exclusively formula'
    if re.search(r'solid', x) and _bm and _formula:
        return 'breastmilk and formula with solids'
    if re.search(r'solid', x) and _bm:
        return 'breastmilk with solids'
    if (_bm and _formula) or re.search(r'\bmixed (fed|feed|feeding|milk)\b', x):
        return 'mixed breastmilk and formula'
    if re.search(r'not\s+exclusive', x) and _bm:
        return 'not exclusively breastfed'
    if re.search(r'exclusive', x) and _bm:
        return 'exclusively breastmilk'
    if _bm:
        return 'breastfed (NOS)'
    if _formula:
        return 'exclusively formula'
    # deficiency / supplementation / fiber diets (common in rodent studies)
    if re.search(r'choline def|methionine choline|\bmcd\b|\bcdaa\b|\bcdahfd\b', x):
        return 'choline-deficient diet'
    if re.search(r'iron def|low iron|iron deficient', x):
        return 'iron-deficient diet'
    if re.search(r'iron suppl|iron fortif|iron enrich|high iron', x):
        return 'iron-supplemented diet'
    if re.search(r'high fib(er|re)', x):
        return 'high-fiber diet'
    if re.search(r'low fib(er|re)', x):
        return 'low-fiber diet'
    # feeding schedule (check time-restricted before calorie-restricted)
    if re.search(r'time restricted|\btrf\b|\btre\b', x):
        return 'time-restricted diet'
    if re.search(r'intermittent fasting|alternate day fasting|\badf\b|\bfasting\b', x):
        return 'intermittent fasting'
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
    if re.search(r'calorie restrict|caloric restrict|energy restrict', x):
        return 'calorie-restricted diet'
    # purified/defined vs standard chow vs generic control
    if re.search(r'purified|\bdefined\b|ain ?76|ain ?93|ain76|ain93', x):
        return 'purified diet'
    if re.search(r'\bchow\b|regular chow|standard chow|\bncd\b', x):
        return 'chow diet'
    if re.search(r'\bcontrol\b|\bnormal\b|\bstandard\b|\bcd\b|\bcon\b|\bctrl\b|\bsc\b|normal diet', x):
        return 'control diet (NOS)'
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
            response = requests.get(
                f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=taxonomy&term={species_name}&retmode=xml"
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