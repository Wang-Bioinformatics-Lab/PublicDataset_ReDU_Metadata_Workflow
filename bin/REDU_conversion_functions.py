import json
import os
import requests
from extend_allowed_terms import adapt_allowed_terms 
import time
from bs4 import BeautifulSoup
from owlready2 import get_ontology
import owlready2
import pandas as pd
import tqdm



def complete_and_fill_REDU_table(df, allowedTerm_dict):
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
    # Convert year to string for comparison
    if 'YearOfAnalysis' in df.columns:
        df['YearOfAnalysis'] = df['YearOfAnalysis'].astype(str)

    # Add missing columns with their respective default missing value from the dictionary
    for key, value in allowedTerm_dict.items():
        if key not in df.columns:
            df[key] = value['missing']

    # Replace values with the respective "missing" value if they're not in the allowed terms or are missing/empty
    for key, value in allowedTerm_dict.items():
        allowed_terms = value['allowed_values']
        missing_value = value['missing']
        if key in df.columns:
            if key not in ["MassiveID", "filename", "AgeInYears"]:
                df[key] = df[key].apply(lambda x: x if x in allowed_terms else missing_value).fillna(missing_value).replace("", missing_value)
            else:
                df[key] = df[key].fillna(missing_value).replace("", missing_value)

    # Ensure the dataframe contains only the columns specified in the dictionary plus 'USI'
    return df[[*allowedTerm_dict.keys()]]


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

def update_unassigned_terms(organism_name, column_key = "Samples_Organism", autoupdated=False, unassigned_file='unassigned_terms.json'):
    if os.path.exists(unassigned_file):
        with open(unassigned_file, 'r') as file:
            unassigned_data = json.load(file)
    else:
        unassigned_data = {column_key: {}}
    
    if organism_name in unassigned_data[column_key]:
        unassigned_data[column_key][organism_name]["count"] += 1
    else:
        unassigned_data[column_key][organism_name] = {"count": 1, "autoupdated": autoupdated}

    with open(unassigned_file, 'w') as file:
        json.dump(unassigned_data, file, indent=4)


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






def get_raw_data_from_factors_json():

    #ST001491
    stdy_info_req = requests.get(
        'https://www.metabolomicsworkbench.org/rest/study/study_id/{}/factors'.format(str('ST001491')))
    
    stdy_info = stdy_info_req.json()

    print(stdy_info)

if __name__ == '__main__':
    get_raw_data_from_factors_json()