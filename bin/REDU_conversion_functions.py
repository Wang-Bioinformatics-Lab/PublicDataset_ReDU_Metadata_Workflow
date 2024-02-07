import json
import os
import requests
from extend_allowed_terms import adapt_allowed_terms 
import time
from bs4 import BeautifulSoup

def get_taxonomy_id_from_name__allowedTerms(organism_name, **kwargs):

    allowedTerm_dict = kwargs['allowedTerm_dict']
    taxonomy_data = allowedTerm_dict["NCBITaxonomy"]["allowed_values"]
    for entry in taxonomy_data:
        parts = entry.split('|')
        if len(parts) == 2:
            ncbi_id, name = parts
            if name.lower() == organism_name.lower():
                return ncbi_id + '|' + str(name)
        else:
            continue

    req_ncbi_name = get_taxonomy_id_from_name(organism_name)
    autoupdate = False


    if req_ncbi_name is not None:

        print([req_ncbi_name + '__AUTOUPDATE'])

        allowedTerm_dict = adapt_allowed_terms(terms_dict = allowedTerm_dict, 
                                               redu_variable = 'NCBITaxonomy', 
                                               term_list = [req_ncbi_name + '__AUTOUPDATE'], 
                                               add_or_remove = 'add', 
                                               load_dict_from_path = 'no', 
                                               save_dict_to_path = '/home/yasin/projects/ReDU-MS2-GNPS2/workflows/PublicDataset_ReDU_Metadata_Workflow/bin/allowed_terms/allowed_terms_autoupdate.json')

        autoupdate = True

    update_unassigned_terms(organism_name, autoupdated=autoupdate)
    return None

def update_unassigned_terms(organism_name, autoupdated=False, unassigned_file='unassigned_terms.json'):
    if os.path.exists(unassigned_file):
        with open(unassigned_file, 'r') as file:
            unassigned_data = json.load(file)
    else:
        unassigned_data = {"Samples_Organism": {}}
    
    if organism_name in unassigned_data["Samples_Organism"]:
        unassigned_data["Samples_Organism"][organism_name]["count"] += 1
    else:
        unassigned_data["Samples_Organism"][organism_name] = {"count": 1, "autoupdated": autoupdated}

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