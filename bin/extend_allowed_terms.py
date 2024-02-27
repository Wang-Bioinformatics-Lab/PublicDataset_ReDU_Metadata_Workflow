import json
import os
import pandas as pd

def adapt_allowed_terms(terms_dict, redu_variable, term_list, add_or_remove, add_to_terms = '', load_dict_from_path = 'no', save_dict_to_path = 'no'):

    # ['MassiveID', 'filename', 'SampleType', 'SampleTypeSub1', 'NCBITaxonomy', 'YearOfAnalysis', 'SampleCollectionMethod', 'SampleExtractionMethod', 'InternalStandardsUsed', 'MassSpectrometer', 'IonizationSourceAndPolarity', 
    # 'ChromatographyAndPhase', 'SubjectIdentifierAsRecorded', 'AgeInYears', 'BiologicalSex', 'UBERONBodyPartName', 'TermsofPosition', 'HealthStatus', 'DOIDCommonName', 'ComorbidityListDOIDIndex', 'SampleCollectionDateandTime', 
    # 'Country', 'HumanPopulationDensity', 'LatitudeandLongitude', 'DepthorAltitudeMeters', 'sample_name', 'UniqueSubjectID', 'LifeStage', 'UBERONOntologyIndex', 'DOIDOntologyIndex']

    if load_dict_from_path != 'no' and os.path.exists(load_dict_from_path):
        with open(load_dict_from_path, 'r', encoding='utf-8') as jsonfile:
            terms_dict = json.load(jsonfile)

    if len(term_list) > 0:
        pass
    else:
        raise ValueError("No values to add/remove provided.")

    #add new terms
    if add_or_remove == 'add':

        if redu_variable not in terms_dict:
            print(f"Warning: Key '{redu_variable}' not found in the dictionary.")
            return

        for value in term_list:
            if value not in terms_dict[redu_variable]["allowed_values"]:
                value_to_add = value + add_to_terms
                if value_to_add not in terms_dict[redu_variable]["allowed_values"]:
                    terms_dict[redu_variable]["allowed_values"].append(value_to_add)
            else:
                print(f"Warning: '{value}' already exists in '{redu_variable}'.")
    

    #remove terms
    if add_or_remove == 'remove':
        if redu_variable not in terms_dict:
            print(f"Warning: Key '{redu_variable}' not found in the dictionary.")
            return

        for value in term_list:
            if value in terms_dict[redu_variable]["allowed_values"]:
                terms_dict[redu_variable]["allowed_values"].remove(value)
            else:
                print(f"Warning: '{value}' not found in '{redu_variable}'.")

    if save_dict_to_path != 'no':
        with open(save_dict_to_path, 'w', encoding='utf-8') as jsonfile:
            json.dump(terms_dict, jsonfile, ensure_ascii=False, indent=4)


    return terms_dict


def add_new_variable(terms_dict, redu_variable):
    """
    Adds a new key to the dictionary. If the key already exists, raises an error.
    The new key will have a nested dictionary with 'missing' and 'allowed_values'.
    'missing' has a fixed value, and 'allowed_values' is an empty list.
    """
    if redu_variable in terms_dict:
        raise ValueError(f"Error: The key '{redu_variable}' already exists in the dictionary.")
    
    terms_dict[redu_variable] = {
        "missing": "ML import: not available",
        "allowed_values": []
    }



def update_dict_fromautoupdate(data):
    """
    Update the given dictionary by removing '__AUTOUPDATE' from allowed values and eliminating duplicates.

    :param data: dict, the original dictionary to be updated
    """
    for key in data:
        updated_values = set()
        for value in data[key]["allowed_values"]:
            updated_values.add(value.replace("__AUTOUPDATE", ""))
        data[key]["allowed_values"] = list(updated_values)

    return data



def update_translation_sheets(terms_dict, path_to_translation_sheet_folder):

    special_keys = ["MassSpectrometer"]

    REDU_vars = list(terms_dict.keys())

    files = os.listdir(path_to_translation_sheet_folder)
    csv_files = [os.path.join(path_to_translation_sheet_folder, file) for file in files if file.endswith('.csv')]


    for csv_file in csv_files:
    # Extract the key name from the file name
        key_name = os.path.splitext(os.path.basename(csv_file))[0]
        
        # Check if this key has a corresponding entry in the dictionary
        if key_name in REDU_vars:
            # Read the CSV file
            try:
                df = pd.read_csv(csv_file, encoding='utf-8')
            except UnicodeDecodeError:
                df = pd.read_csv(csv_file, encoding='latin1')
            
            # Get the list of allowed terms for this key
            allowed_terms = terms_dict[key_name]['allowed_values']
            
            # Check each allowed term
            for term in allowed_terms:
                # If the term is not present in the "REDU" column, add it
                if term not in df['REDU'].values:
                    if '|' in term and key_name in special_keys:
                        term_before_pipe, _ = term.split('|', 1)  # Only split on the first occurrence
                        new_row = {'MWB': term_before_pipe, 'REDU': term}
                    else:
                        new_row = {'MWB': term, 'REDU': term}
                    new_row_df = pd.DataFrame([new_row])
                    df = pd.concat([df, new_row_df], ignore_index=True)

            df = df.drop_duplicates(subset=['MWB'], keep='first')
            updated_file_name = os.path.splitext(csv_file)[0] + '.csv'
            df.to_csv(updated_file_name, index=False)


    

if __name__ == '__main__':

    path_to_allowed_terms = '/home/yasin/projects/ReDU-MS2-GNPS2/workflows/PublicDataset_ReDU_Metadata_Workflow/bin/allowed_terms/allowed_terms.json'
    with open(path_to_allowed_terms, 'r', encoding='utf-8') as jsonfile:
        terms_dict = json.load(jsonfile)

    output = update_translation_sheets(terms_dict,
                                       '/home/yasin/projects/ReDU-MS2-GNPS2/workflows/PublicDataset_ReDU_Metadata_Workflow/bin/translation_sheets')
    
    print(output)
    #pass


    # path_to_allowed_terms = '/home/yasin/projects/ReDU-MS2-GNPS2/workflows/PublicDataset_ReDU_Metadata_Workflow/bin/allowed_terms/allowed_terms_autoupdate.json'

    # with open(path_to_allowed_terms, 'r', encoding='utf-8') as jsonfile:
    #     terms_dict = json.load(jsonfile)


    # new_redu_variable_to_add = ''
    # redu_variable_for_which_terms_should_be_changed = 'NCBITaxonomy'
    # list_of_terms_to_add_or_remove = []

    # #to add a new variable to REDU
    # terms_dict = add_new_variable(terms_dict, new_redu_variable_to_add)

    # #to add terms
    # terms_dict = adapt_allowed_terms(terms_dict, 
    #                                  redu_variable_for_which_terms_should_be_changed, 
    #                                  list_of_terms_to_add_or_remove, 
    #                                  'add')

    # #to remove_terms
    # terms_dict = adapt_allowed_terms(terms_dict, 
    #                                  redu_variable_for_which_terms_should_be_changed, 
    #                                  list_of_terms_to_add_or_remove, 
    #                                  'remove')

    # terms_dict = update_dict_fromautoupdate(terms_dict)

    # with open(path_to_allowed_terms, 'w', encoding='utf-8') as jsonfile:
    #     json.dump(terms_dict, jsonfile, ensure_ascii=False, indent=4)

