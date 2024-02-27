import pandas as pd    
import csv
import json

if __name__ == "__main__":

    csv_file_path = '/home/yasin/yasin/projects/PublicDataset_ReDU_Metadata_Workflow/allowed_terms/allowed_terms.csv'

    json_dict = {}

    with open(csv_file_path, mode='r', encoding='utf-8') as csvfile:

        reader = csv.DictReader(csvfile)

        for row in reader:
            for header, value in row.items():
                if header not in json_dict:
                    json_dict[header] = {"missing": "ML import: not available", "allowed_values": []}
                if value not in json_dict[header]["allowed_values"] and value != '':
                    json_dict[header]["allowed_values"].append(value)

    json_file_path = '/home/yasin/yasin/projects/PublicDataset_ReDU_Metadata_Workflow/allowed_terms/allowed_terms.json'

    with open(json_file_path, 'w', encoding='utf-8') as jsonfile:

        json.dump(json_dict, jsonfile, ensure_ascii=False, indent=4)
