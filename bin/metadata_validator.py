#!/usr/bin/python

import os
import argparse
import csv
import pandas as pd
import json
from vladiate import Vlad
from vladiate.validators import UniqueValidator, SetValidator, Ignore, IntValidator, RangeValidator, NotEmptyValidator
from vladiate.inputs import LocalFile

def rewrite_metadata(metadata_filename):
    """
    Metadata Fields Rewrite
    Fields changed and will need to be rewritten
    """

    metadata_df = pd.read_csv(metadata_filename, sep="\t")

    #Limiting the number of rows
    metadata_df = metadata_df.truncate(after=10000)

    #Rewriting Year of Analysis
    metadata_list = metadata_df.to_dict(orient="records")
    for metadata_obj in metadata_list:
        try:
            metadata_obj["YearOfAnalysis"] = str(int(float(metadata_obj["YearOfAnalysis"])))
        except:
            continue

    metadata_df = pd.DataFrame(metadata_list)
    metadata_df.to_csv(metadata_filename, sep="\t", index=False)


def perform_validation(filename, path_to_allowed_terms):

    with open(path_to_allowed_terms, 'r', encoding='utf-8') as jsonfile:
        terms = json.load(jsonfile)
        
    
    validators = {
        'filename': [
            NotEmptyValidator()
        ],
        "MassiveID" : [
            Ignore()
        ],
        'SampleType' : [
            SetValidator(valid_set= list(set(terms['SampleType']["allowed_values"])))
        ],
        'SampleTypeSub1' : [
            SetValidator(valid_set= list(set(terms['SampleTypeSub1']["allowed_values"])))
        ],
        'NCBITaxonomy' : [
            SetValidator(valid_set= list(set(terms['NCBITaxonomy']["allowed_values"])))
        ],
        'YearOfAnalysis' : [
            IntValidator(),
            RangeValidator(low=1999, high=2030)
        ],
        'SampleCollectionMethod' : [
            SetValidator(valid_set= list(set(terms['SampleCollectionMethod']["allowed_values"]))) 
        ],
        "SampleExtractionMethod" : [
            SetValidator(valid_set= list(set(terms['SampleExtractionMethod']["allowed_values"]))) 
        ],
        'InternalStandardsUsed' : [
            SetValidator(valid_set= list(set(terms['InternalStandardsUsed']["allowed_values"]))) 
        ],
        'MassSpectrometer' : [
            SetValidator(valid_set= list(set(terms['MassSpectrometer']["allowed_values"]))) 
        ],        
        'IonizationSourceAndPolarity' : [
            SetValidator(valid_set= list(set(terms['IonizationSourceAndPolarity']["allowed_values"])))
        ],
        'ChromatographyAndPhase' : [
            SetValidator(valid_set= list(set(terms['ChromatographyAndPhase']["allowed_values"]))) 
        ],
        'BiologicalSex': [
            SetValidator(valid_set= list(set(terms['BiologicalSex']["allowed_values"]))) 
        ],
        'Country': [
            SetValidator(valid_set= list(set(terms['Country']["allowed_values"])))  
        ],
        'HumanPopulationDensity' : [
            SetValidator(valid_set= list(set(terms['HumanPopulationDensity']["allowed_values"])))
        ],
         'LifeStage' : [
            SetValidator(valid_set= list(set(terms['LifeStage']["allowed_values"]))) 
        ],
        'UBERONOntologyIndex' : [
            SetValidator(valid_set= list(set(terms['UBERONOntologyIndex']["allowed_values"]))) 
        ],
        'DOIDOntologyIndex' : [
            SetValidator(valid_set= list(set(terms['DOIDOntologyIndex']["allowed_values"])))  
        ]
    }

    my_validator = Vlad(source=LocalFile(filename),delimiter="\t",ignore_missing_validators=True,validators=validators)
    passes_validation = my_validator.validate()
    
    if passes_validation:
        print ("Success")
    else:
        print("Fail")

    errors_list = []
    for column in my_validator.failures:
        for line_number in my_validator.failures[column]:
            error_dict = {}
            error_dict["header"] = column
            error_dict["line_number"] = line_number + 1 #0 Indexed with 0 being the header row
            error_dict["error_string"] = str(my_validator.failures[column][line_number])
            errors_list.append(error_dict)

            if len(errors_list) > 500:
                break

        if len(errors_list) > 500:
            break

    for missing_field in my_validator.missing_fields:
        error_dict = {}
        error_dict["header"] = "Missing Header"
        error_dict["line_number"] = "N/A"
        error_dict["error_string"] = "Missing column %s" % (missing_field)

        errors_list.append(error_dict)

    # Reading in the validation
    metadata_df = pd.read_csv(filename, sep="\t")
    row_count = len(metadata_df)

    valid_rows = []
    #Read in the good rows
    try:
        no_validation_lines = [int(error["line_number"]) for error in errors_list]
        with open(filename) as csvfile:
            reader = csv.DictReader(csvfile, delimiter="\t")
            for row in reader:
                if row_count in no_validation_lines:
                    continue
                valid_rows.append(row)
    except:
        #raise
        print("error reading file")

    return passes_validation, my_validator.failures, errors_list, valid_rows, row_count

def perform_summary(filename):
    with open(filename, 'r') as csvfile:
        reader = csv.DictReader(csvfile, delimiter="\t")

        summary_dict = {}
        summary_dict["row_count"] = sum([1 for row in reader])

        summary_list = []
        summary_list.append({"type" : "row_count", "value" : summary_dict["row_count"]})

        return summary_dict, summary_list

def main():
    parser = argparse.ArgumentParser(description='Validate Stuff.')
    parser.add_argument('inputmetadata', help='inputmetadata')
    parser.add_argument('path_to_allowed_terms_json', help='Path to json with allowed terms.')
    args = parser.parse_args()


    passes_validation, failures, errors_list, valid_rows, total_rows = perform_validation(args.inputmetadata, args.path_to_allowed_terms_json)
    no_validation_lines = [int(error["line_number"]) for error in errors_list]

    output_list = ["MING", os.path.basename(args.inputmetadata), str(total_rows), str(len(valid_rows))]
    print("\t".join(output_list))


    #with open(args.inputmetadata, 'rb') as csvfile:
        #dialect = csv.Sniffer().sniff(csvfile.read(1024))
        #csvfile.seek(0)
        #reader = csv.DictReader(csvfile, dialect=dialect)
        #for row in reader:
        #    print(row)

if __name__ == "__main__":
    main()
