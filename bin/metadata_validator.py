#!/usr/bin/python

import os
import argparse
import csv
import pandas as pd
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

def perform_validation(filename, path_to_allowed_term_csv):
    
    terms_df = pd.read_csv(path_to_allowed_term_csv, low_memory=False)

    validators = {
        'filename': [
            NotEmptyValidator()
        ],
        "MassiveID" : [
            Ignore()
        ],
        'SampleType' : [
            SetValidator(valid_set= terms_df['SampleType'].dropna().loc[terms_df['SampleType'] != ""].unique().tolist())
        ],
        'SampleTypeSub1' : [
            SetValidator(valid_set= terms_df['SampleTypeSub1'].dropna().loc[terms_df['SampleTypeSub1'] != ""].unique().tolist())
        ],
        'NCBITaxonomy' : [
            SetValidator(valid_set= terms_df['NCBITaxonomy'].dropna().loc[terms_df['NCBITaxonomy'] != ""].unique().tolist())
        ],
        'YearOfAnalysis' : [
            IntValidator(),
            RangeValidator(low=1999, high=2030)
        ],
        'SampleCollectionMethod' : [
            SetValidator(valid_set= terms_df['SampleCollectionMethod'].dropna().loc[terms_df['SampleCollectionMethod'] != ""].unique().tolist())
        ],
        "SampleExtractionMethod" : [
            SetValidator(valid_set= terms_df['SampleExtractionMethod'].dropna().loc[terms_df['SampleExtractionMethod'] != ""].unique().tolist())
        ],
        'InternalStandardsUsed' : [
            SetValidator(valid_set= terms_df['InternalStandardsUsed'].dropna().loc[terms_df['InternalStandardsUsed'] != ""].unique().tolist())
        ],
        'MassSpectrometer' : [
            SetValidator(valid_set= terms_df['MassSpectrometer'].dropna().loc[terms_df['MassSpectrometer'] != ""].unique().tolist())
        ],        
        'IonizationSourceAndPolarity' : [
            SetValidator(valid_set= terms_df['IonizationSourceAndPolarity'].dropna().loc[terms_df['IonizationSourceAndPolarity'] != ""].unique().tolist())
        ],
        'ChromatographyAndPhase' : [
            SetValidator(valid_set= terms_df['ChromatographyAndPhase'].dropna().loc[terms_df['ChromatographyAndPhase'] != ""].unique().tolist())
        ],
        'BiologicalSex': [
            SetValidator(valid_set= terms_df['BiologicalSex'].dropna().loc[terms_df['BiologicalSex'] != ""].unique().tolist())
        ],
        'Country': [
            SetValidator(valid_set= terms_df['Country'].dropna().loc[terms_df['Country'] != ""].unique().tolist())
        ],
        'HumanPopulationDensity' : [
            SetValidator(valid_set= terms_df['HumanPopulationDensity'].dropna().loc[terms_df['HumanPopulationDensity'] != ""].unique().tolist())
        ],
         'LifeStage' : [
            SetValidator(valid_set= terms_df['LifeStage'].dropna().loc[terms_df['LifeStage'] != ""].unique().tolist())
        ],
        'UBERONOntologyIndex' : [
            SetValidator(valid_set= terms_df['UBERONOntologyIndex'].dropna().loc[terms_df['UBERONOntologyIndex'] != ""].unique().tolist())
        ],
        'DOIDOntologyIndex' : [
            SetValidator(valid_set= terms_df['DOIDOntologyIndex'].dropna().loc[terms_df['DOIDOntologyIndex'] != ""].unique().tolist())
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
    args = parser.parse_args()

    path_to_allowed_term_csv = '/home/yasin/yasin/projects/PublicDataset_ReDU_Metadata_Workflow/allowed_terms/allowed_terms.csv'

    passes_validation, failures, errors_list, valid_rows, total_rows = perform_validation(args.inputmetadata, path_to_allowed_term_csv)
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
