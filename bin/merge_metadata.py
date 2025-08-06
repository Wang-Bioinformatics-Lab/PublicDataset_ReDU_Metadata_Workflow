import pandas as pd
import argparse
import json

def main():
    # parsing arguments
    parser = argparse.ArgumentParser(description='Merge GNPS and ReDU metadata')
    parser.add_argument('--gnps_metadata')
    parser.add_argument('--mwb_metadata')
    parser.add_argument('--metabolights_metadata')
    parser.add_argument('--norman_metadata')
    parser.add_argument('--masst_metadata')
    parser.add_argument('--path_to_allowed_term_json')
    parser.add_argument('--output_metadata')
    args = parser.parse_args()

    # columns_to_use = ["filename", "ATTRIBUTE_DatasetAccession", "SampleType", "SampleTypeSub1", 
    #                   "NCBITaxonomy", "NCBIDivision", "NCBIRank", "YearOfAnalysis", "UBERONBodyPartName", "BiologicalSex", "AgeInYears",  "LifeStage", 
    #                   "Country", "HealthStatus", "ChromatographyAndPhase", "IonizationSourceAndPolarity",
    #                   "MassSpectrometer", "SampleExtractionMethod",  "SampleCollectionMethod", "ComorbidityListDOIDIndex", 
    #                   "DOIDCommonName", "DOIDOntologyIndex", "ENVOEnvironmentBiome", "DepthorAltitudeMeters", "HumanPopulationDensity", "InternalStandardsUsed", 
    #                   "LatitudeandLongitude", "SampleCollectionDateandTime", "ENVOEnvironmentMaterial", "ENVOEnvironmentBiomeIndex", "ENVOEnvironmentMaterialIndex",
    #                   "SubjectIdentifierAsRecorded", "TermsofPosition", "UBERONOntologyIndex", "UniqueSubjectID", "USI", "DataSource"]
    

    # read allowed terms
    with open(args.path_to_allowed_term_json, 'r') as f:
        allowed_terms = json.load(f)

    columns_to_use = list(allowed_terms.keys())

    # replace value "MassiveID" with "ATTRIBUTE_DatasetAccession"
    columns_to_use = [col.replace("MassiveID", "ATTRIBUTE_DatasetAccession") for col in columns_to_use]

    # read GNPS metadata
    gnps_df = pd.read_csv(args.gnps_metadata, sep='\t')
    gnps_df["DataSource"] = "GNPS"

    #drop duplicated files
    duplicates = gnps_df.duplicated(subset='USI', keep=False)
    gnps_df = gnps_df[~duplicates]

    # read Workbench metadata
    mwb_df = pd.read_csv(args.mwb_metadata, sep='\t', dtype=str)
    mwb_df["DataSource"] = "Workbench"

    # read MetaboLights metadata
    metabo_df = pd.read_csv(args.metabolights_metadata, sep='\t', dtype=str)
    metabo_df["DataSource"] = "MetaboLights"

    # read NORMAN metadata
    norman_df = pd.read_csv(args.norman_metadata, sep='\t', dtype=str)
    norman_df["DataSource"] = "NORMAN"

    # read MASST metadata
    masst_df = pd.read_csv(args.masst_metadata, sep='\t', dtype=str)
    masst_df["DataSource"] = "GNPS"

    # merge GNPS and ReDU metadata
    merged_df = pd.concat([gnps_df, mwb_df, metabo_df, masst_df, norman_df], ignore_index=True)

    # include only columns from gnps
    merged_df = merged_df[columns_to_use]

    # drop eventual duplicates (could happen with MASST/GNPS and also for a different reason with MWB) MWB has duplicates across all rows GNPS/MASST we should keep the first entry (from GNPS)
    merged_df = merged_df.drop_duplicates(subset='USI')

    # write merged metadata to file
    merged_df.to_csv(args.output_metadata, sep='\t', index=False)

if __name__ == '__main__':
    main()