import os
import pandas as pd
import argparse
import collections
import requests
import glob
from io import StringIO
from subprocess import PIPE, run

ccms_peak_link = "https://gnps-datasetcache.ucsd.edu/datasette/database/filename.csv?_sort=filepath&collection__exact=ccms_peak&_size=max"
gnps_column_names = ['filename', 'ATTRIBUTE_DatasetAccession', 'AgeInYears', 'BiologicalSex', 'ChromatographyAndPhase', 'ComorbidityListDOIDIndex', 'Country', 'DOIDCommonName', 'DOIDOntologyIndex', 'DepthorAltitudeMeters', 'HealthStatus', 'HumanPopulationDensity', 'InternalStandardsUsed', 'IonizationSourceAndPolarity', 'LatitudeandLongitude', 'LifeStage', 'MassSpectrometer', 'NCBITaxonomy', 'SampleCollectionDateandTime', 'SampleCollectionMethod', 'SampleExtractionMethod', 'SampleType', 'SampleTypeSub1', 'SubjectIdentifierAsRecorded', 'TermsofPosition', 'UBERONBodyPartName', 'UBERONOntologyIndex', 'UniqueSubjectID', 'YearOfAnalysis']
gnps_column_names_added = ['USI']

def _make_usi_from_filename(s):
    # Replace initial "f." with "mzspec:"
    if s.startswith("f."):
        s = "mzspec:" + s[2:]
    # Replace the first "/" with ":"
    first_slash_index = s.find('/')
    if first_slash_index != -1:
        s = s[:first_slash_index] + ':' + s[first_slash_index+1:]
    return s




def _match_filenames_and_add_usi(dataset_metadata_df):
    try:            
        # Get the value of the first row of the 'column_name' column
        dataset = dataset_metadata_df['ATTRIBUTE_DatasetAccession'].iloc[0]

        dataset_files_response = requests.get("{}&dataset__exact={}".format(ccms_peak_link, dataset))
        csvStringIO = StringIO(dataset_files_response.text)
        ccms_df = pd.read_csv(csvStringIO)

        ccms_df["query_path"] = ccms_df["filepath"].apply(lambda x: os.path.basename(x))
    except TypeError:
        print("Error")
        return None

    
    metadata_row_list = dataset_metadata_df.to_dict('records')
    output_row_list = []

    for metadata_row in metadata_row_list:
        filename = os.path.basename(metadata_row["filename"])
        filename2 = filename[:-3] + "ML"

        found_file_paths = []

        # Searching the query_path column
        found_file_paths = ccms_df[ccms_df["query_path"] == filename]["filepath"].tolist()
        found_file_paths += ccms_df[ccms_df["query_path"] == filename2]["filepath"].tolist()

        if len(found_file_paths) == 1:
            # We've found a match
            print("Found match", found_file_paths[0])
            metadata_row["filename"] = "f." + found_file_paths[0]
            metadata_row["USI"] = _make_usi_from_filename(metadata_row["filename"])

            output_row_list.append(metadata_row)
        else:
            # Didn't find or is ambiguous
            continue

    return pd.DataFrame(output_row_list)

def main():
    # parsing args
    parser = argparse.ArgumentParser(description='GNPS Name Matcher')
    parser.add_argument('passed_file_names', help='Input TSV file')
    parser.add_argument('metadata_folder')

    args = parser.parse_args()
     
    ccms_filenames = collections.defaultdict(set)
    
    if args.passed_file_names != 'all':
        # Read the TSV file and specify the delimiter as a tab
        df = pd.read_csv(args.passed_file_names, delimiter='\t', header=None, names=['Name'])
        # Extract the names from a specific column (e.g., column 'Name')
        passed_file_names = df['Name'].tolist()
    else:
        passed_file_names = glob.glob(f"{args.metadata_folder}/*.tsv")

    print("echo Iterating though rows now")
    
    all_metadata_list = []

    for file in passed_file_names:
        # print("Length of visit is ",len(visit))
        print("Working on ", file)
        # print("echo Working on")
        # csv_path = os.path.join(current_dir, './data.csv' file)
        dataset_metadata_df = pd.read_csv( file , delimiter='\t')
        
        #Renaming the coloumn, Matching common columns and rearranging them in same order to final file
        dataset_metadata_df = dataset_metadata_df.rename(columns={'MassiveID': 'ATTRIBUTE_DatasetAccession'})
        common_cols = list(set(gnps_column_names).intersection(set(dataset_metadata_df.columns)))
        dataset_metadata_df = dataset_metadata_df.loc[:, common_cols]
        try:
            dataset_metadata_df = dataset_metadata_df[gnps_column_names]
        except KeyError:
            print(f"Skipping file {file} due to a TypeError.")
            continue

        # Matching the metadata
        enriched_metadata_df = _match_filenames_and_add_usi(dataset_metadata_df)
        if enriched_metadata_df is not None:
            print('returning at least some files')
            all_metadata_list.append(enriched_metadata_df)

    # Create a DataFrame from the list with headers
    merged_metadata_df = pd.concat(all_metadata_list)
    merged_metadata_df = merged_metadata_df[gnps_column_names + gnps_column_names_added]

    # Save the DataFrame to a TSV file without column names
    merged_metadata_df.to_csv('gnps_metadata_all.tsv', sep='\t', index=False)


if __name__ == '__main__':
    main()