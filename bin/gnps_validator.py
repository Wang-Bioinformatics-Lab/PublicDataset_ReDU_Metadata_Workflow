import os 
import pandas as pd
import urllib
import requests
import io 
import sys
import argparse
import subprocess
import collections
from subprocess import PIPE, run

def main():
    # parsing arguments
    parser = argparse.ArgumentParser(description='GNPS Validator')
    parser.add_argument('file_paths', help='File Paths')
    parser.add_argument('metadata_folder')
    args = parser.parse_args()

    # Print message to indicate importing is done
    print("echo Importing Done ...")

    # Read the TSV file containing file paths into a pandas DataFrame
    df = pd.read_csv(args.file_paths, delimiter='\t')

    # Get a list of all the values in the first column (file names)
    file_names = df.iloc[:, 0].tolist()

    # List to store names of files that have passed validation
    passed_file_names = []

    # Loop through each file in the list of file names and validate each one
    print("echo Validating Files now ...")
    for file_name in file_names:
        # Call the metadata_validator.py script and pass the file name as an argument
        import metadata_validator

        try:
            passes_validation, failures, errors_list, valid_rows, total_rows = metadata_validator.perform_validation(os.path.join(args.metadata_folder, os.path.basename(file_name)))
        except:
            pass

        if passes_validation:
            passed_file_names.append(file_name)
            

    # Print message to indicate that validation is complete
    print("echo Validation Completed !")

    # Convert the list of passed file names to a pandas DataFrame
    passed = pd.DataFrame(passed_file_names)
    
    # Print message to indicate that a TSV file is being created for path names
    print("echo Now creating tsv file for Path names ...")

    # Write the passed file names DataFrame to a TSV file
    passed.to_csv('passed_file_names.tsv', sep='\t', index=False, header=False)

    # Print message to indicate that the TSV file has been created and the script is complete
    print("echo TSV File created! Have a good day.")


if __name__ == "__main__":
    main()