import os
import pandas as pd
import argparse
import collections
from subprocess import PIPE, run

def main():
    # parsing args
    parser = argparse.ArgumentParser(description='GNPS Name Matcher')
    parser.add_argument('passed_file_names', help='Input TSV file')
    parser.add_argument('metadata_folder')

    args = parser.parse_args()
     
    ccms_peak_link_start = "https://gnps-datasetcache.ucsd.edu/datasette/database/filename.csv?_sort=filepath&collection__exact=ccms_peak&dataset__exact="
    ccms_peak_link_end = "&_size=max"

    ccms_filenames = collections.defaultdict(set)
    
    print("echo CCMS Read Done!")

    
    
    # Read the TSV file and specify the delimiter as a tab
    df = pd.read_csv(args.passed_file_names, delimiter='\t', header=None, names=['Name'])
    # Extract the names from a specific column (e.g., column 'Name')
    passed_file_names = df['Name'].tolist()

    print("echo Iterating though rows now")
    merged_rows = {}
    # unique = 0
    # dupli = 0

    visit = set()
    for file in passed_file_names:
        # print("Length of visit is ",len(visit))
        print("Working on ", file)
        # print("echo Working on")
        # csv_path = os.path.join(current_dir, './data.csv' file)
        df = pd.read_csv( file , delimiter='\t')
        try:            
            # Get the value of the first row of the 'column_name' column
            dataset = df['MassiveID'].iloc[0]

            ccms_df = pd.read_csv(ccms_peak_link_start + dataset + ccms_peak_link_end)
        except TypeError:
            print(f"Skipping file {file} due to a TypeError.")
            continue

        for index, row in ccms_df.iterrows():
                # dataset = row["dataset"]
                filepath = row["filepath"]
                ccms_filenames[dataset].add(filepath)
        
        for index, row in df.iterrows():
            filename = row["filename"]
            filename2 = row["filename"][:-3] + "ML"
            for key in ccms_filenames[dataset]:
                if key.endswith(filename) or key.endswith(filename2):
                    if key in visit:
                        #  dupli += 1
                        merged_rows[key] = []
                    else:
                        #  unique += 1
                        visit.add(key)
                        new_row = [key] + list(row)
                        merged_rows[key] = list(new_row)
        print("Worked on ", file)

    print("echo Merged Rows list complete")
    non_empty_values = [v for v in merged_rows.values() if v]
    print("Length of Entries in Final file -> ",len(non_empty_values))

    # Create a DataFrame from the list with headers
    fnt = pd.DataFrame(non_empty_values)

    # Save the DataFrame to a TSV file without column names
    fnt.to_csv('check.tsv', sep='\t', index=False)



if __name__ == '__main__':
    main()