# Import necessary libraries
import os
import pandas as pd
import urllib
import requests
import io 
import sys
import collections
import argparse


def main():
    # Args parse
    parser = argparse.ArgumentParser(description='Download GNPS files')
    parser.add_argument('output_metadata_folder')
    args = parser.parse_args()

    # Print message to indicate importing is done
    print("echo Importing Done!")

    # Set the URL for GNPS metadata CSV file and read it using pandas
    gnps_metadata_link = "https://gnps-datasetcache.ucsd.edu/database/filename.csv?_sort=filepath&filepath__endswith=gnps_metadata.tsv&_size=max"
    gnps_df = pd.read_csv(gnps_metadata_link)

    # Print message to indicate that the CSV file has been read
    print("echo GNPS CSV Read Done!")

    # Convert the create time column to a datetime format
    gnps_df['create_time'] = pd.to_datetime(gnps_df['create_time'])

    # Sort the DataFrame by the create time column
    gnps_df = gnps_df.sort_values(by='create_time')

    # Group the DataFrame by the dataset column
    groups = gnps_df.groupby('dataset')

    # Select the row with the highest create time value for each group
    selected_rows = []
    for name, group in groups:
        idx = group['create_time'].idxmax()
        selected_rows.append(group.loc[idx])

    # Names of GNPS Files
    gnps_list = [row['filepath'] for row in selected_rows]
    print("Total gnps file are ", len(gnps_list))
    os.system("echo Filepath list generated ... ")

    # Set the download link for GNPS files
    download_link = "https://massive.ucsd.edu/ProteoSAFe/DownloadResultFile?forceDownload=true&file=f."
    print("echo We are downloading now ...")

    # Create a defaultdict to store file paths
    file_paths = collections.defaultdict(list)

    # Download files from the links in the gnps_list and store them locally
    for index, link in enumerate(gnps_list):
        download_url = download_link + link
        r = requests.get(download_url)
        file_name = os.path.join(args.output_metadata_folder, str(index) + "_gnps_metadata.tsv")
        file_paths["sys_name"].append(file_name)
        file_paths["svr_name"].append(link)
        with open(file_name,'wb') as f:
            f.write(r.content) #Downloading

    # Print message to indicate that downloading has been completed successfully
    print("echo Download has been completed successfully!")

    # Print message to indicate that a TSV file for file paths is being created
    print("echo Creating tsv for file paths ...")

    # Write file paths to a TSV file
    # TODO: Use pandas to write the TSV file    
    with open("file_paths.tsv", "w") as f:
        f.write("LocalName\tServerName\n")
        for i in range(len(file_paths["sys_name"])):
            local = file_paths["sys_name"][i]
            server = file_paths["svr_name"][i]
            f.write(f"{local}\t{server}\n")
        f.close()
    
    # Print message to indicate that the TSV file for file paths has been created
    print("echo Tsv file created for path names!!")
    
    # Print message to indicate that the execution of the code is complete
    print("echo EXECUTION COMPLETE ! HAVE A GOOD DAY ")

    
    

if __name__ == "__main__":
  main()