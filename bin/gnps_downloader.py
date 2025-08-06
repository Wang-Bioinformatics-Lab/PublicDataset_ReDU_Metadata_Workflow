# Import necessary libraries
import os
import pandas as pd
import urllib
import requests
import io 
import sys
import collections
import argparse


def download_gnps_metadata(temp_file_path="temp.csv"):
    base_url = "https://datasetcache.gnps2.org/datasette/database/filename.csv"
    params = {
        "_sort": "filepath",
        "filepath__endswith": "gnps_metadata.tsv",
        "_size": 1000,
    }

    all_chunks = []
    next_url = base_url
    first = True

    print("Starting paginated download...")

    while next_url:
        print(f"Fetching: {next_url}")
        try:
            response = requests.get(next_url, params=params if first else {}, timeout=60)
        except Exception as e:
            print(f"Request failed: {e}")
            break

        if response.status_code != 200 or "html" in response.headers.get("Content-Type", ""):
            print("Bad response, likely server error:")
            print(response.text[:500])
            break

        df_chunk = pd.read_csv(io.StringIO(response.text))
        if df_chunk.empty:
            print("No more data.")
            break

        all_chunks.append(df_chunk)

        # Check for next page
        next_link = None
        if "_next=" in response.url:
            next_link = response.url.split("_next=")[-1]
            next_url = f"{base_url}?_sort=filepath&filepath__endswith=gnps_metadata.tsv&_size=1000&_next={next_link}"
        else:
            break

        first = False

    # Combine and save to temp file
    gnps_df = pd.concat(all_chunks, ignore_index=True)
    gnps_df.to_csv(temp_file_path, index=False)
    print(f"Saved metadata to {temp_file_path}")
    return gnps_df

def main():
    # Args parse
    parser = argparse.ArgumentParser(description='Download GNPS files')
    parser.add_argument('output_metadata_folder')
    args = parser.parse_args()

    # Print message to indicate importing is done
    print("echo Importing Done!")


    gnps_df = download_gnps_metadata()

    # Print message to indicate that the CSV file has been read
    print("echo GNPS CSV Read Done!")

    # Convert the create time column to a datetime format
    gnps_df['create_time'] = pd.to_datetime(gnps_df['create_time'], errors='coerce')


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

    # Get already present github files by listing file paths in the output_metadata_folder
    existing_files = os.listdir(args.output_metadata_folder)

    # Get list of basenames of existing files without file extension
    existing_files = [os.path.basename(file).split(".")[0] for file in existing_files]

    # Set the download link for GNPS files
    download_link = "https://massive.ucsd.edu/ProteoSAFe/DownloadResultFile?forceDownload=true&file=f."
    print("echo We are downloading now ...")

    # Create a defaultdict to store file paths
    file_paths = collections.defaultdict(list)

    # Download files from the links in the gnps_list and store them locally
    for index, link in enumerate(gnps_list):

        # get parentfolder of path in link
        msv_id = link.split(os.sep)[0]

        # Check if the files is already present from github
        if msv_id in existing_files:
            print(f"File {msv_id} already present in the folder from github!")
            continue

        download_url = download_link + link
        r = requests.get(download_url, verify=False)
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