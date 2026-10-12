# Import necessary libraries
import os
import pandas as pd
import urllib
import requests
import io 
import sys
import collections
import argparse
import time


def download_gnps_metadata(out_path="gnps_metadata.tsv", page_size=1000, max_pages=None, max_retries=3):
    """
    Downloads all GNPS metadata ending with 'gnps_metadata.tsv' via the Datasette API,
    paginating with _next and writing results incrementally as TSV.

    Each page is retried up to max_retries times before aborting.
    """
    base_url = "https://datasetcache.gnps2.org/datasette/database/filename.json"
    params = {
        "_sort": "filepath",
        "filepath__endswith": "gnps_metadata.tsv",
        "_size": page_size,
        "_shape": "objects",  # rows as list of dicts
    }

    wrote_header = False
    all_rows = []
    page = 0
    next_token = None

    # truncate file
    with open(out_path, "w", encoding="utf-8") as f_out:
        pass

    print("Starting paginated download...", file=sys.stderr)

    while True:
        if next_token:
            params["_next"] = next_token
        elif "_next" in params:
            del params["_next"]

        page += 1
        if max_pages and page > max_pages:
            print("Reached max_pages limit; stopping.", file=sys.stderr)
            break

        attempt = 0
        success = False
        while attempt < max_retries and not success:
            attempt += 1
            try:
                r = requests.get(base_url, params=params, timeout=120)
                if r.status_code != 200:
                    print(f"[page {page}] Attempt {attempt} bad status {r.status_code}", file=sys.stderr)
                    time.sleep(2)
                    continue
                payload = r.json()
                success = True
            except Exception as e:
                print(f"[page {page}] Attempt {attempt} failed: {e}", file=sys.stderr)
                time.sleep(2)

        if not success:
            print(f"[page {page}] ❌ All {max_retries} attempts failed; aborting.", file=sys.stderr)
            break

        rows = payload.get("rows", [])
        next_token = payload.get("next")

        if not rows:
            print(f"[page {page}] No rows; stopping.", file=sys.stderr)
            break

        df_chunk = pd.DataFrame.from_records(rows)

        # write incrementally
        with open(out_path, "a", encoding="utf-8") as f_out:
            df_chunk.to_csv(f_out, sep="\t", index=False, header=not wrote_header)
        wrote_header = True

        all_rows.extend(rows)
        print(f"[page {page}] Wrote {len(rows)} rows "
              f"(next={'yes' if next_token else 'no'})", file=sys.stderr)

        if not next_token:
            break

        time.sleep(0.2)  # be nice to server

    if not all_rows:
        raise RuntimeError("No data downloaded; all pages failed.")

    df = pd.DataFrame.from_records(all_rows)
    print(f"✅ Saved {len(df)} rows to {out_path}", file=sys.stderr)
    return df

def _dataset_path(dataset, filepath):
    """Returns <dataset>/<filepath>, without doubling the prefix if filepath already has it."""
    if filepath.startswith(dataset + "/"):
        return filepath
    return dataset + "/" + filepath


def _download_metadata_file(session, url, max_retries=5, timeout=120):
    """Downloads one file from MassIVE, retrying on throttling/server errors.

    Returns the content, or None if the file could not be downloaded. Error pages are
    never returned, so they cannot end up saved as metadata.
    """
    for attempt in range(1, max_retries + 1):
        try:
            r = session.get(url, verify=False, timeout=timeout)
        except requests.RequestException as e:
            print(f"Attempt {attempt} failed for {url}: {e}")
        else:
            if r.status_code == 200:
                return r.content
            print(f"Attempt {attempt} bad status {r.status_code} for {url}")
            # 4xx other than throttling will not succeed on retry
            if r.status_code < 500 and r.status_code != 429:
                return None
        time.sleep(min(2 ** attempt, 60))
    return None


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

    # Names of GNPS Files, as dataset-relative paths prefixed with the dataset accession
    # (the cache's filepath column no longer carries the MSV prefix, but the MassIVE
    # download URL needs f.<MSV>/<path>)
    gnps_list = [_dataset_path(row['dataset'], row['filepath']) for row in selected_rows]
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
    failed_downloads = []
    session = requests.Session()

    # Download files from the links in the gnps_list and store them locally
    for index, link in enumerate(gnps_list):

        # get parentfolder of path in link
        msv_id = link.split(os.sep)[0]

        # Check if the files is already present from github
        if msv_id in existing_files:
            print(f"File {msv_id} already present in the folder from github!")
            continue

        download_url = download_link + link
        print(download_url)
        content = _download_metadata_file(session, download_url)
        if content is None:
            failed_downloads.append(link)
            continue
        file_name = os.path.join(args.output_metadata_folder, str(index) + "_gnps_metadata.tsv")
        file_paths["sys_name"].append(file_name)
        file_paths["svr_name"].append(link)
        with open(file_name,'wb') as f:
            f.write(content) #Downloading

    if failed_downloads:
        print(f"Failed to download {len(failed_downloads)} of {len(gnps_list)} metadata files:")
        for link in failed_downloads:
            print(f"  {link}")
    if gnps_list and len(failed_downloads) == len(gnps_list):
        raise RuntimeError("All MassIVE metadata downloads failed.")

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