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
from ftplib import FTP_TLS, all_errors as ftp_errors


def download_gnps_metadata(out_path="gnps_metadata.tsv", page_size=1000, max_pages=None, max_retries=5):
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
                    time.sleep(2 ** attempt)
                    continue
                payload = r.json()
                success = True
            except Exception as e:
                print(f"[page {page}] Attempt {attempt} failed: {e}", file=sys.stderr)
                time.sleep(2 ** attempt)

        if not success:
            # A partial listing would silently drop every dataset on the missing pages
            raise RuntimeError(f"[page {page}] All {max_retries} attempts failed; aborting.")

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


class MassiveFTPS:
    """Anonymous explicit-FTPS access to MassIVE (massive-ftp.ucsd.edu).

    Files live under version roots (v01.., x01, z01, ..) and one dataset can be split
    across several of them, so an accession -> roots index is built once.
    Recommended over the throttled ProteoSAFe/DownloadResultFile HTTP endpoint.
    """

    def __init__(self, host="massive-ftp.ucsd.edu", timeout=120):
        self.host = host
        self.timeout = timeout
        self.ftp = None
        self.index = None

    def _connect(self):
        ftp = FTP_TLS(self.host, timeout=self.timeout)
        ftp.login("anonymous", "")
        ftp.prot_p()
        ftp.encoding = "utf-8"
        self.ftp = ftp

    def _build_index(self):
        index = collections.defaultdict(list)
        for root in sorted(self.ftp.nlst()):
            try:
                for entry in self.ftp.nlst(root):
                    index[entry.split('/')[-1]].append(root)
            except ftp_errors as e:
                print(f"Could not list FTPS root {root}: {e}")
        self.index = index
        print(f"FTPS index built for {len(index)} accessions")

    def download(self, dataset_path, retries=3):
        """Returns the content of <dataset>/<path>, or None if it could not be fetched."""
        dataset, _, file_path = dataset_path.partition('/')
        for attempt in range(1, retries + 1):
            try:
                if self.ftp is None:
                    self._connect()
                if self.index is None:
                    self._build_index()
                self.ftp.voidcmd("TYPE I")
                for root in self.index.get(dataset, []):
                    remote = f"{root}/{dataset}/{file_path}"
                    try:
                        self.ftp.size(remote)
                    except ftp_errors:
                        continue
                    buffer = io.BytesIO()
                    self.ftp.retrbinary(f"RETR {remote}", buffer.write, blocksize=1 << 20)
                    return buffer.getvalue()
                return None
            except (OSError, EOFError) + ftp_errors as e:
                print(f"FTPS attempt {attempt} failed for {dataset_path}: {e}")
                self.close()
                time.sleep(min(2 ** attempt, 30))
        return None

    def close(self):
        if self.ftp is not None:
            try:
                self.ftp.close()
            except Exception:
                pass
        self.ftp = None


def main():
    # Args parse
    parser = argparse.ArgumentParser(description='Download GNPS files')
    parser.add_argument('output_metadata_folder')
    args = parser.parse_args()

    # Print message to indicate importing is done
    print("echo Importing Done!")


    # If the dataset cache is unavailable, continue with the GitHub metadata only; the
    # MassIVE rows of earlier runs are carried forward by saveOlderData
    try:
        gnps_df = download_gnps_metadata()
    except Exception as e:
        print(f"WARNING: could not list MassIVE metadata files ({e}); continuing with GitHub metadata only")
        gnps_df = pd.DataFrame(columns=['dataset', 'filepath', 'create_time'])

    # Print message to indicate that the CSV file has been read
    print("echo GNPS CSV Read Done!")

    # Convert the create time column to a datetime format
    gnps_df['create_time'] = pd.to_datetime(gnps_df['create_time'], errors='coerce')


    # Sort the DataFrame by the create time column (missing times first, so they are
    # only selected when a dataset has no file with a valid create time)
    gnps_df = gnps_df.sort_values(by='create_time', na_position='first')

    # Group the DataFrame by the dataset column
    groups = gnps_df.groupby('dataset')

    # Select the row with the highest create time value for each group
    selected_rows = []
    for name, group in groups:
        selected_rows.append(group.iloc[-1])

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

    # Files are fetched over FTPS; the (throttled) HTTP endpoint is only a fallback
    download_link = "https://massive.ucsd.edu/ProteoSAFe/DownloadResultFile?forceDownload=true&file=f."
    ftps = MassiveFTPS()
    downloaded_via = collections.Counter()
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

        print(f"Downloading {link}")
        content = ftps.download(link)
        if content is not None:
            downloaded_via["FTPS"] += 1
        else:
            content = _download_metadata_file(session, download_link + link)
            if content is not None:
                downloaded_via["HTTP"] += 1
        if content is None:
            failed_downloads.append(link)
            continue
        file_name = os.path.join(args.output_metadata_folder, str(index) + "_gnps_metadata.tsv")
        file_paths["sys_name"].append(file_name)
        file_paths["svr_name"].append(link)
        with open(file_name,'wb') as f:
            f.write(content) #Downloading

    ftps.close()
    print(f"Downloaded metadata files: {dict(downloaded_via)}")
    if failed_downloads:
        print(f"WARNING: failed to download {len(failed_downloads)} of {len(gnps_list)} metadata files:")
        for link in failed_downloads:
            print(f"  {link}")

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