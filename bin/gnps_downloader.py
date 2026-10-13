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
import re
import json
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

CACHE_FILENAME_URL = "https://datasetcache.gnps2.org/datasette/database/filename.json"
PROXI_URL = "https://massive.ucsd.edu/ProteoSAFe/proxi/v0.1/datasets/{}"

# ReDU-format sheets that are not named *gnps_metadata.tsv (other names, .txt, .csv,
# .xlsx, .tsv.tsv): tabular files with 'metadata' in their path (name or folder), outside
# ccms_metadata/. Only those whose header is the ReDU template are used.
TEMPLATE_SHEET_EXTENSIONS = ('.tsv', '.txt', '.csv', '.xlsx', '.xls')
OTHER_SHEETS_WHERE = (
    "filepath like '%metadata%' "
    "and filepath not like '%ccms_metadata/%' and filepath not like '%gnps_metadata.tsv' and ("
    + " or ".join(f"filepath like '%{ext}'" for ext in TEMPLATE_SHEET_EXTENSIONS) + ")")
MAX_SHEET_BYTES = 50 * 1024 * 1024
# a sheet counts as ReDU-format when it has a filename column and at least this many
# ReDU template columns (current names or first-generation names)
MIN_TEMPLATE_COLUMNS = 4

# first-generation template column names (prefixes in front of current names, and names
# that differ) -> current names
LEGACY_PREFIXES = ('ATTRIBUTE_Curated_', 'ATTRIBUTE_', 'Subject_', 'Analysis_', 'Study_', 'Curated_')
LEGACY_ALIASES = {
    'sex': 'BiologicalSex', 'health': 'HealthStatus', 'bodypartontologyname': 'UBERONBodyPartName',
    'bmi': 'BodyMassIndex', 'diseasecommonname': 'DOIDCommonName', 'smoking': 'SmokingStatus',
    'sampletype_sub1': 'SampleTypeSub1', 'comorbidity_listdoids': 'ComorbidityListDOIDIndex',
    'massiveid': 'MassiveID', 'datasetaccession': 'MassiveID',
}


def template_column(column, template_keys):
    """The current ReDU template column a sheet column stands for, or None.
    Same rules as normalize_gnps_metadata in read_and_validate_redu_from_github.py."""
    by_lower = {k.lower(): k for k in template_keys}
    name = str(column).strip()
    while True:
        if name in template_keys:
            return name
        if name.lower() in by_lower:
            return by_lower[name.lower()]
        if name.lower() in LEGACY_ALIASES:
            return LEGACY_ALIASES[name.lower()]
        for prefix in LEGACY_PREFIXES:
            if name.lower().startswith(prefix.lower()) and len(name) > len(prefix):
                name = name[len(prefix):]
                break
        else:
            return None


def is_template_sheet(df, template_keys):
    targets = {template_column(c, template_keys) for c in df.columns} - {None}
    return 'filename' in targets and len(targets) >= MIN_TEMPLATE_COLUMNS


def _excel_sheet(content, template_keys=None):
    """The worksheet of an Excel file to use: the first one, unless only another sheet has
    the ReDU template columns (e.g. instructions on sheet 1, the table on sheet 2)."""
    sheets = pd.read_excel(io.BytesIO(content), dtype=str, sheet_name=None)
    if not sheets:
        return None
    for df in sheets.values():
        df.columns = [str(c).strip() for c in df.columns]
    frames = list(sheets.values())
    if template_keys is None or is_template_sheet(frames[0], template_keys):
        return frames[0]
    redu = [df for df in frames[1:] if is_template_sheet(df, template_keys)]
    return redu[0] if len(redu) == 1 else frames[0]


def read_sheet(content, name, template_keys=None):
    """Reads a tabular metadata file (tsv/txt/csv/xlsx) as strings, or returns None."""
    lower = name.lower()
    try:
        if lower.endswith(('.xlsx', '.xls')):
            # needs openpyxl (.xlsx) / xlrd (.xls); without them the file is skipped
            return _excel_sheet(content, template_keys)
        for encoding in ('utf-8', 'latin-1'):
            try:
                text = content.decode(encoding)
                break
            except UnicodeDecodeError:
                continue
        text = text.lstrip('\ufeff')
        separators = [',', '\t', ';'] if lower.endswith('.csv') else ['\t', ',']
        for sep in separators:
            df = pd.read_csv(io.StringIO(text), sep=sep, dtype=str)
            if df.shape[1] > 1:
                df.columns = [str(c).strip() for c in df.columns]
                return df
    except Exception as e:
        print(f"Could not read {name}: {type(e).__name__}: {e}")
    return None


def _canonical_column(df, target, template_keys):
    """Renames the column standing for target (e.g. 'ATTRIBUTE_filename') to target."""
    if target in df.columns:
        return df
    for column in df.columns:
        if template_column(column, template_keys) == target:
            return df.rename(columns={column: target})
    return df


def _sample_key(filename):
    """File name without folders, upper/lower case and data-file extension, so a newer
    version listing x.mzML also describes the x.mzXML of an older version."""
    base = os.path.basename(str(filename).strip().replace('\\', '/')).lower()
    return re.sub(r'\.(mzml|mzxml|raw|d|wiff|cdf|mzdata|mgf)$', '', base)


def combine_versions(sheets, dataset, template_keys):
    """One table from several metadata files of a dataset, given newest first. Rows of the
    newest file are all kept; an older file only adds rows of this dataset for files (same
    file name, ignoring folders, case and .mzML/.mzXML) that no newer file describes."""
    described = set()
    parts = []
    for rank, (df, source) in enumerate(sheets):
        df = _canonical_column(df, 'filename', template_keys)
        df = _canonical_column(df, 'MassiveID', template_keys)
        if 'filename' not in df.columns:
            continue
        if 'MassiveID' not in df.columns:
            df = df.assign(MassiveID=dataset)
        df = df.copy()
        accession = df['MassiveID'].fillna('').astype(str).str.strip()
        accession = accession.where(accession.str.match(r'^R?MSV\d+$'), dataset)
        if rank > 0:
            df['MassiveID'] = accession
        keys = list(zip(accession, df['filename'].map(_sample_key)))
        # rows without a file name (empty template rows) describe nothing; rows an older
        # file has for another dataset are left to that dataset's own metadata
        has_name = ~df['filename'].fillna('').astype(str).str.strip().str.lower().isin(['', 'nan'])
        own = accession == dataset
        new = [rank == 0 or (named and mine and k not in described) for k, named, mine in zip(keys, has_name, own)]
        described.update(k for k, named in zip(keys, has_name) if named)
        kept = df[new]
        if rank > 0:
            print(f"  {source}: {int(sum(new))} of {len(df)} rows describe files no newer version describes")
        if len(kept):
            parts.append(kept)
    if not parts:
        return None
    return pd.concat(parts, ignore_index=True, sort=False)


def list_cache_files(where, page_size=1000, max_retries=5):
    """Files of the dataset cache matching an SQL condition on filepath (all pages)."""
    params = {"_where": where, "_shape": "objects", "_size": page_size,
              "_col": ["filepath", "dataset", "create_time", "size"]}
    rows = []
    next_token = None
    page = 0
    while True:
        page += 1
        if next_token:
            params["_next"] = next_token
        payload = None
        for attempt in range(1, max_retries + 1):
            try:
                r = requests.get(CACHE_FILENAME_URL, params=params, timeout=300)
                if r.status_code == 200:
                    payload = r.json()
                    break
                print(f"[page {page}] Attempt {attempt} bad status {r.status_code}", file=sys.stderr)
            except Exception as e:
                print(f"[page {page}] Attempt {attempt} failed: {e}", file=sys.stderr)
            time.sleep(2 ** attempt)
        if payload is None:
            raise RuntimeError(f"[page {page}] All {max_retries} attempts failed")
        rows.extend(payload.get("rows", []))
        next_token = payload.get("next")
        if not next_token or not payload.get("rows"):
            break
        time.sleep(0.2)
    return pd.DataFrame.from_records(rows, columns=["dataset", "filepath", "create_time", "size"])


# Proteomics datasets are not taken in through other-named sheets. Same signals and
# order as the survey's dataset classification, from the PROXI record and the cache.
_METABOLOMICS_TERMS = re.compile(r'(?i)metabolom|metabolite|lipidom|\blipid|natural product|small[- ]molecule|exposom|molecular network|\bFBMN\b|feature[- ]based|mzmine|spectral librar|secondary metabol|untargeted|\bLC-?MS/?MS\b.*extract|\bmetabolic profil|volatil|\bDOM\b|dissolved organic|xenobiotic|drug metabol|bile acid|glycan|oligosacchar|polyketide|alkaloid|terpen|flavon|phenolic|peptidic natural|siderophore|lipopeptide|\bMASST\b|chemical standard|reference standard|food|\bdiet\b|plant extract|fecal|stool')
_PROTEOMICS_TERMS = re.compile(r'(?i)proteom|peptid(?!ic natural)|phospho|\bTMT\b|iTRAQ|SILAC|ubiquit|\bSUMO|glycoproteo|cross-?link|interactom|immunopeptid|top-?down|\bPTM|tryp(sin|tic)|\bdigest|histone|acetylom|proteoform|\bDIA\b|SWATH|immunoprecipitation|\bIP-MS|\bAP-MS|BioID|kinase|protein|\bHLA\b|\bMHC\b|antibod|\bPRM\b|\bSRM\b|MS/MS of proteins|deep learning.*spectra|PXD\d+')
_PROTEOMICS_MODIFICATIONS = re.compile(r'(?i)UNIMOD|\bMOD:|oxidation|carbamidomethyl|phospho|acetyl|TMT|iTRAQ|deamid|methyl')
PROTEOMICS_FILE_EXTENSIONS = ('.fasta', '.mzid', '.mzidentml', '.pepxml', '.protxml', '.msf', '.pdresult', '.dat', '.blib')


def fetch_proxi(session, accession, retries=3):
    for attempt in range(1, retries + 1):
        try:
            r = session.get(PROXI_URL.format(accession), timeout=(10, 60))
            if r.status_code == 200:
                return json.loads(r.content.decode("utf-8", errors="replace"))
            if r.status_code < 500 and r.status_code != 429:
                return None
        except (requests.RequestException, ValueError) as e:
            print(f"PROXI attempt {attempt} failed for {accession}: {e}")
        time.sleep(min(2 ** attempt, 30))
    return None


def cache_dataset_signals(accession):
    """(number of proteomics identification files, set of cache sample types) of a dataset."""
    def query(extra):
        r = requests.get(CACHE_FILENAME_URL, params={"dataset": accession, "_size": 0, **extra}, timeout=300)
        r.raise_for_status()
        return r.json()
    ext_where = " or ".join(f"filepath like '%{ext}'" for ext in PROTEOMICS_FILE_EXTENSIONS)
    n_proteomics_files = int(query({"_where": f"({ext_where})"}).get("filtered_table_rows_count") or 0)
    facet = query({"_facet": "sample_type"}).get("facet_results", {}).get("sample_type", {})
    sample_types = {str(f.get("value")) for f in facet.get("results", []) if f.get("count")}
    return n_proteomics_files, sample_types


def dataset_class(proxi, n_proteomics_files, sample_types):
    """'metabolomics', 'proteomics' or 'unclear' for a dataset without GNPS/ReDU metadata."""
    if not proxi:
        return 'unclear'
    def text(entries):
        return '; '.join(str(e.get('value') or '') for e in entries or [] if isinstance(e, dict))
    title = str(proxi.get('title') or '')
    keywords = text(proxi.get('keywords'))
    summary = str(proxi.get('summary') or '')[:600]
    modifications = '; '.join(str(m.get('name') or '') + ' ' + str(m.get('value') or '')
                              for m in proxi.get('modifications') or [] if isinstance(m, dict)
                              and m.get('accession') != 'MS:1002864')
    title_keywords = title + ' ; ' + keywords
    dataset_type = re.search(r'DatasetType:\s*([^;]*)', keywords)
    dataset_type = dataset_type.group(1).strip() if dataset_type else ''
    if re.match(r'(?i)metabolomics|other \((.*lipidom|.*metabolom|exposom|xenobiot|pfas|glycomics|np|natural)', dataset_type):
        return 'metabolomics'
    if re.match(r'(?i)proteomics|other \((.*prote|.*peptid|hdx|native|intact|.*protein|msp-ms)', dataset_type):
        return 'proteomics'
    metabolomics_kw = bool(_METABOLOMICS_TERMS.search(title_keywords))
    proteomics = n_proteomics_files > 0 or bool(_PROTEOMICS_MODIFICATIONS.search(modifications)) \
        or bool(_PROTEOMICS_TERMS.search(title_keywords))
    gnps = 'GNPS' in sample_types or bool(re.search(r'(?i)(^|[^a-z])GNPS', title)) \
        or bool(re.search(r'(?i)^(MetaboLights|Metabolomics Workbench) (MTBLS|ST)\d+', title))
    if gnps and not proteomics:
        return 'metabolomics'
    if gnps and proteomics:
        return 'metabolomics' if metabolomics_kw and not (n_proteomics_files or _PROTEOMICS_MODIFICATIONS.search(modifications)) else 'unclear'
    if metabolomics_kw and not proteomics:
        return 'metabolomics'
    if proteomics:
        return 'proteomics' if not metabolomics_kw else 'unclear'
    if _PROTEOMICS_TERMS.search(summary) and not _METABOLOMICS_TERMS.search(summary):
        return 'proteomics'
    if 'MASSIVE' in sample_types and 'GNPS' not in sample_types and not _METABOLOMICS_TERMS.search(summary):
        return 'proteomics'
    if _METABOLOMICS_TERMS.search(summary) and not _PROTEOMICS_TERMS.search(summary) and 'MASSIVE' not in sample_types:
        return 'metabolomics'
    return 'unclear'


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
    parser.add_argument('--allowed_terms', default=os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'data', 'allowed_terms.json'),
                        help='allowed terms json (only its column names are used, to recognise ReDU-format sheets)')
    parser.add_argument('--no_other_sheets', action='store_true', help='only read files named *gnps_metadata.tsv')
    args = parser.parse_args()

    with open(args.allowed_terms) as f:
        template_keys = set(json.load(f).keys())

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


    # All *gnps_metadata.tsv files of a dataset, newest first (by create time; missing
    # times count as oldest): the newest version is read in full, older versions only
    # add files the newer ones do not describe
    versions = collections.OrderedDict()
    for name, group in gnps_df.groupby('dataset'):
        group = group.sort_values(['create_time', 'filepath'], na_position='first')
        versions[name] = [_dataset_path(name, fp) for fp in reversed(group['filepath'].tolist())]

    # Names of GNPS Files, as dataset-relative paths prefixed with the dataset accession
    # (the cache's filepath column no longer carries the MSV prefix, but the MassIVE
    # download URL needs f.<MSV>/<path>)
    gnps_list = [paths[0] for paths in versions.values()]
    print("Total gnps file are ", len(gnps_list), f"({sum(len(v) for v in versions.values())} versions)")
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

    def fetch(link):
        content = ftps.download(link)
        if content is not None:
            downloaded_via["FTPS"] += 1
        else:
            content = _download_metadata_file(session, download_link + link)
            if content is not None:
                downloaded_via["HTTP"] += 1
        if content is None:
            failed_downloads.append(link)
        return content

    def save(msv_id, index, content, links):
        # The dataset accession leads the name so the harmonization step can recover it
        # for files that lack a MassiveID column
        file_name = os.path.join(args.output_metadata_folder, f"{msv_id}_{index}_gnps_metadata.tsv")
        file_paths["sys_name"].append(file_name)
        file_paths["svr_name"].append(';'.join(links))
        with open(file_name, 'wb') as f:
            f.write(content)

    # Download files from the links in the gnps_list and store them locally
    for index, (msv_id, links) in enumerate(versions.items()):

        # Check if the files is already present from github
        if msv_id in existing_files:
            print(f"File {msv_id} already present in the folder from github!")
            continue

        print(f"Downloading {links[0]}")
        content = fetch(links[0])
        if content is None:
            continue
        if len(links) == 1:
            save(msv_id, index, content, links)
            continue

        # older versions: the newest file is kept as downloaded when they add nothing
        sheets = []
        newest = read_sheet(content, links[0])
        if newest is None:
            print(f"WARNING: could not parse {links[0]}; older versions are not combined with it")
            save(msv_id, index, content, links[:1])
            continue
        sheets.append((newest, links[0]))
        for link in links[1:]:
            print(f"Downloading older version {link}")
            older = fetch(link)
            older = read_sheet(older, link) if older is not None else None
            if older is not None:
                sheets.append((older, link))
        combined = combine_versions(sheets, msv_id, template_keys)
        if combined is None or len(combined) <= len(newest):
            save(msv_id, index, content, links[:1])
        else:
            print(f"{msv_id}: {len(combined) - len(newest)} rows from older versions added to {len(newest)}")
            save(msv_id, index, combined.to_csv(sep='\t', index=False).encode('utf-8'), [l for _, l in sheets])

    # ReDU-format sheets with other names, in datasets that have no *gnps_metadata.tsv
    # and no GitHub file
    if not args.no_other_sheets:
        try:
            other = list_cache_files(OTHER_SHEETS_WHERE)
        except Exception as e:
            print(f"WARNING: could not list other metadata sheets ({e}); continuing without them")
            other = pd.DataFrame(columns=["dataset", "filepath", "create_time", "size"])
        other = other[~other['dataset'].isin(set(versions) | set(existing_files))]
        extensions = TEMPLATE_SHEET_EXTENSIONS
        for extension, module in (('.xlsx', 'openpyxl'), ('.xls', 'xlrd')):
            try:
                __import__(module)  # pandas needs it to read this Excel format
            except ImportError:
                print(f"WARNING: {module} is not installed; {extension} sheets are not checked")
                extensions = tuple(e for e in extensions if e != extension)
        other = other[other['filepath'].str.lower().str.endswith(extensions)]
        other = other[pd.to_numeric(other['size'], errors='coerce').fillna(0) <= MAX_SHEET_BYTES]
        other['create_time'] = pd.to_datetime(other['create_time'], errors='coerce')
        print(f"Checking {len(other)} other metadata sheets in {other['dataset'].nunique()} datasets")
        index = len(versions)
        n_sheets = collections.Counter()
        for msv_id, group in other.groupby('dataset'):
            group = group.sort_values(['create_time', 'filepath'], na_position='first')
            sheets = []
            for filepath in reversed(group['filepath'].tolist()):
                link = _dataset_path(msv_id, filepath)
                content = fetch(link)
                df = read_sheet(content, link, template_keys) if content is not None else None
                if df is None:
                    n_sheets['unreadable'] += 1
                elif not is_template_sheet(df, template_keys):
                    n_sheets['not_redu_format'] += 1
                else:
                    sheets.append((df, link))
            if not sheets:
                continue
            try:
                proxi = fetch_proxi(session, msv_id)
                n_proteomics_files, sample_types = cache_dataset_signals(msv_id)
                kind = dataset_class(proxi, n_proteomics_files, sample_types)
            except Exception as e:
                print(f"Could not classify {msv_id} ({e}); its sheets are not used")
                kind = 'unclear'
            if kind != 'metabolomics':
                print(f"Skipping ReDU-format sheets of {msv_id}: dataset classified as {kind}")
                n_sheets[f'redu_format_{kind}_dataset'] += len(sheets)
                continue
            combined = combine_versions(sheets, msv_id, template_keys)
            if combined is None:
                continue
            n_sheets['redu_format_used'] += len(sheets)
            print(f"{msv_id}: {len(combined)} rows from ReDU-format sheets {[l for _, l in sheets]}")
            save(msv_id, index, combined.to_csv(sep='\t', index=False).encode('utf-8'), [l for _, l in sheets])
            index += 1
        print(f"Other metadata sheets: {dict(n_sheets)}")

    ftps.close()
    print(f"Downloaded metadata files: {dict(downloaded_via)}")
    if failed_downloads:
        print(f"WARNING: failed to download {len(failed_downloads)} of {len(failed_downloads) + sum(downloaded_via.values())} metadata files:")
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