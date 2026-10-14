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
from ftplib import FTP_TLS, error_perm, all_errors as ftp_errors


CACHE_FILENAME_URL = "https://datasetcache.gnps2.org/datasette/database/filename.json"
CACHE_SQL_URL = "https://datasetcache.gnps2.org/datasette/database.json"
# columns of the dataset cache's filename table, in table order
CACHE_COLUMNS = ["usi", "filepath", "dataset", "collection", "is_update", "update_name", "create_time", "size",
                 "size_mb", "sample_type", "spectra_ms1", "spectra_ms2", "instrument_vendor", "instrument_model",
                 "file_processed"]
# The cache answers each SQL query within 60 s (datasette sql_time_limit_ms) or not at all, and
# a 'filepath like ...' condition has to look at every one of its ~17 million rows; one query
# over the whole table took 3 s to over 60 s depending on the load of the cache (it failed on
# most attempts on 2026-10-14). Listings are therefore made of queries over rowid ranges of
# this many rows (~0.3 s each when idle); a range that still takes too long is split in two.
ROWID_CHUNK = 1000000
MIN_ROWID_CHUNK = 10000


class CacheQueryTooSlow(Exception):
    """The dataset cache stopped the query at its SQL time limit."""


def _cache_sql(sql, params, max_retries=6, timeout=(30, 180)):
    """Result rows of one read-only SQL query of the dataset cache (datasette JSON API)."""
    for attempt in range(1, max_retries + 1):
        try:
            r = requests.get(CACHE_SQL_URL, params={"sql": sql, "_shape": "objects", **params}, timeout=timeout)
            if r.status_code == 200:
                payload = r.json()
                if payload.get("ok", True) is not False:
                    return payload
                error = str(payload.get("error"))[:200]
            else:
                error = f"bad status {r.status_code}"
                if r.status_code == 400 and "too long" in r.text:
                    raise CacheQueryTooSlow(f"{params}")
                if 400 <= r.status_code < 500 and r.status_code not in (408, 429):
                    # a request error (e.g. invalid SQL) does not go away on retry
                    raise RuntimeError(f"Dataset cache query {params}: {error}: {r.text[:300]}")
        except (requests.RequestException, ValueError) as e:
            error = f"{type(e).__name__}: {e}"
        print(f"Dataset cache query {params}: attempt {attempt} failed ({error})", file=sys.stderr)
        if attempt < max_retries:
            time.sleep(min(10 * 2 ** (attempt - 1), 120))
    raise RuntimeError(f"Dataset cache query {params} failed {max_retries} times")


def list_cache_files(where, columns=("dataset", "filepath", "create_time", "size"), max_retries=6):
    """Files of the dataset cache's filename table matching an SQL condition, with their rowid
    (all of them, or an exception: a partial listing would silently drop the datasets it misses). The table is
    read in rowid ranges, in increasing order, the last range open-ended, so rows added while
    the listing is made are included; a row seen twice (re-inserted) is kept once."""
    max_rowid = _cache_sql("select max(rowid) as max_rowid from filename", {}, max_retries)["rows"][0]["max_rowid"] or 0
    ranges = collections.deque((lo, lo + ROWID_CHUNK if lo + ROWID_CHUNK <= max_rowid else None)
                               for lo in range(0, max_rowid + 1, ROWID_CHUNK))
    query_columns = ["rowid as _rowid", "usi"] + [c for c in columns if c != "usi"]
    sql = "select {} from filename where rowid >= :lo{} and ({})"
    rows = {}
    n_queries = n_split = 0
    while ranges:
        lo, hi = ranges.popleft()
        params = {"lo": lo} if hi is None else {"lo": lo, "hi": hi}
        size = ROWID_CHUNK if hi is None else hi - lo
        try:
            n_queries += 1
            payload = _cache_sql(sql.format(", ".join(query_columns), "" if hi is None else " and rowid < :hi", where),
                                 params, max_retries)
        except CacheQueryTooSlow:
            payload = None
        # too slow, or more rows than one answer can carry (max_returned_rows): split the range
        if payload is None or payload.get("truncated"):
            if size < 2 * MIN_ROWID_CHUNK:
                raise RuntimeError(f"Dataset cache: no complete answer for rowid range {params}")
            mid = lo + size // 2
            ranges.appendleft((mid, hi))
            ranges.appendleft((lo, mid))
            n_split += 1
            continue
        for row in payload.get("rows", []):
            rows[row["usi"]] = row
    print(f"Dataset cache listing ({where}): {len(rows)} files from {n_queries} queries "
          f"({n_split} ranges split)", file=sys.stderr)
    return pd.DataFrame.from_records(list(rows.values()), columns=list(columns) + ["_rowid"])


def download_gnps_metadata(out_path="gnps_metadata.tsv"):
    """All files named *gnps_metadata.tsv in the dataset cache (all columns), sorted by
    filepath, also written to out_path as TSV. Raises if the listing could not be completed."""
    df = list_cache_files("filepath like '%gnps_metadata.tsv'", CACHE_COLUMNS)
    if not len(df):
        raise RuntimeError("the dataset cache lists no *gnps_metadata.tsv files")
    # (the order of the former single query: by filepath, then rowid)
    df = df.sort_values(["filepath", "_rowid"], kind="stable").drop(columns="_rowid").reset_index(drop=True)
    df.to_csv(out_path, sep="\t", index=False)
    print(f"Saved {len(df)} rows to {out_path}", file=sys.stderr)
    return df


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


def _update_date(link):
    """Date of the MassIVE update a file <dataset>/<path> came with (updates/<YYYY-MM-DD>_<user>_<id>/),
    or '' for a file of the original submission, which is older than every update."""
    match = re.match(r'updates/(\d{4}-\d{2}-\d{2})_', link.partition('/')[2])
    return match.group(1) if match else ''


def order_versions(links, modify_time):
    """Metadata files <dataset>/<path> of one dataset, given in dataset cache create_time order
    (newest first), ordered newest first by when they were put on MassIVE. create_time is when
    the cache first saw a file: re-indexing in 2024 gave old files new times, which put an older
    version first (e.g. MSV000090940: its 2023-01-11 update before the 2023-06-09 one). The date
    of the update a file came with decides; files of the same date (or all of the original
    submission) are ordered by their modification time on MassIVE (modify_time(link), asked only
    for them), and keep the create_time order if that is not known for all of them."""
    dates = [_update_date(link) for link in links]
    tied = [link for link, date in zip(links, dates) if dates.count(date) > 1]
    times = {link: modify_time(link) for link in tied}
    if not all(times.values()):
        times = {}
    # (sorted keeps the given order of equal keys, also with reverse=True)
    return sorted(links, key=lambda link: (_update_date(link), times.get(link) or ''), reverse=True)


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
        self.index_complete = False
        # set by download(): MassIVE does not have the file (see there)
        self.missing = False

    def _connect(self):
        ftp = FTP_TLS(self.host, timeout=self.timeout)
        ftp.login("anonymous", "")
        ftp.prot_p()
        ftp.encoding = "utf-8"
        self.ftp = ftp

    def _build_index(self):
        index = collections.defaultdict(list)
        complete = True
        for root in sorted(self.ftp.nlst()):
            try:
                for entry in self.ftp.nlst(root):
                    index[entry.split('/')[-1]].append(root)
            except ftp_errors as e:
                print(f"Could not list FTPS root {root}: {e}")
                complete = False
        self.index = index
        self.index_complete = complete
        print(f"FTPS index built for {len(index)} accessions")

    def _ready(self):
        if self.ftp is None:
            self._connect()
        if self.index is None:
            self._build_index()

    def download(self, dataset_path, retries=3):
        """Returns the content of <dataset>/<path>, or None if it could not be fetched.
        Afterwards self.missing is True if MassIVE does not have the file: every FTPS version
        root of the dataset answered that it does not exist there (550). Such a file is only
        listed by the dataset cache (e.g. MSV000086415/updates/2021-04-01_jzemlin_2ca71553/,
        whose dataset holds only ccms_peak/ on MassIVE), so asking again elsewhere is futile."""
        self.missing = False
        dataset, _, file_path = dataset_path.partition('/')
        for attempt in range(1, retries + 1):
            try:
                self._ready()
                self.ftp.voidcmd("TYPE I")
                roots = self.index.get(dataset, [])
                for root in roots:
                    remote = f"{root}/{dataset}/{file_path}"
                    try:
                        self.ftp.size(remote)
                    except error_perm:
                        continue  # not in this root
                    buffer = io.BytesIO()
                    self.ftp.retrbinary(f"RETR {remote}", buffer.write, blocksize=1 << 20)
                    return buffer.getvalue()
                self.missing = bool(roots) and self.index_complete
                return None
            except (OSError, EOFError) + ftp_errors as e:
                print(f"FTPS attempt {attempt} failed for {dataset_path}: {e}")
                self.close()
                time.sleep(min(2 ** attempt, 30))
        return None

    def modify_time(self, dataset_path, retries=3):
        """Modification time of <dataset>/<path> on MassIVE (MDTM, 'YYYYMMDDHHMMSS' UTC), or None."""
        dataset, _, file_path = dataset_path.partition('/')
        for attempt in range(1, retries + 1):
            try:
                self._ready()
                for root in self.index.get(dataset, []):
                    try:
                        reply = self.ftp.sendcmd(f"MDTM {root}/{dataset}/{file_path}")
                    except error_perm:
                        continue  # not in this root
                    value = reply[4:].strip() if reply.startswith('213') else ''
                    return value if re.fullmatch(r'\d{14}(\.\d+)?', value) else None
                return None
            except (OSError, EOFError) + ftp_errors as e:
                print(f"FTPS attempt {attempt} (MDTM) failed for {dataset_path}: {e}")
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


    # If the MassIVE metadata files cannot be listed completely, none are read (an incomplete
    # listing would silently drop datasets) and the script exits with an error after the
    # GitHub metadata is passed on; the MassIVE rows of earlier runs are carried forward by
    # saveOlderData
    listing_errors = []
    try:
        gnps_df = download_gnps_metadata()
    except Exception as e:
        listing_errors.append(f"could not list the MassIVE *gnps_metadata.tsv files ({type(e).__name__}: {e})")
        print(f"ERROR: {listing_errors[-1]}; no MassIVE metadata is read this run, only the GitHub metadata", file=sys.stderr)
        gnps_df = pd.DataFrame(columns=['dataset', 'filepath', 'create_time'])

    # Print message to indicate that the CSV file has been read
    print("echo GNPS CSV Read Done!")

    # Convert the create time column to a datetime format
    gnps_df['create_time'] = pd.to_datetime(gnps_df['create_time'], errors='coerce')


    # All *gnps_metadata.tsv files of a dataset, by create time (missing times count as
    # oldest); ordered newest first by order_versions() before they are downloaded: the
    # newest version is read in full, older versions only add files the newer ones do not
    # describe
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
    missing_on_massive = []
    session = requests.Session()

    def fetch(link):
        content = ftps.download(link)
        if content is not None:
            downloaded_via["FTPS"] += 1
            return content
        if ftps.missing:
            # the dataset cache lists a file MassIVE does not have (any more); the HTTP
            # endpoint serves the same files, so it is not asked
            print(f"Not on MassIVE (no FTPS version root of the dataset has it; listed only by the dataset cache): {link}")
            missing_on_massive.append(link)
            return None
        content = _download_metadata_file(session, download_link + link)
        if content is not None:
            downloaded_via["HTTP"] += 1
        else:
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

        if len(links) > 1:
            links = order_versions(links, ftps.modify_time)
        print(f"Downloading {links[0]}")
        content = fetch(links[0])
        # a newest file MassIVE does not have is no version: the next one is the newest
        while content is None and ftps.missing and len(links) > 1:
            links = links[1:]
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
    if not args.no_other_sheets and not listing_errors:
        try:
            other = list_cache_files(OTHER_SHEETS_WHERE)
        except Exception as e:
            listing_errors.append(f"could not list the other metadata sheets ({type(e).__name__}: {e})")
            print(f"ERROR: {listing_errors[-1]}; continuing without them", file=sys.stderr)
            other = pd.DataFrame(columns=["dataset", "filepath", "create_time", "size"])
        # MassIVE datasets only: the cache also lists Workbench (ST...) and MetaboLights
        # (MTBLS...) files, which MassIVE does not serve (each would only cost futile requests)
        other = other[other['dataset'].astype(str).str.match(r'^R?MSV\d+$')]
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
            if len(sheets) > 1:
                by_link = {link: (df, link) for df, link in sheets}
                sheets = [by_link[link] for link in order_versions([link for _, link in sheets], ftps.modify_time)]
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
    if missing_on_massive:
        print(f"{len(missing_on_massive)} metadata files listed by the dataset cache are not on MassIVE:")
        for link in missing_on_massive:
            print(f"  {link}")
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
    
    if listing_errors:
        for error in listing_errors:
            print(f"ERROR: {error}", file=sys.stderr)
        sys.exit(2)

    # Print message to indicate that the execution of the code is complete
    print("echo EXECUTION COMPLETE ! HAVE A GOOD DAY ")

    
    

if __name__ == "__main__":
  main()