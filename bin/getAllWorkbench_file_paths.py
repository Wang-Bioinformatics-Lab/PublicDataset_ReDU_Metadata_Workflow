import os
import re
import pandas as pd
import requests
import argparse
import tqdm
from urllib.parse import urlparse, parse_qs, quote_plus
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry


# Shared session: retries throttling/server errors with backoff instead of losing the study
session = requests.Session()
session.mount("https://", HTTPAdapter(max_retries=Retry(total=5, backoff_factor=2,
                                                         status_forcelist=[429, 500, 502, 503, 504],
                                                         allowed_methods=["GET"])))
REQUEST_TIMEOUT = (10, 300)


def _get(url):
    response = session.get(url, timeout=REQUEST_TIMEOUT)
    response.raise_for_status()
    return response


def _get_json(url):
    return _get(url).json()


def listing_query(url):
    """Query parameters of a show_archive_contents_json.php URL. Archive names may contain a
    raw '#' (ST001957_Method#2_HILIC_Pos_Raw.7z), which is no fragment marker here."""
    return parse_qs(urlparse(url, allow_fragments=False).query)


# --- archive listing -------------------------------------------------------------------
ARCHIVE_LISTING_URL = "https://www.metabolomicsworkbench.org/data/show_archive_contents_json.php?STUDY_ID={}"
# the same listings as text (unzip -l / 7z l output), for studies whose JSON listing fails
ARCHIVE_PAGE_URL = "https://www.metabolomicsworkbench.org/data/show_archive_contents.php?STUDY_ID={}"


class ArchiveListingError(Exception):
    """Neither the JSON listing nor the archive contents page gave the study's files."""


_PAGE_SECTION_RE = re.compile(r'<strong>(ST\d+): Archive File:([^<]*)</strong></br><pre>(.*?)</pre>', re.S)
_UNZIP_HEADER_RE = re.compile(r'^\s*Length\s+Date\s+Time\s+Name\s*$')
_UNZIP_ENTRY_RE = re.compile(r'^\s*(\d+)\s+(?:\d\d-\d\d-\d{2,4}|\d{4}-\d\d-\d\d)\s+\d\d:\d\d   (.*)$')
_UNZIP_FOOTER_RE = re.compile(r'^\s*\d+\s+(\d+) files?\s*$')
_7Z_HEADER_RE = re.compile(r'^\s*Date\s+Time\s+Attr\s+Size\s+Compressed\s+Name\s*$')
_7Z_FOOTER_RE = re.compile(r'^.{19}\s+\d+\s+(?:\d+\s+)?(\d+) files?(?:, (\d+) folders?)?\s*$')


def _unzip_listing(lines):
    """[(path, size)] of the files in 'unzip -l' output (folders end in '/')."""
    start = next(i for i, l in enumerate(lines) if _UNZIP_HEADER_RE.match(l)) + 2
    if not lines[start - 1].startswith('---------'):
        raise ValueError('unzip listing without the dashes line')
    files, n_entries = [], 0
    for i in range(start, len(lines)):
        line = lines[i]
        if line.startswith('---------'):
            footer = _UNZIP_FOOTER_RE.match(lines[i + 1]) if i + 1 < len(lines) else None
            if footer is None or int(footer.group(1)) != n_entries:
                raise ValueError(f'unzip listing: {n_entries} entries, footer {lines[i + 1:i + 2]}')
            return files
        m = _UNZIP_ENTRY_RE.match(line)
        if m is None:
            raise ValueError(f'unzip listing: unexpected line {line[:200]!r}')
        n_entries += 1
        if not m.group(2).endswith('/'):
            files.append((m.group(2), m.group(1)))
    raise ValueError('unzip listing without footer')


def _7z_listing(lines):
    """[(path, size)] of the files in '7z l' output (folders have attribute D)."""
    head = next(i for i, l in enumerate(lines) if _7Z_HEADER_RE.match(l))
    dashes = lines[head + 1]
    cols = [m.span() for m in re.finditer(r'-+', dashes)]
    if len(cols) != 5 or not dashes.startswith('-------------------'):
        raise ValueError(f'7z listing: unexpected column line {dashes!r}')
    (attr_a, attr_b), (size_a, size_b), name_a = cols[1], cols[2], cols[4][0]
    files, n_dirs = [], 0
    for i in range(head + 2, len(lines)):
        line = lines[i]
        if line.startswith('-------------------'):
            footer = _7Z_FOOTER_RE.match(lines[i + 1]) if i + 1 < len(lines) else None
            if footer is None or int(footer.group(1)) != len(files) or int(footer.group(2) or 0) != n_dirs:
                raise ValueError(f'7z listing: {len(files)} files and {n_dirs} folders, footer {lines[i + 1:i + 2]}')
            return files
        attr, size, name = line[attr_a:attr_b], line[size_a:size_b].strip(), line[name_a:]
        if len(line) <= name_a or line[name_a - 2:name_a] != '  ' or not re.fullmatch(r'[D.][A-Z.]{4}', attr) \
                or not re.fullmatch(r'\d*', size):
            raise ValueError(f'7z listing: unexpected line {line[:200]!r}')
        if attr.startswith('D'):
            n_dirs += 1
        else:
            files.append((name, size or '0'))
    raise ValueError('7z listing without footer')


def listing_from_archive_page(study_id, page):
    """The files of show_archive_contents.php (one plain 'unzip -l' or '7z l' listing per
    archive) as rows like those of show_archive_contents_json.php: STUDY_ID, ARCHIVE, URL
    (A = archive, F = path), FILENAME (the path in the archive) and FILESIZE; folders are not
    listed (as in the JSON listing). Anything else - nested archives, other layouts, HTML
    entities, counts that do not match the listing's footer - raises ValueError."""
    sections = _PAGE_SECTION_RE.findall(page)
    if not sections:
        raise ValueError('no archive listing on the page')
    rows = []
    for st, archive, text in sections:
        if st != study_id:
            raise ValueError(f'listing of {st} on the page of {study_id}')
        if 'Contents of ' in text or '==> use' in text or re.search(r'&(?:[a-zA-Z]+|#\d+);', text):
            raise ValueError(f'{archive}: nested archive or escaped text, not a plain listing')
        lines = text.split('\n')
        if any(_UNZIP_HEADER_RE.match(l) for l in lines):
            files = _unzip_listing(lines)
            php = 'file_extract.php'
        elif any(_7Z_HEADER_RE.match(l) for l in lines):
            files = _7z_listing(lines)
            php = 'file_extract_7z.php'
        else:
            raise ValueError(f'{archive}: no unzip or 7z listing')
        for name, size in files:
            rows.append({'STUDY_ID': study_id, 'ARCHIVE': archive,
                         'URL': f'https://www.metabolomicsworkbench.org/data/{php}?A={quote_plus(archive)}&F={quote_plus(name)}',
                         'FILENAME': name, 'FILESIZE': size})
    return rows


def fetch_archive_listing(study_id, get=_get):
    """(rows of show_archive_contents_json.php, source). Some large studies always answer
    that request with HTTP 500; their rows then come from the listings of
    show_archive_contents.php (source 'page'). get(url) returns the response and raises on
    HTTP errors. Raises ArchiveListingError when neither gives the files."""
    try:
        return get(ARCHIVE_LISTING_URL.format(study_id)).json(), 'json'
    except KeyboardInterrupt:
        raise
    except Exception as e:
        json_error = f"{type(e).__name__}: {e}"
    try:
        rows = listing_from_archive_page(study_id, get(ARCHIVE_PAGE_URL.format(study_id)).text)
    except KeyboardInterrupt:
        raise
    except Exception as e:
        raise ArchiveListingError(f"show_archive_contents_json.php: {json_error}; "
                                  f"show_archive_contents.php: {type(e).__name__}: {e}")
    print(f"Archive listing of {study_id}: show_archive_contents_json.php failed ({json_error}); "
          f"took the {len(rows)} files of show_archive_contents.php instead")
    return rows, 'page'



def is_apple_double(path):
    # macOS resource-fork copies (__MACOSX/.../._x.mzML) are not data files
    return "__MACOSX/" in path or os.path.basename(path).startswith("._")

# Define a function to extract the first folder with an extension in order to not get the individual files in e.g. agilent "files"
def extract_first_folder_with_extension(url):
    filepath = listing_query(url).get('F', [None])[0]
    parts = filepath.split("/")

    # Check for extensions in each part
    for i, part in enumerate(parts):
        if "." in part:
            # If an extension is found, return the path up to and including that part
            return "/".join(parts[:i + 1])

    # If no extension is found, return the entire filepath
    return filepath

def _get_metabolomicsworkbench_filepaths(study_id):

    try:
        mw_file_list, _ = fetch_archive_listing(study_id)
        workbench_df = pd.DataFrame(mw_file_list)
        if len(workbench_df) == 0:
            return workbench_df
        workbench_df = workbench_df[~workbench_df['FILENAME'].apply(is_apple_double)]

        workbench_df['raw_sample_name'] = workbench_df['URL'].apply(extract_first_folder_with_extension)


        workbench_df["USI_file"] = workbench_df["URL"].apply(
            lambda url: f"mzspec:{study_id}:{listing_query(url).get('F', [None])[0]}"
        )





        workbench_df["USI_sample"] = workbench_df.apply(
            lambda
                row: f"mzspec:{study_id}:{listing_query(row['URL']).get('A', [None])[0]}-{row['raw_sample_name']}",
            axis=1
        )
        
    except KeyboardInterrupt:
        raise
    except ArchiveListingError as e:
        print(f"DROPPED {study_id}: no archive listing ({e})")
        DROPPED_STUDIES.append(study_id)
        workbench_df = pd.DataFrame()
    except Exception as e:
        print(f"Failed to list files for {study_id}: {type(e).__name__}: {e}")
        DROPPED_STUDIES.append(study_id)
        workbench_df = pd.DataFrame()

    return workbench_df


DROPPED_STUDIES = []

# Some archive listings carry raw 7z listing columns in front of the path,
# e.g. '2014-10-05 10:45:02 ..H.A           80               MX183900 .../x.d'
_7Z_LISTING_PREFIX = re.compile(r'^\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2} \S{5} +\d+ +(?:\d+ +)?')


def clean_listing_path(path):
    return _7Z_LISTING_PREFIX.sub('', path)


def process_filename(filename):
    path_parts = filename.split('/')
    for part in path_parts:
        if part.lower().endswith('.d'):
            return '/'.join(path_parts[:path_parts.index(part)+1])
    return filename

def _get_existing_datasets(path_to_file):
    try:
        existing_datasets = pd.read_csv(path_to_file, sep="\t")
        existing_datasets = set(existing_datasets['datasets'])
        return existing_datasets
    except:
        return set()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Give an MWB study ID and get a tsv with file paths present in the study.')
    parser.add_argument("--study_id", "-mwb_id", type=str, help='An MWB study ID such as "ST002050", "ALL" for every study', required=True)
    parser.add_argument("--output_path", type=str, help='Output file path to tsv file.')
    parser.add_argument("--filter_extensions", type=str, help='Filter extensions to ".mzml", ".mzxml", ".cdf", ".raw", ".wiff", and ".d".', default='False')
    parser.add_argument("--existing_datasets", type=str, help="path to a file of datasets already indexed", default="none")


    args = parser.parse_args()

    existing_datasets = _get_existing_datasets(args.existing_datasets)

    if args.study_id == "ALL":
        # Getting all files
        url = "https://www.metabolomicsworkbench.org/rest/study/study_id/ST/available"
        studies_dict = _get_json(url)

        study_list = []
        for key in studies_dict.keys():
            study_dict = studies_dict[key]
            study_list.append(study_dict['study_id'])

        study_list = sorted(set(study_list))

        all_results_list = []
        for study_id in tqdm.tqdm(study_list):
            if study_id in existing_datasets:
                print("Skipping", study_id, "Already indexed")
                continue

            temp_result_df = _get_metabolomicsworkbench_filepaths(study_id=study_id)
            all_results_list.append(temp_result_df)

        result_df = pd.concat(all_results_list, axis=0)
        if DROPPED_STUDIES:
            print(f"{len(DROPPED_STUDIES)} studies without file list (their files are not in the output): "
                  + ', '.join(DROPPED_STUDIES))

    else:
        result_df = _get_metabolomicsworkbench_filepaths(study_id=args.study_id)

    result_df['study_id'] = result_df['STUDY_ID']
    result_df['file_path'] = result_df['FILENAME'].apply(clean_listing_path).apply(process_filename)

    result_df['USI'] = result_df.apply(
        lambda row: f"mzspec:{row['study_id']}:{process_filename(clean_listing_path(listing_query(row['URL']).get('F', [''])[0]))}",
        axis=1)

    if args.filter_extensions == 'True':
        extensions = [".mzml", ".mzxml", ".cdf", ".raw", ".wiff", ".d"]
        # filter on the collapsed path so .d folders listed only by their inner files are kept
        result_df = result_df[result_df['file_path'].str.lower().str.endswith(tuple(extensions))]


    result_df = result_df[['study_id', 'file_path', 'USI']]

    # one row per file; .d folders collapse many listed files into one path
    result_df = result_df.drop_duplicates(subset=['study_id', 'file_path'], keep='first')
     

    result_df.to_csv(args.output_path, sep='\t', index=False, header=True)

    print(f"Output written to {args.output_path}")