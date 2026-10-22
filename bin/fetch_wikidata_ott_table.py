"""Fetches Wikidata's NCBI taxon ID -> Open Tree of Life (OTT) ID mapping, the table OpenTreeOfLifeTaxonomyID comes
from (ncbi_database.py), and writes the pinned copy shipped in data/. Not run by the workflow: Wikidata is edited
every day, so a live query would change the column without any change in the metadata; the workflow reads the
pinned copy (WIKIDATA_OTT_TABLE in nf_workflow.nf), and this script is only run to move to a newer copy.

Wikidata (https://www.wikidata.org/) is CC0 1.0 (public domain dedication). The query (QUERY: taxa, P31 Q16521,
with their NCBI taxonomy ID, P685, and Open Tree of Life ID, P9157, if any) goes to the QLever Wikidata endpoint;
the answer is kept as it is (TSV with header '?ncbi_id\\t?ott_id', values in double quotes), only gzipped.

usage:
    python fetch_wikidata_ott_table.py data/                       fetch now (date of today, UTC)
    python fetch_wikidata_ott_table.py data/ --from_tsv answer.tsv --fetch_date 2026-10-21
                                                                   an answer fetched earlier (e.g. a run's sparql_oto_ncbi.tsv)

writes data/wikidata_ncbi_ott_<date>.tsv.gz (gzip -9 without name and time stamp: the same answer always gives the
same bytes) and data/wikidata_ncbi_ott_<date>.json (endpoint, query, fetch date, sha256 of the answer and of the
.gz, rows, license). To move to it: test the workflow with it, then change WIKIDATA_OTT_TABLE in nf_workflow.nf
and remove the old pair of files.
"""
import argparse
import datetime
import gzip
import hashlib
import json
import os
import tempfile
import urllib.request

ENDPOINT = 'https://qlever.dev/api/wikidata'
QUERY = ('PREFIX wd: <http://www.wikidata.org/entity/> PREFIX wdt: <http://www.wikidata.org/prop/direct/> '
         'PREFIX bd: <http://www.bigdata.com/rdf#> PREFIX wikibase: <http://wikiba.se/ontology#> '
         'SELECT ?ncbi_id ?ott_id WHERE { ?item wdt:P685 ?ncbi_id . OPTIONAL { ?item wdt:P9157 ?ott_id } '
         '?item wdt:P31 wd:Q16521 }')
HEADER = b'?ncbi_id\t?ott_id'
# far fewer rows than the ~620,000 of 2026-10 means a truncated or failed answer
MIN_ROWS = 300000
USER_AGENT = 'ReDU-metadata-workflow (https://github.com/Wang-Bioinformatics-Lab/PublicDataset_ReDU_Metadata_Workflow)'


def fetch(timeout=600):
    """The endpoint's answer to QUERY (bytes)"""
    request = urllib.request.Request(ENDPOINT, data=QUERY.encode('utf-8'), method='POST', headers={
        'Accept': 'text/tab-separated-values', 'Content-type': 'application/sparql-query', 'User-Agent': USER_AGENT})
    with urllib.request.urlopen(request, timeout=timeout) as response:
        return response.read()


def check(data):
    """Number of rows of an answer; stops on an answer that is not the expected table"""
    lines = data.split(b'\n')
    if lines[0].rstrip(b'\r') != HEADER:
        raise SystemExit(f'Unexpected header {lines[0][:200]!r} (expected {HEADER!r})')
    if not data.endswith(b'\n'):
        raise SystemExit('The answer does not end with a complete line')
    rows = len(lines) - 2
    if rows < MIN_ROWS:
        raise SystemExit(f'Only {rows} rows (expected at least {MIN_ROWS})')
    return rows


def _atomic_bytes(path, data):
    fd, tmp = tempfile.mkstemp(dir=os.path.dirname(os.path.abspath(path)), prefix='.' + os.path.basename(path) + '.')
    try:
        with os.fdopen(fd, 'wb') as fh:
            fh.write(data)
            fh.flush()
            os.fsync(fh.fileno())
        os.chmod(tmp, 0o644)
        os.replace(tmp, path)
    except BaseException:
        if os.path.exists(tmp):
            os.unlink(tmp)
        raise


def build(out_dir, data, fetch_date):
    rows = check(data)
    # mtime=0: the same answer always gives the same bytes
    packed = gzip.compress(data, compresslevel=9, mtime=0)
    table = os.path.join(out_dir, f'wikidata_ncbi_ott_{fetch_date}.tsv.gz')
    _atomic_bytes(table, packed)
    meta = {
        'resource': 'Wikidata: NCBI taxonomy ID (P685) -> Open Tree of Life ID (P9157) of taxa (P31 Q16521)',
        'endpoint': ENDPOINT,
        'query': QUERY,
        'fetch_date': fetch_date,
        'answer_sha256': hashlib.sha256(data).hexdigest(),
        'answer_bytes': len(data),
        'table': os.path.basename(table),
        'table_sha256': hashlib.sha256(packed).hexdigest(),
        'rows': rows,
        'columns': ['?ncbi_id', '?ott_id'],
        'license': 'CC0 1.0 (https://creativecommons.org/publicdomain/zero/1.0/); Wikidata, https://www.wikidata.org/',
        'used_by': 'bin/ncbi_database.py (OpenTreeOfLifeTaxonomyID); WIKIDATA_OTT_TABLE in nf_workflow.nf',
        'build': 'bin/fetch_wikidata_ott_table.py (the .gz holds the answer as it is: gzip -9, no name or time stamp)',
        'bump': ('python bin/fetch_wikidata_ott_table.py data/; test the workflow with the new table; change '
                 'WIKIDATA_OTT_TABLE in nf_workflow.nf; remove the old table and its json'),
    }
    _atomic_bytes(os.path.join(out_dir, f'wikidata_ncbi_ott_{fetch_date}.json'),
                  (json.dumps(meta, indent=4, ensure_ascii=False) + '\n').encode('utf-8'))
    print(f'Wikidata NCBI -> OTT ({fetch_date}): {rows} rows, {len(packed) / 1e6:.1f} MB -> {table}')
    return table


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Pinned copy of the Wikidata NCBI -> OTT table for ReDU')
    parser.add_argument('output_folder')
    parser.add_argument('--from_tsv', default=None, help='an answer of the endpoint fetched earlier, instead of fetching now')
    parser.add_argument('--fetch_date', default=None, help='YYYY-MM-DD the answer was fetched (default: today, UTC)')
    args = parser.parse_args()
    date = args.fetch_date or datetime.datetime.now(datetime.timezone.utc).strftime('%Y-%m-%d')
    datetime.date.fromisoformat(date)
    if args.from_tsv:
        with open(args.from_tsv, 'rb') as fh:
            answer = fh.read()
    else:
        answer = fetch()
    build(args.output_folder, answer, date)
