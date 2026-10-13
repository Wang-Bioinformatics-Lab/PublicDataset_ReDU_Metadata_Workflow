"""Sets up the reference data (ontologies, NCBI taxonomy, MASST tables) in data/.

Runs as the first step of the workflow, so no manual data/get_data.sh run is needed.
What is needed is listed in data/reference_data.json; what is present is recorded in
data/.reference_data_state.json (url and server fingerprint of every file).

- 'pinned' files are downloaded only when missing or when their URL in the list changed.
- 'latest' files are downloaded when the server reports a new version (ETag,
  Last-Modified or size differ from the recorded ones); if the server reports none of
  these, they are refreshed once they are older than --max-age-days.
- Files present without a record (e.g. from an earlier manual get_data.sh run) are
  replaced once, so every file comes from a known source; leftovers of that setup
  (taxdmp.zip, *.part) are removed.
- Downloads go to <name>.part and replace the file only after they pass a format check.
  If a download fails, an existing file is kept; a missing optional file becomes empty
  (the column it feeds stays empty); a missing required file stops the workflow.

With --link_dir, symlinks to all resulting files are created there (the workflow
passes them on as process outputs).

    python fetch_reference_data.py data/reference_data.json data/ [--link_dir .]
"""
import argparse
import datetime
import fcntl
import json
import os
import shutil
import sys
import time
import urllib.request
import zipfile

STATE_FILE = '.reference_data_state.json'
LEFTOVERS = ['taxdmp.zip']
USER_AGENT = 'ReDU-metadata-workflow (https://github.com/Wang-Bioinformatics-Lab/PublicDataset_ReDU_Metadata_Workflow)'


def log(message):
    print(f"[reference data] {message}", flush=True)


def now_utc():
    return datetime.datetime.utcnow().replace(microsecond=0).isoformat() + 'Z'


def request(url, method='GET', timeout=120):
    req = urllib.request.Request(url, method=method, headers={'User-Agent': USER_AGENT})
    return urllib.request.urlopen(req, timeout=timeout)


def remote_fingerprint(url):
    """ETag / Last-Modified / Content-Length of the final URL, or None if unavailable."""
    try:
        with request(url, method='HEAD', timeout=60) as response:
            headers = response.headers
            fingerprint = {key: headers.get(key) for key in ('ETag', 'Last-Modified', 'Content-Length') if headers.get(key)}
            return fingerprint or None
    except Exception as e:
        log(f"could not check {url}: {type(e).__name__}: {e}")
        return None


def download(url, path, attempts=3):
    part = path + '.part'
    for attempt in range(1, attempts + 1):
        try:
            with request(url, timeout=600) as response, open(part, 'wb') as handle:
                expected = response.headers.get('Content-Length')
                shutil.copyfileobj(response, handle, length=1 << 20)
            size = os.path.getsize(part)
            if expected and int(expected) != size:
                raise IOError(f"incomplete download: {size} of {expected} bytes")
            return part
        except Exception as e:
            log(f"download of {url} failed (attempt {attempt}/{attempts}): {type(e).__name__}: {e}")
            if os.path.exists(part):
                os.remove(part)
            time.sleep(10 * attempt)
    return None


def looks_valid(path, check):
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        return False
    if check == 'zip':
        try:
            with zipfile.ZipFile(path) as archive:
                return archive.testzip() is None
        except zipfile.BadZipFile:
            return False
    with open(path, 'rb') as handle:
        head = handle.read(65536)
    if check == 'owl':
        return b'<rdf:RDF' in head or b'<Ontology' in head or b'owl:Ontology' in head
    if check == 'obo':
        return b'format-version:' in head or b'[Term]' in head
    # text: not an HTML error page
    return not head.lstrip().lower().startswith((b'<!doctype html', b'<html'))


def outputs_of(entry):
    return entry.get('extract') or [entry['name']]


def install(entry, part, data_dir):
    """Moves a checked download into place (or extracts the listed members)."""
    if entry.get('extract'):
        with zipfile.ZipFile(part) as archive:
            for member in entry['extract']:
                target = os.path.join(data_dir, member)
                with archive.open(member) as source, open(target + '.part', 'wb') as handle:
                    shutil.copyfileobj(source, handle, length=1 << 20)
                os.replace(target + '.part', target)
        os.remove(part)
    else:
        os.replace(part, os.path.join(data_dir, entry['name']))


def outputs_present(entry, data_dir):
    if entry.get('extract'):
        return all(os.path.exists(os.path.join(data_dir, m)) and os.path.getsize(os.path.join(data_dir, m)) > 0
                   for m in entry['extract'])
    return looks_valid(os.path.join(data_dir, entry['name']), entry['check'])


def needs_update(entry, record, data_dir, max_age_days):
    """Returns (reason or None, fingerprint)."""
    if not outputs_present(entry, data_dir):
        return 'missing', None
    if not record:
        return 'no record of where the present file came from', None
    if record.get('url') != entry['url']:
        return 'version changed in reference_data.json', None
    if entry['mode'] == 'pinned':
        return None, record.get('fingerprint')
    fingerprint = remote_fingerprint(entry['url'])
    if fingerprint is None:
        fetched = record.get('fetched_at', '')
        try:
            age = datetime.datetime.utcnow() - datetime.datetime.strptime(fetched, '%Y-%m-%dT%H:%M:%SZ')
        except ValueError:
            age = datetime.timedelta(days=max_age_days + 1)
        return ('older than %d days' % max_age_days if age.days > max_age_days else None), record.get('fingerprint')
    if fingerprint != record.get('fingerprint'):
        return 'new version on the server', fingerprint
    return None, fingerprint


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('manifest')
    parser.add_argument('data_dir')
    parser.add_argument('--link_dir', default=None)
    parser.add_argument('--max_age_days', type=int, default=7)
    args = parser.parse_args()

    with open(args.manifest) as handle:
        entries = json.load(handle)['files']
    data_dir = os.path.abspath(args.data_dir)
    state_path = os.path.join(data_dir, STATE_FILE)

    # one setup at a time (e.g. a manual run next to the workflow)
    lock = open(os.path.join(data_dir, '.reference_data.lock'), 'w')
    fcntl.flock(lock, fcntl.LOCK_EX)

    state = {}
    if os.path.exists(state_path):
        with open(state_path) as handle:
            state = json.load(handle)

    # leftovers of the earlier manual setup or of interrupted downloads
    for name in os.listdir(data_dir):
        if name.endswith('.part') or (name in LEFTOVERS and name not in state):
            log(f"removing leftover {name}")
            os.remove(os.path.join(data_dir, name))

    failed_required = []
    for entry in entries:
        name = entry['name']
        reason, fingerprint = needs_update(entry, state.get(name), data_dir, args.max_age_days)
        if reason is None:
            log(f"{name}: up to date")
            continue
        log(f"{name}: downloading ({reason})")
        if fingerprint is None:
            fingerprint = remote_fingerprint(entry['url'])
        part = download(entry['url'], os.path.join(data_dir, name))
        if part and looks_valid(part, entry['check']):
            install(entry, part, data_dir)
            state[name] = {'url': entry['url'], 'fingerprint': fingerprint, 'fetched_at': now_utc(),
                           'files': outputs_of(entry)}
            with open(state_path + '.part', 'w') as handle:
                json.dump(state, handle, indent=2)
            os.replace(state_path + '.part', state_path)
            log(f"{name}: installed")
            continue
        if part and os.path.exists(part):
            log(f"{name}: download did not pass the format check, discarded")
            os.remove(part)
        if outputs_present(entry, data_dir):
            log(f"WARNING {name}: keeping the existing file")
        elif not entry.get('required', True):
            log(f"WARNING {name}: not available; the column it feeds stays empty")
            for output in outputs_of(entry):
                open(os.path.join(data_dir, output), 'a').close()
        else:
            failed_required.append(name)

    if failed_required:
        log(f"ERROR: required reference data missing and could not be downloaded: {', '.join(failed_required)}")
        sys.exit(1)

    if args.link_dir:
        for entry in entries:
            for output in outputs_of(entry):
                link = os.path.join(args.link_dir, output)
                if os.path.lexists(link):
                    os.remove(link)
                os.symlink(os.path.join(data_dir, output), link)


if __name__ == '__main__':
    main()
