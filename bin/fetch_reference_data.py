"""Sets up the reference data (ontologies, NCBI taxonomy, MASST tables) in data/.

Runs as the first step of the workflow, so no manual data/get_data.sh run is needed.
What is needed is listed in data/reference_data.json: every file is pinned to the exact
version the workflow was tested with (url + sha256). Newer upstream releases are never
picked up on their own; moving to a new version means testing it and then changing the
pin. What is present is recorded in data/.reference_data_state.json.

- Entries with "refresh": "every_run" (the MASST tables, which are data that grow
  continuously) are downloaded on every run instead; they have no checksum, and the
  previous copy is kept if a download fails.
- A file is downloaded only when it is missing or is not the pinned version (different
  url/sha256 in the record, or the file on disk no longer matches its checksum).
- Files present without a record (e.g. from an earlier manual get_data.sh run) are
  replaced once by the pinned version; leftovers of that setup (taxdmp.zip, *.part)
  are removed.
- Downloads go to <name>.part and replace the file only if their sha256 equals the pin
  and they pass a format check. If that fails, an existing file is kept (with a
  warning); a missing optional file becomes empty (the column it feeds stays empty); a
  missing required file stops the workflow.

- A name with a folder (po_imports/ro_import.owl) is installed in that folder of data/:
  the owl:imports of an ontology, which bin/owl_loading.py reads from there.

With --link_dir, symlinks to all resulting files (and folders) are created there (the
workflow passes them on as process outputs).

    python fetch_reference_data.py data/reference_data.json data/ [--link_dir .]
"""
import argparse
import datetime
import fcntl
import hashlib
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


def sha256_of(path):
    digest = hashlib.sha256()
    with open(path, 'rb') as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b''):
            digest.update(chunk)
    return digest.hexdigest()


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


def sha256_of_outputs(entry, data_dir):
    return sha256_of(os.path.join(data_dir, outputs_of(entry)[0]))


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


def needs_update(entry, record, data_dir):
    """Reason why the file has to be (re)downloaded, or None if the pinned version is present."""
    if not outputs_present(entry, data_dir):
        return 'missing'
    if entry.get('refresh') == 'every_run':
        return 'refreshed on every run'
    if not record:
        return 'no record of which version the present file is'
    if record.get('url') != entry['url'] or record.get('sha256') != entry['sha256']:
        return 'pinned version changed in reference_data.json'
    if entry.get('extract'):
        sizes = record.get('sizes', {})
        if any(os.path.getsize(os.path.join(data_dir, m)) != sizes.get(m) for m in entry['extract']):
            return 'extracted files differ from the pinned version'
    elif sha256_of(os.path.join(data_dir, entry['name'])) != entry['sha256']:
        return 'file on disk differs from the pinned version'
    return None


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('manifest')
    parser.add_argument('data_dir')
    parser.add_argument('--link_dir', default=None)
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
    for folder, _, names in os.walk(data_dir):
        for name in names:
            if name.endswith('.part') or (folder == data_dir and name in LEFTOVERS and name not in state):
                log(f"removing leftover {os.path.relpath(os.path.join(folder, name), data_dir)}")
                os.remove(os.path.join(folder, name))

    failed_required = []
    for entry in entries:
        name = entry['name']
        reason = needs_update(entry, state.get(name), data_dir)
        if reason is None:
            log(f"{name}: pinned version present")
            continue
        log(f"{name}: downloading ({reason})")
        os.makedirs(os.path.dirname(os.path.join(data_dir, name)), exist_ok=True)
        part = download(entry['url'], os.path.join(data_dir, name))
        if part and entry.get('sha256') and sha256_of(part) != entry['sha256']:
            log(f"{name}: download does not match the pinned sha256 (the source changed?), discarded")
            os.remove(part)
            part = None
        if part and looks_valid(part, entry['check']):
            install(entry, part, data_dir)
            state[name] = {'url': entry['url'], 'sha256': entry.get('sha256') or sha256_of_outputs(entry, data_dir), 'fetched_at': now_utc(),
                           'files': outputs_of(entry),
                           'sizes': {m: os.path.getsize(os.path.join(data_dir, m)) for m in outputs_of(entry)}}
            with open(state_path + '.part', 'w') as handle:
                json.dump(state, handle, indent=2)
            os.replace(state_path + '.part', state_path)
            log(f"{name}: installed")
            continue
        if part and os.path.exists(part):
            log(f"{name}: download did not pass the format check, discarded")
            os.remove(part)
        if outputs_present(entry, data_dir):
            what = 'a new copy' if entry.get('refresh') == 'every_run' else 'the pinned version'
            log(f"WARNING {name}: could not get {what}, keeping the existing file")
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
        # one link per file, or per folder (po_imports/...)
        for output in dict.fromkeys(output.split('/')[0] for entry in entries for output in outputs_of(entry)):
            link = os.path.join(args.link_dir, output)
            if os.path.lexists(link):
                os.remove(link)
            os.symlink(os.path.join(data_dir, output), link)


if __name__ == '__main__':
    main()
