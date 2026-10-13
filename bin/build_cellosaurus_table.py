"""Builds the compact Cellosaurus table shipped in data/ from cellosaurus.txt.

Cellosaurus (https://www.cellosaurus.org/, SIB Swiss Institute of Bioinformatics) is
distributed under CC BY 4.0. Expasy only serves the current release, so the workflow
does not download it: the table built here is pinned and shipped in data/, together
with a small JSON that records the release, its sha256 and the attribution.

usage:
    python build_cellosaurus_table.py cellosaurus.txt data/ [--ncbi_nodes_dmp nodes.dmp]

writes data/cellosaurus_<version>_cell_lines.tsv.gz (one row per name and synonym of a
cell line) and data/cellosaurus_<version>.json. Columns:
    CellLineIndex  primary accession (CVCL_xxxx)
    CellLine       recommended name (ID line)
    Synonym        the name itself or one of its synonyms (SY line)
    NCBITaxonIDs   species of origin ('|'-joined NCBI taxon IDs, OX lines)
    Species        their names as given on the OX lines ('|'-joined, common names removed)
    Category       CA line (Cancer cell line, Hybridoma, ...)
    BTO            BRENDA tissue ontology cross-references ('|'-joined, DR BTO lines)
    Problematic    first sentence of 'CC   Problematic cell line', '' if none
    Mammalian      1 if every species of origin is a mammal (NCBI Mammalia), 0 if not,
                   '' without --ncbi_nodes_dmp
"""
import argparse
import csv
import gzip
import hashlib
import io
import json
import os
import re
import tempfile

COLUMNS = ['CellLineIndex', 'CellLine', 'Synonym', 'NCBITaxonIDs', 'Species', 'Category', 'BTO', 'Problematic', 'Mammalian']
MAMMALIA = 40674


def parse_cellosaurus(path):
    """(header dict, list of entries) of a cellosaurus.txt file."""
    header = {}
    entries = []
    entry = None
    with open(path, encoding='utf-8') as fh:
        for line in fh:
            if entry is None and not line.startswith('ID   '):
                m = re.match(r'\s*(Version|Last update):\s*(.+?)\s*$', line)
                if m and m.group(1) not in header:
                    header[m.group(1)] = m.group(2)
                continue
            tag, value = line[:2], line[5:].rstrip('\n')
            if tag == 'ID':
                entry = {'ID': value.strip(), 'AC': '', 'SY': [], 'OX': [], 'species': [], 'CA': '', 'BTO': [], 'problem': ''}
            elif tag == 'AC':
                entry['AC'] = value.strip()
            elif tag == 'SY':
                entry['SY'] += [s.strip() for s in value.split('; ') if s.strip()]
            elif tag == 'OX':
                m = re.match(r'NCBI_TaxID=(\d+);(?:\s*!\s*(.*))?', value)
                if m and m.group(1) not in entry['OX']:
                    entry['OX'].append(m.group(1))
                    entry['species'].append(re.sub(r'\s*\([^()]*\)\s*$', '', (m.group(2) or '').strip()))
            elif tag == 'CA':
                entry['CA'] = value.strip()
            elif tag == 'DR' and value.startswith('BTO; '):
                bto = value[5:].strip()
                if bto not in entry['BTO']:
                    entry['BTO'].append(bto)
            elif tag == 'CC' and value.startswith('Problematic cell line:') and not entry['problem']:
                entry['problem'] = value.split(':', 1)[1].strip().split('. ')[0].rstrip('.')
            elif tag == '//':
                entries.append(entry)
                entry = None
    return header, entries


def mammal_flags(taxids, nodes_dmp):
    """taxon ID -> True if it is in NCBI Mammalia."""
    wanted = set(int(t) for t in taxids)
    parent = {}
    with open(nodes_dmp, encoding='utf-8') as fh:
        for line in fh:
            parts = line.split('\t|\t', 2)
            parent[int(parts[0])] = int(parts[1])
    flags = {}
    for t in wanted:
        node, seen = t, set()
        while node in parent and node not in seen and node != MAMMALIA:
            seen.add(node)
            if parent[node] == node:
                break
            node = parent[node]
        flags[str(t)] = node == MAMMALIA
    return flags


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


def build(src, out_dir, nodes_dmp=None):
    header, entries = parse_cellosaurus(src)
    version = header.get('Version')
    if not version or not re.fullmatch(r'[0-9.]+', version):
        raise SystemExit(f'No release version in the header of {src}')
    flags = mammal_flags({t for e in entries for t in e['OX']}, nodes_dmp) if nodes_dmp else None

    buf = io.StringIO()
    writer = csv.writer(buf, delimiter='\t', lineterminator='\n', quoting=csv.QUOTE_MINIMAL)
    writer.writerow(COLUMNS)
    n_rows = 0
    for e in entries:
        if not re.fullmatch(r'CVCL_[0-9A-Z]{4}', e['AC']):
            raise SystemExit(f'Unexpected accession {e["AC"]!r} for {e["ID"]!r}')
        mammal = '' if flags is None else ('1' if e['OX'] and all(flags[t] for t in e['OX']) else '0')
        names = [e['ID']] + [s for s in dict.fromkeys(e['SY']) if s != e['ID']]
        for name in names:
            writer.writerow([e['AC'], e['ID'], name, '|'.join(e['OX']), '|'.join(e['species']), e['CA'], '|'.join(e['BTO']),
                             e['problem'], mammal])
            n_rows += 1
    # mtime=0: the same release always gives the same bytes
    data = gzip.compress(buf.getvalue().encode('utf-8'), compresslevel=9, mtime=0)
    table = os.path.join(out_dir, f'cellosaurus_{version}_cell_lines.tsv.gz')
    _atomic_bytes(table, data)

    sha = hashlib.sha256()
    with open(src, 'rb') as fh:
        for block in iter(lambda: fh.read(1 << 20), b''):
            sha.update(block)
    meta = {
        'resource': 'Cellosaurus',
        'version': version,
        'release_date': header.get('Last update', ''),
        'source_url': 'https://ftp.expasy.org/databases/cellosaurus/cellosaurus.txt',
        'source_sha256': sha.hexdigest(),
        'source_bytes': os.path.getsize(src),
        'table': os.path.basename(table),
        'table_sha256': hashlib.sha256(data).hexdigest(),
        'cell_lines': len(entries),
        'rows': n_rows,
        'columns': COLUMNS,
        'license': 'CC BY 4.0 (https://creativecommons.org/licenses/by/4.0/)',
        'attribution': ('Cellosaurus, SIB Swiss Institute of Bioinformatics (https://www.cellosaurus.org/). '
                        'Bairoch A. The Cellosaurus, a cell-line knowledge resource. J Biomol Tech. '
                        '2018;29(2):25-38. doi:10.7171/jbt.18-2902-002'),
        'changes': ('Reduced to accession, recommended name, names/synonyms, species of origin (IDs and names), category, '
                    'BTO cross-references and the problematic-cell-line note; Mammalian derived from NCBI Taxonomy.'),
        'build': 'bin/build_cellosaurus_table.py',
    }
    _atomic_bytes(os.path.join(out_dir, f'cellosaurus_{version}.json'),
                  (json.dumps(meta, indent=4, ensure_ascii=False) + '\n').encode('utf-8'))
    print(f'Cellosaurus {version}: {len(entries)} cell lines, {n_rows} rows, {len(data) / 1e6:.1f} MB -> {table}')
    return table


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Compact Cellosaurus cell-line table for ReDU')
    parser.add_argument('cellosaurus_txt')
    parser.add_argument('output_folder')
    parser.add_argument('--ncbi_nodes_dmp', default=None, help='NCBI nodes.dmp; adds the Mammalian column')
    args = parser.parse_args()
    build(args.cellosaurus_txt, args.output_folder, args.ncbi_nodes_dmp)
