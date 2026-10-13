"""Label/synonym tables from OBO files (HPO, ECTO, PECO).

These ontologies are read from their .obo release instead of OWL: the OWL releases of
ECTO and HPO import large parts of other ontologies, while the .obo files hold exactly
the native terms with their labels and exact synonyms.

    obo_table(path, prefix, index_column, root=None) -> DataFrame[index_column, Label, Synonym]

Obsolete terms are dropped. With `root`, only the root and its is_a descendants are kept
(e.g. HP:0000118 'Phenotypic abnormality' excludes inheritance modes, onset, frequency).
A missing or empty file gives an empty table, so these ontologies stay optional.
"""
import os
import re

import pandas as pd


def _parse_obo(path):
    terms = {}
    current = None
    in_term = False
    with open(path, encoding='utf-8') as handle:
        for line in handle:
            line = line.rstrip('\n')
            if line.startswith('['):
                in_term = line == '[Term]'
                current = None
                continue
            if not in_term:
                continue
            if line.startswith('id: '):
                current = {'id': line[4:].strip(), 'name': None, 'synonyms': [], 'parents': [], 'obsolete': False}
                terms[current['id']] = current
            elif current is None:
                continue
            elif line.startswith('name: '):
                current['name'] = line[6:].strip()
            elif line.startswith('synonym: '):
                match = re.match(r'synonym: "((?:[^"\\]|\\.)*)" EXACT', line)
                if match:
                    current['synonyms'].append(match.group(1).replace('\\"', '"'))
            elif line.startswith('is_a: '):
                current['parents'].append(line[6:].split('!')[0].strip())
            elif line.startswith('is_obsolete: true'):
                current['obsolete'] = True
    return terms


def obo_table(path, prefix, index_column, root=None):
    columns = [index_column, 'Label', 'Synonym']
    if not path or path == 'none' or not os.path.exists(path) or os.path.getsize(path) == 0:
        return pd.DataFrame(columns=columns)

    terms = _parse_obo(path)
    keep = {t for t, info in terms.items() if t.startswith(prefix) and not info['obsolete'] and info['name']}

    if root:
        children = {}
        for t, info in terms.items():
            for parent in info['parents']:
                children.setdefault(parent, set()).add(t)
        descendants, stack = set(), [root]
        while stack:
            node = stack.pop()
            if node in descendants:
                continue
            descendants.add(node)
            stack.extend(children.get(node, ()))
        keep &= descendants

    rows = []
    for t in sorted(keep):
        info = terms[t]
        for synonym in [info['name']] + [s for s in info['synonyms'] if s != info['name']]:
            rows.append({index_column: t, 'Label': info['name'], 'Synonym': synonym})
    return pd.DataFrame(rows, columns=columns)


def hpo_table(path):
    return obo_table(path, 'HP:', 'HPOOntologyIndex', root='HP:0000118')


def exposure_table(ecto_path, peco_path):
    # ECTO first: for the few labels both ontologies use, the general ECTO term is kept
    return pd.concat([obo_table(ecto_path, 'ECTO:', 'ExposureIndex'),
                      obo_table(peco_path, 'PECO:', 'ExposureIndex')], ignore_index=True)
