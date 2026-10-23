import json
import os
import re
import requests
import time
import functools
from datetime import datetime
from bs4 import BeautifulSoup
from owlready2 import get_ontology
import owlready2
import owlready2.namespace as _owl_ns
import pandas as pd
import tqdm

from blank_qc_names import blank_qc_type_from_filename
from owl_loading import load_ontology
# several organisms in one NCBITaxonomy value ('ID|Name||ID|Name')
from ncbi_database import is_multi, single_taxid

# owlready2 raises TypeError when an OWL file (e.g. new cl.owl via STATO import)
# declares an IRI as both a property and a class/individual. Patch _load_properties
# to skip only those conflicting entries instead of crashing. (Catching the TypeError
# around the original method would also skip loading every property after the
# conflicting one.) Mirrors Ontology._load_properties of owlready2 0.46.
def _patched_load_properties(self):
    # Update props from other ontologies, if needed
    for prop in list(self.world._props.values()):
        if prop.namespace.world is _owl_ns.owl_world: continue
        prop._check_update(self)

    # Loads new props
    for prop_storid in _owl_ns.itertools.chain(self._get_obj_triples_po_s(_owl_ns.rdf_type, _owl_ns.owl_object_property),
                                               self._get_obj_triples_po_s(_owl_ns.rdf_type, _owl_ns.owl_data_property),
                                               self._get_obj_triples_po_s(_owl_ns.rdf_type, _owl_ns.owl_annotation_property)):
        Prop = self.world._get_by_storid(prop_storid)
        if not isinstance(Prop, _owl_ns.PropertyClass):
            print(f"Skipping '{Prop.iri}': it is both a property and a class/an individual")
            continue
        python_name_d = self.world._get_data_triple_sp_od(prop_storid, _owl_ns.owlready_python_name)
        if python_name_d is not None:
            with _owl_ns.LOADING: Prop.python_name = python_name_d[0]
_owl_ns.Ontology._load_properties = _patched_load_properties



def _is_missing_observation(value):
    return value is None or (isinstance(value, float) and pd.isna(value)) or str(value).strip().lower() in ('', 'nan', 'missing value')


def _blank_qc_sample_type(value):
    """a blank or pooled QC SampleType (blank_analysis, blank_extraction, blank_QC, blank_culturemedia, pool_QC)"""
    return not _is_missing_observation(value) and ('blank' in str(value).lower() or str(value).strip() == 'pool_QC')


def merge_repeated_fileobservations(df):
    """Collapses rows describing the same file into one.

    A file listed for several samples or assays cannot be attributed to one of them, so
    columns with conflicting values become 'missing value'. Blank/missing entries are not
    a conflict: if exactly one real value is present it is kept (e.g. a sample listed twice,
    once with 'France' and once with an empty country).
    A file listed both for a blank / pooled QC and for samples of another SampleType whose name
    marks it as a blank / QC (blank_qc_names.py; MTBLS10288 HF5_Neg_QC4.raw: a pooled QC and five
    experimental samples) is the blank / QC: only its blank / QC rows are merged, so it keeps that
    SampleType and, as any blank / QC, none of the samples' subject fields
    (complete_and_fill_REDU_table clears them). Other such files (batch files holding samples,
    blanks and QCs: MTBLS3839 Batch1_MRM3_June22_2020.wiff) keep the merge's 'missing value'.
    """
    if 'SampleType' in df.columns and 'filename' in df.columns and len(df):
        blank_qc = df['SampleType'].map(_blank_qc_sample_type)
        other = ~blank_qc & ~df['SampleType'].map(_is_missing_observation)
        mixed = df['filename'].isin(set(df.loc[blank_qc, 'filename'])) & df['filename'].isin(set(df.loc[other, 'filename']))
        if mixed.any():
            named = df['filename'].map(lambda f: isinstance(f, str) and blank_qc_type_from_filename(f) is not None)
            if (mixed & ~named).any():
                study = df['MassiveID'].iloc[0] if 'MassiveID' in df.columns else ''
                print(f"[repeated file] {study}: {df.loc[mixed & ~named, 'filename'].nunique()} files listed for blanks / QCs "
                      f"and for samples are not named as blanks / QCs: merged as they are")
            mixed &= named
        if mixed.any():
            for name, group in df[mixed].groupby('filename', sort=True):
                study = group['MassiveID'].iloc[0] if 'MassiveID' in group.columns else ''
                print(f"[repeated file] {study} {name}: listed as {dict(group['SampleType'].astype(str).value_counts().sort_index())}; "
                      f"named as a blank / QC: its blank / QC rows are kept")
            df = df[~mixed | blank_qc]

    def process_group(group):

        for col in group.columns:
            present = [value for value in group[col] if not _is_missing_observation(value)]
            if len({str(value) for value in present}) == 1:
                group[col] = present[0]
            elif present or len(set(group[col])) > 1 or all(value is None for value in group[col]):
                group[col] = 'missing value'
            else:
                # If there's no discrepancy, keep the original value
                group[col] = group[col].iloc[0]
        
        return group

    # Group by 'USI' and apply the processing function to each group
    df = df.groupby('filename').apply(process_group)

    # Remove duplicates based on 'filename' column
    df = df.drop_duplicates(subset='filename')

    return df

def find_column_after_target_column(df, target_column='', search_column_prefix='Samples_Unit'):
    if target_column in df.columns:
        target_index = df.columns.get_loc(target_column)
        # Check the next 3 columns after the target column
        for i in range(1, 4):
            if target_index + i < len(df.columns):
                next_column = df.columns[target_index + i]
                # Check if the next column starts with the search_column_prefix
                if next_column.startswith(search_column_prefix):
                    return next_column
    return ''  # Return None if no matching column is found

def age_category(age):
    if age is None:
        return ''
    try:
        age = float(age)
    except ValueError:
        return ''
    if age < 2:
        return 'Infancy (<2 yrs)'
    elif age <= 8:
        return 'Early Childhood (2 yrs < x <=8 yrs)'
    elif age <= 18:
        return 'Adolescence (8 yrs < x <= 18 yrs)'
    elif age <= 45:
        return 'Early Adulthood (18 yrs < x <= 45 yrs)'
    elif age <= 65:
        return 'Middle Adulthood (45 yrs < x <= 65 yrs)'
    elif age > 65:
        return 'Later Adulthood (>65 yrs)'
    else:
        return ''
    

# ---- NCBI lineage --------------------------------------------------------------------
# The lineage of a taxon comes from the NCBI taxonomy pinned in data/reference_data.json
# (nodes.dmp, names.dmp and merged.dmp of the taxdmp archive) once a converter has called
# use_pinned_ncbi_taxonomy() (MWB_to_REDU.py does, from its --path_to_ncbi_*_dmp
# arguments). It is the lineage NCBI's E-utilities efetch returns for that release: the
# scientific names of all ancestors, top-down, without the root and without the taxon itself
# ('cellular organisms; Eukaryota; ...; Homo' for 9606); merged IDs get the lineage of the
# taxon they were merged into, IDs neither in nodes.dmp nor in merged.dmp have none.
#
# Converters that have not called use_pinned_ncbi_taxonomy() still ask E-utilities (with a
# warning in the log); a missing lineage is never cached there, so a failed request is
# retried for the next study instead of being remembered as 'no lineage'.

_PINNED_LINEAGE = None


class PinnedNcbiLineage:
    """efetch-style lineages from nodes.dmp / names.dmp / merged.dmp of the pinned taxdump."""

    def __init__(self, nodes_dmp, names_dmp, merged_dmp):
        import numpy as np
        from ncbi_database import Taxonomy, ROOT
        t0 = time.time()
        self.taxonomy = Taxonomy(nodes_dmp, names_dmp, merged_dmp)
        self.root = ROOT
        # a lineage only holds ancestors, i.e. taxa that are the parent of another taxon
        parents = np.frombuffer(self.taxonomy.parent, dtype=np.int32)
        internal = set(np.unique(parents[parents >= 0]).tolist())
        self.names = self.taxonomy.scientific_names(internal)
        unnamed = sorted(internal - set(self.names) - {self.root})
        if unnamed:
            raise ValueError(f"{names_dmp}: no scientific name for {len(unnamed)} taxa of {nodes_dmp}, "
                             f"e.g. {unnamed[:10]}")
        self.lineages = {}           # resolved taxon ID -> lineage (only found ones are kept)
        self.reported = set()        # IDs already reported as having no lineage
        print(f"NCBI lineage from the pinned taxonomy: {self.taxonomy.n_taxa} taxa in {nodes_dmp}, "
              f"{len(self.taxonomy.merged)} merged IDs, {len(self.names)} ancestor names "
              f"({time.time() - t0:.1f} s)", flush=True)

    def lineage(self, ncbi_id):
        """Lower-cased lineage tuple of a taxon ID, or None if the pinned taxonomy has no such taxon."""
        text = str(ncbi_id).strip()
        try:
            taxid = int(text)
        except ValueError:
            taxid = None
        target = self.taxonomy.resolve(taxid) if taxid is not None else None
        if target is None:
            if text not in self.reported:
                self.reported.add(text)
                print(f"NCBI taxon {text!r} is neither in nodes.dmp nor in merged.dmp of the pinned taxonomy "
                      "(deleted or invalid ID): no lineage", flush=True)
            return None
        cached = self.lineages.get(target)
        if cached is not None:
            return cached
        ancestors = []
        parent = self.taxonomy.parent
        node = parent[target]
        for _ in range(1000):
            if target == self.root or node == self.root:
                break
            ancestors.append(self.names[node])
            node = parent[node]
        else:
            raise ValueError(f"nodes.dmp: no path from taxon {target} to the root")
        # the same text efetch's <Lineage> holds, split the same way ('' for a child of the root)
        lineage = tuple(s.lower() for s in '; '.join(reversed(ancestors)).split('; '))
        self.lineages[target] = lineage
        return lineage


def use_pinned_ncbi_taxonomy(nodes_dmp, names_dmp, merged_dmp):
    """Take NCBI lineages from the pinned taxdump from now on (loaded here, so a missing or
    broken file stops the run before any study is converted)."""
    global _PINNED_LINEAGE
    _PINNED_LINEAGE = PinnedNcbiLineage(nodes_dmp, names_dmp, merged_dmp)
    return _PINNED_LINEAGE


# ---- NCBI synonyms and merged IDs (pinned taxonomy only) ------------------------------
# Opt-in helpers for converters that have called use_pinned_ncbi_taxonomy(): a species name
# that is no allowed scientific name may be a synonym in names.dmp (renamed taxa:
# 'Blautia coccoides' -> 33035 Blautia producta), and a taxon ID may have been merged into
# another (merged.dmp). get_taxonomy_id_from_name__allowedTerms is unchanged.
_NCBI_SYNONYMS = None          # lower-cased synonym -> {taxon IDs}
_NCBI_SYNONYM_STEMS = None     # (genus, epithet stem) of binomial synonyms -> {taxon IDs}
_NCBI_SYNONYM_CLASSES = ('synonym', 'equivalent name')


def _epithet_stem(epithet):
    # Latin gender endings change when a species moves to another genus
    # ('Lingulodinium polyedrum' / 'polyedra')
    return re.sub(r'(?:us|um|a|is|e)$', '', epithet)


def _binomial_key(name):
    m = re.fullmatch(r'([a-z]+) ([a-z][a-z-]+)', name)
    return (m.group(1), _epithet_stem(m.group(2))) if m else None


def _load_ncbi_synonyms():
    global _NCBI_SYNONYMS, _NCBI_SYNONYM_STEMS
    if _NCBI_SYNONYMS is None:
        synonyms, stems = {}, {}
        with open(_PINNED_LINEAGE.taxonomy.names_dmp, encoding='utf-8') as handle:
            for line in handle:
                if 'synonym' not in line and 'equivalent name' not in line:
                    continue
                fields = line.split('\t|\t')
                if len(fields) < 4 or fields[3].rstrip('|\t\r\n').strip() not in _NCBI_SYNONYM_CLASSES:
                    continue
                name = ' '.join(fields[1].lower().split())
                taxid = int(fields[0])
                synonyms.setdefault(name, set()).add(taxid)
                key = _binomial_key(name)
                if key is not None:
                    stems.setdefault(key, set()).add(taxid)
        _NCBI_SYNONYMS, _NCBI_SYNONYM_STEMS = synonyms, stems
    return _NCBI_SYNONYMS, _NCBI_SYNONYM_STEMS


# placeholder taxa whose synonyms are words such as 'none', 'other', 'not specified'
# (32644 unidentified, 12908 unclassified sequences, 28384 other sequences, ...)
_NCBI_PLACEHOLDER_TAXA = {1, 131567, 32644, 12908, 28384, 2787823, 2787854, 48479}


def _allowed_term_of_ids(ids, allowedTerm_dict):
    resolved = {_PINNED_LINEAGE.taxonomy.resolve(i) for i in ids} - {None}
    if len(resolved) != 1 or resolved & _NCBI_PLACEHOLDER_TAXA:
        return None
    by_id, _ = _taxonomy_term_index(allowedTerm_dict['NCBITaxonomy']['allowed_values'])
    return by_id.get(str(resolved.pop()))


def pinned_ncbi_synonym_term(organism_name, allowedTerm_dict):
    """Allowed 'ID|Name' term of the one taxon a names.dmp synonym / equivalent name
    `organism_name` belongs to (also with another Latin gender ending of the epithet,
    'Lingulodinium polyedra' for the synonym 'Lingulodinium polyedrum'), merged IDs
    followed; None without the pinned taxonomy or when the name is ambiguous."""
    if _PINNED_LINEAGE is None or not isinstance(organism_name, str):
        return None
    name = ' '.join(organism_name.lower().split())
    if not name:
        return None
    synonyms, stems = _load_ncbi_synonyms()
    ids = synonyms.get(name)
    if ids:
        return _allowed_term_of_ids(ids, allowedTerm_dict)
    key = _binomial_key(name)
    if key is not None and key in stems:
        return _allowed_term_of_ids(stems[key], allowedTerm_dict)
    return None


def pinned_ncbi_current_term(ncbi_id, allowedTerm_dict):
    """Allowed 'ID|Name' term of a taxon ID merged into another one (merged.dmp), else None."""
    if _PINNED_LINEAGE is None:
        return None
    try:
        taxid = int(str(ncbi_id).strip())
    except ValueError:
        return None
    if taxid not in _PINNED_LINEAGE.taxonomy.merged:
        return None
    return _allowed_term_of_ids({taxid}, allowedTerm_dict)


# ---- Organism texts that no exact lookup resolves (pinned taxonomy only) ---------------
# ncbi_name_fallback() is the last try of a converter for an organism text its exact lookups
# (allowed scientific names, alias sheets, synonyms) do not resolve, used where the converter
# would otherwise give 'missing value' or keep an organism ID it does not trust. Curators misspell
# names and drag strain numbers down with the accessions: 'Arbidopsis thaliana' with NCBITaxon
# 3702, 3703, 4056 ... (MTBLS2522), 'Zymomonas mobilis subsp. mobilis ATCC 31916', 'ATCC 31917' ...
# with 264297, 264298 ... (MTBLS14081). The converter then decides as if the text had been written
# as the name found (Metabolights2REDU.organism_accession_terms: the text wins over a drag-filled or
# contradicting accession, microbes keep their exception).
#   0. Other white space only (a no-break space, blanks at the end, two blanks): the converter's exact
#      lookup of the text with its white space normalised.
#   a. A strain / culture-collection designation at the end ('... ATCC 31916', '... DSM 20016',
#      '... NRRL B-14023', '... (ATCC 25922)', '... DSM 30083 = ATCC 11775', '... strain K-12',
#      '... str. MG1655'; collections: _STRAIN_COLLECTIONS; not a 'strain ...' part that names another
#      organism after it: 'A strain X; B', 'A str. X and B'): the text itself when names.dmp lists it
#      (scientific name, synonym or equivalent name) for exactly one taxon, that strain; else the
#      text without the designation(s), when that has two words at least, by the converter's exact
#      lookup, as a names.dmp scientific name / synonym / equivalent name of exactly one taxon, or
#      by b.
#   b. One wrong character: the text (lower-cased, whitespace-normalised) has two words at least,
#      the first two of letters only (a hyphen inside the epithet), the others of letters or an
#      infraspecific rank ('subsp.', 'var.' ...), no abbreviation ('E. coli'), no 'sp.', 'spp',
#      'cf.', hybrid 'x', 'and', no digits, 10 characters at least; it is no name in names.dmp
#      itself, and the names (any class, 9+ characters) at Damerau-Levenshtein distance 1 (one
#      character inserted, deleted or replaced, or two neighbours swapped) all belong to one taxon,
#      at species rank or below, and include its scientific name, or (an old name misspelt:
#      'Clostridium clostridioforme' for the synonym 'Clostridium clostridiiforme' of Enterocloster
#      clostridioformis) include a synonym or equivalent name of it while no name of any other taxon
#      lies within distance 2 of the text (optimal string alignment over all names of names.dmp: one
#      wrong character next to an old name is no evidence where the name of another taxon is near too,
#      'Clostridium hystolyticum' between the synonym 'histolyticum' and 'Clostridium hydrolyticum').
#      Names of other classes (common names, 'includes', 'in-part') only count toward the one taxon.
#      Nothing when names of two taxa are that close, for one-word texts and for placeholder taxa.
# The names are indexed once per process, when a text first needs it: every name of names.dmp
# (lower-cased, whitespace-normalised, 8+ characters) under two keys, its length with its first and
# with its last k = (length - 2) // 2 characters. A name within distance 1 of a text of length n has
# length n - 1, n or n + 1 and shares one of its two keys with the text (an edit leaves either the
# first or the last k characters as they are), so the texts' keys find all of them. Python's hash of
# the keys only selects candidates (each is compared with the text), so the results do not depend on
# PYTHONHASHSEED. The names within distance 2 (only for the texts b. asks for them) are found by a
# banded dynamic program over all names of length n - 2 to n + 2 at once (numpy), per length a matrix of
# the names' bytes built when first needed; names with other than ASCII characters (the texts are
# ASCII) are compared one by one.
_STRAIN_COLLECTIONS = ('ATCC', 'DSMZ', 'DSM', 'NRRL', 'NCTC', 'JCM', 'CBS', 'CECT', 'LMG', 'NBRC', 'IFO', 'KCTC',
                       'CCUG', 'CIP', 'NCIMB',
                       # algae, cyanobacteria and other collections: Pasteur Culture Collection of Cyanobacteria,
                       # NCIB (the former name of NCIMB), UTEX, CCAP, SAG (Goettingen), NIES, NCMA and its former
                       # name CCMP, Roscoff (RCC), NCCB, BEI Resources, Nestle Culture Collection (NCC)
                       'PCC', 'NCIB', 'UTEX', 'CCAP', 'SAG', 'NIES', 'NCMA', 'CCMP', 'RCC', 'NCCB', 'BEI', 'NCC')
_STRAIN_DESIGNATION = re.compile(
    r'(?:^|[\s,;=/(]+)(?P<designation>(?:strain|str\.)\s+\S.*'
    r'|(?:' + '|'.join(_STRAIN_COLLECTIONS) + r')(?:\s*[®:#-])?\s*(?:[A-Z]{1,4}[\s-]?)?\d[\w./-]*™?(?:\s*\(T\)|T)?'
    r')[\s)]*$', re.I)
# a 'strain ...' part that goes on with another organism ('A strain X; B', 'A str. X and B') is no designation
_ORGANISM_LIST_PART = re.compile(r'[;,+&]|\s(?:and|or|with|x)\s', re.I)
_FALLBACK_WORD = re.compile(r'[a-z]+(?:-[a-z]+)*')
_INFRASPECIFIC_RANKS = frozenset({'subsp.', 'ssp.', 'var.', 'f.', 'forma', 'subvar.', 'pv.', 'bv.', 'sv.', 'cv.'})
# words that make a text no single organism name ('Streptomyces spp', 'A x B', 'A and B', 'leaves of A')
_NO_NAME_WORDS = frozenset({'sp', 'spp', 'ssp', 'subsp', 'var', 'cf', 'aff', 'nov', 'gen', 'sect', 'x',
                            'and', 'or', 'et', 'with', 'of', 'in', 'from', 'on', 'the'})
_NAME_CLASS_CODES = {'scientific name': 0, 'synonym': 1, 'equivalent name': 2}   # every other class: 3
_FALLBACK_MIN_LENGTH = 10
_SPECIES_COLUMN = 7   # ncbi_database.DESIRED_RANKS.index('species')
_NAME_INDEX = None


def _within_one_edit(a, b):
    """0 when a == b, 1 when their Damerau-Levenshtein (optimal string alignment) distance is 1, else None."""
    if a == b:
        return 0
    la, lb = len(a), len(b)
    if la == lb:
        i = 0
        while a[i] == b[i]:
            i += 1
        if a[i + 1:] == b[i + 1:]:
            return 1
        if i + 1 < la and a[i] == b[i + 1] and a[i + 1] == b[i] and a[i + 2:] == b[i + 2:]:
            return 1
        return None
    if la > lb:
        a, b, la, lb = b, a, lb, la
    if lb != la + 1:
        return None
    i = 0
    while i < la and a[i] == b[i]:
        i += 1
    return 1 if a[i:] == b[i + 1:] else None


class _NcbiNameIndex:
    """The names of names.dmp (lower-cased, whitespace-normalised, at least _FALLBACK_MIN_LENGTH - 2
    characters) with their taxon and name class, for names within one edit of a text, and within two (see
    above)."""

    def __init__(self, names_dmp):
        import numpy as np
        from array import array
        t0 = time.time()
        blob = bytearray()
        ends, taxa, classes = array('q'), array('i'), array('b')
        first, last = array('q'), array('q')
        # (8 characters: the names within two edits of a text of 10, the shortest text b. reads; near() compares
        # names of 9+ only)
        shortest = _FALLBACK_MIN_LENGTH - 2
        with open(names_dmp, encoding='utf-8', errors='replace') as handle:
            for line in handle:
                fields = line.split('\t|\t')
                if len(fields) < 4:
                    continue
                name = ' '.join(fields[1].lower().split())
                n = len(name)
                if n < shortest:
                    continue
                k = (n - 2) // 2
                first.append(hash((n, name[:k])))
                last.append(hash((n, name[-k:])))
                blob += name.encode('utf-8')
                ends.append(len(blob))
                taxa.append(int(fields[0]))
                classes.append(_NAME_CLASS_CODES.get(fields[3].rstrip('|\t\r\n').strip(), 3))
        self.blob = bytes(blob)
        del blob
        self.ends = np.frombuffer(ends, dtype=np.int64)
        self.taxa = np.frombuffer(taxa, dtype=np.int32)
        self.classes = np.frombuffer(classes, dtype=np.int8)
        self.keys = []
        for column in (first, last):
            keys = np.frombuffer(column, dtype=np.int64)
            order = np.argsort(keys, kind='stable').astype(np.int32)
            self.keys.append((keys[order], order))
        del first, last
        self.memo = {}
        self.by_length = {}      # byte length -> (rows, matrix of the ASCII names' bytes), built when first needed
        self.non_ascii = None    # (rows, names) of the names with other than ASCII characters
        print(f"NCBI name index for organism texts: {len(self.ends)} names of {names_dmp} "
              f"({time.time() - t0:.1f} s)", flush=True)

    def name(self, i):
        start = int(self.ends[i - 1]) if i else 0
        return self.blob[start:int(self.ends[i])].decode('utf-8')

    def near(self, text):
        """[(taxon ID, name class code, distance)] of the names within one edit of text (normalised)."""
        hit = self.memo.get(text)
        if hit is not None:
            return hit
        import numpy as np
        out, seen = [], set()
        for n in (len(text) - 1, len(text), len(text) + 1):
            if n < _FALLBACK_MIN_LENGTH - 1:
                continue
            k = (n - 2) // 2
            for (keys, rows), part in zip(self.keys, (text[:k], text[-k:])):
                h = hash((n, part))
                lo, hi = np.searchsorted(keys, h, 'left'), np.searchsorted(keys, h, 'right')
                for i in rows[lo:hi].tolist():
                    if i in seen:
                        continue
                    seen.add(i)
                    d = _within_one_edit(text, self.name(i))
                    if d is not None:
                        out.append((int(self.taxa[i]), int(self.classes[i]), d))
        out.sort()
        self.memo[text] = out
        return out

    def _ascii_names_of_length(self, m):
        import numpy as np
        hit = self.by_length.get(m)
        if hit is None:
            starts = np.concatenate(([0], self.ends[:-1]))
            lengths = self.ends - starts
            blob = np.frombuffer(self.blob, dtype=np.uint8)
            rows = np.nonzero(lengths == m)[0]
            matrix = blob[starts[rows][:, None] + np.arange(m)[None, :]] if len(rows) else np.zeros((0, m), np.uint8)
            ascii_rows = (matrix < 128).all(axis=1)
            hit = self.by_length[m] = (rows[ascii_rows], np.ascontiguousarray(matrix[ascii_rows]))
        return hit

    def _non_ascii_names(self):
        import numpy as np
        if self.non_ascii is None:
            blob = np.frombuffer(self.blob, dtype=np.uint8)
            starts = np.concatenate(([0], self.ends[:-1])).astype(np.int64)
            # (every name has 8+ bytes: no empty segment)
            rows = np.nonzero(np.maximum.reduceat(blob, starts) >= 128)[0] if len(starts) else np.zeros(0, np.int64)
            self.non_ascii = (rows, [self.name(int(i)) for i in rows])
        return self.non_ascii

    def within_two(self, text):
        """[taxon ID] of the names within optimal string alignment distance 2 of an ASCII text (normalised)."""
        import numpy as np
        q = np.frombuffer(text.encode('ascii'), dtype=np.uint8)
        n, cap = len(q), 3
        found = []
        for m in range(max(1, n - 2), n + 3):
            rows, a = self._ascii_names_of_length(m)
            if not len(rows):
                continue
            count = len(rows)
            prev2 = None
            prev = np.broadcast_to(np.minimum(np.arange(m + 1), cap).astype(np.int8), (count, m + 1)).copy()
            for i in range(1, n + 1):
                cur = np.full((count, m + 1), cap, dtype=np.int8)
                cur[:, 0] = min(i, cap)
                for j in range(max(1, i - 2), min(m, i + 2) + 1):
                    cost = (a[:, j - 1] != q[i - 1]).astype(np.int8)
                    v = np.minimum(np.minimum(prev[:, j] + 1, cur[:, j - 1] + 1), prev[:, j - 1] + cost)
                    if i > 1 and j > 1:
                        swap = (a[:, j - 1] == q[i - 2]) & (a[:, j - 2] == q[i - 1])
                        v = np.where(swap, np.minimum(v, prev2[:, j - 2] + 1), v)
                    cur[:, j] = np.minimum(v, cap)
                prev2, prev = prev, cur
            found.extend(int(t) for t in self.taxa[rows[prev[:, m] <= 2]])
        rows, names = self._non_ascii_names()
        for i, name in zip(rows.tolist(), names):
            if abs(len(name) - n) <= 2 and _osa_distance(text, name, 2) <= 2:
                found.append(int(self.taxa[i]))
        return found


def _osa_distance(a, b, limit):
    """Optimal string alignment distance of a and b, or limit + 1 when it is larger than limit."""
    if abs(len(a) - len(b)) > limit:
        return limit + 1
    prev2, prev = None, list(range(len(b) + 1))
    for i in range(1, len(a) + 1):
        cur = [i] + [0] * len(b)
        for j in range(1, len(b) + 1):
            v = min(prev[j] + 1, cur[j - 1] + 1, prev[j - 1] + (a[i - 1] != b[j - 1]))
            if i > 1 and j > 1 and a[i - 1] == b[j - 2] and a[i - 2] == b[j - 1]:
                v = min(v, prev2[j - 2] + 1)
            cur[j] = v
        if min(cur) > limit:
            return limit + 1
        prev2, prev = prev, cur
    return min(prev[-1], limit + 1)


def _ncbi_name_index():
    global _NAME_INDEX
    if _NAME_INDEX is None:
        _NAME_INDEX = _NcbiNameIndex(_PINNED_LINEAGE.taxonomy.names_dmp)
    return _NAME_INDEX


def _at_or_below_species(taxid):
    """Whether a taxon (pinned taxonomy) has rank species or lies below a species."""
    tax = _PINNED_LINEAGE.taxonomy
    node = tax.resolve(taxid)
    for _ in range(1000):
        if node is None:
            return False
        column = tax.column[tax.rank[node]]
        if column == _SPECIES_COLUMN:
            return True
        if column is not None:
            return False   # genus or a rank above it
        up = tax.parent[node]
        if up == node or up < 0:
            return False
        node = up
    return False


def _one_taxon(taxids):
    resolved = {_PINNED_LINEAGE.taxonomy.resolve(t) for t in taxids} - {None}
    if len(resolved) != 1 or resolved & _NCBI_PLACEHOLDER_TAXA:
        return None
    return resolved.pop()


def _named_taxon(text):
    """The one taxon names.dmp lists text (normalised) for as scientific name, synonym or equivalent name, else None."""
    return _one_taxon({t for t, c, d in _ncbi_name_index().near(text) if d == 0 and c < 3})


def _typo_shape(text):
    words = text.split()
    if len(words) < 2 or len(text) < _FALLBACK_MIN_LENGTH:
        return False
    if any(w in _NO_NAME_WORDS for w in words):
        return False
    if not (_FALLBACK_WORD.fullmatch(words[0]) and _FALLBACK_WORD.fullmatch(words[1])):
        return False
    return all(_FALLBACK_WORD.fullmatch(w) or w in _INFRASPECIFIC_RANKS for w in words[2:])


def _typo_taxon(text):
    """The taxon of b. (block comment above) for a normalised text, else None."""
    if not _typo_shape(text):
        return None
    near = _ncbi_name_index().near(text)
    if not near or any(d == 0 for t, c, d in near):
        return None
    taxon = _one_taxon({t for t, c, d in near})
    if taxon is None or not _at_or_below_species(taxon):
        return None
    if any(c == 0 for t, c, d in near):
        return taxon
    # an old name misspelt: a synonym / equivalent name of the taxon one edit away, no other taxon within two
    if any(c in (1, 2) for t, c, d in near) and _only_taxon_within_two(text, taxon):
        return taxon
    return None


def _only_taxon_within_two(text, taxon):
    """Whether every name of names.dmp within distance 2 of text is a name of taxon (merged IDs followed)."""
    resolve = _PINNED_LINEAGE.taxonomy.resolve
    others = {resolve(t) for t in _ncbi_name_index().within_two(text)} - {taxon, None}
    return not others


def ncbi_name_fallback(text, term_of_id, exact=None):
    """(allowed 'ID|Name' term, how) for an organism text no exact lookup resolves, by the rules of the
    block comment above, else (None, None); how: 'white space' (the text with its white space normalised, by
    exact), 'strain name' (the text is an NCBI name of that strain), 'strain' (the text without its designation),
    'strain, typo', 'typo'. term_of_id(ID as str): the
    allowed term of a taxon ID or None; exact(name): the converter's own lookup for the text without its
    designation (allowed names, alias sheets) or None. Needs the pinned taxonomy (use_pinned_ncbi_taxonomy)."""
    if _PINNED_LINEAGE is None or not isinstance(text, str):
        return None, None
    spaced = ' '.join(text.split())
    if spaced != text and spaced and exact is not None:
        # only other white space (a no-break space, a blank at the end, two blanks)
        term = exact(spaced)
        if term:
            return term, 'white space'
    norm = spaced.lower()
    if len(norm) < _FALLBACK_MIN_LENGTH or not re.search(r'[a-z]{3}', norm):
        return None, None
    rest, stripped = spaced, False
    for _ in range(5):
        m = _STRAIN_DESIGNATION.search(rest)
        if m is None or m.start() == 0 or _ORGANISM_LIST_PART.search(m.group('designation')):
            break
        rest, stripped = rest[:m.start()].rstrip(' ,;=/('), True
    if stripped:
        taxon = _named_taxon(norm)
        if taxon is not None:
            term = term_of_id(str(taxon))
            if term:
                return term, 'strain name'
        if len(rest.split()) < 2:
            return None, None
        term = exact(rest) if exact is not None else None
        if term:
            return term, 'strain'
        rest_norm = ' '.join(rest.lower().split())
        taxon = _named_taxon(rest_norm)
        if taxon is not None:
            term = term_of_id(str(taxon))
            if term:
                return term, 'strain'
        taxon = _typo_taxon(rest_norm)
        term = term_of_id(str(taxon)) if taxon is not None else None
        return (term, 'strain, typo') if term else (None, None)
    taxon = _typo_taxon(norm)
    term = term_of_id(str(taxon)) if taxon is not None else None
    return (term, 'typo') if term else (None, None)


# ---- 'Genus sp. <designation>' ------------------------------------------------------------------------------------
# NCBI names many unnamed species by the designation of a strain, isolate or clone after 'sp.' ('Synechocystis sp.
# PCC 6803', 'Streptomyces sp. CNQ-509'): a text or 'ID|Name' naming such an allowed taxon names that organism
# (read_and_validate_redu_from_github.complete_and_fill_REDU_table, NCBITaxonomy). Two exceptions, applied by the
# converters to the terms their name and ID lookups give (sp_name_term) and by the validator to its own lookups:
#   * a designation that is a lab's numbering rather than one of its own (sp_generic_label: after an optional
#     'strain', 'isolate', 'no.' or '#', at most two letters or digits, one to three digits, or a Roman numeral up to
#     X: 'Bacillus sp. 1', 'Nostoc sp. A', 'Streptomyces sp. 12', 'Genus sp. III') names the genus only, unless the
#     given ID is that taxon: NCBI lists such names for other submitters' isolates ('707210|Bacillus sp. 1',
#     '492787|Streptomyces sp. 2'); 'Xylaria sp. P1' with the ID of that taxon stays the taxon;
#   * a text 'Genus sp. X' with an ID whose NCBI name is 'Genus sp. Y' of another designation (sp_other_strain) names
#     the genus, the ID being that of another strain ('585110|Pseudomonas sp. 2 1 26': 585110 is 'Pseudomonas sp.
#     2156'), unless NCBI lists the text for that taxon (another collection's designation: 'Synechocystis sp. ATCC
#     27184' for PCC 6803). Two designations are the same strain when they share a variant: the whole designation and
#     its parts at NCBI's separators (' str. ', ' strain ', ' substr. ', ' = ', parentheses), letters and digits
#     only, lower-cased ('oral taxon 186 str. F0373' gives oraltaxon186strf0373, oraltaxon186 and f0373; 'PCC 7120 =
#     FACHB-418' pcc7120fachb418, pcc7120 and fachb418). So 'F0373', 'PCC 7120', 'CNQ525' (NCBI: 'CNQ-525'), '1 1 30'
#     (NCBI: '1_1_30') name NCBI's strain, and '1' / '12', '2' / '2156', 'PCC 680' / 'PCC 6803' do not; a designation
#     inside another is no evidence.
# The genus: the NCBI genus of the term's taxon in the pinned taxonomy when its name is the text's first word
# (pinned_genus_term), else (no pinned taxonomy, a taxon NCBI has put in another genus) the converter's term for
# 'Genus sp.', which the validator reads as the genus; else none.
_SP_WORDS = frozenset({'sp', 'spp'})
_SP_DESIGNATION_PARTS = re.compile(r'(?i)\s+(?:str\.|strain|substr\.|=)\s+|^(?:str\.|strain)\s+|\s*[()]\s*')
_SP_LABEL_PREFIX = re.compile(r'(?i)^(?:(?:strain|isolate)\s+|(?:str|no|nr)\.\s*|#\s*)')
_SP_GENERIC_LABEL = re.compile(r'[a-z0-9]{1,2}|\d{1,3}|iii|vii|viii')
_GENUS_COLUMN = 6   # ncbi_database.DESIRED_RANKS.index('genus')


def sp_designation(name):
    """The designation of a name 'Genus sp. <designation>' / 'Genus spp. <designation>', white space normalised
    ('' for other names and for 'Genus sp.' alone)."""
    words = str(name).split()
    if len(words) < 3 or words[1].rstrip('.').lower() not in _SP_WORDS:
        return ''
    return ' '.join(words[2:])


def sp_designation_variants(name):
    """The variants of the designation of a name 'Genus sp. <designation>' (block comment above): a frozenset, empty
    for other names."""
    designation = sp_designation(name)
    if not designation:
        return frozenset()
    parts = [designation] + _SP_DESIGNATION_PARTS.split(designation)
    return frozenset(v for v in (re.sub(r'[^a-z0-9]', '', p.lower()) for p in parts if p) if v)


def sp_other_strain(text, name):
    """Whether a text 'Genus sp. X' and a name 'Genus sp. Y' name two strains: both have a designation and the two
    share no variant (block comment above)."""
    mine, theirs = sp_designation_variants(text), sp_designation_variants(name)
    return bool(mine and theirs and not mine & theirs)


def sp_generic_label(name):
    """Whether the designation of a name 'Genus sp. <designation>' is a lab's numbering (block comment above)."""
    label = re.sub(r'[^a-z0-9]', '', _SP_LABEL_PREFIX.sub('', sp_designation(name)).lower())
    return bool(label) and _SP_GENERIC_LABEL.fullmatch(label) is not None


def pinned_genus_term(taxid, genus, term_of_id):
    """term_of_id(ID) of the NCBI genus of taxon `taxid` (the taxon itself or its ancestor at rank genus, merged IDs
    followed) in the pinned taxonomy when the genus is named `genus` (case-insensitive), else None (also without the
    pinned taxonomy)."""
    if _PINNED_LINEAGE is None or not str(genus).strip():
        return None
    tax = _PINNED_LINEAGE.taxonomy
    try:
        node = tax.resolve(int(str(taxid).strip()))
    except (TypeError, ValueError):
        return None
    if node is None:
        return None
    found = tax.lineage_ids(node)[_GENUS_COLUMN]
    if found is None or _PINNED_LINEAGE.names.get(found, '').lower() != str(genus).strip().lower():
        return None
    return term_of_id(str(found))


def sp_name_term(text, term, given_id, term_of_id, term_of_name=None, text_term=None):
    """The NCBITaxonomy term for an organism text and the allowed term a converter's lookups found for it (`term`, by
    the given ID `given_id` or by name; None if none) by the two exceptions of the block comment above: for a text
    'Genus sp. <designation>', its genus where the designation is a lab's numbering and the term is not the given ID's
    taxon, or where the term is the given ID's strain of another designation and not `text_term` (the term the text
    names by the converter's name lookups, None if none); else `term` (also a term of several organisms). The genus:
    pinned_genus_term of the term's taxon, else term_of_name('Genus sp.') (the converter's exact name lookup), else
    None. term_of_id(ID): the allowed term of a taxon ID or None."""
    if not term or is_multi(term) or not sp_designation(text):
        return term
    taxid, _, name = str(term).partition('|')
    taxid = taxid.strip()
    given = taxid == str(given_id or '').strip()
    if not (sp_generic_label(text) and not given) and not (given and term != text_term and sp_other_strain(text, name)):
        return term
    genus = str(text).split()[0]
    found = pinned_genus_term(taxid, genus, term_of_id)
    if found is None and term_of_name is not None:
        found = term_of_name(f'{genus} sp.')
    return found or None


# ---- a bare 'Genus sp.' with the ID of a species or strain of that genus -------------------------------------------
# An organism text at genus level (genus_level_sp: 'Genus sp.', 'Genus spp.', the period optional, also NCBI's own
# genus-level names 'Bacillus sp. (in: firmicutes)', 'Saccharimonas sp. (SeqCode)') that comes with the taxon ID of a
# species or strain of that genus names that organism: the ID is more specific than the text and agrees with it (user
# decision 2026-10-18, "if there is a more specific id take the id"; in every source: the validator's 'ID|Name' of GNPS
# and MASST, MetaboLights' organism text and its accession, Workbench's SUBJECT_SPECIES and TAXONOMY_ID). In the pinned
# NCBI taxonomy (bare_sp_id_term), the ID's taxon (merged IDs followed)
#   * has rank species or lies below a species (a strain, a subspecies), and is none of NCBI's placeholders at genus
#     level ('Streptomyces sp.', 'Bacillus sp. (in: firmicutes)', 'uncultured Streptomyces sp.', 'unidentified
#     Pseudomonas': NCBI_COLLECTING_TAXON);
#   * lies in the genus the text names: the ID's NCBI genus has the text's genus name (case-insensitive; the ID tells
#     which genus of a name that several genera share: 'Mallotus sp.' with the capelin's ID is the capelin, not a plant;
#     'Mus sp.' with the house mouse's ID the house mouse, 'Mus' also naming a subgenus); or, where the text's genus
#     name is the scientific name of no taxon of the pinned taxonomy (names.dmp, any rank: a word that names a genus
#     there names that genus, also a genus without species, and one of another rank names that taxon: 'Bacillus sp.'
#     with the ID of Priestia megaterium, a species NCBI moved out of Bacillus, is an ID of another genus; 'Micromonas'
#     is a green alga, not the bacterial synonym 'Micromonas <firmicutes>' of Parvimonas; 'Bacteria sp.', 'Fungi sp.'
#     are no genus texts) nor an allowed name of a taxon of another rank, names.dmp lists that name as a synonym /
#     equivalent name of the ID's NCBI genus (a renamed genus: 'Verrucosispora' of Micromonospora, 'Hypocrea' of
#     Trichoderma) or the binomial of the ID's species under the text's genus name (a species of a genus name NCBI no
#     longer has; the first two words of the names.dmp name, also with another Latin gender ending of the epithet);
#   * for NCBI's qualified name 'Genus sp. (in: <group>)', lies in the NCBI genus of the taxon of that name (the note
#     tells which genus: 'Bacillus sp. (in: walking sticks)' with a bacterium's ID is no such pair).
# Otherwise None, and each source decides as before: an ID of the genus itself, of a taxon between the genus and its
# species (a subgenus, 'unclassified Streptomyces') or of NCBI's 'Genus sp.' gives the genus; an ID of another genus or
# above the genus follows the source's own rules. A drag-filled MetaboLights accession is no ID of the text
# (Metabolights2REDU.organism_accession_terms). 'Genus sp. <strain, isolate, clone>' is no text at genus level
# (sp_name_term above).
_GENUS_LEVEL_SP = re.compile(r'(?i)\S+ spp?\.?(?: \((?:in: [^()]*|SeqCode)\))?')
_QUALIFIED_SP = re.compile(r'(?i)\S+ spp?\.? \((?:in: [^()]*|SeqCode)\)')
# names.dmp synonyms / equivalent names (lower-cased) for bare_sp_id_term, read once where needed: one-word names -> taxon
# IDs, (genus, epithet stem) of the first two words of a name ('Lactobacillus plantarum (Orla-Jensen 1919) Bergey et
# al. 1923 ...') -> taxon IDs; and the one-word scientific names (lower-cased) of all taxa
_GENUS_SYNONYMS = None
# NCBI's taxa that collect the unnamed members of a taxon ('unclassified Sulfitobacter', 'Pseudomonas environmental
# samples') and its placeholders for organisms known by a sequence only ('uncultured Streptomyces sp.', 'unidentified
# Pseudomonas'): no one organism either
NCBI_COLLECTING_TAXON = re.compile(r'(?i)(?:unclassified|uncultured|unidentified) \S.*|(?:\S.* )?environmental samples')


def genus_level_sp(name):
    """Whether an organism name is an answer at genus level: 'Genus sp.' / 'Genus spp.' (the period optional) with nothing
    after it but NCBI's note on the genus name ('Bacillus sp. (in: firmicutes)', 'Saccharimonas sp. (SeqCode)'); white
    space normalised."""
    return isinstance(name, str) and _GENUS_LEVEL_SP.fullmatch(' '.join(name.split())) is not None


def _species_of(tax, node):
    """The taxon itself or its ancestor of rank species where taxon `node` has rank species or lies below a species,
    else None (a genus, a rank above it, a taxon between a genus and its species)."""
    for _ in range(1000):
        column = tax.column[tax.rank[node]]
        if column == _SPECIES_COLUMN:
            return node
        if column is not None:
            return None
        up = tax.parent[node]
        if up == node or up < 0:
            return None
        node = up
    return None


def _genus_synonyms(names_dmp):
    global _GENUS_SYNONYMS
    if _GENUS_SYNONYMS is None:
        t0 = time.time()
        words, binomials, scientific = {}, {}, set()
        with open(names_dmp, encoding='utf-8') as handle:
            for line in handle:
                if 'scientific name' not in line and 'synonym' not in line and 'equivalent name' not in line:
                    continue
                fields = line.split('\t|\t')
                if len(fields) < 4:
                    continue
                name_class = fields[3].rstrip('|\t\r\n').strip()
                if name_class == 'scientific name':
                    name = fields[1].split()
                    if len(name) == 1:
                        scientific.add(name[0].lower())
                    continue
                if name_class not in _NCBI_SYNONYM_CLASSES:
                    continue
                name = fields[1].lower().split()
                if len(name) == 1:
                    words.setdefault(name[0], set()).add(int(fields[0]))
                elif len(name) > 1:
                    key = _binomial_key(' '.join(name[:2]))
                    if key is not None:
                        binomials.setdefault(key, set()).add(int(fields[0]))
        _GENUS_SYNONYMS = (words, binomials, frozenset(scientific))
        print(f"names.dmp names for genus-level organism texts: {len(words)} one-word synonyms, {len(binomials)} "
              f"binomials, {len(scientific)} one-word scientific names ({time.time() - t0:.1f} s)", flush=True)
    return _GENUS_SYNONYMS


def _genus_of_term(tax, term):
    """The NCBI genus (pinned taxonomy, merged IDs followed) of the taxon of an allowed 'ID|Name' term, else None."""
    try:
        node = tax.resolve(int(str(term).split('|', 1)[0].strip()))
    except (TypeError, ValueError):
        return None
    return tax.lineage_ids(node)[_GENUS_COLUMN] if node is not None else None


def bare_sp_id_term(text, given_id, term_of_id, term_of_name, pinned=None):
    """The allowed term of the taxon of `given_id` where `text` is a bare 'Genus sp.' and the ID that of a species or
    strain of that genus (block comment above), else None. term_of_id(ID): the allowed term of a taxon ID or None;
    term_of_name(name): the source's exact name lookup (allowed name, case-insensitive) or None; pinned: the pinned
    taxonomy (PinnedNcbiLineage; default the one use_pinned_ncbi_taxonomy loaded; without one, None)."""
    pinned = pinned if pinned is not None else _PINNED_LINEAGE
    if pinned is None or not genus_level_sp(text):
        return None
    tax = pinned.taxonomy
    try:
        node = tax.resolve(int(str(given_id).strip()))
    except (TypeError, ValueError):
        return None
    species = _species_of(tax, node) if node is not None else None
    if species is None:
        return None
    term = term_of_id(str(node))
    if not term or is_multi(term):
        return None
    name = term.split('|', 1)[-1].strip()
    if genus_level_sp(name) or NCBI_COLLECTING_TAXON.fullmatch(name):
        return None
    genus = tax.lineage_ids(node)[_GENUS_COLUMN]
    if genus is None:
        return None
    text = ' '.join(text.split())
    if _QUALIFIED_SP.fullmatch(text):
        named = term_of_name(text)
        return term if named and not is_multi(named) and _genus_of_term(tax, named) == genus else None
    word = text.split()[0].lower()
    if pinned.names.get(genus, '').lower() == word:
        return term
    # the names.dmp readings: a genus name that is the scientific name of no taxon (that taxon: another genus, a taxon of
    # another rank) and no allowed name of a taxon of another rank than genus
    words, binomials, scientific = _genus_synonyms(tax.names_dmp)
    if word in scientific:
        return None
    named = term_of_name(word)
    if named:
        try:
            named_node = tax.resolve(int(named.split('|', 1)[0].strip()))
        except (TypeError, ValueError):
            named_node = None
        if named_node is None or tax.column[tax.rank[named_node]] != _GENUS_COLUMN:
            return None
    if genus in words.get(word, ()):
        return term
    species_name = (name if species == node else pinned.names.get(species, '')).lower().split()
    key = _binomial_key(f'{word} {species_name[1]}') if len(species_name) == 2 else None
    if key is not None and species in {tax.resolve(i) for i in binomials.get(key, ())}:
        return term
    return None


# NCBI allows 3 requests/s without an API key; going faster returns HTTP 429
_NCBI_MIN_INTERVAL = 0.4
_ncbi_last_request = [0.0]
NCBI_TIMEOUT = (10, 60)
_LIVE_LINEAGE_WARNED = [False]


def _ncbi_throttle():
    wait = _NCBI_MIN_INTERVAL - (time.time() - _ncbi_last_request[0])
    if wait > 0:
        time.sleep(wait)
    _ncbi_last_request[0] = time.time()


def _ncbi_get(url, retries=4):
    """GET an NCBI URL, throttled, retrying on throttling, server and connection errors.

    Returns the response if it has status 200, raises otherwise.
    """
    for attempt in range(1, retries + 1):
        _ncbi_throttle()
        try:
            response = requests.get(url, timeout=NCBI_TIMEOUT)
        except (requests.ConnectionError, requests.Timeout) as e:
            if attempt == retries:
                raise
            print(f"NCBI request failed (attempt {attempt}/{retries}): {e}")
        else:
            if response.status_code == 200:
                return response
            if attempt == retries or (response.status_code != 429 and response.status_code < 500):
                response.raise_for_status()
            print(f"NCBI responded with status {response.status_code} (attempt {attempt}/{retries})")
        time.sleep(min(2 ** attempt, 30))


class _NoLiveLineage(Exception):
    """E-utilities answered without a lineage (not cached by _live_ncbi_lineage)."""


@functools.lru_cache(maxsize=None)
def _live_ncbi_lineage(ncbi_id, retries=3):
    for attempt in range(1, retries + 1):
        response = _ncbi_get(f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=taxonomy&id={ncbi_id}")
        soup = BeautifulSoup(response.text, "xml")
        taxon = soup.find("Taxon")
        lineage = taxon.find("Lineage") if taxon is not None else None
        if lineage is not None:
            return tuple(s.lower() for s in lineage.text.split("; "))
        # NCBI occasionally answers 200 with an error document; retry before giving up
        print(f"No lineage in NCBI response for {ncbi_id} (attempt {attempt}/{retries})")
        time.sleep(2)
    raise _NoLiveLineage(ncbi_id)


def _get_ncbi_lineage(ncbi_id, retries=3):
    """Returns the lowercased lineage of an NCBI taxon id, or None if there is none.

    From the pinned taxonomy after use_pinned_ncbi_taxonomy(); otherwise from E-utilities,
    where found lineages are cached (the same ids recur across studies) but a missing one is
    not, and network errors propagate.
    """
    if _PINNED_LINEAGE is not None:
        return _PINNED_LINEAGE.lineage(ncbi_id)
    if not _LIVE_LINEAGE_WARNED[0]:
        _LIVE_LINEAGE_WARNED[0] = True
        print("WARNING: NCBI lineages come from the live E-utilities, not from the pinned taxonomy "
              "(use_pinned_ncbi_taxonomy was not called)", flush=True)
    try:
        return _live_ncbi_lineage(ncbi_id, retries)
    except _NoLiveLineage:
        return None


def _sample_type_from_lineage(classification, cell_culture_key1='', cell_culture_key2=''):
    cell_culture_key_words = ["cell", "media", "culture"]
    SampleType = None
    SampleTypeSub1 = None

    if "viridiplantae" in classification:
        SampleType = "plant"
        SampleTypeSub1 = "plant_NOS"
        if "Algae" in classification or 'rhodophyta' in classification or "phaeophyceae" in classification:
            SampleType = "algae"
        if "chlorophyta" in classification or "microalgae" in classification or "microalga" in classification:
            SampleType = "microalgae"
        if "streptophyta" in classification:
            SampleTypeSub1 = "plant_angiospermae"
        if "cyanobacteria" in classification:
            SampleTypeSub1 = "marine_cyanobacteria_insitu"
        if "bacillariophyta" in classification:
            SampleTypeSub1 = "marine_diatom"
    elif "metazoa" in classification:
        SampleType = "animal"
        if "mammalia" in classification and any(
                word in cell_culture_key1.lower() for word in cell_culture_key_words) or any(
                word in cell_culture_key2.lower() for word in cell_culture_key_words):
            SampleType = "culture_mammalian"
            SampleTypeSub1 = "culture_mammalian"
        if "amphibia" in classification:
            if "Caudata" in classification or "urodela" in classification or "echinodermata" in classification:
                SampleTypeSub1 = "salamander"
            else:
                SampleTypeSub1 = "frog"
        if "insecta" in classification:
            SampleTypeSub1 = "insect"
        if "porifera" in classification or "mollusca" in classification:
            SampleTypeSub1 = "marine_invertebrates"
        if "cnidaria" in classification:
            SampleTypeSub1 = "marine_coral"
    elif "fungi" in classification:
        SampleType = "culture_fungal"
        SampleTypeSub1 = "culture_fungal"
    elif "bacteria" in classification:
        SampleType = "culture_bacterial"
        SampleTypeSub1 = "culture_bacterial"
    return [SampleType, SampleTypeSub1]


def get_taxonomy_info(ncbi_id, cell_culture_key1 = '', cell_culture_key2 = ''):
    if ncbi_id is not None and ncbi_id != "NA":
        if _PINNED_LINEAGE is not None:
            # errors of the pinned taxonomy itself propagate; an ID without a lineage was
            # reported by the lookup
            classification = _PINNED_LINEAGE.lineage(ncbi_id)
            if classification is None:
                return [None, None]
            try:
                return _sample_type_from_lineage(classification, cell_culture_key1, cell_culture_key2)
            except Exception as e:
                print(f"Could not get taxonomy info for {ncbi_id}: {type(e).__name__}: {e}")
                return [None, None]
        try:
            #try to get taxa via API
            classification = _get_ncbi_lineage(str(ncbi_id))
            if classification is None:
                raise Exception(f"No lineage found for NCBI id {ncbi_id}.")
            return _sample_type_from_lineage(classification, cell_culture_key1, cell_culture_key2)
        except Exception as e:
            print(f"Could not get taxonomy info for {ncbi_id}: {e}")
            return [None, None]
    else:
        return [None, None]


# separators of an organism text that lists organisms ('A; B', 'A, B', 'A + B', 'A & B', 'A and B', 'A/B' between
# letters: not 'Mus musculus/C57' or a strain 'MHOM/IN/82'); the MetaboLights and Workbench converters
ORGANISM_LIST = re.compile(r'[;,+&]|\band\b|(?<=[A-Za-z]{2})\s*/\s*(?=[A-Za-z]{2})', re.I)


# culture types of one organism; a sample of organisms of several of them is a culture of several species
_CULTURE_SAMPLE_TYPES = frozenset({'culture_bacterial', 'culture_fungal', 'culture_mammalian'})
MULTIPLE_SPECIES_CULTURE = 'culture_multiplespecies'


def multi_sample_type(terms, taxonomy_info=None):
    """[SampleType, SampleTypeSub1] of a sample of several organisms (terms: their 'ID|Name' terms, e.g. the parts
    of an NCBITaxonomy value of several organisms), from each organism's get_taxonomy_info: culture_multiplespecies
    when all are cultures (a co-culture of two bacteria as of a bacterium and a fungus: a culture of several
    species, as microbeMASST's communities); else the type all give, with their SampleTypeSub1 when all give the
    same one (else 'missing value'); else [None, None] (an animal and a microbe, or an organism without a type: the
    sample type is not that of one of them)."""
    taxonomy_info = taxonomy_info or get_taxonomy_info
    infos = []
    for term in terms:
        taxid = single_taxid(term)
        info = taxonomy_info(taxid) if taxid else None
        if not info or info[0] is None:
            return [None, None]
        infos.append(info)
    if not infos:
        return [None, None]
    types = {i[0] for i in infos}
    if len(infos) > 1 and types <= _CULTURE_SAMPLE_TYPES:
        return [MULTIPLE_SPECIES_CULTURE, MULTIPLE_SPECIES_CULTURE]
    if len(types) == 1:
        subs = {i[1] for i in infos}
        return [types.pop(), subs.pop() if len(subs) == 1 else 'missing value']
    return [None, None]

# Vendor / marketing words that carry no model information; dropped before matching an
# instrument string against the allowed MassSpectrometer vocabulary. Shared by the MWB
# and MetaboLights converters.
_MS_STOPWORDS = {
    'thermo', 'abi', 'agilent', 'bruker', 'waters', 'sciex', 'leco', 'shimadzu', 'ab',
    'daltonics', 'hybrid', 'series', 'system', 'lc', 'ms', 'lcms', 'msms', 'scientific',
    'fisher', 'the', 'and', 'with', 'gc',
}


def _ms_tokens(s):
    s = str(s).lower().replace('q-exactive', 'q exactive').replace('qexactive', 'q exactive')
    s = re.sub(r'[^a-z0-9+]+', ' ', s)
    return {t for t in s.split() if t and t not in _MS_STOPWORDS}


_YES = {'yes', 'y', 'true'}
_NO = {'no', 'n', 'false'}


def convert_smoking(x, key=None):
    """Normalize a smoking value to the SmokingStatus vocabulary
    (current/former/never smoker, non-smoker). Only unambiguous text is mapped.

    Yes/no answers are only mapped when the key (column / factor name) says what was
    asked: 'Smoker: no' means not a current smoker ('non-smoker'), 'Ever smoked: no'
    means 'never smoker'. Numeric codes (0/1/2) stay missing because their polarity
    differs between studies, and keys about someone else's or passive exposure
    (prenatal, maternal, pregnancy, second-hand) are ignored. Separators are
    normalized to spaces first so word boundaries match (e.g. 'never_used')."""
    x = re.sub(r'[^a-z0-9]+', ' ', str(x).strip().lower()).strip()
    k = re.sub(r'[^a-z0-9]+', ' ', str(key).strip().lower()).strip() if key is not None else ''
    if re.search(r'prenatal|maternal|mother|pregnan|passive|second hand|secondhand|exposure|exposed', k):
        return 'missing value'
    if re.search(r'\bnever\b', x):
        return 'never smoker'
    if re.search(r'\bnon smoker\b|\bnonsmoker\b|\bnon smoking\b|\bnot smoking\b', x):
        return 'non-smoker'
    if re.search(r'\bformer\b|\bex smoker\b|\bprevious\b|\bpast\b|\bquit\b', x):
        return 'former smoker'
    if re.search(r'\bcurrent\b', x) or x == 'smoker':
        return 'current smoker'
    if k and (x in _YES or x in _NO):
        if re.search(r'\bever\b|histor', k):
            return 'never smoker' if x in _NO else 'missing value'
        if re.search(r'smok|tobacco|cigar', k):
            return 'current smoker' if x in _YES else 'non-smoker'
    return 'missing value'


def convert_diet(x):
    """Normalize a free-text diet value (from a diet-named field) to the Diet
    vocabulary. Human patterns and experimental/composition diet types; the input
    is assumed to come from a diet field, so short codes (cd/hf/hsd) decode safely.
    Bespoke study-specific formulas (e.g. 'ser/gly free') stay missing. Order is
    significant: most specific compositions first."""
    x = re.sub(r'[^a-z0-9]+', ' ', str(x).strip().lower()).strip()
    if not x:
        return 'missing value'
    # human dietary patterns
    if re.search(r'\bvegan\b', x):
        return 'vegan'
    if re.search(r'\bpesc[ae]tarian\b', x):
        return 'pescatarian'
    if re.search(r'\bvegetarian\b|lacto ovo|ovo vegetarian', x):
        return 'vegetarian'
    if re.search(r'\bomnivor', x):
        return 'omnivore'
    # infant feeding (breast milk / formula / solids; diet_columns.infant_feeding_term), before the
    # composition diets so that 'formula' is not read as one
    from diet_columns import infant_feeding_term
    infant = infant_feeding_term(x)
    if infant:
        return infant
    # choline-deficient diets (also with high fat: CDAHFD), before the high-fat diets
    if re.search(r'choline (deficien|deplet|free)|\bmcdd?\b|\bcdaa\b|\bcdahfd\b', x):
        return 'choline-deficient diet'
    if re.search(r'low fib(er|re)', x):
        return 'low-fibre diet'
    if re.search(r'high fib(er|re)|fib(er|re) (supplement|enrich|rich)', x):
        return 'fibre-supplemented diet'
    # experimental / composition diets (specific first)
    if re.search(r'\bhfhs\b|\bhchf\b|high fat high (sugar|sucrose|carb)', x):
        return 'high-fat high-sugar diet'
    if re.search(r'\bhfd\b|\bhf\b|high[ _]?fat|highfat', x):
        return 'high-fat diet'
    if re.search(r'\bhsd\b|high[ _]?sugar|high[ _]?sucrose|\bh sucrose\b|high carb', x):
        return 'high-sugar diet'
    if re.search(r'\blfd\b|low[ _]?fat', x):
        return 'low-fat diet'
    if re.search(r'high[ _]?protein|\bh protein\b', x):
        return 'high-protein diet'
    if re.search(r'ketogenic|\bketo\b', x):
        return 'ketogenic diet'
    if re.search(r'mediterranean', x):
        return 'Mediterranean diet'
    if re.search(r'western|regular american', x):
        return 'Western diet'
    if re.search(r'prebiotic', x):
        return 'prebiotic diet'
    if re.search(r'calorie restrict|caloric restrict|calorie restricted|\brestricted\b', x):
        return 'calorie-restricted diet'
    # control diets: purified (AIN-76 / AIN-93), chow, or not said
    if re.search(r'purified|\bain ?(76|93)', x):
        return 'purified diet'
    if re.search(r'\bchow\b|\bncd\b', x):
        return 'chow diet'
    if re.search(r'\bnormal\b|\bcontrol\b|\bstandard\b|\bcd\b|\bcon\b|\bctrl\b|\bsc\b|normal diet', x):
        return 'control diet'
    return 'missing value'


def packyears_to_numeric(x):
    """Cumulative smoking exposure in pack-years (packs per day x years smoked). Binned or
    threshold entries ('<70', '>=70') are rejected: they are not point values, and taking the
    bound would make up a number. Plausible range 0-200 (4 packs a day for 50 years); the same
    range is enforced centrally in complete_and_fill_REDU_table."""
    s = str(x)
    if re.search(r'[<>\u2264\u2265]|less than|greater than|more than|at least', s, re.I):
        return 'missing value'
    nums = re.findall(r'-?\d+\.?\d*', s)
    if not nums:
        return 'missing value'
    try:
        v = float(nums[0])
    except ValueError:
        return 'missing value'
    return v if 0 <= v <= 200 else 'missing value'


def bmi_to_numeric(x):
    """Extract a plausible BMI number (kg/m^2). Rejects z-scores, codes and other
    non-BMI values by requiring a human-plausible range (10-120: clinically routine BMIs above
    60 and documented extremes near 100 are kept); the same range is enforced centrally in
    complete_and_fill_REDU_table for all sources."""
    nums = re.findall(r'-?\d+\.?\d*', str(x))
    if not nums:
        return 'missing value'
    try:
        v = float(nums[0])
    except ValueError:
        return 'missing value'
    return v if 10 <= v <= 120 else 'missing value'


# Food and beverage products -> GFOP food ontology label (FoodOntologyTerm, labels of
# translation_sheets/GFOP_food_ontology.csv). Only names of things that exist as food or drink
# alone: milk, meat, fish, fruit, berry, grape, juice, oil, seed, leaf, tea leaf, chicken,
# yeast, insect ... also name body fluids, tissues or organisms of a study's subjects and are
# not listed (they count only for samples a study declares food, see food_ontology_term); eggs
# count only for samples of poultry (food_term_of_organism). A value matches as a whole or by
# its last words ('Aroeira honey'); None marks names ending in a listed word that are not that
# food. Animal products: not for plant-based analogues ('vegan cheese', _FOOD_ANALOGUE_WORDS).
_ANIMAL_FOOD_PRODUCTS = {
    # dairy
    'cheese': 'cheese', 'blue cheese': 'blue cheese', 'cheddar cheese': 'cheddar cheese', 'goat cheese': 'goat cheese',
    'cream cheese': 'cream cheese', 'feta cheese': 'feta cheese', 'gouda cheese': 'gouda cheese',
    'mozzarella cheese': 'mozzarella cheese', 'parmesan cheese': 'parmesan cheese', 'pecorino cheese': 'pecorino cheese',
    'ricotta cheese': 'ricotta cheese', 'brie cheese': 'brie cheese',
    'yogurt': 'yogurt', 'yoghurt': 'yogurt', 'kefir': 'dairy', 'fermented milk': 'dairy', 'buttermilk': 'dairy',
    'whey': 'dairy', 'butter': 'butter', 'cocoa butter': None, 'shea butter': None, 'body butter': None,
    'apple butter': None, 'mango butter': None, 'ice cream': 'ice cream', 'sour cream': 'sour cream',
    'infant formula': 'formula', 'baby formula': 'formula', 'formula milk': 'formula',
    'milk powder': 'milk', 'powdered milk': 'milk',
    # honey, eggs, meat
    'honey': 'honey', 'clover honey': 'honey_clover', 'sage honey': 'honey_sage', 'chestnut honey': 'honey_chestnut',
    'egg': 'egg', 'whole egg': 'egg', 'egg white': 'egg', 'egg yolk': 'egg', 'egg albumen': 'egg',
    'beef': 'beef', 'ground beef': 'beef', 'pork': 'pork', 'veal': 'veal', 'salami': 'salami', 'ham': 'ham',
    'prosciutto': 'prosciutto', 'bacon': 'bacon', 'sausage': 'meat', 'pork sausage': 'pork sausage',
    'beef jerky': 'beef jerkey',
}
FOOD_PRODUCT_TERMS = dict(_ANIMAL_FOOD_PRODUCTS, **{
    # beverages (wine: see _WINE_CULTIVARS)
    'beer': 'beer', 'wheat beer': 'wheat beer', 'wine': 'grape_fermented', 'red wine': 'red wine',
    'white wine': 'white wine', 'rose wine': 'rose wine, red', 'rosé wine': 'rose wine, red', 'port wine': 'port, red',
    'rice wine': 'grain_fermented', 'sake': 'grain_fermented', 'cider': 'apple_cider', 'pear cider': None,
    'kombucha': 'kombucha', 'coffee': 'coffee', 'coffee bean': 'coffee', 'tea': 'tea', 'hybrid tea': None,
    'green tea': 'Green tea', 'black tea': 'Black tea', 'white tea': 'White tea', 'oolong tea': 'Oolong tea',
    'pu erh tea': "Pu'er tea", "pu'er tea": "Pu'er tea", 'matcha': 'Matcha tea', 'yerba mate': 'yerba mate',
    'fruit juice': 'fruit', 'grape juice': 'grape', 'grape must': 'grape', 'orange juice': 'orange',
    'apple juice': 'apple', 'tomato juice': 'tomato',
    # other products
    'olive oil': 'olive', 'canola oil': 'canola oil', 'sesame oil': 'sesame oil',
    'tofu': 'soybean', 'tempeh': 'soybean', 'tempe': 'soybean', 'miso': 'soybean', 'natto': 'soybean',
    'soy sauce': 'soybean', 'soy milk': 'soybean', 'soymilk': 'soybean', 'sauerkraut': 'saurkraut',
    'tomato sauce': 'tomato sauce', 'ketchup': 'tomato sauce', 'tomato paste': 'tomato paste',
    'strawberry jam': 'strawberry jam', 'peanut butter': 'peanut butter', 'maple syrup': 'maple syrup',
    'molasses': 'molasses', 'balsamic vinegar': 'balsamic vinegar', 'potato chip': 'potato chip',
    'mashed potato': 'mashed potato', 'nori': 'nori',
})
_ANIMAL_FOOD_TERMS = set(_ANIMAL_FOOD_PRODUCTS.values()) - {None}
# words that make an animal product a plant-based analogue ('Beyond beef', 'almond butter')
_FOOD_ANALOGUE_WORDS = {'plant', 'vegan', 'vegetarian', 'veggie', 'vegetable', 'meatless', 'imitation', 'analog',
                        'analogue', 'substitute', 'mock', 'faux', 'beyond', 'impossible', 'soy', 'soya', 'almond',
                        'oat', 'coconut', 'cashew', 'nut', 'nuts', 'peanut', 'rice', 'pea', 'hemp'}
# grape cultivars -> GFOP label of their wine, else its colour. They narrow a wine ('Wine' and
# 'cv. Merlot' -> 'merlot, red'; 'Chardonnay wine' -> 'white wine'); in front of 'wine' only
# they, colours and styles (_WINE_WORDS) make it a grape wine, not 'plum wine' or 'palm wine'
_WINE_CULTIVARS = dict(
    {'cabernet sauvignon': 'cabernet sauvignon, red', 'merlot': 'merlot, red', 'pinot noir': 'pinot noir, red',
     'syrah': 'syrah, red', 'shiraz': 'syrah, red', 'petite sirah': 'petit sirah, red', 'petit sirah': 'petit sirah, red',
     'zinfandel': 'zinfindel, red', 'white zinfandel': 'zinfindel, white', 'sauvignon blanc': 'sauvignon blanc, white',
     'riesling': 'riesling, white'},
    **dict.fromkeys(['aglianico', 'barbera', 'cabernet franc', 'cannonau', 'carignan', 'carmenere', 'corvina', 'dolcetto',
                     'gamay', 'garnacha', 'grenache', 'lagrein', 'malbec', 'monastrell', 'montepulciano', 'mourvedre',
                     'nebbiolo', 'negroamaro', 'nerello mascalese', 'petit verdot', 'pinotage', 'primitivo', 'raboso',
                     'refosco', 'sagrantino', 'sangiovese', 'tannat', 'tempranillo', 'teroldego', 'touriga nacional'],
                    'red wine'),
    **dict.fromkeys(['albarino', 'arneis', 'chardonnay', 'chenin blanc', 'cortese', 'falanghina', 'fiano', 'friulano',
                     'garganega', 'gewurztraminer', 'glera', 'grenache blanc', 'gruner veltliner', 'italian riesling',
                     'pinot bianco', 'pinot blanc', 'pinot grigio', 'pinot gris', 'riesling italico', 'semillon',
                     'timorasso', 'tocai friulano', 'trebbiano', 'verdejo', 'verdicchio', 'vermentino', 'viognier',
                     'welschriesling'], 'white wine'))
_WINE_CULTIVAR_RE = re.compile(r'\b(?:%s)\b' % '|'.join(sorted(map(re.escape, _WINE_CULTIVARS), key=len, reverse=True)))
_WINE_WORDS = {'red', 'white', 'sparkling', 'still', 'dry', 'medium', 'semi', 'sweet', 'dessert', 'table', 'fortified',
               'commercial', 'young', 'aged', 'grape', 'and', 'vitis', 'vinifera', 'l', 'cv'}
_GRAPE_WINE_TERMS = {'grape_fermented', 'red wine', 'white wine', 'rose wine, red', 'port, red'} | set(_WINE_CULTIVARS.values())
# GFOP terms of a product by the genus of the sample's organism; eggs only of poultry ('egg'
# also names eggs of insects, worms and fish)
_FOOD_TERM_OF_GENUS = {
    'yogurt': {'Bos': 'yogurt_cow'}, 'milk': {'Bos': 'milk_cow', 'Capra': 'milk_goat'}, 'cream': {'Bos': 'cream_cow'},
    'meat': {'Bos': 'beef', 'Sus': 'pork'}, 'sausage': {'Sus': 'pork sausage'},
    'egg': {'Gallus': 'egg_chicken', 'Coturnix': 'egg_quail', 'Meleagris': 'egg_poultry', 'Anas': 'egg_poultry',
            'Anser': 'egg_poultry', 'Numida': 'egg_poultry'},
}
# words after a GFOP name that keep it that food, for declared food ('Pork Loin', 'Coffee neat')
_FOOD_PREPARATION_WORDS = {'raw', 'fresh', 'frozen', 'dried', 'cooked', 'boiled', 'fried', 'grilled', 'roasted', 'baked',
                           'steamed', 'smoked', 'canned', 'whole', 'sliced', 'diced', 'ground', 'minced', 'peeled', 'neat',
                           'diluted', 'brewed', 'loin', 'tenderloin', 'chop', 'chops', 'fillet', 'filet', 'breast',
                           'thigh', 'leg', 'wing', 'shoulder', 'belly', 'steak', 'mince'}
# a digest of a food (in vitro digestion, colonic fermentation) is no sample of the food
_FOOD_DIGEST_RE = re.compile(r'(?i)\bdigest(?:a|ate|ates|ed|ion|ions|s)?\b|'
                             r'\b(?:colonic|fa?ecal|gut)\s+(?:fermentation|incubation|inocul)')
# values that declare a sample a food (MetaboLights 'Food or Food Product', Workbench 'Sample
# type: Food'); trailing words dropped before matching
FOOD_SAMPLE_MARKERS = {'food', 'foods', 'food item', 'food sample', 'food or food product', 'foodstuff'}
_FOOD_VALUE_TAIL = re.compile(r'(?:\s+(?:food\s+)?(?:products?|samples?))+$')
# GFOP labels that name an organism or a whole branch, never the food of a sample
_GFOP_NOT_A_FOOD_NAME = {'plant', 'animal', 'fungi', 'algae', 'mineral', 'water', 'supplement', 'invertebrate',
                         'worm', 'reptile', 'lizard', 'insect', 'arthropod', 'gastropod', 'mouse', 'rat', 'cow',
                         'pig', 'horse', 'yeast', 'chlamydomonas reinhardtii'}


def _food_value_key(value):
    """Lower-case words of a value: separators as spaces, 'milk (fermented)' as 'fermented
    milk', 'yogurt food product' as 'yogurt'."""
    v = str(value).strip().lower()
    m = re.fullmatch(r'([^()]+?)\s*\(([^()]+)\)', v)
    if m:
        v = f'{m.group(2)} {m.group(1)}'
    v = re.sub(r'[\s_\-/,;:]+', ' ', v).strip()
    return _FOOD_VALUE_TAIL.sub('', v) if v not in FOOD_SAMPLE_MARKERS else v


def _food_lookup(words, terms, whole=False):
    """The longest last-words phrase in terms ('cocoa butter' before 'butter'; only all
    words if whole), also with a plural last word ('wines', 'blueberries'), or None."""
    for i in range(1 if whole else len(words)):
        tail = words[i:]
        last = tail[-1]
        for w in (last, last[:-3] + 'y' if last.endswith('ies') else None, last[:-1] if last.endswith('s') else None):
            phrase = ' '.join(tail[:-1] + [w]) if w else None
            if phrase in terms:
                return phrase
    return None


def _food_match(key, labels, prefix=False):
    """(phrase, term) of the food a value key names, or None: its last words in
    FOOD_PRODUCT_TERMS, else all of it in labels, else (prefix) its first words in either
    when the rest are _FOOD_PREPARATION_WORDS ('pork loin'). The term is None for analogues
    of animal products ('vegan cheese') and for wine after other words than _WINE_WORDS and
    grape cultivars ('plum wine')."""
    words = key.split(' ')
    phrase = _food_lookup(words, FOOD_PRODUCT_TERMS)
    label = _food_lookup(words, labels, whole=True) if phrase is None else None
    if phrase is not None:
        match = (phrase, FOOD_PRODUCT_TERMS[phrase], words[:len(words) - len(phrase.split(' '))])
    elif label is not None:
        match = (label, labels[label], [])
    else:
        match = next(((p, t[p], words[j:]) for j in range(len(words) - 1, 0, -1)
                      if prefix and _FOOD_PREPARATION_WORDS.issuperset(words[j:]) for t in (FOOD_PRODUCT_TERMS, labels)
                      for p in [_food_lookup(words[:j], t, whole=True)] if p is not None), None)
    if match is None:
        return None
    phrase, term, other = match
    if term in _ANIMAL_FOOD_TERMS and _FOOD_ANALOGUE_WORDS.intersection(other):
        term = None
    elif term in _GRAPE_WINE_TERMS and not _WINE_WORDS.issuperset(
            w.strip('.') for w in _WINE_CULTIVAR_RE.sub(' ', ' '.join(other)).split()):
        term = None
    return phrase, term


def _wine_of_cultivars(term, keys):
    """A grape wine term narrowed by the grape cultivars the keys name: one cultivar -> its
    wine ('cv. Merlot' -> 'merlot, red'), several of one colour -> that colour; never against
    the colour already known."""
    def colour(t):
        return 'red wine' if t == 'red wine' or t.endswith(', red') else 'white wine'
    found = {_WINE_CULTIVARS[m] for k in keys for m in _WINE_CULTIVAR_RE.findall(k)}
    found = {t for t in found if term == 'grape_fermented' or colour(t) == term}
    if len(found) == 1:
        return found.pop()
    colours = {colour(t) for t in found}
    return colours.pop() if term == 'grape_fermented' and len(colours) == 1 else term


def food_sample_marker(value):
    """Whether a value declares its sample a food ('Food or Food Product', 'Food')."""
    return isinstance(value, str) and _food_value_key(value) in FOOD_SAMPLE_MARKERS


def food_term_of_organism(term, organism='', sample_type='', phrase=None):
    """The FoodOntologyTerm for a sample of an organism (its NCBITaxonomy term) and SampleType:
    the product of the genus where GFOP has one ('yogurt' of Bos -> 'yogurt_cow', a 'sausage'
    (phrase) of Sus -> 'pork sausage'); None for eggs of other animals than poultry and for
    grape wine of a sample of another plant (fruit wines of Lycium samples)."""
    # (a sample of several organisms has no genus)
    genus = organism.split('|', 1)[1].strip().split(' ')[0] if isinstance(organism, str) and '|' in organism \
        and not is_multi(organism) else ''
    if term in _GRAPE_WINE_TERMS and sample_type == 'plant' and genus and genus != 'Vitis':
        return None
    by_genus = dict(_FOOD_TERM_OF_GENUS.get(term, {}), **_FOOD_TERM_OF_GENUS.get(phrase, {}))
    if term == 'egg' and genus not in by_genus:
        return None
    return by_genus.get(genus, term)


def food_ontology_term(values, gfop_labels, food_sample=False, organism=None, sample_type='', extra_values=(),
                       detail_values=()):
    """GFOP label (FoodOntologyTerm) of a sample from the values that describe it (most
    reliable first), or None. A value names a product of FOOD_PRODUCT_TERMS; for a sample its
    study declares food (food_sample, or one of the values is a FOOD_SAMPLE_MARKERS word) a
    GFOP label counts too, as the whole value ('Beef', 'Blueberries') or before preparation
    words ('Pork Loin', 'Coffee neat'), except labels of organisms. extra_values (e.g. other
    factors of a food sample) only count with a product name. A later value, also of
    detail_values (the sample's other characteristics and factors), can only narrow the
    product ('Wine', then 'red wine'), as can grape cultivars ('cv. Merlot'). A value naming a
    digest ('In Vitro Digestions') makes it no food sample. With an organism (NCBITaxonomy
    term, may be '') the term must fit it (food_term_of_organism)."""
    if any(isinstance(v, str) and _FOOD_DIGEST_RE.search(v) for v in (*values, *extra_values, *detail_values)):
        return None
    keys = [_food_value_key(v) for v in values if isinstance(v, str) and v.strip()]
    food_sample = food_sample or any(k in FOOD_SAMPLE_MARKERS for k in keys)
    labels = {str(l).lower(): l for l in gfop_labels if str(l).lower() not in _GFOP_NOT_A_FOOD_NAME} if food_sample else {}
    extra = [_food_value_key(v) for v in extra_values if isinstance(v, str) and v.strip()]
    detail = [_food_value_key(v) for v in detail_values if isinstance(v, str) and v.strip()]
    first = None
    for i, key in enumerate(keys + extra + detail):
        if not key or key in FOOD_SAMPLE_MARKERS or (first is None and i >= len(keys) + len(extra)):
            continue
        match = _food_match(key, labels if i < len(keys) else {}, prefix=food_sample and i < len(keys))
        if match is None:
            continue
        if first is None:
            first = match
        elif match[0].endswith(' ' + first[0]) and match[1] is not None:
            first = match
            break
    if first is None or first[1] is None:
        return None
    phrase, term = first
    if term in ('grape_fermented', 'red wine', 'white wine'):
        term = _wine_of_cultivars(term, keys + extra + detail)
    if organism is not None:
        term = food_term_of_organism(term, organism, sample_type, phrase)
    return term if term in gfop_labels else None


def map_instrument_to_allowed(instrument, allowed_values):
    """Map a free-text instrument name (e.g. "Thermo Scientific Q-Exactive") onto an
    allowed MassSpectrometer vocabulary entry ("Q Exactive|MS:1001911") by token-subset
    match: every (non-stopword) token of a vocab label must be present in the instrument
    string. The most specific label wins (most tokens), so "HF" never steals "HF-X".
    Returns the allowed value or None if nothing matches confidently."""
    it = _ms_tokens(instrument)
    if not it:
        return None
    best, best_n = None, 0
    for val in allowed_values:
        lab_tokens = _ms_tokens(str(val).split('|')[0])
        if lab_tokens and lab_tokens.issubset(it) and len(lab_tokens) > best_n:
            best, best_n = val, len(lab_tokens)
    return best


_TAXONOMY_TERM_INDEX = []   # [(allowed-values list, its length, by_id, by_name)]


def _taxonomy_term_index(taxonomy_data):
    """(ID -> first 'ID|Name' term, lower-cased name -> first 'ID|Name' term) over the
    allowed terms with exactly one '|', built once per allowed-terms list (the list is kept
    referenced with the index, so its id cannot be reused)."""
    if _TAXONOMY_TERM_INDEX and _TAXONOMY_TERM_INDEX[0][0] is taxonomy_data \
            and _TAXONOMY_TERM_INDEX[0][1] == len(taxonomy_data):
        return _TAXONOMY_TERM_INDEX[0][2], _TAXONOMY_TERM_INDEX[0][3]
    by_id, by_name = {}, {}
    for entry in taxonomy_data:
        parts = entry.split('|')
        if len(parts) == 2:
            ncbi_id, name = parts
            # the entry itself is 'ID|Name' (the same text as ncbi_id + '|' + name); keeping it
            # instead of a new string saves one string per allowed term (~3M, ~0.2 GB)
            by_id.setdefault(ncbi_id, entry)
            by_name.setdefault(name.lower(), entry)
    _TAXONOMY_TERM_INDEX[:] = [(taxonomy_data, len(taxonomy_data), by_id, by_name)]
    return by_id, by_name


def get_taxonomy_id_from_name__allowedTerms(organism_name, **kwargs):
    """First allowed 'ID|Name' NCBITaxonomy term with the ID popped from kwargs['ncbi_id']
    (a set), else the first one whose name is organism_name (case-insensitive), else None.
    Dictionary lookups; the same result as scanning the allowed terms in order."""

    allowedTerm_dict = kwargs['allowedTerm_dict']
    taxonomy_data = allowedTerm_dict["NCBITaxonomy"]["allowed_values"]
    
    ncbi_id_input = ''
    if 'ncbi_id' in kwargs.keys():
        ncbi_id_input = str(kwargs['ncbi_id'].pop())

    try:
        ncbi_id_numeric = int(ncbi_id_input)
    except ValueError:
        # ncbi_id_input cannot be made numeric, return None or an appropriate error message
        ncbi_id_input = ''

    by_id, by_name = _taxonomy_term_index(taxonomy_data)
    if ncbi_id_input != '':
        hit = by_id.get(str(ncbi_id_input))
        if hit is not None:
            return hit

    organism_name = str(organism_name)
    return by_name.get(organism_name.lower())



def get_uberon_table(owl_path):
    onto = load_ontology(owl_path)

    multi_cellular = onto.search(iri="http://purl.obolibrary.org/obo/UBERON_0010000")[0]
    organ = onto.search(iri="http://purl.obolibrary.org/obo/UBERON_0000062")[0]
    bodily_fluid = onto.search(iri="http://purl.obolibrary.org/obo/UBERON_0006314")[0]

    downstream_classes_multi_cellular = set(multi_cellular.descendants())
    downstream_classes_organ = set(organ.descendants())
    downstream_classes_bodily_fluid = set(bodily_fluid.descendants())
            
    data = []
    for cls in onto.classes():
        label = cls.label.first() if cls.label else None
        synonyms = [synonym for synonym in cls.hasExactSynonym] if hasattr(cls, 'hasExactSynonym') else []

        if not synonyms:
            synonyms = [label] if label else []
            
        is_multi_cellular = cls in downstream_classes_multi_cellular
        is_organ = cls in downstream_classes_organ
        is_bodily_fluid = cls in downstream_classes_bodily_fluid

        for synonym in synonyms:
            data.append({
                "UBERONOntologyIndex": cls.name,
                "Label": label,
                "Synonym": synonym,
                "Is Multicellular": is_multi_cellular,
                "Is Organ": is_organ,
                "Is Fluid": is_bodily_fluid
            })

    df = pd.DataFrame(data)

    df = df[(df["Label"].notna() | df["Synonym"].notna()) & ((df["Label"] != '') | (df["Synonym"] != ''))]
    df = df[df['UBERONOntologyIndex'].str.startswith('UBERON_', na=False)]
    # obsolete classes keep their label as 'obsolete <name>'; they must not be assigned
    df = df[~df['Label'].astype(str).str.lower().str.startswith('obsolete')]

    return df



def get_ontology_table(owl_path, ont_prefix, rm_synonym_info=False, descendant_node=None, index_column_name = 'UBERONOntologyIndex', imports_dir=None):
    onto = load_ontology(owl_path, imports_dir)

    filter_class = None
    descendants = set()
    if descendant_node:
        # Search for the specific class using its node ID
        filter_class = onto.search_one(iri="*" + descendant_node)
        if filter_class:
            # Get all its descendants
            descendants = set(filter_class.descendants())
            # Optionally remove the filter_class itself from the set of descendants to exclude it
            descendants.discard(filter_class)

    data = []
    for cls in tqdm.tqdm(onto.classes(), desc="Processing classes"):
        # Skip the class if it is the filter_class itself and descendant_node is provided
        if descendant_node and cls == filter_class:
            continue

        label = cls.label.first() if cls.label else None
        synonyms = [synonym for synonym in cls.hasExactSynonym] if hasattr(cls, 'hasExactSynonym') else []

        if not synonyms:
            synonyms = [label] if label else []

        # Determine if the current class is a descendant of the specified node
        is_descendant = cls in descendants

        for synonym in synonyms:
            data.append({
                index_column_name: cls.name,
                "Label": label,
                "Synonym": synonym,
                "Is Descendant": is_descendant  # Add this information to each entry
            })

    df = pd.DataFrame(data)
    df = df[(df["Label"].notna() | df["Synonym"].notna()) & ((df["Label"] != '') | (df["Synonym"] != ''))]
    df = df[df[index_column_name].str.startswith(ont_prefix, na=False)]
    # obsolete classes keep their label as 'obsolete <name>'; they must not be assigned
    df = df[~df['Label'].astype(str).str.lower().str.startswith('obsolete')]

    if rm_synonym_info:
        df['Synonym'] = df['Synonym'].str.replace(r' \([^)]*\)', '', regex=True)

    # Optionally, you might want to filter the dataframe to include only the descendants
    if descendant_node:
        df = df[df["Is Descendant"] == True]

    return df


# ---------------------------------------------------------------------------
# SampleCollectionDateandTime normalisation
# ---------------------------------------------------------------------------
# Submitters use a wide diversity of date formats and precisions: year only,
# year-month, full date, date+time, textual months, US M/D/Y, ISO Y-M-D, 2-digit
# years, ISO with timezone, and ranges. normalize_redu_datetime() accepts all of
# these and emits ISO-8601 at the SAME precision that was submitted (it never
# fabricates a day or a time that was not given):
#   2018 | 2018-07 | 2018-07-12 | 2018-07-12T22:00 | 2018-07-12T22:00:00
# and intervals as "<start>/<end>" (e.g. 2022-06-14/2022-06-30). Junk / missing
# markers / genuinely unparseable values return None so the caller can substitute
# the column's missing value.

_DATE_MONTHS = {}
for _i, _names in enumerate(
    [('jan', 'january'), ('feb', 'february'), ('mar', 'march'), ('apr', 'april'),
     ('may',), ('jun', 'june'), ('jul', 'july'), ('aug', 'august'),
     ('sep', 'sept', 'september'), ('oct', 'october'), ('nov', 'november'),
     ('dec', 'december')], start=1):
    for _n in _names:
        _DATE_MONTHS[_n] = _i

_DATE_MISSING_TOKENS = {
    '', 'nan', 'nat', 'none', 'null', 'na', 'n/a', 'missing value', 'not applicable',
    'not collected', 'not specified', 'not_applicable', 'ml import: not available',
    'unknown', 'unspecified', 'not available', 'not known', 'not provided',
    'not recorded', 'tbd',
}


def _date_yy(y):
    """Expand a 2-digit year: 00-68 -> 2000-2068, 69-99 -> 1969-1999."""
    y = int(y)
    if y < 100:
        y = 2000 + y if y <= 68 else 1900 + y
    return y


def _date_valid(y, m, d, H, M, S):
    try:
        datetime(y, m or 1, d or 1, H or 0, M or 0, S or 0)
        return 1900 <= y <= 2099
    except (ValueError, TypeError):
        return False


def _date_fmt(y, m, d, H, M, S, prec):
    if prec == 'year':
        return "%04d" % y
    if prec == 'month':
        return "%04d-%02d" % (y, m)
    if prec == 'day':
        return "%04d-%02d-%02d" % (y, m, d)
    if prec == 'minute':
        return "%04d-%02d-%02dT%02d:%02d" % (y, m, d, H, M)
    if prec == 'second':
        return "%04d-%02d-%02dT%02d:%02d:%02d" % (y, m, d, H, M, S)


def _parse_single_date(s, dayfirst=False):
    s = str(s).strip()
    low = s.lower()
    if low in _DATE_MISSING_TOKENS or not s or set(s) <= set('#?-/ .:'):
        return None

    # optional trailing time (':' or '.' separator), preceded by space/T/_, with an
    # optional ISO timezone designator. The time prefix is required so a bare
    # "04-04-2023" is never mistaken for a tz offset.
    H = M = S = None
    tm = re.search(r'[ T_](\d{1,2})[:.](\d{2})(?:[:.](\d{2}))?(?:\s*(?:Z|[+-]\d{2}:?\d{2}))?\s*$', s)
    if tm:
        H = int(tm.group(1)); M = int(tm.group(2))
        S = int(tm.group(3)) if tm.group(3) else None
        if H > 23 or M > 59 or (S is not None and S > 59):
            H = M = S = None
        else:
            s = s[:tm.start()].strip()

    y = m = d = None
    prec = None
    mo = None
    if re.fullmatch(r'\d{4}', s):
        y = int(s); prec = 'year'
    elif (mo := re.fullmatch(r'(\d{4})[-/.](\d{1,2})[-/.](\d{1,2})', s)):
        y, m, d = int(mo.group(1)), int(mo.group(2)), int(mo.group(3)); prec = 'day'
    elif (mo := re.fullmatch(r'(\d{4})[-/.](\d{1,2})', s)):
        y, m = int(mo.group(1)), int(mo.group(2)); prec = 'month'
    elif (mo := re.fullmatch(r'(\d{1,2})[/-](\d{1,2})[/-](\d{2,4})', s)):
        a, b, c = int(mo.group(1)), int(mo.group(2)), _date_yy(mo.group(3))
        if a > 12 and b <= 12:            # unambiguously D/M/Y
            d, m = a, b
        elif b > 12 and a <= 12:          # unambiguously M/D/Y
            m, d = a, b
        elif dayfirst:                     # ambiguous -> caller's derived preference
            d, m = a, b
        else:                              # ambiguous -> platform default (US M/D/Y)
            m, d = a, b
        y = c; prec = 'day'
    elif (mo := re.fullmatch(r'(\d{1,2})[/-](\d{4})', s)):
        m, y = int(mo.group(1)), int(mo.group(2)); prec = 'month'
    elif (mo := re.fullmatch(r'(\d{1,2})[-\s]([A-Za-z]{3,})[-\s,]+(\d{2,4})', s)):
        mon = _DATE_MONTHS.get(mo.group(2).lower()) or _DATE_MONTHS.get(mo.group(2).lower()[:3])
        if mon:
            d, m, y = int(mo.group(1)), mon, _date_yy(mo.group(3)); prec = 'day'
    elif (mo := re.fullmatch(r'([A-Za-z]{3,})[-\s,]+(\d{1,2})[-\s,]+(\d{2,4})', s)):
        mon = _DATE_MONTHS.get(mo.group(1).lower()) or _DATE_MONTHS.get(mo.group(1).lower()[:3])
        if mon:
            m, d, y = mon, int(mo.group(2)), _date_yy(mo.group(3)); prec = 'day'
    elif (mo := re.fullmatch(r'([A-Za-z]{3,})[-\s,]+(\d{2,4})', s)):
        mon = _DATE_MONTHS.get(mo.group(1).lower()) or _DATE_MONTHS.get(mo.group(1).lower()[:3])
        if mon:
            m, y = mon, _date_yy(mo.group(2)); prec = 'month'

    if y is None:
        return None
    if H is not None and prec == 'day':
        prec = 'second' if S is not None else 'minute'
    if not _date_valid(y, m, d, H or 0, M or 0, S or 0):
        return None
    return _date_fmt(y, m, d, H or 0, M or 0, S or 0, prec)


def normalize_redu_datetime(raw, dayfirst=False):
    """Normalise a submitted collection date/time to precision-preserving ISO-8601,
    or return None if it is missing / junk / unparseable.

    dayfirst disambiguates ONLY genuinely ambiguous slash dates (both fields <= 12);
    a value that itself proves an order (a field > 12) always wins over dayfirst.
    """
    if raw is None:
        return None
    s = str(raw).strip()
    if not s:
        return None
    parts = re.split(r'\s+-\s+|\s+to\s+', s)
    if len(parts) == 1:
        mo = re.fullmatch(r'(\d{1,2}/\d{1,2}/\d{2,4})\s*-\s*(\d{1,2}/\d{1,2}/\d{2,4})', s)
        if mo:
            parts = [mo.group(1), mo.group(2)]
    if len(parts) == 2:
        a, b = _parse_single_date(parts[0], dayfirst), _parse_single_date(parts[1], dayfirst)
        if a and b:
            return a + '/' + b
        return a or b
    return _parse_single_date(s, dayfirst)


if __name__ == '__main__':
    pass