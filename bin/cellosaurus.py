"""Cell lines (Cellosaurus) for ReDU: the shipped table, exact name lookups and the
deterministic matcher used by the converters.

The table is data/cellosaurus_<version>_cell_lines.tsv.gz, built by
build_cellosaurus_table.py (one row per name/synonym of a cell line). CellLine is the
Cellosaurus recommended name, CellLineIndex its accession (CVCL_xxxx).
"""
import functools
import re

import pandas as pd

# values that never name a cell line (missing values, genotypes, treatments, plate wells)
_MISSING_LIKE = {'', '-', '--', 'na', 'n/a', 'n.a.', 'nan', 'none', 'null', 'unknown', 'not applicable',
                 'not available', 'not specified', 'not reported', 'not collected', 'missing value', 'other',
                 'yes', 'no', 'control', 'controls', 'blank', 'qc', 'pool', 'pooled', 'media', 'medium',
                 'see above', 'parental', 'vehicle', 'dmso', 'pbs', 'lps', 'treated', 'untreated', 'normal',
                 'mock', 'sample', 'samples', 'cell', 'cells', 'cell line', 'cell lines', 'cultured cells'}
_GENOTYPE_RE = re.compile(r'^(?:wt|ko|kd|oe|del|dko|tko|het|hom|null|mut|mutant|wild[- ]?type|knock[- ]?(?:out|down)|'
                          r'[+-]/[+-]|ctrl|(?:sh|si)-?rna\S*|(?:sh|si)(?:ctrl|nt|scr|neg|[0-9]+)|empty vector|ev|vector|scramble[d]?|'
                          r'(?:wt|ko|kd|oe)[0-9]*|[a-z0-9]{1,6}[- ]?(?:ko|kd|oe|wt|del)|(?:ko|kd|oe|del)[- ]?[a-z0-9]{1,8})$', re.I)
# 96-well plate positions (A1-H12); 'B16' etc. stay possible line names
_WELL_RE = re.compile(r'^[a-h]\s*(?:0?[1-9]|1[0-2])$', re.I)
# mouse/rat strains: names like 'C57BL/6J' are synonyms of some Cellosaurus lines
_STRAIN_RE = re.compile(r'c57|\bbl/?6|balb|\bfvb\b|\bdba\b|\b129s|\bnod\b|\bscid\b|sprague|wistar|\bnude\b|'
                        r'\bcd-?1\b|swiss|long[- ]evans', re.I)
# generic cell-type abbreviations that name primary cells, not a specific line
_GENERIC_ABBREVIATIONS = {'mef', 'mefs', 'nk', 'pbmc', 'pbmcs', 'hsc', 'hscs', 'msc', 'mscs', 'ipsc', 'ipscs', 'esc',
                          'escs', 'bmdm', 'bmdms', 'rbc', 'rbcs', 'nsc', 'nscs', 'npc', 'npcs', 'treg', 'tregs', 'th17',
                          'cd4', 'cd8', 'huvec', 'huvecs', 'hmec', 'mep', 'ips', 'es', 'hff', 'hffs', 'mel'}
_ORGANISM_WORDS = r'(?:human|mouse|murine|rat|hamster|canine|dog|monkey|porcine|pig|bovine|chicken|primary)'
_VENDOR_RE = re.compile(r'\b(?:atcc|dsmz|ecacc|jcrb|riken|cls|sigma|rrid)\b[\s#:no.]*[a-z]*[- ]?\d+[\w.-]*', re.I)
_CATALOG_RE = re.compile(r'\b(?:crl|htb|ccl|acc|tib|hb)[- ]?\d+(?:\.\d+)?\b', re.I)
_CELL_WORDS_RE = re.compile(r'\bcell[ -]?lines?\b|\bcells?\b|\bline\b', re.I)
_CVCL_RE = re.compile(r'cvcl_[0-9a-z]{4}', re.I)


def load_cellosaurus_table(path):
    """The shipped table as a DataFrame of strings, or None for path 'none'/None."""
    if path is None or str(path).lower() == 'none':
        return None
    return pd.read_csv(path, sep='\t', dtype=str, keep_default_na=False)


def _squash(text):
    return re.sub(r'[^a-z0-9]', '', str(text).lower())


def normalize_value(value):
    """Tier (c) normalisation: no '(...)' / '[...]', vendor and catalogue numbers, the
    words 'cell(s)', 'cell line', 'line' and a leading organism word; then only letters
    and digits ('Hep G2 cells (ATCC HB-8065)' -> 'hepg2')."""
    s = str(value).lower().replace('®', '').replace('™', '')
    s = re.sub(r'\([^()]*\)|\[[^\[\]]*\]', ' ', s)
    s = _VENDOR_RE.sub(' ', s)
    s = _CATALOG_RE.sub(' ', s)
    s = _CELL_WORDS_RE.sub(' ', s)
    s = re.sub(r'^\s*' + _ORGANISM_WORDS + r'\b', ' ', s.strip())
    return _squash(s)


def is_rejected_value(value):
    """True for values that are no cell-line name: missing-like words, genotypes and
    treatments (WT, KO, DEL, shRNA), plate wells (D5, G6), mouse strains, numbers."""
    v = ' '.join(str(value).strip().split())
    low = v.lower()
    if low in _MISSING_LIKE or low.rstrip('s') in _GENERIC_ABBREVIATIONS or low in _GENERIC_ABBREVIATIONS:
        return True
    if _GENOTYPE_RE.match(low) or _WELL_RE.match(low) or _STRAIN_RE.search(low):
        return True
    if re.fullmatch(r'[\d\s.,;:%/+-]*', low):
        return True
    return False


class CellosaurusIndex:
    """Lookups over the shipped table."""

    def __init__(self, table):
        table = table.fillna('')
        first = table.drop_duplicates('CellLineIndex')
        self.info = {}
        for ac, name, taxids, species, category, mammal, problem in zip(
                first['CellLineIndex'], first['CellLine'], first['NCBITaxonIDs'],
                first.get('Species', pd.Series('', index=first.index)), first['Category'],
                first.get('Mammalian', pd.Series('', index=first.index)), first.get('Problematic', pd.Series('', index=first.index))):
            self.info[ac] = {'name': name, 'taxids': tuple(t for t in taxids.split('|') if t),
                             'species': tuple(s for s in species.split('|') if s), 'category': category,
                             'mammal': mammal == '1', 'problematic': problem}
        self.by_name = {}          # recommended name (exact) -> accession
        self.by_id = {}            # lower-case recommended name -> accessions
        self.by_synonym = {}       # lower-case synonym -> accessions (recommended names excluded)
        self.by_norm = {}          # normalised name or synonym -> accessions
        for ac, info in self.info.items():
            self.by_name.setdefault(info['name'], ac)
            self.by_id.setdefault(info['name'].lower(), set()).add(ac)
        for ac, name, syn in zip(table['CellLineIndex'], table['CellLine'], table['Synonym']):
            if syn != name:
                self.by_synonym.setdefault(syn.lower(), set()).add(ac)
            n = _squash(syn)
            if n:
                self.by_norm.setdefault(n, set()).add(ac)
        self.accessions = set(self.info)

    # ---- exact lookups (validation of submitted CellLine values) ----
    def canonical(self, value):
        """(recommended name, accession) for a submitted CellLine value: the recommended
        name itself, its accession (CVCL_xxxx), the recommended name in other case, or a
        synonym of exactly one line; else None."""
        v = ' '.join(str(value).strip().split())
        if not v:
            return None
        if v in self.by_name:
            return v, self.by_name[v]
        m = re.fullmatch(r'(?i)cvcl_([0-9a-z]{4})', v)
        if m:
            ac = 'CVCL_' + m.group(1).upper()
            return (self.info[ac]['name'], ac) if ac in self.info else None
        for lookup in (self.by_id, self.by_synonym):
            hits = lookup.get(v.lower(), set())
            if len(hits) == 1:
                ac = next(iter(hits))
                return self.info[ac]['name'], ac
            if hits:
                return None
        return None

    # ---- matching free text (converters) ----
    def candidates(self, value):
        """[(tier, accessions)] for a free-text value: (a) accession or recommended name,
        (b) synonym, (c) normalised name. Several accessions in one value -> 'multi'."""
        v = ' '.join(str(value).strip().split())
        low = v.lower()
        cvcl = {m.upper() for m in _CVCL_RE.findall(v)}
        if len(cvcl) > 1:
            return [('multi', cvcl)]
        if cvcl:
            ac = 'CVCL_' + next(iter(cvcl))[5:]
            return [('a', {ac})] if ac in self.info else []
        out = []
        if low in self.by_id:
            out.append(('a', set(self.by_id[low])))
        if low in self.by_synonym:
            out.append(('b', set(self.by_synonym[low])))
        n = normalize_value(v)
        if len(n) >= 2 and n in self.by_norm:
            out.append(('c', set(self.by_norm[n])))
        return out

    def species_consistent(self, ac, study_taxid, study_name, lineage_fn=None):
        """True if a line comes from the study organism: the same NCBI taxon, the same
        species name, a subspecies/strain of the line's species, or the genus/species
        just above a line of a subspecies."""
        info = self.info[ac]
        sid = str(study_taxid or '').strip()
        sname = str(study_name or '').strip().lower()
        if sid and sid in info['taxids']:
            return True
        species = [s.lower() for s in info['species']]
        if sname and sname in species:
            return True
        if lineage_fn is None:
            return False
        try:
            if sid:
                lin = lineage_fn(sid) or ()
                if any(s in lin for s in species):
                    return True
            for t in info['taxids']:
                lin = lineage_fn(t) or ()
                if sname and sname in lin[-2:]:
                    return True
        except Exception:
            return False
        return False

    def match(self, value, study_taxid='', study_name='', cultured=False, lineage_fn=None):
        """(accession or None, reason) for a free-text value of a sample or study.

        Tiers (a)/(b) first, (c) only when they give no line of the study organism. On
        ambiguity the study organism decides (also a synonym of the study organism over an
        exact name of another one); more than one candidate left, or a species conflict,
        gives None. Values of <= 3 characters only in cultured-cell studies and only when
        unique. A human line in a mouse study (xenograft) only in cultured-cell studies,
        and only when it is the one candidate."""
        v = ' '.join(str(value).strip().split())
        if is_rejected_value(v):
            return None, 'not a cell line value'
        cands = self.candidates(v)
        if not cands:
            return None, 'no match'
        if cands[0][0] == 'multi':
            return None, 'several cell lines'
        short = len(v) <= 3
        if short and not cultured:
            return None, 'short name outside a cultured-cell study'
        has_species = bool(str(study_taxid or '').strip() or str(study_name or '').strip())
        ab = [(t, a) for t, a in cands if t in ('a', 'b')]
        c = [(t, a) for t, a in cands if t == 'c']
        reason = 'no match'
        for group in ([g for g in (ab,) if g] + [g for g in (c,) if g]):
            pool = set().union(*[a for _, a in group])
            exact = set().union(*[a for t, a in group if t == 'a']) if any(t == 'a' for t, _ in group) else set()
            if short and len(pool) > 1:
                return None, 'short name, not unique'
            if not has_species:
                if cultured and len(pool) == 1:
                    return next(iter(pool)), 'unique, no study organism'
                reason = 'no study organism'
                continue
            same = {ac for ac in pool if self.species_consistent(ac, study_taxid, study_name, lineage_fn)}
            if len(same) == 1:
                return next(iter(same)), 'species'
            if len(same) > 1:
                if len(same & exact) == 1:
                    return next(iter(same & exact)), 'species, exact name'
                return None, 'ambiguous within the study organism'
            # no line of the study organism; a human line in a cultured-cell study of mice
            # (often the study organism of the culture's other component)
            mouse_study = str(study_taxid).strip() == '10090' or str(study_name).strip().lower() == 'mus musculus'
            if cultured and mouse_study and len(pool) == 1 and self.info[next(iter(pool))]['taxids'] == ('9606',):
                return next(iter(pool)), 'human line in a cultured-cell study of mice'
            reason = 'species conflict'
        return None, reason
