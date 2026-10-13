"""Deterministic matching of free-text mass spectrometer names to the REDU
MassSpectrometer vocabulary ('<MS ontology label>|MS:xxxxxxx').

Repositories describe instruments as e.g. 'ABI Sciex 6500+', 'Thermo Orbitrap ID-X
Tribrid' or 'Agilent 6470 QQQ', while the vocabulary holds the PSI-MS labels
('Orbitrap ID-X', '6470 Triple Quadrupole LC/MS'). Names are reduced to their core
tokens (vendor names, generic words and punctuation removed, common abbreviations
expanded) and compared with the core tokens of every label. A term is only returned
when exactly one label matches; ambiguous names (e.g. 'Sciex 6500+', which may be a
QTRAP 6500+ or a Triple Quad 6500+) return None.

    match_mass_spectrometer(name, allowed_ms_terms, accession=None) -> 'Label|MS:xxxx' or None
"""
import re

_VENDOR_PHRASES = [
    'thermo fisher scientific', 'thermo scientific', 'thermo fisher', 'thermofisher', 'thermo finnigan',
    'thermo electron', 'thermo', 'finnigan',
    'ab sciex', 'abi sciex', 'ab sciex', 'applied biosystems', 'absciex', 'sciex', 'abi', 'ab',
    'agilent technologies', 'agilent',
    'waters corporation', 'waters',
    'bruker daltonics', 'bruker daltonik', 'bruker',
    'shimadzu', 'leco', 'jeol', 'varian', 'perkinelmer', 'perkin elmer', 'hitachi',
]

# words that do not identify a model
_GENERIC = {
    'lc', 'gc', 'ms', 'msms', 'mass', 'spectrometer', 'spectrometry', 'system', 'systems',
    'instrument', 'model', 'detector', 'hybrid', 'tribrid', 'uhplc', 'hplc', 'uplc', 'with',
    'and', 'the', 'esi', 'apci', 'msd', 'series', 'high', 'resolution', 'accurate', 'hr', 'hrms',
    'interfaced', 'coupled', 'to', 'a', 'an', 'of', 'mode', 'lcms', 'gcms',
}

# abbreviations -> canonical token sequences
_ALIASES = {
    'qqq': ['triple', 'quadrupole'],
    'tqms': ['triple', 'quadrupole'],
    'tqs': ['tqs'],
    'quad': ['quadrupole'],
    'qtof': ['qtof'],
    'q tof': ['qtof'],
    'triple quad': ['triple', 'quadrupole'],
    'tripletof': ['tripletof'],
    'triple tof': ['tripletof'],
    # 'Q-Exactive' / 'QExactive' -> the 'Q Exactive' labels
    'qexactive': ['q', 'exactive'],
}

_ACCESSION = re.compile(r'MS[:_](\d{7})', re.I)


def _core_tokens(text):
    s = str(text).lower()
    s = s.replace('/', ' ').replace(',', ' ').replace('(', ' ').replace(')', ' ')
    # 'GC-QTOF', 'LC-MS': the separation technique is its own word
    s = re.sub(r'(?<![a-z])(gc|lc)-', r'\1 ', s)
    # join hyphenated model parts: 'id-x' -> 'idx', 'hf-x' -> 'hfx', 'q-tof' -> 'qtof'
    s = re.sub(r'(?<=\w)-(?=\w)', '', s)
    s = re.sub(r'[^a-z0-9+ ]+', ' ', s)
    s = ' ' + re.sub(r'\s+', ' ', s).strip() + ' '
    for phrase in sorted(_VENDOR_PHRASES, key=len, reverse=True):
        s = s.replace(' ' + phrase + ' ', ' ')
    for phrase in sorted(_ALIASES, key=len, reverse=True):
        if ' ' in phrase:
            s = s.replace(' ' + phrase + ' ', ' ' + ' '.join(_ALIASES[phrase]) + ' ')
    tokens = []
    for t in s.split():
        tokens.extend(_ALIASES.get(t, [t]))
    return frozenset(t for t in tokens if t not in _GENERIC)


_index_cache = {}


def _index(allowed_ms_terms):
    key = id(allowed_ms_terms), len(allowed_ms_terms)
    if key in _index_cache:
        return _index_cache[key]
    by_id, by_label, by_core = {}, {}, {}
    for term in allowed_ms_terms:
        label, _, acc = str(term).partition('|')
        by_id[acc.strip().upper()] = term
        by_label.setdefault(label.strip().lower(), set()).add(term)
        if 'instrument model' in label.lower() or label.lower().endswith(' series'):
            continue  # vendor-level / series umbrella terms are only used on an exact label match
        core = _core_tokens(label)
        if core:
            by_core.setdefault(core, set()).add(term)
    _index_cache.clear()
    _index_cache[key] = (by_id, by_label, by_core)
    return _index_cache[key]


def match_mass_spectrometer(name, allowed_ms_terms, accession=None):
    """Returns the allowed term 'Label|MS:xxxxxxx' for an instrument name (and/or
    MS ontology accession), or None if it cannot be matched unambiguously."""
    by_id, by_label, by_core = _index(allowed_ms_terms)

    for candidate in (accession, name):
        if candidate is None:
            continue
        m = _ACCESSION.search(str(candidate))
        if m and ('MS:' + m.group(1)) in by_id:
            return by_id['MS:' + m.group(1)]

    if name is None or str(name).strip().lower() in ('', 'nan', 'none', 'missing value', 'na', 'n/a'):
        return None

    exact = by_label.get(str(name).strip().lower())
    if exact and len(exact) == 1:
        return next(iter(exact))

    core = _core_tokens(name)
    if not core:
        return None
    name_is_gc = bool(re.search(r'(?<![a-z])gc', str(name).lower()))
    same = by_core.get(core)
    if same and not name_is_gc:
        # 'gc' is a generic word, so 'Q Exactive GC Orbitrap' has the core of 'Orbitrap Q-Exactive'
        same = {t for t in same if not re.search(r'(?<![a-z])gc', t.lower())}
    if same:
        return next(iter(same)) if len(same) == 1 else None

    # Name carries extra words (e.g. 'Q Exactive Plus Orbitrap'): take the most specific
    # label whose tokens are all in the name, if it is unique and identifies a model.
    contained = [(len(c), {t for t in terms if name_is_gc or not re.search(r'(?<![a-z])gc', t.lower())})
                 for c, terms in by_core.items() if c < core]
    contained = [(n, terms) for n, terms in contained if terms]
    if not contained:
        return None
    best = max(n for n, _ in contained)
    best_terms = set().union(*[terms for n, terms in contained if n == best])
    if len(best_terms) != 1 or best < 2 and not any(ch.isdigit() for t in core for ch in t):
        return None
    term = next(iter(best_terms))
    # the label must contain the model number when the name has one
    name_numbers = {t for t in core if any(ch.isdigit() for ch in t)}
    label_core = _core_tokens(term.partition('|')[0])
    if name_numbers and not name_numbers <= label_core:
        return None
    return term
