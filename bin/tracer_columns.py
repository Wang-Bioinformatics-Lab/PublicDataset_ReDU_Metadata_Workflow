"""MetabolicLabelingTracer from the sample group texts (SampleGroupFreeText), for every row of the final table
(final_text_columns.py), and the tracer phrases of submitted sheets (read_and_validate_redu_from_github.py).

A tracer phrase is a compound of the pinned tracer table (tracer_compounds.tsv, built by build_tracer_compounds.py
from PubChem) written next to an isotope label: '[U-13C6]glucose', 'U-13C Glucose', '13C6-Glc', 'glucose-13C6',
'13C5,15N2-glutamine', 'glucose-d7', 'Uniformly labeled [13C]glucose', '18O2 air', 'D2O'. A word next to a label that
is no compound of the table is no tracer ('D8 infected', 'Washout D10', 'Cell line D17').
- The label must cover every atom of its element: 'U-' or a count equal to the compound's atoms (hydrogens: those
  on carbon; all where none is). Such a phrase becomes the isotopologue of tracer_isotopologues.tsv. A label
  without a count ('13C-glucose': which carbons is not said), with fewer atoms ('13C2-glucose', 'd9-choline') or
  with positions ('[1,2-13C2]glucose', 'amide-15N', 'methyl-d3') says too little for a SMILES: no value (counted).
- An enrichment only where it is unambiguous: '99 atom %', '99% 13C', '13C 99%' -> 'C{13C:99}' (a bare % can be a
  serum or a dose: '10% FBS').
- Several tracers in a value are joined by ';' in the order written.

The column (sample groups): filled on real samples (not blanks / QCs / standards) whose column is missing, from
the parts of their text with a tracer; in a factor of a dataset where some row has a tracer, an explicit
unlabelled value ('unlabeled', '12C-glucose', 'natural abundance'; under a tracer factor name also 'none' / 'no
tracer' / 'control') gives 'unlabeled'. microbiomeMASST rows are left alone. The factor counts as held where the
row's column has the value and the part says nothing more (a dose, time or other condition keeps it in the text).
"""
import collections
import csv
import os
import re

import pandas as pd

from sample_groups import MISSING, is_real_sample, masst_rows, parts_of

_HERE = os.path.dirname(os.path.abspath(__file__))
_ELEMENT_OF = {'13C': 'C', '15N': 'N', '18O': 'O', '34S': 'S', '2H': 'H'}
_ISO_ORDER = ['13C', '15N', '18O', '34S', '2H']

# one isotope of a label: optional positions ('1,2-', 'amide-', 'methyl-'), optional 'U-', isotope, count
_ISO = r'(?:13C|C13|15N|N15|18O|O18|34S|2H|D)'
_COUNT = r'(?:\d{1,2}(?![0-9]*(?:C|N|O|S|H)(?![a-z])))?'
_ITEM = (rf"(?:(?:\d+(?:'|′)?(?:\s*,\s*\d+(?:'|′)?)*|amide|amino|alpha|methyl|carboxyl|ring|ureido|guanido|guanidino|"
         rf"indole|imidazole|epsilon|alpha-amino)\s*-\s*)?(?:U\s*-?\s*)?{_ISO}{_COUNT}")
_LABEL = rf'(?:\[\s*{_ITEM}(?:\s*[,;/+]?\s*{_ITEM})*\s*\]|\(\s*{_ITEM}(?:\s*[,;/+]?\s*{_ITEM})*\s*\)|{_ITEM}(?:\s*[,/+]?\s*{_ITEM})*)'
_UNIFORM_WORDS = r'(?:uniformly[\s\-]*label+ed|uniformly|universally[\s\-]*label+ed)'
_LABEL_END = re.compile(rf'(?:{_UNIFORM_WORDS}\s*)?({_LABEL})[\s\-_]*(?:label+ed\s+)?$', re.I)
_LABEL_START = re.compile(rf'^[\s\-_]*\(?\s*({_LABEL})\s*\)?', re.I)
_ITEM_RE = re.compile(rf"(?P<pos>(?:\d+(?:'|′)?(?:\s*,\s*\d+(?:'|′)?)*|[a-z][a-z\-]*)\s*-\s*)?(?P<u>U\s*-?\s*)?"
                      rf"(?P<iso>{_ISO})(?P<n>{_COUNT})", re.I)
_ENRICH = [re.compile(r'(\d+(?:\.\d+)?)\s*(?:atom\s*%|%\s*atom|at\.?\s*%)\s*(13C|15N|18O|34S|2H|D)?', re.I),
           re.compile(r'(\d+(?:\.\d+)?)\s*%\s*(13C|15N|18O|34S|2H|D)\b', re.I),
           re.compile(r'\b(13C|15N|18O|34S|2H|D)\s*(\d+(?:\.\d+)?)\s*%', re.I)]
_HEAVY_WATER = re.compile(r'(?<![A-Za-z0-9])(?:D2O|heavy water|deuterium oxide|2H2O)(?![A-Za-z0-9])', re.I)
_UNLABELED = re.compile(r'(?<![a-z0-9])(?:un-?label+ed|non-?label+ed|not label+ed|natural(?:\s+isotop\w*)?\s+abundance|'
                        r'12C(?:\d+)?(?![0-9])|light)(?![a-z])', re.I)
_TRACER_KEY = re.compile(r'(?i)tracer|label|isotop|precursor')
_NO_TRACER_VALUE = re.compile(r'(?i)^\s*(?:none|no|no tracer|no label|control|ctrl|vehicle|unlabell?ed control)\s*$')
# what may stay of a part next to its tracers for the column to hold all it says
_FILLER = re.compile(r'(?i)\b(?:label+ed|label+ing|tracer|tracers|isotope|isotopic|stable|uniformly|universally|with|and|'
                     r'plus|media|medium|enriched|enrichment|atom|substrate|fed|feeding|treatment|treated|incubation)\b|'
                     r'[\s\-_+,;:/()\[\]]+')   # (numbers stay: a dose, time or bare % is something the column cannot hold)


def _load():
    comp, iso = {}, {}
    with open(os.path.join(_HERE, 'tracer_compounds.tsv')) as f:
        for r in csv.DictReader(f, delimiter='\t'):
            comp[r['compound']] = r
    with open(os.path.join(_HERE, 'tracer_isotopologues.tsv')) as f:
        for r in csv.DictReader(f, delimiter='\t'):
            iso[(r['compound'], r['labels'])] = r['smiles']
    names = {n: c for c, r in comp.items() for n in r['names'].split(';')}
    return comp, iso, names


COMPOUNDS, ISOTOPOLOGUES, NAMES = _load()
# longest names first: 'glycerol' before 'gly', 'l-glutamine' before 'glutamine'
_NAME_RE = re.compile(r'(?i)(' + '|'.join(re.escape(n) for n in sorted(NAMES, key=len, reverse=True)) + r')(?![a-z])')


def _iso(token):
    t = token.upper()
    return {'C13': '13C', 'N15': '15N', 'O18': '18O', 'D': '2H'}.get(t, t)


def _parse_label(label, uniform_words=False):
    """[(isotope, count or None, uniform, positional)] of a label text"""
    out = []
    for m in _ITEM_RE.finditer(label):
        pos = (m.group('pos') or '').strip(' -')
        iso = _iso(m.group('iso'))
        n = int(m.group('n')) if m.group('n') else None
        if iso == '2H' and m.group('iso') in ('D', 'd') and n is None:
            continue   # a bare 'D' is no label ('D-glucose', 'Day D')
        uniform = bool(m.group('u')) or uniform_words or pos.upper() == 'U'
        out.append((iso, n, uniform, bool(pos) and pos.upper() != 'U'))
    return out or None


def _labels_for(compound, items):
    """(sorted isotope labels, reason): reason is None when every label covers all atoms of its element"""
    c = COMPOUNDS[compound]
    labels = []
    for iso, n, uniform, positional in items:
        total = int(c['n_' + _ELEMENT_OF[iso]])
        if total == 0:
            return None, 'element_not_in_compound'
        if positional:
            return None, 'positional'
        if n is None and not uniform and total > 1:   # (one atom of the element: nothing to specify)
            return None, 'unspecified'
        if n is not None and n != total:
            return None, 'partial'
        labels.append(iso)
    return sorted(set(labels), key=_ISO_ORDER.index), None


def _enrichment(text, labels):
    """{isotope: %} said unambiguously in text, for the isotopes of labels"""
    out = {}
    for rx in _ENRICH:
        for m in rx.finditer(text):
            a, b = m.group(1), m.group(2)
            pct, iso = (a, b) if rx is not _ENRICH[2] else (b, a)
            iso = _iso(iso) if iso else (labels[0] if len(labels) == 1 and rx is _ENRICH[0] else None)
            if iso in labels and 0 < float(pct) <= 100:
                out.setdefault(iso, pct)
    return out


def tracers_in(text):
    """([tracer value 'SMILES[ El{iso:%}]', ...], Counter of reasons for tracer phrases that give no value, residual
    text: what the value says beyond its tracers)"""
    text = ' '.join(str(text).replace('‐', '-').replace('–', '-').split())
    found, reasons, spans = [], collections.Counter(), []
    if _HEAVY_WATER.search(text):
        m = _HEAVY_WATER.search(text)
        found.append(('water', ['2H'], m.span()))
    for m in _NAME_RE.finditer(text):
        compound = NAMES[m.group(1).lower()]
        left, right = text[:m.start()], text[m.end():]
        lm, rm = _LABEL_END.search(left), _LABEL_START.match(right)
        # the name must not run on from a word ('D-glucose' after 'Ham', 'Cellulose' after 'U')
        if not lm and not rm:
            continue
        if lm and re.search(r'[A-Za-z]$', left) and not re.search(r'(?:13C|15N|18O|34S|2H|\]|\)|\d)$', left, re.I):
            lm = None
        # a label on either side; where both are, the one that gives a value
        tried = []
        for side, lab in (('left', lm), ('right', rm)):
            if not lab:
                continue
            uniform_words = side == 'left' and bool(re.search(_UNIFORM_WORDS, left[lab.start():], re.I))
            items = _parse_label(lab.group(1), uniform_words)
            if items:
                span = (lab.start(), m.end()) if side == 'left' else (m.start(), m.end() + lab.end())
                tried.append(_labels_for(compound, items) + (span,))
        if not tried:
            continue
        labels, reason, span = next((t for t in tried if t[1] is None), tried[0])
        if reason:
            reasons[reason] += 1
            spans.append(span)
            continue
        found.append((compound, labels, span))
    values, seen = [], set()
    for compound, labels, span in sorted(found, key=lambda x: x[2][0]):
        smiles = ISOTOPOLOGUES[(compound, ';'.join(labels))]
        enrich = _enrichment(text, labels)
        braces = ''.join(f" {_ELEMENT_OF[i]}{{{i}:{enrich[i]}}}" for i in labels if i in enrich)
        v = smiles + braces
        if v not in seen:
            seen.add(v)
            values.append(v)
        spans.append(span)
    residual = text
    for a, b in sorted(spans, reverse=True):
        residual = residual[:a] + ' ' + residual[b:]
    for rx in _ENRICH:
        residual = rx.sub(' ', residual)
    return values, reasons, _FILLER.sub(' ', residual).strip()


def tracer_value(text):
    """The MetabolicLabelingTracer value of a tracer phrase ('[U-13C6]glucose 99 atom %'), or None"""
    values, _, _ = tracers_in(text)
    return ';'.join(values) if values else None


def apply(df):
    """Returns (df with MetabolicLabelingTracer filled, per-factor table for sample_groups.remove_held_factors,
    stats)."""
    df = df.copy()
    col = 'MetabolicLabelingTracer'
    if col not in df.columns:
        df[col] = MISSING
    df[col] = df[col].fillna(MISSING).replace('', MISSING)
    ds = df['ATTRIBUTE_DatasetAccession'].astype(str)
    stats = collections.Counter()
    recs, cache = [], {}
    for i, text, m in zip(df.index, df['SampleGroupFreeText'], masst_rows(df)):
        if m:
            continue
        for key, value, _, part in parts_of(text):
            ck = (key, value)
            if ck not in cache:
                vals, reasons, residual = tracers_in(value)
                if not vals and _TRACER_KEY.search(key) is None:
                    # '[13C2] Acetate: 200 uM': the tracer is the factor name
                    kv, kr, _ = tracers_in(key)
                    if kv:
                        vals, reasons, residual = kv, kr, value
                    else:
                        reasons = reasons + kr
                unl = not vals and bool(_UNLABELED.search(value) or (_TRACER_KEY.search(key) and _NO_TRACER_VALUE.match(value)))
                if unl:
                    residual = _FILLER.sub(' ', _NAME_RE.sub(' ', _UNLABELED.sub(' ', _NO_TRACER_VALUE.sub(' ', value)))).strip()
                cache[ck] = (';'.join(vals) if vals else None, unl, reasons, residual)
            tracer, unl, reasons, residual = cache[ck]
            for r, n in reasons.items():
                stats[f'phrases_without_value_{r}'] += n
            if tracer or unl:
                recs.append((i, key, value, part, tracer, unl, residual))
    fac = pd.DataFrame(recs, columns=['row', 'key', 'value', 'part', 'tracer', 'unlabeled', 'residual'])
    fac['dataset'] = ds.reindex(fac['row']).to_numpy()
    fac['real'] = df['SampleType'].map(is_real_sample).reindex(fac['row']).to_numpy().astype(bool) if len(fac) else []
    # unlabelled values only count in a factor of a dataset with a tracer
    with_tracer = set(fac.loc[fac['tracer'].notna(), ['dataset', 'key']].itertuples(index=False, name=None))
    fac = fac[fac['tracer'].notna() | [(d, k) in with_tracer for d, k in zip(fac['dataset'], fac['key'])]].copy()

    new = pd.Series(None, index=df.index, dtype=object)
    for i, g in fac[fac['real']].groupby('row', sort=True):
        tr = [t for t in g['tracer'] if t]
        if tr:
            new[i] = ';'.join(dict.fromkeys(';'.join(tr).split(';')))
        elif g['unlabeled'].any():
            new[i] = 'unlabeled'
    fill = (df[col] == MISSING) & new.notna()
    df.loc[fill, col] = new[fill]
    stats['MetabolicLabelingTracer_filled'] = int(fill.sum())
    stats['MetabolicLabelingTracer_filled_unlabeled'] = int((fill & (new == 'unlabeled')).sum())
    stats['MetabolicLabelingTracer_datasets'] = int(ds[fill & (new != 'unlabeled')].nunique())
    stats['MetabolicLabelingTracer_conflicts_(column_kept)'] = int(((df[col] != MISSING) & ~fill & new.notna()
                                                                   & (df[col] != new)).sum())
    row_v = df[col].astype(str).reindex(fac['row']).to_numpy() if len(fac) else []
    fac['captured'] = [(t is not None and rv == t) or (t is None and u and rv == 'unlabeled')
                       for t, u, rv in zip(fac['tracer'], fac['unlabeled'], row_v)]
    fac['held'] = fac['captured'] & (fac['residual'] == '')
    return df, fac, stats
