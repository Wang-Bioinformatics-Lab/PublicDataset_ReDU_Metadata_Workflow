"""InternalStandardsUsed from the sample group texts (SampleGroupFreeText), and which internal-standard
and exposure factors of the texts InternalStandardsUsed / Exposure hold, for every row of the final
table (final_text_columns.py, after diet_columns.apply).

InternalStandardsUsed: only a factor whose name says it is the internal standard ('Internal
standard: caffeine', 'Internal standards used: ...', 'ISTD: ...'); never the false friends:
Workbench 'IR_IS' (insulin resistant / sensitive), spike-ins of analytes or cells ('Spiked Cell
Number', 'Spiking state'), isotope tracers ('Treatment: deuterated methionine'), calibration-curve
samples, or samples that are themselves the standard solution ('Treatment: Internal standard',
'Group: ISTD', 'Blank_IS': no compound named).
- The value must name an allowed term, case, spaces and '_' aside; several names ('a, b', 'a + b',
  'a and b') give the allowed combination of exactly these names (any order), else nothing.
- An internal standard is a property of the analysis: blanks, QCs and standards carrying the
  factor get it too (they are spiked as well), and a dataset whose standard differs per assay gets
  it per row.
- A filled column is kept. The factor counts as held where the column's names include the
  factor's names; a different value keeps the factor in the text (counted as a conflict).

Exposure: never filled from the text here. The converters decide what the sampled individual was
exposed to, with context the final table lacks (study titles, culture columns, pre-exposure time
points, PECO design descriptors). A factor about a treatment / exposure ('Pesticide Exposure',
'Treatment', 'Compound', 'Drug', ...; on GNPS rows every microbiomeMASST intervention, which
masst_tables.py matches to ECTO 'exposure to <value>') is held where the row's Exposure is the ECTO /
PECO term its value names (label or exact synonym, case and separators aside, also with 'treated
with' / 'exposure' wording); a value that also carries a dose, concentration or time is not held
(the column cannot hold the dose).
"""
import collections
import re

import pandas as pd

from sample_groups import MISSING, is_real_sample, masst_rows, norm, parts_of


def squash(text):
    """lower case, alphanumerics only ('(2,4-dichlorophenoxy)acetate' -> '24dichlorophenoxyacetate');
    optical-rotation signs stay ('(-)-methyl jasmonate' is not 'methyl jasmonate')"""
    t = str(text).lower().replace('(+)', ' plus ').replace('(-)', ' minus ').replace('(±)', ' pm ')
    return re.sub(r'[^a-z0-9]+', '', t)


# ---------------------------------------------------------------- internal standards

# factor names that say the factor is the internal standard
_IS_KEY = re.compile(r'^(?:internal standards?|internal std|internal standards? (?:used|added|compounds?|names?|mix|'
                     r'mixture|solution)|istds?|is|is compounds?)$')
# ... of these only with a name ('IS: none' is no statement about an internal standard)
_IS_KEY_SHORT = {'is', 'istd', 'istds'}
_IS_SPLIT = re.compile(r'\s*(?:;|\s\+\s|,\s+|\s+and\s+|\s*&\s*)\s*')


def is_internal_standard_key(key):
    return bool(_IS_KEY.match(norm(key)))


def _is_key(name):
    return ''.join(str(name).lower().replace('_', ' ').split())


class InternalStandardIndex:
    """Allowed InternalStandardsUsed terms by their name key (case, spaces and '_' aside); the
    combinations ('a;b') by their set of names."""

    def __init__(self, allowed_terms):
        by_key = collections.defaultdict(set)
        for t in allowed_terms:
            by_key[_is_key(t)].add(t)
        # (a key two terms share is ambiguous)
        self.by_key = {k: next(iter(v)) for k, v in by_key.items() if len(v) == 1}
        combos = collections.defaultdict(set)
        for t in allowed_terms:
            if ';' in t:
                combos[frozenset(_is_key(n) for n in t.split(';'))].add(t)
        self.combos = {k: next(iter(v)) for k, v in combos.items() if len(v) == 1}

    def term(self, value, key=''):
        """allowed term or None for one internal-standard factor value"""
        v = ' '.join(str(value).split()).strip(' ;,')
        if not v:
            return None
        hit = self.by_key.get(_is_key(v))
        if hit and not (hit == 'none' and norm(key) in _IS_KEY_SHORT):
            return hit
        names = [p for p in _IS_SPLIT.split(v) if p]
        if len(names) > 1:
            return self.combos.get(frozenset(_is_key(n) for n in names))
        return None

    @staticmethod
    def names(term):
        return {_is_key(n) for n in str(term).split(';')}


# ---------------------------------------------------------------- exposures (held test)

# factor names about a treatment / exposure of the sampled organism
_EXP_KEY = re.compile(r'(?<![a-z])(?:expos\w*|treat\w*|compounds?|chemicals?|drugs?|toxi\w*|pollut\w*|contamina\w*|'
                      r'pesticides?|herbicides?|fungicides?|insecticides?|xenobiotics?|agents?|stressors?|inhibitors?|'
                      r'antibiotics?|medications?|medicines?|substances?|interventions?|administ\w*|dosing|gavage|'
                      r'injections?|challenge|stimul\w*|metals?|chemotherap\w*|therap\w*)(?![a-z])')
# ... but not these (own columns, analysis, vehicles, cultures, the parents, doses and times)
_EXP_KEY_SKIP = re.compile(r'(?<![a-z])(?:diet\w*|feed\w*|food|chow|nutri\w*|tracers?|label\w*|isotop\w*|internal|standards?|'
                           r'spik\w*|calibrat\w*|solvents?|vehicles?|media|medium|culture\w*|cells?|extraction|'
                           r'maternal|mother|paternal|father|parental|prenatal|in utero|offspring|dose|doses|dosage|'
                           r'concentrations?|amount|time|duration|days?|hours?|weeks?|units?|route|frequency)(?![a-z])')
# explicit controls / none values (count as held where the row has an Exposure)
_NONE_VALUE = re.compile(r'^(?:control|controls|control group|ctrl|ctl|con|cont|vehicle|veh|untreated|non treated|'
                         r'not treated|none|no|no treatment|unexposed|non exposed|not exposed|no exposure|sham|mock|'
                         r'placebo|naive|dmso|pbs|saline|water|solvent|solvent control|vehicle control|negative control|'
                         r'(?:control|vehicle) \(?(?:dmso|pbs|saline|water|methanol|ethanol)\)?)$')
# doses, concentrations and times in a value ('cisplatin 10 uM', '5 mg/kg')
_DOSE = re.compile(r'(?<![a-z0-9.])\d+(?:[.,]\d+)?\s*(?:mg/kg(?:/d(?:ay)?)?|g/kg|mg|µg|ug|ng|g/l|mg/l|mg/ml|ug/ml|µg/ml|'
                   r'µm|um|nm|mm|m|mmol/l|µmol/l|umol/l|micromolar|millimolar|nanomolar|ppm|ppb|%|h|hr|hrs|hours?|d|'
                   r'days?|weeks?|wks?|min)(?![a-z0-9])', re.I)
_AFFIX = re.compile(r'^(?:treated with|treatment with|treatment|exposed to|exposure to|exposure|administration of|'
                    r'injected with|injection of|given)\s+|\s+(?:treated|treatment|exposed|exposure|group|administration|'
                    r'injection|injected)$')
# ECTO terms too general to say what the organism was exposed to, or vehicles / materials
_GENERIC = {'exposure to water', 'exposure to solvent', 'exposure to dimethyl sulfoxide', 'exposure to drug',
            'exposure to chemical', 'exposure to chemical entity', 'exposure to molecular entity', 'exposure to lipid',
            'exposure to peptide', 'exposure to amino acid', 'exposure to polysaccharide', 'exposure to acids',
            'exposure to acid', 'exposure to carbohydrate', 'exposure to organic compound', 'exposure to treatment',
            'exposure to material entity', 'exposure to chemical with biological effect', 'exposure to pesticide',
            'exposure to herbicide', 'exposure to fungicide', 'exposure to insecticide', 'exposure to agrochemical',
            'exposure to antibiotic', 'exposure to toxin', 'exposure to metal', 'exposure to salt', 'exposure to oil',
            'exposure to inhibitor', 'exposure to air', 'chemical exposure'}


def is_exposure_key(key):
    k = norm(key)
    return bool(_EXP_KEY.search(k)) and not _EXP_KEY_SKIP.search(k)


def is_none_value(value):
    return bool(_NONE_VALUE.match(norm(value)))


class EctoIndex:
    """ECTO / PECO label and exact-synonym lookups from the prepared exposure table
    (prepare_ontologies.py: ExposureIndex, Label, Synonym), allowed labels only; a name must give
    one term."""

    def __init__(self, exposure_table, allowed_labels):
        allowed = set(allowed_labels) if allowed_labels else None
        exact, squashed = collections.defaultdict(set), collections.defaultdict(set)
        if exposure_table is not None and len(exposure_table):
            e = exposure_table[exposure_table['ExposureIndex'].astype(str).str.startswith(('ECTO:', 'PECO:'))]
            for label, syn in e[['Label', 'Synonym']].astype(str).itertuples(index=False):
                low = label.lower()
                if (allowed is not None and label not in allowed) or low in _GENERIC or low.startswith('obsolete') or \
                        re.search(r'\b(?:maternal|paternal|prenatal|via \w+|diet|food|feeding)\b', low):
                    continue
                for name in (label, syn):
                    exact[' '.join(name.lower().replace('_', ' ').split())].add(label)
                    squashed[squash(name)].add(label)
        self.exact = {k: next(iter(v)) for k, v in exact.items() if len(v) == 1}
        self.squashed = {k: next(iter(v)) for k, v in squashed.items() if len(v) == 1}

    def _lookup(self, name):
        n = ' '.join(name.lower().replace('_', ' ').split())
        for cand in (n, 'exposure to ' + n, n + ' exposure'):
            if cand in self.exact:
                return self.exact[cand]
        for cand in (n, 'exposure to ' + n, n + ' exposure'):
            t = self.squashed.get(squash(cand))
            if t:
                return t
        return None

    def term(self, value):
        """the term the value names without loss (exact, case / separators aside, or with 'treated
        with' / 'exposure' wording), or None (also for values with a dose or several agents)"""
        v = ' '.join(str(value).split()).lower()
        if not v or _NONE_VALUE.match(norm(v)):
            return None
        t = self._lookup(v)
        if t:
            return t
        if re.search(r'[+/&]|\band\b|\bwith\b|\bplus\b|;|,\s', v) and not re.match(r'^\d+,\d', v):
            # several agents ('cisplatin + paclitaxel') or a qualified treatment: no single term
            return None
        undosed = ' '.join(_DOSE.sub(' ', v).split()).strip(' -_,;:()')
        if undosed != v.strip(' -_,;:()'):
            # (the column cannot hold the dose)
            return None
        stripped = _AFFIX.sub('', _AFFIX.sub('', undosed)).strip(' -_,;:()')
        if stripped and len(stripped) > 2 and stripped != undosed:
            return self._lookup(stripped)
        return None


# ---------------------------------------------------------------- the pass


def apply(df, allowed_terms, exposure_table):
    """Returns (df with InternalStandardsUsed updated, per-factor table, stats). allowed_terms:
    allowed_terms.json dict; exposure_table: Exposure_ontology.csv (prepare_ontologies.py). The
    per-factor table has the columns sample_groups.remove_held_factors needs."""
    df = df.copy()
    is_index = InternalStandardIndex(allowed_terms['InternalStandardsUsed']['allowed_values'])
    ecto = EctoIndex(exposure_table, allowed_terms['Exposure']['allowed_values'])
    ds = df['ATTRIBUTE_DatasetAccession'].astype(str)
    stats = collections.Counter()

    # one row per internal-standard / exposure factor part of every row's text
    recs, cache = [], {}
    for i, text, m in zip(df.index, df['SampleGroupFreeText'], masst_rows(df)):
        for key, value, itype, part in parts_of(text, m):
            ck = (key, value, itype is not None)
            if ck not in cache:
                if itype is None and is_internal_standard_key(key):
                    cache[ck] = ('istd', is_index.term(value, key), False)
                elif itype is not None or is_exposure_key(key):
                    cache[ck] = ('exposure', ecto.term(value), is_none_value(value))
                else:
                    cache[ck] = None
            if cache[ck]:
                recs.append((i, key, value, part) + cache[ck])
    fac = pd.DataFrame(recs, columns=['row', 'key', 'value', 'part', 'kind', 'target', 'none_value'])
    fac['dataset'] = ds.reindex(fac['row']).to_numpy()
    fac['real'] = df['SampleType'].map(is_real_sample).reindex(fac['row']).to_numpy().astype(bool)

    # InternalStandardsUsed: on every row carrying the factor (blanks / QCs too); two
    # internal-standard factors on one row: the allowed combination of their names
    old_is = df['InternalStandardsUsed'].astype(str)
    ist = fac[(fac['kind'] == 'istd') & fac['target'].notna()]
    new_is = pd.Series(None, index=df.index, dtype=object)
    for i, terms in ist.groupby('row', sort=True)['target'].agg(lambda x: sorted(set(x))).items():
        new_is[i] = terms[0] if len(terms) == 1 else \
            is_index.combos.get(frozenset(n for t in terms for n in is_index.names(t)))
    fill = (old_is == MISSING) & new_is.notna()
    df.loc[fill, 'InternalStandardsUsed'] = new_is[fill]
    stats['InternalStandardsUsed_filled'] = int(fill.sum())
    stats['InternalStandardsUsed_conflicts_(column_kept)'] = sum(
        o != MISSING and isinstance(n, str) and not is_index.names(n) <= is_index.names(o) for n, o in zip(new_is, old_is))
    unknown = fac.loc[(fac['kind'] == 'istd') & fac['target'].isna(), 'value']
    stats['internal_standard_values_not_in_vocabulary'] = int(unknown.nunique())

    # held: the row's column holds the factor's value
    row_is = df['InternalStandardsUsed'].astype(str).reindex(fac['row']).to_numpy()
    row_e = df['Exposure'].astype(str).reindex(fac['row']).to_numpy()
    captured = []
    for kind, target, ri, re_ in zip(fac['kind'], fac['target'], row_is, row_e):
        if kind == 'istd':
            captured.append(target is not None and ri != MISSING and is_index.names(target) <= is_index.names(ri))
        else:
            captured.append(target is not None and re_ == target)
    fac['captured'] = captured
    # an explicit control / none value only where the row has an Exposure (else the factor stays in the text)
    fac['held'] = fac['captured'] | (fac['none_value'].astype(bool) & (row_e != MISSING))
    return df, fac, stats
