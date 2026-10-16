"""Reference materials and in-house pools named in sample-level metadata, for every source (GNPS
metadata files, MASST tables, MetaboLights, Metabolomics Workbench, NORMAN: through
complete_and_fill_REDU_table) and for the final table (rows kept from earlier runs: final_table_rules.py).

A row is a reference material when a sample-level field of its own names one (rules A-D below). It
then gets SampleType 'reference material', a SampleTypeSub1 as specific as the material is known, and
the organism and body part of the material only (never those of the study's samples):

    material named by its number (table)      'reference material_NIST SRM 1950'   organism, body part of the table
    NIST material, number unknown              'reference material_NIST'             human + the matrix (plasma, serum,
                                                                                     urine, milk, feces) if one is named
                                                                                     (no human where another animal is)
    other reference materials                  'missing value'                       human only if the text or the row's
                                                                                     own species field says so;
                                                                                     body part: the matrix; food: none
    purchased pools                            'missing value'                       human only if the text or the row's
                                                                                     own species field says so;
                                                                                     body part: the matrix

The matrix is a matrix word of the value that names the material, else one the row's own sample-level
fields give as their whole value ('Sample type: plasma', MetaboLights 'Organism part: blood plasma'):
never the body part the row has from elsewhere (Workbench study-level sample type).
'key: value' items (Workbench factors, GNPS group texts) are read by their value ('NIST: No' names no
material); 'no', 'false', '0' and the like name none.

Text is lower-cased, letters and digits are split apart ('NIST1950' -> 'nist 1950') and other
characters become spaces; words are matched whole ('administered', 'mechanistic', 'genistein' hold no
'nist').

Rules (sample-level fields only: file base name, sample / source / extract / assay names, characteristics,
factors, subject-sample-factor items, sheet columns other than SampleType and method columns, the MASST
intervention type; never protocol fields, where 'NIST 14/17/20' are spectral libraries and 'SRM' is
selected reaction monitoring):
  A numbered: 'SRM 1950', 'NIST SRM 1950', 'standard reference material 1950', 'RM 8231', 'RGTM 10162',
    'NIST 1950', 'NIST plasma 1950'; the number must be in the material table (reference_materials.tsv;
    bare numbers are never read: '1950' is mostly a date or a run number)
  B NIST word: 'NIST', 'NISTplasma', 'NIST pool', 'NIST QC', 'NIST reference material'
  C bare 'SRM': only where the study names exactly one material
  D phrases: 'standard / certified / candidate reference material', MASST 'human / food reference material'
  For B, C, D the study may name the material: if its other rows or its study-level text (title,
  description, protocols, mwTab sections) name exactly one material of the table that fits the row's
  matrix (and 'human' / 'food'), the row gets that material; a table number standing alone in the file
  name ('L_8256_NEG') picks one of several.
Not reference materials:
  - purchased pools (UTAK, BioIVT and its former name BioreclamationIVT / 'Biorec', 'pooled commercial
    plasma', 'pool serum from Sigma'; not 'commercial serum-free medium' or 'blood agar'): SampleType
    'reference material' without a NIST label (SampleTypeSub1 'missing value')
  - in-house pools ('long term reference', 'study reference material / sample / pool', 'HHEAR / CHEAR
    reference (material, plasma, ...)'): pool_QC (rows the source labels as samples become pool_QC; their
    organism and body part stay, as for every pooled QC); not the abbreviation 'LTR' (also 'long-term
    responders', letters of HHEAR / CHEAR sample codes such as 'H-9LTR0-PE-00', litres)
Guards:
  - an amount factor ('SRM 1950: 500', 'NIST 8231 = 0', 'SRM 1950 (uL): 0.5', 'Volume SRM1950: 10 uL',
    'SRM 1950: 50%') marks a mixture design: the row is left alone
  - a blank ('blank_analysis', 'blank_extraction', 'blank_culturemedia') keeps its type; so does a row whose
    reference-material value also holds a blank word ('NIST blank') or whose file is named as a blank
    (blank_qc_names); a personal care product keeps its type
  - a vendor name (purchased pool) beats a bare 'SRM' or a reference-material phrase of the row
  - a material named only in the directory of the file is ignored
  - CRM / ERM / BCR tokens are not read (CRM197, BCR-ABL); only 'certified reference material'
  - a row already labelled with the same or a more specific material keeps its label (a second pass
    changes nothing)

The material table is pinned in bin/reference_materials.tsv; its SampleTypeSub1 values are in
data/allowed_terms.json.
"""
import collections
import csv
import os
import re

import pandas as pd

from blank_qc_names import blank_qc_type_from_filename

MISSING = 'missing value'
HUMAN = '9606|Homo sapiens'
REFERENCE_MATERIAL = 'reference material'
SUB1_PREFIX = 'reference material_'
SUB1_NIST = SUB1_PREFIX + 'NIST'
_TABLE = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'reference_materials.tsv')

# an extra evidence column a converter may add to its frame: 'field\x1fvalue' items joined by '\x1e'
EVIDENCE_COLUMN = '_refmat_evidence'
# set on the rows this module labels as reference material (complete_and_fill_REDU_table then keeps
# their organism and body part, which come from the material)
PLACED_COLUMN = '_refmat_placed'
FIELD_SEP, ITEM_SEP = '\x1f', '\x1e'
# the converters' internal SampleGroupFreeText form (MWB_to_REDU _GROUP_SID, _GROUP_KV, _GROUP_PART)
_GROUP_SID, _GROUP_KV, _GROUP_PART = '\x1d', '\x1e', '\x1f'


def _read_table(path=_TABLE):
    out = {}
    with open(path, encoding='utf-8', newline='') as fh:
        rows = [r for r in csv.reader((l for l in fh if not l.startswith('#')), delimiter='\t')]
    header = rows[0]
    for r in rows[1:]:
        if not r or not r[0].strip():
            continue
        rec = dict(zip(header, r + [''] * (len(header) - len(r))))
        out[rec['number'].strip()] = {k: v.strip() for k, v in rec.items()}
    return out


MATERIALS = _read_table()


def material_sub1(number):
    return SUB1_PREFIX + MATERIALS[number]['label']


# every SampleTypeSub1 value this module writes (allowed_terms.json must hold them)
SUB1_VALUES = sorted({SUB1_NIST} | {material_sub1(n) for n in MATERIALS})
_SUB1_TO_NUMBER = {material_sub1(n): n for n in MATERIALS}
BLANK_TYPES = ('blank_analysis', 'blank_extraction', 'blank_culturemedia')
NOT_SAMPLE_TYPES = ('blank_analysis', 'blank_extraction', 'blank_culturemedia', 'blank_QC', 'pool_QC',
                    REFERENCE_MATERIAL, 'personal care product', 'molecular_standards')


def normalize(text):
    t = str(text).lower()
    t = re.sub(r'(?<=[a-z])(?=\d)|(?<=\d)(?=[a-z])', ' ', t)
    return re.sub(r'[^a-z0-9]+', ' ', t).strip()


# fields that describe how the sample was measured, not what it is
PROTOCOL_FIELD = re.compile(
    r'(?i)parameter value|instrument|spectromet|chromatograph|column|ionization|polarity|scan|extraction ?method|'
    r'collection ?method|internal ?standard|protocol|software|librar|database|ms_method|lc_method|scan_range|method|'
    r'massspectrometer|description|title|summary|keywords?|comment\[(?!sample)|sampleextraction|samplecollection|'
    r'massive|file_usi|filepath|^usi$|term (?:source|accession)|^unit')
# fields that hold ReDU's own sample type vocabulary (GNPS 'reference material' is any standard or product)
TYPE_FIELD = re.compile(r'(?i)^(?:attribute_)?(?:curated_)?sample_?type(?:_?sub1)?$|^study_sampletypeasrecorded$')

_NUM_RULES = [
    re.compile(r'\b(?:nist )?(?:standard reference material|srm|s r m)(?: material)? ?(\d{3,5})(?: [a-d])?\b'),
    re.compile(r'\b(?:nist )?(?:reference material|rm|rgtm) ?(\d{4,5})\b'),
    re.compile(r'\bnist(?: plasma| serum| urine| pool| milk)? ?(\d{3,5})(?: [a-d])?\b'),
]
_NIST = re.compile(r'\bnist(?:plasma|serum|urine|pool|qc|ref|milk)?\b')
# the NIST mass spectral library, not a material ('NIST14', 'NIST 2020', 'NIST/EPA/NIH', 'NIST MS library',
# 'NIST MS Search'): removed before the NIST rules (normalized text)
# (an edition number only with a library word: 'NIST_20.mzML' is the 20th NIST injection)
_NIST_LIBRARY = re.compile(r'\bnist (?:0?5|08|11|14|17|20|23|2005|2008|2011|2014|2017|2020|2023) (?:ms |msms )?'
                           r'(?:librar\w*|lib|search|database|match\w*)\b|'
                           r'\bnist (?:ms |msms |mass spectral |spectral )?(?:librar\w*|lib|search|database|epa|nih)\b|'
                           r'\b(?:mainlib|replib)\b')
_SRM_TOKEN = re.compile(r'\bsrm\b')
_PHRASE = re.compile(r'\b(?:(standard|certified|candidate) reference material|(human|food) reference material)s?\b')
# (not the abbreviation 'LTR': also 'long-term responders' (ST003966), letters of HHEAR / CHEAR sample codes
# ('H-9LTR0-PE-00', 'C-1LTR8-S-00') and litres)
_INHOUSE = re.compile(r'\b(?:long term reference|study reference (?:material|sample|pool)|(?:hhear|chear) reference)\b')
_PURCHASED = re.compile(r'\b(?:utak|bioivt|biorec|bioreclamation(?: ?ivt)?)\b|\bcommercial(?:ly)?(?: available)?(?: pooled?)?(?: human)? '
                        r'(?:plasma|serum|urine|pool|pooled|blood)\b(?! free| agar| medium| media| broth)|'
                        r'\bpooled? commercial (?:human )?(?:plasma|serum|urine)\b|'
                        r'\bpool(?:ed)? (?:human )?(?:plasma|serum|urine) from (?:sigma|commercial)\b|'
                        r'\b(?:plasma|serum|urine) pool from (?:sigma|commercial)\b')
# amount factors of a mixture design: 'SRM 1950: 500', 'NIST 8231 = 0', 'RM8231:-', also with a unit ('SRM 1950 (uL): 0.5',
# 'Volume SRM1950: 10 uL', 'SRM 1950: 50%'): a key naming the material, a number (or none) as its value
_AMOUNT = re.compile(r'(?i)(?:^|(?<=[;|,]))[^:=;|,]*?(?<![a-z])(?:nist|srm|rm)[ _-]*(\d{3,5})[a-d]?(?![0-9])[^:=;|,]*[:=]\s*'
                     r'(?:-|n/?a|none|[-+]?\d+(?:\.\d+)?)\s*(?:%|[a-z\u00b5]{1,3}(?:/[a-z\u00b5]{1,3})?)?\s*(?=$|[;|,])')
_BLANK_WORD = re.compile(r'(?<![a-z])(?:blank|blk|blnk)(?![a-z])')
_MATRIX = re.compile(r'(?:\b|(?<=nist))(plasma|serum|urine|milk|feces|fecal|faecal|faeces|stool)\b')
_MATRIX_PART = {'plasma': 'blood plasma', 'serum': 'blood serum', 'urine': 'urine', 'milk': 'milk', 'feces': 'feces',
                'fecal': 'feces', 'faecal': 'feces', 'faeces': 'feces', 'stool': 'feces'}
_PART_BIOFLUIDS = {'blood plasma', 'blood serum', 'urine', 'milk', 'feces'}
_HUMAN_WORD = re.compile(r'\bhuman\b')
# animals other than humans: no human organism for an unnumbered material of their matrix ('NIST mouse plasma',
# 'certified reference material bovine serum')
_NONHUMAN_WORD = re.compile(r'\b(?:bovine|cow|cows|calf|fetal bovine|fbs|mouse|mice|murine|rat|rats|porcine|pig|pigs|swine|'
                            r'canine|dog|dogs|equine|horse|horses|ovine|sheep|caprine|goat|goats|chicken|avian|fish|salmon|'
                            r'shrimp|monkey|monkeys|macaque|primate|baboon|rabbit|rabbits|hamster|guinea pig)\b')
# a value that only names a matrix ('plasma', 'blood serum', 'human urine', 'EDTA plasma'): the row's own matrix
_MATRIX_ONLY = re.compile(r'^(?:(?:human|pooled?|blood|whole|edta|heparin|citrate|frozen|fresh) )*'
                          r'(plasma|serum|urine|milk|feces|fecal|faecal|faeces|stool)(?: samples?)?$')
# a species / organism field of the row itself that says human ('Species: Human', MetaboLights 'Organism: Homo
# sapiens', GNPS NCBITaxonomy '9606|Homo sapiens'): the row's own organism
_SPECIES_FIELD = re.compile(r'(?i)species|organism(?! ?part)|ncbitaxonomy|taxon')
_HUMAN_VALUE = {'human', 'humans', 'homo sapiens', '9606 homo sapiens', 'h sapiens'}
# values that name nothing ('NIST: No', 'Spiked NIST: 0')
_NO_VALUE = {'no', 'n', 'false', 'none', '0', 'na', 'n a', 'not applicable', 'not available', 'absent', 'unspiked'}
_FOOD_WORD = re.compile(r'\bfood\b')
# numbered materials named in study-level text: 'NIST SRM 1950', 'SRM1950', 'NIST serum (909C)',
# 'NIST urine (SRM 3672)', 'RM 8231', 'RGTM 10162'
_TEXT_NUMBER = re.compile(r'(?i)(?<![a-z])(?:nist[ _-]*(?:[a-z]+[ _-]+){0,2}\(?[ _-]*(?:srm|rm|standard reference material|'
                          r'reference material)?|srm|rm|rgtm|(?:standard )?reference material(?: \((?:s?rm)\))?)'
                          r'[ _#-]*(\d{3,5})(?:[ _-]?[a-d])?(?![0-9])')
# a cheap test before the rules (raw text)
_CANDIDATE = re.compile(r'(?i)nist|srm|reference[ _-]?material|rgtm|(?<![a-z])rm[ _-]?\d{4}|utak|bioivt|biorec|sigma|commercial|'
                        r'long[ _-]?term|study[ _-]?reference|hhear|chear')
CANDIDATE_RE = _CANDIDATE


def is_amount_factor(value):
    """whether a value is an amount factor of a material of the table ('SRM 1950: 500', 'NIST 8231 = 0'): a
    mixture design"""
    return any(m.group(1) in MATERIALS for m in _AMOUNT.finditer(str(value)))


def study_text_materials(*texts):
    """catalog numbers of the table named in study-level text (sorted list)"""
    found = set()
    for text in texts:
        if not text:
            continue
        for m in _TEXT_NUMBER.finditer(str(text)):
            if m.group(1) in MATERIALS:
                found.add(m.group(1))
    return sorted(found)


def _value_part(value):
    """the value of a 'key: value' / 'key=value' item (for the blank-word guard: 'Blank:NIST' is a factor named Blank)"""
    return re.sub(r'^[^:=|]{1,40}[:=]', '', str(value)).lower()


def _keyed(name, value):
    """whether a field value is a 'key: value' item (Workbench factor items, GNPS group texts): its key is the
    field name or the field is a group text part"""
    name = str(name or '')
    return name.startswith('group:') or (bool(name) and str(value).lower().startswith(name.lower() + ':'))


def is_matrix_value(value):
    """whether a value (of a 'key: value' item: its value) only names a matrix ('plasma', 'blood serum')"""
    return bool(_MATRIX_ONLY.match(normalize(_value_part(value))))


def is_human_species_item(name, value):
    """whether a field is a species / organism field of the row that says human"""
    return bool(_SPECIES_FIELD.search(str(name or ''))) and normalize(_value_part(value)) in _HUMAN_VALUE


class Evidence:
    """what the sample-level fields of one row say"""
    __slots__ = ('numbers', 'amount', 'nist', 'srm', 'phrases', 'inhouse', 'purchased', 'blank_word', 'file_blank', 'matrix',
                 'human', 'food', 'dir_only', 'name_text', 'own_matrix', 'nonhuman', 'own_human')

    def __init__(self):
        self.numbers, self.phrases, self.matrix, self.own_matrix = set(), set(), set(), set()
        self.amount = self.nist = self.srm = self.inhouse = self.purchased = False
        self.blank_word = self.file_blank = self.human = self.food = self.dir_only = self.nonhuman = self.own_human = False
        self.name_text = ''

    def any_rm(self):
        return bool(self.numbers or self.nist or self.srm or self.phrases)


def row_evidence(fields, directory=''):
    """Evidence of one row from (field name, value) pairs of its sample-level fields (protocol and
    sample-type fields are skipped) and the directory of its file (only to tell 'named only there')."""
    ev = Evidence()
    for name, value in fields:
        if value is None:
            continue
        value = str(value)
        if not value or value.strip().lower() in ('nan', MISSING, 'not applicable'):
            continue
        if name and name != 'file' and (PROTOCOL_FIELD.search(name) or TYPE_FIELD.search(name)):
            continue
        if name == 'file':
            # a table number standing alone in the file name can pick one of the study's materials
            ev.name_text += ' ' + value
            # a file named as a blank ('SRM1950_blank', 'IS blank') is no material
            if (blank_qc_type_from_filename(value) or '').startswith('blank'):
                ev.file_blank = True
        # a 'key: value' item names a material by its value ('NIST: No' names none); the amount test reads the
        # whole item ('SRM 1950: 0' is an amount factor of a mixture design)
        keyed = name != 'file' and _keyed(name, value)
        text = _value_part(value) if keyed else value
        hit = False
        if _CANDIDATE.search(value):
            amounts = {m.group(1) for m in _AMOUNT.finditer(value)}
            if amounts & set(MATERIALS):
                ev.amount = True
                hit = True
        if keyed and normalize(text) in _NO_VALUE:
            continue
        if not _CANDIDATE.search(value):
            # a value that only names a matrix: the row's own matrix
            if name != 'file' and _MATRIX_ONLY.match(normalize(text)):
                ev.own_matrix.update(_MATRIX_PART[w] for w in _MATRIX.findall(normalize(text)))
            if name != 'file' and is_human_species_item(name, value):
                ev.own_human = True
            continue
        if not _CANDIDATE.search(text):
            continue
        t = _NIST_LIBRARY.sub(' ', normalize(text))
        for rx in _NUM_RULES:
            for m in rx.finditer(t):
                if m.group(1) in MATERIALS:
                    ev.numbers.add(m.group(1))
                    hit = True
        if _NIST.search(t):
            ev.nist = hit = True
        if _SRM_TOKEN.search(t):
            ev.srm = hit = True
        for m in _PHRASE.finditer(t):
            ev.phrases.add(m.group(1) or m.group(2))
            hit = True
        if _INHOUSE.search(t):
            ev.inhouse = True
        if _PURCHASED.search(t):
            ev.purchased = True
        if hit or ev.purchased:
            if _BLANK_WORD.search(normalize(_value_part(value))) and hit:
                ev.blank_word = True
            ev.matrix.update(_MATRIX_PART[w] for w in _MATRIX.findall(t))
            if _HUMAN_WORD.search(t):
                ev.human = True
            if _NONHUMAN_WORD.search(t):
                ev.nonhuman = True
            if _FOOD_WORD.search(t):
                ev.food = True
    if not ev.matrix:
        ev.matrix = set(ev.own_matrix)
    if not ev.any_rm() and directory and _CANDIDATE.search(directory):
        d = normalize(directory)
        if _NIST.search(d) or any(rx.search(d) and rx.search(d).group(1) in MATERIALS for rx in _NUM_RULES):
            ev.dir_only = True
    return ev


class Decision:
    __slots__ = ('status', 'sample_type', 'sub1', 'organism', 'body_part', 'material', 'rank')

    def __init__(self, status, sample_type=None, sub1=None, organism=None, body_part=None, material='', rank=0):
        self.status, self.sample_type, self.sub1 = status, sample_type, sub1
        self.organism, self.body_part, self.material, self.rank = organism, body_part, material, rank

    def __repr__(self):
        return (f'Decision({self.status!r}, {self.sample_type!r}, {self.sub1!r}, {self.organism!r}, '
                f'{self.body_part!r}, {self.material!r})')


def _parts(ev, own_part=None):
    """the matrix the row names (row_evidence: the value naming the material, else the row's own matrix items);
    own_part (the body part the row has from elsewhere, e.g. a Workbench study-level sample type) is not used"""
    return set(ev.matrix)


def _fits(number, parts, ev):
    m = MATERIALS[number]
    if parts and m['UBERONBodyPartName'] not in parts:
        return False
    if ev.human and m['NCBITaxonomy'] != HUMAN:
        return False
    if ev.food and not ev.human and m['NCBITaxonomy']:
        return False
    return True


def classify(ev, study_materials=(), own_part=''):
    """Decision for one row (None: nothing to do). study_materials: catalog numbers the study names (its
    other rows and its study-level text); own_part: not used (the matrix comes from the row's own fields)."""
    if ev.amount:
        return Decision('mixture_design')
    if ev.blank_word or (ev.file_blank and (ev.any_rm() or ev.purchased or ev.inhouse)):
        # 'NIST blank', 'SRM1950_Blk', a file named as a blank: a blank of the material's run, not the material
        return Decision('blank_word')
    if ev.purchased and not ev.numbers and not ev.nist:
        # a purchased pool named by its vendor ('BioIVT', 'UTAK', 'Pooled Commercial Plasma') is that product,
        # also where the row is in a group named 'SRM' or 'reference material'
        return _purchased(ev, own_part)
    if ev.numbers:
        if len(ev.numbers) == 1:
            n = next(iter(ev.numbers))
            m = MATERIALS[n]
            return Decision('numbered', REFERENCE_MATERIAL, material_sub1(n), m['NCBITaxonomy'] or MISSING,
                            m['UBERONBodyPartName'] or MISSING, n, 3)
        ms = [MATERIALS[n] for n in sorted(ev.numbers)]
        orgs, parts = {m['NCBITaxonomy'] for m in ms}, {m['UBERONBodyPartName'] for m in ms}
        return Decision('several', REFERENCE_MATERIAL, SUB1_NIST, (orgs.pop() or MISSING) if len(orgs) == 1 else MISSING,
                        (parts.pop() or MISSING) if len(parts) == 1 else MISSING, '+'.join(sorted(ev.numbers)), 2)
    if ev.nist or ev.srm or ev.phrases:
        parts = _parts(ev, own_part)
        cands = {n for n in study_materials if n in MATERIALS and _fits(n, parts, ev)}
        in_name = {n for n in cands if re.search(r'(?<![0-9])' + n + r'(?![0-9])', ev.name_text)}
        pick = None
        if len(cands) > 1 and len(in_name) == 1:
            pick = next(iter(in_name))
        elif len(cands) == 1:
            pick = next(iter(cands))
        if pick:
            m = MATERIALS[pick]
            return Decision('resolved_by_study', REFERENCE_MATERIAL, material_sub1(pick), m['NCBITaxonomy'] or MISSING,
                            m['UBERONBodyPartName'] or MISSING, pick, 3)
        if not ev.nist and not ev.phrases:
            return Decision('srm_unresolved')
        part = next(iter(parts)) if len(parts) == 1 else MISSING
        if ev.nist:
            # every NIST biofluid material is human (not where the text names another animal)
            human = part != MISSING and not (ev.nonhuman and not ev.human)
            return Decision('nist_unspecified', REFERENCE_MATERIAL, SUB1_NIST, HUMAN if human else MISSING,
                            part, 'NIST', 2)
        food = ('food' in ev.phrases or ev.food) and 'human' not in ev.phrases and not ev.human
        if food:
            return Decision('phrase_unspecified', REFERENCE_MATERIAL, MISSING, MISSING, MISSING, '', 1)
        # an unnumbered reference material is human only where its text says so ('human reference material',
        # 'standard reference material human plasma'): reference materials of plasma or serum are also made from
        # other animals ('certified reference material bovine serum')
        human = 'human' in ev.phrases or ev.human or (ev.own_human and not ev.nonhuman)
        return Decision('phrase_unspecified', REFERENCE_MATERIAL, MISSING, HUMAN if human else MISSING, part, '', 1)
    if ev.dir_only:
        return Decision('dir_only')
    if ev.inhouse:
        return Decision('inhouse_pool', 'pool_QC', 'pool_QC')
    if ev.purchased:
        return _purchased(ev, own_part)
    return None


def _purchased(ev, own_part):
    """a purchased pool: no material of the table, so no organism unless its text says human ('commercial human
    plasma'; a 'pool serum from Sigma' of a mouse study may be mouse serum); the body part of its matrix"""
    parts = _parts(ev, own_part)
    part = next(iter(parts)) if len(parts) == 1 else MISSING
    human = ev.human or (ev.own_human and not ev.nonhuman)
    return Decision('purchased_pool', REFERENCE_MATERIAL, MISSING, HUMAN if human else MISSING, part, '', 1)


def current_rank(sample_type, sub1):
    """how specific the reference-material label a row already has is (see Decision.rank)"""
    if sample_type != REFERENCE_MATERIAL:
        return -1
    if sub1 in _SUB1_TO_NUMBER:
        return 3
    if sub1 == SUB1_NIST:
        return 2
    # another SampleTypeSub1 the source gave (GNPS 'reference material_chemicalstandard') is as specific as an
    # unnumbered reference material; so is a reference material without one: an unnumbered or purchased material
    # read again (the final table's pass, which reads fewer fields than the converter) changes nothing on it
    return 1


def applies(decision, sample_type, sub1):
    """whether a row of this SampleType / SampleTypeSub1 takes the decision"""
    if decision is None or decision.sample_type is None:
        return False
    st = '' if sample_type is None else str(sample_type)
    if st in BLANK_TYPES or st == 'personal care product':
        return False
    if decision.sample_type == 'pool_QC':
        # in-house pools: only rows labelled as samples (or not at all) become pooled QCs
        return st not in NOT_SAMPLE_TYPES
    return decision.rank > current_rank(st, '' if sub1 is None else str(sub1))


def _split_evidence(text):
    if not isinstance(text, str) or not text:
        return []
    out = []
    for item in text.split(ITEM_SEP):
        if FIELD_SEP in item:
            k, v = item.split(FIELD_SEP, 1)
        else:
            k, v = '', item
        out.append((k, v))
    return out


def encode_evidence(pairs):
    """(field, value) pairs -> one EVIDENCE_COLUMN cell (values with a candidate word, and values that only name a
    matrix: the row's own matrix)"""
    out = []
    for k, v in pairs:
        if v is None:
            continue
        v = str(v).replace(ITEM_SEP, ' ').replace(FIELD_SEP, ' ')
        if not v.strip() or not (_CANDIDATE.search(v) or is_matrix_value(v) or is_human_species_item(k, v)):
            continue
        out.append(f'{str(k).replace(ITEM_SEP, " ").replace(FIELD_SEP, " ")}{FIELD_SEP}{v}')
    return ITEM_SEP.join(dict.fromkeys(out))


def _frame_fields(df, i, columns, extra):
    fields = []
    fn = df.at[i, 'filename'] if 'filename' in columns else ''
    fn = '' if fn is None else str(fn)
    fields.append(('file', os.path.basename(fn.rstrip('/'))))
    if 'SubjectIdentifierAsRecorded' in columns:
        fields.append(('SubjectIdentifierAsRecorded', df.at[i, 'SubjectIdentifierAsRecorded']))
    if 'SampleGroupFreeText' in columns:
        text = str(df.at[i, 'SampleGroupFreeText'])
        if _GROUP_SID in text or _GROUP_KV in text:
            # the converters' internal form ('<sample>\x1d<key>\x1e<value>' parts joined by '\x1f'): 'key: value'
            for part in text.split(_GROUP_PART):
                key, _, value = part.rsplit(_GROUP_SID, 1)[-1].partition(_GROUP_KV)
                fields.append(('group:' + key.strip()[:40], f'{key.strip()}: {value.strip()}'))
        else:
            for part in text.split('; '):
                fields.append(('group:' + part.split(':', 1)[0].strip()[:40], part))
    if extra is not None:
        fields.extend(_split_evidence(extra.get(i)))
    return fields, os.path.dirname(fn.rstrip('/'))


def decide_frame(df, dataset_column='MassiveID', evidence=None, study_text=None):
    """Decisions for the rows of a frame in ReDU columns (filename, SubjectIdentifierAsRecorded,
    SampleGroupFreeText, SampleType, SampleTypeSub1, UBERONBodyPartName) plus optional extra evidence
    (a Series of EVIDENCE_COLUMN cells, by index). study_text: {dataset: [catalog numbers]} or a list of
    catalog numbers for all rows, from study-level text. Returns {index: Decision} for rows with evidence
    (decisions that take no effect included: mixture designs, unresolved 'SRM', directory-only)."""
    if len(df) == 0:
        return {}
    columns = set(df.columns)
    text = pd.Series('', index=df.index, dtype=object)
    for c in ('filename', 'SubjectIdentifierAsRecorded', 'SampleGroupFreeText'):
        if c in columns:
            text = text + ' ' + df[c].fillna('').astype(str)
    if evidence is not None:
        evidence = evidence.reindex(df.index)
        text = text + ' ' + evidence.fillna('').astype(str)
    cand = df.index[text.str.contains(_CANDIDATE, regex=True)]
    if len(cand) == 0:
        return {}
    extra = evidence.to_dict() if evidence is not None else None
    evs = {}
    for i in cand:
        fields, directory = _frame_fields(df, i, columns, extra)
        ev = row_evidence(fields, directory)
        if ev.any_rm() or ev.inhouse or ev.purchased or ev.dir_only or ev.amount:
            evs[i] = ev
    if not evs:
        return {}
    dataset = df[dataset_column].fillna('').astype(str) if dataset_column in columns else pd.Series('', index=df.index)
    named = collections.defaultdict(set)
    for i, ev in evs.items():
        named[dataset.at[i]].update(ev.numbers)
    if isinstance(study_text, dict):
        for d, nums in study_text.items():
            named[str(d)].update(n for n in nums if n in MATERIALS)
    elif study_text:
        for d in set(dataset[list(evs)]):
            named[d].update(n for n in study_text if n in MATERIALS)
    own_part = df['UBERONBodyPartName'].fillna('').astype(str) if 'UBERONBodyPartName' in columns else None
    out = {}
    for i, ev in evs.items():
        d = classify(ev, sorted(named[dataset.at[i]]), own_part.at[i] if own_part is not None else '')
        if d is not None:
            out[i] = d
    return out


def apply_decisions(df, decisions, stats=None, label=''):
    """Writes the decisions that apply (applies()) into SampleType, SampleTypeSub1, NCBITaxonomy and
    UBERONBodyPartName (reference materials; in-house pools only change their type) and sets PLACED_COLUMN
    on the reference materials. Returns the frame and the indexes changed."""
    stats = stats if stats is not None else collections.Counter()
    if not decisions:
        return df, []
    for c in ('SampleType', 'SampleTypeSub1', 'NCBITaxonomy', 'UBERONBodyPartName'):
        if c not in df.columns:
            df[c] = MISSING
    if PLACED_COLUMN not in df.columns:
        df[PLACED_COLUMN] = False
    changed = []
    for i, d in decisions.items():
        st, s1 = df.at[i, 'SampleType'], df.at[i, 'SampleTypeSub1']
        stats[f'{label}decision {d.status}'] += 1
        if d.sample_type is None:
            continue
        if not applies(d, st, s1):
            stats[f'{label}kept {d.status} (row is {st})'] += 1
            continue
        same = (st == d.sample_type and s1 == d.sub1 and (d.sample_type != REFERENCE_MATERIAL or (
            df.at[i, 'NCBITaxonomy'] == d.organism and df.at[i, 'UBERONBodyPartName'] == d.body_part)))
        if same:
            # (a second pass)
            stats[f'{label}unchanged {d.status}'] += 1
            if d.sample_type == REFERENCE_MATERIAL:
                df.at[i, PLACED_COLUMN] = True
            continue
        stats[f'{label}applied {d.status} from {st}'] += 1
        df.at[i, 'SampleType'] = d.sample_type
        df.at[i, 'SampleTypeSub1'] = d.sub1
        if d.sample_type == REFERENCE_MATERIAL:
            df.at[i, 'NCBITaxonomy'] = d.organism
            df.at[i, 'UBERONBodyPartName'] = d.body_part
            df.at[i, PLACED_COLUMN] = True
        changed.append(i)
    return df, changed


def annotate(df, dataset_column='MassiveID', evidence=None, study_text=None, stats=None, label=''):
    """decide_frame + apply_decisions. evidence: None, a Series, or the name of a column of df (e.g.
    EVIDENCE_COLUMN)."""
    if isinstance(evidence, str):
        evidence = df[evidence] if evidence in df.columns else None
    decisions = decide_frame(df, dataset_column=dataset_column, evidence=evidence, study_text=study_text)
    return apply_decisions(df, decisions, stats=stats, label=label)


def place_labelled(df, stats=None):
    """Rows the source already labels with a material of the table (SampleType 'reference material',
    SampleTypeSub1 'reference material_NIST SRM 1950', ...) that annotate did not place take the
    material's organism and body part too (PLACED_COLUMN set)."""
    if 'SampleType' not in df.columns or 'SampleTypeSub1' not in df.columns:
        return df
    if PLACED_COLUMN not in df.columns:
        df[PLACED_COLUMN] = False
    placed = df[PLACED_COLUMN].astype(str).str.lower().eq('true')
    rows = df.index[(df['SampleType'] == REFERENCE_MATERIAL) & df['SampleTypeSub1'].isin(list(_SUB1_TO_NUMBER)) & ~placed]
    for c in ('NCBITaxonomy', 'UBERONBodyPartName'):
        if c not in df.columns:
            df[c] = MISSING
    for i in rows:
        m = MATERIALS[_SUB1_TO_NUMBER[df.at[i, 'SampleTypeSub1']]]
        df.at[i, 'NCBITaxonomy'] = m['NCBITaxonomy'] or MISSING
        df.at[i, 'UBERONBodyPartName'] = m['UBERONBodyPartName'] or MISSING
        df.at[i, PLACED_COLUMN] = True
    if stats is not None and len(rows):
        stats['labelled by the source, organism / body part of the material'] += len(rows)
    return df


def placed_material(sub1):
    """the table entry of a 'reference material_<label>' SampleTypeSub1, or None"""
    n = _SUB1_TO_NUMBER.get(sub1)
    return MATERIALS[n] if n else None
