"""Diet and FastingStatus from the sample group texts (SampleGroupFreeText), for every row of the
final table (rows of this run and rows carried over from earlier runs alike; final_text_columns.py),
and the mapping microbiomeMASST diet interventions share (masst_tables.py).

Diet holds one value: a regimen term (the Diet terms of allowed_terms.json that are not foods,
ONS-aligned: 'high-fat diet', 'control diet', ...) or, for the food a subject ate or was fed as its
dietary treatment, a GFOP food term (the FoodOntologyTerm vocabulary, which updateAllowedTerms adds
to Diet's allowed values).
- A value names a regimen only through an explicit diet wording ('HFD', 'high fat diet', 'normal
  chow', 'Diet, Mediterranean', 'lithogenic diet', 'intermittent fasting'); bare codes ('CD', 'HF',
  'control', 'normal', 'TRF', 'inulin') only under a factor that is itself about the diet ('Diet',
  'DietGroup', microbiomeMASST 'Intervention: mice diet: ...').
- A food only where the factor is about feeding (a diet / feed / meal / supplement factor, or a
  treatment whose value says fed / feeding / diet / supplement) and the food, without dose and
  wrapper words, is a GFOP label ('yeast feeding' -> yeast, 'Supplement: Cranberry' -> cranberry);
  never from a factor that merely names a food ('Food:', 'Cheese maturation:', 'Milk:').
- One value per row: a regimen other than 'control diet' wins over a food named with it
  ('HFD + apple' -> high-fat diet; the factor stays in the free text); 'control diet' only for an
  unmodified control / chow diet (a control diet with something added is not the control group),
  so 'basal diet + <food>' gives the food; different regimens (or two foods) of one row: no Diet.
- Not mapped: diets that change over time (reverse, switch, washout, pre-/post-diet), time-point
  values, bespoke formulas without a regimen word, and samples whose time-point factor says they
  were taken before the intervention ('Visit: Baseline').
- A yes / no factor named after a diet ('HFD: 1') gives that diet for yes; no gives nothing.
- A food only on a sample with a body part (a 'milk' file without one is more likely the food).
FastingStatus: 'fasted' where a part of a value is only a fasting state ('Challenge: Fasting',
'WT - 16h fast', 'Norm Chow 6h Fast'), both states under a fasting factor ('Nutritional status: Fed
state, normal chow diet'); other values naming a state ('Hsc70 fed', 'WT - Fed', 'ALS: Fast',
'Starving') only where the dataset's same factor names both states ('fed' alone can be a feeding
treatment, as in the Workbench rules). A fasting diet group ('DietGroup: fasting', intermittent
fasting) is a regimen, not the state at sampling.

Guards (subject_guard): a value from the text needs a sample of an animal that eats: SampleType
'animal', NCBITaxonomy filled and in Metazoa (NCBIKingdom of the pinned lineage, humans included),
no FoodOntologyTerm and no food / beverage SampleType (the sample is not itself a food), no blank,
QC, standard, reference material, culture or cell line; human dietary patterns (omnivore,
vegetarian, vegan, pescatarian) only for Homo sapiens. Values already in the columns: a food needs
the same; a regimen or fasting state set by a converter is only cleared where the row is not an
animal sample or is a food / blank / QC / culture / cell line, or its taxon is not Metazoa (an
unresolved NCBITaxonomy on an animal sample keeps it).

The text itself is not changed here: apply() returns, per factor part, whether the columns hold
its value (held / captured), for sample_groups.remove_held_factors; an explicit control / none value
('Supplement: Placebo', 'Treatment: control') counts as held only where the row has a Diet.
"""
import collections
import re

import pandas as pd

from sample_groups import MISSING, is_real_sample, masst_rows, norm, parts_of

HUMAN_PATTERNS = {'omnivore', 'vegetarian', 'vegan', 'pescatarian'}
CONTROL = 'control diet'

# ---------------------------------------------------------------- text helpers


def _rx(pattern):
    return re.compile(r'(?<![a-z0-9])(?:' + pattern + r')(?![a-z0-9])')


# regimen wordings (on norm() text), most specific first. (term, wording, short): short wordings are
# codes or bare words that name a diet only under a diet factor
_FIBRES = (r'inulin|arabinoxylans?|oligofructose|fructo ?oligosaccharides?|galacto ?oligosaccharides?|fos|gos|'
           r'pectins?|resistant starch|beta glucans?|psyllium|guar gum')
_REGIMENS = [
    ('high-fat high-sugar diet', r'high fat high (?:sugar|sucrose|carbohydrate|carb)s?|high (?:sugar|sucrose) high fat|'
                                 r'high fat (?:and|&) (?:high )?(?:sugar|sucrose)|hfhsd?|hchf|hfhc', False),
    ('high-fat high-sugar diet', r'hf hs', True),
    ('Western diet', r'western(?: style)? diet|diet western|westernized diet', False),
    ('Western diet', r'western|wd', True),
    ('high-fat diet', r'(?:very )?high fat(?: diet)?|highfat|hfd|hf diet', False),
    ('high-fat diet', r'hf', True),
    ('low-fat diet', r'low fat(?: diet)?|lfd|lf diet', False),
    ('high-sugar diet', r'high (?:sugar|sucrose|fructose) diet', False),
    ('high-sugar diet', r'high (?:sugar|sucrose|fructose)|hsd', True),
    ('high-protein diet', r'high protein(?: diet)?', False),
    ('ketogenic diet', r'ketogenic(?: diet)?', False),
    ('ketogenic diet', r'keto|kd', True),
    ('Mediterranean diet', r'mediterranean(?! dash| mind)(?: diet)?', False),
    ('lithogenic diet', r'lithogenic diet', False),
    ('lithogenic diet', r'lithogenic', True),
    ('iron-deficient diet', r'(?:iron|fe) (?:deficien(?:t|cy)|deplet(?:ed|ion)) diet|low iron diet', False),
    ('iron-deficient diet', r'(?:iron|fe) (?:deficien(?:t|cy)|deplet(?:ed|ion))|low iron', True),
    ('fibre-supplemented diet', r'(?<!low )(?<!no )(?<!non )(?:high )?fib(?:er|re)s?(?: supplemented| enriched| rich)? ?diet|'
                                r'fib(?:er|re) supplement(?:ed|ation)?|(?<!low )high fib(?:er|re)|'
                                r'(?:' + _FIBRES + r')(?: supplemented| enriched)? diet', False),
    ('fibre-supplemented diet', _FIBRES, True),
    ('prebiotic diet', r'prebiotics? diet', False),
    # (bare 'prebiotic' also under a feeding factor: 'Milk formula intervention: prebiotic'; a
    # synbiotic is the prebiotic with a probiotic, which the term does not hold: _PROBIOTIC_TOO)
    ('prebiotic diet', r'prebiotics?|synbiotics?|symbiotics?', True),
    ('intermittent fasting', r'intermittent fast(?:ing)?', False),
    ('time-restricted feeding', r'time restricted (?:feeding|eating|diet)', False),
    ('time-restricted feeding', r'time restricted|trf', True),
    # ('food restricted' / 'restricted feeding' can be time-restricted feeding)
    ('calorie-restricted diet', r'calori(?:c|e)s? restrict(?:ed|ion)|restrict(?:ed|ion) calori(?:c|e)s?', False),
    ('calorie-restricted diet', r'cr|restricted', True),
    ('omnivore', r'omnivor(?:e|es|ous)', False),
    ('vegetarian', r'(?:lacto )?(?:ovo )?vegetarian', False),
    ('vegan', r'vegan', False),
    ('pescatarian', r'pesc[ae]tarian', False),
    (CONTROL, r'(?:normal|norm|regular|standard|std|control|ctrl|vivarium|viv|rodent|laboratory) chow(?: diet)?|chow(?: diet)?|'
              r'(?:normal|standard|std|regular|control|ctrl|basal|purified control|control purified|'
              r'standard rodent|standard laboratory) diet|ctrldiet|controldiet|ncd|nc diet|'
              r'ain 93(?:g|m)?(?: control)?(?: diet)?|ain 76a?', False),
    (CONTROL, r'control|ctrl|con|normal|standard|std|regular|basal|cd|nd|nc|sc', True),
]
_REGIMENS = [(t, _rx(p), short) for t, p, short in _REGIMENS]
# words that do not add information to a diet value
_FILLER = {'diet', 'diets', 'dietary', 'group', 'groups', 'fed', 'feeding', 'feed', 'mice', 'mouse', 'rat', 'rats',
           'the', 'a', 'an', 'of', 'on', 'state', 'status', 'regimen', 'intervention', 'type', 'ad', 'libitum',
           'lib', 'only', 'alone', 'animals', 'arm', 'subjects', 'was', 'were', 'is', 'are', 'given',
           'received', 'receiving', 'by', 'and', 'with'}
# a control value carries these too ('Normal diet group (Con)')
_CONTROL_FILLER = {'control', 'controls', 'ctrl', 'con', 'normal', 'standard', 'std', 'regular', 'basal', 'cd',
                   'ncd', 'nd', 'sd', 'sc', 'chow'}
# values that change the diet over time or name a time point: no diet from them
_TIME = _rx(r'pre|post|before|after|baseline|washout|wash out|reverse|reversal|reversed|switch|switched|switching|'
            r'off|withdrawal|withdrawn|recovery|day|days|week|weeks|wk|wks|month|months|year|years|'
            r'\d+ ?(?:d|w|wk|wks|mo|m)|on diet')
_NEGATION = _rx(r'no|non|not|without|free of')
# the diet of the parents ('Offspring of High-Fat-Diet'), not of the sampled subject
_PARENTS = _rx(r'offspring|maternal|mother|mothers|dam|dams|paternal|father|fathers|sire|parental|in utero|perinatal')
# a control diet changed itself (then it is not the control group's diet)
_MODIFIED = _rx(r'supplement\w*|enrich\w*|fortif\w*|containing|contains|deficien\w*|deficent|free|depleted|devoid|'
                r'added|addition|rich|low [a-z]+|high [a-z]+|with')
_PROBIOTIC_TOO = _rx(r'synbiotics?|symbiotics?')
# named diets without a Diet term: a value naming one of them and a regimen term lists several diets
# ('anti-inflammatory, western' is not the Western diet); under a diet factor (whose values are
# diets) as such, under other factors only as '<name> diet'
_OTHER_DIET_NAMES = (r'anti ?inflammatory|dash|mind|paleo(?:lithic)?|atkins|nordic|okinawan?|macrobiotic|'
                     r'flexitarian|low carb(?:ohydrate)?s?|high carb(?:ohydrate)?s?|low protein|'
                     r'protein restrict(?:ed|ion)|gluten free|(?:low )?fodmap|elimination')
_OTHER_DIET = _rx(_OTHER_DIET_NAMES)
_OTHER_DIET_WORDED = _rx(r'(?:' + _OTHER_DIET_NAMES + r') diets?')

# ---------------------------------------------------------------- factor keys

# factor names about the diet (bare codes allowed), about feeding (foods allowed), about fasting
_DIET_KEY = _rx(r'diet|diets|dietary|diet group|dietgroup|diet type|diet_type|chow|feeding regime|food group')
# (not 'food' / 'nutrition': 'Food: cheese' can name the analysed food, not what the subject ate)
_FEEDING_KEY = _rx(r'diet|diets|dietary|feed|feeding|meal|meals|supplement|supplements|supplementation|formula|'
                   r'intake|fed')
_FASTING_KEY = _rx(r'fast|fasting|fasted|nutritional status|nutritional state|feeding state|feeding status|'
                   r'prandial state|prandial status|fed state')
# factor names that are never about what the subject ate
_NEVER_KEY = _rx(r'food type|food source|food spoiler|medium|media|culture|additive|additives|substrate|'
                 r'carbon source|isotope|isotopomer|tracer|label|labelled|internal standard|standard|solvent|'
                 r'extraction|extract|matrix|sample type|origin|compound|population|cell|cells')
_FEEDING_WORD = _rx(r'diet|diets|fed|feed|feeding|supplement|supplemented|supplementation|intake|meal|chow')
# value words: the sample is a culture / medium, not a subject; or a labelled tracer
_NEVER_VALUE = _rx(r'medium|media|culture|cultured|broth|agar|in vitro|cells|13c|15n|deuterated|labell?ed')

# ---------------------------------------------------------------- fasting
_FASTED = _rx(r'(?:\d+(?:\.\d+)? ?(?:h|hr|hrs|hour|hours|d|day|days) )?(?:fasted|fasting)|'
              r'\d+(?:\.\d+)? ?(?:h|hr|hrs|hour|hours) fast|\d+(?:\.\d+)?(?:h|hr|hrs) fast|overnight fast\w*')
_FAST_BARE = _rx(r'fast')
# food withheld: fasted where the same factor also names fed samples ('Fed ad libitum' / 'Starving'),
# only as the whole value ('Zinc Starvation', 'N-starvation' are nutrient limitations)
_STARVED = _rx(r'starv(?:ed|ing|ation)|unfed|food deprived|food deprivation')
_NONFASTED = _rx(r'non ?fasted|nonfasted|unfasted|not fasted|fed state|fed ad libitum|ad libitum fed|fed|re ?fed|'
                 r'postprandial|post prandial')
_FAST_NOT = _rx(r'intermittent|time restricted|fasting mimicking|alternate day|fast growing|fast progress\w*|'
                r'fast transit|carnitine|fast food|breakfast|fasting glucose|fasting insulin|fasting blood')
_FAST_DURATION = re.compile(r'\d')

# time-point factors whose value is the sampling before the intervention
_TIME_KEY = re.compile(r'(?i)time|visit|phase|period|\bday\b|\bweek\b|\bmonth\b')
_BASELINE = re.compile(r'^(?:\d+ )?(?:baseline|pre|pre ?intervention|pre ?treatment|pre ?diet|pre ?supplementation|'
                       r'before (?:intervention|treatment|diet)|t ?0|(?:week|wk|day|month) ?0|0 ?(?:h|d|w|wk|weeks?|days?))'
                       r'(?: \w+)?$|.*\bbaseline\b.*')

# explicit control / none values of a factor (count as held by the column, see sample_groups)
_NONE_VALUE = re.compile(r'(?i)^(?:control|controls|ctrl|con|placebo|vehicle|veh|untreated|none|no|sham|mock|naive|'
                         r'no treatment|no supplement|unsupplemented|not supplemented)$')

# ---------------------------------------------------------------- GFOP foods
_DOSE = re.compile(r'(?<![a-z0-9])(?:\d+(?:\.\d+)?\s*(?:%|percent|g|mg|ug|kg|ml|l|mm|um|mmol|umol|ppm)'
                   r'(?:\s*(?:w\s*w|w\s*v|kg|g|day|d|bw|cattle))*|w\s*w|w\s*v|g kg|mg kg|\d+(?:\.\d+)?)(?![a-z0-9])')
_FOOD_WRAPPER = {'diet', 'diets', 'dietary', 'feed', 'feeding', 'fed', 'supplement', 'supplemented', 'supplementation',
                 'intake', 'group', 'groups', 'whole', 'dried', 'fresh', 'raw', 'the', 'a', 'an', 'of', 'with',
                 'enriched', 'rich', 'added', 'addition', 'treated', 'treatment', 'meal', 'per', 'day', 'daily',
                 'mice', 'rats', 'given', 'consumption', 'drinking', 'in', 'to'}
# GFOP classes too general for a diet, and water (in a diet factor mostly the vehicle / control); also
# not among Diet's allowed values (update_allowed_terms_from_ontologies.py)
NOT_DIET_FOODS = {'plant', 'animal', 'fungi', 'algae', 'mineral', 'supplement', 'animal feed', 'water', 'tap water',
                  'filtered water', 'bottled water'}
_SPLIT = re.compile(r'\s*(?:\+|,|;|&|\band\b|\bwith\b|\bplus\b|supplemented with|enriched with)\s*', re.I)


def gfop_index(labels):
    """norm(label) -> GFOP label, for labels whose normalised form is unique"""
    by = collections.defaultdict(set)
    for label in labels:
        key = norm(label)
        if key:
            by[key].add(label)
    return {k: next(iter(v)) for k, v in by.items() if len(v) == 1}


def _singular(word_key):
    if word_key.endswith('ies'):
        return word_key[:-3] + 'y'
    if word_key.endswith('es'):
        return word_key[:-2]
    if word_key.endswith('s'):
        return word_key[:-1]
    return word_key


def food_of(component, gfop):
    """GFOP label of a component that names one food (doses and wrapper words removed), or None"""
    text = _DOSE.sub(' ', component.lower())
    words = [w for w in norm(text).split() if w not in _FOOD_WRAPPER]
    if not words:
        return None
    key = ' '.join(words)
    for k in (key, _singular(key)):
        if k in gfop and k not in NOT_DIET_FOODS:
            return gfop[k]
    return None


# ---------------------------------------------------------------- one factor value


class Mapped:
    """What one factor value says about the diet: regimen terms, foods, fasting state, and whether the
    value is fully held by them (nothing left over: no dose, drug, duration, ...)."""
    __slots__ = ('regimens', 'foods', 'fasting', 'complete', 'none_value', 'why')

    def __init__(self):
        self.regimens, self.foods, self.fasting = [], [], None
        self.complete, self.none_value, self.why = False, False, ''


def _regimen_spans(text, diet_key, also_short=(), feeding_key=False):
    """[(term, start, end)] of the regimen wordings in a norm() text (non-overlapping, specific first);
    short wordings under a diet factor (bare 'prebiotic' and fibres also under a feeding factor), and
    those of the terms in also_short anywhere"""
    found, taken = [], [False] * (len(text) + 1)
    for term, rx, short in _REGIMENS:
        if short and not diet_key and term not in also_short and \
                not (feeding_key and term in ('prebiotic diet', 'fibre-supplemented diet')):
            continue
        for m in rx.finditer(text):
            if any(taken[m.start():m.end()]):
                continue
            for i in range(m.start(), m.end()):
                taken[i] = True
            found.append((term, m.start(), m.end()))
    return found


def _without(text, spans):
    for _, s, e in sorted(spans, key=lambda x: -x[1]):
        text = text[:s] + ' ' + text[e:]
    return text


def map_value(key, value, gfop, masst_type=None):
    """Mapped for one factor ('key: value'); masst_type: the microbiomeMASST intervention type
    ('mice diet', 'human diet', 'mice prebiotic', ...) of 'Intervention: <type>: <value>' parts."""
    out = Mapped()
    k = norm(masst_type if masst_type is not None else key)
    raw = str(value)
    v = norm(raw)
    if not v:
        out.why = 'no value'
        return out
    out.none_value = bool(_NONE_VALUE.match(v))
    if _NEVER_KEY.search(k) or _NEVER_VALUE.search(v):
        out.why = 'not a diet factor (medium / culture / tracer / sample type / ...)'
        return out
    if masst_type is not None:
        # microbiomeMASST types: 'mice diet', 'human diet', 'pig diet', 'IBD diet' (the value is the diet)
        diet_key = feeding_key = bool(re.search(r'\bdiet\b', k))
        fasting_key = False
    else:
        diet_key = bool(_DIET_KEY.search(k))
        feeding_key = bool(_FEEDING_KEY.search(k))
        fasting_key = bool(_FASTING_KEY.search(k))
    # a yes / no factor named after a diet ('HFD: 1', 'High fat diet: yes'); 'no' / '0' is not
    # necessarily the control diet, so it gives nothing and keeps the factor
    named = [s for s in _regimen_spans(k, False) if s[0] != CONTROL]
    if masst_type is None and named and not [w for w in _without(k, named).split() if w not in _FILLER]:
        if v in ('1', 'yes', 'y', 'true'):
            out.regimens, out.complete = [named[0][0]], True
            return out
        if v in ('0', 'no', 'n', 'false'):
            out.none_value, out.why = False, 'no (named diet not given)'
            return out
    undosed = norm(_DOSE.sub(' ', raw.lower()))
    # (not the words of a regimen name: 'AIN-93M' is a control diet, not 93 months)
    if _TIME.search(_without(undosed, _regimen_spans(undosed, diet_key, feeding_key=feeding_key))):
        out.why = 'time point / diet change'
        return out
    if _NEGATION.search(v) and not out.none_value:
        out.why = 'negation'
        return out
    if _PARENTS.search(v) or _PARENTS.search(k):
        out.why = "the parents' diet"
        return out

    # fasting state: 'fasted' where a component of the value is only a fasting state ('Fasting',
    # 'WT - 16h fast') and both states under a fasting factor ('Nutritional status: Fed state, ...');
    # else a candidate ('Hsc70 fed', 'WT - Fed', 'ALS: Fast') that counts only where the dataset's same
    # factor has both states (apply())
    if not _FAST_NOT.search(v) and not (diet_key and not fasting_key):
        # (a fasting diet group, intermittent / therapeutic fasting, is a regimen, not the state at sampling)
        found = set()
        for comp in [norm(c) for c in re.split(r'\s*(?:/|,|;|\+|\s-\s|_)\s*', raw) if norm(c)]:
            if _FASTED.search(comp):
                state, rx = 'fasted', _FASTED
            elif _FAST_BARE.search(comp):
                state, rx = 'fasted', _FAST_BARE
            elif _STARVED.search(comp):
                if [w for w in _STARVED.sub(' ', comp).split() if w not in _FILLER]:
                    continue
                state, rx = 'fasted', _FAST_BARE
                comp = _STARVED.sub(' fast ', comp)
            elif _NONFASTED.search(comp):
                state, rx = 'non-fasted', _NONFASTED
            else:
                continue
            left = _without(rx.sub(' ', comp), _regimen_spans(rx.sub(' ', comp), diet_key))
            alone = rx is not _FAST_BARE and not [w for w in left.split() if w not in _FILLER]
            # 'fed' alone can be a feeding treatment (Workbench rule): non-fasted needs fasted samples too
            found.add(state if (fasting_key or (alone and state == 'fasted')) else state + '?')
        states = {s.rstrip('?') for s in found}
        if len(states) == 1:
            s = next(iter(states))
            out.fasting = s if s in found else s + '?'
        out.fasting = {'fasted?': 'fast?'}.get(out.fasting, out.fasting)

    # the components of the value ('HFD + myriocin', 'Fasted, normal chow diet')
    rest, comp_terms, comp_control, other_diets = [], [], [], []
    # ('High-fat and sucrose diet' is one diet, not two components)
    joined = re.sub(r'(?i)high[\s_-]*fat\s*(?:and|&)\s*(?:high[\s_-]*)?(sugar|sucrose)', r'high fat high \1', raw)
    for comp in [c for c in _SPLIT.split(joined) if norm(c)]:
        c = norm(comp)
        if re.search(r'\bsalt\b', v):
            # 'high salt (HSD)': HSD is the high-salt diet there
            c = re.sub(r'\bhsd\b', ' ', c)
        spans = _regimen_spans(c, diet_key, feeding_key=feeding_key)
        terms = list(dict.fromkeys(t for t, _, _ in spans))
        left = _without(c, spans)
        # abbreviations of a term already named ('High-fat diet group (HF)')
        left = _without(left, _regimen_spans(left, False, also_short=set(terms)))
        other_diets += [m.group(0) for m in (_OTHER_DIET if diet_key else _OTHER_DIET_WORDED).finditer(left)]
        if out.fasting:
            left = _NONFASTED.sub(' ', _FAST_BARE.sub(' ', _STARVED.sub(' ', _FASTED.sub(' ', left))))
        # (doses and other numbers are left over: the column cannot hold them)
        words = [w for w in left.split() if w not in _FILLER] + (['(probiotic)'] if _PROBIOTIC_TOO.search(c) else [])
        if not terms and (feeding_key or _FEEDING_WORD.search(c)):
            food = food_of(comp, gfop)
            if food:
                out.foods.append(food)
                words = ['(dose)'] if _DOSE.search(comp.lower()) else []
        others = [t for t in terms if t != CONTROL]
        if CONTROL in terms and not others:
            comp_control.append(True)
        comp_terms.append(others)
        rest += words
    others = list(dict.fromkeys(t for ts in comp_terms for t in ts))
    if len(others) > 1 or (others and comp_control) or ((others or comp_control) and other_diets):
        # several diets in one value ('Mediterranean, western, bisphenol', 'NC, choline deficient high fat diet',
        # 'anti-inflammatory, western': also where one of them has no Diet term)
        out.why = 'several diets: ' + ' | '.join(others + ([CONTROL] if comp_control else []) +
                                                 [f'{d} (no Diet term)' for d in dict.fromkeys(other_diets)])
        out.fasting = None
        out.foods = []
        return out
    if others:
        out.regimens = others
    elif comp_control:
        # control diet only for an unmodified control diet ('Normal diet group (Con)', also next to a
        # separate treatment: 'AMX + std. diet', 'NCD + sesame oil' (the vehicle of an injection)), not for
        # a control diet changed itself ('supplemented with', 'deficient', a % dose: 'basal diet + 0.3% X')
        rest = [w for w in rest if w not in _CONTROL_FILLER and len(w) > 1]
        if out.foods or _MODIFIED.search(re.sub(r'\bfed with\b', ' ', v)) or re.search(r'%|w ?/ ?w|g ?/ ?kg', raw):
            if not out.foods:
                out.why = 'control diet changed itself: ' + ' '.join(rest)
                rest.append('(changed control diet)')
        else:
            out.regimens = [CONTROL]
    out.foods = list(dict.fromkeys(out.foods))
    if len(out.foods) > 1:
        out.why = 'several foods: ' + ' | '.join(out.foods)
    duration = out.fasting is not None and bool(_FAST_DURATION.search(v))
    n_items = len(out.regimens) + len(out.foods)
    out.complete = (n_items == 1 or (n_items == 0 and out.fasting is not None)) and not rest and not duration
    if not out.regimens and not out.foods and not out.fasting:
        out.why = out.why or ('no regimen / GFOP food / fasting wording' if not rest else 'not mapped: ' + ' '.join(rest))
    elif rest:
        out.why = 'extra: ' + ' '.join(rest)
    elif duration:
        out.why = 'extra: fasting duration'
    return out


def row_diet(mapped, human=False):
    """(Diet or None, note) of one row from the Mapped values of its factors: a regimen other than the
    control diet wins over a food; a food given on the control diet gives the food; several regimens
    or several foods: none"""
    regs = [t for m in mapped for t in m.regimens]
    regs = [t for t in regs if t not in HUMAN_PATTERNS or human]
    foods = list(dict.fromkeys(f for m in mapped for f in m.foods))
    other = list(dict.fromkeys(t for t in regs if t != CONTROL))
    if len(other) > 1 or (other and CONTROL in regs):
        return None, 'conflict: ' + ' | '.join(dict.fromkeys(regs))
    if other:
        return other[0], ''
    if len(foods) > 1:
        return None, 'conflict: ' + ' | '.join(foods)
    if foods:
        return foods[0], ''
    if CONTROL in regs:
        return CONTROL, ''
    return None, ''


def masst_diet_type(itype):
    """a microbiomeMASST diet intervention type ('mice diet', 'human diet', 'pig diet', 'IBD diet'; not the
    'DrugsNDiet' communities)"""
    return isinstance(itype, str) and re.search(r'\bdiet\b', itype.lower()) is not None


def study_level_diet(diets):
    """The Diet that every one of these files / rows (one value each, None / NaN for none) would get, else None.
    microbiomeMASST lists a study's diet for all its files, also for the control arm (MSV000085879
    'ketogenic', MSV000097318 'high-fat'), so a Diet from its diet interventions that every real file of a
    dataset would get, in any wording ('high-fat, dark' / 'high-fat, light'), is not used
    (masst_tables.py, and apply() for the texts of microbiomeMASST rows)."""
    values = list(diets)
    return values[0] if values and isinstance(values[0], str) and all(v == values[0] for v in values) else None


# ---------------------------------------------------------------- the pass


def subject_guard(df):
    """Series of '' (row may hold a Diet / FastingStatus value from the text) or the first reason it
    may not: kind (blank / QC / standard / reference / culture / other SampleType / cell line),
    food (FoodOntologyTerm or a food / beverage SampleType), organism (NCBITaxonomy missing or not
    Metazoa)."""
    st = df['SampleType'].astype(str)
    reason = pd.Series('', index=df.index, dtype=object)
    blankish = ~st.map(is_real_sample)
    food = (df['FoodOntologyTerm'].astype(str) != MISSING) | st.isin(['food', 'beverage']) | \
        df['SampleTypeSub1'].astype(str).str.match(r'(?i)food|beverage|reference material_animalfeed')
    culture = st.str.startswith('culture_')
    taxon_missing = df['NCBITaxonomy'].astype(str).isin([MISSING, '', 'nan'])
    not_metazoa = df['NCBIKingdom'].astype(str) != 'Metazoa'
    cell_line = df['CellLine'].astype(str) != MISSING if 'CellLine' in df.columns else pd.Series(False, index=df.index)
    other_type = st != 'animal'
    for cond, why in ((blankish, 'kind: blank / QC / standard / reference material'),
                      (culture, 'kind: culture (cells / microbes)'),
                      (food, 'food: the sample is a food'),
                      (other_type & ~blankish & ~culture & ~food, 'kind: SampleType not animal'),
                      (cell_line, 'kind: cell line'),
                      (taxon_missing, 'organism: NCBITaxonomy missing'),
                      (not_metazoa & ~taxon_missing, 'organism: not Metazoa')):
        reason = reason.where((reason != '') | ~cond, why)
    return reason


def apply(df, gfop_labels, regimen_terms):
    """Returns (df with Diet and FastingStatus updated, per-factor table, per-row log). df: the final
    table (needs ATTRIBUTE_DatasetAccession, DataSource, SampleType, SampleTypeSub1, FoodOntologyTerm,
    NCBITaxonomy, NCBIKingdom, UBERONBodyPartName, CellLine, Diet, FastingStatus, SampleGroupFreeText).
    The per-factor table has the columns sample_groups.remove_held_factors needs."""
    gfop = gfop_index(gfop_labels)
    regimens = set(regimen_terms)
    allowed_diet = regimens | set(gfop_labels)
    df = df.copy()
    guard = subject_guard(df)
    human = df['NCBITaxonomy'].astype(str).str.startswith('9606|')
    ds = df['ATTRIBUTE_DatasetAccession'].astype(str)

    # 1. every factor of every row
    recs, cache = [], {}
    for i, text, m in zip(df.index, df['SampleGroupFreeText'], masst_rows(df)):
        for key, value, itype, part in parts_of(text, m):
            ck = (key, value, itype)
            if ck not in cache:
                cache[ck] = map_value(key, value, gfop, itype)
            recs.append((i, key, value, itype, part, cache[ck]))
    fac = pd.DataFrame([(i, k, v, t, p, r.regimens, r.foods, r.fasting, r.complete, r.none_value, r.why, r)
                        for i, k, v, t, p, r in recs],
                       columns=['row', 'key', 'value', 'itype', 'part', 'regimens', 'foods', 'fasting', 'complete',
                                'none_value', 'why', 'mapped'])
    fac['dataset'] = ds.reindex(fac['row']).to_numpy()
    real = df['SampleType'].map(is_real_sample)
    fac['real'] = real.reindex(fac['row']).to_numpy().astype(bool)

    # 2. fasting: '?' candidates count only where the dataset's same factor has both states
    states = fac[fac['fasting'].notna()].groupby(['dataset', 'key'])['fasting'].agg(set)

    def resolve(f, dataset, key):
        if not isinstance(f, str) or f in ('fasted', 'non-fasted'):
            return f if isinstance(f, str) else None
        s = states.get((dataset, key), set())
        if not (s & {'fasted', 'fast?'} and s & {'non-fasted', 'non-fasted?'}):
            return None
        return 'fasted' if f == 'fast?' else 'non-fasted'
    fac['fasting'] = [resolve(f, d, k) for f, d, k in zip(fac['fasting'], fac['dataset'], fac['key'])]

    # 3. per row: Diet and FastingStatus from its factors
    new_diet = pd.Series(None, index=df.index, dtype=object)
    new_fast = pd.Series(None, index=df.index, dtype=object)
    diet_note = pd.Series('', index=df.index, dtype=object)
    for i, g in fac.groupby('row', sort=True):
        new_diet[i], diet_note[i] = row_diet(list(g['mapped']), human[i])
        fs = set(g['fasting'].dropna())
        if len(fs) == 1:
            new_fast[i] = next(iter(fs))
        elif len(fs) > 1:
            diet_note[i] += ' fasting conflict'
    # values must be allowed terms
    new_diet = new_diet.where(new_diet.isin(allowed_diet))
    # microbiomeMASST diet interventions (texts of GNPS rows): a Diet that every real row carrying one would
    # get is the study's (study_level_diet; as in masst_tables.py, which leaves out the interventions given
    # for every file)
    masst_diet = fac['real'].astype(bool) & fac['itype'].map(masst_diet_type).astype(bool)
    for _, rows in fac[masst_diet].groupby('dataset', sort=True)['row'].unique().items():
        if study_level_diet(new_diet.loc[rows]):
            diet_note.loc[rows] = diet_note.loc[rows] + ' study-level microbiomeMASST diet'
            new_diet.loc[rows] = None
    # a sample taken before the intervention started ('Visit: Baseline') was not on its group's diet yet
    if len(fac):
        base = fac.loc[fac['key'].str.contains(_TIME_KEY) & fac['value'].map(lambda x: bool(_BASELINE.match(norm(x)))),
                       'row'].unique()
        diet_note[base] = diet_note[base] + ' baseline'
        new_diet[base] = None

    old_diet = df['Diet'].astype(str)
    old_fast = df['FastingStatus'].astype(str)
    log = pd.DataFrame({'dataset': ds, 'source': df['DataSource'], 'guard': guard, 'old_diet': old_diet,
                        'text_diet': new_diet, 'old_fast': old_fast, 'text_fast': new_fast, 'note': diet_note})
    ok = guard == ''
    # a food only on a sample of the subject's body (a milk or yogurt file without a body part, listed
    # under 'human diet: milk', is more likely the food itself)
    no_part = df['UBERONBodyPartName'].astype(str).isin([MISSING, '', 'nan'])
    food_ok = ok & ~no_part
    # the column's own value (a dedicated diet field) is kept where the guards allow it; the text
    # fills only missing values, and only on rows that pass all guards
    old_ok = ~guard.str.startswith(('kind', 'food')) & (guard != 'organism: not Metazoa')
    foods = set(gfop_labels)
    keep = (old_diet.isin(regimens) & old_ok) | (old_diet.isin(foods) & food_ok) | ~old_diet.isin(allowed_diet)
    diet = old_diet.where((old_diet != MISSING) & keep, MISSING)
    diet = diet.where(diet != MISSING, new_diet.where(ok & (new_diet.isin(regimens) | food_ok)).fillna(MISSING))
    fast = old_fast.where((old_fast != MISSING) & old_ok, MISSING)
    fast = fast.where(fast != MISSING, new_fast.where(ok).fillna(MISSING))
    log['diet'] = diet
    log['fast'] = fast
    df['Diet'] = diet
    df['FastingStatus'] = fast

    # 4. held: the columns hold the factor's value (sample_groups.remove_held_factors)
    row_d = diet.reindex(fac['row']).to_numpy()
    row_f = fast.reindex(fac['row']).to_numpy()
    held, captured = [], []
    for r, f, d, fs in zip(fac['mapped'], fac['fasting'], row_d, row_f):
        items = r.regimens + r.foods
        mapped = bool(items or f)
        # (a fasting state the dataset does not confirm is not held either)
        h = mapped and r.complete and (not items or d == items[0]) and (r.fasting is None or (f is not None and fs == f))
        captured.append(h)
        # an explicit control / none value ('Supplement: Placebo') only where the row has a Diet; else the
        # factor stays in the text (sample_groups.remove_held_factors)
        held.append(h or (r.none_value and not mapped and d != MISSING))
    fac['held'] = held
    fac['captured'] = captured
    return df, fac.drop(columns=['mapped']), log
