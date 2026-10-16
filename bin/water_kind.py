"""Kind of an environmental water sample whose own sample text names water but no kind of water (the
'no stated kind' case of get_enviromental_water in MWB_to_REDU.py and Metabolights2REDU.py).

water_kind_of_unstated_water(sample_text, study_text) -> 'water_surface' or None:
- the sample text names a surface fresh-water body ('river water', 'Filtered lake water'), or
- the study text (title, summary / description) names lakes, rivers, ponds, streams, creeks, ditches or
  reservoirs (ST004041: Lake Aydat, 'Water sample'; MTBLS12874: Xinglinwan Reservoir, 'Fresh water'),
and the text used names no other kind of water (sea / ocean / coastal / estuarine / brackish, ground
water, waste water / sewage / effluents, drinking / tap / bottled water, storm or rain water, pore
water): then the water samples of the study can be of another kind, and nothing is set. Sample texts of
culture media, blanks and solvents (culture, media, Milli-Q, deionized, artificial ...) are no
environmental water. Precision first: a study naming both a lake and the sea gets nothing.
"""
import html
import re

# a surface fresh-water body (whole words: not 'downstream', 'bloodstream')
SURFACE_WATER_BODY_RE = re.compile(r'(?i)\b(?:lakes?|rivers?|ponds?|streams?|creeks?|ditch(?:es)?|reservoirs?|surface waters?)\b')
# other kinds of water a text can name
OTHER_WATER_RE = re.compile(
    r'(?i)\b(?:seas?|seawaters?|sea waters?|oceans?|oceanic|marine|coastal|estuar\w*|brackish|saltwater|salt water|'
    r'groundwaters?|ground waters?|aquifers?|well waters?|boreholes?|'
    r'wastewaters?|waste waters?|sewage|sewers?|sewerage|effluents?|influents?|wwtps?|treatment plants?|'
    r'drinking waters?|tap waters?|potable|bottled waters?|mineral waters?|'
    r'stormwaters?|storm waters?|rainwaters?|rain waters?|porewaters?|pore waters?|interstitial waters?)\b')
# sample texts that are no environmental water sample
NOT_ENVIRONMENTAL_WATER_RE = re.compile(
    r'(?i)\b(?:cultures?|cultured|cultivation|media|medium|broth|blanks?|milli-?q|mq|deionized|deionised|distilled|ultrapure|'
    r'ultra-pure|solvents?|hplc|lc-?ms|artificial|artifical|synthetic|h2o)\b')
# study-text phrases with a water-body word that name no water body (removed before the study text is read):
# a supplier, a disease, viral / drug reservoirs, streams of waste, blood, data or gas
STUDY_FALSE_FRIENDS_RE = re.compile(
    r'(?i)\bcharles\s+river\b|\briver\s+blindness\b|'
    r'\b(?:hiv|viral|virus|latent|drug|cellular|tissue|reservoir\s+of)\s+(?:\w+\s+)?reservoirs?\b|\breservoirs?\s+(?:of|for)\s+'
    r'(?:hiv|viruses?|infections?|pathogens?|drugs?|cells?)\b|'
    r'\b(?:waste|blood|data|gas|air|sample|process|product|side|liquid)\s+streams?\b')


# place names with a water-body word ('Salt Lake City', 'Lake Forest', 'Red River Valley', 'Reservoir Road'):
# a capitalised water-body word followed by a capitalised place word (case-sensitive; removed first)
STUDY_PLACE_NAME_RE = re.compile(
    r'\b(?:Lakes?|Rivers?|Ponds?|Streams?|Creeks?|Reservoirs?)\s+(?:City|Town|Township|Village|County|Valley|Forest|Park|'
    r'Heights|Hills|Falls|Campus|Road|Street|Avenue|Drive|Hospital|Medical|University|College|Institute|Center|Centre)\b|'
    r'\b(?:Lake|River|Creek)\s+(?=[A-Z][a-z]+\s+(?:City|County|Valley|University|College|Hospital)\b)')


def _plain(text):
    """text without HTML tags and entities (MetaboLights descriptions), whitespace runs as one space"""
    return ' '.join(html.unescape(re.sub(r'<[^>]+>', ' ', str(text or ''))).split())


def water_kind_of_unstated_water(sample_text, study_text=''):
    """'water_surface' for a water sample of no stated kind (see the module docstring), else None"""
    sample = _plain(sample_text)
    if NOT_ENVIRONMENTAL_WATER_RE.search(sample) or OTHER_WATER_RE.search(sample):
        return None
    if SURFACE_WATER_BODY_RE.search(sample):
        return 'water_surface'
    study = STUDY_FALSE_FRIENDS_RE.sub(' ', STUDY_PLACE_NAME_RE.sub(' ', _plain(study_text)))
    if SURFACE_WATER_BODY_RE.search(study) and not OTHER_WATER_RE.search(study):
        return 'water_surface'
    return None
