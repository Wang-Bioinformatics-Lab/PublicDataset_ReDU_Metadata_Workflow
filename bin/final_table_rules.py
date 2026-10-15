"""Rules applied to every row of the final table, after the NCBI lineage (enrich_ncbi_information) and
before the text columns (final_text_columns.py), so also to the rows saveOlderData keeps from earlier runs:
- personal_care_products: personal care products are SampleType 'personal care product' with the most
  specific product category
- life_stage_humans_only: LifeStage (human age brackets) only on rows of humans
- sex_from_body_part: BiologicalSex of mammal samples with no recorded sex from a body part only one
  sex has (UBERON_sex_specific_parts.csv of prepare_ontologies.py)
Values are only filled or renamed, never removed (but LifeStage on rows not of humans); a second pass
changes nothing.

    python final_table_rules.py merged_with_ncbi.tsv merged_with_rules.tsv \\
        --path_to_sex_specific_parts UBERON_sex_specific_parts.csv
"""
import argparse
import collections

import pandas as pd

from ncbi_database import write_tsv_atomic
# (the food track's marker of a digest of a food: in vitro digestion, colonic fermentation)
from REDU_conversion_functions import _FOOD_DIGEST_RE

MISSING = 'missing value'
# SampleTypeSub1 of personal care products (GNPS template, personalCareProductMASST)
PERSONAL_CARE_SUB1 = 'reference material_personalcareproduct'
# LifeStage (humans only) of a subject that does not lactate: the row of a milk sample then
# describes the infant who drank it
CHILD_LIFE_STAGES = {'Infancy (<2 yrs)', 'Early Childhood (2 yrs < x <=8 yrs)'}
# NCBI taxon of Homo sapiens (NCBITaxonomy is '<taxon id>|<name>')
HUMAN_TAXON_ID = '9606'


def _missing(series):
    return series.fillna('').astype(str).str.strip().str.lower().isin(('', 'nan', MISSING))


def personal_care_products(df):
    """A personal care product is not a reference material: SampleType 'reference material' with
    SampleTypeSub1 reference material_personalcareproduct (GNPS metadata files, MASST rows of earlier
    runs) becomes 'personal care product'; a ProductCategory path ('Personal Care > Cosmetics > Skin
    Care > Sunscreen', earlier runs) keeps its last, most specific category. Returns the frame and a
    Counter."""
    stats = collections.Counter()
    if 'SampleType' in df.columns and 'SampleTypeSub1' in df.columns:
        product = (df['SampleType'] == 'reference material') & (df['SampleTypeSub1'] == PERSONAL_CARE_SUB1)
        df.loc[product, 'SampleType'] = 'personal care product'
        stats['SampleType_personal_care_product'] = int(product.sum())
    if 'ProductCategory' in df.columns:
        path = df['ProductCategory'].fillna('').astype(str)
        is_path = path.str.contains('>', regex=False)
        df.loc[is_path, 'ProductCategory'] = path[is_path].map(
            lambda p: ([x.strip() for x in p.split('>') if x.strip()] or [MISSING])[-1])
        stats['ProductCategory_path_to_last_category'] = int(is_path.sum())
    return df, stats


def life_stage_humans_only(df):
    """LifeStage holds human age brackets ('Infancy (<2 yrs)', ...): 'missing value' on every row
    whose organism is not Homo sapiens, i.e. NCBITaxonomy not taxon 9606 ('9606|Homo sapiens') and
    NCBISpecies (the lineage of ncbi_database.py) not Homo sapiens (its subspecies, e.g.
    '63221|Homo sapiens neanderthalensis'); a row without an organism is not taken for human. GNPS
    sheets give the submitter's LifeStage for any organism (mice, cats, fish, bacteria, offices),
    Workbench a stated one ('newborn', 'adolescent') for any organism. Run before
    sex_from_body_part, which reads LifeStage (a second pass then changes nothing). Returns the
    frame and a Counter."""
    stats = collections.Counter()
    if 'LifeStage' not in df.columns or 'NCBITaxonomy' not in df.columns:
        return df, stats
    human = df['NCBITaxonomy'].fillna('').astype(str).str.split('|', n=1).str[0].str.strip() == HUMAN_TAXON_ID
    if 'NCBISpecies' in df.columns:
        human |= df['NCBISpecies'] == 'Homo sapiens'
    clear = ~human & (df['LifeStage'] != MISSING)
    source = df['DataSource'] if 'DataSource' in df.columns else pd.Series('', index=df.index)
    for s, n in source[clear].value_counts().sort_index().items():
        stats[f'LifeStage_not_human_cleared {s}'] = int(n)
    stats['LifeStage_not_human_cleared'] = int(clear.sum())
    df.loc[clear, 'LifeStage'] = MISSING
    return df, stats


def sex_from_body_part(df, parts):
    """Fills a missing BiologicalSex with the sex of the body part (parts: prepare_ontologies.py's
    UBERON_sex_specific_parts.csv; UBERONBodyPartName = its Label), only on rows of
    - mammals (NCBIClass Mammalia; other animals include hermaphrodites, e.g. oysters and many fish)
      sampled as SampleType 'animal' (not blanks, QCs, standards, foods, cultures), without a
      FoodOntologyTerm (milk analysed as food) or CellLine;
    - milk: not where the row's subject is a child (LifeStage), whose sex the row then gives;
    - not on digests (in vitro digestion, colonic fermentation: _FOOD_DIGEST_RE of
      REDU_conversion_functions.py in SampleGroupFreeText), nor on the rows of the same body part of a
      dataset that has such rows (MTBLS6768: in vitro digests of cow-milk infant formula, one of them a
      file listed for two samples, whose text the merge left empty): no animal was sampled;
    - not in a dataset whose rows of the body part record another sex (the column then describes
      someone else, e.g. the infants of a breast milk study).
    Recorded values are never changed: rows whose recorded sex differs from the body part's are only
    counted. Returns the frame and a Counter."""
    stats = collections.Counter()
    if len(df) == 0 or any(c not in df.columns for c in ('BiologicalSex', 'UBERONBodyPartName', 'NCBIClass', 'SampleType')):
        return df, stats
    by_label = parts.groupby('Label')['BiologicalSex'].agg(lambda s: s.iloc[0] if s.nunique() == 1 else '')
    milk_labels = set(parts.loc[parts['Is Milk'].astype(str).str.lower() == 'true', 'Label'])
    part = df['UBERONBodyPartName'].fillna('').astype(str)
    sex = part.map(by_label).fillna('')
    recorded = df['BiologicalSex'].fillna('').astype(str).str.strip().str.lower()
    unrecorded = _missing(df['BiologicalSex'])
    dataset = df['ATTRIBUTE_DatasetAccession'] if 'ATTRIBUTE_DatasetAccession' in df.columns else pd.Series('', index=df.index)
    other_sex = (sex != '') & ~unrecorded & (recorded != sex)
    for p, n in part[other_sex].value_counts().sort_index().items():
        stats[f'recorded_sex_differs_from_body_part_(kept) {p}'] = int(n)
    conflicting = set(zip(dataset[other_sex], part[other_sex]))
    digest = df['SampleGroupFreeText'].fillna('').astype(str).str.contains(_FOOD_DIGEST_RE) \
        if 'SampleGroupFreeText' in df.columns else pd.Series(False, index=df.index)
    digested = set(zip(dataset[digest], part[digest]))

    candidate = (sex != '') & unrecorded
    reasons = [
        ('not_a_mammal', df['NCBIClass'] != 'Mammalia'),
        ('not_SampleType_animal', df['SampleType'] != 'animal'),
        ('food', ~_missing(df['FoodOntologyTerm']) if 'FoodOntologyTerm' in df.columns else None),
        ('cell_line', ~_missing(df['CellLine']) if 'CellLine' in df.columns else None),
        ('milk_of_a_child_subject', part.isin(milk_labels) & df['LifeStage'].isin(CHILD_LIFE_STAGES)
         if 'LifeStage' in df.columns else None),
        ('dataset_records_another_sex', pd.Series([k in conflicting for k in zip(dataset, part)], index=df.index)
         if conflicting else None),
        ('digest', digest | pd.Series([k in digested for k in zip(dataset, part)], index=df.index) if digested else None),
    ]
    fill = candidate.copy()
    for reason, mask in reasons:
        if mask is None:
            continue
        stats[f'not_filled_{reason}'] = int((fill & mask).sum())
        fill &= ~mask
    df.loc[fill, 'BiologicalSex'] = sex[fill]
    for (p, s), n in pd.DataFrame({'p': part[fill], 's': sex[fill]}).value_counts().sort_index().items():
        stats[f'filled {s} {p}'] = int(n)
    stats['filled_rows'] = int(fill.sum())
    return df, stats


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Rules applied to every row of the final table')
    parser.add_argument('input_tsv')
    parser.add_argument('output_tsv')
    parser.add_argument('--path_to_sex_specific_parts', default=None,
                        help='UBERON_sex_specific_parts.csv of prepare_ontologies.py (BiologicalSex from the body part)')
    args = parser.parse_args()

    # all cells as the previous step wrote them (no NA parsing: e.g. 'NA' stays 'NA')
    df = pd.read_csv(args.input_tsv, sep='\t', dtype=str, keep_default_na=False, na_filter=False, low_memory=False)
    df, stats = personal_care_products(df)
    df, life_stage_stats = life_stage_humans_only(df)
    stats.update(life_stage_stats)
    if args.path_to_sex_specific_parts:
        parts = pd.read_csv(args.path_to_sex_specific_parts, dtype=str, keep_default_na=False)
        df, sex_stats = sex_from_body_part(df, parts)
        stats.update(sex_stats)
    for key, value in sorted(stats.items()):
        print(f"  {key}: {value}", flush=True)
    write_tsv_atomic(df, args.output_tsv)
