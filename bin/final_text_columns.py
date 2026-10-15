"""The last edits of the final table's SampleGroupFreeText and of the columns its factors are moved
into, on every row (this run's rows and the rows save_older_data.py kept from earlier runs alike),
after the NCBI lineage is known (ncbi_database.py) and the rules of final_table_rules.py (personal
care products, BiologicalSex), before the cache columns (add_cache_columns.py).

1. Diet / FastingStatus from the texts (diet_columns.py), InternalStandardsUsed from them
   (istd_exposure_columns.py; Exposure is only checked, the converters fill it);
2. the factors these columns hold are removed (sample_groups.remove_held_factors), per dataset and
   factor and separately for each column's factors;
3. the group-consistency rule (sample_groups.consistent_sample_groups) as the last edit.
A second run on its own output changes nothing.
"""
import argparse
import collections
import json

import pandas as pd

import diet_columns
import istd_exposure_columns
from ncbi_database import write_tsv_atomic
from sample_groups import MISSING, consistent_sample_groups, remove_held_factors

# (rows of earlier runs' tables can lack these columns or carry empty values)
TEXT_COLUMNS = ('Diet', 'FastingStatus', 'SampleGroupFreeText', 'InternalStandardsUsed', 'Exposure')


def finalize(df, allowed_terms, exposure_table):
    """df: the final table as text (all columns str). Returns the updated frame."""
    df = df.copy()
    for c in TEXT_COLUMNS:
        if c not in df.columns:
            df[c] = MISSING
        df[c] = df[c].fillna(MISSING).replace('', MISSING)
    before = {c: df[c].copy() for c in ('Diet', 'FastingStatus', 'InternalStandardsUsed', 'SampleGroupFreeText')}
    foods = allowed_terms['FoodOntologyTerm']['allowed_values']
    regimens = [t for t in allowed_terms['Diet']['allowed_values'] if t not in set(foods)]
    df, diet_fac, diet_log = diet_columns.apply(df, foods, regimens)
    df, ie_fac, ie_stats = istd_exposure_columns.apply(df, allowed_terms, exposure_table)
    texts, removed_diet = remove_held_factors(df['SampleGroupFreeText'], diet_fac)
    texts, removed_ie = remove_held_factors(texts, ie_fac)
    df['SampleGroupFreeText'] = texts
    kinds = collections.Counter(ie_fac.drop_duplicates(['dataset', 'key']).set_index(['dataset', 'key'])['kind']
                                .reindex(sorted(removed_ie)).tolist()) if removed_ie else {}
    print(f"[final text columns] Diet / FastingStatus factors removed from the texts: {len(removed_diet)} "
          f"(dataset x factor); internal-standard / exposure factors removed: {dict(kinds)}; {dict(ie_stats)}; "
          f"guard reasons of rows with a Diet value from the text: "
          f"{diet_log.loc[diet_log['text_diet'].notna(), 'guard'].replace('', 'ok').value_counts().to_dict()}", flush=True)
    df = consistent_sample_groups(df)
    changed = {c: int((df[c] != v).sum()) for c, v in before.items()}
    print(f"[final text columns] rows changed: {changed}", flush=True)
    return df


def main():
    parser = argparse.ArgumentParser(description='Diet / FastingStatus / InternalStandardsUsed from the sample group '
                                                 'texts, and the SampleGroupFreeText rules, for the final table')
    parser.add_argument('input_tsv')
    parser.add_argument('output_tsv')
    parser.add_argument('--path_to_allowed_term_json', required=True)
    parser.add_argument('--path_to_exposure_csv', required=True, help='Exposure_ontology.csv of prepare_ontologies.py')
    args = parser.parse_args()
    with open(args.path_to_allowed_term_json) as f:
        allowed_terms = json.load(f)
    exposure_table = pd.read_csv(args.path_to_exposure_csv, dtype=str, keep_default_na=False)
    df = pd.read_csv(args.input_tsv, sep='\t', dtype=str, keep_default_na=False)
    write_tsv_atomic(finalize(df, allowed_terms, exposure_table), args.output_tsv)


if __name__ == '__main__':
    main()
