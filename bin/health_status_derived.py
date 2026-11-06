"""HealthStatusDerived of the final table: 'healthy' or 'health condition' from HealthStatus, DOIDCommonName,
HPOCommonName and ComorbidityListDOIDIndex, so that healthy subjects and subjects with a disease or condition can be
compared with one column. The last edit of these columns (finalizeTextColumns, after masst_fallback.py), on every row.

- 'health condition': a disease (DOIDCommonName), a condition or phenotype (HPOCommonName), comorbidities, or a
  HealthStatus other than 'healthy' (acute illness, chronic illness, unhealthy (NOS));
- 'healthy': HealthStatus 'healthy' and none of these;
- missing value: no health information, or 'healthy' together with a disease or condition (the row contradicts
  itself). Rows that are no samples of a study's subjects carry no health fields (final_table_rules.py) and so no
  value either.
A second run on its own output changes nothing.
"""
import argparse

import pandas as pd

from ncbi_database import write_tsv_atomic

MISSING = 'missing value'
COLUMN = 'HealthStatusDerived'
_CONDITION_COLUMNS = ('DOIDCommonName', 'HPOCommonName', 'ComorbidityListDOIDIndex')
_ILL = {'acute illness', 'chronic illness', 'unhealthy (NOS)'}


def _filled(s):
    return ~s.fillna('').astype(str).str.strip().str.lower().isin(['', 'nan', MISSING])


def derive(df):
    """(df with HealthStatusDerived, stats)"""
    df = df.copy()
    status = df['HealthStatus'].fillna(MISSING).astype(str) if 'HealthStatus' in df.columns \
        else pd.Series(MISSING, index=df.index)
    condition = pd.Series(False, index=df.index)
    for c in _CONDITION_COLUMNS:
        if c in df.columns:
            condition |= _filled(df[c])
    healthy = status.eq('healthy')
    ill = status.isin(_ILL) | (condition & ~healthy)
    out = pd.Series(MISSING, index=df.index, dtype=object)
    out[ill] = 'health condition'
    out[healthy & ~condition] = 'healthy'
    if COLUMN in df.columns:
        df[COLUMN] = out
    else:   # next to HealthStatus (the order of allowed_terms.json)
        df.insert(df.columns.get_loc('HealthStatus') + 1 if 'HealthStatus' in df.columns else len(df.columns), COLUMN, out)
    stats = {'healthy': int((out == 'healthy').sum()), 'health condition': int((out == 'health condition').sum()),
             'healthy_with_condition_(left_missing)': int((healthy & condition).sum())}
    return df, stats


def main():
    parser = argparse.ArgumentParser(description='HealthStatusDerived of the final table (in place)')
    parser.add_argument('table_tsv')
    args = parser.parse_args()
    df = pd.read_csv(args.table_tsv, sep='\t', dtype=str, keep_default_na=False)
    df, stats = derive(df)
    print(f"[HealthStatusDerived] {stats}", flush=True)
    write_tsv_atomic(df, args.table_tsv)


if __name__ == '__main__':
    main()
