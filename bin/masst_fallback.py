"""When a MASST step fell back this run (MASST_to_REDU or gnpsmatchName_masst, from their status files), the fields
MASST filled into other sources' rows in the earlier table (--old_redu) are taken from that table, after
finalizeTextColumns (final_text_columns.py), in place. Without a fallback the table is not read or written.

MASST fills missing fields of the rows other sources describe (merge_metadata.supplement_from_masst; also the
fields the later steps derive from them: the lineage of a filled NCBITaxonomy, Diet from a filled
SampleGroupFreeText, ...). The published table does not record which cells MASST filled, so a dataset is taken from
the earlier table only where this run's rows can differ from it by MASST alone (MASST is deterministic: its tables
are pinned in data/reference_data.json):
- the dataset has the same files (USI without MS file extension) in both tables;
- every cell of this run that has a value has the same value in the earlier table: this run's rows are the earlier
  rows with cells left 'missing value';
- none of these cells is of a field MASST never gives (SUBMITTER_ONLY_FIELDS of save_older_data.py but
  InternalStandardsUsed, which finalizeTextColumns derives from SampleGroupFreeText, and the identity fields
  SUPPLEMENT_NEVER of merge_metadata.py): such a value was removed upstream.
(The columns add_cache_columns.py sets next are left out. The columns ncbi_database.py computes for every row from
the pinned NCBI taxonomy and Wikidata table, DERIVED, are left out where their source column has a value this run:
there they are this run's own computation, and a difference to the earlier table comes from other reference data,
not from MASST; where the source is missing this run they go with it.)
Then the cells of the dataset that are missing this run take the earlier table's values; any other dataset is left
as this run made it (its metadata changed, or the change cannot be told apart from one). A value removed upstream
in a dataset that has not changed otherwise, in a field MASST can fill, is taken back until MASST works again: the
table does not tell it from a MASST fill.

    python masst_fallback.py merged_with_text_columns.tsv merged_metadata.tsv \\
        --masst_status status_MASST_to_REDU.tsv status_gnpsmatchName_masst.tsv
"""
import argparse
import collections
import os

import numpy as np
import pandas as pd

from merge_metadata import SUPPLEMENT_NEVER, usi_noext
from ncbi_database import DESIRED_RANKS, ISOLATION_SOURCE, ISOLATION_SOURCE_LINEAGE, write_tsv_atomic
from save_older_data import SUBMITTER_ONLY_FIELDS
from source_status import failed_steps, read_step_status

# fields a value of the earlier table in a cell that is missing this run cannot come from MASST
NEVER_FROM_MASST = (set(SUBMITTER_ONLY_FIELDS) - {'InternalStandardsUsed'}) | SUPPLEMENT_NEVER
# the columns add_cache_columns.py sets after this step (rows kept from earlier runs carry them already): not compared
# and not taken
CACHE_COLUMNS = ('classification', 'MS2spectra_count')
# column computed from the pinned reference data in every run -> its source column (ncbi_database.py)
DERIVED = {**{'NCBI' + rank.capitalize(): 'NCBITaxonomy' for rank in DESIRED_RANKS},
           'OpenTreeOfLifeTaxonomyID': 'NCBITaxonomy',
           **{column: ISOLATION_SOURCE for column in ISOLATION_SOURCE_LINEAGE.values()}}
_MISSING = ('', 'nan', 'missing value')


def _missing(values):
    """bool array: the values (object array of str) that are missing"""
    return pd.Series(values, dtype=object).fillna('').astype(str).str.strip().str.lower().isin(_MISSING).to_numpy()


def _dataset_of(usi):
    return usi.str.extract(r'^mzspec:([^:]+):', expand=False).fillna('')


def keep_earlier_masst_fields(current, earlier):
    """current: this run's table, earlier: the earlier table (both all str, keep_default_na=False). Returns (current
    with the cells of the datasets taken from the earlier table filled, a per-column Counter of the filled cells,
    a Counter of the datasets per outcome, {dataset not taken: why} of the datasets with cells to take)."""
    outcome = collections.Counter()
    filled = collections.Counter()
    if not len(current) or not len(earlier) or 'USI' not in earlier.columns:
        return current, filled, outcome, {}
    key = usi_noext(current['USI'].astype(str))
    old_key = usi_noext(earlier['USI'].astype(str))
    earlier = earlier.loc[~old_key.duplicated().to_numpy()]
    old_key = old_key[~old_key.duplicated()]
    dataset = _dataset_of(current['USI'].astype(str)).to_numpy()
    position = pd.Series(np.arange(len(earlier)), index=old_key.to_numpy())
    pos = position.reindex(key.to_numpy()).to_numpy()
    matched = ~np.isnan(pos)
    pos = np.where(matched, pos, 0).astype(int)

    columns = [c for c in current.columns if c != 'USI' and c not in CACHE_COLUMNS]
    conflict = ~matched
    gain = np.zeros(len(current), dtype=bool)
    gains = {}
    # per dataset: the columns of the rows that differ otherwise (rows), for the log
    why = collections.defaultdict(collections.Counter)
    for d, n in collections.Counter(dataset[~matched]).items():
        why[d]['files not in the earlier table'] = n
    # rows whose source of a DERIVED column has a value this run (its derived columns are not compared there)
    own_source = {src: ~_missing(current[src].to_numpy(dtype=object)) for src in set(DERIVED.values()) if src in current.columns}
    for c in columns:
        mine = current[c].to_numpy(dtype=object)
        if c in earlier.columns:
            theirs = earlier[c].to_numpy(dtype=object)[pos]
        else:
            theirs = np.full(len(current), '', dtype=object)
        differ = matched & (mine != theirs)
        if DERIVED.get(c) in own_source:
            differ &= ~own_source[DERIVED[c]]
        if not differ.any():
            continue
        i = np.flatnonzero(differ)
        mine_missing, theirs_missing = _missing(mine[i]), _missing(theirs[i])
        # a value this run that the earlier table lacks or has otherwise; an earlier value MASST cannot have given
        bad = ~mine_missing
        if c in NEVER_FROM_MASST:
            bad |= ~theirs_missing
        conflict[i[bad]] = True
        for d, n in collections.Counter(dataset[i[bad]]).items():
            why[d][c] += n
        take = i[mine_missing & ~theirs_missing & ~bad]
        if len(take):
            gain[take] = True
            gains[c] = (take, theirs[take])
    # per dataset: the same files in both tables, no row with a conflict, at least one cell to take
    rows = pd.DataFrame({'dataset': dataset, 'conflict': conflict, 'gain': gain})
    per_dataset = rows.groupby('dataset').agg(conflict=('conflict', 'any'), gain=('gain', 'any'), n=('gain', 'size'))
    old_dataset = _dataset_of(earlier['USI'].astype(str))
    n_old = old_dataset.value_counts()
    in_current = old_key.isin(set(key)).to_numpy()
    n_old_matched = old_dataset[in_current].value_counts()
    per_dataset['dropped'] = (n_old.reindex(per_dataset.index, fill_value=0)
                              - n_old_matched.reindex(per_dataset.index, fill_value=0)) > 0
    candidates = per_dataset[per_dataset['gain'] & (per_dataset.index != '')]
    taken = candidates.index[~candidates['conflict'] & ~candidates['dropped']]
    outcome['datasets_with_cells_to_take'] = len(candidates)
    outcome['taken_from_earlier_table'] = len(taken)
    outcome['not_taken_rows_differ_or_added'] = int(candidates['conflict'].sum())
    outcome['not_taken_rows_dropped_this_run'] = int((~candidates['conflict'] & candidates['dropped']).sum())
    not_taken = {d: ', '.join(f'{c} {n}' for c, n in why[d].most_common()) or 'earlier files not in this run'
                 for d in candidates.index.difference(taken)}
    in_taken = pd.Series(dataset).isin(set(taken)).to_numpy()
    current = current.copy()
    for c, (take, values) in gains.items():
        keep = in_taken[take]
        if keep.any():
            col = current[c].to_numpy(dtype=object).copy()
            col[take[keep]] = values[keep]
            current[c] = col
            filled[c] = int(keep.sum())
    return current, filled, outcome, not_taken


def main():
    parser = argparse.ArgumentParser(description='The fields MASST filled in the earlier table, when MASST fell back this run')
    parser.add_argument('table', help='this run\'s table (finalizeTextColumns), rewritten in place')
    parser.add_argument('earlier_table', help='the earlier table (--old_redu)')
    parser.add_argument('--masst_status', nargs='*', default=None,
                        help='status files of MASST_to_REDU and gnpsmatchName_masst (status_<step>.tsv)')
    args = parser.parse_args()
    failed = failed_steps(read_step_status(args.masst_status))
    if not failed:
        print('MASST steps ok: the table is left as it is')
        return
    if not os.path.exists(args.earlier_table) or os.path.getsize(args.earlier_table) == 0:
        print(f"MASST fell back ({', '.join(failed)}), no earlier table: the table is left as it is")
        return
    try:
        earlier = pd.read_csv(args.earlier_table, sep='\t', dtype=str, keep_default_na=False)
    except pd.errors.EmptyDataError:
        earlier = pd.DataFrame()
    current = pd.read_csv(args.table, sep='\t', dtype=str, keep_default_na=False)
    earlier = earlier[[c for c in earlier.columns if c == 'USI' or c in current.columns]]
    current, filled, outcome, not_taken = keep_earlier_masst_fields(current, earlier)
    print(f"MASST fell back ({', '.join(failed)}): datasets whose missing cells the earlier table fills: "
          f"{dict(outcome)}; cells taken from it: {sum(filled.values())} {dict(sorted(filled.items()))}", flush=True)
    for d, reason in sorted(not_taken.items()):
        print(f"  not taken from the earlier table: {d} (rows that differ otherwise, per column: {reason})")
    if filled:
        write_tsv_atomic(current, args.table)


if __name__ == '__main__':
    main()
