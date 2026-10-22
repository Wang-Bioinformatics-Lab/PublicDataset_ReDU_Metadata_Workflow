"""Each repository's last update, for source_summary.tsv (sourceStatus in nf_workflow.nf).

Every step with a fallback writes one line per part it fetches (step, DataSource, status, started) to its status
file: ok; failed (it fell back, the repository keeps its rows from earlier runs); or degraded (it fell back to a
mode that still delivers this run's rows: the dataset cache queried per dataset, MassSpectrometer not filled
from the dataset descriptions); and when the step started (UTC). A repository is updated when none of its steps
failed and it delivered rows of its own this run (rows not only MASST describes, saveOlderData's own_this_run);
then its last_update is the start of this run (UTC, passed in as it is), or the start of the earliest of its
steps that started before it (with -resume, a step whose result is from an earlier session); otherwise the
last_update of the summary next to --old_redu (the run that wrote that table). That summary may be missing or
come from a workflow version without last_update: then the date is empty.

    python source_status.py source_counts.tsv previous_source_summary.tsv source_summary.tsv \
        --step_status status_*.tsv --run_start 2026-10-21T08:00:00Z
"""
import argparse
import os

import pandas as pd

STATUS_COLUMNS = ['step', 'DataSource', 'status', 'started']


def read_step_status(paths):
    """The status lines of the steps (STATUS_COLUMNS, tab-separated, no header) as one table; started is empty in
    the lines of a workflow version without it"""
    tables = [pd.read_csv(path, sep='\t', header=None, names=STATUS_COLUMNS, dtype=str, keep_default_na=False)
              for path in paths or [] if os.path.getsize(path)]
    if not tables:
        return pd.DataFrame(columns=STATUS_COLUMNS)
    return pd.concat(tables, ignore_index=True).fillna('')


def failed_steps(step_status):
    """The steps of step_status (read_step_status) that fell back: any status other than ok / degraded; sorted"""
    return sorted(set(step_status.loc[~step_status['status'].isin(['ok', 'degraded']), 'step']))


def last_update(run_start, started):
    """run_start, or the earliest start of the steps (started) before it, as it is given: with -resume, a step whose
    result is from an earlier session started then. run_start when it is no date"""
    started = pd.Series(list(started), dtype=str)
    times = pd.to_datetime(started, utc=True, errors='coerce', format='ISO8601')
    earlier = times < pd.to_datetime(run_start, utc=True, errors='coerce')
    return started[times[earlier].idxmin()] if earlier.any() else run_start


def read_last_updates(path):
    """{DataSource: last_update} of the source summary of an earlier run; {} when there is none, it cannot be
    read or it has no last_update (a summary of an earlier workflow version)"""
    try:
        previous = pd.read_csv(path, sep='\t', dtype=str, keep_default_na=False)
    except Exception as e:
        print(f"No earlier source summary in {path} ({type(e).__name__}: {e})")
        return {}
    if 'DataSource' not in previous.columns or 'last_update' not in previous.columns:
        print(f"No last update per repository in {path}")
        return {}
    return dict(zip(previous['DataSource'], previous['last_update']))


def source_summary(counts, step_status, last_updates, run_start):
    """saveOlderData's rows per source (counts, indexed by DataSource) with last_update and last_run_status of each
    repository that has steps in step_status ('updated', or 'not updated' and why); a repository without rows gets a
    row of zeros. own_this_run is used and left out."""
    summary = counts.copy()
    own = summary.pop('own_this_run') if 'own_this_run' in summary.columns else None
    repositories = sorted(set(step_status['DataSource']))
    summary = summary.reindex(pd.Index(summary.index.union(repositories), name='DataSource'), fill_value=0)
    summary['last_update'] = ''
    summary['last_run_status'] = ''
    for repository in repositories:
        steps = step_status[step_status['DataSource'] == repository]
        # (any other status than ok / degraded counts as failed; sorted: the status files come in any order)
        failed = sorted(set(steps.loc[~steps['status'].isin(['ok', 'degraded']), 'step']))
        degraded = sorted(set(steps.loc[steps['status'] == 'degraded', 'step']))
        no_rows = not failed and own is not None and own.get(repository, 0) == 0
        why = ([f"fell back: {', '.join(failed)}"] if failed else []) + (['no rows of its own this run'] if no_rows else []) \
            + ([f"degraded: {', '.join(degraded)}"] if degraded else [])
        updated = not failed and not no_rows
        summary.loc[repository, 'last_update'] = last_update(run_start, steps['started']) if updated else last_updates.get(repository, '')
        summary.loc[repository, 'last_run_status'] = ('updated' if updated else 'not updated') + (f" ({'; '.join(why)})" if why else '')
    return summary


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Each repository\'s last update for the source summary')
    parser.add_argument('source_counts', help='saveOlderData\'s rows per source (save_older_data.py --summary_output)')
    parser.add_argument('previous_summary', help='source_summary.tsv of the run that wrote --old_redu (may be missing or older)')
    parser.add_argument('output_summary')
    parser.add_argument('--step_status', nargs='*', default=None, help='status files of the steps with a fallback')
    parser.add_argument('--run_start', required=True, help='start of this run (UTC), the last update of the repositories it updates')
    args = parser.parse_args()

    counts = pd.read_csv(args.source_counts, sep='\t', index_col='DataSource', keep_default_na=False)
    summary = source_summary(counts, read_step_status(args.step_status), read_last_updates(args.previous_summary), args.run_start)
    summary.to_csv(args.output_summary, sep='\t')
    print(summary.to_string())
