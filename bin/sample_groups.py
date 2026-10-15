"""SampleGroupFreeText of the final table: the text helpers shared by the passes that move parts
of the text into columns of their own (diet_columns.py: Diet / FastingStatus,
istd_exposure_columns.py: InternalStandardsUsed, and the check against Exposure), the removal of
the factors those columns hold, and the group-consistency rule (the last edit of the text,
final_text_columns.py).

A text is 'Factor: value' parts joined by '; '. On GNPS rows the texts come from microbiomeMASST
('Intervention: <type>: <value>'); there each intervention type is a factor of its own.

Removal of moved factors (remove_held_factors), per dataset (ATTRIBUTE_DatasetAccession) and
factor: a factor is removed from all rows of the dataset only if on every real sample carrying it
its value is held by the column it was moved to (an explicit control / none value: the column has a
value on that row, e.g. 'control diet'), and at least one of its values is held; otherwise it stays
on all rows (values the column cannot hold are not lost, and the groups stay comparable).

Group consistency (consistent_sample_groups), per dataset, over its real samples (is_real_sample).
A sample's label is its parts plus its values of the columns parts are moved into (LABEL_COLUMNS:
a sample whose only factor now sits in Diet still belongs to its group, it just has nothing left
to show; not InternalStandardsUsed, a property of the analysis, not a group):
- blanks, QCs, standards, reference material and personal care products get 'missing value': their
  SampleType says what they are;
- a part (or column value) that every labelled real sample has is a property of the study, not a
  group (MTBLS2776 'Treatment: caffeine', the antibiotics of MSV000083777): removed from all rows
  of the dataset. Repeated until no such part is left, because a sample whose only parts were
  study-level has no label any more; but a removal that would leave fewer than 2 different labels
  where there were more is not made, it would wipe the design (ST001714: 'Parasite Stage' on every
  sample, 'Concentration: 10X IC50' on the treated ones only: both stay);
- the texts are kept where the real samples carry at least 2 different labels and the texts on
  average at least 2 samples per text (as in the converters); otherwise every row of the dataset gets
  'missing value'. Real samples the source gives no group for keep 'missing value'.
"""
import collections
import re

import pandas as pd

MISSING = 'missing value'
_NA_TEXT = {'', 'nan', 'none', 'missing value', 'not applicable', 'not specified'}
_NOT_REAL_TYPES = {'pool_QC', 'reference material', 'personal care product', 'molecular_standards'}
# columns that factors of the texts are moved into and that describe the sampled subject: part of
# a sample's group label
LABEL_COLUMNS = ('Diet', 'FastingStatus', 'Exposure')


def norm(text):
    """lower case, every non-alphanumeric run -> one space"""
    return ' '.join(re.sub(r'[^a-z0-9]+', ' ', str(text).lower()).split())


def is_real_sample(sample_type):
    """not a blank, QC, standard, reference material or personal care product (a missing SampleType
    counts as a sample)"""
    t = '' if pd.isna(sample_type) else str(sample_type).strip()
    return 'blank' not in t.lower() and t not in _NOT_REAL_TYPES


def masst_rows(df):
    """rows whose texts are microbiomeMASST interventions (GNPS rows)"""
    if 'DataSource' not in df.columns:
        return pd.Series(False, index=df.index)
    return df['DataSource'].astype(str).eq('GNPS')


def parts_of(text, masst=False):
    """[(factor, value, intervention type or None, part)] of a SampleGroupFreeText; on microbiomeMASST
    rows (masst) the factor of 'Intervention: <type>: <value>' is 'Intervention: <type>'"""
    out = []
    if not isinstance(text, str) or text in ('', MISSING):
        return out
    for part in text.split('; '):
        name, _, value = part.partition(': ')
        if masst and name == 'Intervention':
            itype, _, ivalue = value.partition(': ')
            out.append((f'Intervention: {itype}', ivalue, itype, part))
        else:
            out.append((name, value, None, part))
    return out


def remove_held_factors(texts, fac):
    """The removal rule. fac: one row per factor part of a row's text, with the columns row (index
    of texts), dataset, key (the factor), part, real, held (the column holds the value, or it is a
    control / none value on a row whose column has a value) and captured (held by the column, not only
    a control / none value).
    Returns (texts without the removed factors, set of removed (dataset, key))."""
    if not len(fac):
        return texts, set()
    real = fac[fac['real']].groupby(['dataset', 'key'], sort=True)['held'].all()
    captured = fac.groupby(['dataset', 'key'], sort=True)['captured'].any()
    # (a factor only blanks / QCs carry is left to consistent_sample_groups, which clears their texts)
    remove = {k for k, ok in real.items() if ok and captured.get(k, False)}
    drop = fac[[(d, k) in remove for d, k in zip(fac['dataset'], fac['key'])]].groupby('row', sort=True)['part'].agg(set)
    texts = texts.astype(object).copy()
    for i, ps in drop.items():
        kept = [p for p in str(texts[i]).split('; ') if p not in ps]
        texts[i] = '; '.join(kept) if kept else MISSING
    return texts, remove


def part_identity(part, masst=False):
    """What a 'Factor: value' part says, case and spacing aside. microbiomeMASST parts (masst:
    GNPS rows) by their value under any type (the type where there is none): MSV000083777 lists
    its antibiotic cocktail under 'mice bile acids' for 120 files and under 'mice drug' for one."""
    p = ' '.join(part.lower().split())
    if masst and p.startswith('intervention: '):
        return 'intervention: ' + p[len('intervention: '):].split(': ', 1)[-1]
    return p


def study_level_parts(labels):
    """(study, kept): the part identities (labels: one set per real sample) that every labelled
    sample has, also those that every sample still labelled without them has; but not the parts
    whose removal would leave fewer than 2 different labels where there were more (kept)"""
    study = set()
    while True:
        labelled = [s - study for s in labels if s - study]
        common = set.intersection(*labelled) if labelled else set()
        if not common:
            return study, set()
        if len({frozenset(s - common) for s in labelled if s - common}) < 2 <= len({frozenset(s) for s in labelled}):
            return study, common
        study |= common


def consistent_sample_groups(df, column='SampleGroupFreeText', dataset_column='ATTRIBUTE_DatasetAccession',
                             label_columns=LABEL_COLUMNS):
    """See the module docstring. Never removes rows; a second pass changes nothing."""
    if column not in df.columns or len(df) == 0:
        return df
    stats = collections.Counter()
    text = df[column].fillna('').astype(str).str.strip()
    text = text.where(~text.str.lower().isin(_NA_TEXT), '')
    real = df['SampleType'].map(is_real_sample)
    masst = masst_rows(df)
    stats['non_real_rows_cleared'] = int(((~real) & (text != '')).sum())
    out = text.where(real, '')
    in_scope = real & df[dataset_column].isin(df.loc[real & (text != ''), dataset_column].unique())
    label_columns = [c for c in label_columns if c in df.columns]
    for _, idx in text[in_scope].groupby(df.loc[in_scope, dataset_column], sort=True).groups.items():
        parts = text.loc[idx].map(lambda s: s.split('; ') if s else [])
        ids = [[part_identity(p, m) for p in ps] for ps, m in zip(parts, masst.loc[idx])]
        # (column values as parts of the label; 'column\t...' never equals a text part)
        values = [[f'{c}\t{v}' for c, v in zip(label_columns, vs) if not pd.isna(v) and str(v).lower() not in _NA_TEXT]
                  for vs in zip(*[df.loc[idx, c] for c in label_columns])] if label_columns else [[]] * len(idx)
        labels = [set(i) | set(v) for i, v in zip(ids, values)]
        study, not_removed = study_level_parts(labels)
        kept = pd.Series(['; '.join(p for p, i in zip(ps, pi) if i not in study) for ps, pi in zip(parts, ids)],
                         index=idx) if study else text.loc[idx]
        stats['study_level_parts_removed'] += len(study & {i for pi in ids for i in pi})
        if not_removed & {i for pi in ids for i in pi}:
            stats['datasets_study_level_parts_kept (else fewer than 2 groups)'] += 1
        n_groups = len({frozenset(s - study) for s in labels if s - study})
        texts = kept[kept != '']
        if n_groups < 2 or texts.nunique() > len(texts) / 2:
            stats['datasets_cleared (fewer than 2 groups)' if n_groups < 2 else 'datasets_cleared (texts > half of the samples)'] += 1
            stats['real_rows_cleared'] += int((kept != '').sum())
            out.loc[idx] = ''
        else:
            stats['datasets_kept'] += 1
            stats['real_rows_without_group_text_kept'] += int((kept == '').sum())
            out.loc[idx] = kept
    print(f"[SampleGroupFreeText] group consistency: {dict(sorted(stats.items()))}")
    df = df.copy()
    df[column] = out.where(out != '', MISSING)
    return df
