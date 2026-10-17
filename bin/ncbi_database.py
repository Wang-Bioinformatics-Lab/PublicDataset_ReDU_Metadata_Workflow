"""Adds the NCBI lineage of each row's NCBITaxonomy (NCBISuperkingdom, NCBIKingdom, NCBIPhylum, NCBIClass,
NCBIOrder, NCBIFamily, NCBIGenus, NCBISpecies), fills SampleType from NCBIKingdom where it is 'missing value',
and adds OpenTreeOfLifeTaxonomyID (Wikidata NCBI -> OTT table sparql_oto_ncbi.tsv in the working directory).

The lineage comes from the NCBI taxonomy pinned in data/reference_data.json (nodes.dmp, names.dmp and
merged.dmp of the taxdmp archive, set up by bin/fetch_reference_data.py), the release every other step uses.
Earlier versions used ete3.NCBITaxa(), i.e. whichever taxonomy copy ete3 had in ~/.etetoolkit (or the newest
NCBI release, downloaded on first use); for the same taxdump the columns are the same.

- Taxon ID of a row: the first number in the first '|' field of NCBITaxonomy ('562|Escherichia coli' -> 562).
- A sample of several organisms has one 'ID|Name' term per organism, joined by '||' (MULTI_SEP; the ReDU
  template's convention: '9606|Homo sapiens||5833|Plasmodium falciparum'). Its lineage columns and
  OpenTreeOfLifeTaxonomyID hold one value per organism in the same order, joined by '||'
  ('Metazoa||missing value'; aligned(); a column without a value for any of the organisms is a single
  'missing value'), and SampleType is filled only when all organisms give the same one.
- Each lineage column holds the scientific name (names.dmp) of the deepest ancestor at that rank, the taxon
  itself included. NCBI renamed the top rank 'superkingdom' to 'domain' in 2025; 'domain' fills
  NCBISuperkingdom.
- Merged IDs (merged.dmp: NCBI merged the taxon into another one) get the lineage of the taxon they were
  merged into, as with ete3. NCBITaxonomy itself is left as it is.
- IDs neither in nodes.dmp nor in merged.dmp (deleted by NCBI, see delnodes.dmp, or never valid) get
  'missing value' in all lineage columns, as with ete3. Both kinds are listed in the log.
- Names are used as they are in names.dmp. (ete3 also cut '"' at either end of a name, which changed 4 names
  in the 2026-10-01 release, e.g. 'Exapion sp. "malvae"' -> 'Exapion sp. "malvae'.)

Only two small arrays over all ~3M taxa are held (parent and rank per taxon ID); names are read for the
taxa that are needed only.

    python ncbi_database.py input.tsv output.tsv --path_to_ncbi_nodes_dmp nodes.dmp \\
        --path_to_ncbi_names_dmp names.dmp --path_to_ncbi_merged_dmp merged.dmp
"""
import argparse
import csv
import os
import re
import tempfile
from array import array

import numpy as np
import pandas as pd


DESIRED_RANKS = ['superkingdom', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']
# NCBI renamed the 'superkingdom' rank to 'domain' (2025); keep the old column filled
RANK_ALIASES = {'domain': 'superkingdom'}
ROOT = 1
MAX_DEPTH = 1000  # no NCBI lineage comes close; stops a loop in a broken nodes.dmp


# ---- several organisms in one NCBITaxonomy value ----------------------------------------------------------
# A sample of several organisms (a co-culture, a host with its symbiont or parasite, a mixture, a pool of
# several species) has one allowed 'ID|Name' term per organism, joined by MULTI_SEP in the order the source
# names them, each taxon once ('6613|Euprymna scolopes||312309|Aliivibrio fischeri ES114'). The derived columns
# (NCBIRank, NCBIDivision, the lineage, OpenTreeOfLifeTaxonomyID) hold one value per organism, aligned with
# the terms. A test that the organism lies in a group (Metazoa, Mammalia, Viridiplantae, ...) holds for such a
# value only when it holds for every part (all_parts); a test that it is one taxon (human) never holds
# (single_taxid is '' for it). These helpers are used by every step that reads NCBITaxonomy.
MULTI_SEP = '||'
MISSING_VALUE = 'missing value'
# a submitted value: '||' between terms, or ', ' / '; ' before the next 'ID|' ('43114|A b, 1008953|C d'); a comma
# inside a name ('457404|Fusobacterium ulcerans, 12_1B') is no separator
_SUBMITTED_SEP = re.compile(r'\s*\|\|\s*|\s*[,;]\s*(?=\d+\s*\|)')
# parts of a submitted value that give no organism ('9606|Homo sapiens||missing value': the human only)
_PLACEHOLDER_PARTS = frozenset({'nan', 'na', 'n/a', 'none', 'null', MISSING_VALUE, 'not applicable', 'not available'})


def _missing_text(value):
    if value is None or (isinstance(value, float) and value != value):
        return True
    return str(value).strip().lower() in ('', 'nan', MISSING_VALUE)


def taxonomy_parts(value):
    """The parts of an NCBITaxonomy value (or of a derived column's aligned value) split at MULTI_SEP,
    stripped, empty parts dropped; [] for a missing value."""
    if _missing_text(value):
        return []
    return [p.strip() for p in str(value).split(MULTI_SEP) if p.strip()]


def submitted_taxonomy_parts(value):
    """The parts of a submitted NCBITaxonomy value: split at '||', and at ',' / ';' directly before the next
    'ID|' (several 'ID|Name' pairs); parts that are a placeholder ('missing value', 'NA' ...) dropped; [] for a
    missing value."""
    if _missing_text(value):
        return []
    parts = [p.strip() for p in _SUBMITTED_SEP.split(str(value)) if p.strip()]
    if len(parts) < 2:
        return parts
    return [p for p in parts if p.lower() not in _PLACEHOLDER_PARTS]


def is_multi(value):
    """Whether an NCBITaxonomy value names several organisms."""
    return isinstance(value, str) and MULTI_SEP in value and len(taxonomy_parts(value)) > 1


def taxon_ids(value):
    """The taxon ID (text before the first '|') of each part of a value."""
    return [p.split('|', 1)[0].strip() for p in taxonomy_parts(value)]


def single_taxid(value):
    """The taxon ID of a value of one organism ('9606|Homo sapiens' -> '9606'); '' for a value of several
    organisms and for a missing one."""
    parts = taxonomy_parts(value)
    return parts[0].split('|', 1)[0].strip() if len(parts) == 1 else ''


def join_parts(terms, key=None):
    """The terms (each 'ID|Name' or itself a MULTI_SEP value) as one NCBITaxonomy value: each taxon once (the
    first; key(term) gives a term's taxon, default its ID), in the given order, joined by MULTI_SEP; a single
    term when one taxon is left, '' for none."""
    out, seen = [], set()
    for term in terms:
        for part in taxonomy_parts(term):
            k = key(part) if key is not None else part.split('|', 1)[0].strip()
            if k in seen:
                continue
            seen.add(k)
            out.append(part)
    return MULTI_SEP.join(out)


def aligned(values):
    """One value per part of an NCBITaxonomy value joined by MULTI_SEP, a part without a value as
    'missing value'; 'missing value' when no part has one."""
    values = [MISSING_VALUE if _missing_text(v) else str(v) for v in values]
    if not values or all(v == MISSING_VALUE for v in values):
        return MISSING_VALUE
    return MULTI_SEP.join(values)


def all_parts(value, predicate):
    """Whether predicate holds for every part of a value (an NCBITaxonomy value or an aligned derived column);
    False for a missing value. For a value of one organism: predicate(value)."""
    parts = taxonomy_parts(value)
    return bool(parts) and all(predicate(p) for p in parts)


def _per_value(series, function, dtype):
    memo, out = {}, []
    for v in series.tolist():
        k = v if isinstance(v, str) else None
        if k not in memo:
            memo[k] = function(k)
        out.append(memo[k])
    return pd.Series(out, index=series.index, dtype=dtype)


def all_parts_series(series, predicate):
    """all_parts per cell of a Series (evaluated once per distinct value; cells that are no text: False)."""
    return _per_value(series, lambda v: all_parts(v, predicate), bool)


def single_taxid_series(series):
    """single_taxid per cell of a Series (evaluated once per distinct value; cells that are no text: '')."""
    return _per_value(series, single_taxid, object)


class Taxonomy:
    """Lineages from NCBI nodes.dmp, names.dmp (scientific names) and merged.dmp."""

    def __init__(self, nodes_dmp, names_dmp, merged_dmp):
        self.names_dmp = names_dmp

        # nodes.dmp: tax_id | parent tax_id | rank | ...  (fields end in '\t|')
        taxids, parents, codes = array('i'), array('i'), array('h')
        rank_code = {}
        with open(nodes_dmp, encoding='utf-8') as handle:
            for line in handle:
                fields = line.split('|', 3)
                taxids.append(int(fields[0]))
                parents.append(int(fields[1]))
                rank = fields[2].strip()
                code = rank_code.get(rank)
                if code is None:
                    code = rank_code[rank] = len(rank_code)
                codes.append(code)
        if not taxids:
            raise ValueError(f"{nodes_dmp} has no taxa")
        taxids = np.frombuffer(taxids, dtype=np.int32)
        parents = np.frombuffer(parents, dtype=np.int32)
        size = int(max(taxids.max(), parents.max())) + 1
        parent = np.full(size, -1, dtype=np.int32)
        parent[taxids] = parents
        rank = np.full(size, -1, dtype=np.int16)
        rank[taxids] = np.frombuffer(codes, dtype=np.int16)
        if (parent[parents] < 0).any():
            raise ValueError(f"{nodes_dmp}: a parent taxon is not in the file")
        if parent[ROOT] != ROOT:
            raise ValueError(f"{nodes_dmp}: taxon {ROOT} is not the root")
        # array.array: plain Python ints when indexed (faster than numpy scalars in the walk below)
        self.parent = array('i', parent.tobytes())
        self.rank = array('h', rank.tobytes())
        # rank code -> lineage column (index in DESIRED_RANKS), or None
        self.column = [None] * len(rank_code)
        for name, code in rank_code.items():
            name = RANK_ALIASES.get(name, name)
            if name in DESIRED_RANKS:
                self.column[code] = DESIRED_RANKS.index(name)
        self.n_taxa = len(taxids)

        # merged.dmp: old tax_id | new tax_id
        self.merged = {}
        with open(merged_dmp, encoding='utf-8') as handle:
            for line in handle:
                fields = line.split('|', 2)
                self.merged[int(fields[0])] = int(fields[1])

    def __contains__(self, taxid):
        return 0 <= taxid < len(self.parent) and self.parent[taxid] >= 0

    def resolve(self, taxid):
        """ID in nodes.dmp for a taxon ID (merged IDs: the taxon they were merged into), or None."""
        if taxid in self:
            return taxid
        new = self.merged.get(taxid)
        return new if new is not None and new in self else None

    def lineage_ids(self, taxid):
        """[taxon ID or None per DESIRED_RANKS]: the deepest ancestor at each rank (taxon itself included)."""
        found = [None] * len(DESIRED_RANKS)
        node = taxid
        for _ in range(MAX_DEPTH):
            col = self.column[self.rank[node]]
            if col is not None and found[col] is None:
                found[col] = node
            if node == ROOT:
                return found
            node = self.parent[node]
        raise ValueError(f"nodes.dmp: no path from taxon {taxid} to the root")

    def scientific_names(self, taxids):
        """{taxon ID: scientific name} for the given IDs."""
        wanted = set(taxids)
        names = {}
        with open(self.names_dmp, encoding='utf-8') as handle:
            for line in handle:
                # names.dmp: tax_id | name_txt | unique name | name class |
                if 'scientific name' not in line:
                    continue
                fields = line.split('\t|\t')
                if len(fields) < 4 or fields[3].rstrip('|\t\r\n').strip() != 'scientific name':
                    continue
                taxid = int(fields[0])
                if taxid in wanted:
                    names[taxid] = fields[1].strip()
        return names


def get_lineage_as_dict(taxonomy, tax_ids):
    """DataFrame: taxid (str) and one column per DESIRED_RANKS (scientific name or None),
    for the IDs found in the taxonomy (directly or via merged.dmp)."""
    lineage_ids, merged, not_found = {}, [], []
    for tid in tax_ids:
        tid = int(tid)
        target = taxonomy.resolve(tid)
        if target is None:
            not_found.append(tid)
            continue
        if target != tid:
            merged.append(f"{tid}->{target}")
        lineage_ids[str(tid)] = taxonomy.lineage_ids(target)  # <- keep keys as strings

    names = taxonomy.scientific_names({t for ids in lineage_ids.values() for t in ids if t is not None})
    unnamed = {t for ids in lineage_ids.values() for t in ids if t is not None and t not in names}
    results = {key: [names.get(t) if t is not None else None for t in ids] for key, ids in lineage_ids.items()}

    print(f"Lineage found for {len(lineage_ids)} of {len(lineage_ids) + len(not_found)} taxon IDs "
          f"({taxonomy.n_taxa} taxa in nodes.dmp)", flush=True)
    if merged:
        print(f"{len(merged)} merged taxon IDs, lineage of the taxon they were merged into (merged.dmp): "
              + ', '.join(merged), flush=True)
    if not_found:
        print(f"{len(not_found)} taxon IDs not in nodes.dmp or merged.dmp (deleted or invalid), "
              "lineage 'missing value': " + ', '.join(map(str, sorted(not_found))), flush=True)
    if unnamed:
        print(f"WARNING: no scientific name in names.dmp for {len(unnamed)} taxa: "
              + ', '.join(map(str, sorted(unnamed))), flush=True)

    df = pd.DataFrame.from_dict(results, orient='index', columns=DESIRED_RANKS)
    df.index.name = 'taxid'
    df.reset_index(inplace=True)
    return df


def annotate_taxonomy(df: pd.DataFrame, nodes_dmp, names_dmp, merged_dmp) -> pd.DataFrame:
    if 'NCBITaxonomy' not in df.columns:
        return df

    print("annotate_taxonomy got a totoal of", len(df), "rows", flush=True)
    taxonomy = Taxonomy(nodes_dmp, names_dmp, merged_dmp)

    _df = df.copy()
    _df.dropna(subset=['NCBITaxonomy'], inplace=True)

    # Extract the FIRST numeric token before '|' and coerce safely
    first_tokens = (
        _df['NCBITaxonomy']
        .astype(str)
        .str.split('|').str[0]
        .str.extract(r'(\d+)', expand=False)
    )
    # values of several organisms: the taxon ID of each part
    multi_values = [v for v in pd.unique(_df['NCBITaxonomy'].astype(str)) if is_multi(v)]
    part_tokens = pd.Series([_part_taxid(p) for v in multi_values for p in taxonomy_parts(v)], dtype=object)
    first_tokens = pd.concat([first_tokens[~_df['NCBITaxonomy'].astype(str).isin(multi_values)], part_tokens])
    tax_ids = (
        pd.to_numeric(first_tokens, errors='coerce')
        .dropna().astype('int64').unique()
    )
    print(f"Found {len(tax_ids)} unique NCBI Taxonomy IDs", flush=True)

    lineage_df = get_lineage_as_dict(taxonomy, tax_ids)

    # Map rank → pretty column name
    rank_to_col = {
        "superkingdom": "NCBISuperkingdom",
        "kingdom": "NCBIKingdom",
        "phylum": "NCBIPhylum",
        "class": "NCBIClass",
        "order": "NCBIOrder",
        "family": "NCBIFamily",
        "genus": "NCBIGenus",
        "species": "NCBISpecies",
    }
    lineage_df = lineage_df.rename(columns=rank_to_col)
    lineage_df.fillna("missing value", inplace=True)
    lineage_df['taxid'] = lineage_df['taxid'].astype(str)

    def safe_extract_taxid(val):
        m = re.search(r'(\d+)', str(val).split('|')[0])
        return m.group(1) if m else "missing value"

    df = df.copy()
    df['taxid'] = df['NCBITaxonomy'].apply(safe_extract_taxid).astype(str)
    df['taxid'] = df['taxid'].fillna("missing value")
    # rows of several organisms: no single taxon (their lineage is set per part below)
    multi = df['NCBITaxonomy'].map(is_multi).astype(bool)
    df.loc[multi, 'taxid'] = "missing value"

    # --- drop any existing taxonomy columns before merging ---
    cols_to_drop = list(rank_to_col.values())
    df = df.drop(columns=[c for c in cols_to_drop if c in df.columns], errors="ignore")

    # merge fresh lineage data
    df = df.merge(lineage_df, on='taxid', how='left')
    df = df.fillna("missing value")

    # Drop helper taxid column
    df.drop(columns=['taxid'], inplace=True)

    # several organisms: each lineage column holds one value per part, in the order of the parts
    if multi.any():
        lineage_of = {t: row for t, row in zip(lineage_df['taxid'], lineage_df[cols_to_drop].itertuples(index=False, name=None))}
        multi = multi.to_numpy()
        for value in pd.unique(df.loc[multi, 'NCBITaxonomy']):
            per_part = [lineage_of.get(_part_taxid(p) or '') for p in taxonomy_parts(value)]
            rows = multi & (df['NCBITaxonomy'] == value).to_numpy()
            for i, column in enumerate(cols_to_drop):
                df.loc[rows, column] = aligned([lin[i] if lin is not None else None for lin in per_part])
        print(f"{int(multi.sum())} rows of several organisms (NCBITaxonomy parts joined by '{MULTI_SEP}'): "
              f"lineage per organism", flush=True)
    return df


def _part_taxid(part):
    """The taxon ID of one 'ID|Name' part (the first number before its first '|'), or None."""
    m = re.search(r'(\d+)', str(part).split('|')[0])
    return m.group(1) if m else None


def update_sampletype(df):

    # Example df
    # df = pd.DataFrame({"NCBIKingdom": [...], "SampleType": [...]})

    # Define a mapping for NCBIKingdom → SampleType
    kingdom_to_type = {
        "Viridiplantae": "plant",
        "Metazoa": "animal",
        "Fungi": "fungi",
        "Bacillati": "bacteria",
        "Pseudomonadati": "bacteria",
        "Fusobacteriati": "bacteria",
        "Thermotogati": "bacteria",
        # archaeal kingdoms
        "Thermoproteati": "archaea",
        "Methanobacteriati": "archaea",
        "Nanobdellati": "archaea",
        "Promethearchaeati": "archaea"
    }

    # Update SampleType only if it's "missing value"
    mask = df["SampleType"] == "missing value"
    # (aligned values of several organisms: one type only when all parts give the same)
    kingdoms = {k: _kingdom_type(k, kingdom_to_type) for k in pd.unique(df.loc[mask, "NCBIKingdom"]) if isinstance(k, str)}
    df.loc[mask, "SampleType"] = df.loc[mask, "NCBIKingdom"].map(kingdoms).fillna("missing value")

    return df


def _kingdom_type(kingdom, kingdom_to_type):
    """SampleType of an NCBIKingdom value; of an aligned value (several organisms) only when every part gives the
    same one, else None."""
    if not isinstance(kingdom, str) or MULTI_SEP not in kingdom:
        return kingdom_to_type.get(kingdom)
    types = {kingdom_to_type.get(p) for p in kingdom.split(MULTI_SEP)}
    return types.pop() if len(types) == 1 else None


def write_tsv_atomic(df, path):
    """Writes df to path via a temporary file next to it, which replaces path only after it was
    read back with the same header and number of rows."""
    fd, tmp = tempfile.mkstemp(dir=os.path.dirname(os.path.abspath(path)), prefix='.' + os.path.basename(path) + '.')
    os.close(fd)
    try:
        df.to_csv(tmp, index=False, sep='\t')
        with open(tmp, newline='', encoding='utf-8') as handle:
            header = next(csv.reader(handle, delimiter='\t'))
        n_rows = len(pd.read_csv(tmp, sep='\t', usecols=[0], dtype=str, na_filter=False, low_memory=False))
        if header != [str(c) for c in df.columns] or n_rows != len(df):
            raise IOError(f"{path}: written file read back with {len(header)} columns and {n_rows} rows, "
                          f"expected {len(df.columns)} and {len(df)}")
        os.chmod(tmp, 0o644)
        os.replace(tmp, path)
    except BaseException:
        if os.path.exists(tmp):
            os.remove(tmp)
        raise
    print(f"Wrote {path}: {n_rows} rows, {len(header)} columns", flush=True)


def multi_ott(df, dt_oto_ncbi):
    """OpenTreeOfLifeTaxonomyID of the rows of several organisms: the OTT id(s) of each part (dt_oto_ncbi: NCBI,
    OpenTreeOfLifeTaxonomyID; one row per NCBI id), aligned with the parts."""
    multi = df['NCBITaxonomy'].map(is_multi).astype(bool)
    if not multi.any():
        return df
    ott = dict(zip(dt_oto_ncbi['NCBI'], dt_oto_ncbi['OpenTreeOfLifeTaxonomyID']))
    values = {v: aligned([ott.get(t) for t in taxon_ids(v)]) for v in pd.unique(df.loc[multi, 'NCBITaxonomy'])}
    df.loc[multi, 'OpenTreeOfLifeTaxonomyID'] = df.loc[multi, 'NCBITaxonomy'].map(values)
    return df


def _keep_previous_ott(df, previous_ott):
    """Without a Wikidata table this run, keeps the OTT ids the rows had (row order is unchanged)."""
    if previous_ott is not None and len(previous_ott) == len(df):
        df['OpenTreeOfLifeTaxonomyID'] = previous_ott.fillna('missing value').to_numpy()
    return df


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Annotate NCBITaxonomy with its NCBI lineage, from the pinned NCBI taxonomy dump.")
    parser.add_argument("input_tsv", type=str)
    parser.add_argument("output_tsv", type=str)
    parser.add_argument("--path_to_ncbi_nodes_dmp", required=True, help="nodes.dmp of the pinned NCBI taxdump")
    parser.add_argument("--path_to_ncbi_names_dmp", required=True, help="names.dmp of the pinned NCBI taxdump")
    parser.add_argument("--path_to_ncbi_merged_dmp", required=True, help="merged.dmp of the pinned NCBI taxdump")
    args = parser.parse_args()

    # all columns as text: numeric parsing turns e.g. 66 into 66.0 for columns with gaps
    df = pd.read_csv(args.input_tsv, sep='\t', dtype=str, low_memory=False)
    # earlier runs' tables already carry the OTT column (saveOlderData keeps all old columns);
    # merging over it would create _x/_y columns and later fail, so it is rebuilt every run
    df = df.drop(columns=[c for c in df.columns if c.startswith('OpenTreeOfLifeTaxonomyID_')])
    previous_ott = df.pop('OpenTreeOfLifeTaxonomyID') if 'OpenTreeOfLifeTaxonomyID' in df.columns else None
    annotated_df = annotate_taxonomy(df, args.path_to_ncbi_nodes_dmp, args.path_to_ncbi_names_dmp,
                                     args.path_to_ncbi_merged_dmp)
    annotated_df = update_sampletype(annotated_df)

    try:
        # read the sparql_oto_ncbi tsv file
        dt_oto_ncbi = pd.read_csv('sparql_oto_ncbi.tsv', sep='\t', dtype=str)
        dt_oto_ncbi = dt_oto_ncbi.rename(columns={'?ncbi_id': 'NCBI', '?ott_id': 'uid'})
        dt_oto_ncbi = dt_oto_ncbi[(dt_oto_ncbi['NCBI'].notna()) & (dt_oto_ncbi['NCBI'] != '') & (dt_oto_ncbi['uid'].notna()) & (dt_oto_ncbi['uid'] != '')]
        # OTT ids are numbers; Wikidata also holds free-text values such as 'mmm nji'
        dt_oto_ncbi = dt_oto_ncbi[dt_oto_ncbi['uid'].str.fullmatch(r'\d+') & dt_oto_ncbi['NCBI'].str.fullmatch(r'\d+')]

        dt_oto_ncbi['OpenTreeOfLifeTaxonomyID'] = 'ott' + dt_oto_ncbi['uid']
        dt_oto_ncbi = dt_oto_ncbi[['NCBI', 'OpenTreeOfLifeTaxonomyID']].drop_duplicates()
        # some NCBI ids have several OTT ids in Wikidata: one row per NCBI id (ids joined with ';'),
        # otherwise the merge below duplicates the samples of those taxa
        dt_oto_ncbi = (dt_oto_ncbi.groupby('NCBI', sort=False)['OpenTreeOfLifeTaxonomyID']
                       .agg(lambda ids: ';'.join(sorted(set(ids)))).reset_index())

        annotated_df['NCBI'] = annotated_df['NCBITaxonomy'].apply(lambda x: x.split('|')[0] if pd.notna(x) and '|' in x and not is_multi(x) else None)

        annotated_df = annotated_df.merge(dt_oto_ncbi, on='NCBI', how='left')
        annotated_df.drop(columns=['NCBI'], inplace=True)

        annotated_df['OpenTreeOfLifeTaxonomyID'] = annotated_df['OpenTreeOfLifeTaxonomyID'].fillna('missing value')
        # several organisms: the OTT id(s) of each part, aligned with the parts
        annotated_df = multi_ott(annotated_df, dt_oto_ncbi)

    except FileNotFoundError:
        print("sparql_oto_ncbi.tsv not found, skipping OpenTreeOfLifeTaxonomyID annotation.", flush=True)
        annotated_df = _keep_previous_ott(annotated_df, previous_ott)
    except (KeyError, pd.errors.EmptyDataError, pd.errors.ParserError) as e:
        # e.g. the SPARQL endpoint returned an error page instead of the TSV
        print(f"WARNING: sparql_oto_ncbi.tsv is not a valid NCBI/OTT table ({type(e).__name__}: {e}), "
              "skipping OpenTreeOfLifeTaxonomyID annotation.", flush=True)
        annotated_df = annotated_df.drop(columns=['NCBI'], errors='ignore')
        annotated_df = _keep_previous_ott(annotated_df, previous_ott)



    write_tsv_atomic(annotated_df, args.output_tsv)
