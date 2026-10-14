"""Adds the NCBI lineage of each row's NCBITaxonomy (NCBISuperkingdom, NCBIKingdom, NCBIPhylum, NCBIClass,
NCBIOrder, NCBIFamily, NCBIGenus, NCBISpecies), fills SampleType from NCBIKingdom where it is 'missing value',
and adds OpenTreeOfLifeTaxonomyID (Wikidata NCBI -> OTT table sparql_oto_ncbi.tsv in the working directory).

The lineage comes from the NCBI taxonomy pinned in data/reference_data.json (nodes.dmp, names.dmp and
merged.dmp of the taxdmp archive, set up by bin/fetch_reference_data.py), the release every other step uses.
Earlier versions used ete3.NCBITaxa(), i.e. whichever taxonomy copy ete3 had in ~/.etetoolkit (or the newest
NCBI release, downloaded on first use); for the same taxdump the columns are the same.

- Taxon ID of a row: the first number in the first '|' field of NCBITaxonomy ('562|Escherichia coli' -> 562).
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

    # --- drop any existing taxonomy columns before merging ---
    cols_to_drop = list(rank_to_col.values())
    df = df.drop(columns=[c for c in cols_to_drop if c in df.columns], errors="ignore")

    # merge fresh lineage data
    df = df.merge(lineage_df, on='taxid', how='left')
    df = df.fillna("missing value")

    # Drop helper taxid column
    df.drop(columns=['taxid'], inplace=True)
    return df


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
    df.loc[mask, "SampleType"] = df.loc[mask, "NCBIKingdom"].map(kingdom_to_type).fillna("missing value")

    return df


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

        annotated_df['NCBI'] = annotated_df['NCBITaxonomy'].apply(lambda x: x.split('|')[0] if pd.notna(x) and '|' in x else None)

        annotated_df = annotated_df.merge(dt_oto_ncbi, on='NCBI', how='left')
        annotated_df.drop(columns=['NCBI'], inplace=True)

        annotated_df['OpenTreeOfLifeTaxonomyID'] = annotated_df['OpenTreeOfLifeTaxonomyID'].fillna('missing value')

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
