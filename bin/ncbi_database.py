import ete3
import pandas as pd
from tqdm import tqdm
import re


### Credit to Michael Strobel (adapted)

def get_lineage_as_dict(ncbi, tax_ids, pbar=True):
    desired_ranks = ['superkingdom','kingdom','phylum','class','order','family','genus','species']
    results = {}

    it = tqdm(tax_ids, desc="Fetching lineage for tax_ids") if pbar else tax_ids
    for tid in it:
        try:
            lineage = ncbi.get_lineage(tid)
            names = ncbi.get_taxid_translator(lineage)
            ranks = ncbi.get_rank(lineage)
            lineage_by_rank = {ranks[t]: names[t] for t in lineage if ranks[t] in desired_ranks}
            row = {rank: lineage_by_rank.get(rank) for rank in desired_ranks}
            results[str(tid)] = row  # <- keep keys as strings
        except Exception as e:
            print(f"Error fetching lineage for taxid {tid}: {e}")

    df = pd.DataFrame.from_dict(results, orient='index')
    df.index.name = 'taxid'
    df.reset_index(inplace=True)
    return df


def annotate_taxonomy(df: pd.DataFrame) -> pd.DataFrame:
    if 'NCBITaxonomy' not in df.columns:
        return df

    print("annotate_taxonomy got a totoal of", len(df), "rows", flush=True)
    ncbi = ete3.NCBITaxa()

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

    lineage_df = get_lineage_as_dict(ncbi, tax_ids, pbar=True)

    # Map ete3 rank → pretty column name
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

    df = df.merge(lineage_df, on='taxid', how='left')
    df = df.fillna("missing value")

    # Drop helper taxid column
    df.drop(columns=['taxid'], inplace=True)
    return df




if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Annotate NCBITaxonomy with taxonomic classifications using ete3.")
    parser.add_argument("input_tsv", type=str)
    parser.add_argument("output_tsv", type=str)
    args = parser.parse_args()

    df = pd.read_csv(args.input_tsv, sep='\t', dtype={'NCBITaxonomy': 'string'}, low_memory=False)
    annotated_df = annotate_taxonomy(df)
    annotated_df.to_csv(args.output_tsv, index=False, sep='\t')
