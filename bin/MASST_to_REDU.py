import argparse
import pandas as pd
import json
import os
from read_and_validate_redu_from_github import complete_and_fill_REDU_table

# NCBI taxids used to tell microbes apart (Bacteria division also holds archaea,
# 'Plants and Fungi' holds both)
TAXID_ARCHAEA = '2157'
TAXID_FUNGI = '4751'


def load_ncbi_parents(path_to_nodes_dmp):
    """taxid -> parent taxid from NCBI nodes.dmp"""
    parents = {}
    with open(path_to_nodes_dmp) as f:
        for line in f:
            taxid, parent = line.split('\t|\t', 2)[:2]
            parents[taxid] = parent
    return parents


def has_ancestor(taxid, ancestor, parents):
    seen = set()
    while taxid and taxid not in seen:
        if taxid == ancestor:
            return True
        seen.add(taxid)
        parent = parents.get(taxid)
        if parent == taxid:
            break
        taxid = parent
    return False


def microbe_sample_type(taxid, division_by_taxid, parents):
    """SampleType/SampleTypeSub1 of a microbeMASST culture, or ('', '') if unclear"""
    division = division_by_taxid.get(taxid)
    if division == 'Bacteria' and parents is not None:
        if has_ancestor(taxid, TAXID_ARCHAEA, parents):
            return 'archaea', ''
        return 'culture_bacterial', 'culture_bacterial'
    if division == 'Plants and Fungi' and parents is not None and has_ancestor(taxid, TAXID_FUNGI, parents):
        return 'culture_fungal', 'culture_fungal'
    # microbeMASST also holds mammalian cell cultures (e.g. human cell lines)
    if division in ('Primates', 'Rodents', 'Mammals'):
        return 'culture_mammalian', 'culture_mammalian'
    return '', ''


def apply_blank_qc_flags(df):
    """Blank/QC files carry no organism: label them as blank/QC and drop the taxonomy"""
    blank = df['Blank'].astype(str).str.strip().str.lower().eq('yes')
    qc = df['QC'].astype(str).str.strip().str.lower().eq('yes') & ~blank
    df.loc[blank, 'SampleType'] = 'blank_analysis'
    df.loc[blank, 'SampleTypeSub1'] = 'blank_analysis'
    df.loc[qc, 'SampleType'] = 'pool_QC'
    df.loc[qc, 'SampleTypeSub1'] = 'pool_QC'
    df.loc[blank | qc, 'NCBITaxonomy'] = ''
    print(f"Blank/QC files: {blank.sum()} blanks, {qc.sum()} QCs")
    return df





if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='GNPS Validator')
    parser.add_argument('path_REDUtable')
    parser.add_argument('output_folder')
    parser.add_argument('--path_microbeMASST')
    parser.add_argument('--path_plantMASST')
    parser.add_argument('--path_ncbiRanksDivisions')
    parser.add_argument('--path_ncbi_nodes', default=None, help='NCBI nodes.dmp, to tell bacteria/archaea and fungi/plants apart')
    parser.add_argument("--AllowedTermJson_path", type=str, help="Path to json with allowed terms")
    args = parser.parse_args()



    df_redu = pd.read_csv(args.path_REDUtable, sep = '\t', dtype=str)
    df_redu['filename'] = df_redu['filename'].str.slice(2)
    # Files already described by GNPS metadata, identified by dataset + path (the same
    # relative path, e.g. ccms_peak/Blank.mzML, exists in many datasets)
    redu_files = set(zip(df_redu['ATTRIBUTE_DatasetAccession'], df_redu['filename']))

    def _masst_files_not_in_redu(df):
        df = df.dropna(subset=['Filepath', 'MassIVE']).copy()
        df['filename'] = df['Filepath'].apply(lambda x: '/'.join(x.split('/')[1:]))
        return df[[(m, f) not in redu_files for m, f in zip(df['MassIVE'], df['filename'])]]

    df_mm = pd.read_csv(args.path_microbeMASST, dtype=str)
    print(df_mm)
    df_mm = _masst_files_not_in_redu(df_mm)

    df_pm = pd.read_csv(args.path_plantMASST, dtype=str)
    print(df_pm)
    df_pm = _masst_files_not_in_redu(df_pm)


    with open(args.AllowedTermJson_path, 'r') as json_file:
        allowed_terms = json.load(json_file)


    NCBIRankDivision_table = pd.read_csv(args.path_ncbiRanksDivisions, index_col = False)
    NCBIRankDivision_table = NCBIRankDivision_table.drop_duplicates(subset=['TaxonID'])

    

    df_list = []

    if len(df_mm) > 0:
        df_mm['NCBITaxonomy'] = df_mm['Taxa_NCBI'] + '|' + df_mm['Taxaname_file']
        df_mm = df_mm.rename(columns={'MassIVE': 'MassiveID'})
        parents = load_ncbi_parents(args.path_ncbi_nodes) if args.path_ncbi_nodes else None
        division_by_taxid = dict(zip(NCBIRankDivision_table['TaxonID'].astype(str), NCBIRankDivision_table['NCBIDivision']))
        types = [microbe_sample_type(str(t).strip(), division_by_taxid, parents) for t in df_mm['Taxa_NCBI']]
        df_mm['SampleType'] = [t[0] for t in types]
        df_mm['SampleTypeSub1'] = [t[1] for t in types]
        df_mm = apply_blank_qc_flags(df_mm)
        df_mm = df_mm[['MassiveID', 'filename', 'SampleType', 'SampleTypeSub1', 'NCBITaxonomy']]
        df_list.append(df_mm)


    if len(df_pm) > 0:
        df_pm['NCBITaxonomy'] = df_pm['Taxa_NCBI'] + '|' + df_pm['Taxaname_file']
        df_pm = df_pm.rename(columns={'MassIVE': 'MassiveID'})
        df_pm['SampleType'] = 'plant'
        df_pm['SampleTypeSub1'] = ''
        df_pm = apply_blank_qc_flags(df_pm)
        df_pm = df_pm[['MassiveID', 'filename', 'SampleType', 'SampleTypeSub1', 'NCBITaxonomy']]
        df_list.append(df_pm)


    if len(df_list) > 0:

        df_massts = pd.concat(df_list, ignore_index=True)
        df_massts_filled = complete_and_fill_REDU_table(df_massts, allowedTerm_dict=allowed_terms, NCBIRankDivision_table=NCBIRankDivision_table)

        #save output to csv
        for massive_id in df_massts_filled['MassiveID'].unique():
            # Filter the DataFrame for the current MassiveID
            df_filtered = df_massts_filled[df_massts_filled['MassiveID'] == massive_id]
            
            # Construct the file name
            file_name = f"{args.output_folder}/{massive_id}.tsv"
            
            # Save the filtered DataFrame to a TSV file
            df_filtered.to_csv(file_name, sep='\t', index=False)


