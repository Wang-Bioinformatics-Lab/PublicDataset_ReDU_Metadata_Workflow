#!/usr/bin/env nextflow
nextflow.enable.dsl=2

TOOL_FOLDER = "$baseDir/bin"
DATA_FOLDER = "$baseDir/data"


//This is default empty file, we will overwrite this in web application
params.old_redu = './data/empty_redu.tsv'


process updateAllowedTerms {
    publishDir "./nf_output", mode: 'copy'

    conda "$TOOL_FOLDER/conda_env.yml"

    input:
    val x

    output:
    file 'allowed_terms.json'

    """
    python $TOOL_FOLDER/update_allowed_terms_from_ontologies.py  \
    $DATA_FOLDER/allowed_terms.json \
    --path_to_ncbi_dump $DATA_FOLDER/names.dmp \
    --path_to_uberon_owl $DATA_FOLDER/uberon.owl \
    --path_to_po_owl $DATA_FOLDER/po.owl \
    --path_to_cl_owl $DATA_FOLDER/cl.owl \
    --path_to_doid_owl $DATA_FOLDER/doid.owl \
    --path_to_ms_owl $DATA_FOLDER/ms.owl \
    --path_to_material_envs_owl $DATA_FOLDER/material-hierarchy.owl \
    --path_to_biome_envs_owl $DATA_FOLDER/biome-hierarchy.owl
    """
}



process downloadMetadata_massive_and_github {
    publishDir "./nf_output", mode: 'copy'

    conda "$TOOL_FOLDER/conda_env.yml"

    input:
    val x

    output:
    file 'file_paths.tsv'
    file 'metadata_folder'

    """
    mkdir metadata_folder
    wget https://github.com/Wang-Bioinformatics-Lab/ReDU_metadata/archive/refs/heads/main.zip
    unzip main
    mv ReDU_metadata-main/metadata/*.tsv metadata_folder/
    python $TOOL_FOLDER/gnps_downloader.py metadata_folder
    """
}


process mwbRun {
    conda "$TOOL_FOLDER/conda_env.yml"

    publishDir "./nf_output", mode: 'copy'

    input:
    path uberon_po_cl_csv_path
    path ENVO_bio_csv
    path ENVO_material_csv
    path ncbi_rank_division
    path allowed_terms

    output:
    file 'REDU_from_MWB_all.tsv'

    """
    python $TOOL_FOLDER/MWB_to_REDU.py \
    --study_id ALL \
    --path_to_csvs $TOOL_FOLDER/translation_sheets \
    --path_to_allowed_term_json ${allowed_terms} \
    --duplicate_raw_file_handling keep_all \
    --path_to_uberon_cl_po_csv ${uberon_po_cl_csv_path} \
    --path_to_envo_biome_csv ${ENVO_bio_csv} \
    --path_to_envo_material_csv ${ENVO_material_csv} \
    --path_ncbi_rank_division ${ncbi_rank_division} \
    --path_to_polarity_info $DATA_FOLDER/MWB_polarity_table.csv
    """
}


process mwbFiles {
    publishDir "./nf_output", mode: 'copy'

    conda "$TOOL_FOLDER/conda_env.yml"

    input:
    val x

    output:
    file 'mwb_files_all.tsv'

    """
    python $TOOL_FOLDER/getAllWorkbench_file_paths.py \
    --study_id ALL \
    --output_path mwb_files_all.tsv \
    --filter_extensions True
    """
}

process formatmwb {
    publishDir "./nf_output", mode: 'copy'

    conda "$TOOL_FOLDER/conda_env.yml"

    input:
    file mwb_metadata
    file mwb_files

    output:
    file 'mwb_redu.tsv'

    """
    python $TOOL_FOLDER/MWB_merge.py  \
    $mwb_metadata \
    $mwb_files \
    mwb_redu.tsv 
    """
}


process mlRun {

    conda "$TOOL_FOLDER/conda_env.yml"

    publishDir "./nf_output", mode: 'copy'

    input:
    path uberon_po_cl_csv_path
    path ENVO_bio_csv
    path ENVO_material_csv
    path ncbi_rank_division
    path allowed_terms

    output:
    file 'Metabolights2REDU_ALL.tsv'


    """
    python $TOOL_FOLDER/Metabolights2REDU.py \
    --study_id ALL  \
    --path_to_translation_sheet_csvs $TOOL_FOLDER/translation_sheets_metabolights \
    --path_to_allowed_term_json ${allowed_terms} \
    --path_to_uberon_cl_po_csv ${uberon_po_cl_csv_path} \
    --path_to_envo_biome_csv ${ENVO_bio_csv} \
    --path_to_envo_material_csv ${ENVO_material_csv} \
    --path_ncbi_rank_division ${ncbi_rank_division}
    """
}

process mlFiles {
    publishDir "./nf_output", mode: 'copy'

    conda "$TOOL_FOLDER/conda_env.yml"

    input:
    val x

    output:
    file 'MetabolightsFilePaths_ALL.tsv'

    """
    python $TOOL_FOLDER/getAllMetabolights_file_paths.py \
    --output_filename MetabolightsFilePaths_ALL.tsv \
    --user_token e6db13e8-bfa7-452c-83dd-92ddf10677c1
    """
}

process formatml {
    publishDir "./nf_output", mode: 'copy'

    conda "$TOOL_FOLDER/conda_env.yml"


    input:
    file ml_metadata
    file ml_files

    output:
    file 'ml_redu.tsv'

    """
    python $TOOL_FOLDER/ML_merge.py \
    $ml_metadata \
    $ml_files \
    ml_redu.tsv 
    """
}

process normanRun {

    conda "$TOOL_FOLDER/conda_env.yml"

    publishDir "./nf_output", mode: 'copy'

    input:
    path uberon_po_cl_csv_path
    path ENVO_bio_csv
    path ENVO_material_csv
    path ncbi_rank_division
    path allowed_terms

    output:
    file 'NORMAN2REDU_ALL.tsv'


    """
    python $TOOL_FOLDER/NORMAN_to_REDU.py \
    --study_id ALL  \
    --path_to_allowed_term_json ${allowed_terms} \
    --path_to_uberon_cl_po_csv ${uberon_po_cl_csv_path} \
    --path_to_envo_biome_csv ${ENVO_bio_csv} \
    --path_to_envo_material_csv ${ENVO_material_csv} \
    --path_ncbi_rank_division ${ncbi_rank_division} \
    --output NORMAN2REDU_ALL.tsv
    """
}

process mergeAllMetadata {
    publishDir "./nf_output", mode: 'copy'

    conda "$TOOL_FOLDER/conda_env.yml"

    input:
    file gnps_metadata
    file mwb_redu
    file metabolights_redu
    file norman_redu
    file masst_metadata

    output:
    path 'merged_metadata_wo_cache_columns.tsv'

    """
    python $TOOL_FOLDER/merge_metadata.py \
    --gnps_metadata $gnps_metadata \
    --mwb_metadata $mwb_redu \
    --metabolights_metadata $metabolights_redu \
    --norman_metadata $norman_redu \
    --masst_metadata $masst_metadata \
    --output_metadata merged_metadata_wo_cache_columns.tsv
    """
}

process prepare_ontologies {
    publishDir "./nf_output", mode: 'copy'

    conda "$TOOL_FOLDER/conda_env.yml"

    input:
    val x

    output:
    path 'UBERON_CL_PO_ontology.csv'
    path 'DOID_ontology.csv'
    path 'ENVO_biome_ontology.csv'
    path 'ENVO_material_ontology.csv'
    path 'NCBI_Rank_Division.csv'

    """
    python $TOOL_FOLDER/prepare_ontologies.py \
    --path_to_uberon_owl $DATA_FOLDER/uberon.owl \
    --path_to_cl_owl $DATA_FOLDER/cl.owl \
    --path_to_po_owl $DATA_FOLDER/po.owl \
    --path_to_doid_owl $DATA_FOLDER/doid.owl \
    --path_to_biome_envs_owl $DATA_FOLDER/biome-hierarchy.owl \
    --path_to_material_envs_owl $DATA_FOLDER/material-hierarchy.owl \
    --path_to_ncbi_nodes_dmp $DATA_FOLDER/nodes.dmp \
    --path_to_ncbi_division_dmp $DATA_FOLDER/division.dmp
    """
}


// This cleans up the metadata from MassIVE into the appropriate CV terms
process gnpsHarmonize {
    publishDir "./nf_output", mode: 'copy'

    conda "$TOOL_FOLDER/conda_env.yml"

    input:
    path metadata_ch
    path UBERON_CL_PO_ontology_csv
    path DOID_ontology_csv
    path ENVO_bio_csv
    path ENVO_material_csv
    path ncbi_rank_division
    path allowed_terms

    output:
    file 'adjusted_metadata_folder'

    """
    mkdir adjusted_metadata_folder
    python $TOOL_FOLDER/read_and_validate_redu_from_github.py \
    ${metadata_ch} \
    adjusted_metadata_folder \
    --AllowedTermJson_path ${allowed_terms} \
    --path_ncbi_rank_division ${ncbi_rank_division} \
    --path_to_uberon_cl_po_csv ${UBERON_CL_PO_ontology_csv} \
    --path_to_doid_csv ${DOID_ontology_csv} \
    --path_to_envo_biome_csv ${ENVO_bio_csv} \
    --path_to_envo_material_csv ${ENVO_material_csv}
    """
}


process gnpsmatchName {
    publishDir "./nf_output", mode: 'copy'

    conda "$TOOL_FOLDER/conda_env.yml"


    input:
    // file 'passed_file_names.tsv'
    file 'adjusted_metadata_folder'
    path allowed_terms

    output:
    file 'gnps_metadata_all.tsv'

    """
    python $TOOL_FOLDER/gnps_name_matcher.py \
    all \
    adjusted_metadata_folder \
    gnps_metadata_all.tsv \
    ${allowed_terms}
    """
}


// Getting microbmemasst data and putting it in a tentaive redu format
process MASST_to_REDU {
    publishDir "./nf_output", mode: 'copy'

    conda "$TOOL_FOLDER/conda_env.yml"

    input:
    path redu_table
    path ncbi_rank_division
    path allowed_terms

    output:
    file 'adjusted_metadata_folder'

    """
    mkdir adjusted_metadata_folder
    python $TOOL_FOLDER/MASST_to_REDU.py \
    ${redu_table} \
    adjusted_metadata_folder \
    --path_microbeMASST $DATA_FOLDER/microbe_masst_table.csv \
    --path_plantMASST $DATA_FOLDER/plant_masst_table.csv \
    --path_ncbiRanksDivisions ${ncbi_rank_division} \
    --AllowedTermJson_path ${allowed_terms}
    """
}

// This is taking metadata from microbeMASST and PlatnMASST and modifying it
process gnpsmatchName_masst {
    publishDir "./nf_output", mode: 'copy'

    conda "$TOOL_FOLDER/conda_env.yml"


    input:
    file 'adjusted_metadata_folder' 
    path allowed_terms

    output:
    file 'masst_metadata_all.tsv'

    """
    python $TOOL_FOLDER/gnps_name_matcher.py \
    all \
    adjusted_metadata_folder \
    masst_metadata_all.tsv \
    ${allowed_terms}
    """
}

process add_cache_columns {
    publishDir "./nf_output", mode: 'copy'

    cache false

    conda "$TOOL_FOLDER/conda_env.yml"

    input:
    file merged_ch

    output:
    file 'merged_metadata.tsv'

    """
    python $TOOL_FOLDER/add_cache_columns.py \
    ${merged_ch} \
    merged_metadata.tsv
    """
}


process saveOlderData {
    publishDir "./nf_output", mode: 'copy'

    conda "$TOOL_FOLDER/conda_env.yml"

    input:
    path merged_ch
    path old_redu

    output:
    path 'merged_with_old.tsv'

    """
    python $TOOL_FOLDER/save_older_data.py \
    ${merged_ch} \
    ${old_redu} \
    merged_with_old.tsv
    """
}


workflow {

    //  Prepare ontologies
    (uberon_cl_co_onto, doid_onto, envo_bio, envo_material, ncbi_rank_division) = prepare_ontologies(1)
    allowed_terms = updateAllowedTerms(1)

    // Massive REDU data, called before GitHub because taking it from MassIVE as the place to keep metadata and not github
    (file_paths_ch, metadata_ch) = downloadMetadata_massive_and_github(1)
    msv_metadata_ch = gnpsHarmonize(metadata_ch, uberon_cl_co_onto, doid_onto, envo_bio, envo_material, ncbi_rank_division, allowed_terms)
    gnps_metadata_ch = gnpsmatchName(msv_metadata_ch, allowed_terms)

    // MicrobeMASST and PlantMASST
    masst_metadata_ch = MASST_to_REDU(gnps_metadata_ch, ncbi_rank_division, allowed_terms)
    masst_metadata_wFiles_ch = gnpsmatchName_masst(masst_metadata_ch, allowed_terms)

    // Metabolomics Workbench
    mwb_metadata_ch = mwbRun(uberon_cl_co_onto, envo_bio, envo_material, ncbi_rank_division, allowed_terms)
    mwb_files_ch = mwbFiles(1)
    mwb_redu_ch = formatmwb(mwb_metadata_ch, mwb_files_ch)

    // Metabolights
    ml_metadata_ch = mlRun(uberon_cl_co_onto, envo_bio, envo_material, ncbi_rank_division, allowed_terms)
    ml_files_ch = mlFiles(1)
    ml_redu_ch = formatml(ml_metadata_ch, ml_files_ch)

    // NORMAN
    norman_metadata_ch = normanRun(uberon_cl_co_onto, envo_bio, envo_material, ncbi_rank_division, allowed_terms)

    // Combine everything
    merged_ch = mergeAllMetadata(gnps_metadata_ch, mwb_redu_ch, ml_redu_ch, norman_metadata_ch, masst_metadata_wFiles_ch)


    // Make sure we dont loose older data
    old_redu_path_ch = Channel.fromPath(params.old_redu)
    merged_with_old_ch = saveOlderData(merged_ch, old_redu_path_ch)

    // Add Cache columns
    cache_enriched_metadata = add_cache_columns(merged_with_old_ch)

}
