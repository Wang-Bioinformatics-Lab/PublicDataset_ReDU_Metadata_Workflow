#!/usr/bin/env nextflow
nextflow.enable.dsl=2

TOOL_FOLDER = "$baseDir/bin"
DATA_FOLDER = "$baseDir/data"


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
    --path_to_ms_owl $DATA_FOLDER/ms.owl
    """
}



process downloadMetadata {
    publishDir "./nf_output", mode: 'copy'

    conda "$TOOL_FOLDER/conda_env.yml"

    input:
    val x

    output:
    file 'file_paths.tsv'
    file 'metadata_folder'

    """
    mkdir metadata_folder
    python $TOOL_FOLDER/gnps_downloader.py metadata_folder
    """
}

process validateMetadata {
    publishDir "./nf_output", mode: 'copy'

    conda "$TOOL_FOLDER/conda_env.yml"

    input:
    file 'file_paths.tsv'
    file 'metadata_folder'
    file 'allowed_terms.json'

    output:
    file 'passed_file_names.tsv'

    """
    python $TOOL_FOLDER/gnps_validator.py \
    file_paths.tsv \
    metadata_folder \
    --AllowedTermJson_path 'allowed_terms.json'
    """
}

process gnpsmatchName {
    publishDir "./nf_output", mode: 'copy'

    conda "$TOOL_FOLDER/conda_env.yml"

    cache false

    input:
    // file 'passed_file_names.tsv'
    val passed_file_names
    file 'metadata_folder' 

    output:
    file 'gnps_metadata_all.tsv'

    """
    python $TOOL_FOLDER/gnps_name_matcher.py \
    ${passed_file_names} \
    metadata_folder
    """

    // """
    // python $TOOL_FOLDER/gnps_name_matcher.py \
    // passed_file_names.tsv \
    // metadata_folder
    // """
}

process mwbRun {
    conda "$TOOL_FOLDER/conda_env.yml"

    publishDir "./nf_output", mode: 'copy'

    input:
    path uberon_po_cl_csv_path
    path allowed_terms

    output:
    file 'REDU_from_MWB_all.tsv'

    """
    python $TOOL_FOLDER/MWB_to_REDU.py \
    --study_id ALL \
    --path_to_csvs $TOOL_FOLDER/translation_sheets \
    --path_to_allowed_term_json ${allowed_terms} \
    --duplicate_raw_file_handling keep_all \
    --path_to_uberon_cl_po_csv ${uberon_po_cl_csv_path}
    """

    // """
    // python $TOOL_FOLDER/MWB_to_REDU.py \
    // --study_id ALL \
    // --path_to_csvs $TOOL_FOLDER/translation_sheets \
    // --path_to_allowed_term_json $TOOL_FOLDER/allowed_terms/allowed_terms.json \
    // --duplicate_raw_file_handling keep_pols_dupl
    // """
}

// process validateMetadata_MWB {
//     publishDir "./nf_output", mode: 'copy'

//     conda "$TOOL_FOLDER/conda_env.yml"

//     input:
//     file 'file_paths.tsv'
//     file 'metadata_folder' 

//     output:
//     file 'passed_file_names.tsv'

//     """
//     python $TOOL_FOLDER/gnps_validator.py \
//     file_paths.tsv \
//     metadata_folder \
//     --AllowedTermJson_path $TOOL_FOLDER/allowed_terms/allowed_terms.json
//     """
// }

process mwbFiles {
    publishDir "./nf_output", mode: 'copy'

    conda "$TOOL_FOLDER/conda_env.yml"

    input:
    val x

    output:
    file 'mwb_files_all.tsv'

    """
    python $TOOL_FOLDER/MWB_to_fileDF.py \
    --study_id ALL \
    --output_path mwb_files_all.tsv
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
    python $TOOL_FOLDER/MWB_merge.py \
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
    path allowed_terms

    output:
    file 'Metabolights2REDU_ALL.tsv'


    """
    python $TOOL_FOLDER/Metabolights2REDU.py \
    --study_id ALL  \
    --path_to_translation_sheet_csvs $TOOL_FOLDER/translation_sheets_metabolights \
    --path_to_allowed_term_json ${allowed_terms} \
    --path_to_uberon_cl_po_csv ${uberon_po_cl_csv_path}
    """

    // """
    // python $TOOL_FOLDER/Metabolights2REDU.py \
    // --study_id ALL  \
    // --path_to_translation_sheet_csvs $TOOL_FOLDER/translation_sheets_metabolights \
    // --path_to_allowed_term_json $TOOL_FOLDER/allowed_terms/allowed_terms.json \
    // --path_to_uberon_owl $TOOL_FOLDER/allowed_terms/uberon-base.owl \
    // --path_to_cl_owl $TOOL_FOLDER/allowed_terms/cl.owl \
    // --path_to_plant_owl $TOOL_FOLDER/allowed_terms/po.owl
    // """
}

process mlFiles {
    publishDir "./nf_output", mode: 'copy'

    conda "$TOOL_FOLDER/conda_env.yml"

    input:
    val x

    output:
    file 'MetabolightsFilePaths_ALL.tsv'

    """
    python $TOOL_FOLDER/GetAllMetabolightsFiles.py \
    --study_id ALL
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


process mergeAllMetadata {
    publishDir "./nf_output", mode: 'copy'

    conda "$TOOL_FOLDER/conda_env.yml"

    cache false

    input:
    file gnps_metadata
    file mwb_redu
    file metabolights_redu

    output:
    file 'merged_metadata.tsv'

    """
    python $TOOL_FOLDER/merge_metadata.py \
    $gnps_metadata \
    $mwb_redu \
    $metabolights_redu \
    merged_metadata.tsv
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

    """
    python $TOOL_FOLDER/prepare_ontologies.py \
    --path_to_uberon_owl $DATA_FOLDER/uberon.owl \
    --path_to_cl_owl $DATA_FOLDER/cl.owl \
    --path_to_po_owl $DATA_FOLDER/po.owl \
    --path_to_doid_owl $DATA_FOLDER/doid.owl
    """
}


process read_and_clean_github_redu_metadata {
    publishDir "./nf_output", mode: 'copy'

    conda "$TOOL_FOLDER/conda_env.yml"

    input:
    path UBERON_CL_PO_ontology_csv
    path DOID_ontology_csv
    path allowed_terms

    output:
    file 'metadata_folder'

    """
    mkdir metadata_folder
    python $TOOL_FOLDER/read_and_validate_redu_from_github.py \
    /home/yasin/projects/ReDU_metadata/metadata \
    metadata_folder \
    --AllowedTermJson_path $TOOL_FOLDER/allowed_terms/allowed_terms.json \
    --path_to_uberon_cl_po_csv ${UBERON_CL_PO_ontology_csv} \
    --path_to_doid_csv ${DOID_ontology_csv}
    """
}

process gnpsmatchName_github {
    publishDir "./nf_output", mode: 'copy'

    conda "$TOOL_FOLDER/conda_env.yml"

    cache false

    input:
    // file 'passed_file_names.tsv'
    val passed_file_names
    file 'metadata_folder' 

    output:
    file 'gnps_metadata_all.tsv'

    """
    python $TOOL_FOLDER/gnps_name_matcher.py \
    ${passed_file_names} \
    metadata_folder
    """
}



workflow {

    //  Prepare ontologies
    (uberon_cl_co_onto, doid_onto) = prepare_ontologies(1)
    allowed_terms = updateAllowedTerms(1)

    // Github REDU data
    // prepared_files_folder = read_and_clean_github_redu_metadata(uberon_cl_co_onto, doid_onto, allowed_terms)
    // gnps_github_metadata_ch = gnpsmatchName_github('all', prepared_files_folder)
    // need to compare with gnps_massive_metadata

    // Massive REDU data
    (file_paths_ch, metadata_ch) = downloadMetadata(1)
    (passed_paths_ch) = validateMetadata(file_paths_ch, metadata_ch, allowed_terms)
    gnps_metadata_ch = gnpsmatchName(passed_paths_ch, metadata_ch)
    
    // Metabolomics Workbench
    mwb_metadata_ch = mwbRun(uberon_cl_co_onto, allowed_terms)
    mwb_files_ch = mwbFiles(1)
    mwb_redu_ch = formatmwb(mwb_metadata_ch, mwb_files_ch)

    // Metabolights
    ml_metadata_ch = mlRun(uberon_cl_co_onto, allowed_terms)
    ml_files_ch = mlFiles(1)
    ml_redu_ch = formatml(ml_metadata_ch, ml_files_ch)

    merged_ch = mergeAllMetadata(gnps_metadata_ch, mwb_redu_ch, ml_redu_ch)
}
