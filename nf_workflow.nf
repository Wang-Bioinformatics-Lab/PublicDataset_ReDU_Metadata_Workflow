#!/usr/bin/env nextflow
nextflow.enable.dsl=2

TOOL_FOLDER = "$baseDir/bin"

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

    output:
    file 'passed_file_names.tsv'

    """
    python $TOOL_FOLDER/gnps_validator.py \
    file_paths.tsv \
    metadata_folder \
    --AllowedTermJson_path $TOOL_FOLDER/allowed_terms/allowed_terms.json
    """
}

process gnpsmatchName {
    publishDir "./nf_output", mode: 'copy'

    conda "$TOOL_FOLDER/conda_env.yml"

    cache false

    input:
    file 'passed_file_names.tsv'
    file 'metadata_folder' 

    output:
    file 'gnps_metadata_all.tsv'

    """
    python $TOOL_FOLDER/gnps_name_matcher.py \
    passed_file_names.tsv \
    metadata_folder
    """
}

process mwbRun {
    conda "$TOOL_FOLDER/conda_env.yml"

    publishDir "./nf_output", mode: 'copy'

    input:
    val x

    output:
    file 'REDU_from_MWB_all.tsv'

    """
    echo "Conda env info:"
    conda info
    echo "Python version:"
    python --version
    echo "Which Python:"
    which python
    echo "Pip list:"
    pip list
    python $TOOL_FOLDER/MWB_to_REDU.py \
    --study_id ALL \
    --path_to_csvs $TOOL_FOLDER/translation_sheets \
    --path_to_allowed_term_json $TOOL_FOLDER/allowed_terms/allowed_terms.json \
    --duplicate_raw_file_handling keep_all \
    --path_to_uberon_owl $TOOL_FOLDER/allowed_terms/uberon-base.owl \
    --path_to_cl_owl $TOOL_FOLDER/allowed_terms/cl.owl \
    --path_to_plant_owl $TOOL_FOLDER/allowed_terms/po.owl
    """

    // """
    // python $TOOL_FOLDER/MWB_to_REDU.py \
    // --study_id ALL \
    // --path_to_csvs $TOOL_FOLDER/translation_sheets \
    // --path_to_allowed_term_json $TOOL_FOLDER/allowed_terms/allowed_terms.json \
    // --duplicate_raw_file_handling keep_pols_dupl
    // """
}

process validateMetadata_MWB {
    publishDir "./nf_output", mode: 'copy'

    conda "$TOOL_FOLDER/conda_env.yml"

    input:
    file 'file_paths.tsv'
    file 'metadata_folder' 

    output:
    file 'passed_file_names.tsv'

    """
    python $TOOL_FOLDER/gnps_validator.py \
    file_paths.tsv \
    metadata_folder \
    --AllowedTermJson_path $TOOL_FOLDER/allowed_terms/allowed_terms.json
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
    val x

    output:
    file 'Metabolights2REDU_ALL.tsv'

    """
    python $TOOL_FOLDER/Metabolights2REDU.py \
    --study_id ALL  \
    --path_to_translation_sheet_csvs $TOOL_FOLDER/translation_sheets_metabolights \
    --path_to_allowed_term_json $TOOL_FOLDER/allowed_terms/allowed_terms.json \
    --path_to_uberon_owl $TOOL_FOLDER/allowed_terms/uberon-base.owl \
    --path_to_cl_owl $TOOL_FOLDER/allowed_terms/cl.owl \
    --path_to_plant_owl $TOOL_FOLDER/allowed_terms/po.owl
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


workflow {
    (file_paths_ch, metadata_ch) = downloadMetadata(1)
    (passed_paths_ch) = validateMetadata(file_paths_ch, metadata_ch)
    gnps_metadata_ch = gnpsmatchName(passed_paths_ch, metadata_ch)
    
    // Metabolomics Workbench
    mwb_metadata_ch = mwbRun(1)
    mwb_files_ch = mwbFiles(1)
    mwb_redu_ch = formatmwb(mwb_metadata_ch, mwb_files_ch)

    // Metabolights
    ml_metadata_ch = mlRun(1)
    ml_files_ch = mlFiles(1)
    ml_redu_ch = formatml(ml_metadata_ch, ml_files_ch)

    merged_ch = mergeAllMetadata(gnps_metadata_ch, mwb_redu_ch, ml_redu_ch)
}
