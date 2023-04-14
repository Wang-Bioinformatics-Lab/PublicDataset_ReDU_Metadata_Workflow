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
    metadata_folder
    """
}

process matchName {
    publishDir "./nf_output", mode: 'copy'

    conda "$TOOL_FOLDER/conda_env.yml"

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
    publishDir "./nf_output", mode: 'copy'

    cache false

    conda "$TOOL_FOLDER/conda_env.yml"

    input:
    val x

    output:
    file 'mwb_metadata_all.tsv'

    """
    python $TOOL_FOLDER/MWB_to_REDU.py \
    --study_id ALL \
    --path_to_csvs $TOOL_FOLDER/translation_sheets
    """
}



workflow {
    (file_paths_ch, metadata_ch) = downloadMetadata(1)
    (passed_paths_ch) = validateMetadata(file_paths_ch, metadata_ch)
    matchName(passed_paths_ch, metadata_ch)

    mwbRun(1)
}
