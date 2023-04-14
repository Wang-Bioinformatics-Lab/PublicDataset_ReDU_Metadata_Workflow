#!/usr/bin/env nextflow
nextflow.enable.dsl=2

TOOL_FOLDER = "$baseDir/bin"

process downloadMetadata {
    publishDir "./nf_output", mode: 'copy'

    conda "$TOOL_FOLDER/conda_env.yml"

    input:
    val x

    output:
    //file 'outputs.tsv.*'
    file 'file_paths.tsv'
    file 'metadata_folder'
    //file 'passed_file_names.tsv'
    //file 'check.tsv'

    """
    mkdir metadata_folder
    python $TOOL_FOLDER/gnps_downloader.py metadata_folder
    """
}

process validateMetadata {
    publishDir "./nf_output", mode: 'copy'

    cache false

    conda "$TOOL_FOLDER/conda_env.yml"

    input:
    file 'file_paths.tsv'
    file 'metadata_folder' 

    //output:
    //file 'outputs.tsv.*'
    
    //file 'passed_file_names.tsv'
    //file 'check.tsv'

    """
    python $TOOL_FOLDER/gnps_validator.py \
    file_paths.tsv \
    metadata_folder
    """
}




//python $TOOL_FOLDER/gnps_validator.py passed_file_names.tsv
//python $TOOL_FOLDER/gnps_name_matcher.py check.tsv

workflow {
    (file_paths_ch, metadata_ch) = downloadMetadata(1)
    validateMetadata(file_paths_ch, metadata_ch)
}
