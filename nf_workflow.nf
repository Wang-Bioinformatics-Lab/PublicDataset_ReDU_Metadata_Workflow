#!/usr/bin/env nextflow
nextflow.enable.dsl=2

TOOL_FOLDER = "$baseDir/bin"

process processData {
    publishDir "./nf_output", mode: 'copy'

    conda "$TOOL_FOLDER/conda_env.yml"

    input:
    val x

    output:
    file 'outputs.tsv.*'
    file 'file_paths.tsv'
    file 'passed_file_names.tsv'
    file 'check.tsv'

    """
    python $TOOL_FOLDER/gnps_downloader.py
    python $TOOL_FOLDER/gnps_validator.py passed_file_names.tsv
    python $TOOL_FOLDER/gnps_name_matcher.py check.tsv
    """
}

workflow {
    processData(1)
}
