#!/usr/bin/env nextflow
nextflow.enable.dsl=2

TOOL_FOLDER = "$baseDir/bin"
DATA_FOLDER = "$baseDir/data"


//This is default empty file, we will overwrite this in web application
params.old_redu = './data/empty_redu.tsv'

// Each source (MassIVE/GitHub, Workbench, MetaboLights, NORMAN, MASST) is retrieved
// independently. If one fails or exceeds this limit, its step emits an empty result and
// saveOlderData keeps that source's rows from earlier runs, so a problem with one
// repository or dataset never blocks the update of the others. The limit stays below
// the 23.5 h celery time limit so the merge steps still run.
params.source_timeout = '20h'


// First step: sets up all reference data in data/ (ontologies, NCBI taxonomy, MASST
// tables) as listed in data/reference_data.json. Files already present in the listed
// version are kept; missing, outdated or unrecorded files (e.g. from an earlier manual
// get_data.sh run) are downloaded. Outputs are links to the files in data/.
process fetchReferenceData {
    conda "$TOOL_FOLDER/conda_env.yml"

    // always check the sources, also with -resume
    cache false

    output:
    path 'uberon.owl', emit: uberon
    path 'po.owl', emit: po
    path 'cl.owl', emit: cl
    path 'doid.owl', emit: doid
    path 'ms.owl', emit: ms
    path 'biome-hierarchy.owl', emit: biome
    path 'material-hierarchy.owl', emit: material
    path 'microbe_masst_table.csv', emit: microbe_masst
    path 'plant_masst_table.csv', emit: plant_masst
    path 'names.dmp', emit: names
    path 'nodes.dmp', emit: nodes
    path 'division.dmp', emit: division
    path 'hp.obo', emit: hp
    path 'ecto.obo', emit: ecto
    path 'peco.obo', emit: peco

    """
    python $TOOL_FOLDER/fetch_reference_data.py $DATA_FOLDER/reference_data.json $DATA_FOLDER --link_dir .
    """
}


process updateAllowedTerms {
    publishDir "./nf_output", mode: 'copy'

    conda "$TOOL_FOLDER/conda_env.yml"

    // Declared as inputs (not plain strings) so that updating them invalidates the cache
    input:
    path allowed_terms_seed, stageAs: 'allowed_terms_seed.json'
    path names_dmp
    path uberon_owl
    path po_owl
    path cl_owl
    path doid_owl
    path ms_owl
    path material_owl
    path biome_owl
    path hp_obo
    path ecto_obo
    path peco_obo

    output:
    file 'allowed_terms.json'

    """
    python $TOOL_FOLDER/update_allowed_terms_from_ontologies.py  \
    ${allowed_terms_seed} \
    --path_to_ncbi_dump ${names_dmp} \
    --path_to_uberon_owl ${uberon_owl} \
    --path_to_po_owl ${po_owl} \
    --path_to_cl_owl ${cl_owl} \
    --path_to_doid_owl ${doid_owl} \
    --path_to_ms_owl ${ms_owl} \
    --path_to_material_envs_owl ${material_owl} \
    --path_to_biome_envs_owl ${biome_owl} \
    --path_to_hp_obo ${hp_obo} \
    --path_to_ecto_obo ${ecto_obo} \
    --path_to_peco_obo ${peco_obo}
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
    if wget -q --tries=3 --timeout=120 https://github.com/Wang-Bioinformatics-Lab/ReDU_metadata/archive/refs/heads/main.zip && unzip -q main.zip; then
        mv ReDU_metadata-main/metadata/*.tsv metadata_folder/
    else
        echo "WARNING: could not fetch ReDU_metadata from GitHub; its rows from earlier runs are kept" >&2
    fi
    timeout ${params.source_timeout} python $TOOL_FOLDER/gnps_downloader.py metadata_folder \
    || echo "WARNING: MassIVE metadata download failed; MassIVE rows from earlier runs are kept" >&2
    [ -f file_paths.tsv ] || printf 'LocalName\tServerName\n' > file_paths.tsv
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
    path doid_csv

    output:
    file 'REDU_from_MWB_all.tsv'

    """
    timeout ${params.source_timeout} python $TOOL_FOLDER/MWB_to_REDU.py \
    --study_id ALL \
    --path_to_csvs $TOOL_FOLDER/translation_sheets \
    --path_to_allowed_term_json ${allowed_terms} \
    --duplicate_raw_file_handling keep_all \
    --path_to_uberon_cl_po_csv ${uberon_po_cl_csv_path} \
    --path_to_envo_biome_csv ${ENVO_bio_csv} \
    --path_to_envo_material_csv ${ENVO_material_csv} \
    --path_ncbi_rank_division ${ncbi_rank_division} \
    --path_to_polarity_info $DATA_FOLDER/MWB_polarity_table.csv \
    --path_to_doid_csv ${doid_csv} \
    || { echo "WARNING: mwbRun failed or timed out; Workbench rows from earlier runs are kept" >&2; : > REDU_from_MWB_all.tsv; }
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
    timeout ${params.source_timeout} python $TOOL_FOLDER/getAllWorkbench_file_paths.py \
    --study_id ALL \
    --output_path mwb_files_all.tsv \
    --filter_extensions True \
    || { echo "WARNING: mwbFiles failed or timed out; Workbench rows from earlier runs are kept" >&2; : > mwb_files_all.tsv; }
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
    timeout ${params.source_timeout} python $TOOL_FOLDER/MWB_merge.py  \
    $mwb_metadata \
    $mwb_files \
    mwb_redu.tsv \
    || { echo "WARNING: formatmwb failed or timed out; Workbench rows from earlier runs are kept" >&2; : > mwb_redu.tsv; }
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
    path doid_csv

    output:
    file 'Metabolights2REDU_ALL.tsv'


    """
    timeout ${params.source_timeout} python $TOOL_FOLDER/Metabolights2REDU.py \
    --study_id ALL  \
    --path_to_translation_sheet_csvs $TOOL_FOLDER/translation_sheets_metabolights \
    --path_to_allowed_term_json ${allowed_terms} \
    --path_to_uberon_cl_po_csv ${uberon_po_cl_csv_path} \
    --path_to_envo_biome_csv ${ENVO_bio_csv} \
    --path_to_envo_material_csv ${ENVO_material_csv} \
    --path_ncbi_rank_division ${ncbi_rank_division} \
    --path_to_doid_csv ${doid_csv} \
    || { echo "WARNING: mlRun failed or timed out; MetaboLights rows from earlier runs are kept" >&2; : > Metabolights2REDU_ALL.tsv; }
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
    timeout ${params.source_timeout} python $TOOL_FOLDER/getAllMetabolights_file_paths.py \
    --output_filename MetabolightsFilePaths_ALL.tsv \
    || { echo "WARNING: mlFiles failed or timed out; MetaboLights rows from earlier runs are kept" >&2; : > MetabolightsFilePaths_ALL.tsv; }
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
    timeout ${params.source_timeout} python $TOOL_FOLDER/ML_merge.py \
    $ml_metadata \
    $ml_files \
    ml_redu.tsv \
    || { echo "WARNING: formatml failed or timed out; MetaboLights rows from earlier runs are kept" >&2; : > ml_redu.tsv; }
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
    timeout ${params.source_timeout} python $TOOL_FOLDER/NORMAN_to_REDU.py \
    --study_id ALL  \
    --path_to_allowed_term_json ${allowed_terms} \
    --path_to_uberon_cl_po_csv ${uberon_po_cl_csv_path} \
    --path_to_envo_biome_csv ${ENVO_bio_csv} \
    --path_to_envo_material_csv ${ENVO_material_csv} \
    --path_ncbi_rank_division ${ncbi_rank_division} \
    --output NORMAN2REDU_ALL.tsv \
    || { echo "WARNING: normanRun failed or timed out; NORMAN rows from earlier runs are kept" >&2; : > NORMAN2REDU_ALL.tsv; }
    """
}

process mergeAllMetadata {
    publishDir "./nf_output", mode: 'copy'

    conda "$TOOL_FOLDER/conda_env.yml"

    input:
    path allowed_terms
    path gnps_metadata
    path mwb_redu
    path metabolights_redu
    path norman_redu
    path masst_metadata

    output:
    path 'merged_metadata_wo_cache_columns.tsv'

    """
    python $TOOL_FOLDER/merge_metadata.py \
    --path_to_allowed_term_json ${allowed_terms} \
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
    path uberon_owl
    path cl_owl
    path po_owl
    path doid_owl
    path biome_owl
    path material_owl
    path nodes_dmp
    path division_dmp
    path hp_obo
    path ecto_obo
    path peco_obo

    output:
    path 'UBERON_CL_PO_ontology.csv'
    path 'DOID_ontology.csv'
    path 'ENVO_biome_ontology.csv'
    path 'ENVO_material_ontology.csv'
    path 'NCBI_Rank_Division.csv'
    path 'HPO_ontology.csv'
    path 'Exposure_ontology.csv'

    """
    python $TOOL_FOLDER/prepare_ontologies.py \
    --path_to_uberon_owl ${uberon_owl} \
    --path_to_cl_owl ${cl_owl} \
    --path_to_po_owl ${po_owl} \
    --path_to_doid_owl ${doid_owl} \
    --path_to_biome_envs_owl ${biome_owl} \
    --path_to_material_envs_owl ${material_owl} \
    --path_to_ncbi_nodes_dmp ${nodes_dmp} \
    --path_to_ncbi_division_dmp ${division_dmp} \
    --path_to_hp_obo ${hp_obo} \
    --path_to_ecto_obo ${ecto_obo} \
    --path_to_peco_obo ${peco_obo}
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
    timeout ${params.source_timeout} python $TOOL_FOLDER/read_and_validate_redu_from_github.py \
    ${metadata_ch} \
    adjusted_metadata_folder \
    --AllowedTermJson_path ${allowed_terms} \
    --path_ncbi_rank_division ${ncbi_rank_division} \
    --path_to_uberon_cl_po_csv ${UBERON_CL_PO_ontology_csv} \
    --path_to_doid_csv ${DOID_ontology_csv} \
    --path_to_envo_biome_csv ${ENVO_bio_csv} \
    --path_to_envo_material_csv ${ENVO_material_csv} \
    || { echo "WARNING: gnpsHarmonize failed or timed out; GNPS rows from earlier runs are kept" >&2; true; }
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
    timeout ${params.source_timeout} python $TOOL_FOLDER/gnps_name_matcher.py \
    all \
    adjusted_metadata_folder \
    gnps_metadata_all.tsv \
    ${allowed_terms} \
    || { echo "WARNING: gnpsmatchName failed or timed out; GNPS rows from earlier runs are kept" >&2; : > gnps_metadata_all.tsv; }
    """
}


// Fills missing MassSpectrometer values from the MassIVE dataset description (only when
// it lists exactly one instrument that maps unambiguously). Optional enrichment: on any
// failure the rows are passed through unchanged.
process gnpsInstruments {
    publishDir "./nf_output", mode: 'copy'

    conda "$TOOL_FOLDER/conda_env.yml"

    input:
    path gnps_metadata
    path allowed_terms

    output:
    file 'gnps_metadata_instruments.tsv'

    """
    timeout 2h python $TOOL_FOLDER/massive_dataset_instruments.py \
    ${gnps_metadata} \
    gnps_metadata_instruments.tsv \
    --path_to_allowed_term_json ${allowed_terms} \
    || { echo "WARNING: gnpsInstruments failed or timed out; passing rows through" >&2; cp ${gnps_metadata} gnps_metadata_instruments.tsv; }
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
    path microbe_masst_table
    path plant_masst_table
    path nodes_dmp

    output:
    file 'adjusted_metadata_folder'

    """
    mkdir adjusted_metadata_folder
    timeout ${params.source_timeout} python $TOOL_FOLDER/MASST_to_REDU.py \
    ${redu_table} \
    adjusted_metadata_folder \
    --path_microbeMASST ${microbe_masst_table} \
    --path_plantMASST ${plant_masst_table} \
    --path_ncbiRanksDivisions ${ncbi_rank_division} \
    --AllowedTermJson_path ${allowed_terms} \
    --path_ncbi_nodes ${nodes_dmp} \
    || { echo "WARNING: MASST_to_REDU failed or timed out; MASST rows from earlier runs are kept" >&2; true; }
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
    timeout ${params.source_timeout} python $TOOL_FOLDER/gnps_name_matcher.py \
    all \
    adjusted_metadata_folder \
    masst_metadata_all.tsv \
    ${allowed_terms} \
    || { echo "WARNING: gnpsmatchName_masst failed or timed out; MASST rows from earlier runs are kept" >&2; : > masst_metadata_all.tsv; }
    """
}

process add_cache_columns {
    // move (a rename within this filesystem) so the web app never reads a half-copied file
    publishDir "./nf_output", mode: 'move'

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
    path ml_files

    output:
    path 'merged_with_old.tsv', emit: merged
    path 'source_summary.tsv'

    """
    python $TOOL_FOLDER/save_older_data.py \
    ${merged_ch} \
    ${old_redu} \
    merged_with_old.tsv \
    --metabolights_files ${ml_files} \
    --summary_output source_summary.tsv
    """
}



process enrich_ncbi_information {
    publishDir "./nf_output", mode: 'copy'

    conda "$TOOL_FOLDER/conda_env_ete3.yml"

    input:
    path merged_ch

    output:
    path 'merged_with_ncbi.tsv'

    //https://qlever.dev/wikidata/adgYjs
    """
    curl -sS --fail --max-time 600 https://qlever.dev/api/wikidata \
    -H "Accept: text/tab-separated-values" \
    -H "Content-type: application/sparql-query" \
    --data "PREFIX wd: <http://www.wikidata.org/entity/> PREFIX wdt: <http://www.wikidata.org/prop/direct/> PREFIX bd: <http://www.bigdata.com/rdf#> PREFIX wikibase: <http://wikiba.se/ontology#> SELECT ?ncbi_id ?ott_id WHERE { ?item wdt:P685 ?ncbi_id . OPTIONAL { ?item wdt:P9157 ?ott_id } ?item wdt:P31 wd:Q16521 }" \
    -o sparql_oto_ncbi.tsv || rm -f sparql_oto_ncbi.tsv

    python $TOOL_FOLDER/ncbi_database.py \
    ${merged_ch}  \
    merged_with_ncbi.tsv \
    || { echo "WARNING: NCBI lineage annotation failed; passing rows through without it" >&2; cp ${merged_ch} merged_with_ncbi.tsv; }
    """
}


// Files tracked in data/ (allowed terms seed); fail early with a clear message if one is missing
def dataFile(name) {
    file("$DATA_FOLDER/$name", checkIfExists: true)
}

workflow {

    //  Reference data, then ontologies
    ref = fetchReferenceData()
    (uberon_cl_co_onto, doid_onto, envo_bio, envo_material, ncbi_rank_division, hpo_onto, exposure_onto) = prepare_ontologies(
        ref.uberon, ref.cl, ref.po, ref.doid, ref.biome, ref.material,
        ref.nodes, ref.division, ref.hp, ref.ecto, ref.peco)
    allowed_terms = updateAllowedTerms(
        dataFile('allowed_terms.json'), ref.names, ref.uberon, ref.po, ref.cl, ref.doid, ref.ms,
        ref.material, ref.biome, ref.hp, ref.ecto, ref.peco)

    // Massive REDU data, called before GitHub because taking it from MassIVE as the place to keep metadata and not github
    (file_paths_ch, metadata_ch) = downloadMetadata_massive_and_github(1)
    msv_metadata_ch = gnpsHarmonize(metadata_ch, uberon_cl_co_onto, doid_onto, envo_bio, envo_material, ncbi_rank_division, allowed_terms)
    gnps_metadata_ch = gnpsmatchName(msv_metadata_ch, allowed_terms)
    gnps_metadata_instruments_ch = gnpsInstruments(gnps_metadata_ch, allowed_terms)

    // MicrobeMASST and PlantMASST
    masst_metadata_ch = MASST_to_REDU(gnps_metadata_ch, ncbi_rank_division, allowed_terms,
        ref.microbe_masst, ref.plant_masst, ref.nodes)
    masst_metadata_wFiles_ch = gnpsmatchName_masst(masst_metadata_ch, allowed_terms)

    // Metabolomics Workbench
    mwb_metadata_ch = mwbRun(uberon_cl_co_onto, envo_bio, envo_material, ncbi_rank_division, allowed_terms, doid_onto)
    mwb_files_ch = mwbFiles(1)
    mwb_redu_ch = formatmwb(mwb_metadata_ch, mwb_files_ch)

    // Metabolights
    ml_metadata_ch = mlRun(uberon_cl_co_onto, envo_bio, envo_material, ncbi_rank_division, allowed_terms, doid_onto)
    ml_files_ch = mlFiles(1)
    ml_redu_ch = formatml(ml_metadata_ch, ml_files_ch)

    // NORMAN
    norman_metadata_ch = normanRun(uberon_cl_co_onto, envo_bio, envo_material, ncbi_rank_division, allowed_terms)

    // Combine everything
    merged_ch = mergeAllMetadata(allowed_terms, gnps_metadata_instruments_ch, mwb_redu_ch, ml_redu_ch, norman_metadata_ch, masst_metadata_wFiles_ch)


    // Make sure we dont loose older data
    old_redu_path_ch = Channel.fromPath(params.old_redu)
    merged_with_old_ch = saveOlderData(merged_ch, old_redu_path_ch, ml_files_ch).merged

    // Enrich NCBI information
    merged_with_ncbi_ch = enrich_ncbi_information(merged_with_old_ch)

    // Add Cache columns
    cache_enriched_metadata = add_cache_columns(merged_with_ncbi_ch)

}
