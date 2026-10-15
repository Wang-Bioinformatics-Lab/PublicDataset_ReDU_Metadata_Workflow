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

// Cellosaurus cell lines (CellLine / CellLineIndex): a pinned table shipped in data/
// (built by bin/build_cellosaurus_table.py; Expasy only serves the current release).
// Release, sha256 and the CC BY 4.0 attribution are in data/cellosaurus_<version>.json.
CELLOSAURUS_TABLE = 'cellosaurus_56.0_cell_lines.tsv.gz'


// First step: sets up all reference data in data/ (ontologies, NCBI taxonomy, MASST
// tables) in the versions pinned in data/reference_data.json (url + sha256), i.e. the
// versions the workflow was tested with. Files already present in the pinned version are
// kept; missing, different or unrecorded files (e.g. from an earlier manual get_data.sh
// run) are replaced by the pinned version. Newer upstream releases are only used after
// their pin is changed. Outputs are links to the files in data/.
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
    path 'food_masst_table.csv', emit: food_masst
    path 'food_masst_table_full.tsv', emit: food_masst_full
    path 'personalCareProduct_masst_table.csv', emit: personal_care_masst
    path 'microbiome_masst_table.tsv', emit: microbiome_masst
    path 'names.dmp', emit: names
    path 'nodes.dmp', emit: nodes
    path 'division.dmp', emit: division
    path 'merged.dmp', emit: merged
    path 'hp.obo', emit: hp
    path 'ecto.obo', emit: ecto
    path 'peco.obo', emit: peco

    """
    python $TOOL_FOLDER/fetch_reference_data.py $DATA_FOLDER/reference_data.json $DATA_FOLDER --link_dir .
    """
}


// One copy of the GNPS2 dataset cache's file table (https://datasetcache.gnps2.org/dataset/uniquemri,
// a ~1 GB export the cache regenerates about every 3 h; ~3 min). The name matching steps take each
// dataset's mz(X)ML files from it instead of one dataset cache query per dataset (~1700 queries,
// 20-90 min per run; a query that timed out dropped all rows of its dataset). It is fetched once the
// MassIVE listing is done (so the dataset cache is asked one request at a time), while gnpsHarmonize
// runs. If the download fails or is incomplete the file is left empty and the name matching steps
// query the cache per dataset, as before.
process fetchDatasetCache {
    conda "$TOOL_FOLDER/conda_env.yml"

    // the table changes every few hours; never reuse one from an earlier run (-resume)
    cache false

    input:
    val after_massive_download

    output:
    path 'uniquemri.csv', emit: table

    """
    for attempt in 1 2 3; do
        problem=""
        if curl -sS --fail --speed-limit 50000 --speed-time 300 --max-time 5400 \
                -D uniquemri.headers -o uniquemri.csv https://datasetcache.gnps2.org/dataset/uniquemri; then
            expected=\$(grep -i '^content-length:' uniquemri.headers | tail -n 1 | tr -dc '0-9' || true)
            actual=\$(stat -c %s uniquemri.csv)
            if [ -n "\$expected" ] && [ "\$expected" != "\$actual" ]; then
                problem="incomplete download: \$actual of \$expected bytes"
            elif [ "\$(head -c 21 uniquemri.csv)" != 'usi,filepath,dataset,' ]; then
                problem="unexpected format"
            elif [ "\$(tail -c 1 uniquemri.csv | od -An -c | tr -d ' ')" != '\\n' ]; then
                problem="the file does not end with a complete line"
            fi
        else
            problem="download failed (curl exit \$?)"
        fi
        if [ -z "\$problem" ]; then break; fi
        echo "dataset cache table, attempt \$attempt: \$problem" >&2
        if [ \$attempt -lt 3 ]; then sleep 60; fi
    done
    if [ -n "\$problem" ]; then
        echo "WARNING: dataset cache table not available (\$problem); the name matching queries the dataset cache per dataset" >&2
        : > uniquemri.csv
    fi
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
    path cellosaurus_tsv

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
    --path_to_peco_obo ${peco_obo} \
    --path_to_cellosaurus_tsv ${cellosaurus_tsv}
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
    path doid_parents_csv
    path hpo_csv
    path exposure_csv
    path cellosaurus_tsv
    path nodes_dmp
    path names_dmp
    path merged_dmp
    val after_listing

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
    --path_to_doid_parents_csv ${doid_parents_csv} \
    --path_to_hpo_csv ${hpo_csv} \
    --path_to_exposure_csv ${exposure_csv} \
    --path_to_cellosaurus_tsv ${cellosaurus_tsv} \
    --path_to_ncbi_nodes_dmp ${nodes_dmp} \
    --path_to_ncbi_names_dmp ${names_dmp} \
    --path_to_ncbi_merged_dmp ${merged_dmp} \
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
    path doid_parents_csv
    path hpo_csv
    path exposure_csv
    // staged for the MetaboLights cell lines; Metabolights2REDU.py does not take
    // --path_to_cellosaurus_tsv yet (its CellLine / CellLineIndex are 'missing value')
    path cellosaurus_tsv
    // pinned NCBI taxdump (fetchReferenceData): lineages for SampleType without E-utilities
    path nodes_dmp
    path names_dmp
    path merged_dmp
    val after_listing

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
    --path_to_doid_parents_csv ${doid_parents_csv} \
    --path_to_hpo_csv ${hpo_csv} \
    --path_to_exposure_csv ${exposure_csv} \
    --path_to_ncbi_nodes_dmp ${nodes_dmp} \
    --path_to_ncbi_names_dmp ${names_dmp} \
    --path_to_ncbi_merged_dmp ${merged_dmp} \
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
    --include_compressed \
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
    path masst_supplement

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
    --masst_supplement $masst_supplement \
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
    path 'DOID_parents.csv'
    path 'UBERON_sex_specific_parts.csv'

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
    path hpo_csv
    path exposure_csv
    path names_dmp
    path cellosaurus_tsv

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
    --path_to_hpo_csv ${hpo_csv} \
    --path_to_exposure_csv ${exposure_csv} \
    --path_to_ncbi_names_dmp ${names_dmp} \
    --path_to_cellosaurus_tsv ${cellosaurus_tsv} \
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
    path files_table

    output:
    file 'gnps_metadata_all.tsv'

    """
    timeout ${params.source_timeout} python $TOOL_FOLDER/gnps_name_matcher.py \
    all \
    adjusted_metadata_folder \
    gnps_metadata_all.tsv \
    ${allowed_terms} \
    --files_table ${files_table} \
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


// MASST tables (microbe, plant, food, personal care products, microbiome) in ReDU format, for
// files no GNPS metadata file describes; a missing or empty table is skipped with a warning
process MASST_to_REDU {
    publishDir "./nf_output", mode: 'copy'

    conda "$TOOL_FOLDER/conda_env.yml"

    input:
    path redu_table
    path ncbi_rank_division
    path allowed_terms
    path microbe_masst_table
    path plant_masst_table
    path food_masst_table
    path food_masst_table_full
    path personal_care_masst_table
    path microbiome_masst_table
    path nodes_dmp
    path names_dmp
    path merged_dmp
    path uberon_cl_po_csv
    path doid_csv
    path exposure_csv
    path files_table

    output:
    path 'adjusted_metadata_folder', emit: folder
    // MASST rows of files GNPS metadata describes: they only fill missing fields in the merge
    path 'masst_supplement.tsv', emit: supplement

    """
    mkdir adjusted_metadata_folder
    : > masst_supplement.tsv
    timeout ${params.source_timeout} python $TOOL_FOLDER/MASST_to_REDU.py \
    ${redu_table} \
    adjusted_metadata_folder \
    --path_microbeMASST ${microbe_masst_table} \
    --path_plantMASST ${plant_masst_table} \
    --path_foodMASST ${food_masst_table} \
    --path_foodMASST_full ${food_masst_table_full} \
    --path_personalCareProductMASST ${personal_care_masst_table} \
    --path_microbiomeMASST ${microbiome_masst_table} \
    --path_ncbiRanksDivisions ${ncbi_rank_division} \
    --AllowedTermJson_path ${allowed_terms} \
    --path_ncbi_nodes ${nodes_dmp} \
    --path_ncbi_names ${names_dmp} \
    --path_ncbi_merged ${merged_dmp} \
    --path_to_uberon_cl_po_csv ${uberon_cl_po_csv} \
    --path_to_doid_csv ${doid_csv} \
    --path_to_exposure_csv ${exposure_csv} \
    --files_table ${files_table} \
    --supplement_output masst_supplement.tsv \
    || { echo "WARNING: MASST_to_REDU failed or timed out; MASST rows from earlier runs are kept" >&2; : > masst_supplement.tsv; }
    """
}

// This is taking metadata from microbeMASST and PlatnMASST and modifying it
process gnpsmatchName_masst {
    publishDir "./nf_output", mode: 'copy'

    conda "$TOOL_FOLDER/conda_env.yml"


    input:
    file 'adjusted_metadata_folder' 
    path allowed_terms
    path files_table

    output:
    file 'masst_metadata_all.tsv'

    """
    timeout ${params.source_timeout} python $TOOL_FOLDER/gnps_name_matcher.py \
    all \
    adjusted_metadata_folder \
    masst_metadata_all.tsv \
    ${allowed_terms} \
    --files_table ${files_table} \
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
    path files_table

    output:
    path 'merged_with_old.tsv', emit: merged
    path 'source_summary.tsv'
    path 'saveolder_dropped.tsv'

    """
    python $TOOL_FOLDER/save_older_data.py \
    ${merged_ch} \
    ${old_redu} \
    merged_with_old.tsv \
    --metabolights_files ${ml_files} \
    --files_table ${files_table} \
    --summary_output source_summary.tsv \
    --dropped_output saveolder_dropped.tsv
    """
}



// Lineage columns (NCBISuperkingdom ... NCBISpecies) from the pinned NCBI taxonomy
// (nodes/names/merged.dmp of reference_data.json), the release all other steps use
process enrich_ncbi_information {
    publishDir "./nf_output", mode: 'copy'

    conda "$TOOL_FOLDER/conda_env.yml"

    input:
    path merged_ch
    path nodes_dmp
    path names_dmp
    path merged_dmp

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
    --path_to_ncbi_nodes_dmp ${nodes_dmp} \
    --path_to_ncbi_names_dmp ${names_dmp} \
    --path_to_ncbi_merged_dmp ${merged_dmp}
    """
}


// Rules for every row of the final table, also the rows kept from earlier runs (final_table_rules.py):
// personal care products, BiologicalSex of mammal samples from a body part only one sex has
process finalTableRules {
    publishDir "./nf_output", mode: 'copy'

    conda "$TOOL_FOLDER/conda_env.yml"

    input:
    path merged_ch
    path sex_specific_parts

    output:
    path 'merged_with_rules.tsv'

    """
    python $TOOL_FOLDER/final_table_rules.py \
    ${merged_ch} \
    merged_with_rules.tsv \
    --path_to_sex_specific_parts ${sex_specific_parts}
    """
}


// Last edits of SampleGroupFreeText and of the columns its factors are moved into, on every row (also the
// rows saveOlderData keeps from earlier runs), after finalTableRules: Diet / FastingStatus /
// InternalStandardsUsed from the texts, removal of the factors these columns hold, the group-consistency rule
process finalizeTextColumns {
    publishDir "./nf_output", mode: 'copy'

    conda "$TOOL_FOLDER/conda_env.yml"

    input:
    path merged_ch
    path allowed_terms
    path exposure_csv

    output:
    path 'merged_with_text_columns.tsv'

    """
    python $TOOL_FOLDER/final_text_columns.py \
    ${merged_ch} \
    merged_with_text_columns.tsv \
    --path_to_allowed_term_json ${allowed_terms} \
    --path_to_exposure_csv ${exposure_csv}
    """
}


// Files tracked in data/ (allowed terms seed, Cellosaurus table); fail early with a clear message if one is missing
def dataFile(name) {
    file("$DATA_FOLDER/$name", checkIfExists: true)
}

workflow {

    //  Reference data, then ontologies
    ref = fetchReferenceData()
    (uberon_cl_co_onto, doid_onto, envo_bio, envo_material, ncbi_rank_division, hpo_onto, exposure_onto, doid_parents, sex_specific_parts) = prepare_ontologies(
        ref.uberon, ref.cl, ref.po, ref.doid, ref.biome, ref.material,
        ref.nodes, ref.division, ref.hp, ref.ecto, ref.peco)
    cellosaurus = dataFile(CELLOSAURUS_TABLE)
    allowed_terms = updateAllowedTerms(
        dataFile('allowed_terms.json'), ref.names, ref.uberon, ref.po, ref.cl, ref.doid, ref.ms,
        ref.material, ref.biome, ref.hp, ref.ecto, ref.peco, cellosaurus)

    // Massive REDU data, called before GitHub because taking it from MassIVE as the place to keep metadata and not github
    (file_paths_ch, metadata_ch) = downloadMetadata_massive_and_github(1)
    // the dataset cache's file table for the name matching steps, after the MassIVE listing
    dataset_cache = fetchDatasetCache(file_paths_ch)
    msv_metadata_ch = gnpsHarmonize(metadata_ch, uberon_cl_co_onto, doid_onto, envo_bio, envo_material, ncbi_rank_division, allowed_terms, hpo_onto, exposure_onto, ref.names, cellosaurus)
    gnps_metadata_ch = gnpsmatchName(msv_metadata_ch, allowed_terms, dataset_cache.table)
    gnps_metadata_instruments_ch = gnpsInstruments(gnps_metadata_ch, allowed_terms)

    // MASST tables: microbeMASST, plantMASST, foodMASST, personalCareProductMASST, microbiomeMASST
    masst_metadata_ch = MASST_to_REDU(gnps_metadata_ch, ncbi_rank_division, allowed_terms,
        ref.microbe_masst, ref.plant_masst, ref.food_masst, ref.food_masst_full, ref.personal_care_masst,
        ref.microbiome_masst, ref.nodes, ref.names, ref.merged, uberon_cl_co_onto, doid_onto, exposure_onto, dataset_cache.table)
    masst_metadata_wFiles_ch = gnpsmatchName_masst(masst_metadata_ch.folder, allowed_terms, dataset_cache.table)

    // Metabolomics Workbench
    // one request stream per repository: the file listing first, then the converter
    mwb_files_ch = mwbFiles(1)
    mwb_metadata_ch = mwbRun(uberon_cl_co_onto, envo_bio, envo_material, ncbi_rank_division, allowed_terms, doid_onto, doid_parents, hpo_onto, exposure_onto, cellosaurus,
        ref.nodes, ref.names, ref.merged, mwb_files_ch)
    mwb_redu_ch = formatmwb(mwb_metadata_ch, mwb_files_ch)

    // Metabolights
    // one request stream per repository: the file listing first, then the converter
    ml_files_ch = mlFiles(1)
    ml_metadata_ch = mlRun(uberon_cl_co_onto, envo_bio, envo_material, ncbi_rank_division, allowed_terms, doid_onto, doid_parents, hpo_onto, exposure_onto, cellosaurus,
        ref.nodes, ref.names, ref.merged, ml_files_ch)
    ml_redu_ch = formatml(ml_metadata_ch, ml_files_ch)

    // NORMAN
    norman_metadata_ch = normanRun(uberon_cl_co_onto, envo_bio, envo_material, ncbi_rank_division, allowed_terms)

    // Combine everything
    merged_ch = mergeAllMetadata(allowed_terms, gnps_metadata_instruments_ch, mwb_redu_ch, ml_redu_ch, norman_metadata_ch, masst_metadata_wFiles_ch, masst_metadata_ch.supplement)


    // Make sure we dont loose older data
    old_redu_path_ch = Channel.fromPath(params.old_redu)
    merged_with_old_ch = saveOlderData(merged_ch, old_redu_path_ch, ml_files_ch, dataset_cache.table).merged

    // Enrich NCBI information
    merged_with_ncbi_ch = enrich_ncbi_information(merged_with_old_ch, ref.nodes, ref.names, ref.merged)

    // Rules for every row (needs the lineage)
    merged_with_rules_ch = finalTableRules(merged_with_ncbi_ch, sex_specific_parts)

    // Diet / FastingStatus / InternalStandardsUsed from the sample group texts; SampleGroupFreeText rules
    merged_with_text_ch = finalizeTextColumns(merged_with_rules_ch, allowed_terms, exposure_onto)

    // Add Cache columns
    cache_enriched_metadata = add_cache_columns(merged_with_text_ch)

}
