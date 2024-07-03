#!/bin/bash

# Get the directory where the script is located
SCRIPT_DIR="$(dirname "$0")"

# Function to download a file if it does not already exist and echo its download status
download_file() {
    local url=$1
    local path=$2
    if [ ! -f "$path" ]; then
        echo "Downloading $(basename "$path")..."
        curl -o "$path" "$url"
        echo "$(basename "$path") has been downloaded to ${SCRIPT_DIR}"
    else
        echo "$(basename "$path") already exists. Skipping download."
    fi
}

# Download commands with check for each file
download_file "http://aber-owl.net/media/ontologies/UBERON/298/uberon.owl" "${SCRIPT_DIR}/uberon.owl"
download_file "http://aber-owl.net/media/ontologies/PO/26/po.owl" "${SCRIPT_DIR}/po.owl"
download_file "http://aber-owl.net/media/ontologies/CL/104/cl.owl" "${SCRIPT_DIR}/cl.owl"
download_file "http://aber-owl.net/media/ontologies/DOID/678/doid.owl" "${SCRIPT_DIR}/doid.owl"
download_file "http://aber-owl.net/media/ontologies/MS/188/ms.owl" "${SCRIPT_DIR}/ms.owl"
download_file "https://raw.githubusercontent.com/EnvironmentOntology/envo/master/subsets/biome-hierarchy.owl" "${SCRIPT_DIR}/biome-hierarchy.owl"
download_file "https://raw.githubusercontent.com/EnvironmentOntology/envo/master/subsets/material-hierarchy.owl" "${SCRIPT_DIR}/material-hierarchy.owl"
download_file "https://raw.githubusercontent.com/robinschmid/microbe_masst/master/data/microbe_masst_table.csv" "${SCRIPT_DIR}/microbe_masst_table.csv"
download_file "https://raw.githubusercontent.com/robinschmid/microbe_masst/master/data/plant_masst_table.csv" "${SCRIPT_DIR}/plant_masst_table.csv"
download_file "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdmp.zip" "${SCRIPT_DIR}/taxdmp.zip"

# Unzip taxdmp.zip if it has been downloaded
if [ -f "${SCRIPT_DIR}/taxdmp.zip" ]; then
    echo "Unzipping taxdmp.zip..."
    unzip -j "${SCRIPT_DIR}/taxdmp.zip" "names.dmp" -d "${SCRIPT_DIR}"
    echo "names.dmp has been extracted."
    unzip -j "${SCRIPT_DIR}/taxdmp.zip" "nodes.dmp" -d "${SCRIPT_DIR}"
    echo "nodes.dmp has been extracted."
    unzip -j "${SCRIPT_DIR}/taxdmp.zip" "division.dmp" -d "${SCRIPT_DIR}"
    echo "division.dmp has been extracted."
    
    # Remove the zip file to clean up
    rm "${SCRIPT_DIR}/taxdmp.zip"
    echo "taxdmp.zip has been removed."
else
    echo "taxdmp.zip not found. Skipping extraction."
fi

echo "All downloads and extraction completed."
