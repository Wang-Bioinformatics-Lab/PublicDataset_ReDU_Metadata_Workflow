**GNPS Metadata Tools** <br>
This repository contains three Python scripts for working with GNPS metadata files: _gnps_downloader.py_, _gnps_validator.py_, and _gnps_name_matcher.py_. These scripts allow you to aggregate lists of GNPS metadata files, validate the files using the GNPS metadata validator, and match the files to their respective CCMS peak dataset names.<br>

**Installation**<br>
To use these scripts, you'll need to have Python 3 installed on your system. <br>
You can download Python from the official Python website: https://www.python.org/downloads/<br>

-> You can install dependencies using pip:<br>
_pip install pandas urllib requests_<br>
<br>
**Usage**<br>
_**gnps_downloader.py**_<br>
This script aggregates a list of GNPS metadata files, sorts the files by their creation time, and downloads the latest GNPS metadata file. The script then appends the file path and file name into a TSV file.<br>

To run the script, use the following command: <br>
_python3 gnps_downloader.py_

_**gnps_validator.py**_<br>
This script runs the downloaded GNPS metadata files against the metadata validator and stores the list of file names that have passed through the validator. The script also rejects files that haven't passed and appends the passed file names into a TSV file.

To run the script, use the following command:<br>
_python3 gnps_validator.py_

_**gnps_name_matcher.py**_<br>
This script matches the GNPS metadata files to their respective CCMS peak dataset names and gives out a TSV file that contains all the names that match unambiguously.

To run the script, use the following command:<br>
_python3 gnps_name_matcher.py_


_**data/allowed_terms.json**_<br>
The terms allowed in REDU are pulled from this json. Terms from controlled ontologies for variables MassSpectrometer, NCBITaxonomy, UBERONBodyPartName, DOIDCommonName, HPOCommonName (HPO) and Exposure (ECTO, PECO) are added to json within the workflow. The reference data (ontologies, NCBI taxonomy, MASST tables) is set up automatically by the first workflow step from data/reference_data.json: each file is pinned to the version the workflow was tested with (url + sha256); files already present in that version are kept, others are downloaded (data/get_data.sh does the same by hand). The NCBI lineage columns (NCBISuperkingdom ... NCBISpecies) are computed from the same pinned NCBI taxonomy. One deliberate exception: OpenTreeOfLifeTaxonomyID comes from a live Wikidata query on every run, like the MASST tables. Newer ontology or NCBI releases are not used automatically: to move to one, test the workflow with it, then change its url and sha256 in reference_data.json. Additional terms can be added to the json, but dont forget to update also in the the google sheet (https://docs.google.com/spreadsheets/d/1v71bnUd8fiXX51zuZIUAvYETWmpwFQj-M3mu4CNsHBU/edit#gid=791995663). 

_**data/cellosaurus_56.0_cell_lines.tsv.gz**_<br>
Cell lines for CellLine / CellLineIndex: a compact table (accession, recommended name, names and synonyms, species of origin, category, BTO cross-references, problematic-cell-line note) built from Cellosaurus release 56.0 (25 June 2026) with bin/build_cellosaurus_table.py. Expasy only serves the current Cellosaurus release, so the table is pinned and shipped in data/ rather than downloaded; data/cellosaurus_56.0.json records the release, the sha256 of the source cellosaurus.txt and of the table. To move to a newer release, download https://ftp.expasy.org/databases/cellosaurus/cellosaurus.txt, run `python bin/build_cellosaurus_table.py cellosaurus.txt data/ --ncbi_nodes_dmp nodes.dmp`, test, and change CELLOSAURUS_TABLE in nf_workflow.nf.
Cellosaurus is made available by the SIB Swiss Institute of Bioinformatics under the Creative Commons Attribution 4.0 International license (CC BY 4.0, https://creativecommons.org/licenses/by/4.0/); https://www.cellosaurus.org/. Bairoch A. The Cellosaurus, a cell-line knowledge resource. J Biomol Tech. 2018;29(2):25-38. The shipped table is a reduced subset of the release (changes listed in data/cellosaurus_56.0.json).
