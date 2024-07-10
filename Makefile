run:
	nextflow run ./nf_workflow.nf -resume -c nextflow.config -with-report report.html -with-trace trace.html

run_hpcc:
	nextflow run ./nf_workflow.nf -resume -c nextflow_hpcc.config

run_docker:
	nextflow run ./nf_workflow.nf -resume -with-docker <CONTAINER NAME>