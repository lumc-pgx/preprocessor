include: "helper.snake"

PARAMS = PreprocessingHelper(config, "preprocessing")

onsuccess: PARAMS.onsuccess()
onerror: PARAMS.onerror()

# main workflow
localrules:
    all, fastq_to_fasta

rule all:
    input:
        PARAMS.outputs
        
include: "rules/source_data.snake"
include: "rules/barcoding.snake"
include: "rules/merge_subreadset.snake"
include: "rules/demultiplex.snake"
include: "rules/consolidate_xml.snake"
include: "rules/laa.snake"
include: "rules/laa_summary.snake"
include: "rules/laa_whitelist.snake"
include: "rules/ccs.snake"
include: "rules/ccs_check.snake"
include: "rules/fastq_to_fasta.snake"
include: "rules/ccs_check_summary.snake"
