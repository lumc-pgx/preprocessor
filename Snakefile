include: "helper.snake"

PARAMS = PreprocessingHelper(config, "preprocessing")

onsuccess: PARAMS.onsuccess()
onerror: PARAMS.onerror()

# main workflow

rule all:
    input:
        PARAMS.outputs
        
include: "rules/source_data.snake"
include: "rules/barcoding.snake"
include: "rules/merge_subreadset.snake"
include: "rules/demultiplex.snake"
include: "rules/consolidate_xml.snake"
include: "rules/ccs.snake"
