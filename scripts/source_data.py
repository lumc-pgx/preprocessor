from snakemake import shell

config = snakemake.config
wildcards = snakemake.wildcards
basecalls = snakemake.input
subreads = snakemake.output[0]
scraps = snakemake.output[1]

if config["SEQUENCING_PLATFORM"] == "RS2":
    shell("cd preprocessor/basecalls && bax2bam -o {wildcards.moviename} {basecalls}")
else:
    subreads_bam = next(f for f in input if f.endswith("subreads.bam"))
    scraps_bam = next(f for f in input if f.endswith("scraps.bam"))
    shell("ln -s {subreads_bam} {subreads}")
    shell("ln -s {scraps_bam} {scraps}")
