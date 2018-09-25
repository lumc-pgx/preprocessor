config = Snakemake.config

if config["SEQUENCING_PLATFORM"] == "RS2":
    shell("cd preprocessor/basecalls && bax2bam -o {wildcards.moviename} {input}")
else:
    subreads = next(f for f in input if f.endswith("subreads.bam"))
    scraps = next(f for f in input if f.endswith("scraps.bam"))
    shell("ln -s {subreads} {output[0]}")
    shell("ln -s {scraps}   {output[1]}")
