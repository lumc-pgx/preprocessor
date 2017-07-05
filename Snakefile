# load config file
configfile: srcdir("config.yaml")

import os
import glob

RS2_BASECALLS = glob.glob(config["SOURCE_DATA_PATH"] + "/*.bax.h5")
SEQUEL_BASECALLS = glob.glob(config["SOURCE_DATA_PATH"] + "/*.bam")


# handlers for workflow exit status
onsuccess:
    print("Pharmacogenomics preprocessing workflow completed successfully")
onerror:
    print("Error encountered while executing workflow")
    shell("cat {log}")


rule all:
    input:
        "barcoded_subreads/m150227_003039_42162_c100786542550000001823153608251501_s1_p0.subreads.bam"


rule source_data:
    """
    Populate basecalls folder with subreads bam file
    """
    input:
        sorted(RS2_BASECALLS) if len(RS2_BASECALLS) == 3 else SEQUEL_BASECALLS
    output:
        "basecalls/{moviename}.subreads.bam",
        "basecalls/{moviename}.scraps.bam"
    params:
        bax2bam = config["SMRTCMD_PATH"] + "/bax2bam"
    run:
        shell("mkdir -p basecalls")
        if len(RS2_BASECALLS) == 3:
            shell("cd basecalls && {params.bax2bam} -o {wildcards.moviename} {input}")
        else:
            subreads = next(f for f in SEQUEL_BASECALLS if f.endswith("subreads.bam"))
            scraps = next(f for f in SEQUEL_BASECALLS if f.endswith("scraps.bam"))
            shell("ln -s {subreads} {output[0]}")
            shell("ln -s {scraps}   {output[1]}")


rule barcoding:
    """
    Create barcoded bam file from subreads bam file
    """
    input:
        subreads = "basecalls/{moviename}.subreads.bam",
        scraps = "basecalls/{moviename}.scraps.bam",
        barcodes = config["BARCODES"]
    output:
        "barcoded_subreads/{moviename}.subreads.bam",
        "barcoded_subreads/{moviename}.scraps.bam",
        "barcoded_subreads/{moviename}.subreadset.xml"
    params:
        bam2bam = config["SMRTCMD_PATH"] + "/bam2bam"
    shell:
        """
        mkdir -p barcoded_subreads &&   
        {params.bam2bam} --barcodes {input.barcodes} -o barcoded_subreads/{wildcards.moviename} {input.subreads} {input.scraps}
        """

