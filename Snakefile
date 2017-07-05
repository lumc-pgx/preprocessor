# load config file
configfile: srcdir("config.yaml")

import os
import glob
from Bio import SeqIO

RS2_BASECALLS = glob.glob(config["SOURCE_DATA_PATH"] + "/*.bax.h5")
SEQUEL_BASECALLS = glob.glob(config["SOURCE_DATA_PATH"] + "/*.bam")
MOVIE = os.path.basename((RS2_BASECALLS + SEQUEL_BASECALLS)[0]).split(".")[0]
BARCODE_IDS = [x.id for x in SeqIO.parse(config["BARCODES"], "fasta")]


def tool_param_string(config_key):
    """
    Return a string containing tool parameters from the config,
    formatted for feeding to a command line tool
    """
    return " ".join(["--"+" ".join([str(x) for x in list(i.items())[0] if x]) for i in config[config_key]])


# handlers for workflow exit status
onsuccess:
    print("Pharmacogenomics preprocessing workflow completed successfully")
onerror:
    print("Error encountered while executing workflow")
    shell("cat {log}")


rule all:
    input:
        expand("LAA/{barcodes}.fastq", barcodes=BARCODE_IDS)


rule source_data:
    """
    Populate basecalls folder with subreads bam file
    Convert RS2 format bax.h5 to sequel bam if required
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
    Create barcoded bam file from subreads bam filea
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


rule demultiplex:
    """
    Create a single bam file for each barcode
    """
    input:
        "barcoded_subreads/{}.subreadset.xml".format(MOVIE)
    output:
        expand("demultiplexed/{bc_id}.xml", bc_id=BARCODE_IDS)
    run:
        shell("mkdir -p demultiplexed")
        for i, bc in enumerate(BARCODE_IDS):
            outfile = output[i]
            shell("dataset filter {input} {outfile} 'bc=[{i},{i}]'")


rule consolidate:
    """
    Create bam file from dataset xml
    """
    input:
        "{stage}/{barcode}.xml"
    output:
        "{stage}/{barcode}.bam"
    shell:
        "dataset consolidate {input} {output} {input}"


rule laa:
    """
    Run LAA
    """
    input:
        bam = "demultiplexed/{barcode}.bam",
        barcodes = config["BARCODES"]
    output:
        results = "LAA/{barcode}.fastq",
        noise = "LAA/{barcode}_chimeras_noise.fastq",
        report = "LAA/{barcode}_summary.csv",
        input_report = "LAA/{barcode}_input.csv",
        subread_report = "LAA/subreads.{barcode}--{barcode}.csv"
    params:
        subread_prefix = "LAA/subreads",
        laa_params = tool_param_string("LAA_PARAMS")
    shell:
        "echo {params.laa_params} && "
        "mkdir -p LAA && "
        "laa {params.laa_params} "
        "--barcodes {input.barcodes} "
        "--resultFile {output.results} "
        "--junkFile {output.noise} "
        "--reportFile {output.report} "
        "--inputReportFile {output.input_report} "
        "--subreadsReportPrefix {params.subread_prefix} "
        "{input.bam}"
        
