# load config file
configfile: srcdir("config.yaml")

# imports
import os
import glob
from Bio import SeqIO

# globals
RS2_BASECALLS = glob.glob(config["SOURCE_DATA_PATH"] + "/*.bax.h5")
SEQUEL_BASECALLS = glob.glob(config["SOURCE_DATA_PATH"] + "/*.bam")
MOVIE = os.path.basename((RS2_BASECALLS + SEQUEL_BASECALLS)[0]).split(".")[0]
BARCODE_IDS = [x.id for x in SeqIO.parse(config["BARCODES"], "fasta")]


# utility functions
def parse_param(config_entry):
    """
    Parse a config entry from the config file
    """
    try:
        return list(config_entry.items())[0]
    except AttributeError:
        return [config_entry]


def tool_param_string(config_dict):
    """
    Convert config settings to a parameter string for feeding to command line tools
    :param config_dict: portion of the config which contains the params
    :return: A string containing the formatted parameters
    """
    return " ".join(["--"+" ".join([str(x) for x in parse_param(i)]) for i in config_dict])


# handlers for workflow exit status
onsuccess:
    print("Pharmacogenomics preprocessing workflow completed successfully")
onerror:
    print("Error encountered while executing workflow")
    shell("cat {log}")


# main workflow
rule all:
    input:
        expand("LAA/{barcodes}.fastq", barcodes=BARCODE_IDS),
        expand("CCS/{barcodes}.bam", barcodes=BARCODE_IDS),
        #expand("ccs_check/{barcodes}/snrs.csv", barcodes=BARCODE_IDS)
        "ccs_check/PB_M13F_01/snrs.csv"


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
        laa = config["SMRTCMD_PATH"] + "/laa",
        subread_prefix = "LAA/subreads",
        laa_params = tool_param_string(config["LAA_PARAMS"])
    shell:
        "echo {params.laa_params} && "
        "mkdir -p LAA && "
        "{params.laa} {params.laa_params} "
        "--barcodes {input.barcodes} "
        "--resultFile {output.results} "
        "--junkFile {output.noise} "
        "--reportFile {output.report} "
        "--inputReportFile {output.input_report} "
        "--subreadsReportPrefix {params.subread_prefix} "
        "{input.bam}"


rule ccs:
    """
    Run CCS
    """
    input:
        "demultiplexed/{barcode}.bam"
    output:
        bam = "CCS/{barcode}.bam",
        report = "CCS/{barcode}.report.csv"
    params:
        ccs = config["SMRTCMD_PATH"] + "/ccs",
        ccs_params = tool_param_string(config["CCS_PARAMS"])
    shell:
        "echo {params.ccs_params} && "
        "mkdir -p CCS && "
        "{params.ccs} {params.ccs_params} --reportFile {output.report} "
        "{input} {output.bam}"


rule ccs_check:
    """
    Run ccs_check to analyze the ccs results
    """
    input:
        bam = os.path.abspath("CCS/{barcode}.bam"),
        genome = os.path.abspath(config["GENOME"])
    output:
        "ccs_check/{barcode}/qv_calibration.csv",
        "ccs_check/{barcode}/snrs.csv",
        "ccs_check/{barcode}/variants.csv",
        "ccs_check/{barcode}/zmws.csv",
        "ccs_check/{barcode}/zscores.csv"
    params:
        ccs_check_dir = config["CCS_CHECK_PATH"],
        ccs_check = config["CCS_CHECK_PATH"] + "/ccscheck"
    shell:
        "mkdir -p ccs_check && cd ccs_check && rm -rf {wildcards.barcode} && "
        "export LD_LIBRARY_PATH={params.ccs_check_dir}:$LD_LIBRARY_PATH && "
        "{params.ccs_check} {input.bam} {wildcards.barcode} {input.genome}"
