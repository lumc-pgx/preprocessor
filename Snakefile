# load config file
configfile: srcdir("config.yaml")

# imports
import os
import glob
from Bio import SeqIO
import yaml
import datetime

# globals
PLATFORM = config["SEQUENCING_PLATFORM"].upper()
assert PLATFORM in ("RS2", "SEQUEL"), "Invalid sequencing platform specified"

if PLATFORM == "RS2":
    BASECALLS = [sorted(glob.glob(pth + "/*.bax.h5")) for pth in config["SOURCE_DATA_PATHS"]]
else:
    BASECALLS = [sorted(glob.glob(pth + "/*.bam")) for pth in config["SOURCE_DATA_PATHS"]]

MOVIES = [os.path.basename(BASECALLS[i][0]).split(".")[0] for i in range(len(BASECALLS))]
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
    print("Preprocessing workflow completed successfully")
onerror:
    print("Error encountered while executing workflow")
    shell("cat {log}")


# main workflow
localrules:
    all, fastq_to_fasta

rule all:
    input:
        expand("ccs_check/{barcodes}/variants.csv", barcodes=BARCODE_IDS), # run CCS and ccs_check
        expand("LAA/{barcodes}.fasta", barcodes=BARCODE_IDS),
        "summary/LAA/laa_summary.csv",
        "config.{}.yaml".format("{:%Y-%m-%d_%H:%M:%S}".format(datetime.datetime.now()))


rule source_data:
    """
    Populate basecalls folder with subreads bam file
    Convert RS2 format bax.h5 to sequel bam if required
    """
    input:
        lambda wildcards: BASECALLS[MOVIES.index(wildcards.moviename)]
    output:
        "basecalls/{moviename}.subreads.bam",
        "basecalls/{moviename}.scraps.bam"
    params:
        bax2bam = config["SMRTCMD_PATH"] + "/bax2bam"
    run:
        shell("mkdir -p basecalls")
        if PLATFORM == "RS2":
            shell("cd basecalls && {params.bax2bam} -o {wildcards.moviename} {input}")
        else:
            subreads = next(f for f in input if f.endswith("subreads.bam"))
            scraps = next(f for f in input if f.endswith("scraps.bam"))
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
        "barcoded_subreads/{moviename}.subreadset.xml",
        "barcoded_subreads/{moviename}.report.json.gz"
    params:
        bam2bam = config["SMRTCMD_PATH"] + "/bam2bam"
    shell:
        "{params.bam2bam} --barcodes {input.barcodes} -o barcoded_subreads/{wildcards.moviename} {input.subreads} {input.scraps}"


rule barcoding_summary:
    """
    Put the barcode summary report into the summary folder
    """
    input:
        "barcoded_subreads/{moviename}.report.json.gz"
    output:
        "summary/barcoding/{moviename}.report.json.gz"
    shell:
        "mkdir -p summary/barcoding && cp {input} {output}"

rule merge:
    """
    Merge barcoded subreads from multiple movies into a single dataset
    """
    input:
        expand("barcoded_subreads/{moviename}.subreadset.xml", moviename=MOVIES)
    output:
        "merged_subreads/merged.subreadset.xml"
    params:
        dataset = config["SMRTCMD_PATH"] + "/dataset"
    shell:
        "{params.dataset} merge {output} {input}"


rule demultiplex:
    """
    Create a single dataset for each barcode
    """
    input:
        "merged_subreads/merged.subreadset.xml"
    output:
        expand("demultiplexed/{bc_id}.xml", bc_id=BARCODE_IDS)
    params:
        dataset = config["SMRTCMD_PATH"] + "/dataset"
    run:
        for i, bc in enumerate(BARCODE_IDS):
            outfile = output[i]
            shell("{params.dataset} filter {input} {outfile} 'bc=[{i},{i}]'")


rule consolidate:
    """
    Create bam file from dataset xml
    """
    input:
        "demultiplexed/{barcode}.xml"
    output:
        bam = "consolidated/{barcode}.bam",
        xml = "consolidated/{barcode}.xml"
    params:
        dataset = config["SMRTCMD_PATH"] + "/dataset"
    shell:
        "{params.dataset} consolidate {input} {output.bam} {output.xml}"


rule laa:
    """
    Run LAA
    """
    input:
        bam = "consolidated/{barcode}.bam",
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
    run:
        shell(
            "{params.laa} {params.laa_params} "
            "--barcodes {input.barcodes} "
            "--resultFile {output.results} "
            "--junkFile {output.noise} "
            "--reportFile {output.report} "
            "--inputReportFile {output.input_report} "
            "--subreadsReportPrefix {params.subread_prefix} "
            "{input.bam}")
        # LAA does not create output files if no amplicons were constructed
        # Create empty files here for those cases
        for output_file in output:
            if not os.path.isfile(output_file):
                shell("touch {output_file}")
        

rule laa_summary:
    """
    Summarize LAA results
    """
    input:
        expand("LAA/{barcodes}_summary.csv", barcodes=BARCODE_IDS)
    output:
        "summary/LAA/laa_summary.csv"
    shell:
        "mkdir -p summary/LAA && head -n1 {input[0]} > {output} && "
        "cat {input} | grep -v 'BarcodeName,FastaName' >> {output}"


rule ccs:
    """
    Run CCS
    """
    input:
        "consolidated/{barcode}.bam"
    output:
        bam = "CCS/{barcode}.bam",
        report = "CCS/{barcode}.report.csv"
    params:
        ccs = config["SMRTCMD_PATH"] + "/ccs",
        ccs_params = tool_param_string(config["CCS_PARAMS"])
    shell:
        "{params.ccs} {params.ccs_params} --reportFile {output.report} "
        "{input} {output.bam}"


rule ccs_check:
    """
    Run ccs_check to analyze the ccs results
    """
    input:
        bam = "CCS/{barcode}.bam",
        genome = config["GENOME"]
    output:
        "ccs_check/{barcode}/qv_calibration.csv",
        "ccs_check/{barcode}/snrs.csv",
        "ccs_check/{barcode}/variants.csv",
        "ccs_check/{barcode}/zmws.csv",
        "ccs_check/{barcode}/zscores.csv"
    params:
        ccs_check_dir = config["CCS_CHECK_PATH"],
        ccs_check = config["CCS_CHECK_PATH"] + "/ccscheck",
        bam_path = os.path.abspath("CCS/{barcode}.bam"),
        genome_path = os.path.abspath(config["GENOME"])
    shell:
        "mkdir -p ccs_check && cd ccs_check && rm -rf {wildcards.barcode} && "
        "export LD_LIBRARY_PATH={params.ccs_check_dir}:$LD_LIBRARY_PATH && "
        "{params.ccs_check} {params.bam_path} {wildcards.barcode} {params.genome_path}"


rule fastq_to_fasta:
    """
    Convert fastq to fasta
    """
    input:
        "{source}.fastq"
    output:
        "{source}.fasta"
    run:
        with open(input[0], "r") as fastq, open(output[0], "w") as fasta:
            SeqIO.write((rec for rec in SeqIO.parse(fastq, "fastq")), fasta, "fasta")


rule write_config:
    output:
        "config.{timestamp}.yaml"
    run:
        with open(output[0], "w") as outfile:
            yaml.dump(config, outfile, default_flow_style=False)
