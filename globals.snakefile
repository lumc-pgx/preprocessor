import os
import glob
import yaml
import datetime

# yaml representer for dumping config
from yaml.representer import Representer
import collections
yaml.add_representer(collections.defaultdict, Representer.represent_dict)

PLATFORM = config.get("SEQUENCING_PLATFORM", "RS2").upper()
assert PLATFORM in ("RS2", "SEQUEL"), "Invalid sequencing platform specified"

if PLATFORM == "RS2":
    BASECALLS = [sorted(glob.glob(pth + "/*.bax.h5")) for pth in config["SOURCE_DATA_PATHS"]]
else:
    BASECALLS = [sorted(glob.glob(pth + "/*.bam")) for pth in config["SOURCE_DATA_PATHS"]]

MOVIES = [os.path.basename(BASECALLS[i][0]).split(".")[0] for i in range(len(BASECALLS))]
with open(config["BARCODES"], "r") as bc_file:
    BARCODE_IDS = [line.strip()[1:] for line in bc_file if line.startswith(">")]


# determine the files to be generated depending on which stages are to be run
TARGET_FILES = []
if config.get("STAGES", {}).get("LAA", False):
    TARGET_FILES.append("summary/LAA/laa_summary.csv")
    TARGET_FILES += expand("LAA/{barcodes}.fasta", barcodes=BARCODE_IDS)

if config.get("STAGES", {}).get("CCS", False):
    TARGET_FILES += expand("CCS/{barcodes}.bam", barcodes=BARCODE_IDS)

if config.get("STAGES", {}).get("CCS_CHECK", False):
    TARGET_FILES += expand("summary/ccs_check/{barcodes}.html", barcodes=BARCODE_IDS)



# utility functions
def parse_param(config_pair):
    """
    Parse a config entry from the config file
    """
    
    # convert boolean to flag
    if config_pair[1] is True:
        return [config_pair[0]]
    
    if config_pair[1] is False:
        return []
    
    return config_pair


def tool_param_string(config_dict):
    """
    Convert config settings to a parameter string for feeding to command line tools
    :param config_dict: portion of the config which contains the params
    :return: A string containing the formatted parameters
    """
    return " ".join(["--"+" ".join([str(x) for x in parse_param([i, config_dict[i]])]) for i in config_dict])


# handlers for workflow exit status
def onsuccess():
    print("{} workflow completed successfully".format(WORKFLOW_NAME))
    config_file = "config.{}.yaml".format("{:%Y-%m-%d_%H:%M:%S}".format(datetime.datetime.now()))
    with open(config_file, "w") as outfile:
        print(yaml.dump(config, default_flow_style=False), file=outfile)

def onerror():
    print("Error encountered while executing workflow")
    shell("cat {log}")

