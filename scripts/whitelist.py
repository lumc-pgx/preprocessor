from collections import defaultdict
from math import floor


def zmw_assignments(filename):
    """
    Parse the summary file from LAA and itentify which amplicon the reads associated with each ZMW map to

    :param filename: path to the LAA summary file to be parsed
    :return: A (default)dict where the keys are the ZMW and the values are a list of indices to the amplicon sequence
    """
    zmws = defaultdict(set)
    with open(filename, "r") as infile:
        for row in infile:
            if row.startswith("SubreadId"):
                continue

            row = [x.strip() for x in row.split(",")]
            zmw = "/".join(row[0].split("/")[:-1])

            assignments = [bool(int(round(floor(float(x) + 0.5)))) for x in row[1:]]

            try:
                which = assignments.index(True)
                zmws[zmw].add(which)
            except ValueError:
                pass

    return zmws


# parse the laa summary file
zmws = zmw_assignments(snakemake.input[0])

# write the whitelist which contains only subreads for zmws which contribute to a single amplicon sequence
with open(snakemake.output[0], "w") as outfile:
    for zmw, count in zmws.items():
        if len(count) <= 1:
            print(zmw, file=outfile)
