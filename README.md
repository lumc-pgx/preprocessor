# Preprocessor

**Initial data processing for PacBio amplicon data.**  

The preprocessor module performs the following operations:  
- Combine raw data from multiple SMRT cells associated with a single experiment.
- Barcoding
- Demultiplexing
- Long amplicon analysis (LAA2)
- Circular consensus sequence (CCS2)
- Comparison of LAA and CCS features  

```plantuml
digraph snakemake_dag {
    graph[bgcolor=white, margin=0];
    node[shape=box, style=rounded, fontname=sans,                 fontsize=10, penwidth=2];
    edge[penwidth=2, color=grey];
        0[label = "barcoding", color = "0.00 0.6 0.85", style="rounded"];
        1[label = "link_laa", color = "0.08 0.6 0.85", style="rounded"];
        2[label = "ccs_check", color = "0.12 0.6 0.85", style="rounded"];
        3[label = "all", color = "0.58 0.6 0.85", style="rounded"];
        4[label = "source_data", color = "0.25 0.6 0.85", style="rounded"];
        5[label = "laa_1", color = "0.54 0.6 0.85", style="rounded"];
        6[label = "fastq_to_fasta", color = "0.38 0.6 0.85", style="rounded"];
        7[label = "demultiplex", color = "0.33 0.6 0.85", style="rounded"];
        8[label = "ccs", color = "0.42 0.6 0.85", style="rounded"];
        9[label = "ccs_check_summary", color = "0.46 0.6 0.85", style="rounded"];
        10[label = "laa_summary", color = "0.50 0.6 0.85", style="rounded"];
        11[label = "laa_whitelist", color = "0.17 0.6 0.85", style="rounded"];
        12[label = "merge", color = "0.29 0.6 0.85", style="rounded"];
        13[label = "laa_2", color = "0.21 0.6 0.85", style="rounded"];
        14[label = "consolidate", color = "0.62 0.6 0.85", style="rounded"];
        4 -> 0
        13 -> 1
        8 -> 2
        6 -> 3
        8 -> 3
        9 -> 3
        10 -> 3
        14 -> 5
        1 -> 6
        12 -> 7
        14 -> 8
        1 -> 9
        2 -> 9
        1 -> 10
        5 -> 11
        0 -> 12
        14 -> 13
        11 -> 13
        7 -> 14
}            
```
  
## Requirements
- [Conda/Miniconda](https://conda.io/miniconda.html)  
- [ccscheck](https://github.com/PacificBiosciences/ccscheck)

## Installation
- Clone the repository
  - `git clone https://git.lumc.nl/PharmacogenomicsPipe/preprocessor.git`

- Change to the preprocessor directory
  - `cd preprocessor`

- Create a conda environment for running the pipeline
  - `conda env create -n preprocessor -f environment.yaml`

- In order to use the pipeline on the cluster, update your .profile to use the drmaa library:
  - `echo "export DRMAA_LIBRARY_PATH=libdrmaa.so.1.0" >> ~/.profile`
  - `source ~/.profile`

## Configuration
Pipeline configuration settings can be altered by editing [config.yaml](config.yaml).  

## Execution
- Activate the conda environment
  - `source activate preprocessor`
- For parallel execution on the cluster
  - `pipe-runner`
- To specify that the pipeline should write output to a location other than the default:
  - `pipe-runner --directory path/to/output/directory`
          
