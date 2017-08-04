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
	0[label = "consolidate", color = "0.00 0.6 0.85", style="rounded"];
	1[label = "fastq_to_fasta", color = "0.31 0.6 0.85", style="rounded"];
	2[label = "ccs_check_summary", color = "0.05 0.6 0.85", style="rounded"];
	3[label = "demultiplex", color = "0.10 0.6 0.85", style="rounded"];
	4[label = "ccs_check", color = "0.26 0.6 0.85", style="rounded"];
	5[label = "ccs", color = "0.21 0.6 0.85", style="rounded"];
	6[label = "laa", color = "0.36 0.6 0.85", style="rounded"];
	7[label = "all", color = "0.41 0.6 0.85", style="rounded"];
	8[label = "laa_summary", color = "0.46 0.6 0.85", style="rounded"];
	9[label = "source_data", color = "0.15 0.6 0.85", style="rounded"];
	10[label = "barcoding", color = "0.51 0.6 0.85", style="rounded"];
	11[label = "merge", color = "0.62 0.6 0.85", style="rounded"];
	3 -> 0
	6 -> 1
	4 -> 2
	6 -> 2
	11 -> 3
	5 -> 4
	0 -> 5
	0 -> 6
	1 -> 7
	8 -> 7
	2 -> 7
	5 -> 7
	6 -> 8
	9 -> 10
	10 -> 11
}
```
  
## Requirements
- [Conda/Miniconda](https://conda.io/miniconda.html)  
- [PacBio SmrtLink](https://github.com/PacificBiosciences/SMRT-Link)
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
          
