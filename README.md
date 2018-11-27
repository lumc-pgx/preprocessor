# Preprocessor

**Initial data processing for PacBio amplicon data.**  

The preprocessor module performs the following operations:  
- Combine raw data from multiple SMRT cells
- Transparently handle data from both RS2 and Sequel platforms
- Barcoding
- Demultiplexing
- CCS

```plantuml
digraph snakemake_dag {
        graph [bb="0,0,92,468",
                bgcolor=white,
                margin=0
        ];
        node [fontname=sans,
                fontsize=10,
                label="\N",
                penwidth=2,
                shape=box,
                style=rounded
        ];
        edge [color=grey,
                penwidth=2
        ];
        0        [color="0.38 0.6 0.85",
                height=0.5,
                label=merge,
                pos="54,306",
                width=0.75];
        5        [color="0.57 0.6 0.85",
                height=0.5,
                label=demultiplex,
                pos="54,234",
                width=1];
        0 -> 5   [pos="e,54,252.1 54,287.7 54,279.98 54,270.71 54,262.11"];
        1        [color="0.48 0.6 0.85",
                height=0.5,
                label=barcoding,
                pos="54,378",
                width=0.88889];
        1 -> 0   [pos="e,54,324.1 54,359.7 54,351.98 54,342.71 54,334.11"];
        2        [color="0.00 0.6 0.85",
                height=0.5,
                label=consolidate,
                pos="54,162",
                width=0.97917];
        3        [color="0.10 0.6 0.85",
                height=0.5,
                label=ccs,
                pos="27,90",
                width=0.75];
        2 -> 3   [pos="e,33.597,108.1 47.326,143.7 44.285,135.81 40.617,126.3 37.239,117.55"];
        4        [color="0.19 0.6 0.85",
                height=0.5,
                label=all,
                pos="54,18",
                width=0.75];
        2 -> 4   [pos="e,57.654,36.092 57.654,143.91 59.676,133.57 61.981,120.09 63,108 64.344,92.057 64.344,87.943 63,72 62.283,63.499 60.931,54.312 59.488,\
46.012"];
        3 -> 4   [pos="e,47.403,36.104 33.674,71.697 36.715,63.813 40.383,54.304 43.761,45.546"];
        5 -> 2   [pos="e,54,180.1 54,215.7 54,207.98 54,198.71 54,190.11"];
        6        [color="0.29 0.6 0.85",
                height=0.5,
                label=source_data,
                pos="54,450",
                width=1.0556];
        6 -> 1   [pos="e,54,396.1 54,431.7 54,423.98 54,414.71 54,406.11"];
}
```

## Requirements
- [Conda/Miniconda](https://conda.io/miniconda.html)  

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
          
