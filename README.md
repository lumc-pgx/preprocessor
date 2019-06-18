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
	graph [bb="0,0,93,468",
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
	0	 [color="0.19 0.6 0.85",
		height=0.5,
		label=source_data,
		pos="38,450",
		width=1.0556];
	2	 [color="0.00 0.6 0.85",
		height=0.5,
		label=barcoding,
		pos="38,378",
		width=0.88889];
	0 -> 2	 [pos="e,38,396.1 38,431.7 38,423.98 38,414.71 38,406.11"];
	1	 [color="0.57 0.6 0.85",
		height=0.5,
		label=all,
		pos="38,18",
		width=0.75];
	4	 [color="0.38 0.6 0.85",
		height=0.5,
		label=merge,
		pos="38,306",
		width=0.75];
	2 -> 4	 [pos="e,38,324.1 38,359.7 38,351.98 38,342.71 38,334.11"];
	3	 [color="0.29 0.6 0.85",
		height=0.5,
		label=consolidate,
		pos="38,162",
		width=0.97917];
	3 -> 1	 [pos="e,34.752,36.112 34.752,143.89 32.954,133.54 30.905,120.06 30,108 28.803,92.045 28.803,87.955 30,72 30.637,63.518 31.838,54.336 33.121,\
46.036"];
	5	 [color="0.48 0.6 0.85",
		height=0.5,
		label=ccs,
		pos="66,90",
		width=0.75];
	3 -> 5	 [pos="e,59.158,108.1 44.921,143.7 48.075,135.81 51.878,126.3 55.381,117.55"];
	6	 [color="0.10 0.6 0.85",
		height=0.5,
		label=demultiplex,
		pos="38,234",
		width=1];
	4 -> 6	 [pos="e,38,252.1 38,287.7 38,279.98 38,270.71 38,262.11"];
	5 -> 1	 [pos="e,44.842,36.104 59.079,71.697 55.925,63.813 52.122,54.304 48.619,45.546"];
	6 -> 3	 [pos="e,38,180.1 38,215.7 38,207.98 38,198.71 38,190.11"];
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
          
