# Preprocessor

*Initial data processing for PacBio amplicon data.*  

The preprocessor module performs the following operations:  
- Combine raw data from multiple SMRT cells associated with a single experiment.
- Barcoding
- Demultiplexing
- Long amplicon analysis (LAA2)
- Circular consensus sequence (CCS2)
- Comparison of LAA and CCS features  

## Requirements
- [Conda/Miniconda](https://conda.io/miniconda.html)  
- [PacBio SmrtLink](https://github.com/PacificBiosciences/SMRT-Link)
- [ccscheck](https://github.com/PacificBiosciences/ccscheck)

## Installation
- Clone the repository
  - `git clone https://wgallard@git.lumc.nl/PharmacogenomicsPipe/preprocessor.git`

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
  - `./run_cluster.sh`

