# Preprocessor

Initial data processing for PacBio data.  
Generates amplicon sequences and associated metrics.

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

## Configuration
Pipeline configuration settings can be altered by editing [config.yaml](config.yaml).  

## Execution
- For parallel execution on the cluster, run `run_cluster.sh`
- For execution in a single thread, run `run_local.sh`
