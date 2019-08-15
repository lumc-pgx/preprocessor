# Preprocessor

**Initial data processing for PacBio amplicon data.**  

The preprocessor module performs the following operations:  
- Combine raw data from multiple SMRT cells
- Transparently handle data from both RS2 and Sequel platforms
- Barcoding
- Demultiplexing
- CCS

![rule graph](static/rulegraph.png)

## Requirements
- [Conda/Miniconda](https://conda.io/miniconda.html)  

## Installation
- Clone the repository
  - `git clone https://github.com/lumc-pgx/preprocessor.git`

- Change to the preprocessor directory
  - `cd preprocessor`

- Create a conda environment for running the pipeline
  - `conda env create -n preprocessor -f environment.yaml`

## Configuration
Pipeline configuration settings can be altered by editing [config.yaml](config.yaml).  

## Execution
- Activate the conda environment
  - `source activate preprocessor`
- Run the pipeline using Snakemake
          
