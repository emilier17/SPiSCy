Thanks for checking out SPiSCy. Get starting by reading the documentation below, or watch the [SPiSCY Tour video](placeholder) for a detailed walkthrough. 


## Overview

SPiSCy (<u>S</u>nakemake <u>Pi</u>line for <u>S</u>pectral <u>Cy</u>tometry) is a bioinformatic pipeline to fully analyze large spectral flow cytometry datasets. SPiSCY aims to be:
- flexible, by having extensive custom configuration
- reproducible, by having detailed logs
- complete, by starting from raw FCS files and performing all steps for a full analysis
- automatic, by using the workflow management system Snakemake


## Features

Starting from FCS files straight from the cytometer, SPiSCY performs:
- extensive data preprocessing (gating, transformation, QC, normalization)
- dimensionality reduction (choice of 5 methods - PCA, KernelPCA, Isomap, FastICA, or direct markers)
- clustering (choice of 6 methods - FlowSOM, PARC, PhenoGraph, BIRCH, CytoVI, and HDBSCAN)
- comparaison between clustering results
- differential analysis (abundance and marker expression)

All steps produce intermediate results and detailed log files that allow you to monitor the progress and success of each step.


SPiSCy can handle different experimental setups via its customizable configuration files:
- Large datasets, up to millions of cells and several markers
- Different types of files, like samples and controls


SPiSCy can be run on Windows or HPC infracstructure. If using HPC, jobs and resource requests are automatically handled.  


## Implementation

SPiSCy is written in R and Python and uses Snakemake to automate the workflow. Apptainer is used as the container system for HPC infrastructure. Docker is used for the Windows compatible version. 


Snakemake works by setting up rules that define how to create output files from input files. For more information about how Snakemake works, please see the official [documentation](https://snakemake.readthedocs.io/en/stable/). 


SPiSCY is organized per official Snakemake recommendations:

```
.
└── spiscy
    ├── apptainers/  <-- containers for HPC environment
    ├── config/  <-- configuration files for custom user settings
    ├── data/  <-- input fcs files and metadata files
    ├── docker/  <-- containers for Windows environment
    ├── results/  <-- created by Snakemake as pipeline progresses through steps
    ├── workflow
    │   ├── profiles/  <-- SLURM configuration file for HPC environment
    │   ├── rules  <-- order of each step of the pipeline
    │   │   ├── clustering.smk
    │   │   ├── preprocessing.smk
    │   │   ├── differential_analysis.smk
    │   │   └── scripts/  <-- R or Python script for each step of the pipeline
    │   └── Snakefile  <-- master Snakemake control file
    ├── master_spiscy.sh  <-- bash script for HPC environment
    ├── run_spiscy.sh  <-- bash script for HPC environment
    └── test_spiscy.sh  <-- bash script for HPC environment
```

## Workflow

This is the SPiSCy rulegraph, or what happens in what order to go from raw FCS files to final results (all). 

![Rulegraph](images/rulegraph.jpg)


## Limits

**Coding knowledge.** While SPiSCy is deisnged to be mostly automatic, you will still need to use the command line to install and run SPiSCy. You will also need to be able to modify yaml files to customize configuration settings. These steps are explained in the [Go to section](#installation-and-setup) and [Go to section](#usage) sections, as well as the [video](placeholder). 


**Biological knowledge.** Technical knowledge about cytometry and biological knowledge about the dataset is required to use SPiSCy. 
Preprocessing steps:
- Gating: selec the cells of interest using typical bi-marker plots
- Transformation: know the typical distributions of the markers in your panel
- Batch correction: evaluate how batch correction affected marker distribution

Clustering steps:
- Clustering: test different dimensionality reduction and clustering parameters to find most ideal settings for your dataset
- Evaluation: annotate clusters based on markers, interpret the biological coherence of identified clusters


**Data and preprocessing quality.** Clustering and differential analysis quality is highly dependent on the quality of the raw data and the configuration choices made during the preprocessing steps. To ensure optimal cytometry quality, please follow the [best practices](https://www.thermofisher.com/ca/en/home/references/newsletters-and-journals/bioprobes-journal-of-cell-biology-applications/bioprobes-79/best-practices-multiparametric-flow-cytometry.html). To ensure preprocessing quality, please verify and validate the intermediate results produced by each preprocessing step. For more detail on how to monitor the pipeline, please see [Go to section](#monitor-the-pipeline-log-files-and-intermediate-results) sections, or the [video](placeholder). 


**Tool assumptions.** Every method used in SPiSCy has its own set of assumptions about the data. If you want to make sure the methods are adapted to your data, please see a link to each tool's documentation in [Go to section](#tools-and-documentation). 



## Requirements

**Computing**
For large datasets (tens of millions of cells), and HPC infrastructure is ideal to run SPiSCY. The total number of cells is the most important variable for memory usage. Here is an example of resources used on HPC to analyze a dataset with 32 million cells and 17 markers (215 FCS files). 


Preprocessing
- Time (all steps): 2h30
- Highest memory usage: 32 GB
- CPUs: 4


Clustering
- Time: depends on method and subsample size. 
    - FlowSOM: 10 min (clustering 13 million cells)
    - BIRCH: 14 min (clustering 13 million cells)
    - PARC
    - PhenoGraph
    - CytoVI
    - HDBSCAN: 1h37 (clustering 855 000 cells)
- Highest memory usage: 64 GB
- CPUs: 8


Differential analysis
- Time: approx 20 min per clustering result
- Highest memory usage: 64 GB
- Number of CPUs: 8


**Differential analysis**
Requires a minimum amount of samples per condition tested. A warning will be emitted if there are too few samples for a particular statistical test. 



## Installation and Setup

### HPC

### Windows



## Usage

### Customize configuration settings

Each rule has its own configuration file that should be customized. Each configuration file has extensive comments to help inform choices. The settings are spoken about in detail in the [video](placeholder). 

### Launch the pipeline

### Monitor the pipeline: log files and intermediate results

Each rule produces 2 log files : stderr (errors and warnings) and sdtout (general information). These files are available under results/logs/. Addtionnally, if using HPC, SLURM specific logs are available under results/logs/slurm/. 

Each rule will also produce intermediate results, such as plots, heatmaps and csv files. To make sure every rule worked correctly, please check the results folder. Interpretating these results is spoken about in detail in the [video](placeholder).


## Tools and Documentation
Preprocessing

Clustering

Differential analysis