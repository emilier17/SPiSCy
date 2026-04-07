<img src="images/banner.png" width="200"/>


## Get Started


Read the documentation below, or watch the [SPiSCy Tour video](placeholder) for a detailed walkthrough.

## Overview

SPiSCy (**S**nakemake **Pi**peline for **S**pectral **Cy**tometry) is a bioinformatic pipeline to fully analyze large spectral flow cytometry datasets. SPiSCy aims to be:
- Flexible, by having extensive custom configuration
- Reproducible, by having detailed logs
- Complete, by starting from raw FCS files and performing all steps for a full analysis
- Automatic, by using the workflow management system Snakemake


## Features

Starting with FCS files from the cytometer, SPiSCy performs:
- Data preprocessing (gating, transformation, QC, normalization)
- Dimensionality reduction (choice of 5 methods - PCA, KernelPCA, Isomap, FastICA, or direct markers)
- Clustering (choice of 6 methods - FlowSOM, PARC, PhenoGraph, BIRCH, CytoVI, or HDBSCAN)
- Comparaison between clustering results
- Differential analysis (abundance and marker expression)

All steps produce intermediate results and detailed log files that allow you to monitor the progress and success of each step.


SPiSCy can handle different experimental setups via its customizable configuration files:
- Large datasets, up to millions of cells and several markers
- Different types of files, like samples and controls


SPiSCy can be run on Windows or HPC infracstructure. If using HPC, jobs and resource requests are automatically handled.  


## Implementation

SPiSCy is written in R and Python and uses Snakemake to automate the workflow. Apptainer is used as the container system for HPC infrastructure. Docker is used for the Windows compatible version. 


Snakemake works by setting up rules that define how to create output files from input files. For more information about how Snakemake works, please see the official [documentation](https://snakemake.readthedocs.io/en/stable/). 


SPiSCy is organized as follows:

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

**Coding knowledge.** While SPiSCy is designed to be mostly automatic, you will still need to use the command line to install and run SPiSCy. You will also need to be able to modify yaml files to customize configuration settings. These steps are explained in the [Installation](#installation-and-setup) and [Usage](#usage) sections, as well as the [video](placeholder). 


**Biological knowledge.** Technical knowledge about cytometry and biological knowledge about the dataset is required to use SPiSCy. 

Preprocessing steps:
- Gating: select the cells of interest using typical bi-marker plots
- Transformation: know the typical distributions of the markers in your panel
- Batch correction: evaluate how batch correction affected marker distribution

Clustering steps:
- Clustering: test different dimensionality reduction and clustering parameters to find most ideal settings for your dataset
- Evaluation: annotate clusters based on markers, interpret the biological coherence of identified clusters


**Data and preprocessing quality.** Clustering and differential analysis quality is highly dependent on the quality of the raw data and the configuration choices made during the preprocessing steps. To ensure optimal cytometry quality, please follow the [best practices](https://www.thermofisher.com/ca/en/home/references/newsletters-and-journals/bioprobes-journal-of-cell-biology-applications/bioprobes-79/best-practices-multiparametric-flow-cytometry.html). To ensure preprocessing quality, please verify and validate the intermediate results produced by each preprocessing step. For more detail on how to monitor the pipeline, please see the [Monitor the pipeline](#monitor-the-pipeline-log-files-and-intermediate-results) section, or the [video](placeholder). 


**Tool assumptions.** Every method used in SPiSCy has its own set of assumptions about the data. If you want to make sure the methods are adapted to your data, please see each tool's documentation in [Tools and documentation](#tools-and-documentation). 



## Requirements

**Computing.** For large datasets (tens of millions of cells), an HPC infrastructure is ideal to run SPiSCy. The total number of cells is the most important variable for memory usage. Here is an example of resources used on HPC to analyze a dataset with 32 million cells and 17 markers (215 FCS files). 


Preprocessing
- Time (all steps): 2h30
- Highest memory usage: 32 GB
- CPUs: 4


Clustering
- Time varies per dimensionality reduction and clustering method. The execution times here include PCA dimensionality reduction, clustering and label propagagtion. 
    - FlowSOM: 10 min (clustering 13 million cells)
    - BIRCH: 14 min (clustering 13 million cells)
    - PARC: 25 min (clustering 950 000 cells)
    - PhenoGraph: 2h08 (clustering 955 000 cells)
    - CytoVI: 4h19 (clustering 855 000 cells)
    - HDBSCAN: 1h37 (clustering 855 000 cells)
- Highest memory usage: 64 GB
- CPUs: 8


Differential analysis
- Time: approx 20 min per clustering result
- Highest memory usage: 64 GB
- Number of CPUs: 8


**Differential analysis.** Requires a minimum amount of samples per condition tested. A warning will be emitted if there are too few samples for a particular statistical test. 



## Installation and Setup

### HPC

1. Clone the repository. 


Depending on your HPC infrastructure, you might be able to directly clone a GitHub repository with:

```
git clone https://github.com/emilier17/spiscy.git 
```

If direct cloning is not allowed, then clone the repository to a personal computer and upload the folder to the HPC using approved methods. 

2. Check software availability. 


Make sure R 4.4.0, Python 3.11, and Apptainer 1.4.5 are available on your HPC. This step depends on your HPC's infrastructure, so follow its documentation. Check if modules are already available:
```
module spider r
module spider python
module spider apptainer
```

If any module is missing, then follow your HPC installation instructions or request access. 


3. Create a Snakemake environment with the Snakemake slurm executor plugin.


Navigate to the spiscy folder. 

```
module load python/3.11
virtualenv --no-download ~/envs/snakemake  # or choose your own path
source ~/envs/snakemake/bin/activate
pip install --no-index snakemake==8.20.6
pip install snakemake-executor-plugin-slurm
```

4. Prepare the data folder.


Place all FCS files in ```data/all_raw``` folder.


Place the following CSV files in ```data``` folder: ```metadata.csv``` and ```marker_info.csv```


```metadata.csv``` must minimally contain a column called ```filename```, which lists all the FCS filenames present in ```data/all_raw```, without the .fcs extension. Example:

<img src="images/metadata_csv.png" width="600"/>


```marker_info.csv``` must contain columns ```marker_name```, ```channel```, and ```marker_class```. In ```marker_class```, indicate if the marker is used for distinguishing stable cell populations (type) or if the marker is used to characterize cells (state). Example:

<img src="images/marker_info_csv.png" width="250"/>


5. Adapt slurm profile and bash scripts. 


In ```run_spiscy.sh```, ```test_spiscy.sh```, and ```master_spiscy.sh```, paste the complete path to the spiscy folder in the indicated spaces. 


In ```run_spiscy.sh``` and ```test_spiscy.sh```, make sure the path to the snakemake environment is correct if you chose a different path in step 3. 


In ```master_spiscy.sh```, adjust SBATCH parameters: account name, requested time, memory, and CPUs.

> [!NOTE]
> Each child job cannot request more CPUs than master_spiscy.sh allocates, so set CPUs to the maximum required by any rule. Check [Requirements](#requirements) for estimated runtimes. 


In ```workflow/profiles/slurm/config.yaml```, indicate the account name. 

6. Download containers and convert to Apptainer.

```
cd apptainers
module load apptainer/1.4.5
apptainer pull spiscy_pybase.sif docker://ghcr.io/emilier17/spiscy_pybase:1.0
apptainer pull spiscy_rbase.sif docker://ghcr.io/emilier17/spiscy_rbase:1.0
```

This might take a few minutes. 

7. Optional cleanup


Delete folders and files for Windows usage.


You may delete the ```docker``` folder and the ```docker-compose.yaml``` file. 

8. Customize configuration files.


See [Customize configuration settings](#customize-configuration-settings) to adapt the configuration settings to your dataset. 


### Windows

1. Clone the repository.


In your terminal, change the current working directory to the location where you want spiscy to be and run:

```
git clone https://github.com/emilier17/spiscy.git 
```

A directory called spiscy will be created in the current working directory. Change the working directory to spiscy. All next steps take place in this working directory. 
```
cd spiscy
```

3. Build the docker images.


This will take several minutes. 
```
docker build -f docker/Dockerfile_rbase -t spiscy_rbase:1.0 .
docker build -f docker/Dockerfile_pybase -t spiscy_pybase:1.0 .
```

4. Prepare the data folder.


See instructions in step 4 in the HPC setup

5. Optional cleanup


Delete folders and files for HPC usage.


You may delete the folders ```apptainers``` and ```workflow/profiles```, as well as files ```master_spiscy.sh```, ```run_spiscy.sh```, and ```test_spiscy.sh```. 

6. Customize config files.


See [Customize configuration settings](#customize-configuration-settings) to adapt the configuration settings to your dataset. 


## Usage

### Customize configuration settings

Each rule has its own configuration file that should be customized. Each configuration file has extensive comments to help inform choices. Check out the [video](placeholder) for a detailed explanation of all configuration settings. 

> [!IMPORTANT]
> If using control files, make sure each control file has a unique identifying word in their filename. Then make sure to write that control ID in ```config/normalization.yaml```. This setting allows to distinguish between sample and control files. 



### Launch the pipeline

#### HPC

1. Navigate to the spiscy folder. 
2. Test the pipeline.


You might need to change permissions to run the file:
```
./test_spiscy.sh
```

3. Launch the pipeline. 


SPiSCy can either be run in a login node with:
```
./run_spiscy.sh
```

or submitted as a job with: 
```
sbatch master_spiscy.sh
```

Either way, SPiSCy will be able to automatically submit its own jobs to SLURM with automatic resource calculations. Running SPiSCy as a job itself allows SPiSCy to continue running if you are disconnected from the login node. 


If the pipeline failed and shut down, then you need to unlock the spiscy folder before relaunching.
1. Make sure spiscy is no longer running.
2. Navigate to the spiscy folder.
3. Activate the Snakemake environment (set in step 3 of installation).
```
source ~/envs/snakemake/bin/activate
```
4. Unlock the directory.
```
snakemake -n --unlock
```
5. Relaunch the pipeline.


#### Windows

Preprocessing and differential analysis steps use the rbase container, while clustering steps use the pybase container. According to which steps you want run, you will need to activate the appropriate container, run the chosen rules, then close the container, and move on to the next steps in the other container. 

1. Navigate to the spiscy folder.
2. Run the appropriate container:
```
docker compose run rbase
```

or

```
docker compose run pybase
```
3. Activate the appropriate conda environment:
```
conda activate rbase
```

or

```
conda activate pybase
```
4. Test the pipeline.
```
snakemake -np
```
5. Run all appropriate rules for the chosen container.


-j sets the number of cores. 

To run all preprocessing steps:
```
snakemake -j 1 -p --until merge_csv
``` 


To run all clustering steps:
```
snakemake -j 1 -p --until split_clusters_by_filename
``` 

To run all differential analysis:
```
snakemake -j 1 -p --until diff_analysis
``` 

6. Stop the container. 

```
exit
docker compose down --remove-orphans
```


### Monitor the pipeline: log files and intermediate results

Each rule produces 2 log files : stderr (errors and warnings) and sdtout (general information). These files are available under results/logs/. Addtionnally, if using HPC, SLURM specific logs are available under results/logs/slurm/. 

Each rule will also produce intermediate results, such as plots, heatmaps and csv files. To make sure every rule worked correctly, please check the corresponding rule folder in the results folder. Check out the [video](placeholder) for more details about interpretating the results.


## Tools and Documentation

This is a list of the main packages and their documentations. 


Preprocessing (R 4.4.0)
- [cytonorm](https://github.com/saeyslab/CytoNorm ) 2.0.8
- [flowcore](https://rdrr.io/bioc/flowCore/) 2.18.0
- [flowstats](https://rdrr.io/github/RGLab/flowStats/) 4.18.0
- [flowvs](https://rdrr.io/bioc/flowVS/) 1.38.0
- [peacoqc](https://rdrr.io/bioc/PeacoQC/) 1.16.0
- [ggplot2](https://ggplot2.tidyverse.org/) 3.5.1
- [umap](https://github.com/tkonopka/umap) 0.2.10.0


Clustering (python 3.11.14)
- [flowsom](https://flowsom.readthedocs.io/en/latest/) 0.2.2
- [matplotlib](https://matplotlib.org/) 3.10.8
- [numpy](https://numpy.org/) 2.3.5
- [pandas](https://pandas.pydata.org/) 2.3.3
- [parc](https://parc.readthedocs.io/en/latest/index.html) 0.40
- [phenograph](https://github.com/dpeerlab/phenograph) 1.5.7
- scikit-learn 1.8.0
    - [hdbscan](https://scikit-learn.org/stable/modules/generated/sklearn.cluster.HDBSCAN.html)
    - [birch](https://scikit-learn.org/stable/modules/generated/sklearn.cluster.Birch.html)
- [seaborn](https://seaborn.pydata.org/) 0.13.2
- scvi-tools 1.4.1
    - [cytovi](https://docs.scvi-tools.org/en/1.4.1/user_guide/models/cytovi.html)
- [umap-learn](https://umap-learn.readthedocs.io/en/latest/) 0.5.11


Differential analysis (R 4.4.0)
- [diffcyt](https://rdrr.io/bioc/diffcyt/) 1.26.0
- [edger](https://bioconductor.org/packages/release/bioc/html/edgeR.html) 4.4.0
- [limma](https://bioconductor.org/packages/release/bioc/html/limma.html) 3.62.1



## Acknowledgments

SPiSCy logo: Éloïse Blais ([portfolio](https://www.artstation.com/esiole))
