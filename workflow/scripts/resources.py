#=====================================================================#
# SPiSCy (Snakemake PIpeline for Spectral CYtometry)
#
# Functions for dynamically calculating memory and time for SLURM. 
# Based on number and size of input fcs files
# 
# Author: Émilie Roy
# Date: Sept 2025
# Version: v.1.0
#======================================================================#

import os
import yaml
import subprocess
import glob

#### Locate necessary config files and load them ####

ROOT_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
CSV_DIR = os.path.join(ROOT_DIR, "results", "csv", "samples_final")
PREDICT_PATH = os.path.join(ROOT_DIR, "config", "predict_cofactors.yaml")
DETECT_BEF_PATH = os.path.join(ROOT_DIR, "config", "detect_batch_effect_before.yaml")
DETECT_AFT_PATH = os.path.join(ROOT_DIR, "config", "detect_batch_effect_after.yaml")
BIRCH_PATH = os.path.join(ROOT_DIR, "config", "run_birch.yaml")
CYTOVI_PATH = os.path.join(ROOT_DIR, "config", "run_cytovi.yaml")
FLOWSOM_PATH = os.path.join(ROOT_DIR, "config", "run_flowsom.yaml")
HDBSCAN_PATH = os.path.join(ROOT_DIR, "config", "run_hdbscan.yaml")
PARC_PATH = os.path.join(ROOT_DIR, "config", "run_parc.yaml")
PHENOGRAPH_PATH = os.path.join(ROOT_DIR, "config", "run_phenograph.yaml")
EVAL_CLUSTER_PATH = os.path.join(ROOT_DIR, "config", "evaluate_clustering.yaml")

with open(PREDICT_PATH) as f:
    PREDICT_CONFIG = yaml.safe_load(f)

with open(DETECT_BEF_PATH) as f:
    DETECT_BEF_CONFIG = yaml.safe_load(f)

with open(DETECT_AFT_PATH) as f:
    DETECT_AFT_CONFIG = yaml.safe_load(f)

with open(BIRCH_PATH) as f:
    BIRCH_CONFIG = yaml.safe_load(f)

with open(CYTOVI_PATH) as f:
    CYTOVI_CONFIG = yaml.safe_load(f)

with open(FLOWSOM_PATH) as f:
    FLOWSOM_CONFIG = yaml.safe_load(f)

with open(HDBSCAN_PATH) as f:
    HDBSCAN_CONFIG = yaml.safe_load(f)

with open(PARC_PATH) as f:
    PARC_CONFIG = yaml.safe_load(f)

with open(PHENOGRAPH_PATH) as f:
    PHENOGRAPH_CONFIG = yaml.safe_load(f)

with open(EVAL_CLUSTER_PATH) as f:
    EVAL_CLUSTER_CONFIG = yaml.safe_load(f)


#### Memory and time requests: preprocessing rules ####

def file_size_mb(file_path):
    """
    Return file size in mb
    """
    return os.path.getsize(file_path) / 1024**2

def mem_per_fcs(wildcards, input, threads=1, attempt=1):
    """
    Returns memory in MB for a job with 1 fcs file as input.
    5x the actual fcs file size, or minimum 3GB
    """
    fcs_size = file_size_mb(input[0])
    mem = max(5 * fcs_size, 3000)
    return int(mem * attempt)

def min_per_mb(wildcards, input, threads=1, attempt=1):
    """
    Returns runtime in minutes for a job with 1 file as input
    1 minute per 10 MB, minimum 5 minutes
    """
    file_size = file_size_mb(input[0])
    mins = max(file_size / 10, 5)
    return int(mins * attempt)

def mem_for_flowset(wildcards, input, threads=1, attempt=1):
    """
    Returns memory in MB for a job with many fcs file as input (flowset)
    9.5x sum size of fcs file, or minimum 32GB
    """
    total_mb = sum(file_size_mb(f) for f in input)
    mem = max(9.5 * total_mb, 32000)
    return int(mem * attempt)

def mem_for_evaluate_normalization(wildcards, input, threads=1, attempt=1):
    """
    Returns memory in MB for rule evaluate_normalization
    2.5x sum size of fcs file, or minimum 16GB
    """
    total_mb = sum(file_size_mb(f) for f in input)
    mem = max(2.5 * total_mb, 16000)
    return int(mem * attempt)


def time_predict_cofactors(wildcards, input, threads=1, attempt=1):
    """
    Returns runtime in minutes for predict_cofactors rule
    Depends:
        number of markers in predict_cofactors.yaml config (10 min per marker)
        sampling size for fcs file in predict_cofactors.yaml (extra time)
    """
    nb_markers = len(PREDICT_CONFIG["markers_to_transform"])
    nb_files = len(input)
    agg_size = int(PREDICT_CONFIG["downsample_size"]) * nb_files
    time_for_agg = agg_size / 8000 
    total_time = (nb_markers * 10) + time_for_agg
    return (total_time)

def min_per_flowset(wildcards, input, threads=1, attempt=1):
    """
    Returns runtime in minutes for a job with many fcs files as input (flowset)
    15 sec per fcs file in flowframe
    """
    nb_files = len(input)
    runtime_min = nb_files * 0.25
    return int(runtime_min * attempt)

def time_detect_before(wildcards, input, threads=1, attempt=1):
    """
    Returns runtime in minutes for detect_batch_effect_before rule
    Standard min_per_flowset time + extra time according to number of sampled cells in detect_batch_effect_before.yaml
    """
    min_time = min_per_flowset(wildcards, input, threads, attempt)
    sample_size = int(DETECT_BEF_CONFIG["sample_size"])
    nb_files = len(input)
    total_cells = nb_files * sample_size
    extra_time = total_cells * 0.0001
    return (min_time + extra_time)

def time_detect_after(wildcards, input, threads=1, attempt=1):
    """
    Returns runtime in minutes for detect_batch_effect_after rule
    Standard min_per_flowset time + extra time according to number of sampled cells in detect_batch_effect_after.yaml
    """
    min_time = min_per_flowset(wildcards, input, threads, attempt)
    sample_size = int(DETECT_AFT_CONFIG["sample_size"])
    nb_files = len(input)
    total_cells = nb_files * sample_size
    extra_time = total_cells * 0.0001
    return (min_time + extra_time)




#### Memory and time requests: clustering rules ####

def count_sample_csv_files():
    """
    Returns number of csv files in results/csv/samples_final
    Useful for calculating runtime of clustering algos
    """
    if not os.path.exists(CSV_DIR):
        raise FileNotFoundError(
            f"Directory not found: {CSV_DIR}"
        )

    csv_files = glob.glob(os.path.join(CSV_DIR, "*.csv"))

    if not csv_files:
        raise ValueError(
            f"No CSV files found in {CSV_DIR}"
        )

    return len(csv_files)

def count_rows_csv(filepath):
    """
    Returns the number of rows in a csv file, minus header row
    """
    result = subprocess.run(
        ["wc", "-l", filepath],
        capture_output=True,
        text=True,
        check=True)
    total_lines = int(result.stdout.strip().split()[0])
    return max(total_lines - 1, 0)

def count_rows_dir(directory):
    """
    Returns sums of all rows of csv files in a directory
    """
    total = 0
    for f in glob.glob(f"{directory}/*.csv"):
        total += count_rows_csv(f)
    return total

def count_columns_csv(filepath):
    """
    Returns number of columns by reading only first line of csv
    """
    with open(filepath, "r") as f:
        header = f.readline().strip()
    return len(header.split(","))

def mem_for_csv(wildcards, input, threads=1, attempt=1):
    """
    Returns memory in MB for a job with a big csv as input (ex final_samples.csv)
    2x sum size of the csv file, or minimum 48GB
    """
    total_mb = sum(file_size_mb(f) for f in input)
    mem = max(2 * total_mb, 48000)
    return int(mem * attempt)

def mem_for_diff_analysis(wildcards, input, threads=1, attempt=1):
    """
    Returns memory in MB for differential_analysis rule
    2.5x sum size of the csv file, or minimum 64GB
    """
    total_mb = sum(file_size_mb(f) for f in input)
    mem = max(2.5 * total_mb, 64000)
    return int(mem * attempt)

def mins_birch(wildcards, input, threads=1, attempt=1):
    """
    Returns runtime in minutes for running birch clustering
    Depends:
        sample size per csv file set in run_birch.yaml (direct clustering)
        number of cells which need cluster label propagation
    """
    sample_size = int(BIRCH_CONFIG["sample_size_per_file"])
    nb_files = count_sample_csv_files()
    nb_cells_to_cluster = sample_size * nb_files
    total_nb_cells = count_rows_csv(input["csv"])
    nb_cells_to_label = total_nb_cells - nb_cells_to_cluster

    cluster_unit = 1_000_000
    label_unit = 1_000_000

    mins_per_cluster_unit = 0.5
    mins_per_label_unit = 5

    clustering_time = (nb_cells_to_cluster/cluster_unit) * mins_per_cluster_unit
    labeling_time = (nb_cells_to_label / label_unit) * mins_per_label_unit

    total_time = clustering_time + labeling_time
    return (total_time)

def mins_cytovi(wildcards, input, threads=1, attempt=1):
    """
    Returns runtime in minutes for running cytovi clustering
    Depends:
        sample size per csv file set in run_cytovi.yaml (direct clustering)
        number of cells which need cluster label propagation
    """
    sample_size = int(CYTOVI_CONFIG["sample_size_per_file"])
    nb_files = count_sample_csv_files()
    nb_cells_to_cluster = sample_size * nb_files
    total_nb_cells = count_rows_csv(input["csv"])
    nb_cells_to_label = total_nb_cells - nb_cells_to_cluster

    cluster_unit = 100_000
    label_unit = 1_000_000

    mins_per_cluster_unit = 8
    mins_per_label_unit = 6

    clustering_time = (nb_cells_to_cluster/cluster_unit) * mins_per_cluster_unit
    labeling_time = (nb_cells_to_label / label_unit) * mins_per_label_unit

    total_time = clustering_time + labeling_time
    return (total_time)

def mins_flowsom(wildcards, input, threads=1, attempt=1):
    """
    Returns runtime in minutes for running flowsom clustering
    Depends:
        sample size per csv file set in run_flowsom.yaml (direct clustering)
        number of cells which need cluster label propagation
    """
    sample_size = int(FLOWSOM_CONFIG["sample_size_per_file"])
    nb_files = count_sample_csv_files()
    nb_cells_to_cluster = sample_size * nb_files
    total_nb_cells = count_rows_csv(input["csv"])
    nb_cells_to_label = total_nb_cells - nb_cells_to_cluster

    cluster_unit = 1_000_000
    label_unit = 1_000_000

    mins_per_cluster_unit = 1
    mins_per_label_unit = 0.5

    clustering_time = (nb_cells_to_cluster/cluster_unit) * mins_per_cluster_unit
    labeling_time = (nb_cells_to_label / label_unit) * mins_per_label_unit

    total_time = clustering_time + labeling_time
    return (total_time)

def mins_hdbscan(wildcards, input, threads=1, attempt=1):
    """
    Returns runtime in minutes for running HDBSCAN clustering
    Depends:
        sample size per csv file set in run_hdbscan.yaml (direct clustering)
        number of cells which need cluster label propagation
    """
    sample_size = int(HDBSCAN_CONFIG["sample_size_per_file"])
    nb_files = count_sample_csv_files()
    nb_cells_to_cluster = sample_size * nb_files
    total_nb_cells = count_rows_csv(input["csv"])
    nb_cells_to_label = total_nb_cells - nb_cells_to_cluster

    cluster_unit = 100_000
    label_unit = 1_000_000

    mins_per_cluster_unit = 8
    mins_per_label_unit = 5

    clustering_time = (nb_cells_to_cluster/cluster_unit) * mins_per_cluster_unit
    labeling_time = (nb_cells_to_label / label_unit) * mins_per_label_unit

    total_time = clustering_time + labeling_time
    return (total_time)

def mins_parc(wildcards, input, threads=1, attempt=1):
    """
    Returns runtime in minutes for running PARC clustering
    Depends:
        sample size per csv file set in run_parc.yaml (direct clustering)
        number of cells which need cluster label propagation
    """
    sample_size = int(PARC_CONFIG["sample_size_per_file"])
    nb_files = count_sample_csv_files()
    nb_cells_to_cluster = sample_size * nb_files
    total_nb_cells = count_rows_csv(input["csv"])
    nb_cells_to_label = total_nb_cells - nb_cells_to_cluster

    cluster_unit = 100_000
    label_unit = 1_000_000

    mins_per_cluster_unit = 3
    mins_per_label_unit = 5

    clustering_time = (nb_cells_to_cluster/cluster_unit) * mins_per_cluster_unit
    labeling_time = (nb_cells_to_label / label_unit) * mins_per_label_unit

    total_time = clustering_time + labeling_time
    return (total_time)

def mins_phenograph(wildcards, input, threads=1, attempt=1):
    """
    Returns runtime in minutes for running Phenograph clustering
    Depends:
        number of cells to directly cluster
        number of cells which need cluster label propagation
    """
    sample_size = int(PHENOGRAPH_CONFIG["sample_size_per_file"])
    nb_files = count_sample_csv_files()
    nb_cells_to_cluster = sample_size * nb_files
    total_nb_cells = count_rows_csv(input["csv"])
    nb_cells_to_label = total_nb_cells - nb_cells_to_cluster

    cluster_unit = 100_000
    label_unit = 1_000_000

    mins_per_cluster_unit = 9
    mins_per_label_unit = 6

    clustering_time = (nb_cells_to_cluster/cluster_unit) * mins_per_cluster_unit
    labeling_time = (nb_cells_to_label / label_unit) * mins_per_label_unit

    total_time = clustering_time + labeling_time
    return (total_time)

def mins_evaluate_clustering(wildcards, input, threads=1, attempt=1):
    """
    Returns runtime in minutes for evaluating clustering
    Depends:
        total number of cells (sum of all rows of csv file)
        number of cells to run umap on (in evaluate_clustering.yaml)
        number of markers (# columns in csv file - 1)
    """

    total_nb_cells = count_rows_csv(input["markers"])
    nb_markers = count_columns_csv(input["markers"]) - 1
    sample_size = int(EVAL_CLUSTER_CONFIG["sample_size_umap"])
    nb_files = count_sample_csv_files()
    nb_cells_to_umap = sample_size * nb_files

    cell_unit = 1_000_000
    umap_unit = 100_000

    mins_per_cell_unit = 1
    mins_per_umap_unit = 2

    marker_multiplier = nb_markers / 20

    cell_time = (total_nb_cells / cell_unit) * mins_per_cell_unit
    umap_time = (nb_cells_to_umap / umap_unit) * mins_per_umap_unit

    total_time = (cell_time * marker_multiplier) + umap_time

    return(total_time)


def mins_compare_clustering(wildcards, input, threads=1, attempt=1):
    """
    Returns runtime in minutes for comparing clustering results
    Depends:
        total number of cells (# rows in final_samples.csv - 1)
        number of markers (# columns in final_samples.csv - 1)
    """

    total_nb_cells = count_rows_csv(input["markers"])
    nb_markers = count_columns_csv(input["markers"]) - 1

    cell_unit = 1_000_000

    mins_per_cell_unit = 0.5
    marker_multiplier = nb_markers / 20

    cell_time = (total_nb_cells / cell_unit) * mins_per_cell_unit

    total_time = (cell_time * marker_multiplier)

    return(total_time)


def mins_split_samples(wildcards, input, threads=1, attempt=1):
    """
    Returns runtime in minutes for splitting final_samples into
    sample level individual csv files
    Depends:
        total number of cells (# rows in final_samples.csv - 1)
    """
    total_nb_cells = count_rows_csv(input["markers"])
    cell_unit = 1_000_000
    mins_per_cell_unit = 1
    total_time = (total_nb_cells / cell_unit) * mins_per_cell_unit
    return(total_time)

def mins_split_clusters(wildcards, input, threads=1, attempt=1):
    """
    Returns runtime in minutes for splitting final_samples into
    sample level individual csv files
    Depends:
        total number of cells (# rows in final_samples.csv - 1)
    """
    total_nb_cells = count_rows_csv(input["markers"])
    cell_unit = 1_000_000
    mins_per_cell_unit = 0.5
    total_time = (total_nb_cells / cell_unit) * mins_per_cell_unit
    return(total_time)


def mins_diff_analysis(wildcards, input, threads=1, attempt=1):
    """
    Returns runtime in minutes for differential analysis
    Depends:
        total number of cells (sum of all rows of csv files)
        number of cells to run umap on (in evaluate_clustering.yaml)
        number of markers (# columns in csv file - 1)
    """
    total_nb_cells = count_rows_dir(input["markers"])
    first_csv = sorted(glob.glob(os.path.join(input["markers"], "*.csv")))[0]
    nb_markers = count_columns_csv(first_csv) - 1
    sample_size = int(EVAL_CLUSTER_CONFIG["sample_size_umap"])
    nb_files = count_sample_csv_files()
    nb_cells_to_umap = sample_size * nb_files

    cell_unit = 1_000_000
    umap_unit = 100_000

    mins_per_cell_unit = 2
    mins_per_umap_unit = 2

    marker_multiplier = nb_markers / 20

    cell_time = (total_nb_cells / cell_unit) * mins_per_cell_unit
    umap_time = (nb_cells_to_umap / umap_unit) * mins_per_umap_unit

    total_time = (cell_time * marker_multiplier) + umap_time

    return(total_time)
