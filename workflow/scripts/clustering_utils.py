#==================================================================================================================#
# SPiSCy (Snakemake PIpeline for Spectral CYtometry)
#
# Phase: clustering
# Step: functions for clustering
#
# Author: Émilie Roy
# Date: Jan 2026
# Version: v.1.0
#
# Output: functions for clustering to import into individual clustering scripts
#==================================================================================================================#



#######################
#### GENERAL UTILS ####
#######################

import os
import sys
import atexit
from datetime import datetime


def snakemake_logs(snakemake_obj):
    """
    Redirect stdout and stderr to Snakemake log files
    """
    log_stdout = open(snakemake_obj.log.stdout, "w")
    log_stderr = open(snakemake_obj.log.stderr, "w")

    sys.stdout = log_stdout
    sys.stderr = log_stderr

    def _cleanup():
        sys.stdout.flush()
        sys.stderr.flush()
        log_stdout.close()
        log_stderr.close()

    atexit.register(_cleanup)

def start_time():
    """
    Get script start time
    Returns start_time as a datetime.now() format
    """
    start_time = datetime.now()
    print("Start:", start_time.strftime("%Y-%m-%d %H:%M:%S"))
    return start_time

def end_time():
    """
    Get script end time
    Returns end_time as a datetime.now() format
    """
    end_time = datetime.now()
    print("End:", end_time.strftime("%Y-%m-%d %H:%M:%S"))
    return end_time

def elapsed_time(start_time, end_time):
    """
    Calculate total script execution time. Prints execution time without miliseconds
    
    :param start_time: datetime format
    :param end_time: datetime format
    """
    elapsed = end_time - start_time
    elapsed_short = str(elapsed).split(".")[0]
    print(f"Total execution time: {elapsed_short}")

def print_pd_dataframe(dataframe, msg):
    """
    Print a preview of a pandas dataframe: columns and head
    
    :param dataframe: pandas dataframe
    :param msg: additional printing information about df
    """

    print(f"[TABLE PREVIEW] {msg}:")
    print(dataframe.head())
    print(f"Columns used: {dataframe.columns}")
    print("\n")

def print_save_msg(filepath, filename):
    """
    Print save message for a file
    
    :param filepath: filepath to where file was saved
    :str filename: name of what was saved
    """
    print(f">> Saved {filename} to {filepath} \n")


def check_pos_int_not_null(value, value_name, config_filepath):
    """
    Verify that a variable set in a config file is a positive non-null integer.
    Stop script and send error if not. 
    
    :param value: variable
    :param value_name: variable name
    :param config_filepath: filepath of config file where the variable is set
    """
    if not isinstance(value, int) or value <= 0:
        sys.exit(f"In {os.path.basename(config_filepath)}: {value_name}={value}. {value_name} must be an integer larger than 0.")

def check_pos_float_not_null(value, value_name, config_filepath):
    """
    Verify that a variable set in a config file is a positive non-null float.
    Stop script and send error if not. 
    
    :param value: variable
    :param value_name: variable name
    :param config_filepath: filepath of config file where the variable is set
    """
    if not isinstance(value, float) or value <= 0:
        sys.exit(f"In {os.path.basename(config_filepath)}: {value_name}={value}. {value_name} must be an integer larger than 0.")

def check_boolean(value, value_name, config_filepath):
    """
    Verify that a variable is boolean, either 'True' or 'False'
    """
    if not isinstance(value, bool):
        sys.exit(f"In {os.path.basename(config_filepath)}: {value_name}={value}. {value_name} must be set to either 'True' or 'False'")

def check_string(value, value_name, config_filepath):
    """
    Verify that a variable is a string
    """
    if not isinstance(value, str):
        sys.exit(f"In {os.path.basename(config_filepath)}: {value_name}={value}. {value_name} must be a valid string")

def check_list_strings(yaml_list, list_name, config_filepath):
    if not isinstance(yaml_list, list):
        raise ValueError(f"In {os.path.basename(config_filepath)}: {list_name} must be a non-empty list.")
    if len(yaml_list) == 0:
        raise ValueError(f"In {os.path.basename(config_filepath)}: {list_name} at least 1 marker must be choosen for clustering.")
    if not all(isinstance(x, str) and x.strip() for x in yaml_list):
        raise ValueError(f"In {os.path.basename(config_filepath)}: {list_name} must contain marker names in ''.")







##########################
#### CLUSTERING UTILS ####
##########################

import pandas as pd
import numpy as np
import anndata as ad
from sklearn.decomposition import PCA, KernelPCA, FastICA
from sklearn.manifold import Isomap


def downsample_pd_dataframe(dataframe, sample_size):
    """
    Samples n rows from each file in the dataframe
    Returns downsampled dataframe
    
    :param dataframe: pandas dataframe with a 'filename' column
    :param sample_size: number of rows to sample per filename
    """
    
    print(f"Sampling {sample_size} cells from each file")

    small_df = (
        dataframe
        .groupby("filename", group_keys=False)
        .apply(lambda x: x.sample(n=min(sample_size, len(x)), random_state=17))
        .reset_index(drop=True)
    )

    unique_filenames = small_df["filename"].unique()
    print(f"Number of files sampled: {len(unique_filenames)}")
    print(f"Total number of cells in dataset after sampling: {small_df.shape[0]}")
    print("\n")

    return small_df

def selective_downsample_pd_dataframe(dataframe, keyword, sample_size):
    """
    Samples n rows from each file in the dataframe, but only if the
    filename contains an exact keyword
    Returns downsampled dataframe
    
    :param dataframe: pandas dataframe with a 'filename' column
    :param keyword: string that marks from which file to sample
    :param sample_size: number of rows to sample per filename
    """
    
    print(f"Sampling at most {sample_size} rows from files containing '{keyword}'")

    filtered_df = dataframe[dataframe["filename"].str.contains(keyword, regex=False)]
    
    if filtered_df.empty:
        sys.exit(f"No files matched the keyword '{keyword}'")

    small_df = (
        filtered_df
        .groupby("filename", group_keys=False)
        .apply(lambda x: x.sample(n=min(sample_size, len(x)), random_state=17))
        .reset_index(drop=True)
    )

    unique_filenames = small_df["filename"].unique()
    print(f"Number of files sampled: {len(unique_filenames)}")
    print(f"Total number of cells in dataset after sampling: {small_df.shape[0]}")
    print("\n")

    return small_df


def fit_transform_dr(method_choice, input_df, clustering_markers=None,markers_to_reduce=None, n_components=None, n_neighbors=None):
    """
    Fit a dimensionality reduction method of choice on dataframe and transform the data
    Method optionsa are: directly use the markers values, PCA, KernelPCA, Isomap, or FastICA
    
    :param method_choice: string from dim_reduction.yaml config file indicating chosen method
    :param input_df: pandas dataframe that must be transformed
    :rest: parameters for each method from dim_reduction.yaml

    Returns transformed data as pandas dataframe, names of columns, and the fitted reducer
    """
    method = method_choice.lower()

    if method == "direct_markers":
        print(f"Dimensionality reduction via {method_choice} with {clustering_markers}", flush=True)
        small_reduced_df = input_df[clustering_markers].copy()
        reduced_columns = clustering_markers
        reducer = None # placeholder (because function returns 3 things)
        return small_reduced_df, reduced_columns, reducer

    elif method == "pca":
        print(f"Dimensionality reduction via {method_choice} on {markers_to_reduce}. n_components={n_components}", flush=True)
        reducer = PCA(n_components=n_components)
        reduced_columns = [f"PC{x}" for x in range(1, n_components + 1)]

    elif method == "kernelpca":
        print(f"Dimensionality reduction via {method_choice} on {markers_to_reduce}. n_components={n_components} and kernel=rbf", flush=True)
        reducer = KernelPCA(n_components=n_components, kernel='rbf', n_jobs=-1)
        reduced_columns = [f"kPCA{x}" for x in range(1, n_components + 1)]

    elif method == "isomap":
        print(f"Dimensionality reduction via {method_choice} on {markers_to_reduce}. n_components={n_components} and n_neighbors={n_neighbors}", flush=True)
        reducer = Isomap(n_components=n_components, n_neighbors=n_neighbors, n_jobs=-1)
        reduced_columns = [f"IM{x}" for x in range(1, n_components + 1)]

    elif method == "fastica":
        print(f"Dimensionality reduction via {method_choice} on {markers_to_reduce}. n_components={n_components}.", flush=True)
        reducer = FastICA(n_components=n_components, random_state=17)
        reduced_columns = [f"FI{x}" for x in range(1, n_components + 1)]
    
    else:
        sys.exit(f"Unrecognized dimensionality reduction method: {method_choice}")
    
    index = input_df.index
    input_dr = input_df[markers_to_reduce].to_numpy(dtype="float32")
    reduced = reducer.fit_transform(input_dr)
    small_reduced_df = pd.DataFrame(reduced, index=index, columns=reduced_columns)

    if method == "pca":
        print(f"Percentage of variance explained by each of the PCA components: {reducer.explained_variance_ratio_}", flush=True)
    
    return small_reduced_df, reduced_columns, reducer


def transform_dr(method_choice, input_df, reducer, chunksize, clustering_markers=None, markers_to_reduce=None, n_components=None, n_neighbors=None):
    """
    Transform a pandas dataframe with a fitted reducted (from fit_transform_dr)
    
    :param method_choice: string from dim_reduction.yaml config file indicating chosen method
    :param input_df: pandas dataframe that must be transformed
    :reducer: fitted reducer produced by fit_transform_dr
    :chunksize: int that indicates how many rows to transform in one go
    """
    method = method_choice.lower()

    if method == "direct_markers":
        print(f"Dimensionality reduction via {method_choice} with {clustering_markers}", flush=True)
        small_reduced_df = input_df[clustering_markers].copy()
        reduced_columns = clustering_markers
        return small_reduced_df, reduced_columns

    if method == "pca":
        print(f"Applying {method_choice} model", flush=True)
        reduced_columns = [f"PC{x}" for x in range(1, n_components + 1)]

    elif method == "kernelpca":
        print(f"Applying {method_choice} model", flush=True)
        reduced_columns = [f"kPCA{x}" for x in range(1, n_components + 1)]

    elif method == "isomap":
        print(f"Applying {method_choice} model", flush=True)
        reduced_columns = [f"IM{x}" for x in range(1, n_components + 1)]

    elif method == "fastica":
        print(f"Applying {method_choice} model", flush=True)
        reduced_columns = [f"FI{x}" for x in range(1, n_components + 1)]

    else:
        sys.exit(f"Unrecognized dimensionality reduction method: {method_choice}")


    input_df_dr = input_df[markers_to_reduce]

    n_rows = input_df_dr.shape[0]
    total_chunks = (n_rows + chunksize - 1) // chunksize
    reduced_chunks = []
    
    for i, start in enumerate(range(0, n_rows, chunksize), start=1):
        end = min(start + chunksize, n_rows)
        chunk_start = datetime.now()
        print(chunk_start.strftime("%Y-%m-%d %H:%M:%S"), f"- Currently processing chunk # {i}/{total_chunks} (rows {start}:{end})", flush=True)
        end = start + chunksize

        chunk = input_df_dr.iloc[start:end]
        index = chunk.index
        chunk_np = chunk.to_numpy(dtype="float32")

        reduced_chunk = reducer.transform(chunk_np)

        chunk_df = pd.DataFrame(
            reduced_chunk,
            index=index,
            columns=reduced_columns
        )

        reduced_chunks.append(chunk_df)

    small_reduced_df = pd.concat(reduced_chunks)

    if method == "pca":
        print(f"Percentage of variance explained by each PCA component: {reducer.explained_variance_ratio_}")

    return small_reduced_df, reduced_columns



def pd_df_to_anndata(X, df, var_names, sample_key=None):
    columns = ["cell_id", "filename"]

    if sample_key is not None:
        columns.append(sample_key)

    adata = ad.AnnData(X=X,
                    obs=df[columns].copy(),
                    var=pd.DataFrame(index=var_names))
    return adata



def attach_direct_cluster_labels(downsampled_df_with_clusters, complete_df):
    """
    Attach cluster labels to complete dataframe
    Add extra column "label_propagation" to distinguish between cells that received
    their cluster label directly from the algorithm (==0) or propagated via NeareastNeighbor (==1)
    
    :param downsampled_df_with_clusters: downsampled_dataframe that was used for direct clustering, with labels attached
    :param complete_df: complete dataframe with full dataset
    """

    # Add the clustering results to the full dataframe
    samples_with_clusters = complete_df.merge(
        downsampled_df_with_clusters[["cell_id", "cluster"]],
        on="cell_id",
        how="left"
    )

    # Initialize label_propagation: assume propagated (1) by default
    samples_with_clusters["label_propagation"] = 1
    samples_with_clusters["label_propagation"] = (samples_with_clusters["label_propagation"].astype("uint8"))

    # Change label_propagation for cells that received a cluster directly from the algorithm (value 0)
    samples_with_clusters.loc[
        samples_with_clusters["cluster"].notna(), 
        "label_propagation"
    ] = 0

    return samples_with_clusters


def get_nb_unlabeled_cells(full_df_with_clusters):
    """
    Returns number of cells without a cluster label in a dataframe
    Prints the number
    
    :param full_df_with_clusters: Description
    """
    nb_unlabeled = full_df_with_clusters['cluster'].isna().sum()
    print(f"There are {nb_unlabeled} cells without a cluster assignment.")
    return nb_unlabeled

def majority_vote_1d(arr):
    """
    Majority vote to assign a cluster label to an unlabeled cell
    The unlabeled cell will get the most popular cluster label of its neighbors
    
    :param arr: array of the neighbors' cluster labels
    """
    values, counts = np.unique(arr, return_counts=True)
    return values[np.argmax(counts)]

def propagate_labels_by_chunk(full_df_with_clusters, chunksize, nearest_neigh_model, cluster_labels, method_choice, reducer=None, clustering_markers=None, markers_to_reduce=None):
    """
    Transform values with the reducer produced by fit_transform_dr
    Propagate cluster labels to cells without a cluster assignment via NearestNeighbors
    Proceed chunk by chunk

    :param full_df_with_clusters: complete dataset, with "cluster" column
    :param chunksize: number of cells to process in a go
    :param nearest_neigh_model: fitted NearestNeighbors model
    :param cluster_labels: cluster labels of training cells
    :param method_choice: dimensionality reduction method
    :param reducer: fitted DR model (PCA, KernelPCA, etc.)
    :param clustering_markers: markers used when method_choice == "direct_markers"
    :param markers_to_reduce: markers used as input for DR
    """

    method = method_choice.lower()

    print(f"Dimensionality reduction with {method} and label propagation by chunk")

    n_rows = full_df_with_clusters.shape[0]
    total_chunks = (n_rows + chunksize - 1) // chunksize

    for i, start in enumerate(range(0, n_rows, chunksize), start=1):
        end = min(start + chunksize, n_rows)
        chunk_start = datetime.now()
        print(chunk_start.strftime("%Y-%m-%d %H:%M:%S"),f"- Currently processing chunk # {i}/{total_chunks} (rows {start}:{end})", flush=True)

        chunk = full_df_with_clusters.iloc[start:end]
        unlabeled_mask = chunk["cluster"].isna()

        if not unlabeled_mask.any():
            continue

        unlabeled_chunk = chunk.loc[unlabeled_mask]

        # DR for the chunk
        if method == "direct_markers":
            X_chunk = unlabeled_chunk[clustering_markers].to_numpy(dtype="float32")
        else:
            X_input = unlabeled_chunk[markers_to_reduce].to_numpy(dtype="float32")
            X_chunk = reducer.transform(X_input)

        # Propagate cluster label by nearest neighbor
        idx = nearest_neigh_model.kneighbors(X_chunk, return_distance=False)
        neighbor_labels = cluster_labels[idx]

        voted_labels = np.apply_along_axis(
            majority_vote_1d,
            axis=1,
            arr=neighbor_labels
        ).ravel()

        full_df_with_clusters.loc[unlabeled_chunk.index, "cluster"] = voted_labels

    return full_df_with_clusters


def propagate_labels_by_chunk_cytovi(full_df_with_clusters, adata_full, chunksize, nearest_neigh_model, cluster_labels):
    """
    Propagated cluster labels to cells without a cluster assignment via NearestNeighbor
    using CytoVI latent sapce
    Proceed chunk by chunk

    Assumes
    - Row order of full_df_with_clusters matches adata_full.obs
    - cluster_labels correspond to training latent space
    
    :param full_df_with_clusters: complete dataset, with "cluster" column
    :param chunksize: number of cells to process in a go
    :param nearest_neigh_model: NearestNeighbor model that has been fit
    :param adata_full: anndata containing latent space representation from CytoVI model
    :param cluster_labels: direct cluster labels that will be used for majority vote
    """
    # checks
    assert full_df_with_clusters.shape[0] == adata_full.n_obs
    assert "X_cytovi" in adata_full.obsm
    assert np.array_equal(
        full_df_with_clusters["cell_id"].values,
        adata_full.obs["cell_id"].values)

    n_rows = full_df_with_clusters.shape[0]
    total_chunks = (n_rows + chunksize - 1) // chunksize

    latent_full = adata_full.obsm["X_cytovi"]

    for i, start in enumerate(range(0, n_rows, chunksize), start=1):
        end = min(start + chunksize, n_rows)
        chunk_start = datetime.now()
        print(chunk_start.strftime("%Y-%m-%d %H:%M:%S"),f"- Currently processing chunk # {i}/{total_chunks} (rows {start}:{end})", flush=True)

        chunk = full_df_with_clusters.iloc[start:end]
        unlabeled_mask = chunk["cluster"].isna()

        if not unlabeled_mask.any():
            continue

        # Direct row-aligned indexing into latent space
        global_indices = np.arange(start, end)[unlabeled_mask.values]

        X_chunk = latent_full[global_indices].astype("float32")

        idx = nearest_neigh_model.kneighbors(X_chunk,return_distance=False)

        neighbor_labels = cluster_labels[idx]
        voted_labels = np.apply_along_axis(
            majority_vote_1d,
            axis=1,
            arr=neighbor_labels
        ).ravel()

        full_df_with_clusters.loc[
            chunk.index[unlabeled_mask],
            "cluster"
        ] = voted_labels

    return full_df_with_clusters


def print_cluster_size_bef_aft_propagation(df_before, df_after, groupby):
    """
    Print cluster sizes before and after propagagtion for verification
    
    :param df_before: pandas dataframe, before propagation
    :param df_after: pandas dataframe, after cluster label propagation
    :str groupby: group the cells by cluster typically (could be metacluster)
    """
    print(f"Number of cells in each {groupby} from the sampled dataset:")
    print(df_before.groupby(groupby, observed=True).size())

    print("Number of cells in each cluster from the complete dataset:")
    print(df_after.groupby(groupby, observed=True).size())

def save_standard_clustering_output(full_df_with_clusters, output_path, other_columns=None):
    """
    Save a standard csv for clustering results.
    Always includes:
      - filename
      - row_id
      - cluster
      - label_propagation

    Optionally includes additional columns.
    """
    base_columns = ["filename", "row_id", "cluster", "label_propagation"]

    if other_columns is None:
        columns = base_columns
    elif isinstance(other_columns, str):
        columns = base_columns + [other_columns]
    else:
        columns = base_columns + list(other_columns)

    clusters_df = full_df_with_clusters[columns]
    clusters_df.to_csv(output_path, index=False)




#####################################
#### CLUSTERING EVALUATION UTILS ####
#####################################

import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import pickle
import umap
from sklearn.metrics import silhouette_score, calinski_harabasz_score, davies_bouldin_score


def count_rows_per_condition(dataframe, groupby, extra_groupby=None):
    """
    Count the number of cells per condition (ex in each batch)
    Returns pandas dataframe with the counts per condition
    
    :param dataframe : pandas dataframe
    :param groupby : str correspond to a column in the dataframe
    :param extra_groupy : any other column
    """
    if extra_groupby is None:
        columns = [groupby]
    elif isinstance(extra_groupby, str):
        columns = [groupby, extra_groupby]
    else:
        raise TypeError("extra_groupby must be None or str")

    counts = (
        dataframe
        .groupby(columns)
        .size()
        .reset_index(name="n_cells")
    )
    return counts

def cluster_cell_counts(dataframe, clustering_key):
    """
    Count the total number of cells in a cluster
    Returns pandas dataframe with counts per cluster
    
    :param dataframe : pandas dataframe
    :param clustering_key : granularity of clustering (cluster or metacluster)
    """
    totals = (
        dataframe
        .groupby(clustering_key)
        .size()
        .reset_index(name="total_cells")
    )
    return totals

def calculate_pct(counts, totals, clustering_key):
    """
    Calculates percentage of counts on totals
    
    :param counts : dataframe produced by count_rows_per_condition
    :param totals : dataframe produced by cluster_cell_counts
    """
    pct = counts.merge(totals, on=clustering_key)
    pct["pct_cells"] = pct["n_cells"] / pct["total_cells"] * 100
    return pct


def batch_contribution_qc(qc_df, clustering_key, out_png, out_csv):
    """
    Calculates the number of cells in each batch contributing to the cluster formation (as percentage)
    Plots a heatmap and saves a csv of counts
    
    :param qc_df : dataframe with all cells and their corresponding "batch" and "cluster" columns
    :param clustering_key : granularity of clustering (cluster or metacluster)
    :param out_png : filepath to save heatmap
    :param out_csv : filepath to save dataframe with results
    """
    # Count cells per batch per cluster
    counts = count_rows_per_condition(qc_df, clustering_key, "batch")
    
    # Count the total number of cells in a clustering level
    totals = cluster_cell_counts(qc_df, clustering_key)

    # Calculate the percentages
    pct = calculate_pct(counts, totals, clustering_key)
    print_pd_dataframe(pct, f"Number of total cells in a {clustering_key} and number of cells per batch contributating to that {clustering_key}")

    # Dataframe for heatmap creation
    heatmap_df = ( # pivot axis to have clusters as rows and batch number as columns
        pct
        .pivot(index=clustering_key, columns="batch", values="pct_cells")
        .fillna(0)
    )

    # make sure that sum of percentages is about 100% for each cluster
    assert np.allclose(heatmap_df.sum(axis=1), 100, atol=1e-6)

    # save heatmap
    sns.heatmap(heatmap_df, cmap="rocket_r")
    plt.title(f"Batch contribution to {clustering_key} formation")
    plt.savefig(out_png)
    plt.close()

    # save dataframe used for plotting as csv
    heatmap_df.to_csv(out_csv)
    print(f">> Saved batch contribution plot to {out_png} and results to {out_csv} \n")



def file_contribution_qc(qc_df, clustering_key, out_png, out_csv):
    """
    Calculates the number of cells in each file contributing to the cluster formation (as percentage)
    Plots a heatmap and saves a csv of counts
    
    :param qc_df : dataframe with all cells and their corresponding "file" and "cluster" columns
    :param clustering_key : granularity of clustering (cluster or metacluster)
    :param out_png : filepath to save heatmap
    :param out_csv : filepath to save dataframe with results
    """
    # Count cells per filename per cluster
    counts = count_rows_per_condition(qc_df, clustering_key, "filename")
    
    # Count the total number of cells in a clustering level
    totals = cluster_cell_counts(qc_df, clustering_key)

    # Calculate the percentages
    pct = calculate_pct(counts, totals, clustering_key)
    print_pd_dataframe(pct, f"Number of total cells in a {clustering_key} and number of cells per file contributating to that {clustering_key}")

    # changing orientation for seaborn heatmap: rows are clusters and columns are filenames
    heatmap_df = (
        pct
        .pivot(
            index=clustering_key,
            columns="filename",
            values="pct_cells"
        )
        .fillna(0)
    )

    # make sure each cluster should sum to about 100%
    assert np.allclose(heatmap_df.sum(axis=1), 100, atol=1e-6)
    
    # save heatmap
    sns.heatmap(heatmap_df, cmap="rocket_r")
    plt.xticks([]) # don't show column names (filenames)
    plt.title(f"File contribution to {clustering_key} formation")
    plt.savefig(out_png)
    plt.close()

    # save dataframe used for plotting as csv
    pct.to_csv(out_csv)

    print(f">> Saved file contribution plot to {out_png} and results to {out_csv} \n")

def sort_and_reset(expr_df, clustering_key):
    """
    Sort dataframe by clustering key and reset index.
    """
    return expr_df.sort_values(clustering_key).reset_index(drop=True)

def get_marker_columns(expr_df, clustering_key, extra_cols=None):
    """
    Return columns representing markers, excluding specified extra columns.
    """
    if extra_cols is None:
        extra_cols = ["filename", "row_id", "label_propagation", "cluster", "metacluster", clustering_key]
    return [col for col in expr_df.columns if col not in extra_cols]

def downsample_dataframe(data, sample_size, clustering_key, random_state=17):
    downsampled_df = (
        data
        .groupby(clustering_key, group_keys=False)
        .apply(lambda x: x.sample(n=min(sample_size, len(x)), random_state=random_state))
        .reset_index(drop=True)
    )
    return downsampled_df

def aggregate_median(expr_df, clustering_key, marker_cols):
    """
    Compute median per marker per cluster.
    """
    return expr_df.groupby(clustering_key)[marker_cols].median()

def cluster_size_qc(qc_df, clusters_df, clustering_key, out_png, out_csv):
    # count total sizes of clustering levels and sort (for neater visualization)
    totals = cluster_cell_counts(qc_df, clustering_key)
    totals_sorted = (totals.sort_values("total_cells", ascending=False))

    # Set cluster number as a string (not a float)
    totals_sorted[clustering_key] = totals_sorted[clustering_key].astype(str)

    # Calculate cluster size as percentage of total cells
    nb_total_cells = clusters_df.shape[0] # number of rows in clusters table 
    totals_sorted.insert(2, "pct_total", (totals_sorted["total_cells"]/nb_total_cells)*100)
    print_pd_dataframe(totals_sorted, f"Number of cells in each {clustering_key} and percentages of all cells")

    # Draw barplot
    nb_clusters = totals.shape[0]
    fig_width = max(10, nb_clusters*0.15)
    fig_height = 10

    plt.figure(figsize=(fig_width, fig_height))
    sns.barplot(totals_sorted, x=clustering_key, y="total_cells", color="firebrick")
    plt.xticks(rotation=90, fontsize=7)
    plt.yticks(fontsize=8)
    plt.ylabel("Number of cells")
    ax2 = plt.twinx() # to add a second y axis range corresponding to number of cells
    sns.barplot(totals_sorted, x=clustering_key, y="pct_total", color="firebrick")
    plt.xticks(rotation=90, fontsize=7)
    plt.yticks(fontsize=8)
    plt.ylabel("Percentage of total number of cells")
    plt.title(f"{clustering_key} sizes")
    plt.savefig(out_png)
    plt.close()

    # save dataframe used for plotting as csv
    totals_sorted.to_csv(out_csv)

    print(f">> Saved {clustering_key} size barplot to {out_png} and results to {out_csv} \n")



def median_expr_qc(expr_df, clustering_key, out_png, out_csv):
    # sort the dataframe per clustering level number
    expr_df_sorted = sort_and_reset(expr_df, clustering_key)

    # get the median per marker for all cells in a level
    marker_cols = get_marker_columns(expr_df_sorted, clustering_key)
    expr_df_sorted_median = aggregate_median(expr_df_sorted, clustering_key, marker_cols)
    print_pd_dataframe(expr_df_sorted_median, f"Median marker expression values per {clustering_key}")

    # draw heatmap
    nb_clusters = expr_df_sorted_median.shape[0]
    nb_markers = expr_df_sorted_median.shape[1]
    fig_width = max(10, nb_markers*0.2)
    fig_height = max(10, nb_clusters*0.15)

    plt.figure(figsize=(fig_width, fig_height))
    sns.heatmap(expr_df_sorted_median,cmap="rocket_r", cbar_kws={"shrink": 0.4})
    plt.yticks(fontsize=8)
    plt.title("Median marker expression in clusters")
    plt.tight_layout()
    plt.savefig(out_png)
    plt.close()

    # save dataframe used for plotting as csv
    expr_df_sorted_median.to_csv(out_csv)

    print(f">> Saved median marker exprssion plot to {out_png} and results to {out_csv} \n")



def sample_expr_qc(expr_df, clustering_key, sample_size, out_png, out_csv):
    # sort dataframe by cluster
    expr_df_sorted = sort_and_reset(expr_df, clustering_key)

    # sample n cells per cluster
    expr_df_sampled = downsample_dataframe(expr_df_sorted, sample_size, clustering_key)

    # sort again to ensure clusters are consecutive after sampling
    expr_df_sampled = expr_df_sampled.sort_values(clustering_key).reset_index(drop=True)

    # save sampled data
    expr_df_sampled.to_csv(out_csv, index=False)

    # extract cluster info
    cluster_series = expr_df_sampled[clustering_key].astype("category")
    cluster_sizes = cluster_series.value_counts(sort=False)
    cluster_order = cluster_sizes.index.tolist()

    # keep only marker columns for heatmap
    marker_cols = get_marker_columns(expr_df_sampled, clustering_key)
    heatmap_df = expr_df_sampled[marker_cols]

    # color palette per cluster
    palette = sns.color_palette("Paired", n_colors=len(cluster_order))
    cluster_colors = dict(zip(cluster_order, palette))

    # figure sizing
    nb_markers = heatmap_df.shape[1]
    nb_cells = heatmap_df.shape[0]
    fig_width = max(10, nb_markers * 0.2)
    fig_height = max(10, nb_cells * 0.02)

    fig, (ax_colors, ax_heatmap) = plt.subplots(
        ncols=2,
        figsize=(fig_width, fig_height),
        gridspec_kw={"width_ratios": [0.4, 9]}
    )

    # annotation axis setup
    ax_colors.set_xlim(0, 1)
    ax_colors.set_ylim(0, nb_cells)
    ax_colors.invert_yaxis()
    ax_colors.axis("off")

    # draw cluster blocks
    y_start = 0
    for cluster in cluster_order:
        height = cluster_sizes[cluster]
        color = cluster_colors[cluster]

        rect = patches.Rectangle(
            (0, y_start),
            width=1,
            height=height,
            facecolor=color,
            edgecolor="black",
            linewidth=0.3
        )
        ax_colors.add_patch(rect)

        ax_colors.text(
            0.5,
            y_start + height / 2,
            str(cluster),
            ha="center",
            va="center",
            fontsize=6
        )

        y_start += height

    # plot heatmap (cell-level)
    sns.heatmap(
        heatmap_df,
        ax=ax_heatmap,
        cmap="rocket_r",
        cbar=True,
        yticklabels=False,
        cbar_kws={"shrink": 0.25}
    )

    plt.title(f"Sampled cells' marker expression per {clustering_key}")
    plt.tight_layout()
    plt.savefig(out_png, dpi=300)
    plt.close()

    print(f">> Saved sample's marker expression plot to {out_png} and results to {out_csv}\n")


def load_model(model_pkl):
    with open(model_pkl, "rb") as file:
        model = pickle.load(file)
    return model


def run_umap_markers(data_df, clustering_markers, n_neighbors=15, min_dist=0.3, metric="euclidean"):
    # Initialize UMAP object
    reducer = umap.UMAP(
        n_neighbors=n_neighbors,
        min_dist=min_dist,
        metric=metric,
        random_state=17)
    
    # Get embedding
    embedding = reducer.fit_transform(data_df[clustering_markers])
    
    # Add embedding to dataframe
    data_df["UMAP1"] = embedding[:,0]
    data_df["UMAP2"] = embedding[:,1]

    return data_df, reducer


def define_centroids(data_df, clustering_markers, clustering_key):
    # Centroids are defined by the cells directly clustered (label_propagation==0)
    # The centroids are defined by the median marker expression
    cluster_centroids = (
        data_df.loc[data_df["label_propagation"]==0]
        .groupby(clustering_key, observed=True)[clustering_markers]
        .median()
    )
    centroid_values = cluster_centroids.values
    return centroid_values, cluster_centroids


def build_legend(data_df, clustering_key):
    # Build colormap (1 color for each cluster)
    nb_clusters = data_df[clustering_key].nunique()
    cmap = plt.get_cmap("tab20", nb_clusters)
    
    # Build legend elements
    handles = []
    for cluster, color in zip(data_df[clustering_key].cat.categories, cmap.colors):
        handles.append(plt.Line2D([], [], color=color, marker='o', linestyle='', markersize=5, label=str(cluster)))
    handles.append(plt.Line2D([], [], color='black', marker='X', linestyle='', markersize=5, label='Centroid'))
    handles.append(plt.Line2D([], [],
            marker='o',
            linestyle='',
            markersize=5,
            markerfacecolor='white',
            markeredgecolor='black',
            label='Algorithm assigned'))
    handles.append(
        plt.Line2D(
            [], [],
            marker='o',
            linestyle='',
            markersize=5,
            markerfacecolor='lightgray',
            markeredgecolor='none',
            label='Inferred'
        )
    )
    
    # Dot colors
    colors = cmap(data_df[clustering_key].cat.codes)

    return handles, colors, cmap

def plot_umap(data_df, clustering_key, colors, cmap, handles, umap_png, centroids_umap, cluster_centroids, plot_centroids=False):
    # Initialize figure
    fig, ax = plt.subplots()

    # Get mask for directly labeled cells and propagated label cells
    mask_direct = data_df["label_propagation"] == 0
    mask_prop = ~mask_direct

    # Plot propagated cells (no outline)
    ax.scatter(
        data_df.loc[mask_prop, "UMAP1"],
        data_df.loc[mask_prop, "UMAP2"],
        c=colors[mask_prop],
        alpha=0.8,
        s=2,
        edgecolor="none"
    )

    # Plot directly clustered cells with black outline
    ax.scatter(
        data_df.loc[mask_direct, "UMAP1"],
        data_df.loc[mask_direct, "UMAP2"],
        c=colors[mask_direct],
        alpha=0.8,
        s=6,
        edgecolor="black",
        linewidth=0.5
    )

    # Optional: overlay centroids
    if plot_centroids != False:
        centroid_colors = cmap(cluster_centroids.index.codes)
        ax.scatter(
            centroids_umap[:,0],
            centroids_umap[:,1],
            c=centroid_colors,
            s=50,
            edgecolor="black",
            marker="X",
            label="Centroid"
        )
    
    # Add legend
    ax.legend(handles=handles, title=clustering_key, bbox_to_anchor=(1, 0.5), loc="center left", fontsize=6, title_fontsize=8)

    ax.set_title(f"{clustering_key} with label propagation")
    ax.set_xlabel("UMAP1")
    ax.set_ylabel("UMAP2")
    plt.tight_layout()
    plt.savefig(umap_png, dpi=300)
    plt.close()

    print(f">> Saved UMAP visualisation to {umap_png} \n")

def get_marker_columns_for_coloring(df, exclude_cols):
    """
    Return all marker columns that can be used to color UMAPs.
    """
    return [c for c in df.columns if c not in exclude_cols]

def plot_umap_marker_grid(data_df,marker_cols,out_png,ncols=5,point_size=2,cmap="mako"):
    """
    Plot a grid of UMAPs colored by marker expression.

    Parameters
    ----------
    data_df : pd.DataFrame that must contain UMAP1 and UMAP2.
    marker_cols : list Marker columns to color by.
    out_png : str Output PNG path.
    """
    assert "UMAP1" in data_df.columns
    assert "UMAP2" in data_df.columns

    n_markers = len(marker_cols)
    nrows = int(np.ceil(n_markers / ncols))

    fig, axes = plt.subplots(
        nrows=nrows,
        ncols=ncols,
        figsize=(ncols * 3, nrows * 3),
        squeeze=False
    )

    for ax, marker in zip(axes.flat, marker_cols):
        sc = ax.scatter(
            data_df["UMAP1"],
            data_df["UMAP2"],
            c=data_df[marker],
            s=point_size,
            cmap=cmap
        )
        ax.set_title(marker, fontsize=8)
        ax.set_xticks([])
        ax.set_yticks([])

        plt.colorbar(sc, ax=ax, fraction=0.046, pad=0.04)

    # Turn off unused axes
    for ax in axes.flat[n_markers:]:
        ax.axis("off")

    plt.tight_layout()
    plt.savefig(out_png, dpi=300)
    plt.close()

def compute_clustering_metrics(X, labels, sample_size):
    sil = silhouette_score(X, labels=labels, metric="euclidean", sample_size=sample_size)
    ch = calinski_harabasz_score(X, labels)
    db = davies_bouldin_score(X, labels=labels)
    return sil, ch, db