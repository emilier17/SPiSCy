#==================================================================================================================#
# SPiSCy (Snakemake PIpeline for Spectral CYtometry)
#
# Phase: clustering
# Step: evaluate the quality of clustering results
# 
# Author: Émilie Roy
# Date: Dec 2025
# Version: v.1.0
#
# Input : standardized clustering results file (clusters.csv)
# Outputs : model quality control plots and graphs (.png)
#==================================================================================================================#

from clustering_utils import (
    snakemake_logs, start_time,  end_time, elapsed_time, print_pd_dataframe,  batch_contribution_qc,
    file_contribution_qc, cluster_size_qc, sample_expr_qc, median_expr_qc, downsample_dataframe,
    load_model, define_centroids, run_umap_markers, build_legend, plot_umap, plot_umap_marker_grid,
    get_marker_columns_for_coloring, check_pos_int_not_null, check_boolean
    )

snakemake_logs(snakemake)
start = start_time()




#### IMPORTS AND SNAKEMAKE VARIABLES #####

# Library
import pandas as pd
import yaml


# Snakemake variables
# Inputs
clusters_csv = snakemake.input[0]
metadata_csv = snakemake.input[1]
marker_csv = snakemake.input[2]
spe_config = snakemake.input[3]
model_pkl = snakemake.input[4]
# Outputs
batch_contribution_png = snakemake.output[0]
batch_contribution_csv = snakemake.output[1]
file_contribution_png = snakemake.output[2]
file_contribution_csv = snakemake.output[3]
cluster_size_png = snakemake.output[4]
cluster_size_csv = snakemake.output[5]
marker_exp_median_png = snakemake.output[6]
marker_exp_median_csv = snakemake.output[7]
marker_exp_sample_png = snakemake.output[8]
marker_exp_sample_csv = snakemake.output[9]
umap_png = snakemake.output[10]
umap_markers_png = snakemake.output[11]
# Params
clustering_level = snakemake.params.clustering_key

print(f"Clustering_level: {clustering_level}")

# Config files
with open(spe_config) as f:
    config_file = yaml.safe_load(f)

sample_size = config_file["sample_size"]
sample_size_umap = config_file["sample_size_umap"]
umap_plot_centroids = config_file["umap_plot_centroids"]

check_pos_int_not_null(sample_size, "sample_size", spe_config)
check_pos_int_not_null(sample_size_umap, "sample_size_umap", spe_config)
check_boolean(umap_plot_centroids, "umap_plot_centroids", spe_config)




#### IMPORT STANDARDIZED CLUSTER AND MARKER FILES, AND METADATA ####

clusters = pd.read_csv(clusters_csv) # standard output file of clustering algorithm: cluster number that each cell belongs to
metadata = pd.read_csv(metadata_csv, usecols=["filename", "batch"])
markers = pd.read_csv(marker_csv) # marker expression values for each cell after all preprocessing steps

clusters[clustering_level] = clusters[clustering_level].astype("category") # set cluster label as category

print_pd_dataframe(clusters, "Clustering results")
print_pd_dataframe(markers, "Marker expression values")


# make sure that the clusters output contains all clustering_levels set in config
if clustering_level not in clusters.columns:
    raise ValueError(
        f"'{clustering_level}' requested in {spe_config}, but not found in clusters.csv"
        f"Available clustering levels: {list(clusters.columns)}"
    )




#### EXTRACT QC INFORMATION ####

## Merging metadata and clusters tables ##

# Add batch number to every cell
qc_df = clusters.merge(
    metadata[["filename", "batch"]],
    on="filename",
    how="left",
    validate="many_to_one"
)

# Make sure every filename is associated to a batch number
if qc_df["batch"].isna().any():
    missing = qc_df.loc[qc_df["batch"].isna(), "filename"].unique()
    raise ValueError(f"Missing batch for filenames: {missing}")

print_pd_dataframe(qc_df, "Metadata and clustering results merge")



## Merging clusters and markers tables ##

# Every cell has its cluster number (+metacluster number if available) and marker expression values
base_cols = ["filename", "row_id", "label_propagation"]
clustering_cols = [c for c in ["cluster", "metacluster"] if c in clusters.columns]

expr_df = markers.merge(
    clusters[base_cols + clustering_cols],
    on = ["filename", "row_id"],
    how = "inner",
    validate = "one_to_one"
).drop(columns=["filename"])

# make sure every cell has been merged
assert expr_df.shape[0] == markers.shape[0]

print_pd_dataframe(expr_df, "Marker and clustering results merge")




#### QC1: BATCH CONTRIBUTION TO CLUSTER FORMATION ####

#-- Heatmap that shows the batch composition of a cluster, in terms of % of cells in a cluster from each batch. 
#-- Allows to assess batch effect presence. 
#-- Each row is a cluster and each column is a batch number.
#-- If a majority of cells in a cluster come from a single batch, then a batch effect might be the cause.

batch_contribution_qc(qc_df,
                      clustering_key=clustering_level,
                      out_png=batch_contribution_png,
                      out_csv=batch_contribution_csv)




#### QC2: FILE CONTRIBUTION TO CLUSTER FORMATION ####

#-- Heatmap that shows the file composition of a cluster, in terms of % of cells in a cluster from each file. 
#-- Allows to assess if cells from different files are contributating equally to cluster formation. 
#-- Each row is a cluster and each column is a file.
#-- If a majority of cells in a cluster come from a single file, then investigating that file could be important.

file_contribution_qc(qc_df,
                     clustering_key=clustering_level,
                     out_png=file_contribution_png,
                     out_csv=file_contribution_csv)





#### QC3: NUMBER OF CELLS IN CLUSTERS ####

#-- Barplot that shows the number of cells in each cluster
#-- Double y axis for amount of cells as a raw number and percentage of total number of cells
#-- Allows to evaluate typical cluster sizes, if any very small clusters are present
#-- Useful if analyis is focused on identifying rare cell types

cluster_size_qc(qc_df=qc_df,
                clusters_df=clusters,
                clustering_key=clustering_level,
                out_png=cluster_size_png,
                out_csv=cluster_size_csv)




#### QC4: MARKER EXPRESSION IN CLUSTERS ####

#-- Heatmap showing marker expression in each cluster, similar to Seurat's DoHeatmap function
#-- Columns are markers and rows are grouped by cluster
#-- Rows: two granularity options. Either, median expression of the marker in the whole cluster, or see marker expression for n randomly sampled cells per cluster
#-- Allows to see if different clusters have similar expression patterns
#-- Allows to see which markers might be most useful for clustering (different expression levels among the cells) 

# Large Granularity: cluster's median expression

median_expr_qc(expr_df,
                clustering_key=clustering_level,
                out_png=marker_exp_median_png,
                out_csv=marker_exp_median_csv)


# Finer Granularity: sampled cells' expression

sample_expr_qc(expr_df,
               clustering_key=clustering_level,
               sample_size=sample_size,
               out_png=marker_exp_sample_png,
               out_csv=marker_exp_sample_csv)





#### QC5: UMAP VISUALIZATION OF CLUSTERS ####

#-- Run UMAP dimension reduction
#-- View 2D coordinates on plot

downsampled_df = downsample_dataframe(data=expr_df,
                                    sample_size=sample_size_umap,
                                    clustering_key=clustering_level)
model = load_model(model_pkl)
clustering_markers = model["clustering_markers"]
centroid_values, cluster_centroids = define_centroids(data_df=expr_df,
                                                    clustering_markers=clustering_markers,
                                                    clustering_key=clustering_level)
markers_umap_df, reducer = run_umap_markers(data_df=downsampled_df,
                                            clustering_markers=clustering_markers)
centroids_umap = reducer.transform(centroid_values)
handles, colors, cmap = build_legend(data_df=downsampled_df,
                                    clustering_key=clustering_level)
plot_umap(data_df=downsampled_df,
        clustering_key=clustering_level,
        colors=colors,
        cmap=cmap,
        handles=handles,
        umap_png=umap_png,
        centroids_umap=centroids_umap,
        cluster_centroids=cluster_centroids,
        plot_centroids=umap_plot_centroids)





#### QC6: MARKER EXPRESSION IN UMAP VISUALIZATION  ####

#-- View every marker's expression in the UMAP visualization of clusters
#-- Marker expression displayed as heatmap

non_marker_columns = ["filename", "row_id", "label_propagation", "cluster", "metacluster", "UMAP1", "UMAP2", clustering_level]
marker_columns = get_marker_columns_for_coloring(expr_df, exclude_cols = non_marker_columns)
plot_umap_marker_grid(markers_umap_df, marker_cols=marker_columns, out_png=umap_markers_png)




end = end_time()
elapsed_time(start, end)