#==================================================================================================================#
# SPiSCy (Snakemake PIpeline for Spectral CYtometry)
#
# Phase: clustering
# Step: perform clustering with HDBSCAN
# Original publicatipon : https://doi.org/10.1145/235968.233324
#
# Author: Émilie Roy
# Date: Jan 2026
# Version: v.1.0
#
# Input : final csv file (.csv)
# Outputs : standard output clusters.csv + model (.pkl)
#==================================================================================================================#

from clustering_utils import (
    snakemake_logs, start_time, check_pos_int_not_null, check_pos_float_not_null, end_time,
    elapsed_time, print_pd_dataframe, print_save_msg, downsample_pd_dataframe, 
    attach_direct_cluster_labels, get_nb_unlabeled_cells, print_cluster_size_bef_aft_propagation,
    save_standard_clustering_output, selective_downsample_pd_dataframe, check_boolean, 
    check_string, check_list_strings, fit_transform_dr, transform_dr 
    )

snakemake_logs(snakemake)
start = start_time()



#### IMPORTS AND SNAKEMAKE VARIABLES #####

# Library
import os
import sys
import pandas as pd
import numpy as np
import yaml
import pickle
from sklearn.cluster import Birch
from datetime import datetime



# Snakemake variables
# Inputs
samples_csv = snakemake.input[0]
clustering_config = snakemake.input[1]
sampling_config = snakemake.input[2]
label_prop_config = snakemake.input[3]
dr_config = snakemake.input[4]
# Outputs
model_pkl = snakemake.output[0]
clusters_csv = snakemake.output[1]


## clustering config file
with open(clustering_config) as f:
    config_file = yaml.safe_load(f)

sample_size = config_file["sample_size_per_file"]
radius_threshold = config_file["radius_threshold"]
branching_factor = config_file["branching_factor"]
nb_clusters = config_file["nb_clusters"]

check_pos_int_not_null(sample_size, "sample_size", clustering_config)
check_pos_float_not_null(radius_threshold, "radius_threshold", clustering_config)
check_pos_int_not_null(branching_factor, "branching_factor", clustering_config)
if isinstance(nb_clusters, str) and nb_clusters != "None":
    sys.exit(f"In {os.path.basename(clustering_config)}: nb_clusters can either be a positive integer or 'None'")
else:
    check_pos_int_not_null(nb_clusters, "nb_clusters", clustering_config)



## sampling config file
with open(sampling_config) as f:
    config_file = yaml.safe_load(f)

selective_sampling = config_file["selective_sampling"]
filename_keyword = config_file["filename_keyword"]

check_boolean(selective_sampling, "selective_sampling", sampling_config)
check_string(filename_keyword, "filename_keyword", sampling_config)


## label propagation config
with open(label_prop_config) as f:
    config_file = yaml.safe_load(f)

chunksize = config_file["chunksize"]

check_pos_int_not_null(chunksize, "chunksize", label_prop_config)



## dimensionality reduction config file
with open(dr_config) as f:
    config_file = yaml.safe_load(f)

# initialize dr variables
clustering_markers = None
markers_to_reduce = None
n_components = None
n_neighbors = None

method_choice = config_file["method_choice"]

if method_choice.lower() == "direct_markers":
    clustering_markers = config_file["direct_markers"]
    check_list_strings(clustering_markers, "direct_markers", dr_config)
    print(clustering_markers)

elif method_choice.lower() == "pca":
    n_components = config_file["PCA"]["n_components"]
    markers_to_reduce = config_file["PCA"]["markers_to_reduce"]
    check_pos_int_not_null(n_components, "n_components", dr_config)
    check_list_strings(markers_to_reduce, "markers_to_reduce", dr_config)

elif method_choice.lower() == "kernelpca":
    n_components = config_file["KernelPCA"]["n_components"]
    markers_to_reduce = config_file["KernelPCA"]["markers_to_reduce"]
    check_pos_int_not_null(n_components, "n_components", dr_config)
    check_list_strings(markers_to_reduce, "markers_to_reduce", dr_config)

elif method_choice.lower() == "isomap":
    n_components = config_file["Isomap"]["n_components"]
    n_neighbors = config_file["Isomap"]["n_neighbors"]
    markers_to_reduce = config_file["Isomap"]["markers_to_reduce"]
    check_pos_int_not_null(n_components, "n_components", dr_config)
    check_pos_int_not_null(n_neighbors, "n_neighbors", dr_config)
    check_list_strings(markers_to_reduce, "markers_to_reduce", dr_config)

elif method_choice.lower() == "fastica":
    n_components = config_file["FastICA"]["n_components"]
    markers_to_reduce = config_file["FastICA"]["markers_to_reduce"]
    check_pos_int_not_null(n_components, "n_components", dr_config)
    check_list_strings(markers_to_reduce, "markers_to_reduce", dr_config)

else:
    sys.exit(f"In {os.path.basename(dr_config)}, unrecognized method_choice='{method_choice}'. Options for dimensionality reduction techniques are 'direct_markers', 'PCA', 'KernelPCA', 'Isomap', or 'FastICA'.")




#### IMPORT PREPROCESSED DATA ####

print("Importing preprocessed data...", flush=True)

# Only import the columns we need
if method_choice.lower() == "direct_markers":
    samples = pd.read_csv(samples_csv, usecols=["filename", "row_id"] + clustering_markers)
else:
    samples = pd.read_csv(samples_csv)
print_pd_dataframe(samples, "Imported preprocessed data")

# Create unique cell identifier by combining <filename> and <row_id>
samples["cell_id"] = samples["filename"] + "__" + samples["row_id"].astype(str)




#### DOWNSAMPLE ####

if selective_sampling == True:
    small_df = selective_downsample_pd_dataframe(dataframe=samples, keyword=filename_keyword, sample_size=sample_size)
else:  
    small_df = downsample_pd_dataframe(dataframe=samples, sample_size=sample_size)





#### DIMENSIONALITY REDUCTION ####

small_reduced_df, reduced_columns, reducer = fit_transform_dr(method_choice=method_choice,
                                                              input_df=small_df,
                                                              clustering_markers=clustering_markers,
                                                              markers_to_reduce=markers_to_reduce,
                                                              n_components=n_components,
                                                              n_neighbors=n_neighbors)

print_pd_dataframe(small_reduced_df, "Reduced data, input for clustering")

assert not small_reduced_df.isna().any().any()
assert np.isfinite(small_reduced_df.values).all()




#### RUN BIRCH ####

clustering_start = datetime.now()
print("Clustering start:", clustering_start.strftime("%Y-%m-%d %H:%M:%S"), flush=True)
print(f"Running Birch with threshold={radius_threshold}, branching_factor={branching_factor}, and n_clusters={nb_clusters}... \n", flush=True)

brc = Birch(threshold=radius_threshold, branching_factor=branching_factor, n_clusters=nb_clusters, compute_labels=True)
brc.fit(small_reduced_df)

clustering_end = datetime.now()
elapsed_clustering = str(clustering_end - clustering_start).split(".")[0]
print(f"Total clustering time: {elapsed_clustering}", flush=True)
print("Number of subclusters:", len(brc.subcluster_centers_), flush=True)

# Save model
if method_choice.lower() == "direct_markers":
    brc_model = {
        "date" : datetime.now(),
        "model" : brc.subcluster_centers_,
        "radius_threshold" : radius_threshold,
        "branching_factor" : branching_factor,
        "nb_clusters" : nb_clusters,
        "dimensionality_reduction" : method_choice,
        "clustering_markers" : clustering_markers}
else:
    brc_model = {
        "date" : datetime.now(),
        "model" : brc.subcluster_centers_,
        "radius_threshold" : radius_threshold,
        "branching_factor" : branching_factor,
        "nb_clusters" : nb_clusters,
        "dimensionality_reduction" : method_choice,
        "clustering_markers" : markers_to_reduce}

with open(model_pkl, "wb") as f:
    pickle.dump(brc_model, f)

print_save_msg(filepath=model_pkl, filename="Birch model")


# Get and attach directly obtained labels 
brc_labels = brc.labels_
small_df["cluster"] = pd.Categorical(brc_labels)
samples_with_clusters = attach_direct_cluster_labels(small_df, samples)
print_pd_dataframe(small_df, "Dataframe with cluster results attached")




#### PREDICT LABELS FOR REST OF DATA #####

nb_unlabeled_bef = get_nb_unlabeled_cells(samples_with_clusters)
print("Assigning unlabeled cells to their nearest cluster...", flush=True)

# Dataframe of unlabeled cells
unlabeled_df = samples_with_clusters[samples_with_clusters["cluster"].isna()]

# Same DR method for unlabeled data
unlabeled_df_reduced, reduced_columns = transform_dr(method_choice=method_choice,
                                input_df=unlabeled_df,
                                reducer=reducer,
                                chunksize=chunksize,
                                clustering_markers=clustering_markers,
                                markers_to_reduce=markers_to_reduce,
                                n_components=n_components,
                                n_neighbors=n_neighbors)

# Run birch (predict mode) to get remaining labels
pred_labels = brc.predict(unlabeled_df_reduced)

# Attach labels to unlabeled cells dataframe
unlabeled_df["cluster"] = pd.Categorical(pred_labels)

# Attach predicted labels to complete dataframe
samples_with_clusters.loc[unlabeled_df.index, "cluster"] = unlabeled_df["cluster"]


nb_unlabeled_aft = get_nb_unlabeled_cells(samples_with_clusters)
assert nb_unlabeled_aft == 0, f"{nb_unlabeled_aft} cells still unlabeled after propagation"
print(f"All {nb_unlabeled_bef} unlabeled cells successfully assigned a cluster", flush=True)

print_cluster_size_bef_aft_propagation(small_df, samples_with_clusters, "cluster")





#### SAVE STANDARD CLUSTERING OUTPUT ####

save_standard_clustering_output(samples_with_clusters, clusters_csv)
print_save_msg(filepath=clusters_csv, filename="standard clustering output")


end = end_time()
elapsed_time(start, end)