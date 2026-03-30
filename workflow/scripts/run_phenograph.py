#==================================================================================================================#
# SPiSCy (Snakemake PIpeline for Spectral CYtometry)
#
# Phase: clustering
# Step: perform clustering with the Phenograph algorithm (KNN + Jaccard Index + community detection)
# Original publication : 10.1016/j.cell.2015.05.047
#
# Author: Émilie Roy
# Date: Jan 2026
# Version: v.1.0
#
# Input : final csv file (.csv)
# Outputs : standard output clusters.csv + phenograph model (.pkl)
#==================================================================================================================#


from clustering_utils import (
    snakemake_logs, start_time, check_pos_int_not_null, end_time, elapsed_time,
    print_pd_dataframe, print_save_msg, downsample_pd_dataframe, selective_downsample_pd_dataframe,
    attach_direct_cluster_labels, get_nb_unlabeled_cells, propagate_labels_by_chunk,
    print_cluster_size_bef_aft_propagation, save_standard_clustering_output, check_boolean,
    check_string, fit_transform_dr, check_list_strings 
    )



snakemake_logs(snakemake)
start = start_time()



#### IMPORTS AND SNAKEMAKE VARIABLES #####

# Library
import sys
import os
import pandas as pd
import numpy as np
import phenograph as ph
import yaml
from sklearn.neighbors import NearestNeighbors
from multiprocessing import freeze_support
import pickle
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


# Config files
with open(clustering_config) as f:
    config_file = yaml.safe_load(f)

k = config_file["k"]
dist_metric = config_file["distance_metric"]
community_detection = config_file["clustering_algo"]
sample_size = config_file["sample_size_per_file"]

if dist_metric not in ['euclidean', 'manhattan', 'correlation', 'cosine']:
    sys.exit(f"In {os.path.basename(clustering_config)}: {dist_metric} is unknown. Distance metric options are 'euclidean', 'manhattan', 'correlation', or 'cosine'")
if community_detection not in ['louvain', 'leiden']:
    sys.exit(f"In {os.path.basename(clustering_config)}: {community_detection} is unknown. Community detection algorithm options are 'louvain' or 'leiden'")
check_pos_int_not_null(k, "k", clustering_config)
check_pos_int_not_null(sample_size, "sample_size", clustering_config)




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

nb_neighbors = config_file["nb_neighbors"]
chunksize = config_file["chunksize"]

check_pos_int_not_null(nb_neighbors, "nb_neighbors", label_prop_config)
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




def main():

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




    #### RUN PHENOGRAPH ####
    
    clustering_start = datetime.now()
    print("Clustering start:", clustering_start.strftime("%Y-%m-%d %H:%M:%S"), flush=True)
    print(f"Running PhenoGraph with k={k}, distance metric={dist_metric}, and clustering_algo={community_detection}... \n", flush=True)
    
    communities, graph, Q = ph.cluster(
        small_reduced_df,
        k = k,
        primary_metric = dist_metric,
        clustering_algo = community_detection,
        seed = 17
    )

    clustering_end = datetime.now()
    elapsed_clustering = str(clustering_end - clustering_start).split(".")[0]
    print(f"Total clustering time: {elapsed_clustering}", flush=True)


    # Save model
    if method_choice.lower() == "direct_markers":
        phenograph_model = {
            "date" : datetime.now(),
            "k" : k,
            "metric" : dist_metric,
            "clustering_aglo" : community_detection,
            "graph" : graph,
            "Q" : Q,
            "communities" : communities,
            "dimensionality_reduction" : method_choice,
            "clustering_markers" : clustering_markers}
    else:
        phenograph_model = {
            "date" : datetime.now(),
            "k" : k,
            "metric" : dist_metric,
            "clustering_aglo" : community_detection,
            "graph" : graph,
            "Q" : Q,
            "communities" : communities,
            "dimensionality_reduction" : method_choice,
            "clustering_markers" : markers_to_reduce}
    
    with open(model_pkl, "wb") as f:
        pickle.dump(phenograph_model, f)

    print_save_msg(filepath=model_pkl, filename="PhenoGraph model")




    #### ATTACH CLUSTERING RESULTS ####

    # Add cluster number (community) to the downsampled dataframe
    small_df["cluster"] = pd.Categorical(communities)

    # Add the clustering results to the main samples dataframe
    samples_with_clusters = attach_direct_cluster_labels(small_df, samples)

    print_pd_dataframe(small_df, "Dataframe with cluster results attached")





    #### ASSIGN UNLABELED (NON-CLUSTERED) CELLS TO A CLUSTER ####

    nb_unlabeled_bef = get_nb_unlabeled_cells(samples_with_clusters)
    print("Assigning unlabeled cells to their nearest cluster...", flush=True)

    # Prepare training data : clustered cells
    X_train = small_reduced_df.values # marker values
    y_train = small_df["cluster"].values # cluster labels

    # Fit NearestNeighbors on cells with a cluster assigned (label)
    nbrs = NearestNeighbors(n_neighbors=nb_neighbors, n_jobs=-1)
    nbrs.fit(X_train)

    # Transform and propagate labels chunk by chunk
    samples_with_clusters = propagate_labels_by_chunk(full_df_with_clusters=samples_with_clusters,
                                                    chunksize=chunksize, 
                                                    nearest_neigh_model=nbrs,
                                                    cluster_labels=y_train,
                                                    clustering_markers=clustering_markers,
                                                    method_choice=method_choice,
                                                    reducer=reducer,
                                                    markers_to_reduce=markers_to_reduce)

    nb_unlabeled_aft = get_nb_unlabeled_cells(samples_with_clusters)
    assert nb_unlabeled_aft == 0, f"{nb_unlabeled_aft} cells still unlabeled after propagation"
    print(f"All {nb_unlabeled_bef} unlabeled cells successfully assigned a cluster")

    print_cluster_size_bef_aft_propagation(small_df, samples_with_clusters, "cluster")



    #### SAVE STANDARD CLUSTERING OUTPUT ####

    save_standard_clustering_output(samples_with_clusters, clusters_csv)
    print_save_msg(filepath=clusters_csv, filename="standard clustering output")



if __name__ == "__main__":
    freeze_support()
    main()
    end = end_time()
    elapsed_time(start, end)
