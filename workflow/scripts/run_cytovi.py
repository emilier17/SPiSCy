#==================================================================================================================#
# SPiSCy (Snakemake PIpeline for Spectral CYtometry)
#
# Phase: clustering
# Step: perform clustering with the CytoVI algorithm from scvi-tools
# Original publication : 10.1038/s41587-021-01206-w.
#
# Author: Émilie Roy
# Date: Jan 2026
# Version: v.1.0
#
# Input : final csv file (.csv)
# Outputs : standard output clusters.csv + cytovi model (.pkl)
#==================================================================================================================#


from clustering_utils import (
    snakemake_logs, start_time, check_pos_int_not_null, end_time, elapsed_time,
    print_pd_dataframe, print_save_msg, downsample_pd_dataframe, attach_direct_cluster_labels,
    get_nb_unlabeled_cells, propagate_labels_by_chunk_cytovi, print_cluster_size_bef_aft_propagation,
    save_standard_clustering_output, pd_df_to_anndata, check_string, 
    selective_downsample_pd_dataframe, check_boolean
    )

snakemake_logs(snakemake)
start = start_time()




#### IMPORTS AND SNAKEMAKE VARIABLES #####

# Library
import sys
import os
import pandas as pd
import numpy as np
import yaml
import pickle
import matplotlib.pyplot as plt
import scanpy as sc
from scvi.external import cytovi
from sklearn.neighbors import NearestNeighbors
from datetime import datetime


# Snakemake variables
# Inputs
samples_csv = snakemake.input[0]
clustering_config = snakemake.input[1]
sampling_config = snakemake.input[2]
label_prop_config = snakemake.input[3]
metadata_csv = snakemake.input[4]
# Outputs
model_pkl = snakemake.output[0]
clusters_csv = snakemake.output[1]
elbo_png = snakemake.output[2]


# Config files
with open(clustering_config) as f:
    config_file = yaml.safe_load(f)

clustering_markers = config_file["clustering_markers"]
sample_size = config_file["sample_size_per_file"]
nb_neighbors_clustering = config_file["nb_neighbors_clustering"]
sample_key = config_file["sample_key"]

# Verify that config choices are valid
if len(clustering_markers) == 0:
    sys.exit(f"In {os.path.basename(spe_config)}: at least 1 marker must be choosen for clustering")
check_pos_int_not_null(sample_size, "sample_size", clustering_config)
check_pos_int_not_null(nb_neighbors_clustering, "nb_neighbors_clustering", clustering_config)
check_string(sample_key, "sample_key", clustering_config)
if not sample_key:
    sample_key = None



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

nb_neighbors_propagation = config_file["nb_neighbors"]
chunksize = config_file["chunksize"]

check_pos_int_not_null(nb_neighbors_propagation, "nb_neighbors", label_prop_config)
check_pos_int_not_null(chunksize, "chunksize", label_prop_config)




#### IMPORT PREPROCESSED DATA ####

print("Importing preprocessed data...", flush=True)

# Only import the columns we need
samples = pd.read_csv(samples_csv, usecols=["filename", "row_id"] + clustering_markers)
print_pd_dataframe(samples, "Imported preprocessed data")

# Create unique cell identifier by combining <filename> and <row_id>
samples["cell_id"] = samples["filename"] + "__" + samples["row_id"].astype(str)

# if sample_key is provided, merge metadata
if sample_key:
    metadata = pd.read_csv(metadata_csv, usecols=[sample_key, "filename"])

    samples = samples.merge(
        metadata[[sample_key, "filename"]],
        on="filename",
        how="left",
        validate="many_to_one"
    )
    
    print_pd_dataframe(samples, f"Metadata merge with sample_key={sample_key}")





#### DOWNSAMPLE ####

if selective_sampling == True:
    small_df = selective_downsample_pd_dataframe(dataframe=samples, keyword=filename_keyword, sample_size=sample_size)
else:  
    small_df = downsample_pd_dataframe(dataframe=samples, sample_size=sample_size)





#### SELECT ONLY CLUSTERING MARKERS ####

small_clustering_df = small_df[clustering_markers]

assert not small_clustering_df.isna().any().any()
assert np.isfinite(small_clustering_df.values).all()

print_pd_dataframe(small_clustering_df, "Input data for clustering")




#### CONVERT TO ANNDATA FOR CYTOVI ####

adata_small = pd_df_to_anndata(X=small_clustering_df.values,
                        df=small_df,
                        var_names=clustering_markers, 
                        sample_key=sample_key)

assert adata_small.n_obs == small_df.shape[0]
assert list(adata_small.var.index) == clustering_markers




#### RUN CYTOVI: GET LATENT REPRESENTATION OF DATA ####

# training (creates latent space)
cytovi.CYTOVI.setup_anndata(adata_small, sample_key=sample_key)
model = cytovi.CYTOVI(adata_small)

clustering_start = datetime.now()
print("Training start:", clustering_start.strftime("%Y-%m-%d %H:%M:%S"), flush=True)
model.train(n_epochs_kl_warmup=50)
clustering_end = datetime.now()
elapsed_clustering = str(clustering_end - clustering_start).split(".")[0]
print(f"Total training time: {elapsed_clustering}", flush=True)
print(model)

# Save model
cytovi_model = {
    "date" : datetime.now(),
    "model" : model,
    "clustering_markers" : clustering_markers}

with open(model_pkl, "wb") as f:
    pickle.dump(cytovi_model, f)

print_save_msg(model_pkl, "CytoVI model")

# save elbo plot
plt.plot(model.history["elbo_train"], label="Train")
plt.plot(model.history["elbo_validation"], label="Validation")
plt.xlabel("Epochs")
plt.ylabel("ELBO")
plt.legend()
plt.title("Training vs Validation ELBO")
plt.savefig(elbo_png)
print_save_msg(elbo_png, "ELBO plot")

# get latent space
latent_small = model.get_latent_representation(adata=adata_small, batch_size=1024)
assert latent_small.shape[0] == adata_small.n_obs
adata_small.obsm["X_cytovi"] = latent_small

# prepare full dataset for mapping onto latent space
full_clustering_df = samples[clustering_markers]

assert not full_clustering_df.isna().any().any()
assert np.isfinite(full_clustering_df.values).all()

adata_full = pd_df_to_anndata(X=full_clustering_df.values,
                                df=samples,
                                var_names=clustering_markers, 
                                sample_key=sample_key)

assert list(adata_full.var.index) == list(adata_small.var.index)

# Map full dataset onto latent space
latent_full = model.get_latent_representation(
    adata=adata_full,
    batch_size=1024
)

assert latent_full.shape[0] == adata_full.n_obs
adata_full.obsm["X_cytovi"] = latent_full



#### PERFORM CLUSTERING ON THE LATENT REPRESENTATION ####

# Compute neighbors
sc.pp.neighbors(adata_small, n_neighbors=nb_neighbors_clustering, use_rep="X_cytovi", transformer="pynndescent")

# Run Leiden community detection
sc.tl.leiden(adata_small, resolution=0.4, key_added="cluster", flavor="igraph")
print(adata_small.obs["cluster"].value_counts(), flush=True)




#### ATTACH CLUSTERING RESULTS ####

# Add cluster number to the downsampled dataframe
small_df["cluster"] = pd.Categorical(adata_small.obs["cluster"].values)

# Add the clustering results to the main samples dataframe
samples_with_clusters = attach_direct_cluster_labels(small_df, samples)

print_pd_dataframe(small_df, "Dataframe with cluster results attached")





#### ASSIGN UNLABELED (NON-CLUSTERED) CELLS TO A CLUSTER ####

nb_unlabeled_bef = get_nb_unlabeled_cells(samples_with_clusters)
print("Assigning unlabeled cells to their nearest cluster...", flush=True)

# Prepare training data : clustered cells
X_train = adata_small.obsm["X_cytovi"] # CytoVI latent space
y_train = small_df["cluster"].values # cluster labels

# Fit NearestNeighbors on cells with a cluster assigned (label)
nbrs = NearestNeighbors(n_neighbors=nb_neighbors_propagation, n_jobs=-1)
nbrs.fit(X_train)

# Propagate labels chunk by chunk (n cells at a time)
samples_with_clusters = propagate_labels_by_chunk_cytovi(full_df_with_clusters = samples_with_clusters,
                                                        adata_full = adata_full,
                                                        chunksize = chunksize,
                                                        nearest_neigh_model = nbrs,
                                                        cluster_labels=y_train)

nb_unlabeled_aft = get_nb_unlabeled_cells(samples_with_clusters)
assert nb_unlabeled_aft == 0, f"{nb_unlabeled_aft} cells still unlabeled after propagation"
print(f"All {nb_unlabeled_bef} unlabeled cells successfully assigned a cluster", flush=True)

print_cluster_size_bef_aft_propagation(small_df, samples_with_clusters, "cluster")




#### SAVE STANDARD CLUSTERING OUTPUT ####

save_standard_clustering_output(samples_with_clusters, clusters_csv)
print_save_msg(filepath=clusters_csv, filename="standard clustering output")


end = end_time()
elapsed_time(start, end)