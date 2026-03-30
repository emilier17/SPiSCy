import pandas as pd
import numpy as np
import phenograph as ph
import yaml
import umap 
import matplotlib.pyplot as plt
import re
import pickle
# Library
import os
import sys
import scvi
from scvi.external import cytovi
import anndata as ad
import scanpy as sc
from sklearn.neighbors import NearestNeighbors



#clusters_csv = "results/clustering/phenograph/clusters.csv"
metadata_csv = "data/metadata.csv"
samples_csv = "results/csv/final_samples.csv"
spe_config = "config/run_cytovi.yaml"
model_pkl = "results/clustering/phenograph/phenograph_model.pkl"
cluster_csv = "results/clustering/cytovi/clusters.csv"


with open(spe_config) as f:
    config_file = yaml.safe_load(f)

clustering_markers = config_file["clustering_markers"]
sample_size = config_file["sample_size_per_file"]



## import

samples = pd.read_csv(samples_csv)

print("[TABLE PREVIEW] Imported preprocessed data:")
print(samples.head())
print("\n")



## downsample

print(f"Sampling {sample_size} from each file")

small_df = (
    samples
    .groupby("filename", group_keys=False)
    .apply(lambda x: x.sample(
        n=min(sample_size, len(x)),
        random_state=17
    ))
    .reset_index(drop=True)
)

unique_filenames = small_df["filename"].unique()
print(f"Number of files sampled: {len(unique_filenames)}")
print(f"Total number of cells in dataset after sampling: {small_df.shape[0]}")
print("\n")



## select markers only

small_clustering_df = small_df[clustering_markers]

assert not small_clustering_df.isna().any().any()
assert np.isfinite(small_clustering_df.values).all()

print("[TABLE PREVIEW] Data that will be used as input to clustering algorithm:")
print(small_clustering_df.head())
print("\n")



## convert to anndata format for CytoVI

X = small_clustering_df.values
adata = ad.AnnData(X=X,
                   obs=small_df[["cell_id", "filename"]].copy(),
                   var=pd.DataFrame(index=clustering_markers)
                   )

assert adata.n_obs == small_df.shape[0]
assert list(adata.var.index) == clustering_markers

print(adata)


## prepare full dataset for inference (same markers, same order)
full_clustering_df = samples[clustering_markers]

assert not full_clustering_df.isna().any().any()
assert np.isfinite(full_clustering_df.values).all()

X_full = full_clustering_df.values

adata_full = ad.AnnData(
    X=X_full,
    obs=samples[["cell_id", "filename"]].copy(),
    var=pd.DataFrame(index=clustering_markers)
)

assert list(adata_full.var.index) == list(adata.var.index)




### LATENT REPRENSETATION

## training (creates latent space)
cytovi.CYTOVI.setup_anndata(adata, sample_key="filename")
model = cytovi.CYTOVI(adata)
model.train(n_epochs_kl_warmup=50)
print(model)

plt.plot(model.history["elbo_train"], label="Train")
plt.plot(model.history["elbo_validation"], label="Validation")
plt.xlabel("Epochs")
plt.ylabel("ELBO")
plt.legend()
plt.title("Training vs Validation ELBO")

## get latent space
latent_small = model.get_latent_representation()
assert latent_small.shape[0] == adata.n_obs
adata.obsm["X_cytovi"] = latent_small


## latent embedding for full dataset (no retraining)
latent_full = model.get_latent_representation(
    adata=adata_full,
    batch_size=1024  # tune based on GPU memory
)

assert latent_full.shape[0] == adata_full.n_obs
adata_full.obsm["X_cytovi"] = latent_full



## CLUSTERING
## compute neighbors
sc.pp.neighbors(adata, n_neighbors=15, use_rep="X_cytovi", transformer="pynndescent")

## running Leiden clustering on the latent space 
sc.tl.leiden(adata, resolution=0.4, key_added="cluster", flavor="igraph")
print(adata.obs["cluster"].value_counts())

sc.tl.umap(adata)
sc.pl.umap(adata, color="cluster")
plt.show()


## propagate cluster labels



## save standard output
cluster_df = adata.obs[["cell_id", "filename", "cluster"]].copy()
cluster_df.to_csv(cluster_csv, index=False)
