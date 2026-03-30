#==================================================================================================================#
# SPiSCy (Snakemake PIpeline for Spectral CYtometry)
#
# Phase: clustering
# Step: perform clustering with the FlowSOM algorithm (self-organized maps)
# Original publicatipon : https://doi.org/10.1093/bioinformatics/btae179
#
# Author: Émilie Roy
# Date: Jan 2026
# Version: v.1.0
#
# Input : final csv file (.csv)
# Outputs : standard output clusters.csv + model (.pkl)
#==================================================================================================================#


# Redirect prints and errors into log files (stdout and stderr)

from clustering_utils import (
    snakemake_logs, start_time, check_pos_int_not_null, end_time, elapsed_time,
    print_pd_dataframe, print_save_msg, downsample_pd_dataframe, get_nb_unlabeled_cells,
    print_cluster_size_bef_aft_propagation, save_standard_clustering_output, pd_df_to_anndata, 
    check_list_strings, fit_transform_dr, transform_dr, check_boolean, check_string,
    selective_downsample_pd_dataframe
    )



snakemake_logs(snakemake)
start = start_time()



#### IMPORTS AND SNAKEMAKE VARIABLES #####

# Library
import os
import pandas as pd
import numpy as np
from datetime import datetime
import yaml
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import pickle
import flowsom as fs


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
summary_pdf = snakemake.output[2]


## clustering config file
with open(clustering_config) as f:
    config_file = yaml.safe_load(f)

sample_size = config_file["sample_size_per_file"]
x_dims = config_file["x_dims"]
y_dims = config_file["y_dims"]
nb_metaclusters = config_file["nb_metaclusters"]

check_pos_int_not_null(sample_size, "sample_size", clustering_config)
check_pos_int_not_null(x_dims, "x_dims", clustering_config)
check_pos_int_not_null(y_dims, "y_dims", clustering_config)
check_pos_int_not_null(nb_metaclusters, "nb_metaclusters", clustering_config)


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




#### CONVERT TO ANNDATA FOR FLOWSOM ####

adata = pd_df_to_anndata(X=small_reduced_df.values,
                         df=small_df,
                         var_names=reduced_columns)

assert adata.n_obs == small_df.shape[0]
assert list(adata.var.index) == reduced_columns





#### RUN FLOWSOM ####

# Run clustering
clustering_start = datetime.now()
print("Clustering start:", clustering_start.strftime("%Y-%m-%d %H:%M:%S"), flush=True)
print(f"Running FlowSOM with xdim={x_dims}, ydim={y_dims}, n_clusters={x_dims*y_dims}, and n_metaclusters={nb_metaclusters}... \n", flush=True)

fsom = fs.FlowSOM(adata,
                cols_to_use=reduced_columns,
                xdim=x_dims,
                ydim=y_dims,
                n_clusters=nb_metaclusters,
                seed=17)

clustering_end = datetime.now()
elapsed_clustering = str(clustering_end - clustering_start).split(".")[0]
print(f"Total clustering time: {elapsed_clustering}", flush=True)



# Save model
if method_choice.lower() == "direct_markers":
    flowsom_model = {
        "date" : datetime.now(),
        "cluster_data" : fsom.get_cluster_data(),
        "x_dims" : x_dims,
        "y_dims" : y_dims,
        "nb_metaclusters" : nb_metaclusters,
        "dimensionality_reduction" : method_choice,
        "clustering_markers" : clustering_markers}
else:
    flowsom_model = {
        "date" : datetime.now(),
        "cluster_data" : fsom.get_cluster_data(),
        "x_dims" : x_dims,
        "y_dims" : y_dims,
        "nb_metaclusters" : nb_metaclusters,
        "dimensionality_reduction" : method_choice,
        "clustering_markers" : markers_to_reduce}

with open(model_pkl, "wb") as f:
    pickle.dump(flowsom_model, f)

print_save_msg(filepath=model_pkl, filename="FlowSOM model")



#### ATTACH CLUSTERING RESULTS ####

# get the clustering results
adata.obs["FlowSOM_clusters"] = fsom.get_cell_data().obs["clustering"]
adata.obs["FlowSOM_metaclusters"] = fsom.get_cell_data().obs["metaclustering"]

# attach results to pandas dataframe
small_df["cluster"] = pd.Categorical(adata.obs["FlowSOM_clusters"].values)
small_df["metacluster"] = pd.Categorical(adata.obs["FlowSOM_metaclusters"].values)

# Add the clustering results to the main samples dataframe
samples_with_clusters = samples.merge(
    small_df[["cell_id", "cluster", "metacluster"]],
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


print_pd_dataframe(small_df, "Dataframe with cluster results attached")






#### ASSIGN UNSAMPLED CELLS TO CLUSTERS ####

nb_unlabeled_bef = get_nb_unlabeled_cells(samples_with_clusters)
print("Assigning unlabeled cells to their nearest cluster...")

# Get unlabeled cells
unlabeled_df = samples_with_clusters[samples_with_clusters["cluster"].isna()]

# apply same DR method
unlabeled_df_reduced, reduced_columns = transform_dr(method_choice=method_choice,
                                input_df=unlabeled_df,
                                reducer=reducer,
                                chunksize=chunksize,
                                clustering_markers=clustering_markers,
                                markers_to_reduce=markers_to_reduce,
                                n_components=n_components,
                                n_neighbors=n_neighbors)


# Convert to anndata format for flowsom
adata_unlabeled = pd_df_to_anndata(X=unlabeled_df_reduced.values,
                                   df=unlabeled_df,
                                   var_names=reduced_columns)


assert adata_unlabeled.n_obs == unlabeled_df.shape[0]
assert list(adata_unlabeled.var.index) == reduced_columns

# Use FlowSOM's function new_data to map new cells to an existing FlowSOM grid
fsom_unlabeled = fsom.new_data(adata_unlabeled)

# get the clustering results and attach to pandas dataframe
prop = fsom_unlabeled.get_cell_data().obs[
    ["clustering", "metaclustering"]
].copy()

prop["cell_id"] = adata_unlabeled.obs["cell_id"].values

samples_with_clusters = samples_with_clusters.merge(
    prop.rename(columns={
        "clustering": "cluster_prop",
        "metaclustering": "metacluster_prop"
    }),
    on="cell_id",
    how="left"
)

# allow new categories (new cluster IDs) to be added
samples_with_clusters["cluster"] = samples_with_clusters["cluster"].cat.add_categories(
    np.setdiff1d(prop["clustering"].unique(), samples_with_clusters["cluster"].cat.categories))

samples_with_clusters["cluster"] = samples_with_clusters["cluster"].fillna(
    samples_with_clusters["cluster_prop"])

samples_with_clusters["metacluster"] = samples_with_clusters["metacluster"].cat.add_categories(
    np.setdiff1d(prop["metaclustering"].unique(), samples_with_clusters["metacluster"].cat.categories))

samples_with_clusters["metacluster"] = samples_with_clusters["metacluster"].fillna(
    samples_with_clusters["metacluster_prop"])

nb_unlabeled_aft = get_nb_unlabeled_cells(samples_with_clusters)
assert nb_unlabeled_aft == 0, f"{nb_unlabeled_aft} cells still unlabeled after propagation"
print(f"All {nb_unlabeled_bef} unlabeled cells successfully assigned a cluster")

print_cluster_size_bef_aft_propagation(small_df, samples_with_clusters, "metacluster")
print_cluster_size_bef_aft_propagation(small_df, samples_with_clusters, "cluster")




#### SAVE STANDARD OUTPUT ####

save_standard_clustering_output(samples_with_clusters, clusters_csv, "metacluster")
print_save_msg(filepath=clusters_csv, filename="standard clustering output")


#### FLOWSOM SPECIFIC VISUALIZATIONS ####

print("Ploting FlowSOM specific visualizations")

# 3 types of plots
plots = [
    ("FlowSOM tree", lambda: fs.pl.plot_stars(
        fsom,
        background_values=fsom.get_cluster_data().obs["metaclustering"]
        )),
    ("FlowSOM nodes (clusters)", lambda: fs.pl.plot_stars(
        fsom,
        background_values=fsom.get_cluster_data().obs.metaclustering,
        view="grid",
        equal_node_size=True,
        equal_background_size=True,
        )),
    ("FlowSOM cluster numbers on tree", lambda: fs.pl.plot_numbers(
        fsom,
        level="clusters",
        text_size=5
        ))
]

# add in an expression plot for every marker
if method_choice=="direct_markers":
    for marker in clustering_markers:
        plots.append(
            (
                f"Marker expression: {marker}",
                lambda m=marker: fs.pl.plot_marker(
                    fsom,
                    marker=np.array([m])
                )
            )
        )


# save each plot to a page in the summary pdf
with PdfPages(summary_pdf) as pdf:
    for title, plot_fn in plots:
        plot_fn()
        fig = plt.gcf()
        fig.suptitle(title, fontsize=8)
        pdf.savefig(fig, bbox_inches="tight")
        plt.close(fig)




end = end_time()
elapsed_time(start, end)
