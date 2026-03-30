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
    snakemake_logs, start_time, end_time, elapsed_time, compute_clustering_metrics, 
    check_pos_int_not_null
    )

snakemake_logs(snakemake)
start = start_time()



##########################################
#### IMPORTS AND SNAKEMAKE VARIABLES #####
##########################################

# Library
import os
import sys
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import yaml
from sklearn.metrics import adjusted_rand_score, adjusted_mutual_info_score
from sklearn.metrics.pairwise import cosine_similarity
from sklearn.decomposition import PCA
from matplotlib.backends.backend_pdf import PdfPages
from itertools import combinations


# Snakemake variables
# Inputs
samples_csv = snakemake.input[0]
spe_config = snakemake.input[1]
clusters_csvs = snakemake.input["clusters_csvs"]
marker_median_csvs = snakemake.input["median_marker_csvs"]

# Outputs
summary_pdf = snakemake.output[0]

# Config files
with open(spe_config) as f:
    config_file = yaml.safe_load(f)

clustering_markers = config_file["clustering_markers"]
sample_size_agreement = config_file["sample_size_agreement"]
sample_size_silhouette = config_file["sample_size_silhouette"]

if len(clustering_markers) == 0:
    sys.exit(f"In {os.path.basename(spe_config)}: at least 1 marker must be chosen")
check_pos_int_not_null(sample_size_agreement, "sample_size_agreement", spe_config)
check_pos_int_not_null(sample_size_silhouette, "sample_size_silhouette", spe_config)





#################
#### IMPORTS ####
#################

# associate each clusters.csv to its clustering method+level, and store in dictionnary
clustering_results = {}
non_cluster_cols = {"filename", "row_id", "label_propagation"}

for csv in clusters_csvs:
    clusters = pd.read_csv(csv)
    method_dir = str(os.path.dirname(csv))
    base_method = method_dir.split("/")[-1]

    # detect clustering columns dynamically
    clustering_columns = [col for col in clusters.columns if col not in non_cluster_cols]

    for level in clustering_columns:
        method_name = f"{base_method}_{level}"
        # keep only relevant columns for this level
        df_level = clusters[["row_id", level]].rename(columns={level: "cluster"})
        clustering_results[method_name] = df_level

print(f"Clustering method + clustering granularity (level) found: {clustering_results.keys()}")





#######################
#### CLUSTER SIZES ####
#######################

#-- Compare number of cluster and cluster size distributions for each method

# number of clusters and cluster sizes
print("Plotting cluster sizes...", flush=True)

cluster_size_rows = []

for method, df in clustering_results.items():
    counts = df["cluster"].value_counts()

    for cluster_id, size in counts.items():
        cluster_size_rows.append({
            "method": method,
            "cluster_nb": cluster_id,
            "cluster_size": size
        })

cluster_sizes_df = pd.DataFrame(cluster_size_rows)






###################
#### AGREEMENT ####
###################

#-- See if different methods assigned the same cells to the same clusters
#-- Two measures:
# Adjusted Rand Index. [-0.5 to 1]. 1 indicates identical assignment 
# Adjusted Mutual Information. [0 to 1]. 1 indicates identical assignment

print("Clustering agreement: calculating ARI and AMI...", flush=True)

# Calculate ARI and AMI score for each combination of methods (pairwise comparaison)
methods = sorted(clustering_results.keys())
method_pairs = combinations(methods, 2)

# initialize empty square matrices
ari_matrix = pd.DataFrame(np.nan, index=methods, columns=methods)
ami_matrix = pd.DataFrame(np.nan, index=methods, columns=methods)

# diagonal = perfect agreement
np.fill_diagonal(ari_matrix.values, 1.0)
np.fill_diagonal(ami_matrix.values, 1.0)

for m1, m2 in method_pairs:
    labels1 = clustering_results[m1]["cluster"].to_numpy()
    labels2 = clustering_results[m2]["cluster"].to_numpy()

    n_cells = min(len(labels1), len(labels2), sample_size_agreement)
    idx = np.random.choice(len(labels1), size=n_cells, replace=False)

    ari = adjusted_rand_score(labels1[idx], labels2[idx])
    ami = adjusted_mutual_info_score(labels1[idx], labels2[idx])

    # fill symmetric positions
    ari_matrix.loc[m1, m2] = ari
    ari_matrix.loc[m2, m1] = ari

    ami_matrix.loc[m1, m2] = ami
    ami_matrix.loc[m2, m1] = ami





############################
#### STRUCTURAL METRICS ####
############################

#-- Compare the clustering structural quality of each method
#-- Three scores (more details https://scikit-learn.org/stable/modules/clustering.html#clustering-performance-evaluation):
# Silhouette Coefficient: how dense and well separated are clusters. [-1 to +1]. Higher is better
# Calinski-Harabasz Index: how dense and well separated are clusters. Higher is better
# Davies-Bouldin: how well separated clusters are. 0 is lowest score. Lower is better

print("Structural metrics: calculating Silhouette Coefficient, Calinski-Harabasz Index, and Davies-Bouldin score...", flush=True)

samples = pd.read_csv(samples_csv)
X = samples[clustering_markers]

rows = []

for method, results in clustering_results.items():
    labels = results["cluster"].to_numpy()
    sil, ch, db = compute_clustering_metrics(X, labels, sample_size=sample_size_silhouette)
    rows.append({
        "method": method,
        "silhouette": sil,
        "calinski_harabasz": ch,
        "davies_bouldin": db,
    })

metrics_df = pd.DataFrame(rows).set_index("method")
raw_metrics_df = metrics_df.copy()

# inverse davies bouldin metric (OG: smaller is better) to match other metrics evaluation (bigger is better)
metrics_df["davies_bouldin_inv"] = 1 / metrics_df["davies_bouldin"]
metrics_df = metrics_df.drop(columns="davies_bouldin")

# min-max scaling for each metric (for heatmap)
scaled = (metrics_df - metrics_df.min()) / (metrics_df.max() - metrics_df.min())
metrics_heatmap = sns.heatmap(scaled)




#################
#### BIOLOGY ####
#################

#-- Do the different clustering algorithms find similar cell populations (based on marker patterns) ?
#-- Two ways to answer this:
        # Measure cosine similarity between median marker values of each cluster. View as heatmap
        # Plot 2 PCA components of median marker values of each cluster. View as scatterplot

print("Biology evaluation: calculating cosine similarity and plotting PCA...", flush=True)

# Import all csvs that contains median marker values for each cluster of each method
all_medians_list = []

for csv in marker_median_csvs:
    medians_df = pd.read_csv(csv)

    method_dir = str(os.path.dirname(csv))
    base_method = method_dir.split("/")[-3]
    method_level = method_dir.split("/")[-1]
    method_name = base_method + "_" + method_level

    medians_df = medians_df.rename({medians_df.columns[0]:"cluster"}, axis="columns")
    medians_df.insert(0, "cluster_id", method_name + "_" + medians_df["cluster"].astype(str))
    medians_df = medians_df.drop(columns="cluster")
    all_medians_list.append(medians_df)

all_medians_df = pd.concat(all_medians_list)
all_medians_df.reset_index(inplace=True, drop=True)


# separate cluster_ids and features
cluster_ids = all_medians_df["cluster_id"]
marker_cols = all_medians_df.columns.drop("cluster_id")
X_all = all_medians_df[marker_cols]
X_cluster = X_all[clustering_markers]

# remove all flowsom_clusters_x, for cosine similarity heatmap
no_flowsom_clusters = all_medians_df[~all_medians_df["cluster_id"].str.startswith("flowsom_cluster_")].copy()
cluster_ids_filtered = no_flowsom_clusters["cluster_id"]
X_all_cos_sim = no_flowsom_clusters[marker_cols]
X_cluster_cos_sim = X_all_cos_sim[clustering_markers]

# Calculate cosine correlation between clusters (don't use flowsom_cluster - too many)
cosine_sims_cluster = pd.DataFrame(
    cosine_similarity(X_cluster_cos_sim),
    index=cluster_ids_filtered,
    columns=cluster_ids_filtered
)
cosine_sims_all = pd.DataFrame(
    cosine_similarity(X_all_cos_sim),
    index=cluster_ids_filtered,
    columns=cluster_ids_filtered
)

# PCA on the median markers and plotting top 2 components
pca = PCA(n_components=2)
X_all_pca = pca.fit_transform(X_all)
X_cluster_pca = pca.fit_transform(X_cluster)

all_medians_df["PC1_all"] = X_all_pca[:,0]
all_medians_df["PC2_all"] = X_all_pca[:,1]
all_medians_df["PC1_cluster"] = X_cluster_pca[:,0]
all_medians_df["PC2_cluster"] = X_cluster_pca[:,1]

# Merge cluster sizes into all_medians_df
cluster_sizes_df["cluster_id"] = cluster_sizes_df["method"] + "_" + cluster_sizes_df["cluster_nb"].astype(str)
all_medians_df = all_medians_df.merge(
    cluster_sizes_df[["cluster_id", "cluster_size"]],
    on="cluster_id",
    how="left"
)

# max and min dot sizes in plots
size_min = 10
size_max = 200
sizes = all_medians_df["cluster_size"]

# min-max linear scaling
scaled_sizes = size_min + (sizes - sizes.min()) / (sizes.max() - sizes.min()) * (size_max - size_min)
all_medians_df["method_level"] = all_medians_df["cluster_id"].apply(lambda x: "_".join(x.split("_")[:2]))





###################################
#### PLOTTING AND GENERATE PDF ####
###################################

with PdfPages(summary_pdf) as pdf:
    #-- Page 1: Heatmap of normalized metrics
    fig, ax = plt.subplots(figsize=(6, 4))
    sns.heatmap(scaled, ax=ax, cmap="viridis", annot=True, fmt=".2f")
    ax.set_title("Normalized Clustering Metrics", fontsize=10)
    ax.set_xlabel("Metric")
    ax.set_ylabel("Clustering Method")
    ax.tick_params(axis="x", rotation=45)
    pdf.savefig(fig, bbox_inches="tight")
    plt.close(fig)

    #-- Page 2: Barplots with raw metric values
    for metric in raw_metrics_df.columns:
        print(metric)
        fig, ax = plt.subplots(figsize=(6, 4))
        sns.barplot(
            x=raw_metrics_df.index,
            y=raw_metrics_df[metric],
            ax=ax,
            palette="Set2"
        )
        ax.set_title(f"{metric.replace('_',' ').title()} per Method", fontsize=10)
        ax.set_ylabel("Score")
        ax.set_xlabel("Method")
        ax.tick_params(axis="x", rotation=45)
        pdf.savefig(fig, bbox_inches="tight")
        plt.close(fig)
    
    #-- Page 3: ARI heatmap
    fig, ax = plt.subplots(figsize=(6, 4))
    sns.heatmap(ari_matrix, ax=ax, cmap="viridis", annot=True, fmt=".2f")
    ax.set_title("Adjusted Rand Index", fontsize=10)
    pdf.savefig(fig, bbox_inches="tight")
    plt.close(fig)

    #-- Page 4: AMI heatmap
    fig, ax = plt.subplots(figsize=(6, 4))
    sns.heatmap(ami_matrix, ax=ax, cmap="viridis", annot=True, fmt=".2f")
    ax.set_title("Adjusted Mutual Information", fontsize=10)
    pdf.savefig(fig, bbox_inches="tight")
    plt.close(fig)
    fig, ax = plt.subplots(figsize=(8,5))

    #-- Page 5: number of clusters and cluster sizes
    sns.boxplot(
        data=cluster_sizes_df,
        x="method",
        y="cluster_size",
        ax=ax,
        showfliers=False,
        color="paleturquoise"
    )

    # jittered points = individual clusters
    sns.stripplot(
        data=cluster_sizes_df,
        x="method",
        y="cluster_size",
        ax=ax,
        color="black",
        alpha=0.4,
        size=3
    )

    ax.set_title("Cluster Size Distribution Across Methods", fontsize=10)
    ax.set_xlabel("Clustering Method")
    ax.set_ylabel("Cluster Size (# cells)")
    ax.set_yscale("log")
    ax.tick_params(axis="x", rotation=45)

    n_clusters_per_method = (
    cluster_sizes_df.groupby("method")["cluster_id"]
    .nunique()
    .reindex(cluster_sizes_df["method"].unique())
    )

    ymax = cluster_sizes_df["cluster_size"].max()
    ax.set_ylim(top=ymax * 2)

    for i, (method, n_clusters) in enumerate(n_clusters_per_method.items()):
        ax.text(
            i,
            ymax * 1.2,
            f"n={n_clusters}",
            ha="center",
            va="bottom",
            fontsize=8
        )

    plt.tight_layout()
    pdf.savefig(fig, bbox_inches="tight")
    plt.close(fig)


    #-- Page 6: cluster expression patterns heatmap
    n = cosine_sims_all.shape[0]
    figsize = (max(8, n * 0.35), max(6, n * 0.35))
    fig, ax = plt.subplots(figsize=figsize)
    sns.heatmap(cosine_sims_all, ax=ax, cmap="viridis", fmt=".2f", cbar_kws={"shrink": 0.5})
    ax.set(xlabel="", ylabel="")
    ax.set_title("Cosine similarity scores between median marker expression (all markers)", fontsize=10)
    pdf.savefig(fig, bbox_inches="tight")
    plt.close(fig)

    #-- Page 7: cluster expression patterns heatmap
    n = cosine_sims_cluster.shape[0]
    figsize = (max(8, n * 0.35), max(6, n * 0.35))
    fig, ax = plt.subplots(figsize=figsize)
    sns.heatmap(cosine_sims_cluster, ax=ax, cmap="viridis", fmt=".2f", cbar_kws={"shrink": 0.5})
    ax.set(xlabel="", ylabel="")
    ax.set_title("Cosine similarity scores between median marker expression (clustering markers)", fontsize=10)
    pdf.savefig(fig, bbox_inches="tight")
    plt.close(fig)


    #-- Page 8 : PCA median marker expression
    fig, ax = plt.subplots(figsize=(7, 5))
    method_levels = all_medians_df["method_level"].unique()
    colors = sns.color_palette("tab10", n_colors=len(method_levels))  # assign a color to each method_level

    # scatter points by method_level
    for method_level, color in zip(method_levels, colors):
        subset = all_medians_df[all_medians_df["method_level"] == method_level]
        subset_sizes = scaled_sizes[subset.index]  # actual sizes

        ax.scatter(
            subset["PC1_all"],
            subset["PC2_all"],
            s=subset_sizes,
            color=color,
            label=method_level,
            alpha=0.7
        )

        # annotate cluster numbers
        for _, row in subset.iterrows():
            cluster_num = row["cluster_id"].split("_")[-1]
            ax.text(
                row["PC1_all"],
                row["PC2_all"],
                cluster_num,
                fontsize=2,
                alpha=0.8,
                ha="center",
                va="center"
            )

    ax.set_xlabel("PC1")
    ax.set_ylabel("PC2")
    ax.set_title("PCA of cluster median marker expression (all markers)")

    ## Legend for method_level (color)
    # create dummy points of fixed size
    color_handles = [plt.scatter([], [], color=color, s=50) for color in colors]
    color_labels = method_levels
    color_legend = ax.legend(color_handles, color_labels, title="Method / Level",
                            bbox_to_anchor=(1.05, 1), loc="upper left", fontsize=6,
                            title_fontsize=8, labelspacing=1)
    ax.add_artist(color_legend)

    ## Legend for cluster size
    size_values = [size_min, round(size_max/2), size_max]  # example scaled sizes
    size_labels = [sizes.min(), round(sizes.max()/2), sizes.max()]
    size_handles = [plt.scatter([], [], s=s, color='gray', alpha=0.7) for s in size_values]
    ax.legend(size_handles, size_labels, title="Cluster size (# cells)",
            bbox_to_anchor=(1.05, 0.5), loc="upper left", fontsize=6,
            title_fontsize=8, labelspacing=2, borderpad=1)

    plt.tight_layout()
    pdf.savefig(fig, bbox_inches="tight")
    plt.close(fig)


    #-- Page 9: PCA median marker expression
    fig, ax = plt.subplots(figsize=(7, 5))
    method_levels = all_medians_df["method_level"].unique()
    colors = sns.color_palette("tab10", n_colors=len(method_levels))  # assign a color to each method_level

    # scatter points by method_level
    for method_level, color in zip(method_levels, colors):
        subset = all_medians_df[all_medians_df["method_level"] == method_level]
        subset_sizes = scaled_sizes[subset.index]  # actual sizes

        ax.scatter(
            subset["PC1_cluster"],
            subset["PC2_cluster"],
            s=subset_sizes,
            color=color,
            label=method_level,
            alpha=0.7
        )

        # annotate cluster numbers
        for _, row in subset.iterrows():
            cluster_num = row["cluster_id"].split("_")[-1]
            ax.text(
                row["PC1_cluster"],
                row["PC2_cluster"],
                cluster_num,
                fontsize=2,
                alpha=0.8,
                ha="center",
                va="center"
            )

    ax.set_xlabel("PC1")
    ax.set_ylabel("PC2")
    ax.set_title("PCA of cluster median marker expression (clustering markers)")

    ## Legend for method_level (color)
    # create dummy points of fixed size
    color_handles = [plt.scatter([], [], color=color, s=50) for color in colors]
    color_labels = method_levels
    color_legend = ax.legend(color_handles, color_labels, title="Method / Level",
                            bbox_to_anchor=(1.05, 1), loc="upper left", fontsize=6,
                            title_fontsize=8, labelspacing=1)
    ax.add_artist(color_legend)

    ## Legend for cluster size
    size_values = [size_min, round(size_max/2), size_max]  # example scaled sizes
    size_labels = [sizes.min(), round(sizes.max()/2), sizes.max()]
    size_handles = [plt.scatter([], [], s=s, color='gray', alpha=0.7) for s in size_values]
    ax.legend(size_handles, size_labels, title="Cluster size (# cells)",
            bbox_to_anchor=(1.05, 0.5), loc="upper left", fontsize=6,
            title_fontsize=8, labelspacing=2, borderpad=1)

    plt.tight_layout()
    pdf.savefig(fig, bbox_inches="tight")
    plt.close(fig)





end = end_time()
elapsed_time(start, end)