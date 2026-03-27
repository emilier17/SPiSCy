from scripts.resources import (
    min_per_mb, mem_for_csv, mins_birch, mins_cytovi, mins_flowsom, mins_hdbscan,
    mins_parc, mins_phenograph, mins_evaluate_clustering, mins_compare_clustering,
    mins_split_samples, mins_split_clusters
    )

container: f"{workflow.basedir}/../apptainers/spiscy_pybase.sif"


rule run_flowsom:
    input:
        csv="results/csv/final_samples.csv",
        clustering_config="config/run_flowsom.yaml",
        sampling_config="config/clustering_sampling.yaml",
        label_prop_config="config/cluster_label_propagation.yaml",
        dr_config="config/dim_reduction.yaml"
    output:
        model="results/clustering/flowsom/flowsom_model.pkl",
        clusters="results/clustering/flowsom/clusters.csv",
        summary="results/clustering/flowsom/summary.pdf"
    log:
        stdout="results/logs/clustering/run_flowsom.stdout",
        stderr="results/logs/clustering/run_flowsom.stderr"
    threads:4
    resources:
        mem_mb=mem_for_csv,
        runtime=mins_flowsom
    script:
        "../scripts/run_flowsom.py"

rule run_phenograph:
    input:
        csv="results/csv/final_samples.csv",
        clustering_config="config/run_phenograph.yaml",
        sampling_config="config/clustering_sampling.yaml",
        label_prop_config="config/cluster_label_propagation.yaml",
        dr_config="config/dim_reduction.yaml"
    output:
        model="results/clustering/phenograph/phenograph_model.pkl",
        clusters="results/clustering/phenograph/clusters.csv"
    log:
        stdout="results/logs/clustering/run_phenograph.stdout",
        stderr="results/logs/clustering/run_phenograph.stderr"
    threads:4
    resources:
        mem_mb=mem_for_csv,
        runtime=mins_phenograph
    script:
        "../scripts/run_phenograph.py"

rule run_parc:
    input:
        csv="results/csv/final_samples.csv",
        clustering_config="config/run_parc.yaml",
        sampling_config="config/clustering_sampling.yaml",
        label_prop_config="config/cluster_label_propagation.yaml",
        dr_config="config/dim_reduction.yaml"
    output:
        model="results/clustering/parc/parc_model.pkl",
        clusters="results/clustering/parc/clusters.csv"
    log:
        stdout="results/logs/clustering/run_parc.stdout",
        stderr="results/logs/clustering/run_parc.stderr"
    threads:4
    resources:
        mem_mb=mem_for_csv,
        runtime=mins_parc
    script:
        "../scripts/run_parc.py"

rule run_cytovi:
    input:
        csv="results/csv/final_samples.csv",
        clustering_config="config/run_cytovi.yaml",
        sampling_config="config/clustering_sampling.yaml",
        label_prop_config="config/cluster_label_propagation.yaml",
        metadata="data/metadata.csv"
    output:
        model="results/clustering/cytovi/cytovi_model.pkl",
        clusters="results/clustering/cytovi/clusters.csv",
        elbo_plot="results/clustering/cytovi/elbo_plot.png"
    log:
        stdout="results/logs/clustering/run_cytovi.stdout",
        stderr="results/logs/clustering/run_cytovi.stderr"
    threads:8
    resources:
        mem_mb=mem_for_csv,
        runtime=mins_cytovi
    script:
        "../scripts/run_cytovi.py"

rule run_hdbscan:
    input:
        csv="results/csv/final_samples.csv",
        clustering_config="config/run_hdbscan.yaml",
        sampling_config="config/clustering_sampling.yaml",
        label_prop_config="config/cluster_label_propagation.yaml",
        dr_config="config/dim_reduction.yaml"
    output:
        model="results/clustering/hdbscan/hdbscan_model.pkl",
        clusters="results/clustering/hdbscan/clusters.csv",
    log:
        stdout="results/logs/clustering/run_hdbscan.stdout",
        stderr="results/logs/clustering/run_hdbscan.stderr"
    threads:8
    resources:
        mem_mb=mem_for_csv,
        runtime=mins_hdbscan
    script:
        "../scripts/run_hdbscan.py"

rule run_birch:
    input:
        csv="results/csv/final_samples.csv",
        clustering_config="config/run_birch.yaml",
        sampling_config="config/clustering_sampling.yaml",
        label_prop_config="config/cluster_label_propagation.yaml",
        dr_config="config/dim_reduction.yaml"
    output:
        model="results/clustering/birch/birch_model.pkl",
        clusters="results/clustering/birch/clusters.csv",
    log:
        stdout="results/logs/clustering/run_birch.stdout",
        stderr="results/logs/clustering/run_birch.stderr"
    threads:4
    resources:
        mem_mb=mem_for_csv,
        runtime=mins_birch
    script:
        "../scripts/run_birch.py"

rule evaluate_clustering:
    input:
        clusters="results/clustering/{algorithm}/clusters.csv",
        metadata="data/metadata.csv",
        markers="results/csv/final_samples.csv",
        specific_config="config/evaluate_clustering.yaml",
        model="results/clustering/{algorithm}/{algorithm}_model.pkl"
    output:
        batch_contribution_png="results/clustering/{algorithm}/QC/{level}/batch_contribution.png",
        batch_contribution_csv="results/clustering/{algorithm}/QC/{level}/batch_contribution.csv",
        file_contribution_png="results/clustering/{algorithm}/QC/{level}/file_contribution.png",
        file_contribution_csv="results/clustering/{algorithm}/QC/{level}/file_contribution.csv",
        cluster_size_png="results/clustering/{algorithm}/QC/{level}/{level}_size.png",
        cluster_size_csv="results/clustering/{algorithm}/QC/{level}/{level}_size.csv",
        marker_exp_median_png="results/clustering/{algorithm}/QC/{level}/marker_expression_median.png",
        marker_exp_median_csv="results/clustering/{algorithm}/QC/{level}/marker_expression_median.csv",
        marker_exp_sample_png="results/clustering/{algorithm}/QC/{level}/marker_expression_sample.png",
        marker_exp_sample_csv="results/clustering/{algorithm}/QC/{level}/marker_expression_sample.csv",
        umap_png="results/clustering/{algorithm}/QC/{level}/umap.png",
        umap_marker_png="results/clustering/{algorithm}/QC/{level}/umap_markers.png"
    params:
        clustering_key="{level}"
    log:
        stdout="results/logs/clustering/evaluate_{algorithm}_{level}.stdout",
        stderr="results/logs/clustering/evaluate_{algorithm}_{level}.stderr"
    threads:1
    resources:
        mem_mb=mem_for_csv,
        runtime=mins_evaluate_clustering
    script:
        "../scripts/evaluate_clustering.py"

if len(clustering_algos) > 1:
    rule compare_clustering:
        input:
            markers="results/csv/final_samples.csv",
            specific_config="config/compare_clustering.yaml",
            clusters_csvs=expand("results/clustering/{algorithm}/clusters.csv", algorithm=clustering_algos),
            median_marker_csvs=median_marker_csvs
        output:
            summary="results/clustering/clustering_comparaison.pdf"
        log:
            stdout="results/logs/clustering/compare_clustering.stdout",
            stderr="results/logs/clustering/compare_clustering.stderr"
        threads:8
        resources:
            mem_mb=mem_for_csv,
            runtime=mins_compare_clustering
        script:
            "../scripts/compare_clustering.py"


rule split_markers_by_filename:
    input:
        markers="results/csv/final_samples.csv"
    output:
        outdir=directory("results/diff_analysis/markers")
    log:
        stdout="results/logs/diff_analysis/input/split_markers.stdout",
        stderr="results/logs/diff_analysis/input/split_markers.stderr"
    threads:1
    resources:
        mem_mb=mem_for_csv,
        runtime=mins_split_samples
    script:
        "../scripts/split_markers_by_filename.py"

rule split_clusters_by_filename:
    input:
        markers="results/clustering/{algorithm}/clusters.csv"
    output:
        outdir=directory("results/diff_analysis/labels/{algorithm}")
    log:
        stdout="results/logs/diff_analysis/labels/split_clusters_{algorithm}.stdout",
        stderr="results/logs/diff_analysis/labels/split_clusters_{algorithm}.stderr"
    threads: 1
    resources:
        mem_mb=mem_for_csv,
        runtime=mins_split_clusters
    script:
        "../scripts/split_clusters_by_filename.py"