from scripts.resources import mem_for_diff_analysis, mins_diff_analysis

container: f"{workflow.basedir}/../apptainers/spiscy_rbase.sif"

rule diff_analysis:
    input:
        markers="results/diff_analysis/markers/",
        labels_dir="results/diff_analysis/labels/{algorithm}",
        spe_config="config/diff_analysis.yaml",
        metadata="data/metadata.csv",
        marker_info="data/marker_info.csv"
    output:
        stats_dir=directory("results/diff_analysis/{algorithm}/{level}/stats"),
        heatmaps_dir=directory("results/diff_analysis/{algorithm}/{level}/heatmaps"),
        limma_voom_plots_dir=directory("results/diff_analysis/{algorithm}/{level}/limma_voom_plots")
    log:
        stdout="results/logs/diff_analysis/diff_analysis_{algorithm}_{level}.stdout",
        stderr="results/logs/diff_analysis/diff_analysis_{algorithm}_{level}.stderr"
    params:
        clustering_key="{level}"
    threads:8
    resources:
        mem_mb=mem_for_diff_analysis,
        runtime=mins_diff_analysis
    script:
        "../scripts/differential_analysis.R"
