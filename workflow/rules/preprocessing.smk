from scripts.resources import mem_per_fcs, min_per_mb, mem_for_flowset, time_predict_cofactors, min_per_flowset, time_detect_before, time_detect_after, mem_for_evaluate_normalization

localrules: merge_counts, cleanup_cell_counts, merge_reports, cleanup_qc_reports, merge_csv

container: f"{workflow.basedir}/../apptainers/spiscy_rbase.sif"


rule prelim_gating:
    input:
        fcs="data/all_raw/{sample}.fcs",
        general_config="config/general_config.yaml",
        specific_config="config/prelim_gating.yaml"
    output:
        fcs_out="results/prelim_gating/all_gated/{sample}.fcs",
        gate_plots="results/prelim_gating/gates/{sample}_gates.png"
    log:
        stdout="results/logs/prelim_gating/{sample}.stdout",
        stderr="results/logs/prelim_gating/{sample}.stderr"
    threads: 1
    resources:
        mem_mb=mem_per_fcs,
        runtime=min_per_mb
    script:
        "../scripts/prelim_gating.R"


skip_FlowVS = transformation_config.get("skip_FlowVS", False)

if not skip_FlowVS:
    rule predict_cofactors:
        input:
            fcs=expand("results/prelim_gating/all_gated/{sample}.fcs", sample=ALL_FILES),
            general_config="config/general_config.yaml",
            specific_config="config/predict_cofactors.yaml"
        output:
            cofactors_csv="results/transformation/cofactors.csv",
            cofactors_graph="results/transformation/cofactor_bartlett.pdf",
            flowvs_cofactors="results/transformation/predicted_cofactors.csv"
        log:
            stdout="results/logs/predict_cofactors/predict_cofactors.stdout",
            stderr="results/logs/predict_cofactors/predict_cofactors.stderr"
        threads: 8
        resources:
            mem_mb=mem_for_flowset,
            runtime=time_predict_cofactors
        script:
            "../scripts/predict_cofactors.R"

rule transformation:
    input:
        fcs=expand("results/prelim_gating/all_gated/{sample}.fcs", sample=ALL_FILES),
        general_config="config/general_config.yaml",
        specific_config="config/transformation.yaml",
        cofactors_csv="results/transformation/cofactors.csv" if skip_FlowVS else rules.predict_cofactors.output.cofactors_csv
    output:
        fcs_out=expand("results/transformation/all_transformed/{sample}.fcs", sample=ALL_FILES),
        plots=expand("results/transformation/before_after_{marker}.png", marker=MARKERS)
    log:
        stdout="results/logs/transformation/transformation.stdout",
        stderr="results/logs/transformation/transformation.stderr"
    threads: 1
    resources:
        mem_mb=mem_for_flowset,
        runtime=min_per_flowset
    script:
        "../scripts/transformation.R"

rule secondary_gating:
    input:
        fcs="results/transformation/all_transformed/{sample}.fcs",
        general_config="config/general_config.yaml",
        specific_config="config/secondary_gating.yaml"
    output:
        fcs_out="results/secondary_gating/all_final_gated/{sample}.fcs",
        cell_counts=temp("results/secondary_gating/cell_counts/{sample}.csv"),
        gate_plot="results/secondary_gating/gates/{sample}_gates.png"
    log:
        stdout="results/logs/secondary_gating/{sample}.stdout",
        stderr="results/logs/secondary_gating/{sample}.stderr"
    threads: 1
    resources:
        mem_mb=mem_per_fcs,
        runtime=min_per_mb
    script:
        "../scripts/secondary_gating.R"

rule merge_counts:
    input:
        expand("results/secondary_gating/cell_counts/{sample}.csv", sample=ALL_FILES)
    output:
        "results/secondary_gating/cell_counts.csv"
    run:
        import pandas as pd
        dfs = [pd.read_csv(f) for f in input]
        merged = pd.concat(dfs, ignore_index=True)
        merged.to_csv(output[0], index=False)
    
rule cleanup_cell_counts:
    input:
        "results/secondary_gating/cell_counts.csv"
    output:
        touch("results/secondary_gating/.cleanup_cell_counts.done")
    run:
        import os, shutil
        shutil.rmtree("results/secondary_gating/cell_counts", ignore_errors=True)

rule quality_control:
    input:
        fcs="results/secondary_gating/all_final_gated/{sample}.fcs",
        general_config="config/general_config.yaml",
        specific_config="config/QC.yaml"
    output:
        fcs_out="results/QC/all_QCed/{sample}.fcs",
        plot="results/QC/PeacoQC_plots/PeacoQC_{sample}.png",
        report="results/QC/PeacoQC_report_{sample}.txt"
    log:
        stdout="results/logs/QC/{sample}.stdout",
        stderr="results/logs/QC/{sample}.stderr"
    shadow: "minimal"
    threads: 1
    resources:
        mem_mb=mem_per_fcs,
        runtime=min_per_mb
    script:
        "../scripts/QC.R"

rule merge_reports:
    input:
        expand("results/QC/PeacoQC_report_{sample}.txt", sample=ALL_FILES)
    output:
        "results/QC/PeacoQC_report.csv"
    run:
        import pandas as pd
        dfs = [pd.read_csv(f, sep="\t") for f in input]
        merged = pd.concat(dfs, ignore_index=True)
        merged.to_csv(output[0], index=False) 

rule cleanup_qc_reports:
    input:
        "results/QC/PeacoQC_report.csv"
    output:
        touch("results/QC/.cleanup_qc_reports.done")
    run:
        import glob, os
        for report in glob.glob("results/QC/PeacoQC_report_*.txt"):
            os.remove(report)

rule detect_batch_effect_before:
    input:
        fcs=expand("results/QC/all_QCed/{sample}.fcs", sample=ALL_FILES),
        metadata="data/metadata.csv",
        general_config="config/general_config.yaml",
        specific_config="config/detect_batch_effect_before.yaml"
    output:
        umaps=expand("results/detect_batch_effect/before_umap_{condition}.png", condition=exp_conditions_before)
    log:
        stdout="results/logs/detect_batch_effect/detect_batch_effect_before.stdout",
        stderr="results/logs/detect_batch_effect/detect_batch_effect_before.stderr"
    threads: 1
    resources:
        mem_mb=mem_for_flowset,
        runtime=time_detect_before
    script:
        "../scripts/detect_batch_effect_before.R"


skip_normalization = normalization_config.get("skip_normalization", False)


if skip_normalization:
    rule create_csv:
            input:
                fcs="results/QC/all_QCed/{sample}.fcs",
                general_config="config/general_config.yaml",
                specific_config="config/create_csv.yaml"
            output:
                csv="results/csv/samples_final/{sample}.csv"
            log:
                stdout="results/logs/create_csv/create_csv_{sample}.stdout",
                stderr="results/logs/create_csv/create_csv_{sample}.stderr"
            threads: 1
            resources:
                mem_mb=mem_per_fcs,
                runtime=min_per_mb
            script:
                "../scripts/create_csv.R"
    
    rule merge_csv:
        input:
            expand("results/csv/samples_final/{sample}.csv", sample=SAMPLE_FILES)
        output:
            "results/csv/final_samples.csv"
        run:
            import os
            import pandas as pd
            dfs = []
            for f in input:
                df = pd.read_csv(f)

                filename = os.path.splitext(os.path.basename(f))[0]
                df.insert(0, "filename", filename)

                # Enforce column order: filename, row_id, rest
                cols = df.columns.tolist()
                cols = ["filename", "row_id"] + [c for c in cols if c not in ("filename", "row_id")]
                df = df[cols]

                dfs.append(df)

            merged = pd.concat(dfs, ignore_index=True)
            assert not merged.duplicated(["filename", "row_id"]).any()
            merged.to_csv(output[0], index=False)


if not skip_normalization:
    rule perform_normalization:
        input:
            fcs_sample=expand("results/QC/all_QCed/{sample}.fcs", sample=SAMPLE_FILES),
            fcs_control=expand("results/QC/all_QCed/{sample}.fcs", sample=CONTROL_FILES),
            general_config="config/general_config.yaml",
            specific_config="config/normalization.yaml",
            metadata="data/metadata.csv"
        output:
            fcs_out=expand("results/normalization/samples_normalized/{sample}.fcs", sample=SAMPLE_FILES),
            norm_model="results/normalization/norm_model.rds",
            marker_ranges_bef="results/normalization/marker_ranges_before.csv",
            marker_ranges_aft="results/normalization/marker_ranges_after.csv",
            batches_used="results/normalization/batches_used.txt"
        log:
            stdout="results/logs/normalization/perform_normalization.stdout",
            stderr="results/logs/normalization/perform_normalization.stderr"
        threads: 1
        resources:
            mem_mb=mem_for_flowset,
            runtime=min_per_flowset
        script:
            "../scripts/perform_normalization.R"
    

    # Load batches used
    if os.path.exists("results/normalization/batches_used.txt"):
        with open("results/normalization/batches_used.txt") as f:
            BATCHES_USED = [line.strip() for line in f if line.strip()]
    else:
        print("Warning could not find results/normalization/batches_used.txt yet. Make sure perform_normalization has run.")
        BATCHES_USED = []  # placeholder before first run


    rule evaluate_normalization:
        input:
            fcs_before_norm=expand("results/QC/all_QCed/{sample}.fcs", sample=SAMPLE_FILES),
            fcs_after_norm=expand("results/normalization/samples_normalized/{sample}.fcs", sample=SAMPLE_FILES),
            batches_used="results/normalization/batches_used.txt",
            norm_model="results/normalization/norm_model.rds",
            specific_config="config/normalization.yaml",
            metadata="data/metadata.csv"
        output:
            test_CV="results/normalization/testCV.png",
            splines_plots=expand("results/normalization/splines/splines_batch{batch_nb}.png", batch_nb=BATCHES_USED),
            before_after_densities=expand("results/normalization/before_after/bef_aft_densities_{marker}.png", marker=norm_markers),
            before_after_densities_legend="results/normalization/before_after/densities_legend.png"
        log:
            stdout="results/logs/normalization/evaluate_normalization.stdout",
            stderr="results/logs/normalization/evaluate_normalization.stderr"
        threads: 1
        resources:
            mem_mb=mem_for_evaluate_normalization,
            runtime=min_per_flowset
        script:
            "../scripts/evaluate_normalization.R"


    rule detect_batch_effect_after:
        input:
            fcs=expand("results/normalization/samples_normalized/{sample}.fcs", sample=SAMPLE_FILES),
            metadata="data/metadata.csv",
            general_config="config/general_config.yaml",
            specific_config="config/detect_batch_effect_after.yaml"
        output:
            umaps=expand("results/detect_batch_effect/after_umap_{condition}.png", condition=exp_conditions_after)
        log:
            stdout="results/logs/detect_batch_effect/detect_batch_effect_after.stdout",
            stderr="results/logs/detect_batch_effect/detect_batch_effect_after.stderr"
        threads: 1
        resources:
            mem_mb=mem_for_flowset,
            runtime=time_detect_after
        script:
            "../scripts/detect_batch_effect_after.R"

    rule create_csv:
        input:
            fcs="results/normalization/samples_normalized/{sample}.fcs",
            general_config="config/general_config.yaml",
            specific_config="config/create_csv.yaml"
        output:
            csv="results/csv/samples_final/{sample}.csv"
        log:
            stdout="results/logs/create_csv/create_csv_{sample}.stdout",
            stderr="results/logs/create_csv/create_csv_{sample}.stderr"
        threads: 1
        resources:
            mem_mb=mem_per_fcs,
            runtime=min_per_mb
        script:
            "../scripts/create_csv.R"
    
    rule merge_csv:
        input:
            expand("results/csv/samples_final/{sample}.csv", sample=SAMPLE_FILES)
        output:
            "results/csv/final_samples.csv"
        run:
            import os
            import pandas as pd
            dfs = []
            for f in input:
                df = pd.read_csv(f)

                filename = os.path.splitext(os.path.basename(f))[0]
                df.insert(0, "filename", filename)

                # Enforce column order: filename, row_id, rest
                cols = df.columns.tolist()
                cols = ["filename", "row_id"] + [c for c in cols if c not in ("filename", "row_id")]
                df = df[cols]

                dfs.append(df)

            merged = pd.concat(dfs, ignore_index=True)
            assert not merged.duplicated(["filename", "row_id"]).any()
            merged.to_csv(output[0], index=False)
