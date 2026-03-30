#===============================================================================================#
# SPiSCy (Snakemake PIpeline for Spectral CYtometry)
#
# Phase: differential analysis
# Step: perform differential analysis (abundance + state) for each clustering algo result. 
# Package: diffcyt
# 
# Author: Émilie Roy
# Date: Feb 2026
# Version: v.1.0
#
# Input : marker values for each sample (.csv) + cluster labels for each cell in the sample (.csv)
# Outputs : heatmaps (.png), stat results (.csv), diagnostic plots for limma and voom (.png)
#===============================================================================================#

# Redirect Rscript messages and errors to log files
if (exists("snakemake")) {
    log_stdout = file(snakemake@log$stdout, open = "wt")
    log_stderr = file(snakemake@log$stderr, open = "wt")
    sink(log_stdout, type = "output")
    sink(log_stderr, type = "message")
    on.exit({
        sink(type = "output")
        sink(type = "message")
        close(log_stdout)
        close(log_stderr)
    }, add=TRUE)
}

cat("Start:", date(), "\n")
start_time = Sys.time()



###############################################
#### LIBRARY AND SNAKEMAKE VARIABLES SETUP ####
###############################################

# Library
set.seed(123)
library(diffcyt)
library(SummarizedExperiment)

# Snakemake inputs
markers_dir = snakemake@input[["markers"]]
labels_dir = snakemake@input[["labels_dir"]]
spe_config_file = snakemake@input[["spe_config"]]
metadata_csv = snakemake@input[["metadata"]]
marker_info_csv = snakemake@input[["marker_info"]]
res_csv = snakemake@output[["res"]]
stats_dir=snakemake@output[["stats_dir"]]
heatmaps_dir=snakemake@output[["heatmaps_dir"]]
limma_voom_plots_dir=snakemake@output[["limma_voom_plots_dir"]]
clustering_level=snakemake@params[["clustering_key"]]

# Load config file (YAML)
spe_config = yaml::read_yaml(spe_config_file)

# Extract and adapt config file variables to work with R
cat_variables = spe_config$categorical_variables
design_cols = spe_config$design_matrix_columns
contrast_vars = spe_config$contrast_matrix_var
formula_vars = spe_config$formula_variables
methods_DA = spe_config$diff_abundance_methods
methods_DS = spe_config$diff_state_methods
heatmap_args = spe_config$heatmap_viz

cat("Clustering key:", clustering_level, "\n")



####################################
#### PREPARE INPUTS FOR DIFFCYT ####
####################################

## List of dataframes containing marker values
markers_csvs = list.files(markers_dir, full.names=TRUE)
markers_list = lapply(markers_csvs, read.csv, stringsAsFactors = FALSE)
names(markers_list) = gsub("\\.csv$", "", basename(markers_csvs))
cat("Found", length(markers_list), "sample files", "\n")

## List of dataframes container cluster labels
labels_csvs = list.files(labels_dir, full.names=TRUE)
cat("Found", length(labels_csvs), "cluster label files found", "\n")

labels_list = lapply(labels_csvs, read.csv, stringsAsFactors = FALSE)
names(labels_list) = gsub("\\.csv$", "", basename(labels_csvs))

## Experiment_info df : metadata about the experiment
exp_info_df = read.csv(metadata_csv, header=TRUE)
names(exp_info_df)[names(exp_info_df) == 'filename'] = 'sample_id'
if (!("sample_id" %in% colnames(exp_info_df))) {
    stop(paste("'filename' column missing in metadata table:", metadata_csv))
  }

# Keep only rows corresponding to the samples in markers_list
exp_info_df_min = exp_info_df[exp_info_df$sample_id %in% names(markers_list), , drop = FALSE]
exp_info_df_min = exp_info_df_min[match(names(markers_list), exp_info_df_min$sample_id), ]

metadata_cols = colnames(exp_info_df_min)

cat("[TABLE PREVIEW] Imported metadata table:", "\n")
print(head(exp_info_df_min))


# Set all categorical variables as factors
cat("Setting these variables to categorical:", cat_variables, "\n")
for (col in cat_variables) {
  if (col %in% colnames(exp_info_df_min)) {
    exp_info_df_min[[col]] = factor(exp_info_df_min[[col]])
  } else {
    stop(paste("Column not found in dataframe:", col))
  }
}
if (length(markers_csvs) != nrow(exp_info_df_min)) {
  stop("Not all input samples have been found in metadata. Is there a filename mismatch?")
}

## Marker_info df : information about markers used in experiment
marker_info_df = read.csv(marker_info_csv, header=TRUE)
if (!("marker_name" %in% colnames(marker_info_df))) {
    stop(paste("'marker_name' column missing in marker_info table:", marker_info_csv))
  }
if (!("marker_class" %in% colnames(marker_info_df))) {
    stop(paste("'marker_class' column missing in marker_info table:", marker_info_csv))
  }

cat("Imported marker_info table:", "\n")
print(marker_info_df)




#############################
#### CHECK ROW ALIGNMENT ####
#############################

# Make sure cells are in the same order for marker data / cluster labels
check_row_alignment = mapply(function(marker_df, label_df, nm) {
  if (!("row_id" %in% colnames(marker_df))) {
    stop(paste("row_id missing in markers for sample:", nm))
  }
  if (!("row_id" %in% colnames(label_df))) {
    stop(paste("row_id missing in labels for sample:", nm))
  }
  
  res = all.equal(marker_df$row_id, label_df$row_id)
  
  if (!isTRUE(res)) {
    stop(paste("Row order mismatch detected in sample:", nm))
  }
  TRUE
}, markers_list, labels_list, names(markers_list))

cat("Row order between labels and cells matches for all samples", "\n")

# Remove row_id (necessary for diffcyt) if row order matches
markers_list = lapply(markers_list, function(df) {
  df[, !(colnames(df) %in% "row_id"), drop = FALSE]
})

labels_list = lapply(labels_list, function(df) {
  df[, !(colnames(df) %in% "row_id"), drop = FALSE]
})




########################
#### READ IN INPUTS ####
########################

prepped_data = prepareData(d_input = markers_list,
                        experiment_info = exp_info_df_min,
                        marker_info = marker_info_df)

# Add in the cluster labels
cluster_vector = unlist(lapply(labels_list, function(df) df[[clustering_level]]))
clusters = factor(cluster_vector, levels = sort(unique(cluster_vector)))
names(clusters) = NULL
rowData(prepped_data)$cluster_id = clusters






############################
#### CALCULATE FEATURES ####
############################

cat("Calculating features...", "\n")

counts_df = calcCounts(prepped_data)
medians_df = calcMedians(prepped_data)
medians_cluster_marker_df = calcMediansByClusterMarker(prepped_data)





##############################
#### SET UP TEST MATRICES ####
##############################

## Design matrix (for methods voom, limma, and edgeR)

# Make sure design columns have been specified and that they exist in metadata
if (!is.null(design_cols) || !length(design_cols) == 0) {
  cat("Generating design matrix...", "\n")
  
  design_cols = unique(design_cols)
  invalid_design_cols = setdiff(design_cols, metadata_cols)
  
  if (length(invalid_design_cols)) {
    stop("Unknown design columns: ", paste(invalid_design_cols, collapse=", "), ". ", "Design columns should match column names of metadata.csv.")
  }
  
  design = createDesignMatrix(exp_info_df_min, cols_design = design_cols)
  cat("Columns of design matrix:", "\n")
  print(colnames(design))
}


## Model formula (for methods GLMM and LMM)

# Make sure design columns have been specified and that they exist in metadata
if (!is.null(formula_vars) || !length(formula_vars) == 0) {
  cat("Generating model formula...", "\n")
  
  vars_fixed = formula_vars$fixed
  vars_random = formula_vars$random
  
  invalid_vars_fixed = setdiff(vars_fixed, metadata_cols)
  invalid_vars_random = setdiff(vars_random, metadata_cols)
  
  if (length(invalid_vars_fixed)) {
    stop("Unknown fixed variables in model formula: ", paste(invalid_vars_fixed, collapse=", "), ". ", "Fixed variables should match column names of metadata.csv.")
  }
  if (length(invalid_vars_random)) {
    stop("Unknown random variable in model formula: ", paste(invalid_vars_random, collapse=", "), ". ", "Random variables should match column names of metadata.csv.")
  }
  
  vars_overlap = intersect(vars_fixed, vars_random)
  if (length(vars_overlap)) {
    stop("In model formula, columns cannot be both fixed and random effects: ", paste(vars_overlap, collapse=", "))
  }
  
  formula = createFormula(experiment_info = exp_info_df_min,
                          cols_fixed = vars_fixed,
                          cols_random = vars_random)

  cat("Created model formula: ", "\n")
  print(formula$formula)  
}


## Contrast vectors
create_contrast_vector = function(design, factor_name) {
  # Find the first column whose name contains the factor name
  matches = grep(paste0("^", factor_name), colnames(design))
  
  if (length(matches) == 0) {
    stop(paste0("No column found containing for contrast matrix '", factor_name, "'"))
  }
  
  # Create one contrast vector per matched column
  contrast_list = lapply(matches, function(idx) {
    vec = rep(0, ncol(design))
    vec[idx] = 1
    vec
  })
  # Name each vector according to the design matrix column
  names(contrast_list) = colnames(design)[matches]
  
  contrast_list
}


# make sure contrast variables have been specified and that they exist in metadata 
if (is.null(contrast_vars) || length(contrast_vars) == 0) {
  stop("No contrast variables defined. At least one variable must be specified to run a differential analysis.")
}

contrast_vars = unique(contrast_vars)
invalid_contrast_var = setdiff(contrast_vars, metadata_cols)

if (length(invalid_contrast_var)) {
  stop("Unknown contrast variables: ", paste(contrast_vars, collapse=", "), ". ", "Contrast variables should match column names of metadata.csv.")
}

# generate contrast vectors based on chosen variable
contrast_vectors = lapply(contrast_vars, function(var) {
  create_contrast_vector(design, var)
})
names(contrast_vectors) = contrast_vars

cat("Contrast vectors:", "\n")
print(contrast_vectors)

# generate contrast matrix for each contrast vector
contrast_mtxs = lapply(contrast_vectors, function(var_list) {
  out = lapply(var_list, createContrast)
  names(out) = names(var_list)
  out
})
names(contrast_mtxs) = names(contrast_vectors)

cat("Contrast matrices:", "\n")
print(str(contrast_mtxs))

if (any(sapply(contrast_mtxs, function(x) is.null(names(x))))) {
  stop("Contrast matrices missing names — structure invalid.")
}


###############################
#### DIFFERENTIAL ANALYSIS ####
###############################

## Make sure appropriate methods are specified and all requisites for each method are defined
allowed_methods_DA = c("voom", "edgeR", "GLMM")
allowed_methods_DS = c("limma", "LMM")

methods_DA = unique(methods_DA)
methods_DS = unique(methods_DS)

invalid_methods_DA = setdiff(methods_DA, allowed_methods_DA)
invalid_methods_DS = setdiff(methods_DS, allowed_methods_DS)

if (length(invalid_methods_DA)) {
  stop("Unknown differential abundance methods: ", paste(invalid_methods_DA, collapse = ", ", "."), "Options are: 'voom', 'edgeR', and 'GLMM'.")
}

if (length(invalid_methods_DS)) {
  stop("Unknown differential state methods: ", paste(invalid_methods_DS, collapse = ", ", "."), "Options are: 'limma' and 'LMM'.")
}

if ("GLMM" %in% methods_DA || "LMM" %in% methods_DS) {
  if (!(exists("formula"))) {
    stop("GLMM and/or LMM requested, but no model formula defined.")
  }
}

if (any(c("edgeR", "voom") %in% methods_DA) || "limma" %in% methods_DS) {
  if (!exists("design")) {
    stop("edgeR, voom, and/or limma requested, but no design matrix defined.")
  }
}


## Helper function to skip method if error and to handle pdf_prefix for voom+limma
run_safe = function(fun, method_name, pdf_prefix = NULL, ...) {
  tryCatch(
    {
      # run the function
      res = fun(...)
      
      # if pdf_prefix is provided and method is voom or limma, rename PDFs
      if (!is.null(pdf_prefix) && method_name %in% c("voom", "limma")) {
        # assume the function saves PDFs in a "path" argument
        args = list(...)
        path = if("path" %in% names(args)) args$path else "."
        
        # define old PDF filenames per method
        old_files = switch(method_name,
                            voom  = c("voom_before.pdf", "voom_after.pdf"),
                            limma = c("SA_plot.pdf"),
                            character(0))
        
        # construct new filenames with prefix
        new_files = paste0(pdf_prefix, "_", old_files)
        
        # rename only if the files exist
        for (i in seq_along(old_files)) {
          old_path = file.path(path, old_files[i])
          new_path = file.path(path, new_files[i])
          if (file.exists(old_path)) file.rename(old_path, new_path)
        }
      }
      
      res
    },
    error = function(e) {
      writeLines(
        paste0("SKIPPING method ", method_name, ": ", conditionMessage(e)),
        con = stderr()
      )
      NULL
    }
  )
}


## Differential cluster abundance (options are edgeR, voom, GLMM)

dir.create(limma_voom_plots_dir, recursive = TRUE, showWarnings = FALSE)

results_DA = lapply(names(contrast_mtxs), function(var_name) {
  
  lapply(names(contrast_mtxs[[var_name]]), function(ct_name) {
    
    contrast = contrast_mtxs[[var_name]][[ct_name]]
    prefix = ct_name

    # print(contrast)
    # print(prefix)
    # print(class(contrast))
    # print(class(counts_df))
    # print(class(design))
    # print(class(formula))
    # print(is(contrast))
    # print(is(counts_df))
    # print(is(design))
    # print(is(formula))
    
    method_results = lapply(methods_DA, function(method) {
      if (method == "voom") {
        run_safe(
          testDA_voom,
          method_name = "voom",
          d_counts = counts_df,
          design = design,
          contrast = contrast,
          plot = TRUE,
          path = limma_voom_plots_dir,
          pdf_prefix = prefix)

        } else if (method == "edgeR") {
          run_safe(
            testDA_edgeR,
            method_name = "edgeR",
            d_counts = counts_df,
            design = design,
            contrast = contrast)

        } else if (method == "GLMM") {
          run_safe(
            testDA_GLMM,
            method_name = "GLMM",
            d_counts = counts_df,
            formula = formula,
            contrast = contrast)

        }
      
    })
    
    names(method_results) = methods_DA
    method_results
    
  }) |> setNames(names(contrast_mtxs[[var_name]]))
  
}) |> setNames(names(contrast_mtxs))

cat("Saved results for differential abundance analysis:", "\n")
print(results_DA)


## Differential cellular states (options are limma or LMM)
results_DS = lapply(names(contrast_mtxs), function(var_name) {
  
  lapply(names(contrast_mtxs[[var_name]]), function(ct_name) {
    
    contrast = contrast_mtxs[[var_name]][[ct_name]]
    prefix = ct_name
    
    method_results = lapply(methods_DS, function(method) {
      if (method == "limma") {
        run_safe(
          testDS_limma,
          method_name="limma",
          pdf_prefix=prefix,
          d_counts=counts_df,
          d_medians=medians_df,
          design=design,
          contrast=contrast,
          plot=TRUE,
          path=limma_voom_plots_dir)

      } else if (method == "LMM") {
        run_safe(
          testDS_LMM,
          method_name="LMM",
          d_counts=counts_df,
          d_medians=medians_df,
          formula=formula,
          contrast=contrast) 
      }
      
    })
    
    names(method_results) = methods_DS
    method_results
    
  }) |> setNames(names(contrast_mtxs[[var_name]]))
  
}) |> setNames(names(contrast_mtxs))

cat("Saved results for differential state analysis:", "\n")
print(results_DS)





###################################
#### SAVE STATS RESULTS TO CSV ####
###################################

dir.create(stats_dir, recursive = TRUE, showWarnings = FALSE)

## Helper function to save results of a method to CSV
save_method_results = function(res, path, prefix, method_name) {
  if (is.null(res)) {
    cat("No results for method ", method_name, "; skipping CSV.", "\n")
    return(NULL)
  }
  
  # rowData contains the statistical results
  df = as.data.frame(rowData(res))
  
  # create output filename
  out_file = file.path(path, paste0(prefix, "_", method_name, ".csv"))
  
  # write CSV
  write.csv(df, out_file, row.names = FALSE)
  
  cat("Saved results for ", method_name, " to ", out_file, "\n")
  out_file
}

## Save each contrast matrice results
for (var_name in names(results_DA)) {
  for (ct_name in names(results_DA[[var_name]])) {
    
    prefix = ct_name
    
    for (method_name in names(results_DA[[var_name]][[ct_name]])) {
      res = results_DA[[var_name]][[ct_name]][[method_name]]
      
      save_method_results(res, path = stats_dir, prefix = prefix, method_name = method_name)
    }
  }
}


for (var_name in names(results_DS)) {
  for (ct_name in names(results_DS[[var_name]])) {
    
    prefix = ct_name
    
    for (method_name in names(results_DS[[var_name]][[ct_name]])) {
      res = results_DS[[var_name]][[ct_name]][[method_name]]
      
      save_method_results(res, path = stats_dir, prefix = prefix, method_name = method_name)
    }
  }
}





##################
#### HEATMAPS ####
##################

dir.create(heatmaps_dir, recursive = TRUE, showWarnings = FALSE)

## Arguments for png creation
nb_type_markers = nrow(marker_info_df[marker_info_df$marker_class=="type",])
nb_files = nrow(exp_info_df_min)
png_width = (nb_type_markers + nb_files) * 100

nb_clusters = length(levels(clusters))
png_height = nb_clusters * 200


## Check arguments for visualization
sign_threshold = heatmap_args$threshold
top_nb = heatmap_args$threshold

if (is.null(sign_threshold) || !is.numeric(sign_threshold) || sign_threshold<=0) {
  stop("`threshold` must be defined and bigger than 0.")
}

if (is.null(top_nb) || !is.numeric(top_nb) || top_nb<=0) {
  stop("`top_nb` must be defined and bigger than 0.")
}


## Differential abundance
for (var_name in names(results_DA)) {
  for (ct_name in names(results_DA[[var_name]])) {
    for (method_name in names(results_DA[[var_name]][[ct_name]])) {
      
      res = results_DA[[var_name]][[ct_name]][[method_name]]
      
      # skip NULL results
      if (is.null(res)) next
      
      if (inherits(res, "SummarizedExperiment")) {
        
        # construct filename
        file_name = paste0(ct_name, "_", method_name, ".png")
        file_path = file.path(heatmaps_dir, file_name)
        
        # open PNG device
        png(filename = file_path, width = png_width, height = png_height, res = 300)
        
        # generate heatmap
        plotHeatmap(analysis_type="DA",
                    threshold=sign_threshold,
                    res=res,
                    d_se=diffcyt_input,
                    d_counts=counts_df,
                    d_medians=medians_df,
                    d_medians_by_cluster_marker = medians_cluster_marker_df)
        dev.off()
      }
    }
  }
}


## Differential state
for (var_name in names(results_DS)) {
  for (ct_name in names(results_DS[[var_name]])) {
    for (method_name in names(results_DS[[var_name]][[ct_name]])) {
      
      res = results_DS[[var_name]][[ct_name]][[method_name]]
      
      # skip NULL results
      if (is.null(res)) next
      
      if (inherits(res, "SummarizedExperiment")) {
        
        # construct filename
        file_name = paste0(ct_name, "_", method_name, ".png")
        file_path = file.path(heatmaps_dir, file_name)
        
        # open PNG device
        png(filename = file_path, width = png_width, height = png_height, res = 300)
        
        # generate heatmap
        plotHeatmap(analysis_type="DS",
                    threshold=sign_threshold,
                    res=res,
                    d_se=diffcyt_input,
                    d_counts=counts_df,
                    d_medians=medians_df,
                    d_medians_by_cluster_marker = medians_cluster_marker_df)
        dev.off()
      }
    }
  }
}



end_time = Sys.time()
elapsed = end_time - start_time
cat("End:", date(), "\n")
cat("Total execution time:", elapsed)