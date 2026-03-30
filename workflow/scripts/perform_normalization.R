#===============================================================================================#
# SPiSCy (Snakemake PIpeline for Spectral CYtometry)
#
# Phase: preprocessing
# Step: correct batch effects with CytoNorm
# Script 7/8
# 
# Author: Émilie Roy
# Date: Sept 2025
# Version: v.1.0
#
# Input : QCed fcs files (.fcs) from QC.R
# Outputs : final fcs files (.fcs) + marker_ranges_before.csv
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



#### LIBRARY AND SNAKEMAKE VARIABLES SETUP ####

# Library
set.seed(123)
library(flowViz)         # flow cytometry plots
library(flowCore)        # manipulating fcs files
library(ggplot2)         # data visualization
library(dplyr)           # sampling
library(gridExtra)       # plotting
library(ggcyto)          # fortify method to make flowsets compatible with ggplot2
library(CytoNorm)        # normalization
library(FlowSOM)         # for AggregateFlowFrames

# Snakemake inputs
fcs_sample = snakemake@input[["fcs_sample"]]
fcs_control = snakemake@input[["fcs_control"]]
spe_config_file = snakemake@input[["specific_config"]]
gen_config_file = snakemake@input[["general_config"]]
metadata_path = snakemake@input[["metadata"]]
output_fcs = snakemake@output[["fcs_out"]]
norm_model_rds = snakemake@output[["norm_model"]]
marker_ranges_after = snakemake@output[["marker_ranges_aft"]]
marker_ranges_before = snakemake@output[["marker_ranges_bef"]]
batches_used = snakemake@output[["batches_used"]]


# Load config file (YAML)
spe_config = yaml::read_yaml(spe_config_file)
gen_config = yaml::read_yaml(gen_config_file)

# Extract and adapt config file variables to work with R
channel_alias_map = gen_config$channel_alias_map
channel_alias_df = data.frame(
  channels = names(channel_alias_map),
  alias = unlist(channel_alias_map),
  stringsAsFactors = FALSE,
  row.names = NULL
)

controls_exist = as.logical(spe_config$controls_exist)
control_ID = as.character(spe_config$control_ID)
markers_to_normalize = spe_config$markers_to_normalize
clustering_markers = spe_config$markers_for_clustering
nb_sampled_cells = as.numeric(spe_config$nb_sampled_cells)
flowsom_params = spe_config$flowsom_parameters
norm_params = spe_config$normalization_parameters




#### HELPER FUNCTIONS ####

get_marker_ranges = function(fs, output_csv = FALSE, csv_path = NULL) {
  # fs must be a flowset
  if (!inherits(fs, "flowSet")) {
  stop("Input must be a flowSet.")
  }
  ffs = flowSet_to_list(fs)


  # initialize results
  marker_ranges = list()

  # loop through flowFrames
  for (ff in ffs) {
    expr = exprs(ff)
    col_names = colnames(expr)

    # loop through each marker
    for (marker in col_names) {
      values = expr[, marker]

      if (!marker %in% names(marker_ranges)) {
        marker_ranges[[marker]] = list(
          min = min(values, na.rm = TRUE),
          max = max(values, na.rm = TRUE)
        )
      } else {
        marker_ranges[[marker]]$min = min(marker_ranges[[marker]]$min, min(values, na.rm = TRUE))
        marker_ranges[[marker]]$max = max(marker_ranges[[marker]]$max, max(values, na.rm = TRUE))
      }
    }
  }

  # convert to dataframe
  marker_summary = data.frame(
    marker = names(marker_ranges),
    min = sapply(marker_ranges, function(x) x$min),
    max = sapply(marker_ranges, function(x) x$max),
    row.names = NULL
  )

  # optionally save as CSV
  if (output_csv) {
    if (is.null(csv_path)) stop("Please provide a csv_file path to save output.")
    write.csv(marker_summary, file = csv_path, row.names = FALSE)
  }

  return(marker_summary)
}




#### IMPORT METADATA AND FCS FILES ####

#-- Prepare the controls (either existing or artificially created) to train the CytoNorm model
#-- Import fcs files into appropriate flowsets (set of flowframes aka fcs file)

#1) Import metadata
metadata = read.csv(metadata_path, header=TRUE)

#1) Create artificial controls if no controls exist
if (!controls_exist) {
  if (nb_sampled_cells == 0 || is.null(nb_sampled_cells)) {
    stop("Must sample at least 1 cell from every fcs file to create aggregated control files (nb_sampled_cells > 1)")
  }
    # setting sample_files
    sample_files = fcs_sample
    fcs_keys = tools::file_path_sans_ext(basename(sample_files))

    # Build lookup between basenames and full paths
    fcs_lookup = setNames(sample_files, fcs_keys)

    # Add matching paths to metadata (only if filename exists in lookup)
    metadata$file_path = fcs_lookup[metadata$filename]

    # Ensure all metadata files matched
    if (any(is.na(metadata$file_path))) {
      warning("Some files in metadata have no matching FCS file: ",
              paste(metadata$filename[is.na(metadata$file_path)], collapse = ", "))
      metadata = metadata[!is.na(metadata$file_path), ]
    }

    # Split files by batch
    batch_groups = split(metadata$file_path, metadata$batch)

    # Sample n cells from every file in a batch to create the artificial control for that batch
    control_files = list()
    for (batch_nb in names(batch_groups)) {
      cat("Creating artificial control for batch:", batch_nb, "\n")
      cat("Sampling", nb_sampled_cells, "from each file below to create an aggregate file", "\n")
      print(batch_groups[[batch_nb]])
      
      agg_batch_nb = AggregateFlowFrames(
        fileNames = batch_groups[[batch_nb]],
        cTotal = length(batch_groups[[batch_nb]]) * nb_sampled_cells
      )
      # change flowFrame identifier
      identifier(agg_batch_nb) = paste0("artificial_control_batch_", batch_nb)

      control_files = append(control_files, agg_batch_nb)
    }

} else if (controls_exist) {
  sample_files = fcs_sample
  control_files = fcs_control

} else {
  stop(paste(
    "Invalid option for 'controls_exist'. Must be either 'True' or 'False'"
  ))
}

cat("Since 'controls_exist' set to", controls_exist, ",", "using following files:", "\n")
cat(length(control_files), "control files:", "\n")
if (controls_exist){
  print(basename(control_files))
} else {
  for (i in seq_along(control_files)) {
    ff_name = identifier(control_files[[i]])
    print(ff_name)
  }
}
cat(length(sample_files), "sample files:", "\n")
print(basename(sample_files))


#2) Read fcs files into two different flowsets
# flowset of technical replicates to train normalization model
if (controls_exist) {
  fs_controls = read.flowSet(
    files = control_files,
    transformation = FALSE,
    truncate_max_range = FALSE,
    channel_alias = channel_alias_df
  )
} else {
  fs_controls = flowSet(control_files)
  sampleNames(fs_controls) = sapply(control_files, identifier) # assigning identifier to sampleName for each flowframe in flowset
}

# flowset of biological samples on which normalization model will be applied
fs_samples = read.flowSet(
  files = sample_files,
  transformation = FALSE,
  truncate_max_range = FALSE,
  channel_alias = channel_alias_df
)
cat("Controls flowset:", "\n")
cat("* if using artificial control samples, it's normal for the flowset to have extra columns compared to the samples flowset", "\n")
print(fs_controls)
cat("Samples flowset", "\n")
print(fs_samples)



#### NORMALIZATION: DETECTING AND CORRECTING BATCH EFFECT ####

#--If the experiment was done in multiple batches, then a potential batch effect might arise. 
#--Instead of clustering based on biological differences, clustering might be performed based on technical differences.
#--That's why its important to asses if a batch effect is present, and if so, correct it. 
#--To model the batch effect, CytoNorm suggests a technical replicate that was processed in every batch.
#--However, it is possible to use an artificial control sample constructed for every batch.
#--Once the batch effect has been modeled, it can be corrected using CytoNorm.

#1) Associating each file to its batch number
# for sample files
sample_basenames = tools::file_path_sans_ext(basename(sampleNames(fs_samples)))
sample_labels = metadata$batch[match(sample_basenames, metadata$filename)]


# for control files
if (controls_exist){
  control_basenames = tools::file_path_sans_ext(basename(sampleNames(fs_controls)))
  control_labels = metadata$batch[match(control_basenames, metadata$filename)]
} else {
  control_basenames = sampleNames(fs_controls)
  control_labels = sub("artificial_control_batch_", "", sampleNames(fs_controls))
}

sample_df = data.frame(sample_basenames, sample_labels)
colnames(sample_df) = c("Sample", "Batch number")
control_df = data.frame(control_basenames, control_labels)
colnames(control_df) = data.frame("Control", "Batch number")
cat("File and batch number associations:", "\n")
print(sample_df)
print(control_df)


#2) Getting marker ranges across all files to avoid extrapolation issues during model training (https://github.com/saeyslab/CytoNorm/issues/28)
marker_ranges = get_marker_ranges(fs=fs_controls, output_csv = TRUE, csv_path = marker_ranges_before)

cluster_ranges = marker_ranges[marker_ranges$marker %in% markers_to_normalize, ]
min_range = min(cluster_ranges$min)
max_range = max(cluster_ranges$max)


#3) Modeling batch effect using control files
cat("Training normalization model...", "\n")
cat("Markers that will be normalized:","\n")
print(markers_to_normalize)
cat("Markers used for FlowSOM clustering in the normalization model:", "\n")
print(clustering_markers)
norm_model = CytoNorm.train(files = fs_controls,
                             labels = control_labels,
                             channels = markers_to_normalize,
                             transformList = NULL,
                             FlowSOM.params = list(nCells = as.numeric(flowsom_params$nb_cells_in_model),
                                                   xdim = as.numeric(flowsom_params$xdim),
                                                   ydim = as.numeric(flowsom_params$ydim),
                                                   nClus = as.numeric(flowsom_params$nb_metaclusters),
                                                   scale = FALSE,
                                                   colsToUse = clustering_markers),
                             normMethod.train = QuantileNorm.train,
                             normParams = list(nQ = as.numeric(norm_params$nQ), limit=c(min_range, max_range)),
                             seed = 1,
                             verbose = FALSE,
                             clean = TRUE,
                             recompute = TRUE)


#5) Applying the normalization model to the sample files
cat("Applying normalization model to samples...", "\n")
CytoNorm.normalize(model = norm_model,
                   files = fs_samples,
                   labels = sample_labels,
                   transformList = NULL,
                   transformList.reverse = NULL,
                   normMethod.normalize = QuantileNorm.normalize,
                   outputDir = dirname(output_fcs[[1]]),
                   prefix = "",
                   clean = TRUE,
                   verbose = FALSE)


#6) Importing the normalized .fcs files into a new flowSet
fs_norm = read.flowSet(path = dirname(output_fcs[[1]]),
                        pattern = ".fcs",
                        transformation = FALSE,
                        truncate_max_range = FALSE,
                        channel_alias = channel_alias_df)


#7) Getting the new ranges after normalization to make sure no abnormally high negative or positive values
marker_ranges = get_marker_ranges(fs=fs_norm, output_csv = TRUE, csv_path = marker_ranges_after)

#8) Extract batche numbers used for next preprocessing step
batches_used = names(norm_model$clusterRes[[1]]$splines)
writeLines(batches_used, file.path(dirname(norm_model_rds), "batches_used.txt"))

#9) Save the normalization model for next script
cat("Saving normalization model:", basename(norm_model_rds), "\n")
saveRDS(object=norm_model, file=norm_model_rds)


end_time = Sys.time()
elapsed = end_time - start_time
cat("End:", date(), "\n")
cat("Total execution time:", elapsed)
