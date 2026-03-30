#=============================================================================================#
# SPiSCy (Snakemake PIpeline for Spectral CYtometry)
#
# Phase: preprocessing
# Step: apply arcsinh transformation cofactors to fluorescent signal for each marker
# Script 3/8
# 
# Author: Émilie Roy
# Date: Sept 2025
# Version: v.1.0
#
# Input : gated fcs file from prelim_gating.R (.fcs) + cofactors.csv
# Outputs : transformed fcs file (.fcs) + plots for validation (.png)
#=============================================================================================#

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
library(flowVS)          # predict arcsinh cofactors
library(gridExtra)       # plotting
library(ggcyto)          # fortify method to make flowsets compatible with ggplot2

# Snakemake inputs
fcs_folder = snakemake@input[["fcs"]]
spe_config_file = snakemake@input[["specific_config"]]
gen_config_file = snakemake@input[["general_config"]]
cofactors_file = snakemake@input[["cofactors_csv"]]
output_fcs = snakemake@output[["fcs_out"]]
plot_files = snakemake@output[["plots"]]

# Load config file (YAML)
spe_config = yaml::read_yaml(spe_config_file)
gen_config = yaml::read_yaml(gen_config_file)

# Extract and adapt config file variables to work with R
markers_to_transform = as.vector(spe_config$markers_to_transform)
names(plot_files) = markers_to_transform

channel_alias_map = gen_config$channel_alias_map
channel_alias_df = data.frame(
  channels = names(channel_alias_map),
  alias = unlist(channel_alias_map),
  stringsAsFactors = FALSE,
  row.names = NULL
)

manual_cofactors = spe_config$arcsinh_cofactors
manual_cofactors = data.frame(
  marker = names(manual_cofactors),
  cofactor = unlist(manual_cofactors),
  stringsAsFactors = FALSE
)

sample_size = as.numeric(spe_config$downsample_size)



#### IMPORT FCS FILES FROM PREVIOUS PREPROCESSING STEP ####

#-- Import fcs files into a flowset (set of flowframes aka fcs file)

fs_gated = read.flowSet(
  files = fcs_folder,
  transformation = FALSE,
  truncate_max_range = FALSE,
  channel_alias = channel_alias_df
)



#### HELPER FUNCTIONS ####

# Downsample amount of cells for each flowframe in a flowset
downsampling_flowset = function(x, samplesize, replace = TRUE, prob = NULL) {
  if (missing(samplesize)) {
    samplesize = min(flowCore::fsApply(x, nrow))
  }
  # extract list of flowFrames
  ff_list = flowCore::fsApply(x, identity, simplify = FALSE)
  # downsample each flowFrame
  ff_list = lapply(ff_list, function(ff) {
    i = sample(nrow(ff), size = min(samplesize, nrow(ff)), replace = replace, prob = prob)
    ff[i, ]
  })
  # give names back to the list
  names(ff_list) = sampleNames(x)
  # rebuild flowSet
  flowCore::flowSet(ff_list)
}




#### ARCSINH TRANSFORMATION OF FLUORESCENCE SIGNAL ####

#-- Fluorescent signal from flow cytometry are usually transformed by arcsinh function.
#-- Transforming signal is necessary to compare values, for clustering and for visualization.
#-- FlowVS calculates the optimal arcsinh cofactor by a Bartlett test

#1) Retrieve predicted cofactors from a previous run, and overwrite them if user provided manual cofactors
if (!file.exists(cofactors_file)) {
    stop("Cofactors file not found. Please run predict_cofactors.R first.")
}
cofactors_csv = read.csv(cofactors_file, header=TRUE, sep=",")

# Merge manual cofactors if provided
if (nrow(manual_cofactors) > 0) {
    cofactors_csv = merge(
        cofactors_csv,
        manual_cofactors,
        by.x = "markers_to_transform",
        by.y = "marker",
        all.x = TRUE
    )
    # overwrite only where manual values exist
    cofactors_csv$cofactors = ifelse(
        is.na(cofactors_csv$cofactor),
        cofactors_csv$cofactors,
        cofactors_csv$cofactor
    )
    cofactors_csv$cofactor = NULL
    write.csv(cofactors_csv, file = cofactors_file, row.names=FALSE)
}

cofactors = cofactors_csv$cofactors


#2) Apply cofactors to all flowFrames in the flowset. Create a new flowset with the transformed data
cofactors_csv = read.csv(cofactors_file, header=TRUE, sep=",")
cofactor_map = setNames(cofactors_csv$cofactors, cofactors_csv$markers_to_transform)
cofactors_ordered = cofactor_map[markers_to_transform]

cat("Transforming following markers with following cofactors:", "\n")
print(cofactors_csv)
fs_transform = transFlowVS(fs_gated, channels = markers_to_transform, cofactors_ordered)

# reset filenames of fs_transforms since FlowVS changes them
filenames = sampleNames(fs_gated)
sampleNames(fs_transform) = filenames


#3) Visualize transformation effect on fluorescence signals before/after (FlowViz)
# transforming flowFrame into dataframe for ggplot2 and adding sampleNames for graph grouping
cat("Downsampling flowframes to have", sample_size, "cells", "\n")
fs_gated_small = downsampling_flowset(x = fs_gated, samplesize = sample_size)
fs_transform_small = downsampling_flowset(x = fs_transform, samplesize = sample_size)
df_before = fortify(fs_gated_small)
df_after = fortify(fs_transform_small)

# plotting functions before/after transformation
plot_before = function(marker) {
  plot = ggplot(df_before, aes(x = .data[[marker]], color = factor(name))) +
    geom_density() +
    theme_minimal() +
    theme(legend.position = "none") +
    ggtitle("Before transformation")
  return(plot)
}
plot_after = function(marker) {
  plot = ggplot(df_after, aes(x = .data[[marker]], color = factor(name))) +
    geom_density() +
    theme_minimal() +
    theme(legend.position = "none") +
    ggtitle("After transformation")
  return(plot)
}

# saving plots for all markers. Each curve is a fcs file (flowframe) from the flowset.
for (marker in markers_to_transform) {
  before = plot_before(marker)
  after = plot_after(marker)

  png(plot_files[[marker]], width=1800,height=900,res=300)
  print(grid.arrange(before, after, ncol=2))
  dev.off()
}

#4) Save each flowframe in the flowset as a fcs file for the next preprocessing step
output_fcs_named = setNames(output_fcs, gsub("\\.fcs$", "", basename(output_fcs)))
samples = gsub("\\.fcs$", "", sampleNames(fs_transform))

# reorder according to flowSet sampleNames
ordered_files = output_fcs_named[samples]

# safety check
if (any(is.na(ordered_files))) {
  stop("Some sampleNames(fs_transform) did not match expected Snakemake outputs.")
}
# write fcs file for every flowframe in flowset
write.flowSet(fs_transform, dirname(output_fcs[1]), basename(ordered_files))


end_time = Sys.time()
elapsed = end_time - start_time
cat("End:", date(), "\n")
cat("Total execution time:", elapsed)