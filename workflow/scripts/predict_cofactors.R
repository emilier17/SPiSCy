#===============================================================================================#
# SPiSCy (Snakemake PIpeline for Spectral CYtometry)
#
# Phase: preprocessing
# Step: predict arcsinh cofactors with FlowVS for transformation of fluorescent signal per marker
# Script 2/8
# 
# Author: Émilie Roy
# Date: Sept 2025
# Version: v.1.0
#
# Input : gated fcs file from prelim_gating.R (.fcs)
# Outputs : cofactors.csv + cofactor_graph.pdf
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
library(flowCore)        # manipulating fcs files
library(ggplot2)         # data visualization
library(dplyr)           # sampling
library(flowVS)          # predict arcsinh cofactors
library(gridExtra)       # plotting

# Snakemake inputs
fcs_folder = snakemake@input[["fcs"]]
spe_config_file = snakemake@input[["specific_config"]]
gen_config_file = snakemake@input[["general_config"]]
cofactor_file = snakemake@output[["cofactors_csv"]]
cofactors_graph = snakemake@output[["cofactors_graph"]]


# Load config file (YAML)
spe_config = yaml::read_yaml(spe_config_file)
gen_config = yaml::read_yaml(gen_config_file)

# Extract and adapt config file variables to work with R
markers_to_transform = as.vector(spe_config$markers_to_transform)

channel_alias_map = gen_config$channel_alias_map
channel_alias_df = data.frame(
  channels = names(channel_alias_map),
  alias   = unlist(channel_alias_map),
  stringsAsFactors = FALSE,
  row.names = NULL
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



#### FLOWVS: PREDICT OPTIMAL ARCSINH TRANSFORMATION COFACTORS ####

#-- Fluorescent signal from flow cytometry are usually transformed by arcsinh function.
#-- Transforming signal is necessary to compare values, for clustering and for visualization.
#-- FlowVS calculates the optimal arcsinh cofactor according to a Bartlett test (minimize channel-specific variance)
#-- In practice, the optimal cofactor found by FlowVS isn't always the most appropriate cofactor that describes the biological distribution of the marker

#1) Select markers to undergo transformation
cat("Predicing optimal cofactor for following channels: ", markers_to_transform, "\n")

#2) Create function: downsample amount of cells per flowFrame to speed up FlowVS calculations
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

#3) Run FlowVS to predict optimal cofactors for arcsinh transformation
cat("Downsampling each flowframe to have", sample_size, "cells", "\n")
gated_fs_small = downsampling_flowset(x = fs_gated, samplesize = sample_size)

# Predict cofactors and save Bartlett test results per marker
pdf(cofactors_graph)
cofactors = estParamFlowVS(gated_fs_small, channels = markers_to_transform)
dev.off()
cat("Saved Bartlett test results: ", cofactors_graph, "\n")

# Save predicted cofactors in a csv
cofactordata = data.frame(markers_to_transform, cofactors)
write.csv(x = cofactordata, file = cofactor_file, row.names = FALSE)
cat("Saved predicted cofactors: ", cofactor_file, "\n")


end_time = Sys.time()
elapsed = end_time - start_time
cat("End:", date(), "\n")
cat("Total execution time:", elapsed)
