#===============================================================================================#
# SPiSCy (Snakemake PIpeline for Spectral CYtometry)
#
# Phase: preprocessing
# Step: visualize potential batch effects related to experimental conditions
# Script 6/8
# 
# Author: Émilie Roy
# Date: Sept 2025
# Version: v.1.0
#
# Input : QCed fcs files (.fcs) from QC.R
# Outputs : umap plots (.png)
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
library(CytoNorm)        # normalization
library(umap)            # umap

# Snakemake inputs
fcs_files = snakemake@input[["fcs"]]
metadata_path = snakemake@input[["metadata"]]
spe_config_file = snakemake@input[["specific_config"]]
gen_config_file = snakemake@input[["general_config"]]
metadata_csv = snakemake@input[["metadata"]]
umap_plots = snakemake@output[["umaps"]]

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

sample_size = as.numeric(spe_config$sample_size)
exp_conditions = as.character(spe_config$experimental_conditions)
umap_markers = as.character(spe_config$umap_markers)
data_to_use = as.character(spe_config$dataset_type)
control_ID = as.character(spe_config$control_ID)



#### HELPER FUNCTIONS ####

# Downsample flowframes in a flowset to contain a certain amount of cells
Downsampling_FlowSet <- function(x, samplesize , replace=TRUE, prob=NULL){
  if(missing(samplesize))
    samplesize <- min(flowCore::fsApply(x,nrow))
  flowCore::fsApply(x, function(ff){
    i <- sample(nrow(ff), size = samplesize, replace=replace, prob)
    ff[i,]
  })
}

# Color a UMAP plot according to certain condition
plot_umap = function(dataframe, condition) {
  plot = ggplot(dataframe, aes(x = UMAP1, y = UMAP2, color = factor(.data[[condition]]))) +
    geom_point(size = 0.01, alpha = 1) +
    labs(title = "Before normalization", color = condition, x = "UMAP1", y = "UMAP2") +
    theme_minimal() +
    theme(
      axis.text = element_text(size = 3),
      axis.title = element_text(size = 3),
      legend.title = element_text(size = 3),
      legend.text = element_text(size = 2),
      legend.spacing.y = unit(1, "cm"),
      legend.key.size = unit(0.15, "cm")
      ) 
  return(plot)
}



#### IMPORT FCS FILES AND METADATA ####

#1) Import metadata
metadata = read.csv(metadata_path, header=TRUE)

#2) Find appropriate fcs files according to user choice
if (length(control_ID) == 0) {
  control_ID <- NULL
}
has_control_id = !is.null(control_ID) && nzchar(control_ID)

if (data_to_use == "controls") {
  if (!has_control_id) {
    stop("To use 'controls', a 'control_ID' must be specified. Control_ID must be a term present in every filename of a control file")
  }
  files <- fcs_files[grepl(control_ID, tools::file_path_sans_ext(basename(fcs_files)), ignore.case = TRUE)]
  if (length(files) == 0) {
    stop("No control files found using specified control_ID")
  }
  
} else if (data_to_use == "samples") {
  if (has_control_id && any(grepl(control_ID, tools::file_path_sans_ext(basename(fcs_files)), ignore.case = TRUE))) {
    control_files <- fcs_files[grepl(control_ID, tools::file_path_sans_ext(basename(fcs_files)), ignore.case = TRUE)]
    files <- setdiff(fcs_files, control_files)
  } else {
    # No control_ID provided or no matching files — treat all as samples
    files <- fcs_files
  }
  
} else if (data_to_use == "both") {
  files <- fcs_files
  
} else {
  stop(paste(
    "Invalid option for 'data_to_use'. Must be one of 'samples', 'controls', or 'both'."
  ))
}

cat("Since files to use set to:", data_to_use, ", using following", length(files), "files:", "\n")
print(basename(files))


#3) Import fcs files into a flowset (set of flowframes aka fcs file)
fs_QCed = read.flowSet(
  files = files,
  transformation = FALSE,
  truncate_max_range = FALSE,
  channel_alias = channel_alias_df
)

#4) Downsample flowframes
cat("Downsampling flowframes to contain", sample_size, "cells", "\n")
cat("Total amount of cells across all files:", length(files)*sample_size, "\n")
fs_small <- Downsampling_FlowSet(x=fs_QCed, samplesize = sample_size)


#5) Convert flowset into a dataframe for umap
fs_small_df <- data.frame(matrix(ncol = length(colnames(fs_small)), nrow = 0)) #create dataframe with same column names as flowset
colnames(fs_small_df) <- colnames(fs_small)

fs_small_df <- fsApply(fs_small, function(ff) {
  # Extract expression matrix
  temp <- as.data.frame(exprs(ff))
  
  # Add filename as an identifier
  temp$filename <- tools::file_path_sans_ext(basename(keyword(ff)$FILENAME))
  
  return(temp)
}) %>% bind_rows()


#6) Add in metadata to dataframe
fs_small_df <- fs_small_df %>%
  left_join(metadata, by = "filename")


#7) Split dataframe into data and labels dataframes
fs_small_df_data = fs_small_df[, umap_markers] # for generating umap
if (!all(exp_conditions %in% colnames(metadata))) {
    stop("One or more experimental conditions are not columns in metadata.")
}
fs_small_df_labels = fs_small_df[, exp_conditions] # for visualization conditions' effects



#### UMAP : VISUALIZE EXPERIMENTAL CONDITIONS BATCH EFFECT ####

#-- Generate a umap based on the user's markers. The resulting umap will then be colored by the different experimental conditions.
#-- If the umap's coloration seems to correlate with an experimental condition (particularly batch number) and NOT a biological difference, then
#-- a batch effect might be present and needs to be corrected (next step)

#1) Create umap
cat("Generating UMAP using following markers:", umap_markers, "\n")
fs_small_umap = umap(fs_small_df_data)

#2) Extract umap coords and add metadata
umap_df <- as.data.frame(fs_small_umap$layout)
colnames(umap_df) <- c("UMAP1", "UMAP2")
umap_df <- cbind(umap_df, fs_small_df_labels)

#3) Save umap for each condition
for(i in seq_along(exp_conditions)){
  condition <- exp_conditions[i]
  output_file <- umap_plots[i]

  cat("Coloring UMAP according to", condition, "\n")
  
  plot = plot_umap(umap_df, condition)

  png(output_file, height=1700, width=2000, res=600)
  print(plot)
  dev.off()
}


end_time = Sys.time()
elapsed = end_time - start_time
cat("End:", date(), "\n")
cat("Total execution time:", elapsed)