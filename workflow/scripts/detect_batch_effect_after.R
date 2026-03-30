#===============================================================================================#
# SPiSCy (Snakemake PIpeline for Spectral CYtometry)
#
# Phase: preprocessing
# Step: visualize batch effect with umap, after normalization
# Script 9/9
# 
# Author: Émilie Roy
# Date: Sept 2025
# Version: v.1.0
#
# Input : normalized fcs files (.fcs) from perform_normalization.R
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
    labs(title = "After normalization", color = condition, x = "UMAP1", y = "UMAP2") +
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

#2) Import fcs files into a flowset (set of flowframes aka fcs file)
fs_normalized = read.flowSet(
  files = fcs_files,
  transformation = FALSE,
  truncate_max_range = FALSE,
  channel_alias = channel_alias_df
)
cat("Using following", length(fcs_files), "files:", "\n")
print(basename(fcs_files))


#3) Downsample flowframes
cat("Downsampling flowframes to contain", sample_size, "cells", "\n")
cat("Total amount of cells across all files:", length(fcs_files)*sample_size, "\n")
fs_small <- Downsampling_FlowSet(x=fs_normalized, samplesize = sample_size)


#4) Convert flowset into a dataframe for umap
fs_small_df <- data.frame(matrix(ncol = length(colnames(fs_small)), nrow = 0)) #create dataframe with same column names as flowset
colnames(fs_small_df) <- colnames(fs_small)

fs_small_df <- fsApply(fs_small, function(ff) {
  # Extract expression matrix
  temp <- as.data.frame(exprs(ff))
  
  # Add filename as an identifier
  temp$filename <- tools::file_path_sans_ext(basename(keyword(ff)$FILENAME))
  
  return(temp)
}) %>% bind_rows()


#5) Add in metadata to dataframe
fs_small_df <- fs_small_df %>%
  left_join(metadata, by = "filename")


#6) Split dataframe into data and labels dataframes
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