#=============================================================================================#
# SPiSCy (Snakemake PIpeline for Spectral CYtometry)
#
# Phase: preprocessing
# Step: convert fcs file to csv file for downstream clustering
# Script 9/9
# 
# Author: Émilie Roy
# Date: Sept 2025
# Version: v.1.0
#
# Input : final fcs file (.fcs)
# Outputs : csv file (.csv) with columns = markers and each row = cell
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
library(flowCore)        # manipulating fcs files
library(dplyr)           # selecting columns

# Snakemake inputs and ouputs
fcs_file = snakemake@input[["fcs"]]
gen_config_file = snakemake@input[["general_config"]]
spe_config_file = snakemake@input[["specific_config"]]
csv_file = snakemake@output[["csv"]]

# Load config files
spe_config = yaml::read_yaml(spe_config_file)
gen_config = yaml::read_yaml(gen_config_file)

# Adapt config file types to work with flowCore
channel_alias_map = gen_config$channel_alias_map
channel_alias_df = data.frame(
  channels = names(channel_alias_map),
  alias   = unlist(channel_alias_map),
  stringsAsFactors = FALSE,
  row.names = NULL
)

markers_to_keep = spe_config$markers_to_export


#### READ FCS FILE AND CONVERT TO CSV ####

ff = read.FCS(
  filename = fcs_file,
  transformation = FALSE,
  truncate_max_range = FALSE,
  channel_alias = channel_alias_df)

# values for all markers in the flowframe
marker_df = as.data.frame(exprs(ff))

# save only the marker data for user chosen markers (config)
if (length(markers_to_keep)!=0){
    cat("Saving values for following markers only:", markers_to_keep, "\n")
    marker_df = marker_df %>% select(all_of(markers_to_keep))
    cat("Dimensions of final csv:", dim(marker_df), "\n")
} else {
    cat("Keeping all columns in final csv. Dimensions of final csv:", dim(marker_df), "\n")
}

# creating a unique id for every cell in the csv : <filename>__<row_nb_in_df>
# this will allow to match a cell's markers' expression values to the cell's clustering result
marker_df$row_id = seq_len(nrow(marker_df))
marker_df = marker_df[, c("row_id", setdiff(colnames(marker_df), "row_id"))]

# export as csv
cat("Exporting data as csv for file:", basename(fcs_file), "\n")
write.csv(marker_df, csv_file, row.names=FALSE)


end_time = Sys.time()
elapsed = end_time - start_time
cat("End:", date(), "\n")
cat("Total execution time:", elapsed)