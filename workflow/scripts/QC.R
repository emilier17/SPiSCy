#=============================================================================================#
# SPiSCy (Snakemake PIpeline for Spectral CYtometry)
#
# Phase: preprocessing
# Step: acquisition quality control
# Script 5/8
# 
# Author: Émilie Roy
# Date: Sept 2025
# Version: v.1.0
#
# Input : final gated fcs files from secondary_gating.R (.fcs)
# Outputs : QCed fcs file (.fcs) + plots for validation (.png) + report (csv)
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
library(flowCore)        # manipulating fcs files
library(PeacoQC)         # quality control

# Snakemake inputs
fcs_file = snakemake@input[["fcs"]]
spe_config_file = snakemake@input[["specific_config"]]
gen_config_file = snakemake@input[["general_config"]]
output_fcs = snakemake@output[["fcs_out"]]
sample_report_output = snakemake@output[["report"]]
sample_plot_output = snakemake@output[["plot"]]

# Load config file (YAML)
spe_config = yaml::read_yaml(spe_config_file)
gen_config = yaml::read_yaml(gen_config_file)

# Extract and adapt config file variables to work with R
channel_alias_map = gen_config$channel_alias_map
channel_alias_df = data.frame(
  channels = names(channel_alias_map),
  alias   = unlist(channel_alias_map),
  stringsAsFactors = FALSE,
  row.names = NULL
)

# Extract user settings for PeacoQC
MAD_value = spe_config$MAD
IT_value = spe_config$IT
criteria = as.character(spe_config$criteria)
markers_to_evaluate = spe_config$markers_to_evaluate

if ( is(MAD_value, 'numeric') == FALSE || is(IT_value, 'numeric') == FALSE ) {
  stop("MAD and IT must be numeric values")
}

#### IMPORT FCS FILE FROM PREVIOUS PREPROCESSING STEP ####

#-- Import fcs file into a flowframe

ff = read.FCS(
  filename = fcs_file,
  transformation = FALSE,
  truncate_max_range = FALSE,
  channel_alias = channel_alias_df)



#### ACQUISITION QUALITY CONTROL ####

#--Remove low quality and abnormal events caused by acquisition issues (clogs, flow rate changes...)
#--PeacoQC detects outliers with IT and MAD

# setup path
QC_folder = dirname(sample_report_output) # so results/QC

#1) Perform PeacoQC
cat("Running PeacoQC with criteria:", criteria, "\n")

if (criteria == "all") {
  if (is.null(MAD_value) || is.null(IT_value)){
    stop("A non-null value must be provided for MAD and IT if criteria selected is 'all'")
  }
  cat("Using MAD =", MAD_value, ",", "IT =", IT_value, "\n")
  peacoqc_res <- PeacoQC(ff = ff,
                       channels = markers_to_evaluate,
                       determine_good_cells = criteria, # cell must pass criteria to be included in final fcs file
                       save_fcs = FALSE,
                       suffix_fcs = NULL,
                       plot = TRUE,
                       report = TRUE,
                       name_directory = QC_folder,
                       IT_limit = as.numeric(IT_value),
                       MAD = as.numeric(MAD_value)
                       )

} else if (criteria == "IT") {
  if (is.null(IT_value)){
    stop("A non-null value must be provided for IT if criteria selected is 'IT'")
  }
  cat("Using IT =", IT_value, "\n")
  peacoqc_res <- PeacoQC(ff = ff,
                       channels = markers_to_evaluate,
                       determine_good_cells = criteria, # cell must pass IT criteria to be included in final fcs file
                       save_fcs = FALSE,
                       suffix_fcs = NULL,
                       plot = TRUE,
                       report = TRUE,
                       name_directory = QC_folder,
                       IT_limit = as.numeric(IT_value)
                       )

} else if (criteria == "MAD") {
  if (is.null(MAD_value)){
    stop("A non-null value must be provided for MAD if criteria selected is 'MAD'")
  }
  cat("Using MAD =", MAD_value, "\n")
  peacoqc_res <- PeacoQC(ff = ff,
                       channels = markers_to_evaluate,
                       determine_good_cells = criteria, # cell must pass MAD criteria to be included in final fcs file
                       save_fcs = FALSE,
                       suffix_fcs = NULL,
                       plot = TRUE,
                       report = TRUE,
                       name_directory = QC_folder,
                       MAD = as.numeric(MAD_value)
                       )

} else {
  stop("Unrecognized argument for PeacoQC. Options are 'all', 'IT' or 'MAD'")
}


#2) Rename PeacoQC_report to expected output (PeacoQC doesn't allow custom naming of report file)
file.rename(file.path(QC_folder, "PeacoQC_report.txt"), sample_report_output)

#3) Save gated flowframe as .fcs for next preprocessing step
cat("Saving fcs:", "\n")
write.FCS(peacoqc_res$FinalFF, filename=output_fcs)



end_time = Sys.time()
elapsed = end_time - start_time
cat("End:", date(), "\n")
cat("Total execution time:", elapsed)

