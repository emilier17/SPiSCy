#===============================================================================================#
# SPiSCy (Snakemake PIpeline for Spectral CYtometry)
#
# Phase: preprocessing
# Step: evaluate the effects of CytoNorm normalization on marker values and distribution
# Script 8/8
# 
# Author: Émilie Roy
# Date: Sept 2025
# Version: v.1.0
#
# Input : normalized fcs files (.fcs) from perform_normalization.R
# Outputs : marker_ranges_after.csv + splines_plots.png + density_plots.png
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
library(ggpubr)          # for ggarrange
library(dplyr)           # sampling
library(gridExtra)       # plotting
library(ggcyto)          # fortify method to make flowsets compatible with ggplot2
library(CytoNorm)        # normalization


# Snakemake inputs
fcs_before_norm = snakemake@input[["fcs_before_norm"]]
fcs_after_norm = snakemake@input[["fcs_after_norm"]]
gen_config_file = snakemake@input[["general_config"]]
spe_config_file = snakemake@input[["specific_config"]]
norm_model_rds = snakemake@input[["norm_model"]]
batches_used = snakemake@input[["batches_used"]]
metadata_path = snakemake@input[["metadata"]]
bef_aft_plots = snakemake@output[["before_after_densities"]]
bef_aft_legend = snakemake@output[["before_after_densities_legend"]]
splines_plots = snakemake@output[["splines_plots"]]
test_CV = snakemake@output[["test_CV"]]

# Load config file (YAML)
spe_config = yaml::read_yaml(spe_config_file)

# Extract and adapt config file variables to work with R
markers_to_normalize = spe_config$markers_to_normalize
controls_exist = as.logical(spe_config$controls_exist)



#### IMPORT METADATA, NORM MODEL AND FCS FILES ####

#1) Import normalization model saved as .rds
norm_model = readRDS(norm_model_rds)
cat("Imported normalization model:", norm_model_rds, "\n")

#2) Import metadata
metadata = read.csv(metadata_path, header=TRUE)


#3) Waterfall plot (testCV) : check assumptions about batch effect on clusters and metaclusters
# CytoNorm assumes that the batch effect will not impact clustering or metaclustering
# TestCV measures the % of cells from each file assigned to each metaclusters,then measures the difference in % between all files (via coefficient of variation - CV)
# The waterfall plot helps decide an optimal amount of metaclusters that minimizes
cat("TestCV running and plotting...", "\n") 
cv_results = testCV(norm_model$fsom, cluster_values = c(3:20), plot = FALSE)
png(test_CV, width = 1800, height = 2500, res = 300, pointsize = 8)
PlotOverviewCV(fsom = norm_model$fsom, cv_res = cv_results, show_cv = 0.8, max_cv = 1.5)
dev.off()


#4) View the splines used by the normalization model for every batch and marker
cat("Saving splines plot for each batch...", "\n")

spline_plots = plotSplines(model=norm_model, channels=markers_to_normalize, groupClusters=TRUE)

# extract batch number from snakemake file path
batch_numbers = gsub(".*splines_batch(.*)\\.png$", "\\1", basename(splines_plots))

for (i in seq_along(batch_numbers)) {
  batch = batch_numbers[i]
  output_file = splines_plots[i]

  cat("batch:", batch, "\n")
  cat("output_file:", output_file, "\n")

  if (!batch %in% names(spline_plots)) {
    warning(paste("Skipping batch", batch, ", not found in model."))
    next
  }

  png(output_file, width = 2000, height = 1200, res = 200)
  print(spline_plots[[batch]])
  dev.off()
}


#5) Visualize marker densities before and after normalization for each marker

# first must separate the files into before/after and associate each file to its batch number
before_basenames = tools::file_path_sans_ext(basename(fcs_before_norm))
before_batches = metadata$batch[match(before_basenames, metadata$filename)]
after_basenames  = tools::file_path_sans_ext(basename(fcs_after_norm))
after_batches = metadata$batch[match(after_basenames, metadata$filename)]

if (!controls_exist) {
    control_indices = grep("artificial_control_batch_", before_basenames)
    before_batches[control_indices] = sub("artificial_control_batch_", "", before_basenames[control_indices])
    after_batches[control_indices] = sub("artificial_control_batch_", "", after_basenames[control_indices])
}

# build the input list necessary for plotDensities
batches = unique(before_batches)
input_list = list()
for (b in batches) {
    before_files = fcs_before_norm[before_batches == b]
    after_files  = fcs_after_norm[after_batches == b]
    
    if (length(before_files) > 0) {
        input_list[[paste0("BATCH", b)]] = before_files
    }
    if (length(after_files) > 0) {
        input_list[[paste0("BATCH", b, "_norm")]] = after_files
    }
}

cat("Input list for plotDensities. Before / after normalization fcs files associated to their batch numbers", "\n")
print(input_list)

# creating density plots for each marker that was normalized
cat("Plotting marker distributions per batch, before / after normalization...")
densities_plot = plotDensities(input=input_list, channels=markers_to_normalize, model=norm_model)

# extracting the legend and plots from the results
plot_list = densities_plot[setdiff(names(densities_plot), "legend")]
legend_plot = densities_plot$legend
plot_names = names(plot_list)

# keep only original and normalized plots
plot_list = plot_list[grep("original$|normalized$", plot_names)]  
plot_names = names(plot_list)

# remove plot titles and reduce axis text size
plot_list = lapply(plot_list, function(p) {
  p + 
    ggtitle(NULL) +
    theme(
      axis.title = element_text(size=5),
      axis.text = element_text(size = 5),     
      plot.margin = margin(6, 2, 2, 2)
    )
})


# save for each marker a png of before / after normalization
for (i in seq_along(markers_to_normalize)){
  marker = markers_to_normalize[i]
  outfile = bef_aft_plots[i]

  cat("Marker:", marker, "\n")
  cat("Saved to file:", outfile, "\n")

  # find the corresponding plots in plot_list
  original_name <- paste(marker, "original")
  normalized_name <- paste(marker, "normalized")

  if (!original_name %in% names(plot_list) ||
      !normalized_name %in% names(plot_list)) {
      stop(paste("Missing plots for marker:", marker))
    }

  before_plot <- plot_list[[original_name]]
  after_plot <- plot_list[[normalized_name]]

  # two column layout
  combined_plots = ggarrange(before_plot, after_plot, ncol=2, nrow=1)

  # save to outfile
  png(outfile, width=1400, height=700, res=400)
  print(combined_plots)
  dev.off()

}


class(legend_plot)

# save plot lgend
png(bef_aft_legend)
grid::grid.newpage()
grid::grid.draw(legend_plot)
dev.off()


end_time = Sys.time()
elapsed = end_time - start_time
cat("End:", date(), "\n")
cat("Total execution time:", elapsed)