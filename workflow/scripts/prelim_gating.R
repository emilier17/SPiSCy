#=============================================================================================#
# SPiSCy (Snakemake PIpeline for Spectral CYtometry)
#
# Phase: preprocessing
# Step: preliminary gating to remove debris and select live single cells
# Script 1/8
# 
# Author: Émilie Roy
# Date: Sept 2025
# Version: v.1.0
#
# Input : raw unmixed fcs file (.fcs)
# Outputs : gated fcs file (.fcs) + plots for validation (.png)
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
library(ggplot2)         # data visualization
library(ggpubr)          # for ggarrange
library(dplyr)           # sampling
library(gridExtra)       # plotting

# Snakemake inputs
fcs_file = snakemake@input[["fcs"]]
spe_config_file = snakemake@input[["specific_config"]]
gen_config_file = snakemake@input[["general_config"]]
output_fcs = snakemake@output[["fcs_out"]]
plot_output = snakemake@output[["gate_plots"]]

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


#### IMPORT FCS FILE INTO A FLOWFRAME ####

#--Import raw fcs file into a flowframe object

ff = read.FCS(
  filename = fcs_file,
  transformation = FALSE,
  truncate_max_range = FALSE,
  channel_alias = channel_alias_df)
#Truncate_max_range=FALSE, because max_range is adapted for standard FC values, not spectral.



#### HELPER FUNCTIONS ####

# Goes through user inputs for gates and makes sure they are properly formated for flowcore gates
create_gate = function(gate_def) {
  type = gate_def$type

  # RectangleGate
  if (type == "rectangle") {
    marker_x = as.character(gate_def$mk1)
    marker_y = if (!is.null(gate_def$mk2)) as.character(gate_def$mk2) else NULL
    args_list = list()
    args_list[[marker_x]] = as.numeric(gate_def$coords$mk1)
    if (!is.null(marker_y) && !is.null(gate_def$coords$mk2)) {
      args_list[[marker_y]] = as.numeric(gate_def$coords$mk2)
    }
    args_list$filterId = gate_def$filterId
    return(do.call(rectangleGate, args_list))
  
  # PolygonGate
  } else if (type == "polygon") {
    marker_x = as.character(gate_def$mk1)
    marker_y = as.character(gate_def$mk2)
    coords_mat = matrix(
      as.numeric(unlist(gate_def$coords)), 
      ncol = 2,
      byrow = TRUE
    )
    colnames(coords_mat) = c(marker_x, marker_y)
    return(polygonGate(filterId = gate_def$filterId, parameters = c(marker_x, marker_y), .gate = coords_mat))
  

   } else {
    stop(paste("Unknown gate type:", type))
  }
}

  #TODO : implement create_gate for all other gate types
   #else if (type == "ellipsoid") {
  #   covar = gate_def$covar
  #   mean_vec = gate_def$mean
  #   distance = gate_def$distance %||% 1  # default Mahalanobis distance
  #   if (!is.matrix(covar) || !is.numeric(mean_vec) || length(mean_vec) != ncol(covar)) {
  #     stop("ellipsoidGate requires a numeric covariance matrix and a mean vector of matching length")
  #   }
  #   return(ellipsoidGate(.gate = covar, mean = mean_vec, distance = distance, filterId = gate_def$filterId))


  # } else if (type == "quad") {
  #   if (is.null(gate_def$coords) || length(gate_def$coords) != 2) {
  #     stop("quadGate requires exactly two coordinates (thresholds) for the two-dimensional space.")
  #   }
  #   if (is.null(gate_def$x) || is.null(gate_def$y)) {
  #     stop("quadGate requires 'x' and 'y' channels to define the quadrant boundaries.")
  #   }
  #   gate_vec = setNames(unlist(gate_def$coords), c(gate_def$x, gate_def$y))
  #   return(quadGate(.gate = gate_vec, filterId = gate_def$filterId))
  
  
  # } else if (type == "polytope") {
  #   if (is.null(gate_def$.gate)) {
  #     stop("polytopeGate requires '.gate' parameter defining vertices in n dimensions.")
  #   }
  #   if (is.null(gate_def$b)) {
  #     stop("polytopeGate requires 'b' parameter representing half-space thresholds.")
  #   }
  #   gate_matrix = if (is.list(gate_def$.gate)) {
  #     do.call(rbind, gate_def$.gate)
  #   } else if (is.matrix(gate_def$.gate)) {
  #     gate_def$.gate
  #   } else {
  #     stop("polytopeGate '.gate' must be a list or matrix of vertices")
  #   }
  #   return(polytopeGate(.gate = gate_matrix, b = gate_def$b, filterId = gate_def$filterId))


# Plotting function (returns plot)
plot_gate = function(gate_name, mk1, mk2, percentile, fcs_after, fcs_before) {
  fcs_clean_df = as.data.frame(exprs(fcs_after))
  fcs_raw_df   = as.data.frame(exprs(fcs_before))
  
  # Subsample for speed
  sample_size = 100000  # max number of cells to visualize on the plot
  if (nrow(fcs_raw_df) > sample_size) {
    fcs_raw_df = sample_n(fcs_raw_df, sample_size)  #sampling
  }

  if (nrow(fcs_clean_df) > sample_size) {
    fcs_clean_df = sample_n(fcs_clean_df, sample_size)  #sampling
  }
  
  # Flag gated and non-gated cells
  fcs_clean_df$gated = TRUE
  fcs_raw_df$gated   = FALSE
  combined_df = rbind(fcs_raw_df, fcs_clean_df)

  # Setup axis names
  mk1_sym = sym(mk1)
  mk2_sym = sym(mk2)

  # Compute axis limits (percentiles)
  mk1_data = combined_df[[mk1]]
  mk2_data = combined_df[[mk2]]
  mk1_lim = quantile(mk1_data, probs = c(percentile, 1 - percentile), na.rm = TRUE)
  mk2_lim = quantile(mk2_data, probs = c(percentile, 1 - percentile), na.rm = TRUE)

  # Plotting
  plot = ggplot(combined_df, aes(x = !!mk1_sym, y = !!mk2_sym, color = as.factor(gated))) +
            geom_point(alpha = 0.5) +
            scale_color_manual(values = c("FALSE" = "gray", "TRUE" = "blue"),
                              labels = c("Excluded", "Included"),
                              name = "Cells") +
            labs(title = paste0(gate_name, " gate"),
                x = mk1,
                y = mk2) +
            coord_cartesian(xlim = mk1_lim, ylim = mk2_lim) +
            theme_minimal()
  
  return (plot)
}



#### PRELIMINARY GATING ####

#--Apply the gating sequence defined in gate_chains to the flowframe

fcs_raw = ff

#1) Selecting gating sequence
gates_sequence = spe_config$gate_chains$default

#2) Applying gating sequence
fcs_clean = fcs_raw
for (gate_name in gates_sequence) {
  gate_def = spe_config$gates[[gate_name]]
  cat("Applying gate:", gate_def$filterId, "|", "type:", gate_def$type, "\n")
  gate_obj = create_gate(gate_def)
  fcs_clean = Subset(fcs_clean, gate_obj)
}

#3) Generate all plots per gate for the fcs file
cat("Generating plots for:", basename(fcs_file), "\n")
all_plots = list()
for (gate_name in gates_sequence) {
  gate_def = spe_config$gates[[gate_name]]
  mk1 = gate_def$mk1
  mk2 = if (!is.null(gate_def$mk2)) gate_def$mk2 else "SSC-H"

  plot = plot_gate(gate_name, mk1, mk2, 0.01, fcs_clean, fcs_raw)
  all_plots[[gate_name]] = plot
}


#4) Create one png with all plots
nrows = ceiling(length(gates_sequence) / 2)
ncols = 2
png_width = 1500
png_height = 600 * nrows

all_plots = lapply(all_plots, function(p) {
  p + 
    theme(
      axis.text = element_text(size = 8),    # make axis and title font size smaller 
      axis.title = element_text(size = 15),
      plot.title = element_text(size = 20),
      plot.margin = margin(6, 2, 2, 2)
    )
})

combined_plots = ggarrange(
  plotlist = all_plots,
  nrow = nrows,
  ncol = ncols,
  labels = NULL)


png(plot_output, width=png_width, height=png_height)
print(combined_plots)
dev.off()


#5) Save gated flowframe as .fcs for next preprocessing step
cat("Saving fcs: ", "\n")
write.FCS(fcs_clean, filename=output_fcs)


end_time = Sys.time()
elapsed = end_time - start_time
cat("End:", date(), "\n")
cat("Total execution time:", elapsed)
