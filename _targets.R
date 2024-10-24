library(targets)
library(dplyr)
library(tidyr)
library(compositions)
# This is an example _targets.R file. Every
# {targets} pipeline needs one.
# Use tar_script() to create _targets.R and tar_edit()
# to open it again for editing.
# Then, run tar_make() to run the pipeline
# and tar_read(data_summary) to view the results.

# Define custom functions and other global objects.
# This is where you write source(\"R/functions.R\")
# if you keep your functions in external scripts.
source("scripts/read_data_KL15_XRF.R")
source("scripts/fit_density_pca.R")
source("scripts/fit_composition_pca.R")

# Set target-specific options such as packages:
# tar_option_set(packages = c("dplyr", "tidyr"))

# # Set the error option to "continue"
# tar_option_set(error = "continue")

# End this file with a list of target objects.
list(
  dir <- paste0(getwd(), "/data/Africa_NE_200/data/"),
  tar_target(data_odp_967_22, {
    tryCatch({
    data_odp_967_22 <- read.table(paste0(dir,'data_odp_967_22_43247_2021_339_MOESM2_ESM_XRF_Ti_Al.txt'), header = TRUE, sep = "\t")
    # get the age and Ti_Al attributes
    data_odp_967_22 <- data_odp_967_22[,c(8, 19)]
    # filter out data with age less than 200k
    data_odp_967_22 <- data_odp_967_22[data_odp_967_22[, 1] < 200,]
}, error = function(e) {
  message("Error in reading the file: ", e)
})
    data_odp_967_22
  }),
  # get the functionality from read_data_KL15_XRF.R
  tar_target(data_kl15_xrf, {
    tryCatch({
      data_kl15_xrf <- read.table(paste0(dir,'data_KL15_XRF.txt'), header = TRUE, sep = "\t")
      data_kl15_xrf <- rename(data_kl15_xrf, depth = User_ID)
    }, error = function(e) {
      message("Error in reading the file: ", e)
    })
  }),
  tar_target(data_kl15_agem, {
    tryCatch({
      data_kl15_agem <- read.table(paste0(dir,'data_KL15-2_smooth_spline_ages.txt'), header = TRUE, sep = "\t")
      data_kl15_agem <- rename(data_kl15_agem, age = best)
    }, error = function(e) {
      message("Error in reading the file: ", e)
    })
  }),
  tar_target(data_kl15_qf, {
    tryCatch({
      data_kl15_qf <- read.table(paste0(dir,'data_KL15_qf.txt'),
                                 header = TRUE, sep = "\t")
    }, error = function(e) {
      message("Error in reading the file: ", e)
    })
  }),
  tar_target(results, {
    results <- read_data_kl15_xrf(data_kl15_xrf, data_kl15_agem)
  }),
  tar_target(data_kl15_itpol, {
    data_kl15_itpol <- results$data_kl15_itpol
  }),
  tar_target(data_kl15, {
    data_kl15 <- results$data_kl15
  }),
  tar_target(missings_depth, {
    missings_depth <- results$missings_depth
  }),
    tar_target(data_kl15_comp, {
    data_comp <- results$data_comp
  }),
    tar_target(data_kl15_comp_clr, {
    data_clr <- results$data_clr
  }),
  tar_target(data_kl15_comp_ilr, {
    data_ilr <- results$data_ilr
  }),
  tar_target(data_kl15_comp_alr, {
    data_alr <- results$data_alr
  }),
  tar_target(simulation_densities_greven, {
    # define simulation parameters
    n_data <- 30
    n_samples <- 40
    x_grid_sim <- seq(0, 1, length.out = 1000)
    lambda_1 <- 0.5
    lambda_2 <- 0.2

    sim_densities_results <-
      simulate_densities_greven(
                                n_data,
                                n_samples,
                                x_grid = x_grid_sim,
                                lambda_1, lambda_2)
  }),
  tar_target(density_pca, {
    x_data <- simulation_densities_greven$x_data
    density_pca <- fit_density_pca(x_data, max_iter = 50)
  }),
  tar_target(plot_pca_results, {
    plot_pca(density_pca$pca, x_grid = density_pca$x_grid)
  }),
  tar_target(simulation_composition_1, {
    n_components <- 13
    n_data <- 30
    n_counts <- 2000
    n_samples <- 30
    lambda_1 <- 0.5
    lambda_2 <- 0.2

    sim_composition_1_results <-
      simulate_composition_1(
                             n_components,
                             n_data,
                             n_counts,
                             n_samples,
                             lambda_1, lambda_2)
  }),
  tar_target(composition_pca_sim_1, {
    x_data <- simulation_composition_1$x_data
    composition_pca_sim_1_results <- fit_compositional_pca_vs1_0(x_data, max_iter = 50)
  })
)
