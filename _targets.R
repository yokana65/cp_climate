library(targets)
library(dplyr)
library(tidyr)
library(compositions)
library(mvtnorm)
library(future.apply)

source("scripts/read_data_KL15_XRF.R")
source("scripts/fit_density_pca.R")
source("scripts/fit_composition_pca.R")
source("scripts/fit_composition_pca_ilr.R")
source("scripts/fit_pca_ilr.R")
source("scripts/fit_pca_clr.R")
source("scripts/fit_pca_ilr_2.R")
source("scripts/conditional_scores_function.R")
source("scripts/gradient.R")
source("scripts/simulations.R")
source("scripts/simulation_2.R")
source("scripts/simulation_3.R")
source("scripts/simulation_4.R")
source("scripts/simulation_5.R")
source("scripts/simulation_6.R")
source("scripts/simulation_7.R")
source("scripts/bootstrap.R")
source("scripts/helper_functions.R")

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
    n_data <- 40
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
    composition_pca_sim_1_results <-
      fit_compositional_pca_vs1_0(x_data, max_iter = 50)
  }),
  tar_target(pca_count_ilr_std, {
    x_data <- data_kl15_comp
    set.seed(10)
    pca_results_ilr_std <-
      fit_compositional_pca_ilr_sc(x_data,
                                   max_iter = 50,
                                   r = 10,
                                   lambda = 1,
                                   eps = 0.01,
                                   sc_factor = 1,
                                   sum_exp = TRUE)
  }),
  tar_target(sim_comp_1, {
    n_observations <- 2000
    eigenvalues <- c(0.6, 0.3, 0.05, 0.05)
    mean <- c(0, 2, 0.5, -2, -0.5)
    n_counts <- 300
    set.seed(13)

    sim_composition_1_results <-
      build_setting_2comp_5parts(n_observations,
                                 eigenvalues,
                                 mean,
                                 n_counts)
  }),
  tar_target(sim_comp_1_smi, {
    n_observations <- 2000
    eigenvalues <- c(0.6, 0.3, 0.05, 0.05)
    mean <- c(0, 2, 0.5, -2, -0.5) # the coordinates are too extreme -> or try less extreme
    n_counts <- 30

    set.seed(11)
    sim_composition_1_results <-
      build_setting_2comp_5parts(n_observations,
                                 eigenvalues,
                                 mean,
                                 n_counts)
  }),
  ##***** Sim Setting 1: build_setting_2comp_5parts***#
  tar_target(sim_comp_1_smi_nSim, {
    n_simulations <- 20
    n_observations <- 100
    eigenvalues <- c(0.6, 0.3, 0.05, 0.05)
    mean <- c(0, 1.5, 0.5, -1.5, -0.5) # the coordinates are too extreme -> or try less extreme
    n_counts <- 30
    sim_composition_1_results <- list(length(n_simulations))
    
    for (i in 1:n_simulations) {
      set.seed(1 * i)
      sim_composition_1_results[[i]] <-
        build_setting_2comp_5parts(n_observations,
                                 eigenvalues,
                                 mean,
                                 n_counts)
    }
    sim_composition_1_results
  }),
  tar_target(sim_comp_1_smi_nSim_300, {
    n_simulations <- 20
    n_observations <- 300
    eigenvalues <- c(0.6, 0.3, 0.05, 0.05)
    mean <- c(0, 1.5, 0.5, -1.5, -0.5) # the coordinates are too extreme -> or try less extreme
    n_counts <- 30
    sim_composition_1_results <- list(length(n_simulations))
    
    for (i in 1:n_simulations) {
      set.seed(1 * i)
      sim_composition_1_results[[i]] <-
        build_setting_2comp_5parts(n_observations,
                                 eigenvalues,
                                 mean,
                                 n_counts)
    }
    sim_composition_1_results
  }),
  tar_target(sim_comp_1_smi_nSim_60, {
    n_simulations <- 20
    n_observations <- 100
    eigenvalues <- c(0.6, 0.3, 0.05, 0.05)
    mean <- c(0, 1.5, 0.5, -1.5, -0.5) # the coordinates are too extreme -> or try less extreme
    n_counts <- 60
    sim_composition_1_results <- list(length(n_simulations))
    
    for (i in 1:n_simulations) {
      set.seed(1 * i)
      sim_composition_1_results[[i]] <-
        build_setting_2comp_5parts(n_observations,
                                 eigenvalues,
                                 mean,
                                 n_counts)
    }
    sim_composition_1_results
  }),
  tar_target(sim_comp_1_smi_nSim_90, {
    n_simulations <- 20
    n_observations <- 100
    eigenvalues <- c(0.6, 0.3, 0.05, 0.05)
    mean <- c(0, 1.5, 0.5, -1.5, -0.5) # the coordinates are too extreme -> or try less extreme
    n_counts <- 90
    sim_composition_1_results <- list(length(n_simulations))
    
    for (i in 1:n_simulations) {
      set.seed(1 * i)
      sim_composition_1_results[[i]] <-
        build_setting_2comp_5parts(n_observations,
                                 eigenvalues,
                                 mean,
                                 n_counts)
    }
    sim_composition_1_results
  }),
  tar_target(sim_comp_1_smi_nSim_90_300, {
    n_simulations <- 20
    n_observations <- 300
    eigenvalues <- c(0.6, 0.3, 0.05, 0.05)
    mean <- c(0, 1.5, 0.5, -1.5, -0.5) 
    n_counts <- 90
    sim_composition_1_results <- list(length(n_simulations))
    
    for (i in 1:n_simulations) {
      set.seed(1 * i)
      sim_composition_1_results[[i]] <-
        build_setting_2comp_5parts(n_observations,
                                 eigenvalues,
                                 mean,
                                 n_counts)
    }
    sim_composition_1_results
  }),
  ##***** Sim Setting 2: build_setting_2comp_5parts_vs2***#####################
  tar_target(sim_comp_2_smi_nSim, {
    n_simulations <- 20
    n_observations <- 100
    eigenvalues <- c(0.6, 0.3)
    mean <- c(0, 1, 0.5, -1, -0.5) 
    n_counts <- 30
    sim_composition_2_results <- list(length(n_simulations))
    
    for (i in 1:n_simulations) {
      set.seed(1 * i)
      sim_composition_2_results[[i]] <-
        build_setting_2comp_5parts_vs2(n_observations,
                                 eigenvalues,
                                 mean,
                                 n_counts)
    }
    sim_composition_2_results
  }),
  tar_target(sim_comp_2_smi_nSim_20, {
    n_simulations <- 100
    n_observations <- 100
    eigenvalues <- c(0.6, 0.3)
    mean <- c(0, 1, 0.5, -1, -0.5) 
    n_counts <- 20
    sim_composition_2_results <- list(length(n_simulations))
    
    for (i in 1:n_simulations) {
      set.seed(1 + i)
      sim_composition_2_results[[i]] <-
        build_setting_2comp_5parts_vs2(n_observations,
                                 eigenvalues,
                                 mean,
                                 n_counts)
    }
    sim_composition_2_results
  }),
  tar_target(sim_comp_2_smi_nSim_40, {
    n_simulations <- 100
    n_observations <- 100
    eigenvalues <- c(0.6, 0.3)
    mean <- c(0, 1, 0.5, -1, -0.5) 
    n_counts <- 40
    sim_composition_2_results <- list(length(n_simulations))
    
    for (i in 1:n_simulations) {
      set.seed(1 + i)
      sim_composition_2_results[[i]] <-
        build_setting_2comp_5parts_vs2(n_observations,
                                 eigenvalues,
                                 mean,
                                 n_counts)
    }
    sim_composition_2_results
  }),
  tar_target(sim_comp_2_smi_nSim_80, {
    n_simulations <- 100
    n_observations <- 100
    eigenvalues <- c(0.6, 0.3)
    mean <- c(0, 1, 0.5, -1, -0.5) 
    n_counts <- 80
    sim_composition_2_results <- list(length(n_simulations))
    
    for (i in 1:n_simulations) {
      set.seed(1 + i)
      sim_composition_2_results[[i]] <-
        build_setting_2comp_5parts_vs2(n_observations,
                                 eigenvalues,
                                 mean,
                                 n_counts)
    }
    sim_composition_2_results
  }),
  tar_target(sim_comp_2_smi_nSim_160, {
    n_simulations <- 100
    n_observations <- 100
    eigenvalues <- c(0.6, 0.3)
    mean <- c(0, 1, 0.5, -1, -0.5) 
    n_counts <- 160
    sim_composition_2_results <- list(length(n_simulations))
    
    for (i in 1:n_simulations) {
      set.seed(1 + i)
      sim_composition_2_results[[i]] <-
        build_setting_2comp_5parts_vs2(n_observations,
                                 eigenvalues,
                                 mean,
                                 n_counts)
    }
    sim_composition_2_results
  }),
  tar_target(sim_comp_2_smi_nSim_60, {
    n_simulations <- 20
    n_observations <- 100
    eigenvalues <- c(0.6, 0.3)
    mean <- c(0, 1, 0.5, -1, -0.5) 
    n_counts <- 60
    sim_composition_2_results <- list(length(n_simulations))
    
    for (i in 1:n_simulations) {
      set.seed(1 * i)
      sim_composition_2_results[[i]] <-
        build_setting_2comp_5parts_vs2(n_observations,
                                 eigenvalues,
                                 mean,
                                 n_counts)
    }
    sim_composition_2_results
  }),
  tar_target(sim_comp_2_smi_nSim_120, {
    n_simulations <- 20
    n_observations <- 100
    eigenvalues <- c(0.6, 0.3)
    mean <- c(0, 1, 0.5, -1, -0.5) 
    n_counts <- 120
    sim_composition_2_results <- list(length(n_simulations))
    
    for (i in 1:n_simulations) {
      set.seed(1 * i)
      sim_composition_2_results[[i]] <-
        build_setting_2comp_5parts_vs2(n_observations,
                                 eigenvalues,
                                 mean,
                                 n_counts)
    }
    sim_composition_2_results
  }),
  tar_target(sim_comp_2_smi_nSim_120_300, {
    n_simulations <- 20
    n_observations <- 300
    eigenvalues <- c(0.6, 0.3)
    mean <- c(0, 1, 0.5, -1, -0.5) 
    n_counts <- 120
    sim_composition_2_results <- list(length(n_simulations))
    
    for (i in 1:n_simulations) {
      set.seed(1 * i)
      sim_composition_2_results[[i]] <-
        build_setting_2comp_5parts_vs2(n_observations,
                                 eigenvalues,
                                 mean,
                                 n_counts)
    }
    sim_composition_2_results
  }),
  tar_target(sim_comp_2_smi_nSim_eb_20, {
    n_simulations <- 100
    n_observations <- 100
    eigenvalues <- c(0.6, 0.3)
    mean <- c(0, 0.05, 0.1, -0.05, -0.1) 
    n_counts <- 20
    sim_composition_2_results <- list(length(n_simulations))
    
    for (i in 1:n_simulations) {
      set.seed(1 * i)
      sim_composition_2_results[[i]] <-
        build_setting_2comp_5parts_vs2(n_observations,
                                 eigenvalues,
                                 mean,
                                 n_counts)
    }
    sim_composition_2_results
  }),
  tar_target(sim_comp_2_smi_nSim_eb_40, {
    n_simulations <- 100
    n_observations <- 100
    eigenvalues <- c(0.6, 0.3)
    mean <- c(0, 0.05, 0.1, -0.05, -0.1)
    n_counts <- 40
    sim_composition_2_results <- list(length(n_simulations))
    
    for (i in 1:n_simulations) {
      set.seed(1 * i)
      sim_composition_2_results[[i]] <-
        build_setting_2comp_5parts_vs2(n_observations,
                                 eigenvalues,
                                 mean,
                                 n_counts)
    }
    sim_composition_2_results
  }),
  tar_target(sim_comp_2_smi_nSim_eb_80, {
    n_simulations <- 100
    n_observations <- 100
    eigenvalues <- c(0.6, 0.3)
    mean <- c(0, 0.05, 0.1, -0.05, -0.1)
    n_counts <- 80
    sim_composition_2_results <- list(length(n_simulations))
    
    for (i in 1:n_simulations) {
      set.seed(1 * i)
      sim_composition_2_results[[i]] <-
        build_setting_2comp_5parts_vs2(n_observations,
                                 eigenvalues,
                                 mean,
                                 n_counts)
    }
    sim_composition_2_results
  }),
  tar_target(sim_comp_2_smi_nSim_eb_160, {
    n_simulations <- 100
    n_observations <- 100
    eigenvalues <- c(0.6, 0.3)
    mean <- c(0, 0.05, 0.1, -0.05, -0.1)
    n_counts <- 160
    sim_composition_2_results <- list(length(n_simulations))
    
    for (i in 1:n_simulations) {
      set.seed(1 * i)
      sim_composition_2_results[[i]] <-
        build_setting_2comp_5parts_vs2(n_observations,
                                 eigenvalues,
                                 mean,
                                 n_counts)
    }
    sim_composition_2_results
  }),
  tar_target(sim_comp_2_smi_nSim_eb_20_pi, {
    n_simulations <- 100
    n_observations <- 100
    eigenvalues <- c(0.6, 0.3)
    mean <- c(0, 0.05, 0.1, -0.05, -0.1)
    n_counts <- 20
    sim_composition_2_results <- list(length(n_simulations))
    
    for (i in 1:n_simulations) {
      set.seed(1 * i)
      sim_composition_2_results[[i]] <-
        build_setting_2comp_5parts_vs7(n_observations,
                                 eigenvalues,
                                 mean,
                                 n_counts)
    }
    sim_composition_2_results
  }),
  tar_target(sim_comp_2_smi_nSim_eb_40_pi, {
    n_simulations <- 100
    n_observations <- 100
    eigenvalues <- c(0.6, 0.3)
    mean <- c(0, 0.05, 0.1, -0.05, -0.1)
    n_counts <- 40
    sim_composition_2_results <- list(length(n_simulations))
    
    for (i in 1:n_simulations) {
      set.seed(1 * i)
      sim_composition_2_results[[i]] <-
        build_setting_2comp_5parts_vs7(n_observations,
                                 eigenvalues,
                                 mean,
                                 n_counts)
    }
    sim_composition_2_results
  }),
  tar_target(sim_comp_2_smi_nSim_eb_80_pi, {
    n_simulations <- 100
    n_observations <- 100
    eigenvalues <- c(0.6, 0.3)
    mean <- c(0, 0.05, 0.1, -0.05, -0.1)
    n_counts <- 80
    sim_composition_2_results <- list(length(n_simulations))
    
    for (i in 1:n_simulations) {
      set.seed(1 * i)
      sim_composition_2_results[[i]] <-
        build_setting_2comp_5parts_vs7(n_observations,
                                 eigenvalues,
                                 mean,
                                 n_counts)
    }
    sim_composition_2_results
  }),
  tar_target(sim_comp_2_smi_nSim_eb_160_pi, {
    n_simulations <- 100
    n_observations <- 100
    eigenvalues <- c(0.6, 0.3)
    mean <- c(0, 0.05, 0.1, -0.05, -0.1)
    n_counts <- 160
    sim_composition_2_results <- list(length(n_simulations))
    
    for (i in 1:n_simulations) {
      set.seed(1 * i)
      sim_composition_2_results[[i]] <-
        build_setting_2comp_5parts_vs7(n_observations,
                                 eigenvalues,
                                 mean,
                                 n_counts)
    }
    sim_composition_2_results
  }),
  ###********Setting 3*******************######################################
  tar_target(sim_3_smi_nSim_20, {
    n_simulations <- 60
    n_observations <- 100
    eigenvalues = c(0.5, 0.3, 0.1, 0.1)
    # mean = c(0, 0.9, 0.3, -0.8, -0.2)
    mean = c(0, 1.5, 0.5, -1.5, -0.5)
    n_counts <- 20
    sim_composition_3_results <- list(length(n_simulations))
    
    for (i in 1:n_simulations) {
      set.seed(1 * i)
      sim_composition_3_results[[i]] <-
        build_setting_4comp_5parts_vs2(n_observations,
                                 eigenvalues,
                                 mean,
                                 n_counts)
    }
    sim_composition_3_results
  }),
  tar_target(sim_3_smi_nSim_30, {
    n_simulations <- 20
    n_observations <- 100
    eigenvalues = c(0.5, 0.3, 0.1, 0.1)
    # mean = c(0, 0.9, 0.3, -0.8, -0.2)
    mean = c(0, 1.5, 0.5, -1.5, -0.5)
    n_counts <- 30
    sim_composition_3_results <- list(length(n_simulations))
    
    for (i in 1:n_simulations) {
      set.seed(1 * i)
      sim_composition_3_results[[i]] <-
        build_setting_4comp_5parts_vs2(n_observations,
                                 eigenvalues,
                                 mean,
                                 n_counts)
    }
    sim_composition_3_results
  }),
  tar_target(sim_3_smi_nSim_40, {
    n_simulations <- 60
    n_observations <- 100
    eigenvalues = c(0.5, 0.3, 0.1, 0.1)
    # mean = c(0, 0.9, 0.3, -0.8, -0.2)
    mean = c(0, 1.5, 0.5, -1.5, -0.5)
    n_counts <- 40
    sim_composition_3_results <- list(length(n_simulations))
    
    for (i in 1:n_simulations) {
      set.seed(1 * i)
      sim_composition_3_results[[i]] <-
        build_setting_4comp_5parts_vs2(n_observations,
                                 eigenvalues,
                                 mean,
                                 n_counts)
    }
    sim_composition_3_results
  }),
  tar_target(sim_3_smi_nSim_80, {
    n_simulations <- 60
    n_observations <- 100
    eigenvalues = c(0.5, 0.3, 0.1, 0.1)
    # mean = c(0, 0.9, 0.3, -0.8, -0.2)
    mean = c(0, 1.5, 0.5, -1.5, -0.5)
    n_counts <- 80
    sim_composition_3_results <- list(length(n_simulations))
    
    for (i in 1:n_simulations) {
      set.seed(1 * i)
      sim_composition_3_results[[i]] <-
        build_setting_4comp_5parts_vs2(n_observations,
                                 eigenvalues,
                                 mean,
                                 n_counts)
    }
    sim_composition_3_results
  }),
  tar_target(sim_3_smi_nSim_160, {
    n_simulations <- 60
    n_observations <- 100
    eigenvalues = c(0.5, 0.3, 0.1, 0.1)
    # mean = c(0, 0.9, 0.3, -0.8, -0.2)
    mean = c(0, 1.5, 0.5, -1.5, -0.5)
    n_counts <- 160
    sim_composition_3_results <- list(length(n_simulations))
    
    for (i in 1:n_simulations) {
      set.seed(1 * i)
      sim_composition_3_results[[i]] <-
        build_setting_4comp_5parts_vs2(n_observations,
                                 eigenvalues,
                                 mean,
                                 n_counts)
    }
    sim_composition_3_results
  }),
  tar_target(sim_3_smi_nSim_60, {
    n_simulations <- 20
    n_observations <- 100
    eigenvalues = c(0.5, 0.3, 0.1, 0.1)
    # mean = c(0, 0.9, 0.3, -0.8, -0.2)
    mean = c(0, 1.5, 0.5, -1.5, -0.5)
    n_counts <- 60
    sim_composition_3_results <- list(length(n_simulations))
    
    for (i in 1:n_simulations) {
      set.seed(1 * i)
      sim_composition_3_results[[i]] <-
        build_setting_4comp_5parts_vs2(n_observations,
                                 eigenvalues,
                                 mean,
                                 n_counts)
    }
    sim_composition_3_results
  }),
  tar_target(sim_3_smi_nSim_120, {
    n_simulations <- 20
    n_observations <- 100
    eigenvalues = c(0.5, 0.3, 0.1, 0.1)
    # mean = c(0, 0.9, 0.3, -0.8, -0.2)
    mean = c(0, 1.5, 0.5, -1.5, -0.5)
    n_counts <- 120
    sim_composition_3_results <- list(length(n_simulations))
    
    for (i in 1:n_simulations) {
      set.seed(1 * i)
      sim_composition_3_results[[i]] <-
        build_setting_4comp_5parts_vs2(n_observations,
                                 eigenvalues,
                                 mean,
                                 n_counts)
    }
    sim_composition_3_results
  }),
  tar_target(sim_3_smi_nSim_400, {
    n_simulations <- 20
    n_observations <- 100
    eigenvalues = c(0.5, 0.3, 0.1, 0.1)
    # mean = c(0, 0.9, 0.3, -0.8, -0.2)
    mean = c(0, 1.5, 0.5, -1.5, -0.5)
    n_counts <- 400
    sim_composition_3_results <- list(length(n_simulations))
    
    for (i in 1:n_simulations) {
      set.seed(1 * i)
      sim_composition_3_results[[i]] <-
        build_setting_4comp_5parts_vs2(n_observations,
                                 eigenvalues,
                                 mean,
                                 n_counts)
    }
    sim_composition_3_results
  }),
  ##***** Sim Setting 4: build_setting_4comp_5parts_vs3***###############
  tar_target(sim_4_smi_nSim_30, {
    n_simulations <- 20
    n_observations <- 100
    eigenvalues = c(0.5, 0.3, 0.15, 0.05)
    mean = c(0, 0.9, 0.3, -0.8, -0.2)
    n_counts <- 30
    sim_composition_4_results <- list(length(n_simulations))
    
    for (i in 1:n_simulations) {
      set.seed(1 * i)
      sim_composition_4_results[[i]] <-
        build_setting_4comp_5parts_vs3(n_observations,
                                 eigenvalues,
                                 mean,
                                 n_counts)
    }
    sim_composition_4_results
  }),
  tar_target(sim_4_smi_nSim_60, {
    n_simulations <- 20
    n_observations <- 100
    eigenvalues = c(0.5, 0.3, 0.15, 0.05)
    mean = c(0, 0.9, 0.3, -0.8, -0.2)
    n_counts <- 60
    sim_composition_4_results <- list(length(n_simulations))
    
    for (i in 1:n_simulations) {
      set.seed(1 * i)
      sim_composition_4_results[[i]] <-
        build_setting_4comp_5parts_vs3(n_observations,
                                 eigenvalues,
                                 mean,
                                 n_counts)
    }
    sim_composition_4_results
  }),
  tar_target(sim_4_smi_nSim_120, {
    n_simulations <- 20
    n_observations <- 100
    eigenvalues = c(0.5, 0.3, 0.15, 0.05)
    mean = c(0, 0.9, 0.3, -0.8, -0.2)
    n_counts <- 120
    sim_composition_4_results <- list(length(n_simulations))
    
    for (i in 1:n_simulations) {
      set.seed(1 * i)
      sim_composition_4_results[[i]] <-
        build_setting_4comp_5parts_vs3(n_observations,
                                 eigenvalues,
                                 mean,
                                 n_counts)
    }
    sim_composition_4_results
  }),
  tar_target(sim_4_smi_nSim_400, {
    n_simulations <- 20
    n_observations <- 100
    eigenvalues = c(0.5, 0.3, 0.15, 0.05)
    mean = c(0, 0.9, 0.3, -0.8, -0.2)
    n_counts <- 400
    sim_composition_4_results <- list(length(n_simulations))
    
    for (i in 1:n_simulations) {
      set.seed(1 * i)
      sim_composition_4_results[[i]] <-
        build_setting_4comp_5parts_vs3(n_observations,
                                 eigenvalues,
                                 mean,
                                 n_counts)
    }
    sim_composition_4_results
  }),
  ###### Sim Setting 5: build_setting_4comp_13parts_vs1 ################
  tar_target(sim_5_01, {
    n_simulations <- 20
    data <- data_kl15
    eigenvalues = c(0.5, 0.2, 0.12, 0.08)
    mean = c(-1.257, -1.741, 1.713, -0.6,
             1.07, -3.624, -1.254, 1.05, -0.9,
             0.505, 4.0, -0.49, 1.532)
    scale = 0.01
    sim_composition_5_results <- list(length(n_simulations))
    
    for (i in 1:n_simulations) {
      set.seed(1 * i)
      sim_composition_5_results[[i]] <-
        build_setting_4comp_13parts_vs1(data,
                                 eigenvalues,
                                 mean,
                                 scale)
    }
    sim_composition_5_results
  }),
  tar_target(sim_5_001, {
    n_simulations <- 20
    data <- data_kl15
    eigenvalues = c(0.5, 0.2, 0.12, 0.08)
    mean = c(-1.257, -1.741, 1.713, -0.6,
             1.07, -3.624, -1.254, 1.05, -0.9,
             0.505, 4.0, -0.49, 1.532)
    scale = 0.001
    sim_composition_5_results <- list(length(n_simulations))
    
    for (i in 1:n_simulations) {
      set.seed(1 * i)
      sim_composition_5_results[[i]] <-
        build_setting_4comp_13parts_vs1(data,
                                 eigenvalues,
                                 mean,
                                 scale)
    }
    sim_composition_5_results
  }),
  tar_target(sim_5_05, {
    n_simulations <- 20
    data <- data_kl15
    eigenvalues = c(0.5, 0.2, 0.12, 0.08)
    mean = c(-1.257, -1.741, 1.713, -0.6,
             1.07, -3.624, -1.254, 1.05, -0.9,
             0.505, 4.0, -0.49, 1.532)
    scale = 0.05
    sim_composition_5_results <- list(length(n_simulations))
    
    for (i in 1:n_simulations) {
      set.seed(1 * i)
      sim_composition_5_results[[i]] <-
        build_setting_4comp_13parts_vs1(data,
                                 eigenvalues,
                                 mean,
                                 scale)
    }
    sim_composition_5_results
  }),
  tar_target(sim_5_0001, {
    n_simulations <- 20
    data <- data_kl15
    eigenvalues = c(0.5, 0.2, 0.12, 0.08)
    mean = c(-1.257, -1.741, 1.713, -0.6,
             1.07, -3.624, -1.254, 1.05, -0.9,
             0.505, 4.0, -0.49, 1.532)
    scale = 0.0001
    sim_composition_5_results <- list(length(n_simulations))
    
    for (i in 1:n_simulations) {
      set.seed(1 * i)
      sim_composition_5_results[[i]] <-
        build_setting_4comp_13parts_vs1(data,
                                 eigenvalues,
                                 mean,
                                 scale)
    }
    sim_composition_5_results
  }),
  ###### Sim Setting 6: build_setting_4comp_13parts_vs2 ################
  tar_target(sim_6_sc1_n50, {
    n_simulations <- 20
    n_observations <- 50
    data <- data_kl15
    eigenvalues = c(0.5, 0.2, 0.12, 0.08)
    mean = c(-1.257, -1.741, 1.713, -0.6,
             1.07, -3.624, -1.254, 1.05, -0.9,
             0.505, 4.0, -0.49, 1.532)
    scale = 1
    sim_composition_6_results <- list(length(n_simulations))
    
    for (i in 1:n_simulations) {
      set.seed(1 * i)
      sim_composition_6_results[[i]] <-
        build_setting_4comp_13parts_vs2(data,
                                 eigenvalues,
                                 mean,
                                 scale,
                                 n_observations)
    }
    sim_composition_6_results
  }),
  tar_target(sim_6_sc01_n50, {
    n_simulations <- 20
    n_observations <- 50
    data <- data_kl15
    eigenvalues = c(0.5, 0.2, 0.12, 0.08)
    mean = c(-1.257, -1.741, 1.713, -0.6,
             1.07, -3.624, -1.254, 1.05, -0.9,
             0.505, 4.0, -0.49, 1.532)
    scale = 0.1
    sim_composition_6_results <- list(length(n_simulations))
    
    for (i in 1:n_simulations) {
      set.seed(1 * i)
      sim_composition_6_results[[i]] <-
        build_setting_4comp_13parts_vs2(data,
                                 eigenvalues,
                                 mean,
                                 scale,
                                 n_observations)
    }
    sim_composition_6_results
  }),
  tar_target(sim_6_sc001_n50, {
    n_simulations <- 20
    n_observations <- 50
    data <- data_kl15
    eigenvalues = c(0.5, 0.2, 0.12, 0.08)
    mean = c(-1.257, -1.741, 1.713, -0.6,
             1.07, -3.624, -1.254, 1.05, -0.9,
             0.505, 4.0, -0.49, 1.532)
    scale = 0.01
    sim_composition_6_results <- list(length(n_simulations))
    
    for (i in 1:n_simulations) {
      set.seed(1 * i)
      sim_composition_6_results[[i]] <-
        build_setting_4comp_13parts_vs2(data,
                                 eigenvalues,
                                 mean,
                                 scale,
                                 n_observations)
    }
    sim_composition_6_results
  }),
  tar_target(sim_6_01_n100, {
    n_simulations <- 20
    n_observations <- 100
    data <- data_kl15
    eigenvalues = c(0.5, 0.2, 0.12, 0.08)
    mean = c(-1.257, -1.741, 1.713, -0.6,
             1.07, -3.624, -1.254, 1.05, -0.9,
             0.505, 4.0, -0.49, 1.532)
    scale = 0.01
    sim_composition_6_results <- list(length(n_simulations))
    
    for (i in 1:n_simulations) {
      set.seed(1 * i)
      sim_composition_6_results[[i]] <-
        build_setting_4comp_13parts_vs2(data,
                                 eigenvalues,
                                 mean,
                                 scale,
                                 n_observations)
    }
    sim_composition_6_results
  }),
  tar_target(sim_6_sc1_n100, {
    n_simulations <- 20
    n_observations <- 100
    data <- data_kl15
    eigenvalues = c(0.5, 0.2, 0.12, 0.08)
    mean = c(-1.257, -1.741, 1.713, -0.6,
             1.07, -3.624, -1.254, 1.05, -0.9,
             0.505, 4.0, -0.49, 1.532)
    scale = 1
    sim_composition_6_results <- list(length(n_simulations))
    
    for (i in 1:n_simulations) {
      set.seed(1 * i)
      sim_composition_6_results[[i]] <-
        build_setting_4comp_13parts_vs2(data,
                                 eigenvalues,
                                 mean,
                                 scale,
                                 n_observations)
    }
    sim_composition_6_results
  }),
  tar_target(sim_6_sc01_n100, {
    n_simulations <- 20
    n_observations <- 100
    data <- data_kl15
    eigenvalues = c(0.5, 0.2, 0.12, 0.08)
    mean = c(-1.257, -1.741, 1.713, -0.6,
             1.07, -3.624, -1.254, 1.05, -0.9,
             0.505, 4.0, -0.49, 1.532)
    scale = 0.1
    sim_composition_6_results <- list(length(n_simulations))
    
    for (i in 1:n_simulations) {
      set.seed(1 * i)
      sim_composition_6_results[[i]] <-
        build_setting_4comp_13parts_vs2(data,
                                 eigenvalues,
                                 mean,
                                 scale,
                                 n_observations)
    }
    sim_composition_6_results
  }),
  tar_target(sim_6_sc001_n100, {
    n_simulations <- 20
    n_observations <- 100
    data <- data_kl15
    eigenvalues = c(0.5, 0.2, 0.12, 0.08)
    mean = c(-1.257, -1.741, 1.713, -0.6,
             1.07, -3.624, -1.254, 1.05, -0.9,
             0.505, 4.0, -0.49, 1.532)
    scale = 0.01
    sim_composition_6_results <- list(length(n_simulations))
    
    for (i in 1:n_simulations) {
      set.seed(1 * i)
      sim_composition_6_results[[i]] <-
        build_setting_4comp_13parts_vs2(data,
                                 eigenvalues,
                                 mean,
                                 scale,
                                 n_observations)
    }
    sim_composition_6_results
  }),
  tar_target(sim_6_sc0001_n100, {
    n_simulations <- 20
    n_observations <- 100
    data <- data_kl15
    eigenvalues = c(0.5, 0.2, 0.12, 0.08)
    mean = c(-1.257, -1.741, 1.713, -0.6,
             1.07, -3.624, -1.254, 1.05, -0.9,
             0.505, 4.0, -0.49, 1.532)
    scale = 0.001
    sim_composition_6_results <- list(length(n_simulations))
    
    for (i in 1:n_simulations) {
      set.seed(1 * i)
      sim_composition_6_results[[i]] <-
        build_setting_4comp_13parts_vs2(data,
                                 eigenvalues,
                                 mean,
                                 scale,
                                 n_observations)
    }
    sim_composition_6_results
  }),
  tar_target(sim_6_01_n200, {
    n_simulations <- 20
    n_observations <- 200
    data <- data_kl15
    eigenvalues = c(0.5, 0.2, 0.12, 0.08)
    mean = c(-1.257, -1.741, 1.713, -0.6,
             1.07, -3.624, -1.254, 1.05, -0.9,
             0.505, 4.0, -0.49, 1.532)
    scale = 0.01
    sim_composition_6_results <- list(length(n_simulations))
    
    for (i in 1:n_simulations) {
      set.seed(1 * i)
      sim_composition_6_results[[i]] <-
        build_setting_4comp_13parts_vs2(data,
                                 eigenvalues,
                                 mean,
                                 scale,
                                 n_observations)
    }
    sim_composition_6_results
  }),
  tar_target(sim_6_01_n400, {
    n_simulations <- 20
    n_observations <- 400
    data <- data_kl15
    eigenvalues = c(0.5, 0.2, 0.12, 0.08)
    mean = c(-1.257, -1.741, 1.713, -0.6,
             1.07, -3.624, -1.254, 1.05, -0.9,
             0.505, 4.0, -0.49, 1.532)
    scale = 0.01
    sim_composition_6_results <- list(length(n_simulations))
    
    for (i in 1:n_simulations) {
      set.seed(1 * i)
      sim_composition_6_results[[i]] <-
        build_setting_4comp_13parts_vs2(data,
                                 eigenvalues,
                                 mean,
                                 scale,
                                 n_observations)
    }
    sim_composition_6_results
  }),
  tar_target(sim_6_01_n800, {
    n_simulations <- 20
    n_observations <- 800
    data <- data_kl15
    eigenvalues = c(0.5, 0.2, 0.12, 0.08)
    mean = c(-1.257, -1.741, 1.713, -0.6,
             1.07, -3.624, -1.254, 1.05, -0.9,
             0.505, 4.0, -0.49, 1.532)
    scale = 0.01
    sim_composition_6_results <- list(length(n_simulations))
    
    for (i in 1:n_simulations) {
      set.seed(1 * i)
      sim_composition_6_results[[i]] <-
        build_setting_4comp_13parts_vs2(data,
                                 eigenvalues,
                                 mean,
                                 scale,
                                 n_observations)
    }
    sim_composition_6_results
  }),
  tar_target(sim_6_sc1_n200, {
    n_simulations <- 20
    n_observations <- 200
    data <- data_kl15
    eigenvalues = c(0.5, 0.2, 0.12, 0.08)
    mean = c(-1.257, -1.741, 1.713, -0.6,
             1.07, -3.624, -1.254, 1.05, -0.9,
             0.505, 4.0, -0.49, 1.532)
    scale = 1
    sim_composition_6_results <- list(length(n_simulations))
    
    for (i in 1:n_simulations) {
      set.seed(1 * i)
      sim_composition_6_results[[i]] <-
        build_setting_4comp_13parts_vs2(data,
                                 eigenvalues,
                                 mean,
                                 scale,
                                 n_observations)
    }
    sim_composition_6_results
  }),
  tar_target(sim_6_sc1_n400, {
    n_simulations <- 20
    n_observations <- 400
    data <- data_kl15
    eigenvalues = c(0.5, 0.2, 0.12, 0.08)
    mean = c(-1.257, -1.741, 1.713, -0.6,
             1.07, -3.624, -1.254, 1.05, -0.9,
             0.505, 4.0, -0.49, 1.532)
    scale = 1
    sim_composition_6_results <- list(length(n_simulations))
    
    for (i in 1:n_simulations) {
      set.seed(1 * i)
      sim_composition_6_results[[i]] <-
        build_setting_4comp_13parts_vs2(data,
                                 eigenvalues,
                                 mean,
                                 scale,
                                 n_observations)
    }
    sim_composition_6_results
  }),
  tar_target(sim_6_sc1_n800, {
    n_simulations <- 20
    n_observations <- 800
    data <- data_kl15
    eigenvalues = c(0.5, 0.2, 0.12, 0.08)
    mean = c(-1.257, -1.741, 1.713, -0.6,
             1.07, -3.624, -1.254, 1.05, -0.9,
             0.505, 4.0, -0.49, 1.532)
    scale = 1
    sim_composition_6_results <- list(length(n_simulations))
    
    for (i in 1:n_simulations) {
      set.seed(1 * i)
      sim_composition_6_results[[i]] <-
        build_setting_4comp_13parts_vs2(data,
                                 eigenvalues,
                                 mean,
                                 scale,
                                 n_observations)
    }
    sim_composition_6_results
  }),
  tar_target(sim_6_sc1_n1600, {
    n_simulations <- 20
    n_observations <- 1600
    data <- data_kl15
    eigenvalues = c(0.5, 0.2, 0.12, 0.08)
    mean = c(-1.257, -1.741, 1.713, -0.6,
             1.07, -3.624, -1.254, 1.05, -0.9,
             0.505, 4.0, -0.49, 1.532)
    scale = 1
    sim_composition_6_results <- list(length(n_simulations))
    
    for (i in 1:n_simulations) {
      set.seed(1 * i)
      sim_composition_6_results[[i]] <-
        build_setting_4comp_13parts_vs2(data,
                                 eigenvalues,
                                 mean,
                                 scale,
                                 n_observations)
    }
    sim_composition_6_results
  }),
  tar_target(sim_6_sc1_n2000, {
    n_simulations <- 20
    n_observations <- 2000
    data <- data_kl15
    eigenvalues = c(0.5, 0.2, 0.12, 0.08)
    mean = c(-1.257, -1.741, 1.713, -0.6,
             1.07, -3.624, -1.254, 1.05, -0.9,
             0.505, 4.0, -0.49, 1.532)
    scale = 1
    sim_composition_6_results <- list(length(n_simulations))
    
    for (i in 1:n_simulations) {
      set.seed(1 * i)
      sim_composition_6_results[[i]] <-
        build_setting_4comp_13parts_vs2(data,
                                 eigenvalues,
                                 mean,
                                 scale,
                                 n_observations)
    }
    sim_composition_6_results
  }),
  ######################################################################
  ##### Simulation runs ################################################
  ##********* Setting 6**********#
  tar_target(pca_6_01_n100, {
    number_simulations <- length(sim_6_01_n100)
    pca_results_list <- list(length(sim_6_01_n100))
    for  (i in 1:length(sim_6_01_n100)) {
      sim <-  sim_6_01_n100[[i]]
      x_data <- sim$x_data
      cat("Iteration sim_6_01_n100", i, "\n")

      set.seed(1 * i)
      result <- fit_pca_ilr_vs_2(
          x_data, 
          sc_factor = 1,
          max_iter = 80,
          eps = 0.04
      )
      pca_results_list[[i]] <- result
    }
    pca_results_list
  }),
  tar_target(pca_6_01_n200, {
    number_simulations <- length(sim_6_01_n200)
    pca_results_list <- list(length(sim_6_01_n200))
    for  (i in 1:length(sim_6_01_n200)) {
      sim <-  sim_6_01_n200[[i]]
      x_data <- sim$x_data
      cat("Iteration sim_6_01_n200", i, "\n")

      set.seed(1 * i)
      result <- fit_pca_ilr_vs_2(
          x_data, 
          sc_factor = 1,
          max_iter = 80,
          eps = 0.03
      )
      pca_results_list[[i]] <- result
    }
    pca_results_list
  }),
  tar_target(pca_6_01_n400, {
    number_simulations <- length(sim_6_01_n400)
    pca_results_list <- list(length(sim_6_01_n400))
    for  (i in 1:length(sim_6_01_n400)) {
      sim <-  sim_6_01_n400[[i]]
      x_data <- sim$x_data
      cat("Iteration sim_6_01_n400", i, "\n")

      set.seed(1 * i)
      result <- fit_pca_ilr_vs_2(
          x_data, 
          sc_factor = 1,
          r = 15,
          max_iter = 80,
          eps = 0.04
      )
      pca_results_list[[i]] <- result
    }
    pca_results_list
  }),
  tar_target(pca_6_01_n800, {
    number_simulations <- length(sim_6_01_n800)
    pca_results_list <- list(length(sim_6_01_n800))
    for  (i in 1:length(sim_6_01_n800)) {
      sim <-  sim_6_01_n800[[i]]
      x_data <- sim$x_data
      cat("Iteration sim_6_01_n800", i, "\n")

      set.seed(1 * i)
      result <- fit_pca_ilr_vs_2(
          x_data, 
          sc_factor = 1,
          r = 15,
          max_iter = 80,
          eps = 0.04
      )
      pca_results_list[[i]] <- result
    }
    pca_results_list
  }),
  tar_target(pca_6_sc1_n200, {
    number_simulations <- length(sim_6_sc1_n200)
    pca_results_list <- list(length(sim_6_sc1_n200))
    for  (i in 1:length(sim_6_sc1_n200)) {
      sim <-  sim_6_sc1_n200[[i]]
      x_data <- sim$x_data
      cat("Iteration sim_6_sc1_n200", i, "\n")

      set.seed(1 * i)
      result <- fit_pca_ilr_vs_2(
          x_data, 
          sc_factor = 1,
          max_iter = 80,
          eps = 0.05
      )
      pca_results_list[[i]] <- result
    }
    pca_results_list
  }),
  tar_target(pca_6_sc1_n400, {
    number_simulations <- length(sim_6_sc1_n400)
    pca_results_list <- list(length(sim_6_sc1_n400))
    for  (i in 1:length(sim_6_sc1_n400)) {
      sim <-  sim_6_sc1_n400[[i]]
      x_data <- sim$x_data
      cat("Iteration sim_6_sc1_n400", i, "\n")

      set.seed(1 * i)
      result <- fit_pca_ilr_vs_2(
          x_data, 
          sc_factor = 1,
          max_iter = 80,
          eps = 0.05
      )
      pca_results_list[[i]] <- result
    }
    pca_results_list
  }),
  tar_target(pca_6_sc1_n800, {
    number_simulations <- length(sim_6_sc1_n800)
    pca_results_list <- list(length(sim_6_sc1_n800))
    for  (i in 1:length(sim_6_sc1_n800)) {
      sim <-  sim_6_sc1_n800[[i]]
      x_data <- sim$x_data
      cat("Iteration sim_6_sc1_n800", i, "\n")

      set.seed(1 * i)
      result <- fit_pca_ilr_vs_2(
          x_data, 
          sc_factor = 1,
          r = 15,
          max_iter = 80,
          eps = 0.05
      )
      pca_results_list[[i]] <- result
    }
    pca_results_list
  }),
  tar_target(pca_6_sc1_n1600, {
    number_simulations <- length(sim_6_sc1_n1600)
    pca_results_list <- list(length(sim_6_sc1_n1600))
    for  (i in 1:length(sim_6_sc1_n1600)) {
      sim <-  sim_6_sc1_n1600[[i]]
      x_data <- sim$x_data
      cat("Iteration sim_6_sc1_n1600", i, "\n")

      set.seed(1 * i)
      result <- fit_pca_ilr_vs_2(
          x_data, 
          sc_factor = 1,
          r = 15,
          max_iter = 80,
          eps = 0.05
      )
      pca_results_list[[i]] <- result
    }
    pca_results_list
  }),
  tar_target(pca_6_sc1_n2000, {
    number_simulations <- length(sim_6_sc1_n2000)
    pca_results_list <- list(length(sim_6_sc1_n2000))
    for  (i in 1:length(sim_6_sc1_n2000)) {
      sim <-  sim_6_sc1_n2000[[i]]
      x_data <- sim$x_data
      cat("Iteration sim_6_sc1_n2000", i, "\n")

      set.seed(1 * i)
      result <- fit_pca_ilr_vs_2(
        x_data, 
        sc_factor = 1,
        r = 15,
        max_iter = 80,
        eps = 0.05
      )
      pca_results_list[[i]] <- result
    }
    pca_results_list
  }),
  tar_target(pca_6_sc01_n2000, {
    number_simulations <- length(sim_6_sc1_n2000)
    pca_results_list <- list(length(sim_6_sc1_n2000))
    for  (i in 1:length(sim_6_sc1_n2000)) {
      sim <-  sim_6_sc1_n2000[[i]]
      x_data <- sim$x_data
      scaled_data <- lapply(x_data, function(x) x * 0.1)
      cat("Iteration sim_6_sc01_n2000", i, "\n")

      set.seed(1 * i)
      result <- fit_pca_ilr_vs_2(
        scaled_data, 
        sc_factor = 1,
        r = 15,
        max_iter = 80,
        eps = 0.05
      )
      pca_results_list[[i]] <- result
    }
    pca_results_list
  }),
  tar_target(pca_6_sc001_n2000, {
    number_simulations <- length(sim_6_sc1_n2000)
    pca_results_list <- list(length(sim_6_sc1_n2000))
    for  (i in 1:length(sim_6_sc1_n2000)) {
      sim <-  sim_6_sc1_n2000[[i]]
      x_data <- sim$x_data
      scaled_data <- lapply(x_data, function(x) x * 0.01)
      cat("Iteration sim_6_sc001_n2000", i, "\n")

      set.seed(1 * i)
      result <- fit_pca_ilr_vs_2(
        scaled_data, 
        sc_factor = 1,
        r = 15,
        max_iter = 80,
        eps = 0.05
      )
      pca_results_list[[i]] <- result
    }
    pca_results_list
  }),
  tar_target(pca_6_sc0001_n2000, {
    number_simulations <- length(sim_6_sc1_n2000)
    pca_results_list <- list(length(sim_6_sc1_n2000))
    for  (i in 1:length(sim_6_sc1_n2000)) {
      sim <-  sim_6_sc1_n2000[[i]]
      x_data <- sim$x_data
      scaled_data <- lapply(x_data, function(x) x * 0.001)
      cat("Iteration sim_6_sc0001_n2000", i, "\n")

      set.seed(1 * i)
      result <- fit_pca_ilr_vs_2(
        scaled_data, 
        sc_factor = 1,
        r = 15,
        max_iter = 80,
        eps = 0.05
      )
      pca_results_list[[i]] <- result
    }
    pca_results_list
  }),
  tar_target(pca_6_sc1_n100, {
    number_simulations <- length(sim_6_sc1_n100)
    pca_results_list <- list(length(sim_6_sc1_n100))
    for  (i in 1:length(sim_6_sc1_n100)) {
      sim <-  sim_6_sc1_n100[[i]]
      x_data <- sim$x_data
      cat("Iteration sim_6_sc1_n100", i, "\n")

      set.seed(1 * i)
      result <- fit_pca_ilr_vs_2(
          x_data, 
          sc_factor = 1,
          r = 15,
          max_iter = 80,
          eps = 0.07
      )
      pca_results_list[[i]] <- result
    }
    pca_results_list
  }),
  tar_target(pca_6_sc01_n100, {
    number_simulations <- length(sim_6_sc1_n100)
    pca_results_list <- list(length(sim_6_sc1_n100))
    for  (i in 1:length(sim_6_sc1_n100)) {
      sim <-  sim_6_sc1_n100[[i]]
      x_data <- sim$x_data
      scaled_data <- lapply(x_data, function(x) x * 0.1)
      cat("Iteration sim_6_sc01_n100", i, "\n")

      set.seed(1 * i)
      result <- fit_pca_ilr_vs_2(
        scaled_data, 
        sc_factor = 1,
        r = 15,
        max_iter = 80,
        eps = 0.05
      )
      pca_results_list[[i]] <- result
    }
    pca_results_list
  }),
  tar_target(pca_6_sc001_n100, {
    number_simulations <- length(sim_6_sc1_n100)
    pca_results_list <- list(length(sim_6_sc1_n100))
    for  (i in 1:length(sim_6_sc1_n100)) {
      sim <-  sim_6_sc1_n100[[i]]
      x_data <- sim$x_data
      scaled_data <- lapply(x_data, function(x) x * 0.01)
      cat("Iteration sim_6_sc001_n100", i, "\n")

      set.seed(1 * i)
      result <- fit_pca_ilr_vs_2(
        scaled_data, 
        sc_factor = 1,
        r = 15,
        max_iter = 80,
        eps = 0.05
      )
      pca_results_list[[i]] <- result
    }
    pca_results_list
  }),
  tar_target(pca_6_sc0001_n100, {
    number_simulations <- length(sim_6_sc1_n100)
    pca_results_list <- list(length(sim_6_sc1_n100))
    for  (i in 1:length(sim_6_sc1_n100)) {
      sim <-  sim_6_sc1_n100[[i]]
      x_data <- sim$x_data
      scaled_data <- lapply(x_data, function(x) x * 0.001)
      cat("Iteration sim_6_sc0001_n100", i, "\n")

      set.seed(1 * i)
      result <- fit_pca_ilr_vs_2(
        scaled_data, 
        sc_factor = 1,
        r = 15,
        max_iter = 80,
        eps = 0.05
      )
      pca_results_list[[i]] <- result
    }
    pca_results_list
  }),
  tar_target(pca_6_sc1_n50_vs2, {
    number_simulations <- length(sim_6_sc1_n50)
    pca_results_list <- list(length(sim_6_sc1_n50))
    for  (i in 1:length(sim_6_sc1_n50)) {
      sim <-  sim_6_sc1_n50[[i]]
      x_data <- sim$x_data
      cat("Iteration sim_6_sc1_n50", i, "\n")
      set.seed(1 * i)
      result <- fit_pca_ilr_vs_2(
        x_data, 
        sc_factor = 1,
        r = 15,
        max_iter = 80,
        eps = 0.07
      )
      pca_results_list[[i]] <- result
    }
    pca_results_list
  }),
  tar_target(pca_6_sc01_n50_vs2, {
    number_simulations <- length(sim_6_sc01_n50)
    pca_results_list <- list(length(sim_6_sc01_n50))
    for  (i in 1:length(sim_6_sc01_n50)) {
      sim <-  sim_6_sc01_n50[[i]]
      x_data <- sim$x_data
      cat("Iteration sim_6_sc01_n50", i, "\n")
      set.seed(1 * i)
      result <- fit_pca_ilr_vs_2(
        x_data, 
        sc_factor = 1,
        r = 15,
        max_iter = 80,
        eps = 0.07
      )
      pca_results_list[[i]] <- result
    }
    pca_results_list
  }),
  tar_target(pca_6_sc001_n50_vs2, {
    number_simulations <- length(sim_6_sc001_n50)
    pca_results_list <- list(length(sim_6_sc001_n50))
    for  (i in 1:length(sim_6_sc001_n50)) {
      sim <-  sim_6_sc001_n50[[i]]
      x_data <- sim$x_data
      cat("Iteration sim_6_sc001_n50", i, "\n")
      set.seed(1 * i)
      result <- fit_pca_ilr_vs_2(
        x_data, 
        sc_factor = 1,
        r = 15,
        max_iter = 80,
        eps = 0.06
      )
      pca_results_list[[i]] <- result
    }
    pca_results_list
  }),
  tar_target(pca_6_sc01_n100_vs2, {
    number_simulations <- length(sim_6_sc01_n100)
    pca_results_list <- list(length(sim_6_sc01_n100))
    for  (i in 1:length(sim_6_sc01_n100)) {
      sim <-  sim_6_sc01_n100[[i]]
      x_data <- sim$x_data
      cat("Iteration sim_6_sc01_n100", i, "\n")
      set.seed(1 * i)
      result <- fit_pca_ilr_vs_2(
        x_data, 
        sc_factor = 1,
        r = 15,
        max_iter = 80,
        eps = 0.06
      )
      pca_results_list[[i]] <- result
    }
    pca_results_list
  }),
  tar_target(pca_6_sc001_n100_vs2, {
    number_simulations <- length(sim_6_sc001_n100)
    pca_results_list <- list(length(sim_6_sc001_n100))
    for  (i in 1:length(sim_6_sc001_n100)) {
      sim <-  sim_6_sc001_n100[[i]]
      x_data <- sim$x_data
      cat("Iteration sim_6_sc001_n100", i, "\n")
      set.seed(1 * i)
      result <- fit_pca_ilr_vs_2(
        x_data, 
        sc_factor = 1,
        r = 15,
        max_iter = 80,
        eps = 0.05
      )
      pca_results_list[[i]] <- result
    }
    pca_results_list
  }),
  tar_target(pca_6_sc0001_n100_vs2, {
    number_simulations <- length(sim_6_sc0001_n100)
    pca_results_list <- list(length(sim_6_sc0001_n100))
    for  (i in 1:length(sim_6_sc0001_n100)) {
      sim <-  sim_6_sc0001_n100[[i]]
      x_data <- sim$x_data
      cat("Iteration sim_6_sc0001_n100", i, "\n")
      set.seed(1 * i)
      result <- fit_pca_ilr_vs_2(
        x_data, 
        sc_factor = 1,
        r = 15,
        max_iter = 80,
        eps = 0.05
      )
      pca_results_list[[i]] <- result
    }
    pca_results_list
  }),
  ##********* Setting 5**********#
  tar_target(pca_5_0001, {
    number_simulations <- length(sim_5_0001)
    pca_results_list <- list(length(sim_5_0001))
    for  (i in 1:length(sim_5_0001)) {
      sim <-  sim_5_0001[[i]]
      x_data <- sim$x_data
      cat("Iteration sim_5_0001", i, "\n")

      set.seed(1 * i)
      result <- fit_pca_ilr_vs_2(
          x_data, 
          sc_factor = 1,
          max_iter = 80,
          eps = 0.03
      )
      pca_results_list[[i]] <- result
    }
    pca_results_list
  }),
  tar_target(pca_5_001, {
    number_simulations <- length(sim_5_001)
    pca_results_list <- list(length(sim_5_001))
    for  (i in 1:length(sim_5_001)) {
      sim <-  sim_5_001[[i]]
      x_data <- sim$x_data
      cat("Iteration sim_5_001", i, "\n")

      set.seed(1 * i)
      result <- fit_pca_ilr_vs_2(
          x_data, 
          sc_factor = 1,
          max_iter = 80,
          eps = 0.03
      )
      pca_results_list[[i]] <- result
    }
    pca_results_list
  }),
  tar_target(pca_5_01, {
    number_simulations <- length(sim_5_01)
    pca_results_list <- list(length(sim_5_01))
    for  (i in 1:length(sim_5_01)) {
      sim <-  sim_5_01[[i]]
      x_data <- sim$x_data
      cat("Iteration sim_5_01", i, "\n")

      set.seed(1 * i)
      result <- fit_pca_ilr_vs_2(
          x_data, 
          sc_factor = 1,
          max_iter = 80,
          eps = 0.03
      )
      pca_results_list[[i]] <- result
    }
    pca_results_list
  }),
  tar_target(pca_5_05, {
    number_simulations <- length(sim_5_05)
    pca_results_list <- list(length(sim_5_05))
    for  (i in 1:length(sim_5_05)) {
      sim <-  sim_5_05[[i]]
      x_data <- sim$x_data
      cat("Iteration sim_5_05", i, "\n")

      set.seed(1 * i)
      result <- fit_pca_ilr_vs_2(
          x_data, 
          sc_factor = 1,
          max_iter = 80,
          eps = 0.04
      )
      pca_results_list[[i]] <- result
    }
    pca_results_list
  }),
  ##********* Setting 4**********#
  tar_target(pca_smi_nSim_4_30, {
    number_simulations <- length(sim_4_smi_nSim_30)
    pca_results_list <- list(length(sim_4_smi_nSim_30))
    for  (i in 1:length(sim_4_smi_nSim_30)) {
      sim <-  sim_4_smi_nSim_30[[i]]
      x_data <- sim$x_data
      cat("Iteration smi_nSim_4_30", i, "\n")

      set.seed(1 * i)
      result <- fit_pca_ilr_vs_4(
          x_data, 
          sc_factor = 1,
          max_iter = 80,
          eps = 0.02
      )
      pca_results_list[[i]] <- result
    }
    pca_results_list
  }),
  tar_target(pca_smi_nSim_4_60, {
    number_simulations <- length(sim_4_smi_nSim_60)
    pca_results_list <- list(length(sim_4_smi_nSim_60))
    for  (i in 1:length(sim_4_smi_nSim_60)) {
      sim <-  sim_4_smi_nSim_60[[i]]
      x_data <- sim$x_data
      cat("Iteration smi_nSim_4_60", i, "\n")

      set.seed(1 * i)
      result <- fit_pca_ilr_vs_4(
          x_data, 
          sc_factor = 1,
          max_iter = 80,
          eps = 0.02
      )
      pca_results_list[[i]] <- result
    }
    pca_results_list
  }),
  tar_target(pca_smi_nSim_4_120, {
    number_simulations <- length(sim_4_smi_nSim_120)
    pca_results_list <- list(length(sim_4_smi_nSim_120))
    for  (i in 1:length(sim_4_smi_nSim_120)) {
      sim <-  sim_4_smi_nSim_120[[i]]
      x_data <- sim$x_data
      cat("Iteration smi_nSim_3_120", i, "\n")

      set.seed(1 * i)
      result <- fit_pca_ilr_vs_4(
          x_data, 
          sc_factor = 1,
          max_iter = 80,
          eps = 0.03
      )
      pca_results_list[[i]] <- result
    }
    pca_results_list
  }),
  tar_target(pca_smi_nSim_4_400, {
    number_simulations <- length(sim_4_smi_nSim_400)
    pca_results_list <- list(length(sim_4_smi_nSim_400))
    for  (i in 1:length(sim_4_smi_nSim_400)) {
      sim <-  sim_4_smi_nSim_400[[i]]
      x_data <- sim$x_data
      cat("Iteration smi_nSim_4_400", i, "\n")

      set.seed(1 * i)
      result <- fit_pca_ilr_vs_4(
          x_data, 
          sc_factor = 1,
          r = 50,
          max_iter = 80,
          eps = 0.04
      )
      pca_results_list[[i]] <- result
    }
    pca_results_list
  }),
  ##********* Setting 3**********#
  tar_target(pca_smi_nSim_3_30, {
    number_simulations <- length(sim_3_smi_nSim_30)
    pca_results_list <- list(length(sim_3_smi_nSim_30))
    for  (i in 1:length(sim_3_smi_nSim_30)) {
      sim <-  sim_3_smi_nSim_30[[i]]
      x_data <- sim$x_data
      cat("Iteration smi_nSim_3_30", i, "\n")

      set.seed(1 * i)
      result <- fit_pca_ilr_vs_4(
          x_data, 
          sc_factor = 1,
          max_iter = 80,
          eps = 0.02
      )
      pca_results_list[[i]] <- result
    }
    pca_results_list
  }),
  tar_target(pca_smi_nSim_3_20, {
    number_simulations <- length(sim_3_smi_nSim_20)
    pca_results_list <- list(length(sim_3_smi_nSim_20))
    for  (i in 1:length(sim_3_smi_nSim_20)) {
      sim <-  sim_3_smi_nSim_20[[i]]
      x_data <- sim$x_data
      cat("Iteration smi_nSim_3_20", i, "\n")

      set.seed(1 * i)
      result <- fit_pca_ilr_vs_4(
          x_data, 
          sc_factor = 1,
          max_iter = 80,
          eps = 0.035
      )
      pca_results_list[[i]] <- result
    }
    pca_results_list
  }),
  tar_target(pca_smi_nSim_3_40, {
    number_simulations <- length(sim_3_smi_nSim_40)
    pca_results_list <- list(length(sim_3_smi_nSim_40))
    for  (i in 1:length(sim_3_smi_nSim_40)) {
      sim <-  sim_3_smi_nSim_40[[i]]
      x_data <- sim$x_data
      cat("Iteration smi_nSim_3_40", i, "\n")

      set.seed(1 * i)
      result <- fit_pca_ilr_vs_4(
          x_data, 
          sc_factor = 1,
          max_iter = 80,
          eps = 0.025
      )
      pca_results_list[[i]] <- result
    }
    pca_results_list
  }),
  tar_target(pca_smi_nSim_3_80, {
    number_simulations <- length(sim_3_smi_nSim_80)
    pca_results_list <- list(length(sim_3_smi_nSim_80))
    for  (i in 1:length(sim_3_smi_nSim_80)) {
      sim <-  sim_3_smi_nSim_80[[i]]
      x_data <- sim$x_data
      cat("Iteration smi_nSim_3_80", i, "\n")

      set.seed(1 * i)
      result <- fit_pca_ilr_vs_4(
          x_data, 
          sc_factor = 1,
          max_iter = 80,
          eps = 0.02
      )
      pca_results_list[[i]] <- result
    }
    pca_results_list
  }),
  tar_target(pca_smi_nSim_3_160, {
    number_simulations <- length(sim_3_smi_nSim_160)
    pca_results_list <- list(length(sim_3_smi_nSim_160))
    for  (i in 1:length(sim_3_smi_nSim_160)) {
      sim <-  sim_3_smi_nSim_160[[i]]
      x_data <- sim$x_data
      cat("Iteration smi_nSim_3_160", i, "\n")

      set.seed(1 * i)
      result <- fit_pca_ilr_vs_4(
          x_data, 
          sc_factor = 1,
          max_iter = 80,
          eps = 0.025
      )
      pca_results_list[[i]] <- result
    }
    pca_results_list
  }),
  tar_target(pca_smi_nSim_3_60, {
    number_simulations <- length(sim_3_smi_nSim_60)
    pca_results_list <- list(length(sim_3_smi_nSim_60))
    for  (i in 1:length(sim_3_smi_nSim_60)) {
      sim <-  sim_3_smi_nSim_60[[i]]
      x_data <- sim$x_data
      cat("Iteration smi_nSim_3_60", i, "\n")

      set.seed(1 * i)
      result <- fit_pca_ilr_vs_4(
          x_data, 
          sc_factor = 1,
          max_iter = 80,
          eps = 0.02
      )
      pca_results_list[[i]] <- result
    }
    pca_results_list
  }),
  tar_target(pca_smi_nSim_3_120, {
    number_simulations <- length(sim_3_smi_nSim_120)
    pca_results_list <- list(length(sim_3_smi_nSim_120))
    for  (i in 1:length(sim_3_smi_nSim_120)) {
      sim <-  sim_3_smi_nSim_120[[i]]
      x_data <- sim$x_data
      cat("Iteration smi_nSim_3_120", i, "\n")

      set.seed(1 * i)
      result <- fit_pca_ilr_vs_4(
          x_data, 
          sc_factor = 1,
          max_iter = 80,
          eps = 0.02
      )
      pca_results_list[[i]] <- result
    }
    pca_results_list
  }),
  tar_target(pca_smi_nSim_3_400, {
    number_simulations <- length(sim_3_smi_nSim_400)
    pca_results_list <- list(length(sim_3_smi_nSim_400))
    for  (i in 1:length(sim_3_smi_nSim_400)) {
      sim <-  sim_3_smi_nSim_400[[i]]
      x_data <- sim$x_data
      cat("Iteration smi_nSim_3_400", i, "\n")

      set.seed(1 * i)
      result <- fit_pca_ilr_vs_4(
          x_data, 
          sc_factor = 1,
          max_iter = 80,
          eps = 0.04
      )
      pca_results_list[[i]] <- result
    }
    pca_results_list
  }),
  ##********* Setting 2**********#
  tar_target(pca_sim1_ilr_StdPara_smi_nSim_2_60, {
    number_simulations <- length(sim_comp_2_smi_nSim_60)
    pca_results_list <- list(length(sim_comp_2_smi_nSim_60))
    for  (i in 1:length(sim_comp_2_smi_nSim_60)) {
      sim <-  sim_comp_2_smi_nSim_60[[i]]
      x_data <- sim$x_data
      cat("Iteration smi_nSim_2_60", i, "\n")

      set.seed(1 * i)
      result <- fit_pca_ilr_vs_4(
          x_data, 
          sc_factor = 1,
          max_iter = 80,
          eps = 0.02
      )
      pca_results_list[[i]] <- result
    }
    pca_results_list
  }),
  tar_target(pca_sim1_ilr_StdPara_smi_nSim_2_20, {
    number_simulations <- length(sim_comp_2_smi_nSim_20)
    pca_results_list <- list(length(sim_comp_2_smi_nSim_20))
    for  (i in 1:length(sim_comp_2_smi_nSim_20)) {
      sim <-  sim_comp_2_smi_nSim_20[[i]]
      x_data <- sim$x_data
      cat("Iteration smi_nSim_2_20", i, "\n")

      set.seed(1 * i)
      result <- fit_pca_ilr_vs_4(
          x_data, 
          sc_factor = 1,
          max_iter = 80,
          eps = 0.03
      )
      pca_results_list[[i]] <- result
    }
    pca_results_list
  }),
  tar_target(pca_sim1_ilr_StdPara_smi_nSim_2_40, {
    number_simulations <- length(sim_comp_2_smi_nSim_40)
    pca_results_list <- list(length(sim_comp_2_smi_nSim_40))
    for  (i in 1:length(sim_comp_2_smi_nSim_40)) {
      sim <-  sim_comp_2_smi_nSim_40[[i]]
      x_data <- sim$x_data
      cat("Iteration smi_nSim_2_40", i, "\n")

      set.seed(1 * i)
      result <- fit_pca_ilr_vs_4(
          x_data, 
          sc_factor = 1,
          max_iter = 80,
          eps = 0.02
      )
      pca_results_list[[i]] <- result
    }
    pca_results_list
  }),
  # tar_target(pca_sim1_2_80_l06, {
  #   number_simulations <- length(sim_comp_2_smi_nSim_80)
  #   pca_results_list <- list(length(sim_comp_2_smi_nSim_80))
  #   for  (i in 1:length(sim_comp_2_smi_nSim_80)) {
  #     sim <-  sim_comp_2_smi_nSim_80[[i]]
  #     x_data <- sim$x_data
  #     cat("pca_sim1_2_80_l06", i, "\n")

  #     set.seed(1 * i)
  #     result <- fit_pca_ilr_vs_4(
  #         x_data, 
  #         sc_factor = 1,
  #         lambda = 0.6,
  #         max_iter = 80,
  #         eps = 0.02
  #     )
  #     pca_results_list[[i]] <- result
  #   }
  #   pca_results_list
  # }),
  tar_target(pca_sim1_2_20_eb, {
    number_simulations <- length(sim_comp_2_smi_nSim_eb_20)
    pca_results_list <- list(length(sim_comp_2_smi_nSim_eb_20))
    for  (i in 1:length(sim_comp_2_smi_nSim_eb_20)) {
      sim <-  sim_comp_2_smi_nSim_eb_20[[i]]
      x_data <- sim$x_data
      cat("sim_comp_2_smi_nSim_eb_20", i, "\n")

      set.seed(1 * i)
      result <- fit_pca_ilr_vs_4(
          x_data, 
          sc_factor = 1,
          lambda = 1,
          max_iter = 80,
          eps = 0.02
      )
      pca_results_list[[i]] <- result
    }
    pca_results_list
  }),
  tar_target(pca_sim1_2_40_eb, {
    number_simulations <- length(sim_comp_2_smi_nSim_eb_40)
    pca_results_list <- list(length(sim_comp_2_smi_nSim_eb_40))
    for  (i in 1:length(sim_comp_2_smi_nSim_eb_40)) {
      sim <-  sim_comp_2_smi_nSim_eb_40[[i]]
      x_data <- sim$x_data
      cat("sim_comp_2_smi_nSim_eb_40", i, "\n")

      set.seed(1 * i)
      result <- fit_pca_ilr_vs_4(
          x_data, 
          sc_factor = 1,
          lambda = 1,
          max_iter = 80,
          eps = 0.02
      )
      pca_results_list[[i]] <- result
    }
    pca_results_list
  }),
  tar_target(pca_sim1_2_80_eb, {
    number_simulations <- length(sim_comp_2_smi_nSim_eb_80)
    pca_results_list <- list(length(sim_comp_2_smi_nSim_eb_80))
    for  (i in 1:length(sim_comp_2_smi_nSim_eb_80)) {
      sim <-  sim_comp_2_smi_nSim_eb_80[[i]]
      x_data <- sim$x_data
      cat("sim_comp_2_smi_nSim_eb_80", i, "\n")

      set.seed(1 * i)
      result <- fit_pca_ilr_vs_4(
          x_data, 
          sc_factor = 1,
          lambda = 1,
          max_iter = 80,
          eps = 0.02
      )
      pca_results_list[[i]] <- result
    }
    pca_results_list
  }),
  tar_target(pca_sim1_2_160_eb, {
    number_simulations <- length(sim_comp_2_smi_nSim_eb_160)
    pca_results_list <- list(length(sim_comp_2_smi_nSim_eb_160))
    for  (i in 1:length(sim_comp_2_smi_nSim_eb_160)) {
      sim <-  sim_comp_2_smi_nSim_eb_160[[i]]
      x_data <- sim$x_data
      cat("sim_comp_2_smi_nSim_eb_16 0", i, "\n")

      set.seed(1 * i)
      result <- fit_pca_ilr_vs_4(
          x_data, 
          sc_factor = 1,
          lambda = 1,
          max_iter = 80,
          eps = 0.02
      )
      pca_results_list[[i]] <- result
    }
    pca_results_list
  }),
  # tar_target(pca_sim1_2_80_l08, {
  #   number_simulations <- length(sim_comp_2_smi_nSim_80)
  #   pca_results_list <- list(length(sim_comp_2_smi_nSim_80))
  #   for  (i in 1:length(sim_comp_2_smi_nSim_80)) {
  #     sim <-  sim_comp_2_smi_nSim_80[[i]]
  #     x_data <- sim$x_data
  #     cat("Iteration pca_sim1_2_80_l08", i, "\n")

  #     set.seed(1 * i)
  #     result <- fit_pca_ilr_vs_4(
  #         x_data, 
  #         sc_factor = 1,
  #         lambda = 0.8,
  #         max_iter = 80,
  #         eps = 0.02
  #     )
  #     pca_results_list[[i]] <- result
  #   }
  #   pca_results_list
  # }),
  # tar_target(pca_sim1_2_80_l12, {
  #   number_simulations <- length(sim_comp_2_smi_nSim_80)
  #   pca_results_list <- list(length(sim_comp_2_smi_nSim_80))
  #   for  (i in 1:length(sim_comp_2_smi_nSim_80)) {
  #     sim <-  sim_comp_2_smi_nSim_80[[i]]
  #     x_data <- sim$x_data
  #     cat("Iteration smi_nSim_2_80_l12", i, "\n")

  #     set.seed(1 * i)
  #     result <- fit_pca_ilr_vs_4(
  #         x_data, 
  #         sc_factor = 1,
  #         lambda = 1.2,
  #         max_iter = 80,
  #         eps = 0.02
  #     )
  #     pca_results_list[[i]] <- result
  #   }
  #   pca_results_list
  # }),
  # tar_target(pca_sim1_2_80_l14, {
  #   number_simulations <- length(sim_comp_2_smi_nSim_80)
  #   pca_results_list <- list(length(sim_comp_2_smi_nSim_80))
  #   for  (i in 1:length(sim_comp_2_smi_nSim_80)) {
  #     sim <-  sim_comp_2_smi_nSim_80[[i]]
  #     x_data <- sim$x_data
  #     cat("Iteration smi_nSim_2_l14", i, "\n")

  #     set.seed(1 * i)
  #     result <- fit_pca_ilr_vs_4(
  #         x_data, 
  #         sc_factor = 1,
  #         lambda = 1.4,
  #         max_iter = 80,
  #         eps = 0.02
  #     )
  #     pca_results_list[[i]] <- result
  #   }
  #   pca_results_list
  # }),
  tar_target(pca_sim1_ilr_StdPara_smi_nSim_2_80, {
    number_simulations <- length(sim_comp_2_smi_nSim_80)
    pca_results_list <- list(length(sim_comp_2_smi_nSim_80))
    for  (i in 1:length(sim_comp_2_smi_nSim_80)) {
      sim <-  sim_comp_2_smi_nSim_80[[i]]
      x_data <- sim$x_data
      cat("Iteration smi_nSim_2_80", i, "\n")

      set.seed(1 * i)
      result <- fit_pca_ilr_vs_4(
          x_data, 
          sc_factor = 1,
          max_iter = 80,
          eps = 0.02
      )
      pca_results_list[[i]] <- result
    }
    pca_results_list
  }),
  tar_target(pca_sim1_ilr_StdPara_smi_nSim_2_160, {
    number_simulations <- length(sim_comp_2_smi_nSim_160)
    pca_results_list <- list(length(sim_comp_2_smi_nSim_160))
    for  (i in 1:length(sim_comp_2_smi_nSim_160)) {
      sim <-  sim_comp_2_smi_nSim_160[[i]]
      x_data <- sim$x_data
      cat("Iteration smi_nSim_2_160", i, "\n")

      set.seed(1 * i)
      result <- fit_pca_ilr_vs_4(
          x_data, 
          sc_factor = 1,
          r = 20,
          max_iter = 80,
          eps = 0.035
      )
      pca_results_list[[i]] <- result
    }
    pca_results_list
  }),
  tar_target(pca_sim1_ilr_StdPara_smi_nSim_2_120, {
    number_simulations <- length(sim_comp_2_smi_nSim_120)
    pca_results_list <- list(length(sim_comp_2_smi_nSim_120))
    for  (i in 1:length(sim_comp_2_smi_nSim_120)) {
      sim <-  sim_comp_2_smi_nSim_120[[i]]
      x_data <- sim$x_data
      cat("Iteration smi_nSim_2_120", i, "\n")

      set.seed(1 * i)
      result <- fit_pca_ilr_vs_4(
          x_data, 
          sc_factor = 1,
          max_iter = 80,
          eps = 0.02
      )
      pca_results_list[[i]] <- result
    }
    pca_results_list
  }),
  tar_target(pca_sim1_ilr_StdPara_smi_nSim_2_120_300, {
    number_simulations <- length(sim_comp_2_smi_nSim_120_300)
    pca_results_list <- list(length(sim_comp_2_smi_nSim_120_300))
    for  (i in 1:length(sim_comp_2_smi_nSim_120_300)) {
      sim <-  sim_comp_2_smi_nSim_120_300[[i]]
      x_data <- sim$x_data
      cat("Iteration smi_nSim_2_120_300", i, "\n")

      set.seed(1 * i)
      result <- fit_pca_ilr_vs_4(
        x_data, 
        sc_factor = 1,
        max_iter = 80,
        eps = 0.02
      )
      pca_results_list[[i]] <- result
    }
    pca_results_list
  }),
  tar_target(pca_sim1_ilr_StdPara_smi_nSim_2, {
    number_simulations <- length(sim_comp_2_smi_nSim)
    pca_results_list <- list(length(sim_comp_2_smi_nSim))
    for  (i in 1:length(sim_comp_2_smi_nSim)) {
      sim <-  sim_comp_2_smi_nSim[[i]]
      x_data <- sim$x_data
      cat("Iteration smi_nSim_2", i, "\n")

      set.seed(1 * i)
      result <- fit_pca_ilr_vs_4(
          x_data, 
          sc_factor = 1,
          max_iter = 80,
          eps = 0.01
      )
      pca_results_list[[i]] <- result
    }
    pca_results_list
  }),
  tar_target(pca_sim1_ilr_StdPara_smi_nSim_2_60_e003, {
    number_simulations <- length(sim_comp_2_smi_nSim_60)
    pca_results_list <- list(length(sim_comp_2_smi_nSim_60))
    for  (i in 1:length(sim_comp_2_smi_nSim_60)) {
      sim <-  sim_comp_2_smi_nSim_60[[i]]
      x_data <- sim$x_data
      cat("Iteration smi_nSim_2_60", i, "\n")

      set.seed(1 * i)
      result <- fit_pca_ilr_vs_4(
          x_data, 
          sc_factor = 1,
          max_iter = 80,
          eps = 0.03
      )
      pca_results_list[[i]] <- result
    }
    pca_results_list
  }),
  tar_target(pca_sim1_ilr_StdPara_smi_nSim_2_120_e003, {
    number_simulations <- length(sim_comp_2_smi_nSim_120)
    pca_results_list <- list(length(sim_comp_2_smi_nSim_120))
    for  (i in 1:length(sim_comp_2_smi_nSim_120)) {
      sim <-  sim_comp_2_smi_nSim_120[[i]]
      x_data <- sim$x_data
      cat("Iteration smi_nSim_2_120", i, "\n")

      set.seed(1 * i)
      result <- fit_pca_ilr_vs_4(
          x_data, 
          sc_factor = 1,
          max_iter = 80,
          eps = 0.03
      )
      pca_results_list[[i]] <- result
    }
    pca_results_list
  }),
  tar_target(pca_sim1_ilr_StdPara_smi_nSim_2_120_300_e003, {
    number_simulations <- length(sim_comp_2_smi_nSim_120_300)
    pca_results_list <- list(length(sim_comp_2_smi_nSim_120_300))
    for  (i in 1:length(sim_comp_2_smi_nSim_120_300)) {
      sim <-  sim_comp_2_smi_nSim_120_300[[i]]
      x_data <- sim$x_data
      cat("Iteration smi_nSim_2_120_300", i, "\n")

      set.seed(1 * i)
      result <- fit_pca_ilr_vs_4(
        x_data, 
        sc_factor = 1,
        max_iter = 80,
        eps = 0.03
      )
      pca_results_list[[i]] <- result
    }
    pca_results_list
  }),
  tar_target(pca_sim1_ilr_StdPara_smi_nSim_2_e003, {
    number_simulations <- length(sim_comp_2_smi_nSim)
    pca_results_list <- list(length(sim_comp_2_smi_nSim))
    for  (i in 1:length(sim_comp_2_smi_nSim)) {
      sim <-  sim_comp_2_smi_nSim[[i]]
      x_data <- sim$x_data
      cat("Iteration smi_nSim_2", i, "\n")

      set.seed(1 * i)
      result <- fit_pca_ilr_vs_4(
          x_data, 
          sc_factor = 1,
          max_iter = 80,
          eps = 0.03
      )
      pca_results_list[[i]] <- result
    }
    pca_results_list
  }),
  ##********* Setting 1**********#
  tar_target(pca_sim1_ilr_StdPara_smi_nSim_1, {
    number_simulations <- length(sim_comp_1_smi_nSim)
    pca_results_list <- list(length(sim_comp_1_smi_nSim))
    for  (l in 1:length(sim_comp_1_smi_nSim)) {
      sim <-  sim_comp_1_smi_nSim[[l]]
      x_data <- sim$x_data
      cat("Iteration", l, "\n")

      set.seed(1 * l)
      result <- fit_pca_ilr_vs_4(
          x_data, 
          sc_factor = 1,
          max_iter = 80,
          eps = 0.02
      )
      pca_results_list[[l]] <- result
    }
    pca_results_list
  }),
  tar_target(pca_sim1_ilr_StdPara_smi_nSim_1_300, {
    number_simulations <- length(sim_comp_1_smi_nSim_300)
    pca_results_list <- list(length(sim_comp_1_smi_nSim_300))
    for  (l in 1:length(sim_comp_1_smi_nSim_300)) {
      sim <-  sim_comp_1_smi_nSim_300[[l]]
      x_data <- sim$x_data
      cat("Iteration nSim_300", l, "\n")

      set.seed(1 * l)
      result <- fit_pca_ilr_vs_4(
          x_data, 
          sc_factor = 1,
          max_iter = 80,
          eps = 0.02
      )
      pca_results_list[[l]] <- result
    }
    pca_results_list
  }),
  tar_target(pca_sim1_ilr_StdPara_smi_nSim_1_60, {
    number_simulations <- length(sim_comp_1_smi_nSim_60)
    pca_results_list <- list(length(sim_comp_1_smi_nSim_60))
    for  (l in 1:length(sim_comp_1_smi_nSim_60)) {
      sim <-  sim_comp_1_smi_nSim_60[[l]]
      x_data <- sim$x_data
      cat("Iteration StdPara_smi_nSim_1_60", l, "\n")

      set.seed(1 * l)
      result <- fit_pca_ilr_vs_4(
        x_data, 
        sc_factor = 1,
        max_iter = 80,
        eps = 0.03
      )
      pca_results_list[[l]] <- result
    }
    pca_results_list
  }),
  tar_target(pca_sim1_ilr_StdPara_smi_nSim_1_90, {
    number_simulations <- length(sim_comp_1_smi_nSim_90)
    pca_results_list <- list(length(sim_comp_1_smi_nSim_90))
    for  (l in 1:length(sim_comp_1_smi_nSim_90)) {
      sim <-  sim_comp_1_smi_nSim_90[[l]]
      x_data <- sim$x_data
      cat("Iteration smi_nSim_1_90", l, "\n")

      set.seed(1 * l)
      result <- fit_pca_ilr_vs_4(
          x_data, 
          sc_factor = 1,
          max_iter = 80,
          eps = 0.03
      )
      pca_results_list[[l]] <- result
    }
    pca_results_list
  }),
  tar_target(pca_sim1_ilr_StdPara_smi_nSim_1_90_300, {
    number_simulations <- length(sim_comp_1_smi_nSim_90_300)
    pca_results_list <- list(length(sim_comp_1_smi_nSim_90_300))
    for  (l in 1:length(sim_comp_1_smi_nSim_90_300)) {
      sim <-  sim_comp_1_smi_nSim_90_300[[l]]
      x_data <- sim$x_data
      cat("Iteration smi_nSim_1_90_300", l, "\n")

      set.seed(1 * l)
      result <- fit_pca_ilr_vs_4(
          x_data, 
          sc_factor = 1,
          max_iter = 80,
          eps = 0.03
      )
      pca_results_list[[l]] <- result
    }
    pca_results_list
  }),
  tar_target(pca_count_ilr_vs1_1perc, {
    x_data <- data_kl15_comp
    # rescale the data to avoid overflow issues in the log likelihood
    x_data <- x_data * 0.01
    set.seed(12)
    pca_results_ilr_std <-
      fit_pca_ilr_vs_2(x_data,
                       max_iter = 50,
                       r = 10,
                       lambda = 1,
                       eps = 0.01,
                       sc_factor = 1,
                       sum_exp = TRUE)
  }),
  tar_target(pca_count_ilr_vs1_01, {
    x_data <- data_kl15_comp
    x_data <- x_data * 0.1
    set.seed(12)
    pca_results_ilr_std <-
      fit_pca_ilr_vs_2(x_data,
                       max_iter = 50,
                       r = 10,
                       lambda = 1,
                       eps = 0.02,
                       sc_factor = 1,
                       sum_exp = TRUE)
  }),
  tar_target(pca_count_ilr_vs1_0001, {
    x_data <- data_kl15_comp
    x_data <- x_data * 0.001
    set.seed(12)
    pca_results_ilr_std <-
      fit_pca_ilr_vs_2(x_data,
                       max_iter = 50,
                       r = 10,
                       lambda = 1,
                       eps = 0.02,
                       sc_factor = 1,
                       sum_exp = TRUE)
  }),
  tar_target(pca_count_ilr_vs1_sc, {
    x_data <- data_kl15_comp
    set.seed(12)
    pca_results_ilr_std <-
      fit_pca_ilr_vs_4(x_data,
                       max_iter = 10,
                       r = 20,
                       lambda = 1,
                       eps = 0.035,
                       sc_factor = 1, 
                       sum_exp = TRUE)
  }),
  tar_target(pca_count_ilr_vs2_sc, {
    x_data <- data_kl15_comp
    set.seed(12)
    pca_results_ilr_std <-
      fit_pca_ilr_vs_4(x_data,
                       max_iter = 50,
                       r = 20,
                       lambda = 0.8,
                       eps = 0.03,
                       sc_factor = 1, 
                       sum_exp = TRUE)
  }),
  tar_target(pca_count_ilr_vs5_sc01, {
    x_data <- data_kl15_comp
    x_data <- x_data * 0.01
    set.seed(12)
    pca_results_ilr_std <-
      fit_pca_vs_5(x_data,
                   max_iter = 50,
                   r = 10,
                   lambda = 1,
                   eps = 0.03,
                   sc_factor = 1,
                   scores = TRUE,
                   sum_exp = TRUE,
                   fix_sign = TRUE)
 }),
  tar_target(pca_count_ilr_vs5, {
    x_data <- data_kl15_comp
    set.seed(12)
    pca_results_ilr_std <-
      fit_pca_vs_5(x_data,
                   max_iter = 50,
                   r = 10,
                   lambda = 1,
                   eps = 0.03,
                   sc_factor = 1,
                   scores = TRUE,
                   sum_exp = TRUE,
                   fix_sign = TRUE)
 }),
  tar_target(pca_count_ilr_vs5_e06, {
    x_data <- data_kl15_comp
    set.seed(12)
    pca_results_ilr_std <-
      fit_pca_vs_5(x_data,
                   max_iter = 50,
                   r = 10,
                   lambda = 1,
                   eps = 0.06,
                   sc_factor = 1,
                   scores = TRUE,
                   sum_exp = TRUE,
                   fix_sign = TRUE)
 }),
  tar_target(pca_count_ilr_vs5_sc01_e06, {
    x_data <- data_kl15_comp
    x_data <- x_data * 0.01
    set.seed(12)
    pca_results_ilr_std <-
      fit_pca_vs_5(x_data,
                   max_iter = 50,
                   r = 10,
                   lambda = 1,
                   eps = 0.06,
                   sc_factor = 1,
                   scores = TRUE,
                   sum_exp = TRUE,
                   fix_sign = TRUE)
 }),
  tar_target(pca_count_ilr_vs5_sc01_noSign, {
    x_data <- data_kl15_comp
    x_data <- x_data * 0.01
    set.seed(12)
    pca_results_ilr_std <-
      fit_pca_vs_5(x_data,
                   max_iter = 50,
                   r = 10,
                   lambda = 1,
                   eps = 0.03,
                   sc_factor = 1,
                   scores = TRUE,
                   sum_exp = TRUE,
                   fix_sign = FALSE)
 }),
  tar_target(pca_count_ilr_vs5_sc01_r30, {
    x_data <- data_kl15_comp
    x_data <- x_data * 0.01
    set.seed(12)
    pca_results_ilr_std <-
      fit_pca_vs_5(x_data,
                   max_iter = 50,
                   r = 30,
                   lambda = 1,
                   eps = 0.03,
                   sc_factor = 1,
                   scores = TRUE,
                   sum_exp = TRUE,
                   fix_sign = TRUE)
 }),
  tar_target(pca_count_ilr_vs6_sc01, {
    x_data <- data_kl15_comp
    x_data <- x_data * 0.01
    set.seed(12)
    pca_results_ilr_std <-
      fit_pca_vs_6(x_data,
                   max_iter = 50,
                   r = 10,
                   lambda = 1,
                   eps = 0.03,
                   sc_factor = 1,
                   scores = TRUE,
                   sum_exp = TRUE,
                   fix_sign = TRUE)
 }),
  tar_target(pca_count_ilr_vs6_acomp, {
    x_data <- data_kl15_comp
    x_data <- acomp(x_data)
    set.seed(12)
    pca_results_ilr_std <-
      fit_pca_vs_6(x_data,
                   max_iter = 50,
                   r = 10,
                   lambda = 1,
                   eps = 0.03,
                   sc_factor = 1,
                   scores = TRUE,
                   sum_exp = TRUE,
                   fix_sign = TRUE)
 }),
  tar_target(bbootstrap_01_1, {
    x_data <- data_kl15_comp
    results_mcem <- pca_count_ilr_vs1_1perc
    reference_pca <- results_mcem$pca    
    set.seed(12)
    results <-
      bbootstrap_pca_ts_par(x_data,
                        block_length = 50,
                        replicates = 100,
                        reference_pca = reference_pca,
                        scale = 0.01,
                        eps = 0.03,
                        workers = 10)
  }),
  tar_target(bbootstrap_01_2, {
    x_data <- data_kl15_comp
    results_mcem <- pca_count_ilr_vs1_1perc
    reference_pca <- results_mcem$pca    
    set.seed(12)
    results <-
      bbootstrap_pca_ts_par_vs2(x_data,
                        block_length = 50,
                        replicates = 100,
                        reference_pca = reference_pca,
                        scale = 0.01,
                        eps = 0.03,
                        workers = 10)
  })
  # tar_target(bbootstrap_01_2, {
  #   x_data <- data_kl15_comp
  #   results_mcem <- pca_count_ilr_vs1_01
  #   reference_pca <- results_mcem$pca
  #   set.seed(123)
  #   results <-
  #     bbootstrap_pca_ts_par(x_data,
  #                       block_length = 50,
  #                       replicates = 900,
  #                       reference_pca = reference_pca,
  #                       scale = 0.01,
  #                       eps = 0.03,
  #                       workers = 10)
  # })
)
