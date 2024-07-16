library(dplyr)
library(tidyr)
library(compositions)

read_data_kl15_xrf <- function(data_kl15_xrf, data_kl15_agem) {
  data_kl15 <- full_join(data_kl15_agem, data_kl15_xrf, by = "depth") 
  col_names <- colnames(data_kl15[c(4,6:33)])
  missings_depth <- data_kl15[is.na(data_kl15[,6]), 1]

  # Remove rows where the sixth column's value is NA => removes 49 NA's
  data_kl15 <- data_kl15[!is.na(data_kl15[,6]),]

  # Divide the age column by 1000
  data_kl15[,4] <- data_kl15[,4] / 1000


  # Interpolating data to evenly spaced time axis.
  min <- min(data_kl15[,"age"])
  max <- max(data_kl15[,"age"])
  # TODO: check Resolution = (Max-Min)/Length => Resolution of around 260 years
  length <- nrow(data_kl15)
  resolution <- (max - min) / length

  # Create a sequence from agemodelmin to agemodelmax with a step size of agemodelres
  data_kl15_itpol <- data.frame(matrix(ncol = 29, nrow = length(seq(min, max, by = resolution))))
  data_kl15_itpol[,1] <- seq(min, max, by = resolution)


  # Interpolate data to evenly spaced time axis, that is fit the row values to the new time axis
  for (i in 2:28) {
    data_kl15_itpol[,i] <- approx(data_kl15[,4], data_kl15[,i+5], xout = data_kl15_itpol[,1])$y
  }

  # get col names from data_kl15 and name data_kl15_itpol respectivaly
  colnames(data_kl15_itpol) <- col_names

  # define a colum "aggregate" that is the sum of all elements in a row except of the first column
  data_kl15_itpol$aggregate <- rowSums(data_kl15_itpol[2:28])

  # Select the variables of interest
  data_select <- data_kl15[, c("depth", "age", "Br_Area", "Rb_Area", "Sr_Area", "Zr_Area", "Ru_Area", "Mg_Area",
  "Al_Area", "Si_Area", "S_Area", "K_Area", "Ca_Area", "Ti_Area", "Fe_Area")]

  # for compositional data, we are interested in the sum of all components
  data_select$aggregate <- rowSums(data_select[3:ncol(data_select)])

  # get the compositional elements
  data_comp <- data_select[, 3:(ncol(data_select)-1)]

  # compute centered log-ratios
  data_clr <- clr(data_comp)

  results <- list(
    data_kl15_itpol = data_kl15_itpol,
    data_kl15 = data_select,
    data_comp = data_comp,
    missings_depth = missings_depth,
    data_clr = data_clr)

  return(results)
}
