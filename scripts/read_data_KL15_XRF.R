library(dplyr)
library(tidyr)

# environment
dir <- paste0(getwd(), "/data/Africa_NE_200/data/")

# TODO: chain all commands with pipelines (readr)

# read the sedimant data and the age model
data_kl15_xrf <- read.table(paste0(dir,'data_KL15_XRF.txt'), header = TRUE, sep = "\t")
data_kl15_agem <- read.table(paste0(dir,'data_KL15-2_smooth_spline_ages.txt'), header = TRUE, sep = "\t")
data_kl15_qf <- read.table(paste0(dir,'data_KL15_qf.txt'), header = TRUE, sep = "\t")

# rename the columns
data_kl15_xrf <- rename(data_kl15_xrf, depth = User_ID)
data_kl15_agem <- rename(data_kl15_agem, age = best)

# TODO: try to implement:

# Perform an outer join
T <- full_join(data_kl15_agem, data_kl15_xrf, by = "common_column") # replace "common_column" with the column you want to join by

# Get the variable names
data_kl15_string <- colnames(T)

# Convert the data frame to a matrix
data_kl15 <- as.matrix(T)

# Remove rows where the first column's value is in the first column of data_kl15_qf
data_kl15 <- data_kl15[!data_kl15[,1] %in% data_kl15_qf[,1],]

# Remove rows where the sixth column's value is NA
data_kl15 <- data_kl15[!is.na(data_kl15[,6]),]

# Divide the fourth column by 1000
data_kl15[,4] <- data_kl15[,4] / 1000

# Create a sequence from agemodelmin to agemodelmax with a step size of agemodelres
data_kl15_age <- data.frame(matrix(ncol = 29, nrow = length(seq(agemodelmin, agemodelmax, by = agemodelres))))
data_kl15_age[,1] <- seq(agemodelmin, agemodelmax, by = agemodelres)

# Interpolate data to evenly spaced time axis
for (i in 2:29) {
  data_kl15_age[,i] <- approx(data_kl15[,4], data_kl15[,i+5], xout = data_kl15_age[,1])$y
}