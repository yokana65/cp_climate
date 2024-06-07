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

# Perform an outer join to merge the age model with the XRF data by "depth"
T <- full_join(data_kl15_agem, data_kl15_xrf, by = "depth") 
col_names <- colnames(T[c(4,6:33)])
missings_depth <- T[is.na(T[,6]), 1]

# remove QF marked measurements:
# data_kl15 <- data_kl15[!data_kl15[,1] %in% data_kl15_qf[,1],]
# TODO: ask markus: it does not remove anything; whats the purpose of this action?

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
