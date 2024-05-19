dir <- paste0(getwd(), "/data/Africa_NE_200/data/")

# ODP 967
data_odp_967_22 <- read.table(paste0(dir,'data_odp_967_22_43247_2021_339_MOESM2_ESM_XRF_Ti_Al.txt'), header = TRUE, sep = "\t")
# get the age and Ba_ppm attributes
data_odp_967_22 <- data_odp_967_22[,c(8, 19)]
# filter out data with age less than 200k
data_odp_967_22 <- data_odp_967_22[data_odp_967_22[, 1] < 200,]

# ODP 721/722
data_odp721_722_terr <- read.table(paste0(dir,'data_odp_721_722_terr.txt'), header = TRUE, sep = "\t")
data_odp721_722_terr <- data_odp721_722_terr[,c(2, 3)]
data_odp721_722_terr <- data_odp721_722_terr[data_odp721_722_terr[, 1] < 200,]

# ODP 709
data_odp_709 <- read.table(paste0(dir,'ODP709_Ca_K_Ti_Fe_ratio_600kyr.txt'), header = TRUE, sep = "\t")
data_odp_709 <- data_odp_709[data_odp_709[, 1] < 200,]

# ICPD Chew Bahir
data_icdp_chb <- read.table(paste0(dir,'data_icdp_chb14_2_cb01_cb03_age_rrmarch2021_mht560.txt'), header = TRUE, sep = "\t")
data_icdp_chb <- data_icdp_chb[data_icdp_chb[, 1] < 200,]

# KL 09
data_kl09 <- read.table(paste0(dir,'data_KL09_41467_2014_BFncomms6076_MOESM1344_ESM.txt'), header = TRUE, sep = "\t")
data_kl09 <- data_kl09[,c(1, 2)]
data_kl09 <- data_kl09[data_kl09[, 1] < 200,]

# KL 11
data_kl11 <- read.table(paste0(dir,'data_KL11_geochem.txt'), header = TRUE, sep = "\t")
data_kl11 <- data_kl11[,c(2, 3)]
data_kl11 <- data_kl11[data_kl11[, 1] < 200,]

# KL 15
data_kl15 <- read.table(paste0(dir,'data_KL15_XRF_tuned.txt'), header = TRUE, sep = "\t")
# ToDo: this dataset is missing but is probably calculated from another dataset
datakl15agem <- read.table(paste0(dir,'data_KL15-2_smooth_spline_ages.txt'), header = TRUE, sep = "\t")
names(data_kl15)[names(data_kl15) == "NEW_depth_cm"] <- "depth"
names(datakl15agem)[names(datakl15agem) == "best"] <- "age"
T <- merge(datakl15agem, data_kl15, by = "row.names", all = TRUE)
data_kl15 <- T[,c(4, 8)]
data_kl15[,1] <- data_kl15[,1] / 1000
data_kl15 <- data_kl15[data_kl15[, 1] < 200,]
data_kl15 <- data_kl15[!is.na(data_kl15[,2]),]

# Lake Tana
data_lake_tana <- read.table(paste0(dir,'data_lake_tana'), header = TRUE, sep = "\t")
data_lake_tana <- data_lake_tana[,c(2, 4)]
data_lake_tana[,1] <- data_lake_tana[,1] / 1000
data_lake_tana <- data_lake_tana[data_lake_tana[, 1] < 200,]