library(targets)
# This is an example _targets.R file. Every
# {targets} pipeline needs one.
# Use tar_script() to create _targets.R and tar_edit()
# to open it again for editing.
# Then, run tar_make() to run the pipeline
# and tar_read(data_summary) to view the results.

# Define custom functions and other global objects.
# This is where you write source(\"R/functions.R\")
# if you keep your functions in external scripts.
# source("scripts/read_data_200kyr_all.r")

# Set target-specific options such as packages:
# tar_option_set(packages = "utils") # nolint

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
  })
 )
