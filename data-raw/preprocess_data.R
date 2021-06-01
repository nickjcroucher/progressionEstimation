require(tidyverse)

process_pneumo_df <- function(fn) {
  input_df <-
    read.csv(fn) %>%
    dplyr::mutate(study = factor(DS)) %>%
    dplyr::mutate(categorisation = factor(Serotype)) %>%
    dplyr::rename(carriage_samples = n.swab) %>%
    dplyr::rename(surveillance_population = N) %>%
    dplyr::rename(time_interval = time.int) %>%
    dplyr::select(study, categorisation, carriage_samples, surveillance_population, time_interval)
  return(input_df)
}

# Load raw data from .csv file
S_pneumoniae_infant_serotype <- process_pneumo_df("data-raw/S_pneumoniae_infant_serotype.csv")
S_pneumoniae_adult_serotype <- read.csv("data-raw/S_pneumoniae_adult_serotype.csv")

# Save the cleaned data in the required R package location
usethis::use_data(S_pneumoniae_infant_serotype, overwrite = TRUE)
usethis::use_data(S_pneumoniae_adult_serotype, overwrite = TRUE)
