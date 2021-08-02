require(tidyverse)

process_pneumo_df <- function(fn) {
  input_df <-
    read.csv(fn) %>%
      dplyr::select(study, type, carriage_samples, surveillance_population, time_interval, carriage, disease)
  return(input_df)
}

make_factors <- function(df) {
  df %<>%
    dplyr::mutate(study = factor(study)) %>%
    dplyr::mutate(type = factor(type))
  return(df)
}

process_strains <- function(df) {
  df %<>%
    dplyr::mutate(strain =
                    dplyr::case_when(
                      !grepl("\\D", as.character(strain)) ~ paste0("GPSC", strain),
                      TRUE ~ strain
                    )
    ) %>%
    dplyr::mutate(strain = factor(strain))
  return(df)
}

# Load raw data from .csv file
S_pneumoniae_infant_serotype <- process_pneumo_df("data-raw/S_pneumoniae_infant_serotype.csv") %>%
                                  make_factors()
S_pneumoniae_adult_serotype <- process_pneumo_df("data-raw/S_pneumoniae_adult_serotype.csv") %>%
                                  make_factors()
S_pneumoniae_infant_strain <- read.csv("data-raw/S_pneumoniae_infant_strain.csv") %>%
                                  make_factors() %>%
                                  process_strains()
S_pneumoniae_mixed_strain <- read.csv("data-raw/S_pneumoniae_mixed_strain.csv") %>%
                                  make_factors() %>%
                                  process_strains()

# Save the cleaned data in the required R package location
usethis::use_data(S_pneumoniae_infant_serotype, overwrite = TRUE)
usethis::use_data(S_pneumoniae_adult_serotype, overwrite = TRUE)
usethis::use_data(S_pneumoniae_infant_strain, overwrite = TRUE)
usethis::use_data(S_pneumoniae_mixed_strain, overwrite = TRUE)

