process_pneumo_df<-function(in.df) {
  in.df %<>%
    dplyr::select(DS, Serotype, carriage, disease, n.swab, N, time.int) %>%
    dplyr::rename(study = DS) %>%
    dplyr::rename(carriage_samples = n.swab) %>%
    dplyr::rename(surveillance_population = N) %>%
    dplyr::rename(time_interval = time.int) %>%
    dplyr::rename(categorisation = Serotype)
}

# Load raw data from .csv file
S_pneumoniae_infant_serotype <- read.csv("data-raw/S_pneumoniae_infant_serotype.csv")
S_pneumoniae_adult_serotype <- read.csv("data-raw/S_pneumoniae_adult_serotype.csv")

# Save the cleaned data in the required R package location
usethis::use_data(S_pneumoniae_infant_serotype, overwrite = TRUE)
usethis::use_data(S_pneumoniae_adult_serotype, overwrite = TRUE)
