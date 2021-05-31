# Load raw data from .csv file
S_pneumoniae_infant_serotype <- read.csv("data-raw/S_pneumoniae_infant_serotype.csv")
S_pneumoniae_adult_serotype <- read.csv("data-raw/S_pneumoniae_adult_serotype.csv")

# Save the cleaned data in the required R package location
usethis::use_data(S_pneumoniae_infant_serotype)
usethis::use_data(S_pneumoniae_adult_serotype)
