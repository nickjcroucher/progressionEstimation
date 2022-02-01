require(tidyverse)
require(magrittr)
require(progressionEstimation)
require(testthat)

testthat::test_that("Input XLSX can be processed",{
  input_xlsx <-
    progressionEstimation::process_input_xlsx("progression_estimation_input_test.xlsx")
  testthat::expect_equal(nrow(input_xlsx),22)
})

testthat::test_that("Input XLSX can be processed and appended to data",{
  input_xlsx <-
    progressionEstimation::process_input_xlsx("progression_estimation_input_test.xlsx")
  combined_data <-
    progressionEstimation::combine_with_existing_datasets(input_xlsx,S_pneumoniae_infant_serotype %>% dplyr::filter(grepl("^South.Africa",study)))
  testthat::expect_gt(nrow(combined_data),nrow(input_xlsx))
})
