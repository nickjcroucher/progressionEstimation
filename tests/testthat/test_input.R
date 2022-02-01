require(tidyverse)
require(magrittr)
require(progressionEstimation)
require(testthat)

testthat::test_that("Input XLSX can be processed",{
  input_xlsx <-
    progressionEstimation::process_input_xlsx("progression_estimation_input_test.xlsx")
  testthat::expect_equal(nrow(input_xlsx),22)
  is_valid <- validate_progression_estimation_dataset(input_xlsx)
  testthat::expect_equal(is_valid,NULL)
})

testthat::test_that("Input XLSX can be processed and appended to data",{
  input_xlsx <-
    progressionEstimation::process_input_xlsx("progression_estimation_input_test.xlsx")
  combined_data <-
    progressionEstimation::combine_with_existing_datasets(input_xlsx,S_pneumoniae_infant_serotype)
  testthat::expect_gt(nrow(combined_data),nrow(input_xlsx))
})
