require(tidyverse)
require(magrittr)
require(progressionEstimation)
require(testthat)

# Generate test data subset
make_test_df <- function() {
  test.df <-
    S_pneumoniae_infant_strain %>%
    dplyr::filter(type %in% c("19A","14")) %>%
    dplyr::filter(study == "South.Africa.post.PCV7" | study == "South.Africa.post.PCV13") %>%
    dplyr::mutate(study = factor(study)) %>%
    dplyr::mutate(strain = factor(strain))
  return(test.df)
}

# Test input data processing
testthat::test_that("Input data can be processed",{
  test.df <- make_test_df()
  test.in <- progressionEstimation::process_input_data(test.df,
                                                       use_strain = TRUE)
  testthat::expect_equal(test.in$j_max,2)
})

# Test simple model fit
testthat::test_that("Type-specific model can still be fitted",{
  test.df <- make_test_df()
  test.in <- progressionEstimation::process_input_data(test.df,
                                                       use_strain = TRUE)
  type_specific.out <- progressionEstimation::fit_progression_rate_model(
    test.in,
    type_specific = TRUE,
    location_adjustment = FALSE,
    strain_as_primary_type = FALSE,
    strain_as_secondary_type = FALSE,
    stat_model = "poisson",
    num_iter = 10000,
    num_chains = 1
  )
  testthat::expect_match(type_specific.out@model_name,"type_specific_poisson")
})

# Test simple model fit to strains
testthat::test_that("Type-specific model can be fitted to strains",{
  test.df <- make_test_df()
  test.in <- progressionEstimation::process_input_data(test.df,
                                                       type = "strain",
                                                       use_strain = FALSE)
  testthat::expect_gt(test.in$j_max,2)
  type_specific.out <- progressionEstimation::fit_progression_rate_model(
    test.in,
    type_specific = TRUE,
    location_adjustment = FALSE,
    strain_as_primary_type = FALSE,
    strain_as_secondary_type = FALSE,
    stat_model = "poisson",
    num_iter = 10000,
    num_chains = 1
  )
  testthat::expect_match(type_specific.out@model_name,"type_specific_poisson")
})

# Test simple model fit to type-strain combinations
testthat::test_that("Type-specific model can be fitted to strain-type combinations",{
  test.df <- make_test_df()
  test.in <- progressionEstimation::process_input_data(test.df,
                                                       combine_strain = TRUE)
  testthat::expect_gt(test.in$j_max,2)
  type_specific.out <- progressionEstimation::fit_progression_rate_model(
    test.in,
    type_specific = TRUE,
    location_adjustment = FALSE,
    strain_as_primary_type = FALSE,
    strain_as_secondary_type = FALSE,
    stat_model = "poisson",
    num_iter = 10000,
    num_chains = 1
  )
  testthat::expect_match(type_specific.out@model_name,"type_specific_poisson")
})

# Test model fit to to types, strains as secondary
testthat::test_that("Combined model can be fitted to type and strain",{
  test.df <- make_test_df()
  test.in <- progressionEstimation::process_input_data(test.df,
                                                       use_strain = TRUE)
  strain_secondary.out <- progressionEstimation::fit_progression_rate_model(
    test.in,
    type_specific = TRUE,
    location_adjustment = FALSE,
    strain_as_primary_type = FALSE,
    strain_as_secondary_type = TRUE,
    stat_model = "poisson",
    num_iter = 10000,
    num_chains = 1
  )
  testthat::expect_match(strain_secondary.out@model_name,"type_specific_type_modified_by_strain_poisson")
})

# Test model fit to to types, strains as secondary with negbin model
testthat::test_that("Combined model can be fitted to type and strain with negbin",{
  test.df <- make_test_df()
  test.in <- progressionEstimation::process_input_data(test.df,
                                                       use_strain = TRUE)
  strain_secondary.out <- progressionEstimation::fit_progression_rate_model(
    test.in,
    type_specific = TRUE,
    location_adjustment = FALSE,
    strain_as_primary_type = FALSE,
    strain_as_secondary_type = TRUE,
    stat_model = "negbin",
    num_iter = 10000,
    num_chains = 1
  )
  testthat::expect_match(strain_secondary.out@model_name,"type_specific_type_modified_by_strain_negbin")
})

# Test model fit to to types, strains as primary
testthat::test_that("Combined model can be fitted to strain and type",{
  test.df <- make_test_df()
  test.in <- progressionEstimation::process_input_data(test.df,
                                                       use_strain = TRUE)
  strain_primary.out <- progressionEstimation::fit_progression_rate_model(
    test.in,
    type_specific = TRUE,
    location_adjustment = FALSE,
    strain_as_primary_type = TRUE,
    strain_as_secondary_type = FALSE,
    stat_model = "poisson",
    num_iter = 10000,
    num_chains = 1
  )
  testthat::expect_match(strain_primary.out@model_name,"type_specific_strain_modified_by_type_poisson")
})

# Test model fit to to types, strains as primary with negbin model
testthat::test_that("Combined model can be fitted to strain and type with negbin",{
  test.df <- make_test_df()
  test.in <- progressionEstimation::process_input_data(test.df,
                                                       use_strain = TRUE)
  strain_primary.out <- progressionEstimation::fit_progression_rate_model(
    test.in,
    type_specific = TRUE,
    location_adjustment = FALSE,
    strain_as_primary_type = TRUE,
    strain_as_secondary_type = FALSE,
    stat_model = "negbin",
    num_iter = 10000,
    num_chains = 1
  )
  testthat::expect_match(strain_primary.out@model_name,"type_specific_strain_modified_by_type_negbin")
})

# Test model fit to to types, strains as secondary, location adjustment
testthat::test_that("Combined model can be fitted to type and strain with location adjustment",{
  test.df <- make_test_df()
  test.in <- progressionEstimation::process_input_data(test.df,
                                                       use_strain = TRUE)
  strain_secondary.out <- progressionEstimation::fit_progression_rate_model(
    test.in,
    type_specific = TRUE,
    location_adjustment = TRUE,
    strain_as_primary_type = FALSE,
    strain_as_secondary_type = TRUE,
    stat_model = "poisson",
    num_iter = 10000,
    num_chains = 1
  )
  testthat::expect_match(strain_secondary.out@model_name,"adjusted_type_specific_type_modified_by_strain_poisson")
})

# Test model fit to to types, strains as secondary with negbin model, location adjustment
testthat::test_that("Combined model can be fitted to type and strain with location adjustment and negbin",{
  test.df <- make_test_df()
  test.in <- progressionEstimation::process_input_data(test.df,
                                                       use_strain = TRUE)
  strain_secondary.out <- progressionEstimation::fit_progression_rate_model(
    test.in,
    type_specific = TRUE,
    location_adjustment = TRUE,
    strain_as_primary_type = FALSE,
    strain_as_secondary_type = TRUE,
    stat_model = "negbin",
    num_iter = 10000,
    num_chains = 1
  )
  testthat::expect_match(strain_secondary.out@model_name,"adjusted_type_specific_type_modified_by_strain_negbin")
})

# Test model fit to to types, strains as primary, location adjustment
testthat::test_that("Combined model can be fitted to strain and type with location adjustment",{
  test.df <- make_test_df()
  test.in <- progressionEstimation::process_input_data(test.df,
                                                       use_strain = TRUE)
  strain_primary.out <- progressionEstimation::fit_progression_rate_model(
    test.in,
    type_specific = TRUE,
    location_adjustment = TRUE,
    strain_as_primary_type = TRUE,
    strain_as_secondary_type = FALSE,
    stat_model = "poisson",
    num_iter = 10000,
    num_chains = 1
  )
  testthat::expect_match(strain_primary.out@model_name,"adjusted_type_specific_strain_modified_by_type_poisson")
})

# Test model fit to to types, strains as primary with negbin model, location adjustment
testthat::test_that("Combined model can be fitted to strain and type with location adjustment and negbin",{
  test.df <- make_test_df()
  test.in <- progressionEstimation::process_input_data(test.df,
                                                       use_strain = TRUE)
  strain_primary.out <- progressionEstimation::fit_progression_rate_model(
    test.in,
    type_specific = TRUE,
    location_adjustment = TRUE,
    strain_as_primary_type = TRUE,
    strain_as_secondary_type = FALSE,
    stat_model = "negbin",
    num_iter = 10000,
    num_chains = 1
  )
  testthat::expect_match(strain_primary.out@model_name,"adjusted_type_specific_strain_modified_by_type_negbin")
})

# Test model fit to to combined type and strain
testthat::test_that("Combined model can be fitted to type and strain",{
  test.df <- make_test_df()
  test.in <- progressionEstimation::process_input_data(test.df,
                                                       combine_strain = TRUE)
  type_specific.out <- progressionEstimation::fit_progression_rate_model(
    test.in,
    type_specific = TRUE,
    location_adjustment = FALSE,
    strain_as_primary_type = FALSE,
    strain_as_secondary_type = FALSE,
    stat_model = "poisson",
    num_iter = 10000,
    num_chains = 1
  )
  combined_strain.processed.out <-
    progressionEstimation::process_progression_rate_model_output(model_output = type_specific.out,
                                                                 combined_strain_type = TRUE,
                                                                 input_df = test.df
    )

  strain_combined_prediction_plot <-
    progressionEstimation::plot_case_carrier_predictions(combined_strain.processed.out,label_col = "combined")
  expect_error(print(strain_combined_prediction_plot), NA)
  strain_combined_type_rates_plot <-
    progressionEstimation::plot_progression_rates(combined_strain.processed.out,type = "combined")
  expect_error(print(strain_combined_type_rates_plot), NA)
})

# Test model outputs can be processed and plotted - strain secondary
testthat::test_that("Model fit to type and strain can be processed",{
  test.df <- make_test_df()
  test.in <- progressionEstimation::process_input_data(test.df,
                                                       use_strain = TRUE)
  strain_secondary.out <- progressionEstimation::fit_progression_rate_model(
    test.in,
    type_specific = TRUE,
    location_adjustment = FALSE,
    strain_as_primary_type = FALSE,
    strain_as_secondary_type = TRUE,
    stat_model = "poisson",
    num_iter = 10000,
    num_chains = 1
  )
  testthat::expect_match(strain_secondary.out@model_name,"type_specific_type_modified_by_strain_poisson")
  strain_secondary.processed.out <-
    progressionEstimation::process_progression_rate_model_output(model_output = strain_secondary.out,
                                                               input_df = test.df,
                                                               strain_as_secondary_type = TRUE
                                                              )
  testthat::expect_gt(mean(strain_secondary.processed.out$surveillance_population),mean(strain_secondary.processed.out$nu))
  strain_secondary_prediction_plot <-
    progressionEstimation::plot_case_carrier_predictions(strain_secondary.processed.out)
  expect_error(print(strain_secondary_prediction_plot), NA)
  strain_secondary_type_rates_plot <-
    progressionEstimation::plot_progression_rates(strain_secondary.processed.out)
  expect_error(print(strain_secondary_type_rates_plot), NA)
  strain_secondary_strain_rates_plot <-
    progressionEstimation::plot_progression_rates(strain_secondary.processed.out,type = "strain")
  expect_error(print(strain_secondary_strain_rates_plot), NA)
})

# Test model outputs can be processed and plotted - strain primary
testthat::test_that("Model fit to strain and type can be processed",{
  test.df <- make_test_df()
  test.in <- progressionEstimation::process_input_data(test.df,
                                                       use_strain = TRUE)
  strain_primary.out <- progressionEstimation::fit_progression_rate_model(
    test.in,
    type_specific = TRUE,
    location_adjustment = FALSE,
    strain_as_primary_type = TRUE,
    strain_as_secondary_type = FALSE,
    stat_model = "poisson",
    num_iter = 10000,
    num_chains = 1
  )
  testthat::expect_match(strain_primary.out@model_name,"type_specific_strain_modified_by_type_poisson")
  strain_primary.processed.out <-
    progressionEstimation::process_progression_rate_model_output(model_output = strain_primary.out,
                                                                 input_df = test.df,
                                                                 strain_as_primary_type = TRUE
    )
  testthat::expect_gt(mean(strain_primary.processed.out$surveillance_population),mean(strain_primary.processed.out$nu))
  strain_primary_prediction_plot <-
    progressionEstimation::plot_case_carrier_predictions(strain_primary.processed.out)
  expect_error(print(strain_primary_prediction_plot), NA)
  strain_primary_type_rates_plot <-
    progressionEstimation::plot_progression_rates(strain_primary.processed.out)
  expect_error(print(strain_primary_type_rates_plot), NA)
  strain_primary_strain_rates_plot <-
    progressionEstimation::plot_progression_rates(strain_primary.processed.out,type = "strain")
  expect_error(print(strain_primary_strain_rates_plot), NA)
})
