require(tidyverse)
require(magrittr)
require(progressionEstimation)
require(testthat)

# Generate test data subset
make_test_df <- function() {
  test.df <-
    S_pneumoniae_infant_serotype %>%
      dplyr::filter(type %in% c("19A","14")) %>%
      dplyr::filter(grepl("pre.PCV$",study)) %>%
      dplyr::mutate(type = factor(type)) %>%
      dplyr::mutate(study = factor(study))
  return(test.df)
}

# Test input data processing
testthat::test_that("Input data can be processed",{
  test.df <- make_test_df()
  test.in <- progressionEstimation::process_input_data(test.df)
  testthat::expect_equal(test.in$j_max,2)
})

# Test simple model fit
testthat::test_that("Null model can be fitted",{
  test.df <- make_test_df()
  test.in <- progressionEstimation::process_input_data(test.df)
  null.out <- progressionEstimation::fit_progression_rate_model(
    test.in,
    type_specific = FALSE,
    location_adjustment = FALSE,
    stat_model = "poisson",
    num_iter = 10000,
    num_chains = 1
  )
  testthat::expect_match(null.out@model_name,"null_poisson")
})

# Test negative binomial model fit
testthat::test_that("Negative binomial model can be fitted",{
  test.df <- make_test_df()
  test.in <- progressionEstimation::process_input_data(test.df)
  null.out <- progressionEstimation::fit_progression_rate_model(
    test.in,
    type_specific = FALSE,
    location_adjustment = FALSE,
    stat_model = "negbin",
    num_iter = 10000,
    num_chains = 1
  )
  testthat::expect_match(null.out@model_name,"null_negbin")
})

# Test type-specific model fit
testthat::test_that("Type-specific model can be fitted",{
  test.df <- make_test_df()
  test.in <- progressionEstimation::process_input_data(test.df)
  type_specific.out <- progressionEstimation::fit_progression_rate_model(
    test.in,
    type_specific = TRUE,
    location_adjustment = FALSE,
    stat_model = "poisson",
    num_iter = 10000,
    num_chains = 1
  )
  testthat::expect_match(type_specific.out@model_name,"type_specific_poisson")
})

# Test location-adjusted model fit
testthat::test_that("Location-adjusted model can be fitted",{
  test.df <- make_test_df()
  test.in <- progressionEstimation::process_input_data(test.df)
  location_adjusted.out <- progressionEstimation::fit_progression_rate_model(
    test.in,
    type_specific = FALSE,
    location_adjustment = TRUE,
    stat_model = "poisson",
    num_iter = 10000,
    num_chains = 1
  )
  testthat::expect_match(location_adjusted.out@model_name,"adjusted_null_poisson")
})

# Test type-specific location-adjusted model fit
testthat::test_that("Type-specific location-adjusted model can be fitted",{
  test.df <- make_test_df()
  test.in <- progressionEstimation::process_input_data(test.df)
  type_specific_location_adjusted.out <- progressionEstimation::fit_progression_rate_model(
    test.in,
    type_specific = TRUE,
    location_adjustment = TRUE,
    stat_model = "poisson",
    num_iter = 10000,
    num_chains = 1
  )
  testthat::expect_match(type_specific_location_adjusted.out@model_name,"adjusted_type_specific_poisson")
})

# Test type-specific location-adjusted negative binomial model fit
testthat::test_that("Type-specific location-adjusted model can be fitted",{
  test.df <- make_test_df()
  test.in <- progressionEstimation::process_input_data(test.df)
  type_specific_location_adjusted.out <-
    progressionEstimation::fit_progression_rate_model(
      test.in,
      type_specific = TRUE,
      location_adjustment = TRUE,
      stat_model = "negbin",
      num_iter = 10000,
      num_chains = 1
  )
  testthat::expect_match(type_specific_location_adjusted.out@model_name,"adjusted_type_specific_negbin")
})

# Test model comparisons with Bayes factors
testthat::test_that("Model comparison with Bayes factors",{
  test.df <- make_test_df()
  test.in <- progressionEstimation::process_input_data(test.df)
  null.out <- progressionEstimation::fit_progression_rate_model(
    test.in,
    type_specific = FALSE,
    location_adjustment = FALSE,
    stat_model = "poisson",
    num_iter = 10000,
    num_chains = 1
  )
  type.out <- progressionEstimation::fit_progression_rate_model(
    test.in,
    type_specific = TRUE,
    location_adjustment = FALSE,
    stat_model = "poisson",
    num_iter = 10000,
    num_chains = 1
  )
  model.comp <-
    progressionEstimation::compare_model_fits_with_bf(list(null.out,type.out))
  testthat::expect_equal(model.comp[1,2],0)
})

# Test model comparisons with LOO-CV
testthat::test_that("Model comparison with LOO-CV",{
  test.df <- make_test_df()
  test.in <- progressionEstimation::process_input_data(test.df)
  null.out <- progressionEstimation::fit_progression_rate_model(
    test.in,
    type_specific = FALSE,
    location_adjustment = FALSE,
    stat_model = "poisson",
    num_iter = 10000,
    num_chains = 1
  )
  type.out <- progressionEstimation::fit_progression_rate_model(
    test.in,
    type_specific = TRUE,
    location_adjustment = FALSE,
    stat_model = "poisson",
    num_iter = 10000,
    num_chains = 1
  )
  model.comp <-
    progressionEstimation::compare_model_fits_with_loo(list(null.out,type.out))
  testthat::expect_equal(model.comp[1,2],0)
})

# Test processing of null model
testthat::test_that("Null model can be processed",{
  test.df <- make_test_df()
  test.in <- progressionEstimation::process_input_data(test.df)
  null.out <- progressionEstimation::fit_progression_rate_model(
    test.in,
    type_specific = FALSE,
    location_adjustment = FALSE,
    stat_model = "poisson",
    num_iter = 10000,
    num_chains = 1
  )
  null.processed.out <-
    progressionEstimation::process_progression_rate_model_output(null.out,test.df)
  testthat::expect_gt(mean(null.processed.out$surveillance_population),mean(null.processed.out$nu))
})

# Test plotting of case-carrier predictions
testthat::test_that("Plotting of case and carrier predictions",{
  test.df <- make_test_df()
  test.in <- progressionEstimation::process_input_data(test.df)
  null.out <- progressionEstimation::fit_progression_rate_model(
    test.in,
    type_specific = FALSE,
    location_adjustment = FALSE,
    stat_model = "poisson",
    num_iter = 10000,
    num_chains = 1
  )
  null.processed.out <-
    progressionEstimation::process_progression_rate_model_output(null.out,test.df)
  predictions_plot <-
    progressionEstimation::plot_case_carrier_predictions(null.processed.out)
  expect_error(print(predictions_plot), NA)
  testthat::expect_match(null.out@model_name,"null_poisson")
})

# Test plotting of progression rates
testthat::test_that("Plotting of progression rates",{
  test.df <- make_test_df()
  test.in <- progressionEstimation::process_input_data(test.df)
  null.out <- progressionEstimation::fit_progression_rate_model(
    test.in,
    type_specific = FALSE,
    location_adjustment = FALSE,
    stat_model = "poisson",
    num_iter = 10000,
    num_chains = 1
  )
  null.processed.out <-
    progressionEstimation::process_progression_rate_model_output(null.out,test.df)
  rates_plot <-
    progressionEstimation::plot_progression_rates(null.processed.out)
  expect_error(print(rates_plot), NA)
})

# Test further plotting of progression rates
testthat::test_that("Further plotting of progression rates",{
  test.df <- make_test_df()
  test.in <- progressionEstimation::process_input_data(test.df)
  null.out <- progressionEstimation::fit_progression_rate_model(
    test.in,
    type_specific = FALSE,
    location_adjustment = FALSE,
    stat_model = "poisson",
    num_iter = 10000,
    num_chains = 1
  )
  null.processed.out <-
    progressionEstimation::process_progression_rate_model_output(null.out,test.df)
  rates_plot <-
    progressionEstimation::plot_progression_rates(null.processed.out,
                                                  use_sample_size = TRUE,
                                                  colour_col = "study"
                                                  )
  expect_error(print(rates_plot), NA)
})

# Test processing and plotting of type-specific model
testthat::test_that("Type-specific model can be processed",{
  test.df <- make_test_df()
  test.in <- progressionEstimation::process_input_data(test.df)
  type_specific.out <- progressionEstimation::fit_progression_rate_model(
    test.in,
    type_specific = TRUE,
    location_adjustment = FALSE,
    stat_model = "poisson",
    num_iter = 10000,
    num_chains = 1
  )
  type_specific.processed.out <-
    progressionEstimation::process_progression_rate_model_output(type_specific.out,test.df)
  testthat::expect_gt(mean(type_specific.processed.out$surveillance_population),mean(type_specific.processed.out$nu))
  type_specific_prediction_plot <-
    progressionEstimation::plot_case_carrier_predictions(type_specific.processed.out)
  expect_error(print(type_specific_prediction_plot), NA)
  type_specific_rates_plot <-
    progressionEstimation::plot_progression_rates(type_specific.processed.out)
  expect_error(print(type_specific_rates_plot), NA)
})

# Test processing and plotting of location-adjusted model
testthat::test_that("Location-adjusted model can be processed",{
  test.df <- make_test_df()
  test.in <- progressionEstimation::process_input_data(test.df)
  location_adjusted.out <- progressionEstimation::fit_progression_rate_model(
    test.in,
    type_specific = FALSE,
    location_adjustment = TRUE,
    stat_model = "poisson",
    num_iter = 10000,
    num_chains = 1
  )
  location_adjusted.processed.out <-
    progressionEstimation::process_progression_rate_model_output(location_adjusted.out,test.df)
  testthat::expect_gt(mean(location_adjusted.processed.out$surveillance_population),
                      mean(location_adjusted.processed.out$nu))
  location_adjusted_prediction_plot <-
    progressionEstimation::plot_case_carrier_predictions(location_adjusted.processed.out)
  expect_error(print(location_adjusted_prediction_plot), NA)
  location_adjusted_rates_plot <-
    progressionEstimation::plot_progression_rates(location_adjusted.processed.out)
  expect_error(print(location_adjusted_rates_plot), NA)
  location_adjusted_scale_factors_plot <-
    progressionEstimation::plot_study_scale_factors(location_adjusted.processed.out)
  expect_error(print(location_adjusted_scale_factors_plot), NA)
})

# Test processing and plotting of type-specific location-adjusted model
testthat::test_that("Location-adjusted model can be processed",{
  test.df <- make_test_df()
  test.in <- progressionEstimation::process_input_data(test.df)
  type_specific_location_adjusted.out <-
    progressionEstimation::fit_progression_rate_model(
      test.in,
      type_specific = TRUE,
      location_adjustment = TRUE,
      stat_model = "negbin",
      num_iter = 10000,
      num_chains = 1
    )
  type_specific_location_adjusted.processed.out <-
    progressionEstimation::process_progression_rate_model_output(type_specific_location_adjusted.out,test.df)
  testthat::expect_gt(mean(type_specific_location_adjusted.processed.out$surveillance_population),
                      mean(type_specific_location_adjusted.processed.out$nu))
  type_specific_location_adjusted_prediction_plot <-
    progressionEstimation::plot_case_carrier_predictions(type_specific_location_adjusted.processed.out)
  expect_error(print(type_specific_location_adjusted_prediction_plot), NA)
  type_specific_location_adjusted_rates_plot <-
    progressionEstimation::plot_progression_rates(type_specific_location_adjusted.processed.out)
  expect_error(print(type_specific_location_adjusted_rates_plot), NA)
  type_specific_location_adjusted_scale_factors_plot <-
    progressionEstimation::plot_study_scale_factors(type_specific_location_adjusted.processed.out)
  expect_error(print(type_specific_location_adjusted_scale_factors_plot), NA)
})

# Test processing and plotting of type-specific location-adjusted model
testthat::test_that("Location-specific rates can be extracted",{
  test.df <- make_test_df()
  test.in <- progressionEstimation::process_input_data(test.df)
  type_specific_location_adjusted.out <-
    progressionEstimation::fit_progression_rate_model(
      test.in,
      type_specific = TRUE,
      location_adjustment = TRUE,
      stat_model = "negbin",
      num_iter = 10000,
      num_chains = 1
    )
  type_specific_location_adjusted.processed.out <-
    progressionEstimation::process_progression_rate_model_output(type_specific_location_adjusted.out,test.df)
  location_specific_output_df <-
    progressionEstimation::get_type_invasiveness_for_study(
      test.df,
      type_specific_location_adjusted.out,
      study = "E.W.pre.PCV"
    )
  testthat::expect_equal(nrow(location_specific_output_df),2)
})

# https://r-pkgs.org/tests.html
# https://stackoverflow.com/questions/31038709/how-to-write-a-test-for-a-ggplot-plot
