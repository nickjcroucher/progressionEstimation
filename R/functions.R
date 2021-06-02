require(tidyverse)
require(magrittr)
require(rstan)
require(bridgesampling)
require(loo)
require(kableExtra)
require(xlsx)

#' Generate data frame from input spreadsheet
#'
#' @param fn Name of spreadsheet file containing data
#' @param use_strain Boolean variable indicating whether the strain column should also be used for subtyping
#'
#' @return Data frame containing data extracted from spreadsheet
#' @export
#'
#' @examples
process_input_xlsx <- function(fn = "progression_estimation_input.xlsx", use_strain = FALSE) {
  max_col_num <- 7
  if (use_strain) {
    max_col_num <- 8
  }
  input_df <-
    xlsx::read.xlsx(
      file = fn,
      sheetIndex = 1,
      colIndex = 1:max_col_num,
      header = TRUE
    ) %>%
    dplyr::mutate(study = factor(study)) %>%
    dplyr::mutate(categorisation = factor(categorisation))
  return(input_df)
}

#' Generate model input from input data frame
#'
#' @param input_df Data frame containing details of case-carrier studies
#' @param subtype Name of column in input data frame used for subtyping information
#'
#' @return A list of lists used as an input to stan models
#' @export
#'
#' @examples
#' process_input_data(S_pneumoniae_infant_serotype, subtype = "categorisation")
process_input_data <- function(input_df, subtype = "categorisation") {
  i_values <- as.integer(input_df$study)
  j_values <- as.integer(input_df[, which(colnames(input_df) == subtype)])
  c_ij <- input_df$carriage
  d_ij <- input_df$disease
  n_i <- input_df$carriage_samples
  N_i <- input_df$surveillance_population
  t_i <- input_df$time_interval
  progression_rate_data <- list(
    i_max = max(i_values),
    j_max = max(j_values),
    n_obs = length(c_ij),
    i_values = i_values,
    j_values = j_values,
    c_ij = c_ij,
    d_ij = d_ij,
    n_i = n_i,
    N_i = N_i,
    t_i = t_i
  )
  return(progression_rate_data)
}

#' Fit a progression rate model to input data
#'
#' @param input_data
#' @param model
#' @param num_chains
#' @param num_iter
#' @param adapt_delta
#'
#' @return
#' @export
#'
#' @examples
fit_progression_rate_model<-function(input_data, model = "Poisson", num_chains = 4, num_iter = 1e4, adapt_delta = 0.8) {
  return(input_data)
}
