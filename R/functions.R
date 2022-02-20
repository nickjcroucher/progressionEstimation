require(tidyverse)
require(magrittr)
require(rstan)
require(bridgesampling)
require(loo)
require(kableExtra)
require(xlsx)
require(cowplot)
require(ggrepel)

#' Generate data frame from input spreadsheet
#'
#' @param fn Name of spreadsheet file containing data
#' @param use_strain Boolean variable indicating whether the strain column should also be used for subtyping
#'
#' @return Data frame containing data extracted from spreadsheet
#' @export
#'
#' @importFrom magrittr %>%
#' @importFrom magrittr %<>%
#' @importFrom stringr str_trim
#'
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
    tidyr::drop_na() %>%
    dplyr::mutate(study = factor(str_trim(study))) %>%
    dplyr::mutate(type = factor(str_trim(type, side = "both"))) %>%
    dplyr::mutate(carriage = as.numeric(str_trim(carriage, side = "both"))) %>%
    dplyr::mutate(disease = as.numeric(str_trim(disease, side = "both"))) %>%
    dplyr::mutate(carriage_samples = as.numeric(str_trim(carriage_samples, side = "both"))) %>%
    dplyr::mutate(surveillance_population = as.numeric(str_trim(surveillance_population, side = "both"))) %>%
    dplyr::mutate(time_interval = as.numeric(str_trim(time_interval, side = "both")))
  if (use_strain) {
    input_df %<>%
      dplyr::mutate(strain = factor(str_trim(strain, side = "both")))
  }
  return(input_df)
}

combine_rows <- function(df, col_name = "type") {
  df %<>%
    dplyr::group_by(study,!!! dplyr::syms(col_name)) %>%
      dplyr::mutate(carriage = sum(carriage)) %>%
      dplyr::mutate(disease = sum(disease)) %>%
      dplyr::ungroup() %>%
      dplyr::distinct()
  return(df)
}

#' Generate model input from input data frame
#'
#' @param input_df Data frame containing details of case-carrier studies
#' @param type Name of column in input data frame used for typing information
#' @param use_strain Boolean specifying whether strain information should be used in addition to type information
#' @param combine_strain Boolean specifying whether strain information should be combined with type information
#' @param condense Boolean specifying whether repeated entries of the same type should be condensed into a single entry
#'
#' @return A list of lists used as an input to stan models
#' @export
#'
#' @importFrom rlang :=
#' @importFrom rlang !!
#'
process_input_data <- function(input_df, type = "type", use_strain = FALSE, combine_strain = FALSE, condense = FALSE) {
  if (!(type %in% colnames(input_df))) {
    stop("Type column not in input data")
  }
  # Process input data
  if (combine_strain | use_strain) {
    input_df %<>%
      tidyr::unite("combined",!!type,strain, sep='_', remove = FALSE) %>%
      dplyr::mutate(combined = factor(combined))
    if (condense) {
      input_df <- combine_rows(input_df, col_name = "combined")
    }
    if (combine_strain) {
      type = "combined"
    } else if (use_strain) {
      input_df %<>% dplyr::select(-combined)
    }
  } else if (condense) {
    input_df <- combine_rows(input_df %>% dplyr::select(study,
                                                        !!type,
                                                        carriage,
                                                        disease,
                                                        carriage_samples,
                                                        surveillance_population,
                                                        time_interval),
                             col_name = type)
  }
  # Convert to factors
  input_df %<>%
    dplyr::mutate(study = factor(study)) %>%
    dplyr::mutate((!!type) := factor(!!! dplyr::syms(type)))
  if ("strain" %in% colnames(input_df)) {
    input_df %<>%
      dplyr::mutate(strain = factor(strain))
  }
  # Calculate input
  input_df <- as.data.frame(input_df)
  i_values <- as.integer(input_df$study)
  j_values <- as.integer(input_df[, which(colnames(input_df) == type)])
  c_ij <- input_df$carriage
  d_ij <- input_df$disease
  n_i <- input_df$carriage_samples
  N_i <- input_df$surveillance_population
  t_i <- input_df$time_interval
  if (use_strain) {
    k_values <- as.integer(input_df$strain)
    progression_rate_data <- list(
      i_max = max(i_values),
      j_max = max(j_values),
      k_max = max(k_values),
      n_obs = length(c_ij),
      i_values = i_values,
      j_values = j_values,
      k_values = k_values,
      c_ij = c_ij,
      d_ij = d_ij,
      n_i = n_i,
      N_i = N_i,
      t_i = t_i
    )
  } else {
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
  }
  return(progression_rate_data)
}

#' Fit a progression rate model to case-carrier data
#'
#' @param input_data List of lists generated by `process_input_data`
#' @param type_specific Boolean specifying whether progression rates vary between types
#' @param location_adjustment Boolean specifying whether progression rates vary between locations
#' @param stat_model Whether progression to disease is "poisson" (Poisson process) or "negbin" (overdispersed negative binomial distribution)
#' @param strain_as_primary_type Whether strain should be used as the primary determinant of progression rate, and the other type used as the secondary determinant
#' @param strain_as_secondary_type Whether strain should be used as the secondary determinant of progression rate, and the other type used as the primary determinant
#' @param model_description Descriptive title of model
#' @param num_chains Number of MCMCs to be run for inference
#' @param num_iter Length of MCMCs
#' @param num_cores Number of threads used to calculate MCMCs
#' @param adapt_delta_value Target average acceptance probability of MCMCs (default = 0.8)
#' @param stepsize_value Initial MCMC step size (default = 1)
#' @param max_treedepth_value Depth of tree explored by MCMC sampler (default = 10)
#'
#' @return A stanfit object
#' @export
#'
fit_progression_rate_model<-function(input_data,
                                     type_specific = TRUE,
                                     location_adjustment = TRUE,
                                     stat_model = "poisson",
                                     strain_as_primary_type = FALSE,
                                     strain_as_secondary_type = FALSE,
                                     model_description = NULL,
                                     num_chains = 4,
                                     num_iter = 1e4,
                                     num_cores = 1,
                                     adapt_delta_value = 0.8,
                                     stepsize_value = 1,
                                     max_treedepth_value = 10) {
  # Validate input
  model_suffix = match.arg(stat_model, c("poisson","negbin"), several.ok = FALSE)
  if ((strain_as_primary_type | strain_as_secondary_type) & !("k_max" %in% names(input_data))) {
    stop("If strain to be used in typing, then needs to be in input data")
  }
  if (strain_as_primary_type & strain_as_secondary_type) {
    stop("Strain can only be primary type or seconday type, not both")
  }
  # Select model based on specifications
  model_prefix = "null"
  if (type_specific & location_adjustment) {
    model_prefix = "adjusted_type_specific"
  } else if (type_specific & !(location_adjustment)) {
    model_prefix = "type_specific"
  } else if (!(type_specific) & location_adjustment) {
    model_prefix = "adjusted_null"
  }
  # Adjust names if modifying serotype by strain
  if (strain_as_primary_type & type_specific) {
    model_prefix = paste0(model_prefix,"_strain_modified_by_type")
  } else if (strain_as_secondary_type & type_specific) {
    model_prefix = paste0(model_prefix,"_type_modified_by_strain")
  }
  # Complete model name
  model_name = paste0(model_prefix,'_',model_suffix)
  # Validate model name
  if (!(model_name %in% names(stanmodels))) {
    stop(paste(model_name,"not in list of valid model names"))
  }
  # Sample from model
  model_output<-
    rstan::sampling(stanmodels[[model_name]],
                    data = input_data,
                    iter = num_iter,
                    cores = num_cores,
                    chains = num_chains,
                    control = list(adapt_delta = adapt_delta_value,
                                   stepsize = stepsize_value,
                                   max_treedepth = max_treedepth_value)
    )
  # Rename if specified
  if (!(is.null(model_description))) {
    model_output@model_name <- model_description
  }
  # Return output
  return(model_output)
}

get_mean<-function(parameter,model) {
  return(as.numeric(rstan::summary(model,pars=c(parameter))$summary[,1]))
}

get_upper<-function(parameter,model) {
  return(as.numeric(rstan::summary(model,pars=c(parameter))$summary[,8]))
}

get_lower<-function(parameter,model) {
  return(as.numeric(rstan::summary(model,pars=c(parameter))$summary[,4]))
}

get_median<-function(parameter,model) {
  return(as.numeric(rstan::summary(model,pars=c(parameter))$summary[,6]))
}

#' Process the model output for downstream analysis
#'
#' @param model_output Stanfit object returned by model fitting
#' @param input_df Data frame used as input to model fitting
#' @param type Name of column used to define types
#' @param strain_as_primary_type Whether strain was used as the primary determinant of progression rate, and the other type used as the secondary determinant
#' @param strain_as_secondary_type Whether strain was used as the secondary determinant of progression rate, and the other type used as the primary determinant
#' @param combined_strain_type Whether strain and type were jointly used to subdivide the population
#' @param condense Whether data should be reclassified by combining repeated entries for the same type to match with model input
#'
#' @return A data frame
#' @export
#'
#' @importFrom stats setNames
#'
process_progression_rate_model_output<-function(model_output,
                                                input_df,
                                                type = "type",
                                                strain_as_primary_type = FALSE,
                                                strain_as_secondary_type = FALSE,
                                                combined_strain_type = FALSE,
                                                condense = FALSE) {
  # Add model name
  input_df %<>%
    dplyr::mutate(model_name = model_output@model_name)
  # Process input data
  if (strain_as_primary_type | strain_as_secondary_type | combined_strain_type) {
    input_df %<>%
      tidyr::unite("combined",!!type,strain, sep='_', remove = FALSE) %>%
      dplyr::mutate(combined = factor(combined))
    if (combined_strain_type) {
      type = "combined"
    }
    input_df %<>%
      dplyr::select(model_name,
                    study,
                    !!type,
                    carriage,
                    disease,
                    carriage_samples,
                    surveillance_population,
                    time_interval,
                    strain)
  }
  # Remove unused rows
  if (condense) {
    if (strain_as_primary_type | strain_as_secondary_type) {
      input_df <- combine_rows(input_df %>%
                                 dplyr::select(model_name,
                                               study,
                                               !!type,
                                               carriage,
                                               disease,
                                               carriage_samples,
                                               surveillance_population,
                                               time_interval,
                                               strain)) %>%
        dplyr::distinct()
    } else {
      input_df <- combine_rows(input_df %>%
                                 dplyr::select(model_name,
                                               study,
                                               !!type,
                                               carriage,
                                               disease,
                                               carriage_samples,
                                               surveillance_population,
                                               time_interval),
                               col_name = type)
    }
  }
  # Extract factor levels
  i_levels = levels(input_df %>% dplyr::pull(study))
  j_levels = levels(input_df %>% dplyr::pull(!!type))
  if (strain_as_primary_type | strain_as_secondary_type) {
    k_levels = levels(input_df %>% dplyr::pull(strain))
  }
  # Carriage prevalence estimates
  carriage_df <- data.frame(
    "rho" = get_median("rho_ij",model_output),
    "rho_lower" = get_lower("rho_ij",model_output),
    "rho_upper" = get_upper("rho_ij",model_output)
  )
  input_df %<>% dplyr::bind_cols(carriage_df)
  # Variation by location
  scale_parameter <- 1
  if ("gamma_i" %in% model_output@model_pars) {
    location_parameters <- data.frame(
      "study" = i_levels,
      "gamma" = get_median("gamma_i",model_output),
      "gamma_lower" = get_lower("gamma_i",model_output),
      "gamma_upper" = get_upper("gamma_i",model_output)
    )
  } else {
    location_parameters <- data.frame(
      "study" = i_levels,
      "gamma" = 1,
      "gamma_lower" = 1,
      "gamma_upper" = 1
    )
  }
  input_df %<>% dplyr::left_join(location_parameters, by = c("study"="study"))
  # Calculate invasiveness values
  nu_name = "nu"
  if ("nu_j" %in% model_output@model_pars) {
    nu_name = "nu_j"
  }
  progression_rates_df <- data.frame(
    "type" = j_levels,
    "nu" = as.numeric(get_median(nu_name,model_output)),
    "nu_lower" = as.numeric(get_lower(nu_name,model_output)),
    "nu_upper" = as.numeric(get_upper(nu_name,model_output))
  )
  input_df %<>% dplyr::left_join(progression_rates_df, by = setNames("type",type))

  if ("nu_k" %in% model_output@model_pars) {
    secondary_progression_rates_df <- data.frame(
      "type" = k_levels,
      "secondary_nu" = get_median("nu_k",model_output),
      "secondary_nu_lower" = get_lower("nu_k",model_output),
      "secondary_nu_upper" = get_upper("nu_k",model_output)
    )
    input_df %<>% dplyr::left_join(secondary_progression_rates_df, by = setNames("type","strain"))
  }

  if ("phi_nb" %in% model_output@model_pars) {
    precision_parameters_df <- data.frame(
      "phi" = get_median("phi_nb",model_output),
      "phi_lower" = get_lower("phi_nb",model_output),
      "phi_upper" = get_upper("phi_nb",model_output)
    )
    input_df %<>% dplyr::bind_cols(precision_parameters_df)
  }

  # Extract predictions and intervals
  input_df %<>%
    dplyr::mutate(carriage_prediction = get_median("c_ij_pred",model_output)) %>%
    dplyr::mutate(carriage_prediction_lower = get_lower("c_ij_pred",model_output)) %>%
    dplyr::mutate(carriage_prediction_upper =  get_upper("c_ij_pred",model_output)) %>%
    dplyr::mutate(disease_prediction = get_median("d_ij_pred",model_output)) %>%
    dplyr::mutate(disease_prediction_lower = get_lower("d_ij_pred",model_output)) %>%
    dplyr::mutate(disease_prediction_upper =  get_upper("d_ij_pred",model_output))

  # Add in absolute deviation
  input_df %<>%
    dplyr::mutate(carriage_abs_dev = abs(carriage - carriage_prediction)) %>%
    dplyr::mutate(disease_abs_dev = abs(disease - disease_prediction))

  return(input_df)
}

#' Function for plotting the observed and predicted case-carrier counts
#'
#' @param model_output_df Data frame include input data and model fit output
#' @param n_label Number of top-ranked observations to label
#' @param label_col Column to use for labels
#' @param legend Boolean specifying whether a legend should be included in the plot
#' @param just_legend Whether to only return the legend
#'
#' @return ggplot2 plot
#' @export
#'
#' @importFrom ggplot2 geom_abline
#' @importFrom ggplot2 geom_point
#'
plot_case_carrier_predictions <- function(model_output_df, n_label = 3, label_col = "type", legend = TRUE, just_legend = FALSE) {

  if (!("carriage_prediction" %in% colnames(model_output_df))) {
    stop("Need to include model output in data frame for plotting")
  }

  if (!is.null(label_col)) {
    carriage_labels <-
      model_output_df %>%
      dplyr::slice_max(carriage_prediction, n = n_label) %>%
      dplyr::select(!!! dplyr::syms(label_col), carriage, carriage_prediction)

    disease_labels <-
      model_output_df %>%
      dplyr::slice_max(disease_prediction, n = n_label) %>%
      dplyr::select(!!! dplyr::syms(label_col), disease, disease_prediction)
  }

  carriage_plot <-
    ggplot(model_output_df,
           aes(x = carriage,
               y = carriage_prediction,
               ymin = carriage_prediction_lower,
               ymax = carriage_prediction_upper)) +
    geom_abline(slope = 1, intercept = 0, lty = 2, colour = "coral") +
    ylab("Predicted carriage isolates") +
    xlab("Observed carriage isolates") +
    theme_bw()

  if (!is.null(label_col)) {
    carriage_plot <-
      carriage_plot +
      ggrepel::geom_text_repel(data = carriage_labels,
                               aes(x = carriage,
                                   y = carriage_prediction,
                                   label = get(label_col)),
                               alpha = 0.9,
                               force = 50,
                               inherit.aes = FALSE)
  }

  disease_plot <-
    ggplot(model_output_df,
           aes(x = disease,
               y = disease_prediction,
               ymin = disease_prediction_lower,
               ymax = disease_prediction_upper)) +
    geom_abline(slope = 1, intercept = 0, lty = 2, colour = "coral") +
    ylab("Predicted disease isolates") +
    xlab("Observed disease isolates") +
    theme_bw()

  if (!is.null(label_col)) {
    disease_plot <-
      disease_plot +
      ggrepel::geom_text_repel(data = disease_labels,
                               aes(x = disease,
                                   y = disease_prediction,
                                   label = get(label_col)),
                               alpha = 0.9,
                               force = 50,
                               inherit.aes = FALSE)
  }

  # Add in function for colouring by location if appropriate
  if (model_output_df$study %>% dplyr::n_distinct() > 1) {
    carriage_plot <- carriage_plot +
      geom_point(aes(color = study)) +
      geom_errorbar(aes(color = study), alpha = 0.75)
    disease_plot <- disease_plot +
      geom_point(aes(color = study)) +
      geom_errorbar(aes(color = study), alpha = 0.75)
    if (just_legend | legend) {
      carriage_plot <- carriage_plot +
        theme(legend.position = "bottom")
      disease_plot <- disease_plot +
        theme(legend.position = "bottom")
    } else {
      carriage_plot <- carriage_plot +
        theme(legend.position = "none")
      disease_plot <- disease_plot +
        theme(legend.position = "none")
    }
  } else {
    carriage_plot <- carriage_plot +
      geom_point(color = "blue") +
      geom_errorbar(color = "blue", alpha = 0.75)
    disease_plot <- disease_plot +
      geom_point(color = "blue") +
      geom_errorbar(color = "blue", alpha = 0.75)
  }

  if (just_legend) {
    return(cowplot::get_legend(carriage_plot))
  } else {
    return(cowplot::plot_grid(plotlist = list(carriage_plot, disease_plot)))
  }
}


#' Plot progression rate estimates
#'
#' @param model_output_df Data frame include input data and model fit output
#' @param type Name of column to be used on x axis
#' @param unit_time String specifying the time unit for the y axis label
#' @param type_name Name of typing scheme for x axis label
#' @param colour_col Name of the column used for colour scheme
#' @param colour_palette Named vector to be used for colour scheme
#' @param use_sample_size Whether to change point style based on sample size
#'
#' @return
#' @export
#'
#' @importFrom ggplot2 element_text
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_errorbar
#' @importFrom ggplot2 scale_y_continuous
#' @importFrom ggplot2 scale_shape_binned
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 theme_bw
#' @importFrom ggplot2 xlab
#' @importFrom ggplot2 ylab
#'
plot_progression_rates <- function(model_output_df, type = "type", unit_time = "unit time", type_name = "type",
                                   colour_col = NULL, colour_palette = NULL, use_sample_size = FALSE) {
  if (!(any(grepl("nu",colnames(model_output_df))))) {
    stop("Need to include model output in data frame for plotting")
  }
  progression_rate_values <- c("nu", "nu_lower", "nu_upper")
  if (type == "strain" & "secondary_nu" %in% colnames(model_output_df)) {
    progression_rate_values <- paste0("secondary_", progression_rate_values)
  }
  y_label_text = paste0("Progression rate (disease per carrier per ",unit_time,")")
  if ("secondary_nu" %in% colnames(model_output_df)) {
    y_label_text = paste0("Progression rate contribution (disease per carrier per ",unit_time,")")
  }

  if (use_sample_size) {
    model_output_df %<>%
      dplyr::group_by(!!! dplyr::syms(type)) %>%
      dplyr::mutate(num_observations = sum(carriage+disease)) %>%
      dplyr::ungroup()
  }

  model_output_df %<>%
    dplyr::group_by(!!! dplyr::syms(type)) %>%
    dplyr::slice_head(n=1) %>%
    dplyr::ungroup()

  if (is.null(colour_col)) {
    if (use_sample_size) {
      base_graph <-
        ggplot(model_output_df,
               aes(x = get(!!type),
                   y = get(progression_rate_values[1]),
                   ymin = get(progression_rate_values[2]),
                   ymax = get(progression_rate_values[3]),
                   shape = num_observations))
    } else {
      base_graph <-
        ggplot(model_output_df,
               aes(x = get(!!type),
                   y = get(progression_rate_values[1]),
                   ymin = get(progression_rate_values[2]),
                   ymax = get(progression_rate_values[3])))
    }
  } else {
    if (use_sample_size) {
      base_graph <-
        ggplot(model_output_df,
               aes(x = get(!!type),
                   y = get(progression_rate_values[1]),
                   ymin = get(progression_rate_values[2]),
                   ymax = get(progression_rate_values[3]),
                   colour = get(colour_col),
                   fill = get(colour_col),
                   shape = num_observations))
    } else {
      base_graph <-
        ggplot(model_output_df,
               aes(x = get(!!type),
                   y = get(progression_rate_values[1]),
                   ymin = get(progression_rate_values[2]),
                   ymax = get(progression_rate_values[3]),
                   colour = get(colour_col),
                   fill = get(colour_col)))
    }
  }

  point_graph <-
    base_graph +
    geom_point() +
    geom_errorbar() +
    ylab(y_label_text) +
    xlab(type_name) +
    scale_y_continuous(trans ="log10") +
    theme_bw()

  if (!is.null(colour_col) & !is.null(colour_palette)) {
    point_graph <- point_graph +
      scale_colour_manual(values = colour_palette,
                          name = colour_col) +
      scale_fill_manual(values = colour_palette,
                        name = colour_col) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
            legend.position = "bottom")
  } else if (!is.null(colour_col)) {
    point_graph <- point_graph +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
            legend.position = "bottom")
  } else {
    point_graph <- point_graph +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
  }

  if (use_sample_size) {
    point_graph <- point_graph +
      scale_shape_binned(name = "Number of\nisolates",
                         breaks = c(5,10,25,50,100))
  }

  return(point_graph)
}

#' Compare model fits using Bayes factors
#'
#' @param model_list List of stan fit objects
#' @param num_iter Number of iterations used for bridge sampling
#'
#' @return Data frame containing Bayes factors
#' @export
#'
compare_model_fits_with_bf <- function(model_list, num_iter = 1e3, num_threads = 1) {
  model_names <- unlist(lapply(model_list, getElement, "model_name"))
  model_lml <- lapply(model_list,
                      bridgesampling::bridge_sampler,
                      maxiter = num_iter,
                      cores = num_threads,
                      silent = TRUE)
  model_lml_values <- -1*unlist(lapply(model_lml, getElement, "logml"))
  lml_value_order <- order(model_lml_values)
  model_lml <- model_lml[lml_value_order]
  model_names <- model_names[lml_value_order]
  bf_values <- c()
  for (x in 1:length(model_names)) {
    bf_values <- c(bf_values,
                   bridgesampling::bf(model_lml[[x]], model_lml[[1]], log = TRUE)$bf
    )
  }
  bf_df <- data.frame(
    "Model" = model_names,
    "log_Bayes_factor" = bf_values
  )
  return(bf_df)
}

#' Compare model fits using leave-one-out cross-validation
#'
#' @param model_list List of stan fit objects
#' @param log_lik_param Name of log likelihood parameter
#' @param use_moments Boolean specifying whether to use moment matching
#' @param num_threads Number of core to use
#'
#' @return Data frame containing cross-validation values
#' @export
#'
compare_model_fits_with_loo <- function(model_list,
                                        log_lik_param = "log_lik",
                                        use_moments = FALSE,
                                        num_threads = 1) {
  model_loo <- lapply(model_list, rstan::loo,
                      pars = log_lik_param,
                      moment_match = use_moments,
                      cores = num_threads)
  loo_comparisons <- loo::loo_compare(model_loo)
  model_names <- unlist(lapply(model_list, getElement, "model_name"))
  rownames(loo_comparisons) <- model_names[as.integer(gsub("model","",rownames(loo_comparisons)))]
  return(loo_comparisons)
}

#' Combine input data frames for meta-analyses
#'
#' @param new_df Data frame with new studies
#' @param old_df Data frame with old studies
#'
#' @return Combined data frame
#' @export
#'
combine_with_existing_datasets <- function(new_df, old_df) {
  new_studies <-
    new_df %>% dplyr::select(study) %>% dplyr::distinct() %>% dplyr::pull()
  old_studies <-
    old_df %>% dplyr::select(study) %>% dplyr::distinct() %>% dplyr::pull()
  if (length(intersect(new_studies,old_studies)) > 0) {
    stop("Names of studies in new data must not be present in old studies")
  }
  combined_df <- dplyr::bind_rows(old_df, new_df)
  return(combined_df)
}

#' Plot study scale factors
#'
#' @param model_output_df Data frame including input data and model fit output
#'
#' @return ggplot2 object
#' @export
#'
plot_study_scale_factors <- function(model_output_df) {
  if (!("carriage_prediction" %in% colnames(model_output_df))) {
    stop("Need to include model output in data frame for plotting")
  }

  model_output_df %<>%
    dplyr::group_by(study) %>%
    dplyr::slice_head(n=1) %>%
    dplyr::ungroup()

  ggplot(model_output_df,
         aes(x = study, y = gamma, ymin = gamma_lower, ymax = gamma_upper)) +
    geom_point() +
    geom_errorbar() +
    ylab(paste0("Study scale factor")) +
    xlab("Study") +
    scale_y_continuous(trans = "log10") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
}

#' Validate input dataset for progrsssion estimation
#'
#' @param df Input dataframe containing case and carrier data
#'
#' @return Nothing
#' @export
#'
validate_progression_estimation_dataset <- function(df) {
  # Check column names
  if (!all(c("study","carriage","disease","carriage_samples","surveillance_population","time_interval") %in% colnames(df))) {
    stop('Columns required in dataset: "study","carriage","disease","carriage_samples","surveillance_population","time_interval"')
  }
  # Check consistency of statistics
  filtered_df <-
    df %>%
    dplyr::select(study,carriage_samples,surveillance_population,time_interval) %>%
    dplyr::distinct()
  n_distinct_values <- NULL
  for (var in c("study","carriage_samples","surveillance_population","time_interval")) {
    n_distinct_values <- c(n_distinct_values, filtered_df %>% dplyr::select(!!var) %>% dplyr::n_distinct())
  }
  if (max(n_distinct_values %>% dplyr::n_distinct()) > n_distinct_values[1]) {
    stop("Each study have a unique name associated with consistent carriage sample, surveillance population and time interval values")
  }
}

#' Generate invasiveness estimates specific to a particular study within a meta-analysis
#'
#' @param df Data frame containing input data used for model fitting
#' @param fit stan object corresponding to model fit
#' @param study Study for which invasiveness should be estimated
#' @param type Column name used to specify types for which invasiveness was estimated
#' @param use_strain_invasiveness Whether invasiveness should be returned for strains, if invasiveness was estimated for type and strain
#'
#' @return Data frame containing adjusted estimates of invasiveness for a study
#' @export
#'
#' @importFrom stats quantile
#'
get_type_invasiveness_for_study <- function(df, fit, study = NULL, type = "type", use_strain_invasiveness = FALSE) {

  # Check study name in list
  if (!(study %in% levels(df$study))) {
    stop(paste(study,"not found in list of studies:",levels(df$study)))
  }

  # Check model has a study adjustment
  if (!("gamma_i") %in% fit@model_pars) {
    stop("Function only necessary if a model features a study-specific adjustment factor")
  }

  # Get level
  i_level <- grep(study,levels(df$study))

  # Extract MCMC values
  invasiveness_param <- "nu"
  if ("nu_j" %in% fit@model_pars) {
    invasiveness_param <- "nu_j"
  }
  if (use_strain_invasiveness) {
    invasiveness_param <- "nu_k"
  }
  overall_invasiveness_mat <-
    as.matrix(fit, pars = c(invasiveness_param))
  study_adjustment <-
    as.matrix(fit, pars = c(paste0("gamma_i[",i_level,"]")))

  # Get product
  study_adjusted_invasiveness_mat <-
    t(t(overall_invasiveness_mat) * as.vector(study_adjustment))

  # Create the output data frame
  study_invasiveness_df <- df[df$study == study,]

  # Invasiveness values
  invasiveness_df <-
    data.frame(
      "nu" = apply(study_adjusted_invasiveness_mat, 2, mean),
      "nu_lower" = apply(study_adjusted_invasiveness_mat, 2, quantile, probs = c(0.025)),
      "nu_upper" = apply(study_adjusted_invasiveness_mat, 2, quantile, probs = c(0.975))
    )
  if (use_strain_invasiveness) {
    invasiveness_df[["strain"]] <- levels(df$strain)
  } else {
    invasiveness_df[[type]] <- levels(df[[type]])
  }
  invasiveness_df[["study"]] <- study

  # Add sample size
  if (use_strain_invasiveness) {
    study_invasiveness_df %<>%
      dplyr::left_join(
        invasiveness_df,
        by = c("study","strain")
      )
  } else {
    study_invasiveness_df %<>%
      dplyr::left_join(
        invasiveness_df,
        by = c("study",get("type"))
      )
  }

  # return
  return(study_invasiveness_df)
}
