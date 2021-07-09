# Load stan package
if (!("progressionEstimation" %in% rownames(installed.packages()))) {
  devtools::install_github("https://github.com/nickjcroucher/progressionEstimation")
}
library(progressionEstimation)
pkgbuild::compile_dll()
roxygen2::roxygenize(package.dir = "..")

# Load dependencies
require(tidyverse)
require(magrittr)
require(rstan)
require(cowplot)
require(ggrepel)
require(xlsx)

# Load spreadsheet
input_df <- progressionEstimation::process_input_xlsx("progression_estimation_input.xlsx", use_strain = FALSE)

# Convert to stan input
input_data <- progressionEstimation::process_input_data(input_df)

# Test whether there are multiple studies in the input data
num_studies <-
  input_df %>%
    dplyr::select(study) %>%
    dplyr::n_distinct()

# Fit model
model_fit <-
  progressionEstimation::fit_progression_rate_model(input_data,
                                                 type_specific = TRUE,
                                                 location_adjustment = (num_studies > 1),
                                                 num_chains = 2,
                                                 num_iter = 1e4)

# Process model output
model_output_df <-
  progressionEstimation::process_progression_rate_model_output(model_fit,
                                                               input_df)

# Plot validation measures
validation_plots <-
  cowplot::plot_grid(
    plotlist = list(
      plot(model_fit, plotfun = "rhat", binwidth = 0.00005),
      rstan::traceplot(fit, pars = "lp__")
    ),
    labels= "AUTO"
  )
ggsave("model_convergence_test_plots.png",
       height = 8,
       width = 12)

# Plot invasiveness
progression_rate_plot <-
  progressionEstimation::plot_progression_rates(model_output_df,
                                                unit_time = "year",
                                                type_name = "type")
ggsave("progression_rate_plot.png",
       height = 8,
       width = 12)

# Save invasiveness values
write.csv(model_output_df,
          file = "progression_rate_analysis_output.csv",
          row.names = F,
          quote = F)

# Plot study adjustments
if (num_studies > 1) {
  scale_factor_plot <-
    progressionEstimation::plot_study_scale_factors(model_output_df)
  ggsave("study_adjustment_factor_plot.png",
         height = 8,
         width = 8)
}
