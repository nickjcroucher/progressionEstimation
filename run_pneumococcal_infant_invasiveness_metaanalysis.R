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
require(gtools)

# Fucntions for ordering serotypes
mixedrank = function(x) order(gtools::mixedorder(as.character(x)))

get_serotype_levels <- function(df) {
  s_levels <-
    df %>%
      dplyr::select(type,classification) %>%
      dplyr::distinct() %>%
      group_by(classification) %>%
      dplyr::arrange(mixedrank(type),
                     .by_group = TRUE) %>%
      ungroup() %>%
      dplyr::select(type) %>%
      dplyr::pull()
  return(s_levels)
}

# Define serotype categories and palettes
vaccine_types <- list(
  "PCV7" = c("4","6B","9V","14","18C","19F","23F"),
  "PCV10" = c("1","5","7F"),
  "PCV13" = c("6A","3","19A"),
  "PCV15" = c("22F","33F"),
  "PCV20" = c("8","10A","11A","12F","15B/C"),
  "PPV23" = c("2","9N","17F","20")
)

vaccine_colours <- c("No PCV" = "blue",
                     "PCV7" = "red",
                     "PCV10" = "orange",
                     "PCV13" = "coral2",
                     "PCV15" = "pink",
                     "PCV20" = "purple",
                     "PPV23" = "black")

# Load spreadsheet
new_S_pneumoniae_infant_serotype <- progressionEstimation::process_input_xlsx("progression_estimation_input_test.xlsx", use_strain = FALSE)
new_studies <- unique(new_S_pneumoniae_infant_serotype$study)

# Merge datasets
S_pneumoniae_infant_serotype <- progressionEstimation::combine_with_existing_datasets(new_S_pneumoniae_infant_serotype,
                                                                                      S_pneumoniae_infant_serotype)

# Convert to stan input
input_data <- progressionEstimation::process_input_data(S_pneumoniae_infant_serotype)

# Fit model
model_fit <-
  progressionEstimation::fit_progression_rate_model(input_data,
                                                    type_specific = TRUE,
                                                    location_adjustment = TRUE,
                                                    num_chains = 2,
                                                    num_iter = 1e4)

# Process model output
model_output_df <-
  progressionEstimation::process_progression_rate_model_output(model_fit,
                                                               S_pneumoniae_infant_serotype) %>%
  dplyr::mutate("classification" = dplyr::case_when(
    type %in% vaccine_types[["PCV7"]] ~ "PCV7",
    type %in% vaccine_types[["PCV10"]] ~ "PCV10",
    type %in% vaccine_types[["PCV13"]] ~ "PCV13",
    type %in% vaccine_types[["PCV15"]] ~ "PCV15",
    type %in% vaccine_types[["PCV20"]] ~ "PCV20",
    type %in% vaccine_types[["PPV23"]] ~ "PPV23",
    TRUE ~ "No PCV"
  )
  ) %>%
  dplyr::mutate(classification = factor(classification, levels = c("PCV7",
                                                                   "PCV10",
                                                                   "PCV13",
                                                                   "PCV15",
                                                                   "PCV20",
                                                                   "PPV23",
                                                                   "No PCV")))

# Plot validation measures
validation_plots <-
  cowplot::plot_grid(
    plotlist = list(
      plot(model_fit, plotfun = "rhat", binwidth = 0.00005),
      rstan::traceplot(model_fit, pars = "lp__")
    ),
    labels= "AUTO"
  )
ggsave("model_convergence_test_plots.png",
       height = 8,
       width = 12)

# Order serotypes
serotype_levels <- get_serotype_levels(model_output_df)

model_output_df %<>%
  dplyr::mutate(type = factor(type, levels = serotype_levels))

# Plot invasiveness
progression_rate_plot <-
  progressionEstimation::plot_progression_rates(model_output_df,
                                                unit_time = "year",
                                                type_name = "Serotype",
                                                colour_col = "classification",
                                                colour_palette = vaccine_colours,
                                                use_sample_size = TRUE)
ggsave("progression_rate_plot.png",
       height = 8,
       width = 12)

# Save invasiveness values
write.csv(model_output_df,
          file = "progression_rate_analysis_output.csv",
          row.names = F,
          quote = F)

# Plot study adjustments
scale_factor_plot <-
    progressionEstimation::plot_study_scale_factors(model_output_df)
  ggsave("study_adjustment_factor_plot.png",
         height = 8,
         width = 8)

# Iterate to plot study specific estimates
for (new_study in new_studies) {
  study_output <- progressionEstimation::get_type_invasiveness_for_study(study = new_study,
                                                                         S_pneumoniae_infant_serotype,
                                                                         model_fit) %>%
    dplyr::mutate("classification" = dplyr::case_when(
      type %in% vaccine_types[["PCV7"]] ~ "PCV7",
      type %in% vaccine_types[["PCV10"]] ~ "PCV10",
      type %in% vaccine_types[["PCV13"]] ~ "PCV13",
      type %in% vaccine_types[["PCV15"]] ~ "PCV15",
      type %in% vaccine_types[["PCV20"]] ~ "PCV20",
      type %in% vaccine_types[["PPV23"]] ~ "PPV23",
      TRUE ~ "No PCV"
    )
    ) %>%
    dplyr::mutate(classification = factor(classification, levels = c("PCV7",
                                                                     "PCV10",
                                                                     "PCV13",
                                                                     "PCV15",
                                                                     "PCV20",
                                                                     "PPV23",
                                                                     "No PCV")))
  serotype_levels <- get_serotype_levels(study_output)

  study_output %<>%
    dplyr::mutate(type = factor(type, levels = serotype_levels))
  study_invasiveness_plot <- progressionEstimation::plot_progression_rates(study_output,
                                                    unit_time = "year",
                                                    type_name = "Serotype",
                                                    colour_col = "classification",
                                                    colour_palette = vaccine_colours,
                                                    use_sample_size = TRUE)
  ggsave(study_invasiveness_plot,
         file = paste0(new_study,"_progression_rate_estimates.png"),
         width = 12,
         height = 8)
  write.csv(study_output,
            file = paste0(new_study,"_progression_rate_estimates.csv"),
            quote = F,
            row.names = F)
}
