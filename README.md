# progressionEstimation
R and stan package for the estimation of microbial progression rates

This package uses Bayesian models implemented in stan to estimate the rates at which microbes progress from carriage to disease using case and carrier data.

[![DOI:10.1101/2021.01.08.425840](http://img.shields.io/badge/DOI-10.1101/2021.09.01.458483-B31B1B.svg)](https://doi.org/10.1101/2021.09.01.458483) [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5964023.svg)](https://doi.org/10.5281/zenodo.5964023) [![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-brightgreen.svg)](https://github.com/nickjcroucher/progressionEstimation/blob/master/LICENSE)    [![R-CMD-check](https://github.com/nickjcroucher/progressionEstimation/actions/workflows/check-standard.yaml/badge.svg)](https://github.com/nickjcroucher/progressionEstimation/actions/workflows/check-standard.yaml)   [![codecov](https://codecov.io/gh/nickjcroucher/progressionEstimation/branch/main/graph/badge.svg?token=CZ63KCRN63)](https://codecov.io/gh/nickjcroucher/progressionEstimation)

## Quick start

The case and carrier data can be input into the model using the spreadsheet `progression_estimation_input.xlsx`, or read into R and converted to a data frame or tibble with the same format. The required columns are:

- **study** - unique name of the case-carrier study

- **type** - the subdivision of the microbial population to which the isolates belong

- **carriage** - the number of isolates of this type recovered from asymptomatic carriers (or analogous samples, e.g. environmental samples)

- **disease** - the number of isolates of this type recovered from disease

- **carriage_samples** - the number of samples taken undertaken in the carriage study, including samples that were negative for the microbe of interest (corresponds to the number of individuals in the carriage study, if one sample was taken per individual)

- **surveillance_population** - the number of individuals of the appropriate demographic in the population under surveillance

- **time_interval** - the duration of the disease surveillance period

A further optional column is:

- **strain** - genetic background of the isolates

Further columns can be added to use as alternative typing schemes; these can be selected in model fits using the `type` option. An example of a correctly-formatted spreadsheet [can be found in the examples directory](vignettes/s_pneumoniae_sweden.xlsx). If you are comfortable using R, then the package can be installed using `devtools::install_github("https://github.com/nickjcroucher/progressionEstimation")`. If you are more of a beginner with R, then you can:

- download this Github repository as a zip file using the button above

- extract the directory

- fill in your data in the spreadsheet

- right-click on `run_progression_rate_analysis.R` and open with R

- paste the code into the R terminal - note the first time, the package will time some time to install and compile

- on publication, share your input data spreadsheet (raise an Issue on this repository or [email me](https://www.imperial.ac.uk/people/n.croucher)) so others can include your data and cite your work

## Meta-analysis with stored data sets

There are four case-carrier datasets currently included in the package:

- **S_pneumoniae_infant_serotype** - serotype information on *Streptococcus pneumoniae* case and carriage data from infants

- **S_pneumoniae_adult_serotype** - serotype information on *Streptococcus pneumoniae* carriage data from infants and case data from adults

- **S_pneumoniae_infant_strain** - serotype and strain information on *Streptococcus pneumoniae* case and carriage data from infants

- **S_pneumoniae_mixed_strain** - serotype and strain information on *Streptococcus pneumoniae* carriage data from infants and case data from infants and adults

## More detailed usage examples

For information on how to run alternative analyses, generate a greater range of plots, and use leave-one-out cross validation and Bayes factors for model comparisons, there are suggested workflows in the `vignettes` directory. The `examples` directory contains detailed analysis of *Streptococcus pneumoniae* invasiveness data.
