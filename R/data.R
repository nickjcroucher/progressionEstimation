#' @title Data from case and carrier studies of *S. pneumoniae* serotypes in infant populations
#'
#' @description A data set with counts of serotype observations in infant carriers and infant disease cases, and the associated study information
#'
#' @format A data frame with 7 variables:
#' \describe{
#'   \item{study}{Unique name of the case-carrier study.}
#'   \item{categorisation}{The label of the subdivision of the microbial population (serotype).}
#'   \item{carriage}{The number of isolates of this subtype recovered from asymptomatic carriers (nasopharyngeal swabs of infants).}
#'   \item{disease}{The number of isolates of this subtype recovered from disease (invasive disease in infants).}
#'   \item{carriage_samples}{The number of individuals in the carriage study, including negative samples.}
#'   \item{surveillance_population}{The number of individuals of the relevant demographic in the population under surveillance.}
#'   \item{time_interval}{The duration of the disease surveillance period.}
#' }
#' @source <https://github.com/nickjcroucher/progressionEstimation>
"S_pneumoniae_infant_serotype"
#' @title Data from case and carrier studies of *S. pneumoniae* serotypes in adult populations
#'
#' @description A data set with counts of serotype observations in infant carriers and adult disease cases, and the associated study information
#'
#' @format A data frame with 7 variables:
#' \describe{
#'   \item{study}{Unique name of the case-carrier study.}
#'   \item{categorisation}{The label of the subdivision of the microbial population (serotype).}
#'   \item{carriage}{The number of isolates of this subtype recovered from asymptomatic carriers (nasopharyngeal swabs of infants).}
#'   \item{disease}{The number of isolates of this subtype recovered from disease (invasive disease in adult).}
#'   \item{carriage_samples}{The number of individuals in the carriage study, including negative samples.}
#'   \item{surveillance_population}{The number of individuals of the relevant demographic in the population under surveillance.}
#'   \item{time_interval}{The duration of the disease surveillance period.}
#' }
#' @source <https://github.com/nickjcroucher/progressionEstimation>
"S_pneumoniae_adult_serotype"

