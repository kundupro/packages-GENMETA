#' Simulated data on study information based.
#'
#' The underlying model is assumed to be logistic. Then, outcomes are simulated from a logistic model (maximal) where the covariates(Age, Height and Weight) are simulated from the same multivariate normal distribution (see \code{\link[GMeta]{reference_data}} ). In this example, three different studies are considered. Study 1, study 2 and, study 3 have information on Age and Height, Height and Weight, and, Weight and Age, respectively.
#' Sample sizes of the three studies are taken to be 300, 500 and 1000, respectively.
#' Estimates of the regression coefficients are obtained by fitting a logistic model using the glm function. For study 2, estimate of the variance-covariance matrix is set to NULL
#' and the sample size of the corresponding study is provided to reflect the scenario where it is difficult for the user to provide variance-covariance matrix.
#' \cr \cr NOTE: This is a simulated data on summary-level information (estimates of regression coefficients, variance-covariance matrices) across studies. In real scenario, the users, usually, have those information. They have to just put that information in study_info argument appropriately.
#'
#' @docType data
#'
#' @usage data(study_info)
#'
#' @format A list of lists containing study information on 3 studies.
#'
#' @references See the simulation section of the paper by Tang, R., Kundu, P. and Chatterjee, N. (2017) Generalized Meta-Analysis for Multivariate Regression Models Across Studies with Disparate Covariate Information. \href{https://arxiv.org/abs/1708.03818}{arXiv:1708.03818v1 [stat.ME]}.
#' @references For reference_data, see \code{\link[GMeta]{reference_data}} .
#' @keywords datasets
#'
#'
#' @examples
#' data(study_info)

"study_info"
