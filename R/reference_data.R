#' Simulated reference data on a set of 50 individuals.
#'
#' This is a data set containing simulated covariates, Age, Height and Weight, from a multivariate normal distribution with zero mean vector and variance-covariance matrix, matrix(c(1,0.3,0.6,0.3,1,0.1,0.6,0.1,1),3,3).
#' In other words, the covariates in the underlying population have unit variance with correlation between Age and Height, Height and Weight, and, Height and Age are 0.3, 0.1 and 0.6, respectively.
#' The first column of the data set contains vector of ones indicating the intercept.
#'
#' @docType data
#'
#' @usage data(reference_data)
#'
#' @format A data matrix containing 50 rows and 4 columns.
#'
#' @keywords datasets
#' @references See the simulation section of the paper by Tang, R., Kundu, P. and Chatterjee, N. (2017) Generalized Meta-Analysis for Multivariate Regression Models Across Studies with Disparate Covariate Information. \href{https://arxiv.org/abs/1708.03818}{arXiv:1708.03818v1 [stat.ME]}.
#'
#' @examples
#' data(reference_data)

"reference_data"
