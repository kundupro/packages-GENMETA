#' Auxiliary for controlling the IRWLS algorithm
#'
#' This is an auxiliary function for the iteratively reweighted least squares algorithm for GMeta.
#' This is used internally by the myoptim function, but can be used by the user to create a control argument in the GMeta function
#' @param epsilon a positive numeric indicating converegence tolerence; the algotithm stops when the absolute difference between the estimates in current and previous step is less than epsilon, i.e, \eqn{|estimate_new - estimate_old| < \epsilon}
#' @param maxit a positive numeric indicating the maximum number of iterations to be used in the algorithm. Default is 1000.
#' @return A list with components named as the arguments.
#' @examples
#' control <- GMeta.control(1e-08, 100)
#' @export
GMeta.control <- function(epsilon = 1e-06, maxit = 1000)
{
  return(list(epsilon = epsilon, maxit = maxit))
}
