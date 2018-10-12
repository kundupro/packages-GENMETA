\name{MetaG.control}
\alias{MetaG.control}
\title{Auxiliary for controlling the IRWLS algorithm}
\usage{
  MetaG.control(epsilon = 1e-06, maxit = 1000)
}
\arguments{
  \item{epsilon}{a positive numeric indicating converegence
  tolerence; the algorithm stops when the absolute
  difference between the estimates in current and previous
  step is less than epsilon, i.e, \eqn{|estimate_{new} -
  estimate_{old}| < \epsilon}}

  \item{maxit}{a positive number indicating the maximum
  number of iterations to be used in the algorithm. Default
  is 1000.}
}
\value{
  A list with components named as the arguments.
}
\description{
  This is an auxiliary function for the iteratively
  reweighted least squares algorithm for GMeta. This is
  used internally by the myoptim function, but can be used
  by the user to create a control argument in the MetaG
  function
}
\examples{
control <- MetaG.control(1e-08, 100)
}

