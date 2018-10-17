#' Implementing Generalized Meta-analysis
#'
#' Generalized Meta-analysis(GENMETA) is an approach for combining information on multivariate regression parameters across multiple different studies which have different, but, possibly overlapping information on subsets of covariates.
#' GENMETA implements the generalized meta-analysis using IRWLS algorithm.
#' @param study_info a list of lists containing information about the studies; the main list contains a list for each study, which must have the fields:
#' \itemize{
#' \item{"Coeff": a named numeric vector containing the estimates of regression parameters (including intercept) where the names identify the covariates. For example, names(study_info$Coeff) <- c("(Intercept)", "Age", "Height", "Weight").}
#' \item{"Covariance": a matrix containing an estimate of variance-covariance matrix of the regression parameters. This can be NULL if the "Sample_size" is provided.}
#' \item{"Sample_size": a numeric containing sample size of the study. This can be NULL if the "Covariance" is provided.}
#' }
#' @param ref_dat a data matrix containing all the distinct covariates across studies from a reference set of individuals. This is used for estimating joint distribution of the covariates. The data matrix must have the vector of ones as its first column. The column names of the data matrix should match the names of the covariates from the studies.
#' @param model a description of the type of regression model; this is a character string naming the regression model. The current version is for "logistic" and "linear".
#' @param variable_intercepts an optional logical (applicable only when the model is "logistic"); if TRUE, the intercepts of the true models for each of the studies are assumed to be different. Default is FALSE.
#' @param initial_val an optional numeric vector containing initial values for the maximal model parameters which is needed for the IRWLS algorithm. Default is set to the one obtained from standard meta-analysis of the available estimates for each of the parameters across studies.
#' @param control an optional list containing the epsilon (positive numeric) and maxiter (positive number) needed for convergence of the algorithm. Default epsilon and maximum iterations are 1e-06 and 1000, respectively. For creating a control argument for GENMETA, see \code{\link[GENMETA]{GENMETA.control}}.
#' @details Generalized Meta-analysis (GENMETA) is a tool that allows researchers to quickly build models for multivariate meta-analysis in the presence of disparate covariate information across studies. It is implemented based on mainly two input arguments:
#' \itemize{
#' \item{Information on the model parameters from each of the studies.}
#' \item{Reference data for estimation of the joint distribution of all the distinct covariates across studies.}}
#' The software provides flexibility to the users to choose the intercepts to be different (when the model is logistic) across studies through the input argument, variable_intercepts.
#' It also allows estimation of the regression parameters, only from the sample sizes of the studies when it is difficult to obtain estimate of the variance-covariance matrices.
#' \cr \cr \bold{Note}: GENMETA will not work if both the estimates of the covariance matrix and the sample size are NULL.
#' @details When the model is "linear", it is assumed that the outcome is standardized to have unit variance.
#' For more details on the IRWLS, see References.
#' @return An object of class "GENMETA" is a list containing GENMETA estimate, its variance-covariance matrix and estimates the residual variance in the case of "linear" model .
#' \item{Est.coeff}{a numeric vector containing the estimated regression coefficients of the maximal model using optimal weighting matrix.}
#' \item{Est.var.cov}{a matrix containing estimate of variance-covariance matrix of the corresponding GENMETA estimator.}
#' \item{Res.var}{a numeric containing the residual variance of the maximal model when it is linear. It is calculated from the formula : \eqn{1 - \hat{\beta}_{GENMETA}^Tvar(X)\hat{\beta}_{GENMETA}} which is derived by assuming the outcomes to have unit variance. \eqn{var(X)} is calculated from reference data. Res.var is NA when the model is "logistic".}
#' \item{iter}{a numeric containing the number of iterations used in the algorithm}
#' \item{call}{the matched call}
#' \cr The function \code{\link[GENMETA]{GENMETA.summary}} prints a summary of the results obtained from GENMETA.
#' \cr The function \code{\link[GENMETA]{GENMETA.plot}} plots the estimate of the parameter from each of the studies, the summary measure(GENMETA estimate) and their confidence intervals. 
#' @keywords Generalized Meta Analysis
#' @references Tang, R., Kundu, P. and Chatterjee, N. (2017) Generalized Meta-Analysis for Multivariate Regression Models Across Studies with Disparate Covariate Information. \href{https://arxiv.org/abs/1708.03818}{arXiv:1708.03818v1 [stat.ME]}.
#' @seealso \code{\link[GENMETA]{GENMETA.summary}}, \code{\link[GENMETA]{GENMETA.plot}}.
#' @examples
#' # This example shows how to create the inputs GENMETA and then implement generalized meta-analysis
#' #########################
#' ### Basic setting #######
#' #########################
#'d.X = 3 # number of covariates.
#' mu = matrix(rep(0,d.X), nrow=d.X) # mean vector of the covariates.

#' r1 = 0.3 # correlation coefficient of the covariates.
#' r2 = 0.6
#' r3 = 0.1
#' Sigma = matrix(
#'   c(1, r1, r2,
#'     r1, 1, r3,
#'     r2, r3, 1),
#'   nrow=d.X,
#'   ncol=d.X) # covariance matrix of the covariates.


#' beta.star = matrix(c(-1.2, log(1.3), log(1.3), log(1.3)),nrow = d.X+1) # beta.star
#' #beta.star = matrix(c(-1.2, 0.26, 0.26, 0.26),nrow = d.X+1) # beta.star
#' #beta.star = matrix(c(-3, 1, 2, 3),nrow = d.X+1) # beta.star

#' n1 = 300 # sample size of the 1st data set.
#' n2 = 500 # 2nd
#' n3 = 1000 # 3rd

#' n = 50
#' sim=1
#' set.seed(sim)
# Generate the reference data set
#' X.rf = MASS::mvrnorm(n = n, mu, Sigma)


# Generate data set 1. m1 means model 1.
#' X.m1 = MASS::mvrnorm(n = n1, mu, Sigma) # Generate the covariates.
#' X.m1.1 = cbind(rep(1, n1), X.m1) # Add a column of 1's to X.m1.
#' p.m1 = 1/(1+exp(-X.m1.1%*%beta.star)) # the vector of probabilities
#' Y.m1 = rbinom(n1, size=1, p.m1) # the Bernoulli responses
#' # print(p.m1[1])
#' # print(mean(Y.m1))
#' # print(mean(p.m1))

#' # Generate data set 2. m1 means model 2.
#' X.m2 = MASS::mvrnorm(n = n2, mu, Sigma)
#' X.m2.1 = cbind(rep(1, n2), X.m2)
#' p.m2 = 1/(1+exp(-X.m2.1%*%beta.star))
#' Y.m2 = rbinom(n2, size=1, p.m2)

# Generate data set 3. m1 means model 3.
#' X.m3 = MASS::mvrnorm(n = n3, mu, Sigma)
#' X.m3.1 = cbind(rep(1, n3), X.m3)
#' p.m3 = 1/(1+exp(-X.m3.1%*%beta.star))
#' Y.m3 = rbinom(n3, size=1, p.m3)

#' #####
#' ### Create data sets in the format of data frame.
#' #####
#' data.m1 = data.frame(Y=Y.m1, X.m1)
#' data.m2 = data.frame(Y=Y.m2, X.m2)
#' data.m3 = data.frame(Y=Y.m3, X.m3)
# str(data.m1)



#' #####
#' ### Apply logistic regression with reduced models to the data sets
#' #####
#' logit.m1 <- glm(Y ~ X1 + X2, data = data.m1, family = "binomial")
#' # print(logit.m1)
#' if(logit.m1$converged == FALSE)
#' {
#'   print("glm for logit.m1 is not convergent.")
#'   next
#' }

#' logit.m2 <- glm(Y ~ X2 + X3, data = data.m2, family = "binomial")
#' # print(logit.m2)
#' if(logit.m2$converged == FALSE)
#' {
#'   print("glm for logit.m2 is not convergent.")
#'   next
#' }

#' logit.m3 <- glm(Y ~ X1 + X3, data = data.m3, family = "binomial")
#' # print(logit.m3)
#' if(logit.m3$converged == FALSE)
#' {
#'   print("glm for logit.m3 is not convergent.")
#'   next
#' }


#' #####
#' ### Obtain the estimators of the parameters in the reduced models.
#' #####
#' theta.m1 = logit.m1$coefficients
#' theta.m2 = logit.m2$coefficients
#' theta.m3 = logit.m3$coefficients



#' #####
#' ### Find the covariance matrix estimators for the reduced models
#' #####



#' #####
#' # Basic notations for inputs
#' #####

#' K = 3 # Number of data sets

#' A1 = c(1, 2) # index set A1, the indexes of the covariates of data set 1.
#' A2 = c(2, 3) # index set A2
#' A3 = c(1, 3) # index set A3


#' X.m1.used = cbind(rep(1, n1), X.m1[, A1, drop=FALSE])
#' X.m2.used = cbind(rep(1, n2), X.m2[, A2, drop=FALSE])
#' X.m3.used = cbind(rep(1, n3), X.m3[, A3, drop=FALSE])
#' # str(X.m1.used)
#' # str(X.m2.used)
#' # str(X.m3.used)

#' ##### Find Sigma.m1

#' T.1 = matrix(rep(0, (length(A1)+1)^2), nrow=length(A1)+1)
#' T.2 = T.1

#' for (i in 1:n1)
#' {
#' a = as.vector(exp(-X.m1.used[i, , drop=FALSE]%*%theta.m1))
#'   T.1 = T.1 + (a/(1+a)^2) * (t(X.m1.used[i, , drop=FALSE])%*%X.m1.used[i, , drop=FALSE])
#' }

#' for (i in 1:n1)
#' {
#'   a = as.vector(1/( 1 + exp(-X.m1.used[i, , drop=FALSE]%*%theta.m1)))
#'   T.2 = T.2 + (Y.m1[i]-a)^2 * (t(X.m1.used[i, , drop=FALSE])%*%X.m1.used[i, , drop=FALSE])
#' }

#' Sigma.m1 = solve(T.1)%*%T.2%*%solve(T.1) # This is actually Sigma.m1.n1.

#' ##### Find Sigma.m2

#' T.1 = matrix(rep(0, (length(A2)+1)^2), nrow=length(A2)+1)
#' T.2 = T.1

#' for (i in 1:n2)
#' {
#'   a = as.vector(exp(-X.m2.used[i, , drop=FALSE]%*%theta.m2))
#'   T.1 = T.1 + (a/(1+a)^2) * (t(X.m2.used[i, , drop=FALSE])%*%X.m2.used[i, , drop=FALSE])
#' }

#' for (i in 1:n2)
#' {
#'   a = as.vector(1/( 1 + exp(-X.m2.used[i, , drop=FALSE]%*%theta.m2)))
#'   T.2 = T.2 + (Y.m2[i]-a)^2 * (t(X.m2.used[i, , drop=FALSE])%*%X.m2.used[i, , drop=FALSE])
#' }

#' Sigma.m2 = solve(T.1)%*%T.2%*%solve(T.1)


#' ##### Find Sigma.m3

#' T.1 = matrix(rep(0, (length(A3)+1)^2), nrow=length(A3)+1)
#' T.2 = T.1

#' for (i in 1:n3)
#' {
#'   a = as.vector(exp(-X.m3.used[i, , drop=FALSE]%*%theta.m3))
#'   T.1 = T.1 + (a/(1+a)^2) * (t(X.m3.used[i, , drop=FALSE])%*%X.m3.used[i, , drop=FALSE])
#' }

#' for (i in 1:n3)
#' {
#'   a = as.vector(1/( 1 + exp(-X.m3.used[i, , drop=FALSE]%*%theta.m3)))
#'   T.2 = T.2 + (Y.m3[i]-a)^2 * (t(X.m3.used[i, , drop=FALSE])%*%X.m3.used[i, , drop=FALSE])
#' }

#' Sigma.m3 = solve(T.1)%*%T.2%*%solve(T.1)







#' names(theta.m1)=c("(Intercept)","Age","Height")
#' names(theta.m2)=c("(Intercept)","Height", "Weight")
#' names(theta.m3)=c("(Intercept)","Age", "Weight")

###now put in the GENMETA example

#' study1 = list(Coeff=theta.m1,Covariance=NULL,Sample_size=n1)
#' study2 = list(Coeff=theta.m2,Covariance=NULL,Sample_size=n2)
#' study3 = list(Coeff=theta.m3,Covariance=NULL,Sample_size=n3)

#' studies = list(study1,study2,study3)
#' model = "logistic"

#' reference = cbind(rep(1,n), X.rf)
#' colnames(reference) = c("(Intercept)","Age","Height", "Weight")
#' same.inter = GENMETA(studies, reference, model, initial_val = c(-1.2, log(1.3), log(1.3), log(1.3)))
#' diff.inter = GENMETA(studies, reference, model, variable_intercepts=TRUE)


#' @author Prosenjit Kundu, Runlong Tang and Nilanjan Chatterjee.
#' @importFrom magic adiag
#' @import MASS
#' @import stats
#' @import graphics
#' @importFrom Matrix rankMatrix
#' @export
#library(magic)
#library(MASS)

#Definition of GENMETA function

#GENMETA <- function(study_info, ref_dat, model, variable_intercepts=FALSE, control = list(epsilon = 1e-06, maxit = 1000))

GENMETA <- function(study_info, ref_dat, model, variable_intercepts=FALSE, initial_val=NULL, control = list(epsilon = 1e-06, maxit = 1000))
{
  #print("Computing in progress!")
  # optional_arguments <- list(...)
  # if(is.null(optional_arguments$control))
  # {
  #   control <- list("threshold" = 1e-06, "maxit" = 1000)
  #   threshold <- control$threshold
  #   maxit <- control$maxit
  # }
  # else{
  #   threshold <- optional_arguments$control[[1]]
  #   maxit <- optional_arguments$control[[2]]
  # }
 # if(is.null(optional_arguments$variable_intercepts))
 # {
 #   different_intercept <- FALSE
 # }
 #  else{
 #    different_intercept <- optional_arguments$variable_intercepts
 #  }
  # if(missing(control))
  # {
  #   control <- GENMETA.control()
  #   threshold <- control[[1]]
  #   maxit <- control[[2]]
  # }
  # else{
  #   threshold <- control[[1]]
  #   maxit <- control[[2]]
  # }
  call_MetaG <- match.call()
  ## Checking if the control argument is missing
  if(!missing(control))
  {
    control = control
  }
  threshold <- control[[1]]
  maxit <- control[[2]]

   different_intercept <- variable_intercepts
   study_estimates <- study_info

    error_1 <- 0
    error_2 <- 0
    error_3 <- 0
    error_4 <- 0
    error_5 <- 0
    no_of_studies <- length(study_estimates)
    #try(if(no_of_studies != length(study_estimates)) stop("number of studies does not match with length of the list study"))
    temp <- c()
    indicator_missing_covariance_sample_size = 0
    #indicator_missing_dispersion_study = 0
    missing_study_sample_size <- c()
    missing_covariance_study_indices <- c()
    #missing_dispersion_study_indices <- c()

    for(i in 1 : no_of_studies)
    {
        if(is.null(study_estimates[[i]][[3]]) == T)
        missing_study_sample_size = c(missing_study_sample_size, i)
    }
    # if(model == "linear")
    # {
    #   for(i in 1 : no_of_studies)
    #   {
    #     if(is.null(study_estimates[[i]][[4]]) == T)
    #       missing_dispersion_study_indices = c(missing_dispersion_study_indices, i)
    #   }
    # }


    for(i in 1:no_of_studies)
    {
        temp <- union(temp,names(study_estimates[[i]][[1]]))
    }
    for(i in 1:no_of_studies)
    {
        if(is.null(study_estimates[[i]][[2]]) == T && is.null(study_estimates[[i]][[3]]) == T)
        indicator_missing_covariance_sample_size = 1
    }
    # Condition needed for user input of dispersion parameters, for now we assume var(Y) = 1 and therefore estimate them. Later we can add them to increase flexibility.
    # for(i in 1:no_of_studies)
    # {
    #   if(study_estimates[[i]][[4]] == "NULL")
    #     indicator_missing_dispersion_study = 1
    # }

    ## All sanity checks...

    if(indicator_missing_covariance_sample_size == 1){
        print("Error: All the studies should have either an estimate for the var-cov matrix or the sample size. Atleast one of them is missing(NULL) in atleast one of the studies")
        error_1 <-1
    }
    if(indicator_missing_covariance_sample_size == 0)
    {
        for(i in missing_study_sample_size)
        {
            if(length(study_estimates[[i]][[1]]) != ncol(study_estimates[[i]][[2]]))
            {
                print("Error: length of the study parameter(effect sizes) vector does not match with the dimension of its variance covariance matrix")
                error_2 <- 1
            }
            if(sum(names(study_estimates[[i]][[1]]) != colnames(study_estimates[[i]][[2]])) > 0)
            {
                print("Error: names of the variables corresponding to the study specific parameter vector is not same(also not in the same order) as in its variance covariance matrix")
                error_2 <- 1
            }
        }
    }

    if(ncol(ref_dat) < length(temp))
    {
        print("Error: number of covariates in the reference data does not match with that of the maximal model")
        error_3 <- 1
    }

    if(sum(is.na(match(temp, colnames(ref_dat)))) > 0)
    {
        print("Error: names of covariates in the reference data does not match with that of the study specific covariates")
        error_4 <- 1
    }
    if(model == "linear" & variable_intercepts == "TRUE")
    {
      print("Error: when the model is linear, the current version works only when intercepts are assumed same across studies ")
      error_5 <- 1
    }
    ## Needed for calculating initial value
    names_wo_intercept <- c()
    for(i in 1:no_of_studies)
    {
      names_wo_intercept <- union(names_wo_intercept,names(study_estimates[[i]][[1]][-1]))
    }
    weight_sum <- 0
    sum <- 0
    estimates_in_which_studies_indices <- list()
    for(k in 1: length(names_wo_intercept))
    {
      temp_estimates_in_which <- c()
      for(j in 1: no_of_studies)
      {
        if(names_wo_intercept[k] %in% names(study_estimates[[j]][[1]]) == T)
          temp_estimates_in_which <- c(temp_estimates_in_which, j)
      }
      estimates_in_which_studies_indices[[k]] <- temp_estimates_in_which
    }
    ## End of the need for calculating initial value

    ## If the given inputs pass all the sanity checks...
    if(error_1 == 0 && error_2 == 0 && error_3 == 0 && error_4 == 0 && error_5 == 0)
    {
      
      if(length(initial_val) == 0)
      {
        initial_val <- c()
        # Calculating initial value when the different_intercept is TRUE
        if(different_intercept == TRUE)
        {
          initial_val <- c(initial_val, unlist(lapply(lapply(study_estimates, `[[`, 1), `[[`, 1)))
          
          for(k in 1: length(names_wo_intercept))
          {
            for(j in estimates_in_which_studies_indices[[k]])
            {
              if(is.null(study_estimates[[j]][[2]]) == F)
              {
                index_cov <- which(names(study_estimates[[j]][[1]]) %in% names_wo_intercept[k]  == T)
                weight <- 1/study_estimates[[j]][[2]][index_cov, index_cov]
                weight_sum <- weight_sum + weight
                sum <- sum + study_estimates[[j]][[1]][index_cov]*weight
              }
              if(is.null(study_estimates[[j]][[2]]) == T)
              {
                index_cov <- which(names(study_estimates[[j]][[1]]) %in% names_wo_intercept[k]  == T)
                weight <- study_estimates[[j]][[3]]
                weight_sum <- weight_sum + weight
                sum <- sum + study_estimates[[j]][[1]][index_cov]*weight
              }
            }
            initial_val[k+length(study_estimates)] <- sum/weight_sum
            weight_sum <- 0
            sum <- 0
          }
        }
        # End of Calculating initial value when the different_intercept is TRUE
        
        # Calculating initial value when the different_intercept is FALSE
        if(different_intercept == F)
        {
          for(k in 1: length(names_wo_intercept))
          {
            for(j in estimates_in_which_studies_indices[[k]])
            {
              if(is.null(study_estimates[[j]][[2]]) == F)
              {
                index_cov <- which(names(study_estimates[[j]][[1]]) %in% names_wo_intercept[k]  == T)
                weight <- 1/study_estimates[[j]][[2]][index_cov, index_cov]
                weight_sum <- weight_sum + weight
                sum <- sum + study_estimates[[j]][[1]][index_cov]*weight
              }
              if(is.null(study_estimates[[j]][[2]]) == T)
              {
                index_cov <- which(names(study_estimates[[j]][[1]]) %in% names_wo_intercept[k]  == T)
                weight <- study_estimates[[j]][[3]]
                weight_sum <- weight_sum + weight
                sum <- sum + study_estimates[[j]][[1]][index_cov]*weight
              }
            }
            initial_val[k+1] <- sum/weight_sum
            weight_sum <- 0
            sum <- 0
          }
          for(j in 1:no_of_studies)
          {
            if(is.null(study_estimates[[j]][[2]]) == F)
            {
              weight <- 1/study_estimates[[j]][[2]][1, 1]
              weight_sum <- weight_sum + weight
              sum <- sum + study_estimates[[j]][[1]][1]*weight
            }
            if(is.null(study_estimates[[j]][[2]]) == T)
            {
              weight <- study_estimates[[j]][[3]]
              weight_sum <- weight_sum + weight
              sum <- sum + study_estimates[[j]][[1]][1]*weight
            }
          }
          initial_val[1] <- sum/weight_sum
        }
        # End of Calculating initial value when the different_intercept is FALSE
        
        # End of initial_val calculation
      }
      


        ## Creating X_rbind and X_bdiag matrices
        for(i in 1 : no_of_studies)
        {
            if(is.null(study_estimates[[i]][[2]]) == T)
            missing_covariance_study_indices = c(missing_covariance_study_indices, i)
        }

        ## Defining X_rbind here
        X_rbind <- c()
        for(k in 1 : no_of_studies)
        {
            X_rbind <- rbind(X_rbind,ref_dat)
        }
        ## Defining X_abdiag  here
        X_bdiag_list <- list()
        for(k in 1 : no_of_studies)
        {
            col_ind <-  which(colnames(ref_dat) %in% names(study_estimates[[k]][[1]]) == TRUE)
            X_bdiag_list[[k]] <- as.matrix(ref_dat[,col_ind])
        }
        X_bdiag <- Reduce(magic::adiag,X_bdiag_list)
        ## End of Creating X_rbind and X_bdiag matrices

        model_optim <- model
        study_optim <- study_estimates
        eps_outer = 0


        ## Iterated GMM algorithm when the model is linear...
        if(model_optim == "linear")
        {
            # Calculating dispersion parameters by assuming the outcomes to have unit variance...
            study_indices <- seq(1,no_of_studies,1)
            if(length(missing_covariance_study_indices) > 0)
            non_missing_covariance_study_indices <- study_indices[-which(study_indices %in% missing_covariance_study_indices)]
            if(length(missing_covariance_study_indices) == 0)
            non_missing_covariance_study_indices = study_indices
            disp <- rep(NA, no_of_studies)
            for(j in study_indices)
            {
                col_ind <-  which(colnames(ref_dat) %in% names(study_estimates[[j]][[1]]) == TRUE)
                disp[j] <- (1- ((nrow(ref_dat)-1)/nrow(ref_dat)*var(ref_dat[,col_ind] %*% study_estimates[[j]][[1]])))
            }
            # End of Calculating dispersion parameters by assuming the outcomes to have unit variance

            C_init = diag(ncol(X_bdiag))
            Gamma_hat <- matrix(NA, nrow(C_init), ncol(ref_dat))
            lambda_ref <- list()
            k_gamma_hat <- 1

            # Calculating Gamma matrix...
            for(k in 1:no_of_studies)
            {
                col_ind <-  which(colnames(ref_dat) %in% names(study_estimates[[k]][[1]]) == TRUE)
                Gamma_hat[k_gamma_hat:(k_gamma_hat + length(col_ind) -1), ] <- (t(ref_dat[, col_ind]) %*% ref_dat)/(disp[[k]] * nrow(ref_dat))
                k_gamma_hat <- k_gamma_hat + length(col_ind)
            }
            # End of Calculating Gamma matrix

            # Calculating Lambda matrix for those studies where the var-cov matrices are provided by the user...
            for(k in non_missing_covariance_study_indices)
            {
                col_ind <-  which(colnames(ref_dat) %in% names(study_estimates[[k]][[1]]) == TRUE)
                temp_lambda_ref <- t(ref_dat[,col_ind]) %*% ref_dat[,col_ind]
                lambda_ref[[k]] <- (temp_lambda_ref %*% study_estimates[[k]][[2]] %*% temp_lambda_ref)/nrow(ref_dat)
            }
            A_n1 <- matrix(NA, nrow(C_init), ncol(ref_dat))
            B_n1 <- matrix(NA, nrow(C_init), 1)
            k_A = 1
            for(k in 1: no_of_studies)
            {
                col_ind <-  which(colnames(ref_dat) %in% names(study_estimates[[k]][[1]]) == TRUE)
                A_n1[k_A:(k_A + length(col_ind) - 1), ] <- (t(ref_dat[, col_ind]) %*% ref_dat) * (1/disp[k])
                B_n1[k_A:(k_A + length(col_ind) - 1), ] <- (t(ref_dat[, col_ind]) %*% ref_dat[, col_ind] %*% study_estimates[[k]][[1]]) * (1/disp[k])
                k_A <- k_A + length(col_ind)
            }
            beta_old_first <- solve(t(A_n1) %*% C_init %*% A_n1, tol = 1e-60)
            beta_old_sec <- t(A_n1) %*% C_init %*% B_n1
            beta_old_identity <- beta_old_first %*% beta_old_sec
            beta_old <- beta_old_first %*% beta_old_sec

            A_n2 <- 1/(disp^2)
            B_n2 <- rep(NA, no_of_studies)
            proceed <- TRUE
            U <- matrix(NA, nrow(C_init), nrow(ref_dat))
            no_of_iter <- 0
            while(proceed)
            {
                no_of_iter <- no_of_iter + 1
                disp_max_old <- (1- ((nrow(ref_dat)-1)/nrow(ref_dat)*var(ref_dat %*% beta_old)))
                #print(disp_max_old)
                for(k in missing_covariance_study_indices)
                {
                    col_ind <-  which(colnames(ref_dat) %in% names(study_estimates[[k]][[1]]) == TRUE)
                    temp_lambda_ref_1 <- as.numeric(disp_max_old) * t(ref_dat[, col_ind]) %*% ref_dat[, col_ind]
                    temp_lambda_W <- diag((as.vector(ref_dat %*% beta_old) - as.vector(ref_dat[, col_ind] %*% study_estimates[[k]][[1]])))
                    temp_lambda_ref_2 <- t(ref_dat[, col_ind]) %*% temp_lambda_W %*% ref_dat[, col_ind]
                    lambda_ref[[k]] <- (temp_lambda_ref_1 + temp_lambda_ref_2)/(study_estimates[[k]][[3]] * (disp[k]^2))
                }
                Lambda_ref <- Reduce(magic::adiag,lambda_ref)
                k_U <- 1
                for(k in 1:no_of_studies)
                {
                    col_ind <-  which(colnames(ref_dat) %in% names(study_estimates[[k]][[1]]) == TRUE)
                    temp_delta_k <- diag(as.vector(ref_dat %*% beta_old) - as.vector(ref_dat[, col_ind] %*% study_estimates[[k]][[1]]))/(disp[k]^2)
                    U[k_U:(k_U + length(col_ind) - 1), ] <- t(ref_dat[, col_ind]) %*% temp_delta_k
                    k_U <- k_U + length(col_ind)
                }
                Delta_hat <- (U %*% t(U))/(nrow(ref_dat))
                C_new <- solve(Lambda_ref + Delta_hat, tol = 1e-60)
                asy_var_C_identity <- (solve(t(Gamma_hat) %*% Gamma_hat, tol = 1e-60) %*% (t(Gamma_hat) %*% C_new %*% Gamma_hat) %*% solve(t(Gamma_hat) %*% Gamma_hat, tol = 1e-60))/(nrow(ref_dat))
                beta_new_first <- solve(t(A_n1) %*% C_new %*% A_n1, tol = 1e-60)
                beta_new_sec <- t(A_n1) %*% C_new %*% B_n1
                beta_new <- beta_new_first %*% beta_new_sec
                eps_outer = sqrt(sum((beta_new - beta_old)^2))
                beta_old <- beta_new
                if(eps_outer < threshold)
                proceed <- FALSE

            }
            asy_var_opt <- solve(t(Gamma_hat) %*% C_new %*% Gamma_hat, tol = 1e-60)/nrow(ref_dat)
            beta_old_identity <- as.vector(beta_old_identity)
            names(beta_old_identity) <-  colnames(ref_dat)
            if(is.null(asy_var_C_identity) == FALSE)
            {
              colnames(asy_var_C_identity) <- colnames(ref_dat)
              rownames(asy_var_C_identity) <- colnames(ref_dat)
            }
            if(is.null(asy_var_opt) == FALSE)
            {
              colnames(asy_var_opt) <- colnames(ref_dat)
              rownames(asy_var_opt) <- colnames(ref_dat)
            }
            
            beta_old <- as.vector(beta_old)
            names(beta_old) <- colnames(ref_dat)
            
            linear_result <- list("Est.coeff" = beta_old, "Est.var.cov" = asy_var_opt, "Res.var" = disp_max_old, "iter" = no_of_iter, "call" = call_MetaG)
            class(linear_result) <- "GENMETA"
            return(linear_result)

        }


        ## When the model is logistic, calls the myoptim function which implements the NR method...
        if(model_optim == "logistic")
        {
                C_init = diag(ncol(X_bdiag))
                no_of_iter_outer = 0
                total_iter = 0
                #print(initial_val)
                output_identity <- myoptim(no_of_studies, study_optim, ref_dat, X_rbind, X_bdiag_list, C_init, initial_val, threshold, model_optim, missing_covariance_study_indices, different_intercept, no_of_iter_outer)
                beta_identity <- output_identity$beta_optim
                beta_initial <- output_identity$beta_optim
                asy_var_beta_converged_identity <- output_identity$Asy_var_optim
                total_iter_identity <- output_identity$iter_IRWLS
                C_iter <- output_identity$C_optim
                proceed <- TRUE
                mark_failure = 0
                
                if(sum(is.na(beta_identity)) > 0 || output_identity$Status == 0)
                {
                  #print("Jacobian is computationally singular or the algo didn't converge")
                  beta_initial = rep(NA, ncol(ref_dat))
                  asy_var_beta_converged = NULL
                }else{
                  
                  while(proceed)
                  {
                    no_of_iter_outer <- no_of_iter_outer + 1
                    output_optim <- myoptim(no_of_studies, study_optim, ref_dat, X_rbind, X_bdiag_list, C_iter, beta_initial, threshold, model_optim, missing_covariance_study_indices, different_intercept, no_of_iter_outer)
                    beta_iter_old <- output_optim$beta_optim
                    C_iter <- output_optim$C_optim
                    if(sum(is.na(beta_iter_old)) > 0 || output_optim$Status == 0)
                    {
                      mark_failure = 1
                      break;
                    }
                    eps_outer = sqrt(sum((beta_iter_old - beta_initial)^2))
                    total_iter <- output_optim$iter_IRWLS + no_of_iter_outer
                    #print(eps_outer)
                    beta_initial <- beta_iter_old
                    if(eps_outer < threshold)
                      proceed <- FALSE
                  }
                  
                  total_iter <- total_iter + total_iter_identity
                  if(mark_failure == 1)
                  {
                    asy_var_beta_converged <- NULL
                  }else{
                    asy_var_beta_converged <- output_optim$Asy_var_optim
                    if(different_intercept == T)
                    {
                      temp_name <- c()
                      for(i in 1 : no_of_studies)
                      {
                        temp_name <- c(temp_name, paste0("(Intercept_Study_",i,")"))
                      }
                      temp_name <- c(temp_name, colnames(ref_dat[,-1]))
                      beta_identity <- as.vector(beta_identity)
                      #print(beta_initial)
                      names(beta_identity) <- temp_name
                      beta_initial <- as.vector(beta_initial)
                      names(beta_initial) <- temp_name
                      #print(beta_initial)
                      
                      if(is.null(asy_var_beta_converged_identity) == FALSE)
                      {
                        colnames(asy_var_beta_converged_identity) <- temp_name
                        rownames(asy_var_beta_converged_identity) <- temp_name
                      }
                      if(is.null(asy_var_beta_converged) == FALSE)
                      {
                        colnames(asy_var_beta_converged) <- temp_name
                        rownames(asy_var_beta_converged) <- temp_name
                      }
                      
                    }
                    
                    if(different_intercept == F)
                    {
                      beta_identity <- as.vector(beta_identity)
                      names(beta_identity) <- colnames(ref_dat)
                      beta_initial <- as.vector(beta_initial)
                      names(beta_initial) <- colnames(ref_dat)
                      if(is.null(asy_var_beta_converged_identity) == FALSE)
                      {
                        colnames(asy_var_beta_converged_identity) <- colnames(ref_dat)
                        rownames(asy_var_beta_converged_identity) <- colnames(ref_dat)
                      }
                      if(is.null(asy_var_beta_converged) == FALSE)
                      {
                        colnames(asy_var_beta_converged) <- colnames(ref_dat)
                        rownames(asy_var_beta_converged) <- colnames(ref_dat)
                      }
                      
                    }
                    
                  }
                  
                  
                }
                
                
               
               

                

                logistic_result <- list("Est.coeff" = beta_initial, "Est.var.cov" = asy_var_beta_converged, "Res.var" = NA, "iter" = total_iter, "call" = call_MetaG)
                class(logistic_result) <- "GENMETA"
                return(logistic_result)

        }

    }


}

