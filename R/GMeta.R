#' Implementing Generalized Meta-analysis
#'
#' Generalized Meta-analysis(GMeta) is an approach for combining information on multivariate regression parameters across multiple different studies which have different, but, possibly overlapping information on subsets of covariates.
#' GMeta implements the generalized meta-analysis using IRWLS algorithm.
#' @param study_info a list of lists containing information about the studies; the main list contains a list for each study, which must have the fields:
#' \itemize{
#' \item{"Coeff": a named numeric vector containing the estimates of regression parameters (including intercept) where the names identify the covariates. For example, names(study_info$Coeff) <- c("(Intercept)", "Age", "Height", "Weight").}
#' \item{"Covariance": a matrix containing an estimate of variance-covariance matrix of the regression parameters. This can be NULL if the "Sample_size" is provided.}
#' \item{"Sample_size": a numeric containing sample size of the study. This can be NULL if the "Covariance" is provided.}
#' }
#' @param ref_dat a data matrix containing all the distinct covariates across studies from a reference set of individuals. This is used for estimating joint distribution of the covariates. The data matrix must have the vector of ones as its first column. The column names of the data matrix should match the names of the covariates from the studies.
#' @param model a description of the type of regression model; this is a character string naming the regression model. The current version is for "logistic" and "linear".
#' @param variable_intercepts an optional logical (applicable only when the model is "logistic"); if TRUE, the intercepts of the true models for each of the studies are assumed to be different. Default is FALSE.
#' @param control an optional list containing the epsilon (numeric) and maxiter (numeric) needed for convergence of the algorithm. Default epsilon and maximum iterations are 1e-06 and 1000, respectively. For creating a control argument for GMeta, see \code{\link[GMeta]{GMeta.control}}.
#' @details Generalized Meta-analysis (GMeta) is a tool that allows researchers to quickly build models for multivariate meta-analysis in the presence of disparate covariate information across studies. It is implemented based on mainly two input arguments:
#' \itemize{
#' \item{Information on the model parameters from each of the studies.}
#' \item{Reference data for estimation of the joint distribution of all the distinct covariates across studies.}}
#' The software provides flexibility to the users to choose the intercepts to be different (when the model is logistic) across studies through the input argument, variable_intercepts.
#' It also allows estimation of the regression parameters, only from the sample sizes of the studies when it is difficult to obtain estimate of the variance-covariance matrices.
#' \cr \cr \bold{Note}: GMeta will not work if both the estimates of the covariance matrix and the sample size are NULL.
#' @details When the model is "linear", it is assumed that the outcome is standardized to have unit variance.
#' For more details on the IRWLS, see References.
#' @return An object of class "GMeta" is a list containing GMeta estimate, its variance-covariance matrix and estimates the residual variance in the case of "linear" model .
#' \item{Est.coeff}{a numeric vector containing the estimated regression coefficients of the maximal model using optimal weighting matrix.}
#' \item{Est.var.cov}{a matrix containing estimate of variance-covariance matrix of the corresponding GMeta estimator.}
#' \item{Res.var}{a numeric containing the residual variance of the maximal model when it is linear. It is calculated from the formula : \eqn{1 - \hat{\beta}_{GMeta}^Tvar(X)\hat{\beta}_{GMeta}} which is derived by assuming the outcomes to have unit variance. \eqn{var(X)} is calculated from reference data. Res.var is NA when the model is "logistic".}
#' \item{iter}{a numeric containing the number of iterations used in the algorithm}
#' \item{call}{the matched call}
#' \cr The function \code{\link[GMeta]{GMeta.summary}} prints a summary of the results obtained from GMeta.
#' \cr The function \code{\link[GMeta]{GMeta.plot}} plots the estimate of the parameter from each of the studies, the summary measure(GMeta estimate) and their confidence intervals. 
#' @keywords Generalized Meta Analysis
#' @references Tang, R., Kundu, P. and Chatterjee, N. (2017) Generalized Meta-Analysis for Multivariate Regression Models Across Studies with Disparate Covariate Information. \href{https://arxiv.org/abs/1708.03818}{arXiv:1708.03818v1 [stat.ME]}.
#' @seealso \code{\link[GMeta]{GMeta.summary}}, \code{\link[GMeta]{study_info}}, \code{\link[GMeta]{reference_data}}, \code{\link[GMeta]{GMeta.plot}}.
#' @examples
#' # This example shows the GMeta implementation on a simulated data set for logistic regression
#' data(reference_data)
#' head(reference_data)
#' data(study_info)
#' head(study_info)
#' model <- "logistic"
#' # When the true intercept parameters of the studies are different.
#' result_diff <- GMeta(study_info, reference_data, model, variable_intercepts = TRUE)
#' print(result_diff)
#' # When the true intercept parameters of the studies are same.
#' result_same <- GMeta(study_info, reference_data, model)
#' print(result_same)

#' @author Prosenjit Kundu, Runlong Tang and Nilanjan Chatterjee.
#' @import magic
#' @import MASS
#' @import stats
#' @import graphics
#' @export
#library(magic)
#library(MASS)

#Definition of Gmeta function

#GMeta <- function(study_info, ref_dat, model, variable_intercepts=FALSE, control = list(epsilon = 1e-06, maxit = 1000))

GMeta <- function(study_info, ref_dat, model, variable_intercepts=FALSE, control = list(...))
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
  #   control <- GMeta.control()
  #   threshold <- control[[1]]
  #   maxit <- control[[2]]
  # }
  # else{
  #   threshold <- control[[1]]
  #   maxit <- control[[2]]
  # }
  call_gmeta <- match.call()
  ## Checking if the control argument is missing
  if(missing(control))
  {
    control <- GMeta.control()
  }else{
    control <- do.call("GMeta.control", control)
  }
  threshold <- control[[1]]
  maxit <- control[[2]]

   different_intercept <- variable_intercepts
   study_estimates <- study_info

    error_1 <- 0
    error_2 <- 0
    error_3 <- 0
    error_4 <- 0
    #error_5 <- 0
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
        print("number of covariates in the reference data does not match with that of the full model")
        error_3 <- 1
    }

    if(sum(temp != colnames(ref_dat)) > 0)
    {
        print("names of covariates in the reference data does not match with that of the study specific covariates")
        error_4 <- 1
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
    if(error_1 == 0 && error_2 == 0 && error_3 == 0 && error_4 == 0)
    {
      # Calculating initial value when the different_intercept is TRUE
      initial_val <- c()
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


        ## Creating X_rbind and X_bdiag matrices
        for(i in 1 : no_of_studies)
        {
            if(is.null(study_estimates[[i]][[2]]) == T)
            missing_covariance_study_indices = c(missing_covariance_study_indices, i)
        }


        X_rbind <- c()
        for(k in 1 : no_of_studies)
        {
            X_rbind <- rbind(X_rbind,ref_dat)
        }
        X_bdiag_list <- list()
        for(k in 1 : no_of_studies)
        {
            col_ind <-  which(colnames(ref_dat) %in% names(study_estimates[[k]][[1]]) == TRUE)
            X_bdiag_list[[k]] <- ref_dat[,col_ind]
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
            non_missing_covariance_study_indices <- study_indices[-which(study_indices %in% missing_covariance_study_indices)]
            disp <- rep(NA, no_of_studies)
            for(j in study_indices)
            {
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
                C_new <- solve(Lambda_ref + Delta_hat, tol = 1e-06)
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
            colnames(asy_var_C_identity) <- colnames(ref_dat)
            rownames(asy_var_C_identity) <- colnames(ref_dat)
            beta_old <- as.vector(beta_old)
            names(beta_old) <- colnames(ref_dat)
            colnames(asy_var_opt) <- colnames(ref_dat)
            rownames(asy_var_opt) <- colnames(ref_dat)
            linear_result <- list("Est.coeff" = beta_old, "Est.var.cov" = asy_var_opt, "Res.var" = disp_max_old, "iter" = no_of_iter, "call" = call_gmeta)
            class(linear_result) <- "GMeta"
            return(linear_result)

        }


        ## When the model is logistic, calls the myoptim function which implements the NR method...
        if(model_optim == "logistic")
        {
                C_init = diag(ncol(X_bdiag))
                no_of_iter_outer = 0
                total_iter = 0
                output_identity <- myoptim(no_of_studies, study_optim, ref_dat, X_rbind, X_bdiag_list, C_init, initial_val, threshold, model_optim, missing_covariance_study_indices, different_intercept, no_of_iter_outer)
                beta_identity <- output_identity$beta_optim
                beta_initial <- output_identity$beta_optim
                asy_var_beta_converged_identity <- output_identity$Asy_var_optim
                total_iter_identity <- output_identity$iter_IRWLS
                C_iter <- output_identity$C_optim
                proceed <- TRUE
                while (proceed)
                {
                    #Define C with beta_initial
                    no_of_iter_outer <- no_of_iter_outer + 1
                    if(sum(is.na(beta_initial)) > 0)
                    {
                        break;
                    }
                    output_optim <- myoptim(no_of_studies, study_optim, ref_dat, X_rbind, X_bdiag, C_iter, beta_initial, threshold, model_optim, missing_covariance_study_indices, different_intercept, no_of_iter_outer)
                    beta_iter_old <- output_optim$beta_optim
                    C_iter <- output_optim$C_optim
                    eps_outer = sqrt(sum((beta_iter_old - beta_initial)^2))
                    total_iter <- output_optim$iter_IRWLS + no_of_iter_outer
                    #print(eps_outer)
                    beta_initial <- beta_iter_old
                    if(eps_outer < threshold)
                    proceed <- FALSE
                }
                total_iter <- total_iter + total_iter_identity
                if(sum(is.na(beta_initial)) > 0)
                {
                    asy_var_beta_converged <- NULL
                }else{
                    asy_var_beta_converged <- output_optim$Asy_var_optim
                }

                if(different_intercept == T)
                {
                  temp_name <- c()
                  for(i in 1 : no_of_studies)
                  {
                    temp_name <- c(temp_name, paste0("(Intercept_Study_",i,")"))
                  }
                  temp_name <- c(temp_name, colnames(ref_dat[,-1]))
                  beta_identity <- as.vector(beta_identity)
                  names(beta_identity) <- temp_name
                  beta_initial <- as.vector(beta_initial)
                  names(beta_initial) <- temp_name
                  colnames(asy_var_beta_converged_identity) <- temp_name
                  rownames(asy_var_beta_converged_identity) <- temp_name
                  colnames(asy_var_beta_converged) <- temp_name
                  rownames(asy_var_beta_converged) <- temp_name
                }

                if(different_intercept == F)
                {
                  beta_identity <- as.vector(beta_identity)
                  names(beta_identity) <- colnames(ref_dat)
                  beta_initial <- as.vector(beta_initial)
                  names(beta_initial) <- colnames(ref_dat)
                  colnames(asy_var_beta_converged_identity) <- colnames(ref_dat)
                  rownames(asy_var_beta_converged_identity) <- colnames(ref_dat)
                  colnames(asy_var_beta_converged) <- colnames(ref_dat)
                  rownames(asy_var_beta_converged) <- colnames(ref_dat)
                }

                logistic_result <- list("Est.coeff" = beta_initial, "Est.var.cov" = asy_var_beta_converged, "Res.var" = NA, "iter" = total_iter, "call" = call_gmeta)
                class(logistic_result) <- "GMeta"
                return(logistic_result)

        }

    }


}

