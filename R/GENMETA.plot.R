#' Generalized Meta-analysis(forest plot)
#'
#' This function plots the confidence intervals with boxes as the study specific estimates and diamond as the GENMETA estimate. For the current version, it assumes that the estimate of the variance-covariance matrix in each of the studies is provided.
#' It is demonstrated using a different dataset, "study_info_plot", which meets the assumption.
#' @param x an object of class "GENMETA"
#' @param study_info_plot a list of lists containing information about the studies(similar to the study_info argument used in GENMETA function.)
#' @import graphics
#' @examples
#' # This example shows how to obtain the forest plot of GENMETA object.
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


#' beta.star = matrix(c(-1.2, log(1.3), log(1.3), log(1.3)),nrow = d.X+1) 
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

#' study1 = list(Coeff=theta.m1,Covariance=Sigma.m1,Sample_size=n1)
#' study2 = list(Coeff=theta.m2,Covariance=Sigma.m2,Sample_size=n2)
#' study3 = list(Coeff=theta.m3,Covariance=Sigma.m3,Sample_size=n3)

#' studies = list(study1,study2,study3)
#' model = "logistic"

#' reference = cbind(rep(1,n), X.rf)
#' colnames(reference) = c("(Intercept)","Age","Height", "Weight")
#' result_diff <- GENMETA(studies, reference, model, variable_intercepts = TRUE)
#' GENMETA.plot(result_diff, studies)
#' @export

GENMETA.plot <- function(x, study_info_plot)
{
  no_of_studies <- length(study_info_plot)
  row_names <- c()
  tot_var_names <- c()
  for(i in 1: no_of_studies)
  {
    tot_var_names <- union(tot_var_names, names(study_info_plot[[i]][[1]]))
  }
  tot_var_names_wo_intercept <- tot_var_names[-1]

  y_indices <- seq((no_of_studies+1), 1, -1)

  for(i in 1:no_of_studies)
  {
    row_names <- c(row_names, paste0("Study",i))
  }
  row_names <- c(row_names, "GENMETA")


  #op <- par(ask=TRUE)
  for(i in 1:length(tot_var_names_wo_intercept))
  {

    data_plot <- matrix(NA, (no_of_studies+1), 4)
    data_plot[,1] <- y_indices

    for(j in 1: no_of_studies)
    {
      if(tot_var_names_wo_intercept[i] %in% names(study_info_plot[[j]][[1]]) == T)
      {
        index <- which(names(study_info_plot[[j]][[1]]) %in% tot_var_names_wo_intercept[i] == T)
        data_plot[j,2] <- as.vector(study_info_plot[[j]][[1]])[index]
        data_plot[j,3] <- as.vector(study_info_plot[[j]][[1]])[index] - 1.96*sqrt(as.vector(diag(study_info_plot[[j]][[2]]))[index])
        data_plot[j,4] <- as.vector(study_info_plot[[j]][[1]])[index] + 1.96*sqrt(as.vector(diag(study_info_plot[[j]][[2]]))[index])
      }

    }
    index <- which(names(x[[1]]) %in% tot_var_names_wo_intercept[i] == T)
    data_plot[(no_of_studies+1),2] <- as.vector(x[[1]])[index]
    data_plot[(no_of_studies+1),3] <- as.vector(x[[1]])[index] - 1.96*sqrt(as.vector(diag(x[[2]])))[index]
    data_plot[(no_of_studies+1),4] <- as.vector(x[[1]])[index] + 1.96*sqrt(as.vector(diag(x[[2]])))[index]



    data_plot <- as.data.frame(data_plot)
    colnames(data_plot) <- c("Y", "coef","low","high")
    rownames(data_plot) <- row_names
    #row_names <- list(list(row_names))
    #Enter the forestplot command
    #forestplot(row_names, data_plot$coef, data_plot$low, data_plot$high, zero = 1, cex  = 2, lineheight = "auto", xlab = "log odds ratio")
    #plot.new()

    plot(data_plot[,2], data_plot[,1], axes=F, xlab = "", ylab = "", pch = 15, xlim = c(-1.5,1.5))
    for(k in 1:no_of_studies)
    {
      segments(data_plot[k,3], data_plot[k,1], x1=data_plot[k,4], y1=data_plot[k,1], col = "black")
      #points(data_plot[k,2], data_plot[k,1], pch = 18, cex=2)
    }
    segments(data_plot[(no_of_studies+1),3], data_plot[(no_of_studies+1),1], x1=data_plot[(no_of_studies+1),4], y1=data_plot[(no_of_studies+1),1], col = "black")
    points(data_plot[(no_of_studies+1),2], data_plot[(no_of_studies+1),1], pch = 18, cex=2)
    mtext("Log Odds Ratio", side=1, las=0, col="black", outer = F, line = 2.5)
    abline(v=0, lty = 3)
    axis(1, at = seq(-2,2,0.5), labels = seq(-2,2,0.5))
    axis(2, at = data_plot[,1], labels = rownames(data_plot), las = 2)
    title(main = tot_var_names_wo_intercept[i])

    #legend("topright", legend = c("Study specific estimates", "GENMETA estimate", "95% CI"), lty = c(NA,NA,1), pch = c(15,18,NA))
    readline(prompt="Press [enter] for the forestplot of next covariate:")
  }
}
