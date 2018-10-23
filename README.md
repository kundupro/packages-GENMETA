# GENMETA : An R package
Generalized Meta-Analysis(GENMETA) is an approach for combining information on multivariate regression parameters across multiple different studies which have different, but, possibly overlapping information on subsets of covariates. GENMETA implements the generalized meta-analysis using iteratively reweighted least squares (IRWLS) algorithm. For details of the method, please see the Reference section below.  
This file provides guidelines for implementing generalized meta-analysis(GENMETA)  


NOTE: If the user wish to directly download the source files from Github, please go to the link https://github.com/28pro92/GENMETA and see the README for further instructions.   
If the user wish to directly install the package from R console, please read the following instructions:  

# Installation procedure:  
install.packages("devtools", dependencies=TRUE)  
library(devtools)  

install_github("28pro92/packages-GENMETA")  
library(GENMETA)  


# Example describing how to create input aruments for GENMETA
#--------- The following code is for the simulation results in Table 3 shown in Reference----------------#  


# Setting-I (Ideal setting)  
library(MASS)  
library(stats)  
library(magic)  
#--#######################  
#--### Basic setting######  
#--#######################  

#---number of covariates in the full model---#  
d.X = 3  

#---mean vector of the covariates---#  
mu = matrix(rep(0,d.X), nrow=d.X)   

#--- correlation coefficient of the covariates---#  
r1 = 0.3   
r2 = 0.6  
r3 = 0.1  

#--- covariance matrix of the covariates---#  
Sigma = matrix(   
  c(1, r1, r2,    
    r1, 1, r3,   
    r2, r3, 1),   
  nrow=d.X, ncol=d.X)  


#-- True parameter---#  
beta.star = matrix(c(-1.2, log(1.3), log(1.3), log(1.3)),nrow = d.X+1)   


#-----Sample sizes of the three studies-----#  
n1 = 300 # sample size of the 1st data set.  
n2 = 500 # 2nd  
n3 = 1000 # 3rd  

#---- Sample size of the reference dataset---#  
n = 50  


no.of.simulations = 1000  

#----- Defining a matrix to store the results for each simulation----#  
sim.matrix = matrix(NA, 1000, 8)  

start.time = Sys.time()  

#------ Repeat for each simulation ----#  

for(sim in 1:no.of.simulations)  
{  
  set.seed(sim)  


  #--################################################  
  #--##### Generating the reference and studies #####    
  #--################################################  

  #---Generate the reference data set---#  

  X.rf = mvrnorm(n = n, mu, Sigma)  
  
  
  #----Generate data set 1. m1 means model 1.---#  

  X.m1 = mvrnorm(n = n1, mu, Sigma) # Generate the covariates.  
  X.m1.1 = cbind(rep(1, n1), X.m1) # Add a column of 1's to X.m1.  
  p.m1 = 1/(1+exp(-X.m1.1%\*%beta.star)) # the vector of probabilities  
  Y.m1 = rbinom(n1, size=1, p.m1) # the Bernoulli responses  
  
  
  #----Generate data set 2. m2 means model 2.---#  

  X.m2 = mvrnorm(n = n2, mu, Sigma)  
  X.m2.1 = cbind(rep(1, n2), X.m2)  
  p.m2 = 1/(1+exp(-X.m2.1%\*%beta.star))  
  Y.m2 = rbinom(n2, size=1, p.m2)  
  
  #----Generate data set 3. m3 means model 3.---#  

  X.m3 = mvrnorm(n = n3, mu, Sigma)  
  X.m3.1 = cbind(rep(1, n3), X.m3)  
  p.m3 = 1/(1+exp(-X.m3.1%\*%beta.star))  
  Y.m3 = rbinom(n3, size=1, p.m3)  
  
  #--#######################################################  
  #--### Create data sets in the format of data frame.######  
  #--#######################################################  
  data.m1 = data.frame(Y=Y.m1, X.m1)  
  data.m2 = data.frame(Y=Y.m2, X.m2)  
  data.m3 = data.frame(Y=Y.m3, X.m3)  
  
  
  
  
  #--####################################################################  
  #--### Fit logistic regression with reduced models to the data sets ###  
  #--####################################################################  

  logit.m1 <- glm(Y ~ X1 + X2, data = data.m1, family = "binomial")  
  if(logit.m1$converged == FALSE)  
  {  
    print("glm for logit.m1 is not convergent.")  
    next  
  }  
  
  logit.m2 <- glm(Y ~ X2 + X3, data = data.m2, family = "binomial")  

  if(logit.m2$converged == FALSE)  
  {  
    print("glm for logit.m2 is not convergent.")  
    next  
  }  
  
  logit.m3 <- glm(Y ~ X1 + X3, data = data.m3, family = "binomial")  
  
  if(logit.m3$converged == FALSE)  
  {  
    print("glm for logit.m3 is not convergent.")  
    next  
 }  
  
  
  #--####################################################################  
  #--### Obtain the estimators of the parameters in the reduced models.##   
  #--####################################################################  

  theta.m1 = logit.m1$coefficients  
  theta.m2 = logit.m2$coefficients  
  theta.m3 = logit.m3$coefficients  
  
  
  #--####################################################################  
  #--#### The following is used to calculate robust variance estimates###  
  #--#### for each of three parameter vectors in the reduced models. ####  
  #--#### This is needed if the user inputs the variance-covariance #####  
  #--### matrices from each of studies. Later, we show that it can be ##   
  #--### computed from reference data set if these are not provided ####    
  #--###################################################################  
  
  #--####################################################################  
  #--### Find the covariance matrix estimators for the reduced models####  
  #--####################################################################  
  
  
  #--#############################  
  #-- Basic notations for inputs#  
  #--#############################  
  
#-- Number of data sets---#   

K = 3   
  
#-- index set A1, the indexes of the covariates of data set 1.--#  

  A1 = c(1, 2)  
  A2 = c(2, 3) # index set A2  
  A3 = c(1, 3) # index set A3  
  
#--- X.m1.used denotes the design matrix in model 1----#  
  
  X.m1.used = cbind(rep(1, n1), X.m1\[, A1, drop=F\])  
  X.m2.used = cbind(rep(1, n2), X.m2\[, A2, drop=F\])  
  X.m3.used = cbind(rep(1, n3), X.m3\[, A3, drop=F\])  
  
  
  

#----- Computing robust estimate of Sigma.m1 -------#  
  
  T.1 = matrix(rep(0, (length(A1)+1)^2), nrow=length(A1)+1)  
  T.2 = T.1  
  
  for (i in 1:n1)  
  {  
    a = as.vector(exp(-X.m1.used\[i, , drop=F\]%\*%theta.m1))  
    T.1 = T.1 + (a/(1+a)^2) \* (t(X.m1.used\[i, , drop=F\]) %*% X.m1.used\[i, , drop=F\])  
  }  
  
  for (i in 1:n1)  
  {  
    a = as.vector(1/( 1 + exp(-X.m1.used\[i, , drop=F\]%\*%theta.m1)))  
    T.2 = T.2 + (Y.m1\[i\]-a)^2 \* (t(X.m1.used\[i, , drop=F\])%\*%X.m1.used\[i, , drop=F\])  
  }  
  
  Sigma.m1 = solve(T.1)%\*%T.2%\*%solve(T.1)   
  
  
#----- Computing robust estimate of Sigma.m2 -------#  

  
  T.1 = matrix(rep(0, (length(A2)+1)^2), nrow=length(A2)+1)  
  T.2 = T.1  
  
  for (i in 1:n2)  
  {  
    a = as.vector(exp(-X.m2.used\[i, , drop=F\]%\*%theta.m2))  
    T.1 = T.1 + (a/(1+a)^2) \* (t(X.m2.used\[i, , drop=F\])%\*%X.m2.used\[i, , drop=F\])  
  }  
  
  for (i in 1:n2)  
  {  
    a = as.vector(1/( 1 + exp(-X.m2.used\[i, , drop=F\]%\*%theta.m2)))  
    T.2 = T.2 + (Y.m2\[i\]-a)^2 \* (t(X.m2.used\[i, , drop=F\])%\*%X.m2.used\[i, , drop=F\])  
  }  
  
  Sigma.m2 = solve(T.1)%\*%T.2%\*%solve(T.1)  
  
  
 



 #----- Computing robust estimate of Sigma.m3 -------#   

  
  T.1 = matrix(rep(0, (length(A3)+1)^2), nrow=length(A3)+1)  
  T.2 = T.1  
  
  for (i in 1:n3)  
  {  
    a = as.vector(exp(-X.m3.used\[i, , drop=F\]%\*%theta.m3))  
    T.1 = T.1 + (a/(1+a)^2) \* (t(X.m3.used\[i, , drop=F\])%\*%X.m3.used\[i, , drop=F\])  
  }  
  
  for (i in 1:n3)  
  {  
    a = as.vector(1/( 1 + exp(-X.m3.used\[i, , drop=F\]%\*%theta.m3)))  
    T.2 = T.2 + (Y.m3\[i\]-a)^2 \* (t(X.m3.used\[i, , drop=F\])%\*%X.m3.used\[i, , drop=F\])  
  }  
  
  Sigma.m3 = solve(T.1)%\*%T.2%\*%solve(T.1)  
  
  
  
  
  names(theta.m1)=c("(Intercept)","Age","Height")  
  names(theta.m2)=c("(Intercept)","Height", "Weight")  
  names(theta.m3)=c("(Intercept)","Age", "Weight")  
  
  #---- now put in the GENMETA example ----#  

 #--### The results shown in Table 2 is obtained by considering #####  
 #--#### study-specific variance-covariance matrices. This is shown ##  
 #--#### below. ######################################################  
 #--##################################################################   
  study1 = list(Coeff=theta.m1,Covariance=Sigma.m1,Sample_size=n1)  
  study2 = list(Coeff=theta.m2,Covariance=Sigma.m2,Sample_size=n2)  
  study3 = list(Coeff=theta.m3,Covariance=Sigma.m3,Sample_size=n3)  

#--##### If the study-specific variance-covariance matrices are not ###  
#--##### available, then it is set to NULL. The following lines #######  
#--##### demonstrate it…                                       #######  
#--##### study1 = list(Coeff=theta.m1,Covariance=NULL,Sample_size=n1)##  
#--##### study2 = list(Coeff=theta.m2,Covariance=NULL,Sample_size=n2)##  
#--##### study3 = list(Coeff=theta.m3,Covariance=NULL,Sample_size=n3)##  

 #---- Creating the study list for GENMETA input ---#  
  
 studies = list(study1,study2,study3)  
  model = "logistic"  
  
  #------ Creating the reference data for GENMETA input -----#  

  reference = cbind(rep(1,n), X.rf)  
  colnames(reference) = c("(Intercept)","Age","Height", "Weight")  

#-------Calling the GENMETA function-------#  

  result.same = GENMETA(studies, reference, model, initial_val = c(-1.2, log(1.3), log(1.3), log(1.3)))  
  
  
#------Checking the conditions where the optim did not converge--#   

 if(sum(is.na(result.same$Est.coeff)) == 0 && is.null(result.same$Est.var.cov) != TRUE)  
  {  
    sim.matrix\[sim, ] = c(result.same$Est.coeff, diag(result.same$Est.var.cov))  
  }  
  
  print(sim)  
}  

stop.time = Sys.time()  

elapsed.time = stop.time - start.time  


#-----Computing the 95% confidence-intervals---#  
Lower.CI = na.omit(sim.matrix)\[,1:4\] - 1.96\*sqrt(na.omit(sim.matrix)\[,5:8\])  
Upper.CI = na.omit(sim.matrix)\[,1:4\] + 1.96\*sqrt(na.omit(sim.matrix)\[,5:8\])  

#------Computing the coverage rate---------#    

Coverage.rate = matrix(0,dim(na.omit(sim.matrix))\[1\], 4)  
for(k in 1:dim(na.omit(sim.matrix))\[1\])  
{  
  if(beta.star\[1,1\] >= Lower.CI\[k,1\] & beta.star\[1,1\] <= Upper.CI\[k,1\])  
    Coverage.rate\[k,1\] = 1  
  if(beta.star\[2,1\] >= Lower.CI\[k,2\] & beta.star\[2,1\] <= Upper.CI\[k,2\])  
    Coverage.rate\[k,2\] = 1  
  if(beta.star\[3,1\] >= Lower.CI\[k,3\] & beta.star\[3,1\] <= Upper.CI\[k,3\])  
    Coverage.rate\[k,3\] = 1  
  if(beta.star\[4,1\] >= Lower.CI\[k,4\] & beta.star\[4,1\] <= Upper.CI\[k,4\])  
    Coverage.rate\[k,4\] = 1  
}  

square = function(x)  
{  
  sum(x^2)  
}  

#----Computing the average length of CI----#  
  
Average.length = apply((Upper.CI-Lower.CI), 2, mean)  

#------Computing root mean square error -----#  

RMSE = sqrt(apply(t(t(na.omit(sim.matrix)\[,1:4\]) - as.numeric(beta.star)), 2, square)/dim(na.omit(sim.matrix))\[1\])  

#-----Combining the results------#  

result = data.frame(cbind(apply(na.omit(sim.matrix), 2, mean)\[1:4\], beta.star, (apply(na.omit(sim.matrix), 2, mean)\[1:4\]- beta.star), sqrt(apply(na.omit(sim.matrix), 2, mean)\[5:8\]), sqrt(apply(na.omit(sim.matrix)\[,1:4\], 2, var))), apply(Coverage.rate,2, mean), Average.length, RMSE)  

colnames(result) = c("Coeff","True", "Bias", "ESD", "SD", "Coverage.rate", "AL", "RMSE")  

#-------Writing into a file-----#  
write.csv(result, file = "Simulation_1.csv", row.names = F)  












# Setting-II : Means are different across studies
#--############################################################################  
#--######## For this and all other settings, the codes in the basic setting and for generating ####  
#--######## the data sets change. Rest remain the same. ##############################  
#--############################################################################  


#--#######################  
#--### Basic setting######   
#--#######################  

#-- number of covariates.  
d.X = 3   

#-- Different mean vectors of the covariates.  
mu.1 = matrix(rep(0,d.X), nrow=d.X)  
mu.2 = matrix(rep(1,d.X), nrow=d.X)   
mu.3 = matrix(rep(0.5,d.X), nrow=d.X)   

#-- correlation coefficient(low) of the covariates.  
r1.1 = 0.2   
r2.1 = 0.4  
r3.1 = 0  

#-- covariance matrix of the covariates.  
Sigma.1 = matrix(   
  c(1, r1.1, r2.1,    
    r1.1, 1, r3.1,   
    r2.1, r3.1, 1),   
  nrow=d.X,   
  ncol=d.X)    


#-- correlation coefficient(base/medium) of the covariates.  
r1.2 = 0.3    
r2.2 = 0.6  
r3.2 = 0.1  

#-- covariance matrix of the covariates.  
Sigma.2 = matrix(   
  c(1, r1.2, r2.2,    
    r1.2, 1, r3.2,   
    r2.2, r3.2, 1),   
  nrow=d.X,   
  ncol=d.X)  

#-- correlation coefficient(high) of the covariates.  
r1.3 = 0.4   
r2.3 = 0.8  
r3.3 = 0.2  

#-- covariance matrix of the covariates.  
Sigma.3 = matrix(   
  c(1, r1.3, r2.3,    
    r1.3, 1, r3.3,   
    r2.3, r3.3, 1),   
  nrow=d.X,   
  ncol=d.X) # covariance matrix of the covariates.  

#--- True parameter  
beta.star = matrix(c(-1.2, log(1.3), log(1.3), log(1.3)),nrow = d.X+1)   


#-----Sample sizes of the three studies-----#  
n1 = 300 # sample size of the 1st data set.  
n2 = 500 # 2nd  
n3 = 1000 # 3rd  

#---- Sample size of the reference dataset---#  
n = 50  


no.of.simulations = 1000  

#----- Defining a matrix to store the results for each simulation----#  
sim.matrix = matrix(NA, 1000, 8)  

start.time = Sys.time()  

#------ Repeat for each simulation ----#  

for(sim in 1:no.of.simulations)  
{  
  set.seed(sim)  


 #--#################################################  
 #--###### Generating the reference and studies #####  
 #--#################################################  

  #---Generate the reference data set---#  

  X.rf = mvrnorm(n = n, mu.1, Sigma.2)  
  
  
  #----Generate data set 1. m1 means model 1.---#  

  X.m1 = mvrnorm(n = n1, mu.1, Sigma.2) # Generate the covariates.  
  X.m1.1 = cbind(rep(1, n1), X.m1) # Add a column of 1's to X.m1.  
  p.m1 = 1/(1+exp(-X.m1.1%\*%beta.star)) # the vector of probabilities  
  Y.m1 = rbinom(n1, size=1, p.m1) # the Bernoulli responses  
 
  
  #----Generate data set 2. m2 means model 2.---#  

  X.m2 = mvrnorm(n = n2, mu.2, Sigma.2)  
  X.m2.1 = cbind(rep(1, n2), X.m2)  
  p.m2 = 1/(1+exp(-X.m2.1%\*%beta.star))  
  Y.m2 = rbinom(n2, size=1, p.m2)  
  
  #----Generate data set 3. m3 means model 3.---#  

  X.m3 = mvrnorm(n = n3, mu.3, Sigma.2)  
  X.m3.1 = cbind(rep(1, n3), X.m3)  
  p.m3 = 1/(1+exp(-X.m3.1%\*%beta.star))  
  Y.m3 = rbinom(n3, size=1, p.m3)  
  
  #--#######################################################  
  #--### Create data sets in the format of data frame.######  
  #--#######################################################  
  data.m1 = data.frame(Y=Y.m1, X.m1)  
  data.m2 = data.frame(Y=Y.m2, X.m2)  
  data.m3 = data.frame(Y=Y.m3, X.m3)  
 
  
  
  
  #--####################################################################  
  #--### Fit logistic regression with reduced models to the data sets ###  
  #--####################################################################  

  logit.m1 <- glm(Y ~ X1 + X2, data = data.m1, family = "binomial")  
  if(logit.m1$converged == FALSE)  
  {  
    print("glm for logit.m1 is not convergent.")  
    next  
  }  
  
  logit.m2 <- glm(Y ~ X2 + X3, data = data.m2, family = "binomial")  

  if(logit.m2$converged == FALSE)  
  {  
    print("glm for logit.m2 is not convergent.")  
    next  
  }  
  
  logit.m3 <- glm(Y ~ X1 + X3, data = data.m3, family = "binomial")  
  
  if(logit.m3$converged == FALSE)  
  {  
    print("glm for logit.m3 is not convergent.")  
    next  
 }  
  
  
 #--####################################################################  
  #--### Obtain the estimators of the parameters in the reduced models.##   
  #--####################################################################  

  theta.m1 = logit.m1$coefficients  
  theta.m2 = logit.m2$coefficients  
  theta.m3 = logit.m3$coefficients  
  
  
  #--####################################################################  
  #--#### The following is used to calculate robust variance estimates###  
  #--#### for each of three parameter vectors in the reduced models. ####  
  #--#### This is needed if the user inputs the variance-covariance #####  
  #--### matrices from each of studies. Later, we show that it can be ##   
  #--### computed from reference data set if these are not provided ####    
  #--###################################################################  
  
  #--####################################################################  
  #--### Find the covariance matrix estimators for the reduced models####  
  #--####################################################################  
  
  
  #--#############################  
  #-- Basic notations for inputs#  
  #--#############################  
  
#-- Number of data sets---#   

K = 3   
  
#-- index set A1, the indexes of the covariates of data set 1.--#  

  A1 = c(1, 2)  
  A2 = c(2, 3) # index set A2  
  A3 = c(1, 3) # index set A3  
  
#--- X.m1.used denotes the design matrix in model 1----#  
  
  X.m1.used = cbind(rep(1, n1), X.m1\[, A1, drop=F\])  
  X.m2.used = cbind(rep(1, n2), X.m2\[, A2, drop=F\])  
  X.m3.used = cbind(rep(1, n3), X.m3\[, A3, drop=F\])  
  
  
  

#----- Computing robust estimate of Sigma.m1 -------#  
  
  T.1 = matrix(rep(0, (length(A1)+1)^2), nrow=length(A1)+1)  
  T.2 = T.1  
  
  for (i in 1:n1)  
  {  
    a = as.vector(exp(-X.m1.used\[i, , drop=F\]%\*%theta.m1))  
    T.1 = T.1 + (a/(1+a)^2) \* (t(X.m1.used\[i, , drop=F\])%\*%X.m1.used\[i, , drop=F\])  
  }  
  
  for (i in 1:n1)  
  {  
    a = as.vector(1/( 1 + exp(-X.m1.used\[i, , drop=F\]%\*%theta.m1)))  
    T.2 = T.2 + (Y.m1\[i\]-a)^2 \* (t(X.m1.used\[i, , drop=F\])%\*%X.m1.used\[i, , drop=F\])  
  }  
  
  Sigma.m1 = solve(T.1)%\*%T.2%\*%solve(T.1)   
  
  
#----- Computing robust estimate of Sigma.m2 -------#  

  
  T.1 = matrix(rep(0, (length(A2)+1)^2), nrow=length(A2)+1)  
  T.2 = T.1  
  
  for (i in 1:n2)  
  {  
    a = as.vector(exp(-X.m2.used\[i, , drop=F\]%\*%theta.m2))  
    T.1 = T.1 + (a/(1+a)^2) \* (t(X.m2.used\[i, , drop=F\])%\*%X.m2.used\[i, , drop=F\])  
  }  
  
  for (i in 1:n2)  
  {  
    a = as.vector(1/( 1 + exp(-X.m2.used\[i, , drop=F\]%\*%theta.m2)))  
    T.2 = T.2 + (Y.m2\[i\]-a)^2 \* (t(X.m2.used\[i, , drop=F\])%*%X.m2.used\[i, , drop=F\])  
  }  
  
  Sigma.m2 = solve(T.1)%\*%T.2%\*%solve(T.1)  
  
  
 



 #----- Computing robust estimate of Sigma.m3 -------#   

  
  T.1 = matrix(rep(0, (length(A3)+1)^2), nrow=length(A3)+1)  
  T.2 = T.1  
  
  for (i in 1:n3)  
  {  
    a = as.vector(exp(-X.m3.used\[i, , drop=F\]%\*%theta.m3))  
    T.1 = T.1 + (a/(1+a)^2) \* (t(X.m3.used\[i, , drop=F\])%\*%X.m3.used\[i, , drop=F\])  
  }  
  
  for (i in 1:n3)  
  {  
    a = as.vector(1/( 1 + exp(-X.m3.used\[i, , drop=F\]%\*%theta.m3)))  
    T.2 = T.2 + (Y.m3\[i\]-a)^2 \* (t(X.m3.used\[i, , drop=F\])%\*%X.m3.used\[i, , drop=F\])  
  }  
  
  Sigma.m3 = solve(T.1)%\*%T.2%\*%solve(T.1)  
  
  
  
  
  names(theta.m1)=c("(Intercept)","Age","Height")  
  names(theta.m2)=c("(Intercept)","Height", "Weight")  
  names(theta.m3)=c("(Intercept)","Age", "Weight")  
  
  #---- now put in the GENMETA example ----#  

 #--### The results shown in Table 2 is obtained by considering #####  
 #--#### study-specific variance-covariance matrices. This is shown ##  
 #--#### below. ######################################################  
 #--##################################################################   
  study1 = list(Coeff=theta.m1,Covariance=Sigma.m1,Sample_size=n1)  
  study2 = list(Coeff=theta.m2,Covariance=Sigma.m2,Sample_size=n2)  
  study3 = list(Coeff=theta.m3,Covariance=Sigma.m3,Sample_size=n3)  

#--##### If the study-specific variance-covariance matrices are not ###  
#--##### available, then it is set to NULL. The following lines #######  
#--##### demonstrate it…                                       #######  
#--##### study1 = list(Coeff=theta.m1,Covariance=NULL,Sample_size=n1)##  
#--##### study2 = list(Coeff=theta.m2,Covariance=NULL,Sample_size=n2)##  
#--##### study3 = list(Coeff=theta.m3,Covariance=NULL,Sample_size=n3)##  

 #---- Creating the study list for GENMETA input ---#  
  
 studies = list(study1,study2,study3)  
  model = "logistic"  
  
  #------ Creating the reference data for GENMETA input -----#  

  reference = cbind(rep(1,n), X.rf)  
  colnames(reference) = c("(Intercept)","Age","Height", "Weight")  

#-------Calling the GENMETA function-------#  

  result.same = GENMETA(studies, reference, model, initial_val = c(-1.2, log(1.3), log(1.3), log(1.3)))  
  
  
#------Checking the conditions where the optim did not converge--#   

 if(sum(is.na(result.same$Est.coeff)) == 0 && is.null(result.same$Est.var.cov) != TRUE)  
  {  
    sim.matrix\[sim, \] = c(result.same$Est.coeff, diag(result.same$Est.var.cov))  
  }  
  
  print(sim)  
}  

stop.time = Sys.time()  

elapsed.time = stop.time - start.time  


#-----Computing the 95% confidence-intervals---#  
Lower.CI = na.omit(sim.matrix)\[,1:4\] - 1.96\*sqrt(na.omit(sim.matrix)\[,5:8\])  
Upper.CI = na.omit(sim.matrix)\[,1:4\] + 1.96\*sqrt(na.omit(sim.matrix)\[,5:8\])  

#------Computing the coverage rate---------#    

Coverage.rate = matrix(0,dim(na.omit(sim.matrix))\[1\], 4)  
for(k in 1:dim(na.omit(sim.matrix))\[1\])  
{  
  if(beta.star\[1,1\] >= Lower.CI\[k,1\] & beta.star\[1,1\] <= Upper.CI\[k,1\])  
    Coverage.rate\[k,1\] = 1  
  if(beta.star\[2,1\] >= Lower.CI\[k,2\] & beta.star\[2,1\] <= Upper.CI\[k,2\])  
    Coverage.rate\[k,2\] = 1  
  if(beta.star\[3,1\] >= Lower.CI\[k,3\] & beta.star\[3,1\] <= Upper.CI\[k,3\])  
    Coverage.rate\[k,3\] = 1  
  if(beta.star\[4,1\] >= Lower.CI\[k,4\] & beta.star\[4,1\] <= Upper.CI\[k,4\])  
    Coverage.rate\[k,4\] = 1  
}  

square = function(x)  
{  
  sum(x^2)  
}  

#----Computing the average length of CI----#  
  
Average.length = apply((Upper.CI-Lower.CI), 2, mean)  

#------Computing root mean square error -----#  

RMSE = sqrt(apply(t(t(na.omit(sim.matrix)\[,1:4\]) - as.numeric(beta.star)), 2, square)/dim(na.omit(sim.matrix))\[1\])  

#-----Combining the results------#  

result = data.frame(cbind(apply(na.omit(sim.matrix), 2, mean)\[1:4\], beta.star, (apply(na.omit(sim.matrix), 2, mean)\[1:4\]- beta.star), sqrt(apply(na.omit(sim.matrix), 2, mean)\[5:8\]), sqrt(apply(na.omit(sim.matrix)\[,1:4\], 2, var))), apply(Coverage.rate,2, mean), Average.length, RMSE)  

colnames(result) = c("Coeff","True", "Bias", "ESD", "SD", "Coverage.rate", "AL", "RMSE")  

#-------Writing into a file-----#

write.csv(result, file = "Simulation_9.csv", row.names = F)



# Reference
Tang, R., Kundu, P. and Chatterjee, N. (2017) Generalized Meta-Analysis for Multivariate Regression Models Across Studies with Disparate Covariate Information. https://arxiv.org/abs/1708.03818



