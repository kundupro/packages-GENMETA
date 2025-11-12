myoptim <- function(no_of_studies, study_optim, ref_dat_optim, X_rbind, X_bdiag_list, C, initial_val, threshold, model_optim, missing_covariance_study_indices, different_intercept, no_of_iter_outer)
{
  beta_old <- as.vector(initial_val)
  eps_inner <- 0
  study <- study_optim
  ref_dat <- ref_dat_optim
  threshold_optim <- threshold
  status = 1
  
  #print(beta_old)
  
  X_abdiag = Reduce(magic::adiag,X_bdiag_list)

  ## Logistic regression with same intercepts
  if(different_intercept == FALSE)
  {
    iter = 0
    continue <- TRUE
    r_first_U <- c()
    r_second_U <- c()
    U <- matrix(NA, nrow(C), nrow(ref_dat))
    Gamma_hat <- matrix(NA, nrow(C), ncol(ref_dat))
    while(continue)
    {
      r = c()
      r_first <- c()
      r_second <- c()
      Dn_1 <- matrix(NA, ncol(ref_dat), nrow(C))
      Dn_2 <- c()
      W_star_first_list <- list()
      V <- c()

      k_j  <- 1

      W_temp_1 <- as.vector((1/(1 + exp(-ref_dat %*% beta_old)))*(1/(1 + exp(ref_dat %*% beta_old))))
      #print(class(W_temp_1))
      W <- diag(W_temp_1)

      #W <- diag(rep(1/W_temp_1, no_of_studies))
      #print(is.nan(W))
      e1 <- as.vector(1/(1 + exp(-ref_dat %*% beta_old)))
      e2 <- as.vector((exp(-ref_dat %*% beta_old) - 1)/(exp(-ref_dat %*% beta_old) + 1))
      e3 <- as.vector(1/(1 + exp(ref_dat %*% beta_old)))
      l <- e1*e2*e3
      #print(class(l))
      nan_indices <- which(l %in% NaN == TRUE)
      l[nan_indices] <- 0
      L <- diag(l)



      #print(sum(is.nan(W)))

      for(k in 1 : no_of_studies)
      {
        r_first[[k]] <- e1
        col_ind <-  which(colnames(ref_dat) %in% names(study[[k]][[1]]) == TRUE)
        r_second[[k]] <-  as.vector(1/(1 + exp(-ref_dat[,col_ind] %*% study[[k]][[1]])))
        #print(length(t(ref_dat[,col_ind]) %*% (r_first[[k]] - r_second[[k]])))
        Dn_2 <- c(Dn_2, as.vector(t(ref_dat[,col_ind]) %*% (r_first[[k]] - r_second[[k]])))
        ncol_study <- length(col_ind)
        #print(ncol_study)
        #print(class((r_first[[k]] - r_second[[k]])))
        Dn_1[,k_j: (k_j + ncol_study -1)] <- t(ref_dat) %*% W %*% ref_dat[, col_ind]
        #V <- c(V, t(r_first[[k]] - r_second[[k]]) %*% ref_dat[, col_ind] %*% t(ref_dat[, col_ind]) %*% L)
        W_star_first_list[[k]] <- W %*% ref_dat[,col_ind]
        k_j <- k_j + ncol_study
        r = c(r, (r_first[[k]] - r_second[[k]]))
      }

      
      #print(dim(Dn_1))
      #print(length(Dn_2))
      Dn <- Dn_1 %*% C %*% Dn_2
      
      L = diag(rep(l, no_of_studies))
      
      V = as.vector(t(r) %*% X_abdiag %*% C %*% t(X_abdiag) %*% L)
      #print(class(V))
      W_star_second <- diag(V)

      W_star_first <- Reduce(magic::adiag,W_star_first_list) %*% C %*% t(Reduce(magic::adiag,W_star_first_list))
      #print(class(W_star_first))
      W_star <- W_star_first + W_star_second

      # print(W_star_first)
      # 
      # if(is.symmetric.matrix(W_star_first) != TRUE)
      #   print("W_star_first Not symmetric")
      # 
      # if(is.symmetric.matrix(W_star_first) != TRUE)
      #   print("W_star_second Not symmetric")
      # 
      # if(is.symmetric.matrix(W_star) != TRUE)
      #   print("W_star Not symmetric")
      # 
      # if(is.negative.definite(W_star, tol=1e-8) == TRUE)
      #   print("W_star Negative definite")
      #Define Jacobian here
      J_n <- t(X_rbind) %*% W_star %*% X_rbind
      
      #Scaling Jacobian
      max_absolute_J_n = abs(max(diag(J_n)))
      J_n_init = J_n/max_absolute_J_n
      eps = max(svd(J_n_init)$d)/pracma::cond(J_n_init)
      ill_conditioned_J_n = FALSE
      # if(is.symmetric.matrix(J_n) != TRUE)
      # {
      #   print("J_n Not symmetric")
      #   print(class(W_star))
      #   print(class(X_rbind))
      #   if(is.symmetric.matrix(W_star) != TRUE)
      #     print("W_star Not symmetric")
      # }

      #print(class(J_n))
      #print(cond(J_n))
      #print(is.nan(J_n_beta))
      #print(det(J_n))
      if(pracma::cond(J_n_init) > 1000)
        {
          ill_conditioned_J_n = TRUE
          #perturb_seq = seq(eps, mean(diag(J_n_init)), 1e-08)
          #perturb_seq = seq(0, mean(diag(J_n_init)), 0.01)
          well_condition_status = TRUE 
          cond_max_iter = 2000
          perturb_seq_index = 1
          while(well_condition_status & perturb_seq_index <= cond_max_iter)
            {
              J_n_init = J_n_init + diag(eps, nrow(J_n_init))
              #J_n_init = J_n_init + diag(perturb_seq[perturb_seq_index], nrow(J_n_init))
              #J_n = J_n + perturb_seq[perturb_seq_index]* diag(diag(J_n))
              if(pracma::cond(J_n_init) <= 1000)
                well_condition_status = FALSE
              perturb_seq_index = perturb_seq_index + 1
              eps = eps + 1e-08
            }
        }else{
              J_init = J_init
        }

        if(ill_conditioned_J_n == FALSE)
          {
            beta_new <- beta_old - (solve(J_n_init, tol=1e-60) %*% (Dn/max_absolute_J_n))
          }else if(ill_conditioned_J_n == TRUE & well_condition_status == FALSE)
          {
            beta_new <- beta_old - (solve(J_n_init, tol=1e-60) %*% (Dn/max_absolute_J_n))
          }else{
            beta_new <- beta_old - (MASS::ginv(J_n) %*% Dn)
          }
      #---- lines 129 to 134 updated by adding 0.01*mean of diagonal of J_n to J_n
      #if (abs(det(J_n))>1e-20){
        #print("original")
       # J_n_inv <- solve(J_n, tol = 1e-60)
      #}else{
       #   J_n_inv <- solve(J_n+(mean(diag(J_n))*0.01)*diag(1,dim(J_n)[1]), tol = 1e-60)
      #}
      #if(det(J_n) == 0)
      #{
        #beta_old <- rep(NA, ncol(ref_dat))
        #print(det(J_n))
        #print("The Jacobian is singular")
        #break;
      #}
      ######beta_new <- beta_old - (solve(J_n, tol=1e-60) %*% (Dn/max_absolute_J_n))
      #print(D_n_beta_t)
      #print(beta_old)
      #print(beta_new)
      eps_inner <- sqrt(sum((beta_new - beta_old)^2))
      #print(eps_inner)
      beta_old <- beta_new
      iter = iter + 1
      #print("Number of iterations \n")
      print(iter)
      if(eps_inner < threshold_optim || iter > 500)
      {
         continue <- FALSE
         if(iter >= 500 && eps_inner >= threshold_optim)
         {
           status = 0
         }
      }


    }
    no_of_iter_outer <- no_of_iter_outer + 1

    if(sum(is.na(beta_old)) == 0)
    {
      #Define the indices of the studies where the estimate of the variance-covariance matix or sample size is missing(NULL)
      study_indices <- seq(1,no_of_studies,1)
      non_missing_covariance_study_indices <- study_indices[-which(study_indices %in% missing_covariance_study_indices)]
      #print(non_missing_covariance_study_indices)
      #print(missing_covariance_study_indices)

      #Define lambda_ref here
      lambda_ref <- list()
      if(length(missing_covariance_study_indices) > 0)
      {
        for(k in missing_covariance_study_indices)
        {
          col_ind <-  which(colnames(ref_dat) %in% names(study[[k]][[1]]) == TRUE)
          #Define W_lambda_ref_logistic
          w_lambda_ref_logistic_vec <- (((1/(1 + exp(ref_dat[,col_ind] %*% study[[k]][[1]])))^2)*(1/(1 + exp(-ref_dat %*% beta_old)))) + (((1/(1 + exp(-ref_dat[,col_ind] %*% study[[k]][[1]])))^2)*(1/(1 + exp(ref_dat %*% beta_old))))
          #print(w_lambda_ref_logistic_vec )
          W_lambda_ref_logistic <- diag(as.vector(w_lambda_ref_logistic_vec))
          #print(class(W_lambda_ref_logistic))
          lambda_ref[[k]] <- (t(ref_dat[,col_ind]) %*% W_lambda_ref_logistic %*% ref_dat[,col_ind])/(study[[k]][[3]])
          #print(dim(lambda_ref[[k]]))
        }

        for(k in non_missing_covariance_study_indices)
        {
          col_ind <-  which(colnames(ref_dat) %in% names(study[[k]][[1]]) == TRUE)
          #Define W_k here
          temp_weight_logistic <- as.vector((1/(1 + exp(-ref_dat[, col_ind] %*% study[[k]][[1]])))*(1/(1 + exp(ref_dat[,col_ind] %*% study[[k]][[1]]))))
          #temp_weight_logistic <- (exp(ref_dat[,col_ind] %*% study[[k]][[1]]))/((1 + exp(ref_dat[,col_ind] %*% study[[k]][[1]]))^2)
          W_k <- t(ref_dat[,col_ind]) %*% diag(temp_weight_logistic) %*% ref_dat[,col_ind]
          lambda_ref[[k]] <- (W_k %*% study[[k]][[2]] %*% t(W_k))/(nrow(ref_dat))
          #print(dim(lambda_ref[[k]]))
        }

      }
      ## When all the estimates for var-cov are provided
      if(length(missing_covariance_study_indices) == 0)
      {
        #print("all the estimates for var-cov are provided")
        for(k in 1 : no_of_studies)
        {
          col_ind <-  which(colnames(ref_dat) %in% names(study[[k]][[1]]) == TRUE)
          #Define W_k here
          temp_weight_logistic <- as.vector((1/(1 + exp(-ref_dat[, col_ind] %*% study[[k]][[1]])))*(1/(1 + exp(ref_dat[,col_ind] %*% study[[k]][[1]]))))
          #temp_weight_logistic <- (exp(ref_dat[,col_ind] %*% study[[k]][[1]]))/((1 + exp(ref_dat[,col_ind] %*% study[[k]][[1]]))^2)
          W_k <- t(ref_dat[,col_ind]) %*% diag(temp_weight_logistic) %*% ref_dat[,col_ind]
          lambda_ref[[k]] <- (W_k %*% study[[k]][[2]] %*% t(W_k))/(nrow(ref_dat))
          #print(dim(lambda_ref[[k]]))
        }

      }

      Lambda_ref <- Reduce(magic::adiag,lambda_ref)
      #print(dim(Lambda_ref))
      k_U = 1
      r_first_1 <- as.vector(1/(1 + exp(-ref_dat %*% beta_old)))
      W_Gamma_temp_1 <- as.vector((1/(1 + exp(-ref_dat %*% beta_old)))*(1/(1 + exp(ref_dat %*% beta_old))))
      for(k in 1: no_of_studies)
      {
        col_ind <-  which(colnames(ref_dat) %in% names(study[[k]][[1]]) == TRUE)
        r_first_U[[k]] <- r_first_1
        r_second_U[[k]] <- as.vector(1/(1 + exp(-ref_dat[,col_ind] %*% study[[k]][[1]])))
        U[k_U:(k_U + length(col_ind)-1), ] <- t(ref_dat[ ,col_ind]) %*% diag(r_first_U[[k]]-r_second_U[[k]])
        Gamma_hat[k_U:(k_U + length(col_ind)-1), ] <- t(ref_dat[, col_ind]) %*% diag(W_Gamma_temp_1) %*% ref_dat
        k_U <- k_U + length(col_ind)
      }

      # Defining delta_hat here...
      Delta_hat <- (U %*% t(U))/(nrow(ref_dat))

      # Defining optimal C here...
      #Scaling inv_C
      inv_C <- Lambda_ref + Delta_hat
      max_absolute_inv_C = abs(max(diag(inv_C)))
      inv_C_init = inv_C/max_absolute_inv_C
      eps_inv_C = max(svd(inv_C_init)$d) / pracma::cond(inv_C_init)
      ill_conditioned_inv_C = FALSE
      
      if(pracma::cond(inv_C_init) > 1000){
        ill_conditioned_inv_C = TRUE
        #perturb_seq = seq(eps, mean(diag(J_n_init)), 1e-08)
        #perturb_seq = seq(0, mean(diag(J_n_init)), 0.01)
        well_condition_status = TRUE 
        cond_max_iter = 2000
        perturb_seq_index = 1
        while(well_condition_status & perturb_seq_index <= cond_max_iter)
        {
          inv_C_init = inv_C_init + diag(eps, nrow(inv_C_init))
          #J_n_init = J_n_init + diag(perturb_seq[perturb_seq_index], nrow(J_n_init))
          #J_n = J_n + perturb_seq[perturb_seq_index]* diag(diag(J_n))
          if(pracma::cond(inv_C_init) <= 1000)
            well_condition_status = FALSE
          perturb_seq_index = perturb_seq_index + 1
          eps = eps + 1e-08
        }
      }else{
        inv_C_init = inv_C_init
      }
      
      if(ill_conditioned_inv_C == FALSE)
      {
        C_beta = solve(inv_C_init, tol=1e-60)
      }else if(ill_conditioned_inv_C == TRUE & well_condition_status == FALSE)
      {
        C_beta = solve(inv_C_init, tol=1e-60)
      }else{
        C_beta = MASS::ginv(inv_C_init)
      }
      
      
      #C_beta <- solve(inv_C, tol = 1e-60)
      
      Gamma_hat <- Gamma_hat/nrow(ref_dat)
      
      #info <- (t(Gamma_hat) %*% C_beta %*% Gamma_hat)
      info <- ((t(Gamma_hat) %*% C_beta %*% Gamma_hat))/max_absolute_inv_C
      
      if(det(info) == 0)
      {
        asy_var_opt = NULL
      }else{
        if(no_of_iter_outer == 1)
        {
          if(det(t(Gamma_hat) %*% Gamma_hat) == 0)
          {
            asy_var_opt = NULL
          }else{
            asy_var_opt <- (solve(t(Gamma_hat) %*% Gamma_hat, tol = 1e-60) %*% (t(Gamma_hat) %*% (Lambda_ref + Delta_hat) %*% Gamma_hat) %*% solve(t(Gamma_hat) %*% Gamma_hat, tol = 1e-60))/(nrow(ref_dat))
          }
        }
          
        # Defining the asymptotic variance with optimal C here...
        asy_var_opt <- solve(info, tol = 1e-60)/nrow(ref_dat)
        #print(asy_var_opt)
      }
    }

    if(sum(is.na(beta_old)) > 0)
    {
      asy_var_opt = NULL
      C_beta = NULL
    }

    #print(asy_var_opt)
    # Returning objects for the inner loop(NR method)
    if(status == 0)
      print("WARNING:THE ALGORITHM DID NOT CONVERGE")
    return(list("beta_optim" = beta_old, "C_optim" = C_beta, "Asy_var_optim" = asy_var_opt,  "iter_IRWLS" = iter - 1, "Status" = status))
  }








   ## Logistic regression with different intercepts
  if(different_intercept == TRUE)
  {
    W_Gamma <- list()
    X_rbind_star <- matrix(0, nrow(ref_dat)*no_of_studies, length(beta_old))
    k_X_rbind_star <- 1
    for(k in 1:no_of_studies)
    {
      X_rbind_star[(k_X_rbind_star:(k_X_rbind_star + nrow(ref_dat) -1)),k] <- 1
      X_rbind_star[(k_X_rbind_star:(k_X_rbind_star + nrow(ref_dat) -1)), ((no_of_studies + 1):length(beta_old))] <- ref_dat[,-1]
      k_X_rbind_star <- k_X_rbind_star + nrow(ref_dat)
    }
    #print(dim(X_rbind_star))
    iter = 0
    continue <- TRUE
    r_first_U <- c()
    r_second_U <- c()
    U <- matrix(NA, nrow(C), nrow(ref_dat))
    Gamma_hat <- matrix(NA, nrow(C), length(beta_old))
    while(continue)
    {
      #print(beta_old)
      beta_k <- list()
      for(k in 1:no_of_studies)
      {
        beta_k[[k]] <- beta_old[c(k,((no_of_studies + 1):length(beta_old)))]
      }
      #print(beta_k)
      r_first <- c()
      r_second <- c()
      r <- c()
      Dn_1 <- matrix(NA, length(beta_old), nrow(C))
      Dn_2 <- c()
      W <- list()
      L <- list()
      W_star_first_list <- list()
      V <- c()

      k_j  <- 1



      for(k in 1:no_of_studies)
      {
        W_temp_1 <- as.vector((1/(1 + exp(-ref_dat %*% beta_k[[k]])))*(1/(1 + exp(ref_dat %*% beta_k[[k]]))))
        W[[k]] <- diag(W_temp_1)
      }

      for(k in 1:no_of_studies)
      {
        e1 <- as.vector(1/(1 + exp(-ref_dat %*% beta_k[[k]])))
        e2 <- as.vector((exp(-ref_dat %*% beta_k[[k]]) - 1)/(exp(-ref_dat %*% beta_k[[k]]) + 1))
        e3 <- as.vector(1/(1 + exp(ref_dat %*% beta_k[[k]])))
        l <- e1*e2*e3
        nan_indices <- which(l %in% NaN == TRUE)
        #print(length(nan_indices))
        l[nan_indices] <- 0
        L[[k]] <- diag(l)
      }




      #print(sum(is.nan(W)))
      k_j_star <- 1
      for(k in 1 : no_of_studies)
      {
        r_first[[k]] <- as.vector(1/(1 + exp(-ref_dat %*% beta_k[[k]])))
        col_ind <-  which(colnames(ref_dat) %in% names(study[[k]][[1]]) == TRUE)
        r_second[[k]] <-  as.vector(1/(1 + exp(-ref_dat[,col_ind] %*% study[[k]][[1]])))
        #print(length(t(ref_dat[,col_ind]) %*% (r_first[[k]] - r_second[[k]])))
        Dn_2 <- c(Dn_2, as.vector(t(ref_dat[,col_ind]) %*% (r_first[[k]] - r_second[[k]])))
        ncol_study <- length(col_ind)
        #print(ncol_study)
        #print(class((r_first[[k]] - r_second[[k]])))
        Dn_1[,k_j: (k_j + ncol_study -1)] <- t(X_rbind_star[k_j_star:(k_j_star + nrow(ref_dat) - 1), ]) %*% W[[k]] %*% ref_dat[, col_ind]
        #V <- c(V, t(r_first[[k]] - r_second[[k]]) %*% ref_dat[, col_ind] %*% t(ref_dat[, col_ind]) %*% L[[k]])
        W_star_first_list[[k]] <- W[[k]] %*% ref_dat[,col_ind]
        k_j <- k_j + ncol_study
        k_j_star <- k_j_star + nrow(ref_dat)
        r = c(r, (r_first[[k]] - r_second[[k]]))
      }

      #print(dim(Dn_1))
      #print(length(Dn_2))
      Dn <- Dn_1 %*% C %*% Dn_2
      #print(length(Dn))
      V = as.vector(t(r) %*% X_abdiag %*% C %*% t(X_abdiag) %*% Reduce(magic::adiag,L))
      W_star_second <- diag(V)

      W_star_first <- Reduce(magic::adiag,W_star_first_list) %*% C %*% t(Reduce(magic::adiag,W_star_first_list))
      W_star <- W_star_first + W_star_second

      #Define Jacobian here
      J_n <- t(X_rbind_star) %*% W_star %*% X_rbind_star
      
      #Scaling Jacobian
      max_absolute_J_n = abs(max(diag(J_n)))
      J_n_init = J_n/max_absolute_J_n
      eps = max(svd(J_n_init)$d) / pracma::cond(J_n_init)
      ill_conditioned_J_n = FALSE
      # if(is.symmetric.matrix(J_n) != TRUE)
      # {
      #   print("J_n Not symmetric")
      #   print(class(W_star))
      #   print(class(X_rbind))
      #   if(is.symmetric.matrix(W_star) != TRUE)
      #     print("W_star Not symmetric")
      # }

      #print(class(J_n))
      #print(cond(J_n))
      #print(is.nan(J_n_beta))
      #print(det(J_n))
      if(pracma::cond(J_n_init) > 1000)
        {
          ill_conditioned_J_n = TRUE
          #perturb_seq = seq(eps, mean(diag(J_n_init)), 1e-08)
          #perturb_seq = seq(0, mean(diag(J_n_init)), 0.01)
          well_condition_status = TRUE 
          cond_max_iter = 2000
          perturb_seq_index = 1
          while(well_condition_status & perturb_seq_index <= cond_max_iter)
            {
              J_n_init = J_n_init + diag(eps, nrow(J_n_init))
              #J_n_init = J_n_init + diag(perturb_seq[perturb_seq_index], nrow(J_n_init))
              #J_n = J_n + perturb_seq[perturb_seq_index]* diag(diag(J_n))
              if(pracma::cond(J_n_init) <= 1000)
                well_condition_status = FALSE
              perturb_seq_index = perturb_seq_index + 1
              eps = eps + 1e-08
            }
        }else{
              J_init = J_init
        }

        if(ill_conditioned_J_n == FALSE)
          {
            beta_new <- beta_old - (solve(J_n_init, tol=1e-60) %*% (Dn/max_absolute_J_n))
          }else if(ill_conditioned_J_n == TRUE & well_condition_status == FALSE)
          {
            beta_new <- beta_old - (solve(J_n_init, tol=1e-60) %*% (Dn/max_absolute_J_n))
          }else{
            beta_new <- beta_old - (MASS::ginv(J_n) %*% Dn)
          }
      #--lines 401 to 406 upadted before based on adding 0.01 to mean diagonal
      #if (abs(det(J_n))>1e-20){
        #print("original")
       # J_n_inv <- solve(J_n, tol = 1e-60)
      #}else{
      #    J_n_inv <- solve(J_n+(mean(diag(J_n))*0.01)*diag(1,dim(J_n)[1]), tol = 1e-60)
      #}
      #if(det(J_n) == 0)
      #{ beta_old <- rep(NA, (ncol(ref_dat) - 1 + no_of_studies))
      #break;
      #}
      #######beta_new <- beta_old - (solve(J_n, tol=1e-60) %*% (Dn*max_absolute_J_n))
      #print(D_n_beta_t)
      #print(beta_old)
      #print(beta_new)
      eps_inner <- sqrt(sum((beta_new - beta_old)^2))
      #print(eps_inner)
      beta_old <- beta_new
      iter = iter + 1
      #print("Number of iterations \n")
      #print(iter)
      if(eps_inner < threshold_optim || iter > 500)
       {
        
        continue <- FALSE
        
        if(iter >= 500 && eps_inner >= threshold_optim)
        {
          status = 0
        }
      } 


    }

    no_of_iter_outer <- no_of_iter_outer + 1

    if(sum(is.na(beta_old)) == 0)
    {
      for(k in 1:no_of_studies)
      {
        beta_k[[k]] <- beta_old[c(k,((no_of_studies + 1):length(beta_old)))]
      }

      #Define the indices of the studies where the estimate of the variance-covariance matix or sample size is missing(NULL)
      study_indices <- seq(1,no_of_studies,1)
      non_missing_covariance_study_indices <- study_indices[-which(study_indices %in% missing_covariance_study_indices)]
      #print(non_missing_covariance_study_indices)
      #print(missing_covariance_study_indices)

      #Define lambda_ref here
      lambda_ref <- list()
      if(length(missing_covariance_study_indices) > 0)
      {
        for(k in missing_covariance_study_indices)
        {
          col_ind <-  which(colnames(ref_dat) %in% names(study[[k]][[1]]) == TRUE)
          #Define W_lambda_ref_logistic
          w_lambda_ref_logistic_vec <- (((1/(1 + exp(ref_dat[,col_ind] %*% study[[k]][[1]])))^2)*(1/(1 + exp(-ref_dat %*% beta_k[[k]])))) + (((1/(1 + exp(-ref_dat[,col_ind] %*% study[[k]][[1]])))^2)*(1/(1 + exp(ref_dat %*% beta_k[[k]]))))
          #print(w_lambda_ref_logistic_vec )
          W_lambda_ref_logistic <- diag(as.vector(w_lambda_ref_logistic_vec))
          #print(class(W_lambda_ref_logistic))
          lambda_ref[[k]] <- (t(ref_dat[,col_ind]) %*% W_lambda_ref_logistic %*% ref_dat[,col_ind])/(study[[k]][[3]])
          #print(dim(lambda_ref[[k]]))
        }

        for(k in non_missing_covariance_study_indices)
        {
          col_ind <-  which(colnames(ref_dat) %in% names(study[[k]][[1]]) == TRUE)
          #Define W_k here
          temp_weight_logistic <- as.vector((1/(1 + exp(-ref_dat[, col_ind] %*% study[[k]][[1]])))*(1/(1 + exp(ref_dat[,col_ind] %*% study[[k]][[1]]))))
          #temp_weight_logistic <- (exp(ref_dat[,col_ind] %*% study[[k]][[1]]))/((1 + exp(ref_dat[,col_ind] %*% study[[k]][[1]]))^2)
          W_k <- t(ref_dat[,col_ind]) %*% diag(temp_weight_logistic) %*% ref_dat[,col_ind]
          lambda_ref[[k]] <- (W_k %*% study[[k]][[2]] %*% t(W_k))/(nrow(ref_dat))
          #print(dim(lambda_ref[[k]]))
        }

      }
      ## When all the estimates for var-cov are provided
      if(length(missing_covariance_study_indices) == 0)
      {
        for(k in 1 : no_of_studies)
        {
          col_ind <-  which(colnames(ref_dat) %in% names(study[[k]][[1]]) == TRUE)
          #Define W_k here
          temp_weight_logistic <- as.vector((1/(1 + exp(-ref_dat[, col_ind] %*% study[[k]][[1]])))*(1/(1 + exp(ref_dat[,col_ind] %*% study[[k]][[1]]))))
          #temp_weight_logistic <- (exp(ref_dat[,col_ind] %*% study[[k]][[1]]))/((1 + exp(ref_dat[,col_ind] %*% study[[k]][[1]]))^2)
          W_k <- t(ref_dat[,col_ind]) %*% diag(temp_weight_logistic) %*% ref_dat[,col_ind]
          lambda_ref[[k]] <- (W_k %*% study[[k]][[2]] %*% t(W_k))/(nrow(ref_dat))
          #print(dim(lambda_ref[[k]]))
        }

      }

      Lambda_ref <- Reduce(magic::adiag,lambda_ref)
      #print(dim(Lambda_ref))

      ## Defining Gamma_hat...
      k_U = 1
      k_j_star <- 1
      #beta_k <- list()

      #r_first_1 <- as.vector(1/(1 + exp(-ref_dat %*% beta_old)))
      #W_Gamma_temp_1 <- as.vector((1/(1 + exp(-ref_dat %*% beta_old)))*(1/(1 + exp(ref_dat %*% beta_old))))
      for(k in 1: no_of_studies)
      {
        col_ind <-  which(colnames(ref_dat) %in% names(study[[k]][[1]]) == TRUE)
        r_first_U[[k]] <- as.vector(1/(1 + exp(-ref_dat %*% beta_k[[k]])))
        r_second_U[[k]] <- as.vector(1/(1 + exp(-ref_dat[,col_ind] %*% study[[k]][[1]])))
        U[k_U:(k_U + length(col_ind)-1), ] <- t(ref_dat[ ,col_ind]) %*% diag(r_first_U[[k]]-r_second_U[[k]])
        W_Gamma_temp_1 <- as.vector((1/(1 + exp(-ref_dat %*% beta_k[[k]])))*(1/(1 + exp(ref_dat %*% beta_k[[k]]))))
        Gamma_hat[k_U:(k_U + length(col_ind)-1), ] <- t(ref_dat[, col_ind]) %*% diag(W_Gamma_temp_1) %*% X_rbind_star[k_j_star:(k_j_star + nrow(ref_dat) - 1), ]
        k_U <- k_U + length(col_ind)
        k_j_star <- k_j_star + nrow(ref_dat)
      }

      Delta_hat <- (U %*% t(U))/(nrow(ref_dat))
      inv_C <- Lambda_ref + Delta_hat
      
      # if(pracma::cond(inv_C) > 1000)
      # {
      #   perturb_seq_C = seq(0, mean(diag(inv_C)), 0.01)
      #   well_condition_status_C = TRUE 
      #   perturb_seq_index_C = 1
      #   while(well_condition_status_C)
      #   {
      #     inv_C = inv_C + perturb_seq_C[perturb_seq_index_C]* diag(diag(inv_C))
      #     if(pracma::cond(inv_C) <= 1000)
      #       well_condition_status_C = FALSE
      #     perturb_seq_index_C = perturb_seq_index_C + 1
      #   }
      # }
      
      #Scaling inv_C
      inv_C <- Lambda_ref + Delta_hat
      max_absolute_inv_C = abs(max(diag(inv_C)))
      inv_C_init = inv_C/max_absolute_inv_C
      eps_inv_C = max(svd(inv_C_init)$d) / pracma::cond(inv_C_init)
      ill_conditioned_inv_C = FALSE
      
      
      
      if(pracma::cond(inv_C_init) > 1000){
        ill_conditioned_inv_C = TRUE
        #perturb_seq = seq(eps, mean(diag(J_n_init)), 1e-08)
        #perturb_seq = seq(0, mean(diag(J_n_init)), 0.01)
        well_condition_status = TRUE 
        cond_max_iter = 2000
        perturb_seq_index = 1
        while(well_condition_status & perturb_seq_index <= cond_max_iter)
        {
          inv_C_init = inv_C_init + diag(eps, nrow(inv_C_init))
          #J_n_init = J_n_init + diag(perturb_seq[perturb_seq_index], nrow(J_n_init))
          #J_n = J_n + perturb_seq[perturb_seq_index]* diag(diag(J_n))
          if(pracma::cond(inv_C_init) <= 1000)
            well_condition_status = FALSE
          perturb_seq_index = perturb_seq_index + 1
          eps = eps + 1e-08
        }
      }else{
        inv_C_init = inv_C_init
      }
      
      if(ill_conditioned_inv_C == FALSE)
      {
        C_beta = solve(inv_C_init, tol=1e-60)
      }else if(ill_conditioned_inv_C == TRUE & well_condition_status == FALSE)
      {
        C_beta = solve(inv_C_init, tol=1e-60)
      }else{
        C_beta = MASS::ginv(inv_C_init)
      }
      
      
      #C_beta <- solve(inv_C, tol = 1e-60)
      
      Gamma_hat <- Gamma_hat/nrow(ref_dat)

      #info <- (t(Gamma_hat) %*% C_beta %*% Gamma_hat)
      info <- ((t(Gamma_hat) %*% C_beta %*% Gamma_hat))/max_absolute_inv_C
      
      if(det(info) == 0)
      {
        asy_var_opt = NULL
      }else{
        if(no_of_iter_outer == 1)
          asy_var_opt <- (solve(t(Gamma_hat) %*% Gamma_hat, tol = 1e-60) %*% (t(Gamma_hat) %*% (Lambda_ref + Delta_hat) %*% Gamma_hat) %*% solve(t(Gamma_hat) %*% Gamma_hat, tol = 1e-60))/(nrow(ref_dat))
        asy_var_opt <- solve(info, tol = 1e-60)/nrow(ref_dat)
      }

    }

    #print(beta_old)
    if(sum(is.na(beta_old)) > 0)
    {
      asy_var_opt = NULL
      C_beta = NULL
    }
    #print(class(info))
    #asy_var_beta <- (solve(info, tol = 1e-30))/nrow(ref_dat)
    #asy_var_beta <- diag(3)
    #print(asy_var_opt)
    if(status == 0)
      print("WARNING:THE ALGORITHM DID NOT CONVERGE")
    return(list("beta_optim" = beta_old, "C_optim" = C_beta, "Asy_var_optim" = asy_var_opt,  "iter_IRWLS" = iter - 1, "Status" = status))
  }

}
