################################################################################
## SHARED FUNTIONS FOR SIMULATION (updated 8 FEB 2022)
################################################################################

library(Matrix)
library(mvtnorm)
library(reshape2)
library(dplyr)
library(ggplot2)
library(foreach)
library(cowplot)
library(msos)
library(lavaan)
library(semTools)
library(nlme)
library(parallel)


################################################################################
## FOR OU COVARIANCE MATRIX
################################################################################

# Calculate Kronecker sum
kronsum <- function(A, B){
  dimA <- nrow(A)
  dimB <- nrow(B)
  kronecker(A, diag(dimB)) + kronecker(diag(dimA), B)
}

# for marginal covariance (if eta(t_0) is unknown)
Gamma2 <- function(theta_boup, sigma_boup, t, s){ # t <= s
  if(t > s){
    return('error: t <= s is required')
  }
  theta_ks <- kronsum(theta_boup, theta_boup)
  theta_ks.inv <- solve(theta_ks)
  
  term1 <- theta_ks.inv
  term2 <- expm(t*theta_ks - kronsum(theta_boup*t, theta_boup*s))
  term3 <- matrix(sigma_boup %*% t(sigma_boup), ncol = 1)
  
  out <- matrix(term1 %*% term2 %*% term3,
                nrow = nrow(theta_boup), ncol = nrow(theta_boup))
  
  return(t(as.matrix(out,
                     row = nrow(theta_boup), ncol = nrow(theta_boup))))
}


# Calculate sigma using identifiability constraint that the OUP has a stationary
#  variance of 1
calc_sigma <- function(theta){
  p <- nrow(theta)
  theta_ks <- kronecker(theta, diag(p)) + kronecker(diag(p), theta)
  theta_ks.inv <- solve(theta_ks)
  sigma_vec <- c(diag(p))
  X <- theta_ks.inv[sigma_vec == 1, sigma_vec == 1]
  Y <- matrix(1, nrow = p, ncol = 1)
  sigma2_solutions <- c(solve(X, Y))
  return(sqrt(sigma2_solutions))
}


################################################################################
## FOR RE-SCALING OUP PARAMETERS TO SATISFY THE IDENTIFIABILITY CONSTRAINT
################################################################################

# Function for the covariance of the OU process (just the diagonal block)
construct_V <- function(theta, sigma) {
  kron_sum_theta = kronsum(theta, theta)  
  vec_sigmasq = as.vector(sigma %*% t(sigma))
  matrix(solve(kron_sum_theta, vec_sigmasq), nrow = nrow(theta))
}

# Project theta to identifiability space of theta*
update_theta <- function(theta, constants) {
  constants %*% theta %*% solve(constants)
}

# Project sigma to identifiability space of sigma*
update_sigma <- function(sigma, constants) {
  constants %*% sigma
}

# Calculate the constants that rescale the original OUP to one that satisfies
#   the identifiability constraint of constant stationary variance = 1
# Goal: want (c1, c2) s.t. cov2cor(Psi(theta, sigma)) = Psi(theta*, sigma*)
# Input: original OU parameters, guess of rescaling constants
# Output: rescaling constants (assuming > 0) and optimization convergence code
calc_constants <- function(guess_constants_vec, theta, sigma){
  log_constants_vec <- log(guess_constants_vec)
  obj_fun <- function(log_constants_vec, theta, sigma){
    constants <- diag(x = exp(log_constants_vec),
                      nrow = length(log_constants_vec))
    target_correlation <- cov2cor(construct_V(theta, sigma))
    theta_star <- update_theta(theta, constants)
    sigma_star <- update_sigma(sigma, constants)
    # cat('  ..... current sigma* is', diag(sigma_star), '\n')
    diag(sigma_star) <- pmax(diag(sigma_star), rep(0.001, nrow(sigma_star)))
    rescaled_covariance <- construct_V(theta_star, sigma_star)
    cur_diff <- sum(abs(target_correlation - rescaled_covariance))
    return(cur_diff)
  }
  if(nrow(theta)==1){
    opt <- optimize(interval = c(log(1e-4), log(99999)),
                    f = obj_fun, 
                    theta = theta, sigma = sigma)
    return(list(c_vec = exp(opt$minimum), code = 0))
  }else{
    opt <- optim(par = rep(0, length(log_constants_vec)), fn = obj_fun, theta = theta,
                 sigma = sigma)
    return(list(c_vec = exp(opt$par), code = opt$convergence))
  }

}
# example: calc_constants(c(1, 1), theta, sigma)


################################################################################
# ANALYTICAL GRADIENT (for measurement submodel block update)
################################################################################

# gradient w.r.t. Lambda
grad_lambda <- function(lambda_vec, Yi, Sigma_u, Theta, Psi){ #lambda_vec
  ni <- length(Yi)/k
  ones_ni <- matrix(1, nrow = ni, ncol = ni)
  lambda_vec[1] <- exp(lambda_vec[1])
  Lambda[nonzero==1] <- lambda_vec
  A <- kronecker(diag(ni), Lambda) %*% Psi %*% kronecker(diag(ni), t(Lambda))
  BC <- kronecker(ones_ni, Sigma_u) + kronecker(diag(ni), Theta) 
  Sigma_star <- A + BC
  Sigma_star <- (Sigma_star + t(Sigma_star)) / 2 # make exactly symmetric
  Sigma_star.chol <- chol(Sigma_star)
  Sigma_star.inv <- chol2inv(Sigma_star.chol)
  j <- 1
  grad_vec <- rep(0, length(lambda_vec))
  
  for (i in 1:length(Lambda)){
    if (nonzero[i] == 1){
      J.k <- matrix(0, nrow = k, ncol = p)
      J.k[i] <- 1
      Ci <- kronecker(diag(ni), Lambda) %*% Psi %*% kronecker(diag(ni), t(J.k))
      Ci <- Ci + t(Ci)
      inside <- -Sigma_star.inv %*% Ci %*% Sigma_star.inv
      term2 <- as.numeric(t(Yi) %*% inside %*% Yi)
      term1 <- sum(diag(Sigma_star.inv %*% Ci))
      
      grad_vec[j] <- term1 + term2
      
      j <- j + 1
    }
  }
  ddlog <- rep(1, length(lambda_vec))
  ddlog[1] <- lambda_vec[1]
  return(-0.5*grad_vec*ddlog)
}

# gradient w.r.t. Sigma_u
grad_sigma_u <- function(R.sigma_u_vec, Yi, Lambda, Theta, Psi){
  ni <- length(Yi)/k
  ones_ni <- matrix(1, nrow = ni, ncol = ni)
  
  R.Sigma_u <- matrix(0, nrow = k, ncol = k)
  R.Sigma_u[upper.tri(R.Sigma_u, diag = TRUE)] <- R.sigma_u_vec
  diag(R.Sigma_u) <- exp(diag(R.Sigma_u))
  Sigma_u <- crossprod(R.Sigma_u) # parameterized using cholesky decomp 
  
  A <- kronecker(diag(ni), Lambda) %*% Psi %*% kronecker(diag(ni), t(Lambda))
  BC <- kronecker(ones_ni, Sigma_u) + kronecker(diag(ni), Theta) 
  Sigma_star <- A + BC
  Sigma_star <- (Sigma_star + t(Sigma_star)) / 2 # make exactly symmetric
  Sigma_star.chol <- chol(Sigma_star)
  Sigma_star.inv <- chol2inv(Sigma_star.chol)

  # grad_vec <- rep(0, length(R.sigma_u_vec)) # for non-diag Sigma_u
  grad_vec <- rep(0, k) # for diag Sigma_u
  
  # # for non-diag Sigma_u
  # j <- 1
  # for (i in 1:length(R.Sigma_u)){
  #   if (upper.tri(R.Sigma_u, diag = TRUE)[i] == 1){
  #     J.k <- matrix(0, nrow = k, ncol = k)
  #     J.k[i] <- 1
  #     ddR <- t(J.k) %*% R.Sigma_u + t(R.Sigma_u) %*% J.k
  #     Ci <- kronecker(ones_ni, ddR)
  # 
  #     term1 <- sum(diag(Sigma_star.inv %*% Ci))
  # 
  #     inside <- -Sigma_star.inv %*% Ci %*% Sigma_star.inv
  #     term2 <- as.numeric(t(Yi) %*% inside %*% Yi)
  # 
  #     grad_vec[j] <- term1 + term2
  #     j <- j + 1
  #   }
  # }
  
  # for diagonal Sigma_u
  for (i in 1:k){
      J.k <- matrix(0, nrow = k, ncol = k)
      J.k[i,i] <- 1
      ddR <- t(J.k) %*% R.Sigma_u + t(R.Sigma_u) %*% J.k
      Ci <- kronecker(ones_ni, ddR)

      term1 <- sum(diag(Sigma_star.inv %*% Ci))

      inside <- -Sigma_star.inv %*% Ci %*% Sigma_star.inv
      term2 <- as.numeric(t(Yi) %*% inside %*% Yi)
      
      grad_vec[i] <- term1 + term2
  }
  
  D <- matrix(0, nrow = k, ncol = k) # Note: switch 0 to 1 to estimate correlations
  diag(D) <- diag(R.Sigma_u)
  # D.vec <- D[upper.tri(D, diag = TRUE)] # for non-diag Sigma_u
  D.vec <- diag(D) # for diag Sigma_u
  return(-0.5 * grad_vec*D.vec)
}

# for Sigma_epsilon (called theta here)
grad_theta <- function(log.theta_vec, Yi, Lambda, Sigma_u, Psi){
  ni <- length(Yi)/k
  ones_ni <- matrix(1, nrow = ni, ncol = ni)
  
  Theta <- diag(exp(log.theta_vec)^2) # parameterized as log(sigma_e)
  
  A <- kronecker(diag(ni), Lambda) %*% Psi %*% kronecker(diag(ni), t(Lambda))
  BC <- kronecker(ones_ni, Sigma_u) + kronecker(diag(ni), Theta) 
  Sigma_star <- A + BC
  Sigma_star <- (Sigma_star + t(Sigma_star)) / 2 # make exactly symmetric
  Sigma_star.chol <- chol(Sigma_star)
  Sigma_star.inv <- chol2inv(Sigma_star.chol)

  grad_vec <- rep(0, length(log.theta_vec))
  
  for (i in 1:length(log.theta_vec)){
    J.k <- matrix(0, nrow = k, ncol = k)
    J.k[i,i] <- 2*exp(log.theta_vec[i])
    # J.k[i,i] <- 2*(log.theta_vec[i]) 
    
    Ci <- kronecker(diag(ni), J.k)
    term1 <- sum(diag(Sigma_star.inv %*% Ci))
    inside <- -Sigma_star.inv %*% Ci %*% Sigma_star.inv
    term2 <- as.numeric(t(Yi) %*% inside %*% Yi)
    grad_vec[i] <- term1 + term2
  }
  
  return(-0.5*grad_vec*exp(log.theta_vec))
  # return(-0.5*grad_vec) 
  
}

# gradient of the log-likelihood for i w.r.t. all meas submod parameters
lli_grad <- function(m_params, Yi, Psi_i){
  lambda_vec <- m_params[1:k]
  R.sigma_u_vec <- rep(0, choose(k, 2)+k)

  m = diag(k)
  upperm = upper.tri(m, diag = T)
  m[upperm] = 1:(choose(k, 2) + k)
  
  R.sigma_u_vec[diag(m)] <- m_params[(k+1):(2*k)]
  log.theta_vec <- m_params[(2*k+1):(3*k)]

  # set up Lambda
  cur_Lambda <- matrix(0, nrow = k, ncol = p)
  cur_Lambda[nonzero==1] <- lambda_vec
  cur_Lambda[1,1] <- exp(cur_Lambda[1,1])
  # set up Sigma_u
  R.Sigma_u <- matrix(0, nrow = k, ncol = k)
  R.Sigma_u[upper.tri(R.Sigma_u, diag = TRUE)] <- R.sigma_u_vec
  diag(R.Sigma_u) <- exp(diag(R.Sigma_u))
  cur_Sigma_u <- crossprod(R.Sigma_u) # parameterized using cholesky decomp
  # set up Sigma_e (called Theta here)
  cur_Theta <- diag((exp(log.theta_vec))^2) # parameterized as log(sigma_e)
  # cur_Theta <- diag(((log.theta_vec))^2) 
  
  ni = length(Yi)/k
  Sigma_star_inv <- calc_Sigma_star_inv(ni, cur_Lambda, cur_Sigma_u,
                                        cur_Theta, Psi_i)
  
  c(c(grad_lambda_cpp_slow(k, ni, Yi, cur_Lambda, Psi_i, Sigma_star_inv,
                      nonzero)),
    c(grad_sigma_u_cpp_slow(k, ni, Yi, cur_Sigma_u, Sigma_star_inv)),
    c(grad_sigma_e_cpp_slow(k, ni, Yi, cur_Theta, Sigma_star_inv)))
}

# gradient for negative log likelihood w.r.t FA parameters
ll_grad <- function(m_params, all_person_Y_centered, Psi_list){
  # define non-zero indices of covariance matrices (for Sigma_u and Sigma_e)
  m <- diag(k)
  upperm <- upper.tri(m, diag = T)
  m[upperm] = 1:(choose(k, 2) + k)
  sigma_u_inds <- diag(m)+k
  sigma_e_inds <- (max(sigma_u_inds)+1):(max(sigma_u_inds)+k)
  new_grad_inds <- c(1:k, sigma_u_inds, sigma_e_inds)
  
  ########
  lambda_vec <- m_params[1:k]
  R.sigma_u_vec <- rep(0, choose(k, 2)+k)
  m = diag(k)
  upperm = upper.tri(m, diag = T)
  m[upperm] = 1:(choose(k, 2) + k)
  R.sigma_u_vec[diag(m)] <- m_params[(k+1):(2*k)]
  log.theta_vec <- m_params[(2*k+1):(3*k)]
  
  # set up Lambda
  cur_Lambda <- matrix(0, nrow = k, ncol = p)
  cur_Lambda[nonzero==1] <- lambda_vec
  cur_Lambda[1,1] <- exp(cur_Lambda[1,1])
  # set up Sigma_u
  R.Sigma_u <- matrix(0, nrow = k, ncol = k)
  R.Sigma_u[upper.tri(R.Sigma_u, diag = TRUE)] <- R.sigma_u_vec
  diag(R.Sigma_u) <- exp(diag(R.Sigma_u))
  cur_Sigma_u <- crossprod(R.Sigma_u) # parameterized using cholesky decomp
  # set up Sigma_e
  cur_Sigma_e <- diag((exp(log.theta_vec))^2) # parameterized as log(sigma_e)
  ########
  
  # cat('using', n_cores, 'cores \n')
  # cl <- parallel::makeCluster(n_cores)
  # doParallel::registerDoParallel(cl)
  clusterExport(cl, c("measurement_times", "k", "lli_grad", "p", "nonzero"))
  grad_list <- foreach(i = 1:n_subj, .packages = 'mrabbott') %dopar% {
    ni <- length(measurement_times[i,!is.na(measurement_times[i,])])
    Y_i <- all_person_Y_centered[i, 1:(k*ni)]
    Psi_i <- Psi_list[[i]]
    Sigma_star_inv <- calc_Sigma_star_inv(ni, cur_Lambda, cur_Sigma_u,
                                          cur_Sigma_e, Psi_i)
    I_SYYt = diag(k*ni) - Sigma_star_inv %*% Y_i %*% t(Y_i)
    Sigma_term <-  I_SYYt %*% Sigma_star_inv
    fa_grads(k, ni, Y_i, cur_Lambda, cur_Sigma_u, cur_Sigma_e, Psi_i,
             Sigma_star_inv, nonzero, I_SYYt, Sigma_term)
  }
  # stopImplicitCluster()
  # parallel::stopCluster(cl)
  
  grad_mat <- matrix(unlist(grad_list), ncol = length(new_grad_inds),
                     byrow = TRUE)
  out <- -colSums(grad_mat)
  return(out)
}

################################################################################
#### Define FABOUP likelihood w.r.t FA parameters
################################################################################

# for single subject
FA_negllk_i <- function(m_params_vec, i){
  ni <- length(measurement_times[i,!is.na(measurement_times[i,])])
  # update Lambda
  lambda_vec <- m_params_vec[1:k]
  lambda_vec[1] <- exp(m_params_vec[1])
  cur_Lambda <- matrix(0, nrow = k, ncol = p)
  cur_Lambda[nonzero == 1] <- lambda_vec
  cur_Lambda_mat <- kronecker(diag(ni), cur_Lambda)
  # update Sigma_u
  log_sigma_u_vec <- m_params_vec[(k+1):(2*k)]
  # add lower bound to prevent est from getting stuck at -Inf
  log_sigma_u_vec <- pmax(log_sigma_u_vec, log(1e-4))
  cur_Sigma_u <- diag((exp(log_sigma_u_vec))^2)
  cur_Sigma_u_mat <- kronecker(matrix(1, nrow = ni, ncol = ni), cur_Sigma_u)
  # update Theta
  log_theta_vec <- m_params_vec[(2*k+1):(3*k)]
  # add lower bound to prevent est from getting stuck at -Inf
  log_theta_vec <- pmax(log_theta_vec, log(1e-4))
  cur_Theta <- diag((exp(log_theta_vec))^2)
  cur_Theta_mat <- kronecker(diag(ni), cur_Theta)
  # update covariance matrix for Yi
  cur_Sigma_star <- cur_Lambda_mat %*% Psi_list[[i]]
  cur_Sigma_star <- cur_Sigma_star %*% t(cur_Lambda_mat)
  cur_Sigma_star <- cur_Sigma_star + cur_Sigma_u_mat + cur_Theta_mat
  cur_Sigma_star <- (cur_Sigma_star + t(cur_Sigma_star))/2
  
  # return negative log likelihood
  return(0.5 * (logdet(cur_Sigma_star)
                + tr(solve(cur_Sigma_star, empirical_CovY_list[[i]]))))
}

# negative log likelihood function (w.r.t FA parameters) for all subjects
FA_negllk <- function(m_params_vec, Psi_list){
  # likelihood
  out <- (sum(unlist(lapply(1:n_subj,
                            FUN = function(i){ FA_negllk_i(m_params_vec, i) }))))
  # gradient
  attr(out, 'gradient') <- ll_grad(m_params_vec, all_person_Y_centered,
                                   Psi_list)
  
  return(out)
}

################################################################################
#### Define FABOUP likelihood w.r.t OUP parameters
################################################################################

### 1. Define likelihood for Y w.r.t. to BOUP parameters
# NOTE: sigma is estimated on log scale
BOUP_neg_loglik_i <- function(params, c_vec, i) {
  # cat('person #:', i, '\n')
  theta_param_length = length(as.vector(theta_boup))
  theta_current = matrix(params[1:theta_param_length], nrow = nrow(theta_boup))
  sigma_current = diag(x = exp(params[(theta_param_length+1):length(params)]),
                       nrow = length(exp(params[(theta_param_length+1):length(params)])))
  
  # Rescale BOUP parameters to satisfy the identifiability constraint
  theta_current <- update_theta(theta_current,
                                diag(x = c_vec, nrow = length(c_vec)))
  sigma_current <- update_sigma(sigma_current,
                                diag(x = c_vec, nrow = length(c_vec)))
  
  meas_times_i <- measurement_times[i,!is.na(measurement_times[i,])]
  ni <- length(meas_times_i)
  
  theta_t = t(theta_current)
  sigma2_vec = matrix(sigma_current %*% t(sigma_current), ncol = 1)
  kron_sum_theta = kronsum(theta_current, theta_current)
  cur_Omega_i <- calc_precision_cpp(kron_sum_theta = kron_sum_theta,
                                    theta = theta_current, theta_t = theta_t,
                                    sigma2_vec = sigma2_vec, ni = ni,
                                    times = meas_times_i)
  cur_Psi_i <- solve(cur_Omega_i)
  cur_Psi_i <- (cur_Psi_i + t(cur_Psi_i))/2
  
  # update Cov(Y_i)
  Lambda_mat <- kronecker(diag(ni), Lambda) %*% cur_Psi_i %*% kronecker(diag(ni), t(Lambda))
  Sigma_u_mat <- kronecker(matrix(1, nrow = ni, ncol = ni), Sigma_u)
  Sigma_e_mat <- kronecker(diag(ni), Sigma_e)
  CovY_i <- Lambda_mat + Sigma_u_mat + Sigma_e_mat # Cov(Y_i) = Sigma^*_i
  CovY_i <- (CovY_i + t(CovY_i))/2
  
  # -loglik(Y)
  0.5 * (logdet(CovY_i) + tr(solve(CovY_i, empirical_CovY_list[[i]])))
}

BOUP_negloglik_n <- function(params, c_vec) {
  theta_param_length <- length(as.vector(theta_boup))
  theta_current <- matrix(params[1:theta_param_length], nrow = nrow(theta_boup))
  
  # check that theta corresponds to mean-reverting process
  theta_ok <- as.numeric(all(Re(eigen(theta_current)$values) > 0))

  if (theta_ok == 0){
    cat('warning: eigenvalues of theta do not have positive real parts \n')
    return(999999)
  }else{
    sigma_current = diag(x = exp(params[(theta_param_length+1):length(params)]),
                         nrow = length(exp(params[(theta_param_length+1):length(params)])))
    current_c_fit <- calc_constants(c_vec, theta_current, sigma_current)
    current_c_vec <- current_c_fit$c_vec
    
    # return negative log likelihood
    return(sum(unlist(lapply(1:n_subj, FUN = function(i){
      BOUP_neg_loglik_i(params, current_c_vec, i) }))))
  }
}


################################################################################
#### Define FABOUP likelihood w.r.t ALL parameters
################################################################################

# for single subject
FABOUP_negllk_i <- function(m_params_vec, i){
  
  meas_times_i <- measurement_times[i,!is.na(measurement_times[i,])]
  ni <- length(meas_times_i)
  # update Lambda
  lambda_vec <- m_params_vec[1:k]
  lambda_vec[1] <- exp(m_params_vec[1])
  cur_Lambda <- matrix(0, nrow = k, ncol = p)
  cur_Lambda[nonzero == 1] <- lambda_vec
  cur_Lambda_mat <- kronecker(diag(ni), cur_Lambda)
  # update Sigma_u
  log_sigma_u_vec <- m_params_vec[(k+1):(2*k)]
  cur_Sigma_u <- diag((exp(log_sigma_u_vec))^2)
  cur_Sigma_u_mat <- kronecker(matrix(1, nrow = ni, ncol = ni), cur_Sigma_u)
  # update Sigma_e (called Theta here)
  log_theta_vec <- m_params_vec[(2*k+1):(3*k)]
  cur_Theta <- diag((exp(log_theta_vec))^2)
  cur_Theta_mat <- kronecker(diag(ni), cur_Theta)
  # update theta_OU
  theta_param_length = length(as.vector(theta_boup))
  theta_current = matrix(m_params_vec[(3*k+1):(3*k+theta_param_length)], nrow = nrow(theta_boup))
  # update sigma_OU
  sigma_current = diag(x = exp(m_params_vec[(3*k+theta_param_length+1):length(m_params_vec)]),
                       nrow = p)
  # update OUP covariance matrix
  theta_t = t(theta_current)
  sigma2_vec = matrix(sigma_current %*% t(sigma_current), ncol = 1)
  kron_sum_theta = kronsum(theta_current, theta_current)
  ni = length(c(meas_times_i))
  cur_Omega_i <- calc_precision_cpp(kron_sum_theta = kron_sum_theta,
                                    theta = theta_current, theta_t = theta_t,
                                    sigma2_vec = sigma2_vec,
                                    ni = ni, times = c(meas_times_i))
  # for debugging
  # if (det(cur_Omega_i) == 0){
  #   print(theta_current)
  #   print(sigma_current)
  #   print(cur_Omega_i)
  #   write.csv(data.frame(omega_by.cols = as.vector(cur_Omega_i)),
  #             paste0(wd, 'data/omega_i.csv'))
  # }
  # cat('    eigenvalues of precision matrix are', eigen(cur_Omega_i)$values, '      \n \n ')
  # if (!all(eigen(cur_Omega_i)$values > 0)){
  #   cat('~ ~ ~ ~ using generalize inverse... ~ ~ ~ ~ \n \n \n ')
  #   cur_Psi_i <- MASS::ginv(cur_Omega_i)
  # }else{
  #   
  # }
  
  cur_Psi_i <- solve(cur_Omega_i)
  cur_Psi_i <- (cur_Psi_i + t(cur_Psi_i))/2

  # update covariance matrix for Yi
  cur_Sigma_star <- cur_Lambda_mat %*% cur_Psi_i
  cur_Sigma_star <- cur_Sigma_star %*% t(cur_Lambda_mat)
  cur_Sigma_star <- cur_Sigma_star + cur_Sigma_u_mat + cur_Theta_mat
  cur_Sigma_star <- (cur_Sigma_star + t(cur_Sigma_star))/2
  
  return(0.5 * (logdet(cur_Sigma_star)
                + tr(solve(cur_Sigma_star, empirical_CovY_list[[i]]))))
  
}

# for all i = 1, ..., N subjects
FABOUP_negloglik_n <- function(m_params_vec, calc_sigma = TRUE) {
  theta_param_length = length(as.vector(theta_boup))
  theta_current = matrix(m_params_vec[(3*k+1):(3*k+theta_param_length)],
                         nrow = nrow(theta_boup))
  
  if (calc_sigma == TRUE){
    sigma_current = calc_sigma(theta_current)
    m_params_long <- c(m_params_vec, log(sigma_current))
  }else{
    m_params_long <- m_params_vec
  }
  
  # check mean reverting OUP
  theta_ok <- as.numeric(all(Re(eigen(theta_current)$values) > 0))
  
  if (theta_ok == 0){ 
    cat('warning: eigenvalues of theta do not have positive real parts \n')
    return(999999)
  }else{
    sum(unlist(lapply(1:n_subj, FUN = function(i){ FABOUP_negllk_i(m_params_long, i) })))
  }
}


