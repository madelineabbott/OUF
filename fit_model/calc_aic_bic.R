################################################################################
###               Calculate AIC and BIC for our OU factor model              ###
################################################################################

# We account for the fact that our identifiability constraints reduce the
#  number of free parameters in our model by p (the dimension of the OU process)

################################################################################
# Define the -log-likelihood function for calculate AIC and BIC
################################################################################

# This version of the likelihood include the constant term (which is omitted
#  from the versions used during estimation)

# -log-lik for single subject
exact_negllk_i <- function(m_params_vec, i){
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
  # update Sigma_e
  log_sigma_e_vec <- m_params_vec[(2*k+1):(3*k)]
  cur_Sigma_e <- diag((exp(log_sigma_e_vec))^2)
  cur_Sigma_e_mat <- kronecker(diag(ni), cur_Sigma_e)
  # update theta_ou
  theta_param_length = length(as.vector(theta_ou))
  theta_current = matrix(m_params_vec[(3*k+1):(3*k+theta_param_length)], nrow = nrow(theta_ou))
  # update sigma_ou
  sigma_current = diag(x = exp(m_params_vec[(3*k+theta_param_length+1):length(m_params_vec)]),
                       nrow = p)
  # update OU covariance matrix
  theta_t = t(theta_current)
  sigma2_vec = matrix(sigma_current %*% t(sigma_current), ncol = 1)
  kron_sum_theta = kronsum(theta_current, theta_current)
  ni = length(c(meas_times_i))
  cur_Omega_i <- calc_precision_cpp(kron_sum_theta = kron_sum_theta,
                                    theta = theta_current, theta_t = theta_t,
                                    sigma2_vec = sigma2_vec,
                                    ni = ni, times = c(meas_times_i))
  cur_Psi_i <- solve(cur_Omega_i)
  cur_Psi_i <- (cur_Psi_i + t(cur_Psi_i))/2

  # update covariance matrix for Yi
  cur_Sigma_star <- cur_Lambda_mat %*% cur_Psi_i
  cur_Sigma_star <- cur_Sigma_star %*% t(cur_Lambda_mat)
  cur_Sigma_star <- cur_Sigma_star + cur_Sigma_u_mat + cur_Sigma_e_mat
  cur_Sigma_star <- (cur_Sigma_star + t(cur_Sigma_star))/2
  
  A <- 0.5 * (logdet(cur_Sigma_star)
              + tr(solve(cur_Sigma_star, empirical_CovY_list[[i]])))
  constant_term <- k * ni / 2 * log(2*pi)
  
  return(A + constant_term)
  
}


# -log-lik for all i = 1, ..., N subjects
exact_negloglik_n <- function(m_params_vec, calc_sigma = TRUE) {

  theta_param_length = length(as.vector(theta_ou))
  theta_current = matrix(m_params_vec[(3*k+1):(3*k+theta_param_length)],
                         nrow = nrow(theta_ou))
  
  if (calc_sigma == TRUE){
    sigma_current = calc_sigma(theta_current)
    m_params_long <- c(m_params_vec, log(sigma_current))
  }else{
    m_params_long <- m_params_vec
  }
  
  theta_ok <- as.numeric(all(Re(eigen(theta_current)$values) > 0))
  
  if (theta_ok == 0){ 
    cat('OH NO! eigenvalues of theta do not have positive real parts \n')
    return(999999)
  }else{
    sum(unlist(lapply(1:n_subj, FUN = function(i){ exact_negllk_i(m_params_long, i) })))
  }
}

################################################################################
# Calculate AIC and BIC
################################################################################

# Load saved parameter estimates
cur_results <- try(read.csv(paste0(res_wd, 'data/ouf_fit_', 'v', v, '_p', p,
                                   '_g', g, '.csv')),
                   silent = TRUE)

# Initialize parameter estimates at final values
final_est <- cur_results %>%
  filter(faboup_iter == max(cur_results$faboup_iter)) %>%
  mutate(param_hat = case_when(variable %in% c('lambda_1',
                                               paste0(paste0('sigma_ou_', 1:p), 1:p)) ~ log(est),
                               variable %in% c(paste0('sigma2_u_', 1:k),
                                               paste0('sigma2_e_', 1:k)) ~ log(sqrt(est)),
                               TRUE ~ est))

init <- final_est$param_hat[which(final_est$variable %in% c(paste0('lambda_', 1:k),
                                                            paste0('sigma2_u_', 1:k),
                                                            paste0('sigma2_e_', 1:k),
                                                            paste0(paste0('theta_ou_', 1:p),
                                                                   sort(rep(1:p, p))),
                                                            paste0(paste0('sigma_ou_', 1:p), 1:p)))]

negloglik <- exact_negloglik_n(m_params_vec = init, calc_sigma = FALSE)



# Calculate model selection metrics (lower is better)

# AIC = 2k - 2log(L) 
cur_aic <- 2*(length(init)-p) + 2*negloglik # num. of FREE parameters - constraints
# BIC = log(n)k - 2log(L)
cur_bic <- 2*log(n_subj)*(length(init)-p) + 2*negloglik # num. of FREE parameters - constraints

gof_res <- data.frame(mod_p = p,
                      gof_metric = c('negloglik', 'AIC', 'BIC'),
                      metric_value = c(negloglik, cur_aic, cur_bic))

write.csv(gof_res,
          paste0(res_wd, 'data/model_selection/compare_fit_v', v, '_p', p,
                 '_g', g, '.csv'),
          row.names = FALSE)




