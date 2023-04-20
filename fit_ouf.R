################################################################################
##        Fit dynamic Ornstein-Uhlenbeck factor model to simulated data       ##
################################################################################

# To fit this model, will need to install C++ code (contains gradient formulas)
#devtools::install_github("https://github.com/madelineabbott/FABOUP_fast.git", ref="main")

# Set working directory
wd <- ''

# Source some functions used for estimation
source(file = paste0(wd, 'functions_ouf.R'))
# Source Cpp code for calculating OU process precision matrix
Rcpp::sourceCpp(paste0(wd, 'ou_precision_matrix.cpp')) 

# Load simulated data: contains 4 longitudinal outcomes, Y1, ..., Y4
sim_dat <- read.csv(paste0(wd, 'sim_dat.csv'))

# Assume that we want to fit an OUF model with 2 latent factors
k <- 4 # set number of longitudinal outcomes
p <- 2 # set number of latent factors
cap_theta <- 7 # upper bound for theta, based on gap between meas. occasions
# We also assume that each of the measured longitudinal outcomes load onto only
#  one of the latent factors.  Define indicator matrix of same size as loadings
#  matrix that defines the location of the non-zero loadings
nonzero <- matrix(c(1, 1, 0, 0, 0, 0, 1, 1), ncol = p)

# Parameters can either be initialized empirically or with user-defined values.
#  Specify the approach to be used here:
init_method <- 'empirical' # other option is 'user_def'

# Number of cores--used for parallel calculation of gradients
n_cores <- 4


set.seed(22210)

################################################################################
## Initialize parameters
################################################################################

# Parameters can either be initialized using user-defined values or empirically.

if (init_method == 'empirical'){
  # To initialize parameter EMPIRICALLY, source the R file below. (This may take a
  #  few minutes to run.  Warnings may be printed but can generally be ignored.)
  source(paste0(wd, 'init_ouf.R')) 
  init_params # list of initialize parameters
}else if (init_method == 'user_def'){
  # To initialize parameters with USER-DEFINED VALUES, update the values in the
  #  list below.
  init_params = list(lambda_vec = c(1, 1, 1, 1), # nonzero elements of Lambda
                     sigma2_u_vec = c(0.5, 0.5, 0.5, 0.5), # diag of Sigma_u
                     sigma2_e_vec = c(0.5, 0.5, 0.5, 0.5), # diag of Sigma_e
                     theta_ou = matrix(c(2, 1, 1, 2), 2, 2), # OU parameter
                     sigma_ou = diag(c(1, 1))) # OU parameter
  init_params
}

# Once initial parameters have been defined, we'll update the matrices used when
#  fitting the model.
theta_ou <- init_params$theta_ou
sigma_ou <- init_params$sigma_ou

# First, update OU process covariance matrix using initial OU parameters
Psi_list <- lapply(1:n_subj, FUN = function(i){
  meas_times_i <- measurement_times[i,!is.na(measurement_times[i,])]
  theta_t = t(theta_ou)
  sigma2_vec = matrix(sigma_ou %*% t(sigma_ou), ncol = 1)
  kron_sum_theta = kronsum(theta_ou, theta_ou)
  cur_Omega_i <- calc_precision_cpp(kron_sum_theta = kron_sum_theta,
                                    theta = theta_ou, theta_t = theta_t,
                                    sigma2_vec = sigma2_vec,
                                    ni = length(meas_times_i),
                                    times = meas_times_i)
  cur_Psi_i <- solve(cur_Omega_i)
  cur_Psi_i <- (cur_Psi_i + t(cur_Psi_i))/2
  cur_Psi_i
})

################################################################################
## Define functions used to update each block of parameters
################################################################################

# Block update for measurement submodel parameters given fixed structural
#  submodel parameters
fitFA <- function(initial, gradtol, steptol){
  # initialize parameter values
  lambda_vec_0 <- initial[1:k]
  sigma_u_vec_0 <- initial[(k+1):(2*k)]
  theta_vec_0 <- initial[(2*k+1):(3*k)]
  m_params_vec <- c(log(lambda_vec_0[1]), lambda_vec_0[2:k], # lambda
                    log(sqrt(sigma_u_vec_0)),
                    log(sqrt(theta_vec_0)))
  
  start_time <- Sys.time()
  fit <- nlm(f = FA_negllk, p = m_params_vec, Psi_list = Psi_list,
             check.analyticals = FALSE, stepmax = 10, print.level = 0,
             gradtol = gradtol, steptol = steptol, iterlim = 200)
  end_time <- Sys.time()
  run_time <- as.numeric(end_time) - as.numeric(start_time)
  
  # convert from model-scale to original-scale
  params_hat <- c(exp(fit$estimate[1]), fit$estimate[2:k],
                  exp(fit$estimate[(k+1):(2*k)])^2,
                  exp(fit$estimate[(2*k+1):(3*k)])^2)
  cur_gradient <- c(fit$gradient)
  
  # collect things we'd like to save from this block update
  current_fit <- data.frame(variable = c(paste0('lambda_', 1:k),
                                         paste0('sigma2_u_', 1:k),
                                         paste0('sigma2_e_', 1:k)),
                            est = c(params_hat), # updated parameter estimates
                            grad = c(cur_gradient), # final value of gradient
                            obj = fit$minimum, # objective function
                            iter = fit$iterations, # number of iteration required
                            convg = fit$code, # convergence status
                            run_time = run_time) # computation time
  return(current_fit)
}

# Block update for structural submodel parameters given fixed measurement
#  submodel parameters
fitBOUP <- function(initial, c_vec, gradtol, steptol){
  theta_param_length <- length(as.vector(theta_ou))
  # initial parameters should be on model scale (sigmas need to be on log scale)
  start <- c(initial[1:theta_param_length],
             log(initial[(theta_param_length+1):(theta_param_length+nrow(theta_ou))])) 
  
  start_time <- Sys.time()
  # define lower bounds for theta parameters
  lower_bnds <- matrix(-Inf, nrow = nrow(theta_ou), ncol = ncol(theta_ou))
  diag(lower_bnds) <- 1e-4
  lower_bnds <- c(lower_bnds)
  # define upper bounds for theta parameters
  # note: upper bounds for theta depend on gap times between meas. occasions
  upper_bnds <- matrix(Inf, nrow = nrow(theta_ou), ncol = ncol(theta_ou))
  diag(upper_bnds) <- cap_theta
  upper_bnds <- c(upper_bnds)
  
  fit2 <- nlminb(start = start, objective = BOUP_negloglik_n,
                 c_vec = c_vec, control = list(iter.max = 200),
                 lower = c(lower_bnds, rep(-Inf, nrow(theta_ou))),
                 upper = c(upper_bnds, rep(Inf, nrow(theta_ou))))
  end_time <- Sys.time()
  run_time <- as.numeric(end_time) - as.numeric(start_time)
  
  # collect things we'd like to save from this block update
  current_fit <- data.frame(variable = c(paste0('lambda_', 1:k),
                                         paste0('sigma2_u_', 1:k),
                                         paste0('sigma2_e_', 1:k),
                                         paste0(paste0('theta_ou_', 1:p),
                                                sort(rep(1:p, p))),
                                         paste0(paste0('sigma_', 1:p), 1:p)),
                            # updated parameter estimates
                            est = c(c(Lambda[nonzero==1]),
                                    diag(Sigma_u),
                                    diag(Sigma_e),
                                    fit2$par[1:theta_param_length],
                                    exp(fit2$par[(theta_param_length+1):(theta_param_length+nrow(theta_ou))])),
                            # gradient
                            grad = c(rep(NA, 3*k), rep(NA, theta_param_length + nrow(theta_ou))), 
                            # final value of objective function
                            obj = fit2$objective, 
                            # total number of iteration required 
                            iter = fit2$iterations, 
                            # convergence status 
                            convg = fit2$convergence, 
                            # computation time
                            run_time = run_time)
  
  return(current_fit)
}

################################################################################
## Set up simulated data 
################################################################################

# This set up is only needed if user-defined values are specified.  If empirical
#  initialization is used, then all of this set up is done in the init_ouf.R file.

if (init_method == 'user_def'){
  # Define sample size and max number of measurement occasions
  n_subj <- length(unique(sim_dat$user.id))
  measurement_n <- as.numeric(table(sim_dat$user.id))
  max_ni <- max(measurement_n)
  
  ### set up measurement occasions
  measurement_times <- matrix(nrow = n_subj, ncol = max_ni)
  for (i in 1:n_subj){
    ni <- measurement_n[i]
    measurement_times[i,1:ni] <- sim_dat$end_hrs[sim_dat$user.id == i]
  }
  
  # observed longitudinal outcome
  all_person_Y <- matrix(NA, nrow = n_subj, ncol = k*max_ni)
  for (i in 1:n_subj){
    ni <- measurement_n[i]
    nrow_CovY <- ni * k
    obs_Y <- sim_dat %>%
      filter(user.id == i) %>%
      dplyr::select(all_of(emotion_items))
    all_person_Y[i, 1:nrow_CovY] <- matrix(c(t(obs_Y)), nrow = 1)
  }
  
  # calculate sample covariance matrix for each subject
  #   center observed data
  colMeansY_matrix = matrix(colMeans(all_person_Y, na.rm = T),
                            nrow = nrow(all_person_Y), ncol = ncol(all_person_Y),
                            byrow =  T)
  all_person_Y_centered = all_person_Y - colMeansY_matrix
  
  empirical_CovY_list <- lapply(1:n_subj, FUN = function(i){
    person_data_centered <- all_person_Y_centered[i,!is.na(all_person_Y_centered[i,])]
    person_data_centered %*% t(person_data_centered)
  })
  
  
  ### 2. Set dimensions of OUP parameter values & factor model parameters
  nonzero <- nonzero_for_fitting
  theta_ou <- diag(x = 1, nrow = p)
  sigma_ou <- diag(x = 1, nrow = p)
  Lambda <- matrix(0, nrow = k, ncol = p)
  Lambda[nonzero == 1] <- 1
  Sigma_u <- diag(k)
  Sigma_e <- diag(k)
  c_star <- rep(1, p)
}

################################################################################
## Estimate parameters using block coordinate descent algorithm
################################################################################

# Rename initial parameter estimates
prelim_lambda_vec <- init_params$lambda_vec
prelim_sigma_u_vec <- init_params$sigma2_u_vec
prelim_theta_vec <- init_params$sigma2_e_vec
theta_ou <- init_params$theta_ou
sigma_ou <- init_params$sigma_ou


res <- matrix(nrow = 0, ncol = 8,
              dimnames = list(NULL, c('variable', 'est', 'obj',
                                      'iter', 'convg', 'run_time',
                                      'faboup_iter', 'faboup_convg')))
max_r <- 80 # maximum number of iterations across blocks
for (r in 1:max_r){
  cat('  block-wise update: iteration number', r, '\n')
  # initialize
  if(r == 1){
    Lambda <- matrix(0, nrow = k, ncol = p)
    Lambda[nonzero==1] <- prelim_lambda_vec
    Sigma_u <- diag(prelim_sigma_u_vec)
    Sigma_e <- diag(prelim_theta_vec)
    
    # save initial estimates
    res <- data.frame(variable = c(paste0('lambda_', 1:k),
                                   paste0('sigma2_u_', 1:k),
                                   paste0('sigma2_e_', 1:k),
                                   paste0(paste0('theta_ou_', 1:p),
                                          sort(rep(1:p, p))),
                                   paste0(paste0('sigma_', 1:p), 1:p)),
                      # current parameter values
                      est = c(c(Lambda[nonzero==1]), diag(Sigma_u),
                              diag(Sigma_e), c(theta_ou), diag(sigma_ou)), 
                      grad = NA, # gradient
                      obj = NA, # final value of objective function
                      iter = NA, # total number of iterations required 
                      convg = NA, # convergence status
                      run_time = NA,
                      faboup_iter = 0, # covergence of block update
                      faboup_convg = 0) # covergence of entire coord descent alg.
    
    # Save results
    # write.csv(res %>% mutate(update_type = 'block',
    #                          max_faboup_iter = max_r,
    #                          faboup_convg = 0),
    #           paste0(wd, 'OUF_fit.csv'))
  }
  
  #### Update measurement submodel parameters ####
  # unknown meas. submod. parameters: initial (on original scale, not log scale)
  # fixed struct. submod. parameters: Psi_list, theta_ou, sigma_ou
  FA_initial <- c(Lambda[nonzero == 1], diag(Sigma_u), diag(Sigma_e))
  cl <- parallel::makeCluster(n_cores)
  doParallel::registerDoParallel(cl)
  updateFA <- fitFA(initial = FA_initial,
                    gradtol = max(1e-4/(10^r), 1e-8),
                    steptol = max(1e-4/(10^r), 1e-8))
  doParallel::stopImplicitCluster()
  parallel::stopCluster(cl)
  res <- rbind(res, head(updateFA, 3*k) %>% mutate(faboup_iter = r,
                                                   faboup_convg = 0))

  # Update measurement submodel parameters
  Lambda[nonzero == 1] <- updateFA$est[1:k]
  Sigma_u <- diag(updateFA$est[(k+1):(2*k)])
  Sigma_e <- diag(updateFA$est[(2*k+1):(3*k)])
  
  #### Update structural submodel parameters ####
  # unknown struct. submod. parameters: initial (on original scale, not log scale)
  # fixed meas. submod. parameters: Lambda, Sigma_u, Sigma_e
  BOUP_initial <- c(theta_ou, diag(sigma_ou))
  updateBOUP <- fitBOUP(initial = BOUP_initial,
                        c_vec = c_vec,
                        gradtol = max(1e-4/(10^r), 1e-8),
                        steptol = max(1e-4/(10^r), 1e-8))
  res <- rbind(res, tail(updateBOUP, length(BOUP_initial))
               %>% mutate(faboup_iter = r, faboup_convg = 0))

  # Update structural submodel parameters
  theta_ou <- matrix(updateBOUP$est[(3*k+1):(3*k+length(theta_ou))],
                       nrow(theta_ou), nrow(theta_ou))
  sigma_ou <- diag(x = updateBOUP$est[(1+3*k+length(theta_ou)):(p+3*k+length(theta_ou))],
                     nrow = p)
  
  # Identifiability assumption: calculate scaling constants and rescale updated
  #  structural submodel parameters (OU process parameters)
  c_fit <- calc_constants(c_vec, theta_ou, sigma_ou)
  c_vec <- c_fit$c_vec
  # These are the rescaled values that satisfy the identifiability constraint
  theta_ou <- update_theta(theta_ou, diag(x = c_vec, nrow = length(c_vec)))
  sigma_ou <- update_sigma(sigma_ou, diag(x = c_vec, nrow = length(c_vec)))
  
  # Optional: save results
  # write.csv(res %>% mutate(update_type = 'block',
  #                          max_faboup_iter = max_r,
  #                          faboup_convg = 0),
  #           paste0(wd, 'OUF_fit.csv'))
  
  # Update covariance matrix for latent factors
  Psi_list <- lapply(1:n_subj, FUN = function(i){
    meas_times_i <- measurement_times[i,!is.na(measurement_times[i,])]
    theta_t = t(theta_ou)
    sigma2_vec = matrix(sigma_ou %*% t(sigma_ou), ncol = 1)
    kron_sum_theta = kronsum(theta_ou, theta_ou)
    ni = length(c(meas_times_i))
    Psi_i <- solve(calc_precision_cpp(kron_sum_theta = kron_sum_theta,
                                      theta = theta_ou, theta_t = theta_t,
                                      sigma2_vec = sigma2_vec,
                                      ni = ni, times = c(meas_times_i)))
    Psi_i <- (Psi_i + t(Psi_i))/2
    Psi_i
    
  })
  # Psi_i should have diagonal elements of 1 (so is a correlation matrix) as a
  #  result of the identifiability assumption
  
  # Check for convergence across iterations of block updates
  max_abs_reldiff = 1
  llk_diff = 100
  if (r > 1){
    # convergence criteria 1: maximum relative change in estimates
    cur_estimates <- res$est[which(res$faboup_iter == r)]
    prev_estimates <- res$est[which(res$faboup_iter == r-1)]
    max_abs_reldiff <- max(abs((cur_estimates - prev_estimates)/prev_estimates))
    # convergence criteria 2: change in log likelihood
    cur_llk <- tail(res$obj[which(res$faboup_iter == r)],1)
    prev_llk <- tail(res$obj[which(res$faboup_iter == r-1)],1)
    llk_diff = prev_llk - cur_llk # positive change in llk_diff is good
  }
  if(max_abs_reldiff < 1e-6){
    cat('  ---- convergence achieved! (criteria # 1: small parameter change) ---- \n')
    res <- res %>% mutate(faboup_convg = case_when(faboup_iter == r ~ 1,
                                                   TRUE ~ 0))
    # write.csv(res, paste0(wd, 'OUF_fit.csv'))
    break()
  }else if (llk_diff > 0 & llk_diff < 1e-6){
    cat('  ---- convergence achieved! (criteria # 2: small llk change) ---- \n')
    res <- res %>% mutate(faboup_convg = case_when(faboup_iter == r ~ 2,
                                                   TRUE ~ 0))
    # write.csv(res, paste0(wd, 'OUF_fit.csv'))
    break()
  }
}


################################################################################
## Estimate standard errors based on Hessian from joint likelihood
################################################################################

# Extract final parameter estimates from the block coordinate descent algorithm
final_est <- res %>%
  filter(faboup_iter == max(res$faboup_iter)) %>%
  dplyr::select(c(variable, est)) %>%
  mutate(param_hat = case_when(variable %in% c('lambda_1',
                                               paste0(paste0('sigma_', 1:p), 1:p)) ~ log(est),
                               variable %in% c(paste0('sigma2_u_', 1:k),
                                               paste0('sigma2_e_', 1:k)) ~ log(sqrt(est)),
                               TRUE ~ est))

param_hat <- final_est$param_hat[which(final_est$variable %in% c(paste0('lambda_', 1:k),
                                                      paste0('sigma2_u_', 1:k),
                                                      paste0('sigma2_e_', 1:k),
                                                      paste0(paste0('theta_ou_', 1:p),
                                                             sort(rep(1:p, p)))))]

# Calculate numerical approximation to the Hessian of the joint log likelihood
start.time <- Sys.time()
fit_hess <- optimHess(par = param_hat, fn = FABOUP_negloglik_n)
end.time <- Sys.time()
run.time_h <- as.numeric(end.time) - as.numeric(start.time)

new_se <- sqrt(diag(solve(fit_hess)))


# Results (on transformed scale)
res_se <- data.frame(variable = final_est$variable[1:length(param_hat)],
                     # estimates (on transformed scale)
                     m_est = c(param_hat[1:length(param_hat)]),
                     # standard errors (on transformed scale)
                     m_se = new_se,
                     run_time_sec = run.time_h, 
                     update_type = 'joint')

# Transform from scale used for estimation (e.g., log(lambda_1)) to model scale
res_se <- res_se %>%
  mutate(est = case_when(variable %in% c('lambda_1',
                                         'sigma_ou_11',
                                         'sigma_ou_22') ~ exp(m_est),
                         variable %in% c('sigma2_u_1', 'sigma2_u_2',
                                         'sigma2_u_3', 'sigma2_u_4',
                                         'sigma2_e_1', 'sigma2_e_2',
                                         'sigma2_e_3', 'sigma2_e_4') ~ exp(m_est)^2,
                         TRUE ~ m_est)) %>%
  # Apply Delta rule where needed for standard errors
  mutate(se = case_when(variable %in% c('lambda_1',
                                        'sigma_ou_11',
                                        'sigma_ou_22') ~ sqrt(m_se^2 * exp(2 * m_est)),
                        variable %in% c('sigma2_u_1', 'sigma2_u_2',
                                        'sigma2_u_3', 'sigma2_u_4',
                                        'sigma2_e_1', 'sigma2_e_2',
                                        'sigma2_e_3', 'sigma2_e_4') ~ sqrt(m_se^2 * 4 * exp(4 * m_est)),
                        TRUE ~ m_se)) %>%
  dplyr::select(c(variable, est, se, run_time_sec))

# calculate sigma_ou estimates (standard errors can be estimated w/ bootstrapping)
theta_ou_hat <- matrix(res_se$est[which(res_se$variable %in%
                                          paste0(paste0('theta_ou_', 1:p),
                                                 sort(rep(1:p, p))))], p, p)
sigma_ou_hat <- calc_sigma(theta_ou_hat)
res_se <- rbind(res_se, data.frame(variable = paste0('sigma_ou_', 1:p, 1:p),
                                   est = sigma_ou_hat,
                                   se = rep(NA, p),
                                   run_time_sec = rep(NA, p)))

# Because we have assumed only one non-zero loading per row in Lambda, we can
#  flip the correlation of the latent variables and the signs of the loadings
#  but still have the same model.  If we interpret the latent factors as
#  positive and negative affect, we'd expect them to be negatively correlated--
#  a positive correlation was estimated so we flip the signs for interpretation.

res_se <- res_se %>%
  mutate(est = case_when(variable %in% c('lambda_3', 'lambda_4',
                                         'theta_ou_12', 'theta_ou_21') ~ -est,
                         TRUE ~ est))


cat('---- FINAL ESTIMATES ---- \n')
res_se
# write.csv(res_se, paste0(wd, 'OUF_fit_SEs.csv'))

################################################################################
## Compare estimates to true values
################################################################################

# for comparison, true parameters values are: 
true_values <- data.frame(variable = c(paste0('lambda_', 1:k),
                                       paste0('sigma2_u_', 1:k),
                                       paste0('sigma2_e_', 1:k),
                                       paste0(paste0('theta_ou_', 1:p),
                                              sort(rep(1:p, p))),
                                       paste0(paste0('sigma_ou_', 1:p),
                                              1:p)),
                      true_value = c(1.1577, 1.7366, -0.3948, 1.9742, # loadings
                                    1.1, 1.3, 1.4, 0.9, # diag(Sigma_u)
                                    0.6,  0.5,  0.4,  0.7, # diag(Sigma_e)
                                    1,  3.9095, 0.6139,  5, # theta_ou
                                    1.0365,  2.0261)) # diag(sigma_ou)

res_se_with_truth <- left_join(res_se, true_values)

ggplot(data = res_se_with_truth) +
  geom_point(aes(x = variable, y = est, color = 'Estimate')) +
  geom_errorbar(aes(x = variable, ymin = est -1.96*se, ymax = est + 1.96*se,
                    color = 'Estimate')) +
  geom_point(aes(x = variable, y = true_value, color = 'true value'), ) +
  theme_bw() + labs(x = 'Parameter', y = 'Estimate (95% CI)') +
  scale_color_manual(values = c('black', 'orange')) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))





