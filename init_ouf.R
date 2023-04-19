################################################################################
##      Use empirical approach to INITIALIZE parameter estimates before       ##
##    fitting the dynamic Ornstein-Uhlenbeck factor model to simulated data   ##
################################################################################


################################################################################
## Restructure data
################################################################################

# Define sample size and number of measurements per individual
n_subj <- length(unique(sim_dat$user.id))
measurement_n <- as.numeric(table(sim_dat$user.id))
max_ni <- max(measurement_n)

# Set up measurement times for each individual
measurement_times <- matrix(nrow = n_subj, ncol = max_ni)
for (i in 1:n_subj){
  ni <- measurement_n[i]
  measurement_times[i,1:ni] <- sim_dat$meas_time[sim_dat$user.id == i]
}

# Set up matrix of observed longitudinal outcomes
all_person_Y <- matrix(NA, nrow = n_subj, ncol = k*max_ni)
for (i in 1:n_subj){
  ni <- measurement_n[i]
  nrow_CovY <- ni * k
  obs_Y <- sim_dat %>%
    filter(user.id == i) %>%
    dplyr::select(c('Y1', 'Y2', 'Y3', 'Y4'))
  all_person_Y[i, 1:nrow_CovY] <- matrix(c(t(obs_Y)), nrow = 1)
}

# Calculate sample covariance matrix for each subject
#   center observed data
colMeansY_matrix = matrix(colMeans(all_person_Y, na.rm = T),
                          nrow = nrow(all_person_Y), ncol = ncol(all_person_Y),
                          byrow =  T)
all_person_Y_centered = all_person_Y - colMeansY_matrix

empirical_CovY_list <- lapply(1:n_subj, FUN = function(i){
  person_data_centered <- all_person_Y_centered[i,!is.na(all_person_Y_centered[i,])]
  person_data_centered %*% t(person_data_centered)
})


# Set dimensions of OUP parameter values & factor model parameters
theta_boup <- diag(x = 1, nrow = p)
sigma_boup <- diag(x = 1, nrow = p)
Lambda <- matrix(0, nrow = k, ncol = p)
Lambda[nonzero == 1] <- 1
Sigma_u <- diag(k)
Sigma_e <- diag(k)
c_star <- rep(1, p)


################################################################################
## Estimate initial parameter values
################################################################################

cat('---- Initializing Lambda ---- \n')

# Start by initializing the loadings: fit a standard cross-sectional factor
#  model that ignores the longitudinal aspect of our data
init_fac_mod <- 'X1  =~ Y1 + Y2
                   X2 =~ Y3 + Y4'  # structure of initial factor model
fac_mod <- lavaan::cfa(init_fac_mod, data = sim_dat, std.lv = TRUE)
std_loadings <- lavInspect(fac_mod, "std")$lambda # extract loadings


# Next, initialize the elements of Sigma_u and Sigma_e

cat('---- Initializing Sigma_u and Sigma_e ---- \n')

# Compute the factor scores from the previously fitted model
eta_pred <- lavPredict(fac_mod, type = "lv", method = "regression")
# And add factor scores to long version of data
colnames(eta_pred) <- paste0('eta', 1:p)
sim_dat <- cbind(sim_dat, eta_pred)

# Then, fit linear mixed models (without intercepts)
lmm1 <- lme(Y1 ~ eta1 - 1, random = ~ 1 | user.id,
            data = sim_dat, method = 'REML')
lmm2 <- lme(Y2 ~ eta1 - 1, random = ~ 1 | user.id,
            data = sim_dat, method = 'REML')
lmm3 <- lme(Y3 ~ eta2 - 1, random = ~ 1 | user.id,
            data = sim_dat, method = 'REML')
lmm4 <- lme(Y4 ~ eta2 - 1, random = ~ 1 | user.id,
            data = sim_dat, method = 'REML')


# Combine parameter estimates into vectors of initial values
prelim_lambda_vec <- c(fixef(lmm1), fixef(lmm2), fixef(lmm3), fixef(lmm4))
prelim_sigma_u_vec <- c(as.numeric(VarCorr(lmm1)[1]),
                        as.numeric(VarCorr(lmm2)[1]),
                        as.numeric(VarCorr(lmm3)[1]),
                        as.numeric(VarCorr(lmm4)[1]))
prelim_theta_vec <- c(as.numeric(VarCorr(lmm1)[2]),
                      as.numeric(VarCorr(lmm2)[2]),
                      as.numeric(VarCorr(lmm3)[2]),
                      as.numeric(VarCorr(lmm4)[2]))
# Since we model many of these parameters on the log scale, we would like to
#  prevent the initial values from being too close to 0 to avoid logged values
#  getting stuck near -infty
prelim_lambda_vec <- pmax(abs(prelim_lambda_vec), 0.1)*sign(prelim_lambda_vec)
prelim_sigma_u_vec <- pmax(prelim_sigma_u_vec, 0.1)
prelim_theta_vec <- pmax(prelim_theta_vec, 0.1)


# Finally, we initialize the OU process parameters

cat('---- Initializing OU process parameters ---- \n')

# Select the predicted factor scores and convert to wide format
boup_traj <- sim_dat %>%
  dplyr::select(c(user.id, paste0('eta', 1:p)))
boup_traj_wide <- matrix(NA, nrow = nrow(all_person_Y),
                         ncol = p*max_ni)
for (i in 1:n_subj){
  cur_dat <- t(boup_traj %>% filter(user.id == i) %>% dplyr::select(-c(user.id)))
  cur_mat <- matrix(as.matrix(cur_dat), ncol = 1)
  boup_traj_wide[i, 1:(p*measurement_n[i])] <- t(cur_mat)
}

# Then calculate the empirical covariance matrix for observed etas for each subj.
colMeans_matrix = matrix(colMeans(boup_traj_wide, na.rm = T),
                         nrow = nrow(boup_traj_wide),
                         ncol = ncol(boup_traj_wide), byrow =  T)
all_data_centered = boup_traj_wide - colMeans_matrix
empirical_Cov_list <- lapply(1:n_subj, FUN = function(i){
  person_data_centered <- all_data_centered[i,!is.na(all_data_centered[i,])]
  person_data_centered %*% t(person_data_centered)
})

# To estimate the OU process parameters from the factor scores, we define the OU
#  process likelihood, but include an additional white noise error term to soak
#  up some of the extra noise the predicted factor scores
# Note that sigma is estimated on the log scale
# This likelihood is for a single individual
ETA_neg_loglik_i <- function(params, i) {
  #cat('person #:', i, '\n')
  theta_current = matrix(params[1:(p^2)], nrow = nrow(theta_boup))
  sigma_current = diag(x = exp(params[(p^2+1):(length(params)-1)]),
                       nrow = length(exp(params[(p^2+1):(length(params)-1)])))
  gamma_current = exp(params[length(params)])^2 # noise term (on log sqrt scale)
  meas_times_i <- measurement_times[i,!is.na(measurement_times[i,])]
  ni <- length(meas_times_i)
  
  theta_t = t(theta_current)
  sigma2_vec = matrix(sigma_current %*% t(sigma_current), ncol = 1)
  kron_sum_theta = kronsum(theta_current, theta_current)
  
  # calculate precision matrix for OU process
  cur_Omega_i <- calc_precision_cpp(kron_sum_theta = kron_sum_theta,
                                    theta = theta_current, theta_t = theta_t,
                                    sigma2_vec = sigma2_vec, ni = ni,
                                    times = meas_times_i)
  # cur_Psi_i <- solve(cur_Omega_i)
  # Sometimes the precision matrix cannot be inverted, so in that case add a
  #  small amount to diagonal (since this is just for initialization of params)
  cur_Psi_i <- try(solve(cur_Omega_i), silent = TRUE)
  if (class(cur_Psi_i)[[1]] == 'try-error'){
    cat('      Warning: Omega is not invertible for i =', i, '\n')
    diag(cur_Omega_i) <- diag(cur_Omega_i) + 0.001
    cur_Psi_i <- try(solve(cur_Omega_i), silent = TRUE)
    cat('      Adding 0.001 to the diagonal to try to make Omega intervible \n')
  }
  cur_Psi_i <- (cur_Psi_i + t(cur_Psi_i))/2

  # Add white noise (gamma)
  cur_Psi_i <- cur_Psi_i + diag(gamma_current, nrow(cur_Psi_i))
  0.5 * (logdet(cur_Psi_i) + tr(solve(cur_Psi_i, empirical_Cov_list[[i]])))
}

# This is the likelihood for many individuals
ETA_negloglik_n <- function(params) {
  # cat('current params: ', params, '\n')
  theta_current = matrix(params[1:(p^2)], nrow = nrow(theta_boup))
  
  # check that theta corresponds to mean reverting OU process
  # note: conditions checked here are slightly stronger than needed for mean-reverting
  check_diag_pos <- as.numeric(all(diag(theta_current) > 0))
  check_det_pos <- as.numeric(det(theta_current) > 0)
  theta_ok <-  check_diag_pos * check_det_pos
  
  if (theta_ok == 0){ 
    cat('Warning: eigenvalues of theta do not have positive real parts \n')
    return(999999)
  }else{
    nllk <- sum(unlist(lapply(1:n_subj, FUN = function(i){ ETA_neg_loglik_i(params, i) })))
    return(nllk)
  }
}

# Finally, estimate the OU process parameters from the factor scores
theta_start <- matrix(1, p, p); diag(theta_start) <- 2
sigma_start <- rep(log(2), p)
gamma_start <- log(sqrt(2))
start <- c(theta_start, sigma_start, gamma_start)


start_time <- Sys.time()
lower_bnds <- matrix(-Inf, p, p); diag(lower_bnds) <- 1e-4
fit1 <- nlminb(start = start, objective = ETA_negloglik_n,
               control = list(iter.max = 200),
               lower = c(lower_bnds, rep(-Inf, p+1))) 
end_time <- Sys.time()
run_time <- as.numeric(end_time) - as.numeric(start_time)

# Update initial estimates of the OU process parameters
theta_boup <- matrix(fit1$par[1:length(theta_boup)], nrow = p, ncol = p)
sigma_boup <- diag(x = exp(fit1$par[(length(theta_boup) + 1):(length(theta_boup) + nrow(sigma_boup))]),
                   nrow = nrow(sigma_boup))

# Set bounds on initial values because sometimes empirical initialization
# can be unstable
diag(theta_boup) <- pmin(diag(theta_boup), 7)
lowerbd_sigma <- 0.001
diag(sigma_boup) <- pmax(diag(sigma_boup), 0.001)

# Finally, we may need to rescale the OU parameters slightly to ensure that they
#  still correspond to a mean-reverting OU process.
# If new theta does not correspond to a mean-reverting OUP, then make sure to
#    scale down the off-diagonals too to maintain a positive determinant.
# First calculate product of off diagonals
off_diag_prod <- prod(theta_boup[row(theta_boup) != col(theta_boup)])
# Then calculate product of diagonals.
diag_prod <- prod(diag(theta_boup))
# If determinant is not positive, then shrink the offdiagonal elements
if (off_diag_prod >= diag_prod){
  scale_factor <- diag_prod / off_diag_prod
  # Rescale diagonal elements (99% is just to ensure determinant is not exactly 0)
  new_offdiag <- theta_boup[row(theta_boup) != col(theta_boup)] * scale_factor^{1/nrow(theta_boup)}*0.99
  theta_boup[row(theta_boup) != col(theta_boup)] <- new_offdiag
}


# Calculate scaling constants: these are used to impose our identifiability
#  constraint that the factors have a stationary variance of 1
c_start <- rep(1, nrow(theta_boup))
c_fit <- try(calc_constants(c_start, theta_boup, sigma_boup))
if (class(c_fit) == 'try-error'){
  cat('   --- Warning: empirical initialization of OUP parameters did not work --- \n   using default starting values instead \n')
  # default to a standard theta and sigma pair...
  standard_theta <- matrix(0.5, nrow(theta_boup), ncol(theta_boup))
  diag(standard_theta) <- 1
  theta_boup <- standard_theta * sign(theta_boup + 1e-9)
  sigma_boup <- diag(x = 1, nrow = nrow(sigma_boup))
  c_fit <- calc_constants(c_start, theta_boup, sigma_boup)
}
c_vec <- c_fit$c_vec
# These are the rescaled values that satisfy the identifiability constraint
theta_boup <- update_theta(theta_boup, diag(x = c_vec, nrow = length(c_vec)))
sigma_boup <- update_sigma(sigma_boup, diag(x = c_vec, nrow = length(c_vec)))

# To improve initialization of BOUP parameters, compare loglikelihood for the
#  empirical initial values and default initial values, then pick the best one
standard_theta <- matrix(0.5, nrow(theta_boup), ncol(theta_boup))
diag(standard_theta) <- 1
default_theta_boup <- standard_theta * sign(theta_boup + 1e-9)
default_sigma_boup <- diag(x = 1, nrow = nrow(sigma_boup))

# Log likelihood is calculated using initial measurement submod parameters
Lambda <- matrix(0, nrow = k, ncol = p)
Lambda[nonzero==1] <- prelim_lambda_vec
Sigma_u <- diag(prelim_sigma_u_vec)
Sigma_e <- diag(prelim_theta_vec)

empirical_negllk <- BOUP_negloglik_n(params = c(theta_boup, log(diag(sigma_boup))),
                                     c_vec = c_vec)
c_vec_default <- calc_constants(c_start, default_theta_boup, default_sigma_boup)$c_vec
default_negllk <- BOUP_negloglik_n(params = c(default_theta_boup, log(diag(default_sigma_boup))),
                                   c_vec = c_vec_default)
if (empirical_negllk > default_negllk){
  print('default')
  theta_boup <- default_theta_boup
  sigma_boup <- default_sigma_boup
}else{
  print('empirical')
}



################################################################################
## Collect initial parameter estimates
################################################################################

names(prelim_lambda_vec) <- NULL
init_params = list(lambda_vec = prelim_lambda_vec, # nonzero elements of Lambda
                   sigma2_u_vec = prelim_sigma_u_vec, # diag of Sigma_u
                   sigma2_e_vec = prelim_theta_vec, # diag of Sigma_e
                   theta_ou = theta_boup, # OU parameters
                   sigma_ou = sigma_boup) 

cat('Initial parameter estimates are stored in the list called "init_params"')






