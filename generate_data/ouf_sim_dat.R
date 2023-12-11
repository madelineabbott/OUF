################################################################################
##                               Simulate data                                ##
################################################################################

# dat_wd <- [directory containing subfolder called 'data' for saving datasets]
# fun_wd <- [directory containing the R and C++ function files]

# Save true parameter values?
save_true <- TRUE
# Save longitudinal data?
save_data <- TRUE
# Which set of true parameter values?
v <- 1
# Random seed?
g <- 1

# load useful functions
source(file = paste0(fun_wd, "functions_ouf.R"))

# load C++ code for calculating OU precision matrix
Rcpp::sourceCpp(paste0(fun_wd, "ou_precision_matrix.cpp"))

################################################################################
#### STEP 1: DEFINE TRUE PARAMETER VALUES
################################################################################

cat('---- Defining parameter values ---- \n')

# Set parameter values for data-generating structural and measurement submodels.
# We consider models with one, two, or three latent factors.

if (v == 1){ # values similar to those from motivating data
  # OUP parameters
  theta_ou <- matrix(c(1, 4, 0.6, 5),
                       nrow = 2, ncol = 2)
  sigma_ou <- diag(c(1, 2))
  
  # Factor model parameters
  k <- 4 # number of longitudinal outcomes (e.g. EMA items)
  nonzero <- matrix(c(1, 1, 0, 0, 0, 0, 1, 1), ncol = 2) # locations of non-0 loadings
  # Loading matrix
  Lambda <- matrix(c(1.2, 1.8, 0, 0, 0, 0, -0.4, 2), ncol = 2) 
  # Covariance matrix for random intercept
  Sigma_u <- diag(c(1.1, 1.3, 1.4, 0.9)) 
  # Covariance matrix for measurement error
  Sigma_e <- diag(c(0.6, 0.5, 0.4, 0.7)) 
  
} else if (v == 2){ # values s.t. theta* = theta (no need for rescaling to correlation scale)
  # OUP parameters
  theta_ou <- matrix(c(1, 1.8, 0.4, 3),
                       nrow = 2, ncol = 2)
  sigma_ou <- diag(calc_sigma(theta_ou))
  
  # Factor model parameters
  k <- 4 # number of longitudinal outcomes (e.g. EMA items)
  nonzero <- matrix(c(1, 1, 0, 0, 0, 0, 1, 1), ncol = 2) # locations of non-0 loadings
  # Loading matrix
  Lambda <- matrix(c(1.2, 1.8, 0, 0, 0, 0, -0.4, 2), ncol = 2) 
  # Covariance matrix for random intercept
  Sigma_u <- diag(c(1.1, 1.3, 1.4, 0.9)) 
  # Covariance matrix for measurement error
  Sigma_e <- diag(c(0.6, 0.5, 0.4, 0.7)) 
  
} else if (v == 3){ # values s.t. theta* != theta
  # OUP parameters
  theta_ou <- matrix(c(1, 2, 0.5, 5),
                       nrow = 2, ncol = 2)
  sigma_ou <- diag(c(2, 3))
  
  # Factor model parameters
  k <- 4 # number of longitudinal outcomes (e.g. EMA items)
  nonzero <- matrix(c(1, 1, 0, 0, 0, 0, 1, 1), ncol = 2) # locations of non-0 loadings
  # Loading matrix
  Lambda <- matrix(c(1.2, 1.8, 0, 0, 0, 0, -0.4, 2), ncol = 2) 
  # Covariance matrix for random intercept
  Sigma_u <- diag(c(1.1, 1.3, 1.4, 0.9)) 
  # Covariance matrix for measurement error
  Sigma_e <- diag(c(0.6, 0.5, 0.4, 0.7)) 
  
} else if (v == 4){ # true model only has one latent factor
  # OUP parameters
  theta_ou <- matrix(c(0.8), nrow = 1, ncol = 1)
  sigma_ou <- matrix(c(1), nrow = 1, ncol = 1)
  
  # Factor model parameters
  k <- 4 # number of longitudinal outcomes (e.g. EMA items)
  nonzero <- matrix(c(1, 1, 1, 1), ncol = 1) # locations of non-0 loadings
  # Loading matrix
  Lambda <- matrix(c(1.2, 1.8, -0.4, 2), ncol = 1) 
  # Covariance matrix for random intercept
  Sigma_u <- diag(c(1.1, 1.3, 1.4, 0.9)) 
  # Covariance matrix for measurement error
  Sigma_e <- diag(c(0.6, 0.5, 0.4, 0.7)) 
  
} else if (v == 5){ # true model has three latent factors ("low signal")
  # OUP parameters
  theta_ou <- matrix(c(1, 1.8, 0.9,
                         0.4, 3, 1,
                         0.6, 0.9, 1.2), nrow = 3, ncol = 3)
  sigma_ou <- diag(c(1.2, 0.8, 0.4))
  
  # Factor model parameters
  k <- 4 # number of longitudinal outcomes (e.g. EMA items)
  nonzero <- matrix(c(1, 1, 0, 0,
                      0, 0, 1, 0,
                      0, 0, 0, 1), ncol = 3) # locations of non-0 loadings
  # Loading matrix
  Lambda <- matrix(0, nrow = k, ncol = ncol(theta_ou))
  Lambda[nonzero==1] <- c(1.2, 1.8, -0.4, 2)
  # Covariance matrix for random intercept
  Sigma_u <- diag(c(1.1, 1.3, 1.4, 0.9)) 
  # Covariance matrix for measurement error
  Sigma_e <- diag(c(0.6, 0.5, 0.4, 0.7)) 
  
}else if (v == 6){ # true model has three latent factor ("high signal")
  # OUP parameters
  theta_ou <- matrix(c(2, 0.8, 0.7,
                         0.2, 1.1, 0.5,
                         0.4, 0.5, 1.2), nrow = 3, ncol = 3)
  sigma_ou <- diag(c(1.2, 0.8, 0.4))
  
  # Factor model parameters
  k <- 4 # number of longitudinal outcomes (e.g. EMA items)
  nonzero <- matrix(c(1, 1, 0, 0,
                      0, 0, 1, 0,
                      0, 0, 0, 1), ncol = 3) # locations of non-0 loadings
  # Loading matrix
  Lambda <- matrix(0, nrow = k, ncol = ncol(theta_ou))
  Lambda[nonzero==1] <- c(1.2, 1.8, -0.4, 2)
  # Covariance matrix for random intercept
  Sigma_u <- diag(c(1.1, 1.3, 1.4, 0.9)) 
  # Covariance matrix for measurement error
  Sigma_e <- diag(c(0.6, 0.5, 0.4, 0.7)) 
  
}else if (v == 7){ # true model has two latent factors ("low signal")
  # OUP parameters
  theta_ou <- matrix(c(1, 2, 1.5, 5),
                       nrow = 2, ncol = 2)
  # eigen(theta_ou) # check eigenvalues
  sigma_ou <- diag(c(2, 3))
  
  # Factor model parameters
  k <- 4 # number of longitudinal outcomes (e.g. EMA items)
  nonzero <- matrix(c(1, 1, 0, 0, 0, 0, 1, 1), ncol = 2) # locations of non-0 loadings
  # Loading matrix
  Lambda <- matrix(c(1.2, 1.8, 0, 0, 0, 0, -0.4, 2), ncol = 2) 
  # Covariance matrix for random intercept
  Sigma_u <- diag(c(1.1, 1.3, 1.4, 0.9)) 
  # Covariance matrix for measurement error
  Sigma_e <- diag(c(0.6, 0.5, 0.4, 0.7)) 
  
}else if (v == 8){ # true model has two latent factors ("high signal")
  # OUP parameters
  theta_ou <- matrix(c(2, 0.4, 0.5, 4),
                       nrow = 2, ncol = 2)
  # eigen(theta_ou) # check eigenvalues
  sigma_ou <- diag(c(2, 1))
  cov2cor(construct_V(theta_ou, sigma_ou))
  
  # Factor model parameters
  k <- 4 # number of longitudinal outcomes (e.g. EMA items)
  nonzero <- matrix(c(1, 1, 0, 0, 0, 0, 1, 1), ncol = 2) # locations of non-0 loadings
  # Loading matrix
  Lambda <- matrix(c(1.2, 1.8, 0, 0, 0, 0, -0.4, 2), ncol = 2) 
  # Covariance matrix for random intercept
  Sigma_u <- diag(c(1.1, 1.3, 1.4, 0.9)) 
  # Covariance matrix for measurement error
  Sigma_e <- diag(c(0.6, 0.5, 0.4, 0.7)) 
}else if (v == 9){ # for comparison with Tran et al. (2021)
  # OUP parameters
  theta_ou <- matrix(c(1, 4, 0.6, 5),
                       nrow = 2, ncol = 2)
  sigma_ou <- diag(c(1, 2))
  
  # Factor model parameters
  k <- 4 # number of longitudinal outcomes (e.g. EMA items)
  nonzero <- matrix(c(1, 1, 0, 0, 0, 0, 1, 1), ncol = 2) # indices of non-0 loadings
  # Loading matrix
  Lambda <- matrix(c(1.2, 1.8, 0, 0, 0, 0, 0.4, 2), ncol = 2) # all positive
  # Covariance matrix for random intercept
  Sigma_u <- diag(c(1.1, 1.3, 1.4, 0.9)) 
  # Covariance matrix for measurement error
  Sigma_e <- diag(c(0.6, 0.5, 0.4, 0.7)) 
}

################################################################################
#### STEP 2: GENERATE LONGITUDINAL DATA
################################################################################

cat('---- Generating data ---- \n')

set.seed(g)

# Set sample size and maximum number of measurement occasions
n_subj <- 200
max_ni <- 20 

# Set true number of latent factors used for data generation
true_p <- nrow(theta_ou) 

# Generate measurement times (unbalanced data)
measurement_times <- matrix(nrow = n_subj, ncol = max_ni)
measurement_n <- rep(NA, n_subj)
for (i in 1:n_subj){
  ni <- sample(10:max_ni, size = 1)
  measurement_n[i] <- ni
  measurement_times[i,1:ni] <- cumsum(runif(0.1, 2, n = ni))
}

# Calculate OU covariance matrix of latent factors for each subject
Psi_list <- lapply(1:n_subj, FUN = function(i){
  meas_times_i <- measurement_times[i,!is.na(measurement_times[i,])]
  
  theta_t = t(theta_ou)
  sigma2_vec = matrix(sigma_ou %*% t(sigma_ou), ncol = 1)
  kron_sum_theta = kronsum(theta_ou, theta_ou)
  ni = length(c(meas_times_i))
  Omega_i <- calc_precision_cpp(kron_sum_theta = kron_sum_theta, theta = theta_ou,
                                theta_t = theta_t, sigma2_vec = sigma2_vec,
                                ni = ni, times = c(meas_times_i))
  Psi_i <- solve(Omega_i)
  Psi_i <- (Psi_i + t(Psi_i))/2
  Psi_i
})

# Generate observed values of the latent process (eta) for each subject
all_data <- matrix(nrow = n_subj, ncol = true_p*max_ni)
for (i in 1:n_subj) {
  Psi_i <- Psi_list[[i]]
  person_data <- rmvnorm(n = 1, rep(0,nrow(Psi_i)), Psi_i)
  all_data[i,1:nrow(Psi_i)] <- person_data
}

# Calculate empirical covariance matrix for observed etas for each subj.
colMeans_matrix = matrix(colMeans(all_data, na.rm = T), nrow = nrow(all_data),
                         ncol = ncol(all_data), byrow =  T)
all_data_centered = all_data - colMeans_matrix # center data

empirical_Cov_list <- lapply(1:n_subj, FUN = function(i){
  person_data_centered <- all_data_centered[i,!is.na(all_data_centered[i,])]
  person_data_centered %*% t(person_data_centered)
})


# Generate observed Y values (longitudinal outcome)
# First calculate Y covariance matrix for each subject
CovY_list <- lapply(1:n_subj, FUN = function(i){
  meas_times_i <- measurement_times[i,!is.na(measurement_times[i,])]
  ni <- length(meas_times_i)
  Psi_i <- Psi_list[[i]]
  Lambda_mat <- kronecker(diag(ni), Lambda) %*% Psi_i %*% kronecker(diag(ni), t(Lambda))
  Sigma_u_mat <- kronecker(matrix(1, nrow = ni, ncol = ni), Sigma_u)
  Sigma_e_mat <- kronecker(diag(ni), Sigma_e)
  
  Lambda_mat + Sigma_u_mat + Sigma_e_mat # Cov(Y_i)
})

# And then sample Y_i ~ N(0, Sigma*_i)
all_person_Y <- matrix(NA, nrow = n_subj, ncol = k*max_ni)
for (i in 1:n_subj){
  CovY_i <- CovY_list[[i]]
  all_person_Y[i, 1:nrow(CovY_i)] <- rmvnorm(n = 1,
                                             mean = rep(0, nrow(CovY_i)),
                                             sigma = CovY_i)
}

# Calculate sample covariance matrix of observed long. outcome Y each subject
colMeansY_matrix = matrix(colMeans(all_person_Y, na.rm = T),
                          nrow = nrow(all_person_Y), ncol = ncol(all_person_Y),
                          byrow =  T)
all_person_Y_centered = all_person_Y - colMeansY_matrix # center obs. Y data

empirical_CovY_list <- lapply(1:n_subj, FUN = function(i){
  person_data_centered <- all_person_Y_centered[i,!is.na(all_person_Y_centered[i,])]
  person_data_centered %*% t(person_data_centered)
})


# Rescale theta_ou, sigma_ou, and Lambda (these are the true values targeted by
#  the estimation algorithm)
d <- diag(x = diag(construct_V(theta_ou, sigma_ou))^(1/2),
          nrow(theta_ou), ncol(theta_ou))
d.inv <- solve(d)
# True loadings matrix targeted by estimation algorithm
Lambda_star <- Lambda %*% d

# Calculate scaling constant used for OU process parameters
c_star_fit <- calc_constants(rep(1, nrow(theta_ou)), theta_ou, sigma_ou)
cat('  convergence code for calculating constant matrix:', c_star_fit$code, '\n')
cat('  new constants are:', c_star_fit$c_vec, '\n')
c_star <- c_star_fit$c_vec

# Rescale theta_ou and sigma_ou s.t. the OU process has stationary var of 1
theta_ou_star <- update_theta(theta_ou,
                                diag(x = c_star, nrow = length(c_star)))
sigma_ou_star <- update_sigma(sigma_ou,
                                diag(x = c_star, nrow = length(c_star)))
Sigma_u_true <- Sigma_u # doesn't require rescaling
Sigma_e_true <- Sigma_e # doesn't require rescaling

# Save true parameters
true_res <- data.frame(variable = c(paste0('lambda_', 1:k),
                                    paste0('sigma2_u_', 1:k),
                                    paste0('sigma2_e_', 1:k),
                                    paste0(paste0('theta_ou_', 1:true_p),
                                           sort(rep(1:true_p, true_p))),
                                    paste0(paste0('sigma_ou_', 1:true_p), 1:true_p)),
                       value = c(c(Lambda_star[nonzero==1]),
                                 diag(Sigma_u_true),
                                 diag(Sigma_e_true),
                                 c(theta_ou_star),
                                 diag(sigma_ou_star)))

if (save_true == T){
  write.csv(true_res, paste0(dat_wd, 'ouf_true_', 'v', v, '.csv'))
}


################################################################################
#### STEP 3: Convert longitudinal data to long format
################################################################################

# Condense observed longitudinal outcomes over time
Y_wide <- all_person_Y 
colnames(Y_wide) <- rep(paste0('Y', 1:4), max_ni) # four longitudinal outcomes
rownames(Y_wide) <- 1:n_subj

Y1 <- Y_wide[,which(colnames(Y_wide) == 'Y1')]
Y2 <- Y_wide[,which(colnames(Y_wide) == 'Y2')]
Y3 <- Y_wide[,which(colnames(Y_wide) == 'Y3')]
Y4 <- Y_wide[,which(colnames(Y_wide) == 'Y4')]

Y_long <- data.frame(Y1 = matrix(as.matrix(t(Y1)), ncol = 1)[,1],
                     Y2 = matrix(as.matrix(t(Y2)), ncol = 1)[,1],
                     Y3 = matrix(as.matrix(t(Y3)), ncol = 1)[,1],
                     Y4 = matrix(as.matrix(t(Y4)), ncol = 1)[,1])
Y_long <- Y_long %>%
  filter(complete.cases(.)) %>%
  mutate(id = rep(1:n_subj, measurement_n))

# Add in measurement times
all_meas_times <- matrix(t(measurement_times), ncol = 1)
Y_long <- Y_long %>%
  mutate(meas_time = all_meas_times[!is.na(all_meas_times)])

# Rename a few columns for consistency
Y_long <- Y_long %>%
  mutate(user.id = id) %>%
  dplyr::select(c(user.id, meas_time, Y1, Y2, Y3, Y4))

# Save data
if (save_data == T){
  write.csv(Y_long, file = paste0(dat_wd, 'data/sim_dat_v', v,
                                  '_g', g, '.csv'), row.names = F)
}

