# helper functions for implementing scenario 2 estimators

is_binary <- function(x) {
  unique_values <- unique(x)
  if (length(unique_values) <= 2 && all(unique_values %in% c(0, 1))) {
    return(TRUE)
  } else if (length(unique_values) == 1 && unique_values %in% c(0, 1)) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}
mu_hat_est <-function(data,
                      sl.lib = c('SL.mean','SL.glm','SL.glm.interaction'),
                      cov_names = c('X1','X2','X3')) {
  #' Function for modeling outcome/conditional mean
  #'
  #' INPUTS:
  #' - data: a dataframe generated from ---
  #' - sl.lib: libraries for SuperLearner
  #' - cov_names: vector of covariate names to pass into SuperLearner
  #' - 
  #' 
  #' OUTPUTS:
  #' 
  
  # First, make sure to only fit on validation data
  data_fit <- data %>% filter(val_idx==1)
  
  Y = data_fit$Y
  X = data_fit %>% select(all_of(c(cov_names,'A')))
  Xmain <- data %>% select(all_of(c(cov_names,'A')))
  
  # Check is y is binary or not
  fam <- 'gaussian'
  if (is_binary(Y)) {fam <- 'binomial'}
  
  # Fit the outcome model
  mu_mod <- SuperLearner::SuperLearner(Y=Y,
                                       X=X,
                                       SL.library=sl.lib,
                                       family=fam)
  
  # Get predicted values under treatment
  XA1 <- Xmain %>% mutate(A=1)
  mu_hat1 <- predict(mu_mod,
                        XA1,
                        onlySL=TRUE)
  
  # Get predicted values under no treatment
  XA0 <- Xmain %>% mutate(A=0)
  mu_hat0 <- predict(mu_mod,
                     XA0,
                     onlySL=TRUE)
  
  return(list(mu_hat1=mu_hat1, mu_hat0 = mu_hat0))
}

pi_hat_est <- function(data,
                       sl.lib = c('SL.mean','SL.glm','SL.glm.interaction'),
                       cov_names = c('X1','X2','X3')) {
  #' Function for modeling propensity scores
  #'
  #'
  
  # Only fit with val data
  data_fit <- data %>% filter(val_idx==1)
  A <- data_fit$A ; X <- data_fit %>% select(all_of(cov_names))
  Xmain <- data %>% select(all_of(cov_names))
  
  # Fit the ps model on the validation data
  pi_mod <- SuperLearner::SuperLearner(Y=A,
                                       X=X,
                                       SL.library =sl.lib,
                                       family='binomial')
  
   # Get predicted values over the whole dataset
  pi_hat <- predict(pi_mod,Xmain,onlySL=TRUE)
  
  return( pi_hat )
  
}

kappa_hat_est <- function(data,
                          sl.lib = c('SL.mean','SL.glm','SL.glm.interaction'),
                          cov_names = c('X1','X2','X3')) {
  #' Function for modeling participation probabilities P(S=1|X)
  #'
  #'

  # Fit the ps model on the validation data
  S <- data$val_idx ; X <- data %>% select(all_of(cov_names))
  kappa_mod <- SuperLearner::SuperLearner(Y=S,
                                       X=X,
                                       SL.library =sl.lib,
                                       family='binomial')
  
  # Get predicted values over the whole dataset
  kappa_hat <- predict(kappa_mod,onlySL=TRUE)
  
  return(kappa_hat)
  
}



generalizability_eif <- function(mu_hat1, mu_hat0,
                                 pi_hat,
                                 kappa_hat,
                                 data) {
  #' Function for obtaining the generalizability efficient influence function
  #' terms
  #' 
  #' INPUTS:
  #' - mu_hat, pi_hat and kappa_hat contain vectors of the fitted
  #'   values for all 3 nuisance models
  #' - S: a binary vector denoting whether each observation is a member
  #'   of the validation data
  #' - n_main and n_val: integers denoting the sample size of the main and val
  #'   data (these sets are nested, i.e. size of whole sample is n_main)
  #'
  #' OUTPUTS: 
  #' - 
  #' 
  
  S <- data$val_idx ; A <- data$A ; Y <- data$Y
  nmain <- length(S) ; nval <- sum(S)
  
  psihat1 <-  (S*A*(Y-mu_hat1)/(kappa_hat*pi_hat) +
                           mu_hat1)
  
  psihat0 <- (S*(1-A)*(Y-mu_hat0)/(kappa_hat*(1-pi_hat)) +
                            mu_hat0)
  
  return(psihat1 - psihat0)
  
}

generalizability_est <- function(data, 
                                 rho=NA) {
  
  # Get nuisance model estimates
  # Outcome model
  mu_hat_ests <- mu_hat_est(data)
  mu_hat1 <- mu_hat_ests$mu_hat1$pred
  mu_hat0 <- mu_hat_ests$mu_hat0$pred
  
  # Propensity score model
  pi_hat <- pi_hat_est(data)$pred
  
  # Participation model
  # Trying setting where we always estimated participation probability, 
  # regardless of study design
  # kappa_hat <- rho
  # if (is.na(rho)) {
    kappa_hat <- kappa_hat_est(data)$pred
  # }
  
  # Form EIF
  genz_eif <- generalizability_eif(mu_hat1, mu_hat0,
                                   pi_hat,
                                   kappa_hat,
                                   data)
  
  # Use EIF to get ATE est and variance est
  # 7/6/2023: changed denom
  ate_est <- mean(genz_eif)
  var_est <- var(genz_eif)/length(genz_eif) # sqrt(length(genz_eif))
  
  return(list(ATE=ate_est,
              EIF=genz_eif))
  
}

cv_est_generalizability <- function(data,
                                    rho=NA) {
  #'
  #'
  #'
  
  n <- nrow(data)
  
  # First, get ATE estimate pre-variance reduction
  ate_eif_genz <- generalizability_est(data,rho)
  tau_hat_val <- ate_eif_genz$ATE
  
  # Next, get error-prone ATE estimates (val data -> main)
  data_ep <- data %>% mutate(A=Astar)
  ate_eif_ep_val_genz <- generalizability_est(data_ep,rho)
  tau_hat_ep_val <- ate_eif_ep_val_genz$ATE
  
  # Next, get error-prone ATE estimates over the main data
  data_main <- data %>% mutate(A=Astar) # use ep A vals 
  ate_eif_ep_main <- get_tau_hat(data_main)
  tau_hat_ep_main <- ate_eif_ep_main$estimates$RD['Estimate']
  
  # Now, construct CV estimator. First, get Gamma
  eif_val <- ate_eif_genz$EIF
  eif_val_ep <- ate_eif_ep_val_genz$EIF
  eif_main_ep <- ate_eif_ep_main$obs_est$aipw_eif1 -
                 ate_eif_ep_main$obs_est$aipw_eif0
  gamma_hat <- 1/n * (cov(eif_val,eif_main_ep) -
                      cov(eif_val,eif_val_ep))

  # print(gamma_hat)
  
  # Next, get V
  V_hat <- 1/n * var(eif_val_ep + eif_main_ep)
  # print(gamma_hat/V_hat)
  # Finally, form the CV estimator
  tau_cv <- tau_hat_val - gamma_hat/V_hat * (tau_hat_ep_main - tau_hat_ep_val)
  
  return(list(tau_cv=tau_cv,
              gamma_hat=gamma_hat,
              V_hat=V_hat))
}
