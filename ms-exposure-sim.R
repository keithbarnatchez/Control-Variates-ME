# SIMULATION PARAMETERS
# covariances
# ------------------------------------------------------------------------------

.libPaths("~/R_packages")
source('ms-exposure-funcs.R')
source('generalizability_funcs.R')
library(foreach)
library(doParallel)
library(AIPW)
library(mvtnorm)
library(SuperLearner)
library(reshape2)
library(tidyverse)
# library(ggridges)
library(mice)
# ------------------------------------------------------------------------------
# Set up relevant paths for output
flnm <- gsub(':','-',Sys.time()) # filename suffixed with date/time
flnm <- paste('sim_results_',gsub(' ','_',flnm),'.csv',sep = '')

simdir <- 'output/sim_results/' # directory
fullpath <- paste(simdir,flnm,sep='') # construct path to final file name 
# ------------------------------------------------------------------------------
# SIMULATION PARAMETERS

# COVARIATE DISTRIBUTION
basemu <- c(1,1,1) # baseline mean values
rho12 <- 0.25
rho13 <- 0.5
rho23 <- -0.4
sig <- gen_sigma(rho12,rho13,rho23)

# COEFFICIENTS FOR ALL MODELS
alphas <- c(0.1,-0.5,0.3,0.85) # propensity score coeffs
betas <- c(0,1,-3,.5) # outcome model main effectd from x
delta <- c(.1,.1,.1) # outcome ME coeffs
tau <- 1 # baseline tmt effect without interaction 
gammas <- c(0.2,0.4,-0.6) # interaction terms 
etas <- c(0.1,-0.2,0.6) # selection prob coefficients (can be turned off)

# OUTCOME DISTRIBUTION
sig_y <- 1 # variance in outcome across population
outcome_model <- 'normal'

# SAMPLE SIZE
n <- 5000 ; rho <- 0.1

# SENSITIVITY AND SPECIFICITY OF M.E. PROCESS
sens <- 0.70 ; spec <- 0.95
error_systematic <- FALSE

sens_grid <-   c(0.8, 0.85,0.9,0.95)  # sensitivities
n_grid <-  c(5000)  # sample size of main data
rho_grid <- c(0.1, 0.2, 0.3, 0.4, 0.5) # val. data proportion
err_grid <- c(TRUE,FALSE) # measurement errors systematic vs purely random
eta_grid <- c(TRUE,FALSE) # when true, selection into val. data depends on X

# Create the grid
param_combos <- expand.grid(sens_grid=sens_grid,n_grid=n_grid,rho_grid=rho_grid,
                            err_grid=err_grid,eta_grid=eta_grid)
# ------------------------------------------------------------------------------
# RUN SIMULATION

res_list <- list()
# set up range for n1
nsim <- 250 ; nboot <- 100
for (i in 1:nrow(param_combos)) {
  
  curr_etas <- rep(NA,3)
  if (param_combos$eta_grid[i] == TRUE) {
    curr_etas <- etas
  }
  
  # set up parameters for this iteration
  params <- list(n=param_combos$n_grid[i], rho=param_combos$rho_grid[i],
                 sig=sig,sig_y=sig_y,
                 basemu=basemu,
                 alphas=alphas,betas=betas,tau=tau,gammas=gammas,
                 etas=curr_etas,
                 sens=param_combos$sens_grid[i],spec=spec,
                 outcome_model=outcome_model,
                 error_systematic=param_combos$err_grid[i])
  
  # run simulation
  res <- main_sim(nsim, params, nboot)
  
  # add parameters to result data frame
  params_df_base <- as.data.frame(t(unlist(params)))
  params_df <- params_df_base
  for (r in 2:nrow(res)) {
    params_df <- rbind(params_df,params_df_base)
  }
  res <- cbind(params_df, res) # merge the param values for this grid point on
  rownames(res) <- NULL
  
  # append results to list
  res_list[[i]] <- res
  
}

# combine all results into a single data frame
final_res <- do.call(rbind, res_list)

final_res <- final_res %>% mutate(cicovmime = 1-as.numeric( cilowmime>1 | cihimime<1  ),
                                  cicovtilde = 1-as.numeric( cilowtilde>1 | cihitilde<1  ),
                                  cicovoracle = 1-as.numeric( ciloworacle>1 | cihioracle<1  ),
                                  cicovcvgen = 1-as.numeric( cilowcvgen>1 | cihicvgen<1  ))

# write the output
write.csv(final_res,file=fullpath)

