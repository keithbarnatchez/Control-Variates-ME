# genz_est_test.R
# testing code for implementation of generalizability estimator

library(AIPW)
library(SuperLearner)
library(reshape2)
library(tidyverse)
library(ggridges)
library(mice)
rm(list=ls())

source('ms-exposure-funcs.R')
source('generalizability_funcs.R')
# ------------------------------------------------------------------------------
# Set up relevant paths for output
flnm <- gsub(':','-',Sys.time()) # filename suffixed with date/time
flnm <- paste('sim_results_',gsub(' ','_',flnm),'.csv',sep = '')

simdir <- 'output/sim_results/genz' # directory
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

# OUTCOME DISTRIBUTION
sig_y <- 5 # variance in outcome across popuation
outcome_model <- 'normal'

# SAMPLE SIZE
n <- 3000 ; rho <- 0.1

# SENSITIVITY AND SPECIFICITY OF M.E. PROCESS
sens <- .9 ; spec <- .9
error_systematic <- FALSE

sens_grid <- c(0.8,0.85,0.9,0.95) 
n_grid <- c(1000, 2000, 3000, 5000)
rho_grid <- c(0.05, 0.1, 0.15, 0.2)
etas <- rep(NA,3) #c(0.1,-0.2,0.6)


nsim <- 50
tmt_ests <- matrix(NA, nrow=nsim, ncol=8)
for (s in 1:nsim) {
  print(s)
  # sample data
  data <- gen_data(n, rho,
                 sig,sig_y,
                 basemu,
                 alphas,betas,tau,gammas,
                 etas,
                 sens, spec,
                 outcome_model,
                 error_systematic)
  
  # Non-cv generalizability est
  # tmt_ests[s,1] <- generalizability_est(data,rho)$ATE
  
  # # Error prone via val data
  # data_val <- data %>% mutate(A=Astar)
  # tmt_ests[s,2] <- generalizability_est(data_val,rho)$ATE
  # 
  # # Error prone with all data
  # data_main <- data %>% mutate(A=Astar)
  # tmt_ests[s,3] <- get_tau_hat(data_main)$estimates$RD['Estimate']
  
  # CV generalizability
  temp <- cv_est_generalizability(data)
  
  tmt_ests[s,4] <- temp$tau_cv
  
  tmt_ests[s,5] <- temp$gamma_hat
  
  tmt_ests[s,6] <- temp$V_hat
  
  # CV under setting 1
  tmt_ests[s,7] <- get_ATE_cv(data)$tau_cv
  
  # of the MIME
  data_mime <- data
  data_mime$A[!data$val_idx] <- NA # treat as missing
  tmt_ests[s,8] <- mime(data_mime)$ATE
}

# simmys <- 1e5
# for (s in 1:simmys) {
#   print(s)
#   data <- gen_data(n, rho,
#                    sig,sig_y,
#                    basemu,
#                    alphas,betas,tau,gammas,
#                    sens, spec,
#                    outcome_model,
#                    error_systematic)
#   if (any(is.na(data$A))) {
#     print(s)
#     stop('check now')
#   }
# }

# write.csv(tmt_ests,file=fullpath)


