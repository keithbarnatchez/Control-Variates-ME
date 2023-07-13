library(AIPW)
library(SuperLearner)
library(reshape2)
library(tidyverse)
library(ggridges)
library(mice)
rm(list=ls())

source('ms_functions.R')

# covariances
rho12 <- c(0.25,-0.25)
rho13 <- c(0.5,-0.5)
rho23 <- c(0.6,-0.6)
sigs <- gen_sigmas(rho12,rho13,rho23)

# coefficients
alphas <- c(0.1,-0.5,0.3,0.85) # propensity score coeffs
betas <- c(0,1,-3,.5) # outcome model main effectd from x
tau <- c(1,1) # baseline tmt effect without interaction 
gammas <- c(0.2,0.4,-0.6) # interaction terms 

sig_y <- c(1,1) # variance in outcome across popuation
shift <- 1 # shift in means across populations
basemu <- c(1,1,1) # baseline mean values
n1 <- 300 ; n2 <- 3000
outcome_model <- 'normal'

sig_u <- 1
#######################################################
res_list <- list()
# set up range for n1
nsim <- 100 ; nboot <- 50

# set up parameters for this iteration
params <- list(n1=n1,n2=n2,
               sig1=sigs[[1]],sig2=sigs[[2]],sig_y=sig_y,
               basemu=basemu,shift=shift,
               alphas=alphas,betas=betas,tau=tau,gammas=gammas,
               sig_u = sig_u,
               outcome_model=outcome_model)


thedata <- gen_data(n1, n2, # sample sizes
                      sigs[[1]], sigs[[2]],sig_y, # covariance matrices and outcome variance
                      basemu,shift, # baseline mean and shift for other sample
                      alphas,betas,tau,gammas, sig_u,
                      outcome_model='normal')

thedatamat <- rbind(thedata[[1]],thedata[[2]])
thedatamat$S <- c(rep(1,n1),rep(0,n2))
thedatamat <- thedatamat %>% mutate(X3 = ifelse(S==1,X3,NA))

nsim <- 50
thetaus <- rep(NA, nsim)
for (s in 1:nsim) {
  thedata <- gen_data(n1, n2, # sample sizes
                      sigs[[1]], sigs[[2]],sig_y, # covariance matrices and outcome variance
                      basemu,shift, # baseline mean and shift for other sample
                      alphas,betas,tau,gammas, sig_u,
                      outcome_model='normal')
  
  thetau1s <- get_tau(thedata[[1]])
  
  # thedatamat <- rbind(thedata[[1]],thedata[[2]])
  # thedatamat$S <- c(rep(1,n1),rep(0,n2))
  # thedatamat <- thedatamat %>% mutate(X3 = ifelse(S==1,X3,NA))
  # thetaus[s] <- mime(thedatamat)
}
