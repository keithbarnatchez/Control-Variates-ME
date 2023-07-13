# testing code for multi study data
#
#
#
#
#

library(AIPW)
library(SuperLearner)
library(reshape2)
library(tidyverse)
rm(list=ls())

expit <- function(o) {
  return( exp(o)/(1+exp(o)))
}

# Covariance matrices
gen_sigmas <- function(rho12,rho13,rho23) {
  #' Given vectors of covariance terms in each population,
  #' generates covariance matrices (with 1s on the diagonals)

  sig1 <-  matrix(c(1,rho12[1],rho13[1],
                    rho12[1],1,rho23[1],
                    rho13[1],rho23[1],1),nrow=3,byrow=F)
  sig2 <-  matrix(c(1,rho12[2],rho13[2],
                    rho12[2],1,rho23[2],
                    rho13[2],rho23[2],1),nrow=3,byrow=F)
  
  return(list(sig1,sig2))
}

gen_data <- function(n1, n2, # sample sizes
                     sig1, sig2,sig_y, # covariance matrices and outcome variance
                     basemu,shift, # baseline mean and shift for other sample
                     alphas,betas,tau,gammas) { # coefficients
  
  # generate covariates
  Xmat1 <- mvtnorm::rmvnorm(n1,mean = basemu,sigma=sig1)
  Xmat2 <- mvtnorm::rmvnorm(n2,mean = (basemu+shift),sigma=sig2)
  
  # Generate tmt model
  probs1 <- expit(cbind(1,Xmat1)%*%alphas)
  probs2 <- expit(cbind(1,Xmat2)%*%alphas)
  A1 <- rbinom(n1,1,prob=probs1) ; A2 <- rbinom(n2,1,prob=probs2)
  
  # Outcome model
  y1 <- tau*A1 + cbind(1,Xmat1)%*%betas + A1*Xmat1%*%gammas + rnorm(n1,sd=1)
  y2 <- tau*A2 + cbind(1,Xmat2)%*%betas + A2*Xmat2%*%gammas + rnorm(n2,sd=1)
  
  return(
    list(
      df1=data.frame(Y=y1,A=A1,X1=Xmat1[,1],X2=Xmat1[,2],X3=Xmat1[,3]),
      df2=data.frame(Y=y2,A=A2,X1=Xmat2[,1],X2=Xmat2[,2],X3=Xmat2[,3])
    )
  )
  
}

get_tau <- function(df) {
  #' Estimates ATE in df via AIPW
  #'
  #' Returns the ATE estimate
  sl.lib <- c("SL.mean","SL.glm")
  res <- AIPW$new(Y=df$Y,
                      A=df$A,
                      W=subset(df,select=c("X1","X2","X3")),
                      Q.SL.library = sl.lib,
                      g.SL.library = sl.lib,
                      k_split = 1,
                  verbose=FALSE)$fit()$summary()
  
  return(res$estimates$RD['Estimate'])
}

get_psi <- function(dfs) {
  #' Constructs psi1 and psi2 from dfs 1 and 2 by performing an error-prone
  #' AIPW regression (only using X1 and X2)
  #'
  
  # Construct error prone estimates in both datasets
  sl.lib <- c("SL.mean","SL.glm")
  psi1 <- AIPW$new(Y=dfs[[1]]$Y,
                  A=dfs[[1]]$A,
                  W=subset(dfs[[1]],select=c("X1","X2")),
                  Q.SL.library = sl.lib,
                  g.SL.library = sl.lib,
                  k_split = 1,
                  verbose=FALSE)$fit()$summary()$estimates$RD['Estimate']
  
  psi2 <- AIPW$new(Y=dfs[[2]]$Y,
                  A=dfs[[2]]$A,
                  W=subset(dfs[[2]],select=c("X1","X2")),
                  Q.SL.library = sl.lib,
                  g.SL.library = sl.lib,
                  k_split = 1,
                  verbose=FALSE)$fit()$summary()$estimates$RD['Estimate']
  
  return(c(psi1, psi2))
  
}

get_vcov <- function(dfs, nboot=25) {
  
  tau_m <- get_tau(dfs[[1]])
  psi_m <- get_psi(dfs)
  
  psi <- matrix(NA,nrow=nboot,ncol=2)
  tau <- rep(NA,nboot)
  for (b in 1:nboot) {
 
    # first sample with replacement
    dfstemp <- dfs
    idx1 <- sample(1:nrow(dfstemp[[1]]),replace=TRUE)
    idx2 <- sample(1:nrow(dfstemp[[2]]),replace=TRUE)
    dfstemp[[1]] <- dfstemp[[1]][idx1,]
    dfstemp[[2]] <- dfstemp[[2]][idx2,]
    
    # get psi and tau
    tautemp <- get_tau(dfstemp[[1]])
    psitemp <- get_psi(dfstemp)
    
    psi[b,] <- psitemp ; tau[b] <- tautemp
  }
  print(psi)
  gamma_hat <- (1/(nboot-1))*sum( (tau-tau_m)*(psi[,1]-psi[,2]-psi_m[1]+psi_m[2]) )
  v_hat <- (1/(nboot-1))*sum( (psi[,1]-psi[,2]-psi_m[1]+psi_m[2])^2 )

  return(list(gamma_hat=gamma_hat,
              v_hat=v_hat))
}

control_variates <- function(dfs,nboot=100) {
  Sighat <- get_vcov(dfs,nboot)
  tau1hat <- get_tau(dfs[[1]])
  psi <- get_psi(dfs)
  
  ghat <- Sighat$gamma_hat
  vhatinv <- 1/Sighat$v_hat
  
  tau_tilde <- tau1hat - ghat*vhatinv*(psi[1]-psi[2])
  
  return(tau_tilde)
}

main_sim <- function(nsim,
                     params,
                     nboot=50) {
  
  # extract names from params
  list2env(params,envir = environment() )
  print(shift)
  tau1hat <- rep(NA,nsim)
  tau1tilde <- rep(NA,nsim)
  psi1 <- rep(NA,nsim) ; psi2 <- rep(NA,nsim)
  for (s in 1:nsim) {
    print(s)
    # sample data
    dfs <- gen_data(n1,n2,
                    sig1,sig2,sig_y,
                    basemu,shift,
                    alphas,betas,tau,gammas)
    
    # get estimates
    tau1hat[s] <- get_tau(dfs[[1]])
    tau1tilde[s] <- control_variates(dfs,nboot)
    psis <- get_psi(dfs)
    psi1[s] <- psis[1] ; psi2[s] <- psis[2]
  }
  return(data.frame(tau1hat=tau1hat,
                    tau1tilde=tau1tilde,
                    psi1=psi1,psi2=psi2))
}
#######################################################
# SIMULATION PARAMETERS
# covariances
rho12 <- c(0.25,0.25)
rho13 <- c(0.7,0.7)
rho23 <- c(-0.2,-0.2)
sigs <- gen_sigmas(rho12,rho13,rho23)

# coefficients
alphas <- c(0.1,-0.5,0.3,0.45) # propensity score coeffs
betas <- c(0,1,-3,.5) # outcome model main effectd from x
tau <- 1 # baseline tmt effect without interaction 
gammas <- c(0.2,0.4,-0.6) # interaction terms 

sig_y <- c(1,2) # variance in outcome across popuation
shift <- 0 # shift in means across populations
basemu <- c(1,1,1) # baseline mean values
n1 <- 150 ; n2 <- 1500
#######################################################
# SIMULATION 1: IDEAL SCENARIO
nsim <- 100
params <- list(n1=n1,n2=n2,
               sig1=sigs[[1]],sig2=sigs[[2]],sig_y=sig_y,
               basemu=basemu,shift=shift,
               alphas=alphas,betas=betas,gammas=gammas)
res_none <- main_sim(nsim,params,nboot=50)
#######################################################
# Scenario 2: small shift

params$shift <- c(0.025,-0.025,0.05)
res_low <- main_sim(nsim,params,nboot=50)
#######################################################
# Scenario 3: moderate shift
params$shift <- c(0.05,-0.05,0.1)
res_mod <- main_sim(nsim,params,nboot=50)
#######################################################
# Scenario 4: big shift
params$shift <- c(0.1,-0.1,0.2)
res_high <- main_sim(nsim,params,nboot=50)
#######################################################

# Plots

convert_res <- function(resdf) {
  nn <- nrow(resdf)
  newvec <- c(resdf[,1],resdf[,2])
  labs <- c(rep('Main',nn),rep('C.V.',nn))
  
  df <- data.frame(tau_est=newvec,method=labs)
  return(df)
}

# Make plotting dfs
res_none_long <- convert_res(res_none)
res_none_long$scenario <- 'No shift'

res_low_long <- convert_res(res_low)
res_low_long$scenario <- 'Low shift'

res_high_long <- convert_res(res_high)
res_high_long$scenario <- 'High shift'

res_full <- rbind(res_high_long,res_low_long,res_none_long)

res_full %>% ggplot(aes(x=factor(method),y=tau_est, fill=factor(method))) +
  geom_boxplot() + geom_hline(yintercept = 1,color='red') + theme_bw() + 
  facet_wrap(~factor(scenario)) + theme(legend.position="bottom")

# Do same but with psi

convert_psi <- function(resdf) {
  nn <- nrow(resdf)
  newvec <- c(resdf[,3],resdf[,4])
  labs <- c(rep('C.V. 1',nn),rep('C.V. 2',nn))
  
  df <- data.frame(psi_hat=newvec,CV=labs)
  return(df)
}

# Make plotting dfs
psi_none_long <- convert_psi(res_none)
psi_none_long$scenario <- 'No shift'

psi_low_long <- convert_psi(res_low)
psi_low_long$scenario <- 'Low shift'

psi_high_long <- convert_psi(res_high)
psi_high_long$scenario <- 'High shift'

psi_full <- rbind(psi_high_long,psi_low_long,psi_none_long)

psi_full %>% ggplot(aes(x=factor(CV),y=psi_hat, fill=factor(CV))) +
  geom_boxplot() + geom_hline(yintercept = 1,color='red') + theme_bw() + 
  facet_wrap(~factor(scenario)) + theme(legend.position="bottom")


