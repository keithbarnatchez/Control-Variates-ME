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

make_outcomes <- function(mu1,mu2,outcome_model) {
  n1 <- length(mu1) ; n2 <- length(mu2)
  if (outcome_model == 'normal') {
    y1 <- mu1 + rnorm(n1,sd=1)
    y2 <- mu2 + rnorm(n2,sd=1)
  }
  if (outcome_model == 'binary') {
    y1 <- rbinom(n1, size=1, prob=expit(mu1))
    y2 <- rbinom(n2, size=1, prob=expit(mu2))
  }
  if (outcome_model == 'exponential') {
    y1 <- rexp(n1,rate=1/abs(mu1))
    y2 <- rexp(n2,rate=1/abs(mu2))
  }
  return(list(y1=y1,y2=y2))
}

gen_data <- function(n1, n2, # sample sizes
                     sig1, sig2,sig_y, # covariance matrices and outcome variance
                     basemu,shift, # baseline mean and shift for other sample
                     alphas,betas,tau,gammas, sig_u,
                     outcome_model='normal') { # coefficients
  
  # generate covariates
  Xmat1 <- mvtnorm::rmvnorm(n1,mean = basemu,sigma=sig1)
  Xmat2 <- mvtnorm::rmvnorm(n2,mean = (basemu+shift),sigma=sig2)
  
  # Generate tmt model
  probs1 <- expit(cbind(1,Xmat1)%*%alphas)
  probs2 <- expit(cbind(1,Xmat2)%*%alphas)
  A1 <- rbinom(n1,1,prob=probs1) ; A2 <- rbinom(n2,1,prob=probs2)
  
  # Outcome model
  mu1 <- tau[1]*A1 + cbind(1,Xmat1)%*%betas + A1*Xmat1%*%gammas
  mu2 <- tau[1]*A2 + cbind(1,Xmat2)%*%betas + A2*Xmat2%*%gammas
  ylist <- make_outcomes(mu1,mu2,outcome_model)
  y1 <- ylist[[1]] ; y2 <- ylist[[2]]
  
  # Meas. error model
  W1 <- Xmat1[,3] + rnorm(n1,mean = 0,sd = sig_u)
  W2 <- Xmat2[,3] + rnorm(n2,mean = 0,sd = sig_u)
  
  return(
    list(
      df1=data.frame(Y=y1,A=A1,X1=Xmat1[,1],X2=Xmat1[,2],X3=Xmat1[,3],W=W1),
      df2=data.frame(Y=y2,A=A2,X1=Xmat2[,1],X2=Xmat2[,2],X3=Xmat2[,3],W=W2)
    )
  )
  
}

get_tau <- function(df) {
  #' Estimates ATE in df via AIPW
  #'
  #' Returns the ATE estimate
  sl.lib <- c("SL.mean","SL.glm","SL.glm.interaction")
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
  
  sl.lib <- c("SL.mean","SL.glm","SL.glm.interaction")
  psi1 <- AIPW$new(Y=dfs[[1]]$Y,
                   A=dfs[[1]]$A,
                   W=subset(dfs[[1]],select=c("X1","X2","W")),
                   Q.SL.library = sl.lib,
                   g.SL.library = sl.lib,
                   k_split = 1,
                   verbose=FALSE)$fit()$summary()$estimates$RD['Estimate']
  
  psi2 <- AIPW$new(Y=dfs[[2]]$Y,
                   A=dfs[[2]]$A,
                   W=subset(dfs[[2]],select=c("X1","X2","W")),
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
  vhatinv <- solve(Sighat$v_hat)
  
  tau_tilde <- tau1hat - ghat%*%vhatinv*(psi[1]-psi[2])
  
  return(tau_tilde)
}

mime <- function(data,m=20) {
  #' Performs multiple imputation for ME correction, following Webb-Vargas et 
  #' al. (2015) with the one modification that we assume there is an internal,
  #' not external, validation sample. Makes use of the 'mice' package to perform
  #' multiple imputations
  
  # Run the (congenial) MI procedure
  imps <- mice(data,method='pmm',m=m)
  
  # Loop through imputed datasets, estimate PS 
  vals <- sapply(1:m, function(d, imps, ...) {
    
    # use AIPW on d-th dataset
    curr_data <- complete(imps,d)
    
    sl.lib <- c("SL.mean","SL.glm","SL.glm.interaction")
    res <- AIPW$new(Y=curr_data$Y,
                    A=curr_data$A,
                    W=subset(curr_data,select=c("X1","X2","X3")),
                    Q.SL.library = sl.lib,
                    g.SL.library = sl.lib,
                    k_split = 1,
                    verbose=FALSE)$fit()$summary()
    
    return(c(tau = res$estimates$RD['Estimate'], se = 1))
    
  }, imps = imps)
  
  # grab ests and their SEs
  betas <- vals[1,]
  ses <- vals[2,]
  
  # compute SE est 
  ovr_mean <- mean(betas) 
  bw_var   <- sum( (betas - ovr_mean)^2 )/(length(betas)-1) # var bw ests
  wi_var   <- mean(ses) # within variance
  SE_est <- wi_var + ((1 + (1/length(betas))) * bw_var )
  CI <- c(-qnorm(.975)*SE_est + ovr_mean,
          qnorm(.975)*SE_est + ovr_mean)
  
  return(list(ATE=ovr_mean,
              CI=CI))
  
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
  tau1mime <- rep(NA, nsim) ; tau1rcal <- rep(NA,nsim)
  tau1naive <- rep(NA, nsim) ; tau1oracle <- rep(NA, nsim)
  for (s in 1:nsim) {
    print(s)
    # sample data
    dfs <- gen_data(n1,n2,
                    sig1,sig2,sig_y,
                    basemu,shift,
                    alphas,betas,tau,gammas,sig_u,
                    outcome_model)
    
    # get estimates
    tau1hat[s] <- get_tau(dfs[[1]])
    tau1tilde[s] <- control_variates(dfs,nboot)
    psis <- get_psi(dfs)
    psi1[s] <- psis[1] ; psi2[s] <- psis[2]
    
    # regression calibration
    tau1rcal[s] <- rcal(dfs)
    
    # naive (use W as if true meas)
    tau1naive[s] <- AIPW$new(Y=dfs[[1]]$Y,
                             A=dfs[[1]]$A,
                             W=subset(dfs[[1]],select=c("X1","X2","W")),
                             Q.SL.library = c("SL.mean","SL.glm","SL.glm.interaction"),
                             g.SL.library = c("SL.mean","SL.glm","SL.glm.interaction"),
                             k_split = 1,
                             verbose=FALSE)$fit()$summary()$estimates$RD['Estimate']
    
    # ideal (access to true meas)
    tau1oracle[s] <- AIPW$new(Y=dfs[[2]]$Y,
                              A=dfs[[2]]$A,
                              W=subset(dfs[[2]],select=c("X1","X2","X3")),
                              Q.SL.library = c("SL.mean","SL.glm","SL.glm.interaction"),
                              g.SL.library = c("SL.mean","SL.glm","SL.glm.interaction"),
                              k_split = 1,
                              verbose=FALSE)$fit()$summary()$estimates$RD['Estimate']
    
    # do multiple imputation method 
    dfs[[2]]$X3 <- NA # treat as missing
    data <- rbind(dfs[[1]],dfs[[2]])
    mime_res <- mime(data)
    tau1mime[s] <- mime_res$ATE
    
    
  }
  return(data.frame(tau1hat=tau1hat,
                    tau1tilde=tau1tilde,
                    psi1=psi1,psi2=psi2,
                    tau1mime=tau1mime,
                    tau1rcal=tau1rcal,
                    tau1naive=tau1naive,
                    tau1oracle=tau1oracle))
}

rcal <- function(dfs,
                 sl.lib= c("SL.mean","SL.glm","SL.glm.interaction")) {
  #  fit calibration model
  dfcal <- dfs[[1]]
  dfcal <- dfcal %>% select(-Y)
  calmod <- lm(X3 ~ ., data=dfcal)
  
  # impute with calibration data
  maindataforpred <- dfs[[2]] %>% select(-c(Y,X3))
  Xpred <- predict(calmod,newdata=maindataforpred)
  dfmain <- dfs[[2]] %>% mutate(X3=Xpred)
  
  # fit final model with imputed X
  tau1hat <- AIPW$new(Y=dfmain$Y,
                      A=dfmain$A,
                      W=subset(dfmain,select=c("X1","X2","X3")),
                      Q.SL.library = sl.lib,
                      g.SL.library = sl.lib,
                      k_split = 1,
                      verbose=FALSE)$fit()$summary()$estimates$RD['Estimate']
  
  return(tau1hat)
  
}