# Helper functions for simulation study (DR extension)

# minor utility funcs
expit <- function(o) {
  return( exp(o)/(1+exp(o)))
}

demean <- function(vec) {
  return(vec-mean(vec))
}

#-------------------------------------------------------------------------------
# DATA FUNCTIONS

# Covariance matrices
gen_sigma <- function(rho12,rho13,rho23) {
  #' Given covariance terms in each population,
  #' generates covariance matrices (with 1s on the diagonals)

  sig <-  matrix(c(1,rho12,rho13,
                   rho12,1,rho23,
                   rho13,rho23,1),nrow=3,byrow=F)

  return(sig)
}

make_outcomes <- function(mu,outcome_model) {
  n <- length(mu) ;
  if (outcome_model == 'normal') {
    y <- mu + rnorm(n,sd=1)
  }
  if (outcome_model == 'binary') {
    y <- rbinom(n, size=1, prob=expit(mu))
  }
  if (outcome_model == 'exponential') {
    y <- rexp(n,rate=1/abs(mu))
  }
  return(y)
}

mis_class_model <- function(A, sens, spec) {
  Astar <- A ; n <- length(A)
  for (i in 1:n) {
    if (A[i]==1) {
      if (runif(1) < (1-sens[i])) {Astar[i] <- 1-A[i]}
    }
    if (A[i]==0) {
      if (runif(1) < (1-spec[i])) {Astar[i] <- 1-A[i]}
    }
  }
  return(Astar)
}

sampling_vector <- function(Xmat,rho,etas) {
  
  # First, get sampling probabilities. We set these so that on average a subject
  # has prob. of rho of being sampled
  if (any(is.na(etas))) {
    final_probs = rep(rho,nrow(Xmat))
  } else {
    init_probs <- expit(Xmat %*% etas)
   final_probs <- rho * (init_probs / mean(init_probs))
  }

  # Then, get the final sample membership vector
  samp_idx <- rbinom(length(final_probs),size=1,prob = final_probs)
  return(list(samp_idx=samp_idx,
              val_probs=final_probs)) 
}

gen_data <- function(n, rho, # sample sizes
                     sig,sig_y, # covariance matrices and outcome variance
                     basemu, # baseline mean and shift for other sample
                     alphas,betas,tau,gammas,
                     etas,
                     sens, spec,
                     outcome_model='normal',
                     error_systematic=FALSE) { # coefficients

  # generate covariates
  Xmat <- mvtnorm::rmvnorm(n,mean = basemu,sigma=sig)

  # Generate tmt model
  probs <- expit(cbind(1,Xmat)%*%alphas)
  A <- rbinom(n,1,prob=probs)

  # Outcome model
  mu <- tau[1]*A + cbind(1,Xmat)%*%betas + A*Xmat%*%gammas
  y <- make_outcomes( mu,outcome_model )

  # Misclassification model
  mults <- 0
  if (error_systematic==TRUE) {
    mults <- expit( cbind(1,Xmat)%*%betas )/10
  }
  sens <- rep(sens,n) + mults ; spec <- rep(spec,n) + mults
  Astar <- mis_class_model(A,sens,spec)

  # Sample validation data
  val_idx_res <- sampling_vector(Xmat,rho,etas)
  val_idx <- val_idx_res$samp_idx

  df <- data.frame(Y=y,Astar=Astar,A=A, X1=Xmat[,1], X2=Xmat[,2], X3=Xmat[,3],
                   val_idx=val_idx,val_probs=val_idx_res$val_probs)

  return(df)

}
#-------------------------------------------------------------------------------
# ESTIMATION FUNCTIONS

get_tau_hat <- function(df,
                        sl.lib = c("SL.mean","SL.glm","SL.glm.interaction")) {
  #' Estimates ATE in df via AIPW
  #'
  #' Returns the ATE estimate
  #'

  res <- AIPW$new(Y=df$Y,
                  A=df$A,
                  W=subset(df,select=c("X1","X2","X3")),
                  Q.SL.library = sl.lib,
                  g.SL.library = sl.lib,
                  k_split = 1,
                  verbose=FALSE)$fit()$summary()

  return(res)
}

get_psi <- function(df) {
  #' Constructs psi1 and psi2 from dfs 1 and 2 by performing an error-prone
  #' AIPW regression (only using X1 and X2)
  #'
  #' Recall psi1 and psi2 in our examples are the error-prone treatment effect
  #' estimators

  dfs <- list()
  dfs[[1]] <- df # full data
  dfs[[2]] <- df %>% filter(val_idx==1) # validation data

  sl.lib <- c("SL.mean","SL.glm","SL.glm.interaction")
  psi1 <- AIPW$new(Y=dfs[[1]]$Y,
                   A=dfs[[1]]$Astar,
                   W=subset(dfs[[1]],select=c("X1","X2","X3")),
                   Q.SL.library = sl.lib,
                   g.SL.library = sl.lib,
                   k_split = 1,
                   verbose=FALSE)$fit()$summary() #$estimates$RD['Estimate']

  psi2 <- AIPW$new(Y=dfs[[2]]$Y,
                   A=dfs[[2]]$Astar,
                   W=subset(dfs[[2]],select=c("X1","X2","X3")),
                   Q.SL.library = sl.lib,
                   g.SL.library = sl.lib,
                   k_split = 1,
                   verbose=FALSE)$fit()$summary() #$estimates$RD['Estimate']

  res_list = list(psi1=psi1,psi2=psi2)

  return(res_list)

}

get_ATE_cv <- function(df) {
  #' Given a dataframe df (containing a validation data indicator), implements
  #' the control variates method to obtain ATE estimate
  #'
  #' INPUTS:
  #' - a dataframe from the gen_data() function
  #' OUTPUTS:
  #' - a list containing the ATE estimate, variance estimate and 95% CI

  # Keep track of val data
  df_val <- df %>% filter(val_idx==1)

  # Step 1: estimate ATE in validation data
  tau_val_mod <- get_tau_hat(df_val)
  tau_hat_val <- tau_val_mod$estimates$RD['Estimate']
  v_hat <- var(tau_val_mod$obs_est$aipw_eif1 -
                 tau_val_mod$obs_est$aipw_eif0)/nrow(df_val)

  # Step 2: estimate control variates
  cv_mods <- get_psi(df)
  tau_ep_val <- cv_mods$psi1$estimates$RD['Estimate']
  tau_ep_main <- cv_mods$psi2$estimates$RD['Estimate']

  # Step 3: estimate Gamma and V
  Sig_hat <- get_vcov_asym(tau_val_mod,cv_mods,df$val_idx)
  gamma_hat <- Sig_hat$gamma_hat ; V_hat <- Sig_hat$V_hat

  # Step 4: subtract off (may need to flip the subtraction sign)
  tau_cv <- tau_hat_val - (gamma_hat/V_hat)*(tau_ep_main-tau_ep_val)

  # Get variance estimate and 95% CI
  var_hat <- v_hat - gamma_hat^2/V_hat
  ci_low <- tau_cv - qnorm(0.975)*sqrt(var_hat)
  ci_high <- tau_cv + qnorm(0.975)*sqrt(var_hat)


  # Return the ATE est and associated variance
  return(list(tau_cv=tau_cv,
              var_hat=var_hat,
              CI=c(ci_low,ci_high))
  )

}

get_vcov_boot <- function(df, nboot=25) {
  #' Code for estimating the Gamma and V terms from the asymptotic covariance
  #' matrix (see e.g. thm 1) via bootstrapping
  #'
  #' Note: this code isn't being used in simulations at the moment since we
  #' have faster alternatives

  tau_m <- get_tau(dfs[[1]])
  psi_m <- get_psi(df)

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

mime <- function(data,m=20,
                 congenial=T) {
  #' Performs multiple imputation for ME correction, following Webb-Vargas et
  #' al. (2015) with the one modification that we assume there is an internal,
  #' not external, validation sample. Makes use of the 'mice' package to perform
  #' multiple imputations

  if (congenial==F) { # if want to do uncongenial mime, i.e. without Y
    Yextract <- data$Y
    data <- data %>% select(-Y)
  }

  # Run the (congenial) MI procedure
  imps <- mice(data,method='pmm',m=m,printFlag = F)

  # Loop through imputed datasets, estimate PS
  vals <- sapply(1:m, function(d, imps, ...) {

    # use AIPW on d-th dataset
    curr_data <- complete(imps,d)
    if (congenial==F) {
      curr_data$Y <- Yextract # add Y back in
    }

    sl.lib <- c("SL.mean","SL.glm","SL.glm.interaction")
    res <- AIPW$new(Y=curr_data$Y,
                    A=curr_data$A,
                    W=subset(curr_data,select=c("X1","X2","X3")),
                    Q.SL.library = sl.lib,
                    g.SL.library = sl.lib,
                    k_split = 1,
                    verbose=FALSE)$fit()$summary()

    return(c(tau = res$estimates$RD['Estimate'],
             varval = res$estimates$RD['SE']^2))

  }, imps = imps)

  # grab ests and their SEs
  betas <- vals[1,]
  vars <- vals[2,]

  # compute SE est
  ovr_mean <- mean(betas)
  bw_var   <- sum( (betas - ovr_mean)^2 )/(length(betas)-1) # var bw ests
  wi_var   <- mean(vars) # within variance
  SE_est <- sqrt(wi_var + ((1 + (1/length(betas))) * bw_var ))
  CI <- c(-qnorm(.975)*SE_est + ovr_mean,
          qnorm(.975)*SE_est + ovr_mean)

  return(list(ATE=ovr_mean,
              var_est=SE_est^2,
              CI=CI))

}

main_sim <- function(nsim,
                     params,
                     nboot=500) {
  #' Main function for implementing simulation exercise
  #'
  #' INPUTS:
  #' - nsim: desired number of iterations
  #' - params: a list containing all relevant parameters needed for generating
  #'           the data (see the main simulation execution file for details)
  #' - nboot: when bootstrapping, the # of bootstrap iterations
  #'
  #' OUTPUTS:
  #' - a dataframe containing the results of each iteration

  # extract names from params
  list2env(params,envir = environment() )

  #register parallel backend
  registerDoParallel(cores=16)

  results <-foreach(s = 1:nsim, .combine = rbind) %dopar% {
    
    # sample data
    df <- gen_data(n, rho,
                   sig,sig_y,
                   basemu,
                   alphas,betas,tau,gammas,
                   etas,
                   sens, spec,
                   outcome_model,
                   error_systematic)
  

    df_val <- df %>% filter(val_idx==1)
    
    # Val data only
    tau1mod <- get_tau_hat(df_val)
    tau1hat <- tau1mod$estimates$RD['Estimate']
    var1hat  <- tau1mod$estimates$RD['SE']^2
    cilowhat <- tau1hat - 1.96*sqrt(var1hat)
    cihihat <- tau1hat + 1.96*sqrt(var1hat)

    # naive (use W as if true meas)
    naivemod <- AIPW$new(Y=df$Y,
                         A=df$Astar,
                         W=subset(df,select=c("X1","X2","X3")),
                         Q.SL.library = c("SL.mean","SL.glm","SL.glm.interaction"),
                         g.SL.library = c("SL.mean","SL.glm","SL.glm.interaction"),
                         k_split = 1,
                         verbose=FALSE)$fit()$summary()
    tau1naive <- naivemod$estimates$RD['Estimate']
    var1naive <- naivemod$estimates$RD['SE']^2
    cilownaive <- tau1naive - 1.96*sqrt(var1naive)
    cihinaive <- tau1naive + 1.96*sqrt(var1naive)

    # ideal (access to true meas)
    oraclemod <- AIPW$new(Y=df$Y,
                          A=df$A,
                          W=subset(df,select=c("X1","X2","X3")),
                          Q.SL.library = c("SL.mean","SL.glm","SL.glm.interaction"),
                          g.SL.library = c("SL.mean","SL.glm","SL.glm.interaction"),
                          k_split = 1,
                          verbose=FALSE)$fit()$summary()
    tau1oracle <- oraclemod$estimates$RD['Estimate']
    var1oracle <- oraclemod$estimates$RD['SE']^2
    ciloworacle <- tau1oracle - 1.96*sqrt(var1oracle)
    cihioracle <- tau1oracle + 1.96*sqrt(var1oracle)
    
    # Control variates: m mis=specified
    rhotemp=NA
    if (any(is.na(etas))) {
      rhotemp = rho
    }
    tau1cvgen_main_m <- cv_est_generalizability(df,rhotemp,
                                                m_lib='SL.glm')
    tau1cvgen_m <- tau1cvgen_main_m$tau_cv
    var1cvgen_m <- tau1cvgen_main_m$var_hat
    cilowcvgen_m <- tau1cvgen_m - 1.96*sqrt(var1cvgen_m)
    cihicvgen_m <- tau1cvgen_m + 1.96*sqrt(var1cvgen_m)
    tau1gen_m <- tau1cvgen_main_m$tau_hat_val
    var1gen_m <- tau1cvgen_main_m$v_hat_val
    
    # Control variates: g mis=specified
    tau1cvgen_main_g <- cv_est_generalizability(df,rhotemp,
                                                g_lib='SL.mean')
    tau1cvgen_g <- tau1cvgen_main_g$tau_cv
    var1cvgen_g <- tau1cvgen_main_g$var_hat
    cilowcvgen_g <- tau1cvgen_g - 1.96*sqrt(var1cvgen_g)
    cihicvgen_g <- tau1cvgen_g + 1.96*sqrt(var1cvgen_g)
    tau1gen_g <- tau1cvgen_main_g$tau_hat_val
    var1gen_g <- tau1cvgen_main_g$v_hat_val

    # do multiple imputation method (congenial and uncongenial)
    data <- df
    # data$A[!df$val_idx] <- NA # treat as missing

    # Congenial
    # tmp <-  try(mime(data,congenial = T,m=10))
    # if (inherits(tmp,'try-error')) {
    #   mime_res <- list(ATE=NA,var_est=NA)
    # } else {
    #   mime_res <- tmp
    # }
    tau1mime <- 1#mime_res$ATE
    var1mime <- 1#mime_res$var_est
    cilowmime <-1# tau1mime - 1.96*sqrt(var1mime)
    cihimime <- 1#tau1mime + 1.96*sqrt(var1mime)
    
    
    return(
      data.frame(tau1hat=tau1hat,
                      tau1mime=tau1mime,
                      tau1naive=tau1naive,
                      tau1oracle=tau1oracle,
                      tau1cv_m=tau1cvgen_m,
                      tau1cv_g=tau1cvgen_g,
                      var1cv_m=var1cvgen_m,
                      var1cv_g=var1cvgen_g,
                      var1hat=var1hat,var1oracle=var1oracle,
                      var1naive=var1naive,var1mime=var1mime,
                      cilowmime=cilowmime,cihimime=cihimime,
                      cilownaive=cilownaive,cihinaive=cihinaive,
                      ciloworacle=ciloworacle,cihioracle=cihioracle,
                      cilowgen_m=cilowcvgen_m,cihigen_m=cihicvgen_m,
                      cilowgen_g=cilowcvgen_g,cihigen_g=cihicvgen_g)
      )
  }

  # stop parallel backend
  stopImplicitCluster()

  return( results )
}
