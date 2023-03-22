rm(list = ls())

## Preliminaries

library(parallel)
library(abind)
library(splines)
library(SuperLearner)

# Code for generating and fitting data
source("~/Github/ME-Causal-Review/R/error_outcome/calibrate.R")
source("~/Github/ME-Causal-Review/R/error_outcome/me_erf.R")

### Test RC/DR estimator
n.sim <- 100

# model arguments
a.vals <- seq(6, 14, by = 0.02)
n <- 4000
sigma <- 2
omega <- 1
tau <- 0.5
start <- Sys.time()

out <- mclapply(1:n.sim, function(i, ...) {
  
  # covariates
  x1 <- stats::rnorm(n, 0, 1)
  x2 <- stats::rnorm(n, 0, 1)
  x3 <- stats::rnorm(n, 0, 1)
  x4 <- stats::rnorm(n, 0, 1)
  x <- cbind(x1, x2, x3, x4)
  
  mu_gps <- 10 + 0.5*x[,1] - 0.5*x[,2] - 0.5*x[,3] + 0.5*x[,4]
  
  a <- rnorm(n, mu_gps, sigma)
  
  mu_out <- 2 - 0.75*x[,1] - 0.25*x[,2] + 0.25*x[,3] + 0.75*x[,4] +
    0.25*(a - 10) - 0.75*cos(pi*(a - 6)/4) - 0.25*(a - 10)*x[,1]
  
  eta_out <- 0.5 + 0.25*x[,1] - 0.75*x[,2] + 0.75*x[,3] - 0.25*x[,4] + 0.25*(a - 6)
  
  y_true <- rnorm(n, mu_out, omega)
  y_me <- y_true + rnorm(n, eta_out, tau)
  
  rho <- plogis(-2 + 0.5*x[,1] - 0.5*x[,4])
  val <- rbinom(n, 1, rho)
  
  true_erc <- sapply(a.vals, function(a.tmp, ...) {
    
    mean(2 - 0.75*x[val == 0,1] - 0.25*x[val == 0,2] + 0.25*x[val == 0,3] + 0.75*x[val == 0,4] +
           0.25*(a.tmp - 10) - 0.75*cos(pi*(a.tmp - 6)/4) - 0.25*(a.tmp - 10)*x[val == 0,1])
    
  })
  
  A0 <- a[val == 0]
  A1 <- a[val == 1]
  X0 <- cbind(1, x[val == 0,])
  X1 <- cbind(1, x[val == 1,])
  Y0 <- y_me[val == 0]
  Y1_me <- y_me[val == 1]
  Y1_true <- y_true[val == 1]
  
  adjust <- out_me(A0 = A0, A1 = A1, X0 = X0, X1 = X1,
                   Y0 = Y0, Y1_me = Y1_me, Y1_true = Y1_true, 
                   a.vals = a.vals, bw = 0.5)
  naive <- out_naive(A0, X0, Y0, a.vals = a.vals, bw = 0.5)
  
  return(list(est = t(data.frame(true_erc = true_erc, naive_est = naive$estimate, adjust_est = adjust$estimate)),
              se = t(data.frame(adjust_se = adjust$se, naive_se = naive$se))))
  
}, mc.cores = 8, mc.preschedule = TRUE)

stop <- Sys.time()
stop - start

est <- abind(lapply(out, function(lst, ...) lst$est), along = 3)
se <- abind(lapply(out, function(lst, ...) lst$se), along = 3)
mu.mat <- matrix(rep(rowMeans(est[1,,]), n.sim), nrow = length(a.vals), ncol = n.sim)

# Estimate
out_est <- t(apply(est, 1, rowMeans, na.rm = T))
colnames(out_est) <- a.vals
rownames(out_est) <- c("True ERF","Naive","Bias Corrected")

# Mean Absolute Error
out_bias <- t(apply(est[2:3,,], 1, function(x) rowMeans(abs(x - mu.mat), na.rm = T)))
colnames(out_bias) <- a.vals
rownames(out_bias) <- c("Naive","Corrected")

# Mean Squared Error
out_mse <- t(apply(est[2:3,,], 1, function(x) rowMeans((x - mu.mat)^2, na.rm = T)))
colnames(out_mse) <- a.vals
rownames(out_mse) <- c("Naive","Corrected")

# Coverage Probability
cp <- list(as.matrix((est[2,,] - 1.96*se[1,,]) < mu.mat & (est[2,,] + 1.96*se[1,,]) > mu.mat),
           as.matrix((est[3,,] - 1.96*se[2,,]) < mu.mat & (est[3,,] + 1.96*se[2,,]) > mu.mat))
out_cp <- do.call(rbind, lapply(cp, rowMeans, na.rm = T))
colnames(out_cp) <- a.vals
rownames(out_cp) <- c("Naive","Corrected")

# Plots
plot(a.vals, out_est[1,], type = "l", col = "black", lwd = 2, ylim  = c(0,5),
     main = "Exposure Response Curve", xlab = "Exposure", ylab = "Rate of Event")
lines(a.vals, out_est[2,], type = "l", col = "blue", lwd = 2, lty = 1)
lines(a.vals, out_est[3,], type = "l", col = "red", lwd = 2, lty = 1)
grid()
legend("topleft", legend=c("True ERF","Naive ERF","Corrected ERF"),
       col= c("black", "blue", "red"), lty = c(1,1,1), lwd=2, cex=0.8)
