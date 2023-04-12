out_me <- function(A0, A1, X0, X1, Y0, Y1_me, Y1_true, 
                   a.vals = seq(min(A0), max(A0), length = 100), bw = 1, 
                   family = gaussian(), se.fit = TRUE, ...) {
  
  n1 <- nrow(X1)
  n0 <- nrow(X0)
  
  S <- rep(c(0,1), times = c(n0, n1))
  X <- rbind(X0, X1)
  A <- c(A0, A1)
  
  n <- n1 + n0
  m <- ncol(X)
  
  ## weight model
  
  # standardize a
  astar_0 <- c(A0 - mean(A0))/var(A0)
  astar2_0 <- c((A0 - mean(A0))^2/var(A0) - 1)
  astar_1 <- c(A1 - mean(A1))/var(A1)
  astar2_1 <- c((A1 - mean(A1))^2/var(A1) - 1)
  
  astar <- c(astar_0, astar_1)
  astar2 <- c(astar2_0, astar2_1)
  
  cmat <- cbind((1 - S)*X*astar, (1 - S)*astar2,
                S*X*astar, S*astar2, (1-S)*X, S*X)
  target <- c(rep(0, (2*m + 2)), n0*colMeans(X0), n1*colMeans(X0))
  
  mod <- calibrate(cmat = cmat, target = target)
  ipw <- mod$weights
  
  Y <- c(Y0, Y1_me - Y1_true)
  psi <- Y
  
  ## kernel weighted least squares
  out <- sapply(a.vals, kern_ipw, a = A, s = S, x = X, psi = psi, ipw = ipw, 
                cmat = cmat, astar = astar, astar2 = astar2, bw = bw, se.fit = se.fit)
  
  ## spline method
  # out0 <- spl_ipw(a = A0, psi = psi0, a.vals = a.vals, df = 6, se.fit = TRUE)
  # out1 <- spl_ipw(a = A1, psi = psi1, a.vals = a.vals, df = 6, se.fit = TRUE)
  
  estimate <- out[1,]
  se <- sqrt(out[2,])
  
  return(list(estimate = estimate, se = se, ipw0 = ipw[S == 0], ipw1 = ipw[S == 1]))
  
}

out_naive <- function(A0, X0, Y0, a.vals = seq(min(A0), max(A0), length = 100), bw = 1,...) {
  
  n0 <- nrow(X0)
  m <- ncol(X0)
  
  ## weight model
  
  # standardize a
  astar <- c(A0 - mean(A0))/var(A0)
  astar2 <- c((A0 - mean(A0))^2/var(A0) - 1)
  
  cmat <- cbind(X0*astar, astar2)
  target <- rep(0, m + 1)
  
  mod <- calibrate(cmat = cmat, target = target)
  ipw <- mod$weights
  
  ## kernel weighted least squares
  psi0 <- c(ipw*Y0)
  out <- sapply(a.vals, kern_naive, a = A0, psi = psi0, bw = bw, se.fit = TRUE)
  
  estimate <- out[1,]
  se <- sqrt(out[2,])
  
  return(list(estimate = estimate, se = se, ipw = ipw))
  
}

esteq <- function(p, s, x, psi, g.std, astar, astar2, ipw, eta, theta) {
  
  eq1 <- (1 - s)*p*x*astar
  eq2 <- (1 - s)*p*astar2
  eq3 <- s*p*x*astar
  eq4 <- s*p*astar2
  
  eq5 <- (1 - s)*(p*x - theta)
  eq6 <- s*(p*x - theta)
  eq7 <- (1 - s)*(x - theta)
  
  eq8 <- 2*p*(psi - eta)*g.std
  
  eq <- c(eq1, eq2, eq3, eq4, eq5, eq6, eq7, eq8) 
  return(eq)
  
}

# Kernel weighted least squares
kern_ipw <- function(a.new, a, s, x, psi, astar, astar2, cmat, ipw, bw = 1, se.fit = FALSE) {
  
  n <- length(a)
  m <- ncol(x)
  
  # Gaussian Kernel
  a.std <- (a - a.new) / bw
  k.std <- dnorm(a.std) / bw
  g.std <- cbind((1 - s)*a.std, s*a.std, (1 - s), s)
  
  b <- lm(psi ~ -1 + g.std, weights = k.std*ipw)$coefficients
  mu <- unname(b[3] - b[4])
  
  if (se.fit) {
    
    U <- matrix(0, ncol = 5*m + 2, nrow = 5*m + 2)
    V <- matrix(0, ncol = 5*m + 6, nrow = 4)
    meat <- matrix(0, ncol = 5*m + 6, nrow = 5*m + 6)
    theta <- colMeans(x[s == 0, ])
    eta <- c(g.std%*%b)
    
    for (i in 1:n) {
      
      U[1:(4*m + 2),1:(4*m + 2)] <- U[1:(4*m + 2),1:(4*m + 2)] - ipw[i] * tcrossprod(cmat[i,])
      U[(2*m + 3):(3*m + 2),(4*m + 3):(5*m + 2)] <- U[(2*m + 3):(3*m + 2),(4*m + 3):(5*m + 2)] - diag(1 - s[i], m, m)
      U[(3*m + 3):(4*m + 2),(4*m + 3):(5*m + 2)] <- U[(3*m + 3):(4*m + 2),(4*m + 3):(5*m + 2)] - diag(s[i], m, m)
      U[(4*m + 3):(5*m + 2),(4*m + 3):(5*m + 2)] <- U[(4*m + 3):(5*m + 2),(4*m + 3):(5*m + 2)] - diag(1 - s[i], m, m)
      
      V[,1:(4*m + 2)] <- V[,1:(4*m + 2)] - k.std[i]*ipw[i]*(psi[i] - eta[i])*tcrossprod(g.std[i,],cmat[i,])
      V[,(5*m + 3):(5*m + 6)] <- V[,(5*m + 3):(5*m + 6)] - k.std[i]*ipw[i]*tcrossprod(g.std[i,])
      
      meat <- meat + tcrossprod(esteq(p = ipw[i], s = s[i], x = x[i,], psi = psi[i], g.std = g.std[i,],
                                      astar = astar[i], astar2 = astar2[i], eta = eta[i], theta = theta))
      
    }
    
    invbread <- matrix(0, nrow = 5*m + 6, ncol = 5*m + 6)
    invbread[1:(5*m + 2),1:(5*m + 2)] <- U
    invbread[(5*m + 3):(5*m + 6), ] <- V
    
    bread <- try(solve(invbread), silent = TRUE)
    
    if (inherits(bread, "try-error")) {
      
      sandwich <- NA
      variance <- NA
      
    } else {
      
      sandwich <- bread %*% meat %*% t(bread)
      variance <- sandwich[5*m + 5, 5*m + 5] + sandwich[5*m + 6, 5*m + 6] -
        2*sandwich[5*m + 5, 5*m + 6]
      
    }
    
    return(c(mu = mu, sig2 = variance))
    
    
  } else
    return(mu)
  
}

# Kernel weighted least squares
kern_naive <- function(a.new, a, psi, bw = 1, se.fit = FALSE) {
  
  n <- length(a)
  
  # Gaussian Kernel
  a.std <- (a - a.new) / bw
  k.std <- dnorm(a.std) / bw
  g.std <- cbind(1, a.std)
  
  b <- lm(psi ~ -1 + g.std, weights = k.std)$coefficients
  mu <- unname(b[1])
  
  if (se.fit) {
    
    eta <- c(g.std %*% b)
    U <- solve(crossprod(g.std, k.std*g.std))
    V <- cbind(k.std * (psi - eta),
               a.std * k.std * (psi - eta))
    Sig <- U %*% crossprod(V) %*% U
    
    return(c(mu = mu, sig2 = unname(Sig[1,1])))
    
    
  } else
    return(mu)
  
}