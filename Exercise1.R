# Here is the solution to exercise 1
# First one has to download the data from the 'Data.R' file

# QUESTION 1a -------------------------------------------------------------
y <- estimation_data$y
mu <- mean(y)
std <- sd(y)
p1 <- 0.5
sigma_1 <- 0.5 * std
sigma_2 <- std

nll <- function(par, y) {
  mu       <- par[1]
  sigma_1  <- par[2]
  sigma_2  <- par[3]
  p1       <- par[4]
  f <- p1 * dnorm(y, mean = mu, sd = sigma_1) +
    (1 - p1) * dnorm(y, mean = mu, sd = sigma_2)
  -sum(log(f))
}

start <- c(mu, sigma_1, sigma_2, p1)

mle <- optim(par = start, fn = nll, y = y,
             method = "L-BFGS-B", #deze moeten we ff want anders kan je die bounds er niet op zetten
             lower = c(-Inf, 0, 0, 0),      
             upper = c( Inf,  Inf,  Inf, 1)
             )

# calculate P(s_t=2|y_t=0)
mu_hat  <- mle$par[1]
s1_hat  <- mle$par[2]
s2_hat  <- mle$par[3]
p1_hat  <- mle$par[4]

y0 <- 0
num <- (1 - p1_hat) * dnorm(y0, mean = mu_hat, sd = s2_hat)
den <- p1_hat * dnorm(y0, mean = mu_hat, sd = s1_hat) + num
prob <- as.numeric(num / den)

print(s1_hat)
print(s2_hat)
print(prob)

# QUESTION 1b -------------------------------------------------------------
# function for Hamilton filter
update_step <- function(y, p11, p22, mu, sigma1, sigma2, xi) {
  # transition matrix
  P <- matrix(c(p11, 1 - p22,
                1 - p11, p22), nrow = 2)
  
  # initialization
  n <- length(y)
  xi_t  <- matrix(NA_real_, nrow = n, ncol = 2) 
  ll_total  <- 0                             
  
  for (t in seq_len(n)) {
    # prediction step
    xi_pred <- as.numeric(P %*% xi)
    
    # densities
    eta1 <- dnorm(y[t], mean = mu, sd = sigma1)
    eta2 <- dnorm(y[t], mean = mu, sd = sigma2)
    
    # predictive density for log-likelihood (use xi_pred, not xi)
    ft <- (eta1 * xi_pred[1]) + (eta2 * xi_pred[2])
    ll_total <- ll_total + log(ft)
    
    # update (filter)
    num <- c(eta1, eta2) * xi_pred
    den <- sum(num)
    xi  <- num / den
    
    xi_t[t, ] <- xi
  }
  list(logLik = ll_total, xi_t = xi_t)
}

# initialize
p11 <- 0.5
p22 <- 0.5

# three initializations as per assignment
xi1 <- c(1, 0)
xi2 <- c(0, 1)
xi3 <- c(0.5, 0.5)

# negative log-likelihood calculation
nll_ms <- function(par, y, xi0) {
  p11 <- par[1]; p22 <- par[2]; mu <- par[3]; s1 <- par[4]; s2 <- par[5]
  res <- update_step(y, p11, p22, mu, s1, s2, xi0)  # your filter
  -res$logLik
}

# bounds
lower <- c(0,    0,    -Inf, 1e-6, 1e-6)
upper <- c(1,    1,     Inf,  Inf,  Inf)

# start params
start <- c(p11, p22, mu, sigma_1, sigma_2)

fit1 <- optim(start, nll_ms, y = y, xi0 = xi1, method = "L-BFGS-B", lower = lower, upper = upper)
fit2 <- optim(start, nll_ms, y = y, xi0 = xi2, method = "L-BFGS-B", lower = lower, upper = upper)
fit3 <- optim(start, nll_ms, y = y, xi0 = xi3, method = "L-BFGS-B", lower = lower, upper = upper)

print(fit1$par)
print(fit1$value)
print(fit2$par)
print(fit2$value)
print(fit3$par)
print(fit3$value)
# QUESTION 1c -------------------------------------------------------------
# Extract estimated p11 and p22 with initialisation in state 1
p11_hat_1 <- fit1$par[1]
p22_hat_1 <- fit1$par[2]

print(p11_hat_1)
print(p22_hat_1)
# Calculate steady state probability of being in state 2: (1-p11)/(2-p11-p22)
steady_state_2 <- (1-p11_hat_1) / (2-p11_hat_1-p22_hat_1)
print(steady_state_2)
# Expected duration of a high volatility episode in quarters (to get in years, divide by 4)
expected_duration_2 <- 1 / (1 - p22_hat_1)
print(expected_duration_2)
# 1-step ahead OOS forecast



# QUESTION 1d -------------------------------------------------------------
library(mvtnorm)

# Hamilton filter multivariate
Hamilton_filter_mv <- function(p11, p22, mu_list, Sigma_list, xi0_in, Y) {
  Y   <- as.matrix(Y)
  Tn  <- nrow(Y)
  K   <- 2
  
  P <- matrix(c(p11, 1 - p22,
                1 - p11, p22), nrow = 2, byrow = TRUE)
  
  predictedxi <- matrix(NA_real_, nrow = K, ncol = Tn + 1)
  filteredxi  <- matrix(NA_real_, nrow = K, ncol = Tn)
  likelihood  <- matrix(NA_real_, nrow = K, ncol = Tn)
  loglik      <- 0
  
  # initial prediction
  predictedxi[, 1] <- P %*% xi0_in
  
  for (t in 1:Tn) {
    # state-conditional densities
    for (k in 1:K) {
      likelihood[k, t] <- dmvnorm(Y[t, ], mean = mu_list[[k]], sigma = Sigma_list[[k]])
    }
    
    num <- predictedxi[, t] * likelihood[, t]
    den <- sum(num)
    
    filteredxi[, t] <- num / den
    loglik <- loglik + log(den)
    
    predictedxi[, t + 1] <- P %*% filteredxi[, t]
  }
  
  predictedxi <- predictedxi[, 1:Tn, drop = FALSE]
  
  list(filteredxi  = filteredxi,
       predictedxi = predictedxi,
       loglik      = loglik,
       P           = P)
}

# Hamilton smoother multivariate
Hamilton_smoother_mv <- function(p11, p22, mu_list, Sigma_list, xi0_in, Y) {
  hf <- Hamilton_filter_mv(p11, p22, mu_list, Sigma_list, xi0_in, Y)
  filteredxi  <- hf$filteredxi  
  predictedxi <- hf$predictedxi  
  loglik      <- hf$loglik
  P           <- hf$P
  
  Tn <- ncol(filteredxi)
  K  <- nrow(filteredxi)
  
  smoothedxi <- matrix(NA_real_, nrow = K, ncol = Tn)
  smoothedxi[, Tn] <- filteredxi[, Tn]
  
  for (t in (Tn - 1):1) {
    tmp <- t(P) %*% (smoothedxi[, t + 1] / predictedxi[, t + 1])
    smoothedxi[, t] <- filteredxi[, t] * as.numeric(tmp)
  }
  
  tmp0    <- t(P) %*% (smoothedxi[, 1] / predictedxi[, 1])
  xi0_out <- xi0_in * as.numeric(tmp0)
  
  # cross terms Pstar(i,j,t) = Pr(S_{t-1}=i, S_t=j | Y)
  Pstar <- array(NA_real_, dim = c(K, K, Tn))
  
  # t = 1 uses xi0_in
  denom1 <- as.numeric(predictedxi[, 1] %*% c(1, 1))
  num1   <- smoothedxi[, 1, drop = FALSE] %*% t(xi0_in)
  Pstar[, , 1] <- P * (num1 / denom1)
  
  for (t in 2:Tn) {
    denom <- as.numeric(predictedxi[, t] %*% c(1, 1))
    num   <- smoothedxi[, t, drop = FALSE] %*% 
      t(filteredxi[, t - 1, drop = FALSE])
    Pstar[, , t] <- P * (num / denom)
  }
  
  list(smoothedxi = smoothedxi,
       xi0_out    = xi0_out,
       Pstar      = Pstar,
       loglik     = loglik)
}

EM_step_mv <- function(p11, p22, mu1, mu2, Sigma1, Sigma2, xi0_in, Y) {
  Y   <- as.matrix(Y)
  Tn  <- nrow(Y)
  
  mu_list    <- list(mu1, mu2)
  Sigma_list <- list(Sigma1, Sigma2)
  
  sm <- Hamilton_smoother_mv(p11, p22, mu_list, Sigma_list, xi0_in, Y)
  smoothedxi <- sm$smoothedxi   
  xi0_out    <- sm$xi0_out
  Pstar      <- sm$Pstar
  loglik     <- sm$loglik
  
  # extract p*_ij(t)
  p11star <- Pstar[1, 1, ]
  p12star <- Pstar[1, 2, ]
  p21star <- Pstar[2, 1, ]
  p22star <- Pstar[2, 2, ]
  
  p1star  <- p11star + p12star     # Pr(S_{t-1} = 1 | Y)
  p2star  <- p21star + p22star     # Pr(S_{t-1} = 2 | Y)
  
  # M-step: transition probabilities (raw, from theory)
  p11_out_raw <- sum(p11star) / (xi0_out[1] + sum(p1star[1:(Tn - 1)]))
  p22_out_raw <- sum(p22star) / (xi0_out[2] + sum(p2star[1:(Tn - 1)]))
  
  
  # Enforce probability bounds (consistent with Markov chain theory)
  p11_out <- min(max(p11_out_raw, 0), 1)
  p22_out <- min(max(p22_out_raw, 0), 1)
  
  # M-step: means (2D vectors)
  mu1_out <- colSums(Y * p1star) / sum(p1star)
  mu2_out <- colSums(Y * p2star) / sum(p2star)
  
  # M-step: covariance matrices (2x2)
  Y1c <- sweep(Y, 2, mu1_out, "-")
  Y2c <- sweep(Y, 2, mu2_out, "-")
  
  Sigma1_out <- t(Y1c) %*% (Y1c * p1star) / sum(p1star)
  Sigma2_out <- t(Y2c) %*% (Y2c * p2star) / sum(p2star)
  
  # Check covariance matrices are valid (theoretical requirement)
  det1 <- det(Sigma1_out)
  det2 <- det(Sigma2_out)
  if (!is.finite(det1) || !is.finite(det2) || det1 <= 0 || det2 <= 0) {
    cat("Invalid covariance matrix in EM_step_mv\n")
    cat("det(Sigma1_out) =", det1, "\n")
    cat("det(Sigma2_out) =", det2, "\n")
    cat("Sigma1_out =\n"); print(Sigma1_out)
    cat("Sigma2_out =\n"); print(Sigma2_out)
    stop("Covariance matrix not positive definite")
  }
  
  # Ensure xi0_out is a proper probability vector
  if (any(xi0_out < 0) || sum(xi0_out) <= 0) {
    cat("xi0_out not a valid probability vector:\n")
    print(xi0_out)
    stop("xi0_out invalid")
  }
  xi0_out <- xi0_out / sum(xi0_out)
  
  list(
    p11    = p11_out,
    p22    = p22_out,
    mu1    = mu1_out,
    mu2    = mu2_out,
    Sigma1 = Sigma1_out,
    Sigma2 = Sigma2_out,
    xi0    = xi0_out,
    loglik = loglik
  )
}

# EM loop
em_ms_mv <- function(Y,
                     p11_init, p22_init,
                     mu1_init, mu2_init,
                     Sigma1_init, Sigma2_init,
                     xi0_init = c(0.5, 0.5),
                     maxiter = 200, tol = 1e-6) {
  
  p11    <- p11_init
  p22    <- p22_init
  mu1    <- mu1_init
  mu2    <- mu2_init
  Sigma1 <- Sigma1_init
  Sigma2 <- Sigma2_init
  xi0    <- xi0_init
  
  loglik_old <- -Inf
  # check outcome for iterations
  for (iter in 1:maxiter) {
    cat("EM iteration:", iter, "\n")
    step <- EM_step_mv(p11, p22, mu1, mu2, Sigma1, Sigma2, xi0, Y)
    
    p11    <- step$p11
    p22    <- step$p22
    mu1    <- step$mu1
    mu2    <- step$mu2
    Sigma1 <- step$Sigma1
    Sigma2 <- step$Sigma2
    xi0    <- step$xi0
    loglik <- step$loglik
    
    cat("  loglik =", loglik, "\n")
    cat("  p11 =", p11, " p22 =", p22, "\n")
    
    if (!is.finite(loglik)) {
      cat("Non-finite loglik at iteration", iter, "\n")
      stop("loglik is NA/NaN/Inf")
    }
    
    if (abs(loglik - loglik_old) < tol) {
      cat("Converged at iteration", iter, "\n")
      break
    }
    loglik_old <- loglik
  }
  
  list(
    p11    = p11,
    p22    = p22,
    mu1    = mu1,
    mu2    = mu2,
    Sigma1 = Sigma1,
    Sigma2 = Sigma2,
    xi0    = xi0,
    loglik = loglik,
    iter   = iter
  )
}

c <- estimation_data$c
Y <- cbind(y, c)

# initialize
p11 <- 0.8
p22 <- 0.8

# Starting values 
mu <- colMeans(Y)
Sigma_1 <- cov(Y)
Sigma_2 <- 0.5 * cov(Y)

fit <- em_ms_mv(
  Y           = Y,
  p11_init    = p11,
  p22_init    = p22,
  mu1_init    = mu,
  mu2_init    = mu,
  Sigma1_init = Sigma_1,
  Sigma2_init = Sigma_2,
  xi0_init    = c(0.5, 0.5)
)

fit$p11
fit$p22
fit$mu1
fit$mu2
fit$Sigma1
fit$Sigma2
fit$loglik
