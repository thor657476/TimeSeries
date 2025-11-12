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

mle <- optim(par = start, fn = nll, y = y)

# calculate P(s_t=2|y_t=0)
mu_hat  <- mle$par[1]
s1_hat  <- mle$par[2]
s2_hat  <- mle$par[3]
p1_hat  <- mle$par[4]

y0 <- 0
num <- (1 - p1_hat) * dnorm(y0, mean = mu_hat, sd = s2_hat)
den <- p1_hat * dnorm(y0, mean = mu_hat, sd = s1_hat) + num
prob <- as.numeric(num / den)




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

# start params
start <- c(p11, p22, mu, sigma_1, sigma_2)

fit1 <- optim(start, nll_ms, y = y, xi0 = xi1)
fit2 <- optim(start, nll_ms, y = y, xi0 = xi2)
fit3 <- optim(start, nll_ms, y = y, xi0 = xi3)