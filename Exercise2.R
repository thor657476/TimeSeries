# Here is the solution to exercise 2
# First one has to download the data from the 'Data.R' file

# QUESTION 2a -------------------------------------------------------------
y <- estimation_data$y
mu <- mean(y)
y_1 <- y - mu

# estimate the AR(1) model using OLS
y_1lagged <- y_1[1 : (length(y_1) - 1)]
y_1t <- y_1[2 : length(y_1)]

#zero intercept
ar_1 <-lm(y_1t ~ 0 + y_1lagged)
summary(ar_1)
sigma_ups_hat <- summary(ar_1)$sigma^2
print(sigma_ups_hat)

# QUESTION 2b -------------------------------------------------------------
# Kalman filter 
kalman_step <- function(y_1, F, Q, R, xi0, P0) {
  T <- length(y_1)
  
  # storage
  predictedxi <- rep(0, T)  
  predictedP  <- rep(0, T)  
  xi          <- rep(0, T)  
  P           <- rep(0, T)  
  
  # first prediction based on xi0 and P0 
  predictedxi[1] <- xi0
  predictedP[1]  <- P0
  
  # first updating step
  xi[1] <- predictedxi[1] + (predictedP[1]*1 *(y[1] -1*predictedxi[1]))/ (predictedP[1 + R])
  P[1]  <- predictedP[1]  - (predictedP[1]*1*predictedP[1])/(predictedP[1] + R)
  
  # log-likelihood
  ll_total <- -0.5 * (log(2 * pi) + log(F1) + v1^2 / F1)
  
  # Kalman filter for t >= 2
  for (t in 2:T) {
    # prediction step
    predictedxi[t] <- F * xi[t - 1]
    predictedP[t]  <- F^2 * P[t - 1] + Q
    
    # prediction error
    vt <- y_1[t] - predictedxi[t]
    Ft <- predictedP[t] + R      # H = 1
    Kt <- predictedP[t] / Ft
    
    # update step
    xi[t] <- predictedxi[t] + Kt * vt
    P[t]  <- predictedP[t]  - Kt * predictedP[t]
    
    # log-likelihood
    ll_total <- ll_total - 0.5 * (log(2 * pi) + log(Ft) + vt^2 / Ft)
  }
  
  list(
    logLik = ll_total,
    xi_t_t = xi,    #test      
    P_t_t  = P,           
    xi_pred = predictedxi,
    P_pred  = predictedP 
  )
}

nll_kalman <- function(par, y_1, xi0, P0) {
  F <- par[1]  
  Q <- par[2]
  R <- par[3]
  
  res <- kalman_step(y_1, F, Q, R, xi0, P0)
  -res$logLik
}
start_par <- c(0.1, 4, 8)  # (phi_start, Q_start, R_start)
P0 <- 8 /( 1- (0.1)^2)
xi0 <- 0
lower_par <- c(-Inf, 1e-6, 1e-6)
upper_par <- c( Inf,  Inf,  Inf)

fit_kalman <- optim(
  par     = start_par,
  fn      = nll_kalman,
  y_1     = y_1,
  xi0     = xi0,
  P0      = P0,
  method  = "L-BFGS-B",
  lower   = lower_par,
  upper   = upper_par
)

F_hat <- fit_kalman$par[1]
Q_hat   <- fit_kalman$par[2]
R_hat   <- fit_kalman$par[3]

kf_res <- kalman_step(y_1, F_hat, Q_hat, R_hat, xi0, P0)

xi_t_t <- kf_res$xi_t_t    # xi_{T|T}
P_t_t  <- kf_res$P_t_t     # P_{T|T}

print(F_hat) # phi hat
print(Q_hat) # sigma_eps^2
print(R_hat) # sigma_nu^2

print(xi_t_t)
print(P_t_t)
# QUESTION 2c -------------------------------------------------------------
c <- estimation_data$c
mu_c <- mean(c)
y_2 <- c - mu_c

Y <- cbind(y_1, y_2)
