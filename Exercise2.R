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
y_t <- y - mu

kalman_loglik_ar1_noise <- function(theta, y) {
  F  <- theta[1]
  Q  <- exp(theta[2])       
  R  <- exp(theta[3])       
  
  n  <- length(y)
  
  # Initialization
  mu_pred <- 0
  P_pred <- Q / (1 - F^2)
  loglik <- 0
  
  for (t in 1:n) {
    # 1-step ahead forecast of y_t
    y_hat <- mu_pred
    S_t   <- P_pred + R      
    v_t   <- y[t] - y_hat   
    
    # add log-density
    loglik <- loglik - 0.5 * (log(2 * pi) + log(S_t) + (v_t^2 / S_t))
    
    # Kalman gain
    K_t <- P_pred / S_t
    
    # Updated state mu_t|t and variance P_t|t
    mu_filt <- mu_pred + K_t * v_t
    P_filt  <- (1 - K_t) * P_pred
    
    # Predict next state mu_{t+1|t}, P_{t+1|t}
    mu_pred <- F * mu_filt
    P_pred  <- F^2 * P_filt + Q
  }
  
  return(-loglik)
}

start_par <- c(F = 0.1, logQ = log(4), logR = log(8))

opt_2b <- optim(
  par    = start_par,
  fn     = kalman_loglik_ar1_noise,
  y      = y_t,
  method = "BFGS",
  hessian = TRUE
)
opt_2b

# Extract estimated parameters
F_hat  <- opt_2b$par[1]
Q_hat  <- exp(opt_2b$par[2])
R_hat  <- exp(opt_2b$par[3])
F_hat
Q_hat
R_hat

# Get xi_T|T and P_T|T
kalman_filter_ar1_noise <- function(F, Q, R, y) {
  n <- length(y)
  
  mu_pred <- 0
  P_pred <- Q / (1 - F^2)
  
  mu_filt <- 0
  P_filt  <- 0
  
  for (t in 1:n) {
    # Forecast
    y_hat <- mu_pred
    S_t   <- P_pred + R
    v_t   <- y[t] - y_hat
    
    # Kalman gain
    K_t <- P_pred / S_t
    
    # Update
    mu_filt <- mu_pred + K_t * v_t
    P_filt  <- (1 - K_t) * P_pred
    
    # Predict next
    mu_pred <- F * mu_filt
    P_pred  <- F^2 * P_filt + Q
  }
  
  # Return last filtered state
  list(mu_TT = mu_filt, P_TT = P_filt)
}

kf_2b <- kalman_filter_ar1_noise(F_hat, Q_hat, R_hat, y_t)

xi_T_T_hat <- kf_2b$mu_TT
P_T_T_hat  <- kf_2b$P_TT
xi_T_T_hat
P_T_T_hat


# QUESTION 2c -------------------------------------------------------------
c <- estimation_data$c
mu_c <- mean(c)
y_2 <- c - mu_c

Y <- cbind(y_1, y_2)
