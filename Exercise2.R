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
y1_t <- estimation_data$y - mean(estimation_data$y)
y2_t <- estimation_data$c - mean(estimation_data$c)

# Y vector
Y <- rbind(y1_t, y2_t)

# Kalman filter
kalman_filter_factor <- function(theta, Y, return_states = FALSE) {
  h1 <- theta[1]
  h2 <- theta[2]
  f0 <- theta[3]
  f1 <- theta[4]
  f2 <- theta[5]
  r1 <- exp(theta[6])
  r2 <- exp(theta[7])
  q1 <- exp(theta[8])
  q2 <- exp(theta[9])
  
  H <- matrix(c(h1, 1, 0, h2, 0, 1), nrow = 3, ncol = 2)
  F_mat <- diag(c(f0, f1, f2))
  R <- diag(c(r1, r2))
  Q <- diag(c(1, q1, q2))
  n <- ncol(Y)
  
  # Diffuse initialization
  xi_pred <- matrix(0, nrow = 3, ncol = 1)
  P_pred <- diag(3) * 1e6
  loglik <- 0
  
  # Store filtered states
  if (return_states) {
    xi_filt_all <- matrix(NA, nrow = 3, ncol = n)
    xi_pred_all <- matrix(NA, nrow = 3, ncol = n)
    P_filt_all <- array(NA, dim = c(3, 3, n))
  }
  
  for (t in 1:n) {
    # Prediction error
    y_pred <- t(H) %*% xi_pred
    v_t <- Y[, t, drop = FALSE] - y_pred
    S_t <- t(H) %*% P_pred %*% H + R
    
    # Check for near-singularity (anders foutcode)
    if (rcond(S_t) < 1e-14) {
      return(1e10)
    }
    
    tryCatch({
      S_inv <- solve(S_t)
      loglik_contrib <- -0.5 * (2 * log(2 * pi) + log(det(S_t)) + as.numeric(t(v_t) %*% S_inv %*% v_t))
      loglik <- loglik + loglik_contrib
      
      K_t <- P_pred %*% H %*% S_inv
      
      # Update
      xi_filt <- xi_pred + K_t %*% v_t
      P_filt <- P_pred - K_t %*% t(H) %*% P_pred
      
      # Ensure P_filt is symmetric and positive definite
      P_filt <- (P_filt + t(P_filt)) / 2
      
      if (return_states) {
        xi_pred_all[, t] <- xi_pred
        xi_filt_all[, t] <- xi_filt
        P_filt_all[, , t] <- P_filt
      }
      
      # Predict
      xi_pred <- F_mat %*% xi_filt
      P_pred <- F_mat %*% P_filt %*% t(F_mat) + Q
      P_pred <- (P_pred + t(P_pred)) / 2
      
    }, error = function(e) {
      return(1e10)
    })
  }
  
  if (return_states) {
    return(list(loglik = loglik, 
                xi_filt = xi_filt_all,
                xi_pred = xi_pred_all,
                P_filt = P_filt_all,
                xi_T_T = xi_filt,
                P_T_T = P_filt))
  } else {
    return(-loglik)  
  }
}

start_par_2c <- c(
  h1 = 0.5,
  h2 = 0.5,
  f0 = 1,     
  f1 = 1,
  f2 = 1,
  log_r1 = log(1),
  log_r2 = log(1),
  log_q1 = log(1),
  log_q2 = log(1)
)

opt_2c <- optim(
  par = start_par_2c,
  fn = kalman_filter_factor,
  Y = Y,
  return_states = FALSE,
  method = "BFGS",
  hessian = TRUE,
  control = list(
    maxit = 2000,
    trace = 1,  
    REPORT = 50,
    reltol = 1e-8
  )
)

# Extract estimated parameters
h1_hat <- opt_2c$par[1]
h2_hat <- opt_2c$par[2]
f0_hat <- opt_2c$par[3]
f1_hat <- opt_2c$par[4]
f2_hat <- opt_2c$par[5]
r1_hat <- exp(opt_2c$par[6])
r2_hat <- exp(opt_2c$par[7])
q1_hat <- exp(opt_2c$par[8])
q2_hat <- exp(opt_2c$par[9])

# Print results
# Ik weet niet welke idioot dit zo heeft geprint maar serieus onnodig wel (kusjes Thor)
cat("h1 =", round(h1_hat, 4), "\n")
cat("h2 =", round(h2_hat, 4), "\n")
cat("f0 =", round(f0_hat, 4), "\n")
cat("f1 =", round(f1_hat, 4), "\n")
cat("f2 =", round(f2_hat, 4), "\n")
cat("r1 =", round(r1_hat, 4), "\n")
cat("r2 =", round(r2_hat, 4), "\n")
cat("q1 =", round(q1_hat, 4), "\n")
cat("q2 =", round(q2_hat, 4), "\n")
cat("Log-likelihood =", round(-opt_2c$value, 4), "\n")

# Get final filtered state
final_state <- kalman_filter_factor(opt_2c$par, Y, return_states = TRUE)

# Out-of-sample forecasting for K=12 observations
y1_full <- df$y - mean(estimation_data$y)
y2_full <- df$c - mean(estimation_data$c)
Y_full <- rbind(y1_full, y2_full)
K <- 12
n_est <- ncol(Y)
forecast_errors_sq <- matrix(0, nrow = 2, ncol = K)

# Get estimated parameters
H <- matrix(c(h1_hat, 1, 0,
              h2_hat, 0, 1), nrow = 3, ncol = 2)
F_mat <- diag(c(f0_hat, f1_hat, f2_hat))
R <- diag(c(r1_hat, r2_hat))
Q <- diag(c(1, q1_hat, q2_hat))

# Initialize with final filtered state from estimation sample
xi_pred <- F_mat %*% final_state$xi_T_T
P_pred <- F_mat %*% final_state$P_T_T %*% t(F_mat) + Q

for (k in 1:K) {
  t <- n_est + k
  
  # One-step-ahead forecast
  y_forecast <- t(H) %*% xi_pred
  
  # Forecast error
  y_actual <- Y_full[, t, drop = FALSE]
  forecast_error <- y_actual - y_forecast
  forecast_errors_sq[, k] <- forecast_error^2
  
  # Update filter with actual observation
  v_t <- y_actual - t(H) %*% xi_pred
  S_t <- t(H) %*% P_pred %*% H + R
  K_t <- P_pred %*% H %*% solve(S_t)
  
  xi_filt <- xi_pred + K_t %*% v_t
  P_filt <- P_pred - K_t %*% t(H) %*% P_pred
  P_filt <- (P_filt + t(P_filt)) / 2
  
  # Predict next
  xi_pred <- F_mat %*% xi_filt
  P_pred <- F_mat %*% P_filt %*% t(F_mat) + Q
  P_pred <- (P_pred + t(P_pred)) / 2
}

# Mean squared forecast error
msfe_2c <- mean(forecast_errors_sq)
cat("\nMean Squared Forecast Error (MSFE) =", round(msfe_2c, 6), "\n")

#Resultaten Bram
#h1_hat  = 2.6886
#h2_hat  = 1.8869
#r1_hat  = 0.0008
#r2_hat  = 3.4522
#loglik  = -1309.5926
#OOS MSFE (K=12, both series) = 8.0258


# QUESTION 2d -------------------------------------------------------------

# ff lekker kalman filtertje eroverheen
kalman_filter_2d <- function(F_mat, Q, R, Y) {
  d <- nrow(Y)
  Tn <- ncol(Y)
  
  # diffuse initialisation
  m_pred <- matrix(0, nrow = d, ncol = 1)      # xi_{1|0}
  P_pred <- diag(d) * 1e6                      # P_{1|0}
  
  # storage
  m_pred_all <- matrix(NA_real_, nrow = d, ncol = Tn)     # xi_{t|t-1}
  P_pred_all <- array(NA_real_, dim = c(d, d, Tn))        # P_{t|t-1}
  m_filt_all <- matrix(NA_real_, nrow = d, ncol = Tn)     # xi_{t|t}
  P_filt_all <- array(NA_real_, dim = c(d, d, Tn))        # P_{t|t}
  
  loglik <- 0
  
  for (t in 1:Tn) {
    y_t <- Y[, t, drop = FALSE]
    
    # H = I_d
    y_pred <- m_pred
    v_t    <- y_t - y_pred
    S_t    <- P_pred + R   # HPH'+R, but H=I
    
    S_inv <- solve(S_t)
    
    # log-likelihood contribution
    loglik <- loglik - 0.5 * (
      d * log(2 * pi) + log(det(S_t)) + t(v_t) %*% S_inv %*% v_t
    )
    
    # Kalman gain and update
    K_t    <- P_pred %*% S_inv
    m_filt <- m_pred + K_t %*% v_t
    P_filt <- (diag(d) - K_t) %*% P_pred
    P_filt <- 0.5 * (P_filt + t(P_filt))  # symmetrise
    
    # store
    m_pred_all[, t] <- m_pred
    P_pred_all[, , t] <- P_pred
    m_filt_all[, t]   <- m_filt
    P_filt_all[, , t] <- P_filt
    
    # predict next
    m_pred <- F_mat %*% m_filt
    P_pred <- F_mat %*% P_filt %*% t(F_mat) + Q
    P_pred <- 0.5 * (P_pred + t(P_pred))
  }
  
  list(
    loglik     = as.numeric(loglik),
    m_pred_all = m_pred_all,
    P_pred_all = P_pred_all,
    m_filt_all = m_filt_all,
    P_filt_all = P_filt_all
  )
}

# smoothen we die zooi ook nog ff
kalman_smoother_2d <- function(F_mat, Q, kf) {
  m_pred_all <- kf$m_pred_all
  P_pred_all <- kf$P_pred_all
  m_filt_all <- kf$m_filt_all
  P_filt_all <- kf$P_filt_all
  
  d  <- nrow(m_filt_all)
  Tn <- ncol(m_filt_all)
  
  m_smooth <- matrix(NA_real_, nrow = d, ncol = Tn)
  P_smooth <- array(NA_real_, dim = c(d, d, Tn))
  
  # initial (last) smoothed equals filtered
  m_smooth[, Tn]   <- m_filt_all[, Tn]
  P_smooth[, , Tn] <- P_filt_all[, , Tn]
  
  for (t in (Tn-1):1) {
    P_pred_tp1 <- P_pred_all[, , t + 1]
    P_filt_t   <- P_filt_all[, , t]
    
    J_t <- P_filt_t %*% t(F_mat) %*% solve(P_pred_tp1)  # smoother gain
    
    m_smooth[, t] <- m_filt_all[, t] + 
      J_t %*% (m_smooth[, t + 1, drop = FALSE] - m_pred_all[, t + 1, drop = FALSE])
    
    P_smooth[, , t] <- P_filt_t + 
      J_t %*% (P_smooth[, , t + 1] - P_pred_tp1) %*% t(J_t)
    
    P_smooth[, , t] <- 0.5 * (P_smooth[, , t] + t(P_smooth[, , t]))
  }
  
  list(m_smooth = m_smooth, P_smooth = P_smooth)
}

# ff dat klote algoritme runnen
EM_step_2d <- function(F_mat, Q, R, Y) {
  d  <- nrow(Y)
  Tn <- ncol(Y)
  
  # E-step: filter & smoother
  kf <- kalman_filter_2d(F_mat, Q, R, Y)
  ks <- kalman_smoother_2d(F_mat, Q, kf)
  
  m_pred_all <- kf$m_pred_all    # xi_{t|t-1}
  P_pred_all <- kf$P_pred_all
  m_filt_all <- kf$m_filt_all    # xi_{t|t}
  P_filt_all <- kf$P_filt_all
  m_smooth   <- ks$m_smooth      # xi_{t|T}
  P_smooth   <- ks$P_smooth      # P_{t|T}
  
  loglik <- kf$loglik
  
  # Sufficient statistics
  S_xx      <- matrix(0, d, d)   # sum E[x_t x_t']
  S_xx_lag0 <- matrix(0, d, d)   # sum_{t=1}^{T-1} E[x_t x_t']
  S_xx_lag1 <- matrix(0, d, d)   # sum_{t=2}^{T} E[x_t x_{t-1}']
  
  # E[x_t x_t'] terms
  for (t in 1:Tn) {
    mt <- matrix(m_smooth[, t], ncol = 1)
    S_xx <- S_xx + (P_smooth[, , t] + mt %*% t(mt))
    
    if (t < Tn) {
      S_xx_lag0 <- S_xx_lag0 + (P_smooth[, , t] + mt %*% t(mt))
    }
  }
  
  # Cross terms E[x_t x_{t-1}'] using slide formula (17)
  for (t in 2:Tn) {
    mt   <- matrix(m_smooth[, t], ncol = 1)
    mtm1 <- matrix(m_smooth[, t - 1], ncol = 1)
    
    # Cov(x_t, x_{t-1} | Y) = P_t|T P_t|t-1^{-1} F P_{t-1|t-1}
    P_tT     <- P_smooth[, , t]
    P_t_tmin1 <- P_pred_all[, , t]    # P_{t|t-1}
    P_tm1_tm1 <- P_filt_all[, , t - 1]
    
    cross_cov <- P_tT %*% solve(P_t_tmin1) %*% F_mat %*% P_tm1_tm1
    
    Exx_tm1 <- mt %*% t(mtm1) + cross_cov
    
    S_xx_lag1 <- S_xx_lag1 + Exx_tm1
  }
  
  # M-step: update F
  F_new <- S_xx_lag1 %*% solve(S_xx_lag0)
  
  # M-step: update Q
  S11 <- matrix(0, d, d)
  for (t in 2:Tn) {
    mt <- matrix(m_smooth[, t], ncol = 1)
    S11 <- S11 + (P_smooth[, , t] + mt %*% t(mt))
  }
  S10 <- S_xx_lag1
  S00 <- S_xx_lag0
  
  Q_new <- (S11 - F_new %*% t(S10) - S10 %*% t(F_new) + F_new %*% S00 %*% t(F_new)) / (Tn - 1)
  Q_new <- 0.5 * (Q_new + t(Q_new))   # symmetrise
  
  # M-step: update R (H = I)
  R_new <- matrix(0, d, d)
  for (t in 1:Tn) {
    mt  <- m_smooth[, t]
    res <- matrix(Y[, t] - mt, ncol = 1)
    R_new <- R_new + (res %*% t(res) + P_smooth[, , t])
  }
  R_new <- R_new / Tn
  R_new <- 0.5 * (R_new + t(R_new))
  
  list(F = F_new, Q = Q_new, R = R_new, loglik = loglik)
}


# habba
d  <- nrow(Y)
Tn <- ncol(Y)

# Starting values from assignment
F_em <- 0.8 * diag(d)
R_em <- 2   * diag(d)
Q_em <- 3   * diag(d)

maxiter <- 200
for (iter in 1:maxiter) {
  
  step <- EM_step_2d(F_em, Q_em, R_em, Y)
  
  F_em   <- step$F
  Q_em   <- step$Q
  R_em   <- step$R
  loglik <- step$loglik
  
  if (iter %in% c(20, 200)) {
    cat("  --- Snapshot at iteration", iter, "---\n")
    cat("  F:\n"); print(F_em)
    cat("  Q:\n"); print(Q_em)
    cat("  R:\n"); print(R_em)
    cat("  loglik =", loglik, "\n")
  }
}

