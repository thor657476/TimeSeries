# Here is the solution to exercise 1
# First one has to download the data from the 'Data.R' file

### 1a
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

