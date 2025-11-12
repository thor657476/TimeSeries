# Here is the solution to exercise 1
# First one has to download the data from the 'Data.R' file

### 1a

y <- estimation_data$y
mu <- mean(y)
std <- sd(y)
p1 <- 0.5
sigma_1 <- 0.5 * std
sigma_2 <- std

nll <- function(mu, sigma_1, sigma_2, p1, y) {
  f <- p1 * dnorm(y, mean = mu, sd = sigma_1) +
    (1 - p1) * dnorm(y, mean = mu, sd = sigma_2)
  -sum(log(f))
}

mle <- optim(par = c(mu, sigma_1, sigma_2, p1), fn = nll, y = y)

# print estimates
mle$par

