##############################################################################################
# Method of finding estimators I - method of moments
#   - Exact sampling distribution
#
# LNp.8-8
# Suppose that X1, ..., Xn are i.i.d. ~ P(lambda).
#                            _
# The estimator of lambda is X
# E(lambda) = lambda
# Var(lambda) = lambda/n
# The sampling distribution of lambda is 1/nP(n*lambda)
# The standard error of lambda is (lambda/n)^(1/2)
##############################################################################################
x <- c(31, 29, 19, 18, 31, 28, 34, 27, 34, 30, 16, 18, 26, 27, 27, 18, 24, 22, 28, 24, 21, 17, 24)
hist(x)

# The estimator of lambda
lambda_estimator <- function(x) { mean(x)}

# The estimate of lambda is 24.91304
(lambda_estimate <- lambda_estimator(x))

# The variance of lambda is 1.083176
# (assume 24.91304 is the real lambda)
lambda_estimate/length(x)

# The standard error lambda estimator is 1.040757
# (assume 24.91304 is the real lambda)
(lambda_estimate/length(x))^(1/2)

# The distribution of lambda estimator
# distribution of 1/n * P(n * lambda).
(parameter <- length(x) * lambda_estimate)
x=c(400:800)
data <- data.frame(x=x, y=x/23, prob=dpois(x, 573))
ggplot(data, aes(x=y, y=prob)) + geom_bar(stat="identity")

# The sum of probability is 1
sum(data$prob)

# The probabily within 2 estimation error is 95.28%
sum(data[data$y>24.91-2.08 & data$y<24.91+2.08,]$prob)

# This is what 1/n * P(n * lambda) means in LNp.8-10
rpois(1, 573) / 23

##############################################################################################
# Method of finding estimators I - method of moments
#   - Asymptotical method
#
# LNp.8-10
# Suppose that X1, ..., Xn are i.i.d. ~ P(lambda).
#
# By CLT, sampling distribution of lambda is approximately normal when n is large enough.
#                            _
# The estimator of lambda is X
# E(lambda) = lambda
# Var(lambda) = lambda/n
# The sampling distribution of lambda is N(lambda, lambda/n)
# The standard error of lambda is (lambda/n)^(1/2)
##############################################################################################
x <- c(31, 29, 19, 18, 31, 28, 34, 27, 34, 30, 16, 18, 26, 27, 27, 18, 24, 22, 28, 24, 21, 17, 24)

# The estimator of lambda
lambda_estimator <- function(x) { mean(x)}

# The estimate of lambda is 24.91304
(lambda_estimate <- lambda_estimator(x))

# The variance of lambda is 1.083176
# (assume 24.91304 is the real lambda)
lambda_estimate/length(x)

# The standard error lambda estimator is 1.040757
# (assume 24.91304 is the real lambda)
(lambda_estimate/length(x))^(1/2)

# The distribution of lambda estimator
lambda_dist2 <- function(x) {dnorm(x, 24.91304, 1.083176)}
(parameter <- length(x) * lambda_estimate)
x=c(400:800)
data <- data.frame(x=x, y=x/23, prob=dpois(x, 573))
ggplot(data, aes(x=y, y=prob)) + geom_bar(stat="identity") + stat_function(fun=lambda_dist2)

# The sum of probability is 1
integrate(lambda_dist2, 20, 30)

# The probabily within 2 estimation error is 94.51%
integrate(lambda_dist2, 24.91-2.08, 24.91+2.08)
# According to Chebyshev's inequality: P(|X-u|>k*sd) <= 1/(k^2)
# The probability of more than 2 standard deviation is smaller than 0.25
# In other words, the probability within 2 standard deviation is at least 0.75
# In this case, by fitting a normal model, same criteria says at least 0.9451 is
# guarantee.
1/4

# This is what 1/n * P(n * lambda) means in LNp.8-10
rpois(1, 573) / 23

##############################################################################################
# Method of finding estimators I - method of moments
#   - Simulation method (bootstrap)
#
# LNp.8-13~14
# Suppose that X1, ..., Xn are i.i.d. ~ Gamma(alpha, lambda).
#                            _
# The estimator of lambda is X/(var^2)
#                            _
# The estimator of alpha is (X^2)/(var^2)
##############################################################################################
# The estimate of lambda and alpha can be obtained (via X = 0.224, variance = 0.1338)
(lambda = 0.224 / 0.1338)
(alpha = (0.224^2) / 0.1338)

# Use the estimate lambda and alpha, to generate Gamma random numbers
# And, calculate the lambda/alpha variance
sim_lambda <- NULL
sim_alpha <- NULL
for( i in c(1:1000) ) {
    (sim_x <- rgamma(227, shape=alpha, rate=lambda))
    (sim_mean <- mean(sim_x))
    (sim_var <- mean(sim_x^2 - sim_mean^2))
    sim_lambda <- c(sim_lambda, sim_mean / sim_var)
    sim_alpha <- c(sim_alpha, (sim_mean^2) / sim_var)
}

# The variance of lambda is 0.12, standard deviation is 0.34
hist(sim_lambda)
(sim_lambda_var <- mean(sim_lambda^2 - mean(sim_lambda)^2))
(sim_lambda_sd <- sim_lambda_var^(1/2))

# The variance of alpha is 0.004, standard deviation is 0.06
hist(sim_alpha)
(sim_alpha_var <- mean(sim_alpha^2 - mean(sim_alpha)^2))
(sim_alpha_sd <- sim_alpha_var^(1/2))
