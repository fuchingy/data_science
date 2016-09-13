##############################################################################################
# Method of finding estimators I - method of moments
#
# LNp.8-8
# Suppose that X1, ..., Xn are i.i.d. ~ Poisson(lambda).
#
# The mean of Poisson(lambda): lambda
# The variance of Poisson(lambda): lambda
##############################################################################################
x <- c(31, 29, 19, 18, 31, 28, 34, 27, 34, 30, 16, 18, 26, 27, 27, 18, 24, 22, 28, 24, 21, 17, 24)
hist(x)

# First sample moment
(first_sample_moment <- sum(x)/length(x))

# Second sample moment
(second_sample_moment <- sum(x^2)/length(x))

# Since the lambda is the mean, whose estimator is exactly the first sample moment, so the
# best estimator for lambda is the first sample moment.
# Put in the observed data, the estimate is 24.91304
(lambda_estimate <- first_sample_moment)

##############################################################################################
# Estimate the variance of the estimate
#   - Exact sampling distribution
#   - Asymptotical method
#   - Simulation method (bootstrap)
#
# Continue with the previous example.
##############################################################################################
# 1. Exact sampling distribution
# By the calculation (LNp.8-10):
#   The sampling distribution of lambda:  1/nP(n*lambda)
#   The standard error of lambda       :  (lambda/n)^(1/2)
sqrt((lambda_estimate/length(x))) # 1.040757

# Plot the sampling distribution of lambda estimator
# Note that 1/n * P(n * lambda) actually means rpois(1, 573) / 23 (LNp.8-10)
(parameter <- length(x) * lambda_estimate)
x=c(400:800)
data <- data.frame(x=x, y=x/23, prob=dpois(x, 573))
ggplot(data, aes(x=y, y=prob)) + geom_bar(stat="identity")

# The probabily within 2 estimation error is 95.28%
sum(data[data$y>24.91-2.08 & data$y<24.91+2.08,]$prob)

# 2. Asymptotical method
# By LLN,
#   E(sample mean)   : population mean (LNp.1~6.95)
#   Var(sample mean) : population variance/n (LNp.1~6.95)
#   Var(sample mean) : sampling variance/n (LNp.1~6.96)
(lambda_estimated <- mean(x))
# The variance of the sample mean is variance/n
var(x)/length(x)
# The standard error of the sample mean is sqrt(variance/n)
(std_error <- sqrt(var(x)/length(x))) # 1.143659
#sqrt(lambda_estimated/length(x)) # 1.143659

# Plot the sampling distribution of lambda estimator
ggplot(data, aes(x=y, y=prob)) + geom_bar(stat="identity") + stat_function(fun=dnorm, args=list(mean=lambda_estimated, sd=std_error), colour="blue")

# 3. Simulation method
# Use random number generator to generate sampling distribution.
# The estimated parameters is from other approachs, such as moment method, MLE.
# Here, I use the estimated parameter from moment method
r_mean <- NULL
for( i in c(1:1000) ) {
    r_mean <- c(r_mean, mean(rpois(23, lambda=lambda_estimated)))
}
hist(r_mean)
var(r_mean)

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
# First sample moment
first_sample_moment <- 0.224

# Second sample moment
second_sample_moment <- 0.183976

# Estimated variance
(estimated_variance <- second_sample_moment - first_sample_moment^2)

# The estimate of lambda and alpha can be obtained (via X = 0.224, variance = 0.1338)
(lambda = first_sample_moment / estimated_variance)
(alpha = (first_sample_moment^2) / estimated_variance)

##############################################################################################
# Estimate the variance of the estimate
#   - Simulation method (bootstrap)
#
#   Exact sampling distribution is difficult to get in this example.
#   Asymptotical method is also difficult, since involving dividend: lambda=sample_mean/variance.
#
# Continue with the previous example.
##############################################################################################

# 2. Simulation method
# Use random number generator to generate sampling distribution.
# The estimated parameters is from other approachs, such as moment method, MLE.
# Here, I use the estimated parameter from moment method
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
(sim_lambda_sd <- sqrt(sim_lambda_var))

# The variance of alpha is 0.004, standard deviation is 0.06
hist(sim_alpha)
(sim_alpha_var <- mean(sim_alpha^2 - mean(sim_alpha)^2))
(sim_alpha_sd <- sqrt(sim_alpha_var))

##############################################################################################
# Method of finding estimators II - Maximum Likelihood Estimator (MLE)
#   - ?
#
# LNp.8-16
# Suppose that X1, ..., Xn are i.i.d. ~ P(lambda).
#                            _
##############################################################################################
dbinom(7, 10, 0.1)
dbinom(7, 10, 0.5)
dbinom(7, 10, 0.7)
dbinom(7, 10, 0.9)
