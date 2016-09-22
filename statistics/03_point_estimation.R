library(ggplot2)
library(stats4) # for mle
library(numDeriv) # For grad

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
#
# Use the sample Poisson example in LNp.8-8 to demo:
# Suppose that X1, ..., Xn are i.i.d. ~ Poisson(lambda).
#
# Instead of conducting mathmetical deduction, R mle()/mle() function is used instead.
# The log likelihood function is used here.
##############################################################################################
x <- c(31, 29, 19, 18, 31, 28, 34, 27, 34, 30, 16, 18, 26, 27, 27, 18, 24, 22, 28, 24, 21, 17, 24)
hist(x)

# Use mle()
LL <- function(lambda) {
    R = dpois(x, lambda)
    -sum(log(R))
}
mle(LL, start = list(lambda = 10))

# Use nlm()
LL <- function(lambda, x) {
    sum(-dpois(x, lambda, log = TRUE))
}
out <- nlm(LL, 10, x, hessian = TRUE) # Trun on hessian for Fisher information
out$estimate

# {TBD}
(fish <- out$hessian)
solve(fish)

##############################################################################################
# Asymptotic theory for method of MLE
#
# Use the sample Poisson example in LNp.8-8 to demo:
# Suppose that X1, ..., Xn are i.i.d. ~ Poisson(lambda).
#
##############################################################################################
x <- c(31, 29, 19, 18, 31, 28, 34, 27, 34, 30, 16, 18, 26, 27, 27, 18, 24, 22, 28, 24, 21, 17, 24)

# The log likelihood function for 1 random variable
LL_1 <- function(x, lambda) {
    dpois(x, lambda, log = TRUE)
}

# Observe how log likelihood function varies with different parameters, given different data
# Same figure in LNp.8-33.
lambda <- c(0:90)
data <- data.frame(lambda=lambda, prob=LL_1(31, lambda))
p <- ggplot(data, aes(x=lambda, y=prob)) + geom_line(colour="blue")

data <- data.frame(lambda=lambda, prob=LL_1(29, lambda))
p <- p + geom_line(data=data, aes(x=lambda, y=prob), colour="green")

data <- data.frame(lambda=lambda, prob=LL_1(19, lambda))
p <- p + geom_line(data=data, aes(x=lambda, y=prob), colour="purple")

data <- data.frame(lambda=lambda, prob=LL_1(18, lambda))
p <- p + geom_line(data=data, aes(x=lambda, y=prob), colour="yellow")

data <- data.frame(lambda=lambda, prob=LL_1(28, lambda))
p <- p + geom_line(data=data, aes(x=lambda, y=prob), colour="orange") + stat_function(fun=LL_joint)
p

##############################################################################################
# Insect counts data, Bliss and Fisher (1953) LNp.67
#
# Use the sample Poisson example in LNp.8-8 to demo:
# Suppose that X1, ..., Xn are i.i.d. ~ Poisson(lambda).
#
##############################################################################################
# Create the sample data
data <- data.frame(bug_per_leaf=c(0:8), obs_cnt=c(70, 38, 17, 10, 9, 3, 2, 1, 0))
data$prob <- data$obs_cnt / sum(data$obs_cnt)
data

# Histgram
ggplot(data, aes(x=bug_per_leaf, y=obs_cnt)) + geom_bar(stat="identity")
ggplot(data, aes(x=bug_per_leaf, y=prob)) + geom_bar(stat="identity")

# Fit Poisson statistical model
# The lambda estimator is sample mean
(sample_mean <- sum(data$bug_per_leaf * data$prob))
# Create estimation data following fitted Poisson model
data$pois_eta_cnt <- 150 * dpois(data$bug_per_leaf, lambda=sample_mean)

# Observe how good Poisson fit
ggplot(data, aes(x=bug_per_leaf, y=obs_cnt)) + geom_bar(stat="identity") + geom_line(aes(y=pois_eta_cnt), colour="blue")

# Fit Generalized Negative Binomial (GNB) statistical model
# Since I can't find built-in GNB model, I create a function as its pmf
dgnb <- function(x, m, k) {
    ((1+m/k)^(-k)) * (gamma(k+x)/factorial(x)/gamma(k)) * (m/(m+k))^x
}

# MLE is choosen as the estimators of m and k.
# The joint log likelihood function is created.
#   Since the joint log likelihood function needs observation data, instead of
# frequency, I recreate the observation data as x.
x <- NULL
for( row_i in c(1:nrow(data)) ) {
    x <- c(x, rep(data[row_i,]$bug_per_leaf, data[row_i,]$obs_cnt))
}
x

# Create joint log likelihood function
LL <- function(m, k) {
    R = dgnb(x, m, k)
    -sum(log(R))
}

# Use mle() to estimate m and k
mle(LL, start = list(m = 10, k = 10))

# Create estimation data following fitted Generalized Negative Binomial model
data$gnb_eta_cnt <- 150 * dgnb(data$bug_per_leaf, m=1.14881, k=1.024593)
data

# Observe how good Poisson and GNB fit
ggplot(data, aes(x=bug_per_leaf, y=obs_cnt)) + geom_bar(stat="identity") +
    geom_line(aes(y=pois_eta_cnt), colour="blue") +
    geom_line(aes(y=gnb_eta_cnt), colour="red")
