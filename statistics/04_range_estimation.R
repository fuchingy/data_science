library(ggplot2)  # For ggplot
library(stats4)   # For mle
library(numDeriv) # For grad
library(plyr)     # For ddply

##############################################################################################
# Interval estimation method I - Exact approach
# Pivotal quantative
#
# LNp.8-78
# Suppose that X1, ..., Xn are i.i.d. ~ Normal(mean, variance).
#
# The estimator of mean:            sample_mean
# The standard error of mean      : sqrt(sample_variance/n)
# The pivotable quantitive of mean: sqrt(n) * (sample_mean - mean) / sample_variance ~ t(n-1)
##############################################################################################
# The pivotal quantative of mean
# This is t(n-1) distribution
mean.pvtq <- function(x, mean) {
    n <- length(x)
    sample_mean <- mean(x)
    sample_var <- var(x)
    sqrt(n) * (sample_mean - mean) / sample_var
}

# Generate Normal random number to produce many pivot quantative to observe if it fits T
# distribution
n <- 50
data <- NULL
for( mean in c(50, 51, 52, 53, 54) ) {
    pvtq <- NULL
    for( i in c(1:1000)){
        x <- rnorm(n, mean=mean, sd=12)
        pvtq <- c(pvtq, mean.pvtq(x, mean))
    }
    this_data <- data.frame(mean=mean, pvtq=pvtq)
    data <- rbind(data, this_data)
}
data$req_x <- data$pvtq + data$mean
head(data)
ddply(data, ~mean, summarize, count=length(pvtq))

# Plot pivotal quantative of all different mean LNp.8-78
ggplot(data, aes(x=as.factor(mean), y=pvtq)) + geom_jitter()

# Plot the band concept LNp.8-78
ggplot(data, aes(x=as.factor(mean), y=req_x)) + geom_jitter()

# Plot pivotal quantative of mean
ggplot(data[data$mean==50,], aes(x=pvtq)) + geom_histogram(aes(y = (..count..)/sum(..count..))) + stat_function(fun=dt, args=list(df=50-1))

# Plot T distribution
ggplot(data.frame(x=seq(-5, 5, 0.01)), aes(x=x, y=dt(x, df=n-1))) + geom_line() #+ stat_function(fun=dt, args=list(df=n-1))

# For get the quantile for 0.05 of T distribution
(a <- qt(0.025, df=n-1))
(b <- qt(0.975, df=n-1))

# Calculate the 95% confidence interval of mean
x <- rnorm(n, mean=90, sd=12)
(lower_bound <- -1 * a * var(x) / sqrt(n)) + mean(x)
(upper_bound <- -1 * b * var(x) / sqrt(n)) + mean(x)

# Comparing with the mean estimate and its standard error
mean(x)
sqrt(var(x)/n)

##############################################################################################
# Interval estimation method II - Asymptotic approach
# Asymptotic pivotal quantative
#
# Use the sample Poisson example in LNp.8-8 to demo:
# Suppose that X1, ..., Xn are i.i.d. ~ Poisson(lambda).
#
# Instead of conducting mathmetical deduction, R nlm() function is used.
# The log likelihood function is used here.
##############################################################################################
x <- c(31, 29, 19, 18, 31, 28, 34, 27, 34, 30, 16, 18, 26, 27, 27, 18, 24, 22, 28, 24, 21, 17, 24)
hist(x)

# Use nlm()
LL <- function(lambda, x) {
    sum(-dpois(x, lambda, log = TRUE))
}
out <- nlm(LL, 10, x, hessian = TRUE) # Trun on hessian for Fisher information
out$estimate

# For get the quantile for 0.05 of N(0,1) distribution
(z <- qnorm(0.975))

# The asymptotic 100(1-a)% confidence interval for theta is
# estimated_theta +- z / sqrt(n * Fisher_info)
# Using nlm(), hessian is the Fisher information of n observation data
(fish <- out$hessian[1])

(lower_bound <- out$estimate - z / sqrt(fish))
(upper_bound <- out$estimate + z / sqrt(fish))

##############################################################################################
# Interval estimation method III - Bootstrap approach
#
# Use the example in LNp.8-89 to demo:
# (X1, X2, X3) ~ multinomial(n, p1, p2, p3)
#              ~ multinomial(n, (1-the)^2, 2the(1-the), the^2)
##############################################################################################
# Estimator of theta (LNp.8-24)
estimator_theta <- function (x2, x3, n) {
    (2 * x3 + x2) / (2*n)
}

# Estimate of theta
(estimated_theta <- 0.4247)

# Generate random theta to know the distribution.
# Step1: Pretend the estimated one is the real theta
# Step2: Generate random number using the estimated theta
# Step3: Calculate simulated theta
(p1 <- (1-estimated_theta)^2)
(p2 <- 2 * estimated_theta * (1-estimated_theta))
(p3 <- estimated_theta^2)
randx <- rmultinom(1000, 1029, c(p1, p2, p3))
(eta_theta <- sqrt(randx[3,]/1029))

# Plot the simulated theta and calculate the quantile
hist(eta_theta)
sd(eta_theta)
quantile(eta_theta, 0.025)
quantile(eta_theta, 0.975)

# Calculate and plot the delta
delta <- eta_theta - estimated_theta
hist(delta)

# The simulated 100(1-a)% confidence interval for theta is
# estimated_theta - delta_975 <= theta <= estimated_theta - delta_025
quantile(delta, 0.025)
quantile(delta, 0.975)

(lower_bound <- estimated_theta - quantile(delta, 0.975))
(upper_bound <- estimated_theta - quantile(delta, 0.025))



##############################################################################################
# Interval estimation method III - Bootstrap approach
#
# Use the example in LNp.8-89 to demo. (Also LNp.8-13, LNp.8-28)
##############################################################################################
# Assume already has estimate using MLE (LNp.8-28)
lambda <- 1.96
alpha <- 0.441

# Generate random alpha to know the distribution.
# Step1: Pretend the estimated alpha and lambda are the real theta
# Step2: Generate random number using the estimated theta
# Step3: Calculate simulated theta
sim_lambda <- NULL
sim_alpha <- NULL
for( i in c(1:1000) ) {
    (sim_x <- rgamma(227, shape=alpha, rate=lambda))
    (sim_mean <- mean(sim_x))
    (sim_var <- mean(sim_x^2 - sim_mean^2))
    sim_lambda <- c(sim_lambda, sim_mean / sim_var)
    sim_alpha <- c(sim_alpha, (sim_mean^2) / sim_var)
}

# Plot the simulated alpha and calculate the quantile
hist(sim_alpha)
sd(sim_alpha)
quantile(sim_alpha, 0.05)
quantile(sim_alpha, 0.95)

# Calculate and plot the delta
delta <- sim_alpha - alpha
hist(delta)

# The simulated 100(1-a)% confidence interval for alpha is
# estimated_alpha - delta_95 <= theta <= estimated_alpha - delta_05
quantile(delta, 0.05)
quantile(delta, 0.95)

(lower_bound <- alpha - quantile(delta, 0.95))
(upper_bound <- alpha - quantile(delta, 0.05))



# Plot the simulated lambda and calculate the quantile
hist(sim_lambda)
sd(sim_lambda)
quantile(sim_lambda, 0.05)
quantile(sim_lambda, 0.95)

# Calculate and plot the lambda
delta <- sim_lambda - lambda
hist(delta)

# The simulated 100(1-a)% confidence interval for lambda is
# estimated_lambda - delta_95 <= theta <= estimated_lambda - delta_05
quantile(delta, 0.05)
quantile(delta, 0.95)

(lower_bound <- lambda - quantile(delta, 0.95))
(upper_bound <- lambda - quantile(delta, 0.05))



