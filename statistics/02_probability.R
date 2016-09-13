library(ggplot2)

##############################################################################################
# Sample mean/variance (LNp.1-6-82)
##############################################################################################
# Create population
population <- rnorm(100000, mean=10, sd = 2)

# Create sample data from the population
n <- 10
observe_data <- sample( population, n)

# Calculate sample mean
# sample_mean = 1/n * (x1 + x2 + ... + xn)
# R's mean() function also does the job
(sample_mean <- sum(observe_data)/n)
mean(observe_data)

# Calculate sample variance
# sample_var = 1/(n-1) * [ (x1-mean)^2 + (x2-mean)^2 + ... + (xn-mean)^2 ]
# R's var() function also does the job
(sample_var <- sum( (observe_data-sample_mean)^2 ) / (n-1))
var(observe_data)

# Note that the lower part is n-1, instead of n. This is because sample variance
# is of n-1 degree of freedom. This is exactly what R function var() assumes.
# Here are 2 ways to calculate variance
sum( (observe_data-sample_mean)^2 ) / (n)
mean((observe_data)^2) - mean(observe_data)^2
# Here is the sample variance
var(observe_data)

##############################################################################################
# LLN (LNp.1-6-95~96)
#
#   Let X1, X2, ..., Xn be a sequence of independent random variables with same mean and
# variance. The mean and variance can be estimated by sample mean and sample variance. Note
# that, those random variables are not necessary to be i.i.d.
##############################################################################################
# Create population
population <- rnorm(100000, mean=10, sd = 2)

# Create sample data from the population
n <- 50
observe_data <- sample( population, n)

mean(observe_data)
sqrt(var(observe_data))

##############################################################################################
# Monte Carlo integration (LNp.1-6-95)
##############################################################################################
g <- function(x) {dnorm(x, 0, sd=1)}
integrate(g, 0, 1)

x <- runif(500, min=0, max=1)
data <- data.frame(x=x, gx=g(x))
mean(data$gx)

ggplot(data, aes(x=x, y=gx)) + geom_point() + geom_bar(stat="identity") + geom_bar(stat="identity", width=0.005)

##############################################################################################
# CLT (LNp.1-6-97)
##############################################################################################
# Let X1, ..., Xn be i.i.d. ~ Exponential(1)
# E(Xi) = 1
# Var(Xi) = 1

# Plot the sample mean given n=5. A normal distribution is shown.
n <- 5
sample_mean <- NULL
for( i in c(1:1000) ) {
    sample_mean <- c(sample_mean, mean(rexp(5, 1)))
}
hist(sample_mean)

# Plot the sample mean given n=10. A normal distribution is shown.
n <- 10
sample_mean <- NULL
for( i in c(1:1000) ) {
    sample_mean <- c(sample_mean, mean(rexp(10, 1)))
}
hist(sample_mean)

# Plot the sample mean given n=20. A normal distribution is shown.
n <- 20
sample_mean <- NULL
for( i in c(1:1000) ) {
    sample_mean <- c(sample_mean, mean(rexp(20, 1)))
}
hist(sample_mean)

# Plot the sample mean given n=100. A normal distribution is shown.
n <- 100
sample_mean <- NULL
for( i in c(1:1000) ) {
    sample_mean <- c(sample_mean, mean(rexp(100, 1)))
}
hist(sample_mean)

##############################################################################################
# Normal dist and Binomial dist approximation
##############################################################################################
# Observe how large n can result in approximation of normal pdf and binomial pmf. (LNp.6-36)
# Note that approximation via pdf/pmf is in-correct, as pointed out in later example
data <- data.frame( x=c(0:50))
ggplot(data, aes(x=x, y=dbinom(x, size=5, prob=0.7))) + geom_bar(stat="identity") + stat_function(fun=dnorm, args=list(mean=3.5, sd=1.05^(1/2))) + ylim(0, 0.4)
ggplot(data, aes(x=x, y=dbinom(x, size=10, prob=0.7))) + geom_bar(stat="identity") + stat_function(fun=dnorm, args=list(mean=7, sd=2.1^(1/2))) + ylim(0, 0.4)
ggplot(data, aes(x=x, y=dbinom(x, size=20, prob=0.7))) + geom_bar(stat="identity") + stat_function(fun=dnorm, args=list(mean=14, sd=4.2^(1/2))) + ylim(0, 0.4)
ggplot(data, aes(x=x, y=dbinom(x, size=50, prob=0.7))) + geom_bar(stat="identity") + stat_function(fun=dnorm, args=list(mean=35, sd=10.5^(1/2))) + ylim(0, 0.4)
