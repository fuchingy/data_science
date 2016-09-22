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
for( mean in c(50, 70, 90, 110, 130) ) {
    pvtq <- NULL
    for( i in c(1:1000)){
        x <- rnorm(n, mean=mean, sd=12)
        pvtq <- c(pvtq, mean.pvtq(x, mean))
    }
    this_data <- data.frame(mean=mean, pvtq=pvtq)
    data <- rbind(data, this_data)
}
data$req_x <- data$pvtq / data$mean
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
