##############################################################################################
# Sample mean/variance (LNp.1-6-82,96)
##############################################################################################
omega <- rnorm(100000, mean=10, sd = 2)

(sample_mean <- mean(sample( omega, 1)))
(sample_mean <- mean(sample( omega, 5)))
(sample_mean <- mean(sample( omega, 10)))
(sample_mean <- mean(sample( omega, 50)))
(sample_mean <- mean(sample( omega, 100)))
(sample_mean <- mean(sample( omega, 1000)))
(sample_mean <- mean(sample( omega, 10000)))

(sample_variance <- mean((sample( omega, 1) - sample_mean)^2))
(sample_variance <- mean((sample( omega, 5) - sample_mean)^2))
(sample_variance <- mean((sample( omega, 10) - sample_mean)^2))
(sample_variance <- mean((sample( omega, 50) - sample_mean)^2))
(sample_variance <- mean((sample( omega, 100) - sample_mean)^2))
(sample_variance <- mean((sample( omega, 1000) - sample_mean)^2))
(sample_variance <- mean((sample( omega, 10000) - sample_mean)^2))

##############################################################################################
# LLN & CLT (LNp.1-6-97)
##############################################################################################
rand <- function(x) {rgamma(rnorm(rexp(x)), shape=2.5) + rbeta(rnorm(x), shape1=100, shape2=2)}
omega <- rand(100000)
hist(omega)

# Plot the sample mean. A normal distribution is shown.
sample_mean <- NULL
for( i in c(1:1000) ) {
    sample_mean <- c(sample_mean, mean(sample( omega, 5)))
}
hist(sample_mean)

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
# Normal dist and Binomial dist approximation
##############################################################################################
# Observe how large n can result in approximation of normal pdf and binomial pmf. (LNp.6-36)
# Note that approximation via pdf/pmf is in-correct, as pointed out in later example
data <- data.frame( x=c(0:50))
ggplot(data, aes(x=x, y=dbinom(x, size=5, prob=0.7))) + geom_bar(stat="identity") + stat_function(fun=dnorm, args=list(mean=3.5, sd=1.05^(1/2))) + ylim(0, 0.4)
ggplot(data, aes(x=x, y=dbinom(x, size=10, prob=0.7))) + geom_bar(stat="identity") + stat_function(fun=dnorm, args=list(mean=7, sd=2.1^(1/2))) + ylim(0, 0.4)
ggplot(data, aes(x=x, y=dbinom(x, size=20, prob=0.7))) + geom_bar(stat="identity") + stat_function(fun=dnorm, args=list(mean=14, sd=4.2^(1/2))) + ylim(0, 0.4)
ggplot(data, aes(x=x, y=dbinom(x, size=50, prob=0.7))) + geom_bar(stat="identity") + stat_function(fun=dnorm, args=list(mean=35, sd=10.5^(1/2))) + ylim(0, 0.4)


##############################################################################################
# Method of moments LNp.8-8
##############################################################################################
x <- c(31, 29, 19, 18, 31, 28, 34, 27, 34, 30, 16, 18, 26, 27, 27, 18, 24, 22, 28, 24, 21, 17, 24)
lambda <- mean(x)

hist(x)
plot(dpois(c(1:50), lambda=lambda))
factor(x)
data.frame(x=x,
