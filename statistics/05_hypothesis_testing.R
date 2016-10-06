library(ggplot2)
##############################################################################################
# LNp.9-8
#
# Two coins:
#   Coin 0: Binomial(10, 0.5)
#   Coin 1: Binomial(10, 0.7)
#
# Omega = { Binomial(10, 0.5), Binomial(10, 0.7) }
#
# Hypothesis:
#   H0: p=0.5
#   HA: p=0.7
##############################################################################################
x <- c(0:10)
x <- rep(c(0:10), 2)
data <- data.frame(x=x, theta=c(rep(0.5, 11), rep(0.7, 11)))
data$prob <- dbinom(x, size=10, prob=data$theta)
data$type <- ""
head(data)

# The test statistics is defined as the sum of X1~X10
# The critical value is set to 6.5
# Rejection Region (RR): test statistics > 6.5
# Accept Region (AR)   : test statistics < 6.5
sample_space <- c(0:10)
critical_value <- 6.5
reject_region <- subset(sample_space, sample_space > critical_value)
accept_region <- setdiff(sample_space, reject_region)

# The alpha for Binomial(10, 0.5) is set to 0.171875
(alpha <- sum(dbinom(reject_region, size=10, prob=0.5)))
# The beta for Binomial(10, 0.7) is set to 0.3503893
(beta <- sum(dbinom(accept_region, size=10, prob=0.7)))
# The power
(power <- 1-beta)

# Divide the H0/HA according to the test statistics
data[data$theta==0.5 & data$x %in% reject_region,]$type <- "type-I (alpha)"
data[data$theta==0.5 & data$x %in% accept_region,]$type <- "H0"
data[data$theta==0.7 & data$x %in% reject_region,]$type <- "H1 (power)"
data[data$theta==0.7 & data$x %in% accept_region,]$type <- "type-II (beta)"

# Plot H0/HA and RR/AR
ggplot(data, aes(x=x)) + geom_bar(aes(y=prob, fill=type), stat="identity") + facet_grid(theta ~.)

# Plot power function
# {TBD: the code is cumbersome.}
prob <- seq(0.5, 1.0, 0.01)
power_vect <- NULL
for( p in prob ) {
    power <- 1 - sum(dbinom(accept_region, size=10, prob=p))
    power_vect <- c(power_vect, power)
}
power_vect
ggplot(data=data.frame(prob=prob, power=power_vect), aes(x=prob, y=power)) + geom_line()

# LNp.9-13
# If we observe 10 successes, then the p-value is 0.0010
1-pbinom(9, 10, 0.5)
# If we observe 9 successes, then the p-value is 0.0107
1-pbinom(8, 10, 0.5)
# If we observe 8 successes, then the p-value is 0.0547
1-pbinom(7, 10, 0.5)

##############################################################################################
# Likelihood ratio test (LNp.9-16)
#
# Two coins:
#   Coin 0: Binomial(10, 0.5)
#   Coin 1: Binomial(10, 0.7)
#
# Omega = { Binomial(10, 0.5), Binomial(10, 0.7) }
#
# Hypothesis:
#   H0: p=0.5
#   HA: p=0.7
#
#   likelihood ratio: H0/HA < c
#
##############################################################################################
data <- data.frame(x=c(0:10), f0=dbinom(x, 10, 0.5), fa=dbinom(x, 10, 0.7))
data$likeratio <- data$f0/data$fa
data$likeratio_sub <- data$f0 - data$fa

ggplot(data, aes(x=x, y=likeratio)) + geom_line()
ggplot(data, aes(x=x, y=likeratio_sub)) + geom_line()

##############################################################################################
# Generalized likelihood ratio test (LNp.9-31)
#
# Omega = { Binomial(10, [0, 1]) }
#
# Hypothesis:
#   H0: p=[0, 0.5]
#   HA: p=[0.5, 1]
##############################################################################################
max_h0 <- function(x) {
    max(dbinom(x, size=10, prob=seq(0, 0.5, 0.01)))
}
max_all <- function(x) {
    max(dbinom(x, size=10, prob=seq(0, 1, 0.01)))
}
glr <- function(x) {
    max_h0(x)/max_all(x)
}

max_h0(2)
max_all(2)
glr(2)

# Plot GLR
# {TBD: the code is cumbersome.}
x_vect <- c(0:10)
glr_vect <- NULL
for( x in x_vect ) {
    glr_vect <- c(glr_vect, glr(x))
}
glr_vect
ggplot(data=data.frame(x=x_vect, y=glr_vect), aes(x=x, y=glr_vect)) + geom_line()

##############################################################################################
# Generalized likelihood ratio test (LNp.9-31)
#
# Omega = { Normal(mean=[-10, 20], sd=5) }
#
# Hypothesis:
#   H0: mean=[0]
#   HA: mean=[-10, 20]-[0]
##############################################################################################
max_h0 <- function(x) {
    max(dnorm(x, mean=0, sd=5))
}
max_all <- function(x) {
    max(dnorm(x, mean=seq(-10, 20, 0.01), sd=5))
}
glr <- function(x) {
    max_h0(x)/max_all(x)
}

max_h0(2)
max_all(2)
glr(2)

# Plot GLR
# {TBD: the code is cumbersome.}
x_vect <- c(-20:20)
glr_vect <- NULL
for( x in x_vect ) {
    glr_vect <- c(glr_vect, glr(x))
}
glr_vect
ggplot(data=data.frame(x=x_vect, y=glr_vect), aes(x=x, y=glr_vect)) + geom_line()

##############################################################################################
# Hardy-Weihberg model (LNp.9-43)
#
# Hypothesis:
#   H0: p1, p2, p3 are specified by the Hardy-Weihberg model
#   HA: p1, p2, p3 do not have the specified form
##############################################################################################
obs <- c(342, 500, 187)
exp <- c(340.6, 502.8, 185.6)

# Set alpha=0.05
# Reject H0 if the value of chisq statistic exceeds 3.84 (95%-quantile)
qchisq(0.95, df=1)

# Perform chisq test
results <- chisq.test(obs, p=exp, rescale.p=TRUE, simulate.p.value=TRUE)

# The statistics (0.0319) is smaller than 3.84 (95%-quantile), thus accept H0
results$statistic

# p-Value is 0.985, which is very hard to reject H0
results$p.value

##############################################################################################
# Fisher's reexamination of Mendel's data (LNp.9-45)
#
# Hypothesis:
#   H0: p1, p2, p3 are specified by the Hardy-Weihberg model
#   HA: p1, p2, p3 do not have the specified form
##############################################################################################
obs <- c(315, 108, 102, 31)
exp <- c(312.75, 104.25, 104.25, 34.75)

# Set alpha=0.05
# Reject H0 if the value of chisq statistic exceeds 7.81 (95%-quantile)
qchisq(0.95, df=3)

# Perform chisq test
results <- chisq.test(obs, p=exp, rescale.p=TRUE, simulate.p.value=TRUE)

# The Peason's statistics (0.604) is smaller than 7.81 (95%-quantile), thus accept H0
results$statistic

# p-Value is 0.903, which is very hard to reject H0
results$p.value
