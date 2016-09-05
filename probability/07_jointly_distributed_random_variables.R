library(ggplot2)
library(plyr)
require(rgl)      # For persp3d
library(distr)    # For DiscreteDistribution
library(numDeriv) # For grad

##############################################################################################
# Discrete joint distribution:
#
#   Using LNp.7-3 as an example, which is a dependent jointly distribution
#   This is an example in data frame
#
##############################################################################################
data <- data.frame( x1=c(0:3,0:3), x2=c(0, 0, 0, 0, 1, 1, 1, 1), joint_prob=c(1/8, 2/8, 1/8, 0, 0, 1/8, 2/8, 1/8))

# calculate joint cdf (Should have a better approach)
cum <- NULL
for(x2 in c(0:1)){
    for(x1 in c(0:3)){
        cum <- rbind(cum, sum(data[data$x1<=x1&data$x2<=x2,]$joint_prob))
    }
}
cum
data$joint_cum_prob <- cum

# alternative approach to calculate cdf
# out_df <- data.frame()
# for (i in c(min(x$x1):max(x$x1))) {
#     for (j in c(min(x$x2):max(x$x2))) {
#         out <- x %>%
#             filter(x1 <= i, x2 <= j) %>%
#             summarise(sum(joint_prob)) %>%
#             as.numeric()
#         out_df <- rbind(out_df, c(i, j, out))
#     }
# }

# plot jointly pmf
ggplot(data, aes(x=x1, y=x2)) + geom_point() + geom_tile(aes(fill=joint_prob))

# plot jointly cdf
ggplot(data, aes(x=x1, y=x2)) + geom_point() + geom_tile(aes(fill=joint_cum_prob))
ggplot(data, aes(x=x1, y=x2)) + geom_point() + geom_contour(aes(z=joint_cum_prob))

# calculate x1 marginal dist by add-up rows
ddply(data, .(x1), summarize, test=sum(joint_prob))
# calculate x2 marginal dist by add-up cols
ddply(data, .(x2), summarize, test=sum(joint_prob))

View(data)

##############################################################################################
# Discrete joint distribution:
#
#   Using LNp.7-3 as an example, which is a dependent jointly distribution
#   This is an example in matrix
#
##############################################################################################
x1 <- c(0, 1, 2, 3)
x2 <- c(0, 1)
joint_prob <- matrix( c(1/8, 0, 2/8, 1/8, 1/8, 2/8, 0, 1/8), nrow=2, ncol=4)

# calculate joint cdf
# (TBD)

# plot jointly pmf
persp3d(x1, x2, joint_prob, col="blue")

# plot jointly cdf
# (TBD)

# calculate x1 marginal dist by add-up rows
colSums (joint_prob, na.rm = FALSE, dims = 1)
# calculate x2 marginal dist by add-up columns
rowSums (joint_prob, na.rm = FALSE, dims = 1)

##############################################################################################
# Independent joint distribution:
#
#   This is a bivariant example: 2 random variables: x and y.
#   Both x and y are exponential distribution.
#   x and y are assumed to be independent, so their joint pdf and cdf are the multiplication
#   of marginal pdf and cdf, respectively. LNp.7-21
#
##############################################################################################
x <- c(0:10)
y <- c(0:10)
rate <- 0.5
joint <- data.frame( x=x, y=y, x_pdf=dexp(x, rate), x_cdf=pexp(x, rate), y_pdf=dexp(y, rate), y_cdf=pexp(y, rate))

# calculate joint pdf (by multiplication, since x, y are assumed to be independent. LNp.7-21)
joint_pdf <- outer(x, y, function(x, y){
    dexp(x, rate) * dexp(y, rate)
})

# calculate joint cdf (by multiplication, since x, y are assumed to be independent. LNp.7-21)
joint_cdf <- outer(x, y, function(x, y){
    pexp(x, rate) * pexp(y, rate)
})

# plot marginal pdf
ggplot(joint, aes(x=x)) + geom_area(stat = "function", fun = dexp, args=list(rate=rate), fill = "black")
ggplot(joint, aes(x=y)) + geom_area(stat = "function", fun = dexp, args=list(rate=rate), fill = "black")

# plot marginal cdf
ggplot(joint, aes(x=x)) + stat_function(fun=pexp, args=list(rate=rate))
ggplot(joint, aes(x=y)) + stat_function(fun=pexp, args=list(rate=rate))

# plot joint pdf
persp3d(x, y, joint_pdf, col="blue")

# Generate random number with such jointly distribution, and visualize with scatter plot
data <- data.frame(x=rexp(5000, rate), y=rexp(5000, rate))
ggplot(data, aes(x=x, y=y)) + geom_point()

# plot joint cdf
persp3d(x, y, joint_cdf, col="blue")

# If two random variables x, y are independent, for any x, the y distribution is identical LNp.7-27
# This can be proven by observing the following: after normalization, the probability for every y
# is the same
joint_pdf[,1] / sum(joint_pdf[,1])
joint_pdf[,2] / sum(joint_pdf[,2])
joint_pdf[,3] / sum(joint_pdf[,3])
joint_pdf[,4] / sum(joint_pdf[,4])
joint_pdf[,5] / sum(joint_pdf[,5])
joint_pdf[,6] / sum(joint_pdf[,6])
joint_pdf[,7] / sum(joint_pdf[,7])
joint_pdf[,8] / sum(joint_pdf[,8])
joint_pdf[,9] / sum(joint_pdf[,9])
joint_pdf[,10] / sum(joint_pdf[,10])
joint_pdf[,11] / sum(joint_pdf[,11])

# However, this is not easy to understand visually, without log(y)
joint_pdf_df <- as.data.frame(joint_pdf)
colnames(joint_pdf_df) <- c("X1", "X2", "X3", "X4", "X5", "X6", "X7", "X8", "X9", "X10", "X11")
ggplot(joint_pdf_df, aes(x=c(1:11), y=X1)) + geom_line() + geom_line(aes(y=X2)) + geom_line(aes(y=X3)) + geom_line(aes(y=X11))
ggplot(joint_pdf_df, aes(x=c(1:11), y=X1)) + geom_line() + geom_line(aes(y=X2)) + geom_line(aes(y=X3)) + geom_line(aes(y=X11)) + scale_y_log10()


# Other bivariant jointly distribution
# I find it interesting to look at
data <- data.frame(x=rnorm(5000), y=rexp(5000, rate))
ggplot(data, aes(x=x, y=y)) + geom_point()

##############################################################################################
# Multinomial distribution:
#
##############################################################################################
# Generate random number based on multinomial distribution
# For a basic experiment with 6 outcomes, each with 1/6 probability.
# For total 10 trials, generate the random outcome based on such multinomial distribution
trial <- 10
die_sum <- 6
die_prob <- c(1/6, 1/6, 1/6, 1/6, 1/6, 1/6)
rmultinom(trial, size = die_sum, prob = die_prob)

# For a basic experiment with 2 outcomes (x1, x2), each with 1/2 probability.
# For total 3 trials, calculate probability for all outcome combinations.
#
# All possible outcomes are listed as follows:
#
#   x1 x2 prob  comment
#    0  3 0.125 In 3 trials, the probability that, x2 occurs 3 times, is 0.125
#    1  2 0.375 In 3 trials, the probability that, x1 occurs 1 time and x2 occurs 2 times, is 0.375
#    2  1 0.375 In 3 trials, the probability that, x1 occurs 2 time and x2 occurs 1 times, is 0.375
#    3  0 0.125 In 3 trials, the probability that, x1 occurs 3 times, is 0.125
data <- data.frame( x1=c(0:3), x2=c(3:0))
data <- ddply(data, .(x1, x2), summarize, joint_prob = dmultinom(c(x1, x2), size = NULL, prob = c(1/2, 1/2)) )
data

# Plot the joint pmf
# The points are only available at diagonal, since x1 + x2 = 3
ggplot(data, aes(x=x1, y=x2)) + geom_point() + geom_tile(aes(fill=joint_prob))

# This is why such bivariant distribution can be discussed on only 1 variable
# LNp.7-16
ggplot(data, aes(x=x1, y=joint_prob)) + geom_point() + geom_line()


# For a basic experiment with 3 outcomes (x1, x2, x3), each with 1/3 probability.
# For total 3 trials, calculate probability for all outcome combinations.
#
# All possible outcomes are listed as follows:
#
# The number of all outcomes is 10, which can be calculated by (n+r-1)!/(n-1)!r!
# The probability of each outcome is a partition problem.
# For example, p(x1=0, x2=1, x3=1) = ( 3! / 0! * 1! * 1! ) / 27 = 6/27 = 0.22222222
#
#   x1 x2 x3 joint_prob number
#   0  0  3  0.03703704  1
#   0  1  2  0.11111111  3
#   0  2  1  0.11111111  3
#   0  3  0  0.03703704  1
#   1  0  2  0.11111111  3
#   1  1  1  0.22222222  6
#   1  2  0  0.11111111  3
#   2  0  1  0.11111111  3
#   2  1  0  0.11111111  3
#   3  0  0  0.03703704  1
# SUM        1          27 = (3 x 3 x 3)
#
# This is in fact corresponds to (x1+x2+x3)^3
data <- data.frame( x1=c(0, 0, 0, 0, 1, 1, 1, 2, 2, 3), x2=c(0, 1, 2, 3, 0, 1, 2, 0, 1, 0), x3=c(3, 2, 1, 0, 2, 1, 0, 1, 0, 0))
data
data <- ddply(data, .(x1, x2, x3), summarize, joint_prob = dmultinom(c(x1, x2, x3), size = NULL, prob = c(1/3, 1/3, 1/3)) )
sum(data$joint_prob)

# Plot the joint pmf
# The points are only available at diagonal, since x1 + x2 = 3
ggplot(data, aes(x=x1, y=x2)) + geom_point() + geom_tile(aes(fill=joint_prob))

# If a balanced (6-slided) die is rolled 12 times, what's the probability that each face
# appears twice?
dmultinom(c(2, 2, 2, 2, 2, 2), size = NULL, prob = c(1/6, 1/6, 1/6, 1/6, 1/6, 1/6))

##############################################################################################
# Joint distribution transformation - method of event
#
#   For 2 independent Poisson distribution X ~ Poisson(l1) and Y ~ Poisson(l2),
#   the tranformation of X+Y=Z will be Z ~ Poisson(l1+l2)
#
##############################################################################################
x <- c(-5:20)
l1 <- 5
x.pois <- data.frame( x=x, pmf=dpois(x, l1), cdf=ppois(x, l1))
x.pois
ggplot(x.pois, aes(x=x, y=pmf)) + geom_bar(stat="identity")

y <- c(-5:20)
l2 <- 10
y.pois <- data.frame( y=y, pmf=dpois(x, l2), cdf=ppois(x, l2))
y.pois
ggplot(y.pois, aes(x=x, y=pmf)) + geom_bar(stat="identity")

x.pois
y.pois

z <- c(-5:20)
l3 <- l1 + l2
z.pois <- data.frame( z=z, pmf=dpois(x, l3), cdf=ppois(x, l3))
z.pois
ggplot(z.pois, aes(x=z, y=pmf)) + geom_bar(stat="identity")

# convolution of Poisson as an example
dpois(0, l1) * dpois(4, l2) +
    dpois(1, l1) * dpois(3, l2) +
    dpois(2, l1) * dpois(2, l2) +
    dpois(3, l1) * dpois(1, l2) +
    dpois(4, l1) * dpois(0, l2)

dpois(4, l3)

##############################################################################################
# Order Statistics
#
# Practice with X(1) and X(n).
# Other jointly distribution is skipped.
##############################################################################################
# Given two random variables following i.i.d
fx1 <- function(x1) {dexp(x1, rate=0.5)}
Fx1 <- function(x1) {pexp(x1, rate=0.5)}
fx2 <- function(x2) {dexp(x2, rate=0.5)}
Fx2 <- function(x2) {pexp(x2, rate=0.5)}

# 1. According to theorem, cdf X(1)=1-[1-F(x)]^n pdf X(1)=nf(x)[1-F(x)]^(n-1)
fmin <- function(x, n=2) {n*dexp(x, rate=0.5) * (1-pexp(x, rate=0.5))^(n-1)}
Fmin <- function(x, n=2) {1-(1-pexp(x, rate=0.5))^n}

# Calculate the probability of min <= 3
# One can use the derived function directly
Fmin(3)
# Or, conduct by reasoning
1-((1-Fx1(3)) * (1 - Fx2(3)))

# 2. According to theorem, cdf X(n)=[F(x)]^n pdf X(n)=nf(x)[F(x)]^(n-1)
fmax <- function(x, n=2) {n*dexp(x, rate=0.5) * (pexp(x, rate=0.5))^(n-1)}
Fmax <- function(x, n=2) {pexp(x, rate=0.5)^n}

# Calculate the probability of max <= 3
# One can use the derived function directly
Fmax(3)
# Or, conduct by reasoning
Fx1(3)*Fx2(3)

# LNp.7-47 example
# n light bulbs are placed in service at time t=0, and allowed to burn continuously. Denote
# their lifetime by X1,...,Xn, and suppose that they are i.i.d. with cdf F. If burned out
# bulbs are not replaced, then the probability that the room is still lighted after two
# month is?
fmax <- function(x, n=5) {n*dexp(x, rate=1) * (pexp(x, rate=1))^(n-1)}
Fmax <- function(x, n=5) {pexp(x, rate=1)^n}
1-Fmax(2)

##############################################################################################
# Conditional distribution LNp.7-55
##############################################################################################
# The joint pdf of x, y is as follows, given ( 0 <= x,y < Inf )
fxy <- function(x, y) {
    sapply(x, function(z, y) {
        2/((1+z+y)^3)
        }, y=y)
}

# The marginal pdf of x is obtained from integrate y part of fxy
# One can verify if fx equals to 1/(1+x)^2
fx <- function(x) {
    sapply(x, function(y) integrate(fxy, lower=0, upper=Inf, y=y)$val)
}

fx(10)
1/(1+10)^2

# Integrate fx again, to prove the total prob of fxy is 1.
integrate(fx, 0, Inf)

# The conditional prob fy|x is fxy/fx
# One can verify if this equals to 2(1+x)^2 / (1+x+y)^3
fy_x <- function(x,y) {
    fxy(x, y) / fx(x)
}

fy_x(1,1)
2*(1+1)^2 / (1+1+1)^3
