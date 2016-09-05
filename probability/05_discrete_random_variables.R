library(ggplot2)
library(plyr)
require(rgl)      # For persp3d
library(distr)    # For DiscreteDistribution
library(numDeriv) # For grad

##############################################################################################
# Univariant discrete random variable (LNp.5-3~LNp.5-11)
#
# DiscreteDistribution is used to create pmf/cdf/... from given random variable results.
##############################################################################################
data <- data.frame(trial1=c("h","h","h", "t", "h", "t", "t", "t"),
                   trial2=c("h","h","t","h","t","h","t","t"),
                   trial3=c("h","t","h","h","t","t","h","t"))

# Generate the outcome of the three random variables x1, x2, x3
data <- ddply(data, .(trial1, trial2, trial3), transform,
              x1=(trial1=="h")+(trial2=="h")+(trial3=="h"),
              x2=sum(trial1=="h"),
              x3=(trial1=="h")+(trial2=="h")+(trial3=="h")-(trial1=="t")-(trial2=="t")-(trial3=="t"))

# Calculate the probability and cumulative probability for x1, x2, x3
x1 <- ddply(data, .(x1), summarize, freq=length(x1), prob=freq/nrow(data))
x2 <- ddply(data, .(x2), summarize, freq=length(x2), prob=freq/nrow(data))
x3 <- ddply(data, .(x3), summarize, freq=length(x3), prob=freq/nrow(data))
x1$acc_prob <- cumsum(x1$prob)
x2$acc_prob <- cumsum(x2$prob)
x3$acc_prob <- cumsum(x3$prob)

# Create pmf and cdf for x1 (x2, x3 are skipped)
x1.dist <- DiscreteDistribution (supp = x1$x1, prob = x1$prob)
x1.pmf <- d(x1.dist)  ## Probability mass function
x1.cdf <- p(x1.dist)  ## Cumulative distribution function

# One can actually use call such function
x1.pmf(2)
x1.pmf(3)
x1.pmf(2.5) # No probability exists at 2.5
x1.cdf(2.5) # probability exists at 2.5 due to cdf natural

# cdf jumps at points with weight
x1.pmf(0) == x1.cdf(0)
x1.pmf(0) == x1.cdf(0.1)
x1.pmf(0) == x1.cdf(0.9)
x1.pmf(0) == x1.cdf(1)
x1.pmf(0) + x1.pmf(1) == x1.cdf(1)
# 1 = pmf(0) + pmf(1) + pmf(2) + pmf(3)
# cdf(2) = pmf(0) + pmf(1) + pmf(2)
x1.pmf(3) == 1 - x1.cdf(2)

# Plot pmf/cdf
plot(x1.dist)

# FIXME: one way to plot the pmf/cdf is using ggplot stat_function, with the pmf/cdf created
# above. However, the plotting is wrong. I think this is a bug.
ggplot(x1, aes(x=x1)) + stat_function(fun=x1.pmf)
ggplot(x1, aes(x=x1)) + stat_function(fun=x1.cdf)

##############################################################################################
# Univariant discrete random variable transformation (LNp.5-12)
##############################################################################################
data <- data.frame(trial1=c("h","h","h", "t", "h", "t", "t", "t"),
                   trial2=c("h","h","t","h","t","h","t","t"),
                   trial3=c("h","t","h","h","t","t","h","t"))

# Generate the outcome of the three random variables x1
data <- ddply(data, .(trial1, trial2, trial3), transform,
              x1=(trial1=="h")+(trial2=="h")+(trial3=="h"))

# Calculate the probability and cumulative probability for x1
x1 <- ddply(data, .(x1), summarize, freq=length(x1), prob=freq/nrow(data))
x1$acc_prob <- cumsum(x1$prob)
x1.dist <- DiscreteDistribution (supp = x1$x1, prob = x1$prob)
x1.pmf <- d(x1.dist)
plot(x1.dist)

# Assume a transformation Y1=g(X1): Y1=X1^2
# Approach 1: map Y1 to the original sample space to get its pmf fy(Y=y)
data$y1 <- data$x1^2
y1 <- ddply(data, .(y1), summarize, freq=length(y1), prob=freq/nrow(data))
y1$acc_prob <- cumsum(y1$prob)
y1.dist <- DiscreteDistribution (supp = y1$y1, prob = y1$prob)
y1.pmf <- d(y1.dist)
plot(y1.dist)
# Approach 2: replace fx(X=x) with fx(X=g'(Y1))
x1.pmf(sqrt(4)) + x1.pmf(-sqrt(4))

y1.pmf(4)

##############################################################################################
# Univariant discrete random variable expectation (LNp.5-17)
##############################################################################################
x <- data.frame(x=c(0,1,2,3,4), fx=c(5/210, 50/210, 100/210, 50/210, 5/210))


# calculate mean: x * f(x)
(mean = sum(x$x * x$fx))

# calculate variance: (x-u)^2 * f(x)
(var = sum( (x$x - mean)^2 * x$fx ))

# calculate moment: x^n * f(x)
# here, n is set to 2: x^2 * f(x)
n <- 2
(mean_of_sqrt = sum((x$x)^n * x$fx))

# variance can be calculated alternative from mean
# var = E(x^2) - (E(x))^2
var
mean_of_sqrt - (mean^2)

# Calculate MSE (Mean Squre Error), given c=3
c <- 3
# calculate by creating another random variable t
x$t <- (x$x - c)^2
sum(x$t * x$fx)
# calculate by equation
bias <- c-mean
var + bias^2

##############################################################################################
#
# Common distribution model:
#
#   Binomial:          binom
#   Negative Binomial: nbinom
#   Geometric:         geom
#   Poisson:           pois
#   Hypergeometric:    hyper
#   Uniform:           unif
#   Exponential:       exp
#   Gamma:             gamma
#   Beta:              beta
#   Normal:            norm
#   Weibull:           weibull
#   Cauchy:            cauchy
#   Chi-square:        chisq (only listed here)
#   t:                 t     (only listed here)
#   F:                 F     (only listed here)
##############################################################################################


##############################################################################################
# Binomial distribution:
#
#   每次嘗試的成功機率為prob，在size次嘗試中，成功x次的機率
#
#   pmf:      dbinom(x, size, prob, log = FALSE)
#   cdf:      pbinom(q, size, prob, lower.tail = TRUE, log.p = FALSE)
#   quantile: qbinom(p, size, prob, lower.tail = TRUE, log.p = FALSE)
#   random:   rbinom(n, size, prob)
#
##############################################################################################
x <- c(1:9)
size <- 9
prob <- 0.3038

dbinom(0, 9, 0.3038) + dbinom(1, 9, 0.3038)
pbinom(1, 9, 0.3038)
qbinom(0.18, 9, 0.3038)
hist(rbinom(10000, 2, 0.5))

binomi <- data.frame( x=x, pmf=dbinom(x, size, prob), cdf=pbinom(x, size, prob))
binomi.summary <- ddply(binomi, .(), summarize, mean=sum(x * pmf), variance=sum(pmf * ((x-mean)^2)))

binomi
binomi.summary

ggplot(binomi, aes(x=x, y=pmf)) + geom_bar(stat="identity")
ggplot(binomi, aes(x=x, y=cdf)) + geom_bar(stat="identity")

# What is the probability that South gets no Aces on at least k=5 of n=9 hands?
prob <- 0.3038
sum(dbinom(c(5:9), 9, prob))

# Given probability = 0.9, how many trails expected to have Aces at n=9 hands?
qbinom(0.9, 9, prob)

##############################################################################################
# Negative Binomial distribution:
#
#   每次成功機率為prob，要求size次成功，總共需要失敗x次的機率
#
#   pmf:      dnbinom(x, size, prob, mu, log = FALSE)
#   cdf:      pnbinom(q, size, prob, mu, lower.tail = TRUE, log.p = FALSE)
#   quantile: qnbinom(p, size, prob, mu, lower.tail = TRUE, log.p = FALSE)
#   random:   rnbinom(n, size, prob, mu)
#
##############################################################################################
x <- c(1:25)
size <- 3
prob <- 1/3

nbinom <- data.frame( x=x, pmf=dnbinom(x, size, prob), cdf=pnbinom(x, size, prob))
nbinom.summary <- ddply(nbinom, .(), summarize, mean=sum(x * pmf), variance=sum(pmf * ((x-mean)^2)))

nbinom
nbinom.summary

ggplot(nbinom, aes(x=x, y=pmf)) + geom_bar(stat="identity")
ggplot(nbinom, aes(x=x, y=cdf)) + geom_bar(stat="identity")

# Equivalence of Binomial and Negative binomial (LNp.5-27)
# Negative binomial: Each time hiring probability is 1/3, the probability of exact 10 trials
#                    required to hire 3 engineers. (10th trial is the last one for the 3rd
#                    engineer.)
# Binomial:          Each time hiring probability is 1/3, the probability of hiring 2 engineers
#                    in exact 9 trials. Then, times by 1/3.
dnbinom(7, 3, 1/3)
dbinom(2, 9, 1/3) * (1/3)
##############################################################################################
# Geometric distribution:
#
#   每次成功機率為prob，要求1次成功，總共需要失敗x次的機率
#
#   pmf:      dgeom(x, prob, log = FALSE)
#   cdf:      pgeom(q, prob, lower.tail = TRUE, log.p = FALSE)
#   quantile: qgeom(p, prob, lower.tail = TRUE, log.p = FALSE)
#   random:   rgeom(n, prob)
#
##############################################################################################
x <- c(1:25)
prob <- 0.2

geom <- data.frame( x=x, pmf=dgeom(x, prob), cdf=pgeom(x, prob))
geom.summary <- ddply(geom, .(), summarize, mean=sum(x * pmf), variance=sum(pmf * ((x-mean)^2)))

geom
geom.summary

ggplot(geom, aes(x=x, y=pmf)) + geom_bar(stat="identity")
ggplot(geom, aes(x=x, y=cdf)) + geom_bar(stat="identity")

# Equivalence of Negative binomial and Geometric (LNp.5-29)
# Geometric is a special case of Negative binomial when success number is 1 (size=1)
dnbinom(9, 1, 1/3)
dgeom(9, 1/3)
##############################################################################################
# Poisson distribution:
#
#   每次嘗試的成功機率為prob，在size次嘗試中，成功x次的機率
#   size很大，size遠大於x，prob很小
#
#   pmf:      dpois(x, lambda, log = FALSE)
#   cdf:      ppois(q, lambda, lower.tail = TRUE, log.p = FALSE)
#   quantile: qpois(p, lambda, lower.tail = TRUE, log.p = FALSE)
#   random:   rpois(n, lambda)
#
##############################################################################################
x <- c(1:20)
size <- 2500
prob <- 0.001
lambda <- size * prob

pois <- data.frame( x=x, pmf=dpois(x, lambda), cdf=ppois(x, lambda))
pois.summary <- ddply(pois, .(), summarize, mean=sum(x * pmf), variance=sum(pmf * ((x-mean)^2)))

pois
pois.summary

ggplot(pois, aes(x=x, y=pmf)) + geom_bar(stat="identity")
ggplot(pois, aes(x=x, y=cdf)) + geom_bar(stat="identity")

# A professor hits the wrong key with probability p=0.001 each time he types a letter. Assume
# independence for the occurrence of errors between different letter typings. What's the
# probability that 5 or more errors in n=2500 letters. (LNp.5-33)
1-sum( dpois(c(0:4), 2500 * 0.001) )

# Traffic accident occurs according to a Poisson process at a rate of 5.5 per month.
#  What is the probability of 3 or more accidents occur in a 2 month periods?
1-sum(dpois(c(0:2), lambda = 5.5 * 2))

##############################################################################################
# Hypergeometric distribution:
#
#   盒子中共有m個白球、n個黑球，從盒子中要抽k球(抽球後不放回)，其中x個白球的機率
#
#   pmf:      dhyper(x, m, n, k, log = FALSE)
#   cdf:      phyper(q, m, n, k, lower.tail = TRUE, log.p = FALSE)
#   quantile: qhyper(p, m, n, k, lower.tail = TRUE, log.p = FALSE)
#   random:   rhyper(nn, m, n, k)
#
##############################################################################################
x <- c(1:10)
m <- 10
n <- 7
k <- 8

hyper <- data.frame( x=x, pmf=dhyper(x, m, n, k), cdf=phyper(x, m, n, k, lambda))
hyper.summary <- ddply(hyper, .(), summarize, mean=sum(x * pmf), variance=sum(pmf * ((x-mean)^2)))

hyper
hyper.summary

ggplot(hyper, aes(x=x, y=pmf)) + geom_bar(stat="identity")
ggplot(hyper, aes(x=x, y=cdf)) + geom_bar(stat="identity")

