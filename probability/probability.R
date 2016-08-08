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
dbinom(2, 9, 1/3) / 3
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

##############################################################################################
# Continuous random variable distribution: pdf and cdf
#
#  cdf = pdf的積分
#  pdf = cdf的微分
#
#   pdf:      dnorm(x, mean = 0, sd = 1, log = FALSE)
#   cdf:      pnorm(q, mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE)
#   quantile: qnorm(p, mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE)
#   random:   rnorm(n, mean = 0, sd = 1)
#
##############################################################################################
x <- c(-5:5)
mean <- 0
sd <- 1

norm <- data.frame( x=x, pdf=dnorm(x, mean, sd), cdf=pnorm(x, mean, sd))
norm

ggplot(norm, aes(x=x)) + stat_function(fun=dnorm, args=list(mean=mean, sd=sd)) + geom_point(aes(x=2.5, y=0.3))
ggplot(norm, aes(x=x)) + stat_function(fun=pnorm, args=list(mean=mean, sd=sd))

# pdf由-Inf積分到2.5，等於cdf(2.5)
integrand <- function(x){ dnorm(x, mean=0, sd=1) }
integrate(integrand, lower=-Inf, upper=2.5)
pnorm(2.5, mean=0, sd=1)

# cdf在2.5的微分，等於pdf(2.5)
f <- function(x) {pnorm(x, mean=0, sd=1)}
grad(f, 2.5)
dnorm(2.5, mean=0, sd=1)

##############################################################################################
# Continuous random variable transformation
#
# Method of cdf:
#  Fy(y) = Fx(g-1(y))     If g is strictly increasing
#  Fy(y) = 1 - Fx(g-1(y)) If g is strictly decreasing
#
# Method of pdf:
#  fy(y) = fx(g-1(y)) * |dg-1(y) / dy|
##############################################################################################
# Assume X is a random variable following normal distribution
# X's pdf/cdf is fx/Fx, respectively
fx <- function(x) {dnorm(x)}
Fx <- function(x) {pnorm(x)}

# Generate X's random number, for verification
data <- data.frame( x=rnorm(1000))

# Take a look at X's pdf and cdf
# blue part is the stat of pdf/cdf
ggplot(data, aes(x=x)) + geom_histogram(aes(y = ..density..)) + stat_function(fun=fx, colour="blue")
tmp_x <- ddply(data, .(x), summarize, freq=length(x), prob=freq/nrow(data))
tmp_x$acc_prob <- cumsum(tmp_x$prob)
ggplot(data, aes(x=x)) + geom_point(data=tmp_x, aes(x=x, y=acc_prob)) + stat_function(fun=Fx, colour="blue")

# Assume Y is a transformation of X: Y=g(X)=3X+2   g-1(y)=(Y-2)/3
# According to LNp.6-8~6-10, Y's pdf/cdf is as follows
g <- function(x) {3*x+2}
g_inv <- function(y) {(y-2)/3}
fy <- function(y) {fx(g_inv(y))/3}
Fy <- function(y) {Fx(g_inv(y))}
plot(g) # g is a strictly monotone function (increasing)

# Generate Y data according to g(X)
data$y <- g(data$x)

# Bi-plot X and Y pdf
# The red part (Y) is a direct transformation result of g(X), as the comparison with the fy
# The blue part is a mathematical deduction result of fy
# The red part and blue part should be similar
fy_2 <- function(x) {grad(Fy, x)}
ggplot(data, aes(x=x)) + geom_histogram(aes(y = ..density..)) + stat_function(fun=fx) + geom_histogram(aes(x=y, y = ..density..), colour="red") + geom_density(aes(x=y, y = ..density..), colour="red") + stat_function(fun=fy, colour="blue")
# The Y pdf can be also obtained from differential of Y cdf
# The purple part should be identical to previous plot's blue part
ggplot(data, aes(x=x)) + geom_histogram(aes(y = ..density..)) + stat_function(fun=fx) + geom_histogram(aes(x=y, y = ..density..), colour="red") + geom_density(aes(x=y, y = ..density..), colour="red") + stat_function(fun=fy_2, colour="purple")


# Bi-plot X and Y cdf
# The red part is Y's cumulative density based on Y's raw data
# The blue part is a mathematical deduction result of fy
# The red part and blue part should be similar
tmp_y <- ddply(data, .(y), summarize, freq=length(y), prob=freq/nrow(data))
tmp_y$acc_prob <- cumsum(tmp_y$prob)
ggplot(data, aes(x=x)) + geom_point(data=tmp_x, aes(x=x, y=acc_prob)) + stat_function(fun=Fx) + geom_point(data=tmp_y, aes(x=y, y=acc_prob), color="red") + stat_function(fun=Fy, colour="blue")

##############################################################################################
# Continuous random variable distribution: transformation
#
# 1. Let the computer generates random number following normal distribution (data$x).
#    Transform the data with this normal distribution's cdf will generate random number
#    following uniform distribution.
#
#    U ~ uniform(0, 1)
#    X's cdf is F(X)
#    X = F-1(U)
#
# 2. Let the computer generates random number between 0~1 following uniform distribution
#    (data$x). Transform the data with the inverse cdf of a distribution can generate random
#    number of that distribution.
#
##############################################################################################
# Generate random number following normal distribution
data <- data.frame(x = rexp(1000))
hist(data$x)

# Transform the random number with cdf of normal distribution, to generate random number
# following uniform distribution
data$unif_rand <- pexp(data$x)
hist(data$unif_rand)

# Generate random number following uniform distribution
data <- data.frame(x = runif(1000))
hist(data$x)

# Transform the random number with inverse cdf of normal distribution, to generate random
# number following normal distribution
data$norm_rand <- qnorm(data$x)
hist(data$norm_rand)

# Transform the random number with inverse cdf of exponential distribution to generate random
# number following exponential distribution
data$exp_rand <- qexp(data$x)
hist(data$exp_rand)

# Alternatively, one can try all kinds of transformation one by one
data <- data.frame(x = rexp(1000))
hist(data$x)

data$y <- qnorm(pexp(data$x)) # cdf -> cdf-1
hist(data$y)

data$x_log <- log(data$x)
hist(data$x_log)
data$x_sqrt <- data$x ^ (0.3)
hist(data$x_sqrt)

##############################################################################################
# Uniform distribution:
#
#   介於min到max之間的值的機率是一樣的
#
#   pdf:      dunif(x, min = 0, max = 1, log = FALSE)
#   cdf:      punif(q, min = 0, max = 1, lower.tail = TRUE, log.p = FALSE)
#   quantile: qunif(p, min = 0, max = 1, lower.tail = TRUE, log.p = FALSE)
#   random:   runif(n, min = 0, max = 1)
#
##############################################################################################
x <- c(1,2,3,4,5,6,7,7.5,8,9,10)
min <- 3
max <- 7

unif <- data.frame( x=x, pdf=dunif(x, min, max), cdf=punif(x, min, max))
unif

x_fx <- function(x, min=3, max=7){ x*dunif(x, min=min, max=max) }
(mean = integrate(x_fx, lower=-Inf, upper=Inf))
mean

bias_fx <- function(x, min=3, max=7, mean=5){ (x-mean)^2*dunif(x, min=min, max=max) }
(var = integrate(bias_fx, lower=-Inf, upper=Inf))

ggplot(unif, aes(x=x)) + stat_function(fun=dunif, args=list(min=min, max=max))
ggplot(unif, aes(x=x)) + geom_area(stat = "function", fun = dunif, args=list(min=min, max=max), fill = "black") + geom_point(aes(x=5, y=0), colour="blue", size=4)
ggplot(unif, aes(x=x)) + stat_function(fun=punif, args=list(min=min, max=max))

##############################################################################################
# Exponential distribution:
#
#   單位時間內成功次數為rate，成功第1次時，總共需要多少時間的機率值
#
#   pdf:      dexp(x, rate = 1, log = FALSE)
#   cdf:      pexp(q, rate = 1, lower.tail = TRUE, log.p = FALSE)
#   quantile: qexp(p, rate = 1, lower.tail = TRUE, log.p = FALSE)
#   random:   rexp(n, rate = 1)
#
##############################################################################################
x <- c(0:10)
rate <- 0.5

exp <- data.frame( x=x, pdf=dexp(x, rate), cdf=pexp(x, rate))
exp

x_fx <- function(x, rate=2){ x*dexp(x, rate=rate) }
(mean = integrate(x_fx, lower=-Inf, upper=Inf))
mean

bias_fx <- function(x, rate=2, mean=0.5){ (x-mean)^2*dexp(x, rate=rate) }
(var = integrate(bias_fx, lower=-Inf, upper=Inf))


ggplot(exp, aes(x=x)) + stat_function(fun=dexp, args=list(rate=rate))
ggplot(exp, aes(x=x)) + geom_area(stat = "function", fun = dexp, args=list(rate=rate), fill = "black") + geom_point(aes(x=0.5, y=0), colour="blue", size=4)
ggplot(exp, aes(x=x)) + stat_function(fun=pexp, args=list(rate=rate))

# LNp6-17, the lambda(rate) discussion
ggplot(exp, aes(x=x)) + stat_function(fun=dexp, args=list(rate=0.5)) + stat_function(fun=dexp, args=list(rate=1), colour="blue") + stat_function(fun=dexp, args=list(rate=2), colour="red")

##############################################################################################
# Gamma distribution:
#
#   單位時間內成功次數為rate，要成功第shape次，總共需要多少時間的機率值
#
#   pdf:      dgamma(x, shape, rate = 1, scale = 1/rate, log = FALSE)
#   cdf:      pgamma(q, shape, rate = 1, scale = 1/rate, lower.tail = TRUE, log.p = FALSE)
#   quantile: qgamma(p, shape, rate = 1, scale = 1/rate, lower.tail = TRUE, log.p = FALSE)
#   random:   rgamma(n, shape, rate = 1, scale = 1/rate)
#
##############################################################################################
x <- c(0:25)
shape <- 2
rate <- 0.5

gamma <- data.frame( x=x, pdf=dgamma(x, shape, rate), cdf=pgamma(x, shape, rate))
gamma

x_fx <- function(x, shape=2, rate=0.5){ x*dgamma(x, shape=shape, rate=rate) }
(mean = integrate(x_fx, lower=-Inf, upper=Inf))
mean

bias_fx <- function(x, shape=2, rate=0.5, mean=4){ (x-mean)^2*dgamma(x, shape=shape, rate=rate) }
(var = integrate(bias_fx, lower=-Inf, upper=Inf))

ggplot(gamma, aes(x=x)) + stat_function(fun=dgamma, args=list(shape=shape, rate=rate))
ggplot(gamma, aes(x=x)) + geom_area(stat = "function", fun = dgamma, args=list(shape=shape, rate=rate), fill = "black")
ggplot(gamma, aes(x=x)) + stat_function(fun=pgamma, args=list(shape=shape, rate=rate))

# LNp6-26, the alpha/lambda(rate) discussion
ggplot(gamma, aes(x=x)) + stat_function(fun=dgamma, args=list(shape=1, rate=1)) + stat_function(fun=dgamma, args=list(shape=2, rate=1), colour="blue") + stat_function(fun=dgamma, args=list(shape=4, rate=1), colour="red")
ggplot(gamma, aes(x=x)) + stat_function(fun=dgamma, args=list(shape=2, rate=2)) + stat_function(fun=dgamma, args=list(shape=2, rate=1), colour="blue") + stat_function(fun=dgamma, args=list(shape=2, rate=0.5), colour="red")

##############################################################################################
# Beta distribution:
#
#   From mathmetical equation beta function
#
#   pdf:      dbeta(x, shape1, shape2, ncp = 0, log = FALSE)
#   cdf:      pbeta(q, shape1, shape2, ncp = 0, lower.tail = TRUE, log.p = FALSE)
#   quantile: qbeta(p, shape1, shape2, ncp = 0, lower.tail = TRUE, log.p = FALSE)
#   random:   rbeta(n, shape1, shape2, ncp = 0)
#
##############################################################################################
x <- c(0:1)
shape1 <- 10
shape2 <- 10

beta <- data.frame( x=x, pdf=dbeta(x, shape1, shape2), cdf=pbeta(x, shape1, shape2))
beta

x_fx <- function(x, shape1=10, shape2=10){ x*dbeta(x, shape1=shape1, shape2=shape2) }
(mean = integrate(x_fx, lower=-Inf, upper=Inf))
mean

bias_fx <- function(x, shape1=10, shape2=10, mean=0.5){ (x-mean)^2*dbeta(x, shape1=shape1, shape2=shape2) }
(var = integrate(bias_fx, lower=-Inf, upper=Inf))

ggplot(beta, aes(x=x)) + stat_function(fun=dbeta, args=list(shape1=shape1, shape2=shape2))
ggplot(beta, aes(x=x)) + geom_area(stat = "function", fun = dbeta, args=list(shape1=shape1, shape2=shape2), fill = "black")
ggplot(beta, aes(x=x)) + stat_function(fun=pbeta, args=list(shape1=shape1, shape2=shape2))

# LNp6-29, the alpha/beta discussion
ggplot(beta, aes(x=x)) + stat_function(fun=dbeta, args=list(shape1=10, shape2=10)) + stat_function(fun=dbeta, args=list(shape1=3, shape2=3), colour="blue") + stat_function(fun=dbeta, args=list(shape1=1, shape2=1), colour="red")
ggplot(beta, aes(x=x)) + stat_function(fun=dbeta, args=list(shape1=6, shape2=6*19)) + stat_function(fun=dbeta, args=list(shape1=2/3, shape2=2/3*19), colour="blue")
ggplot(beta, aes(x=x)) + stat_function(fun=dbeta, args=list(shape1=6*19, shape2=6)) + stat_function(fun=dbeta, args=list(shape1=2/3*19, shape2=2/3), colour="blue")

##############################################################################################
# Normal distribution:
#
#   ?
#
#   pdf:      dnorm(x, mean = 0, sd = 1, log = FALSE)
#   cdf:      pnorm(q, mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE)
#   quantile: qnorm(p, mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE)
#   random:   rnorm(n, mean = 0, sd = 1)
#
##############################################################################################
x <- c(-10:20)
mean <- 0
sd <- 1

norm <- data.frame( x=x, pdf=dnorm(x, mean, sd), cdf=pnorm(x, mean, sd))
norm

x_fx <- function(x, mean=0, sd=1){ x*dnorm(x, mean=mean, sd=sd) }
(mean = integrate(x_fx, lower=-Inf, upper=Inf))
mean

bias_fx <- function(x, mean=0, sd=1){ (x-mean)^2*dnorm(x, mean=mean, sd=sd) }
(var = integrate(bias_fx, lower=-Inf, upper=Inf))

# The mean of a distribution can be also obtained by cdf (LNp.6-15)
f <- function(x, mean=0, sd=1){ 1-pnorm(x, mean=mean, sd=sd) }
area_i <- integrate(f, lower=0, upper=Inf)
f <- function(x, mean=0, sd=1){ pnorm(x, mean=mean, sd=sd) }
area_ii <- integrate(f, lower=Inf, upper=0)
(mean <- area_i$value - area_ii$value)

ggplot(norm, aes(x=x)) + stat_function(fun=dnorm, args=list(mean=mean, sd=sd))
ggplot(norm, aes(x=x)) + geom_area(stat = "function", fun = dnorm, args=list(mean=mean, sd=sd), fill = "black")
ggplot(norm, aes(x=x)) + stat_function(fun=pnorm, args=list(mean=0, sd=1))

# LNp6-31, the bell-shaped discussion
ggplot(norm, aes(x=x)) + stat_function(fun=dnorm, args=list(mean=0, sd=0.5)) + stat_function(fun=dnorm, args=list(mean=0, sd=1), colour="blue") + stat_function(fun=dnorm, args=list(mean=0, sd=2), colour="red")

# LNp6-33, transformation
mean <- 5
sd <- 0.5
a=2
b=3
ggplot(norm, aes(x=x)) + stat_function(fun=dnorm, args=list(mean=mean, sd=sd)) +
    stat_function(fun=dnorm, args=list(mean=a*mean, sd=(a^2)*sd), colour="red") +   # Y=aX
    stat_function(fun=dnorm, args=list(mean=mean+b, sd=sd), colour="blue") +         # Y=X+b
    stat_function(fun=dnorm, args=list(mean=a*mean+b, sd=(a^2)*sd), colour="green")   # Y=aX+b

# LNp6-33, standard normal distribution
ggplot(, aes(x=c(-5:5))) + stat_function(fun=dnorm, args=list(mean=0, sd=1))

##############################################################################################
# Normal dist and Binomial dist approximation
#
##############################################################################################
# Observe how large n can result in approximation of normal pdf and binomial pmf. (LNp.6-36)
# Note that approximation via pdf/pmf is in-correct, as pointed out in later example
data <- data.frame( x=c(0:50))
ggplot(data, aes(x=x, y=dbinom(x, size=5, prob=0.7))) + geom_bar(stat="identity") + stat_function(fun=dnorm, args=list(mean=3.5, sd=1.05^(1/2))) + ylim(0, 0.4)
ggplot(data, aes(x=x, y=dbinom(x, size=10, prob=0.7))) + geom_bar(stat="identity") + stat_function(fun=dnorm, args=list(mean=7, sd=2.1^(1/2))) + ylim(0, 0.4)
ggplot(data, aes(x=x, y=dbinom(x, size=20, prob=0.7))) + geom_bar(stat="identity") + stat_function(fun=dnorm, args=list(mean=14, sd=4.2^(1/2))) + ylim(0, 0.4)
ggplot(data, aes(x=x, y=dbinom(x, size=50, prob=0.7))) + geom_bar(stat="identity") + stat_function(fun=dnorm, args=list(mean=35, sd=10.5^(1/2))) + ylim(0, 0.4)

# An example that approximation via pdf/pmf is in-correct. (LNp.6-37)
# Assume a transformation of Y applies to both binomial(size=50, prob=0.7) and
#   normal(mean=35, sd=10.5^(1/2)), with Y=g(x)=X/2 g-1(y)=2Y

# According to LNp.6-8~6-10, the transformed normal's pdf is as follows:
#  fy(y) = fx(g-1(y)) * |dg-1(y) / dy|
#  Fy(y) = Fx(g-1(y))     If g is strictly increasing
g <- function(x) {3*x+2}
g_inv <- function(y) {2*y}
norm_fx <- function(x) {dnorm(x, mean=35, sd=10.5^(1/2))}
norm_fy <- function(y) {norm_fx(g_inv(y))*2}
norm_Fx <- function(x) {pnorm(x, mean=35, sd=10.5^(1/2))}
norm_Fy <- function(y) {norm_Fx(g_inv(y))}

# The transformed binomial's pmf is as follows:
bino_px <- function(x) {dbinom(x, size=50, prob=0.7)}
bino_py <- function(y) {bino_px(g_inv(y))}
bino_Fx <- function(x) {pbinom(x, size=50, prob=0.7)}
bino_Fy <- function(y) {bino_Fx(g_inv(y))}

# Biplot transformed binomail pmf and normal pdf, to observe the difference
data <- data.frame( x=seq(0,50), by=0.5)
ggplot(data, aes(x=x, y=bino_py(x))) + geom_bar(stat="identity") + stat_function(fun=norm_fy, color="red")

# Biplot transformed binomail cdf and normal cdf, to observe the difference
ggplot(data, aes(x=x, y=bino_Fy(x))) + geom_bar(stat="identity") + stat_function(fun=norm_Fy, color="red")


# Continuity correction (LNp.6-37~6-38)
# pmf of binomial(50, 0.4), mean=20, variance=12
# cdf of normal(20, 12)
dbinom(18, size=50, prob=0.4)
pnorm(18.5, mean=20, sd=12^(1/2)) - pnorm(17.5, mean=20, sd=12^(1/2))

sum(dbinom(c(30:100), size=50, prob=0.4))
1-pnorm(29.5, mean=20, sd=12^(1/2))

##############################################################################################
# Weibull distribution:
#
#   建構在exponential(1)的基礎上，做平移、scale、指數的transformation
#   不具有memoryless
#
#   pdf:      dweibull(x, shape, scale = 1, log = FALSE)
#   cdf:      pweibull(q, shape, scale = 1, lower.tail = TRUE, log.p = FALSE)
#   quantile: qweibull(p, shape, scale = 1, lower.tail = TRUE, log.p = FALSE)
#   random:   rweibull(n, shape, scale = 1)
#
##############################################################################################
x <- c(0:5)
shape <- 2  # b
scale <- 2    # a

weibull <- data.frame( x=x, pdf=dweibull(x, shape, scale), cdf=pweibull(x, shape, scale))
weibull

x_fx <- function(x, shape=2, scale=2){ x*dweibull(x, shape=shape, scale=scale) }
(mean = integrate(x_fx, lower=-Inf, upper=Inf))
mean

bias_fx <- function(x, shape=2, scale=2, mean=1.772){ (x-mean)^2*dweibull(x, shape=shape, scale=scale) }
(var = integrate(bias_fx, lower=-Inf, upper=Inf))

ggplot(weibull, aes(x=x)) + stat_function(fun=dweibull, args=list(shape=shape, scale=scale))
ggplot(weibull, aes(x=x)) + geom_area(stat = "function", fun = dweibull, args=list(shape=shape, scale=scale), fill = "black")
ggplot(weibull, aes(x=x)) + stat_function(fun=pweibull, args=list(shape=shape, scale=scale))

# LNp6-41, the a/b discussion
ggplot(weibull, aes(x=x)) + stat_function(fun=dweibull, args=list(scale=0.5, shape=2), colour="red") +
    stat_function(fun=dweibull, args=list(scale=1.0, shape=2), colour="green") +
    stat_function(fun=dweibull, args=list(scale=1.5, shape=3), colour="blue") +
    stat_function(fun=dweibull, args=list(scale=3.0, shape=4), colour="purple")

# My observation, fixed scale (a)
ggplot(weibull, aes(x=x)) + stat_function(fun=dweibull, args=list(scale=3.0, shape=2), colour="red") +
    stat_function(fun=dweibull, args=list(scale=3.0, shape=3), colour="green") +
    stat_function(fun=dweibull, args=list(scale=3.0, shape=4), colour="blue") +
    stat_function(fun=dweibull, args=list(scale=3.0, shape=5), colour="purple")

# My observation, fixed shape (b)
ggplot(weibull, aes(x=x)) + stat_function(fun=dweibull, args=list(scale=0.5, shape=2), colour="red") +
    stat_function(fun=dweibull, args=list(scale=1.0, shape=2), colour="green") +
    stat_function(fun=dweibull, args=list(scale=1.5, shape=2), colour="blue") +
    stat_function(fun=dweibull, args=list(scale=2.0, shape=2), colour="purple")

##############################################################################################
# Cauchy distribution:
#
#   與Normal一樣有location和scale，但是有長尾特性
#
#   pdf:      dcauchy(x, location = 0, scale = 1, log = FALSE)
#   cdf:      pcauchy(q, location = 0, scale = 1, lower.tail = TRUE, log.p = FALSE)
#   quantile: qcauchy(p, location = 0, scale = 1, lower.tail = TRUE, log.p = FALSE)
#   random:   rcauchy(n, location = 0, scale = 1)
#
##############################################################################################
x <- c(-15:15)
location <- 5
scale <- 0.5

cauchy <- data.frame( x=x, pdf=dcauchy(x, location, scale), cdf=pcauchy(x, location, scale))
cauchy

x_fx <- function(x, location=5, scale=0.5){ x*dcauchy(x, location=location, scale=scale) }
(mean = integrate(x_fx, lower=-Inf, upper=Inf))
mean

bias_fx <- function(x, location=5, scale=0.5, mean=5){ (x-mean)^2*dcauchy(x, location=location, scale=scale) }
(var = integrate(bias_fx, lower=-Inf, upper=Inf))

ggplot(cauchy, aes(x=x)) + stat_function(fun=dcauchy, args=list(location=location, scale=scale))
ggplot(cauchy, aes(x=x)) + geom_area(stat = "function", fun = dcauchy, args=list(location=location, scale=scale), fill = "black")
ggplot(cauchy, aes(x=x)) + stat_function(fun=pcauchy, args=list(location=location, scale=scale))

# LNp6-43, the comparison with standard normal distribution
ggplot(cauchy, aes(x=x)) + stat_function(fun=dcauchy, args=list(location=0, scale=1)) + stat_function(fun=dnorm, args=list(mean=0, sd=1), colour="red")

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

# However, this is not easy to understand visually
joint_pdf_df <- as.data.frame(joint_pdf)
colnames(joint_pdf_df) <- c("X1", "X2", "X3", "X4", "X5", "X6", "X7", "X8", "X9", "X10", "X11")
ggplot(joint_pdf_df, aes(x=c(1:11), y=X1)) + geom_line() + geom_line(aes(y=X2)) + geom_line(aes(y=X3)) + geom_line(aes(y=X11))


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
data <- ddply(data, .(x1, x2, x3), summarize, joint_prob = dmultinom(c(x1, x2, x3), size = NULL, prob = c(1/3, 1/3, 1/3)) )
sum(data$joint_prob)

# Plot the joint pmf
# The points are only available at diagonal, since x1 + x2 = 3
ggplot(data, aes(x=x1, y=x2)) + geom_point() + geom_tile(aes(fill=joint_prob))


##############################################################################################
# Joint distribution transformation - method of event
#
#   For 2 independent Poisson distribution X ~ Poisson(l1) and Y ~ Poisson(l2),
#   the tranformation of X+Y=Z will be Z ~ Poisson(l1+l2)
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
l1 <- 5
x.pois <- data.frame( x=x, pmf=dpois(x, l1), cdf=ppois(x, l1))
x.pois
ggplot(x.pois, aes(x=x, y=pmf)) + geom_bar(stat="identity")

y <- c(1:20)
l2 <- 10
y.pois <- data.frame( y=y, pmf=dpois(x, l2), cdf=ppois(x, l2))
y.pois
ggplot(y.pois, aes(x=x, y=pmf)) + geom_bar(stat="identity")

x.pois
y.pois

z <- c(1:20)
l3 <- l1 + l2
z.pois <- data.frame( z=z, pmf=dpois(x, l3), cdf=ppois(x, l3))
z.pois
ggplot(z.pois, aes(x=z, y=pmf)) + geom_bar(stat="identity")

##############################################################################################
# Misc:
#
#   Create customized discrete distribution function
#
##############################################################################################
# Define full suite of functions (d*, p*, q*, r*) for your distribution
D <- DiscreteDistribution (supp = c(0, 1, 2) , prob = c(0.5, .25, .25))

dD <- d(D)  ## Density function
pD <- p(D)  ## Distribution function
qD <- q(D)  ## Quantile function
rD <- r(D)  ## Random number generation

# Take them for a spin
dD(-1:3)
# [1] 0.00 0.50 0.25 0.25 0.00
pD(-1:3)
# [1] 0.00 0.50 0.75 1.00 1.00
qD(seq(0,1,by=0.1))
# [1] 0 0 0 0 0 0 1 1 2 2 2
rD(20)
# [1] 0 0 2 2 1 0 0 1 0 1 0 2 0 0 0 0 1 2 1 0

test <- data.frame(x=c(0:10))
ggplot(test, aes(x=x)) + stat_function(fun=dD)
