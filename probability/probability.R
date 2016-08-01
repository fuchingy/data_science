library(ggplot2)
library(plyr)
require(rgl)    # for persp3d
library(distr)  # DiscreteDistribution

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
x1.pmf(sqrt(4)) + x1.pmf(-sqrt(4)) == y1.pmf(4)

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
max <- 9

unif <- data.frame( x=x, pdf=dunif(x, min, max), cdf=punif(x, min, max))
unif

ggplot(unif, aes(x=x)) + stat_function(fun=dunif, args=list(min=min, max=max))
ggplot(unif, aes(x=x)) + geom_area(stat = "function", fun = dunif, args=list(min=min, max=max), fill = "black")
ggplot(unif, aes(x=x)) + stat_function(fun=punif, args=list(min=min, max=max))

# {TBD} integrate example
integrand <- function(x){x*punif(x, min, max)}
integrate(integrand,1, 10)
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

ggplot(exp, aes(x=x)) + stat_function(fun=dexp, args=list(rate=rate))
ggplot(exp, aes(x=x)) + geom_area(stat = "function", fun = dexp, args=list(rate=rate), fill = "black")
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
shape <- 1
rate <- 1

gamma <- data.frame( x=x, pdf=dgamma(x, shape, rate), cdf=pgamma(x, shape, rate))
gamma

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
x <- c(0:20)
mean <- 5
sd <- 0.5

norm <- data.frame( x=x, pdf=dnorm(x, mean, sd), cdf=pnorm(x, mean, sd))
norm

ggplot(norm, aes(x=x)) + stat_function(fun=dnorm, args=list(mean=mean, sd=sd))
ggplot(norm, aes(x=x)) + geom_area(stat = "function", fun = dnorm, args=list(mean=mean, sd=sd), fill = "black")
ggplot(norm, aes(x=x)) + stat_function(fun=pnorm, args=list(mean=mean, sd=sd))

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
scale <- 0.5    # a

weibull <- data.frame( x=x, pdf=dweibull(x, shape, scale), cdf=pweibull(x, shape, scale))
weibull

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
