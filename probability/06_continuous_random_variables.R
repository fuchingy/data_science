library(ggplot2)
library(plyr)
require(rgl)      # For persp3d
library(distr)    # For DiscreteDistribution
library(numDeriv) # For grad

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
