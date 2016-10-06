# When sample size is small, the sampling distribution of mean is Tn-1 dist
#   (X - u) / (S/sqrt(n)) ~ Tn-1
# The confidence interval: X +- tn-1 * S / sqrt(n)
#                        : Est +/- t-quantile *std error(Est)
#                        : t=(X'-mu)/(s/sqrt(n))
#                        : X' +/- t_(n-1)*s/sqrt(n)
# When sample size is large enough, the sampling distribution of mean is Normal
#   (X - u) / (S/sqrt(n)) ~ N(0, 1)
# The confidence interval: X +- z * S / sqrt(n)
#                        : Est +/- qnorm *std error(Est)
#                        : Z=(X'-mu)/(sigma/sqrt(n))


myplot <- function(df){
    d <- data.frame(y = c(dnorm(xvals), dt(xvals, df)),
                    x = xvals,
                    dist = factor(rep(c("Normal", "T"), c(k,k))))
    g <- ggplot(d, aes(x = x, y = y))
    g <- g + geom_line(size = 2, aes(colour = dist))
    print(g)
}

myplot(2)


# Note that with 2 degrees of freedom, you only have 3 data points.





library(ggplot2)
data(sleep)
head(sleep)

ggplot(sleep, aes(x=group, y=extra)) + geom_point(aes(colour=ID), size=5) + geom_line(aes(x=group))


g1 <- sleep$extra[1:10]; g2 <- sleep$extra[11:20]
(difference <- g2 - g1)

(mn <- mean(difference))

(s <- sd(difference))
n <- 10

# T confidence interval
# sqrt(n) * (Xn - u ) / Sn ~ Tn-1
mn + c(-1, 1) * qt(.975, n-1) * s / sqrt(n)
t.test(difference)
t.test(g2, g1, paired=TRUE)
t.test(extra ~ I(relevel(group, 2)), paired=TRUE, data=sleep)


library(datasets); data(ChickWeight); library(reshape2)
head(ChickWeight)
class(ChickWeight$Diet)
unique(ChickWeight$Diet)
unique(ChickWeight$Chick)

(wideCW <- dcast(ChickWeight, Diet + Chick ~ Time, value.var="weight"))
head(wideCW)

names(wideCW)[-(1:2)] <- paste("time", names(wideCW)[-(1:2)], sep="")
library(dplyr)
wideCW <- mutate( wideCW, gain=time21-time0)
head(wideCW)

ggplot(ChickWeight, aes(x=Time,y=weight, colour=as.factor(Diet))) + geom_line() + facet_grid(~Chick)


(wideCW14 <- subset(wideCW, Diet %in% c(1,4)))
rbind(
    t.test(gain ~ Diet, paired=FALSE, var.equal=TRUE, data=wideCW14)$conf,
    t.test(gain ~ Diet, paired=FALSE, var.equal=FALSE, data=wideCW14)$conf
)
