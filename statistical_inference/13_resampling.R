library(UsingR)
data(father.son)
x <- father.son$sheight
n <- length(x)
B <- 10000 # 10000 bootstrap resamples
hist(x)

# We collect 1078 data (x)
# Then, we resample 10780000 data from the collected 1078 data, and arrange the
# data as the following matrix.
#
# resamples
#              1078 (n)
# 10000 (B)    ...
#
resamples <- matrix(sample(x, n*B, replace=TRUE), B, n)
head(resamples)
hist(resamples)

# Assume we are interested in median
resampledMedians <- apply(resamples, 1, median)
hist(resampledMedians)

# The estimated standard error of the median is
sd(resampledMedians)

# The 95% confidence interval of the estimated median is
# The instructor says this should be improved by BCa interval (correction for bias)
# BCA: Bias-corrected and accelerated interval
quantile(resampledMedians, c(0.025, 0.975))

# Histgram of bootstrap resamples
g = ggplot(data.frame(medians=resampledMedians), aes(x=medians))
g = g + geom_histogram(color="black", fill="lightblue", binwidth=0.05)
g

# Permutation tests
head(InsectSprays)
ggplot(InsectSprays, aes(factor(spray), count)) + geom_boxplot(aes(fill=factor(spray)))

# Subset only data with spary "B" and "C"
subdata <- InsectSprays[InsectSprays$spray %in% c("B", "C"),]
subdata

# Calculate the difference of mean of B and mean of C.
# y:     11 17 21 11 16 14 17 17 19 21  7 13  0  1  7  2  3  1  2  1  3  0  1  4
# group:  B  B  B  B  B  B  B  B  B  B  B  B  C  C  C  C  C  C  C  C  C  C  C  C
#
# observeStat: 13.25
y <- subdata$count
group <- as.character(subdata$spray)
testStat <- function(w, g) mean(w[g == "B"]) - mean(w[g == "C"])
observeStat <- testStat(y, group)
observeStat

# By null hypothesis: group label is unrelated to the outcome
permutations <- sapply(1 : 10000, function(i) testStat(y, sample(group)))
permutations


# p-Value closed to 0: reject null hypothesis
# => group level is related to the outcome
mean(permutations > observeStat)


g = ggplot(data.frame(permutations=permutations), aes(x=permutations))
g = g + geom_histogram(color="black", fill="lightblue", binwidth=1)
g = g + geom_vline(xintercept = 13)
g
