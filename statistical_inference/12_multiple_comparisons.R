##############################################################################################
# Case study I: no true positives
##############################################################################################

# For each experiment, randomly generate 20 x and y, and calculate p-Value of linear
# regression.
# Then, do such experiment 1000 times, and collect 1000 p-Value
set.seed(1010093)
pValues <- rep(NA, 1000)
for(i in 1:1000){
    y <- rnorm(20)
    x <- rnorm(20)
    pValues[i] <- summary(lm(y ~ x))$coeff[2,4] # p-Value between y and x
}
pValues
hist(pValues)

# Even though we do know that x and y are independent, there are still 51 cases showing
# significant results; has relationship. (alpha=0.05)
sum(pValues < 0.05)

# Controls false positive rate
# Controls FWER
sum(p.adjust(pValues, method="bonferroni") < 0.05)

# Controls FDR
sum(p.adjust(pValues, method="BH") < 0.05)

##############################################################################################
# Case study II: 50% true positives
##############################################################################################

# For first 500 experiments, randomly generate 20 x and y, and calculate p-Value of linear
# regression.
# For second 500 experiments, randomly generate 20 x and y=2*x, and calculate p-Value of linear
# regression.
set.seed(1010093)
pValues2 <- rep(NA, 1000)
for(i in 1:1000){
    y <- rnorm(20)
    # First 500 beta=0, last 500 beta=2
    if( i <= 500) { y <- rnorm(20)} else {y <-rnorm(20, mean=2*x)}
    pValues2[i] <- summary(lm(y ~ x))$coeff[2,4] # P-value between y and x
}
pValues2

trueStatus <- rep(c("zero", "not zero"), each=500)
tail(trueStatus)

#  The first 500 should have no correlation. However, 24 of them shows significant results.
#  The second 500 should have strong correlation, and the results support it.
#
#        trueStatus
#        not zero zero
#  FALSE        0  476
#  TRUE       500   24
#
table(pValues2 < 0.05, trueStatus)

# The false positive rate
24/500

# Controls FWER
#
# 1. No false positives (due to lower alpha)
# 2. 23 misidentified truly significant
#
#        trueStatus
#        not zero zero
#  FALSE       23  500
#  TRUE       477    0
table(p.adjust(pValues2, method="bonferroni") < 0.05, trueStatus)

# Controls BH
#
# 1. All the significant results were correctly identified
# 2. 13 results were incorrectly identified
# 3. Results in between no correction and FWER
#
#        trueStatus
#        not zero zero
#  FALSE        0  487
#  TRUE       500   13
table(p.adjust(pValues2, method="BH") < 0.05, trueStatus)

# Understanding p adjustment
par(mfrow=c(1,2))
plot(pValues2, p.adjust(pValues2, method="bonferroni"), pch=19)
plot(pValues2, p.adjust(pValues2, method="BH"), pch=19)




