#testing gamma rejection sampling function

source("gammaSamp.R")

alpha = c(1.5, 5.3, 10.9, 50.2, 75.2)
beta=2
n=1000

s = matrix(0,n,length(alpha))
probs = rep(0, length(alpha))
for (i in 1:length(alpha)) {
  cat(i, "\n")
  s1 = gammaSamp(n, alpha[i], beta)
  s[,i] = s1$samp
  probs[i] = length(s1$samp)/s1$n.trial
}

library(xtable)
xtable(cbind(alpha,probs))
