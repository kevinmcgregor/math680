#Testing R functions for assignment 4 question 3

source("normCauchyPrior.R")

n = 100
#Number of replications
n.rep = 1000
#Observed data point
x=1

#Running replictions of sampling algorithm and looking at estimated expectations. Use two different methods
#to estimate expected value.
cursamp = rep(0, n)
expect_vals = rep(0, n.rep)
imp_expect_vals = rep(0, n.rep)
for (i in 1:n.rep) {
  cursamp = normCauchyPrior(n, x)$samp
  #Estimate expected value directly from sample of posterior dist
  expect_vals[i] = mean(cursamp)
  #Estimate expected value using importance sampling
  imp_expect_vals[i] = importNormCauchy(n, x)
}

hist(expect_vals)
mean(expect_vals)
var(expect_vals)

hist(imp_expect_vals)
mean(imp_expect_vals)
var(imp_expect_vals)



