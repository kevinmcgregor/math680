#Testing function for question 1
source("~/Documents/mcgill/math680/assignment5/bern_mcmc.R")

#########################
#Part (c)
set.seed(2384)
n = 100
p = 3
Ti = rep(20,n)
X = list()
for (i in 1:n) {
  X[[i]] = cbind(1,rnorm(Ti[i],0.1,1),rnorm(Ti[i],1,sqrt(2)))
}

beta = c(0,1,-0.2)
vu = 2.89
u = rnorm(n, 0, vu)

y = list()
for (i in 1:n) {
  y[[i]] = rbinom(Ti[i], 1, i.logit(X[[i]]%*%beta+u[i]))
}

#Testing the mcmc function
test_post = bern_mcmc(X,y,10000,Ti,1.01,1.01,100)

#Estimated standard error for the betas
beta_std_err = c(sd(test_post$beta.post[,1]),sd(test_post$beta.post[,2]),sd(test_post$beta.post[,3]))

#Diagnostic plots
par(mfrow=c(2,2), mar=rep(2,4))
plot(cumsum(test_post$beta.post[,1])/(1:length(test_post$beta.post[,1])), type="l")
plot(cumsum(test_post$beta.post[,2])/(1:length(test_post$beta.post[,2])), type="l")
plot(cumsum(test_post$beta.post[,3])/(1:length(test_post$beta.post[,3])), type="l")
plot(cumsum(test_post$vu.post)/(1:length(test_post$vu.post)), type="l")

save(test_post, beta_std_err, file="~/Documents/mcgill/math680/assignment5/dat_a5_q1c.RData")

############################################
#Question 1 (d)
set.seed(233)
n2 = 100
p2 = 2
Ti2 = rep(15, n2)

X2 = list()
for (i in 1:n2) {
  X2[[i]] = cbind(1,rnorm(Ti2[i],2,2),rnorm(Ti2[i],-0.5,sqrt(3)))
}

beta2 = c(1,2,-2)
vu2 = 1
u2 = rnorm(n2, -1, vu)

y2 = list()
for (i in 1:n2) {
  y2[[i]] = rbinom(Ti2[i], 1, i.logit(X2[[i]]%*%beta2+u2[i]))
}

#Testing the mcmc function
test_post2 = bern_mcmc(X2,y2,10000,Ti2,1.01,1.01,100)

#Diagnostic plots
par(mfrow=c(2,2), mar=rep(2,4))
plot(cumsum(test_post2$beta.post[,1])/(1:length(test_post2$beta.post[,1])), type="l")
plot(cumsum(test_post2$beta.post[,2])/(1:length(test_post2$beta.post[,2])), type="l")
plot(cumsum(test_post2$beta.post[,3])/(1:length(test_post2$beta.post[,3])), type="l")
plot(cumsum(test_post2$vu.post)/(1:length(test_post2$vu.post)), type="l")

save(test_post2, file="~/Documents/mcgill/math680/assignment5/dat_a5_q1d.RData")








