#Function to use metropolis-hastings to sample from beta dist

#b is half the width of the uniform trial distribution (must be between 0 and 1)
sampleBeta = function(n, alpha, beta, b) {
  if (b<=0 | b>=1) stop("b must be between 0 and 1")
  if (alpha<=0) stop("alpha must be positive")
  if (beta<=0) stop("beta must be positive")
  
  #Initial value
  val = 0.5
  
  samp=rep(0,n)
  for (i in 1:n) {
    #Sample from trial distribution
    new.val = runif(1, val-b, val+b)
    #Generate unif(0,1)
    U = runif(1)
    #Get acceptance prob
    if (U<acceptProb(val, new.val, alpha, beta)) {
      val = new.val
    }
    
    samp[i] = val  
  }
  
  return(samp)
}

#Helper function to get acceptance probability
acceptProb = function(x, y, alpha, beta) {
  tmp = ifelse( y>0 & y<1, (y/x)^(alpha-1)*((1-y)/(1-x))^(beta-1), 0 )
  return(min(1,tmp))
}

#Testing the function
test1 = sampleBeta(10000, 2, 2, 0.5)
mean(test1)
acf(test1, main="Autocorrelation: Beta MH sampling alpha=2, beta=2", lag.max = 1000, col="darkgreen")

test2 = sampleBeta(10000, 5, 2, 0.5)
mean(test2)
acf(test2, main="Autocorrelation: Beta MH sampling alpha=5, beta=2", lag.max = 1000, col="darkgreen")

par(mfrow=c(2,1))
acf(test1, main="Autocorrelation: Beta MH sampling alpha=2, beta=2", lag.max = 1000, col="darkgreen")
acf(test2, main="Autocorrelation: Beta MH sampling alpha=5, beta=2", lag.max = 1000, col="darkgreen")

save(test1, test2, file="~/Documents/mcgill/math680/assignment5/dat_a5_q3.R")

#Checking convergence
par(mfrow=c(2,1), mar=c(3,3,3,3))
plot(cumsum(test1)/(1:length(test1)), type="l"); abline(h=0.5, col="red", lty=2, lwd=1.5)
plot(cumsum(test2)/(1:length(test2)), type="l"); abline(h=5/7, col="red", lty=2, lwd=1.5)



