#Testing for Assign 4 Q3

ratio = function(m,s) {
  curve(exp(-1/2*(x-2)^2)*s*(1+((x-m)/s)^2)/(1+x^2), xlim=c(0,5))
}

