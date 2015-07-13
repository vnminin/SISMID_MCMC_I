## This script illustrates the Importance-Sampling algorithm
## Author: Vladimir N. Minin
## last update: 07/12/2015

## define a threshold value and number of Monte Carlo samples
my.const = 4.5
sim.size = 100000

## true probability of interest
(true.prob = pnorm(my.const,lower.tail=FALSE))


## naive Monte Carlo estimate
(naive.est = mean(ifelse(rnorm(sim.size) > my.const,1,0)))


## importance sampling estimate
shift.exp = rexp(sim.size,rate=1)+my.const
is.est = mean(ifelse(shift.exp > my.const,1,0)*dnorm(shift.exp)*exp(shift.exp-my.const))

print(true.prob)
print(naive.est)
print(is.est)

## these functions can be used to estimate naive and importance sampling
## variances exactly using numerical integration

naive.var = function(my.const, sample.size){
  tail.prob = pnorm(my.const,lower.tail=FALSE)
  return(tail.prob*(1-tail.prob)/sample.size)
}

integrand = function(x,c){
  return(dnorm(x)^2/exp(c-x))
}


is.var = function(my.const, sample.size){
  tail.prob = pnorm(my.const,lower.tail=FALSE)

  return((integrate(integrand, lower = my.const, upper = Inf,c=my.const)$value - tail.prob^2)/sample.size)
}

## Computing Monte Carlo variances

var1 = naive.var(my.const, sim.size)
var2 = is.var(my.const, sim.size)

## Computing Monte Carlo confidence intervals

## Naive:
c(naive.est-1.96*sqrt(var1), naive.est+1.96*sqrt(var1))

## Importance Sampling:
c(is.est-1.96*sqrt(var2), is.est+1.96*sqrt(var2))

## True prob:
true.prob
