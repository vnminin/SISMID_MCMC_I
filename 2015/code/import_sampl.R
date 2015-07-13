## This script illustrates the Importance-Sampling algorithm
## Author: Vladimir N. Minin
## last update: 07/12/2015

## define a threshold value and number of Monte Carlo samples
my_const = 4.5
sim_size = 100000

## true probability of interest
(true_prob = pnorm(my_const,lower.tail=FALSE))


## naive Monte Carlo estimate
(naive_est = mean(ifelse(rnorm(sim_size) > my_const,1,0)))


## importance sampling estimate
shift_exp = rexp(sim_size,rate=1)+my_const
is_est = mean(ifelse(shift_exp > my_const,1,0)*dnorm(shift_exp)*exp(shift_exp-my_const))

print(true_prob)
print(naive_est)
print(is_est)

## these functions can be used to estimate naive and importance sampling
## variances exactly using numerical integration

naive_var = function(my_const, sample_size){
  tail_prob = pnorm(my_const,lower.tail=FALSE)
  return(tail_prob*(1-tail_prob)/sample_size)
}

integrand = function(x,c){
  return(dnorm(x)^2/exp(c-x))
}


is_var = function(my_const, sample_size){
  tail_prob = pnorm(my_const,lower.tail=FALSE)

  return((integrate(integrand, lower = my_const, upper = Inf,c=my_const)$value - tail_prob^2)/sample_size)
}

## Computing Monte Carlo variances

var1 = naive_var(my_const, sim_size)
var2 = is_var(my_const, sim_size)

## Computing Monte Carlo confidence intervals

## Naive:
c(naive_est-1.96*sqrt(var1), naive_est+1.96*sqrt(var1))

## Importance Sampling:
c(is_est-1.96*sqrt(var2), is_est+1.96*sqrt(var2))

## True prob:
true_prob
