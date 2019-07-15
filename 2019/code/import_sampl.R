## This script illustrates the Importance-Sampling algorithm
## Author: Vladimir N. Minin
## last update: 07/14/2019

## define a threshold value and number of Monte Carlo samples
my_const = 4.5
sim_size = 10000000

## true probability of interest
(true_prob = pnorm(my_const,lower.tail=FALSE))


## naive Monte Carlo estimate
naive_samples = ifelse(rnorm(sim_size) > my_const,1,0)
(naive_est = mean(naive_samples))


## importance sampling estimate
shift_exp = rexp(sim_size,rate=1)+my_const
is_samples = ifelse(shift_exp > my_const,1,0)*dnorm(shift_exp)*exp(shift_exp-my_const)
(is_est = mean(is_samples))

print(true_prob)
print(naive_est)
print(is_est)


## Computing Monte Carlo variances

(naive_var_mc = var(naive_samples)/sim_size)
(is_var_mc = var(is_samples)/sim_size)

## Computing Monte Carlo confidence intervals

## Naive (not very useful, usually variance estimate is degenerate):
c(naive_est-1.96*sqrt(naive_var_mc), naive_est+1.96*sqrt(naive_var_mc))

## Importance Sampling:
c(is_est-1.96*sqrt(is_var_mc), is_est+1.96*sqrt(is_var_mc))
