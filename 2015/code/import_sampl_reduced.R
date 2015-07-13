## This script illustrates the Importance-Sampling algorithm
## Author: Vladimir N. Minin
## last update: 07/12/2015

## define a threshold value and number of Monte Carlo samples
my.const = 4.5
sim.size = 10000

## true probability of interest
(true.prob = pnorm(my.const,lower.tail=FALSE))


## naive Monte Carlo estimate
(naive.est = mean(ifelse(rnorm(sim.size) > my.const,1,0)))

## Your task: create a naive and an importance sampling estimate of the normal tail
## probability. To generate realizations from the shifted exponential
## use `rexp()` to generate regular exponentials and then add my.const to 
## them. Also, remember that you don't have to code the formula for the 
## normal density, because it is available via `dnorm()'.
## If you finish early, get Monte Carlo errors for naive and important sampling schemes. 