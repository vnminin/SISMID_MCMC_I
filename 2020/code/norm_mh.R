## This script illustrates the Metropolis-Hastings algorithm for
## approximating the standard normal distribution
## Author: Vladimir N. Minin
## last update: 07/19/2020

#' Generate the next state the Metropolis-Hastings (M-H) chain targeting N(0,1) distribution, using a uniform random walk proposal
#' 
#' @param cur_value Current state of M-H Markov chain 
#' @param tuning_par Tuning parameter of the uniform random walk, where a random increment ~ Unif[-tuning_par, tuning_par] is added to the current state of the Markov chain (detla in the notes)
#'
#' @return Numeric vector of length 2. The first element of the vector contains the next state of the M-H Markov chain. The second element contains 0 if the proposed values was rejected and 1 otherwise.
#'
#' @examples
#' unif_rw_next(0.4, 5.0)
unif_rw_next = function(cur_value, tuning_par){

  return_value = c(cur_value, 0)
  
  prop_value = cur_value + runif(1, min=-tuning_par, max=tuning_par)
  
  if (log(runif(1)) < (cur_value^2 - prop_value^2)/2){
    return_value = c(prop_value, 1)
  }

  return(return_value)      
}

mcmc_size = 10000
start_value = 5.0


mcmc_out_small_jumps = matrix(0, nrow=(mcmc_size), ncol=2)
colnames(mcmc_out_small_jumps) = c("state", "acc_status")

cur_state = c(start_value,1)
mcmc_out_small_jumps[1,] = cur_state

for (i in 2:mcmc_size){
  
  mcmc_out_small_jumps[i,] = unif_rw_next(mcmc_out_small_jumps[i-1,1], 5.0)
}

## acceptance probability
mean(mcmc_out_small_jumps[,"acc_status"])

## mean of the target distribution
mean(mcmc_out_small_jumps[,"state"])

par(mfrow=c(2,2))

## trace plot
plot(1:mcmc_size, mcmc_out_small_jumps[,"state"], type="l", xlab="Iteration", ylab="MCMC State", main="Trace Plot (small jumps)")

## histogram
hist(mcmc_out_small_jumps[,"state"], xlab="MCMC State", main="Target Distribution Histogram")


mcmc_out_large_jumps = matrix(0, nrow=(mcmc_size), ncol=2)
colnames(mcmc_out_large_jumps) = c("state", "acc_status")

cur_state = c(start_value,1)
mcmc_out_large_jumps[1,] = cur_state

for (i in 2:mcmc_size){
  
  cur_state = unif_rw_next(mcmc_out_large_jumps[i-1,1], 0.05)
  
  mcmc_out_large_jumps[i,] = cur_state
}

## acceptance probability
mean(mcmc_out_large_jumps[,"acc_status"])

## mean of the target distribution
mean(mcmc_out_large_jumps[,"state"])

## trace plot
plot(1:mcmc_size, mcmc_out_large_jumps[,"state"], type="l", xlab="Iteration", ylab="MCMC State", main="Trace Plot (large jumps)")

## histogram
hist(mcmc_out_large_jumps[,"state"], xlab="MCMC State", main="Target Distribution Histogram")
