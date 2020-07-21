## This script illustrates the Metropolis-Hastings algorithm for
## approximating the standard normal distribution
## Author: Vladimir N. Minin
## last update: 07/19/2020

## Your task: Add code the following function. 

#' Generate the next state the Metropolis-Hastings (M-H) chain targeting N(0,1) distribution, using a uniform random walk proposal
#' 
#' @param cur_value Current state of M-H Markov chain 
#' @param tuning_par Tuning parameter of the uniform random walk, where a random increment ~ Unif[-tuning_par, tuning_par] is added to the current state of the Markov chain (delta in the notes)
#'
#' @return Numeric vector of length 2. The first element of the vector contains the next state of the M-H Markov chain. The second element contains 0 if the proposed values was rejected and 1 otherwise.
#'
#' @examples
#' unif_rw_next(0.4, 5.0)
unif_rw_next = function(cur_value, tuning_par){

  return_value = c(cur_value, 0)
  
  return(return_value)      
}

mcmc_size = 10000
start_value = 3.0


mcmc_out = matrix(0, nrow=(mcmc_size), ncol=2)
colnames(mcmc_out) = c("state", "acc_status")

cur_value = c(start_value,1)
mcmc_out[1,] = cur_value

## Your task: add code to the below for loop to 
## fill in the mcmc.out matrix defined above with 
## the first column recording the state of the random 
## walk and the second column recording the acceptance status
## of each Metropolis-Hastings move. Don't forget to play 
## with the tuning parameter \delta.

for (i in 2:mcmc_size){
  mcmc_out[i,] = # using row mcmc_out[i-1,] and unif_rw_next
}

## acceptance probability
mean(mcmc_out[,"acc_status"])

## mean of the target distribution
mean(mcmc_out[,"state"])

## trace plot
plot(1:mcmc_size, mcmc_out[,"state"], type="l", xlab="Iteration", ylab="MCMC State")

## histogram
hist(mcmc_out[,"state"], xlab="MCMC State", main="Target Distribution Histogram")

