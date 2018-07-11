## This script illustrates the Metropolis-Hastings algorithm for
## approximating the standard normal distribution
## Author: Vladimir N. Minin
## last update: 07/11/2018

## Your task: Add code the following function. The input of  
## the function is the current value of the random walk and 
## the tuning parameter (\delta in the lab notes). The output
## should be a vector with the first component being the next 
## value of the random walk and the second component being 0 if 
## the proposed value was rejected and 1 if accepted.

unif_rw_next = function(cur_value, tuning_par){

  return_value = c(cur_value, 0)
  
  return(return_value)      
}

mcmc_size = 10000
start.value = 3.0


mcmc_out = matrix(0, nrow=(mcmc_size), ncol=2)
colnames(mcmc_out) = c("state", "acc.status")

cur_value = c(start_value,1)

## Your task: add code to the below for loop to 
## fill in the mcmc.out matrix defined above with 
## the first column recording the state of the random 
## walk and the second column recording the acceptance status
## of each Metropolis-Hastings move. Don't forget to play 
## with the tuning parameter \delta.

for (i in 2:mcmc_size){
  

}

## acceptance probability
mean(mcmc_out[,"acc.status"])

## mean of the target distribution
mean(mcmc_out[,"state"])

## trace plot
plot(1:mcmc_size, mcmc_out[,"state"], type="l", xlab="Iteration", ylab="MCMC State")

## histogram
hist(mcmc_out[,"state"], xlab="MCMC State", main="Target Distribution Histogram")

