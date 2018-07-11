## This script illustrates the Metropolis-Hastings algorithm for
## approximating the posterior distribution of the time of infection
## in a simple SIS model
## Author: Vladimir N. Minin
## last update: 07/11/2018

sis_log_like = function(inf_time, inf_rate, clear_rate, total_time){
  return(log(inf_rate) - inf_rate*inf_time - clear_rate*(total_time-inf_time))
}

# finish this function
sis_proposal = function(cur_inf_time, total_time, win_half_len){
  
  new_inf_time = 0 ## replace with the proper proposal
  
  return(new_inf_time)
}



inf_time_mcmc = function(start_inf_time, inf_rate, clear_rate, total_time, win_half_len, chain_len){
  
  result_mat = matrix(0, chain_len, 3)
  colnames(result_mat) = c("inf_time", "log_like", "acc_ind")
  
  cur_inf_time= start_inf_time
  result_mat[1,1] = start_inf_time
  result_mat[1,2] = sis_log_like(start_inf_time, inf_rate, clear_rate, total_time)
  
  for (i in 2:chain_len){
    ## 1. Generate a new value of the infection time using 
    ##  the function sis_proposal()
    ## 2. Decide whether to accept or reject the proposed value by computing
    ## the Metropolis-Hastings ratio
    ## 3. Save the current or proposed value in result_mat[i,1]
    ##    Save the complete-data log-likelihood evaluated either at the current
    ##    or proposed value of the infection time in result_mat[i,2]
    ##    Save the indicator of the acceptance in result_mat[i,3]
  }
  
  return(result_mat)
}


## run the above functions
test_sample = inf_time_mcmc(start_inf_time=0.1, inf_rate=0.1, clear_rate=0.2, total_time=1.0, win_half_len=0.2, chain_len=10000)

summary(test_sample[1000:10000,])

hist(test_sample[1000:10000,1])

