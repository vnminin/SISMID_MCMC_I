## This script illustrates the Metropolis-Hastings algorithm for
## approximating the posterior distribution of the time of infection
## in a simple SIS model
## Author: Vladimir N. Minin
## last update: 07/19/2020

#' Compute log-likelihood of the disease trajectory of an individual who starts susceptible and become infected during the observation time period
#' 
#' @param inf_time Time of infection (in the notes denoted by I)
#' @param inf_rate Infection rate (lambda_1 in the notes)
#' @param clear_rate Clearance rate (lambda_2 in the notes)
#' @param total_time Length of observation time period (T in the notes)
#'
#' @return Numberic value of the log-likelihood: log(lambda_1) - lambda_1*I - lambda_2*(T-I)
#'
#' @examples
#' sis_log_like(0.4, 0.1, 2.2, 1.0)
sis_log_like = function(inf_time, inf_rate, clear_rate, total_time){
 return(log(inf_rate) - inf_rate*inf_time - clear_rate*(total_time-inf_time))
}

#' Compute log-likelihood of the disease trajectory of an individual who starts suscepible and become infected during the observation time period
#' 
#' @param cur_inf_time Current infection time (in the notes denoted by t_c)
#' @param inf_rate Infection rate (lambda_1 in the notes)
#' @param total_time Length of observation time period (T in the notes)
#' @param win_half_len half the length of the uniform proposal distribution interval (delta in the notes)
#'
#' @return Numberic value of the proposed infection time (t_p in the notes)
#'
#' @examples
#' sis_proposal(0.4, 1.0, 1.5)
sis_proposal = function(cur_inf_time, total_time, win_half_len){
   
  new_inf_time = cur_inf_time + runif(1,-win_half_len, win_half_len)
   
  if (new_inf_time < 0){
     new_inf_time = - new_inf_time
  }else{
    if (new_inf_time > total_time){
      new_inf_time = 2*total_time - new_inf_time 
    }
  }
   
  return(new_inf_time)
}

#' Run MCMC to approximate the posterior distribution of the infection time
#' 
#' @param start_inf_time Initial time of infection at MCMC iteration 1
#' @param inf_rate Infection rate (lambda_1 in the notes)
#' @param clear_rate Clearance reate (lambda_2 in the notes)
#' @param total_time Length of observation time period (T in the notes)
#' @param win_half_len half the length of the uniform proposal distribution interval (delta in the notes)
#' @param chain_len Number of MCMC iterations
#'
#' @return Numeric matrix with rows corresponding to MCMC iterations and columns to values of "inf_time", "log_like", and "acc_ind," where "acc_ind" is an indicator of whether the proposed infection time was accepted or rejected 
#'
#' @examples
#' inf_time_mcmc(start_inf_time=0.1, inf_rate=0.2, clear_rate=2, total_time=1.0, win_half_len=0.2, chain_len=10000)
inf_time_mcmc = function(start_inf_time, inf_rate, clear_rate, total_time, win_half_len, chain_len){
  
  result_mat = matrix(0, chain_len, 3)
  colnames(result_mat) = c("inf_time", "log_like", "acc_ind")
  
  cur_inf_time= start_inf_time
  result_mat[1,1] = start_inf_time
  result_mat[1,2] = sis_log_like(start_inf_time, inf_rate, clear_rate, total_time)
  
  for (i in 2:chain_len){
    prop_inf_time = sis_proposal(cur_inf_time, total_time, win_half_len)
    
    mh_log_ratio = sis_log_like(prop_inf_time, inf_rate, clear_rate, total_time) - 
      sis_log_like(cur_inf_time, inf_rate, clear_rate, total_time)
    
    
    if (log(runif(1)) < mh_log_ratio){
      result_mat[i,1] = prop_inf_time
      result_mat[i,2] = sis_log_like(prop_inf_time, inf_rate, clear_rate, total_time)
      result_mat[i,3] = 1
      cur_inf_time = prop_inf_time
    }else{
      result_mat[i,1] = cur_inf_time
      result_mat[i,2] = sis_log_like(cur_inf_time, inf_rate, clear_rate, total_time)
    }
  }
  
  return(result_mat)
}

## run the above functions

test_sample = inf_time_mcmc(start_inf_time=0.1, inf_rate=0.2, clear_rate=2, total_time=1.0, win_half_len=0.2, chain_len=10000)

summary(test_sample[1000:10000,])

hist(test_sample[1000:10000,1])
