##################################################################################
## This script provides routines for MCMC sampling from the posterior distribution 
## of the beta binomial hierarchical model parameters:
##                        x_i | \theta_i ~ Bin(n_i, \theta_i)
##                theta_i | \alpha, beta ~ Beta(\alpha, \beta) 
##                                \alpha ~ Exp(0.1)
##                                 \beta ~ Exp(0.1)
## Author: Vladimir N. Minin
## last update: 07/16/2017
##################################################################################


propose_pospar = function(cur_prospar, prop_unif, tuning_const){
  return(cur_prospar*exp(tuning_const*(prop_unif-0.5)))
}

log_prop_dens = function(prop_pospar){
  return(-log(prop_pospar))
}

log_posterior = function(my_alpha, my_beta, theta_vec, my_data, sample_size, prior_inten_alpha, prior_inten_beta){

  return(sample_size*(lgamma(my_alpha+my_beta)-lgamma(my_alpha)-lgamma(my_beta))+sum((my_data$x+my_alpha-1)*log(theta_vec))+sum((my_data$n-my_data$x+my_beta-1)*log(1-theta_vec))-prior_inten_alpha*my_alpha -prior_inten_beta*my_beta)  
}

## relevant terms of the posterior distribution that depend on alpha
## transformed to the log scale

log_alpha_posterior = function(my_alpha, my_beta, theta_vec, my_data, sample_size, prior_inten){
  ## contribution of alpha densities
  beta_dens_bit = sample_size*(lgamma(my_alpha+my_beta)-lgamma(my_alpha))+sum((my_data$x+my_alpha-1)*log(theta_vec))
  
  ## contribution of the prior
  prior_bit = -prior_inten*my_alpha

  return(beta_dens_bit+prior_bit)
}

log_beta_posterior = function(my_alpha, my_beta, theta_vec, my_data, sample_size, prior_inten){
  ## contribution of beta densities
  beta_dens_bit = sample_size*(lgamma(my_alpha+my_beta)-lgamma(my_beta))+ sum((my_data$n-my_data$x+my_beta-1)*log(1-theta_vec))
  
  ## contribution of the prior
  prior_bit = -prior_inten*my_beta

  return(beta_dens_bit+prior_bit)
}

mcmc_sampler = function(my_data, init_alpha, init_beta, prior_inten_alpha, prior_inten_beta,
  tuning_alpha, tuning_beta, mcmc_size, mcmc_burnin,mcmc_subsample){

  data_sample_size = dim(my_data)[1]

  ## prepare an MCMC output matrix
  mcmc_out = matrix(0, nrow=(mcmc_size-mcmc_burnin)/mcmc_subsample, ncol=data_sample_size+6)
  colnames(mcmc_out) = c("iter", "posterior", "alpha", "beta", "alpha.acc", "beta.acc",
            paste(rep("theta",data_sample_size),c(1:data_sample_size)))

  cur_alpha = init_alpha
  cur_beta = init_beta
  step_count = mcmc_subsample-1
  
  for (i in 1:mcmc_size){
    cur_alpha_acc = 0
    cur_beta_acc = 0
    
    ## Gibbs step for the components of theta
    ## TO DO: IMPLEMENT GIBSS STEP HERE
    cur_theta = rep(0.5, data_sample_size)

    
    ## Metropolis step for alpha
    ## TO DO: IMPLEMENT M-H STEP HERE FOR ALPHA

    ## Metropolis step for beta
    ## TO DO: IMPLEMENT M-H STEP HERE FOR BETA
    

    if (i > mcmc_burnin){
      
      step_count = step_count+1
      
      if(step_count==mcmc_subsample){
        my_index = (i-mcmc_burnin+mcmc_subsample-1)/mcmc_subsample
        mcmc_out[my_index,1]=i
        mcmc_out[my_index,2]=log_posterior(cur_alpha, cur_beta, cur_theta, my_data, data_sample_size, prior_inten_alpha,prior_inten_beta)
        mcmc_out[my_index,3]=cur_alpha
        mcmc_out[my_index,4]=cur_beta
        mcmc_out[my_index,5]=cur_alpha_acc
        mcmc_out[my_index,6]=cur_beta_acc
        mcmc_out[my_index,7:(data_sample_size+6)] = cur_theta
        
        step_count=0
      }
    }
  }

  return(mcmc_out)
}

## run the sampler

rat_data = read.table("https://raw.githubusercontent.com/vnminin/SISMID_MCMC_I/master/2016/code/rat_tumor.txt", header=TRUE)


rat_results = mcmc_sampler(my_data=rat_data,
             init_alpha=1.0,
             init_beta=1.0,
             prior_inten_alpha=0.1,
             prior_inten_beta=0.1,
             tuning_alpha=0.7,
             tuning_beta=0.7,
             mcmc_size=110000,
             mcmc_burnin=10000,
             mcmc_subsample=10)


print(mean(rat_results[,"alpha.acc"]))
print(mean(rat_results[,"beta.acc"]))


layout(matrix(c(1,1,1,1,2,2,3,3,4,5,5,6), 3, 4, byrow = TRUE),widths=c(1,0.5,0.5,1))

par(mar=c(5,5,4,1),cex.lab=2,cex.axis=1.5)
boxplot(as.data.frame(rat_results[,-c(1:6)]),
        ylab="Bin Success Prob",las=2,cex.axis=1.3)

hist(rat_results[,"alpha"],xlab=expression(alpha),main="")
box()
hist(rat_results[,"beta"],xlab=expression(beta), main="")
box()

par(mar=c(7,5,2,2))
plot(rat_results[,"posterior"],type="l",ylab="Log-Posterior",xlab="Iteration")
plot(rat_results[,"alpha"],type="l",ylab=expression(alpha),xlab="Iteration")
plot(rat_results[,"beta"],type="l",ylab=expression(beta),xlab="Iteration")
