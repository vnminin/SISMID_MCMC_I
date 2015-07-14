##################################################################################
## This script provides routines for MCMC sampling from the posterior distribution 
## of the beta binomial hierarchical model parameters:
##                        x_i | \theta_i ~ Bin(n_i, \theta_i)
##                theta_i | \alpha, beta ~ Beta(\alpha, \beta) 
##                                \alpha ~ Exp(0.1)
##                                 \beta ~ Exp(0.1)
## Author: Vladimir N. Minin
## last update: 07/11/2014
##################################################################################


propose.pospar = function(cur.pospar, prop.unif, tuning.const){
  return(cur.pospar*exp(tuning.const*(prop.unif-0.5)))
}

log.prop.dens = function(prop.pospar){
  return(-log(prop.pospar))
}

log.posterior = function(my.alpha, my.beta, theta.vec, my.data, sample.size, prior.inten.alpha, prior.inten.beta){

  return(sample.size*(lgamma(my.alpha+my.beta)-lgamma(my.alpha)-lgamma(my.beta))+sum((my.data$x+my.alpha-1)*log(theta.vec))+sum((my.data$n-my.data$x+my.beta-1)*log(1-theta.vec))-prior.inten.alpha*my.alpha -prior.inten.beta*my.beta)  
}

## relevant terms of the posterior distribution that depend on alpha
## transformed to the log scale

log.alpha.posterior = function(my.alpha, my.beta, theta.vec, my.data, sample.size, prior.inten){
  ## contribution of alpha densities
  beta.dens.bit = sample.size*(lgamma(my.alpha+my.beta)-lgamma(my.alpha))+sum((my.data$x+my.alpha-1)*log(theta.vec))
  
  ## contribution of the prior
  prior.bit = -prior.inten*my.alpha

  return(beta.dens.bit+prior.bit)
}

log.beta.posterior = function(my.alpha, my.beta, theta.vec, my.data, sample.size, prior.inten){
  ## contribution of beta densities
  beta.dens.bit = sample.size*(lgamma(my.alpha+my.beta)-lgamma(my.beta))+ sum((my.data$n-my.data$x+my.beta-1)*log(1-theta.vec))
  
  ## contribution of the prior
  prior.bit = -prior.inten*my.beta

  return(beta.dens.bit+prior.bit)
}

mcmc.sampler = function(my.data, init.alpha, init.beta, prior.inten.alpha, prior.inten.beta,
  tuning.alpha, tuning.beta, mcmc.size, mcmc.burnin,mcmc.subsample){

  data.sample.size = dim(my.data)[1]

  ## prepare an MCMC output matrix
  mcmc.out = matrix(0, nrow=(mcmc.size-mcmc.burnin)/mcmc.subsample, ncol=data.sample.size+6)
  colnames(mcmc.out) = c("iter", "posterior", "alpha", "beta", "alpha.acc", "beta.acc",
            paste(rep("theta",data.sample.size),c(1:data.sample.size)))

  cur.alpha = init.alpha
  cur.beta = init.beta
  step.count = mcmc.subsample-1
  
  for (i in 1:mcmc.size){
    cur.alpha.acc = 0
    cur.beta.acc = 0
    
    ## Gibbs step for the components of theta
    ## TO DO: IMPLEMENT GIBSS STEP HERE
    cur.theta = rep(0.5, data.sample.size)

    
    ## Metropolis step for alpha
    ## TO DO: IMPLEMENT M-H STEP HERE FOR ALPHA

    ## Metropolis step for beta
    ## TO DO: IMPLEMENT M-H STEP HERE FOR BETA
    

    if (i > mcmc.burnin){
      
      step.count = step.count+1
      
      if(step.count==mcmc.subsample){
        my.index = (i-mcmc.burnin+mcmc.subsample-1)/mcmc.subsample
        mcmc.out[my.index,1]=i
        mcmc.out[my.index,2]=log.posterior(cur.alpha, cur.beta, cur.theta, my.data, data.sample.size, prior.inten.alpha,prior.inten.beta)
        mcmc.out[my.index,3]=cur.alpha
        mcmc.out[my.index,4]=cur.beta
        mcmc.out[my.index,5]=cur.alpha.acc
        mcmc.out[my.index,6]=cur.beta.acc
        mcmc.out[my.index,7:(data.sample.size+6)] = cur.theta
        
        step.count=0
      }
    }
  }

  return(mcmc.out)
}

## run the sampler

rat.data = read.table("http://www.stat.washington.edu/vminin/SISMID/2012/code/rat_tumor.txt", header=TRUE)
num.col = dim(rat.data)[1]


rat.results = mcmc.sampler(my.data=rat.data,
             init.alpha=1.0,
             init.beta=1.0,
             prior.inten.alpha=0.1,
             prior.inten.beta=0.1,
             tuning.alpha=0.7,
             tuning.beta=0.7,
             mcmc.size=110000,
             mcmc.burnin=10000,
             mcmc.subsample=10)


print(mean(rat.results[,"alpha.acc"]))
print(mean(rat.results[,"beta.acc"]))


layout(matrix(c(1,1,1,1,2,2,3,3,4,5,5,6), 3, 4, byrow = TRUE),widths=c(1,0.5,0.5,1))

par(mar=c(5,4,4,1))
boxplot(as.data.frame(rat.results[,-c(1:6)]),
        ylab="Bin Success Prob",las=2,cex.axis=0.8)

hist(rat.results[,"alpha"],xlab=expression(alpha),main="")
box()
hist(rat.results[,"beta"],xlab=expression(beta), main="")
box()

par(mar=c(7,4,2,2))
plot(rat.results[,"posterior"],type="l",ylab="Log-Posterior",xlab="Iteration")
plot(rat.results[,"alpha"],type="l",ylab=expression(alpha),xlab="Iteration")
plot(rat.results[,"beta"],type="l",ylab=expression(beta),xlab="Iteration")
