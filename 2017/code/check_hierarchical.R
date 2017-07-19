check_hierarchical = function(mcmc.sample,mcmc.burnin){
# July 8, 2017/KA
#
# This function draws samples from the posterior predictive 
# distribution of chain frequencies (n1, n11, n111, n12)
# in the hierarchical chain binomial model with household-specific 
# escape probabilities.
#
# The function plots the the samples from the (marginal) 
# distribution of frequencies (n1,n11) for comparison with the
# actually observed value.
#
# The fuction output is the posterior expectations of the four
# chain frequencies in a sample of 334 households
#
# INPUT mcmc.sample = the output from program chain_hierarchical.R
#       mcmc.burnin = the number of samples to be discarded from the MCMC sample
#
# OUTPUT check_hierarchical = the posterior predictive expectations of the 
#                             four frequencies

# The number of MCMC samples
mcmc.size = length(mcmc.sample$al)

# Discard the burn-in samples
al = mcmc.sample$al[(mcmc.burnin+1):mcmc.size]
be = mcmc.sample$be[(mcmc.burnin+1):mcmc.size]

# The number of retained samples
mcmc.subsamplesize = mcmc.size-mcmc.burnin

# Reserve space for the predictive samples of the chain frequencies
Npred = matrix(0,mcmc.subsamplesize,4)

# Reserve space for one chain realisation in each of the 334 households  
L    = matrix(0,334,4) 
 
# Iterate over the MCMC samples
for (i in 1:mcmc.subsamplesize){

  # Sample household-specific escape probabilities and
  # the chains in 334 households
  for (j in 1:334){
    q     = rbeta(1,al[i],be[i])

    P1    = q*q
    P11   = 2*(1-q)*q^2
    P111  = 2*q*(1-q)^2
    P12   = (1-q)^2

    # A realisation of the chain
    L[j,] = rmultinom(1,1,c(P1,P11,P111,P12))
    }

    # Store the i'th MCMC sample of the chain frequencies
    Npred[i,] = apply(L,2,sum)
  } 

# Plot the sample points from the posterior predictive distribution 
# of frequencies n1 and n11
plot(0:60,0:60,type="n",xlab='n1',ylab='n11',cex.lab=2)
points(jitter(Npred[,1],2.5),jitter(Npred[,2],2.5),pch='.',cex=3)

# The actually observed value (n1,n11)
points(34,25,col='red')
abline(v=34,col='red',lwd=0.25)
abline(h=25,col='red',lwd=0.25)

# Return the posterior predictive expected frequencies
check_hierarchical = return(round(apply(Npred,2,mean)))

}

# Run the program to check the fit of the hierarchical model
# N.B. This requires that you have already realized a posterior
#      sample using chain_hierarchical.R

#mcmc.sample = chain_hierarchical(mcmc.size=10000)
check_hierarchical(mcmc.sample,mcmc.burnin=2000)



