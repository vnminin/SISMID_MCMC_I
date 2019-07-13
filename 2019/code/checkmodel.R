  checkmodel = function(mcmc.sample,mcmc.burnin){
# Last updated July 6, 2019
#
# This function draws numerical samples from the joint posterior 
# predictive distribution of frequencies n1, n11, n111 and n12.
# Each individual sample is based on one realisation (sample) from the posterior
# distribution of the escape probability q. The escape probability is assumed to 
# be the same across all 334 households.
#
# The function outputs the mean of the predictive distribution, based on the 
# sample. In addition, the joint predictive distribution of frequencies n1 and n11
# is plotted. 
#
# Ref: the lecture notes.
#
# INPUT mcmc.sample = the output (= MCMC samples) from program mychainGibbs.R 
#                     (or chainGibbs.R)
#       mcmc.burnin = the number of MCMC samples to be discarded
#
# OUTPUT expfreq = the posterior predictive expectations of the 
#                  four chain frequencies (n1,n11,n111,n12)

# The number of MCMC samples in the input
mcmc.size = length(mcmc.sample$q)

# Discard the burn-in samples and retain the rest of the samples for 
# the escape probability q
q = mcmc.sample$q[(mcmc.burnin+1):mcmc.size]

# The number of samples retained (= the size of the sub-sample)
mcmc.subsamplesize = mcmc.size-mcmc.burnin

# Reserve space for the predictive frequencies
Npred = matrix(0,mcmc.subsamplesize,4)

# Calculate the posterior predictive chain probabilities for
# each retained MCMC sample. This produces vectors of length 
# 'mcmc.subsamplesize': 
P1   =  q^2
P11  =  2*(1-q)*q^2   
P111 =  2*(1-q)^2*q
P12  = (1-q)^2

# Sample from the posterior predictive distribution of chain frequencies
# with the total  number of chains as 334.
for (k in 1:mcmc.subsamplesize){
  Npred[k,] = rmultinom(1,334,c(P1[k],P11[k],P111[k],P12[k]))
  } 

# Plot the posterior predictive distribution for frequencies n1 and n11
par(mfrow=c(1,1))
plot(0:60,0:60,type="n",xlab='n1',ylab='n11',cex.lab=2)           
points(jitter(Npred[,1],2.5),jitter(Npred[,2],2.5),pch='.',cex=3) 
 
# The actually observed value
points(34,25,col='red')
abline(v=34,col='red',lwd=0.25)
abline(h=25,col='red',lwd=0.25)

# The posterior predictive expected frequnecies of the four chain types
expfreq = round(334*c(mean(P1),mean(P11),mean(P111),mean(P12)))

# Return the posterior predictive expected frequencies
print("The posterior predictive expected frequencies of the four chain types")
expfreq

return(expfreq)

}

# Draw a sample from the posterior distribution of q 
mcmc.sample = chainGibbs(5000,1,1)

## or
# mcmc.sample = mychainGibbs(34,25,275,5000,1,1)

# Check the model fit
checkmodel(mcmc.sample,mcmc.burnin=500)






