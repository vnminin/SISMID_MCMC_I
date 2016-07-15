  checkmodel_reduced = function(mcmc.sample,mcmc.burnin){
# Last updated July 12, 2016/KA
#
# This function draws numerical samples from the joint posterior 
# predictive distribution of frequencies n1, n11, n111 and n12.
# Each individual sample is based on one realisation from the posterior
# distribution of the escape probability q. The escape probability is
# common across all households.
#
# The function outputs the mean of the predictive distribution, based on the 
# sample. It also plots the samples from the predictive distribution of 
# frequencies n1 and n11. 
#
# Ref: the lecture notes.
#
# INPUT mcmc.sample = the output (= MCMC samples) from program mychainGibbs.R 
#                     (or chainGibbs.R)
#       mcmc.burnin = the number of MCMC samples to be discarded
#
# OUTPUT expfreq = the posterior predictive expectations of the 
#                     four chain frequencies (n1,n11,n111,n12)

# The number of MCMC samples in the input
mcmc.size = length(mcmc.sample$q)

# Step (a): Discard the burn-in samples and retain the rest of the samples for 
#           the escape probability q
q = mcmc.sample$q[(mcmc.burnin+1):mcmc.size]

# The number of samples retained (= the size of the sub-sample)
mcmc.subsamplesize = mcmc.size-mcmc.burnin

# Step (b): Reserve space for the predictive frequencies
Npred = matrix(0,mcmc.subsamplesize,4)

# Step (c): Calculate the posterior predictive chain probabilities for
# each retained MCMC sample. Note that the following produces 4 vectors 
# of length 'mcmc.subsamplesize', i.e., the 4 chain probabilities for each MCMC sample. 
P1   =  q^2
P11  =  2*(1-q)*q^2   
P111 =  2*(1-q)^2*q
P12  = (1-q)^2

# Step (d): Sample from the posterior predictive distribution of chain frequencies
# with the total  number of chains as 334.
for (k in 1:mcmc.subsamplesize){
  Npred[k,] = ...  # YOUR TASK TO WRITE THE RIGHT-HAND SIDE OF THE EXPRESSION.
                   # Use the 'rmultinom' function.
  } 

# Step (f): Plot the posterior predictive distribution for frequencies n1 and n11
par(mfrow=c(1,1))
plot(0:60,0:60,type="n",xlab='n1',ylab='n11',cex.lab=2)
points(jitter(Npred[,1],2.5),jitter(Npred[,2],2.5),pch='.',cex=3)

# The actually observed values
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

# Running the functions: draw a posterior sample and check the model fit
mcmc.sample = chainGibbs(5000,1,1)

## or
#mcmc.sample = mychainGibbs(34,25,275,5000,1,1)

checkmodel_reduced(mcmc.sample,mcmc.burnin=500)
