#################################################################################
##This script provides a function to draw samples from the joint posterior 
## distribution of the escape probability 'q' and the frequency
## n111 of chain 1->1->1, based on the Providence measles 
## outbreak data. The algorithm is based on Gibbs 
## sampling (see the lecture notes).
##
## INPUT mcmc.size  = the number of MCMC iterations
##       alpha,beta = parameters of the Beta prior distribution
##                    of the escape probability
##
## OUTPUT chainGibbs = the MCMC sample of q and n111 as a list
##                     of two entries, each a vector of length mcmc.size
## 
## Updated June 17, 2010/KA
## Last updated July 9, 2018/ MEH
#################################################################################
chainGibbs = function(mcmc.size,alpha,beta){
# Reserve space
q    = rep(0,mcmc.size)
n111 = rep(0,mcmc.size)

# Initialize the model unknowns
q[1]    = 0.5;                         # just an initial guess
n111[1] = round(275*2*q[1]/(2*q[1]+1)) # the (rounded) expected value of n_111, 
                                       # given q =0.5 

# The observations (cf. the lecture)
n1  = 34   # frequency of chain 1
n11 = 25   # frequency of chain 1->1
N3  = 275  # frequency of chains with outbreak size 3 
           # (the total frequency of chains 1->2 and 1->1->1)

# Draw MCMC samples
for (i in 2:mcmc.size){

   # Draw an iterate of q from its full conditional distribution
   q[i] = rbeta(1,2*n1+2*n11+n111[i-1]+alpha,n11+2*N3+beta)
  
   # Draw an iterate of n111 from its full conditional distribution
   n111[i] = rbinom(1,N3,2*q[i]/(2*q[i]+1))
   
   }

# The output: the MCMC samples
chainGibbs = list(q=q,n111=n111)

}

## Run the sampler

## mcmc-size = 5000, alpha = 1, beta = 1

test=chainGibbs(5000,1,1)

#pdf(file="/Users/betz/Documents/TexWork/MCMC/betz/Bayesintro/chaingibbs1.pdf", height=4.5, width=8.9)
# Set up to have 2 plots in one figure
par(mfrow=c(1,2),oma=c(0,0,0,0))

# Assume a burn-in of 500 iterations, so just plot those greater than 500. 

hist(test$q[501:5000],main="",xlab="q")
hist(test$n111[501:5000],main="",xlab="n111")

# or
hist(test[[1]][501:5000],main="",xlab="q")
hist(test[[2]][501:5000],main="",xlab="n111")

# Summary of the results 
summary(test$q[501:5000])
summary(test$n111[501:5000])

#dev.off()