#################################################################################
##This script provides a routine to draws samples from the joint posterior 
## distribution of the escape probability 'q' and the frequency
## n111 of chain 1->1->1, based on a general data set for households of size
## three  where data on n1, n11, and N3 for available for the 
## outbreak. The algorithm is based on Gibbs 
## sampling (see the lecture notes).
##
## INPUT mcmc.size  = the number of MCMC iterations
##       alpha,beta = parameters of the Beta prior distribution
##                    of the escape probability    
##              n1  = frequency of chain 1
##              n11 =  frequency of chain 1->1
##              N3  =  frequency of chains with outbreak size 3 
##                      (the total frequency of chains 1->2 and 1->1->1)
##
## OUTPUT mychainGibbs = the MCMC sample of q and n111 as a list
##                     of two entries, each a vector of length mcmc.size
## 
## Last updated July 9, 2018/ MEH
#################################################################################
mychainGibbs = function(n1,n11,N3,mcmc.size,alpha,beta){
# Reserve space
q    = rep(0,mcmc.size)
n111 = rep(0,mcmc.size)

# Initialize the model unknowns
q[1]    = 0.5;                         # just an initial guess
n111[1] = round(N3*2*q[1]/(2*q[1]+1)) # the (rounded) expected value of n_111, 
                                       # given q =0.5 

# Draw MCMC samples
for (i in 2:mcmc.size){

   # Draw an iterate of q from its full conditional distribution
   q[i] = rbeta(1,2*n1+2*n11+n111[i-1]+alpha,n11+2*N3+beta)
  
   # Draw an iterate of n111 from its full conditional distribution
   n111[i] = rbinom(1,N3,2*q[i]/(2*q[i]+1))
   
   }

# The output: the MCMC samples
mychainGibbs = list(q=q,n111=n111)

}

## Run the sampler

## mcmc-size = 5000, alpha = 1, beta = 1

test1=mychainGibbs(34,25,275,5000,1,1)

# Assume a burn-in of 500 iterations. 

# Summary of the results 
summary(test1$q[501:5000])
summary(test1$n111[501:5000])

#pdf(file="/Users/betz/Documents/TexWork/MCMC/betz/Bayesintro/chaingibbs1.pdf", height=4.5, width=8.9)
# Set up to have 2 plots in one figure
par(mfrow=c(1,2),oma=c(0,0,0,0))

hist(test1$q[501:5000],main="",xlab="q")
hist(test1$n111[501:5000],main="",xlab="n111")

# or
hist(test1[[1]][501:5000],main="",xlab="q")
hist(test1[[2]][501:5000],main="",xlab="n111")


#dev.off()