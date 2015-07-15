#########################################################################
# This script provides routines to draw samples from the joint posterior 
# distribution of parameters 'tildeq' and 'z' in a hierarchical chain 
# binomial model of the measles outbreak data. The model unknowns include
# the household-specific escape probabilities. 
#
# (n1^j,n11^j,n111^j,n12^j)| q_j ~ Multinom(q_j^2,2q_j^2p_j,2q_jp_j^2,p_j^2)
#  q_j| tildeq, z                ~ Beta(tildeq/z,(1-tildeq)/z)
#  tildeq                        ~ Uniform(0,1)
#  z                             ~ Gamma(1.5,1.5)
#
# Last updated July 11, 2015/KA
##########################################################################

propose.par = function(cur.par,prop.unif,delta.w){
# This function samples a new (candidate) value symmetrically 
# around the current value of the parameter iterate
   return(cur.par + (0.5-prop.unif)*delta.w)
 }

log.likelihood = function(my.tildeq,my.z,q.vec){
# This function calculates the log-likelihood function of the two parameters,
# based on the (current iterates) of the household-specific 
# escape probabilities
#
# INPUT my.tildeq = parameter tildeq
#       my.z.     = parameter z
#       q.vec     = a vector of household-specific escape probabilities

 return(sum(log(dbeta(q.vec,my.tildeq/my.z,(1-my.tildeq)/my.z))))

}

chain_hierarchical = function(mcmc.size){
# This is the main sampling routine. The data, the parameters of the 
# prior distributions and the 'tuning parameters' (= the proposal window widths)
# are hard-wired in this code. 
# 
# In addition to the two model parameters ('tildeq' and z), the 
# algorithm draws samples of the 334 escape probabilities and 
# the unknown realisations of chains in the 275 households with 
# the final number infected as 3.
#
# INPUT mcmc.size = the number of MCMC iterations to be sampled
#
# OUTPUT chain_hierarchical = the MCMC sample of parameters 'tildeq' and 'z'
#                             of the Beta distribution for household specific
#                             escape probabilities

# Parameters of the Gamma(nu1,nu2) prior for parameter z
nu1 = 1.5
nu2 = 1.5 

# Tuning parameters: the widths of the proposal windows
delta.tildeq = 0.08 # for parameter q
delta.z      = 0.3  # for parameter z

# Step (a): Reserve space for the MCMC samples
tildeq  = rep(0,mcmc.size)
z       = rep(0,mcmc.size)
q       = matrix(0,nrow=mcmc.size,ncol=334)
n111    = matrix(0,nrow=mcmc.size,ncol=275)

# Step (b): Initialize
cur.tildeq = 0.5
cur.z      = 1.0

tildeq[1]  = cur.tildeq
z[1]       = cur.z
q[1,]      = rep(0.5,334)
n111[1,]   = rbinom(1,275,2*0.5/(2*0.5+1))

# The observations
n1  = 34   # frequency of chain 1
n11 = 25   # frequency of chain 1->1
N3  = 275  # frequency of chains of size 3 (i.e., 1->2 or 1->1->1)

# Draw MCMC samples
for (i in 2:mcmc.size){

 ####################################################
 # Step (c): Update the household-specific escape probabilities
 #           through Gibbs steps
 ####################################################
 alpha =  cur.tildeq/cur.z     # tildeq[i-1]/z[i-1]
 beta  =  (1-cur.tildeq)/cur.z #(1-tildeq[i-1])/z[i-1]

 # households with chain 1
 for (j in 1:n1){
   q[i,j] = rbeta(1,2+ alpha,beta)
   }

 # households with chain 1->1
 for (j in 1:n11){
   q[i,j+n1] = rbeta(1,2+alpha,1+beta)
   }

 # households with chain 1->1->1 or 1->2
 for (j in 1:275){
   q[i,j+n1+n11] = rbeta(1,n111[i-1,j]+alpha,2+beta)
   }


 #############################################################
 # Step (d) Update the household-specific 'frequencies' of chain 1->1->1
 #          through Gibbs steps
 #############################################################
 for (j in 1:N3){ 
 
  # Draw a new iterate of n111 from its full conditional distribution
  n111[i,j] = rbinom(1,1,2*q[i,j+n1+n11]/(2*q[i,j+n1+n11]+1))
  }

 ################################################################
 # Step (e) Update parameter tildeq in a Metropolis-Hastings step
 ################################################################

 # Propose a new value for parameter tildeq
 new.tildeq = propose.par(cur.tildeq,runif(1),delta.tildeq)

 # If the proposed value if within the range of the parameter 
 if ((new.tildeq >0) & (new.tildeq < 1)){

   # The log acceptance ratio = the log likelihood ratio   
   #(a uniform prior is assumed for tildeq and a symmetric proposal has been made)
   log.acc.tildeq = log.likelihood(new.tildeq,cur.z,q[i,]) - log.likelihood(cur.tildeq,cur.z,q[i,]) 

   # Metropolis step for tildeq: determine whether the proposal is accepted
     if (log(runif(1)) < log.acc.tildeq){
     cur.tildeq = new.tildeq
     }
   }

 ############################################################
 # Step (f): Update parameter z in a Metropolis-Hastings step
 ############################################################

 # Propose a new value for parameter z
 new.z  = propose.par(cur.z,runif(1),delta.z)

 # If the proposed value is withing the range of the parameter
 if (new.z > 0) { 
  
   # The log likelihood ratio
   log.like.z = log.likelihood(cur.tildeq,new.z,q[i,]) - log.likelihood(cur.tildeq,cur.z,q[i,]) 
 
   # The log acceptance ratio = the log-posterior ratio (since a symmetric proposals are being made)
   log.acc.z = log.like.z + log(dgamma(new.z,nu1,nu2)) - log(dgamma(cur.z,nu1,nu2))
 
   # Metropolis step for z: determine whether the proposal is accepted
   if (log(runif(1))< log.acc.z){
      cur.z = new.z
      }
   }   

 ####################################################### 
 # Store the current iterates of parameters tildeq and z
 #######################################################
 tildeq[i] = cur.tildeq
 z[i]      = cur.z 
 
}

# The output: the MCMC samples
return(list(tildeq=tildeq,z=z))


}


# Running the program to draw a posterior sample of size 2000 
mcmc.sample = chain_hierarchical(mcmc.size=2000)

# Plot the posterior of tilde q
hist(mcmc.sample$tildeq[500:2000],xlab='tilde q',xlim=c(0.1,0.35))

# The prior and posterior predictive distributions of escape probability
par(mfrow=c(1,2))

tildeq.prior = runif(10000)
z.prior      = rgamma(10000,1.5,1.5)
qprior       = rbeta(10000,tildeq.prior/z.prior,(1-tildeq.prior)/z.prior)
hist(qprior,main="prior predictive q")

tildeq.sample = mcmc.sample$tildeq[501:2000]
z.sample      = mcmc.sample$z[501:2000]

qpost = rbeta(10000,tildeq.sample/z.sample,(1-tildeq.sample)/z.sample)
hist(qpost,main="posterior predictive q")

par(mfrow=c(1,1))
