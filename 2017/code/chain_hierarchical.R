#########################################################################
# This script provides routines to draw samples from the joint posterior 
# distribution of the two parameters ('alpha' and 'beta') in the hierarchical 
# chain binomial model of the measles outbreak data. The model unknowns 
# also include the household-specific escape probabilities and the unknown
# chain identity in households with outbreak size 3.  
#
#  (n1^j,n11^j,n111^j,n12^j)| q_j ~ Multinom(q_j^2,2q_j^2p_j,2q_jp_j^2,p_j^2)
#  q_j| alpha, beta               ~ Beta(alpha,beta)
#  (al,be)                        ~(alpha+beta)^(-5/2) (cf. the lecture notes)
#
# Last updated July 14, 2017/KA
##########################################################################

propose.par = function(cur.par,prop.unif,delta.w){
# This function samples a new candidate value symmetrically 
# around the current value of the parameter iterate
  
  return(cur.par + (0.5-prop.unif)*delta.w)
  }

log.likelihood = function(my.al,my.be,q.vec){
# This function calculates the value of the log-likelihood 
# function of the two parameters (alpha and beta), based on 
# the (current iterates) of the household-specific escape probabilities
#
# INPUT my.al = parameter alpha
#       my.be = parameter beta
#       q.vec = a vector of household-specific escape probabilities

 return(sum(log(dbeta(q.vec,my.al,my.be))))
 }

chain_hierarchical = function(mcmc.size){
# This is the main sampling routine. The data, the parameters of the 
# prior distributon and the 'tuning parameters' (= the proposal window widths)
# are hard-wired in this code. In addition to the two model parameters (alpha and beta), 
# the algorithm draws samples of the 334 escape probabilities and 
# the unknown realisations of transmission chains in the 275 households with 
# the final number infected as 3.
#
# INPUT mcmc.size = the number of MCMC iterations to be sampled
#
# OUTPUT chain_hierarchical = the MCMC sample of parameters 'alpha' and 'beta'
#                             of the Beta distribution for household specific
#                             escape probabilities

# Tuning parameters: the widths of the proposal windows
delta.alpha = 0.25 # for parameter alpha
delta.beta  = 0.40 # for parameter beta

# Step (a): Reserve space for the MCMC samples
al      = rep(0,mcmc.size)
be      = rep(0,mcmc.size)
q       = matrix(0,nrow=mcmc.size,ncol=334)
n111    = matrix(0,nrow=mcmc.size,ncol=275)

# Step (b): Initialize
cur.alpha = 1
cur.beta  = 1

al[1]    = cur.alpha
be[1]    = cur.beta
q[1,]    = rep(0.5,334)
n111[1,] = rbinom(275,1,2*0.5/(2*0.5+1))

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
 # households with chain 1
 for (j in 1:n1){
   q[i,j] = rbeta(1,2+ cur.alpha,cur.beta)
   }

 # households with chain 1->1
 for (j in 1:n11){
   q[i,j+n1] = rbeta(1,2+cur.alpha,1+cur.beta)
   }

 # households with chain 1->1->1 or 1->2
 # N.B. In household j, the (current iterate of the) number 
 #      of escapes is equal to n111[i-1,j]
 for (j in 1:275){
   q[i,j+n1+n11] = rbeta(1,n111[i-1,j]+cur.alpha,2+cur.beta)
   }


 #############################################################
 # Step (d) Update the household-specific 'frequencies' of chain 1->1->1
 #          through Gibbs steps
 #############################################################
 for (j in 1:N3){ 
 
  # Draw a new iterate of n111 from its full conditional distribution
  n111[i,j] = rbinom(1,1,2*q[i,j+n1+n11]/(2*q[i,j+n1+n11]+1))
  }

 #################################################################
 # Step (e) Update parameter 'alpha' in a Metropolis-Hastings step
 #################################################################

 # Propose a new value for parameter alpha
 new.alpha = propose.par(cur.alpha,runif(1),delta.alpha)

 # If the proposed value if within the range of the parameter 
 if (new.alpha >0){

   # The log acceptance ratio = the log likelihood ratio   
   log.likeratio.alpha = log.likelihood(new.alpha,cur.beta,q[i,]) - log.likelihood(cur.alpha,cur.beta,q[i,]) 

   # The log acceptance ratio = the log-posterior ratio (since a symmetric proposals are being made)
   log.acc.alpha = log.likeratio.alpha -(5/2)*log(new.alpha+cur.beta) + (5/2)*log(cur.alpha+cur.beta)
 
   # Metropolis step for alpha: determine whether the proposal is accepted
     if (log(runif(1)) < log.acc.alpha){
     cur.alpha = new.alpha
     }
   }

 #################################################################
 # Step (f): Update parameter 'beta' in a Metropolis-Hastings step
 #################################################################

 # Propose a new value for parameter beta
 new.beta  = propose.par(cur.beta,runif(1),delta.beta)

 # If the proposed value is withing the range of the parameter
 if (new.beta > 0) { 
  
   # The log likelihood ratio
   log.likeratio.beta = log.likelihood(cur.alpha,new.beta,q[i,]) - log.likelihood(cur.alpha,cur.beta,q[i,]) 
 
   # The log acceptance ratio = the log-posterior ratio (since a symmetric proposals are being made)
   log.acc.beta = log.likeratio.beta -(5/2)*log(cur.alpha+new.beta) + (5/2)*log(cur.alpha+cur.beta)
 
   # Metropolis step for beta: determine whether the proposal is accepted
   if (log(runif(1))< log.acc.beta){
      cur.beta = new.beta
      }
   }   

 ######################################################### 
 # Store the current iterates of parameters alpha and beta
 #########################################################
 al[i] = cur.alpha
 be[i] = cur.beta 
 
}

# The output: the MCMC samples
return(list(al=al,be=be))

}

# Call the main routine
mcmc.size = 10000
mcmc.sample = chain_hierarchical(mcmc.size)
mcmc.al = mcmc.sample$al
mcmc.be = mcmc.sample$be

# Discard burn-in samples
burn.in = 2000
mcmc.al = mcmc.al[(burn.in+1):mcmc.size]
mcmc.be = mcmc.be[(burn.in+1):mcmc.size]

# Plot the marginal posterior distributions of parameters alpha and beta
par(mfrow=c(1,2))
hist(mcmc.al,xlab='alpha',main='')
hist(mcmc.be,xlab='beta',main='')

# Plot the joint posterior distribution of alpha and beta
par(mfrow=c(1,1))
plot(mcmc.al,mcmc.be,xlab='alpha',ylab='beta')

# Plot the posterior distribution of the expected escape probability
hist(mcmc.al/(mcmc.al+ mcmc.be),breaks=20,xlab='expected escape probability',main='')

# Plot the posterior predictive distribution of the escape probability
qpost = rbeta((mcmc.size-burn.in),mcmc.al,mcmc.be)
hist(qpost,main="posterior predictive distribution of the escape probability",cex.main=1,xlab="predictive q",breaks=20)

