MH_SIS = function(complete_data,M){
# July 14, 2017/KA
#
# This function realizes a numerical MCMC sample from the joint posterior 
# distribution of the three model parameters ('la', 'mu' and 'initprob').
# The model is a binary Markov process (a simple "SIS model"), observed
# completely from time 0 to time T.
#
# INPUT complete_data = complete data (simulated by program simulateSIS_N)
#       M             = the number of MCMC iterations
#
# OUTPUT MH_SIS = a list of 3 vectors of length M, each containing the 
#                 MCMC samples for one of the model parameters. 
#                 The parameters are
#                    la       = force of infection
#                    mu       = rate of clearing infection
#                    initprob = probability of being infected at time 0 
# 
# N.B. The prior distributions of the rate parameters 'la' and 'mu' 
#      are Gamma(v1,v2) with a (very) small value of v1=v2. The values of 
#      v1 and v2 can be changed in the code (see the code in subroutine 
#      update_parameters). 
#
#      The prior of the 'initprob' parameter is uniform on (0,1).
#
# This function calls subroutine update_parameters.R.

# Initialize parameter iterates
la       = 0.5
mu       = 0.5
initprob = 0.5
par      = list(la,mu,initprob)

# Reserve space to store the iterates
laM       = rep(0,M)
muM       = rep(0,M)
initprobM = rep(0,M)

# The number of individuals in the data
N = length(complete_data)

# Metropolis-Hastings iterations
for (i in 1:M){

  #
  #print(i)

  # Update the parameters
  par = update_parameters(complete_data,par)

  # Store the i'th MCMC iterates 
  laM[i]         = par[[1]]
  muM[i]         = par[[2]]
  initprobM[i]   = par[[3]]

  }

# The output
MH_SIS = list(laM,muM,initprobM)

}



############################################
likelihoodSIS = function(data,la,mu,initprob,index1,index2){
# June 11, 2010/KA
#
# This function calculates the log-likelihood value in the simple
# continuous-time binary Markov process (the simple SIS model). 
# In other words, the function calculates the (log) density of the 
# complete data, using the event-history likelihood in the SIS model.
#
# Ref: the lecture on the analysis of recurrent infections.
#
# N.B. There is a separate Bernoulli model for the initial state
#      at time 0, parameter 'initprob' being the probability of
#      being in state '1' (being infected) at time 0
#
# INPUT data     = complete data
#       la       = the force of infection
#       mu       = the rate of clearing infection
#       initprob = the probability of being infected at time 0
#       index1, index2 = the likelihood contibutions are calculated
#                        for individuals with indices in the range 
#                        [index1,index2]

# Initialize the log-likelihood
loglike = 0

# Add individual log-likelihood contributions
for (i in index1:index2){

 # The sojourn times (lengths of episodes) for individual 'i'
 t  = diff(data[[i]][[1]])# the sojourn times
 le = length(t)           # the number of episodes (including the last, 
                          # censored one)

 # Add the likelihood contribution from the initial state
 if (data[[i]][[2]][1] == 0){
    loglike = loglike + log(1-initprob)
    ins = 1 # the first episode of being susceptible
    inf = 2 # the first episode of being infected
    }
 if (data[[i]][[2]][1] == 1){
    loglike = loglike + log(initprob)
    ins = 2 # susceptible # the first episode of being susceptible
    inf = 1 # infected    # the first episode of being infected
    }

 # Total time spent susceptible
 tsusc = 0 # Default
 if (le >= ins){
    tsusc = sum(t[seq(ins,le,2)])
    }

 # Total times spent infected
 tinf = 0 # Default
 if (le >= inf){
    tinf = sum(t[seq(inf,le,2)])
    }

 # The number of infections
 noinf = 0 # Default
 if (le >=2){
    noinf = sum(data[[i]][[2]][2:le] == 1)
    }

 # The number of clearances of infection
 noclear = 0 # Default
 if (le >=2){
     noclear = sum(data[[i]][[2]][2:le] == 0)
    }

 # Add the log-likelihood contribution from the state 
 # transitions of the individual
 loglike = loglike + noinf*log(la) + noclear*log(mu) -la*tsusc - mu*tinf;
 
 }

# The output
likelihoodSIS = loglike

}


####################################################
update_parameters = function(data,par){
# June 11, 2010/KA
#
# This function realizes one MCMC iterate from the joint posterior 
# distribution of the three model parameters ('la', 'mu' and 'initprob') 
# in a simple binary Markov process (a simple "SIS model") when complete 
# data have been observed.
# 
# INPUT data = the complete data
#       par  = the current parameter vector (la,mu,initprob) as a list
#
# OUTPUT = the updated parameter vector
#
# N.B. The prior distributions of the rate parameters 'la' and 'mu' 
#      are Gamma(v1,v2) with a (very) small value of v1=v2. The values of 
#      v1 and v2 can be changed in the code (see below)
#
#      The prior of the 'initprob' parameter is uniform on (0,1).
#
# This function calls subroutine likelihoodSIS.R.
# This function is a subroutine for functions MH_SIS.R and augmentSIS.R.

# The current MCMC iterate of the parameters
la       = par[[1]]
mu       = par[[2]] 
initprob = par[[3]]

# The parameters of the Gamma(v1,v2) prior for the rate parameters
# The prior mean is v1/v2 and the prior variance v1/(v2*v2)
v1 = 0.00001  # 1 
v2 = 0.00001  # 20

# The widths of proposal windows
Deltala   = 0.1
Deltamu   = 0.1
Deltaprob = 0.1

# The number of individuals in the data
N = length(data)

# Initialize the log-likelihood
likelihoodold = likelihoodSIS(data,la,mu,initprob,1,N)

#######################
# Update parameter 'la'
#######################

# Propose a new value for 'lambda'
proposedla = la + (0.5-runif(1,0,1))*Deltala

# The Metropolis-Hasting step
if (proposedla > 0){

     # The log-likelihood with the proposed value parameter 'la'
     likelihoodnew = likelihoodSIS(data,proposedla,mu,initprob,1,N)
     
     # The log-likelihood ratio
     logratio = likelihoodnew - likelihoodold
    
     # The acceptance ratio (= the log-posterior ratio)
     AR = logratio + log(dgamma(proposedla,v1,v2)) - log(dgamma(la,v1,v2))

     # Determine whether the proposal is accepted or not
     if (log(runif(1,0,1)) < AR){
        la            = proposedla    # accepted
        likelihoodold = likelihoodnew # remember the log-likelihood
     }

    }

#######################
# Update parameter 'mu'
#######################

# Propose a new value for 'mu'
proposedmu = mu + (0.5-runif(1,0,1))*Deltamu

# The Metropolis-Hastings step
if (proposedmu > 0){

     # The log-likelihood with the proposed value of parameter 'mu'
     likelihoodnew = likelihoodSIS(data,la,proposedmu,initprob,1,N)
     
     # The log-likelihood ratio
     logratio = likelihoodnew - likelihoodold

     # The acceptance ratio (= the log-posterior ratio)
     AR = logratio + log(dgamma(proposedmu,v1,v2)) - log(dgamma(mu,v1,v2))
    
     if (log(runif(1,0,1)) < AR){
        mu            = proposedmu
        likelihoodold = likelihoodnew
     }

    }

#############################
# Update parameter 'initprob'
#############################

# Propose a new value for parameter 'initprob'
proposedprob = initprob + (0.5-runif(1,0,1))*Deltaprob

# The Metropolis-Hastings step
if ((proposedprob > 0) & (proposedprob < 1)){

     # The log-likelihood with the proposed value of parameter 'initprob'
     likelihoodnew = likelihoodSIS(data,la,mu,proposedprob,1,N)

     # The log-posterior ratio (assuming a uniform prior)
     logratio      = likelihoodnew - likelihoodold     

     # The acceptance ratio (a uniform prior is assumed here)
     AR = logratio
    
     if (log(runif(1,0,1)) < AR){
        initprob = proposedprob
     }

    }

# The output 
output      = list(1:3)
output[[1]] = la
output[[2]] = mu
output[[3]] = initprob

# 
update_parameters = output

}

# (1) Simulate data
complete_data = simulateSIS_N(N=100,la=0.45,mu=0.67,initprob=0.40,T=12)

# (2) Realise MCMC samples from the posterior distribution
# of the two model parameters
mcmc.size = 1500
par = MH_SIS(complete_data,mcmc.size)

# Plot the sample paths
par(mfrow=c(1,2))
plot(par[[1]],type='l',xlab="iteration",ylab="la")
plot(par[[2]],type='l',xlab="iteration",ylab="mu")

# Calculate the 90% marginal posterior intervals
burn.in=500
no.retained.samples = mcmc.size-burn.in

# lambda
la.retain = par[[1]][(burn.in+1):mcmc.size]
la.retain = sort(la.retain)
mean(la.retain) # the marginal posterior mean
la.retain[round(0.05*no.retained.samples)] # 5% quantile of the marginal posterior
la.retain[round(0.95*no.retained.samples)] # 95% quantile of the marginal posterior

# mu
mu.retain = par[[2]][(burn.in+1):mcmc.size]
mu.retain = sort(mu.retain)
mean(mu.retain) # the marginal posterior mean
mu.retain[round(0.05*no.retained.samples)] # 5% quantile of the marginal posterior
mu.retain[round(0.95*no.retained.samples)] # 95% quantile of the marginal posterior

# Plot the joint posterior of parameters 'la' and 'mu'
par(mfrow=c(1,1))
la.retain = par[[1]][(burn.in+1):mcmc.size]
mu.retain = par[[2]][(burn.in+1):mcmc.size]
plot(la.retain,mu.retain,type='p',xlab='la',ylab='mu')