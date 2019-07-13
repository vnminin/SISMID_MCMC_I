MH_SIS = function(complete_data,M){
# July 9, 2019
#
# This function realizes an MCMC sample of size M from the joint posterior 
# distribution of the three model parameters ('la', 'mu' and 'initprob') in
# a simple SIS model, observed completely from time 0 to time T.
#
# N.B. The individual SIS processes are assumed to be independent. In particular
#      there is no transmission in this model. Note also that because of we assume
#      the complete data are available there is no need for data augmentation.
#
# INPUT complete_data = complete data (need first to be simulated using program simulateSIS_N)
#       M             = the number of MCMC iterations
#
# OUTPUT MH_SIS = a list of length M, each entry of the list containing one
#                 MCMC iterate for each of the three model parameters: 
#                    la       = force of infection
#                    mu       = rate of clearing infection
#                    initprob = probability of being infected at time 0 
# 
# N.B. The prior distributions of the rate parameters 'la' and 'mu' 
#      are Gamma(v1,v2) with a (very) small value of v1=v2. In the standard
#      R parameterization of the Gamma distibution this means that the priors are
#      flat (uninformative) in the range where the likelihood gives support to the
#      parameters.  
#
#      The values of v1 and v2 can be changed in the code (see the code in subroutine 
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
  if (i%%100 == 0) print(i)

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


####################################################
update_parameters = function(data,par){
# Junly 10, 2018/KA
#
# This function realizes an MCMC iterate from the joint posterior 
# distribution of the three model parameters ('la', 'mu' and 'initprob') 
# in a simple binary Markov process (a simple "SIS model") when complete 
# data have been observed.
# 
# INPUT data = the complete data
#       par  = the current parameter vector (la,mu,initprob) as a list
#
# OUTPUT = the updated parameter vector, i.e. the next MCMC iterate
#
# N.B. The prior distributions of the rate parameters 'la' and 'mu' 
#      are Gamma(v1,v2) with a (very) small value of v1=v2. The values of 
#      v1 and v2 can be changed in the code (see below)
#
#      The prior of the 'initprob' parameter is uniform on (0,1).
#
# This function is a subroutine of function MH_SIS.R.
# This function calls subroutine likelihoodSIS.R.

# The current MCMC iterate of the parameters
# N.B. These are also the default output values, changed only if
#      the corresponding Metropolis-Hastings steps below lead to
#      accepting the respective proposals.
la       = par[[1]]
mu       = par[[2]] 
initprob = par[[3]]

# The parameters of the Gamma(v1,v2) prior for the rate parameters
# The prior mean is v1/v2 and the prior variance v1/(v2*v2)
v1 = 0.00001  # You can test with different value of v1 and v2 
v2 = 0.00001  # 

# The widths of proposal windows (tuning parameters)
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

# Metropolis-Hasting step
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
# N.B. There is a separate Bernoulli ("coin-tossing") model for each individual's 
#      initial state at time 0, parameter 'initprob' being the probability of
#      being in state '1' (i.e. infected) at time 0
#
# INPUT data     = complete data
#       la       = the force of infection
#       mu       = the rate of clearing infection
#       initprob = the probability of being infected at time 0
#       index1, index2 = the likelihood contibutions are calculated
#                        for individuals with indices in the range 
#                        [index1,index2]. For the entire data, we use index1=1 and index2=N.
#
# This function is a subroutine of function update_parameters.

# Initialize the log-likelihood
loglike = 0

# Add individual log-likelihood contributions
for (i in index1:index2){

 # The sojourn times (lengths of episodes) for individual 'i'
 t  = diff(data[[i]][[1]])# the sojourn times
 le = length(t)           # the number of episodes (including the last and 
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


##########################################################################
##########################################################################

# (1) Simulate data for N individuals following the SIS dynamics from time
#     0 to 12, with acquisition and clearance rates la and mu, respectively.
#     The probability of infection at time 0 is initprob.
complete_data = simulateSIS_N(N=100,la=0.45,mu=0.67,initprob=0.40,T=12)

# (2) Realize MCMC samples from the posterior distribution of the model parameters
#     of the three model parameters, based on the simulated set of complete data
mcmc.size = 1500
par = MH_SIS(complete_data,mcmc.size)

# Plot the sample paths
par(mfrow=c(1,3))
plot(par[[1]],type='l',xlab="iteration",ylab="la")
plot(par[[2]],type='l',xlab="iteration",ylab="mu")
plot(par[[3]],type='l',xlab="iteration",ylab="initprob")


# Calculate the 90% marginal posterior intervals
burn.in=500
no.retained.samples = mcmc.size-burn.in

# lambda
la.retain = par[[1]][(burn.in+1):mcmc.size]
la.retain = sort(la.retain)
round(mean(la.retain),2)                            # the marginal posterior mean
round(la.retain[round(0.05*no.retained.samples)],2) # 5% quantile of the marginal posterior
round(la.retain[round(0.95*no.retained.samples)],2) # 95% quantile of the marginal posterior

# mu
mu.retain = par[[2]][(burn.in+1):mcmc.size]
mu.retain = sort(mu.retain)
round(mean(mu.retain),2)                            # the marginal posterior mean
round(mu.retain[round(0.05*no.retained.samples)],2) # 5% quantile of the marginal posterior
round(mu.retain[round(0.95*no.retained.samples)],2) # 95% quantile of the marginal posterior

# initprob
initprob.retain = par[[3]][(burn.in+1):mcmc.size]
initprob.retain = sort(initprob.retain)
round(mean(initprob.retain),2)                            # the marginal posterior mean
round(initprob.retain[round(0.05*no.retained.samples)],2) # 5% quantile of the marginal posterior
round(initprob.retain[round(0.95*no.retained.samples)],2) # 95% quantile of the marginal posterior


# Plot the joint posterior of parameters 'la' and 'mu'
par(mfrow=c(1,1))
la.retain = par[[1]][(burn.in+1):mcmc.size]
mu.retain = par[[2]][(burn.in+1):mcmc.size]
plot(la.retain,mu.retain,type='p',xlab='la',ylab='mu')