# This source code contains the main program (sampleSIR.R) and
# all required subroutines.
# Last updated July 9, 2017/KA

######################
# (1) Program readdata
######################
readdata = function(){
# This function outputs the Abakaliki smallpox data.
# The data points are removal times of the 30 cases that
# occurred during the outbreak. Days are counted from 
# the start of the outbreak which is assumed to coincide
# the the index case's infection.

# The observed data as time intervals (in days) between the removals
return(c(14,27,34,36,39,39,39,40,44,49,52,54,54,56,56,61,64,65,69,69,70,71,72,74,74,75,80,80,85,90))
}

############################
# (2) Program initializedata
############################
initializedata = function(remtimes){
#
# This function initializes the complete data for the
# general epidemic model with observed removal times. 
# The complete data contains infection and removal times 
# of the 'n' individuals that were infected during the oubreak.
#
# For each individual, the initial value of the infection 
# time is taken to be halfway between time 0 and the time 
# of removal for the individual.
#
# N.B. There are M-n susceptible individuals that do not
#      acquire infection during the outbreak. These data are
#      not coded in the complete data matrix.
#
# INPUT remtimes        = the removal times as a vector (of length n)
# OUTPUT initializedata = a complete data matrix (nx2 matrix)

# Reserve space the complete data matrix
n            = length(remtimes) # The total number of infected during the outbreak
completedata = matrix(0,n,2)
completedata = as.data.frame(completedata)
colnames(completedata) = c("infectiontime","removaltime")

# The data for the index case is on the first row:
completedata[1,] = c(0,remtimes[1]) # infection at time 0, removal as in the input

# Initialize the infection times for cases 2,...,n
inftimes = rep(0,n)
for (k in 2:n){  
   inftimes[k]      = inftimes[k-1] + (remtimes[k-1]-inftimes[k-1])/2
   completedata[k,] = c(inftimes[k],remtimes[k])
   }

# The output: the complete data matrix
initializedata = completedata;

}



###########################
# (3) Program loglikelihood
###########################
loglikelihood = function(A,M,bet,gam){
#
# This function calculates the value of the log-likelihood.
# The likelihood is based on *individual* histories
# of infection and removal times in the general epidemic model.
# See the lecture notes and the instuctions to this exercise.
#
# INPUT A       = the complete data matrix: a nx2 matrix, including
#                 the infection and removal times for the n infected
#                 individuals in the data
#       M       = the total number of individuals in the community
#       bet,gam = the current iterates of the model parameters
#
# OUTPUT = the value of the log-likelihood.
#
# N.B.     If the complete data matrix is such that
#          the number of infectious or susceptibles is
#          zero during the outbreak, value -Inf is returned.
#          This notifies an impossible realisation in the model.
#
# N.B.     The likelihood function is based on *individual* event-histories.
#          This agrees with the actual augmentation routine (update_inftimes.R),
#          in which unknown *individual* infection times are sampled 
#          conditionally on each individual's removal time.
#
# This subroutine is called by function update_inftimes.R.

# The number of infected individuals during the entire outbreak, including the index case 
n = dim(A)[1]

# 'Recode' the infection and removal times
inftimes = cbind(rep(1,n),A[,1])    # 1 = infection
remtimes = cbind(rep(0,n),A[,2])    # 0 = removal
times    = rbind(inftimes,remtimes) 

# Sort the infection and removal times in ascending order
c     = sort(times[,2],index.return=TRUE);
times = times[c$ix,];

# The indicator vector for infection events
indicator = times[,1]

# The  indices of the infection times
infection = which(indicator==1)

# Reserve space: the numbers of infected at each of the n-1 infections
no_infected        = rep(0,n-1)

#
for (i in 2:n){

# The number of infections before the infection in qustion
ninf = sum(indicator[1:(infection[i]-1)] == 1) 

# The number of removals before the infection in question
nrem = sum(indicator[1:(infection[i]-1)] == 0)

# The number of infected individuals at the time of the infection in question
no_infected[i-1] = ninf-nrem

}

########################################################################
# The log-likeihood contributions: cf. the intructions to this practical
########################################################################

# The log-likelihood contribution from the removals
c1 = (n-1)*log(bet)

# The log-likelihood contribution from the infections
c2 = n*log(gam) + sum(log(no_infected))

# Calculate the contribution of the integral term 
c3 = -gam*totaltime_infected(A) - (bet/M)*totaltime_infpressure(A,M)

# Output: the log-likelihood
return(c1+c2+c3)
 
}



###################################
# (4) Program totaltime_infpressure
###################################
totaltime_infpressure = function(completedata,M){
#
# This function calculates the values of integral 
# int_0^T S(u)I(u)du. This is the total time of
# infectious pressure in the study population. It can also
# interpreted as the total time spent susceptible, weighted
# by the number of infectives.
#
# N.B. This function takes into account also the contribution
#      from the 90 individuals that did not acquire infection
#      during the outbreak.
#
# INPUT completedata = the complete data matrix
#       M            = the total number of individuals
#
# OUTPUT = the value of the integral
#
# This subroutine is called by function update_beta.R.
# This subroutine calls function totaltime_infected.R.

# The number of individuals with infection during the outbreak
Ninf = dim(completedata)[1] 

# The number of individuals the remain suseptible
Nsusc = M-Ninf    

# The infection and removal times 
inftimes = completedata$infectiontime
remtimes = completedata$removaltime

# Initialise the total time of infectious pressure
totaltime = 0

# Calculate the contribution to the total time of infectious pressure
# from those that acquire infection during the outbreak
for (k in 1:Ninf){
 for (j in 1:Ninf){
   t1        = min(remtimes[k],inftimes[j])
   t2        = min(inftimes[k],inftimes[j])
   totaltime = totaltime + (t1-t2)
 }
}


# Add the contribution of all those that remain susceptible throughout
# the outbreak
inftime   = totaltime_infected(completedata) # The total time spent infective
totaltime = totaltime + Nsusc*inftime        #

# The output
return(totaltime)
}

#################################
# (5) Program totaltime_infected
#################################
totaltime_infected = function(completedata){
#
# This function calculates the value of the integral int_0^T I(u)du.
# This is the total time spent infected during the oubreak
# in the study population.
#
# INPUT  completedata  = the complete data matrix
#
# OUTPUT    = the total time spent infected in the study population
#
# This subroutine is called by function update_gamma.R.

# The total sum of the durations of infected periods in the individuals 
return(sum(completedata$removaltime-completedata$infectiontime))

}


#############################
# (6) Program update_inftimes
#############################
update_inftimes = function(completedata,M,bet,gam){
#
# This function updates the unknown infection times
# in the general epidemic model with observed removal times.
# Each of the unknown times is updated separately using
# a Metropolis-Hastings step.
#
# INPUT completedata  = the complete data matrix
#       M             = the total number of individuals
#       bet,gam       = the current iterates of the model parameters
#
# OUTPUT update_times = the updated complete data matrix
#
# This function calls subroutine loglikelihood.R.

# The number of infection times, including the index case
n = length(completedata$infectiontime)

# Iterate over all unknown infection times
for (i in 2:n){

 # The value of the current log-likelihood
 loglikelihood.old = loglikelihood(completedata,M,bet,gam);

 # Propose a new infection time for individual i, uniformly
 # between 0 and the removal time of the individual
 infectiontime.new = runif(1,0,completedata[i,]$removaltime);

 # The modified data matrix: 
 # replace the current infection time with the proposal
 completedata.new                   = completedata;
 completedata.new[i,]$infectiontime = infectiontime.new;
 
 # The value of the log-likelihood under the modified process
 loglikelihood.new = loglikelihood(completedata.new,M,bet,gam);

 # Using the Metrolis-Hastings method, decide if the
 # proposal is accepted 
 if ( log(runif(1,0,1)) < (loglikelihood.new-loglikelihood.old)){
    completedata = completedata.new; # accept 
    }

 # N.B. If the proposal was not accepted, the old complete data
 # matrix will be returned
    
 }

# Return the updated complete data matrix
return(completedata);

}


##########################
# (7) Program update_beta
##########################
update_beta = function(completedata,M){
#
# This function updates parameter 'beta' from its 
# full conditional posterior distribution under
# complete data in the general epidemic model.
#
# INPUT completedata = the complete data matrix
#       M            = the total number of individuals
# 
# OUPUT update_beta  = a Gibbs iterate of parameter 'beta' 
#
# This function calls subroutine totaltime_infpressure.R.
# This function is called by function sampleSIR.R.

# The number of infections, excluding the index case
N.infections = dim(completedata)[1] - 1

# The parameters defining the Beta prior of parameter 'beta'
nu.beta = 0.0001 # 0.0001, 10 
la.beta = 0.0001 # 0.0001, 100

# The parameters of the full conditional (Beta) distribution of parameter 'beta'
al = N.infections + nu.beta;
be = (1/M)*totaltime_infpressure(completedata,M) + la.beta;

# Draw and return one Gamma(al,be) random variate
return(rgamma(1,al,be))

}

##########################
# (8) Program update_gamma
##########################
update_gamma = function(completedata){
#
# This function updates parameter 'gamma' from its
# full conditional posterior distribution under 
# complete data in the general epidemic model.
#
# INPUT completedata  = the complete data matrix
#
# OUTPUT update_gamma = a Gibbs update of parameter 'gamma'
#
# This function calls routines totaltime_infected.R.
# This function is called by functions sampleSIR.R 

# The number of removals during the outbreak
n = dim(completedata)[1]

# The parameters defining the Beta prior of parameter 'gam'
nu.gam = 0.0001  # 0.0001, 10
la.gam = 0.0001 # 0.0001, 100

# The parameters of the full conditional distribution of parameter 'gam'
al = n + nu.gam
be = totaltime_infected(completedata) + la.gam

# Draw and return one Gamma(al,be) random variate
return(rgamma(1,al,be))

}


########################
# (9) Program sampleSIR
########################
sampleSIR = function(remtimes,M,mcmc.size){
#
# This function draws samples from the joint posterior 
# distribution of the two model parameters ('beta' and 'gamma') and 
# the unobserved infection times in the general epidemic model with
# observed removal times.
# 
# The time origin t=0 is taken to be the time when the index
# case becomes infected. The index case's removal time is the
# first one in the input vector 'remtimes'.
#
# The analysis is based on augmenting the set of model unknowns 
# with the unobserved infection times. The infection and removal times
# are implemented as *individual-based* event histories. This affects 
# the choice of the form of the likelihood: see the lecture notes 
# and the instructions to this exercise.
#
# INPUT  remtimes  = the removal times (in ascending order)
#        M         = the total number of individuals in the community
#        mcmc.size = the number of MCMC iterations to be performed
#
# OUTPUT sampleSIR = a (mcmc.size)x2 matrix, the columns including the MCMC 
#                    samples of parameters 'beta' and 'gamma'
#
# This function calls subroutines initializedata.R, update_beta.R, 
# update_gamma.R, and update_inftimes.R.

# (a) Reserve space for the MCMC output
mcmc.samples           = matrix(0,mcmc.size,2); 
mcmc.samples           = as.data.frame(mcmc.samples)
colnames(mcmc.samples) = c("beta","gamma")

# (b) Initialize the model unknowns
bet          = 0.01; # infection rate
gam          = 0.01; # recovery rate
completedata = initializedata(remtimes) # initialize the complete data matrix
                                        # with 2 columns: infection times and
                                        # removal times

# Store the initial values
mcmc.samples[1,]$beta  = bet
mcmc.samples[1,]$gamma = gam

# The MCMC iterations
for (i in 2:mcmc.size){

   # Update the parameters and the unknown infection times 
   bet          = update_beta(completedata,M);            # (c) parameter beta
   gam          = update_gamma(completedata);             # (d) parameter gamma
   completedata = update_inftimes(completedata,M,bet,gam) # (e) unobserved infection times

   # Store the current iterates for 'beta' and 'gamma'
   mcmc.samples[i,]$beta  = bet
   mcmc.samples[i,]$gamma = gam

   }

# Return the MCMC sample of the two model parameters
return(mcmc.samples)

}



#######################################################
#######################################################

# Read in the data and draw the MCMC samples
remtimes    = readdata()
mcmc.size   = 600
mcmc.sample = sampleSIR(remtimes,M=120,mcmc.size)

# Plot the sample paths
par(mfrow=c(1,2))
plot(mcmc.sample$beta,type="l",xlab="iteration",ylab="beta")
plot(mcmc.sample$gam,type="l",xlab="iteration",ylab="gamma")

# Plot the histograms of the posterior distributions after 
# discarding a number of burn-in samples
burnin = 200
retain.sample = mcmc.sample[burnin:mcmc.size,]
hist(retain.sample$beta,xlab="beta (per day)",main=" ")
hist(retain.sample$gam,xlab="gamma (per day)",main=" ")

# Plot the posterior distribution of the 'reproduction number'
par(mfrow=c(1,1))
hist(retain.sample$bet/retain.sample$gam,main="")

