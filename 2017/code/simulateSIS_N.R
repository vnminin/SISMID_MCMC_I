simulateSIS_N = function(N,la,mu,initprob,T){
# July 12, 2014/KA
#
# This function samples realizations of a binary Markov
# process, i.e., it simulates complete data in an alternating Poisson process.
# for a fiven number of independent individuals. The processes are simulated 
# from time 0 to T, at which they are censored.
#
# INPUT N        = the number of individuals
#       la       = the rate (force) of infection per time unit, i.e., 
#                  the hazard of transition 0 -> 1
#       mu       = the rate of clearing infection per time unit, i.e.,
#                  the hazard of  transition 1 -> 0
#       initprob = the probability of being infected at time 0
#       T        = the maximum time, i.e., the time of censoring the 
#                  individual processes
#
# OUTPUT simulateSIS_N = The simulated complete data as a list of length N. 
#                        Each of the N arguments is a list of two arguments: 
#                        list(t,state) where t = the event times and 
#                        state = the corresponding  epidemiological states 
#                        (0 = susceptible, or 1 = infected) for one individual.
#
# This function calls subroutine simulateSIS.R to draw a realization for
# one individual at a time.

# Initialize the list structure
data = list(1:N)

# Simulate data for N independent individuals
for (i in 1:N){

  # Simulate data for individual 'i'
  dataind   = simulateSIS(la,mu,initprob,T)

  # Store in the output list
  data[[i]] = dataind

  }

# The output
simulateSIS_N = data

}


##########################################
simulateSIS = function(la,mu,initprob,T){
# June 11, 2010/KA
#
# This function samples one realization of a binary Markov
# process, i.e., it samples complete data for one (independent)
# individual in an alternating Poisson process. SIS model (without transmission). 
# The process is simulated from time 0 to time T at which it is censored.
#
# INPUT la       = the rate (force) of infection (per time unit)
#       mu       = the rate of clearing infection (per time unit)
#       initprob = the probability of being infected at time 0
#       T        = the maximum time, i.e., the time of censoring the process
#
# OUTPUT simulateSIS = a list with two arguments:
#                      [[1]] = the event times (times of transitions), 
#                              including times 0 and T
#                      [[2]] = the epidemiological states after the event times: 
#                              an alternating vector of 0's (susceptible) and 1's (infected)

# Reserve space (max. 100 event times assumed)
t     = rep(-1,100) # event times
state = rep(-1,100) # states (0 = susceptible, or 1 = infected)

# Sample the first state (at time 0)
t[1]     = 0                       # time
state[1] = runif(1,0,1) < initprob # state

# Simulate the process until time T
k = 1
while (t[k] < T){

 if (state[k] == 0){
   t[k+1]     = t[k] + rexp(1,la); # the time of infection
   state[k+1] = 1
   }
 else {
   t[k+1]     = t[k] + rexp(1,mu); # the time of clearance
   state[k+1] = 0
   }

 #
 k = k+1
 }

# The censoring time
t[k]     = T   # time
state[k] = -1  # no more transitions

# Retain only the event times and the states (i.e. remove -1's).
# The length of the 'time vector' is the number of episodes + 1.
# The length of the 'state vector' is the number of episodes.
t     = t[t>-1]
state = state[state>-1]

# The output 
simulateSIS = list(t,state)

}

complete_data = simulateSIS_N(N=100,la=0.45,mu=0.67,initprob=0.40,T=12)

