## This script illustrates the ergodic theorem using the Ehrenfest model of diffusion
## Author: Vladimir N. Minin
## last update: 07/12/2015

## this function randomly draws a new state of the Ehrenfest model

next_state = function(cur_state, num_mol){
  return.value=NULL
  
  if (runif(1)<cur_state/num_mol){
    return_value = cur_state-1
  }else{
    return_value = cur_state+1
  }
  
  return(return_value)
}

## set the number of molecules and the number of iterations
my_num_mol = 100
sim_size = 10000

## initialize the chain by drawing the initial state uniformly at random from all possible states
my_draws = numeric(sim_size)
my_draws[1] = 100#sample(c(0:my_num_mol), size=1)

## run the Markov chain

for (i in 2:sim_size){
  my_draws[i] = next_state(my_draws[i-1], my_num_mol)
}

## plot the chain

plot(1:sim_size, my_draws, type="l", ylab="Ehhenfest State", xlab="Time Step")

## get the "time averages"

mean(my_draws)
var(my_draws)

## get the "space averages"
my_num_mol/2
my_num_mol/4
