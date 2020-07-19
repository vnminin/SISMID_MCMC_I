## This script illustrates the ergodic theorem using the Ehrenfest model of diffusion
## Author: Vladimir N. Minin
## last update: 07/18/2020

#' Generate the next state of the Ehrenfest diffusion model
#' 
#' @param cur_state Current state of the Ehrenfest diffusion model
#' @param num_mol Total number of moleculers in the Ehrenfest diffusion model
#'
#' @return an integer that is either (cur_state+1) or (cur_state-1) 
#'
#' @examples
#' next_state(32, 1000)
next_state = function(cur.state, num.mol){
  return.value=NULL
  
  ## this function randomly draws a new state of the Ehrenfest model
  ## One of your tasks is to finish writing this function
  return(return.value)
}

## set the number of molecules and the number of iterations
my_num_mol = 100
sim_size = 10000

## initialize the chain by drawing the initial state uniformly at random from all possible states. R function `sample()' will be handy.
my_draws = numeric(sim_size)
#my_draws[1] = FINISH AND UNCOMMENT THIS LINE

## run the Markov chain

for (i in 2:sim_size){
  ## use next.state function to evolve the Markov chain one step at a time
  
  #my_draws[i] = FINISH AND UNCOMMENT THIS LINE
}

## plot the chain

plot(1:sim_size, my_draws, type="l", ylab="Ehhenfest State", xlab="Time Step")

## get the "time averages"

mean(my_draws)
var(my_draws)

## get the "space averages," theoretical binomial 
## mean and variance
my_num_mol/2
my_num_mol/4
