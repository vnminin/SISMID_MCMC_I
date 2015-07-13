## This script illustrates the ergodic theorem using the Ehrenfest model of diffusion
## Author: Vladimir N. Minin
## last update: 07/12/2015

## this function randomly draws a new state of the Ehrenfest model

next.state = function(cur.state, num.mol){
  return.value=NULL
  
  if (runif(1)<cur.state/num.mol){
    return.value = cur.state-1
  }else{
    return.value = cur.state+1
  }
  
  return(return.value)
}

## set the number of molecules and the number of iterations
my.num.mol = 100
sim.size = 1000

## initialize the chain by drawing the initial state uniformly at random from all possible states
my.draws = numeric(sim.size)
my.draws[1] = sample(c(0:my.num.mol), size=1)

## run the Markov chain

for (i in 2:sim.size){
  my.draws[i] = next.state(my.draws[i-1], my.num.mol)
}

## plot the chain

plot(1:sim.size, my.draws, type="l", ylab="Ehhenfest State", xlab="Time Step")

## get the "time averages"

mean(my.draws)
var(my.draws)

## get the "space averages"
my.num.mol/2
my.num.mol/4
