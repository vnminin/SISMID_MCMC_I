## This script illustrates some diagnostic tools available in package coda
## Author: Vladimir N. Minin
## last update: 07/11/18

## first we need to load coda and mcmcse packages (you need to install them first)
library(coda)
library(mcmcse)

## now let's load our M-H and Gibbs sampler examples
source("https://raw.githubusercontent.com/vnminin/SISMID_MCMC_I/master/2016/code/chainGibbs.R")
dev.off()

## run one M-H example chain
gibbs.chain1 = chainGibbs(5000, 1,1)

## convert the output into coda format
coda.gibbs.chain = mcmc(cbind(gibbs.chain1$q[101:5000],gibbs.chain1$n111[101:5000]))

summary(coda.gibbs.chain)

## we can also compute Monte Carlo error for the quantiles
mcse.q.mat(coda.gibbs.chain, 0.025)
mcse.q.mat(coda.gibbs.chain, 0.975)

## plot traceplots and (hopefully) posterior densities
plot(coda.gibbs.chain)

## look at what coda can do
help(package=coda)

## look at the menu options
##codamenu()

## all commands are availabe outside of the menu

## plot autocorrelations plots
autocorr.plot(coda.gibbs.chain)

## calculate effective sample size
effectiveSize(coda.gibbs.chain)


## Run 50 chains with overdispersed starting values
coda.gibbs.chains = list()

for (i in 1:50){
  gibbs.chain = chainGibbsUserStart(1000,1,1,sample(0:275,size=1))
  coda.gibbs.chains[[i]] = mcmc(cbind(gibbs.chain$q,gibbs.chain$n111))
}

coda.gibbs.list = mcmc.list(coda.gibbs.chains)

plot(coda.gibbs.list)

## compute and plot Gelman-Rubin potential reduction factor
gelman.diag(coda.gibbs.list)

gelman.plot(coda.gibbs.list)
