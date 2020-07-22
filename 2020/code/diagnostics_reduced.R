## This script illustrates some diagnostic tools available in package coda
## Author: Vladimir N. Minin
## last update: 07/19/20

## first we need to load coda and mcmcse packages (you need to install them first)
library(coda)
library(mcmcse)

## now let's load our M-H and Gibbs sampler examples
source("https://raw.githubusercontent.com/vnminin/SISMID_MCMC_I/master/2016/code/chainGibbs.R")
dev.off()

## run one M-H example chain
gibbs_chain1 = chainGibbs(5000, 1,1)

## convert the output into coda format
coda_gibbs_chain = mcmc(cbind(gibbs_chain1$q[101:5000],gibbs_chain1$n111[101:5000]))

summary(coda_gibbs_chain)

## we can also compute Monte Carlo error for the quantiles
mcse.q.mat(coda_gibbs_chain, 0.025)
mcse.q.mat(coda_gibbs_chain, 0.975)

## plot traceplots and (hopefully) posterior densities
plot(coda_gibbs_chain)

## look at what coda can do
help(package=coda)

## look at the menu options
##codamenu()

## all commands are availabe outside of the menu

## plot autocorrelations plots
autocorr.plot(coda_gibbs_chain)

## calculate effective sample size
effectiveSize(coda_gibbs_chain)


## Run 50 chains with overdispersed starting values
coda_gibbs_chains = list()

for (i in 1:50){
  gibbs_chain = chainGibbsUserStart(1000,1,1,sample(0:275,size=1))
  coda_gibbs_chains[[i]] = mcmc(cbind(gibbs_chain$q,gibbs_chain$n111))
}

coda_gibbs_list = mcmc.list(coda_gibbs_chains)

plot(coda_gibbs_list)

## compute and plot Gelman-Rubin potential reduction factor
gelman.diag(coda_gibbs_list)

gelman.plot(coda_gibbs_list)


## Exercise: perform diagnostics for SIR data augmentatkon MCMC

source("https://raw.githubusercontent.com/vnminin/SISMID_MCMC_I/master/2019/code/SIRaugmentation.R")
