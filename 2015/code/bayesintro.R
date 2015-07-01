## R program to compute the posterior distribution of the transmission probability
## Beta prior distribution for the binomial likelihood
## Author: Betz Halloran
## Last updated June 28, 2015

#Generate a sequence from 0 to 1 in increments of .01
x=seq(0,1,.01)

x

#Check the length of the vector x

length(x)

#Number of trials (number exposed)
n=c(5,20,50,1000)

n

#Number of successes (infections)
y=c(2,8,20,400)

y

# Check the means of the four trials
y/n

#Set up noninformative Beta prior information
my.alpha=1
my.beta=1
#Set up a matrix to hold results 
posta=matrix(0,nrow=4, ncol= length(x))

posta

# plot the posterior densities using different amounts of data
#open pdf (or ps)  file graphics device

#pdf(file="/Users/betz/Documents/TexWork/MCMC/betz/Bayesintro/betaunif1.pdf", height=6.5, width=8.9)
#set up to have 4 plots in one figure 2 by 2
par(mfrow=c(2,2))
# loop through 4 graphs
for (i in 1:4){
posta[i,] =dbeta(x,my.alpha+y[i],my.beta+n[i]-y[i])
plot(x,posta[i,], type="l", ylab="Posterior Density ", xlab="p ")
}
# close graphics device if using pdf (or ps)
#dev.off()

# return to 1  plot if need be 
par(mfrow=c(1,1))

#Posterior summaries  using closed form distribution

#prior mean
priormean=my.alpha/(my.alpha+my.beta)
priormean

# posterior mean 

postmean= (my.alpha+y)/(my.alpha+my.beta+n)
postmean

#posterior median
# use qbeta to get the values at the given quantiles
postmedian=qbeta(0.5, my.alpha+y, my.beta+n-y)
postmedian

#median, 95%  equitailed posterior or credible interval 
# use qbeta to get the values at the given quantiles
# set up matrix
post95int=matrix(0,4,3)
for (i in 1:4){
 post95int[i,] = qbeta(c(0.5, 0.025,0.975), my.alpha+y[i], my.beta+n[i]-y[i])
}
post95int

round(post95int,3)

# Drawing random samples from the posterior distributions to imitate results 
# of an MCMC output

#Set the number of samples
nsamp=5000
#Set up matrix
post=matrix(0,4,nsamp)

##pdf(file="/Users/betz/Documents/TexWork/MCMC/betz/Bayesintro/betaunif1r.pdf", height=6.5, width=8.9)

par(mfrow=c(2,2))
for (i in 1:4){
post[i,]=rbeta(nsamp,my.alpha+y[i],my.beta+n[i]-y[i])
hist(post[i,], xlim=c(0,1), xlab = "p ", ylab = "Frequency")
}

#dev.off()
par(mfrow=c(1,1))

#Posterior summaries using random samples

# posterior mean

postmeanr= apply(post,1,mean)
postmeanr

#Compare with analytic posterior mean 

postmean

# Get summary of simulated posterior distributions
# Row by row

summary(post[1,])
summary(post[2,])
summary(post[3,])
summary(post[4,])

# Get all summaries of all four rows at the same time

apply(post, 1, summary)


## Now use informative conjugate priors

alpha1=c(1,2,50,3,12,60,4,16)
beta1=c(1,2,50,2,8,40,1,4)

priorsum=alpha1+beta1
priormean=alpha1/(alpha1+beta1)
priorsum
priormean

n1=50
y1=20

post95int2=matrix(0,length(alpha1),3)
for (i in 1:length(alpha1)){
 post95int2[i,] = qbeta(c(0.5,0.025,0.975), alpha1[i]+y1, beta1[i]+n1-y1)
}
round(post95int2,3)

par(mfrow=c(1,1))


## Here is a simple function to compute the posterior distribution 
## of the transmission probability. It plots the posterior 
## distribution and returns the posterior median and
## 95% equitailed credible interval
## n = number of observations (exposed),
## y = number of successes (infections)
## alpha, beta, parameters of the prior Beta distribution

mybeta=function(n,y,alpha,beta){
          x=seq(0,1,.01)
          postbeta =dbeta(x,alpha+y,beta+n-y)
          plot(x,postbeta, type="l", ylab="Posterior Density ", xlab="p ")
          mybeta=qbeta(c(0.5, 0.025,0.975), alpha+y, beta+n-y)
}

test1=mybeta(50,20,1,1)
test1

### Alternatively, you can return the values at the end.

mybeta2=function(n,y,alpha,beta){
          x=seq(0,1,.01)
          postbeta=dbeta(x,alpha+y,beta+n-y)
          plot(x,postbeta, type="l", ylab="Posterior Density ", xlab="p ")
          return(qbeta(c(0.5, 0.025,0.975), alpha+y, beta+n-y))
}

test2=mybeta2(50,20,1,1)
test2

### Alternatively, you can return a list at the end.

mybeta3=function(n,y,alpha,beta){
          x=seq(0,1,.01)
          postbeta=dbeta(x,alpha+y,beta+n-y)
          plot(x,postbeta, type="l", ylab="Posterior Density ", xlab="p ")
          mybeta=list(answer=qbeta(c(0.5, 0.025,0.975), alpha+y, beta+n-y))
}

test3=mybeta3(50,20,1,1)
test3
test3$answer