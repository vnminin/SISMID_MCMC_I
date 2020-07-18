##Introduction to R and Bayes programming 
## Last updated July 12, 2020 / MEH
## SISMID 2020 MCMC I

## The R commands can be cut and paste from this document.

##Simple Beta posterior distribution of the transmission probability 

## R program to compute the posterior distribution of the transmission probability
## Beta prior distribution of the binomial likelihood 

## We want to evaluate the density of the posterior of p along the interval [0,1] 
## To start, generate a sequence from 0 to 1 in increments of .01 that will supply
## the values where we will evaluate the posterior density 

x = seq(0,1, by = .01)
x

## PART 1: Noninformative Beta prior
## Observed data 
## Generate a vector of the observed number of trials in the four experiments 
n=c(5,20,50,1000)
n

## Generate a vector of the number of successes (infections) in the four experiments
y=c(2,8,20,400)
y

## Observed fractions of successs are the same
y/n

##Set up noninformative Beta prior distributions 
my.alpha = 1
my.beta = 1 
my.alpha
my.beta

##Set up a matrix with 4 rows and the number of columns that is the length of the 
## x vector where the values of the posterior densities  will be evaluated. This 
## matrix will hold the values for the four posterior densities.
## The value 0 is a place holder. Other entries could be used.

posta = matrix(0, nrow=4, ncol = length(x))

##plot the four posterior densities using different amounts of data from 
## the four experiments with the noninformative Beta(1,1) prior distribution

## open pdf (or ps) file graphics device 
#pdf(file="/Users/betz/Documents/Bayesintro/betaunif1.pdf", height=6.5, width = 8.9)
## set up to have 4 plots in one figure with 2 by 2
 
par(mfrow=c(2,2))
# loop through the for graphs. Use a for loop 

for (i in 1:4){
  posta[i,] = dbeta(x,my.alpha+y[i],my.beta+n[i]-y[i])
  plot(x,posta[i,], type = "l", ylab ="Posterior Density", xlab="p")
}

##
## With little data (upper left) the posterior distribution is quite spread out, reflecting
## the uncertainty. As the amount of data increases, the posterior distribution becomes more 
## peaked around the data estimate of 0.4.
## Note also that the y-axis changes in the four graphs. The densities all have the same area. 
##
## close graphics device if using pdf (or ps)
#dev.off() 
#return to 1 plot if need be 
par(mfrow=c(1,1))

## Summaries of the posterior distribution 

## Posterior summaries using closed form distributions

#prior mean 
priormean = my.alpha/(my.alpha + my.beta)
priormean

## The prior mean of the noninformative prior is 0.5, the mid-point of the interval 0 to 1. 

# posterior mean 

postmean= (my.alpha+y)/(my.alpha+my.beta+n)
postmean
round(postmean,4)

## Here we have used the closed form for the Beta distribution mean to get the posterior
## mean. As the amount of data increases, the posterior mean approacheds the mean of the 
## which is 0.40. 

#posterior median
# use qbeta to get the values at the given quantiles of the posterior distribution
postmedian=qbeta(0.5, my.alpha+y, my.beta+n-y)
round(postmedian,4)

## The posterior medians also grow closer to 0.4 as the amount of data increases

# Now we want to get an interval around the median for the posterior distributions
# median, 95%  equitailed posterior or credible interval 
# use qbeta to get the values at the given quantiles 0.5, 0.025, 0.975
# set up matrix
post95int=matrix(0,4,3)
for (i in 1:4){
  post95int[i,] = qbeta(c(0.5, 0.025,0.975), my.alpha+y[i], my.beta+n[i]-y[i])
}
round(post95int,3)

## Again one sees quantitatively that the 95% posterior credible interval narrows with
## data


# Drawing random samples from the posterior distributions to imitate results 
# of an MCMC output; use rbeta() command for random samples from a beta distribution

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

## The histograms with random draws from the four posterior Beta distributions look very similar ## to the histograms using the closed form for the Beta distribution. 
## These histograms are an APPROXIMATION to the target posterior Beta distribution. 

#dev.off()
par(mfrow=c(1,1))

#Posterior summaries using random samples 
# We can get summaries of the approximations to the Beta posterior distributions using the random samples drawn from the Beta posterior distributions.

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

## END PART 1. 

## PART 2: Informative Beta priors
## Now use informative conjugate priors
## We look at how different amounts of information in the conjugate Beta prior influences 
## the posterior distribution

alpha1=c(1,2,50,3,12,60,4,16)
beta1=c(1,2,50,2,8,40,1,4)
alpha1
beta1

## The prior sum is the total amount of data in the prior, the prior mean in the mean of the
## Beta prior 

priorsum=alpha1+beta1
priormean=alpha1/(alpha1+beta1)
priorsum
priormean

# The prior sum ranges from 2 (noninformative Beta prior) to 100
# The prior mean ranges from 0.5 to 0.8. 

# The data is going to be held constant while examining the influence of changing 
# the prior information. The data mean is 0.4.

n1=50
y1=20

## Now generate the median and 95% posterior interval for the 8 different 
## posterior distributions

post95int2=matrix(0,length(alpha1),3)
for (i in 1:length(alpha1)){
 post95int2[i,] = qbeta(c(0.5,0.025,0.975), alpha1[i]+y1, beta1[i]+n1-y1)
}
round(post95int2,3)

## Examine how the medians and width of the intervals depends both on the prior sums and the prior means. 

par(mfrow=c(1,1))

## END PART 2. 

## If you have not written functions in R before, the following will hopefully
## be usefule to you. If you know how to write an R function, you can skip this. 
## Here is a simple function to compute the posterior distribution 
## of the transmission probability. It plots the posterior 
## distribution and returns the posterior median and
## 95% equitailed credible interval
## n = number of observations (exposed),
## y = number of successes (infections)
## alpha, beta, parameters of the prior Beta distribution


### Return the values at the end.

mybeta2=function(n,y,alpha,beta){
          x=seq(0,1,by=.01)
          postbeta=dbeta(x,alpha+y,beta+n-y)
          plot(x,postbeta, type="l", ylab="Posterior Density ", xlab="p ")
          mybeta=qbeta(c(0.5, 0.025,0.975), alpha+y, beta+n-y)
          return(mybeta)
}

test2=mybeta2(50,20,1,1)
test2

### Alternatively, you can return a list at the end.

mybeta3=function(n,y,alpha,beta){
          x=seq(0,1,by=.01)
          postbeta=dbeta(x,alpha+y,beta+n-y)
          plot(x,postbeta, type="l", ylab="Posterior Density ", xlab="p ")
          mybeta=list(answer=qbeta(c(0.5, 0.025,0.975), alpha+y, beta+n-y))
          return(mybeta)
}

test3=mybeta3(50,20,1,1)
test3
test3$answer

### Alternatively, you can return a list with more than one object at the end.

mybeta4=function(n,y,alpha,beta){
  x=seq(0,1,by =.01)
  postbeta=dbeta(x,alpha+y,beta+n-y)
  plot(x,postbeta, type="l", ylab="Posterior Density ", xlab="p ")
  answer1=qbeta(.5,alpha+y, beta+n-y)
  answer2=qbeta(.025,alpha+y, beta+n-y)
  answer3=qbeta(.975,alpha+y, beta+n-y)
  mybeta=list(answer1,answer2,answer3)
}

test4=mybeta4(50,20,1,1)
print(test4)
##extract the first element in the list test4
test4[[1]]

