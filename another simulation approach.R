# Survival3Notes in survival simulation paper folder on one drive slide 20

library(survival)
lambda1=0.058
lambda2=0.116
n=80
timeev1=rexp(n,lambda1)
timeev2=rexp(n,lambda2)
timeev=c(timeev1,timeev2)
timecensor=runif((2*n),1,10)
time=apply(cbind(timeev,timecensor),1,min)
stat=(time==timeev)+0
x=c(rep(0,n),rep(1,n))
fitcox=coxph(Surv(time,stat)~x)
fitkm=survfit(Surv(time,stat)~x)
hr=exp(fitcox$coef)
est=summary(fitkm,times=5)$surv
pvwald=summary(fitcox)$waldtest[3]
pvlrt=summary(fitcox)$logtest[3]
pvscore=summary(fitcox)$sctest[3]
out=data.frame(hr,estlowrisk=est[1],esthighrisk=est[2],pvwald,pvlrt,pvscore)
out

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library(survival)
fsim=function(i,n,lambda1,lambda2,f,a)
{
  timeev1=rexp(n,lambda1)
  timeev2=rexp(n,lambda2)
  timeev=c(timeev1,timeev2)
 
  timeev=c(timeev1,timeev2)
  timecensor=runif((2*n),f,a)
  time=apply(cbind(timeev,timecensor),1,min)
  stat=(time==timeev)+0
  x=c(rep(0,n),rep(1,n))
  fitcox=coxph(Surv(time,stat)~x)
  fitkm=survfit(Surv(time,stat)~x)
  hr=exp(fitcox$coef)
  est=summary(fitkm,times=5)$surv
  pvwald=summary(fitcox)$waldtest[3]
  pvlrt=summary(fitcox)$logtest[3]
  pvscore=summary(fitcox)$sctest[3]
  
  numev<-sum(fitkm$n.event)
  
  out=data.frame(hr,numev,
                 estlowrisk=est[1],
                 esthighrisk=est[2],
                 pvwald,pvlrt,pvscore)
  print(i)
  return(out)
}

fsim(i=1, n=80, lambda1=0.058, lambda2= 0.116, f=1, a=9) 

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


set.seed(123)
nsim=1000 ## try it first with 10
a=sapply(c(1:nsim),fsim,n=80,lambda1=0.058, lambda2=0.116,f=1,a=10)
#a # when you try it with 10, check to see how it looks
bcr=data.frame(apply(a,1,unlist))
head(bcr)


sum(bcr$pvwald<=0.05)/nsim
 
sum(bcr$pvlrt<=0.05)/nsim

sum(bcr$pvscore<=0.05)/nsim
 
mean(bcr$hr)
 
mean(bcr$estlowrisk)
 
mean(bcr$esthighrisk)
 
mean(bcr$numev)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


 
# The accrual should not be longer than 6 years.
# The added follow-up at the end of accrual can be as long as 2 years but no longer.
# A difference of 10% at 5 years is considered clinically important.
# The accrual rate is about 50/year.
alpha=0.05
beta=0.2
zalpha=qnorm(1-alpha/2)
zbeta=qnorm(1-beta)
sigma=1/2
S0=0.65
S1=0.75
HR=log(S0)/log(S1)
HR
 
HR=1.5
sqrtnev=(zalpha+zbeta)/(sigma*log(HR))
sqrtnev^2

ac=5
fup=2#3
lambda0=0.058
lambda1=0.0856255  # guess work
Pev0=1-(exp(-lambda0*fup)-exp(-lambda0*(ac+fup)))/(lambda0*ac)
Pev1=1-(exp(-lambda1*fup)-exp(-lambda1*(ac+fup)))/(lambda1*ac)
Pev=mean(c(Pev0,Pev1))
N=191/Pev  #705.3517
N


##########################









