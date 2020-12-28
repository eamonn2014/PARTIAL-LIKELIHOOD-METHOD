
events <- function(alpha=.05, beta=0.1, hr=1.2) {

A <- alpha # alpha
B <- beta  # beta
b <- hr

f<- ((qnorm(1-A/2))+qnorm(1-B))^2
d <- (4*f)/ (log(b)^2)
 
# N <- d/(1-(s1+s2)/2)
# N # total sample
return(d)
}


ev <- ceiling(events())

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

coxdataf <- function(n, allocation=.5, hr=1.2, baseline=.4, followup=50) { 
  
  #n=1000; allocation =.5; hr=2; baseline=.4
  
  trt <- sample(0:1, n,  rep=TRUE, prob=c(1-allocation, allocation))
  
  cens <- 15*runif(n)
  
  h <- baseline*exp(log(hr)*(trt==1))  # hazard function h(t)
  
  dt <- -log(runif(n))/h
  
  label(dt) <- 'Follow-up Time'
  
  e1 <- ifelse(dt <= cens,1,0)
  
  e <- ifelse(e1 <= followup,1,0)
  
  dt <- pmin(dt, cens)
  
  units(dt) <- "Year"
  
  d <<- data.frame(cbind(dt,e,trt=trt))  ##why the << required to circumvent error?
  
  dd <<- datadist(d)
  options(datadist='dd')
  
  foo <-d
  # S <- Surv(dt,e)
  f <- cph(Surv(dt,e) ~  trt, x=TRUE, y=TRUE, data=d )
  teststat<--2*(as.numeric(f$loglik)[1]-as.numeric(f$loglik[2])) #  null-model
  
  p <- pchisq(teststat,df=df,lower.tail=FALSE)
  
  return(p)
  
  #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  
  return(list(p =p, d=d))
  
}


dummy <- coxdata(n=ev, allocation =.5, hr=1.2, baseline=.4, foolup)

mean(replicate(1000, coxdata(n=ev, allocation =.5, hr=1.2, baseline=.4))<0.05)


powerSurvEpi:: powerCT.default0(k = 1, m = 1265, RR = 01.2, alpha = 0.05)
powerSurvEpi:: ssizeCT.default(.9, k=1, pE=.65, pC=.7, RR, alpha = 0.05)



n=800; lambdaT=14.4; lambdaC=12; beta1=1.2; per=0.7 ;

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# https://stattools.crab.org/R/Help%20Documents/Survival_Converter_HelpDoc.html
su <- function(find, median_survival, survival_time_t, hazard_rate, time) {
  
  result = list(median = NA,
                survival = NA,
                hr = NA)
  
  if (find == "Median") {
    result$median = -(log(0.5)) / hazard_rate
    result$median = round(result$median, digits = 2)
  } else if (find == "Survival") {
    result$survival = exp(-time * hazard_rate)
    result$survival = round(result$survival, digits = 2)
  } else if (find == "Hazard_Rate") {
    result$hr = -(log(survival_time_t)) / time
    result$hr = round(result$hr, digits = 3)
  } else if (find == "Hazard_Rate_Median") {
    result$hr = -(log(0.5)) / median_survival
    result$hr = round(result$hr, digits = 3)
  }
  
  return( (result ))
}

su("Median", hazard_rate=14)
su("Survival", time=5,hazard_rate= 1/14)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~start here

# given hr, time  and haz rate provide survival and survival based on hr at time
hr. <- 1.2 
fup.<- 20
hazard_rate <- 1/10   #  number of events per unit time
time <- fup.    
s1.  <- exp(-time * hazard_rate)   
s2.  <- s1.^hr.
s1.
s2.

# s1.=.04
# s2.=.04^.6


pE = 1 - (s1. + s2.) / 2  # stata

 


# s1.% of patients survive for fup. months.
# With a 20% increase in the hazard of the experimental group, we expect roughly 100% x
# s2.to survive for fup. months. This means that, to observe all patients die, 
# we would need to continue our study for more than fup. months. 

pE = 1 - (s1. + s2.) / 2  # stata
# -log(per2) *lambdaT # time
#overall probability of a subject failing in a study at end of study T

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
powerSurvEpi:: ssizeCT.default(power=.9, k=1, pE=1-s1., pC=1-s2., RR=1.2, alpha = 0.05)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

rosner <- function(  s1, s2, alpha, beta, k, hr  ) {
  
  f<- ((qnorm(1-alpha/2))+qnorm(1-beta))^2
  
  ev <- 1/k*( (k*hr+1) / (hr-1)  )^2   * f 
  
  n <- 2* (ev /  (s1 + s2))
  
  return(list(n=n, ev=ev))
  
}

nn<-rosner(s1= 1-s1., s2 = 1-s2., hr=hr., k=1, alpha=.05, beta=.1)   #s1 and
nn



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ENTER 1/HAZARD TO RWEIBULL!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



# this should give me same power as rosner. 

spow <- function(n=nn/pE, lambdaT= 0, lambdaC=1/10, beta1= hr., fup= fup.) {  # need to enter 1/hazard here
  
  x1 = sample(0:1, n,replace=TRUE)#
  
  beta <- log(beta1)
  
  # true event time
  T = rweibull(n, shape=1, scale=lambdaT)*exp(-beta*x1)
  C = rweibull(n, shape=1, scale=lambdaC)                  # censoring time
  time = pmin(T,C)                                         # observed time is min of censored and true
  event = 1*(time==T )  #    event = time==T
  
  # censor at end of follow up #NEWWWWWWWWWWWWWWWWWWWWWWW
 event <- ifelse((event==1 & time>fup) ,0, event)         # set to 1 if event is observed followup over fup
  time <-  ifelse( time>fup, fup, time)
  
  f1 <- coxph(Surv(time, event)~ x1 , method="breslow")
  p <- anova(f1)$`Pr(>|Chi|)`[2] # P-value from likelihood ratio test
  
  #df <- length(coef(f1)) 
  
  #teststat <- -2*(  as.numeric(f1$loglik)[1] - as.numeric(f1$loglik[2])  ) # null-model
  
  #p <- pchisq(teststat, df=df, lower.tail=FALSE)
  
  return(p)
}

# -log(per2) *lambdaT # time
 

 
#mean(replicate(1000, spow(n=2000,  lambdaT=hazard_rate,  lambdaC=17.6,   beta1=hr.,  fup=fup.))    <0.05) #use 1/14.4?
mean(replicate(1000, spow(n=nn[1][[1]]/pE,    lambdaT=1/hazard_rate,  lambdaC=10,   beta1=hr.,  fup=fup.))    <0.05) #use 1/14.4?
#mean(replicate(1000, spow(n=1000,  lambdaT=hazard_rate,  lambdaC=17.6,   beta1=hr.,  fup=fup.))    <0.05) 
#mean(replicate(1000, spow(n=800,   lambdaT=hazard_rate,  lambdaC=17.6,   beta1=hr.,  fup=fup.))    <0.05) 
#mean(replicate(1000, spow(n=600,   lambdaT=hazard_rate,  lambdaC=17.6,   beta1=hr.,  fup=fup.))    <0.05) 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

n=nn[1][[1]]/pE; lambdaT=1/hazard_rate;  lambdaC=10;   beta1=1.2;  fup=fup.

x1 = sample(0:1, n,replace=TRUE)#

beta <- log(beta1)

# true event time
T = rweibull(n, shape=1, scale=lambdaT)*exp(-beta*x1)
C = rweibull(n, shape=1, scale=lambdaC)                  # censoring time
time = pmin(T,C)                                         # observed time is min of censored and true
event = 1*(time==T )  #    event = time==T

summary(time)
summary(event)

# censor at end of follow up #NEWWWWWWWWWWWWWWWWWWWWWWW
event <- ifelse((event==1 & time>fup) ,0, event)         # set to 1 if event is observed followup over fup
time <-  ifelse( time>fup, fup, time)

summary(time)
summary(event)

f1 <- coxph(Surv(time, event)~ x1 , method="breslow")
f1
survfit <- survfit(Surv(time,event) ~ x1)
plot(survfit, ylab="Survival probability", xlab="Time", col=c('blue','red'))























hazard_rate <- 1/14.4
time=20
s1. <- exp(-time * hazard_rate)
s2. <- s1.^1.2
s1.
s2.

rosner <- function(  s1, s2, alpha, beta, k, hr  ) {
  
  f<- ((qnorm(1-alpha/2))+qnorm(1-beta))^2
  
  ev <- 1/k*( (k*hr+1) / (hr-1)  )^2   * f 
  
  n <- 2* (ev /  (s1 + s2))
  
  return(n)
  
}

rosner(s1= 1-s1., s2 = 1-s2., hr=1.2, k=1, alpha=.05, beta=.1)   #s1 and





#s1= .3707; s2 = .4890; hr=.7; k=1; alpha=.05; beta=.2

rosner <- function(  s1, s2, alpha, beta, k, hr  ) {
  
  f<- ((qnorm(1-alpha/2))+qnorm(1-beta))^2
  
  ev <- 1/k*( (k*hr+1) / (hr-1)  )^2   * f 
  
  n <- 2* (ev /  (s1 + s2))
  
  return(n)
  
}


#p838 fundamentals of biostatistics
#rosner(s1= .3707, s2 = .4890, hr=.7, k=1, alpha=.05, beta=.2)


rosner(s1= .7, s2 = .65, hr=1.2, k=1, alpha=.05, beta=.2)   #s1 and s2 


# https://stattools.crab.org/R/Help%20Documents/Survival_Converter_HelpDoc.html#statistical-code
# lets look at the curves, to get totals
hazard_rate <- 1/14.4
time=50
s1. <- exp(-time * hazard_rate)
s2. <- s1.^1.2
s1.
s2.


rosner(s1= 1-s1., s2 = 1-s2., hr=1.2, k=1, alpha=.05, beta=.2)   

# need more if time of trial less
time=40
s1. <- exp(-time * hazard_rate)
s2. <- s1.^1.2
rosner(s1= 1-s1., s2 = 1-s2., hr=1.2, k=1, alpha=.05, beta=.2)   


rosner(s1= 1-s1., s2 = 1-s2., hr=1.2, k=1, alpha=.05, beta=.2)   

# need more if time of trial less
time=20
s1. <- exp(-time * hazard_rate)
s2. <- s1.^1.2
rosner(s1= 1-s1., s2 = 1-s2., hr=1.2, k=1, alpha=.05, beta=.2)   




