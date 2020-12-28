 #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 # median = log2/ hazard
 #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  rm(list=ls())
  library(survival)
  
  set.seed(876)
  n = 100
  beta1 = 0       # NO TRT EFFECT
  lambdaT = .002  # baseline hazard
  lambdaC = .004  # hazard of censoring
  x1 = sample(0:1, n,replace=TRUE)
 
  # true event time
  T = rweibull(n, shape=1, scale=lambdaT*exp(-beta1*x1)) 
  C = rweibull(n, shape=1, scale=lambdaC)                # censoring time
  time = pmin(T,C)                                       # observed time is min of censored and true
  event = time==T                                        # set to 1 if event is observed

  coxph(Surv(time, event)~ x1 , method="breslow")
  
  survfit <- survfit(Surv(time,event) ~ x1)
  plot(survfit, ylab="Survival probability", xlab="Time")
  
  summary(survfit)
  print(survfit)
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  rm(list=ls())
  set.seed(876)
  n = 100
  beta1 = 0;      # no trt effect
  lambdaT = 14.4  # baseline hazard
  lambdaC = 17.6  # hazard of censoring
  x1 = sample(0:1, n,replace=TRUE)
  
  # true event time
  T = rweibull(n, shape=1, scale=lambdaT)
  C = rweibull(n, shape=1, scale=lambdaC)   #censoring time
  time = pmin(T,C)  # observed time is min of censored and true
  event = time==T   # set to 1 if event is observed
  
  library(survival)
  f <- coxph(Surv(time, event)~ 1 , method="breslow")
  mfit<-survfit(f)
  mfit
  
  survfit <- survfit(Surv(time,event) ~ 1)
  plot(survfit, ylab="Survival probability", xlab="Time")
  
  # theoretical median, self learning txt p301
  -log(.5)*lambdaT   # qweibull(.5,1,15.6)
  -log(.7)*lambdaT
  
  ## surv rates 5 10 year survival 70% and 50%
  ## detect 20% increase in haz rate, what does this inply for 5 n 10 year surv rates
  ## S1(t) = So(t)expBx =0.7^1.2=65.2 0.5^1.2=34.5%
  
  
  mfit$surv
  survfit <- survfit(Surv(time,event) ~ 1)
  plot(survfit, ylab="Survival probability", xlab="Time")
  
  
  
  (f<- survreg(Surv(time, event) ~ 1 , dist = "w", control = list(maxiter=90) ))
  hr.aft <- exp(-coef(f))^exp(coef(f)["shape"])
  
  summary(survfit)
  quantile(survfit, probs = c(.3,0.5, 0.7), conf.int = FALSE)
  
  
  1 - pweibull(10.8131,  shape = 1, scale = 15.6)
  1 - pweibull(5.564129, shape = 1, scale = 15.6)
  
  qweibull(.5, shape=1,scale=14.4) # what time give 50%
  qweibull(.3, shape=1,scale=14.4) # what time give 30%
  
  # what scale satisfies   
  1 - pweibull(10 ,  shape = 1, scale = 14.4)
  1 - pweibull(5 ,   shape = 1, scale = 14.4)
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  rm(list=ls())
  set.seed(876)
  n = 99 
  beta1 = log(1.2); 
  lambdaT = 14.4  # baseline hazard
  lambdaC = 17.6  # hazard of censoring
  x1 = sample(0:1, n,replace=TRUE)
  
  # true event time
  T = rweibull(n, shape=1, scale=lambdaT)*exp(-beta1*x1)
  C = rweibull(n, shape=1, scale=lambdaC)   #censoring time
  time = pmin(T,C)  # observed time is min of censored and true
  event = time==T   # set to 1 if event is observed
  
  library(survival)
  f <- coxph(Surv(time, event)~ x1 , method="breslow")
  f
  mfit<-survfit(f)
  mfit
  
  survfit <- survfit(Surv(time,event) ~ x1)
  plot(survfit, ylab="Survival probability", xlab="Time")
  
  # theoretical median, self learning txt p301
  -log(.5)*lambdaT   # qweibull(.5,1,15.6)
  -log(.7)*lambdaT
  
  ## surv rates 5 year  survival 70%  
  ## surv rates 10 year survival 50%
  
  ## detect 20% increase in haz rate, what does this inply for 5 n 10 year surv rates
  ## S1(t) = So(t)expBx =0.7^1.2=65.2 
  ###                    0.5^1.2=43.5%
  ## so 5 yr drops to 65.2 and 10 yr drops to 43.5
  
  
  #mfit$surv
  survfit <- survfit(Surv(time,event) ~ x1)
  plot(survfit, ylab="Survival probability", xlab="Time")
  #survfit
  #summary(survfit)
  quantile(survfit, probs = c(.3,.348,0.5, .565,0.7), conf.int = FALSE)
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  summary(survfit, times=5 )
  summary(survfit, times=10)
  
  
  
  
  
    ## plotting harrells Weibull2 output,not right?
  surv=c(0.7, 0.5)
  times= c(11.9, 23.1)
 
  
  z1 <- -logb(surv[1])
  z2 <- -logb(surv[2])
  t1 <- times[1]
  t2 <- times[2]
  gamma <- logb(z2/z1)/logb(t2/t1)
  alpha <- z1/(t1^gamma)
  
  #eB = (mo/m1)^shape therneau page 62
  
  # 
  # x<-curve(dweibull(x, shape=1, scale = lambdaT), from=0, to=40)
  # x$y <- x$y/max(x$y)
  # plot(x, type = "l", lty = 1)
  # 
  # y <- x$y^exp(beta1)
  # lines(y~x$x)
  # 
  
  #10.3.1
 # scale <- sigma <- 1/alpha
  #mu= -log(lambda)
  
# St <- exp( -exp(-mu/scale)  *  time^(1/scale)  )
 
 

  #https://www.itl.nist.gov/div898/handbook/eda/section3/eda3668.htm
  x <- curve(dweibull(x, shape=gamma, scale = 1/alpha), from=0, to=40)  # note this..correct?
  y <- x$y^exp(beta1)
  lines(y~x$x)
  
  
  
 #  y <- x$y
 #  x <- x$x
 #  x <- x[!is.na(x) & !is.infinite(x)]
 #  y <- y[!is.na(y) & !is.infinite(y)]
 #  
 #  
 # # x<-x[-1]
 #  
 #  y <-  (y) / (max(y))
 #  plot(y~x, type = "l", lty = 1)
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
 