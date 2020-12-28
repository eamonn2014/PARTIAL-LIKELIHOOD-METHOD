

  n <- c(5:40)  # length of DNA
  k <- 4  # bases
  
  factorial(n)/ factorial(n-k)
  
  
  
  # -log(runif(n))/h
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  # hazard will scale it
  
  n <- 100
  h <- 0.5
  i <- runif(n)
  dt <- -log(i)/h
  
  
  
 
  ta <- log(i)  #-infinty to 0
  summary(ta)
  # negate to give 0 to infinty
  (t <- -ta)
  # divide by h, to scale this 
  dt1 <- t/h
  
  summary(t)
  summary(dt1)
  
    plot(t ,i)
    
    plot(dt1 ,i)



#### just understanding
require(survival)

  n <- 10
  h <- 2               # divisor
  i <- runif(n)        # 0-1
  t <- -log(i)/h       # log to give -inf to 0 and negtae , to give 0 to infinity
 
par(mfrow = c(1,2))  
  plot(t ,i)           # basic plots
  event <- rep(1,n) 
  survfit <- survfit(Surv(t,event) ~ 1)
  plot(survfit, ylab="Survival probability", xlab="Time" , conf.int = FALSE)
par(mfrow = c(1,1))


 


#######################bartlett

n <- 100
h <- 1
t <- -log(runif(n))/h

#event <- 1*(t<5)
obstime <- t
#obstime[obstime>=5] <- 5

event=rep(1, n)

survfit <- survfit(Surv(obstime,event) ~ 1)
plot(survfit, ylab="Survival probability", xlab="Time")


plot(survfit, fun="cumhaz", ylab="Cumulative hazard", xlab="Time")
abline(0,h, col='red')

##############################   # sim trt groups

rm(list=ls())
set.seed(2234)

n=100; allocation =.5; hr=2; baseline=.4

trt <- sample(0:1, n,  rep=TRUE, prob=c(1-allocation, allocation))

cens <- 15*runif(n)

h <- baseline*exp(log(hr)*(trt==1))  # hazard function h(t)

dt <- -log(runif(n))/h

label(dt) <- 'Follow-up Time'

e <- ifelse(dt <= cens,1,0)

dt <- pmin(dt, cens)

units(dt) <- "Year"

d <<- data.frame(cbind(dt,e,trt=trt))  ##why the << required to circumvent error?

dd <<- datadist(d)
options(datadist='dd')


# Harrell rms
foo <-d
# S <- Surv(dt,e)
f <- cph(Surv(dt,e) ~  trt, x=TRUE, y=TRUE, data=d )
summary(f)







# cumulative hazards

survfit <- survfit(Surv(dt,e) ~ trt, d)
plot(survfit, fun="cumhaz", ylab="Cumulative hazard", xlab="Time")
# abline(0,h, col='red')
# abline(0,h*hr, col='blue')

abline(0,baseline, col='red')
abline(0,baseline*hr, col='blue')

#############################
 ###why are these so different?


survplot(f,logt=TRUE, loglog=TRUE, 
         col=c("red", "lightblue"))

f1 <- survfit(Surv(dt,e) ~ trt, data = d)  # Harrell
 

 ggsurvplot(f1, fun = "cloglog",  
                 main = "Complementary log−log" ,
                 legend.title = "Trt.")
 

 # hmm?
 b <- basehaz(coxph(Surv(dt,e)~trt),centered=FALSE)
 
 b$trt <- d$trt
 
 plot(b[,1],b[,2])
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 
 # cfit <- cph(Surv(dt, e) ~ trt, data = d, surv=TRUE, x=TRUE, y=TRUE)
 # survest(cfit,newdata=expand.grid(trt=d$trt) , 
 #         times=seq(0,1000,by=1)
 # )$surv
 # 
 # 
 # km <- survfit(Surv(dt, e)~trt, data=d)
 # #urvest <- stepfun(km$time, c(1, km$surv))
 #  trt <- c(  rep(0, km$strata[1])   ,  rep(1, km$strata[2])   )
 # ti <- cbind(km$n.event, km$time, km$surv, trt)
 # names(ti) <- c("e","dt","s","trt")
 # ti <- as.data.frame(ti)
 # 
 #  require(pec)
 # f <- cph(Surv(dt,e) ~  trt, x=TRUE, y=TRUE, data=d , surv=TRUE)
 # predictSurvProb(f, times=dt)
 # 
 # # as follows
 # x$trt=c(0,1)
 # SS <- predictSurvProb(f,newdata=x$trt,times=dt)
 # 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 # getting probability thereneus
 
 fit <- coxph(Surv(dt, e) ~ trt, data = d) 
 temp <- data.frame(x=levels(d$trt))
 expfit <- survfit(fit, temp)
 plot(expfit, xscale=365.25, ylab="Expected")
 plot(expfit )
 
 time <-    expfit$time
 censor <-  expfit$n.censor
 surv <-    unique( expfit$surv,    MARGIN=2 )
 cumhaz <-  unique( expfit$cumhaz,  MARGIN=2 )
 std.err <- unique( expfit$std.err, MARGIN=2 )
 lower <-   unique( expfit$lower,   MARGIN=2 )
 upper <-   unique( expfit$upper,   MARGIN=2 )
 
 dat <- cbind(time, censor, surv1=surv[,1] , surv2=surv[,2], cumhaz1=cumhaz[,1], cumhaz2=cumhaz[,2])
 
 dat <- as.data.frame(dat)
 
 # lets calc log(-log(S))
 
 dat$surv1_cll <- log(-log(dat["surv1"]))
 dat$surv2_cll <- log(-log(dat["surv2"]))
 
 
 # here we manage to duplicate frank harrells plot!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 par(mfrow=c(1,2))
 y <- as.numeric(unlist(dat[,7]))
 x <- log(as.numeric(dat$time ))
 y2 <- as.numeric(unlist(dat[,8]))
 
 
  plot(y~x , type="l",col="red" , xlim=c(-2,3))
  lines(x,y2,col="green")
 
  fit <- cph(Surv(dt, e) ~ trt, data = d, x=TRUE, y=TRUE) 
  survplot(fit,logt=TRUE, loglog=TRUE, 
           col=c("red", "lightblue"))
  par(mfrow=c(1,1))
  ###########################################################
  
  # how is this created?
  
  require("survival")
  require("survminer")
  fit<- survfit(Surv(dt, e) ~ trt, data = d)
  ggsurvplot(fit, data = d, fun = "cloglog")
  
  
  
  library(survival)
  library(MASS)
  plot(survfit(Surv(dt, e) ~ trt),   fun="cloglog")
  
  
  myfun=function(p){return(log(-log(p)))}
  plot(survfit(Surv(dt, e) ~ trt),   fun=myfun, log="x")
  
  plot(survfit(Surv(dt, e) ~ trt)  )
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
 
 fit <- cph( Surv(dt,e)~trt, x=TRUE, y=TRUE, surv=TRUE)
 ss <- survfit(fit)
 
 par(mfrow=c(2,2))
 plot(ss, fun = "cumhaz")
 plot(dat[,6]~ dat[,1], ylim=c(0,12))

 #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 
 survplot(fit,logt=TRUE, loglog=TRUE, 
          col=c("red", "lightblue"))
 
 f1 <- survfit(Surv(dt,e) ~ trt, data = d)  # Harrell
 
 
 ggsurvplot(f1, fun = "cloglog",  
            main = "Complementary log−log" ,
            legend.title = "Trt.")
 
 par(mfrow=c(1,1))
 
 
 
 ##############################################################################################################
 
 #https://stats.stackexchange.com/questions/105881/how-to-simulate-survival-times-using-true-base-line-hazard-function
 
 tdom <- seq(0, 5, by=0.01)
 haz <- rep(0, length(tdom))
 haz[tdom <= 1] <- exp(-0.3*tdom[tdom <= 1])
 haz[tdom > 1 & tdom <= 2.5] <- exp(-0.3)
 haz[tdom > 2.5] <- exp(0.3*(tdom[tdom > 2.5] - 3.5))
 
 
 #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 
 haz <- exp(-0.3*tdom )
 haz <- rep(1,length(tdom))
 haz <- exp(0.3*tdom )
 
 #az <- rweibull(length(tdom), shape=.1, scale=5)

 cumhaz <- cumsum(haz*0.01)
 Surv <- exp(-cumhaz)
 par(mfrow=c(3,1))
 plot(tdom, haz,    type='l', xlab='Time domain', ylab='Hazard')
 plot(tdom, cumhaz, type='l', xlab='Time domain', ylab='Cumulative hazard')
 plot(tdom, Surv,   type='l', xlab='Time domain', ylab='Survival')

 #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 
 
  # generate 100 random samples:
  u <- runif(100)
  failtimes <- tdom[colSums(outer(Surv, u, `>`))]
  
  dev.off()
  library(survival)
  plot(survfit(Surv(failtimes)~1))
 
 
 
 
 
  # Weibull parametrisation
  y<-rweibull(1000, shape=2, scale=5)
  
  survreg(Surv(y)~1, dist="weibull")
  # survreg scale parameter maps to 1/shape, linear predictor to log(scale)
  
 #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  
  require(simsurv)
  N <- 50 # total number of patients
  # Define covariate data
  covs <- data.frame(id = 1:N,
                     trt = rbinom(N, 1, 0.5))
  # Define true coefficient (log hazard ratio)
  pars <- c(trt = -0.5)
  # Simulate the event times
  times <- simsurv(dist = 'weibull',
                   lambdas = 0.1,
                   gammas = 1.5,
                   x = covs,
                   betas = pars) 
 
 
  h <- -log(times$eventtime)
 
 