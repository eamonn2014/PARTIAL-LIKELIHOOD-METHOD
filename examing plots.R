


##############################   # sim trt groups

rm(list=ls())
set.seed(2234)
require(rms)
require(survminer)
require(survival)

n=800; allocation =.5; hr=2; baseline=.4

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
# HARRELS COX
f <- cph(Surv(dt,e) ~  trt, x=TRUE, y=TRUE, data=d )
summary(f)


###########KM#######################################################
survfit <- f1 <- survfit(Surv(dt,e) ~ trt, data = d)
A <- ggsurvplot(f1, main = "Kaplan-Meier Curve" )


####################################################################
# cumulative hazards. no issue here

plot(survfit, fun="cumhaz", ylab="Cumulative hazard", xlab="Time", log=TRUE)
 


plot(survfit, fun="cumhaz", ylab="Cumulative hazard", xlab="Time", log=FALSE)
abline(0,baseline, col='red')
abline(0,baseline*hr, col='blue')














#############################
###why are these so different?

# 
# survplot(f,logt=TRUE, loglog=TRUE, 
#          col=c("red", "lightblue"))
# 
# f1 <- survfit(Surv(dt,e) ~ trt, data = d)  # Harrell
# 
# 
# ggsurvplot(f1, fun = "cloglog",  
#            main = "Complementary log−log" ,
#            legend.title = "Trt.")


# hmm?
# b <- basehaz(coxph(Surv(dt,e)~trt),centered=FALSE)
# 
# b$trt <- d$trt
# 
# plot(b[,1],b[,2])
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
# getting probability Terry Thereneau

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
summary(dat)

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

 
library(MASS)
plot(survfit(Surv(dt, e) ~ trt),   fun="cloglog")


myfun=function(p){return(log(-log(p)))}
plot(survfit(Surv(dt, e) ~ trt),   fun=myfun, log="x")

plot(survfit(Surv(dt, e) ~ trt)  )


###########################################################p95 moore book manually create cloglog
# I can create the clog log 2nd version!!!!
fit<- survfit(Surv(dt, e) ~ trt, data = d)

  time.LA <- fit$time
  surv.LA <- fit$surv
  cloglog.LA <- log(-log(surv.LA))
  logtime.LA <- log(time.LA)
  plot(cloglog.LA ~ logtime.LA, type="s", col="blue", lwd=2, xlim=c(-5, 2))
 
  
  dat2 <- data.frame(cbind(time.LA, surv.LA , cloglog.LA, logtime.LA )) 
  
  summary(dat2)
  summary(dat)
  
  # result.surv.LA <- survfit(Surv(dt) ~ e, subset=(trt == 0))
  #  time.LA <- result.surv.LA$time
  #  surv.LA <- result.surv.LA$surv
  #  cloglog.LA <- log(-log(surv.LA))
  #  logtime.LA <- log(time.LA)
  #  result.surv.M <- survfit(Surv(dt) ~ e, subset=(trt == 1))
  #  time.M <- result.surv.M$time
  #  surv.M <- result.surv.M$surv
  #  cloglog.M <- log(-log(surv.M))
  #  logtime.M <- log(time.M)
  #  plot(cloglog.LA ~ logtime.LA, type="s", col="blue", lwd=2,  xlim=c(-5, 2))
  #  lines(cloglog.M ~ logtime.M, col="red", lwd=2, type="s")
  #  legend("bottomright", legend=c("0",  "1"), col=c("blue","red"), lwd=2)
  # 
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
###########################################################
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

  # time the same
  summary(x)
  summary(logtime.LA)

  








##############################################################




















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