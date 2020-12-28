## breaking down the clog log plots
 

rm(list=ls())
set.seed(2234)
require(rms)
require(survminer)
require(survival)

n=20; allocation =.5; hr=2; baseline=.4

trt <- sample(0:1, n,  rep=TRUE, prob=c(1-allocation, allocation))
cens <- 15*runif(n)
h <- baseline*exp(log(hr)*(trt==1))  # hazard function h(t)
dt <- -log(runif(n))/h
label(dt) <- 'Follow-up Time'
e <- ifelse(dt <= cens,1,0)
dt <- pmin(dt, cens)
units(dt) <- "Year"

d <- NULL
d <- data.frame(cbind(dt,e,trt=trt))  ##why the << required to circumvent error?
dd <- datadist(d)
options(datadist='dd')
 
# HARRELS COX
f <- cph(Surv(dt,e) ~  trt, x=TRUE, y=TRUE, data=d , surv=TRUE)
summary(f)

fit = survfit(Surv(dt,e) ~  trt, data=d) 
summary(fit)

fit = survfit(coxph(Surv(dt,e)~1, data=d[d$trt==1,]), type="aalen")
summary(fit)
-log(fit$surv) 
fit = survfit(coxph(Surv(dt,e)~1, data=d[d$trt!=1,]), type="aalen")
summary(fit)
-log(fit$surv) 
 

####################
#https://sas-and-r.blogspot.com/2010/05/example-739-nelson-aalen-estimate-of.html
####################

library(survival)
time =  c(0.5, 1,1,2,2,3,4,5,6,7,8,9,10,12,14,14,17,20, 21)
event = c(FALSE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, 
          TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, FALSE, TRUE, FALSE, 
          TRUE, FALSE)
ds = data.frame(time, event)
fit = survfit(Surv(time, event) ~ 1, data=ds)


summary(fit)

calcna = function(time, event) {
  na.fit = survfit(coxph(Surv(time,event)~1), type="aalen")
  jumps = c(0, na.fit$time, max(time))
  # need to be careful at the beginning and end
  surv = c(1, na.fit$surv, na.fit$surv[length(na.fit$surv)])
  
  # apply appropriate transformation
  neglogsurv = -log(surv)   
  
  # create placeholder of correct length
  naest = numeric(length(time))  
  for (i in 2:length(jumps)) {
    naest[which(time>=jumps[i-1] & time<=jumps[i])] = 
      neglogsurv[i-1]   # snag the appropriate value
  }
  return(naest)
}


newna = calcna(time, event)
cbind(time, newna)



# s1 <- survfit( coxph( Surv( time,event) ~ 1), type = "aalen")
# 
#  s1$surv <- (-log(s1$surv))
# s1$upper<- (-log(s1$upper))
# s1$lower <- (-log(s1$lower))
# #s1$std.err<- (-log(s1$std.err))
# 
# plot(s1,first=0)
















###########KM#######################################################
f1 <- survfit(Surv(dt,e) ~ trt, data = d)
A <- ggsurvplot(f1, main = "Kaplan-Meier Curve" , fun='pct')
A

####################################################################################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# getting probability Terry Thereneau
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####################################################################################


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

fit <- cph(Surv(dt, e) ~ trt, data = d, x=TRUE, y=TRUE, surv=TRUE)
  
survplot(fit, type="kaplan-meier", conf.int = TRUE, 
         col.fill=c("orange1","steelblue1"), grid=TRUE, what='survival')


plot(cox.zph(fit, transform="identity" ))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~




















#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
S <- Surv(dt,e)
f <- npsurv(S ~ trt)
survplot(f)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

allsurv <- function(fit){
  ggsurvplot(
    fit,
    pval = TRUE,
    pval.coord = c(200, 0.10), 
    conf.int = TRUE,
    xlab = "Days",
    ggtheme = theme_light(), 
    surv.median.line = "hv", 
    legend.labs = c("0","1"),
    legend.title = "",
    palette = c("#8C3F4D","#3E606F"))
}

fit_all <- function(x){
  surv_fit(Surv(dt, e) ~ trt, data = x)
}

allsurv(fit_all(d))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


fit <- coxph(Surv(dt, e) ~ trt, data = d) 
#temp <- data.frame(x=levels(d$trt))
temp <- data.frame((d$trt))

expfit <- survfit(fit, temp)
#plot(expfit, xscale=365.25, ylab="Expected")
plot(expfit )      # this will be based on cox 

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

########################################################################################################
#######################################################################################################

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



as.numeric(unlist(dat$surv1_cll))
dat2$cloglog.LA


 










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


# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# ###########################################################
# par(mfrow=c(1,2))
# 
# y <- as.numeric(unlist(dat[,7]))
# x <- log(as.numeric(dat$time ))
# y2 <- as.numeric(unlist(dat[,8]))
# 
# 
# plot(y~x , type="l",col="red" , xlim=c(-2,3))
# lines(x,y2,col="green")
# 
# fit <- cph(Surv(dt, e) ~ trt, data = d, x=TRUE, y=TRUE) 
# survplot(fit,logt=TRUE, loglog=TRUE, 
#          col=c("red", "lightblue"))
# par(mfrow=c(1,1))
# ###########################################################
# 
# # time the same
# summary(x)
# summary(logtime.LA)
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# ##############################################################
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# fit <- cph( Surv(dt,e)~trt, x=TRUE, y=TRUE, surv=TRUE)
# ss <- survfit(fit)
# 
# par(mfrow=c(2,2))
# plot(ss, fun = "cumhaz")
# plot(dat[,6]~ dat[,1], ylim=c(0,12))
# 
# #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 
# survplot(fit,logt=TRUE, loglog=TRUE, 
#          col=c("red", "lightblue"))
# 
# f1 <- survfit(Surv(dt,e) ~ trt, data = d)  # Harrell
# 
# 
# ggsurvplot(f1, fun = "cloglog",  
#            main = "Complementary log−log" ,
#            legend.title = "Trt.")
# 
# par(mfrow=c(1,1))