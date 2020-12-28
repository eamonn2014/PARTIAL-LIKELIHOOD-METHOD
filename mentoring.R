 
# describe hr
require(Hmisc)
require(rms)

rm(list=ls())
d <- "id event surv dead x
1	1	5	1	1
2	2	8	1	1
3	3	10	1	1
4	4	13	1	1
5	5	18	1	1
6	6	23	1	0
7	7	24	1	1
8	8	25	1	1
9	9	26	1	1
10	10	31	1	1
11	11	35	1	1
12	12	40	1	1
13	13	41	1	1
14	14	47	1	0
15	15	48	1	1
16	16	50	1	1
17	17	59	1	1
18	18	61	1	1
19	19	68	1	1
20	20	69	1	0
21	NA	70	0	0
22	21	71	1	1
23	NA	71	0	0
24	NA	76	0	1
25	NA	100	0	0
26	NA	101	0	0
27	NA	105	0	1
28	NA	107	0	1
29	NA	109	0	1
30	22	113	1	1
31	NA	116	0	1
32	23	118	1	1
33	24	143	1	1
34	25	148	1	0
35	NA	154	0	1
36	NA	162	0	1
37	26	181	1	0
38	NA	188	0	1
39	NA	198	0	0
40	NA	208	0	0
41	NA	212	0	0
42	NA	212	0	1
43	NA	217	0	1
44	NA	224	0	0
45	NA	225	1 1
"

d <- read.table(textConnection(d), header = TRUE)
closeAllConnections()
d

d <- sapply(d, function(x) as.numeric(x))
d <- as.data.frame(d)

require(survival)
f <- coxph(Surv(surv, dead) ~ x, data = d)
summary(f)
f$loglik                                       # null and maximised
 
 
# change the last event no difference
d[nrow(d),"dead"] <-0

f <- coxph(Surv(surv, dead) ~ x, data = d)
summary(f)
f$loglik

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# actual time not important
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
guess=0 #0.9093 
d$expB <- exp(guess*d$x)  # always 1
# maximum likelihood
d$part3 <- log(rev(cumsum(rev(d$expB))))  # log(45) log(44) log(43)....
d$part1 <- d$dead  
d$part2 <- guess*d$x 
d$likelihoodi <- d$part1*(d$part2 - d$part3)
sum(d$likelihoodi)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~plot
f1 <- Surv(d$surv, d$dead)

plot(survfit(f1 ~ d$x), col=c("purple", "red"), fun="cloglog", xlab="Time", ylab="log(-log(Survival)")

hazard.ratio.plot(d$x, f1, e=20, legendloc='ll')
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~log likelihood function
 
# function to LL

loglike <- function(x, dead=d$dead, indep=d$x ) {
  
  expB <- exp(x*indep)
  part1 <- dead  
  part2 <- x*indep 
  part3 <- log(rev(cumsum(rev(expB))))
  likelihoodi <- part1*(part2 - part3)
  L <-sum(likelihoodi)
  return(c(L))
  
}

loglike(0 )  # try a log hr and see if we can minimise
loglike(0.9093  )
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


hr <- optimize(loglike, lower = -2, upper = 2, dead=d$dead, indep=d$x, tol = 0.00000001, maximum=T)
hr$maximum  # here we get the log hr that maximises the data!! 

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~lets simulate data~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# require(rms)
# 
# set.seed(731)
# 
# n <- 5000
# hr <- 2
# 
# trt <- sample(0:1, n,  rep=TRUE, prob=c(.6, .4))
# 
# cens <- 15*runif(n)
# 
# h <- exp(log(hr)*(trt==1))
# 
# dt <- -log(runif(n))/h
# label(dt) <- 'Follow-up Time'
# 
# e <- ifelse(dt <= cens,1,0)
# 
# dt <- pmin(dt, cens)
# 
# units(dt) <- "Year"
# 
# d <- data.frame(cbind(dt,e,trt=trt))
# 
# dd <- datadist(d)
# options(datadist='dd')
# 
# S <- Surv(dt,e)
# (f <- cph(S ~  trt, x=TRUE, y=TRUE,d))
# (f1 <-  summary(f))
# 
# f$loglik # null and maximised
# 
# d <- plyr::arrange(d,dt)  # sort by time
# 
# 
# survplot(f,logt=TRUE, loglog=TRUE, col=c("purple", "orange"))   #Check for Weibull-ness (linearity)
# 
# plot(survfit(S~ d$trt), col=c("purple", "red"), fun="cloglog", xlab="Time", ylab="log(-log(Survival)")
# 
# hazard.ratio.plot(d$trt, S, e=20, legendloc='ll', xlab='Time')
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# reference rms help, function......................
# expect no ties with continuous time

coxdata <- function(n, allocation, hr, baseline) { 
  
  #n=1000; allocation =.5; hr=2; baseline=.4
  
  trt <- sample(0:1, n,  rep=TRUE, prob=c(1-allocation, allocation))
  
  cens <- 15*runif(n)
  
  h <- baseline*exp(log(hr)*(trt==1))
  
  dt <- -log(runif(n))/h
  label(dt) <- 'Follow-up Time'
  
  e <- ifelse(dt <= cens,1,0)
  
  dt <- pmin(dt, cens)
  
  units(dt) <- "Year"
  
  d <<- data.frame(cbind(dt,e,trt=trt))  ##why the << required to circumvent error?
  
  dd <<- datadist(d)
  options(datadist='dd')
  
  foo <-d
  # S <- Surv(dt,e)
  f <- cph(Surv(dt,e) ~  trt, x=TRUE, y=TRUE, data=d )
  f0 <- f$coefficients[[1]] #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  LL1 <- f$loglik[2]
  LL0 <- f$loglik[1]
  
  sf <- summary(f)
  
  f1 <- survfit(Surv(dt,e) ~ trt, data = d)
  
  np <- npsurv(Surv(dt,e) ~ trt,d)
  
  S <- Surv(d$dt, d$e)
  
  #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  
  d <- plyr::arrange(d,dt)
  d$dt <- dd<- NULL
  d$dt <- sort(2*rexp(nrow(d)))#
  
  dx <<- datadist(d)
  options(datadist='dx')
  f0a <- cph(Surv(dt,e) ~  trt, x=TRUE, y=TRUE, data=d )
  f0a <- f0a$coefficients[[1]]
  f2 <- survfit(Surv(dt ,e)  ~ trt, data = d)
  
  #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  
  return(list(f=f, d=d, f1=f1, sf=sf, np=np, LL1=LL1, LL0=LL0, S=S,                  
              
              f0=f0,  f2=f2, f0a=f0a, foo=foo))
  
}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

res <- coxdata(n=100, allocation =.5, hr=2, baseline=.4)

# library(survival)
# library(survminer)
# d <- res$foo
# 
# KM_fit <- survfit(Surv(dt, e) ~ trt ,data = d)
# NA_fit <- survfit(Surv(dt, e) ~ trt, data = d,   type = "fleming-harrington")
# plot(KM_fit, conf.int = FALSE, xlab = "Time", ylab = "Survival Probability")
# lines(NA_fit, conf.int = FALSE, col = 'red')
# 
# 
# with(KM_fit, cbind(time, n.risk, n.event, n, surv,std.err,cumhaz))
# print(KM_fit, print.rmean = TRUE)
# dat_km <- fortify(KM_fit)
# head(dat_km)
# 
# fortify(KM_fit, fun = "cumhaz")
# 
# 
# print(NA_fit, print.rmean = TRUE)
# dat_na <- fortify(NA_fit)
# head(dat_na)
# 
# 
# 











# 
# #  d$id <- 1
# #   d$dt<-d$dt*100
# 
# su_obj <- Surv(d$dt, d$e)
# str(su_obj)
# 
# 
# 
# fit_km <- survfit(su_obj ~ 1, data = d)
# print(fit_km, print.rmean = TRUE)
# 
# dat_km <- fortify(fit_km)
# head(dat_km)
# 
# 
# 
# 
# 
# fit_fh <- survfit(su_obj ~ 1, data = d, type = "fleming-harrington", conf.type = "log-log")
# dat_fh <- fortify(fit_fh)
# ## for the Nelson-Aalen estimator of the cumulative hazard
# #dat_fh <- fortify(fit_fh, fun = "cumhaz")
# head(dat_fh)
# head(dat_km)
# 
# 
# ggplotly(
#   ggplot() +
#     geom_step(data = dat_km, aes(x = time, y = surv, colour = "K-M")) +
#     geom_step(data = dat_fh, aes(x = time, y = surv, colour = "N-A")) +
#     #  geom_step(data = dat_lt, aes(x = cuts[-length(cuts)], y = surv, colour = "LT")) +
#     labs(x = "Time", y = "Survival", colour = "Estimator") +
#     theme_classic()
# )
# 
# 
# 
# 
# 
# 
# 

















#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#















# 
# 
# 
# 
# (f <- res$f)
# summary(f)
# d <- res$d
# head(d)
# 
# S<- Surv(d$dt, d$e)
# 
# survplot(f,logt=TRUE, loglog=TRUE, col=c("orange", "purple"))
# 
# hazard.ratio.plot(d$trt, S, e=20, legendloc='ll')
# 
# f$loglik
# #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 
# 
# # app function
# coxdata <- function(n, allocation, hr, baseline) { 
#   
#   #n=1000; allocation =.5; hr=2; baseline=.4
#   
#   trt <- sample(0:1, n,  rep=TRUE, prob=c(1-allocation, allocation))
#   
#   cens <- 15*runif(n)
#   
#   h <- baseline*exp(log(hr)*(trt==1))
#   
#   dt <- -log(runif(n))/h
#   label(dt) <- 'Follow-up Time'
#   
#   e <- ifelse(dt <= cens,1,0)
#   
#   dt <- pmin(dt, cens)
#   
#   units(dt) <- "Year"
#   
#   d <<- data.frame(cbind(dt,e,trt=trt))  ##why the << required to circumvent error?
#   
#   dd <<- datadist(d)
#   options(datadist='dd')
#   
#   # S <- Surv(dt,e)
#   f <- cph(Surv(dt,e) ~  trt, x=TRUE, y=TRUE, data=d )
#   
#   LL1 <- f$loglik[2]
#   LL0 <- f$loglik[1]
#   
#   sf <- summary(f)
#   
#   f1 <- survfit(Surv(dt,e) ~ trt, data = d)
#   
#   np <- npsurv(Surv(dt,e) ~ trt,d)
#   
#   S <- Surv(d$dt, d$e)
#   
#   return(list(f=f, d=d, f1=f1, sf=sf, np=np, LL1=LL1, LL0=LL0, S=S))
#   
# }
# # function that we will use to maximum the partial log likelihood
# 
# ## guess, data, dead var, trt var.
# loglike2 <- function(x, dat, dead, indep , time) {
#   
#   dd <- dat          #make another data object
#   dd$dead <- dead    # take the key variables to run Cox PH
#   dd$indep <- indep
#   dd$time <- time
#   
#   ## run the analysis to get hr and loglike
#   ddd <<- datadist(dd)
#   options(datadist='ddd')
#   
#   S <- Surv(time,dead)  # run Cox PH
#   f <- cph(S ~  indep, x=TRUE, y=TRUE,dd)
#   
#   #~~~~~~~~~~extract hr and loglikelihood at null and maximised log likelihood
#   dd$hr <- exp(f$coefficients)
#   dd$lognull <- f$loglik[[1]]
#   dd$lognmax <- f$loglik[[2]]
#   
#   #~~~~~~~~~using our guess x calculate log likelihood by jand
#   dd$expB <- exp(x*dd$indep)
#   dd$part1 <- dd$dead  
#   dd$part2 <- x*dd$indep 
#   
#   dd <- plyr::arrange(dd,time)
#   
#   dd$part3 <- log(rev(cumsum(rev(dd$expB))))
#   dd$likelihoodi <- dd$part1*(dd$part2 - dd$part3)
#   dd$guess <- exp(x)
#   dd$L <- sum(dd$likelihoodi)
#   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#   dd <- as.data.frame(dd)
#   dd$dead <- dd$indep <- dd$part1 <- dd$part2 <- dd$part3 <- NULL
#   # return(head(dd))
#   return(dd)
#   
# }
# 
# # this will work in app
# dummy <- coxdata(n=100, allocation =.5, hr=2, baseline=.4)  # app function
# 
# loglike2(0, dat=dummy$d, dead=dummy$d$e, indep=dummy$d$trt, time=dummy$d$dt)  # try
# 
# loglike2(dummy$f$coefficients[[1]], dat=dummy$d, dead=dummy$d$e, indep=dummy$d$trt, time=dummy$d$dt)  # try
# 
# 
# 
# 
#  







