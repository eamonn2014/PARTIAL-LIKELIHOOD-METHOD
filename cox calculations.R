 

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
45	NA	225	0	1
"

d <- read.table(textConnection(d), header = TRUE)
closeAllConnections()
d

d <- sapply(d, function(x) as.numeric(x))
d <- as.data.frame(d)

require(survival)
f <- coxph(Surv(surv, dead) ~ x, data = d)
summary(f)
f$loglik
 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
guess=0
d$expB <- exp(guess*d$x)
# maximum likelihood
d$part3 <- log(rev(cumsum(rev(d$expB))))
d$part1 <- d$dead  
d$part2 <- guess*d$x 
d$likelihoodi <- d$part1*(d$part2 - d$part3)
sum(d$likelihoodi)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~plot
f1 <- Surv(d$surv, d$dead)

plot(survfit(f1 ~ d$x), col=c("purple", "red"), fun="cloglog", xlab="Time", ylab="log(-log(Survival)")

hazard.ratio.plot(d$x, f1, e=20, legendloc='ll')
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~log likelihood function



# res <- coxdata(n=100, allocation=.5, hr=2)
# f <- res$f
# f$loglik
# 
# d <- res$d
# loglike(0 , dead=d$e, indep=d$trt )  
# loglike(-0.025 , dead=d$e, indep=d$trt )  




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

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

 
 hr <- optimize(loglike, lower = -2, upper = 2, dead=d$dead, indep=d$x, tol = 0.00000001, maximum=T)
 hr$maximum  # here we get the log hr that maximises the data!! 
 
 #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~lets simulate data~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 
    require(rms)
 
    set.seed(731)
    
    n <- 5000
    hr <- 2

    trt <- sample(0:1, n,  rep=TRUE, prob=c(.6, .4))
    
    cens <- 15*runif(n)
    
    h <- exp(log(hr)*(trt==1))
    
    dt <- -log(runif(n))/h
    label(dt) <- 'Follow-up Time'
    
    e <- ifelse(dt <= cens,1,0)
    
    dt <- pmin(dt, cens)
    
    units(dt) <- "Year"
    
    d <- data.frame(cbind(dt,e,trt=trt))
 
    dd <- datadist(d)
    options(datadist='dd')
    
    S <- Surv(dt,e)
    (f <- cph(S ~  trt, x=TRUE, y=TRUE,d))
    (f1 <-  summary(f))
    
    f$loglik # null and maximised
    
    d <- plyr::arrange(d,dt)  # sort by time
    
    
    survplot(f,logt=TRUE, loglog=TRUE, col=c("purple", "orange"))   #Check for Weibull-ness (linearity)
    
    plot(survfit(S~ d$trt), col=c("purple", "red"), fun="cloglog", xlab="Time", ylab="log(-log(Survival)")
     
    hazard.ratio.plot(d$trt, S, e=20, legendloc='ll', xlab='Time')
   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    # reference rms help, function......................
    # expect no ties with continuous time
    
    coxdata <- function (n, allocation, hr) { 
      
      trt <- sample(0:1, n,  rep=TRUE, prob=c(1-allocation, allocation))
      
      cens <- 15*runif(n)
      
      h <- exp(log(hr))*(trt==1)
      
      dt <- -log(runif(n))/h
      
      label(dt) <- 'Follow-up Time'
      
      e <- ifelse(dt <= cens,1,0)
      
      dt <- pmin(dt, cens)
      
      units(dt) <- "Year"
      
      d <- data.frame(cbind(dt, e, trt=trt))
      
      d <- plyr::arrange(d, dt)
      
      dd <<- datadist(d)
      options(datadist='dd')
      
      S <- Surv(dt,e)
      
      f <- cph(S ~  trt, x=TRUE, y=TRUE,d)
      
  
      return(list(f=f,d=d))
      
 }
    
    res <- coxdata(n=100, allocation=.5, hr=2)
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
  
    (f <- res$f)
    summary(f)
     d <- res$d
    head(d)
    
    S<- Surv(d$dt, d$e)
    
    survplot(f,logt=TRUE, loglog=TRUE, col=c("orange", "purple"))
    
    hazard.ratio.plot(d$trt, S, e=20, legendloc='ll')
     
    f$loglik
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     
    
    # app function
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
      
      # S <- Surv(dt,e)
      f <- cph(Surv(dt,e) ~  trt, x=TRUE, y=TRUE, data=d )
      
      LL1 <- f$loglik[2]
      LL0 <- f$loglik[1]
      
      sf <- summary(f)
      
      f1 <- survfit(Surv(dt,e) ~ trt, data = d)
      
      np <- npsurv(Surv(dt,e) ~ trt,d)
      
      S <- Surv(d$dt, d$e)
      
      return(list(f=f, d=d, f1=f1, sf=sf, np=np, LL1=LL1, LL0=LL0, S=S))
      
    }
  # function that we will use to maximum the partial log likelihood
   
  ## guess, data, dead var, trt var.
  loglike2 <- function(x, dat, dead, indep , time) {
    
    dd <- dat          #make another data object
    dd$dead <- dead    # take the key variables to run Cox PH
    dd$indep <- indep
    dd$time <- time
    
    ## run the analysis to get hr and loglike
    ddd <<- datadist(dd)
    options(datadist='ddd')
    
    S <- Surv(time,dead)  # run Cox PH
    f <- cph(S ~  indep, x=TRUE, y=TRUE,dd)
    
    #~~~~~~~~~~extract hr and loglikelihood at null and maximised log likelihood
    dd$hr <- exp(f$coefficients)
    dd$lognull <- f$loglik[[1]]
    dd$lognmax <- f$loglik[[2]]
    
    #~~~~~~~~~using our guess x calculate log likelihood by jand
    dd$expB <- exp(x*dd$indep)
    dd$part1 <- dd$dead  
    dd$part2 <- x*dd$indep 
    
    dd <- plyr::arrange(dd,time)
    
    dd$part3 <- log(rev(cumsum(rev(dd$expB))))
    dd$likelihoodi <- dd$part1*(dd$part2 - dd$part3)
    dd$guess <- exp(x)
    dd$L <- sum(dd$likelihoodi)
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    dd <- as.data.frame(dd)
    dd$dead <- dd$indep <- dd$part1 <- dd$part2 <- dd$part3 <- NULL
   # return(head(dd))
    return(dd)
    
  }
  
  # this will work in app
  dummy <- coxdata(n=100, allocation =.5, hr=2, baseline=.4)  # app function
  
  loglike2(0, dat=dummy$d, dead=dummy$d$e, indep=dummy$d$trt, time=dummy$d$dt)  # try
  
  loglike2(dummy$f$coefficients[[1]], dat=dummy$d, dead=dummy$d$e, indep=dummy$d$trt, time=dummy$d$dt)  # try
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # take the data drom coxdata

  loglike2(0,      dat=d, dead=d$e, indep=d$trt, time=d$dt)  # try a log hr and see if we can minimise
  
  
  loglike2(f$coefficients[[1]],  dat=d, dead=d$e, indep=d$trt, time=d$dt)  # try a log hr and see if we can minimise
  
  
  # this will work in app
  dummy <- coxdata(n=100, allocation =.5, hr=2, baseline=.4)  # app function
  
  loglike2(0,      dat=dummy$d, dead=dummy$d$e, indep=dummy$d$trt, time=dummy$d$dt)  # try
  
  
  
  
  
  
  hr <- optimize(loglike2, lower = -2, upper = 1, dead=d$e, dat=d, 
                 indep=d$trt, tol = 0.0000000001, maximum=T)
  hr$maximum
  
  
  
  
  
  
  
 
  ## add cols to dataset
   
  d$beta <- exp(hr$maximum*d$trt)
 # d$likelihood <- 
  
  part1 <- d$e  
  part2 <- hr$maximum*d$trt 
  part3 <- log(rev(cumsum(rev(d$beta))))
  d$likelihoodi <- part1*(part2 - part3)
  d$Partial.LHood <- sum(d$likelihoodi)
  d$hr <- exp(hr$maximum)
  d
  exp(f$coefficients)

  require('survival')
  require('GGally')
  f1 <- survfit(Surv(dt,e) ~ trt, data = d)
  ggsurv(f1, CI=TRUE)

  
  library(survminer)
  library(plotly)
 
  p1 <- ggsurvplot(f1, main = "Kaplan-Meier Curve for the NCCTG Lung Cancer Data")
  plotly::ggplotly(p1[[1]])

  
  S <- Surv(dt,e)
  f <- coxph(S ~  trt, d)
 








