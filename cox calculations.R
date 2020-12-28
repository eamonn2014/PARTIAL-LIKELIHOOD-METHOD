
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# https://stats.stackexchange.com/questions/60238/intuition-for-cumulative-hazard-function-survival-analysis
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# deaths are ages
dx <-  c(3184L, 268L, 145L, 81L, 64L, 81L, 101L, 50L, 72L, 76L, 50L, 
         62L, 65L, 95L, 86L, 120L, 86L, 110L, 144L, 147L, 206L, 244L, 
         175L, 227L, 182L, 227L, 205L, 196L, 202L, 154L, 218L, 279L, 193L, 
         223L, 227L, 300L, 226L, 256L, 259L, 282L, 303L, 373L, 412L, 297L, 
         436L, 402L, 356L, 485L, 495L, 597L, 645L, 535L, 646L, 851L, 689L, 
         823L, 927L, 878L, 1036L, 1070L, 971L, 1225L, 1298L, 1539L, 1544L, 
         1673L, 1700L, 1909L, 2253L, 2388L, 2578L, 2353L, 2824L, 2909L, 
         2994L, 2970L, 2929L, 3401L, 3267L, 3411L, 3532L, 3090L, 3163L, 
         3060L, 2870L, 2650L, 2405L, 2143L, 1872L, 1601L, 1340L, 1095L, 
         872L, 677L, 512L, 376L, 268L, 186L, 125L, 81L, 51L, 31L, 18L, 
         11L, 6L, 3L, 2L)

x <- 0:(length(dx)-1) # age vector

plot((dx/sum(dx))/(1-cumsum(dx/sum(dx))), t="l", xlab="age", ylab="h(t)", 
     main="h(t)", log="y")


plot(cumsum((dx/sum(dx))/(1-cumsum(dx/sum(dx)))), t="l", xlab="age", ylab="H(t)", 
     main="H(t)")


#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

# understanding simple example

require(survival)

  n <- 100
  h <- .5               # divisor, hazard , think of it as scaling...
  i <- runif(n)        # 0-1
  t <- -log(i)/h       # log to give -inf to 0 and negate , to give 0 to infinity

par(mfrow = c(1,3))   # let's plot

  plot(t ,i)          # basic plot
  
  event <- rep(1,n)   # everybody =1, has an event fo this exercise
  survfit <- survfit(Surv(t,event) ~ 1)
  plot(survfit, ylab="Survival probability", xlab="Time" , conf.int = FALSE)

  # https://thestatsgeek.com/2014/03/28/interpreting-changes-in-hazard-and-hazard-ratios/
  plot(exp(-h*t)~t)  #Survival
  
par(mfrow = c(1,1))

##?
b <- basehaz(coxph(Surv(t,event)~1))
plot(sort(t), b[,1] )
 

fit <- coxph( Surv(t,event)~1)

ss <- survfit(fit)

plot(ss, fun = "cumhaz")

plot(Surv(t,event), fun="cumhaz", ylab="Cumulative hazard", xlab="Time")





#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
# https://stackoverflow.com/questions/32019678/obtain-the-baseline-hazard-function-survival-function-from-an-extended-cox-model
model_1<-coxph(Surv(t,event) ~ 1)
mfit<-survfit(model_1)
mfit$surv

all.equal( -log(mfit$surv), mfit$cumhaz )

plot(sort(t),-log(mfit$surv))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
# https://www.itl.nist.gov/div898/handbook/apr/section2/apr222.htm
revrank = c(length(event):1)
haz = event/revrank
plot(sort(t),haz)
cumhaz = cumsum(haz)
plot(sort(t),cumhaz)
 
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
 



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##

## Select failing cases for plotting.
# df = data.frame(time, fail, cumhaz)
# z = subset(df, fail==1)


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~










#Estimate hazard function from right-censored data.
require(muhaz)
haz = muhaz(t,event)
plot(haz)
plot(haz$haz.est)
summary(haz)




 
 
 
 
 muhaz.fit = muhaz(t, event, max.time = max(t))
 # Obtain the estimate of the hazard rate
 muhazrate = muhaz.fit$haz.est
 
 
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
 
 
 
 # Baseline hazard plot (H. Seltman, Feb. 2011)
 # Makes a plot of the baseline hazard based on a coxph model.
 # Works by using basehaz(), R's cumulative hazard function, and then
 # using lowess() smoothing of the simple linear slope estimates.
 #
 # phobj is the result of a coxph() call
 # stratum chooses the stratum for stratified models
 # f is the lowess smoothing parameter
 # add=TRUE is useful for comparing smoothing values (see examples)
 # col sets the curve color 
 # lty sets the curve line type
 #
 # Return value: NULL (invisibly)
 # Side effect: Creates a plot (or adds to an existing plot)
 #
 
 # http://www.stat.cmu.edu/~hseltman/files/baseHazPlot.R
 baseHazPlot = function(phobj, stratum=NULL, f=2/3, 
                        xlab="time", ylab="hazard", add=FALSE, col=1, lty=1, ...) {
   if (!require(survival)) stop("need 'survival' package")
   if (!is(phobj,"coxph")) stop("phobj must be a coxph() result")
   tmp = basehaz(phobj)
   strata = levels(tmp$strata)
   if (is.null(strata)) {
     if (!is.null(stratum)) stop("model was not stratified")
   } else {
     if (is.null(stratum)) stop("specify a stratum number")
     if (!is.numeric(stratum) || stratum!=round(stratum) ||
         stratum<1 || stratum>length(strata))
       stop("bad stratum selection")
     tmp=tmp[tmp$strata==strata[stratum],]
   }
   n = nrow(tmp)
   haztime = apply(cbind(tmp$time[1:(n-1)],tmp$time[2:n]),1,mean)
   hazdiff = apply(cbind(tmp$hazard[1:(n-1)],tmp$hazard[2:n]),1,diff)
   timediff = apply(cbind(tmp$time[1:(n-1)],tmp$time[2:n]),1,diff)
   hazslope =  hazdiff/timediff
   sm = lowess(haztime, hazslope, f=f)
   if (add==FALSE) {
     plot(range(sm$x), range(0,sm$y), type="n", xlab=xlab, ylab=ylab, ...)
   }
   lines(sm$x, sm$y, col=col, lty=lty)
   invisible(NULL)
 }
 
 ## Examples ##
 
   baseHazPlot(coxph(Surv(t,event)~1), main=" hazard")
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 




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
f$loglik

 





require(survival)
f <- coxph(Surv(surv, dead) ~ x, data = d[-45,])
summary(f)
f$loglik

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


(.packages())
if (!require(pacman)) install.packages("pacman")
pacman::p_load(tidyverse, survival, ggfortify, survminer, plotly, gridExtra, 
               Epi, KMsurv, gnm, cmprsk, mstate, flexsurv, splines, epitools, 
               eha, shiny, ctqr, scales)


require(ggfortify)

orca <- read.table("http://www.stats4life.se/data/oralca.txt")


orca <- mutate(orca, all = event != "Alive")
table(orca$all)

fit_km <- survfit(Surv(time, all) ~ 1, data = orca)
print(fit_km, print.rmean = TRUE)
 
dat_km <- fortify(fit_km)
head(dat_km)

# detach("package:Hmisc", unload=TRUE)
# detach("package:rms", unload=TRUE)
# sort((.packages()))
# 
# detach("package:plyr", unload=TRUE)
# #orca$e <- ifelse(orca$event %in% "Alive",0,1)



cuts <- seq(0, 23, 1)
lifetab_dat <- orca %>%
  mutate(time_cat = cut(time, cuts)) %>%
  group_by(time_cat) %>%
  summarise(nlost = sum(all == 0),
            nevent = sum(all == 1))
 
dat_lt <- with(lifetab_dat, lifetab(tis = cuts, ninit = nrow(orca), 
                                    nlost = nlost, nevent = nevent))
round(dat_lt, 4)
 
su_obj <- Surv(orca$time, orca$all)
str(su_obj)

 

fit_fh <- survfit(su_obj ~ 1, data = orca, type = "fleming-harrington", conf.type = "log-log")
dat_fh <- fortify(fit_fh)
## for the Nelson-Aalen estimator of the cumulative hazard
#dat_fh <- fortify(fit_fh, fun = "cumhaz")
head(dat_fh)
head(dat_km)


ggplotly(
  ggplot() +
    geom_step(data = dat_km, aes(x = time, y = surv, colour = "K-M")) +
    geom_step(data = dat_fh, aes(x = time, y = surv, colour = "N-A")) +
    geom_step(data = dat_lt, aes(x = cuts[-length(cuts)], y = surv, colour = "LT")) +
    labs(x = "Time", y = "Survival", colour = "Estimator") +
    theme_classic()
)











 














 
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
    
    library(survival)
    library(survminer)
    d <- res$foo
    
    KM_fit <- survfit(Surv(dt, e) ~ trt ,data = d)
    NA_fit <- survfit(Surv(dt, e) ~ trt, data = d,   type = "fleming-harrington")
    plot(KM_fit, conf.int = FALSE, xlab = "Time", ylab = "Survival Probability")
    lines(NA_fit, conf.int = FALSE, col = 'red')
    
    
    with(KM_fit, cbind(time, n.risk, n.event, n, surv,std.err,cumhaz))
    print(KM_fit, print.rmean = TRUE)
    dat_km <- fortify(KM_fit)
    head(dat_km)
    
    fortify(KM_fit, fun = "cumhaz")
    
    
    print(NA_fit, print.rmean = TRUE)
    dat_na <- fortify(NA_fit)
    head(dat_na)
    
     
    
    
    
    
    
    
    
    
    
    
    
    
    
  #  d$id <- 1
 #   d$dt<-d$dt*100
    
    su_obj <- Surv(d$dt, d$e)
    str(su_obj)
    
   
     
    fit_km <- survfit(su_obj ~ 1, data = d)
    print(fit_km, print.rmean = TRUE)
    
    dat_km <- fortify(fit_km)
    head(dat_km)
    
    
    
    
    
    fit_fh <- survfit(su_obj ~ 1, data = d, type = "fleming-harrington", conf.type = "log-log")
    dat_fh <- fortify(fit_fh)
    ## for the Nelson-Aalen estimator of the cumulative hazard
    #dat_fh <- fortify(fit_fh, fun = "cumhaz")
    head(dat_fh)
    head(dat_km)
    
    
    ggplotly(
      ggplot() +
        geom_step(data = dat_km, aes(x = time, y = surv, colour = "K-M")) +
        geom_step(data = dat_fh, aes(x = time, y = surv, colour = "N-A")) +
      #  geom_step(data = dat_lt, aes(x = cuts[-length(cuts)], y = surv, colour = "LT")) +
        labs(x = "Time", y = "Survival", colour = "Estimator") +
        theme_classic()
    )
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
  
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
 








