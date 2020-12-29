# master R code
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# load the required packages
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# https://thestatsgeek.com/2014/03/28/interpreting-changes-in-hazard-and-hazard-ratios/

set.seed(333) # reproducible
 
library(ggplot2)
library(dplyr)
library(directlabels)
library(Hmisc)
library(ggplot2)
library(tidyverse)
library(plotly)
library(survminer)
library(rms)
# library(scales) # For the trans_format function
# library(shinyalert)
library(DT)
library(survival)
options(max.print=1000000)    

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# function to format 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# https://stackoverflow.com/questions/3245862/format-numbers-to-significant-figures-nicely-in-r
formatz <- function(x){
  
  if (!is.na(x)  ) {
    
    formatC(signif(x,digits=5), digits=5,format="fg", flag="#",big.mark=",")
    
  }
  
}

formatz0 <- function(x){
  sprintf(x, fmt = '%s')  
}
formatz1 <- function(x){
  sprintf(x, fmt = '%#.1f')  
}
formatz2 <- function(x){
  sprintf(x, fmt = '%#.2f')  
}
formatz00 <- function(x){
  round(x,0) 
}
formatz4 <- function(x){
  sprintf(x, fmt = '%#.4f')  
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
logit <- function(p) log(1/(1/p-1))
expit <- function(x) 1/(1/exp(x) + 1)
inv_logit <- function(logit) exp(logit) / (1 + exp(logit))
is.even <- function(x){ x %% 2 == 0 } # function to identify odd maybe useful
 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
# function to create data and analyse
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
coxdata <- function(n, allocation, hr, baseline) { 
  
  #n=1000; allocation =.5; hr=2; baseline=.4
  
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
  d$dt <- sort(2*rexp(nrow(d)))# new times
  
  dx <<- datadist(d)
  options(datadist='dx')
  f0a <- cph(Surv(dt,e) ~  trt, x=TRUE, y=TRUE, data=d )
  f0a <- f0a$coefficients[[1]]
  f2 <- survfit(Surv(dt ,e)  ~ trt, data = d)
  
  #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  
  return(list(f=f, d=d, f1=f1, sf=sf, np=np, LL1=LL1, LL0=LL0, S=S,                  
              
              f0=f0,  f2=f2, f0a=f0a, foo=foo))
  
}

# dummy <- coxdata(n=1000, allocation =.5, hr=2, baseline=.4)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# function that can use the information from coxdata function to learn about the 
# behind the scenes cal in Cox PH, we allow a guess at HR and see the log likelihood
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## enter hr guess, data, dead var, trt var and treat var

loglike2 <- function(x, dat, dead, indep , time) {
  
  dd <- dat          # make another data object
  dd$dead <- dead    # take the key variables to run Cox PH
  dd$indep <- indep
  dd$time <- time
  
  ## run the analysis to get hr and log like
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
  
  dd$guess <- exp(x)
  dd$likelihoodi <- dd$part1*(dd$part2 - dd$part3)
  dd$L <- sum(dd$likelihoodi)
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  dd <- as.data.frame(dd)
  dd$dead <- dd$indep <- dd$part1 <- dd$part2 <- dd$part3 <- NULL
  
  return(dd)
  
}
 
                                                    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                                    
                                                              n          ="800" 
                                                              allocation ="0.5" 
                                                              baseline   ="0.4"
                                                              hr         ="0.75" 
                                                    
                                                    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                                   #"Change in hazard 
                                                              base ="0.03"
                                                              cens="0.02"
                                                              hr2 ="1.2"
                                                              per="0.70"   # Enter a survival probability
                                                              per2="0.50"  # Enter another survival probability' 
                                                             
                                                    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                                    #Power 
                                                               ss=c("0.7, 0.5")    #Enter two survival probs
                                                               ss2=c("11.9, 23.1") #Enter two survival times
                                                               tt=c("500,500")     # Enter ctrl, intervention n 
                                                               hrx ="1.2"          #HR
                                                               af="3"              #Enter accrual time
                                                               af2="160"           #follow up
                                                               sim="500"           #Number of simulations
                                                               t2=".1"             #Intervention non compliance

                                                    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                                    #  Parametric Survival Distributions 
                                                             
                                                            shape ="1"
                                                            scale ="0.03"
                                                             
                                                    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # This is where a new sample is instigated and inputs converted to numeric
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
    
    n <- as.numeric(n)
    allocation <-as.numeric(allocation)
    hr <- as.numeric(hr)
    baseline <- as.numeric(baseline)
    
    ###############################################
    
    base <- as.numeric(base)
    censx <- as.numeric(cens)
    hr2 <-  as.numeric(hr2)
    per <-  as.numeric(per)
    per2 <-  as.numeric(per2)
    ###############################################
     
    ss <- as.numeric(unlist(strsplit(ss,",")))
    
    ss2x <- as.numeric(unlist(strsplit(ss2,",")))
    
    tt <- as.numeric(unlist(strsplit(tt,",")))
     
    af <- as.numeric(unlist(strsplit(af,",")))  
    
    af2 <- as.numeric(unlist(strsplit(af2,","))) 
    
    hrx <- as.numeric(unlist(strsplit(hrx,",")))
    
    nonc <- as.numeric(unlist(strsplit(t2,",")))
    
    sim <- as.numeric(unlist(strsplit(sim,",")))
   
      
      ss1=ss[1]
      ss2=ss[2]
      prob1=ss2x[1]
      prob2=ss2x[2]
      nc=tt[1]
      ni=tt[2]
      AA=af[1]
      FF=af2[1]
      hrx=hrx[1]
      nonc=nonc[1]
      sim=sim[1]
    
    t1=ss1
    t2=ss2
    
    s1=prob1
    s2=prob2
    
    N1=nc
    N2=ni
    
    start=AA
    fini=FF
    
    hrx=hrx
    drop=nonc
    sim=sim
    
    
    library(Hmisc)
    
    Weib.p <- Weibull2(c(s1,s2),c(t1,t2))
    
    rcens <- function(n) runif(n, start, fini)#
    
    ff.dropout <- Quantile2(Weib.p,hratio=function(x) hrx,
                            dropout=function(x) drop)
    
    #plot(ff.dropout)
    
    rcontrol <-      function(n) ff.dropout(n, what='control')
    rintervention <- function(n) ff.dropout(n, what='intervention')
    
    x<-spower(rcontrol, rintervention, rcens, pr=FALSE,
              nc=N1, 
              ni=N2,
              test=logrank, nsim=sim, alpha=0.025, cox=TRUE)
    
   
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # get one simulation realisation and plot it for information
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    ##lifted from spower function
    
    yc <- rcontrol(N1)
    yi <- rintervention(N2)
    cens <- rcens(N1+ N2)
    group <- c(rep(0, N1), rep(1, N2))
    y <- c(yc, yi)
    maxfail <- 0
    maxcens <- 0
    maxfail <- max(maxfail, max(y))
    maxcens <- max(maxcens, max(cens))
    S <- cbind(pmin(y, cens), 1 * (y <= cens))
    nexceed <- 0 + (logrank(S, group) > .025)
    fit <- coxph.fit(as.matrix(group), S, strata = NULL, 
                     offset = NULL, init = NULL, control = coxph.control(iter.max = 10, 
                                                                         eps = 0.0001), method = "efron", 
                     rownames = NULL)
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    d <- cbind(S, group)
    d <- data.frame(d)
    names(d) <- c("dt","e", "trt")
    dd <<- datadist(d)
    options(datadist='dd')
    
    fit <- cph(Surv(dt, e) ~ trt, data = d, x=TRUE, y=TRUE, surv=TRUE)
    
    
    # survplot(fit, type="kaplan-meier", conf.int = TRUE, 
    #          col.fill=c("firebrick1","cyan2"), grid=TRUE, what='survival')
    
    
    f1 <- survfit(Surv(dt,e) ~ trt, data = d)
    
    
    fit
    f1
    d
    plot(ff.dropout)
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  
  
    
    p1 <- ggsurvplot(f1, main = "Kaplan-Meier Curve",
                     
                     
                     legend.title = "Trt."
                     
                     #,xlab=paste0("Time : HR=",round(exp(fit$coefficients),4))
                     ,xlab= "Time")
                     
   p1
   
   print(x, digits=4)
  
   
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Here we plot weibull distributions
   
    
    s1=ss1
    s2=ss2
    
    t1=prob1
    t2=prob2
    
    hr=hrx
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    Weib.p <- Weibull2(c(t1,t2),c(s1,s2))
    #Weib.p
    ff <- Quantile2(Weib.p,hratio=function(x) hr ) # we get weibull parameters 
    
    ##pull out intervention survival probs
    fff <-      attributes(ff)  #lets get the data
    time <-     fff$plot.info$`I Survival`$Time
    survival <- fff$plot.info$`I Survival`$Survival
    #plot(survival~time, lty=1)   # interventions
    
    #pull out control survival data
    timec <-     fff$plot.info$`C Survival`$Time
    survivalc <- fff$plot.info$`C Survival`$Survival
    
    ###~~~~~~~~~~~~~~~~~~~~~~~~# what time is survival at t in intevention
    
    s50i <- which.min(abs(survival-s2)) # what index is time ~ .5
    s70i <- which.min(abs(survival-s1))
    
    Weib.i <- Weibull2(c(  time[s70i], time[s50i]),c(s1, s2))
    #Weib.i
    ffi <-   Quantile2(Weib.i,hratio=function(x) 1)  # use hr of 1 here so no intervention effect
    fff <-   attributes(ffi)  #lets get the data
    timei <- fff$plot.info$`I Survival`$Time
    survivali <- fff$plot.info$`I Survival`$Survival
    
    #par(mfrow=c(2,2))
    #plot(survivalc~timec,    type = "l", lty = 1 , main ="Control arm") # ok
    
    time <-       fff$plot.info$`C Survival`$Time
    survival   <- fff$plot.info$`C Survival`$Survival
    #plot(survival~time,    type = "l", lty = 2,  main ="intervention arm") 
    #plot(ff)
    #plot(ffi)
    #par(mfrow=c(1,1))
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # control, pull out weibull parameters
    p <- lapply(Weib.p, unlist)
    t1x <- (p$alpha)
    t2x <- (p$gamma)
    
    # intervention, pull out weibull parameters
    i <- lapply(Weib.i, unlist)
    t3x <- (i$alpha)
    t4x <- (i$gamma)
    
    A <- expression( paste("control ",      alpha) )
    B <- expression( paste("control ",      gamma) )
    C <- expression( paste("intervention ", alpha) )
    D <- expression( paste("intervention ", gamma) )
    
    #dweibull(x, shape=gamma, scale = 1/alpha), from=0, to=40)
    FF <- expression( paste("Note the Weibull parameterisation shape= ",      gamma," scale=1/ ",      alpha) ) 
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~,
    plot(survivalc~timec,    type = "l", lty = 2,  ylab="Probability of survival",
         
         main =paste("Weibull distibutions, intervention HR =",hr) , col='blue', xlab= "Time",
         sub=FF)
    
    lines(survivali~timei, type = "l", lty = 1, col='red')  
    
    jump <- .85
    jump0 <-.9
    jump1 <-.76
    jump2 <-.72
    
    text(quantile(prob=jump0,c(timec,timei)),  0.96, c(round(t1x ,5)), cex = 1)
    text(quantile(prob=jump1,c(timec,timei)), 0.96, A,                cex = 1)
    
    text(quantile(prob=jump0,c(timec,timei)),  0.78, c(round(t2x ,5)), cex = 1)
    text(quantile(prob=jump1,c(timec,timei)), 0.78, B,                cex = 1)
    
    text(quantile(prob=jump0,c(timec,timei)),  0.84, c(round(t3x ,5)), cex = 1)
    text(quantile(prob=jump1,c(timec,timei)), 0.84, C,                cex = 1)
    
    text(quantile(prob=jump0,c(timec,timei)),  0.90, c(round(t4x ,5)), cex = 1)
    text(quantile(prob=jump1,c(timec,timei)), 0.90, D,                 cex = 1)
    
    s1i=survival[which.min(abs(time-t1))]
    s2i=survivali[which.min(abs(timei-t2))]
    
    text(quantile(prob=jump2,c(timec,timei)), 0.72, paste0("At time ",t1,":"),   cex = 1)
    text(quantile(prob=jump,c(timec,timei)), 0.66, paste0("Control survival prob ",s1," "),   cex = 1)
    text(quantile(prob=jump,c(timec,timei)), 0.60, paste0("Interv. survival prob ",round(s1i,2),""),   cex = 1)
    
    text(quantile(prob=jump2,c(timec,timei)), 0.54, paste0("At time ",t2,":"),   cex = 1)
    text(quantile(prob=jump,c(timec,timei)), 0.48, paste0("Control survival prob ",s2," "),   cex = 1)
    text(quantile(prob=jump,c(timec,timei)), 0.42, paste0("Interv. survival prob ",round(s2i,2)," "),   cex = 1)
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # add arrows to explain  
    arrows(                                   
      t1,                                  # x start  
      s1 ,                                 # surv prob at t1 in control
      t1 ,                                 # x finish
      s1i ,                                # surv prob at t1 in intervention     
      col="black", lty=1 )       
    
    arrows(                                   
      t2,                                    # x start  
      s2 ,                                   # surv prob at t2 in control
      t2 ,                                   # x finish
      s2i,                                   # surv prob at t2 in intervention     
      col="black", lty=1 )          
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 
    
  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # The change in  hazard tab, data generation
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
    
    n =        n
    lambdaT =  base  # event hazard
    lambdaC =  censx  # censoring haz
    beta1 =    hr2   # hr
    per =      per   # survival probability
    per2 =     per2  # survival probability
    
    #  n=800; lambdaT=14.4; lambdaC=12; beta1=1.2; per=0.7 ;  # use this for on the fly
    
    beta1 <- log(as.numeric(beta1))  # log hr
    
    x1 = sample(0:1, n,replace=TRUE)  # assign trt randomly
    
    # create distributions
    T = rweibull(n, shape=1, scale=1/lambdaT)*exp(-beta1*x1)   # say if lambda is entered as 14 big hazard need to use 1/14 here?
    C = rweibull(n, shape=1, scale=1/lambdaC)                  #censoring time
    time = pmin(T,C)  # observed time is min of censored and true, pmin useful here
    event = time==T   # set to 1 if event is observed
    
    # run cox regression
    require(survival)
    f <- coxph(Surv(time, event)~ x1 , method="breslow")
    survfit <- survfit(Surv(time,event) ~ x1)
    # f
    #plot(survfit, ylab="Survival probability", xlab="Time", col=c('blue','red'))
    #run weibull regression
    w <- survreg(formula = Surv(time, event) ~ x1, dist = "w", control = list(maxiter = 90) )
    
    # grab the hr estimates from weibull
    hrx <-  (c(w$coefficient[2],  confint(w)[2,]))
    hrx <- exp(hrx)      #exp(-coef(f))^exp(coef(f)["shape"])
    # grab the hr estimates from Cox
    hrc <- exp(c(f$coefficient,  confint(f) ) )
    
    # capture for later
    ss <- Surv(time, event)
    
    #  for practice run parametric regression using harrell package
    dd <<- datadist( x1=x1)
    options(datadist='dd')
    f.exp  <- psm(ss  ~ x1, dist ='exponential')
    fw    <-  psm(ss  ~ x1,  dist ='weibull')
    
    d <- table(x1,event)
    ev <- d[,"TRUE"]  # number of events
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # new, plotting the hazard, Dirk F moore page 38
    
    library(muhaz)
    
    # plotH
    result.simple <- muhaz(time, event, max.time=max(time),
                           bw.grid=2.25, bw.method="global", b.cor="none")
    
    
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # PlotH1
    result.pe5 <- pehaz(time, event, width=5, max.time=max(time))  # in the text these are months, so will need to adjust? 
    result.pe1 <- pehaz(time, event, width=1, max.time=max(time))
    
    result.smooth <- muhaz(time, event, bw.smooth=max(time),
                           b.cor="left", max.time=max(time))
    
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #One use of smoothing the hazard function is to obtain a smooth estimate of the
    # survival function
    
    haz <- result.smooth$haz.est
    times <- result.smooth$est.grid
    survx <- exp(-cumsum(haz[1:(length(haz)-1)]*diff(times)))
    
    # We may compare our
    # smoothed survival estimate to the Kaplan-Meier estimate as follows:
    result.km <- survfit(Surv(time, event) ~ 1,
                         conf.type="none")
    
     s=survfit 
    
  
  # right plot....................., 
  
    
    # n =        sample$n
    # lambdaT =  sample$base  # event hazard
    # lambdaC =  sample$censx  # censoring haz
    # beta1 =    sample$hr2   # hr
    # per =      sample$per   # survival probability
    # per2 =     sample$per2  # survival probability
    
    
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    yo <- abs(100*((beta1/1)-1))
    
    wordd <- ifelse(beta1 < 1,"will reduce", 
                    ifelse(beta1 > 1, "will increase","will not change")) 
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # plot weibull density
    x <- curve(dweibull(x, shape=1, scale = 1/lambdaT), from=0, to=max(s$time))  
    x$y <- x$y/max(x$y)    #scale the weibull to max of 1 
    plot(x, type = "l", lty = 1, col='blue' , ylab="Survival probability", xlab="Time" , ylim=c(0,1))
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    text(x = max(s$time)*.58, y = .99,                # Text with different color & size
         paste0(" For the blue reference curve:"),
         col = "#1b98e0",
         cex = 1)
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    text(x = max(s$time)*.655, y = .97,                # Text with different color & size
         paste0(" At ",  formatz1(-log(per)* 1/lambdaT), " months the survival probability is ",per* 100,"%"),
         col = "#1b98e0",
         cex = 1)
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    text(x = max(s$time)*.66, y = .95,                # Text with different color & size
         paste0(" At ",  formatz1(-log(per2)* 1/lambdaT), " months the survival probability is ",per2* 100,"%"),
         col = "#1b98e0",
         cex = 1)
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    text(max(s$time)*.57, 0.92, expression( paste(
      "Using S(1)t = S(0)t"^{exp(beta.x)}
    )), cex = 1)
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    text(x = max(s$time)*.7, y = .88,                # Text with different color & size
         paste0(" Postulating treatment " ,wordd," the hazard by ",yo,"%"),
         col = "red",
         cex = 1)
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    text(x = max(s$time)*.69, y = .86,                # Text with different color & size
         paste0(" At ",  formatz1(-log(per)*  1/lambdaT), " months the survival probability becomes ",formatz00(100*(per)^(beta1)) ,"%"),
         col = "red",
         cex = 1)
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    text(x = max(s$time)*.695, y = .84,                # Text with different color & size
         paste0(" At ",  formatz1(-log(per2)*  1/lambdaT), " months the survival probability becomes ",formatz00(100*(per2)^(beta1)) ,"%"),
         col = "red",
         cex = 1)
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    # with the effect of trt and plot
    y <- x$y^(beta1)   # add another line based on S1(t) = S(0) ^exp(B)
    lines(y~x$x, col='red')     
    
    # add arrows to explain  
    arrows(                                  # x start
      -log(per2) *1/lambdaT,                   # time  
      per2 ,#                                # surv prob  
      -log(per2) *1/lambdaT ,  
      per2^beta1 ,
      col="black", lty=1 )       
    
    arrows(                                  # x start
      -log(per) *1/lambdaT,    
      per, # 
      -log(per) *1/lambdaT ,  
      per^beta1 ,
      col="black", lty=1 )    
    
     
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # left plot
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  
    
    plot(s, ylab="Survival probability", xlab="Time", col=c('blue','red'))  # plot the survfit object
    
    text(x = max(s$time)*.58, y = .99,                # Text with different color & size
         paste0("HR from Weibull=",formatz2(1/hr[1])," 95%CI ( ",formatz2(1/hr[3]),", ",formatz2(1/hr[2])," )"),
         col = "#1b98e0",
         cex = 1)
    text(x = max(s$time)*.58, y = .97,                # Text with different color & size
         paste0("HR from Cox PH=",formatz2(hrc[1])," 95%CI ( ",formatz2(hrc[2]),", ",formatz2(hrc[3])," )"),
         col = "#1b98e0",
         cex = 1)
    text(x = max(s$time)*.54, y = .95,                # Text with different color & size
         paste0( "Actual N= ",sum(N), ", Actual events= ",sum(ev)),
         col = "#1b98e0",
         cex = 1)
    
    
    ## power 
    A <- 0.05 # alpha
    B <- 0.1  # beta
    
    s1=per
    s2=per^(beta1)
    
    # (-log(per)* lambdaT)  # time of s1 s2
    
    f <- ((qnorm(1-A/2))+qnorm(1-B))^2
    d1 <- (4*f)/ (log(beta1)^2)
    
    N2 <- d1/(1-(s1+s2)/2)
    
    
    text(x = max(s$time)*.667, y = .93,                # Text with different color & size
         paste0( "Alpha 2 sided ",A," Power ",1-B," events required= ",ceiling(d1), ", N= ",ceiling(N2)),
         col = "#1b98e0",
         cex = 1)
    
    
    
    
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
  
  # estimating hazard plot
     
    H=result.simple
    plot(H)
 
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
  
 
    
    H= result.pe5
    H1= result.pe1
    HS= result.smooth
    
    
    plot(H,  col="green")
    lines(H1, col='blue')
    lines(HS, col='red')
    
    
     
    
    H2= result.km
     
    haz <- result.smooth$haz.est
    times <- result.smooth$est.grid
    surv <- exp(-cumsum(haz[1:(length(haz)-1)]*diff(times)))
    
    plot(H2, conf.int=F, mark="|", xlab="Time",  # confint true throws an error
         #xlim=c(0,30), 
         ylab="Survival probability")
    lines(surv ~ times[1:(length(times) - 1)], col='blue')
    
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    alpha<- as.numeric(shape)
    
    lambda <- as.numeric(scale)
    
    weibSurv <- function(t, shape, scale) pweibull(t, shape=shape,
                                                   scale=scale, lower.tail=F)
    curve(weibSurv(x, shape=alpha, scale=1/lambda), from=0, to=2/lambda,
          ylim=c(0,1), ylab='Survival probability', xlab='Time')
    
    text(x = 2/lambda*.85, y = .95,                # Text with different color & size
         paste0(" Shape ", alpha,  ""),
         col = "blue",
         cex = 1.2)
    
    text(x = 2/lambda*.85, y = .9,                # Text with different color & size
         paste0("  Scale ", lambda ,""),
         col = "blue",
         cex = 1.2)
    
    
    text(x = 2/lambda*.85, y = .85,                # Text with different color & size
         paste0(" Mean ", formatz2(gamma(1 + 1/alpha)/lambda),  ""),
         col = "blue",
         cex = 1.2)
    
    
    
    text(x = 2/lambda*.85, y = .8,                # Text with different color & size
         paste0( " Median ",  formatz2((log(2)^(1/alpha))/lambda ),""),
         col = "blue",
         cex = 1.2)
    
    
    
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~hazard
    
    weibHaz <- function(x, shape, scale) 
      dweibull(x, shape=shape, scale=scale)/
      pweibull(x, shape=shape, scale=scale,lower.tail=F)
    
    curve(weibHaz(x, shape=alpha, scale=1/lambda), from=0, to=2/lambda,
          ylab='Hazard', xlab='Time', col="red")
    
 
  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # text below right plot in changing hazard 
 
    
    yo <- abs(100*((beta1/1)-1))
    
    wordd <- ifelse(beta1 < 1,"will reduce", 
                    ifelse(beta1 > 1, "will increase ","will not change")) 
    
    c(paste0(" 
    With a shape parameter of 1 and a scale parameter of ",lambdaT," a survival curve is created based on a Weibull distribution.
    Here we are measuring time in months.
    The time to reach median survival is equal to -log(0.5) x ",lambdaT,". 
             This equates to a median survival of ", formatz2(-log(.5)*lambdaT)," months. 
             
             Replacing 0.5 with a desired survival percentile will return the associated time. 
             Enter survival probabilities in the two boxes on the left. 
             
             The time at which the survival probability is ",per* 100,"% is ",  formatz2(-log(per)* lambdaT), " months.
             The time at which the survival probability is ",per2*100,"% is ",  formatz2(-log(per2)*lambdaT), " months"
    ))
    
    
    
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
  # another piece of text below right plot in changing hazard 
  
    # 
    # n=          sample$n
    # 
    # lambdaT =  sample$base
    # lambdaC =  sample$censx
    # beta1 =   sample$hr2
    # per = sample$per
    # per2 = sample$per2
    
    yo <- abs(100*((beta1/1)-1))
    
    wordd <- ifelse(beta1 < 1,"will reduce", 
                    ifelse(beta1 > 1, "will increase ","will not change")) 
    
    
    c(paste0("Now if we are postulating that a new treatment " ,wordd," the hazard by ",yo,"% we can use 
               the fact S1(t) = S0(t)^exp Bx, with exp Bx=", (beta1)," the 
              probability of survival at ",  formatz2(-log(per)*  lambdaT), " months now becomes ",formatz2((per)^(beta1)) ,".
          The probability of survival at ",  formatz2(-log(per2)* lambdaT), " months now becomes ",formatz2((per2)^(beta1)) ,"
          . See the arrows in the plot showing the changes in survival. 
             " ) )
     
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
   
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # GENERATE THE DATA Execute analysis
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   
    
    res <- coxdata(n, allocation, hr, baseline)
    
     d=res$d
     f=res$f
     f1=res$f1
     sf=res$sf
     np=res$np 
     LL1=res$LL1
     LL0=res$LL0
     S=res$S
       f0a=res$f0a
      f0=res$f0
      f2=res$f2     
                  
   
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # MAIN PLOT! updated with log transformation  option
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
    
    f <- f1  # Get the  obj
    
    sf  <- sf
    X <- as.numeric(as.character(sf[2,c("Effect")]))
    Y <- as.numeric(as.character(sf[2,c("Lower 0.95")]))
    Z <- as.numeric(as.character(sf[2,c("Upper 0.95")]))
   
    
    p1 <- ggsurvplot(f, main = "Kaplan-Meier Curve",
                     
                     #  text = paste("Time: ", time),
                     
                     #palette = c("orange", "purple")
                     legend.title = "Trt."
                     # , xlab= paste0("Time" )
                     ,xlab=paste0("Time : HR=",formatz4(exp(f0)) )
                     # , xlab= paste0("Time: HR = ", formatz2(X),", 95%CI( ",formatz2(Y),", ",formatz2(Z)," )" )
                     #,#conf.int = TRUE,
                     # ggtheme = theme_bw() # Change ggplot2
                     #   ,text = paste0("wt: ", round(wt), "</br></br>qsec: ", round(qsec))
    )
    ggplotly(p1[[1]] )
    
 
    
    f <- f1  # Get the survfit obj
    
    p1 <-  ggsurvplot(f, fun = "event",   main = "Cumulative proportion",  
                      legend.title = "Trt.")#
    #  palette = c("orange", "purple"))
    ggplotly(p1[[1]])
    
 
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
    
    f <- f1  # Get the survfit obj
    
    p1 <- ggsurvplot(f, fun = "cumhaz",  main = "Cumulative Hazard"  ,
                     legend.title = "Trt.")  
    
    p1$plot <- p1$plot + ggplot2::geom_abline(intercept = 0, slope = baseline,    linetype="dotted", col='red') 
    p1$plot <- p1$plot + ggplot2::geom_abline(intercept = 0, slope = baseline*hr, linetype="dotted", col='blue')
    
    
    ggplotly(p1[[1]])
    
   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     
    f <- f1  # Get the survfit obj
    
    p1 <- ggsurvplot(f, fun = "cloglog",  
                     main = "Complementary log−log" ,
                     legend.title = "Trt.")
    # legend.labs = c("0", "1")) 
    #palette = c("jco"))
    ggplotly(p1[[1]])
    
  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     
    f <- np  # Get the  data
    survdiffplot(f, col='darkgreen' , xlab= "Time")
    
    
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     
    d <- d
    
    dd <<- datadist(d)
    options(datadist='dd')
    
    fit <- cph(Surv(dt, e) ~ trt, data = d, x=TRUE, y=TRUE, surv=TRUE)
    
    survplot(fit, type="kaplan-meier", conf.int = TRUE, 
             col.fill=c("firebrick1","cyan2"), grid=TRUE, what='survival')
    
 
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 
    
    S <- Surv(d$dt, d$e)
    
    hazard.ratio.plot(d$trt, S, e=20, legendloc='ll', xlab='Time', antilog=FALSE, col='blue', smooth=TRUE,
                      ylim=c(-4,4), ylab=c("Log hazard ratio"))
    
   
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
    
    fit <- cph(Surv(dt, e) ~ trt, data = d, x=TRUE, y=TRUE, surv=TRUE)
    
    plot(cox.zph(fit, transform="identity" ), ylim=c(-4,4), ylab=c("Log hazard ratio"))
    
   
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
    
    #plot(survfit(S~ d$trt), col=c("purple", "orange"), fun="cloglog", xlab="Time", ylab="log(-log(Survival)" , lwd=3)
    survplot(f,logt=TRUE, loglog=TRUE, 
             col=c("red", "lightblue")
             
             # col=c("orange", "purple")
    )   #Check for Weibull-ness (linearity)
     
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    p1 <- ggsurvplot(f2, main = "Kaplan-Meier Curve", legend.title = "Trt.",
                     # palette = c("orange", "purple")  ,
                     xlab=paste0("Time : HR=", formatz4(exp(f0a)))
                     # ggtheme = theme_bw() # Change ggplot2 theme
    )
    ggplotly(p1[[1]])
    
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   
    
    d <- plyr::arrange(d, dt)
    
    DT::datatable(d, rownames=FALSE,
                  plugins = 'natural',
                  colnames=c('Time' = 'dt', 'Event or censored' = 'e', 
                             'Treatment'='trt'),
                  
                  options = list(
                    #  dom = 't',
                    columnDefs = list(list(type = 'natural', targets = c(1,2)))
                  )
    ) %>%
      
      formatRound(
        columns= c("Time" ), 
        digits=4 )
 
  
  #~~~~~~~~~~~~~~~~~~~~~~KM table
   
    
    library(survminer)
    
    require(ggfortify)
     
    KM_fit <- survfit(Surv(dt, e) ~ trt ,data = d)
    
    KM <- fortify(KM_fit) # fortify(KM_fit, fun = "cumhaz")'
    
    DT::datatable(KM, rownames=FALSE,
                  plugins = 'natural',
                  #   colnames=c('Time' = 'dt', 'Event or censored' = 'e', 
                  #              'Treatment'='trt'),
                  
                  options = list(
                    #  dom = 't',
                    columnDefs = list(list(type = 'natural', targets = c(1,2)))
                  )
    ) %>%
      
      formatRound(
        columns= c("time" ,"surv","std.err","upper","lower"), 
        digits=4 )
 
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
    require(ggfortify)
     
    KM_fit <- survfit(Surv(dt, e) ~ trt ,data = d)
    
    KM <-  fortify(KM_fit, fun = "cumhaz")
    
    DT::datatable(KM, rownames=FALSE,
                  plugins = 'natural',
                  #   colnames=c('Time' = 'dt', 'Event or censored' = 'e', 
                  #              'Treatment'='trt'),
                  
                  options = list(
                    #  dom = 't',
                    columnDefs = list(list(type = 'natural', targets = c(1,2)))
                  )
    ) %>%
      
      formatRound(
        columns= c("time" ,"surv","std.err","upper","lower"), 
        digits=4 )
 
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~maximum likelihood
  
 # Get the  data
    
    y <- as.numeric(as.character(res$f$coefficients))
    guess=y  # we use the actual model hr estimate
    
    d <- plyr::arrange(d, dt)
    
    
    d$expB <- exp(guess*d$trt)
    d$part2 <- guess*d$trt
    d$part3 <- log(rev(cumsum(rev(d$expB))))
    
    d$likelihoodi <- d$e*(d$part2 - d$part3)
    d$LL <- sum(d$likelihoodi)
    
    datatable(d, rownames=FALSE,
              plugins = 'natural',
              colnames=c('Time' = 'dt', 
                         'Event or censored' = 'e', 
                         'Treatment'='trt',
                         'HR'='expB',
                         'Individual likelihoods' ='likelihoodi',
                         'logHR x trt'='part2',
                         'Log(exp HR) of each risk set'='part3',
                         'Sum the Individual likelihoods to give log likelihood' ='LL'
              ),
              
              options = list(
                #  dom = 't',
                columnDefs = list(list(type = 'natural', targets = c(1,2)))
              )
    ) %>%
      
      formatRound(
        columns= c("Time","HR", 
                   #"A",
                   "logHR x trt","Log(exp HR) of each risk set",'Individual likelihoods'), 
        digits=4 ) %>%
      formatRound(
        columns= c("Sum the Individual likelihoods to give log likelihood"), 
        digits=1 )
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
     
    dummy = res
    
    g <- log(as.numeric(res$f$coefficients))
    
    foo <- loglike2(g, dat=dummy$d, dead=dummy$d$e, indep=dummy$d$trt, time=dummy$d$dt)  # try
    foo$time =NULL
    foo$expB = NULL
    
    datatable(foo, rownames=FALSE,
              plugins = 'natural',
              colnames=c('Time' = 'dt', 
                         'Event or censored' = 'e', 
                         'Treat.'='trt',
                         'Model Hazard Ratio'='hr',
                         'Null Log Likelihood'='lognull',
                         'Maximised Log Likelihood'='lognmax',
                         'HR guess' ='guess',
                         'Individual likelihoods' ='likelihoodi',
                         'Sum the Individual likelihoods to give log likelihood based on guess' ='L'
              ),
              
              options = list(
                #  dom = 't',
                columnDefs = list(list(type = 'natural', targets = c(1,2)))
              )
    ) %>%
      
      formatRound(
        columns= c("Time","Model Hazard Ratio",  
                   'Individual likelihoods'
        ),
        digits=4 ) %>%
      formatRound(
        columns= c( 'Null Log Likelihood',
                    'Maximised Log Likelihood',
                    'Sum the Individual likelihoods to give log likelihood based on guess'), 
        digits=0 )
    
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
    
    # allocation=sample$allocation
    # hr=        sample$hr
    # baseline=  sample$baseline
    
    res <- coxdata(n=10, allocation=allocation, hr=hr, baseline =baseline)
    
    d <- res$d
    d <- plyr::arrange(d, dt)  # sort by time
    
    # Calculate Li for everyone
    d$Numerator   <- exp(res$f$coefficients[[1]] * d$trt)
    d$Denominator <- (rev(cumsum(rev(d$Numerator))))
    d$Li          <- d$Numerator/d$Denominator
    
    # all censored contribute 1 (on multiplicative scale)
    d$Li<- ifelse(d$e %in% 1,d$Li,1)
    
    # get the product of all and log answer
    d$LL <- log(prod(d$Li))  
    d
    # model LL, prove we have ecalc correctly
    res$f$loglik
    
    
    datatable(d, rownames=FALSE,
              plugins = 'natural',
              colnames=c('Time' = 'dt', 
                         'Event or censored' = 'e', 
                         'Treat.'='trt',
                         'Num.'='Numerator',
                         'Den.'='Denominator',
                         'Individual likelihoods'='Li',
                         'log of product of Individual likelihoods to give maximizes log likelihood' ='LL'
                         
              ),
              
              options = list(
                #  dom = 't',
                columnDefs = list(list(type = 'natural', targets = c(1,2)))
              )
    ) %>%
      
      formatRound(
        columns= c("Time","Num.","Den.", 
                   'Individual likelihoods', 'log of product of Individual likelihoods to give maximizes log likelihood'
        ),
        digits=4 ) %>%
      formatRound(
        columns= c( 'Event or censored',
                    'Treat.'), 
        digits=0 )
  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   
    X <- as.numeric(as.character(sf[2,c("Effect")]))
    Y <- as.numeric(as.character(sf[2,c("Lower 0.95")]))
    Z <- as.numeric(as.character(sf[2,c("Upper 0.95")]))
    
    Xp  <- X/(X+1)
    Yp  <- Y/(Y+1)
    Zp  <- Z/(Z+1)
    
    wordup <- ifelse(X>1,"higher", "")
    
    wordup2 <- ifelse(X>1,"increase", "reduction")
    
    paste0( "From the Cox model the estimated hazard ratio is "
            , formatz2(X),", 95%CI ( ",formatz2(Y),", ",formatz2(Z),
            " ) comparing treatment 1 to 0. 
            
             A hazard ratio of  ", formatz2(X)," means that, in each unit of time, someone 
            treated in group 1 has ", formatz00(abs(X/1-1)*100),"% ", wordup ," of the chance of experiencing the event of interest
            in the following unit of time as they would were they taking treatment 0.
            
            There is an estimated ", formatz00(abs(X/1-1)*100),"% ", wordup2 ," in the hazard of the outcome by the Cox model. 
            
            Equivalently, the hazard ratio is equal to the odds that a patient in treatment group 1 experiences the event of interest before a
            a patient in treatment group 0.
            
           Therefore we can reformulate the hazard ratio, possibly more intuitively, as
            the probability that a patient in treatment 
            group 1 experiences the event before a patient in treatment group 0, which is: "
            , formatz2(Xp),", 95%CI ( ",formatz2(Yp),", ",formatz2(Zp),").        ")
    
  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # frank Harrell rms page 479
  
    
    d <-  d
    
    trt <-  d$trt
    e   <-  d$e
    dt  <-  d$dt
    d   <-  d$d
    
    limx <- quantile(dt, prob=.99)
    limx <- max(dt)*.8
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    S <- Surv(dt,e)
    f <- npsurv(S ~ trt)
    
    for (meth in c('exact','breslow','efron')) {
      
      g <- cph(S  ~ trt, method=meth, surv=TRUE, x=TRUE, y=TRUE)
      
    }
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    f.exp  <- psm(S  ~ trt, dist ='exponential')
    fw    <-  psm(S  ~ trt,  dist ='weibull')
    phform <- pphsm(fw)
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    co <- gray(c(0,.8))
    co <- c("red", "lightblue")
    survplot(f, lty=c(1,1)   , lwd=c(1,3), col=co,           label.curves=FALSE, conf='none')
    survplot(g, lty=c(3,3)   , lwd=c(1,3), col=co, add=TRUE, label.curves=FALSE, conf.type='none')
    
    legend(c(limx,160),c(.38,.99),
           c('Nonparametric estimates', 'Cox-Breslow estimates'),
           lty=c(1,3), bty='n',    cex=1.0) # col=co 
    
    legend(c(limx,160),c(.18,.89), 
           c('Trt 0','Trt 1'), lwd=c(1,3), col=co, bty='n',  cex=1.0)
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

 