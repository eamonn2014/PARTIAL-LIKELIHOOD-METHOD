 
# library(survival)
# tm <- c(0, # birth
#         1/365, # first day of life
#         7/365, # seventh day of life
#         28/365, # fourth week of life
#         1:110) # subsequent years
# 
# 
# 
# 
# hazMale <- survexp.us[,"male","2004"] # 2004 males
# hazFemale <- survexp.us[,"female","2004"] # 2004 females
# 
# 
# plot(hazMale ~ tm)
# 
# 

# 
# x <- 1:9; 
# names(x) <- x
# # Multiplication & Power Tables
# x %o% x
# 
# 
# y <- 2:8; 
# names(y) <- paste(y,":", sep = "")
# outer(y, x, "^")
#https://www.statology.org/plot-exponential-distribution-in-r/
   #~~~~~~~~~~~~~~~~~~~~~~#GOOD#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    rm(list=ls())
    
    n <- 100
    # pdf, two versions of the same thing
    lambda <- .23                              # this will be a constant hazard
    end <- ceiling(-(log(1-.999)/ lambda))     # for plotting 99th percentile
    s <- seq(0,end, length.out = end+1)        # reasonable x axis
    log(2)/lambda                              # true median survival
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    # 1 PDF. lets plot using the built in exponential function
    z <- curve(dexp(x, rate = lambda), from=0, to=end, n=end+1, col='blue') 
      
    
    z$x  # n=155 0    to  154
    z$y  # n=155 0.03 to  0.0002955839
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # 2 now program PDF up, this is of length 155
    x <- lambda*exp(-lambda*s)                    # pdf for exponential distribution
    plot(x, type = "l", lty = 1, col='red')           # plot constant hazard effect
    x
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # let me scale my programmed plot:  pdf/ hazard = survival function
    plot(x/lambda~s, type = "l", lty = 1, col='red', 
         main=paste0("s(t), rate = ",lambda,", median = ",round(log(2)/lambda,2),""))    # scaling the y axis to max 1
    abline(v=-log(.5)/lambda )
    abline(h=.5 )
    
    # just for interest
    x                           # just for interest
   # lambda*exp(-lambda*end)     # just for interest, last point
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # let's use the weibull function to plot exponential PDF
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
    plot(dweibull(s, shape=1, scale = 1/lambda), type = "l", lty = 1, col='red', main=paste0("rate = ",lambda,"")) # note scale use reciprocal of rate!
    
    # a function to plot survival function
    weibSurv <- function(t, shape, scale) pweibull(t, shape=shape, scale=scale, lower.tail=F)
    
    curve(weibSurv(x, shape=1, scale=1/lambda), from=0, to=end, n=end+1,
          ylim=c(0,1), ylab='Survival probability', xlab='Time', main=paste0("S(t), rate =",lambda,""))
    abline(v=-log(.5)/lambda )
    abline(h=.5 )
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~hazard = PDF/ Survival function~~~~~~~~~~~~~~~~~
    weibHaz <- function(x, shape, scale) 
       dweibull(x, shape=shape, scale=scale)/
       pweibull(x, shape=shape, scale=scale,lower.tail=F)
    
    curve(weibHaz(x, shape=1, scale=1/lambda), from=0, to=end, main=paste0("H(t), rate = ",lambda,""),
          ylab='Hazard', xlab='Time', col="red")
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Survival function~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    zp <- function(x, shape, scale)
    pweibull(x, shape=shape, scale=scale,lower.tail=F)
    
    curve(zp(x, shape=1, scale=1/lambda), from=0, to=end, main=paste0("S(t), rate = ",lambda,""),
          ylab='Survival probability', xlab='Time', col="red")
    abline(v=-log(.5)/lambda )
    abline(h=.5 )
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~CDF
    zd <- function(x, shape, scale)
       dweibull(x, shape=shape, scale=scale)
    
    curve(zd(x, shape=1, scale=1/lambda), from=0, to=end, main=paste0("rate =",lambda,""),
          ylab='cdf', xlab='Time', col="red")
    
    #~~~~~~~~~ HAZARD CDF/PDF
     D <-   dweibull(x, shape=1, scale=1/lambda)    # pdf
     P <-   pweibull(x, shape=1, scale=1/lambda, lower.tail=F)  # survival
     D/P 
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
     
    #multiplying the scale parameter by 7 is equivalent to dividing the day input by 7, 
    #thereby translating any day input into its equivalent week input (which is what the data is fitted on)
    
# https://stats.stackexchange.com/questions/105881/how-to-simulate-survival-times-using-true-base-line-hazard-function    
# The way to do this is by generating uniform RVs u as quantiles and finding S^−1(u). 
# This can be done analytically, or below I have an example of how to do it numerically with a 
# pseudocontinuous or discrete time using colSums(outer(x,y,'<')) which beats quantile by many flops
     


    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    h <-  rep(lambda, length(s))
    # Cumulative hazard
    ch <- cumsum(h)   # multiply by lambdaT to scale up
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    plot(ch~s)  # cum hazard
    # survival
    S <- c(1,exp(-ch ))   # put a 1 first
    S
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    s1<- c(s, length(s))
    s1
    # compare to P
    P
    
    z <- plot(S~s1, type = "l", lty = 1, col='darkgreen', lwd=2, main=paste0("S(t) rate = ",lambda,"")) ## bit off?
    abline(v=-log(.5)/lambda )
    abline(h=.5 )
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    require(survival)
    options(digits=6)
    # generate n random samples:

    u <- runif(n)
    failtimes <- s[colSums(outer(S, u, `>`))]
    km <- (survfit(Surv(failtimes)~1))
    km
    plot(km,  main=paste0("S(t), rate = ",lambda,""))
    
    curve(weibSurv(x, shape=1, scale=1/lambda), from=0, to=end, n=end+1, col='red', lwd=2, lty=2,
          
          ylim=c(0,1), ylab='Survival probability', xlab='Time', 
          add=TRUE)
    abline(v=-log(.5)/lambda )
    
    abline(h=.5 )
    
    
    
   
    
    
    
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # another go EXPONENTIAL ONLY
    
    rm(list=ls())
    require(survival)
    # a function to plot survival function
    weibSurv <- function(t, shape, scale) pweibull(t, shape=shape, scale=scale, lower.tail=F)
    
    un =100                  ## number of fail times
    # pdf, two versions of the same thing
    lambda <- .003                            # rate
    end <- ceiling(-(log(1-.999)/ lambda))   # for plotting x axis uper limit 
    s <- seq(0,end, length.out = end+1)
    s1 <- c(s, length(s))
    haz <- rep(lambda, length(s)) #
    cumhaz <- c(0,cumsum(haz) )   #  # get that first point of 0
    Surv <- exp(-cumhaz ) 
    
    par(mfrow=c(1,5))
      
       plot(s, haz,    type='l', xlab='Time', ylab='Hazard', main=paste0("h(t),rate = ",lambda,""))
      
       plot(s1, cumhaz, type='l', xlab='Time', ylab='Cumulative hazard', main=paste0("H(t) , rate = ",lambda,""))
       
       plot(s1, Surv,   type='l', xlab='Time', ylab='Survival', 
            main=paste0("S(t), rate = ",lambda,"\nMed Surv=", round(-log(.5)/lambda,2) ))
       abline(v=-log(0.5)/lambda )
       abline(h=0.5 )
       
       plot(s1, log(Surv),   type='l', xlab='Time', ylab='log Survival', main=paste0("S(t), rate = ",lambda,""))
       abline(v=-log(0.5)/lambda )
       abline(h= log(.5))
  
    # generate n random samples:
       
       ###############################
       # A get n random probabilities! (u)
       # B get the Survival times.
       # C check each probability for all survival prob, is the survival prob > probability 
       # D Sum each C
       # select the survival corresponding to each D
       ###############################
       
       u <- runif(un)
       failtimes <- s[colSums(outer(Surv, u, `>`))]
       km <- (survfit(Surv(failtimes)~1))
       km
       plot(km  , main=paste0("rate =",lambda,"\nN=",un),  ylab='Survival probability', xlab='Time')
       curve(weibSurv(x, shape=1, scale=1/lambda), from=0, to=end, n=end+1, 
              col='red', lwd=2, lty=2,
             ylim=c(0,1), add=TRUE)
       abline(v=-log(.5)/lambda )
       abline(h=.5 )
       -log(.5)/lambda 
       
    par(mfrow=c(1,1))
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~repeat but use weibull
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~weibull
    rm(list=ls())
    require(survival)
    
    # a function to plot survival function
    weibSurv <- function(t, shape, scale) pweibull(t, shape=shape, scale=scale, lower.tail=F)
     
    n         <- 100 # use this to simulate fail times
    # Weibull parameters
    lambda    <- .03 # scale   
    k         <- 2#.5  # shape < 1 decreasing with time.., > 1 increasing with time, 1 constant
    
    #~~~~~~~~~~~~~~~~~
    p         <- 0.5
    (med.surv <- 1/lambda* (-log(1-p))^(1/k) )   # true median survival
   
    p <-.95   # 95th percentile
    end <- (med.surv2 <- 1/lambda* (-log(1-p))^(1/k) )  # upper a axis limit
    #~~~~~~~~~~~~~~~~~~
    
    
    x <- seq(0, ceiling(end), length.out = end+1)   # time domain
    S =  pweibull(x, shape=k, scale=1/lambda, lower.tail=F)   # pdf
    # plot(S)
    
    u <- runif(n)
    failtimes <- x[colSums(outer(S, u, `>`))]
    km <- (survfit(Surv(failtimes)~1))
    km
    
    
    #~~~~~~~~~~~~~~~~~~~~~start plot
    par(mfrow=c(1,4))
    
    #~~~~~~~~~~~~
    plot(km  , 
         main=paste0("Shape = ",k,", rate = ",lambda,"\nTrue med surv = ",
                     round(med.surv,1) ,", N=", (n)),  
         ylab='Survival probability', xlab='Time')
    
    curve(weibSurv(x, shape=k, scale=1/lambda), from=0, to=end, n=end+1, 
          col='red', lwd=2, lty=2,
          ylim=c(0,1), add=TRUE)
    abline(v= med.surv , lty=2)
    abline(h=.5 , lty=2 )
    #~~~~~~~~~~~~~~~~~~ log scale
    
    plot(km  , log = TRUE,
         main=paste0("Shape =",k,", rate = ",lambda,"\nTrue med surv = ",
                     round(med.surv,1) ,", N=", (n)),  
         ylab='Survival probability (log scale)', xlab='Time')
   
    curve(weibSurv(x, shape=k, scale=1/lambda), from=0, to=end, n=end+1, log=TRUE,
          col='red', lwd=2, lty=2, #yaxt='n',
           ylim=c(0,1), add=TRUE)
     
    abline(v=  (med.surv) , lty=2)
    abline(h= .5, lty=2 )
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # hazard
    cdf <-   dweibull(x, shape=k, scale=1/lambda)                # pdf
    pdf <-   pweibull(x, shape=k, scale=1/lambda, lower.tail=F)  # survival function
    haz <- cdf/pdf 
    
    plot(x, haz,    type='l', xlab='Time', 
         ylab='Hazard', 
         main=paste0("h(t) Shape = ",k,", rate = ",lambda,""))
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    
    plot(x[-1], cumsum(haz[-1]),  type='l', xlab='Time', 
         ylab='Cumulative Hazard', 
         main=paste0("H(t) Shape = ",k,", rate = ",lambda,""))
    
    
    par(mfrow=c(1,1))
    #~~~~~~~~~~~~~~~~~~~~~end plot
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~repeat but use weibull
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~weibull
 
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    # # 1 PDF. lets plot using the built in exponential function
    # z <- curve(dexp(x, rate = lambda), from=0, to=end, n=end+1, col='blue') 
    # 
    # 
    # z$x  # n=155 0    to  154
    # z$y  # n=155 0.03 to  0.0002955839
    # 
    # # 2 now program CDF up, this is of length 155
    # x <- lambda*exp(-lambda*s)                    # pdf for exponential distribution
    # plot(x, type = "l", lty = 1, col='red')           # plot constant hazard effect
    # x
    # 
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
  
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~weibull
 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#https://stats.stackexchange.com/questions/105881/how-to-simulate-survival-times-using-true-base-line-hazard-function

   # t <- 2
   # h <- 0.03
   # 
   # tdom <- seq(0, t,   by=0.01)
   # haz <-  rep(0, length(tdom))
   # 
   # haz[tdom >= 0] <- h
   # cumhaz <- cumsum(haz)
   # Surv   <- exp(-cumhaz)
   # 
   # par(mfrow=c(1,4))
   #    plot(tdom, haz,    type='l',    xlab='Time', ylab='Hazard')
   #    plot(tdom, cumhaz, type='l',    xlab='Time', ylab='Cumulative hazard')
   #    plot(tdom, Surv,   type='l',    xlab='Time', ylab='Survival')
   #    u <- runif(100)   # generate 100 random samples:
   #    failtimes <- tdom[colSums(outer(Surv, u, `>`))]
   #    plot(survfit(Surv(failtimes)~1))
   # par(mfrow=c(1,1))
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 
# h <-  runif(1,1,10)  # 1/sample(1:9999,1) # (runif(1,1,10))               # sample a random constant hazard
# h <-   (0.03)
# h                                   # print h
# tdom <- seq(0, exp(h)/100, by=0.001)    # time, exp(h)/2 is just a convenient scaling factor, nothing more
# tdom <- seq(0, 80, by=0.1)    # time, exp(h)/2 is just a convenient scaling factor, nothing more
# 
# haz <-  rep(h, length(tdom))
# 
# cumhaz <- cumsum(haz)
# Surv <-  exp(-cumhaz)
# 
# par(mfrow=c(1,3))
# plot(tdom, haz,    type='l', xlab='Time', ylab='Hazard',            main=paste0("Hazard = ", round(exp(-h),5)))
# plot(tdom, cumhaz, type='l', xlab='Time', ylab='Cumulative hazard', main="Cumulative hazard")
# plot(tdom, Surv,   type='l', xlab='Time', ylab='Survival',          main="Survival Function", ylim=c(0,1))
# par(mfrow=c(1,1)) 


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# https://stats.stackexchange.com/questions/60238/intuition-for-cumulative-hazard-function-survival-analysis
    rm(list=ls())
    require(survival)
    options(scipen=999)
# data from here@
# https://www.ons.gov.uk/peoplepopulationandcommunity/birthsdeathsandmarriages/lifeexpectancies/datasets/nationallifetablesunitedkingdomreferencetables
# 2017-2019 pre covid UK deaths aged 1:100 males
dx <- c( 426.7,   24.2,   13.1,   10.0,    9.7,   8.5,   8.8,   6.8,
   6.7,   5.9,   7.4,   8.5,   10.4,   12.6,   12.1,   17.3,
   22.6,   31.4,   39.8,   44.5,   50.4,   50.9,   50.0,   50.0,
   54.7,   59.2,   57.5,   61.4,   68.7,   72.3,   76.1,   82.3,
   82.3,   91.5,   94.1,   105.5,   112.5,   128.3,   123.5,   136.8,
   149.3,   162.8,   176.9,   194.6,    203.0,   227.2,   240.2,   259.7,
   273.9,   305.4,   322.6,   343.1,   370.4,   389.5,   421.1,   445.7,
   502.2,   542.8,   589.1,   634.4,   690.5,   750.0,   825.2,   902.3,
   960.5,   1043.5,   1142.9,   1221.3,    1312.3,   1418.8,   1474.7,   1605.3,
   1731.7,   1933.7,   2078.2,   2255.3,    2442.8,   2606.9,   2807.8,   2968.1,
   3170.9,   3338.1,   3473.9,   3676.0,    3820.8,   3905.8,   3991.4,   3956.3,
   3921.6,   3818.9,   3456.3,   3260.5,   2944.2,   2581.5,   2243.5,   1874.0,
   1521.2,   1150.9,   858.3,   656.0,    430.2)
  
   x <- 0:(length(dx)-1) # age vector ( this could be time )
   n=100
   haz <-    (dx/sum(dx)) / (1-cumsum(dx/sum(dx)))
   cumhaz <- cumsum(haz)
   Surv <-   exp(-cumhaz)

   plot(haz,   t='l', log="y",
        xlab="age", ylab="Annual ‘hazard’ or ‘force of mortality h(t)", 
        main="Annual risk of death in males from all causes for UK 2017-19")  
   #lines(dd$x,  exp(predict(f)), col = 'red',  lwd=2, lty=2)
   
   
# subset data and run regression for linear estimate~~~~~~~~~~~~~~~~~
   
   require(rms)
   z <- 30:90   # these ages only
   dd<- data.frame(cbind(y= haz[z],x=z))
   dx <- datadist(dd)
   options(datadist='dx')
   
   f <- ols(log(y) ~ x , dd)
   y <- haz[1:100]
   x <- 1: 100
   # slope
   exp(f$coefficients[2])

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   plot(x,  y, t='l', log="y",
        xlab="age", ylab="Annual ‘hazard’ or ‘force of mortality h(t)", 
        main="Annual risk of death in males from all causes for UK 2017-19")  
   lines(dd$x,  exp(predict(f)), col = 'red',  lwd=2, lty=2)

par(mfrow=c(1,4))
#plot(haz,    t="l", xlab="age", ylab="Annual ‘hazard’ or ‘force of mortality h(t)", 
#main="Annual risk of death \nmales from all causes \nfor UK 2017-19", log="y")

   plot(x,  y, t='l', log="y",
        xlab="age", ylab="Annual ‘hazard’ or ‘force of mortality h(t)", 
        main="Annual risk of death in \nmales from all causes \nfor UK 2017-19")  
   lines(dd$x,  exp(predict(f)), col = 'red',  lwd=2, lty=2)
   
   plot(cumhaz, t="l", xlab="age", ylab="H(t)", main="H(t)")
   
   plot( Surv,  t="l", xlab="age", ylab="S(t)", main="S(t)")
   
   u <- runif(n)
   failtimes <- x[colSums(outer(Surv, u, `>`))]
   km <- (survfit(Surv(failtimes)~1))
   km
   plot(km, main=paste0("Survival, N=",length(u)),  
   ylab='Survival probability', xlab='age')

par(mfrow=c(1,1))



# subset data and run regression
# require(rms)
#  
# dd<- data.frame(cbind(y= haz[26:90],x=26: 90))
# dx <- datadist(dd)
# options(datadist='dx')
# 
# f <- ols(log(y) ~ x , dd)
# y <- haz[1:100]
# x <- 1: 100
# 
# 
#  plot(x,  log(y), t='l' ) 
#  lines(dd$x, predict(f), col = 'red',  lwd=2, lty=2)
# 
# 
# plot(x,  y, t='l', log="y")  
# lines(dd$x,  exp(predict(f)), col = 'red',  lwd=2, lty=2)

 




# Annual risk of death from all causes for uk,
# 2017–2019, known as the annual ‘hazard’ or ‘force of mortality’


# curve(weibSurv(x, shape=1, scale=1/lambda), from=0, to=end, n=end+1, 
#       main=paste0("rate =",lambda,""), col='red', lwd=2, lty=2,
#       ylim=c(0,1), ylab='Survival probability', xlab='Time', add=TRUE)
# abline(v=-log(.5)/lambda )
# abline(h=.5 )
# -log(.5)/lambda 
 




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# # Moore's book
# par(mfrow=c(1,2)) 
# 
# weibSurv <- function(t, shape, scale) { pweibull(t, shape=shape, scale=scale, lower.tail=F) }
# 
# curve(weibSurv(x, shape=1, scale=1/lambda), from=0, to=80,  # enter reciprocal for 0.03
#       ylim=c(0,1), ylab='Survival probability', xlab='Time')
# 
# # hazard = pdf / survival
# weibHaz <- function(x, shape, scale) { dweibull(x, shape=shape,  scale=scale)/
#                                        pweibull(x, shape=shape,  scale=scale, lower.tail=F) }
# 
# curve(weibHaz(x, shape=1, scale=1/lambda), from=0, to=80, 
#       ylab='Hazard', xlab='Time', col="red")
# par(mfrow=c(1,1)) 
# 
# 



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

inverse = function(fn, min_x, max_x){
   # Returns the inverse of a function for a given range.
   # E.g. inverse(sin, 0, pi/2)(sin(pi/4)) equals pi/4 because 0 <= pi/4 <= pi/2
   fn_inv = function(y){
      uniroot((function(x){fn(x) - y}), lower=min_x, upper=max_x)[1]$root
   }
   return(Vectorize(fn_inv))
}

integrate_from_0 = function(fn, t){
   int_fn = function(t) integrate(fn, 0, t)
   result = sapply(t, int_fn)
   value  = unlist(result["value",])
   msg    = unlist(result["message",])
   value[which(msg != "OK")] = NA
   return(value)
}

random_survival_times = function(hazard_fn, n, max_time=10000){
   # Given a hazard function, returns n random time-to-event observations.
   cumulative_density_fn = function(t) 1 - exp(-integrate_from_0(hazard_fn, t))
   inverse_cumulative_density_fn = inverse(cumulative_density_fn, 0, max_time)
   return(inverse_cumulative_density_fn(runif(n)))
}

h = lambda ; n = 100

hazard_fn = function(t) rep(h, length(t))

survival_times = random_survival_times(hazard_fn, n)
#And let's check that the Kaplan-Meier curve for these survival times approximates, 
#as expected, the curve P(t) = exp(-0.001t):

plot(survfit(Surv(survival_times, rep(1, length(survival_times)))~1), xlab="Time", 
     ylab="Survival Probability",
     main=paste0("Sampled and Expected Survival Curves for h(t) = ",h,""))

neg_exp_fn = function(x){exp(-h * x)}
curve(expr=neg_exp_fn, from=0, to=max(survival_times), add=TRUE, col="red", lwd=2)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


km <-  survfit(Surv(survival_times, rep(1, length(survival_times)))~1)
plot(km)                                    # basic KM plot

 
 -log(0.5)/h # true median survival in reference group
 
# print observed median survival
 
# plot(km$cumhaz)
# 
# prop.fail <- km$n.event/km$n.risk
# # Work out the length of time over which these failure occur:
# time <- km$time
# time0 <- c(0, time[-length(time)])  # ad 0 and remove last 
# 
# 
# #Divide prop.fail by the time interval over which those failures occur
# #(that is, time - time0) to get the probability of failing
# #per unit time, i.e. the instantaneous hazard:
# haz <- prop.fail/(time - time0)
# # Plot the result:
# plot(time, haz, #ylim = c(0,0.9), 
#      type = "s", xlab = "Days to relapse", ylab =
#         "h(t)")
# lines(lowess(time[-1], haz[-1], f = 0.10))
# 
# 
# # Tidier plot:
# plot(time, haz, type = "n", xlab = "Days to relapse", ylab = "h(t)", 
#      ylim =      c(0,1))
# lines(lowess(time[-1], haz[-1], f = 0.10))
# abline(h=h, lty=2)








# https://dk81.github.io/dkmathstats_site/rvisual-cont-prob-dists.html

# Plotting Three Weibull Distributions:
# https://dk81.github.io/dkmathstats_site/rvisual-cont-prob-dists.html
library(gridExtra)
library(grid)
library(ggplot2)
xx= 1.3  # line size

x_lower_wei <- 0
x_upper_wei <- 80

pa<-ggplot(data.frame(x = c(x_lower_wei , x_upper_wei)), aes(x = x)) + 
   xlim(c(x_lower_wei , x_upper_wei)) + 
   stat_function(fun = pweibull, args = list(shape = 1.0, scale = 1/.03, lower.tail=F), size = xx, aes(colour = "1.0 & 0.03")) + 
   stat_function(fun = pweibull, args = list(shape = 1.0, scale = 1/.04, lower.tail=F), size = xx, aes(colour = "1.0 & 0.04")) + 
   stat_function(fun = pweibull, args = list(shape = 0.75, scale = 1/.03, lower.tail=F), size = xx, aes(colour = "0.75 & 0.03")) + 
   stat_function(fun = pweibull, args = list(shape = 1.5, scale = 1/.03, lower.tail=F), size = xx, aes(colour = "1.5 & 0.03")) + 
   scale_color_manual("Shape & Scale \n Parameters", values = c("blue", "black", "red", "darkgreen")) +
   labs(x = "\n t", y = "S(t) \n", 
        title = "Survival Distribution Plots") + 
   theme(plot.title = element_text(hjust = 0.5), 
         axis.title.x = element_text(face="bold", colour="blue", size = 12),
         axis.title.y = element_text(face="bold", colour="blue", size = 12),
         legend.title = element_text(face="bold", size = 10),
         legend.position = "right") + theme_bw() # + theme(legend.position = "none")  


 
# see pweibull help, cum haz  H(t) = - log(1 - F(t)) is

p <- {function(x, shape, scale) 
   -pweibull(x, shape=shape, scale=scale, lower.tail=FALSE, log=TRUE)}

pb<-ggplot(data.frame(x = c(x_lower_wei , x_upper_wei)), aes(x = x)) + 
   xlim(c(x_lower_wei , x_upper_wei)) + 
   stat_function(fun = p, args = list(shape = 1.0, scale = 1/.03 ), size = xx, aes(colour = "1.0 & 0.03")) + 
   stat_function(fun = p, args = list(shape = 1.0, scale = 1/.04 ), size = xx, aes(colour = "1.0 & 0.04")) + 
   stat_function(fun = p, args = list(shape = 0.75, scale =1/.03 ), size = xx, aes(colour = "0.75 & 0.03")) + 
   stat_function(fun = p, args = list(shape = 1.5, scale = 1/.03 ), size = xx, aes(colour = "1.5 & 0.03")) + 
   scale_color_manual("Shape & Scale \n Parameters", values = c("blue", "black", "red", "darkgreen")) +
   labs(x = "\n t", y = "H(t) \n", 
        title = "Cumulative hazard distribution Plots") + 
   theme(plot.title = element_text(hjust = 0.5), 
         axis.title.x = element_text(face="bold", colour="blue", size = 12),
         axis.title.y = element_text(face="bold", colour="blue", size = 12),
         legend.title = element_text(face="bold", size = 10),
         legend.position = "right") + theme_bw() #+ theme(legend.position = "none")  
  # theme(legend.position="bottom") 


# function to calculate hazard:
weibHaz <- {function(x, shape, scale) dweibull(x, shape=shape,
                   scale=scale)/pweibull(x, shape=shape, scale=scale, lower.tail=F)}


pc<-ggplot(data.frame(x = c(x_lower_wei , x_upper_wei)), aes(x = x )) + 
   xlim(c(x_lower_wei , x_upper_wei)) + 
   stat_function(fun = weibHaz, args = list(shape = 1.0, scale = 1/.03), size = xx, aes(colour = "1.0 & 0.03")) + 
   stat_function(fun = weibHaz, args = list(shape = 1.0, scale = 1/.04), size = xx,aes(colour = "1.0 & 0.04")) + 
   stat_function(fun = weibHaz, args = list(shape = 0.75, scale = 1/.03), size = xx,aes(colour = "0.75 & 0.03")) + 
   stat_function(fun = weibHaz, args = list(shape = 1.5, scale = 1/.03), size = xx,aes(colour = "1.5 & 0.03")) + 
   scale_color_manual("Shape & Scale \n Parameters", values = c("blue", "black", "red", "darkgreen")) +
   labs(x = "\n t", y = "h(t) \n", 
        title = "Hazard Distribution Plots") + 
  
   theme(plot.title = element_text(hjust = 0.5), 
         axis.title.x = element_text(face="bold", colour="blue", size = 12),
         axis.title.y = element_text(face="bold", colour="blue", size = 12),
         legend.title = element_text(face="bold", size = 10),
         legend.position = "right") +   theme_bw()   

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



# f2 <- npsurv(S~1)
# n  <- f2$numevents
# x  <- f2$exposure
# f2$numevents\f2$exposure
# 
# 
# 
# grid.arrange(pa,  pc, pb, ncol=3,
#              
#              top = textGrob(paste0( "Weibull"),
#                             gp=gpar(fontsize=20,font=3)) )
# 
# 
# 
# 










## Moore's book p.17

windows()
par(mfrow=c(1,2))

      up= 80
      scale. = 0.03
      A= 1.50
      B= 1.00
      C= 0.75
      D= 0.04
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~Hazard   CDF/PDF 
weibHaz <- {function(x, shape, scale) dweibull(x, shape=shape,
                  scale=scale)/pweibull(x, shape=shape, scale=scale, lower.tail=F)}

      curve(weibHaz(x, shape=A, scale=1/scale.), from=0, to=up, 
         ylab='Hazard', xlab='Time', col="red", 
         main=paste0("The hazard is monotone increasing for ",
         expression(alpha)," >1, monotone decreasing for ",
                     expression(alpha),"<1 \nand constant when equal to 1"))
 
       
      
      curve(weibHaz(x, shape=B, scale=1/scale.), from=0, to=up, 
         ylab='Hazard', xlab='Time', add=T, col="black")
   
      curve(weibHaz(x, shape=C, scale=1/scale.), from=0, to=up, 
         ylab='Hazard', xlab='Time', add=T, col="blue")
      
      curve(weibHaz(x, shape=B, scale=1/D), from=0, to=up, 
            ylab='Hazard', xlab='Time', add=T, col="darkgreen",  lty=2)
   
      text(45, 0.065, bquote("" ~ alpha == .(A) ~ ""  ~ lambda == .(scale.) ~ ""), 
           col="red", cex=1.3)
      
      text(45, 0.034, bquote("" ~ alpha == .(B) ~ ""  ~ lambda == .(scale.) ~ ""), 
           col="black", cex=1.3)
      
      text(45, 0.015, bquote("" ~ alpha == .(C) ~ "" ~ lambda == .(scale.) ~ ""), 
           col="blue", cex=1.3)
   
      text(45, 0.044, bquote("" ~ alpha == .(B) ~ ""  ~ lambda == .(D) ~ ""), 
           col="darkgreen", cex=1.3)
    
#~~~~~~~~~~~~~~~~~~~~~~~~~Survival      PDF
weibHaz <- {function(x, shape, scale)  pweibull(x, shape=shape, scale=scale, lower.tail=F)}

   curve(weibHaz(x, shape=A, scale=1/scale.), from=0, to=up, 
         ylab='Survival', xlab='Time', col="red")
   
   curve(weibHaz(x, shape=B, scale=1/scale.), from=0, to=up, 
         ylab='Survival', xlab='Time', add=T, col="black")
   
   curve(weibHaz(x, shape=C, scale=1/scale.), from=0, to=up, 
         ylab='Survival', xlab='Time', add=T, col="blue")
   
   curve(weibHaz(x, shape=B, scale=1/D), from=0, to=up, 
         ylab='Hazard', xlab='Time', add=T, col="darkgreen",  lty=2)
   
   text(16, 0.065, bquote("" ~ alpha == .(A) ~ ""  ~ lambda == .(scale.) ~ ""), 
        col="red", cex=1.3)
   
   text(15, 0.23, bquote("" ~ alpha == .(B) ~ ""  ~ lambda == .(scale.) ~ ""), 
        col="black", cex=1.3)
   
   text(17, 0.15, bquote("" ~ alpha == .(C) ~ "" ~ lambda == .(scale.) ~ ""), 
        col="blue", cex=1.3)
   
   text(15, 0.31, bquote("" ~ alpha == .(B) ~ ""  ~ lambda == .(D) ~ ""), 
        col="darkgreen", cex=1.3)

par(mfrow=c(1,1))
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
weibHaz <- {function(x, shape, scale)  pweibull(x, shape=shape, scale=scale, lower.tail=F, log=TRUE)}

curve(weibHaz(x, shape=A, scale=1/scale.), from=0, to=up, 
      ylab='Survival', xlab='Time', col="red",
      main=paste0("When alpha is 1 the log scale plot will be a straight line"))  

 curve(weibHaz(x, shape=B, scale=1/scale.), from=0, to=up, 
      ylab='Survival', xlab='Time', add=T, col="black")

curve(weibHaz(x, shape=C, scale=1/scale.), from=0, to=up, 
      ylab='Survival', xlab='Time', add=T, col="blue")

curve(weibHaz(x, shape=B, scale=1/D), from=0, to=up, 
      ylab='Hazard', xlab='Time', add=T, col="darkgreen",  lty=2)

text(16, 0.065-2.4, bquote("" ~ alpha == .(A) ~ ""  ~ lambda == .(scale.) ~ ""), 
     col="red", cex=1.3)
text(15, 0.23-2.85, bquote("" ~ alpha == .(B) ~ ""  ~ lambda == .(scale.) ~ ""), 
     col="black", cex=1.3)

text(16.5, 0.155-3.1, bquote("" ~ alpha == .(C) ~ "" ~ lambda == .(scale.) ~ ""), 
     col="blue", cex=1.3)

text(15, 0.31-3.6, bquote("" ~ alpha == .(B) ~ ""  ~ lambda == .(D) ~ ""), 
     col="darkgreen", cex=1.3)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~






# Inverse cumulative hazard function

#########stata and bender p1716 Generating survival times to simulate Cox proportional
# hazards models

# Times = H^-1[-log(U)]
#  
# 
# uniformly distributed random numbers can be transformed into survival
# times following a specific Cox model
# It is just required to insert the inverse of an appropriate cumulative baseline
# hazard function into equation


#Inverse cumulative hazard function

   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   # In addition to f(), F(), S(), h(), and H() -all different ways of summarizing the same
   # information-a sixth way, Q(), is of use for to create artificial survival-time data sets.
   # if U ∼ U[0; 1], then (1 − U) ∼ U[0; 1] so dispense with 1- in weibull formula
   # uniformly distributed random numbers can be transformed into survival times
   n=100
   # function to plot true survival function
   weibSurv <- function(t, shape, scale) pweibull(t, shape=shape, scale=scale, lower.tail=F)
   
   par(mfrow=c(1,3))  
   
      # simulate weibull shape parameter p, scale = h
      p=1/3    # shape, if 1 this is exponential
      h=0.05   # scale
      
      #~~~~~~~~~~~~~~~~~~~~~~~~~~
      t <- h^-1*(-log(runif(n)))^p # weibull
      
      survfit <- survfit(Surv(t) ~ 1)
      plot(survfit, ylab="Survival probability", xlab="Time", main=paste0("N = ",n,", Weibull, shape = ",1/p,", rate = ",h))
      
      # plot true survival curve for weibull
      curve(weibSurv(x, shape=1/p, scale=1/h), from=0, to=max(t), n=length(t), 
            col='red', lwd=2, lty=2,
            ylim=c(0,1), add=TRUE)
      
      #~~~~~~~~~~~~~~~~~~~~~~~~~~
      t <- -log(runif(n))/h # exponential rate = h 
      
      survfit <- survfit(Surv(t) ~ 1)
      plot(survfit, ylab="Survival probability", xlab="Time", main=paste0("N = ",n,", Exponential, rate = ",h))
    
      # plot true survival curve, constant hazard h
      curve(weibSurv(x, shape=1, scale=1/h), from=0, to=max(t), n=length(t), 
            col='red', lwd=2, lty=2,
            ylim=c(0,1), add=TRUE)


   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   
   require(flexsurv)
   gSurv <- function(t, shape, scale) pgompertz(t, shape=p, rate = h, lower.tail=F)
   
   t <- (1/p)*log((-log(runif(n)))*(p/h)+1) # gompertz see bender table1
   
   survfit <- survfit(Surv(t) ~ 1)
   plot(survfit, ylab="Survival probability", xlab="Time", main=paste0("N = ",n,", Gompertz, shape = ",p,", scale = ",h))
   
   # plot true survival curve, 
   curve(gSurv(x, shape=p, scale=1/h), from=0, to=max(t), n=length(t), 
         col='red', lwd=2, lty=2,
         ylim=c(0,1), add=TRUE)
   
   par(mfrow=c(1,1))
   
   
   
   
   
   #  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   #  # clear environment
   #  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   # 
   #  rm(list=ls())       
   #  options(digits=5)
   #  set.seed(777)
   #  
   #  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   #  # libraries
   #  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   #  
   #  library(survival)   # no need to call as part of base R
   #  library(ggfortify)  # better plotting functions
   #  library(rms)
   #  # beta<- 0.002
   #  # plot(exp(-1*c(seq(.1,10, .1))))
   #  # 
   #  # plot(exp(-2*c(seq(.1,10, .1))))
   #  
   # 
   #  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   #  # simulate censored survival data with one categorical predictor
   #  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   #  
   #  n  = 4000        # total sample size
   #  HR = 2          # true hazard ratio
   #  beta1 = log(HR) # log hazard ratio 
   #  lambdaT = 1/.03  # baseline hazard 
   #  lambdaC = 1/.05  # hazard of censoring
   #  
   #  x1 <- sample(0:1, n, replace=TRUE)                          # treatment groups random assignment
   #  Tx <- rweibull(n, shape=1, scale=lambdaT*exp(-beta1*x1))    # exponential, true event time for each trt. grp.
   #  Cx <- rweibull(n, shape=1, scale=lambdaC)                   # censoring ,
   #  time <- pmin(Tx,Cx)                                         # observed time is min of censored and true
   #  event <- time==Tx                                           # set to 1 if event is observed
   # 
   #  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   #  # KM plots, create survival object and plot 3 different ways
   #  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   #  S <- Surv(time, event)
   #  
   #  
   #  km1 <- survfit(S  ~ x1, conf.type="none")    # survival object
   #  plot(km1)                                    # basic KM plot
   #  
   #  km <- survfit(S  ~ x1 )                      # survival object
   #  autoplot(km)                                 # KM plot...many more options 
   # 
   #  np <- npsurv(S  ~ x1)                        # yet another way, plotly 
   #  survplotp(np, time.inc=lambdaT, times=c(.000693,.000693*2))
   # 
   #  
   #  # nice surv diff plot
   #  survdiffplot(km, col='darkgreen' , xlab= "Time")
   #  
   #  
   #  # check proportional hazards
   # 
   #  dd <- datadist(time, event, x1)
   #  options(datadist='dd')
   #  fit <- cph(S  ~ x1, x=TRUE, y=TRUE, surv=TRUE)  # rms package
   #  
   #  #~formal
   #  cox.zph(fit)
   #  
   #  #~eyeball
   #  plot(cox.zph(fit, transform="identity" ), ylim=c(-4,4), ylab=c("Log hazard ratio"))
   #  
   #  #~eyeball
   #  hazard.ratio.plot(as.numeric(x1), S, e=20, legendloc='ll', xlab='Time', antilog=FALSE, col='blue', smooth=TRUE,
   #                    ylim=c(-4,4), ylab=c("Log hazard ratio"))
   #  abline(h=beta1,lty=2, col='purple')
   #  
   #  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   #  # median survival, of course ratio of median survival between groups is the HR
   #  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   #  # We know in truth median survival :
   #  # The probability of surviving x months or longer after starting treatment was 0.5
   # 
   #  (m0 = -log(0.5)*lambdaT*exp(-beta1*0))   # true median survival in reference group
   #  (m1 = -log(0.5)*lambdaT*exp(-beta1*1))   # true median survival in treated group
   # 
   #  m0/m1
   #  # print observed median survival
   #  km
   #  
   # #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   # # print all data and survival probabilities and extract to data frame
   # #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   # 
   # (res <- fortify(km))  
   # 
   # (summary(km1))  # print events only
   # 
   # summary(km, times = c(59)) # also enter a time of interest if so wish
   # 
   # #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   # # log rank test 
   # #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   # 
   #  survdiff(Surv(time, event) ~ x1)  
   # 
   # #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   # # Cox ph 
   # #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   #   
   #  f <- coxph(Surv(time, event)~ x1 , method="breslow")
   #  summary(f)
   #  exp(confint(f))
   #  
   # #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   # # hazards
   #  
   # plot(km$cumhaz)
   # 
   # prop.fail <- km$n.event/km$n.risk
   # # Work out the length of time over which these failure occur:
   # time <- km$time
   # time0 <- c(0, time[-length(time)])  # ad 0 and remove last 
   # 
   # 
   # #Divide prop.fail by the time interval over which those failures occur
   # #(that is, time - time0) to get the probability of failing
   # #per unit time, i.e. the instantaneous hazard:
   # haz <- prop.fail/(time - time0)
   # # Plot the result:
   # plot(time, haz, #ylim = c(0,0.9), 
   #      type = "s", xlab = "Days to relapse", ylab =
   #          "h(t)")
   # lines(lowess(time[-1], haz[-1], f = 0.10))
   # 
   # 
   # # Tidier plot:
   #   plot(time, haz, type = "n", xlab = "Days to relapse", ylab = "h(t)", 
   #     ylim =      c(0,.1))
   # lines(lowess(time[-1], haz[-1], f = 0.10))
   # 
   #  
   # 
   # 
   # # master R code
   # #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   # # load the required packages
   # #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   # 
   # # https://thestatsgeek.com/2014/03/28/interpreting-changes-in-hazard-and-hazard-ratios/
   # 
   # set.seed(333) # reproducible
   # 
   # library(ggplot2)
   # library(dplyr)
   # library(directlabels)
   # library(Hmisc)
   # library(ggplot2)
   # library(tidyverse)
   # library(plotly)
   # library(survminer)
   # library(rms)
   # # library(scales) # For the trans_format function
   # # library(shinyalert)
   # library(DT)
   # library(survival)
   # options(max.print=1000000)    
   # 
   # #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   # # function to format 
   # #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   # 
   # # https://stackoverflow.com/questions/3245862/format-numbers-to-significant-figures-nicely-in-r
   # formatz <- function(x){
   #   
   #   if (!is.na(x)  ) {
   #     
   #     formatC(signif(x,digits=5), digits=5,format="fg", flag="#",big.mark=",")
   #     
   #   }
   #   
   # }
   # 
   # formatz0 <- function(x){
   #   sprintf(x, fmt = '%s')  
   # }
   # formatz1 <- function(x){
   #   sprintf(x, fmt = '%#.1f')  
   # }
   # formatz2 <- function(x){
   #   sprintf(x, fmt = '%#.2f')  
   # }
   # formatz00 <- function(x){
   #   round(x,0) 
   # }
   # formatz4 <- function(x){
   #   sprintf(x, fmt = '%#.4f')  
   # }
   # 
   # #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   # logit <- function(p) log(1/(1/p-1))
   # expit <- function(x) 1/(1/exp(x) + 1)
   # inv_logit <- function(logit) exp(logit) / (1 + exp(logit))
   # is.even <- function(x){ x %% 2 == 0 } # function to identify odd maybe useful
   # 
   # #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
   # # function to create data and analyse
   # #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   # #
   # coxdata <- function(n, allocation, hr, baseline) { 
   #   
   #   #n=1000; allocation =.5; hr=2; baseline=.4
   #   
   #   trt <- sample(0:1, n,  rep=TRUE, prob=c(1-allocation, allocation))
   #   
   #   cens <- 15*runif(n)
   #   
   #   h <- baseline*exp(log(hr)*(trt==1))  # hazard function h(t)
   #   
   #   dt <- -log(runif(n))/h
   #   
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
   #   foo <-d
   #   # S <- Surv(dt,e)
   #   f <- cph(Surv(dt,e) ~  trt, x=TRUE, y=TRUE, data=d )
   #   f0 <- f$coefficients[[1]] #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
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
   #   #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
   #   
   #   d <- plyr::arrange(d,dt)
   #   d$dt <- dd<- NULL
   #   d$dt <- sort(2*rexp(nrow(d)))# new times
   #   
   #   dx <<- datadist(d)
   #   options(datadist='dx')
   #   f0a <- cph(Surv(dt,e) ~  trt, x=TRUE, y=TRUE, data=d )
   #   f0a <- f0a$coefficients[[1]]
   #   f2 <- survfit(Surv(dt ,e)  ~ trt, data = d)
   #   
   #   #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
   #   
   #   return(list(f=f, d=d, f1=f1, sf=sf, np=np, LL1=LL1, LL0=LL0, S=S,                  
   #               
   #               f0=f0,  f2=f2, f0a=f0a, foo=foo))
   #   
   # }
   # 
   # # dummy <- coxdata(n=1000, allocation =.5, hr=2, baseline=.4)
   # 
   # #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   # # function that can use the information from coxdata function to learn about the 
   # # behind the scenes cal in Cox PH, we allow a guess at HR and see the log likelihood
   # #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   # 
   # ## enter hr guess, data, dead var, trt var and treat var
   # 
   # loglike2 <- function(x, dat, dead, indep , time) {
   #   
   #   dd <- dat          # make another data object
   #   dd$dead <- dead    # take the key variables to run Cox PH
   #   dd$indep <- indep
   #   dd$time <- time
   #   
   #   ## run the analysis to get hr and log like
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
   #   
   #   dd$guess <- exp(x)
   #   dd$likelihoodi <- dd$part1*(dd$part2 - dd$part3)
   #   dd$L <- sum(dd$likelihoodi)
   #   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   #   dd <- as.data.frame(dd)
   #   dd$dead <- dd$indep <- dd$part1 <- dd$part2 <- dd$part3 <- NULL
   #   
   #   return(dd)
   #   
   # }
   # 
   # #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   # 
   # n          ="800" 
   # allocation ="0.5" 
   # baseline   ="0.4"
   # hr         ="0.75" 
   # 
   # #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   # #"Change in hazard 
   # base ="0.03"
   # cens="0.02"
   # hr2 ="1.2"
   # per="0.70"   # Enter a survival probability
   # per2="0.50"  # Enter another survival probability' 
   # 
   # #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   # #Power 
   # ss=c("0.7, 0.5")    #Enter two survival probs
   # ss2=c("11.9, 23.1") #Enter two survival times
   # tt=c("500,500")     # Enter ctrl, intervention n 
   # hrx ="1.2"          #HR
   # af="3"              #Enter accrual time
   # af2="160"           #follow up
   # sim="500"           #Number of simulations
   # t2=".1"             #Intervention non compliance
   # 
   # #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   # #  Parametric Survival Distributions 
   # 
   # shape ="1"
   # scale ="0.03"
   # 
   # #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   # #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   # # This is where a new sample is instigated and inputs converted to numeric
   # #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   # 
   # 
   # n <- as.numeric(n)
   # allocation <-as.numeric(allocation)
   # hr <- as.numeric(hr)
   # baseline <- as.numeric(baseline)
   # 
   # ###############################################
   # 
   # base <- as.numeric(base)
   # censx <- as.numeric(cens)
   # hr2 <-  as.numeric(hr2)
   # per <-  as.numeric(per)
   # per2 <-  as.numeric(per2)
   # ###############################################
   # 
   # ss <- as.numeric(unlist(strsplit(ss,",")))
   # 
   # ss2x <- as.numeric(unlist(strsplit(ss2,",")))
   # 
   # tt <- as.numeric(unlist(strsplit(tt,",")))
   # 
   # af <- as.numeric(unlist(strsplit(af,",")))  
   # 
   # af2 <- as.numeric(unlist(strsplit(af2,","))) 
   # 
   # hrx <- as.numeric(unlist(strsplit(hrx,",")))
   # 
   # nonc <- as.numeric(unlist(strsplit(t2,",")))
   # 
   # sim <- as.numeric(unlist(strsplit(sim,",")))
   # 
   # 
   # ss1=ss[1]
   # ss2=ss[2]
   # prob1=ss2x[1]
   # prob2=ss2x[2]
   # nc=tt[1]
   # ni=tt[2]
   # AA=af[1]
   # FF=af2[1]
   # hrx=hrx[1]
   # nonc=nonc[1]
   # sim=sim[1]
   # 
   # t1=ss1
   # t2=ss2
   # 
   # s1=prob1
   # s2=prob2
   # 
   # N1=nc
   # N2=ni
   # 
   # start=AA
   # fini=FF
   # 
   # hrx=hrx
   # drop=nonc
   # sim=sim
   # 
   # 
   # library(Hmisc)
   # 
   # Weib.p <- Weibull2(c(s1,s2),c(t1,t2))
   # 
   # rcens <- function(n) runif(n, start, fini)#
   # 
   # ff.dropout <- Quantile2(Weib.p,hratio=function(x) hrx,
   #                         dropout=function(x) drop)
   # 
   # #plot(ff.dropout)
   # 
   # rcontrol <-      function(n) ff.dropout(n, what='control')
   # rintervention <- function(n) ff.dropout(n, what='intervention')
   # 
   # x<-spower(rcontrol, rintervention, rcens, pr=FALSE,
   #           nc=N1, 
   #           ni=N2,
   #           test=logrank, nsim=sim, alpha=0.025, cox=TRUE)
   # 
   # 
   # 
   # #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   # # get one simulation realisation and plot it for information
   # #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   # #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
   # ##lifted from spower function
   # 
   # yc <- rcontrol(N1)
   # yi <- rintervention(N2)
   # cens <- rcens(N1+ N2)
   # group <- c(rep(0, N1), rep(1, N2))
   # y <- c(yc, yi)
   # maxfail <- 0
   # maxcens <- 0
   # maxfail <- max(maxfail, max(y))
   # maxcens <- max(maxcens, max(cens))
   # S <- cbind(pmin(y, cens), 1 * (y <= cens))
   # nexceed <- 0 + (logrank(S, group) > .025)
   # fit <- coxph.fit(as.matrix(group), S, strata = NULL, 
   #                  offset = NULL, init = NULL, control = coxph.control(iter.max = 10, 
   #                                                                      eps = 0.0001), method = "efron", 
   #                  rownames = NULL)
   # 
   # #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
   # d <- cbind(S, group)
   # d <- data.frame(d)
   # names(d) <- c("dt","e", "trt")
   # dd <<- datadist(d)
   # options(datadist='dd')
   # 
   # fit <- cph(Surv(dt, e) ~ trt, data = d, x=TRUE, y=TRUE, surv=TRUE)
   # 
   # 
   # # survplot(fit, type="kaplan-meier", conf.int = TRUE, 
   # #          col.fill=c("firebrick1","cyan2"), grid=TRUE, what='survival')
   # 
   # 
   # f1 <- survfit(Surv(dt,e) ~ trt, data = d)
   # 
   # 
   # fit
   # f1
   # d
   # plot(ff.dropout)
   # #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   # 
   # 
   # 
   # 
   # p1 <- ggsurvplot(f1, main = "Kaplan-Meier Curve",
   #                  
   #                  
   #                  legend.title = "Trt."
   #                  
   #                  #,xlab=paste0("Time : HR=",round(exp(fit$coefficients),4))
   #                  ,xlab= "Time")
   # 
   # p1
   # 
   # print(x, digits=4)
   # 
   # 
   # 
   # #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   # 
   # #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   # # Here we plot weibull distributions
   # 
   # 
   # s1=ss1
   # s2=ss2
   # 
   # t1=prob1
   # t2=prob2
   # 
   # hr=hrx
   # 
   # #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   # 
   # Weib.p <- Weibull2(c(t1,t2),c(s1,s2))
   # #Weib.p
   # ff <- Quantile2(Weib.p,hratio=function(x) hr ) # we get weibull parameters 
   # 
   # ##pull out intervention survival probs
   # fff <-      attributes(ff)  #lets get the data
   # time <-     fff$plot.info$`I Survival`$Time
   # survival <- fff$plot.info$`I Survival`$Survival
   # #plot(survival~time, lty=1)   # interventions
   # 
   # #pull out control survival data
   # timec <-     fff$plot.info$`C Survival`$Time
   # survivalc <- fff$plot.info$`C Survival`$Survival
   # 
   # ###~~~~~~~~~~~~~~~~~~~~~~~~# what time is survival at t in intevention
   # 
   # s50i <- which.min(abs(survival-s2)) # what index is time ~ .5
   # s70i <- which.min(abs(survival-s1))
   # 
   # Weib.i <- Weibull2(c(  time[s70i], time[s50i]),c(s1, s2))
   # #Weib.i
   # ffi <-   Quantile2(Weib.i,hratio=function(x) 1)  # use hr of 1 here so no intervention effect
   # fff <-   attributes(ffi)  #lets get the data
   # timei <- fff$plot.info$`I Survival`$Time
   # survivali <- fff$plot.info$`I Survival`$Survival
   # 
   # #par(mfrow=c(2,2))
   # #plot(survivalc~timec,    type = "l", lty = 1 , main ="Control arm") # ok
   # 
   # time <-       fff$plot.info$`C Survival`$Time
   # survival   <- fff$plot.info$`C Survival`$Survival
   # #plot(survival~time,    type = "l", lty = 2,  main ="intervention arm") 
   # #plot(ff)
   # #plot(ffi)
   # #par(mfrow=c(1,1))
   # 
   # #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   # # control, pull out weibull parameters
   # p <- lapply(Weib.p, unlist)
   # t1x <- (p$alpha)
   # t2x <- (p$gamma)
   # 
   # # intervention, pull out weibull parameters
   # i <- lapply(Weib.i, unlist)
   # t3x <- (i$alpha)
   # t4x <- (i$gamma)
   # 
   # A <- expression( paste("control ",      alpha) )
   # B <- expression( paste("control ",      gamma) )
   # C <- expression( paste("intervention ", alpha) )
   # D <- expression( paste("intervention ", gamma) )
   # 
   # #dweibull(x, shape=gamma, scale = 1/alpha), from=0, to=40)
   # FF <- expression( paste("Note the Weibull parameterisation shape= ",      gamma," scale=1/ ",      alpha) ) 
   # #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~,
   # plot(survivalc~timec,    type = "l", lty = 2,  ylab="Probability of survival",
   #      
   #      main =paste("Weibull distibutions, intervention HR =",hr) , col='blue', xlab= "Time",
   #      sub=FF)
   # 
   # lines(survivali~timei, type = "l", lty = 1, col='red')  
   # 
   # jump <- .85
   # jump0 <-.9
   # jump1 <-.76
   # jump2 <-.72
   # 
   # text(quantile(prob=jump0,c(timec,timei)),  0.96, c(round(t1x ,5)), cex = 1)
   # text(quantile(prob=jump1,c(timec,timei)), 0.96, A,                cex = 1)
   # 
   # text(quantile(prob=jump0,c(timec,timei)),  0.78, c(round(t2x ,5)), cex = 1)
   # text(quantile(prob=jump1,c(timec,timei)), 0.78, B,                cex = 1)
   # 
   # text(quantile(prob=jump0,c(timec,timei)),  0.84, c(round(t3x ,5)), cex = 1)
   # text(quantile(prob=jump1,c(timec,timei)), 0.84, C,                cex = 1)
   # 
   # text(quantile(prob=jump0,c(timec,timei)),  0.90, c(round(t4x ,5)), cex = 1)
   # text(quantile(prob=jump1,c(timec,timei)), 0.90, D,                 cex = 1)
   # 
   # s1i=survival[which.min(abs(time-t1))]
   # s2i=survivali[which.min(abs(timei-t2))]
   # 
   # text(quantile(prob=jump2,c(timec,timei)), 0.72, paste0("At time ",t1,":"),   cex = 1)
   # text(quantile(prob=jump,c(timec,timei)), 0.66, paste0("Control survival prob ",s1," "),   cex = 1)
   # text(quantile(prob=jump,c(timec,timei)), 0.60, paste0("Interv. survival prob ",round(s1i,2),""),   cex = 1)
   # 
   # text(quantile(prob=jump2,c(timec,timei)), 0.54, paste0("At time ",t2,":"),   cex = 1)
   # text(quantile(prob=jump,c(timec,timei)), 0.48, paste0("Control survival prob ",s2," "),   cex = 1)
   # text(quantile(prob=jump,c(timec,timei)), 0.42, paste0("Interv. survival prob ",round(s2i,2)," "),   cex = 1)
   # 
   # #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   # # add arrows to explain  
   # arrows(                                   
   #   t1,                                  # x start  
   #   s1 ,                                 # surv prob at t1 in control
   #   t1 ,                                 # x finish
   #   s1i ,                                # surv prob at t1 in intervention     
   #   col="black", lty=1 )       
   # 
   # arrows(                                   
   #   t2,                                    # x start  
   #   s2 ,                                   # surv prob at t2 in control
   #   t2 ,                                   # x finish
   #   s2i,                                   # surv prob at t2 in intervention     
   #   col="black", lty=1 )          
   # 
   # #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   # #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   # 
   # 
   # 
   # 
   # #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   # #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   # #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   # #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   # # The change in  hazard tab, data generation
   # #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   # 
   # 
   # n =        n
   # lambdaT =  base  # event hazard
   # lambdaC =  censx  # censoring haz
   # beta1 =    hr2   # hr
   # per =      per   # survival probability
   # per2 =     per2  # survival probability
   # 
   # #  n=800; lambdaT=14.4; lambdaC=12; beta1=1.2; per=0.7 ;  # use this for on the fly
   # 
   # beta1 <- log(as.numeric(beta1))  # log hr
   # 
   # x1 = sample(0:1, n,replace=TRUE)  # assign trt randomly
   # 
   # # create distributions
   # T = rweibull(n, shape=1, scale=1/lambdaT)*exp(-beta1*x1)   # say if lambda is entered as 14 big hazard need to use 1/14 here?
   # C = rweibull(n, shape=1, scale=1/lambdaC)                  #censoring time
   # time = pmin(T,C)  # observed time is min of censored and true, pmin useful here
   # event = time==T   # set to 1 if event is observed
   # 
   # # run cox regression
   # require(survival)
   # f <- coxph(Surv(time, event)~ x1 , method="breslow")
   # survfit <- survfit(Surv(time,event) ~ x1)
   # # f
   # #plot(survfit, ylab="Survival probability", xlab="Time", col=c('blue','red'))
   # #run weibull regression
   # w <- survreg(formula = Surv(time, event) ~ x1, dist = "w", control = list(maxiter = 90) )
   # 
   # # grab the hr estimates from weibull
   # hrx <-  (c(w$coefficient[2],  confint(w)[2,]))
   # hrx <- exp(hrx)      #exp(-coef(f))^exp(coef(f)["shape"])
   # # grab the hr estimates from Cox
   # hrc <- exp(c(f$coefficient,  confint(f) ) )
   # 
   # # capture for later
   # ss <- Surv(time, event)
   # 
   # #  for practice run parametric regression using harrell package
   # dd <<- datadist( x1=x1)
   # options(datadist='dd')
   # f.exp  <- psm(ss  ~ x1, dist ='exponential')
   # fw    <-  psm(ss  ~ x1,  dist ='weibull')
   # 
   # d <- table(x1,event)
   # ev <- d[,"TRUE"]  # number of events
   # 
   # #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   # # new, plotting the hazard, Dirk F moore page 38
   # 
   # library(muhaz)
   # 
   # # plotH
   # result.simple <- muhaz(time, event, max.time=max(time),
   #                        bw.grid=2.25, bw.method="global", b.cor="none")
   # 
   # 
   # 
   # #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   # # PlotH1
   # result.pe5 <- pehaz(time, event, width=5, max.time=max(time))  # in the text these are months, so will need to adjust? 
   # result.pe1 <- pehaz(time, event, width=1, max.time=max(time))
   # 
   # result.smooth <- muhaz(time, event, bw.smooth=max(time),
   #                        b.cor="left", max.time=max(time))
   # 
   # 
   # #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   # #One use of smoothing the hazard function is to obtain a smooth estimate of the
   # # survival function
   # 
   # haz <- result.smooth$haz.est
   # times <- result.smooth$est.grid
   # survx <- exp(-cumsum(haz[1:(length(haz)-1)]*diff(times)))
   # 
   # # We may compare our
   # # smoothed survival estimate to the Kaplan-Meier estimate as follows:
   # result.km <- survfit(Surv(time, event) ~ 1,
   #                      conf.type="none")
   # 
   # s=survfit 
   # 
   # 
   # # right plot....................., 
   # 
   # 
   # # n =        sample$n
   # # lambdaT =  sample$base  # event hazard
   # # lambdaC =  sample$censx  # censoring haz
   # # beta1 =    sample$hr2   # hr
   # # per =      sample$per   # survival probability
   # # per2 =     sample$per2  # survival probability
   # 
   # 
   # 
   # #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   # yo <- abs(100*((beta1/1)-1))
   # 
   # wordd <- ifelse(beta1 < 1,"will reduce", 
   #                 ifelse(beta1 > 1, "will increase","will not change")) 
   # #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   # # plot weibull density
   # x <- curve(dweibull(x, shape=1, scale = 1/lambdaT), from=0, to=max(s$time))  
   # x$y <- x$y/max(x$y)    #scale the weibull to max of 1 
   # plot(x, type = "l", lty = 1, col='blue' , ylab="Survival probability", xlab="Time" , ylim=c(0,1))
   # #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   # text(x = max(s$time)*.58, y = .99,                # Text with different color & size
   #      paste0(" For the blue reference curve:"),
   #      col = "#1b98e0",
   #      cex = 1)
   # #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   # text(x = max(s$time)*.655, y = .97,                # Text with different color & size
   #      paste0(" At ",  formatz1(-log(per)* 1/lambdaT), " months the survival probability is ",per* 100,"%"),
   #      col = "#1b98e0",
   #      cex = 1)
   # #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   # text(x = max(s$time)*.66, y = .95,                # Text with different color & size
   #      paste0(" At ",  formatz1(-log(per2)* 1/lambdaT), " months the survival probability is ",per2* 100,"%"),
   #      col = "#1b98e0",
   #      cex = 1)
   # #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   # 
   # text(max(s$time)*.57, 0.92, expression( paste(
   #   "Using S(1)t = S(0)t"^{exp(beta.x)}
   # )), cex = 1)
   # 
   # #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   # text(x = max(s$time)*.7, y = .88,                # Text with different color & size
   #      paste0(" Postulating treatment " ,wordd," the hazard by ",yo,"%"),
   #      col = "red",
   #      cex = 1)
   # #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   # text(x = max(s$time)*.69, y = .86,                # Text with different color & size
   #      paste0(" At ",  formatz1(-log(per)*  1/lambdaT), " months the survival probability becomes ",formatz00(100*(per)^(beta1)) ,"%"),
   #      col = "red",
   #      cex = 1)
   # #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   # text(x = max(s$time)*.695, y = .84,                # Text with different color & size
   #      paste0(" At ",  formatz1(-log(per2)*  1/lambdaT), " months the survival probability becomes ",formatz00(100*(per2)^(beta1)) ,"%"),
   #      col = "red",
   #      cex = 1)
   # #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   # 
   # # with the effect of trt and plot
   # y <- x$y^(beta1)   # add another line based on S1(t) = S(0) ^exp(B)
   # lines(y~x$x, col='red')     
   # 
   # # add arrows to explain  
   # arrows(                                  # x start
   #   -log(per2) *1/lambdaT,                   # time  
   #   per2 ,#                                # surv prob  
   #   -log(per2) *1/lambdaT ,  
   #   per2^beta1 ,
   #   col="black", lty=1 )       
   # 
   # arrows(                                  # x start
   #   -log(per) *1/lambdaT,    
   #   per, # 
   #   -log(per) *1/lambdaT ,  
   #   per^beta1 ,
   #   col="black", lty=1 )    
   # 
   # 
   # 
   # #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   # # left plot
   # #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   # 
   # 
   # 
   # plot(s, ylab="Survival probability", xlab="Time", col=c('blue','red'))  # plot the survfit object
   # 
   # text(x = max(s$time)*.58, y = .99,                # Text with different color & size
   #      paste0("HR from Weibull=",formatz2(1/hr[1])," 95%CI ( ",formatz2(1/hr[3]),", ",formatz2(1/hr[2])," )"),
   #      col = "#1b98e0",
   #      cex = 1)
   # text(x = max(s$time)*.58, y = .97,                # Text with different color & size
   #      paste0("HR from Cox PH=",formatz2(hrc[1])," 95%CI ( ",formatz2(hrc[2]),", ",formatz2(hrc[3])," )"),
   #      col = "#1b98e0",
   #      cex = 1)
   # text(x = max(s$time)*.54, y = .95,                # Text with different color & size
   #      paste0( "Actual N= ",sum(N), ", Actual events= ",sum(ev)),
   #      col = "#1b98e0",
   #      cex = 1)
   # 
   # 
   # ## power 
   # A <- 0.05 # alpha
   # B <- 0.1  # beta
   # 
   # s1=per
   # s2=per^(beta1)
   # 
   # # (-log(per)* lambdaT)  # time of s1 s2
   # 
   # f <- ((qnorm(1-A/2))+qnorm(1-B))^2
   # d1 <- (4*f)/ (log(beta1)^2)
   # 
   # N2 <- d1/(1-(s1+s2)/2)
   # 
   # 
   # text(x = max(s$time)*.667, y = .93,                # Text with different color & size
   #      paste0( "Alpha 2 sided ",A," Power ",1-B," events required= ",ceiling(d1), ", N= ",ceiling(N2)),
   #      col = "#1b98e0",
   #      cex = 1)
   # 
   # 
   # 
   # 
   # #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
   # 
   # # estimating hazard plot
   # 
   # H=result.simple
   # plot(H)
   # 
   # #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
   # 
   # 
   # 
   # H= result.pe5
   # H1= result.pe1
   # HS= result.smooth
   # 
   # 
   # plot(H,  col="green")
   # lines(H1, col='blue')
   # lines(HS, col='red')
   # 
   # 
   # 
   # 
   # H2= result.km
   # 
   # haz <- result.smooth$haz.est
   # times <- result.smooth$est.grid
   # surv <- exp(-cumsum(haz[1:(length(haz)-1)]*diff(times)))
   # 
   # plot(H2, conf.int=F, mark="|", xlab="Time",  # confint true throws an error
   #      #xlim=c(0,30), 
   #      ylab="Survival probability")
   # lines(surv ~ times[1:(length(times) - 1)], col='blue')
   # 
   # 
   # #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   # alpha<- as.numeric(shape)
   # 
   # lambda <- as.numeric(scale)
   # 
   # weibSurv <- function(t, shape, scale) pweibull(t, shape=shape,
   #                                                scale=scale, lower.tail=F)
   # curve(weibSurv(x, shape=alpha, scale=1/lambda), from=0, to=2/lambda,
   #       ylim=c(0,1), ylab='Survival probability', xlab='Time')
   # 
   # text(x = 2/lambda*.85, y = .95,                # Text with different color & size
   #      paste0(" Shape ", alpha,  ""),
   #      col = "blue",
   #      cex = 1.2)
   # 
   # text(x = 2/lambda*.85, y = .9,                # Text with different color & size
   #      paste0("  Scale ", lambda ,""),
   #      col = "blue",
   #      cex = 1.2)
   # 
   # 
   # text(x = 2/lambda*.85, y = .85,                # Text with different color & size
   #      paste0(" Mean ", formatz2(gamma(1 + 1/alpha)/lambda),  ""),
   #      col = "blue",
   #      cex = 1.2)
   # 
   # 
   # 
   # text(x = 2/lambda*.85, y = .8,                # Text with different color & size
   #      paste0( " Median ",  formatz2((log(2)^(1/alpha))/lambda ),""),
   #      col = "blue",
   #      cex = 1.2)
   # 
   # 
   # 
   # #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~hazard
   # 
   # weibHaz <- function(x, shape, scale) 
   #   dweibull(x, shape=shape, scale=scale)/
   #   pweibull(x, shape=shape, scale=scale,lower.tail=F)
   # 
   # curve(weibHaz(x, shape=alpha, scale=1/lambda), from=0, to=2/lambda,
   #       ylab='Hazard', xlab='Time', col="red")
   # 
   # 
   # 
   # 
   # #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
   # 
   # 
   # #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   # # text below right plot in changing hazard 
   # 
   # 
   # yo <- abs(100*((beta1/1)-1))
   # 
   # wordd <- ifelse(beta1 < 1,"will reduce", 
   #                 ifelse(beta1 > 1, "will increase ","will not change")) 
   # 
   # c(paste0(" 
   #  With a shape parameter of 1 and a scale parameter of ",lambdaT," a survival curve is created based on a Weibull distribution.
   #  Here we are measuring time in months.
   #  The time to reach median survival is equal to -log(0.5) x ",lambdaT,". 
   #           This equates to a median survival of ", formatz2(-log(.5)*lambdaT)," months. 
   #           
   #           Replacing 0.5 with a desired survival percentile will return the associated time. 
   #           Enter survival probabilities in the two boxes on the left. 
   #           
   #           The time at which the survival probability is ",per* 100,"% is ",  formatz2(-log(per)* lambdaT), " months.
   #           The time at which the survival probability is ",per2*100,"% is ",  formatz2(-log(per2)*lambdaT), " months"
   # ))
   # 
   # 
   # 
   # #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
   # # another piece of text below right plot in changing hazard 
   # 
   # # 
   # # n=          sample$n
   # # 
   # # lambdaT =  sample$base
   # # lambdaC =  sample$censx
   # # beta1 =   sample$hr2
   # # per = sample$per
   # # per2 = sample$per2
   # 
   # yo <- abs(100*((beta1/1)-1))
   # 
   # wordd <- ifelse(beta1 < 1,"will reduce", 
   #                 ifelse(beta1 > 1, "will increase ","will not change")) 
   # 
   # 
   # c(paste0("Now if we are postulating that a new treatment " ,wordd," the hazard by ",yo,"% we can use 
   #             the fact S1(t) = S0(t)^exp Bx, with exp Bx=", (beta1)," the 
   #            probability of survival at ",  formatz2(-log(per)*  lambdaT), " months now becomes ",formatz2((per)^(beta1)) ,".
   #        The probability of survival at ",  formatz2(-log(per2)* lambdaT), " months now becomes ",formatz2((per2)^(beta1)) ,"
   #        . See the arrows in the plot showing the changes in survival. 
   #           " ) )
   # 
   # #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
   # 
   # #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   # # GENERATE THE DATA Execute analysis
   # #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   # 
   # 
   # res <- coxdata(n, allocation, hr, baseline)
   # 
   # d=res$d
   # f=res$f
   # f1=res$f1
   # sf=res$sf
   # np=res$np 
   # LL1=res$LL1
   # LL0=res$LL0
   # S=res$S
   # f0a=res$f0a
   # f0=res$f0
   # f2=res$f2     
   # 
   # 
   # 
   # #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   # 
   # #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   # 
   # 
   # #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   # # MAIN PLOT! updated with log transformation  option
   # #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   # 
   # 
   # f <- f1  # Get the  obj
   # 
   # sf  <- sf
   # X <- as.numeric(as.character(sf[2,c("Effect")]))
   # Y <- as.numeric(as.character(sf[2,c("Lower 0.95")]))
   # Z <- as.numeric(as.character(sf[2,c("Upper 0.95")]))
   # 
   # 
   # p1 <- ggsurvplot(f, main = "Kaplan-Meier Curve",
   #                  
   #                  #  text = paste("Time: ", time),
   #                  
   #                  #palette = c("orange", "purple")
   #                  legend.title = "Trt."
   #                  # , xlab= paste0("Time" )
   #                  ,xlab=paste0("Time : HR=",formatz4(exp(f0)) )
   #                  # , xlab= paste0("Time: HR = ", formatz2(X),", 95%CI( ",formatz2(Y),", ",formatz2(Z)," )" )
   #                  #,#conf.int = TRUE,
   #                  # ggtheme = theme_bw() # Change ggplot2
   #                  #   ,text = paste0("wt: ", round(wt), "</br></br>qsec: ", round(qsec))
   # )
   # ggplotly(p1[[1]] )
   # 
   # 
   # 
   # f <- f1  # Get the survfit obj
   # 
   # p1 <-  ggsurvplot(f, fun = "event",   main = "Cumulative proportion",  
   #                   legend.title = "Trt.")#
   # #  palette = c("orange", "purple"))
   # ggplotly(p1[[1]])
   # 
   # 
   # #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   # 
   # 
   # f <- f1  # Get the survfit obj
   # 
   # p1 <- ggsurvplot(f, fun = "cumhaz",  main = "Cumulative Hazard"  ,
   #                  legend.title = "Trt.")  
   # 
   # p1$plot <- p1$plot + ggplot2::geom_abline(intercept = 0, slope = baseline,    linetype="dotted", col='red') 
   # p1$plot <- p1$plot + ggplot2::geom_abline(intercept = 0, slope = baseline*hr, linetype="dotted", col='blue')
   # 
   # 
   # ggplotly(p1[[1]])
   # 
   # #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   # 
   # f <- f1  # Get the survfit obj
   # 
   # p1 <- ggsurvplot(f, fun = "cloglog",  
   #                  main = "Complementary log−log" ,
   #                  legend.title = "Trt.")
   # # legend.labs = c("0", "1")) 
   # #palette = c("jco"))
   # ggplotly(p1[[1]])
   # 
   # 
   # 
   # #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   # 
   # f <- np  # Get the  data
   # survdiffplot(f, col='darkgreen' , xlab= "Time")
   # 
   # 
   # #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   # 
   # d <- d
   # 
   # dd <<- datadist(d)
   # options(datadist='dd')
   # 
   # fit <- cph(Surv(dt, e) ~ trt, data = d, x=TRUE, y=TRUE, surv=TRUE)
   # 
   # survplot(fit, type="kaplan-meier", conf.int = TRUE, 
   #          col.fill=c("firebrick1","cyan2"), grid=TRUE, what='survival')
   # 
   # 
   # #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   # 
   # 
   # S <- Surv(d$dt, d$e)
   # 
   # hazard.ratio.plot(d$trt, S, e=20, legendloc='ll', xlab='Time', antilog=FALSE, col='blue', smooth=TRUE,
   #                   ylim=c(-4,4), ylab=c("Log hazard ratio"))
   # 
   # 
   # 
   # #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   # 
   # 
   # fit <- cph(Surv(dt, e) ~ trt, data = d, x=TRUE, y=TRUE, surv=TRUE)
   # 
   # plot(cox.zph(fit, transform="identity" ), ylim=c(-4,4), ylab=c("Log hazard ratio"))
   # 
   # 
   # #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   # 
   # 
   # #plot(survfit(S~ d$trt), col=c("purple", "orange"), fun="cloglog", xlab="Time", ylab="log(-log(Survival)" , lwd=3)
   # survplot(f,logt=TRUE, loglog=TRUE, 
   #          col=c("red", "lightblue")
   #          
   #          # col=c("orange", "purple")
   # )   #Check for Weibull-ness (linearity)
   # 
   # #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   # #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
   # p1 <- ggsurvplot(f2, main = "Kaplan-Meier Curve", legend.title = "Trt.",
   #                  # palette = c("orange", "purple")  ,
   #                  xlab=paste0("Time : HR=", formatz4(exp(f0a)))
   #                  # ggtheme = theme_bw() # Change ggplot2 theme
   # )
   # ggplotly(p1[[1]])
   # 
   # 
   # #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   # 
   # 
   # d <- plyr::arrange(d, dt)
   # 
   # DT::datatable(d, rownames=FALSE,
   #               plugins = 'natural',
   #               colnames=c('Time' = 'dt', 'Event or censored' = 'e', 
   #                          'Treatment'='trt'),
   #               
   #               options = list(
   #                 #  dom = 't',
   #                 columnDefs = list(list(type = 'natural', targets = c(1,2)))
   #               )
   # ) %>%
   #   
   #   formatRound(
   #     columns= c("Time" ), 
   #     digits=4 )
   # 
   # 
   # #~~~~~~~~~~~~~~~~~~~~~~KM table
   # 
   # 
   # library(survminer)
   # 
   # require(ggfortify)
   # 
   # KM_fit <- survfit(Surv(dt, e) ~ trt ,data = d)
   # 
   # KM <- fortify(KM_fit) # fortify(KM_fit, fun = "cumhaz")'
   # 
   # DT::datatable(KM, rownames=FALSE,
   #               plugins = 'natural',
   #               #   colnames=c('Time' = 'dt', 'Event or censored' = 'e', 
   #               #              'Treatment'='trt'),
   #               
   #               options = list(
   #                 #  dom = 't',
   #                 columnDefs = list(list(type = 'natural', targets = c(1,2)))
   #               )
   # ) %>%
   #   
   #   formatRound(
   #     columns= c("time" ,"surv","std.err","upper","lower"), 
   #     digits=4 )
   # 
   # #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   # 
   # require(ggfortify)
   # 
   # KM_fit <- survfit(Surv(dt, e) ~ trt ,data = d)
   # 
   # KM <-  fortify(KM_fit, fun = "cumhaz")
   # 
   # DT::datatable(KM, rownames=FALSE,
   #               plugins = 'natural',
   #               #   colnames=c('Time' = 'dt', 'Event or censored' = 'e', 
   #               #              'Treatment'='trt'),
   #               
   #               options = list(
   #                 #  dom = 't',
   #                 columnDefs = list(list(type = 'natural', targets = c(1,2)))
   #               )
   # ) %>%
   #   
   #   formatRound(
   #     columns= c("time" ,"surv","std.err","upper","lower"), 
   #     digits=4 )
   # 
   # 
   # #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~maximum likelihood
   # 
   # # Get the  data
   # 
   # y <- as.numeric(as.character(res$f$coefficients))
   # guess=y  # we use the actual model hr estimate
   # 
   # d <- plyr::arrange(d, dt)
   # 
   # 
   # d$expB <- exp(guess*d$trt)
   # d$part2 <- guess*d$trt
   # d$part3 <- log(rev(cumsum(rev(d$expB))))
   # 
   # d$likelihoodi <- d$e*(d$part2 - d$part3)
   # d$LL <- sum(d$likelihoodi)
   # 
   # datatable(d, rownames=FALSE,
   #           plugins = 'natural',
   #           colnames=c('Time' = 'dt', 
   #                      'Event or censored' = 'e', 
   #                      'Treatment'='trt',
   #                      'HR'='expB',
   #                      'Individual likelihoods' ='likelihoodi',
   #                      'logHR x trt'='part2',
   #                      'Log(exp HR) of each risk set'='part3',
   #                      'Sum the Individual likelihoods to give log likelihood' ='LL'
   #           ),
   #           
   #           options = list(
   #             #  dom = 't',
   #             columnDefs = list(list(type = 'natural', targets = c(1,2)))
   #           )
   # ) %>%
   #   
   #   formatRound(
   #     columns= c("Time","HR", 
   #                #"A",
   #                "logHR x trt","Log(exp HR) of each risk set",'Individual likelihoods'), 
   #     digits=4 ) %>%
   #   formatRound(
   #     columns= c("Sum the Individual likelihoods to give log likelihood"), 
   #     digits=1 )
   # 
   # #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   # 
   # 
   # dummy = res
   # 
   # g <- log(as.numeric(res$f$coefficients))
   # 
   # foo <- loglike2(g, dat=dummy$d, dead=dummy$d$e, indep=dummy$d$trt, time=dummy$d$dt)  # try
   # foo$time =NULL
   # foo$expB = NULL
   # 
   # datatable(foo, rownames=FALSE,
   #           plugins = 'natural',
   #           colnames=c('Time' = 'dt', 
   #                      'Event or censored' = 'e', 
   #                      'Treat.'='trt',
   #                      'Model Hazard Ratio'='hr',
   #                      'Null Log Likelihood'='lognull',
   #                      'Maximised Log Likelihood'='lognmax',
   #                      'HR guess' ='guess',
   #                      'Individual likelihoods' ='likelihoodi',
   #                      'Sum the Individual likelihoods to give log likelihood based on guess' ='L'
   #           ),
   #           
   #           options = list(
   #             #  dom = 't',
   #             columnDefs = list(list(type = 'natural', targets = c(1,2)))
   #           )
   # ) %>%
   #   
   #   formatRound(
   #     columns= c("Time","Model Hazard Ratio",  
   #                'Individual likelihoods'
   #     ),
   #     digits=4 ) %>%
   #   formatRound(
   #     columns= c( 'Null Log Likelihood',
   #                 'Maximised Log Likelihood',
   #                 'Sum the Individual likelihoods to give log likelihood based on guess'), 
   #     digits=0 )
   # 
   # 
   # #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   # 
   # 
   # # allocation=sample$allocation
   # # hr=        sample$hr
   # # baseline=  sample$baseline
   # 
   # res <- coxdata(n=10, allocation=allocation, hr=hr, baseline =baseline)
   # 
   # d <- res$d
   # d <- plyr::arrange(d, dt)  # sort by time
   # 
   # # Calculate Li for everyone
   # d$Numerator   <- exp(res$f$coefficients[[1]] * d$trt)
   # d$Denominator <- (rev(cumsum(rev(d$Numerator))))
   # d$Li          <- d$Numerator/d$Denominator
   # 
   # # all censored contribute 1 (on multiplicative scale)
   # d$Li<- ifelse(d$e %in% 1,d$Li,1)
   # 
   # # get the product of all and log answer
   # d$LL <- log(prod(d$Li))  
   # d
   # # model LL, prove we have ecalc correctly
   # res$f$loglik
   # 
   # 
   # datatable(d, rownames=FALSE,
   #           plugins = 'natural',
   #           colnames=c('Time' = 'dt', 
   #                      'Event or censored' = 'e', 
   #                      'Treat.'='trt',
   #                      'Num.'='Numerator',
   #                      'Den.'='Denominator',
   #                      'Individual likelihoods'='Li',
   #                      'log of product of Individual likelihoods to give maximizes log likelihood' ='LL'
   #                      
   #           ),
   #           
   #           options = list(
   #             #  dom = 't',
   #             columnDefs = list(list(type = 'natural', targets = c(1,2)))
   #           )
   # ) %>%
   #   
   #   formatRound(
   #     columns= c("Time","Num.","Den.", 
   #                'Individual likelihoods', 'log of product of Individual likelihoods to give maximizes log likelihood'
   #     ),
   #     digits=4 ) %>%
   #   formatRound(
   #     columns= c( 'Event or censored',
   #                 'Treat.'), 
   #     digits=0 )
   # 
   # 
   # #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   # 
   # X <- as.numeric(as.character(sf[2,c("Effect")]))
   # Y <- as.numeric(as.character(sf[2,c("Lower 0.95")]))
   # Z <- as.numeric(as.character(sf[2,c("Upper 0.95")]))
   # 
   # Xp  <- X/(X+1)
   # Yp  <- Y/(Y+1)
   # Zp  <- Z/(Z+1)
   # 
   # wordup <- ifelse(X>1,"higher", "")
   # 
   # wordup2 <- ifelse(X>1,"increase", "reduction")
   # 
   # paste0( "From the Cox model the estimated hazard ratio is "
   #         , formatz2(X),", 95%CI ( ",formatz2(Y),", ",formatz2(Z),
   #         " ) comparing treatment 1 to 0. 
   #          
   #           A hazard ratio of  ", formatz2(X)," means that, in each unit of time, someone 
   #          treated in group 1 has ", formatz00(abs(X/1-1)*100),"% ", wordup ," of the chance of experiencing the event of interest
   #          in the following unit of time as they would were they taking treatment 0.
   #          
   #          There is an estimated ", formatz00(abs(X/1-1)*100),"% ", wordup2 ," in the hazard of the outcome by the Cox model. 
   #          
   #          Equivalently, the hazard ratio is equal to the odds that a patient in treatment group 1 experiences the event of interest before a
   #          a patient in treatment group 0.
   #          
   #         Therefore we can reformulate the hazard ratio, possibly more intuitively, as
   #          the probability that a patient in treatment 
   #          group 1 experiences the event before a patient in treatment group 0, which is: "
   #         , formatz2(Xp),", 95%CI ( ",formatz2(Yp),", ",formatz2(Zp),").        ")
   # 
   # 
   # 
   # #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   # # frank Harrell rms page 479
   # 
   # 
   # d <-  d
   # 
   # trt <-  d$trt
   # e   <-  d$e
   # dt  <-  d$dt
   # d   <-  d$d
   # 
   # limx <- quantile(dt, prob=.99)
   # limx <- max(dt)*.8
   # #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   # S <- Surv(dt,e)
   # f <- npsurv(S ~ trt)
   # 
   # for (meth in c('exact','breslow','efron')) {
   #   
   #   g <- cph(S  ~ trt, method=meth, surv=TRUE, x=TRUE, y=TRUE)
   #   
   # }
   # #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   # f.exp  <- psm(S  ~ trt, dist ='exponential')
   # fw    <-  psm(S  ~ trt,  dist ='weibull')
   # phform <- pphsm(fw)
   # #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   # 
   # co <- gray(c(0,.8))
   # co <- c("red", "lightblue")
   # survplot(f, lty=c(1,1)   , lwd=c(1,3), col=co,           label.curves=FALSE, conf='none')
   # survplot(g, lty=c(3,3)   , lwd=c(1,3), col=co, add=TRUE, label.curves=FALSE, conf.type='none')
   # 
   # legend(c(limx,160),c(.38,.99),
   #        c('Nonparametric estimates', 'Cox-Breslow estimates'),
   #        lty=c(1,3), bty='n',    cex=1.0) # col=co 
   # 
   # legend(c(limx,160),c(.18,.89), 
   #        c('Trt 0','Trt 1'), lwd=c(1,3), col=co, bty='n',  cex=1.0)
   # #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   # 
   # 
   # #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   # 
   # #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   # 
   # 
   # 
