#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@'
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@'
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@'
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Demonstrate  it is only ranks that count, mention ties!
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   rm(list=ls())
   require(surival)
   require(rms)
   
   # FUNCTION SIMULATE AND TO DO PLL CALC
   coxdata <- function(n, allocation, hr, baseline) { 
      
      #n=1000; allocation =.5; hr=2; baseline=.4
      
      trt <- sample(0:1, n,  rep=TRUE, prob=c(1-allocation, allocation))
      
      cens <- 15*runif(n)
      
      h <- baseline*exp(log(hr)*(trt==1))  # hazard function h(t)
      
      dt <- -log(runif(n))/h   # <- exemplify this
      
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
      f0 <- f$coefficients[[1]]  
      LL1 <- f$loglik[2]
      LL0 <- f$loglik[1]
      
      sf <- summary(f)
      
      f1 <- survfit(Surv(dt,e) ~ trt, data = d)
      
      np <- npsurv(Surv(dt,e) ~ trt,d)
      
      S <- Surv(d$dt, d$e)
      
      do<-d
      do <- plyr::arrange(do,dt) # collect this dataset

      d <- plyr::arrange(d,dt)
      d$dt <- dd<- NULL
      d$dt <- sort(2*rexp(nrow(d))) # new times here!
      
      dx <<- datadist(d)
      options(datadist='dx')
      f0a <- cph(Surv(dt,e) ~  trt, x=TRUE, y=TRUE, data=d )
      f0a <- f0a$coefficients[[1]]
      f2 <- survfit(Surv(dt ,e)  ~ trt, data = d)
      
      #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
      
      return(list(f=f, d=d, do=do, f1=f1, sf=sf, np=np, LL1=LL1, LL0=LL0, S=S,                  
                  
                  f0=f0,  f2=f2, f0a=f0a, foo=foo))
      
   }
   
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@'
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# quickly show partial log likelihood calc (without logging til end)
# in reality you don't know HR, it is found by iteration, maximizing LL
# create small dataset, 
# we are doing it this way to exemplfy
# large datasets will cause numerical problems
# so we should calculate on log space really

## show the cox phworkings

   set.seed(8713)
   res <- coxdata(n=10, allocation=.5, hr=2, baseline = .4)
   f <- cph(Surv(dt,e) ~ trt, x=TRUE, y=TRUE, data=res$d )
   f
   d <- res$d
   d <- plyr::arrange(d, dt)  # sort by time
   
   # Calculate Li for everyone
   d$numerator   <- exp(res$f$coefficients[[1]] * d$trt)
   
  #d$numerator   <- exp(0 * d$trt)   # show this
   d$denominator <- (rev(cumsum(rev(d$numerator))))
   d$Li          <- d$numerator/d$denominator
   
   # all censored contribute 1 (on multiplicative scale)
   d$Li2<- ifelse(d$e %in% 1,d$Li,1)
   
   # get the product of all and log answer
   d$LL <- log(prod(d$Li2))  
   print(d,digits=6)
   
   # model LL, prove we have ecalc correctly
   res$f$loglik
   f <- cph(Surv(dt,e) ~  trt, x=TRUE, y=TRUE, data=res$d )

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@'
   # time itself does not matter
   res <- coxdata(n=100, allocation=.5, hr=2, baseline = .4)
   
   #~~~~~~~~~~~~~~~~~~~~~~~~~
   fx <-  res$f2 # Get the  obj #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
   f0a <- res$f0a               #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
   p1 <- ggsurvplot(fx, main = "Kaplan-Meier Curve", legend.title = "Trt.",
                    # palette = c("orange", "purple")  ,
                    xlab=paste0("Time : HR=",  (exp(f0a)))
                    # ggtheme = theme_bw() # Change ggplot2 theme
   )
   A<-ggplotly(p1[[1]])
   #~~~~~~~~~~~~~~~~~~~~~~~~~~
   
   fx <-  res$f1 # Get the  obj #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
   f0a <- res$f0                #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
   p2 <- ggsurvplot(fx, main = "Kaplan-Meier Curve", legend.title = "Trt.",
                    # palette = c("orange", "purple")  ,
                    xlab=paste0("Time : HR=",  (exp(f0a)))
                    # ggtheme = theme_bw() # Change ggplot2 theme
   )
   B <- ggplotly(p2[[1]])
   #~~~~~~~~~~~~~~~~~~~~~~~~~
   
   require(gridExtra)
   subplot(A,B,nrows=1, shareX=TRUE , titleX=TRUE)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# tail(res$d)
# tail(res$do)

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@'
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@'
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@'
#~~~~~~~~~~~~~~~~~~~~~~PDF CDF S(t))
rm(list=ls())
require(ggplot2)
set.seed(100)
Time <- data.frame(Time=rnorm(1000, 500, 150))
rand1 <- 400

# pdf
p <- ggplot(Time, aes(Time)) +
   geom_density(fill="grey") #+
# subset region and plot
d <- ggplot_build(p)$data[[1]]

p <- p + geom_area(data = subset(d, x > rand1), aes(x=x, y=y), fill="lightgrey") +
   geom_segment(x=rand1, xend=rand1, 
                y=0, yend=approx(x = d$x, y = d$y, xout = rand1)$y,
                colour="black", size=1)  + ylab("Density (Relative frequency)") +
   
   annotate("text", x=250,  y=0.0002, label= "F(t)", size=4) + 
   annotate("text", x =650, y=0.0002, label = "S(t) = 1 - F(t)", size=4) +
   annotate("text", x =660, y=0.002,  label = "f(t)", size=4) +
   ggtitle("f(t) (death density) as a function of time. The cumulative \nproportion of the population that has died up to time t equals F(t). \nThe prop of the pop that has survived to time t is S(t) = 1- F(t)") +
   theme_bw()

# cdf
P2 <-  ggplot(Time, aes(Time)) + stat_ecdf(geom = "step") +
   ggtitle("F(t). Cumulative proportion of deaths with time")  +
   theme_bw()

# S(t)
pg <- ggplot_build(P2)$data[[1]]
P3 <- ggplot(pg, aes(x = x, y = 1-y )) + geom_step() +
   ggtitle("S(t) The proportion of the population survived to time t \n S(t) = 1 - F(t) = P(T>t)")  + 
   xlab("Time") +
   ylab("Survival probability") +
   theme_bw()


require(gridExtra)
plot1 <- p
plot2 <- P2
plot3 <- P3
grid.arrange(plot1, plot2, plot3, ncol=3)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@'
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@'
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@'
# ANother PDF CDF S(t) not used
   rm(list=ls())
   require(ggplot2)
   set.seed(1088)
   
   n=20
   rand1 <- 400
   Time=rnorm(n, 500, 150)
  # Time=rexp(n, rate=.003)
  
   Time1 <- data.frame(Time=Time)
   
   #mean1 <- mean(amount_spent1$amount_spent)
   rand1 <- 400
   
   
   p <- ggplot(Time1, aes(Time)) +
      geom_density(fill="grey")# +
    #  geom_vline(xintercept=mean1) 
   
   # subset region and plot
   d <- ggplot_build(p)$data[[1]]
   
   # p <- p + geom_area(data = subset(d, x > rand1), aes(x=x, y=y), fill="red") +
   #    geom_segment(x=rand1, xend=rand1, 
   #                 y=0, yend=approx(x = d$x, y = d$y, xout = rand1)$y,
   #                 colour="blue", size=3)
 
   
  
   dput(   sort(as.vector(unlist(Time1))))

 #  pdf
   # p <- ggplot(amount_spent1, aes(Time1)) +
   #    geom_density(fill="grey") #+
   # # subset region and plot
   # d <- ggplot_build(p)$data[[1]]
   #
    p <- p + geom_area(data = subset(d, x > rand1), aes(x=x, y=y), fill="lightgrey") +
      geom_segment(x=rand1, xend=rand1,
                    y=0, yend=approx(x = d$x, y = d$y, xout = rand1)$y,
                    colour="black", size=1)  + ylab("Density (Relative frequency)") +
   
      annotate("text", x=350,  y=0.0002, label= "F(t)", size=4) +
       annotate("text", x =650, y=0.0002, label = "S(t) = 1 - F(t)", size=4) +
       annotate("text", x =600, y=0.002,  label = "f(t)", size=4) +
       ggtitle("f(t) (death density) as a function of time. The cumulative \nproportion of the population that has died up to time t equals F(t). \nThe prop of the pop that has survived to time t is S(t) = 1- F(t)") +
       theme_bw()
   
  
   

# cdf
   P2 <-  ggplot(Time1, aes(Time)) + stat_ecdf(geom = "step") +
      ggtitle("F(t). Cumulative proportion of deaths with time")  +
      theme_bw()

# S(t)
   pg <- ggplot_build(P2)$data[[1]]
   P3 <- ggplot(pg, aes(x = x, y = 1-y )) + geom_step() +
      ggtitle("S(t) The proportion of the population survived to time t \n S(t) = 1 - F(t) = P(T>t)")  + 
      xlab("Time") +
      ylab("Survival probability") +
      theme_bw()


   require(gridExtra)
   plot1 <- p
   plot2 <- P2
   plot3 <- P3
   grid.arrange(plot1, plot2, plot3, ncol=3)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@'
   #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@'
   #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@'
    
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
 

   plot(haz,   t='l', log="y",
        xlab="age", ylab="Annual ‘hazard’ or ‘force of mortality h(t)", 
        main="Annual risk of death in males from all causes for UK 2017-19")  
   #lines(dd$x,  exp(predict(f)), col = 'red',  lwd=2, lty=2)
   
   
   cumhaz <- cumsum(haz)
   Surv <-   exp(-cumhaz)
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

 
grid.arrange(pa,  pc, pb, ncol=3,

             top = textGrob(paste0( "Weibull"),
                            gp=gpar(fontsize=20,font=3)) )

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~









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
# PLOT iNVESRS SAMPLinG MEt#hos AND TRUE DISTRIBUTIONS


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

   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   #Inverse cumulative hazard function
   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   # In addition to f(), F(), S(), h(), and H() -all different ways of summarizing the same
   # information-a sixth way, Q(), is of use for to create artificial survival-time data sets.
   # if U ∼ U[0; 1], then (1 − U) ∼ U[0; 1] so dispense with 1- in weibull formula
   # uniformly distributed random numbers can be transformed into survival times

   n=100
   # functions to plot true survival function curves red dash
      weibSurv <- function(t, shape, scale)   pweibull(t, shape=shape, scale=scale, lower.tail=F)
      gSurv     <- function(t, shape, scale) pgompertz(t, shape=p,     rate = h,    lower.tail=F)
   
   par(mfrow=c(1,3))  
   
      # simulate weibull shape parameter p, scale = h
      p=1/3    # shape, if 1 this is exponential
      h=0.05   # scale
      
      #~~~~~~~~~~~~~~~~~~~~~~~~~~
      t <- h^-1*(-log(runif(n)))^p # weibull random times 
      
      survfit <- survfit(Surv(t) ~ 1)
      plot(survfit, ylab="Survival probability", xlab="Time", 
           main=paste0("N = ",n,", Weibull, shape = ",1/p,", rate = ",h))
      
      # plot true survival curve for weibull
      curve(weibSurv(x, shape=1/p, scale=1/h), from=0, to=max(t), n=length(t), 
            col='red', lwd=2, lty=2,
            ylim=c(0,1), add=TRUE)
      
      #~~~~~~~~~~~~~~~~~~~~~~~~~~
      t <- -log(runif(n))/h # exponential rate = h random times
      
      survfit <- survfit(Surv(t) ~ 1)
      plot(survfit, ylab="Survival probability", xlab="Time", 
           main=paste0("N = ",n,", Exponential, rate = ",h))
    
      # plot true survival curve, constant hazard h
      curve(weibSurv(x, shape=1, scale=1/h), from=0, to=max(t), n=length(t), 
            col='red', lwd=2, lty=2,
            ylim=c(0,1), add=TRUE)

      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      require(flexsurv) # for gompertz true distribution
   
      t <- (1/p)*log((-log(runif(n)))*(p/h)+1) # gompertz, see bender table 1 where -log(runif(n)) is t
      
      survfit <- survfit(Surv(t) ~ 1)
      plot(survfit, ylab="Survival probability", xlab="Time", 
           main=paste0("N = ",n,", Gompertz, shape = ",round(p,2),", scale = ",h))
      
      # plot true survival curve, 
      curve(gSurv(x, shape=p, scale=1/h), from=0, to=max(t), n=length(t), 
            col='red', lwd=2, lty=2,
            ylim=c(0,1), add=TRUE)
   
   par(mfrow=c(1,1))
   
   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   
   
   
    