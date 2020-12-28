## harrells survival functions

s1<-.7
s2<-.5
t1<-17
t2<-23
hr<-1.2


N1<-700
N2<-700
sim <-300
drop <- .1
start <- 36
fini <- 60

library(Hmisc)
Weib.p <- Weibull2(c(t1,t2),c(s1,s2))

ff <- Quantile2(Weib.p,hratio=function(x) hr,  mplot=200)
plot(ff )


rcontrol      <- function(n) ff(n, what='control')
rintervention <- function(n) ff(n, what='intervention')
# 
# Here, let us accrue patients over three years, and
# follow them for an additional seven years. 
# 
# For now we shall assume that the accrual
# will follow a uniform distribution. This leads to the censoring distribution being
# uniform, with a minimum of five years (for a patient entering at the end of the
#                                        accrual period) and a maximum of eight years (for a patient entering at the start of
#                                                                                      the trial). 
 

rcens <- function(n) runif(n, start, fini)#


ff.dropout <- Quantile2(Weib.p,hratio=function(x) hr,
                        dropout=function(x) drop)

plot(ff.dropout)

 rcontrol <-      function(n) ff.dropout(n, what='control')
 rintervention <- function(n) ff.dropout(n, what='intervention')

 x<-spower(rcontrol, rintervention, rcens, 
           nc=N1, 
           ni=N2,
           test=logrank, nsim=sim, alpha=0.025, cox=TRUE)

x






#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
##from spower, trying to get one data set

yc <- rcontrol(N1)
yi <- rintervention(N2)
cens <- rcens(N1+ N2)
group <- c(rep(1, N1), rep(2, N2))
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
p1 <- ggsurvplot(f1, main = "Kaplan-Meier Curve",
                 
                 #  text = paste("Time: ", time),
                 
                 #palette = c("orange", "purple")
                 legend.title = "Trt."
                 # , xlab= paste0("Time" )
                 ,xlab=paste0("Time : HR=",round(exp(fit$coefficients),4))
                 # , xlab= paste0("Time: HR = ", formatz2(X),", 95%CI( ",formatz2(Y),", ",formatz2(Z)," )" )
                 #,#conf.int = TRUE,
                 # ggtheme = theme_bw() # Change ggplot2
                 #   ,text = paste0("wt: ", round(wt), "</br></br>qsec: ", round(qsec))
)
p1
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# see if we can get the weibull parameters for intervention !!

# i. with 2 times and 2 survival prob get WEibull parameters using Weibull2 function
# ii. use Quantile2 function and get curve also get a curve for intervention using hr input
# iii. pull out intervention curve data
# iv. pull out control curve data
# v. Get weibull parameters of intervention, by examining intervention curve and getting 2 times and 2 probs
# vi feed them intervention weibull parameters into weibull2 


hr=1.2
# these are the control time and p(survival)
s1<-.7
s2<-.5
t1<-17
t2<-23
 
# hr=.6
# t1<-27
# t2<-43




Weib.p <- Weibull2(c(t1,t2),c(s1,s2))
Weib.p
ff <- Quantile2(Weib.p,hratio=function(x) hr ) # we get weibull parameters 
 

##pull out intervention survival probs
fff <-      attributes(ff)  #lets get the data
time <-     fff$plot.info$`I Survival`$Time
survival <- fff$plot.info$`I Survival`$Survival
#plot(survival~time, lty=1)   # interventions

#pull out control survival data
timec <-     fff$plot.info$`C Survival`$Time
survivalc <- fff$plot.info$`C Survival`$Survival

###~~~~~~~~~~~~~~~~~~~~~~~~# what time is survival ~.5 intevention

s50i <- which.min(abs(survival-s2)) # what index is time ~ .5
s70i <- which.min(abs(survival-s1))

Weib.i <- Weibull2(c(  time[s70i], time[s50i]),c(s1, s2))
Weib.i
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

text(quantile(prob=jump0,c(timec,timei)),  0.96, c(round(t1x ,5)), cex = 1)
text(quantile(prob=jump1,c(timec,timei)), 0.96, A,                cex = 1)

text(quantile(prob=jump0,c(timec,timei)),  0.90, c(round(t2x ,5)), cex = 1)
text(quantile(prob=jump1,c(timec,timei)), 0.90, B,                cex = 1)

text(quantile(prob=jump0,c(timec,timei)),  0.84, c(round(t3x ,5)), cex = 1)
text(quantile(prob=jump1,c(timec,timei)), 0.84, C,                cex = 1)

text(quantile(prob=jump0,c(timec,timei)),  0.78, c(round(t4x ,5)), cex = 1)
text(quantile(prob=jump1,c(timec,timei)), 0.78, D,                 cex = 1)


s1i=survival[which.min(abs(time-t1))]
s2i=survivali[which.min(abs(timei-t2))]

 
text(quantile(prob=jump,c(timec,timei)), 0.66, paste0("Control survival prob ",s1," at time ",t1,""),   cex = 1)
text(quantile(prob=jump,c(timec,timei)), 0.60, paste0("Interv. survival prob ",round(s1i,2)," at time ",t1,""),   cex = 1)

text(quantile(prob=jump,c(timec,timei)), 0.54, paste0("Control survival prob ",s2," at time ",t2,""),   cex = 1)
text(quantile(prob=jump,c(timec,timei)), 0.48, paste0("Interv. survival prob ",round(s2i,2)," at time ",t2,""),   cex = 1)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# add arrows to explain  
arrows(                                   
    t1,                                  # x start  
    s1 ,#                                # surv prob at t1 in control
    t1 ,                                 # x finish
    s1i ,                                # surv prob at t1 in intervention     
  col="black", lty=1 )       

arrows(                                   
  t2,                                    # x start  
  s2 ,#                                  # surv prob at t2 in control
  t2 ,                                   # x finish
  s2i,                                   # surv prob at t2 in intervention     
  col="black", lty=1 )          

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


 
