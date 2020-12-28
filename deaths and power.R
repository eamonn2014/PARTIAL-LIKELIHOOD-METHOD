



##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    n=800; lambdaT=.03; lambdaC=0.02; beta1=2; per=0.7 ;per2=.652
   
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~use log beta
    
    # sample size, alpha beta and hr and then s1 s2 needed 
    s1=per
    s2=per2
    
     
    A <- 0.05 # alpha
    B <- 0.2  # beta
    b <- beta1
    
    
    f<- ((qnorm(1-A/2))+qnorm(1-B))^2
    d <- (4*f)/ (log(b)^2)
    d  # deaths
    
    N <- d/(1-(s1+s2)/2)
    N # total sample
    
    ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    require(survival)
    lambdaT=.03; lambdaC=0.09
    n <- ceiling(200)
    x1 = sample(0:1, n,replace=TRUE)#
    
    
    #lambdaC <- lambdaT^1.0
    
    beta=log(2)
    
    # true event time
    T = rweibull(n, shape=1, scale=1/lambdaT)*exp(-beta*x1)
    C = rweibull(n, shape=1, scale=1/lambdaT)                  # censoring time
   #  C=1e6
    time = pmin(T,C)                                         # observed time is min of censored and true
    event = time==T                                          # set to 1 if event is observed
    table(event)[1]/table(event)[2]
    table(event) 

    
    
    
    
    
    Surv(time, event)
    
    
    
    
    
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # now we simulate using hazard and N and hr

 # N=2915; lambdaT=14.4; lambdaC=12; beta1=1.2; per=0.7 ;per2=.652

  spow <- function(n=800, lambdaT=0.03, lambdaC=1, beta1=2) {
  
    x1 = sample(0:1, n,replace=TRUE)#
    
    beta=log(beta1)
    
    # true event time
    T = rweibull(n, shape=1, scale=1/lambdaT)*exp(-beta*x1)
    C = rweibull(n, shape=1, scale=1/lambdaC)                  # censoring time
    # C=1e6
    time = pmin(T,C)                                         # observed time is min of censored and true
    event = time==T                                          # set to 1 if event is observed
    table(event)
    f1 <-  coxph(Surv(time, event)~ x1 , method="breslow")
     
    df <- length(coef(f1)) 
   
    teststat<--2*(as.numeric(f1$loglik)[1]-as.numeric(f1$loglik[2])) #  null-model
   
    p <- pchisq(teststat,df=df,lower.tail=FALSE)
  
    return(p)
  }
  
   mean(replicate(1000, spow(n=ceiling(d)))<0.05)
  
  
  
 
 








lrtest(f,f1)
}


survfit <- survfit(Surv(time,event) ~ x1)
f$stats
f
confint(f)



#yr
s1=.7
s2=.652

# 10yr
s1=.5
s2=.435
  

A <- 0.05 # alpha
B <- 0.2  # beta
Theta <- 1.2


f<- ((qnorm(1-A/2))+qnorm(1-B))^2
d <- (4*f)/ (log(Theta)^2)
d  # deaths

N <- d/(1-(s1+s2)/2)
N # total sample


 


lambdaT=14
s1=.7                  # original survival
beta=1.2               # change in hazard
(-log(per)* lambdaT)   # time
s1^beta                # new survival


###

library(survival)
sim = function(n_trt, n_ctrl, m_trt, m_ctrl, t_censoring) {
  # Simulate some data
  t_trt = rexp(n_trt, log(2)/m_trt)
  
  t_ctrl = rexp(n_ctrl, log(2)/m_ctrl)
  t = c(t_trt, t_ctrl)
  
  # Add censoring
  died = (t <= t_censoring)
  t_after_censoring = pmin(t, t_censoring)
  
  # Fit Cox model
  gr = c(rep(1, n_trt), rep(0, n_ctrl))
  l = coxph(Surv(t_after_censoring, died) ~ gr)
  anova(l)$`Pr(>|Chi|)`[2] # P-value from likelihood ratio test
}
set.seed(123) # Reproducible pseudorandom numbers
n_sim = 1000  # Number of simulations (increase this for more accuracy)
p_med = replicate(n_sim, sim(n_trt = 155, n_ctrl = 78,
                             m_trt = 7, m_ctrl = 4.4,
                             t_censoring = 24))
p_hr = replicate(n_sim, sim(n_trt = 84, n_ctrl = 42,
                            m_trt = 4.4/0.54, m_ctrl = 4.4,
                            t_censoring = 24))
mean(p_med <= 0.05) # Estimated power using median
#> [1] 0.901
mean(p_hr <= 0.05)  # Estimated power using control median + HR
#> [1] 0.88