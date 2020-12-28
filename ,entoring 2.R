
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


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



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


set.seed(123)
allocation=.5
hr=        2
baseline=  .4

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
