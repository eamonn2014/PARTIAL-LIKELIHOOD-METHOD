


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
  
}#

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# quickly show partial log likelihood calc (without logging til end)
# in reality you don't know HR, it is found by iteration, maximizing LL

set.seed(1234)

# create small dataset, 
# we are doing it this way to exemplfy
# large datasets will cause numerical problems
# so we should calculate on log space really

res <- coxdata(n=10, allocation=.5, hr=2, baseline = .4)

d <- res$d
d <- plyr::arrange(d, dt)  # sort by time

# Calculate Li for everyone
d$numerator   <- exp(res$f$coefficients[[1]] * d$trt)
d$denominator <- (rev(cumsum(rev(d$numerator))))
d$Li          <- d$numerator/d$denominator

# all censored contribute 1 (on multiplicative scale)
d$Li2<- ifelse(d$e %in% 1,d$Li,1)

# get the product of all and log answer
d$LL <- log(prod(d$Li2))  
d
# model LL, prove we have ecalc correctly
res$f$loglik
