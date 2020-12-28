

# https://www.stata.com/meeting/uk14/abstracts/materials/uk14_crowther.pdf
# slide 9

n=1000
u <- runif(n)
trt <- sample(0:1, n, replace=TRUE)
lambda=.1
gamma=1.2
loghr = .7
s = (-log(u)/lambda*exp(loghr*trt))^(1/gamma)
dead=rep(1,n)

require(survival)
f <- coxph(Surv(s, dead) ~ 1 )
summary(f)
f$loglik                                       # null and maximised


f1 <- Surv(s, dead)

plot(survfit(f1 ~ 1), col=c("purple", "red"), xlab="Time" )
plot(survfit(f1 ~ 1), col=c("purple", "red"), fun="cloglog", xlab="Time", ylab="log(-log(Survival)")