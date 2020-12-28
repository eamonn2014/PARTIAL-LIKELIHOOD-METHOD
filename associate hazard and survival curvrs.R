# moore's book page 17 sample size
# In R the functions “dweibull” and “pweibull” compute the p.d.f. and c.d.f.,
# respectively, of the Weibull distribution. These functions use the arguments “shape”
# and “scale” to represent the parameters ˛ and 1=, respectively. To obtain the
# survival function, we can specify “lower.tail = F” as an option in the “pweibull”
# function. For example, we can plot the Weibull survival function with ˛ D 1 and
#  D 0:03


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
alpha=1
lambda=0.03
upper <- 100
par(mfrow=c(1,2))
weibSurv <- function(t, shape, scale) pweibull(t, shape=shape,
                                               scale=scale, lower.tail=F)
curve(weibSurv(x, shape=alpha, scale=1/lambda), from=0, to=upper,
      ylim=c(0,1), ylab='Survival probability', xlab='Time')

text(x = 80*.695, y = .84,                # Text with different color & size
     paste0(" alpha ", alpha, "; lambda ", lambda ,""),
     col = "red",
     cex = .7)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
weibHaz <- function(x, shape, scale) dweibull(x, shape=shape, scale=scale)/
                                     pweibull(x, shape=shape, scale=scale,lower.tail=F)

curve(weibHaz(x, shape=alpha, scale=1/lambda), from=0, to=upper,
      ylab='Hazard', xlab='Time', col="red")


par(mfrow=c(1,1))
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

tt.weib <- rweibull(1000, shape=alpha, scale=1/lambda)
mean(tt.weib)
median(tt.weib)
 
#theoretical 
 gamma(1 + 1/alpha)/lambda # mean
 
(log(2)^(1/alpha))/lambda # median
 