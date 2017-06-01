## Prediction interval example for Gamma GAM

library(mgcv)

## simulate some data...
f <- function(x) (0.2 * x^11 * (10 * (1 - x))^6 + 10 *
                    (10 * x)^3 * (1 - x)^10)/2
x <- runif(200)
fx <- f(x)
Ey <- exp(fx);scale <- .5 ## mean and GLM scale parameter
## Note that `shape' and `scale' in `rgamma' are almost
## opposite terminology to that used with GLM/GAM...
set.seed(8)
y <- rgamma(Ey*0,shape=1/scale,scale=Ey*scale)

## fit smooth model to x, y data...
model <- gam(y~s(x,k=20),family=Gamma(link=log),method="REML")

## extract parameter estimates and cov matrix...
beta <- coef(model)
Vb <- vcov(model)

## simulate replicate beta vectors from posterior...
Cv <- chol(Vb)
n.rep=10000
nb <- length(beta)
br <- t(Cv) %*% matrix(rnorm(n.rep*nb),nb,n.rep) + beta

## turn these into replicate linear predictors...
xp <- 0:200/200
Xp <- predict(model,newdata=data.frame(x=xp),type="lpmatrix")
lp <- Xp%*%br
fv <- exp(lp) ## ... finally, replicate expected value vectors

## now simulate from Gamma deviates with mean as in fv
## and estimated scale...

yr <- matrix(rgamma(fv*0,shape=1/model$scale,scale=fv*scale),nrow(fv),ncol(fv))

plot(rep(xp,n.rep),yr,pch=".") ## plotting replicates
points(x,y,pch=19,cex=.5) ## and original data

## compute 95% prediction interval...
PI <- apply(yr,1,quantile,prob=c(.025,0.975))
## and plot it...
lines(xp,PI[1,],col=2,lwd=2);lines(xp,PI[2,],col=2,lwd=2)

## Confidence interval for comparison...
pred <- predict(model,newdata=data.frame(x=xp),se=TRUE)
lines(xp,exp(pred$fit),col=3,lwd=2)
u.ci <- exp(pred$fit + 2*pred$se.fit)
l.ci <- exp(pred$fit - 2*pred$se.fit)
lines(xp,u.ci,col=3,lwd=2);lines(xp,l.ci,col=3,lwd=2)

