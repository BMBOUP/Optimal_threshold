rm(list=ls())

set.seed(1234)

## {{{ technicalities
library(survival)
library(timereg)
library(prodlim)

load("GSE14814.randomise.RData")
d <- GSE14814.randomise
## }}}

## {{{ sample size
n <- 1000
timepoints <- c(2,6) # give in increasing order !!
## }}}

## {{{ estimate the model from the data
estima <- comp.risk(Event(OS_time,OS_status)~Treat+const(Marquer)+
                      const(Treat*Marquer),
                    data=d,
                    cause=1,
                    times=timepoints,
                    resample.iid=1,
                    n.sim=1000,
                    model="logistic") 
## }}}

## {{{ generate biomarker
Y <- sample(x=d$Marquer, size=n, replace=TRUE)
# This is equivalent to generate the values by inverting the empirical cumulative distribution function.
# Indeed, in both case, each new value is drwan among values of d$Marquer, with probability 1*/nrow(d)
## }}}

## {{{ generate treatment
treat <- rbinom(n,size=1,prob=0.5) 
## }}}

## {{{ generate time to event

## {{{ first generate the risk given marker and treatment at each time point
RiskGivenYTreat <- function(x){
  predict(estima,newdata = data.frame(Marquer=x[1],Treat=x[2]),times=timepoints)$P1
} 
AllRisks <- t(apply(cbind(Y,treat),1,RiskGivenYTreat))
colnames(AllRisks) <- paste0("Risk at t=",timepoints)
head(AllRisks)
## }}}

## {{{ Second, create a function which linearly extrapolates between the time points, to get the risk for any time 
ExtrapolatedRisk <- function(t,allRisk,timepoints){
  if(t<timepoints[1]){
    out <- t*allRisk[1]/timepoints[1]
  }else{
    i <- sum(t>=timepoints)
    if(i<length(timepoints)){
      a <-  allRisk[i]
      b <- (allRisk[i+1]-allRisk[i])/(timepoints[i+1]-timepoints[i])
      if(b==0){b <- 0.001} # we need to force the risk to be increasing, to avoid computational issue
      out <- a + b*(t-timepoints[i])
    }else{
      a <-  allRisk[i-1]
      b <- (allRisk[i]-allRisk[i-1])/(timepoints[i]-timepoints[i-1])
      out <- a + b*(t-timepoints[i-1])
    }
  }
  min(out,1)
}

# example of the curve of predicted risk given marker and treatment
whichSubject <- 30 
ttt <- seq(from=0,to=30,by=0.01)
yyy <- sapply(ttt,function(x){ExtrapolatedRisk(x,AllRisks[whichSubject,],timepoints)})
plot(ttt,yyy,type="l",ylab="predicted risk", xlab="time", axes=FALSE)
points(timepoints,AllRisks[whichSubject,], col="red", pch=16)
axis(1)
axis(2,las=1)
## }}}

## {{{ Generate time to event
Temps <- rep(NA,n)
U <- runif(n)
for(i in 1:n){
  Temps[i] <- uniroot(f=function(t)(ExtrapolatedRisk(t,AllRisks[i,],timepoints)-U[i]),
                      lower=0,
                      upper=1000,
                      tol=0.0001)$root
}
# note: actually we do not really need to use uniroot in that case, because we can find a simple formula
# for the root of the equation (as it is piecewise linear).
## }}}


## }}}


## {{{ generate censoring time
AllTimeCens <- unique(sort(d$OS_time[d$OS_status==0]))
FitCens <- prodlim(Hist(OS_time,OS_status)~1,data=GSE14814.randomise,reverse=TRUE)
survCens <- predict(FitCens,times=AllTimeCens)
Jumps <- -diff(c(1,survCens,0))
cbind(c(AllTimeCens,max(d$OS_time)+1),Jumps)
length(AllTimeCens)
length(Jumps)
sum(Jumps)
MyTcens <- sample(x=c(AllTimeCens,max(d$OS_time)+1),size=n,replace=TRUE,prob=Jumps)
# the idea is simply that inverting the estimated empirical cumulative distribution function (ecdf) is equivalent to draw censored time, among those observed int the real dataset, whith probability equal to the jump of the Kaplan-Meier estimator of the censoring distribution.
# becase Kaplan-Meier does not always reach 0, I added the time max(d$OS_time)+1 (beyond the last observed time in the data set), which is drawn with prob 1-sum(jumps). 


## Just to Check that it works
plot(FitCens,confint=FALSE)
for(i in 1:10){
  Tcens <- sample(x=c(AllTimeCens,max(d$OS_time)+1),size=n,replace=TRUE,prob=Jumps)
  di <- data.frame(time=Tcens,status=rep(0,n))    
  plot(prodlim(Hist(time,status)~1,data=di,reverse=TRUE),add=TRUE,col="red")
}
plot(FitCens,add=TRUE,confint=FALSE)
## }}}


## {{{ Final simulated data
dataSimul <- data.frame(time=pmin(Temps,MyTcens),
                        status=as.numeric(Temps<=MyTcens),
                        treat=treat,
                        Y=Y)
table(dataSimul$status)
table(d$OS_status)
head(dataSimul)
## }}}


## {{{ Check that estimates from simulated data match with those of the data set

## {{{ plot
plot(prodlim(Hist(time,status)~1,data=dataSimul),col="red",confint=FALSE)
plot(prodlim(Hist(OS_time,OS_status)~1,data=d),col="green", add=TRUE,confint=FALSE)
abline(v=timepoints)
points(timepoints,1-colMeans(AllRisks),cex=3,pch=16)
# Here the two Kaplan-Meier curve should (approximately) be equal
# at the times corresponding to timepoints 
## }}}

## {{{ check values
estima1 <- comp.risk(Event(time,status)~treat+const(Y)+
                       const(treat*Y),
                     data=dataSimul,
                     cause=1,
                     times=timepoints,
                     resample.iid=1,
                     n.sim=5,
                     model="logistic")
estima$cum
estima1$cum
## }}}

## }}}
