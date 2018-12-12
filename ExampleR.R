## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ---- warning=FALSE------------------------------------------------------
library(timereg)
library(survival)
library(prodlim)
library(boot)

# source the required R functions
#source("")
## ------------------------------------------------------------------------

load("JBR.10_data.RData")
head(JBR.10_data)

## ------------------------------------------------------------------------
load("csl_data.RData")
head(csl_data)

## ------------------------------------------------------------------------
Check_effect <- function(time,event,treatment,Marker){
  data <- cbind.data.frame(time=time,event=event,treatment=treatment,Marker=Marker)
estim <- timereg::comp.risk(Event(time,event)~treatment+Marker+
              treatment*Marker,
            cause=1,
            data=data,
            resample.iid=1,
            n.sim=1000,
            model="logistic")
 par(mfrow=c(2,2))
plot(estim)
}

Check_effect(JBR.10_data$OS_time,JBR.10_data$OS_status,JBR.10_data$Treat,
            JBR.10_data$Marquer)

## ------------------------------------------------------------------------
source("Marginal_effect_timepoint.R")

## ------------------------------------------------------------------------
Marginal_effecttimepoint(JBR.10_data$OS_time,JBR.10_data$OS_status,JBR.10_data$Treat,2)

Marginal_effecttimepoint(csl_data$years,csl_data$dc,csl_data$tment,5)

## ------------------------------------------------------------------------
source("get_estimates.R")
## ------------------------------------------------------------------------
risk1 <- get_estimates(csl_data$years,csl_data$dc,csl_data$tment,
                      csl_data$pro,varying=TRUE,5)

risk1$threshold
risk1$confint
risk1$Pneg

risk2 <- get_estimates(JBR.10_data$OS_time,JBR.10_data$OS_status,JBR.10_data$Treat,JBR.10_data$Marquer,varying=TRUE,2)
risk2$threshold
risk2$confint
risk2$Pneg

## ------------------------------------------------------------------------
source("plot_time_dependent_predictiveness.R")
## ------------------------------------------------------------------------
pdf("test1.pdf")
plotime_predictiveness_curve(Object=risk1,timepoint=5,Marker=csl_data$pro)
plotime_predictiveness_curve(Object=risk2,timepoint=2,Marker=JBR.10_data$Marquer)
dev.off()

## ------------------------------------------------------------------------
## function to bootstraps 
timepoint <- 2
marker_value <- function(formula,data,indices){
    d=data[indices,]
    estim <- comp.risk(formula,data=d, cause=1,
                       resample.iid=1,
                       times=timepoint,
                       n.sim =1000,
                       model="logistic")
     out <- -estim$cum[3]/estim$gamma[2]
 if (out >= min(d$Marquer) && out <= max(d$Marquer))
           return(out)
                 else 
          return(NA)
}

## ------------------------------------------------------------------------
## theoretical optimal threshol 
 
result1<-boot(data=JBR.10_data,statistic=marker_value,R=999,formula=Event(OS_time,OS_status)~Treat+const(Marquer)+const(Treat*Marquer))
ci1 <- boot.ci(result1,conf=0.95, type="norm")
 
paste(result1$t0,"(", ci1$normal[2],",",ci1$normal[3],")")

## ------------------------------------------------------------------------
## theoretical Optimal threshol 
csl_data$Marquer <- csl_data$pro
result2 <- boot(data=csl_data,statistic=marker_value,R=999,
                formula=Event(years,dc)~tment+const(pro)+
                  const(tment*pro))
ci2 <- boot.ci(result2,conf=0.95, type="norm")
 
paste(result2$t0,"(", ci2$normal[2],",",ci2$normal[3],")")

## ------------------------------------------------------------------------
## Extension data

set.seed(1234)
d <- JBR.10_data

## {{{ sample size
n <- 5000
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
FitCens <- prodlim(Hist(OS_time,OS_status)~1,data=d,reverse=TRUE)
survCens <- predict(FitCens,times=AllTimeCens)
Jumps <- -diff(c(1,survCens,0))
head(cbind(c(AllTimeCens,max(d$OS_time)+1),Jumps))
MyTcens <- sample(x=c(AllTimeCens,max(d$OS_time)+1),size=n,replace=TRUE,prob=Jumps)
# the idea is simply that inverting the estimated empirical cumulative distribution function (ecdf) is equivalent to draw censored time, among those observed int the real dataset, whith probability equal to the jump of the Kaplan-Meier estimator of the censoring distribution.
# becase Kaplan-Meier does not always reach 0, I added the time max(d$OS_time)+1 (beyond the last observed time in the data set), which is drawn with prob 1-sum(jumps). 


## Just to Check that it works
pdf("p3.pdf")
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

