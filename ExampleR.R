## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ---- warning=FALSE------------------------------------------------------
library(timereg)
library(survival)
library(prodlim)
library(boot)

# source the required R functions
## ------------------------------------------------------------------------
source("Marginal_effect_timepoint.R")
source("get_estimates.R")
source("plot_time_dependent_predictiveness.R")
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


## ------------------------------------------------------------------------
Marginal_effecttimepoint(JBR.10_data$OS_time,JBR.10_data$OS_status,JBR.10_data$Treat,2)

Marginal_effecttimepoint(csl_data$years,csl_data$dc,csl_data$tment,5)

## ------------------------------------------------------------------------

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

## ------------------------------------------------------------------------
pdf("test1.pdf")
plottime_predictiveness_curve(Object=risk1,timepoint=5,Marker=csl_data$pro)
plottime_predictiveness_curve(Object=risk2,timepoint=2,Marker=JBR.10_data$Marquer)
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
## theoretical optimal threshold 
 
result1<-boot(data=JBR.10_data,statistic=marker_value,R=999,formula=Event(OS_time,OS_status)~Treat+const(Marquer)+const(Treat*Marquer))
ci1 <- boot.ci(result1,conf=0.95, type="norm")
 
paste(result1$t0,"(", ci1$normal[2],",",ci1$normal[3],")")

## ------------------------------------------------------------------------
## theoretical Optimal threshold 
csl_data$Marquer <- csl_data$pro
result2 <- boot(data=csl_data,statistic=marker_value,R=999,
                formula=Event(years,dc)~tment+const(pro)+
                  const(tment*pro))
ci2 <- boot.ci(result2,conf=0.95, type="norm")
 
paste(result2$t0,"(", ci2$normal[2],",",ci2$normal[3],")")

## ------------------------------------------------------------------------

