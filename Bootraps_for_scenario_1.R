########## Simulation for scenario 1
rm(list=ls())
library(timereg)
library(boot)

n <- 5000     #  sample size
myseed <- 123456  # seed 

beta1 <- -0.29  # coefficient of the treatment
beta2 <- 0.6    # coefficient of the marker
beta3 <- -1.5   # interaction between treatment and marker

a <- 10  # constant using in the model

#{{{ function to get optimal threshold
marker_value <- function(formula,data,indices){
  d=data[indices,]
  estim <- comp.risk(formula,data=d, cause=1,
                     resample.iid=1,
                     times=timepoint,
                     n.sim = 1000,
                     model="logistic")
  return(-estim$gamma[1]/estim$gamma[3])
  
}

#}}}
infc_2 <- NA
supc_2 <- NA
biasc_2 <- NA
Standc_2 <- NA
Valc_2 <- NA
tCI.total <- 0
timepoint <- 2.43
true.mean <- (-beta1/beta3) # true threshold

for(j in 1:1000)  {
  # {{{ Simulation (une base)
  set.seed(myseed+j)
  
  
  # {{{
  
  Y <- rnorm(n,0,1) # marker
  C <- rexp(n,0.2) # censured

  Treat <- rbinom(n,size=1,prob=0.5)  # le traitement avec P(trt=1)=0.5 
  
  #}}}
  
  #{{{ 
  
  U <- runif(n)
  Tsur <- a*(U/(1-U))*exp(-(beta1*Treat+beta2*Y+beta3*Treat*Y))

  ##}}}
  
  delta <- as.numeric(Tsur<=C)   
  Tsuiv <- pmin(Tsur,C)   
  
  donne<- as.data.frame(cbind(time=Tsuiv,status=delta,Y=Y,Treat=Treat)) 
  
 
  
  foo <- boot(data=donne,statistic = marker_value,R=1000,
              formula=Event(time,status)~const(Treat)+const(Y)+
                const(Treat*Y))
  ci <- boot.ci(foo,conf=0.95, type="norm")
  tCI <- ci$normal[2:3]
  
  infc_2[j] <- ci$normal[2]
  supc_2[j] <- ci$normal[3]
  biasc_2[j] <- mean(foo$t[,1])-true.mean
  Valc_2[j] <- mean(foo$t[,1])
  Standc_2[j] <- sd(foo$t[,1])
  
  if (true.mean > min(tCI) && true.mean < max(tCI)){ 
    tCI.total <- tCI.total + 1
  }
  
} 

tauxc_2 <- (tCI.total/1000)
save(infc_2,file="infc_2.RData")
save(supc_2,file="supc_2.RData")
save(biasc_2,file="biasc_2.RData")
save(Standc_2,file="Standc_2.RData")
save(tauxc_2,file="tauxc_2.RData")
save(Valc_2,file="Valc_2.RData") 