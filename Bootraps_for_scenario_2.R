########## Simulation of scenarion 2 
rm(list=ls())
library(timereg)
library(boot)


n <- 5000     #  sample size
myseed <- 123456  # seed

beta2 <- 0.6     # coefficient of the marker
beta3 <- -1.5  # interaction between marker and treatment
a <- 10    # constant using in the model 

##{{{ function to get optimal threshold
marker_value <- function(formula,data,indices){
  d=data[indices,]
  estim <- comp.risk(formula,data=d, cause=1,
                     resample.iid=1,
                     times=timepoint,
                     n.sim = 1000,
                     model="logistic")
  return(-estim$cum[3]/estim$gamma[2])
  
}

# coefficient of the treatment
beta1<-function(x,a=-1.24,b=-0.31,seuil=3)   
{
  return(ifelse(x<=seuil,a,b))
}


##}}}





out1_1 <- NA
out2_1 <- NA
inf_1 <- NA
sup_1 <- NA
bias_1 <- NA
Stand_1 <- NA
tCI.total <- 0
timepoint <- 1
true.mean <- (-beta1(timepoint)/beta3) # true threshold
for(j in 1:1000)  {
  
  set.seed(myseed+j)
  
#{{{
  Y <- rnorm(n,0,1) # marker
  C <- rexp(n,0.2) # censored

  Treat <- rbinom(n,size=1,prob=0.5)  # treatment
  
# }}}
 #{{{ 
  f <- function(x,i){
    prob <-exp(log(x/a)+beta1(x)*Treat[i]+beta2*Y[i]+beta3*Treat[i]*Y[i])/(1+exp(log(x/a)+beta1(x)*Treat[i]+beta2*Y[i]+beta3*Treat[i]*Y[i]))
    return(prob)
  }
  
  
  U <- runif(n)
  Tsur <- c()
  for(i in 1:n){
    Tsur[i] <- uniroot(function(t)(f(t,i)-U[i]),c(0,15),f.lower=-1,f.upper=1,tol=0.0001)$root
  }

#}}}
  
  delta <- as.numeric(Tsur<=C)   
  Tsuiv <- pmin(Tsur,C)   
  
  donne<- as.data.frame(cbind(time=Tsuiv,status=delta,Y=Y,Treat=Treat)) 

  foo <- boot(data=donne,statistic = marker_value,R=1000,
              formula=Event(time,status)~Treat+const(Y)+
                const(Treat*Y))
  ci <- boot.ci(foo,conf=0.95, type="norm")
  tCI <- ci$normal[2:3]
  
  inf_1[j] <- ci$normal[2]
  sup_1[j] <- ci$normal[3]
  bias_1[j] <- mean(foo$t[,1])-true.mean
  Stand_1[j] <- sd(foo$t[,1])
  
  if (true.mean > min(tCI) && true.mean < max(tCI)){ 
    tCI.total <- tCI.total + 1
  }
  
} 
taux_1 <- (tCI.total/1000)
save(inf_1,file="inf_1.RData")
save(sup_1,file="sup_1.RData")
save(bias_1,file="bias_1.RData")
save(Stand_1,file="Stand_1.RData")
save(taux_1,file="taux_1.RData")

