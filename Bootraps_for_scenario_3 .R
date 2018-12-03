########## Simulation for scenario 3 ######
rm(list=ls())
library(timereg)
library(boot)
# {{{ Parametres des simulations 


n <- 5000     #  sample size
myseed <- 123456  # seed

# coefficients 
beta2 <- 0.6        # coefficient of the marker
beta3 <- -2        # interaction between marker and treatment

a <- 10  # a constant using in the model 



##{{{ function to get optimal threshold in the boot
marker_value <- function(formula,data,indices){
  d=data[indices,]
  estim <- comp.risk(formula,data=d, cause=1,
                     resample.iid=1,
                     times=timepoint,
                     n.sim = 1000,
                     model="logistic")
  return(-estim$cum[3]/estim$gamma[2])
  
}






infd_1 <- NA
supd_1 <- NA
Vald_1 <- NA
biasd_1 <- NA
Standd_1 <- NA
tCI.total <- 0
timepoint <- 1
true.mean <- (-log(1/(1+timepoint))/beta3) # true threshold 

for(j in 1:1000)  {
  # {{{ Simulation (une base)
  set.seed(myseed+j)

# {{{Etape 3

Y <- rnorm(n,0,1) # Marker
C <- rexp(n,0.2) # Censored

Treat <- rbinom(n,size=1,prob=0.5) # treatment  

# }}} function of linear prediction 

f <- function(x,i){
  prob <-exp(log(x/a)+log(1/(x+1))*Treat[i]+beta2*Y[i]+beta3*Treat[i]*Y[i])/(1+exp(log(x/a)+log(1/(x+1))*Treat[i]+beta2*Y[i]+beta3*Treat[i]*Y[i]))
  return(prob)
}

###{{{ 

U <- runif(n)
Tsur <- c()
for(i in 1:n){
  Tsur[i] <- uniroot(function(t)(f(t,i)-U[i]),c(0,15),f.lower=-1,f.upper=1,tol=0.0001)$root
}

##}}}

delta <- as.numeric(Tsur<=C)   
Tsuiv <- pmin(Tsur,C)   

donne<- as.data.frame(cbind(time=Tsuiv,status=delta,Y=Y,Treat=Treat)) 

foo <- boot(data=donne,statistic = marker_value,R=1000,
            formula=Event(time,status)~Treat+const(Y)+
              const(Treat*Y))
ci <- boot.ci(foo,conf=0.95, type="norm")
tCI <- ci$normal[2:3]

infd_1[j] <- ci$normal[2]
supd_1[j] <- ci$normal[3]
biasd_1[j] <- mean(foo$t[,1])-true.mean
Vald_1[j] <- mean(foo$t[,1])
Standd_1[j] <- sd(foo$t[,1])
# coverage probability
if (true.mean > min(tCI) && true.mean < max(tCI)){ 
  tCI.total <- tCI.total + 1
}

} 


tauxd_1 <- (tCI.total/1000)
save(outd1_1,file="outd1_1.RData")
save(outd2_1,file="outd2_1.RData")
save(infd_1,file="infd_1.RData")
save(supd_1,file="supd_1.RData")
save(biasd_1,file="biasd_1.RData")
save(Standd_1,file="Standd_1.RData")
save(tauxd_1,file="tauxd_1.RData")
