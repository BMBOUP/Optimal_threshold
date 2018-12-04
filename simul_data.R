# {{{ Simulate right censored survival data with two covariates treatment and biomarker which follows respectively 
# bernouilli distribution with probability 1/2 and standard normal distribution. The model based is the time-dependent
# logistic model: logit(P(D(t)=1|Treat,Y)=log(t/a)+\beta1(t)Treat+beta2(t)Y+beta3(t)Treat*Y with different scenario
# scenario   : scenario=1 : constant effect 
#              logit(P(D(t)=1|Treat,Y)=log(t/a)-0.29*Treat+0.6*Y-1.5*Treat*Y 
#            : scenario=2 :  increasing effect of treatment with time.
#              logit(P(D(t)=1|Treat,Y)=log(t/a)+beta1(t)*Treat+0.6*Y-1.5*Treat*Y where beta1(t)=-1.24 if t<3 and -0.31 otherwise
#            : scenario=3 : decreasing treatment effect in the model is  with time  
#              logit(P(D(t)=1|Treat,Y)=log(t/a)+log(1/(t+1))*Treat+0.6*Y-2*Treat*Y 
# }}}
#{{{
#inputs description:
# n          : sample size
# scenario   : 1,2,3

# {{{ outputs description:
# donne      : a dataframe to one of the scenario.
# }}}


simuldata <- function(n,scenario) {
  #  simulations parameters 
  myseed <- 123456 # seed
  
  # coefficients of time-dependent logistic model 
  beta2 <- 0.6     # coefficients of biomarkerY
  beta3 <- -1.5  # interaction between marker and treatment 
  a <- 10       # nombre de subdivision
  
  set.seed(myseed) 
  
  Y <- rnorm(n,0,1) # biomarker value
  C <- rexp(n,0.2)  # censoring distribution
  
  Treat <- rbinom(n,size=1,prob=0.5) # treatment indicator
  
  # }}}
  if(scenario==1){
    beta1 <- -0.29 
    U <- runif(n)
    Tsur <- a*(U/(1-U))*exp(-(beta1*Treat+beta2*Y+beta3*Treat*Y))
  }
  
  if(scenario==2) {
    # piecewice constant treatment effect 
    beta1<-function(x,a=-1.24,b=-0.31,seuil=3) # effect of treatment   
    {
      return(ifelse(x<=seuil,a,b))
    }
    
    f <- function(x,i){
      prob <-exp(log(x/a)+beta1(x)*Treat[i]+beta2*Y[i]+beta3*Treat[i]*Y[i])/(1+exp(log(x/a)+beta1(x)*Treat[i]+beta2*Y[i]+beta3*Treat[i]*Y[i]))
      return(prob)
    }
    
    
    U <- runif(n)
    Tsur <- c()
    for(i in 1:n){
      Tsur[i] <- uniroot(function(t)(f(t,i)-U[i]),c(0,15),f.lower=-1,f.upper=1,tol=0.0001)$root
    }
    
  }
  
  if(scenario==3){
    beta3 <- -2
    f <- function(x,i){
      prob <-exp(log(x/a)+log(1/(x+1))*Treat[i]+beta2*Y[i]+beta3*Treat[i]*Y[i])/(1+exp(log(x/a)+log(1/(x+1))*Treat[i]+beta2*Y[i]+beta3*Treat[i]*Y[i]))
      return(prob)
    }
    
    U <- runif(n)
    Tsur <- c()
    for(i in 1:n){
      Tsur[i] <- uniroot(function(t)(f(t,i)-U[i]),c(0,15),f.lower=-1,f.upper=1,tol=0.0001)$root
    }
    
  }
  
  delta <- as.numeric(Tsur<=C)   
  Tsuiv <- pmin(Tsur,C)   
  
  donne<- as.data.frame(cbind(time=Tsuiv,status=delta,Y=Y,Treat=Treat)) 
  return(donne)
}
