
simuldata <- function(n,scenario) {
  # Parametres des simulations 
  myseed <- 123456 # graine
  
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
