# This function estimate the optimal threshold with confidence interval, the percentage of subject
# and risk given biomarker and treatment at given prediction time.
# {{{ input description :
# time           : vector of observed failure times
# event          : vector of indicator of status (0=censored, 1=uncensored)
# Treat          : vector of treatment indicator 1 if treated and 0 others 
# Marker         : vector of continuous biomarker values
# varying        : boolean indicator of time-varying treatment effect
# timepoint      : prediction time

# }}}
# {{{ output description : a list of
# threshold         : the optimal theshold of the biomarker (or cut-off)
# confint           : the confidence intervale low and upper
# Pneg              : the pourcentage of subjects who avoid or take the treatment à time t
# risk_0            : P(D(t)=1|Treat=0,Marker)
# risk_1            : P(D(t)=1|Treat=1,Marker)
# }}}
get_risk <- function(time,event,Treat,Marker,varying,timepoint)
{
  if (length(event)!=length(time) | length(Marker)!=length(time)|length(Treat)!=length(time) | length(event)!=length(time))
  {stop("lengths of vector time, event,treatment and Marker have to be equal\n") }
  if (missing(timepoint)){
    stop("Choose one time for horizon prediction \n") }
  
  data <- cbind.data.frame(time=time,event=event,Treat=Treat,Marker=Marker)
  
  if(varying==FALSE || missing(varying)) {
    estim <- timereg::comp.risk(Event(time,event)~const(Treat)+const(Marker)+
                                  const(Treat*Marker),
                                cause=1,
                                data=data,
                                resample.iid=1,
                                times=timepoint,
                                n.sim=1000,
                                model="logistic")
  }
    if(varying==TRUE){
      estim <- timereg::comp.risk(Event(time,event)~Treat+const(Marker)+
                                    const(Treat*Marker),
                                  cause=1,
                                  data=data,
                                  times=timepoint,
                                  resample.iid=1,
                                  n.sim=1000,
                                  model="logistic")
    }
  predittreat1 <- timereg::predict(estim,newdata=data.frame(Marker=sort(Marker),Treat=1),
                                   n.sim=1,times=timepoint)$P1
  predittreat0<- timereg::predict(estim,newdata=data.frame(Marker=sort(Marker),Treat=0),
                                  n.sim=1,times=timepoint)$P1
  Delta <- predittreat0-predittreat1
  Pneg <- length(Delta[Delta<0])/length(Delta)
  
  f <- function(data,indice){
    vec=data[indice]
    return(quantile(vec,Pneg,names=FALSE))
  }
  
  result <- boot::boot(data=data$Marker,statistic=f,R=999)
  ci <- boot::boot.ci(result,conf=0.95, type="norm")
  risk <- list(threshold=round(result$t0,2),
                 confint=paste("(", round(ci$normal[2],2),",",round(ci$normal[3],2),")"),
                    Pneg=(Pneg*100), 
                  risk_0=predittreat0,
                  risk_1=predittreat1)  
  return(risk)
}
