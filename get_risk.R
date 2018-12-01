
##{{{
## function to get risk on treated and untreated subjects using time dependent logistic model


get_risk <- function(time,event,Treat,Marker,varying,timepoint)
{
  if (length(event)!=length(time) | length(Marker)!=length(time)|length(Treat)!=length(time) | length(event)!=length(time)){
    stop("lengths of vector time, event,treatment and Marker have to be equal\n") }
  if (missing(timepoint)){
    stop("Choose one time for horizon prediction \n") }
  
  data <- cbind.data.frame(time=time,event=event,Treat=Treat,Marker=Marker)
  
  if(varying==TRUE){
    
    estim <- comp.risk(Event(time,event)~Treat+const(Marker)+
                         const(Treat*Marker),
                       cause=1,
                       data=data,
                       times=timepoint,
                       resample.iid=1,
                       n.sim=1000,
                       model="logistic")
  }
  
  if(varying==FALSE || missing(varying)) {
    estim <- comp.risk(Event(time,event)~const(Treat)+const(Marker)+
                         const(Treat*Marker),
                       cause=1,
                       data=data,
                       resample.iid=1,
                       times=timepoint,
                       n.sim=1000,
                       model="logistic")
    
    
  }
  
  predittreat1 <- predict(estim,newdata=data.frame(Marker=sort(Marker),Treat=1),
                          n.sim=1,times=timepoint)$P1
  predittreat0<- predict(estim,newdata=data.frame(Marker=sort(Marker),Treat=0),
                         n.sim=1,times=timepoint)$P1
  
  risk <- cbind.data.frame(risk_0=predittreat0,risk_1=predittreat1)  
  
  return(risk)
}

## }}}
