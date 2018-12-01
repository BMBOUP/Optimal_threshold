


##{{{ check_effect is a function to check the effect of covariate if they are time-varying or not
Check_effect <- function(time,event,treatment,Marker,data){
  estim <- comp.risk(Event(time,event)~treatment+Marker+
                       treatment*Marker,
                     cause=1,
                     data=data,
                     resample.iid=1,
                     n.sim=1000,
                     model="logistic")
  par(mfrow=c(2,2))
  plot(estim)
}
