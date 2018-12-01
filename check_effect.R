

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
