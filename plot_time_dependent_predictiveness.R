
# {{{ input description :

# fitted         : vector of observed failure times
# delta          : vector of indicator of status (0 for censoring, 1 for type of event one, 2 for type of event two and so on...)
# marker         : vector ofmarker values


## function to display time-dependent marker-by-treatment predictiveness curves
plotime_predictiveness_curve <- function(fitted,timepoint,Marker){
  xlim <- c(0,100)
  breaks = seq(xlim[1], xlim[2], length.out = 5)
  k <- length(fitt$risk_0)
  x1<- sort(Marker)
  y1 <- (((1:k)*100)/k)
  z1 <- fitt$risk_0
  z2 <- fitt$risk_1
  plot(y1,z1,type='l',col="red",xlab="(%)population \n marker Value",ylab="Risk given marker (%)",
       main=paste("Time dependent marker-by-treatment \n predictiveness curve at times t=",timepoint,"years"),ylim=c(0,1),axes=FALSE)
  lines(y1,z2,type='l',col="blue",lwd=2)
  axis(1)
  axis(2,at=seq(from=0,to=1,by=0.2),labels=seq(from=0,to=1,by=0.2)*100,las=2)
  axis(1,at= breaks, label = round(quantile(x1, prob = breaks/100), 1),pos=-0.26)   
}
##}}}
