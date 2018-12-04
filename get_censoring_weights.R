# {{{ inputs description :
# time           : vector of observed failure times
# event          : vector of indicator of status (0 for censoring, 1 for type of event)
# timepoint : prediction time 
# }}}
# {{{ outputs description :
# weights     : inverse of probability of censored weighting
# }}}

Marginal_effecttimepoint<-function(time,event,Treat,timepoint) {
  w <- rep(0,length(time))
  get.censoring.weights <- function(timepoint,time,event){
  untreated <- time >= timepoint 
  censure <- time < timepoint & event==0 
  new.times <- time
  Fitcens <- prodlim::prodlim(Hist(time,event)~1,reverse=TRUE)
  new.times [untreated] <- timepoint
  ## increasing time 
  times_sorted <- sort(new.times,decreasing=FALSE,method='shell')
  order <- rank(new.times)
  weights <- prodlim::predict(Fitcens,times=times_sorted)[order]
  weights <- 1/weights
  weights[censure] <- 0 
   return(weights)
    
         }
  w[Treat==1] <- get.censoring.weights(timepoint,time[Treat==1],event[Treat==1])
  rho_1 <- sum(w*I(time<timepoint)*Treat)/sum(w*Treat)
  w[Treat==0] <- get.censoring.weights(timepoint,time[Treat==0],event[Treat==0])

  rho_0 <- sum(w*I(time<timepoint)*(1-Treat))/sum(w*(1-Treat))
  return(rho_0-rho_1)
  }
##}}}
