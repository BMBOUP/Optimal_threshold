# This function returns the marginal treatment effect rho_0(t)-rho_1(t) at time t.
# if positive the standard decision is to treat everybody
# if negative the standard decision is to avoid treatmentfor everybody
# {{{ inputs description :
# time           : vector of observed failure times
# event          : vector of indicator of status (0 for censoring, 1 for type of event)
# Treat          : 1=treated, 0=untreated
# timepoint      : prediction time 
# }}}
# {{{ outputs description :
# rho_0(t)-rho_1(t)   : the marginal effect of treatment at time t=timepoint  without the biomarker
# }}}


Marginal_effecttimepoint<-function(time,event,Treat,timepoint) {
  if (length(event)!=length(time) | length(Treat)!=length(time) | length(event)!=length(Treat))
  {stop("lengths of vector time, event,treatment have to be equal\n") }
  if (missing(timepoint)){
    stop("Choose one time for horizon prediction \n") 
  }
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
  w <- rep(0,length(time))
  w[Treat==1] <- get.censoring.weights(timepoint,
                                        time[Treat==1],
                                       event[Treat==1])
  rho_1 <- sum(w*I(time<timepoint)*Treat)/sum(w*Treat)
  
  w <- rep(0,length(time))
  w[Treat==0] <- get.censoring.weights(timepoint,
                                       time[Treat==0],
                                      event[Treat==0])
  
  rho_0 <- sum(w*I(time<timepoint)*(1-Treat))/sum(w*(1-Treat)) 
  return(rho_0-rho_1)
}
##}}}
