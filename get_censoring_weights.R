

##{{{ This function give the inverse  inverse probability of censoring weighting
## timepoint : is the prediction time
## time and event are  observed time and status of event. 
library(prodlim)
get.censoring.weights <- function(timepoint,time,event){
  untreated <- time >= timepoint 
  censure <- time < timepoint & event==0 
  new.times <- time
  Fitcens <- prodlim(Hist(time,event)~1,reverse=TRUE)
  new.times [untreated] <- timepoint
  ## increasing time 
  times_sorted <- sort(new.times,decreasing=FALSE,method='shell')
  order <- rank(new.times)
  weights <- predict(Fitcens,times=times_sorted)[order]
  weights <- 1/weights
  weights[censure] <- 0 
  return(weights)
}
##}}}