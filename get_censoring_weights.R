# {{{ inputs description :
# time      :
# event     :
# timepoint : prediction time 

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
