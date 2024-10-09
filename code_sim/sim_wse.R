
wse.calc <- function(bs.in, uncensored.in, times.in, ref.times.in){
  
  mat_ref_times <- matrix(rep(t(ref.times.in), length(times.in)), nrow = length(times.in), byrow=T)
  mat_times <- matrix(rep(times.in, length(ref.times.in)), ncol = length(ref.times.in), byrow=F)
  
  loc <- mat_times <= mat_ref_times
  pre_wgt <- colSums(loc)
  wgt <- c(pre_wgt[1],diff(pre_wgt, lag = 1))
  
  wse <- sqrt(sum((((bs.in - uncensored.in)**2)*wgt)/sum(wgt)))
  
  return(wse)
  
}