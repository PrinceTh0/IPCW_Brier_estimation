
#predict survival based on a fitted mboost model
predict.mboost <- function(mboost_model, newdata, times){
  mb_survFit <- survFit(object = mboost_model, newdata = newdata)
  mb_surv <- mb_survFit$surv
  mb_time <- mb_survFit$time
  
  mb_surv <- t(mb_surv)
  mb_surv <- cbind(c(rep(1, nrow(mb_surv))), mb_surv)

  min_event <- sapply(times, function(t){t < mb_time})
  
  m_mb_out <- mb_surv[,apply(min_event,2,"which.max")] # for each ref t, look for smallest event time satisfying the condition ref_t < event_t
  m_mb_out[,colSums(min_event) == 0] <- mb_surv[,ncol(mb_surv)]
  
  return(m_mb_out)
}


#predict survival based on a fitted ranger survival forest model
predict.rangersf <- function(model, newdata, times){
  sf_pred <- predict(object = model, data = newdata)
  sf_surv <- sf_pred$survival
  sf_time <- sf_pred$unique.death.times
  
  sf_surv <- cbind(c(rep(1, nrow(sf_surv))), sf_surv)
  
  min_event <- sapply(times, function(t){t < sf_time})
  
  m_sf_out <- sf_surv[,apply(min_event,2,"which.max")] # for each ref t, look for smallest event time satisfying the condition ref_t < event_t
  m_sf_out[,colSums(min_event) == 0] <- sf_surv[,ncol(sf_surv)]
  
  return(m_sf_out)
}


#predict survival based on a fitted xgboost survival model (discrete timescale with paired event times)
predict.xgb <- function(model,status,observed.times,training.data,new.data,times){
  
  #risk scores exp(X'beta), training data
  hr.train <- predict(object = model, newdata = training.data)
  #risk scores, new data
  hr.new <- predict(object = model, newdata = new.data)
  
  #baseline hazard (Breslow)
  surv.dt <- data.table("time" = observed.times, "state" = status)
  bh.dt <- data.table("time" = times)
  nfail.dt <- data.table("time"=rle(observed.times)$values, "nfail" = rle(observed.times)$lengths)
  
  bh.dt <- merge.data.table(bh.dt,nfail.dt,by="time",all.x = T)
  bh.dt[time == 0, nfail := 0][, nfail := nafill(nfail, "const", 0)][, natrisk := rev(cumsum(rev(nfail)))]
  bh.dt <- merge.data.table(bh.dt, surv.dt[state == 1, .(nevent = .N), by = time], by = "time", all.x = T)
  bh.dt[, nevent := nafill(nevent, "const", 0)][, na := nevent/natrisk][, cumna := cumsum(na)]
  
  hr.dt <- data.table("time"=observed.times,"risk"=hr.train)[, .(risk.weight.t = sum(risk)), by = time]
  bh.dt <- merge.data.table(bh.dt, hr.dt, by = "time", all.x = T)
  bh.dt[, risk.weight.t := nafill(risk.weight.t, "const", 0)]
  bh.dt[, risk.weight := rev(cumsum(rev(risk.weight.t)))][, bh.weight := nevent/risk.weight][, cumbh.weight := cumsum(bh.weight)]
  
  ht <- outer(hr.new, bh.dt$bh.weight, "*")
  
  #cumulative hazard
  Ht <- t(apply(ht, 1, function(x){cumsum(x)}))
  
  #survival function
  St <- exp(-Ht)
  
  return(list("S_hat"=St))
}


#predict survival curves based on a marginal model (KM)
predict.km <- function(model, times){
  
  survfit.model <- survfit(model)
  survfit.dt <- data.table("surv"=survfit.model$surv, "times"=survfit.model$time)
  times.dt <- data.table("times"=times)
  
  surv.dt <- merge.data.table(times.dt, survfit.dt, by = "times", all.x = TRUE)
  surv.dt[times == 0, surv := 1]
  surv.dt[, surv := nafill(surv, "locf")]
  
  return(surv.dt$surv)
}
