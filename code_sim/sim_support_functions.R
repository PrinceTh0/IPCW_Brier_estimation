
#predict survival based on a fitted Cox proportional hazards model
predict.cox <- function(cox_model, newdata, times){
  
  # NULL-model reduces to step function
  if(is.null(cox_model$coefficients)){
    beta <- matrix(0,1,1)
    X <- matrix(0, nrow = nrow(cox_model$x), ncol = 1)
    X_nd <- matrix(0, nrow = nrow(newdata))
  } else{
    beta <- cox_model$coefficients
    X <- cox_model$x
    X_nd <- as.matrix(newdata[, colnames(cox_model$x), with = F])
  }
  
  Y <- cox_model$y
  
  Xbeta <- matrix(X%*%beta, dimnames = list(NULL, "Xbeta"))
  
  Xbeta_nd <- matrix(X_nd%*%beta, dimnames = list(NULL, "Xbeta_nd"))
  
  model_data <- as.data.table(cbind(Y, exp(Xbeta)))
  
  setorder(model_data, time)
  model_data <- model_data[, lapply(.SD, sum), by = time]
  
  # vectorized approach
  atrisk <- rev(cumsum(rev(model_data[,Xbeta])))
  n_fail <- model_data[, status]
  haz <- (n_fail/atrisk)
  
  m_haz <- matrix(rep(haz, length(times)), ncol = length(times), byrow = F)
  m_times <- matrix(rep(times, length(model_data$time)), ncol = length(times), byrow = T)
  m_tevent <- matrix(rep(model_data$time, length(times)), ncol = length(times), byrow = F)
  
  idx <- m_tevent <= m_times
  cum_haz <- colSums(m_haz * idx)
  
  # Breslow estimator
  base_surv <- exp(-cum_haz)
  
  m_Xbeta_nd <- matrix(rep(Xbeta_nd, length(times)), ncol = length(times), byrow = F)
  m_base_surv <- matrix(rep(base_surv, nrow(m_Xbeta_nd)), ncol = length(times), byrow = T)
  
  surv <- m_base_surv ** exp(m_Xbeta_nd)
  
  return("surv" = surv)
}

#predict survival based on a fitted random survival forest model
predict.rsf <- function(rsf_model, newdata, times){
  ptemp <- predict(rsf_model, newdata)$survival
  pos <- sindex(jump.times = rsf_model$time.interest, eval.times = times)
  p <- cbind(1,ptemp)[, pos + 1]
  if (nrow(p) != nrow(newdata) || ncol(p) != length(times)){
    stop("Prediction failed")
  }
  return(p)
}


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


#predict survival based on a fitted elastic net model
predict.glmnet <- function(model, trainingData, testData, times, obs_time){
  
  X_modelData <- data.matrix(trainingData[,grep("(X.*[[:digit:]])",names(trainingData)), with = F])
  X_newData <- data.matrix(testData[,grep("(X.*[[:digit:]])",names(testData)), with = F])
  
  if(obs_time == T){
    Y_modelData <- trainingData[, c("Ttilde","event")]
  } else if(obs_time == F){
    Y_modelData <- trainingData[, c("Ttrue","event_Ttrue")]
  }
  
  
  beta <- model$glmnet.fit$beta[,model$lambda == model$lambda.min]
  
  lp_modelData <- X_modelData%*%beta
  lp_newData <- X_newData%*%beta
  
  atrisk <- rev(cumsum(rev(exp(lp_modelData))))
  n_fail <- Y_modelData$event
  haz <- (n_fail/atrisk)
  
  m_haz <- matrix(rep(haz, length(times)), ncol = length(times), byrow = F)
  if(obs_time == T){
    m_times <- matrix(rep(times, length(Y_modelData$Ttilde)), ncol = length(times), byrow = T)
    m_tevent <- matrix(rep(Y_modelData$Ttilde, length(times)), ncol = length(times), byrow = F)
  } else if(obs_time == F){
    m_times <- matrix(rep(times, length(Y_modelData$Ttrue)), ncol = length(times), byrow = T)
    m_tevent <- matrix(rep(Y_modelData$Ttrue, length(times)), ncol = length(times), byrow = F)
  }
  
  idx <- m_tevent <= m_times
  cum_haz <- colSums(m_haz * idx)
  
  # Breslow estimator
  base_surv <- exp(-cum_haz)
  
  m_lp_newData <- matrix(rep(exp(lp_newData), length(times)), ncol = length(times), byrow = F)
  m_base_surv <- matrix(rep(base_surv, nrow(m_lp_newData)), ncol = length(times), byrow = T)
  
  surv <- m_base_surv ** exp(m_lp_newData)
  
  return("surv" = surv)
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


#predict survival based on a fitted xgboost survival model
predict.xgb <- function(model,status,observed.times,training.data,new.data,times){
  
  #risk scores exp(X'beta), training data
  hr.train <- predict(object = model, newdata = training.data)
  hr.train[is.infinite(hr.train)] <- max(hr.train[!is.infinite(hr.train)])
  #risk scores, new data
  hr.new <- predict(object = model, newdata = new.data)
  hr.new[is.infinite(hr.new)] <- max(hr.new[!is.infinite(hr.new)])
  
  #baseline hazard (Breslow)
  surv.dt <- data.table("time" = observed.times, "state" = status)
  bh.dt <- data.table("time" = times)
  fail.cens.dt <- data.table("time"=rle(observed.times)$values, "fail.cens" = rle(observed.times)$lengths)
  hr.dt <- data.table("time"=observed.times,"risk"=hr.train)
  
  fail.cens.dt <- merge.data.table(fail.cens.dt, surv.dt[, .(nfail = sum(state)), by = time], by = "time") # in case of tied event/censoring times
  fail.cens.dt <- merge.data.table(fail.cens.dt, hr.dt[, .(risk.weight.t = sum(risk)), by = time], by = "time")
  
  indx <- sapply(times, function(x){x < fail.cens.dt$time})
  time2 <- times[apply(indx,1,"which.min")]
  time2[rowSums(indx) == ncol(indx)] <- 9999
  fail.cens.dt$time <- time2
  fail.cens.dt <- fail.cens.dt[, .(fail.cens = sum(fail.cens), nfail = sum(nfail), risk.weight.t = sum(risk.weight.t)), by = time][, natrisk := rev(cumsum(rev(fail.cens)))][, risk.weight := rev(cumsum(rev(risk.weight.t)))]
  
  bh.dt <- merge.data.table(bh.dt,fail.cens.dt,by="time",all.x = T)
  bh.dt[, fail.cens := nafill(fail.cens, "const", 0)][, nfail := nafill(nfail, "const", 0)][, risk.weight.t := nafill(risk.weight.t, "const", 0)]
  bh.dt[, natrisk := nafill(natrisk, "nocb")][, risk.weight := nafill(risk.weight, "nocb")]
  if(sum(is.na(bh.dt$natrisk)>0)){
    bh.dt[, natrisk := nafill(natrisk, "const", bh.dt[which.max(is.na(natrisk))-1]$natrisk - bh.dt[which.max(is.na(natrisk))-1]$nfail)]
  }
  if(sum(is.na(bh.dt$risk.weight)>0)){
    bh.dt[, risk.weight := nafill(risk.weight, "const", bh.dt[which.max(is.na(risk.weight))-1]$risk.weight - bh.dt[which.max(is.na(risk.weight))-1]$risk.weight.t)]
  }
  bh.dt[, bh.weight := nfail/risk.weight]
  bh.dt[is.nan(bh.weight)] <- 0
  
  ht <- outer(hr.new, bh.dt$bh.weight, "*")
  
  #cumulative hazard
  Ht <- t(apply(ht, 1, function(x){cumsum(x)}))
  
  #survival function
  St <- exp(-Ht)
  
  return(list("S_hat"=St))
}

