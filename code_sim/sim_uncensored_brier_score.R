
uncensored_brier_score <- function(validation, training, s, ref_times, surv_formula_true, event_estimator){
  
  trainingData <- copy(training)
  evalData <- copy(validation)
  
  #order by true event times
  setorder(trainingData, Ttrue)
  setorder(evalData, Ttrue)
  
  if(event_estimator == "xgboost"){
    xgb.trainingData <- trainingData[, c("Ttrue",names(trainingData)[grep("(X.*[[:digit:]])",names(trainingData))]), with = F]
    xgb.trainingData.dummy <- model.Matrix(object = Ttrue ~ ., data = xgb.trainingData, sparse = T)[,-1]
    xgb.trainingData.data <- xgb.DMatrix(xgb.trainingData.dummy, label = xgb.trainingData$Ttrue)
    
    xgb.evalData <- evalData[, c("Ttrue",names(evalData)[grep("(X.*[[:digit:]])",names(evalData))]), with = F]
    xgb.evalData.dummy <- model.Matrix(object = Ttrue ~ ., data = xgb.evalData, sparse = T)[,-1]
    xgb.evalData.data <- xgb.DMatrix(xgb.evalData.dummy, label = xgb.evalData$Ttrue)
  }
  
  #specify estimation method
  if(event_estimator == "cox"){
    surv_model_true <- coxph(surv_formula_true, data = trainingData, x = T, y = T)
    S_hat_true <- predict.cox(surv_model_true, evalData, ref_times)
  } else if(event_estimator == "rsf"){
    require(randomForestSRC)
    surv_model_true <- rfsrc(formula = surv_formula_true, data = trainingData, ntime = NULL)
    S_hat_true <- predict.rsf(surv_model_true, evalData, ref_times)
  } else if(event_estimator == "mboost"){
    require(mboost)
    surv_model_true <- mboost(formula = surv_formula_true, data = trainingData, family = CoxPH())
    S_hat_true <- predict.mboost(surv_model_true, evalData, ref_times)
  } else if(event_estimator == "glmnet"){
    require(glmnet)
    surv_model_true <- cv.glmnet(y = Surv(trainingData$Ttrue, trainingData$event_Ttrue), x = data.matrix(trainingData[,grep("(X.*[[:digit:]])",names(trainingData)), with = F]),family = "cox")
    S_hat_true <- predict.glmnet(model = surv_model_true, trainingData = trainingData, testData = evalData, times = ref_times, obs_time = F)
  } else if(event_estimator == "rangersf"){
    require(ranger)
    
    rf_params_trainingData <- rf.tuning(trainingData, surv_formula_true, length(grep("(X.*[[:digit:]])",names(trainingData))))
    
    surv_model_true <- ranger(formula = surv_formula_true, data = trainingData, num.trees = 1000, mtry = rf_params_trainingData$mtry, min.node.size = rf_params_trainingData$min.node.size, max.depth = rf_params_trainingData$max.depth, sample.fraction = rf_params_trainingData$sample.fraction, num.threads = 1)
    S_hat_true <- predict.rangersf(surv_model_true, evalData, ref_times)
    
  } else if(event_estimator == "xgboost"){
    require(xgboost)
    
    xgb_params_trainingData <- xgb.tuning(xgb.trainingData.data)
    
    surv_model_true <- xgboost(params = xgb_params_trainingData$pars_optim, data = xgb.trainingData.data, nrounds = xgb_params_trainingData$nround_optim, verbose = F)
    S_hat_true <- predict.xgb(model = surv_model_true, status = trainingData$event_Ttrue, observed.times = trainingData$Ttrue, training.data = xgb.trainingData.data, new.data = xgb.evalData.data, times = ref_times)$S_hat
  }
  
  #true survival function
  S_Weibull <- function(t, lambda, sigma = s){
    StX <- exp(-lambda*t**(1/sigma))
    return(StX)
  }
  
  m_t <- matrix(rep(t(ref_times), nrow(evalData)), nrow = nrow(evalData), byrow=T)
  m_lambda <- matrix(rep(evalData$lambda, length(ref_times)), ncol = length(ref_times), byrow = F)
  
  m_Ttrue <- matrix(rep(evalData$Ttrue, length(ref_times)), ncol = length(ref_times), byrow = F)
  ind <- m_Ttrue > m_t
  
  S <- S_Weibull(m_t, m_lambda, s)
  
  #plot(ref_times, colMeans(S), type = 'l')
  
  #mean squared error
  MSE <- function(n, S, S_hat_true){
    se <- (S-S_hat_true)**2
    mse <- (1/n)*colSums(se)
    return(mse)
  }
  
  mse <- MSE(nrow(evalData), S, S_hat_true)
  
  #plot(ref_times, mse, type = 'l')
  
  #expected uncensored brier score
  BS <- function(mse, n, S){
    cnst <- (1/n)*colSums(S*(1-S))
    bs <- mse + cnst
    return(bs)
  }
  
  bs_expected <- BS(mse, nrow(evalData), S)
  
  #brier score with step functions
  BS_step <- function(ind, S_hat_true, n){
    se <- (ind-S_hat_true)**2
    bs <- (1/n)*colSums(se)
    return(bs)
  }
  
  #bs_step <- BS_step(ind, S_hat_true, nrow(evalData))
  
  return(list("mse" = mse,"bs" = bs_expected))
}
