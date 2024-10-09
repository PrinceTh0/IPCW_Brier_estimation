bootstrap.bs <- function(data, training, test, times, method, ipcw.weights, censoring.model, tuning){
  
  source("seer_support_functions.R")
  source("seer_ipcw_brier_score.R")
  source("seer_random_forest_tuning.R")
  source("seer_xgboost_tuning.R")
  
  require(data.table)
  require(MatrixModels)
  require(mboost)
  require(pec)
  require(ranger)
  require(survival)
  require(xgboost)
  
  #set.seed(1234)
  
  num_threads <- 1
  
  b.sample.data <- sample(nrow(data), replace = T)
  b.sample.training <- sample(nrow(training), replace = T)
  b.sample.test <- sample(nrow(test), replace = T)
  data.b <- data[b.sample.data]
  training.b <- training[b.sample.training]
  test.b <- test[b.sample.test]
  
  ##order by observed times
  setorder(data.b, time.discrete, -state)
  setorder(training.b, time.discrete, -state)
  setorder(test.b, time.discrete, -state)
  
  ##data prep for xgboost
  if(method == "xgboost"){
    
    #change data format
    xgb.data <- copy(data.b)
    xgb.data[, outcome := ifelse(state == 1, time.discrete, time.discrete * (-1))]
    xgb.data[, time.discrete := NULL][, state := NULL]
    xgb.data.dummy <- model.Matrix(object = outcome ~ ., data = xgb.data, sparse = T)[,-1]
    xgb.data.data <- xgb.DMatrix(xgb.data.dummy, label = xgb.data$outcome)
    xgb.data.data.cens <- xgb.DMatrix(xgb.data.dummy, label = xgb.data$outcome*(-1))
    
    xgb.training <- copy(training.b)
    xgb.training[, outcome := ifelse(state == 1, time.discrete, time.discrete * (-1))]
    xgb.training[, time.discrete := NULL][, state := NULL]
    xgb.training.dummy <- model.Matrix(object = outcome ~ ., data = xgb.training, sparse = T)[,-1]
    xgb.training.data <- xgb.DMatrix(xgb.training.dummy, label = xgb.training$outcome)
    xgb.training.data.cens <- xgb.DMatrix(xgb.training.dummy, label = xgb.training$outcome*(-1))
    
    xgb.test <- copy(test.b)
    xgb.test[, outcome := ifelse(state == 1, time.discrete, time.discrete * (-1))]
    xgb.test[, time.discrete := NULL][, state := NULL]
    xgb.test.dummy <- model.Matrix(object = outcome ~ ., data = xgb.test, sparse = T)[,-1]
    xgb.test.data <- xgb.DMatrix(xgb.test.dummy, label = xgb.test$outcome)
    xgb.test.data.cens <- xgb.DMatrix(xgb.test.dummy, label = xgb.test$outcome*(-1))
    
    #model tuning
    if(tuning == TRUE){
      
      xgb_params_training <- xgb.tuning(xgb.training.data)
      
      if(censoring.model == "full"){
        xgb_params_cens_data <- xgb.tuning(xgb.data.data.cens)
        xgb_params_cens_training <- xgb.tuning(xgb.training.data.cens)
        xgb_params_cens_test <- xgb.tuning(xgb.test.data.cens)
      }
      
    }
  }
  
  surv_formula <- as.formula("Surv(time.discrete, state) ~ .")
  
  cens_formula <- as.formula("Surv(time.discrete, (1-state)) ~ .")
  cens_formula_marginal <- as.formula(Surv(time.discrete, (1-state)) ~ 1)
  
  if(tuning == TRUE & method == "rangersf"){
    
    rf_params_training <- rf.tuning(training.b, surv_formula)
    
    if(censoring.model == "full"){
      rf_params_cens_data <- rf.tuning(data.b, cens_formula)
      rf_params_cens_training <- rf.tuning(training.b, cens_formula)
      rf_params_cens_test <- rf.tuning(test.b, cens_formula)
    }
  }
  
  
  ##fit survival model
  
  if(method == "cox"){
    surv_model <- coxph(formula = surv_formula, data = training.b, x = T, y = T)
  }else if(method == "rangersf"){
    if(tuning == TRUE){
      surv_model <- ranger(formula = surv_formula, data = training.b, mtry = rf_params_training$mtry, min.node.size = rf_params_training$min.node.size, max.depth = rf_params_training$max.depth, sample.fraction = rf_params_training$sample.fraction, num.threads = num_threads)
    }else if(tuning == FALSE){
      surv_model <- ranger(formula = surv_formula, data = training.b, num.threads = num_threads)
    }
  }else if(method == "xgboost"){
    if(tuning == TRUE){
      surv_model <- xgboost(params = xgb_params_training$pars_optim, data = xgb.training.data, nrounds = xgb_params_training$nround_optim, verbose = F)
    }else if(tuning == FALSE){
      surv_model <- xgboost(params = list(booster = "gbtree", objective = "survival:cox", eval_metric = "cox-nloglik", nthread = 1), data = xgb.training.data, nrounds = 500, verbose = F)
    }
  }
  
  ##estimate survival distribution on the test set
  
  if(method == "cox"){
    S_hat <- predictSurvProb(surv_model, test.b, times)
  }else if(method == "rangersf"){
    S_hat <- predict.rangersf(surv_model, test.b, times)
  }else if(method == "xgboost"){
    S_hat <- predict.xgb(model = surv_model, status = training.b$state, observed.times = training.b$time.discrete, training.data = xgb.training.data, new.data = xgb.test.data, times = times)$S_hat
  }
  
  ##fit censoring model
  
  if(censoring.model == "marginal"){
    
    cens_model_data <- coxph(formula = cens_formula_marginal, data = data.b, x = T, y = T)
    cens_model_training <- coxph(formula = cens_formula_marginal, data = training.b, x = T, y = T)
    cens_model_test <- coxph(formula = cens_formula_marginal, data = test.b, x = T, y = T)
    
    G_hat_t_data <- matrix(rep(predict.km(cens_model_data, times), nrow(data.b)), ncol = length(times), byrow = TRUE)
    G_hat_t_training <- matrix(rep(predict.km(cens_model_training, times), nrow(training.b)), ncol = length(times), byrow = TRUE)
    G_hat_t_test <- matrix(rep(predict.km(cens_model_test, times), nrow(test.b)), ncol = length(times), byrow = TRUE)
    
    
  }else if(censoring.model == "full"){
    
    if(method == "cox"){
      cens_model_data <- coxph(formula = cens_formula, data = data.b, x = T, y = T)
      cens_model_training <- coxph(formula = cens_formula, data = training.b, x = T, y = T)
      cens_model_test <- coxph(formula = cens_formula, data = test.b, x = T, y = T)
    }else if(method == "rangersf"){
      if(tuning == TRUE){
        cens_model_data <- ranger(formula = cens_formula, data = data.b, num.trees = 500, mtry = rf_params_cens_data$mtry, min.node.size = rf_params_cens_data$min.node.size, max.depth = rf_params_cens_data$max.depth, sample.fraction = rf_params_cens_data$sample.fraction, num.threads = num_threads)
        cens_model_training <- ranger(formula = cens_formula, data = training.b, num.trees = 500, mtry = rf_params_cens_training$mtry, min.node.size = rf_params_cens_training$min.node.size, max.depth = rf_params_cens_training$max.depth, sample.fraction = rf_params_cens_training$sample.fraction, num.threads = num_threads)
        cens_model_test <- ranger(formula = cens_formula, data = test.b, num.trees = 500, mtry = rf_params_cens_test$mtry, min.node.size = rf_params_cens_test$min.node.size, max.depth = rf_params_cens_test$max.depth, sample.fraction = rf_params_cens_test$sample.fraction, num.threads = num_threads)
      }else if(tuning == FALSE){
        cens_model_data <- ranger(formula = cens_formula, data = data.b, num.threads = num_threads)
        cens_model_training <- ranger(formula = cens_formula, data = training.b, num.threads = num_threads)
        cens_model_test <- ranger(formula = cens_formula, data = test.b, num.threads = num_threads)
      }
    }else if(method == "xgboost"){
      if(tuning == TRUE){
        cens_model_data <- xgboost(params = xgb_params_cens_data$pars_optim, data = xgb.data.data.cens, nrounds = xgb_params_cens_data$nround_optim, verbose = F)
        cens_model_training <- xgboost(params = xgb_params_cens_training$pars_optim, data = xgb.training.data.cens, nrounds = xgb_params_cens_training$nround_optim, verbose = F)
        cens_model_test <- xgboost(params = xgb_params_cens_test$pars_optim, data = xgb.test.data.cens, nrounds = xgb_params_cens_test$nround_optim, verbose = F)
      }else if(tuning == FALSE){
        cens_model_data <- xgboost(params = list(booster = "gbtree", objective = "survival:cox", eval_metric = "cox-nloglik", nthread = 1), data = xgb.data.data.cens, nrounds = 500, verbose = F)
        cens_model_training <- xgboost(params = list(booster = "gbtree", objective = "survival:cox", eval_metric = "cox-nloglik", nthread = 1), data = xgb.training.data.cens, nrounds = 500, verbose = F)
        cens_model_test <- xgboost(params = list(booster = "gbtree", objective = "survival:cox", eval_metric = "cox-nloglik", nthread = 1), data = xgb.test.data.cens, nrounds = 500, verbose = F)
      }
    }
    
    ##estimate censoring distribution on test set
    if(method == "cox"){
      G_hat_t_data <- predictSurvProb(cens_model_data, test.b, times)
      G_hat_t_training <- predictSurvProb(cens_model_training, test.b, times)
      G_hat_t_test <- predictSurvProb(cens_model_test, test.b, times)
    }else if(method == "rangersf"){
      G_hat_t_data <- predict.rangersf(cens_model_data, test.b, times)
      G_hat_t_training <- predict.rangersf(cens_model_training, test.b, times)
      G_hat_t_test <- predict.rangersf(cens_model_test, test.b, times)
    }else if(method == "xgboost"){
      G_hat_t_data <- predict.xgb(model = cens_model_data, status = 1-data.b$state, observed.times = data.b$time.discrete, training.data = xgb.data.data.cens, new.data = xgb.test.data.cens, times = times)$S_hat
      G_hat_t_training <- predict.xgb(model = cens_model_training, status = 1-training.b$state, observed.times = training.b$time.discrete, training.data = xgb.training.data.cens, new.data = xgb.test.data.cens, times = times)$S_hat
      G_hat_t_test <- predict.xgb(model = cens_model_test, status = 1-test.b$state, observed.times = test.b$time.discrete, training.data = xgb.test.data.cens, new.data = xgb.test.data.cens, times = times)$S_hat
    }
  }
  
  ##calculate IPCW brier score
  obs_times <- test.b$time.discrete
  delta <- test.b$state
  
  m_obs_times <- matrix(rep(obs_times, length(times)),ncol=length(times),byrow=F)
  m_delta <- matrix(rep(delta, length(times)),ncol=length(times),byrow=F)
  m_ref_times <- matrix(rep(t(times), length(obs_times)), nrow = length(obs_times), byrow=T)
  
  indiv_times <- cbind(1:nrow(test), obs_times)
  G_hat_event_t_min_1_data <- G_hat_t_data[indiv_times]
  G_hat_event_t_min_1_training <- G_hat_t_training[indiv_times]
  G_hat_event_t_min_1_test <- G_hat_t_test[indiv_times]
  
  m_G_hat_event_t_min_1_data <- matrix(rep(G_hat_event_t_min_1_data, length(times)),ncol=length(times),byrow=F)
  m_G_hat_event_t_min_1_training <- matrix(rep(G_hat_event_t_min_1_training, length(times)),ncol=length(times),byrow=F)
  m_G_hat_event_t_min_1_test <- matrix(rep(G_hat_event_t_min_1_test, length(times)),ncol=length(times),byrow=F)
  
  #set minimum for G_hat (equivalent to setting maximum weight)
  G_hat_t_data <- ifelse(G_hat_t_data < ipcw.weights, ipcw.weights, G_hat_t_data)
  G_hat_t_training <- ifelse(G_hat_t_training < ipcw.weights, ipcw.weights, G_hat_t_training)
  G_hat_t_test <- ifelse(G_hat_t_test < ipcw.weights, ipcw.weights, G_hat_t_test)
  
  m_G_hat_event_t_min_1_data <- ifelse(m_G_hat_event_t_min_1_data < ipcw.weights, ipcw.weights, m_G_hat_event_t_min_1_data)
  m_G_hat_event_t_min_1_training <- ifelse(m_G_hat_event_t_min_1_training < ipcw.weights, ipcw.weights, m_G_hat_event_t_min_1_training)
  m_G_hat_event_t_min_1_test <- ifelse(m_G_hat_event_t_min_1_test < ipcw.weights, ipcw.weights, m_G_hat_event_t_min_1_test)
  
  bs.data <- brier_score(m_delta, m_obs_times, m_ref_times, m_G_hat_event_t_min_1_data, G_hat_t_data, S_hat, n_test, max_weight = 100, brier_general = T)
  bs.training <- brier_score(m_delta, m_obs_times, m_ref_times, m_G_hat_event_t_min_1_training, G_hat_t_training, S_hat, n_test, max_weight = 100, brier_general = T)
  bs.test <- brier_score(m_delta, m_obs_times, m_ref_times, m_G_hat_event_t_min_1_test, G_hat_t_test, S_hat, n_test, max_weight = 100, brier_general = T)  
  
  return(list("method"=method,
              "weights"=ipcw.weights,
              "censoring.model"=censoring.model,
              "bs.data"=bs.data,
              "bs.training"=bs.training,
              "bs.test"=bs.test,
              "samples.data"=b.sample.data,
              "samples.training"=b.sample.training,
              "samples.test"=b.sample.test))
  
}

#res <- bootstrap.bs(data = data, training = training, test = test, times = times, method = "rangersf", ipcw.weights = 0.2, censoring.model = "full", tuning = FALSE)