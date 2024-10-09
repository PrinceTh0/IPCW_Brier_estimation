
sim <- function(n_training, n_test, validation, n_cov_informative, n_cov_noninformative, s, s_c, rate, offset, admin_cens, ref_times, structure, event_estimator, cens_estimator, cens_distr, min_G_hat){

  ##packages
  require(data.table)
  require(evd)
  require(glmnet)
  require(MASS)
  require(MatrixModels)
  require(mboost)
  require(prodlim)
  require(randomForestSRC)
  require(ranger)
  require(survival)
  require(xgboost)
  
  ##dependencies
  #code for generating time-to-event dataset
  source('sim_weibull.R')
  #estimation/prediction functions
  source('sim_support_functions.R')
  #IPCW Brier score
  source('sim_ipcw_brier_score.R')
  #uncensored Brier score
  source('sim_uncensored_brier_score.R')
  #wse evaluation metric
  source('sim_wse.R')
  #xgboost tuning
  source('sim_xgboost_tuning.R')
  #random forest tuning
  source('sim_random_forest_tuning.R')
  
  num_threads <- 1
  
  n_data <- n_training+n_test
  n_validation <- nrow(validation)
  n_predictors <- n_cov_informative+n_cov_noninformative
  
  ##table with simulation parameter settings
  settings <- data.table("n sample" = n_data,
                         "n training" = n_training,
                         "n test" = n_test,
                         "n validation" = n_validation,
                         "informative covariates" = n_cov_informative,
                         "non-informative covariates" = n_cov_noninformative,
                         "s" = s,
                         "s censoring" = s_c,
                         "rate" = rate,
                         "offset" = offset,
                         "administrative censoring" = admin_cens,
                         "structure" = structure,
                         "event estimator" = event_estimator,
                         "censoring estimator" = cens_estimator,
                         "censoring model" = cens_distr)
  
    ##generate time-to-event dataset
    sample <- get.WeibullData(n_data, n_cov_informative, n_cov_noninformative, s, s_c, rate, offset, admin_cens, structure, cens_distr, seed = FALSE)
   
    #split into training & test data
    split <- n_training/(n_training+n_test)
    data_sample <- sample(nrow(sample$data), nrow(sample$data)*split)
    
    data <- as.data.table(sample$data)
    training <- as.data.table(sample$data[data_sample,])
    test <- as.data.table(sample$data[-data_sample,])
    R2 <- as.vector(sample$R2)
    R2c <- as.vector(sample$R2c)
    
    ##order validation set by true times
    setorder(validation, Ttrue)

    ##characteristics of generated time-to-event data
    sample_table <- data.table('n_data' = n_data,
                               'n_training' = n_training,
                               'n_test' = n_test,
                               'n_validation' = n_validation,
                               '%cens_data' = (n_data-sum(data$event))/n_data,
                               '%cens_training' = (n_training-sum(training$event))/n_training,
                               '%cens_test' = (n_test-sum(test$event))/n_test,
                               #'min_time' = data[, min(Ttilde)],
                               #'q1_time' = data[, quantile(Ttilde, 0.25)],
                               #'med_time' = data[, quantile(Ttilde, 0.5)],
                               #'q3_time' = data[, quantile(Ttilde, 0.75)],
                               #'max_time' = data[, max(Ttilde)],
                               'R2' = R2,
                               'R2c' = R2c)

    ##order by observed times
    setorder(data, Ttilde)
    setorder(training, Ttilde)
    setorder(test, Ttilde)
    
    ##prepare data for glmnet
    if(event_estimator == "glmnet" | cens_estimator == "glmnet"){
      trainingCens <- copy(training)
      trainingCens[, event := 1-event]
      testCens <- copy(test)
      testCens[, event := 1-event]
      dataCens <- copy(data)
      dataCens[, event := 1-event]
    }
    
    ##prepare data for xgboost
    if(event_estimator == "xgboost" | cens_estimator == "xgboost"){
      xgb.training <- training[, c("Ttilde","event",names(training)[grep("(X.*[[:digit:]])",names(training))]), with = F]
      xgb.training[, outcome := ifelse(event == 1, Ttilde, Ttilde * (-1))]
      xgb.training[, Ttilde := NULL][, event := NULL]
      xgb.training.dummy <- model.Matrix(object = outcome ~ ., data = xgb.training, sparse = T)[,-1]
      xgb.training.data <- xgb.DMatrix(xgb.training.dummy, label = xgb.training$outcome)
      xgb.training.data.cens <- xgb.DMatrix(xgb.training.dummy, label = xgb.training$outcome*(-1))
      
      xgb.test <- test[, c("Ttilde","event",names(test)[grep("(X.*[[:digit:]])",names(test))]), with = F]
      xgb.test[, outcome := ifelse(event == 1, Ttilde, Ttilde * (-1))]
      xgb.test[, Ttilde := NULL][, event := NULL]
      xgb.test.dummy <- model.Matrix(object = outcome ~ ., data = xgb.test, sparse = T)[,-1]
      xgb.test.data <- xgb.DMatrix(xgb.test.dummy, label = xgb.test$outcome)
      xgb.test.data.cens <- xgb.DMatrix(xgb.test.dummy, label = xgb.test$outcome*(-1))
      
      xgb.data <- data[, c("Ttilde","event",names(data)[grep("(X.*[[:digit:]])",names(data))]), with = F]
      xgb.data[, outcome := ifelse(event == 1, Ttilde, Ttilde * (-1))]
      xgb.data[, Ttilde := NULL][, event := NULL]
      xgb.data.dummy <- model.Matrix(object = outcome ~ ., data = xgb.data, sparse = T)[,-1]
      xgb.data.data <- xgb.DMatrix(xgb.data.dummy, label = xgb.data$outcome)
      xgb.data.data.cens <- xgb.DMatrix(xgb.data.dummy, label = xgb.data$outcome*(-1))
      
      #xgboost model tuning
      xgb_params_training <- xgb.tuning(xgb.training.data)
      
      if(cens_estimator != "km"){
        xgb_params_cens_training <- xgb.tuning(xgb.training.data.cens)
        xgb_params_cens_test <- xgb.tuning(xgb.test.data.cens)
        xgb_params_cens_data <- xgb.tuning(xgb.data.data.cens)
      }
      
    }
    
    #survival model with predictors
    surv_formula <- as.formula(paste("Surv(Ttilde, event)~", paste(names(data)[grep("(X.*[[:digit:]])",names(data))], collapse="+")))
    #censoring model with predictors
    cens_formula <- as.formula(paste("Surv(Ttilde, (1-event))~", paste(names(data)[grep("(X.*[[:digit:]])",names(data))], collapse="+")))
    #marginal censoring model without predictors
    cens_formula_marginal <- as.formula(Surv(Ttilde, 1-event) ~ 1)
    
    #random forest model tuning
    if(event_estimator == "rangersf" | cens_estimator == "rangersf"){
      
      rf_params_training <- rf.tuning(training, surv_formula, n_predictors)
      
      if(cens_estimator != "km"){
        rf_params_cens_training <- rf.tuning(training, cens_formula, n_predictors)
        rf_params_cens_test <- rf.tuning(test, cens_formula, n_predictors)
        rf_params_cens_data <- rf.tuning(data, cens_formula, n_predictors)
      }
    }
    
    ##fit survival model
    
    if(event_estimator == "cox"){
      require(survival)
      surv_model <- coxph(formula = surv_formula, data = training, x = TRUE, y = TRUE)
    } else if(event_estimator == "rsf"){
      require(randomForestSRC)
      surv_model <- rfsrc(formula = surv_formula, data = training, ntime = NULL)
    } else if(event_estimator == "mboost"){
      require(mboost)
      surv_model <- mboost(formula = surv_formula, data = training, family = CoxPH())
    } else if(event_estimator == "glmnet"){
      require(glmnet)
      surv_model <- cv.glmnet(y = Surv(training$Ttilde, training$event), x = data.matrix(training[,grep("(X.*[[:digit:]])",names(training)), with = F]),family = "cox")
    } else if(event_estimator == "rangersf"){
      require(ranger)
      surv_model <- ranger(formula = surv_formula, data = training, num.trees = 1000, mtry = rf_params_training$mtry, min.node.size = rf_params_training$min.node.size, max.depth = rf_params_training$max.depth, sample.fraction = rf_params_training$sample.fraction, num.threads = num_threads)
    } else if(event_estimator == "xgboost"){
      require(xgboost)
      surv_model <- xgboost(params = xgb_params_training$pars_optim, data = xgb.training.data, nrounds = xgb_params_training$nround_optim, verbose = FALSE)
    }
    
    ##estimate survival distribution on the test set
    if(event_estimator == "cox"){
      S_hat <- predict.cox(surv_model, test, ref_times)
    } else if(event_estimator == "rsf"){
      S_hat <- predict.rsf(surv_model, test, ref_times)
    } else if(event_estimator == "mboost"){
      S_hat <- predict.mboost(surv_model, test, ref_times)
    } else if(event_estimator == "glmnet"){
      S_hat <- predict.glmnet(surv_model, training, test, ref_times, T)
    } else if(event_estimator == "rangersf"){
      S_hat <- predict.rangersf(surv_model, test, ref_times)
    } else if(event_estimator == "xgboost"){
      S_hat <- predict.xgb(model = surv_model, status = training$event, observed.times = training$Ttilde, training.data = xgb.training.data, new.data = xgb.test.data, times = ref_times)$S_hat
    }
    
    ##calculate uncensored brier score & mse on validation data
    surv_formula_true <- as.formula(paste("Surv(Ttrue, event_Ttrue)~", paste(names(training)[grep("(X.*[[:digit:]])",names(training))], collapse="+")))
    reference_score <- uncensored_brier_score(validation, training, s, ref_times, surv_formula_true, event_estimator)

    ##fit censoring model on the training, test and combined dataset
    out_bs <- list()

    for(G_fit in c('training','test','all')){

      message(paste0("Fitting G on: ",G_fit," data"))
      
      fit_model_G <- function(cens_estimator, cens_formula, cens_formula_marginal, fitdata){
        if(cens_estimator == "km"){
          cens_model <- coxph(formula = cens_formula_marginal, data = fitdata, x = TRUE, y = TRUE)
        } else if(cens_estimator == "cox"){
          cens_model <- coxph(formula = cens_formula, data = fitdata, x = T, y = T)
        } else if(cens_estimator == "rsf"){
          cens_model <- rfsrc(formula = cens_formula, data = fitdata, ntime = NULL)
        } else if(cens_estimator == "mboost"){
          cens_model <- mboost(formula = cens_formula, data = fitdata, family = CoxPH())
        } else if(cens_estimator == "glmnet"){
          cens_model <- cv.glmnet(y = Surv(fitdata$Ttilde, 1-fitdata$event), x = data.matrix(fitdata[,grep("(X.*[[:digit:]])",names(fitdata)), with = FALSE]),family = "cox")
        } 
        return(cens_model)
      }
      
      
      if(cens_estimator == "xgboost"){
        if(G_fit == 'training'){
          cens_model <- xgboost(params = xgb_params_cens_training$pars_optim, data = xgb.training.data.cens, nrounds = xgb_params_cens_training$nround_optim, verbose = FALSE)
        } else if(G_fit == 'test'){
          cens_model <- xgboost(params = xgb_params_cens_test$pars_optim, data = xgb.test.data.cens, nrounds = xgb_params_cens_test$nround_optim, verbose = FALSE)
        } else if(G_fit == 'all'){
          cens_model <- xgboost(params = xgb_params_cens_data$pars_optim, data = xgb.data.data.cens, nrounds = xgb_params_cens_data$nround_optim, verbose = FALSE)
        }
      } else if(cens_estimator == "rangersf"){
        if(G_fit == 'training'){
          cens_model <- ranger(formula = cens_formula, data = training, num.trees = 1000, mtry = rf_params_cens_training$mtry, min.node.size = rf_params_cens_training$min.node.size, max.depth = rf_params_cens_training$max.depth, sample.fraction = rf_params_cens_training$sample.fraction, num.threads = num_threads)
        } else if(G_fit == 'test'){
          cens_model <- ranger(formula = cens_formula, data = test, num.trees = 1000, mtry = rf_params_cens_test$mtry, min.node.size = rf_params_cens_test$min.node.size, max.depth = rf_params_cens_test$max.depth, sample.fraction = rf_params_cens_test$sample.fraction, num.threads = num_threads)
        } else if(G_fit == 'all'){
          cens_model <- ranger(formula = cens_formula, data = data, num.trees = 1000, mtry = rf_params_cens_data$mtry, min.node.size = rf_params_cens_data$min.node.size, max.depth = rf_params_cens_data$max.depth, sample.fraction = rf_params_cens_data$sample.fraction, num.threads = num_threads)
        }
      } else{
        if(G_fit == 'training'){
          cens_model <- fit_model_G(cens_estimator,cens_formula,cens_formula_marginal,training)
        } else if(G_fit == 'test'){
          cens_model <- fit_model_G(cens_estimator,cens_formula,cens_formula_marginal,test)
        } else if(G_fit == 'all'){
          cens_model <- fit_model_G(cens_estimator,cens_formula,cens_formula_marginal,data)
        }
      }
      
      test_Tmin <- c(0,test$Ttilde[-length(test$Ttilde)])
      
      ##estimate censoring distribution (evaluted at t and Ttilde-) on the test set
      if(cens_estimator == "km"){
        G_hat_t <- predict.cox(cens_model, test, ref_times)
        G_hat_event <- predict.cox(cens_model, test, test_Tmin)
      } else if(cens_estimator == "cox"){
        G_hat_t <- predict.cox(cens_model, test, ref_times)
        G_hat_event <- predict.cox(cens_model, test, test_Tmin)
      } else if(cens_estimator == "rsf"){
        G_hat_t <- predict.rsf(cens_model, test, ref_times)
        G_hat_event <- predict.rsf(cens_model, test, test_Tmin)
      } else if(cens_estimator == "mboost"){
        G_hat_t <- predict.mboost(cens_model, test, ref_times)
        G_hat_event <- predict.mboost(cens_model, test, test_Tmin)
      } else if(cens_estimator == "glmnet"){
        if(G_fit == "training"){
          G_hat_t <- predict.glmnet(cens_model, trainingCens, testCens, ref_times, obs_time = TRUE)
          G_hat_event <- predict.glmnet(cens_model, trainingCens, testCens, test_Tmin, obs_time = TRUE)
        } else if(G_fit == "test"){
          G_hat_t <- predict.glmnet(cens_model, testCens, testCens, ref_times, obs_time = TRUE)
          G_hat_event <- predict.glmnet(cens_model, testCens, testCens, test_Tmin, obs_time = TRUE)
        } else if(G_fit == "all"){
          G_hat_t <- predict.glmnet(cens_model, dataCens, testCens, ref_times, obs_time = TRUE)
          G_hat_event <- predict.glmnet(cens_model, dataCens, testCens, test_Tmin, obs_time = TRUE)
        }
      } else if(cens_estimator == "rangersf"){
        G_hat_t <- predict.rangersf(cens_model, test, ref_times)
        G_hat_event <- predict.rangersf(cens_model, test, test_Tmin)
      } else if(cens_estimator == "xgboost"){
        if(G_fit == "training"){
          G_hat_t <- predict.xgb(model = cens_model, status = 1-training$event, observed.times = training$Ttilde, training.data = xgb.training.data.cens, new.data = xgb.test.data.cens, times = ref_times)$S_hat
          G_hat_event <- predict.xgb(model = cens_model, status = 1-training$event, observed.times = training$Ttilde, training.data = xgb.training.data.cens, new.data = xgb.test.data.cens, times = test_Tmin)$S_hat
        } else if(G_fit == "test"){
          G_hat_t <- predict.xgb(model = cens_model, status = 1-test$event, observed.times = test$Ttilde, training.data = xgb.test.data.cens, new.data = xgb.test.data.cens, times = ref_times)$S_hat
          G_hat_event <- predict.xgb(model = cens_model, status = 1-test$event, observed.times = test$Ttilde, training.data = xgb.test.data.cens, new.data = xgb.test.data.cens, times = test_Tmin)$S_hat
        } else if(G_fit == "all"){
          G_hat_t <- predict.xgb(model = cens_model, status = 1-data$event, observed.times = data$Ttilde, training.data = xgb.data.data.cens, new.data = xgb.test.data.cens, times = ref_times)$S_hat
          G_hat_event <- predict.xgb(model = cens_model, status = 1-data$event, observed.times = data$Ttilde, training.data = xgb.data.data.cens, new.data = xgb.test.data.cens, times = test_Tmin)$S_hat
        }
      }
      
      
      ##calculate IPCW brier score
      obs_times <- test$Ttilde #observed event/censoring times
      delta <- test$event #event indicator
      
      m_obs_times <- matrix(rep(obs_times, length(ref_times)),ncol=length(ref_times),byrow=FALSE)
      m_delta <- matrix(rep(delta, length(ref_times)),ncol=length(ref_times),byrow=FALSE)
      m_ref_times <- matrix(rep(t(ref_times), length(obs_times)), nrow = length(obs_times), byrow=TRUE)
      
      #indiv_times <- as.matrix(merge.data.table(data.table("obs_times"=obs_times)[, idx1 := 1:length(obs_times)],data.table("obs_times"=unique(obs_times))[, idx2 := 1:length(unique(obs_times))], by = "obs_times", all.x = T)[, obs_times := NULL])
      #G_hat_event_t_min_1 <- c(1, G_hat_event[indiv_times][-nrow(indiv_times)])
      
      #G_hat_event_t_min_1 <- c(1, diag(G_hat_event)[-length(diag(G_hat_event))])
      #m_G_hat_event_t_min_1 <- matrix(rep(G_hat_event_t_min_1, length(ref_times)),ncol=length(ref_times),byrow=FALSE)
      
      m_G_hat_event_t_min_1 <- matrix(rep(diag(G_hat_event), length(ref_times)),ncol=length(ref_times),byrow=FALSE)
      
      #set minimum for G_hat (equivalent to setting maximum weight)
      G_hat_t <- ifelse(G_hat_t < min_G_hat, min_G_hat, G_hat_t)
      m_G_hat_event_t_min_1 <- ifelse(m_G_hat_event_t_min_1 < min_G_hat, min_G_hat, m_G_hat_event_t_min_1)
      
      bs <- brier_score(m_delta, m_obs_times, m_ref_times, m_G_hat_event_t_min_1, G_hat_t, S_hat, n_test, max_weight = 100, brier_general = TRUE)
      
      if(G_fit == 'training'){
        out_bs$training <- data.table(matrix(bs,dimnames=list(NULL,c("bs"))))
      } else if(G_fit == 'test'){
        out_bs$test <- data.table(matrix(bs,dimnames=list(NULL,c("bs"))))
      } else if(G_fit == 'all'){
        out_bs$all <- data.table(matrix(bs,dimnames=list(NULL,c("bs"))))
      }
    }
    
    ##calculate weighted squared error of IPCW and uncensored Brier score
    times <- rbind(training[, dataset := "train"][, c("Ttrue","Ttilde", "dataset")], 
                   test[, dataset := "test"][, c("Ttrue","Ttilde", "dataset")])
    
    out_obs_times <- times[, c("Ttilde","dataset")][order(-dataset,Ttilde)]
    out_true_times <- times[, c("Ttrue", "dataset")][order(-dataset,Ttrue)]

    out_uncensored_brier_score <- data.table(matrix(reference_score$bs, ncol = 1, dimnames = list(NULL, c("uncensored"))))
    out_mse <- data.table(matrix(reference_score$mse, ncol = 1, dimnames = list(NULL, c("mse"))))
    
    wse <- list()
    
    wse$train <- wse.calc(bs.in = out_bs$training$bs, uncensored.in = out_uncensored_brier_score$uncensored, times.in = validation$Ttrue, ref.times.in = ref_times)
    wse$test <- wse.calc(bs.in = out_bs$test$bs, uncensored.in = out_uncensored_brier_score$uncensored, times.in = validation$Ttrue, ref.times.in = ref_times)
    wse$all <- wse.calc(bs.in = out_bs$all$bs, uncensored.in = out_uncensored_brier_score$uncensored, times.in = validation$Ttrue, ref.times.in = ref_times)
    
    
    return(list("time" = ref_times,
                "obs.times" = out_obs_times,
                "true.times" = out_true_times,
                "sample.parameters" = sample_table,
                "ipcw.brier.score" = out_bs,
                "uncensored.brier.score" = out_uncensored_brier_score,
                "mse" = out_mse,
                "wse.score" = wse,
                "n.obs" = n_data,
                "settings" = settings))
}

#res <- sim(n_training=10000, n_test=10000, validation=validation, n_cov_informative=10, n_cov_noninformative=0, s=0.578, s_c=0.579, rate=0.864, offset=-0.002, admin_cens=FALSE, ref_times = seq(0,9.99,0.01), structure = "compoundSymmetry", event_estimator = "cox", cens_estimator = "cox", cens_distr = "Weibull", min_G_hat = 0.2)
