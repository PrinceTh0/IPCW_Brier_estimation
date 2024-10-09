xgb.tuning <- function(xgb.data.matrix){
  
  require(ParBayesianOptimization)
 
  folds <- list(
    fold1 = as.integer(seq(1,nrow(xgb.data.matrix),by = 5)),
    fold2 = as.integer(seq(2,nrow(xgb.data.matrix),by = 5)),
    fold3 = as.integer(seq(3,nrow(xgb.data.matrix),by = 5)),
    fold4 = as.integer(seq(4,nrow(xgb.data.matrix),by = 5)),
    fold5 = as.integer(seq(5,nrow(xgb.data.matrix),by = 5))
  )
  
  #set bounds for parameter search
  bounds <- list(eta = c(0.01, 0.3),
                 max_depth = c(1L, 10L),
                 min_child_weight = c(1, 100),
                 subsample = c(0.1, 1),
                 colsample_bytree = c(0.1, 1),
                 gamma = c(0, 10),
                 lambda = c(1, 10),
                 alpha = c(1, 10)
  )
  
  #define scoring function to optimize
  scoring.function <- function(eta, max_depth, min_child_weight, subsample, colsample_bytree, gamma, lambda, alpha) {
    
    pars <- list(
      eta = eta,
      max_depth = max_depth,
      min_child_weight = min_child_weight,
      subsample = subsample,
      colsample_bytree = colsample_bytree,
      gamma = gamma,
      lambda = lambda,
      alpha = alpha,
      
      booster = "gbtree",
      objective = "survival:cox",
      eval_metric = "cox-nloglik",
      nthread = 1
    )
    
    xgbcv1 <- xgb.cv(
      params = pars,
      data = xgb.data.matrix,
      nrounds = 1000,
      folds = folds,
      early_stopping_rounds = 5,
      maximize = FALSE,
      verbose = 0
    )
    
    return(list(Score = -1*min(xgbcv1$evaluation_log$test_cox_nloglik_mean),
                nrounds = xgbcv1$best_iteration
                )
           )
  }
  
  #run Bayesian optimization
  bayes_tune_out <- bayesOpt(FUN = scoring.function, bounds = bounds, initPoints = length(bounds) + 2, iters.n = 8, plotProgress = FALSE, verbose = 0)
  
  pars_optim <- append(list(booster = "gbtree",
                          objective = "survival:cox",
                          eval_metric = "cox-nloglik",
                          nthread = 1),
                 getBestPars(bayes_tune_out, N = 1))

  #find best number of rounds
  xgbcv2 <- xgb.cv(params = pars_optim,
                   data = xgb.data.matrix,
                   nrounds = 1000,
                   folds = folds,
                   prediction = TRUE,
                   early_stopping_rounds = 5,
                   maximize = FALSE,
                   verbose = 0
                   )
  
  nround_optim = xgbcv2$best_iteration
  
  return(list("pars_optim" = pars_optim, "nround_optim" = nround_optim))
  
}
