
rf.tuning <- function(dataset, formula){
  
  require(ParBayesianOptimization)
  
  n_pred <- ncol(dataset)-2
  
  #set bounds for parameter search
  bounds <- list(mtry = c(1, n_pred),
                 max.depth = c(1, 10),
                 min.node.size = c(1, 10),
                 sample.fraction = c(0.1, 1)
  )
  
  #define scoring function to optimize
  scoring.function <- function(mtry, max.depth, min.node.size, sample.fraction){
    
    rf <- ranger(formula = formula, data = dataset, num.trees = 500, mtry = mtry, min.node.size = min.node.size, max.depth = max.depth, sample.fraction = sample.fraction, num.threads = 1)
    
    return(list(Score = -1*rf$prediction.error)
    )
  }
  
  #run Bayesian optimization
  bayes_tune_out <- bayesOpt(FUN = scoring.function, bounds = bounds, initPoints = length(bounds) + 2, iters.n = 8, plotProgress = FALSE, verbose = 0)
  
  pars_optim <- getBestPars(bayes_tune_out, N = 1)
  
  return(pars_optim)
  
}
