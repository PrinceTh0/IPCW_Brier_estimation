
rm(list = ls())

require("batchtools")
require("data.table")

#data directory
data.dir <- "..."

#output directory
out.dir <- "..."

#display progress?
#options(batchtools.progress = FALSE)

#load preprocessed SEER data
load(file.path(data.dir,"breast04.RData"))

#training/test split
set.seed(1234)
split <- 0.5
data <- as.data.table(breast04)
data_sample <- sample(nrow(data), nrow(data)*split)

training <- as.data.table(breast04[data_sample,])
test <- as.data.table(breast04[-data_sample,])

#order by observed times
setorder(data, time.discrete, -state)
setorder(training, time.discrete, -state)
setorder(test, time.discrete, -state)

times <- c(0,unique(data$time.discrete)[order(unique(data$time.discrete))])


bootstrap.slurm <- function(i, data, training, test, times, method, ipcw.weights, censoring.model, tuning){
  
  #dependencies
  source('seer_bootstrap.R')
  
  set.seed(1234*i)
  
  #catch errors in case of model convergence issues etc.
  sim_result <- tryCatch(
    
    expr = {
      bootstrap.bs(data=data, training=training, test=test, times=times, method=method, ipcw.weights=ipcw.weights, censoring.model=censoring.model, tuning=tuning)
    },
    
    error = function(e){
      message("\nAn error occurred:\n")
      message(paste0(e,"\n"))
      return(NULL)
    })
  
  if(is.null(sim_result)){
    next
  }
  
  return(sim_result)
}


##set parameters
ipcw.weights <- 0.2 #min value of 1/IPCW weights
n.samples    <- 10 #number of bootstrap samples
censoring.model <- "full" #marginal censoring model (KM estimator) or full censoring model
tuning <- FALSE #tuning of random forest & xgboost hyperparameters

#submit jobs to cluster
bootstrap.slurm.run <- function(method){
  
  reg <- makeRegistry(file.dir = paste0(out.dir,"/registry_",method,"_",1/ipcw.weights,"_",censoring.model,"_tuning",tuning), seed = 1)
  ids <- batchMap(fun = bootstrap.slurm, args = list(i=1:n.samples, method = method), more.args = list(data = data, training = training, test = test, times = times, ipcw.weights = ipcw.weights, censoring.model = censoring.model, tuning = tuning), reg = reg)
  submitJobs(ids, res = list(walltime = 12*60*60, memory = 64000, ncpus = 1, partition = "batch"), reg = reg)
  
}


methods <- c("cox","xgboost","rangersf")

for(method in methods){
  bootstrap.slurm.run(method = method)
}
