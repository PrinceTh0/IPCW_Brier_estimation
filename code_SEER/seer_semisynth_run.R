################################################################################
#                                                                              #
# SEER SEMISYNTH RUN                                                           #
#                                                                              #
# Generates the semi-synthetic data and submits jobs to batchtools             #
#                                                                              #
#                                                                              #
################################################################################

rm(list = ls())

require("batchtools")
require("data.table")
require("survival")

#data directory
data.dir <- "..."

#output directory
out.dir <- "..."

#display progress?
#options(batchtools.progress = FALSE)

#load preprocessed SEER data
load(file.path(data.dir,"breast04semis.RData"))
breast04semis <- as.data.table(breast04semis)

censoring_formula_km <- as.formula(Surv(time.discrete, 1-state) ~ 1)

#stratify by erstatus
erstatus1_data <- breast04semis[erstatus == 1]
erstatus2_data <- breast04semis[erstatus == 2]

#fit Kaplan-Meier censoring model (stratified by erstatus)
erstatus1_surv_long <- erstatus1_data[,c("time.discrete","state")]
setorder(erstatus1_surv_long,time.discrete, state)

erstatus1_surv_wide <- dcast(erstatus1_surv_long, fun.aggregate = length, formula = time.discrete ~ state, value.var = "state")
setnames(erstatus1_surv_wide,c("0","1"),c("censored","event"))

erstatus1_surv_wide[, atrisk := rev(cumsum(rev(erstatus1_surv_wide[,censored]))) + rev(cumsum(rev(erstatus1_surv_wide[,event])))]

erstatus1_surv_wide[, Skm := cumprod((1 - (erstatus1_surv_wide[,event]/erstatus1_surv_wide[,atrisk])))]
erstatus1_surv_wide[, Ckm := cumprod((1 - (erstatus1_surv_wide[,censored]/erstatus1_surv_wide[,atrisk])))]

erstatus1_surv <- rbind(data.table("time.discrete"=0,"censored"=0,"event"=0,"atrisk"=max(erstatus1_surv_wide[,atrisk]),"Skm"=1,"Ckm"=1),erstatus1_surv_wide)


erstatus2_surv_long <- erstatus2_data[,c("time.discrete","state")]
setorder(erstatus2_surv_long,time.discrete, state)

erstatus2_surv_wide <- dcast(erstatus2_surv_long, fun.aggregate = length, formula = time.discrete ~ state, value.var = "state")
setnames(erstatus2_surv_wide,c("0","1"),c("censored","event"))

erstatus2_surv_wide[, atrisk := rev(cumsum(rev(erstatus2_surv_wide[,censored]))) + rev(cumsum(rev(erstatus2_surv_wide[,event])))]

erstatus2_surv_wide[, Skm := cumprod((1 - (erstatus2_surv_wide[,event]/erstatus2_surv_wide[,atrisk])))]
erstatus2_surv_wide[, Ckm := cumprod((1 - (erstatus2_surv_wide[,censored]/erstatus2_surv_wide[,atrisk])))]

erstatus2_surv <- rbind(data.table("time.discrete"=0,"censored"=0,"event"=0,"atrisk"=max(erstatus2_surv_wide[,atrisk]),"Skm"=1,"Ckm"=1), erstatus2_surv_wide)

#sample new observations from estimated censoring distribution
set.seed(1234)

samples_unif_er1 <- runif(nrow(erstatus1_data[state == 1]))
samples_unif_er1 <- samples_unif_er1[order(samples_unif_er1, decreasing = TRUE)]
samples_unif_er2 <- runif(nrow(erstatus2_data[state == 1]))
samples_unif_er2 <- samples_unif_er2[order(samples_unif_er2, decreasing = TRUE)]

erstatus1_censoringt <- erstatus1_surv[apply(sapply(erstatus1_surv[,Ckm] , function(p){p < samples_unif_er1}),1,"which.max"), time.discrete]
erstatus2_censoringt <- erstatus1_surv[apply(sapply(erstatus2_surv[,Ckm] , function(p){p < samples_unif_er2}),1,"which.max"), time.discrete]

erstatus1_censoringt <- sample(erstatus1_censoringt)
erstatus2_censoringt <- sample(erstatus2_censoringt)

erstatus1_uncensored <- erstatus1_data[state==1]
erstatus1_uncensored <- erstatus1_uncensored[sample(nrow(erstatus1_uncensored))]
erstatus2_uncensored <- erstatus2_data[state==1]
erstatus2_uncensored <- erstatus2_uncensored[sample(nrow(erstatus2_uncensored))]

erstatus1_synthcensored <- copy(erstatus1_uncensored)
erstatus1_synthcensored[,time.discrete := pmin(erstatus1_censoringt,erstatus1_uncensored[,time.discrete])][,state := ifelse(erstatus1_censoringt < erstatus1_uncensored[,time.discrete],0,1)]

erstatus2_synthcensored <- copy(erstatus2_uncensored)
erstatus2_synthcensored[,time.discrete := pmin(erstatus2_censoringt,erstatus2_uncensored[,time.discrete])][,state := ifelse(erstatus2_censoringt < erstatus2_uncensored[,time.discrete],0,1)]

data <- rbind(erstatus1_synthcensored,erstatus2_synthcensored)

#training/test split
split <- 0.5
data_sample <- sample(nrow(data), nrow(data)*split)

training <- as.data.table(data[data_sample,])
test <- as.data.table(data[-data_sample,])

#order by observed times
setorder(data, time.discrete, -state)
setorder(training, time.discrete, -state)
setorder(test, time.discrete, -state)

#sequence of reference times
times <- seq(0,max(data$time.discrete))


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
  submitJobs(ids, res = list(walltime = 1*60*60, memory = 32000, ncpus = 1, partition = "batch"), reg = reg)
  
}


methods <- c("cox","xgboost","rangersf")

for(method in methods){
  bootstrap.slurm.run(method = method)
}
