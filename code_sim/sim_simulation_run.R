################################################################################
# SIMULATION RUN                                                               #
#                                                                              #
# Runs multiple simulations in parallel using SLURM commands.                  # 
#                                                                              #
# default scenario:                                                            #
# n_training = 500                                                             #
# n_test = 500                                                                 #
# n_validation = 2 * 10**5                                                     #
# n_cov_informative = 10                                                       #
# n_cov_nninformative = 0                                                      #
# s = 0.578 -> R2 = 0.5                                                        #
# s_c = 0.579 -> R2 = 0.5                                                      #
# rate = 0.864 -> censoring rate = 0.5 (if cens_distr = exponential)           #
# offset = -0.002 -> censoring rate = 0.5 (if cens_distr = Weibull)            #
# admin_cens = FALSE                                                           #
# structure = 'compoundSymmetry'                                               #
# cens_distr = exponential                                                     #
#                                                                              #
################################################################################

rm(list = ls())

require("batchtools")

#display progress?
#options(batchtools.progress = FALSE)

simulation.slurm <- function(i,n_training,n_test,validation,n_cov_informative,n_cov_noninformative,s,s_c,rate,offset,admin_cens,ref_times,structure,event_estimator,cens_estimator,cens_distr,min_G_hat){
  
  #dependencies
  source('sim_simulation.R')
  
  ref_times  <- seq(0, 9.99, 0.01)
  set.seed(1234*i)
  
  #catch errors in case of model convergence issues etc.
  sim_result <- tryCatch(
    
    expr = {
      sim(n_training=n_training,n_test=n_test,validation=validation,n_cov_informative=n_cov_informative,n_cov_noninformative=n_cov_noninformative,s=s,s_c=s_c,rate=rate,offset=offset,admin_cens=admin_cens,ref_times=ref_times,structure=structure,event_estimator=event_estimator,cens_estimator=cens_estimator,cens_distr=cens_distr,min_G_hat=min_G_hat)
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

##set global parameters
admin_cens.set           <- FALSE
structure.set            <- 'compoundSymmetry' # covariance structure of predictors
iter.set                 <- 1000 # number of simulation runs
min_G_hat.set            <- 0.2 # min value of 1/IPCW weights



################################################################################################################################
################################################################################################################################
### objective 1: assessing the effect of the overall sample size, relative size of training and test set, and the censoring rate
################################################################################################################################
################################################################################################################################

#submit jobs to cluster
simulation.slurm.run <- function(parameter,n_training,n_test,rate,offset){
  
  if(parameter == "n.per.ds"){
    poi <- "n.per.ds"
    poi.value <- n_training
  } else{
    parameters <- as.list(environment())
    poi <- names(parameters[which(names(parameters) == parameter)])
    poi.value <- get(poi)
  }
  
  validation <- get.WeibullData(n=n_validation.set, n_cov_informative=n_cov_informative.set, n_cov_noninformative=n_cov_noninformative.set, s=s.set, s_c=s_c.set, rate=rate, offset=offset, admin_cens=admin_cens.set, structure=structure.set, cens_distr=cens_distr.set, seed = T)$data
  
  reg <- makeRegistry(file.dir = paste0(out.dir,"/registry_",poi,"_",poi.value), seed = 1)
  ids <- batchMap(fun = simulation.slurm, args = list(i=1:iter.set, n_training=n_training, n_test=n_test, n_cov_informative=n_cov_informative.set, n_cov_noninformative=n_cov_noninformative.set, s=s.set, s_c=s_c.set, rate=rate, offset=offset, event_estimator=event_estimator.set, cens_estimator=cens_estimator.set, cens_distr=cens_distr.set, min_G_hat=min_G_hat.set), 
                  more.args= list(validation=validation,admin_cens=admin_cens.set,structure=structure.set), reg = reg)
  submitJobs(ids, res = list(walltime = 01*60*60, memory = 12000, ncpus = 1, partition = "batch"), reg = reg)
  
}

######################
#1a: exponential model
######################

#output directory
out.dir <- "..."

#settings:
n_training.set           <- 500
n_test.set               <- 500
n_validation.set         <- 2*10**5
n_cov_informative.set    <- 10
n_cov_noninformative.set <- 0
s.set                    <- 0.578
s_c.set                  <- 0.579
rate.set                 <- 0.864
offset.set               <- -0.002

cens_distr.set <- 'exponential'
event_estimator.set <- 'cox'
cens_estimator.set  <- 'km'

source("sim_weibull.R")


n.per.ds.values <- c(25,50,100,250,500)

for(n.per.ds in n.per.ds.values){
  simulation.slurm.run(parameter="n.per.ds", n_training=n.per.ds, n_test=n.per.ds, rate=rate.set, offset = offset.set)
}

n_training.values <- c(25,50,100,250,500)

for(n_training in n_training.values){
  simulation.slurm.run(parameter="n_training", n_training=n_training, n_test=n_test.set, rate=rate.set, offset = offset.set)
}

n_test.values <- c(25,50,100,250,500)

for(n_test in n_test.values){
  simulation.slurm.run(parameter="n_test", n_training=n_training.set, n_test=n_test, rate=rate.set, offset = offset.set)
}

rate.values <- c(0.219,0.864,3.013) #exponential censoring model

for(rate in rate.values){
  simulation.slurm.run(parameter="rate", n_training=n_training.set, n_test=n_test.set, rate=rate, offset = offset.set)
}


##################
#1b: Weibull model
##################

#output directory
out.dir <- "..."

#settings:
n_training.set           <- 500
n_test.set               <- 500
n_validation.set         <- 2*10**5
n_cov_informative.set    <- 10
n_cov_noninformative.set <- 0
s.set                    <- 0.578
s_c.set                  <- 0.579
rate.set                 <- 0.864
offset.set               <- -0.002

cens_distr.set <- 'Weibull'
event_estimator.set <- 'cox'
cens_estimator.set  <- 'cox'

source("sim_weibull.R")


n.per.ds.values <- c(25,50,100,250,500)

for(n.per.ds in n.per.ds.values){
  simulation.slurm.run(parameter="n.per.ds", n_training=n.per.ds, n_test=n.per.ds, rate=rate.set, offset = offset.set)
}

n_training.values <- c(25,50,100,250,500)

for(n_training in n_training.values){
  simulation.slurm.run(parameter="n_training", n_training=n_training, n_test=n_test.set, rate=rate.set, offset = offset.set)
}

n_test.values <- c(25,50,100,250,500)

for(n_test in n_test.values){
  simulation.slurm.run(parameter="n_test", n_training=n_training.set, n_test=n_test, rate=rate.set, offset = offset.set)
}

offset.values <- c(-0.803,-0.002,0.8) #Weibull censoring model

for(offset in offset.values){
  simulation.slurm.run(parameter="offset", n_training=n_training.set, n_test=n_test.set, rate=rate.set, offset = offset)
}



################################################################################
################################################################################
### objective 2: assessing the effect of model misspecification
################################################################################
################################################################################

#submit jobs to cluster
simulation.slurm.run.1 <- function(n.per.ds,cens_estimator){
  
  validation <- get.WeibullData(n=n_validation.set, n_cov_informative=n_cov_informative.set, n_cov_noninformative=n_cov_noninformative.set, s=s.set, s_c=s_c.set, rate=rate.set, offset=offset.set, admin_cens=admin_cens.set, structure=structure.set, cens_distr=cens_distr.set, seed = T)$data
  
  reg <- makeRegistry(file.dir = paste0(out.dir,"/registry_","n.per.ds","_",n.per.ds,"_","cens_estimator","_",cens_estimator), seed = 1)
  ids <- batchMap(fun = simulation.slurm, args = list(i=1:iter.set, n_training=n.per.ds, n_test=n.per.ds, n_cov_informative=n_cov_informative.set, n_cov_noninformative=n_cov_noninformative.set, s=s.set, s_c=s_c.set, rate=rate.set, offset=offset.set, event_estimator=event_estimator.set, cens_estimator=cens_estimator, cens_distr=cens_distr.set, min_G_hat=min_G_hat.set), 
                  more.args= list(validation=validation,admin_cens=admin_cens.set,structure=structure.set), reg = reg)
  submitJobs(ids, res = list(walltime = 6*60*60, memory = 16000, ncpus = 1, partition = "batch"), reg = reg)
  
}

#output directory
out.dir <- "..."

#settings:
n_training.set           <- 500
n_test.set               <- 500
n_validation.set         <- 2*10**5
n_cov_informative.set    <- 10
n_cov_noninformative.set <- 0
s.set                    <- 0.578
s_c.set                  <- 0.579
rate.set                 <- 0.864
offset.set               <- -0.002

event_estimator.set <- 'cox'
cens_distr.set <- 'Weibull'

source("sim_weibull.R")

n.per.ds.values <- c(25,50,100,250,500)

for(n.per.ds in n.per.ds.values){
  simulation.slurm.run.1(n.per.ds=n.per.ds, cens_estimator='cox')
}

for(n.per.ds in n.per.ds.values){
  simulation.slurm.run.1(n.per.ds=n.per.ds, cens_estimator='km')
}

#######################################################################################################
#######################################################################################################
### objective 3: assessing the performance of machine learning methods on low and high-dimensional data
#######################################################################################################
#######################################################################################################

#submit jobs to cluster
simulation.slurm.run.2 <- function(parameter, estimator){
  
  parameters <- as.list(environment())
  poi <- names(parameters[which(names(parameters) == parameter)])
  poi.value <- get(poi)
  
  validation <- get.WeibullData(n=n_validation.set, n_cov_informative=n_cov_informative.set, n_cov_noninformative=n_cov_noninformative.set, s=s.set, s_c=s_c.set, rate=rate.set, offset=offset.set, admin_cens=admin_cens.set, structure=structure.set, cens_distr=cens_distr.set, seed = T)$data
  
  reg <- makeRegistry(file.dir = paste0(out.dir,"/registry_",poi,"_",poi.value), seed = 1)
  ids <- batchMap(fun = simulation.slurm, args = list(i=1:iter.set, n_training=n_training.set, n_test=n_test.set, n_cov_informative=n_cov_informative.set, n_cov_noninformative=n_cov_noninformative.set, s=s.set, s_c=s_c.set, rate=rate.set, offset=offset.set, event_estimator=estimator, cens_estimator=estimator, cens_distr=cens_distr.set, min_G_hat=min_G_hat.set), 
                  more.args= list(validation=validation,admin_cens=admin_cens.set,structure=structure.set), reg = reg)
  submitJobs(ids, res = list(walltime = 12*60*60, memory = 64000, ncpus = 1, partition = "batch"), reg = reg)
  
}

############################
#3a: low-dimensional setting
############################

#output directory
out.dir <- "..."

#settings:
n_training.set           <- 500
n_test.set               <- 500
n_validation.set         <- 2*10**5
n_cov_informative.set    <- 10
s.set                    <- 0.578
s_c.set                  <- 0.579
rate.set                 <- 0.864
offset.set               <- -0.002

n_cov_noninformative.set <- 0
cens_distr.set <- "Weibull"

source("sim_weibull.R")


estimator.values <- c("rangersf","xgboost")

for(estimator in estimator.values){
  simulation.slurm.run.2(parameter="estimator", estimator)
}

#############################
#3b: high-dimensional setting
#############################

#output directory
out.dir <- "..."

#settings:
n_training.set           <- 500
n_test.set               <- 500
n_validation.set         <- 2*10**5
n_cov_informative.set    <- 10
s.set                    <- 0.578
s_c.set                  <- 0.579
rate.set                 <- 0.864
offset.set               <- -0.002

n_cov_noninformative.set <- 50
cens_distr.set <- "Weibull"

source("sim_weibull.R")


estimator.values <- c("rangersf","xgboost")

for(estimator in estimator.values){
  simulation.slurm.run.2(parameter="estimator", estimator)
}