rm(list = ls())

require("batchtools")

#options(batchtools.progress = FALSE)

tuning.slurm <- function(n,n_cov_informative,n_cov_noninformative,s,s_c,rate,offset,admin_cens,structure,cens_distr){
  
  source('sim_weibull.R')
  
  require("data.table")
  require("evd")
  require("MASS")
  
  #set.seed(1234)
  
  Weibull.data <- get.WeibullData(n,n_cov_informative,n_cov_noninformative,s,s_c,rate,offset,admin_cens,structure,cens_distr,seed=T)
  
  out <- list("dt" = data.table("s" = s, "s_c" = s_c, "offset" = offset, "rate" = rate, "R2" = Weibull.data$R2, "R2c" = Weibull.data$R2c, "cens" = 1-sum(Weibull.data$data$event)/nrow(Weibull.data$data)))

  return(out)
}

# parameters for data generation
n.set                    <- 10**6
n_cov_informative.set    <- 10
n_cov_noninformative.set <- 0
#s.set                    <- 1
#s_c.set                  <- 1
#rate.set                 <- 1
#offset.set               <- 0
admin_cens.set           <- F
structure.set            <- "compoundSymmetry"
#cens_distr.set           <- "exponential"

#output directory
out.dir <- file.path("~/scratch/survival.simulation.V14/tuning/weibull.noisy/output")


tune <- function(s, s_c, rate, offset, cens_distr){
  
  parameters <- as.list(environment())
  #parameter.fixed <- names(which(lengths(parameters)==1))
  parameter.dynamic <- names(which(lengths(parameters)>1))
  #parameter.fixed.value <- get(parameter.fixed)

  reg <- makeRegistry(file.dir = paste0(out.dir,"/registry_",parameter.dynamic,"_",cens_distr), seed = 1)
  ids <- batchMap(fun = tuning.slurm, args = list(n=n.set,n_cov_informative=n_cov_informative.set,n_cov_noninformative=n_cov_noninformative.set,s=s,s_c=s_c,rate=rate,offset=offset,admin_cens=admin_cens.set,structure=structure.set,cens_distr=cens_distr), reg = reg)
  submitJobs(ids, res = list(walltime = 0.3*60*60, memory = 12000, ncpus = 1, partition = "batch"), reg = reg)
}

tune(s = seq(0.001,3,0.001), s_c = 1, rate = 1, offset = 0, cens_distr = "Weibull")
#tune(s = 0.578, s_c = seq(0.001,3,0.001), rate = 1, offset = 0, cens_distr = "Weibull")

#tune(s = 0.578, s_c = 0.579, rate = seq(0.001,3.500,0.001), offset = 0, cens_distr = "exponential")
#tune(s = 0.578, s_c = 0.579, rate = 1, offset = seq(-1.500,2.500,0.001), cens_distr = "Weibull")

#absolute beta values censoring = event

#R2/R2c 50% -> s~0.578/0.579

#exponential:
#censoring rate 80% -> rate~0.219
#censoring rate 50% -> rate~0.864
#censoring rate 20% -> rate~3.013

#Weibull:
#censoring rate 80% -> offset~-0.803
#censoring rate 50% -> offset~-0.002
#censoring rate 20% -> offset~0.8