
library(data.table)
library(ggplot2)

dir.res <- "..." #results dírectory

rangersf_tuning <- readRDS(file.path(dir.res,"registry_rangersf_5_full_tuningTRUE"))
xgboost_tuning <- readRDS(file.path(dir.res,"registry_xgboost_5_full_tuningTRUE"))

cox_notuning <- readRDS(file.path(dir.res,"registry_cox_5_full_tuningFALSE"))
rangersf_notuning <- readRDS(file.path(dir.res,"registry_rangersf_5_full_tuningFALSE"))
xgboost_notuning <- readRDS(file.path(dir.res,"registry_xgboost_5_full_tuningFALSE"))


create.bs.dt <- function(data.list, tunedYN){
  
  times <- seq(0,156,1)
  
  dt.training <- data.table("times"=times)
  dt.test <- data.table("times"=times)
  dt.all <- data.table("times"=times)
  
  dt.training.long <- melt.data.table(cbind(dt.training,data.list$bs.training), id.vars = "times", variable.name = "sample", value.name = "bs")
  dt.training.long[, "bs.mean" := mean(bs), by = times][, "bs.sd" := sd(bs), by = times]
  dt.training.long[, dataset := "training"][, sample := gsub("bootstrap_", "training_", sample)]
  
  dt.test.long <- melt.data.table(cbind(dt.test,data.list$bs.test), id.vars = "times", variable.name = "sample", value.name = "bs")
  dt.test.long[, "bs.mean" := mean(bs), by = times][, "bs.sd" := sd(bs), by = times]
  dt.test.long[, dataset := "test"][, sample := gsub("bootstrap_", "test_", sample)]
  
  dt.all.long <- melt.data.table(cbind(dt.all,data.list$bs.all), id.vars = "times", variable.name = "sample", value.name = "bs")
  dt.all.long[, "bs.mean" := mean(bs), by = times][, "bs.sd" := sd(bs), by = times]
  dt.all.long[, dataset := "all"][, sample := gsub("bootstrap_", "all_", sample)]
  
  dt.long <- rbind(dt.training.long, dt.test.long, dt.all.long)
  dt.long[, method := data.list$method]
  dt.long[, tuned := ifelse(tunedYN, "Yes", "No")]
  
  return(dt.long)
}

#cox.dt.long.t <- create.bs.dt(cox_tuning)
rangersf.dt.long.t <- create.bs.dt(rangersf_tuning, TRUE)
xgboost.dt.long.t <- create.bs.dt(xgboost_tuning, TRUE)

cox.dt.long.nt <- create.bs.dt(cox_notuning, FALSE)
rangersf.dt.long.nt <- create.bs.dt(rangersf_notuning, FALSE)
xgboost.dt.long.nt <- create.bs.dt(xgboost_notuning, FALSE)

#list for download
seer_plot_data_list <- list("rangersf.dt.long.t" = rangersf.dt.long.t, 
                    "xgboost.dt.long.t" = xgboost.dt.long.t,
                    "cox.dt.long.nt" = cox.dt.long.nt,
                    "rangersf.dt.long.nt" = rangersf.dt.long.nt,
                    "xgboost.dt.long.nt" = xgboost.dt.long.nt)

saveRDS(seer_plot_data_list, file = ".../seer_plot_data_list.rds")
