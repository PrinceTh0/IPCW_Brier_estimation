
#path to slurm registries
registry.path <- "..."

#path to output files
data.path <- "..."

collect.data <- function(registry.path,registry.name){
  
  require(data.table)
  
  output.list <- readRDS(paste0(registry.path,"/",registry.name,"/result.list.rds"))
  
  ref.times <- setnames(as.data.table(output.list[[1]]$time), "ref_times")
  settings <- output.list[[1]]$settings
  
  l.output <- length(output.list)
  
  brier.uncensored <- data.table()
  
  brier.training <- data.table()
  brier.test <- data.table()
  brier.all <- data.table()
 
  obs.times.training <- data.table()
  obs.times.test <- data.table()
  obs.times.all <- data.table()
  
  true.times.training <- data.table()
  true.times.test <- data.table()
  true.times.all <- data.table()
  
  wse.training <- c()
  wse.test <- c()
  wse.all <- c()
 
  mse <- data.table()
  
  sample.parameters <- data.table()
  
  for(l in 1:length(output.list)){
    
    brier.uncensored <- cbind(brier.uncensored, setnames(output.list[[l]]$uncensored.brier.score, paste0("uncensored_",l)))
    
    brier.training <- cbind(brier.training, setnames(output.list[[l]]$ipcw.brier.score$training, paste0("training_",l)))
    brier.test <- cbind(brier.test, setnames(output.list[[l]]$ipcw.brier.score$test, paste0("test_",l)))
    brier.all <- cbind(brier.all, setnames(output.list[[l]]$ipcw.brier.score$all, paste0("all_",l)))
    
    obs.times.training <- cbind(obs.times.training, setnames(output.list[[l]]$obs.times[dataset == "train"][order(Ttilde),.(Ttilde)], paste0("training_",l)))
    obs.times.test <- cbind(obs.times.test, setnames(output.list[[l]]$obs.times[dataset == "test"][order(Ttilde),.(Ttilde)], paste0("test_",l)))
    obs.times.all <- cbind(obs.times.all, setnames(output.list[[l]]$obs.times[order(Ttilde),.(Ttilde)], paste0("all_",l)))
    
    true.times.training <- cbind(true.times.training, setnames(output.list[[l]]$true.times[dataset == "train"][order(Ttrue),.(Ttrue)], paste0("training_",l)))
    true.times.test <- cbind(true.times.test, setnames(output.list[[l]]$true.times[dataset == "test"][order(Ttrue),.(Ttrue)], paste0("test_",l)))
    true.times.all <- cbind(true.times.all, setnames(output.list[[l]]$true.times[order(Ttrue),.(Ttrue)], paste0("all_",l)))
    
    wse.training <- c(wse.training, output.list[[l]]$wse.score$train)
    wse.test <- c(wse.test, output.list[[l]]$wse.score$test)
    wse.all <- c(wse.all, output.list[[l]]$wse.score$all)
    
    mse <- cbind(mse, setnames(output.list[[l]]$mse, paste0("mse_",l)))
    
    sample.parameters <- rbind(sample.parameters, output.list[[l]]$sample.parameters)
  }
  
  return(list("brier.uncensored" = brier.uncensored,
              "brier.training" = brier.training,
              "brier.test" = brier.test,
              "brier.all" = brier.all,
              "obs.times.training" = obs.times.training,
              "obs.times.test" = obs.times.test,
              "obs.times.all" = obs.times.all,
              "true.times.training" = obs.times.training,
              "true.times.test" = obs.times.test,
              "true.times.all" = obs.times.all,
              "wse.training" = wse.training,
              "wse.test" = wse.test,
              "wse.all" = wse.all,
              "mse" = mse,
              "sample.parameters" = sample.parameters,
              "settings" = settings,
              "ref.times" = ref.times,
              "length.output" = l.output))
}

run.collect <- function(sim.parameter, registry.path, data.path){
  
  registry.list <- list.files(registry.path)
  registries <- registry.list[grepl(paste0("_",sim.parameter,"_"),registry.list)]
  
  if(!dir.exists(data.path)){dir.create(data.path)}
  
  for(registry in registries){
    coll.data <- collect.data(registry.path, registry)
    saveRDS(coll.data, file=file.path(data.path,registry))
  }
}

##objective 1
#run.collect("n.per.ds",registry.path,data.path)
#run.collect("rate",registry.path,data.path)
#run.collect("offset",registry.path,data.path)
#run.collect("n_training",registry.path,data.path)
#run.collect("n_test",registry.path,data.path)

#objective 2
#run.collect("cens_estimator",registry.path,data.path)

#objective 3
run.collect("estimator",registry.path,data.path)


