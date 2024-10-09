
#path to slurm registries
registry.path <- "..."

#path to reformatted files
data.path <- "..."

collect.data <- function(registry.path,registry.name){
  
  require(data.table)
  
  output.list <- readRDS(paste0(registry.path,"/",registry.name,"/result.list.rds"))
  
  method <- output.list[[1]]$method
  weights <- output.list[[1]]$weights
  censoring.model <- output.list[[1]]$censoring.model
  n.bootstrap <- length(output.list)
  
  bs.training <- data.table()
  bs.test <- data.table()
  bs.all <- data.table()
  samples.training <- data.table()
  samples.test <- data.table()
  samples.all <- data.table()
  
  for(l in 1:length(output.list)){
    
    bs.training <- cbind(bs.training, setnames(as.data.table(output.list[[l]]$bs.training), paste0("bootstrap_",l)))
    bs.test <- cbind(bs.test, setnames(as.data.table(output.list[[l]]$bs.test), paste0("bootstrap_",l)))
    bs.all <- cbind(bs.all, setnames(as.data.table(output.list[[l]]$bs.data), paste0("bootstrap_",l)))
    
    samples.training <- cbind(samples.training, setnames(as.data.table(output.list[[l]]$samples.training), paste0("bootstrap_",l)))
    samples.test <- cbind(samples.test, setnames(as.data.table(output.list[[l]]$samples.test), paste0("bootstrap_",l)))
    samples.all <- cbind(samples.all, setnames(as.data.table(output.list[[l]]$samples.data), paste0("bootstrap_",l)))
    
  }
  
  return(list("method"=method,
              "weights"=weights,
              "censoring.model"=censoring.model,
              "bs.training"=bs.training,
              "bs.test"=bs.test,
              "bs.all"=bs.all,
              "samples.training"=samples.training,
              "samples.test"=samples.test,
              "samples.all"=samples.all))
}


run.collect <- function(method, registry.path, data.path){
  
  registry.list <- list.files(registry.path)
  registries <- registry.list[grepl(paste0("_",method),registry.list)]
  
  if(!dir.exists(data.path)){dir.create(data.path)}
  
  for(registry in registries){
    coll.data <- collect.data(registry.path, registry)
    saveRDS(coll.data, file=file.path(data.path,registry))
  }
}

run.collect("cox",registry.path,data.path)
run.collect("rangersf",registry.path,data.path)
run.collect("xgboost",registry.path,data.path)