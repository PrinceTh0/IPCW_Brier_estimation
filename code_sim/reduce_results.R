
require(batchtools)

#directory with registries to reduce
reg.dir <- "..."

reg.files <- list.files(reg.dir)

for(reg.file in reg.files){
  
  reg <- loadRegistry(file.dir = file.path(reg.dir,reg.file), writeable = T)
  result.list <- reduceResultsList(reg = reg)
  saveRDS(result.list, file = file.path(reg.dir,reg.file,"result.list.rds"))
  
}
