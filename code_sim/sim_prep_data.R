
#load packages
library(ggplot2)
library(ggpubr)
library(xtable)

#path to data
data.path.o1a <- ".../results/o1a"
data.path.o1b <- ".../results/o1b"
data.path.o2 <- ".../results/o2"
data.path.o3a <- ".../results/o3a"
data.path.o3b <- ".../results/o3b"


prep.plot <- function(sim.parameter, data.path, sim.parameter.type, other = NULL){
  
  require(data.table)
  require(ggplot2)
  
  file.list <- list.files(data.path)
  if(is.null(other)){
    files <- file.list[grepl(paste0("_",sim.parameter,"_"),file.list)]
  }else if(!is.null(other)){
    files <- file.list[grepl(paste0(other),file.list[grepl(paste0("_",sim.parameter,"_"),file.list)])]
  }
  
  wse.dt <- data.table()
  wse.wide.dt <- data.table()
  bs.dt <- data.table()
  times.obs.dt <- data.table()
  
  plot.list <- list()
  
  for(file in files){
    
    coll.data <- readRDS(file = file.path(data.path, file))
    
    if(is.null(other)){
      sim.parameter.value <- gsub(paste0(sim.parameter,"_"),"",gsub("registry_","",file))
    }else if(!is.null(other)){
      sim.parameter.value <- gsub(other,"",gsub(paste0(sim.parameter,"_"),"",gsub("registry_","",file)))
    }
    
    wse.training <- data.table("wse"=coll.data$wse.training,"dataset"="training")
    wse.test <- data.table("wse"=coll.data$wse.test,"dataset"="test")
    wse.all <- data.table("wse"=coll.data$wse.all,"dataset"="all")
    
    wse <- rbind(wse.training,wse.test,wse.all)
    
    if(sim.parameter.type == "numeric"){
      wse[, (sim.parameter) := as.numeric(sim.parameter.value)]
    }else{
      wse[, (sim.parameter) := sim.parameter.value]
    }
    wse.dt <- rbind(wse.dt, wse)
    
    wse.wide <- data.table("training"=coll.data$wse.training,
                              "test"=coll.data$wse.test,
                              "all"=coll.data$wse.all)
    
    if(sim.parameter.type == "numeric"){
      wse.wide[, (sim.parameter) := as.numeric(sim.parameter.value)]
    }else{
      wse.wide[, (sim.parameter) := sim.parameter.value]
    }
    
    #p <- ggplot(wse.wide) +
    #  geom_point(aes(x=test, y=training)) +
    #  lims(x=c(0,max(c(wse.wide$training,wse.wide$test))), y=c(0,max(c(wse.wide$training,wse.wide$test)))) +
    #  ggtitle(label = paste0(sym(sim.parameter)," = ", sim.parameter.value)) +
    #  coord_fixed() +
    #  theme_classic()
    
    #append(plot.list, p)
    
    wse.wide.dt <- rbind(wse.wide.dt, wse.wide)
    
    bs.obs.training <- copy(coll.data$brier.training)
    bs.obs.test <- copy(coll.data$brier.test)
    bs.obs.all <- copy(coll.data$brier.all)
    
    bs.obs.training[, training_bs_mean := apply(.SD, 1, function(x){mean(x, na.rm = T)})]
    bs.obs.test[, test_bs_mean := apply(.SD, 1, function(x){mean(x, na.rm = T)})]
    bs.obs.all[, all_bs_mean := apply(.SD, 1, function(x){mean(x, na.rm = T)})]
    
    bs.obs.training[, training_bs_sd := apply(.SD, 1, function(x){sd(x, na.rm = T)})]
    bs.obs.test[, test_bs_sd := apply(.SD, 1, function(x){sd(x, na.rm = T)})]
    bs.obs.all[, all_bs_sd := apply(.SD, 1, function(x){sd(x, na.rm = T)})]
    
    bs.obs <- cbind(setnames(coll.data$ref.times, "ref_times"), bs.obs.training, bs.obs.test, bs.obs.all)
    
    bs.uncensored <- copy(coll.data$brier.uncensored)
    bs.uncensored <- bs.uncensored[, uncensored_bs_mean := apply(.SD, 1, function(x){mean(x, na.rm = T)})]
    bs.uncensored <- bs.uncensored[, uncensored_bs_sd := apply(.SD, 1, function(x){sd(x, na.rm = T)})][,.(uncensored_bs_mean, uncensored_bs_sd)]
    
    bs.mean <- cbind(bs.obs, bs.uncensored)
    bs.mean <- melt(setnames(bs.mean[,c("ref_times","training_bs_mean","test_bs_mean","all_bs_mean","uncensored_bs_mean")], c("ref_times","training","test","all","uncensored")), id.vars = "ref_times", variable.name = "dataset", value.name = "bs")
    
    bs.sd <- cbind(bs.obs, bs.uncensored)
    bs.sd <- melt(setnames(bs.sd[,c("ref_times","training_bs_sd","test_bs_sd","all_bs_sd","uncensored_bs_sd")], c("ref_times","training","test","all","uncensored")), id.vars = "ref_times", variable.name = "dataset", value.name = "bs_sd")
    
    bs <- merge.data.table(bs.mean, bs.sd, by = c("ref_times","dataset"))
    
    if(sim.parameter.type == "numeric"){
      bs[, (sim.parameter) := as.numeric(sim.parameter.value)]
    }else{
      bs[, (sim.parameter) := sim.parameter.value]
    }
    bs.dt <- rbind(bs.dt, bs)
    
    times.obs.training <- melt.data.table(coll.data$obs.times.training, id.vars = NULL, variable.name = "dataset", value.name = "Ttilde")[, dataset := "training"][order(Ttilde)]
    times.obs.test <- melt.data.table(coll.data$obs.times.test, id.vars = NULL, variable.name = "dataset", value.name = "Ttilde")[, dataset := "test"][order(Ttilde)]
    
    times.obs <- rbind(times.obs.training, times.obs.test)
    times.obs[, dataset := as.factor(dataset)]
    if(sim.parameter.type == "numeric"){
      times.obs[, (sim.parameter) := as.numeric(sim.parameter.value)]
    }else{
      times.obs[, (sim.parameter) := sim.parameter.value]
    }
    
    times.obs[, q75 := quantile(Ttilde, 0.75), by = c(sim.parameter, "dataset")]
    times.obs[, q95 := quantile(Ttilde, 0.95), by = c(sim.parameter, "dataset")]
    times.obs[, q99 := quantile(Ttilde, 0.99), by = c(sim.parameter, "dataset")]
    
    times.obs <- unique(times.obs[, -c("Ttilde")])
    times.obs <- melt(times.obs, by = c(sim.parameter,"dataset"), id.vars = c(sim.parameter,"dataset"))
    times.obs.dt <- rbind(times.obs.dt, times.obs)
    
    #bs.dt <- merge.data.table(bs.dt, times.obs.dt, by = c(sim.parameter,"dataset"), all.x = T)
    
  }
  
  #newlevels <- unique(wse.dt[[sim.parameter]])[order(unique(wse.dt[[sim.parameter]]))]
  wse.dt[, (sim.parameter) := lapply(.SD, function(x){as.factor(x)}), .SD = (sim.parameter)]
  wse.wide.dt[, (sim.parameter) := lapply(.SD, function(x){as.factor(x)}), .SD = (sim.parameter)]
  bs.dt[, (sim.parameter) := lapply(.SD, function(x){as.factor(x)}), .SD = (sim.parameter)]
  times.obs.dt[, (sim.parameter) := lapply(.SD, function(x){as.factor(x)}), .SD = (sim.parameter)]
  
  #plot.wse <- ggplot(wse.dt) +
  #  geom_violin(aes(x=factor(dataset,ordered=T,levels=c("training","test","all")),y=wse, fill = dataset), show.legend = F, draw_quantiles = T, size = 0.2) +
  #  labs(x="fit dataset", y = "wse") +
  #  scale_fill_manual(values = c("#0072b2","#d55e00","#f0e442")) +
  #  lims(y=c(0,0.5)) +
  #  facet_wrap(sym(sim.parameter), scales = "free_x", labeller = label_both) +
  #  theme_classic()
  
  #plot.wse.scatter <- ggplot(wse.wide.dt) +
  #  geom_point(aes(x=test, y=training, color = wse.wide.dt[[sim.parameter]]), show.legend = F) +
  #  facet_wrap(sym(sim.parameter), scales = "free", labeller = label_both) +
  #  #lims(x = c(0,max(wse.wide.dt.min.max$max_test,wse.wide.dt.min.max$max_training)), y = c(0,max(wse.wide.dt.min.max$max_test,wse.wide.dt.min.max$max_training))) +
  #  #lims(x=c(0,0.5), y=c(0,0.5)) +
  #  theme_classic()
  
  #plot.bs <- ggplot(bs.dt, aes(x = ref_times)) +
  #  geom_line(aes(y = bs, color = dataset)) +
  #  geom_ribbon(aes(y = bs, ymin = pmax(0,(bs - bs_sd)), ymax = bs + bs_sd, fill = dataset), alpha = 0.3) +
  #  scale_color_manual(values = c("#f0e442","#d55e00","#0072b2","#009e73")) +
  #  scale_fill_manual(values = c("#f0e442","#d55e00","#0072b2","#009e73")) +
  #  labs(x = "time", y = "brier score") +
  #  facet_wrap(sym(sim.parameter), scales = "free_x", labeller = label_both) +
  #  geom_vline(data = times.obs.dt, aes(xintercept = value, linetype = variable, color = dataset)) +
  #  scale_linetype_manual(values = c("solid","dotdash","dotted"), name = "time quantiles") +
  #  lims(x=c(0,7),y=c(0,NA)) +
  #  theme_classic()
  
  return(list(#"wse" = plot.wse, 
              #"bs" = plot.bs, 
              #"wse.scatter" = plot.wse.scatter, 
              "wse.dt" = wse.dt, 
              "wse.wide.dt" = wse.wide.dt, 
              "bs.dt" = bs.dt
              ))
  
}


wse.dt.n.per.ds.o1a <- prep.plot("n.per.ds", data.path.o1a, "numeric")$wse.dt
wse.dt.n.training.o1a <- prep.plot("n.training", data.path.o1a, "numeric")$wse.dt
wse.dt.n.test.o1a <- prep.plot("n.test", data.path.o1a, "numeric")$wse.dt
wse.dt.rate.o1a <- prep.plot("rate", data.path.o1a, "numeric")$wse.dt

wse.dt.n.per.ds.o1b <- prep.plot("n.per.ds", data.path.o1b, "numeric")$wse.dt
wse.dt.n.training.o1b <- prep.plot("n.training", data.path.o1b, "numeric")$wse.dt
wse.dt.n.test.o1b <- prep.plot("n.test", data.path.o1b, "numeric")$wse.dt
wse.dt.offset.o1b <- prep.plot("offset", data.path.o1b, "numeric")$wse.dt

wse.dt.n.per.ds.cox <- prep.plot("n.per.ds", data.path.o2, "numeric",  "_cens_estimator_cox")$wse.dt
wse.dt.n.per.ds.km <- prep.plot("n.per.ds", data.path.o2, "numeric",  "_cens_estimator_km")$wse.dt

wse.dt.estimator.o3a <- prep.plot("estimator", data.path.o3a, "character")$wse.dt

wse.dt.estimator.o3b <- prep.plot("estimator", data.path.o3b, "character")$wse.dt


#list for download
sim_o1a_plot_data_list <- list("wse.dt.n.per.ds" = wse.dt.n.per.ds.o1a, 
                            "wse.dt.n.training" = wse.dt.n.training.o1a,
                            "wse.dt.n.test" = wse.dt.n.test.o1a,
                            "wse.dt.rate" = wse.dt.rate.o1a)

sim_o1b_plot_data_list <- list("wse.dt.n.per.ds" = wse.dt.n.per.ds.o1b, 
                               "wse.dt.n.training" = wse.dt.n.training.o1b,
                               "wse.dt.n.test" = wse.dt.n.test.o1b,
                               "wse.dt.offset" = wse.dt.offset.o1b)

sim_o2_plot_data_list <- list("wse.dt.n.per.ds.cox" = wse.dt.n.per.ds.cox, 
                               "wse.dt.n.per.ds.km" = wse.dt.n.per.ds.km)

sim_o3a_plot_data_list <- list("wse.dt.estimator" = wse.dt.estimator.o3a)

sim_o3b_plot_data_list <- list("wse.dt.estimator" = wse.dt.estimator.o3b)


saveRDS(sim_o1a_plot_data_list, file = ".../plot_o1a.rds")
saveRDS(sim_o1b_plot_data_list, file = ".../plot_o1b.rds")
saveRDS(sim_o2_plot_data_list, file = ".../plot_o2.rds")
saveRDS(sim_o3a_plot_data_list, file = ".../plot_o3a.rds")
saveRDS(sim_o3b_plot_data_list, file = ".../plot_o3b.rds")


##check values of simulated parameters
check.settings <- function(sim.parameter, data.path, sim.parameter.type, other = NULL){
  
  require(data.table)
  require(ggplot2)
  require(ggpubr)
  
  file.list <- list.files(data.path)
  if(is.null(other)){
    files <- file.list[grepl(paste0("_",sim.parameter,"_"),file.list)]
  }else if(!is.null(other)){
    files <- file.list[grepl(paste0(other),file.list[grepl(paste0("_",sim.parameter,"_"),file.list)])]
  }
  
  parameter.mean.dt <- data.table()
  parameter.sd.dt <- data.table()
  parameter.plot.dt <- data.table()
  n.samples.dt <- data.table()
  
  variables <- c("%cens_data","%cens_training","%cens_test","R2")
  
  for(file in files){
    
    coll.data <- readRDS(file = file.path(data.path, file))
    
    if(is.null(other)){
      sim.parameter.value <- gsub(paste0(sim.parameter,"_"),"",gsub("registry_","",file))
    }else if(!is.null(other)){
      sim.parameter.value <- gsub(other,"",gsub(paste0(sim.parameter,"_"),"",gsub("registry_","",file)))
    }
    
    if(sim.parameter.type == "numeric"){
      n.samples <- data.table("sim.parameter" = sim.parameter, "sim.parameter.value" = as.numeric(sim.parameter.value), "n.samples" = coll.data$length.output)
    }else{
      n.samples <- data.table("sim.parameter" = sim.parameter, "sim.parameter.value" = sim.parameter.value, "n.samples" = coll.data$length.output)
    }
    n.samples.dt <- rbind(n.samples, n.samples.dt)
    
    parameter.mean <- coll.data$sample.parameters[, lapply(.SD, mean), .SD = variables][, setnames(.SD, variables, paste0(variables,"_mean"))]
    if(sim.parameter.type == "numeric"){
      parameter.mean[, (sim.parameter) := as.numeric(sim.parameter.value)]
    }else{
      parameter.mean[, (sim.parameter) := sim.parameter.value]
    }
    parameter.mean.dt <- rbind(parameter.mean.dt,parameter.mean)
    
    parameter.sd <- coll.data$sample.parameters[, lapply(.SD, sd), .SD = variables][, setnames(.SD, variables, paste0(variables,"_sd"))]
    if(sim.parameter.type == "numeric"){
      parameter.sd[, (sim.parameter) := as.numeric(sim.parameter.value)]
    }else{
      parameter.sd[, (sim.parameter) := sim.parameter.value]
    }
    parameter.sd.dt <- rbind(parameter.sd.dt,parameter.sd)
    
    parameter.plot <- coll.data$sample.parameters[, ..variables]
    if(sim.parameter.type == "numeric"){
      parameter.plot[, (sim.parameter) := as.numeric(sim.parameter.value)]
    }else{
      parameter.plot[, (sim.parameter) := sim.parameter.value]
    }
    parameter.plot.dt <- rbind(parameter.plot.dt, parameter.plot)
    
  }
  
  parameter.plot.dt[, (sim.parameter) := lapply(.SD, function(x){as.factor(x)}), .SD = (sim.parameter)]
  
  for(var in variables){
    
    y <- sym(var)
    
    plot.temp <- ggplot(parameter.plot.dt) +
      geom_boxplot(aes(y=!!y)) +
      #geom_boxplot(aes(x=parameter.plot.dt[[sim.parameter]],y=`%cens_training`)) +
      facet_wrap(sym(sim.parameter), scales = "free_x", labeller = label_both) +
      labs(x = NULL, y = var) +
      scale_x_discrete(breaks = NULL) +
      scale_y_continuous(limits = c(0,1), breaks = seq(0,1,0.1), labels = seq(0,1,0.1)) +
      theme_classic()
    
    assign(paste0("plot.",gsub("%","",var)), plot.temp)
    
  }
  
  plot.grid <- ggarrange(plot.cens_data,plot.cens_training,plot.cens_test,plot.R2)
  
  return(list("settings.mean" = parameter.mean.dt,
              "settings.sd" = parameter.sd.dt,
              "plot.grid" = plot.grid,
              "plot.cens.data" = plot.cens_data,
              "plot.cens.training" = plot.cens_training,
              "plot.cens.test" = plot.cens_test,
              "plot.R2" = plot.R2,
              "n.samples" = n.samples.dt))
}

#check.n.per.ds <- check.settings("n.per.ds",data.path.o1a,"numeric")
#check.n.training <- check.settings("n_training",data.path.o1a,"numeric")
#check.n.test <- check.settings("n_test",data.path.o1a,"numeric")
#check.rate <- check.settings("rate",data.path.o1a,"numeric")

#check.n.per.ds <- check.settings("n.per.ds",data.path.o1b,"numeric")
#check.n.training <- check.settings("n_training",data.path.o1b,"numeric")
#check.n.test <- check.settings("n_test",data.path.o1b,"numeric")
#check.offset <- check.settings("offset",data.path.o1b,"numeric")

#check.estimator <- check.settings("estimator",data.path.o3a,"character")
#check.estimator <- check.settings("estimator",data.path.o3b,"character")

#check.n.per.ds.cox <- check.settings("n.per.ds",data.path.o2,"numeric", "_cens_estimator_cox")
#check.n.per.ds.km <- check.settings("n.per.ds",data.path.o2,"numeric", "_cens_estimator_km")