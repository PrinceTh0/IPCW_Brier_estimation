library(data.table)
library(ggplot2)
library(ggpubr)
library(extrafont)
library(xtable)

loadfonts(device = c("postscript"))

data.path <- "..." #path to input
plot.path <- "..." #path to output

plot_o1a_data <- readRDS(file = paste0(data.path,"/plot_o1a.rds"))
plot_o1b_data <- readRDS(file = paste0(data.path,"/plot_o1b.rds"))
plot_o2_data <- readRDS(file = paste0(data.path,"/plot_o2.rds"))
plot_o3a_data <- readRDS(file = paste0(data.path,"/plot_o3a.rds"))
plot_o3b_data <- readRDS(file = paste0(data.path,"/plot_o3b.rds"))


for(l in 1:length(plot_o1a_data)){
  assign(paste0("o1a.",names(plot_o1a_data)[l]),plot_o1a_data[[l]])
}

for(l in 1:length(plot_o1b_data)){
  assign(paste0("o1b.",names(plot_o1b_data)[l]),plot_o1b_data[[l]])
}

for(l in 1:length(plot_o2_data)){
  assign(paste0("o2.",names(plot_o2_data)[l]),plot_o2_data[[l]])
}

for(l in 1:length(plot_o3a_data)){
  assign(paste0("o3a.",names(plot_o3a_data)[l]),plot_o3a_data[[l]])
}

for(l in 1:length(plot_o3b_data)){
  assign(paste0("o3b.",names(plot_o3b_data)[l]),plot_o3b_data[[l]])
}


##set theme
theme.settings <- theme(line = element_line(linewidth = 1),
                        text = element_text(size = 10, family = "serif"),
                        rect = element_rect(linewidth = 1),
                        plot.title = element_text(size = rel(3), hjust = 0.5),
                        legend.title = element_text(size = rel(2.75)),
                        legend.text = element_text(size = rel(2.75)),
                        legend.key.spacing.x = unit(10, "pt"),
                        axis.title = element_text(size = rel(2.75)),
                        axis.text = element_text(size = rel(2.75)),
                        strip.text = element_text(size = rel(3)),
                        strip.background = element_blank(),
                        panel.grid.major.y = element_line(color = "grey", linewidth = 0.2))


##o1a

# n per ds
o1a.wse.dt.n.per.ds[, n.total := as.factor(2*as.numeric(levels(n.per.ds))[n.per.ds])]

#o1a.wse.table.n.per.ds <- o1a.wse.dt.n.per.ds[, .(wse.median = quantile(wse, 0.5, na.rm = TRUE), wse.iqr = IQR(wse, na.rm = TRUE)), by = .(n.total, dataset)][, n.total := as.numeric(levels(n.total))[n.total]][dataset == "all", dataset := "combined"]
#setnames(o1a.wse.table.n.per.ds, names(o1a.wse.table.n.per.ds),c("n","dataset","median","IQR"))
#setorder(o1a.wse.table.n.per.ds,n,dataset)
#print(xtable(o1a.wse.table.n.per.ds, caption = c("numerical results for the Weibull-Cox scenario (ii) with varying n"), align = c("c","l","l","c","c"), digits = c(0,0,0,6,6)), include.rownames = FALSE)

o1a.p.n.per.ds <- ggplot(o1a.wse.dt.n.per.ds) +
  geom_violin(aes(x=n.total,y=log(wse), fill = dataset), show.legend = T, draw_quantiles = c(0.5)) +
  scale_fill_manual(values = c("#00788C","#B07C4F","#76A165"), labels = c("combined", "test", "training"), name = "Fitting dataset:") +
  labs(x="Total size of the dataset", y="log(RWSE)") +
  #lims(y=c(0,0.5)) +
  theme_classic() +
  theme.settings + 
  theme(plot.margin = unit(c(0.7,0.7,0.7,0.7), "cm"))


#n.training

#o1a.wse.table.n.training <- o1a.wse.dt.n.training[, .(wse.median = quantile(wse, 0.5, na.rm = TRUE), wse.iqr = IQR(wse, na.rm = TRUE)), by = .(n.training, dataset)][, n.training := as.numeric(levels(n.training))[n.training]][dataset == "all", dataset := "combined"]
#setnames(o1a.wse.table.n.training, names(o1a.wse.table.n.training),c("n.training","dataset","median","IQR"))
#setorder(o1a.wse.table.n.training,n.training,dataset)
#print(xtable(wse.table.dt.o1.b, caption = c("numerical results for the Weibull-Cox scenario (ii) with varying n.training"), align = c("c","l","l","c","c"), digits = c(0,0,0,6,6)), include.rownames = FALSE)


o1a.p.n.training <- ggplot(o1a.wse.dt.n.training) +
  geom_violin(aes(x=n.training,y=log(wse), fill = dataset), show.legend = T, draw_quantiles = c(0.5)) +
  scale_fill_manual(values = c("#00788C","#B07C4F","#76A165"), labels = c("combined", "test", "training"), name = "Fitting dataset:") +
  labs(x="Size of the training set", y="log(RWSE)") +
  #lims(y=c(0,0.5)) +
  theme_classic() +
  theme.settings + 
  theme(plot.margin = unit(c(0.7,0.7,0.7,0.7), "cm"))


#n.test

#o1a.wse.table.n.test <- o1a.wse.dt.n.test[, .(wse.median = quantile(wse, 0.5, na.rm = TRUE), wse.iqr = IQR(wse, na.rm = TRUE)), by = .(n.test, dataset)][, n.test := as.numeric(levels(n.test))[n.test]][dataset == "all", dataset := "combined"]
#setnames(o1a.wse.table.n.test, names(o1a.wse.table.n.test),c("n.test","dataset","median","IQR"))
#setorder(o1a.wse.table.n.test,n.test,dataset)
#print(xtable(o1a.wse.table.n.test, caption = c("numerical results for the Weibull-Cox scenario (ii) with varying n.test"), align = c("c","l","l","c","c"), digits = c(0,0,0,6,6)), include.rownames = FALSE)

o1a.p.n.test <- ggplot(o1a.wse.dt.n.test) +
  geom_violin(aes(x=n.test,y=log(wse), fill = dataset), show.legend = T, draw_quantiles = c(0.5)) +
  scale_fill_manual(values = c("#00788C","#B07C4F","#76A165"), labels = c("combined", "test", "training"), name = "Fitting dataset:") +
  labs(x="Size of the test set", y="log(RWSE)") +
  #lims(y=c(0,0.15)) +
  theme_classic() +
  theme.settings + 
  theme(plot.margin = unit(c(0.7,0.7,0.7,0.7), "cm"))


## censoring rate

#o1a.wse.table.rate <- o1a.wse.dt.rate[, .(wse.median = quantile(wse, 0.5, na.rm = TRUE), wse.iqr = IQR(wse, na.rm = TRUE)), by = .(rate, dataset)][, rate := as.numeric(levels(rate))[rate]][order(rate)][dataset == "all", dataset := "combined"]
#o1a.wse.table.rate[rate == "0.219", rate := "0.2"][rate == "0.864", rate := "0.5"][rate == "3.013", rate := "0.8"][, rate := as.factor(rate)]
#setnames(o1a.wse.table.rate, names(o1a.wse.table.rate),c("censoring rate","dataset","median","IQR"))
#setorder(o1a.wse.table.rate,`censoring rate`,dataset)
#print(xtable(o1a.wse.table.rate, caption = c("numerical results for the baseline scenario (i) with varying censoring rate"), align = c("c","l","l","c","c"), digits = c(0,0,0,6,6)), include.rownames = FALSE)

o1a.p.censoring.rate <- ggplot(o1a.wse.dt.rate[rate == "0.219", cens.rate := "0.2"][rate == "0.864", cens.rate := "0.5"][rate == "3.013", cens.rate := "0.8"][, cens.rate := as.factor(cens.rate)]) +
  geom_violin(aes(x=cens.rate,y=log(wse), fill = dataset), show.legend = T, draw_quantiles = c(0.5)) +
  scale_fill_manual(values = c("#00788C","#B07C4F","#76A165"), labels = c("combined", "test", "training"), name = "Fitting dataset:") +
  labs(x="Censoring rate", y="log(RWSE)") +
  #lims(y=c(0,0.1)) +
  theme_classic() +
  theme.settings + 
  theme(plot.margin = unit(c(0.7,0.7,0.7,0.7), "cm"))

plot_o1a <- ggarrange(o1a.p.n.per.ds, o1a.p.n.training, o1a.p.n.test, o1a.p.censoring.rate, common.legend = TRUE, legend = "bottom", labels = c("(a)","(b)","(c)","(d)"), hjust = 0, font.label = list(size = 20, face = "bold"))

ggsave(filename = "plot_o1a.eps", plot = plot_o1a, device = cairo_ps, path = plot.path, units = "in", width = 14, height = 8, dpi = 600)



##o1b

# n per ds
o1b.wse.dt.n.per.ds[, n.total := as.factor(2*as.numeric(levels(n.per.ds))[n.per.ds])]

#o1b.wse.table.n.per.ds <- o1b.wse.dt.n.per.ds[, .(wse.median = quantile(wse, 0.5, na.rm = TRUE), wse.iqr = IQR(wse, na.rm = TRUE)), by = .(n.total, dataset)][, n.total := as.numeric(levels(n.total))[n.total]][dataset == "all", dataset := "combined"]
#setnames(o1b.wse.table.n.per.ds, names(o1b.wse.table.n.per.ds),c("n","dataset","median","IQR"))
#setorder(o1b.wse.table.n.per.ds,n,dataset)
#print(xtable(o1b.wse.table.n.per.ds, caption = c("numerical results for the Weibull-Cox scenario (ii) with varying n"), align = c("c","l","l","c","c"), digits = c(0,0,0,6,6)), include.rownames = FALSE)

o1b.p.n.per.ds <- ggplot(o1b.wse.dt.n.per.ds) +
  geom_violin(aes(x=n.total,y=log(wse), fill = dataset), show.legend = T, draw_quantiles = c(0.5)) +
  scale_fill_manual(values = c("#00788C","#B07C4F","#76A165"), labels = c("combined", "test", "training"), name = "Fitting dataset:") +
  labs(x="Total size of the dataset", y="log(RWSE)") +
  #lims(y=c(0,0.5)) +
  theme_classic() +
  theme.settings + 
  theme(plot.margin = unit(c(0.7,0.7,0.7,0.7), "cm"))


#n.training

#o1b.wse.table.n.training <- o1b.wse.dt.n.training[, .(wse.median = quantile(wse, 0.5, na.rm = TRUE), wse.iqr = IQR(wse, na.rm = TRUE)), by = .(n.training, dataset)][, n.training := as.numeric(levels(n.training))[n.training]][dataset == "all", dataset := "combined"]
#setnames(o1b.wse.table.n.training, names(o1b.wse.table.n.training),c("n.training","dataset","median","IQR"))
#setorder(o1b.wse.table.n.training,n.training,dataset)
#print(xtable(o1b.wse.table.n.training, caption = c("numerical results for the Weibull-Cox scenario (ii) with varying n.training"), align = c("c","l","l","c","c"), digits = c(0,0,0,6,6)), include.rownames = FALSE)


o1b.p.n.training <- ggplot(o1b.wse.dt.n.training) +
  geom_violin(aes(x=n.training,y=log(wse), fill = dataset), show.legend = T, draw_quantiles = c(0.5)) +
  scale_fill_manual(values = c("#00788C","#B07C4F","#76A165"), labels = c("combined", "test", "training"), name = "Fitting dataset:") +
  labs(x="Size of the training set", y="log(RWSE)") +
  #lims(y=c(0,0.5)) +
  theme_classic() +
  theme.settings + 
  theme(plot.margin = unit(c(0.7,0.7,0.7,0.7), "cm"))


#n.test

#o1b.wse.table.n.test <- o1b.wse.dt.n.test[, .(wse.median = quantile(wse, 0.5, na.rm = TRUE), wse.iqr = IQR(wse, na.rm = TRUE)), by = .(n.test, dataset)][, n.test := as.numeric(levels(n.test))[n.test]][dataset == "all", dataset := "combined"]
#setnames(o1b.wse.table.n.test, names(o1b.wse.table.n.test),c("n.test","dataset","median","IQR"))
#setorder(o1b.wse.table.n.test,n.test,dataset)
#print(xtable(o1b.wse.table.n.test, caption = c("numerical results for the Weibull-Cox scenario (ii) with varying n.test"), align = c("c","l","l","c","c"), digits = c(0,0,0,6,6)), include.rownames = FALSE)

o1b.p.n.test <- ggplot(o1b.wse.dt.n.test) +
  geom_violin(aes(x=n.test,y=log(wse), fill = dataset), show.legend = T, draw_quantiles = c(0.5)) +
  scale_fill_manual(values = c("#00788C","#B07C4F","#76A165"), labels = c("combined", "test", "training"), name = "Fitting dataset:") +
  labs(x="Size of the test set", y="log(RWSE)") +
  #lims(y=c(0,0.15)) +
  theme_classic() +
  theme.settings + 
  theme(plot.margin = unit(c(0.7,0.7,0.7,0.7), "cm"))


# censoring rate (Weibull)

#o1b.wse.table.offset <- o1b.wse.dt.offset[, .(wse.median = quantile(wse, 0.5, na.rm = TRUE), wse.iqr = IQR(wse, na.rm = TRUE)), by = .(offset, dataset)][, offset := as.numeric(levels(offset))[offset]][order(offset)][dataset == "all", dataset := "combined"]
#o1b.wse.table.offset[offset == "0.8", offset := "0.2"][offset == "-0.002", offset := "0.5"][offset == "-0.803", offset := "0.8"][, offset := as.factor(offset)]
#setnames(o1b.wse.table.offset, names(o1b.wse.table.offset),c("censoring rate","dataset","median","IQR"))
#setorder(o1b.wse.table.offset,`censoring rate`,dataset)
#print(xtable(o1b.wse.table.offset, caption = c("numerical results for the Weibull-Cox scenario (ii) with varying censoring rate"), align = c("c","l","l","c","c"), digits = c(0,0,0,6,6)), include.rownames = FALSE)

o1b.p.offset <- ggplot(o1b.wse.dt.offset[offset == "0.8", cens.rate := "0.2"][offset == "-0.002", cens.rate := "0.5"][offset == "-0.803", cens.rate := "0.8"][, cens.rate := as.factor(cens.rate)]) +
  geom_violin(aes(x=cens.rate,y=log(wse), fill = dataset), show.legend = T, draw_quantiles = c(0.5)) +
  scale_fill_manual(values = c("#00788C","#B07C4F","#76A165"), labels = c("combined", "test", "training"), name = "Fitting dataset:") +
  labs(x="Censoring rate", y="log(RWSE)") +
  #lims(y=c(0,0.1)) +
  theme_classic() +
  theme.settings + 
  theme(plot.margin = unit(c(0.7,0.7,0.7,0.7), "cm"))

plot_o1b <- ggarrange(o1b.p.n.per.ds, o1b.p.n.training, o1b.p.n.test, o1b.p.offset, common.legend = TRUE, legend = "bottom", labels = c("(a)","(b)","(c)","(d)"), hjust = 0, font.label = list(size = 20, face = "bold"))

ggsave(filename = "plot_o1b.eps", plot = plot_o1b, device = cairo_ps, path = plot.path, units = "in", width = 14, height = 8, dpi = 600)



##o2

theme.margins <- theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
                       axis.title.x = element_text(margin = margin(t = 15, r = 0, b = 0, l = 0)))

o2.wse.dt.n.per.ds.cox[, estimator := "cox"]
o2.wse.dt.n.per.ds.km[, estimator := "km"]

miss.wse.dt <- rbind(o2.wse.dt.n.per.ds.cox, o2.wse.dt.n.per.ds.km)
miss.wse.dt[, n := as.factor(2*as.numeric(as.character(n.per.ds)))]

#miss.wse.table <- miss.wse.dt[, .(wse.median = quantile(wse, 0.5, na.rm = TRUE), wse.iqr = IQR(wse, na.rm = TRUE)), by = .(n, dataset, estimator)][, n := as.numeric(levels(n))[n]][order(n,dataset)][dataset == "all", dataset := "combined"]
#setnames(miss.wse.table, names(miss.wse.table),c("n","dataset","estimator","median","IQR"))
#print(xtable(miss.wse.table, caption = c("numerical results for the misspecification scenario (iii)"), align = c("c","l","l","l","c","c"), digits = c(0,0,0,0,6,6)), include.rownames = FALSE)

plot_o2 <- ggplot(miss.wse.dt) +
  geom_violin(aes(x=n, y=log(wse), fill = dataset, color = estimator, linetype = estimator), show.legend = T, draw_quantiles = c(0.5)) +
  scale_fill_manual(values = c("#00788C","#B07C4F","#76A165"), labels = c("combined", "test", "training"), name = "Fitting dataset:") +
  scale_linetype_manual(values = c("solid", "dotted"), labels = c("Cox", "Kaplan-Meier"), name = "Estimator:") +
  scale_color_manual(values = c("#000000", "#ff0003"), labels = c("Cox", "Kaplan-Meier"), name = "Estimator:") +
  labs(x="Total size of the dataset", y="log(RWSE)", color = "Estimator:", linetype = "Estimator:") +
  #lims(y=c(0,0.6)) +
  theme_classic() +
  theme.settings +
  theme(legend.position = "bottom", legend.direction = "horizontal", legend.box = "vertical") +
  theme.margins

ggsave(filename = "plot_o2.eps", plot = plot_o2, device = cairo_ps, path = plot.path, units = "in", width = 14, height = 8, dpi = 600)



##o3a and o3b
theme.x.blank <- theme(axis.title.x = element_blank(),
                       axis.ticks.x = element_blank(),
                       axis.text.x = element_blank())

##o3a

#o3a.wse.table.estimator <- o3a.wse.dt.estimator[, .(wse.median = quantile(wse, 0.5, na.rm = TRUE), wse.iqr = IQR(wse, na.rm = TRUE)), by = .(estimator, dataset)][dataset == "all", dataset := "combined"]
#setnames(o3a.wse.table.estimator, names(o3a.wse.table.estimator),c("estimator","dataset","median","IQR"))
#print(xtable(o3a.wse.table.estimator, caption = c("numerical results for the high-noise ML scenario (v)"), align = c("c","l","l","c","c"), digits = c(0,0,0,6,6)), include.rownames = FALSE)

o3a.p.estimator.cox <- ggplot(o3a.wse.dt.estimator[estimator == "cox", estimator.type := "Cox"][estimator.type == "Cox"]) +
  geom_violin(aes(x=dataset,y=wse, fill = dataset), show.legend = T, draw_quantiles = c(0.5)) +
  scale_fill_manual(values = c("#00788C","#B07C4F","#76A165"), labels = c("combined", "test", "training"), name = "Fitting dataset:") +
  labs(x="Fitting dataset", y="RWSE", title="Cox") +
  lims(y=c(0,0.03)) +
  #lims(y=c(0,0.03)) +
  theme_classic() +
  theme.settings +
  theme.x.blank +
  theme(plot.title = element_text(size=rel(3))) + 
  theme(plot.margin = unit(c(0.9,0.9,0.9,0.9), "cm"))

o3a.p.estimator.glmnet <- ggplot(o3a.wse.dt.estimator[estimator == "glmnet", estimator.type := "Lasso"][estimator.type == "Lasso"]) +
  geom_violin(aes(x=dataset,y=wse, fill = dataset), show.legend = T, draw_quantiles = c(0.5)) +
  scale_fill_manual(values = c("#00788C","#B07C4F","#76A165"), labels = c("combined", "test", "training"), name = "Fitting dataset:") +
  labs(x="Fitting dataset", y="RWSE", title="Lasso") +
  lims(y=c(0,0.125)) +
  #lims(y=c(0,0.08)) +
  theme_classic() +
  theme.settings +
  theme.x.blank +
  theme(plot.title = element_text(size=rel(3))) + 
  theme(plot.margin = unit(c(0.9,0.9,0.9,0.9), "cm"))

o3a.p.estimator.rangersf <- ggplot(o3a.wse.dt.estimator[estimator == "rangersf", estimator.type := "Random forest"][estimator.type == "Random forest"]) +
  geom_violin(aes(x=dataset,y=wse, fill = dataset), show.legend = T, draw_quantiles = c(0.5)) +
  scale_fill_manual(values = c("#00788C","#B07C4F","#76A165"), labels = c("combined", "test", "training"), name = "Fitting dataset:") +
  labs(x="Fitting dataset", y="RWSE", title="Random forest") +
  lims(y=c(0,0.05)) +
  #lims(y=c(0,0.05)) +
  theme_classic() +
  theme.settings +
  theme.x.blank +
  theme(plot.title = element_text(size=rel(3))) + 
  theme(plot.margin = unit(c(0.9,0.9,0.9,0.9), "cm"))

o3a.p.estimator.xgboost <- ggplot(o3a.wse.dt.estimator[estimator == "xgboost", estimator.type := "XGBoost"][estimator.type == "XGBoost"]) +
  geom_violin(aes(x=dataset,y=wse, fill = dataset), show.legend = T, draw_quantiles = c(0.5)) +
  scale_fill_manual(values = c("#00788C","#B07C4F","#76A165"), labels = c("combined", "test", "training"), name = "Fitting dataset:") +
  labs(x=NULL, y="RWSE", title="XGBoost") +
  lims(y=c(0,0.4)) +
  #lims(y=c(0,0.06)) +
  theme_classic() +
  theme.settings +
  theme.x.blank +
  theme(plot.title = element_text(size=rel(3))) + 
  theme(plot.margin = unit(c(0.9,0.9,0.9,0.9), "cm"))

plot_o3a <- ggarrange(o3a.p.estimator.cox, o3a.p.estimator.glmnet, o3a.p.estimator.rangersf, o3a.p.estimator.xgboost, common.legend = TRUE, legend = "bottom", labels = c("(a)","(b)","(c)","(d)"), hjust = 0, font.label = list(size = 20, face = "bold"))

ggsave(filename = "plot_o3a.eps", plot = plot_o3a, device = cairo_ps, path = plot.path, units = "in", width = 14, height = 8, dpi = 600)



##o3b

#o3b.wse.table.estimator <- o3b.wse.dt.estimator[, .(wse.median = quantile(wse, 0.5, na.rm = TRUE), wse.iqr = IQR(wse, na.rm = TRUE)), by = .(estimator, dataset)][dataset == "all", dataset := "combined"]
#setnames(o3b.wse.table.estimator, names(o3b.wse.table.estimator),c("estimator","dataset","median","IQR"))
#print(xtable(o3b.wse.table.estimator, caption = c("numerical results for the high-noise ML scenario (v)"), align = c("c","l","l","c","c"), digits = c(0,0,0,6,6)), include.rownames = FALSE)

o3b.p.estimator.cox <- ggplot(o3b.wse.dt.estimator[estimator == "cox", estimator.type := "Cox"][estimator.type == "Cox"]) +
  geom_violin(aes(x=dataset,y=wse, fill = dataset), show.legend = T, draw_quantiles = c(0.5)) +
  scale_fill_manual(values = c("#00788C","#B07C4F","#76A165"), labels = c("combined", "test", "training"), name = "Fitting dataset:") +
  labs(x="Fitting dataset", y="RWSE", title="Cox") +
  #lims(y=c(0,0.03)) +
  lims(y=c(0,0.03)) +
  theme_classic() +
  theme.settings +
  theme.x.blank +
  theme(plot.title = element_text(size=rel(3))) + 
  theme(plot.margin = unit(c(0.9,0.9,0.9,0.9), "cm"))

o3b.p.estimator.glmnet <- ggplot(o3b.wse.dt.estimator[estimator == "glmnet", estimator.type := "Lasso"][estimator.type == "Lasso"]) +
  geom_violin(aes(x=dataset,y=wse, fill = dataset), show.legend = T, draw_quantiles = c(0.5)) +
  scale_fill_manual(values = c("#00788C","#B07C4F","#76A165"), labels = c("combined", "test", "training"), name = "Fitting dataset:") +
  labs(x="Fitting dataset", y="RWSE", title="Lasso") +
  #lims(y=c(0,0.125)) +
  lims(y=c(0,0.08)) +
  theme_classic() +
  theme.settings +
  theme.x.blank +
  theme(plot.title = element_text(size=rel(3))) + 
  theme(plot.margin = unit(c(0.9,0.9,0.9,0.9), "cm"))

o3b.p.estimator.rangersf <- ggplot(o3b.wse.dt.estimator[estimator == "rangersf", estimator.type := "Random forest"][estimator.type == "Random forest"]) +
  geom_violin(aes(x=dataset,y=wse, fill = dataset), show.legend = T, draw_quantiles = c(0.5)) +
  scale_fill_manual(values = c("#00788C","#B07C4F","#76A165"), labels = c("combined", "test", "training"), name = "Fitting dataset:") +
  labs(x="Fitting dataset", y="RWSE", title="Random forest") +
  #lims(y=c(0,0.05)) +
  lims(y=c(0,0.05)) +
  theme_classic() +
  theme.settings +
  theme.x.blank +
  theme(plot.title = element_text(size=rel(3))) + 
  theme(plot.margin = unit(c(0.9,0.9,0.9,0.9), "cm"))

o3b.p.estimator.xgboost <- ggplot(o3b.wse.dt.estimator[estimator == "xgboost", estimator.type := "XGBoost"][estimator.type == "XGBoost"]) +
  geom_violin(aes(x=dataset,y=wse, fill = dataset), show.legend = T, draw_quantiles = c(0.5)) +
  scale_fill_manual(values = c("#00788C","#B07C4F","#76A165"), labels = c("combined", "test", "training"), name = "Fitting dataset:") +
  labs(x=NULL, y="RWSE", title="XGBoost") +
  #lims(y=c(0,0.4)) +
  lims(y=c(0,0.06)) +
  theme_classic() +
  theme.settings +
  theme.x.blank +
  theme(plot.title = element_text(size=rel(3))) + 
  theme(plot.margin = unit(c(0.9,0.9,0.9,0.9), "cm"))

plot_o3b <- ggarrange(o3b.p.estimator.cox, o3b.p.estimator.glmnet, o3b.p.estimator.rangersf, o3b.p.estimator.xgboost, common.legend = TRUE, legend = "bottom", labels = c("(a)","(b)","(c)","(d)"), hjust = 0, font.label = list(size = 20, face = "bold"))

ggsave(filename = "plot_o3b.eps", plot = plot_o3b, device = cairo_ps, path = plot.path, units = "in", width = 14, height = 8, dpi = 600)
