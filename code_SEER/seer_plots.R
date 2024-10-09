library(data.table)
library(ggplot2)
library(extrafont)

loadfonts(device = c("postscript"))

path.data <- "..." 
path.plots <- "..."

seer_res_data <- readRDS(file = paste0(path.data,"/seer_plot_data_list.rds")) #path to "seer_plot_data_list.rds" data file

for(l in 1:length(seer_res_data)){
  assign(names(seer_res_data)[l],seer_res_data[[l]])
}

### individual plots
xgboost.dt.long <- rbind(xgboost.dt.long.nt, xgboost.dt.long.t)
xgboost.dt.long[, method := "XGBoost"]
rangersf.dt.long <- rbind(rangersf.dt.long.nt, rangersf.dt.long.t)
rangersf.dt.long[, method := "Random Forest"]

theme.settings <- theme(line = element_line(linewidth = 1),
                        text = element_text(size = 10, family = "serif"),
                        #text = element_text(size = 10),
                        rect = element_rect(linewidth = 1),
                        plot.title = element_text(size = rel(3.5), hjust = 0.5),
                        legend.title = element_text(size = rel(2.75)),
                        legend.text = element_text(size = rel(2.75)),
                        legend.position = "bottom",
                        axis.title = element_text(size = rel(2.75)),
                        axis.text = element_text(size = rel(2.75)),
                        axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
                        axis.title.x = element_text(margin = margin(t = 15, r = 0, b = 0, l = 0)),
                        #theme(axis.text.y = element_text(hjust = 0)),
                        strip.text = element_text(size = rel(2.75)),
                        strip.background = element_blank(),
                        panel.grid.major.y = element_line(color = "grey", linewidth = 0.2))


##Cox PH
plot_cox <- ggplot(cox.dt.long.nt) +
  #geom_line(aes(x = times, y = bs, group = sample, linetype = dataset, color = dataset), alpha = 0.45, size = 0.2, show.legend = FALSE) +
  #geom_line(aes(x = times, y = bs.mean, linetype = dataset, color = dataset), linewidth = 1.3) +
  #geom_ribbon(aes(x = times, ymin = bs.mean-bs.sd, ymax = bs.mean+bs.sd, group = sample, fill = dataset), alpha = 0.02, show.legend = FALSE) +
  geom_line(aes(x = times, y = bs.mean, linetype = dataset, color = dataset, linewidth = dataset)) +
  scale_linetype_manual(name = "Dataset:", values = c("solid","dashed","dotted"), labels = c("combined","test","training")) +
  scale_color_manual(name = "Dataset:", values = c("#00788C","#B07C4F","#76A165"), labels = c("combined","test","training")) +
  scale_fill_manual(name = "Dataset:", values = c("#00788C","#B07C4F","#76A165"), labels = c("combined","test","training")) +
  scale_linewidth_manual(name = "Dataset:", values = c(0.7, 1.1, 1.7), labels = c("combined","test","training")) +
  lims(x = c(0,114), y = c(0,0.10)) +
  labs(x = "Time [months]", y = "Brier Score", title = "Cox proportional hazards model") + # with title
  #labs(x = "Time", y = "Brier score") + # no title
  theme_classic() +
  theme.settings

ggsave(filename = "SEER_cox.eps", plot = plot_cox, device = cairo_ps, path = path.plots, units = "in", width = 14, height = 8, dpi = 600)


## xgboost
plot_xgboost <- ggplot(xgboost.dt.long) +
  #geom_line(aes(x = times, y = bs, group = sample, linetype = dataset, color = dataset), alpha = 0.45, size = 0.2, show.legend = FALSE) +
  #geom_line(aes(x = times, y = bs.mean, linetype = dataset, color = dataset), size = 1.3) +
  #geom_ribbon(aes(x = times, ymin = bs.mean-bs.sd, ymax = bs.mean+bs.sd, group = sample, fill = dataset), alpha = 0.02, show.legend = FALSE) +
  geom_line(aes(x = times, y = bs.mean, linetype = dataset, color = dataset, linewidth = dataset)) +
  scale_linetype_manual(name = "Dataset:", values = c("solid","dashed","dotted"), labels = c("combined","test","training")) +
  scale_color_manual(name = "Dataset:", values = c("#00788C","#B07C4F","#76A165"), labels = c("combined","test","training")) +
  scale_fill_manual(name = "Dataset:", values = c("#00788C","#B07C4F","#76A165"), labels = c("combined","test","training")) +
  scale_linewidth_manual(name = "Dataset:", values = c(0.7, 1.1, 1.7), labels = c("combined","test","training")) +
  lims(x = c(0,114), y = c(0,0.15)) +
  labs(x = "Time [months]", y = "Brier score", title = "XGBoost") +
  facet_wrap(~ tuned, labeller = 
               labeller(
                 tuned = ~ paste("Tuned: ", .),
                 .multi_line = FALSE
               )
  ) +
  theme_classic() +
  theme.settings

ggsave(filename = "SEER_xgboost.eps", plot = plot_xgboost, device = cairo_ps, path = path.plots, units = "in", width = 14, height = 8, dpi = 600)


##random forest
plot_random_forest <- ggplot(rangersf.dt.long) +
  #geom_line(aes(x = times, y = bs, group = sample, linetype = dataset, color = dataset), alpha = 0.45, size = 0.2, show.legend = FALSE) +
  #geom_line(aes(x = times, y = bs.mean, linetype = dataset, color = dataset), size = 1.3) +
  #geom_ribbon(aes(x = times, ymin = bs.mean-bs.sd, ymax = bs.mean+bs.sd, group = sample, fill = dataset), alpha = 0.02, show.legend = FALSE) +
  geom_line(aes(x = times, y = bs.mean, linetype = dataset, color = dataset, linewidth = dataset)) +
  scale_linetype_manual(name = "Dataset:", values = c("solid","dashed","dotted"), labels = c("combined","test","training")) +
  scale_color_manual(name = "Dataset:", values = c("#00788C","#B07C4F","#76A165"), labels = c("combined","test","training")) +
  scale_fill_manual(name = "Dataset:", values = c("#00788C","#B07C4F","#76A165"), labels = c("combined","test","training")) +
  scale_linewidth_manual(name = "Dataset:", values = c(0.7, 1.1, 1.7), labels = c("combined","test","training")) +
  lims(x = c(0,114), y = c(0,0.15)) +
  labs(x = "Time [months]", y = "Brier score", title = "Random forest") +
  facet_wrap(~ tuned, labeller = 
               labeller(
                 tuned = ~ paste("Tuned: ", .),
                 .multi_line = FALSE
               )
  ) +
  theme_classic() +
  theme.settings

ggsave(filename = "SEER_random_forest.eps", plot = plot_random_forest, device = cairo_ps, path = path.plots, units = "in", width = 14, height = 8, dpi = 600)
