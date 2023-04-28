library(tidyverse)
library(viridis)
library(ggpubr)
library(rstatix)

dir = "/proj/yunligrp/users/tianyou/gRNA/result/WT_count/comparison"

auc = read_csv(file.path(dir, "WT_spearman_combined_summary.csv"))
auc_long = auc %>% 
  mutate(fold = row_number()) %>%
  pivot_longer(cols = CNN:CNN_no_filter357, 
               names_to = "Method", values_to = "Spearman Correlation")

## boxplot ####
# WT_pval_all = auc_long %>%
#   t_test(`Spearman Correlation` ~ Method, paired = TRUE,
#          ref.group = "CNN") %>%
#   mutate(p.adj.combined = ifelse(p.adj.signif == "ns", "n.s.", 
#                                  paste0(format(p.adj, digits = 2), p.adj.signif)))
xlabs = c("CNN", "XGBoost", "CNN_no7", "CNN_no5", "CNN_no3", "CNN_no57", "CNN_no37", "CNN_no35", "CNN_no357")
means = aggregate(`Spearman Correlation` ~  Method, auc_long, mean)

fig0 = auc_long %>%
  ggplot(aes(x = factor(Method, levels = c("CNN", "XGBoost", "CNN_no_filter7",
                                           "CNN_no_filter5", "CNN_no_filter3",
                                           "CNN_no_filter57", "CNN_no_filter37", 
                                           "CNN_no_filter35", "CNN_no_filter357")), y = `Spearman Correlation`)) +
  geom_boxplot(aes(fill = Method)) +
  geom_point(aes(group = Method), size=4, alpha=0.9) +
  geom_text(data = means, aes(label = round(`Spearman Correlation`, digits = 3), 
                              y = `Spearman Correlation` + 0.04), size = 7) +
  scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  coord_cartesian(ylim = c(0.5,1)) +
  labs(x = "Method") + 
  scale_x_discrete(labels = xlabs) +
  # stat_pvalue_manual(WT_pval_all %>% add_xy_position(x="Method",dodge=0.75), 
  #                    label = "p.adj.combined", tip.length = 0, size = 6,
  #                    step.increase = 0.01) +
  theme(panel.background = element_rect(color=NA, fill = "white"),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        axis.line = element_line(size = 1, linetype = "solid", colour = "black"),
        text = element_text(size = 22),
        legend.position="none") +
  guides(fill = guide_legend(ncol = 2))
  
ggsave(file.path(dir, "WT_boxplot_full.png"), fig0, width = 18, height = 8)


## smaller one with only four comparisons ####
WT_pval = auc_long %>%
  filter(Method %in% c("CNN", "XGBoost", "CNN_no_filter57", "CNN_no_filter357")) %>%
  t_test(`Spearman Correlation` ~ Method, paired = TRUE,
         ref.group = "CNN") %>%
  mutate(p.adj.combined = ifelse(p.adj.signif == "ns", "n.s.",
                                 paste0(format(p.adj, digits = 2), p.adj.signif)))

fig1 = auc_long %>%
  filter(Method %in% c("CNN", "XGBoost", "CNN_no_filter57", "CNN_no_filter357")) %>%
  ggplot(aes(x = factor(Method, levels = c("CNN", "XGBoost", "CNN_no_filter57", "CNN_no_filter357")), 
             y = `Spearman Correlation`)) +
  geom_boxplot(aes(fill = Method)) +
  geom_point(aes(group = Method), size=4, alpha=0.9) +
  scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  coord_cartesian(ylim = c(0.5,1)) +
  labs(x = "Method") + 
  scale_x_discrete(labels = c("CNN", "XGBoost", "CNN_no57", "CNN_no357")) +
  stat_pvalue_manual(WT_pval %>% add_xy_position(x="Method",dodge=0.75),
                     label = "p.adj.combined", tip.length = 0.02, size = 6,
                     step.increase = 0.03) +
  theme(panel.background = element_rect(color=NA, fill = "white"),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        axis.line = element_line(size = 1, linetype = "solid", colour = "black"),
        text = element_text(size = 22),
        legend.position="none") +
  guides(fill = guide_legend(ncol = 2))

ggsave(file.path(dir, "WT_boxplot_four.png"), fig1, width = 8, height = 8)



## scatter plot ####
pred = read_csv(file.path(dir, "WT_predicted_combined.csv")) %>%
  select(grna, true, predict, xgboost_predict)
colnames(pred) = c("grna", "true", "CNN", "XGBoost")
fig = pred %>%
  select(-grna) %>%
  pivot_longer(cols = CNN:XGBoost, names_to = "Method", values_to = "predicted") %>%
  ggplot(aes(x = true, y = predicted)) +
  geom_point(aes(color = Method)) +
  geom_abline(slope = 1, intercept = 0, color = "darkslategray", size = 1.5) +
  scale_x_continuous(expand = c(0, 0), limits = c(0, 7000)) + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, 7000)) +
  scale_color_manual(values = c("blueviolet", "darkorange"),
                     breaks = c("CNN" ,"XGBoost")) +
  theme(panel.background = element_rect(color=NA, fill = "white"),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        axis.line = element_line(size = 1, linetype = "solid", colour = "black"),
        text = element_text(size = 24),
        legend.position = c(1,0.005),
        legend.justification = c("right", "bottom"),
        legend.text = element_text(size = 24),
        legend.background = element_rect(color = NA, fill = "aliceblue")) +
  labs(x = "True Counts in Wild-type Cells",
       y = "Predicted Counts in Wild-type Cells") +
  guides(colour = guide_legend(override.aes = list(size=8)))
ggsave("WT_scatter.png", fig, width = 8, height = 8)


