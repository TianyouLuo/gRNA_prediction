library(tidyverse)
library(viridis)
library(ggpubr)
library(rstatix)

dir = "/proj/yunligrp/users/tianyou/gRNA/result/mhc"
auc = read_csv(file.path(dir, "mhc_auc_all_summary.csv"))
colnames(auc) = c("method", "training_model", "K562", "NPC", "iPSC")

####### Figures separated by training models ######
# auc_long = auc %>% 
#   pivot_longer(cols = K562:iPSC, names_to = "test_data", values_to = "AUC")
# model_list = c("K562", "NPC", "iPSC", "genome-wide")
# fig_list = list()
# for (model in model_list){
#   fig = auc_long %>%
#     filter(training_model == model) %>%
#     ggplot(aes(x = test_data, y = AUC, fill = method)) +
#     geom_boxplot() +
#     geom_point(aes(group = method), size=4, alpha=0.9, position=position_dodge(width = 0.75)) +
#     scale_fill_viridis(discrete = TRUE, alpha=0.6) +
#     coord_cartesian(ylim = c(0.5, 0.85)) +
#     labs(x = "Test Data", title = paste0("Models trained on ", model)) +
#     theme(panel.background = element_rect(color=NA, fill = "white"),
#           panel.border = element_blank(),
#           panel.grid.major = element_blank(),
#           axis.line = element_line(size = 1, linetype = "solid", colour = "black"),
#           text = element_text(size = 20),
#           legend.text = element_text(size = 24),
#           legend.background = element_rect(color = NA, fill = "aliceblue"),
#           plot.title = element_text(vjust = -6, hjust = 0.05))
#   ggsave(file.path(dir, paste0(model, "_model.png")), fig, width = 8, height = 8)
#   fig_list[[model]] = fig
# }
# 
# fig_all = ggarrange(fig_list$K562, fig_list$NPC, fig_list$iPSC, ncol = 3,
#           common.legend = TRUE, legend = "bottom")
# ggsave(file.path(dir, "all_models.png"), fig_all, width = 18, height = 9)



####### Figures separated by testing set ######
auc_long = auc %>% 
  mutate(training_model = ifelse(training_model == "sc", "genome-wide", training_model)) %>%
  pivot_longer(cols = K562:iPSC, names_to = "test_data", values_to = "AUC")

test_list = c("K562", "NPC", "iPSC")
model_list = c("K562", "NPC", "iPSC", "genome-wide")
fig_list = list()
for (test in test_list){
  auc_pval_all = auc_long %>%
    group_by(method) %>%
    filter(test_data == test) %>%
    t_test(AUC ~ training_model, paired = TRUE,
           ref.group = test, p.adjust.method = "BH") %>%
    filter(p.adj <= 0.05) %>%
    mutate(p.adj.combined = ifelse(p.adj.signif == "ns", "n.s.", 
                                   paste0(format(p.adj, digits = 2), p.adj.signif)))
  
  fig = auc_long %>%
    filter(test_data == test) %>%
    ggplot(aes(x = training_model, y = AUC, fill = method)) +
    geom_boxplot() +
    geom_point(aes(group = method), size=4, alpha=0.9, position=position_dodge(width = 0.75)) +
    scale_fill_viridis(discrete = TRUE, alpha=0.6) +
    stat_pvalue_manual(auc_pval_all %>% add_xy_position(x="training_model",group = "method",dodge=0.75), 
                       label = "p.adj.combined", tip.length = 0.03, size = 6,
                       step.increase = 0.1, step.group.by = "method") +
    coord_cartesian(ylim = c(0.5, 0.85)) +
    labs(x = "Training model", title = paste0(test," as test set")) +
    theme(panel.background = element_rect(color=NA, fill = "white"),
          panel.border = element_blank(),
          panel.grid.major = element_blank(),
          axis.line = element_line(size = 1, linetype = "solid", colour = "black"),
          text = element_text(size = 20),
          legend.text = element_text(size = 24),
          legend.background = element_rect(color = NA, fill = "aliceblue"),
          plot.title = element_text(vjust = -6, hjust = 0.05)) +
    guides(fill = guide_legend(override.aes = list(shape = NA)))  ## remove the dot in legend
  ggsave(file.path(dir, paste0(test, "_as_test_set.png")), fig, width = 8, height = 8)
  fig_list[[test]] = fig
}

fig_all = ggarrange(fig_list$K562, fig_list$NPC, fig_list$iPSC, ncol = 3,
                    common.legend = TRUE, legend = "bottom")
ggsave(file.path(dir, "all_models.png"), fig_all, width = 18, height = 9)
