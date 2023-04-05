library(tidyverse)
library(data.table)
library(reticulate)
library(viridis)
library(abind)

CNN_dir = "/proj/yunligrp/users/tianyou/gRNA/result/WT_count/shap"
dat_dir = "/proj/yunligrp/users/tianyou/gRNA/data/WT_fivefold/"

np <- import("numpy")
shap_mean_table = array(dim = c(5, 4, 20))
shap_table = list()

preprocess_seq = function(seq){
  length = nchar(seq[1])
  DATA_X = array(NA, dim = c(length(seq), 4, length))
  print(dim(DATA_X))
  for (l in 1:length(seq)) {
    for (i in 1:length) {
      if (i <= nchar(seq[l])) {
        if (substring(seq[l], i, i) %in% c("A", "a")) {
          DATA_X[l, 1, i] = 1
        } else if (substring(seq[l], i, i) %in% c("C", "c")) {
          DATA_X[l, 2, i] = 1
        } else if (substring(seq[l], i, i) %in% c("G", "g")) {
          DATA_X[l, 3, i] = 1
        } else if (substring(seq[l], i, i) %in% c("T", "t")) {
          DATA_X[l, 4, i] = 1
        } else {
          stop("Non-ATGC character " + substring(seq[l], i, i))
        }
      }
    }
  }
  return(DATA_X)
}

plot_mean_sd = function(mean_dat, sd_dat){
  sd_mat_long = as_tibble(sd_dat, rownames = "BP") %>%
    pivot_longer(`1`:`20`, names_to = "position", values_to = "sd")
  
  fig = as_tibble(mean_dat, rownames = "BP") %>%
    pivot_longer(`1`:`20`, names_to = "position", values_to = "mean") %>%
    left_join(sd_mat_long, by = c("BP", "position")) %>%
    mutate(position = as.factor(as.numeric(position))) %>%
    ggplot(aes(x = position, y = mean, fill = BP)) +
    geom_bar(stat = "identity", position=position_dodge(width = 0.75)) +
    geom_errorbar(aes(x = position, ymin = mean - sd, ymax = mean + sd), 
                  position=position_dodge(width = 0.75), width = 0.5) +
    scale_fill_viridis(discrete = TRUE, alpha=0.6) +
    theme(panel.background = element_rect(color=NA, fill = "white"),
          panel.border = element_blank(),
          panel.grid.major = element_blank(),
          axis.line = element_line(size = 1, linetype = "solid", colour = "black"),
          axis.text = element_text(size = 20),
          text = element_text(size = 24),
          legend.text = element_text(size = 24),
          # legend.background = element_rect(color = NA, fill = "aliceblue"),
          legend.position = c(0.99,0.01),
          legend.justification = c("right", "bottom"),
          legend.title = element_text(vjust = -6, hjust = 0.05)) + 
    labs(y = "Mean SHAP value", x = "Position", fill = "") +
    guides(fill = guide_legend(ncol = 2))
  
  return(fig)
}

plot_mean = function(mean_dat){
  fig = as_tibble(mean_dat, rownames = "BP") %>%
    pivot_longer(`1`:`20`, names_to = "position", values_to = "mean") %>%
    # left_join(sd_mat_long, by = c("BP", "position")) %>%
    mutate(position = as.factor(as.numeric(position))) %>%
    ggplot(aes(x = position, y = mean, fill = BP)) +
    geom_bar(stat = "identity", position=position_dodge(width = 0.75)) +
    # geom_errorbar(aes(x = position, ymin = mean - sd, ymax = mean + sd), 
    #               position=position_dodge(width = 0.75), width = 0.5) +
    scale_fill_viridis(discrete = TRUE, alpha=0.6) +
    theme(panel.background = element_rect(color=NA, fill = "white"),
          panel.border = element_blank(),
          panel.grid.major = element_blank(),
          axis.line = element_line(size = 1, linetype = "solid", colour = "black"),
          axis.text = element_text(size = 20),
          text = element_text(size = 24),
          legend.text = element_text(size = 24),
          # legend.background = element_rect(color = NA, fill = "aliceblue"),
          legend.position = c(0.99,0.01),
          legend.justification = c("right", "bottom"),
          legend.title = element_text(vjust = -6, hjust = 0.05)) + 
    labs(y = "Mean SHAP value", x = "Position", fill = "") +
    guides(fill = guide_legend(ncol = 2))
  
  return(fig)
}

for (fold in 1:5){
  dat = read_csv(file.path(dat_dir, paste0("Maria-gRNAs-k562-validation-screen-full-test-fold",
                                           fold, ".csv")))
  seq = dat$protospacer
  seq_onehot = preprocess_seq(seq)
  shap = np$load(file.path(CNN_dir, paste0("WT_continuous_seqshap-fold", fold, "-test-Dec04.npy")))
  
  shap_mean = apply(shap * seq_onehot, c(2,3), mean, na.rm = TRUE)
  shap_sd = apply(shap * seq_onehot, c(2,3), sd, na.rm = TRUE)
  
  colnames(shap_mean) = c(1:20)
  rownames(shap_mean) = c("A","C","G","T")
  colnames(shap_sd) = c(1:20)
  rownames(shap_sd) = c("A","C","G","T")
  shap_mean_table[fold,,] = shap_mean
  shap_table[[fold]] = shap * seq_onehot

  shap_fig = plot_mean(shap_mean)
  # shap_fig = plot_mean_sd(shap_mean, shap_sd)
  ggsave(file.path(CNN_dir, paste0("WT_seq_meanSHAP_fold",fold,".png")), shap_fig, width = 16, height = 8)
}

fold_mean = apply(shap_mean_table, c(2,3), mean)
fold_sd = apply(shap_mean_table, c(2,3), sd)

colnames(fold_mean) = c(1:20)
rownames(fold_mean) = c("A","C","G","T")
colnames(fold_sd) = c(1:20)
rownames(fold_sd) = c("A","C","G","T")

fold_fig = plot_mean_sd(fold_mean, fold_sd)
ggsave(file.path(CNN_dir, paste0("WT_seq_meanSHAP_allfolds.png")), fold_fig, width = 16, height = 8)


shap_merged = abind(shap_table, along = 1)
all_mean = apply(shap_merged, c(2,3), mean, na.rm = TRUE)
all_sd = apply(shap_merged, c(2,3), sd, na.rm = TRUE)
colnames(all_mean) = c(1:20)
rownames(all_mean) = c("A","C","G","T")
colnames(all_sd) = c(1:20)
rownames(all_sd) = c("A","C","G","T")
all_fig = plot_mean_sd(all_mean, all_sd)
ggsave(file.path(CNN_dir, paste0("WT_seq_SHAP_allmerged.png")), all_fig, width = 16, height = 8)

