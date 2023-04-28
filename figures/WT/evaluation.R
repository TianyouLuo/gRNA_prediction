library(tidyverse)
library(data.table)

CNN_dir = "/proj/yunligrp/users/tianyou/gRNA/result/WT_count"
xgboost_dir = "/proj/milovelab/mu/dukeproj/data/dat_WTcount/result/"
savedir = "/proj/yunligrp/users/tianyou/gRNA/result/WT_count/comparison"
result_mat = matrix(NA, nrow = 5, ncol = 9)
all_pred = list()

for (fold in 1:5){
  CNN_result0 = fread(file.path(CNN_dir, paste0("Maria_plasmid_full-MSE-seq-fold", 
                                               fold, "-Dec04.csv"))) %>% select(-V1)
  CNN_result1 = fread(file.path(CNN_dir, paste0("Maria_plasmid_full-MSE-seq-no5fil-fold", 
                                               fold, "-Feb28.csv"))) %>%
    rename(CNN_no5_predict = predict) %>%
    select(-V1)
  CNN_result2 = fread(file.path(CNN_dir, paste0("Maria_plasmid_full-MSE-seq-no7fil-fold", 
                                               fold, "-Feb28.csv"))) %>%
    rename(CNN_no7_predict = predict) %>%
    select(-V1)
  CNN_result3 = fread(file.path(CNN_dir, paste0("Maria_plasmid_full-MSE-seq-no3fil-fold", 
                                                fold, "-Feb28.csv"))) %>%
    rename(CNN_no3_predict = predict) %>%
    select(-V1)
  CNN_result4 = fread(file.path(CNN_dir, paste0("Maria_plasmid_full-MSE-seq-no57fil-fold", 
                                                fold, "-Feb28.csv"))) %>%
    rename(CNN_no57_predict = predict) %>%
    select(-V1)
  CNN_result5 = fread(file.path(CNN_dir, paste0("Maria_plasmid_full-MSE-seq-no37fil-fold", 
                                                fold, "-Feb28.csv"))) %>%
    rename(CNN_no37_predict = predict) %>%
    select(-V1)
  CNN_result6 = fread(file.path(CNN_dir, paste0("Maria_plasmid_full-MSE-seq-no35fil-fold", 
                                                fold, "-Feb28.csv"))) %>%
    rename(CNN_no35_predict = predict) %>%
    select(-V1)
  CNN_result7 = fread(file.path(CNN_dir, paste0("Maria_plasmid_full-MSE-seq-no357fil-fold", 
                                                fold, "-Feb28.csv"))) %>%
    rename(CNN_no357_predict = predict) %>%
    select(-V1)
  xgb_result = fread(file.path(xgboost_dir, 
                               paste0("Maria-gRNAs-k562-validation-screen-full-test-fold",
                                      fold,".csv")),
                     col.names = c("xgboost_predict"))
  CNN_result = CNN_result0 %>%
    full_join(CNN_result1, by = c("grna", "true")) %>%
    full_join(CNN_result2, by = c("grna", "true")) %>%
    full_join(CNN_result3, by = c("grna", "true")) %>%
    full_join(CNN_result4, by = c("grna", "true")) %>%
    full_join(CNN_result5, by = c("grna", "true")) %>%
    full_join(CNN_result6, by = c("grna", "true")) %>%
    full_join(CNN_result7, by = c("grna", "true"))
  if (dim(xgb_result)[1] != dim(CNN_result)[1]){
    warning(paste0("Model trained on fold ", fold, " has inconsistent dimensions."))
  }
  
  test_result = cbind(CNN_result, xgb_result)
  result_mat[fold, 1] = cor(test_result$true, test_result$predict, method = "spearman")
  result_mat[fold, 2] = cor(test_result$true, test_result$xgboost_predict, method = "spearman")
  result_mat[fold, 3] = cor(test_result$true, test_result$CNN_no5_predict, method = "spearman")
  result_mat[fold, 4] = cor(test_result$true, test_result$CNN_no7_predict, method = "spearman")
  result_mat[fold, 5] = cor(test_result$true, test_result$CNN_no3_predict, method = "spearman")
  result_mat[fold, 6] = cor(test_result$true, test_result$CNN_no57_predict, method = "spearman")
  result_mat[fold, 7] = cor(test_result$true, test_result$CNN_no37_predict, method = "spearman")
  result_mat[fold, 8] = cor(test_result$true, test_result$CNN_no35_predict, method = "spearman")
  result_mat[fold, 9] = cor(test_result$true, test_result$CNN_no357_predict, method = "spearman")
  all_pred[[fold]] = test_result
}

result_mat = as.data.table(result_mat)
colnames(result_mat) = c("CNN","XGBoost", "CNN_no_filter5", "CNN_no_filter7", "CNN_no_filter3",
                         "CNN_no_filter57", "CNN_no_filter37", "CNN_no_filter35", "CNN_no_filter357")
write_csv(result_mat, file.path(savedir, "WT_spearman_combined_summary.csv"))

all_pred = bind_rows(all_pred)
write_csv(all_pred, file.path(savedir, "WT_predicted_combined.csv"))

