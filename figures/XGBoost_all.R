library(readr)
library(tidyverse)
library(stringr)
library(Biostrings)
library(pROC)
library(ggplot2)
library(ggsci, help, pos = 2, lib.loc = NULL)
library(patchwork, help, pos = 2, lib.loc = NULL)
library("wesanderson")
ipsc <- matrix(data = 0, nrow = 5, ncol = 4)|> as.data.frame()
colnames(ipsc) = c("ipsc","k562","npc","genome_wide")
for(j in 1:5){
      gRNA <- read_csv(paste0("data/mhc/ipsc/ipsc-pfdr0.05-pfdr0.2-binary-fold", j, "-test.csv"))
      obs = gRNA$significant
      pred_ipsc <- read_csv(paste0("data/mhc/ipsc/result/ipsc-pfdr0.05-pfdr0.2-binary-fold", j, "-test.csv"),col_names=F)
      pred_k562 <- read_csv(paste0("data/mhc/ipsc/result/ipsc-by-k562-pfdr0.05-pfdr0.2-binary-fold", j, "-test.csv"),col_names=F)
      pred_npc <- read_csv(paste0("data/mhc/ipsc/result/ipsc-by-npc-pfdr0.05-pfdr0.2-binary-fold", j, "-test.csv"),col_names=F)
      pred_sc <- read_csv(paste0("data/mhc/ipsc/result/ipsc-by-sc-pfdr0.05-pfdr0.2-binary-fold", j, "-test.csv"),col_names=F)
      # protospacer <- gRNA$protospacer
      ipsc[j,] <- c(auc(obs,pred_ipsc$X1),auc(obs,pred_k562$X1),auc(obs,pred_npc$X1), auc(obs,pred_sc$X1))
}
summary(ipsc)
data_long <- gather(ipsc, model, auc, ipsc:genome_wide, factor_key=TRUE)
p1 <- data_long  |> ggplot(aes(x= model, y=auc, fill=model))+ geom_boxplot() + theme_minimal() + labs(title = "ipsc as test set")


k562 <- matrix(data = 0, nrow = 5, ncol = 4)|> as.data.frame()
colnames(k562) = c("k562","ipsc","npc","genome_wide")
for(j in 1:5){
      gRNA <- read_csv(paste0("data/mhc/k562/k562-pfdr0.05-pfdr0.2-binary-fold", j, "-test.csv"))
      obs = gRNA$significant
      pred_k562 <- read_csv(paste0("data/mhc/k562/result/k562-pfdr0.05-pfdr0.2-binary-fold", j, "-test.csv"),col_names=F)
      pred_ipsc <- read_csv(paste0("data/mhc/k562/result/k562-by-ipsc-pfdr0.05-pfdr0.2-binary-fold", j, "-test.csv"),col_names=F)
      pred_npc <- read_csv(paste0("data/mhc/k562/result/k562-by-npc-pfdr0.05-pfdr0.2-binary-fold", j, "-test.csv"),col_names=F)
      pred_sc <- read_csv(paste0("data/mhc/k562/result/k562-by-sc-pfdr0.05-pfdr0.2-binary-fold", j, "-test.csv"),col_names=F)
      # protospacer <- gRNA$protospacer
      k562[j,] <- c(auc(obs,pred_k562$X1),auc(obs,pred_ipsc$X1),auc(obs,pred_npc$X1), auc(obs,pred_sc$X1))
}
summary(k562)
data_long <- gather(k562, model, auc, k562:genome_wide, factor_key=TRUE)
p2 <- data_long  |> ggplot(aes(x= model, y=auc, fill=model))+ geom_boxplot() + theme_minimal() + labs(title = "k562 as test set")

npc <- matrix(data = 0, nrow = 5, ncol = 4)|> as.data.frame()
colnames(npc) = c("npc","k562","ipsc","genome_wide")
for(j in 1:5){
      gRNA <- read_csv(paste0("data/mhc/npc/npc-pfdr0.05-pfdr0.2-binary-fold", j, "-test.csv"))
      obs = gRNA$significant
      pred_npc <- read_csv(paste0("data/mhc/npc/result/npc-pfdr0.05-pfdr0.2-binary-fold", j, "-test.csv"),col_names=F)
      pred_k562 <- read_csv(paste0("data/mhc/npc/result/npc-by-k562-pfdr0.05-pfdr0.2-binary-fold", j, "-test.csv"),col_names=F)
      pred_ipsc <- read_csv(paste0("data/mhc/npc/result/npc-by-ipsc-pfdr0.05-pfdr0.2-binary-fold", j, "-test.csv"),col_names=F)
      pred_sc <- read_csv(paste0("data/mhc/npc/result/npc-by-sc-pfdr0.05-pfdr0.2-binary-fold", j, "-test.csv"),col_names=F)
      # protospacer <- gRNA$protospacer
      npc[j,] <- c(auc(obs,pred_npc$X1),auc(obs,pred_k562$X1),auc(obs,pred_npc$X1), auc(obs,pred_sc$X1))
}
summary(npc)
data_long <- gather(npc, model, auc, npc:genome_wide, factor_key=TRUE)
p3 <- data_long  |> ggplot(aes(x= model, y=auc, fill=model))+ geom_boxplot() + theme_minimal() + labs(title = "npc as test set")


library(patchwork, help, pos = 2, lib.loc = NULL)
png(file="plots/cross-apply.png", height = 500, width = 1000)
p1 + p2 + p3
dev.off()

auc_seq_anno =c()
auc_seq =c()
auc_anno =c()
auc_top_seqanno = c()
for(j in 1:5){
      gRNA <- read_csv(paste0("data/dat_discovery/wgCERES-gRNAs-k562-discovery-screen-enh_baseMean125-binary-", j, "-test.csv"))
      obs = gRNA$significance
      pred_seq_anno <- read_csv(paste0("data/dat_discovery/result/wgCERES-gRNAs-k562-discovery-screen-enh_baseMean125-binary-", j, "-test.csv"),col_names=F)
      pred_top_seqanno <- read_csv(paste0("data/dat_discovery/result/wgCERES-gRNAs-k562-discovery-screen-enh_baseMean125-binary-", j, "-test-top-seqanno.csv"),col_names=F)
      pred_seq <- read_csv(paste0("data/dat_discovery/result/wgCERES-gRNAs-k562-discovery-screen-enh_baseMean125-binary-", j, "-test-seq.csv"),col_names=F)
      pred_anno <- read_csv(paste0("data/dat_discovery/result/wgCERES-gRNAs-k562-discovery-screen-enh_baseMean125-binary-", j, "-test-anno.csv"),col_names=F)
      auc_seq_anno[j] = auc(obs,pred_seq_anno$X1)  
      auc_top_seqanno[j] = auc(obs,pred_top_seqanno$X1)  
      auc_seq[j] = auc(obs,pred_seq$X1) 
      auc_anno[j] = auc(obs,pred_anno$X1) 
}
enhancer <- data.frame(seq_anno = auc_seq_anno, seq_top_anno = auc_top_seqanno, seq = auc_seq, anno = auc_anno, row.names = NULL, check.rows = FALSE, check.names = TRUE, stringsAsFactors = default.stringsAsFactors())
p1 <- enhancer  |> pivot_longer(everything(), names_to = "inputs", values_to = "auc")  |> 
      ggplot(aes(x = inputs, y= auc))+
      geom_boxplot()

for(j in 1:5){
      gRNA <- read_csv(paste0("data/dat_discovery/wgCERES-gRNAs-k562-discovery-screen-pro_baseMean125-binary-", j, "-test.csv"))
      obs = gRNA$significance
      pred_seq_anno <- read_csv(paste0("data/dat_discovery/result/wgCERES-gRNAs-k562-discovery-screen-pro_baseMean125-binary-", j, "-test.csv"),col_names=F)
      pred_top_seqanno <- read_csv(paste0("data/dat_discovery/result/wgCERES-gRNAs-k562-discovery-screen-pro_baseMean125-binary-", j, "-test-top-seqanno2.csv"),col_names=F)
      pred_seq <- read_csv(paste0("data/dat_discovery/result/wgCERES-gRNAs-k562-discovery-screen-pro_baseMean125-binary-", j, "-test-seq.csv"),col_names=F)
      pred_anno <- read_csv(paste0("data/dat_discovery/result/wgCERES-gRNAs-k562-discovery-screen-pro_baseMean125-binary-", j, "-test-anno.csv"),col_names=F)
      auc_seq_anno[j] = auc(obs,pred_seq_anno$X1)  
      auc_top_seqanno[j] = auc(obs,pred_top_seqanno$X1)  
      auc_seq[j] = auc(obs,pred_seq$X1) 
      auc_anno[j] = auc(obs,pred_anno$X1) 
}
promoter <- data.frame(seq_anno = auc_seq_anno, seq_top_anno = auc_top_seqanno, seq = auc_seq, anno = auc_anno, row.names = NULL, check.rows = FALSE, check.names = TRUE, stringsAsFactors = default.stringsAsFactors())
p2 <- promoter  |> pivot_longer(everything(), names_to = "inputs", values_to = "auc")  |> 
      ggplot(aes(x = inputs, y= auc))+
      geom_boxplot()
library(patchwork)
pdf("plots/cell_fitness.pdf",height = 185/25.4, width = 180/25.4)
p1 + p2
dev.off()

for(j in 1:5){
      gRNA <- read_csv(paste0("data/scOct4th/sc_fdr_filtered-binary-", j, "-test.csv"))
      obs = gRNA$significance
      pred_seq_anno <- read_csv(paste0("data/scOct4th/result/sc_fdr_filtered-binary-", j, "-test.csv"),col_names=F)
      pred_seq <- read_csv(paste0("data/scOct4th/result/sc_fdr_filtered-binary-", j, "-test-seq.csv"),col_names=F)
      pred_anno <- read_csv(paste0("data/scOct4th/result/sc_fdr_filtered-binary-", j, "-test-anno.csv"),col_names=F)
      # protospacer <- gRNA$protospacer
      auc_seq_anno[j] = auc(obs,pred_seq_anno$X1)   
      auc_seq[j] = auc(obs,pred_seq$X1) 
      auc_anno[j] = auc(obs,pred_anno$X1) 
}

## single cell- top prediction density
library(ensembldb)
library(EnsDb.Hsapiens.v86) # this pkg is about 75 Mb
edb <- EnsDb.Hsapiens.v86
g <- genes(edb)
g <- keepStandardChromosomes(g, pruning.mode="coarse")
lens <- seqlengths(seqinfo(g))
names(lens) <- paste0("chr",c(1,10:19,2,20:22,3:9,"MT","X","Y"))
for(j in 1:5){
  gRNA <- read_csv(paste0("data/scOct4th/sc_fdr_filtered-binary-", j, "-test.csv"))
  obs = gRNA$significance
  pred_seq_anno <- read_csv(paste0("data/scOct4th/result/sc_fdr_filtered-binary-", j, "-test.csv"),col_names=F)
  cmb <- data.frame(gRNA,pred = pred_seq_anno$X1, row.names = NULL, check.rows = FALSE, check.names = TRUE, stringsAsFactors = default.stringsAsFactors())
  cmb_gr <- makeGRangesFromDataFrame(cmb, keep.extra.columns=TRUE, seqnames.field = "grna.chrom",
                                     start.field = "grna.start", end.field = "grna.end")
  query <- GenomicRanges::tileGenome(lens[unique(cmb$dhs_chrom)],
                                     tilewidth = 500000,
                                     cut.last.tile.in.chrom = TRUE
  )
  end(query) = end(query) + 500000
  query$id = seq_along(1:length(query))
  inner_rng <- join_overlap_inner(query, cmb_gr) %>% 
    group_by(id) %>% 
    mutate(top1 = ifelse(pred == max(pred),1,0) |> as.factor(), n = n(),
           top2= ifelse(pred %in% sort(pred, na.last = TRUE, decreasing = TRUE)[1:2],1,0) |> as.factor())
  cmb_filter<- mcols(inner_rng) %>% as.data.frame()
  hist(cmb_filter$n)
  wilcox.test(cmb_filter$pval_fdr_corrected[cmb_filter$top1==1],cmb_filter$pval_fdr_corrected[cmb_filter$top1==0],alternative = "less")
  wilcox.test(cmb_filter$pval_fdr_corrected[cmb_filter$top2==1],cmb_filter$pval_fdr_corrected[cmb_filter$top2==0],alternative = "less")
  wilcox.test(abs(cmb_filter$avg_logFC[cmb_filter$top1==1]),abs(cmb_filter$avg_logFC[cmb_filter$top1==0]),alternative = "greater")
  p1 <- ggplot(cmb_filter)+
    geom_density(aes(pval_fdr_corrected,fill=top1,col=top1),alpha=0.4) +
    labs(x = "p_fdr_corrected")+
    theme_bw()+
    scale_fill_igv()+
    scale_color_igv()
  p2 <- ggplot(mcols(inner_rng) %>% as.data.frame())+
    geom_density(aes(avg_logFC,fill=top1,col=top1),alpha=0.4) +
    labs(x = "avg_lfc")+
    theme_bw()+
    scale_fill_igv()+
    scale_color_igv()
  p1+p2
}

## WTcount
cor =c()
for(j in 1:5){
      gRNA <- read_csv(paste0("/proj/milovelab/mu/dukeproj/data/dat_WTcount/Maria-gRNAs-k562-validation-screen-full-test-fold",j,".csv"))
      obs = gRNA$avgWTcount_K562
      pred <- read_csv(paste0("/proj/milovelab/mu/dukeproj/data/dat_WTcount/result/Maria-gRNAs-k562-validation-screen-full-test-fold",j,".csv"),col_names=F)
      # protospacer <- gRNA$protospacer
      cor[j] = cor(obs,pred$X1, method = "spearman")   
}

## check ranking in cell fitness
cor_list = list()
mypal = pal_npg("nrc", alpha = 0.8)(9)
for(i in 1:5){
      test = read_csv(paste0("data/dat_discovery/wgCERES-gRNAs-k562-discovery-screen-enh_baseMean125-binary-", i, "-test-forrank.csv"))
      pred_xgb <- read_csv(paste0("/proj/milovelab/mu/dukeproj/data/dat_discovery/enhancer/result/wgCERES-gRNAs-k562-discovery-screen-enh_baseMean125-binary-",i,"-test-top-seqanno.csv"),col_names=F)
      pred_cnn <- read_csv(paste0("/proj/yunligrp/users/tianyou/gRNA/result/binary_fivefold/gRNA_binary-enh-BCE-seq-topannot-fold",i,"-Mar11.csv"),col_names=T)
      cmb <- data.frame(test,XGBoost = pred_xgb$X1, CNN = pred_cnn$predict, fold = i, row.names = NULL, check.rows = FALSE, check.names = TRUE, stringsAsFactors = default.stringsAsFactors())
      # cmb_filter <- cmb  |> filter(nsig1==1)  |> 
      #                   select(DHS, padj, pred) |> 
      #                   group_by(DHS) |> mutate(rank_obs = dense_rank(padj),rank_pred = dense_rank(pred))
      cmb_filter <- cmb  |> filter(nsig1==1)  |> 
                        select(DHS, log2FoldChange, padj, XGBoost,CNN, fold) |> 
                        group_by(DHS) |> mutate(
                              top1_XGBoost = ifelse(XGBoost == max(XGBoost),1,0) |> as.factor(),
                              top1_CNN = ifelse(CNN == max(CNN),1,0) |> as.factor(),
                              top2_XGBoost = ifelse(XGBoost %in% sort(XGBoost, na.last = TRUE, decreasing = TRUE)[1:2],1,0) |> as.factor(),
                              top2_CNN = ifelse(CNN %in% sort(CNN, na.last = TRUE, decreasing = TRUE)[1:2],1,0) |> as.factor()) 
      cor_list[[i]] = cmb_filter
}         
cmb_filter <- do.call(rbind, cor_list)  
             
wilcox.test(cmb_filter$padj[cmb_filter$top1_XGBoost==1],cmb_filter$padj[cmb_filter$top1_XGBoost==0],alternative = "less")
summary(cmb_filter$padj[cmb_filter$top1_XGBoost==1])
summary(cmb_filter$padj[cmb_filter$top1_XGBoost==0])
wilcox.test(cmb_filter$padj[cmb_filter$top2_XGBoost==1],cmb_filter$padj[cmb_filter$top2_XGBoost==0],alternative = "less")
summary(cmb_filter$padj[cmb_filter$top2_XGBoost==1])
summary(cmb_filter$padj[cmb_filter$top2_XGBoost==0])
wilcox.test(cmb_filter$padj[cmb_filter$top1_CNN==1],cmb_filter$padj[cmb_filter$top1_CNN==0],alternative = "less")
summary(cmb_filter$padj[cmb_filter$top1_CNN==1])
summary(cmb_filter$padj[cmb_filter$top1_CNN==0])
wilcox.test(cmb_filter$padj[cmb_filter$top2_CNN==1],cmb_filter$padj[cmb_filter$top2_CNN==0],alternative = "less")
summary(cmb_filter$padj[cmb_filter$top2_CNN==1])
summary(cmb_filter$padj[cmb_filter$top2_CNN==0])
wilcox.test(abs(cmb_filter$log2FoldChange[cmb_filter$top1_XGBoost==1]),abs(cmb_filter$log2FoldChange[cmb_filter$top1_XGBoost==0]),alternative = "greater")
summary(abs(cmb_filter$log2FoldChange[cmb_filter$top1_XGBoost==1]))
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.00006 0.01921 0.05502 0.30194 0.47010 3.80941 
summary(abs(cmb_filter$log2FoldChange[cmb_filter$top1_XGBoost==0]))
#     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.000002 0.011316 0.025637 0.081822 0.053211 4.271780 
wilcox.test(abs(cmb_filter$log2FoldChange[cmb_filter$top2_XGBoost==1]),abs(cmb_filter$log2FoldChange[cmb_filter$top2_XGBoost==0]),alternative = "greater")
wilcox.test(abs(cmb_filter$log2FoldChange[cmb_filter$top1_CNN==1]),abs(cmb_filter$log2FoldChange[cmb_filter$top1_CNN==0]),alternative = "greater")
summary(abs(cmb_filter$log2FoldChange[cmb_filter$top1_CNN==1]))
#     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.000011 0.020419 0.064314 0.320838 0.486941 3.809409 
summary(abs(cmb_filter$log2FoldChange[cmb_filter$top1_CNN==0]))
#  Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.000002 0.011227 0.025545 0.079574 0.052574 4.271780 
wilcox.test(abs(cmb_filter$log2FoldChange[cmb_filter$top2_CNN==1]),abs(cmb_filter$log2FoldChange[cmb_filter$top2_CNN==0]),alternative = "greater")
            

p1 <- ggplot(cmb_filter)+
      geom_density(aes(padj,fill=top1_XGBoost,col=top1_XGBoost),alpha=0.8) +
      labs(x = "adjusted p-value", fill = "Top 1 (XGBoost)", col = "Top 1 (XGBoost)")+
      scale_fill_manual(values = wes_palette("Moonrise3"), labels = c("Non-top 1", "Top 1"))+
      scale_colour_manual(values = wes_palette("Moonrise3"), labels = c("Non-top 1", "Top 1"))+
      theme_classic()+
      theme(panel.background = element_rect(fill = "white", colour = "grey50"),
            legend.position = c(.05, .99),
            text = element_text(size = 20),
            legend.text = element_text(size = 24),
            legend.justification = c("left", "top"),
            legend.box.just = "right",
            legend.margin = margin(6, 6, 6, 6))
p3 <- ggplot(cmb_filter)+
      geom_density(aes(padj,fill=top1_CNN,col=top1_CNN),alpha=0.8) +
      labs(x = "adjusted p-value", fill = "Top 1 (CNN)", col = "Top 1 (CNN)")+
      scale_fill_manual(values = wes_palette("GrandBudapest2"), labels = c("Non-top 1", "Top 1"))+
      scale_colour_manual(values = wes_palette("GrandBudapest2"), labels = c("Non-top 1", "Top 1"))+
      theme_classic()+
      theme(panel.background = element_rect(fill = "white", colour = "grey50"),
            legend.position = c(.05, .99),
            text = element_text(size = 20),
            legend.text = element_text(size = 24),
            legend.justification = c("left", "top"),
            legend.box.just = "right",
            legend.margin = margin(6, 6, 6, 6))
p2 <- ggplot(cmb_filter)+
      geom_density(aes(log2FoldChange,fill=top1_XGBoost,col=top1_XGBoost),alpha=0.8) +
      labs(x = "log2FC", fill = "Top 1 (XGBoost)", col = "Top 1 (XGBoost)")+
      scale_fill_manual(values = wes_palette("Moonrise3"), labels = c("Non-top 1", "Top 1"))+
      scale_colour_manual(values = wes_palette("Moonrise3"), labels = c("Non-top 1", "Top 1"))+
      theme_classic()+
      theme(panel.background = element_rect(fill = "white", colour = "grey50"),
            legend.position = c(.05, .99),
            text = element_text(size = 20),
            legend.text = element_text(size = 24),
            legend.justification = c("left", "top"),
            legend.box.just = "right",
            legend.margin = margin(6, 6, 6, 6))
p4 <- ggplot(cmb_filter)+
      geom_density(aes(log2FoldChange,fill=top1_CNN,col=top1_CNN),alpha=0.8) +
      labs(x = "log2FC",  fill = "Top 1 (CNN)", col = "Top 1 (CNN)")+
      scale_fill_manual(values = wes_palette("GrandBudapest2"), labels = c("Non-top 1", "Top 1"))+
      scale_colour_manual(values = wes_palette("GrandBudapest2"), labels = c("Non-top 1", "Top 1"))+
      theme_classic()+
      theme(panel.background = element_rect(fill = "white", colour = "grey50"),
            legend.position = c(.05, .99),
            text = element_text(size = 20),
            legend.text = element_text(size = 24),
            legend.justification = c("left", "top"),
            legend.box.just = "right",
            legend.margin = margin(6, 6, 6, 6))

jpeg(file="/proj/milovelab/mu/dukeproj/data/dat_discovery/result/plots/top1_comparison_enh.jpg",width = 10, height = 10,units = "in",res=450)
(p3 +p1) /(p4+p2) & theme(text = element_text(size = 28),
      plot.tag = element_text(size = 24)
      )
dev.off()
(p1 + p2) / (p3 +p4) + plot_annotation(tag_levels = 'A') &
theme(text = element_text(size = 24),
      plot.tag = element_text(size = 16)
      )

p5 <- ggplot(cmb_filter)+
      geom_density(aes(padj,fill=top2_XGBoost,col=top2_XGBoost),alpha=0.8) +
      labs(x = "adjusted p-value", fill = "Top 2 (XGBoost)", col = "Top 2 (XGBoost)")+
      scale_fill_manual(values = wes_palette("Moonrise3"), labels = c("Non-top 2", "Top 2"))+
      scale_colour_manual(values = wes_palette("Moonrise3"), labels = c("Non-top 2", "Top 2"))+
      theme_classic()+
      theme(panel.background = element_rect(fill = "white", colour = "grey50"),
            legend.position = c(.05, .99),
            text = element_text(size = 20),
            legend.text = element_text(size = 24),
            legend.justification = c("left", "top"),
            legend.box.just = "right",
            legend.margin = margin(6, 6, 6, 6))
p6 <- ggplot(cmb_filter)+
      geom_density(aes(padj,fill=top2_CNN,col=top2_CNN),alpha=0.8) +
      labs(x = "adjusted p-value", fill = "Top 2 (CNN)", col = "Top 2 (CNN)")+
      scale_fill_manual(values = wes_palette("GrandBudapest2"), labels = c("Non-top 2", "Top 2"))+
      scale_colour_manual(values = wes_palette("GrandBudapest2"), labels = c("Non-top 2", "Top 2"))+
      theme_classic()+
      theme(panel.background = element_rect(fill = "white", colour = "grey50"),
            legend.position = c(.05, .99),
            text = element_text(size = 20),
            legend.text = element_text(size = 24),
            legend.justification = c("left", "top"),
            legend.box.just = "right",
            legend.margin = margin(6, 6, 6, 6))
p7 <- ggplot(cmb_filter)+
      geom_density(aes(log2FoldChange,fill=top2_XGBoost,col=top2_XGBoost),alpha=0.8) +
      labs(x = "log2FC", fill = "Top 2 (XGBoost)", col = "Top 2 (XGBoost)")+
      scale_fill_manual(values = wes_palette("Moonrise3"), labels = c("Non-top 2", "Top 2"))+
      scale_colour_manual(values = wes_palette("Moonrise3"), labels = c("Non-top 2", "Top 2"))+
      theme_classic()+
      theme(panel.background = element_rect(fill = "white", colour = "grey50"),
            legend.position = c(.05, .99),
            text = element_text(size = 20),
            legend.text = element_text(size = 24),
            legend.justification = c("left", "top"),
            legend.box.just = "right",
            legend.margin = margin(6, 6, 6, 6))
p8 <- ggplot(cmb_filter)+
      geom_density(aes(log2FoldChange,fill=top2_CNN,col=top2_CNN),alpha=0.8) +
      labs(x = "log2FC",  fill = "Top 2 (CNN)", col = "Top 2 (CNN)")+
      scale_fill_manual(values = wes_palette("GrandBudapest2"), labels = c("Non-top 2", "Top 2"))+
      scale_colour_manual(values = wes_palette("GrandBudapest2"), labels = c("Non-top 2", "Top 2"))+
      theme_classic()+
      theme(panel.background = element_rect(fill = "white", colour = "grey50"),
            legend.position = c(.05, .99),
            text = element_text(size = 20),
            legend.text = element_text(size = 24),
            legend.justification = c("left", "top"),
            legend.box.just = "right",
            legend.margin = margin(6, 6, 6, 6))
jpeg(file="/proj/milovelab/mu/dukeproj/data/dat_discovery/result/plots/top2_comparison_enh.jpg",width = 10, height = 10,units = "in",res=450)
(p6 +p5) /(p8+p7) & theme(text = element_text(size = 28),
      plot.tag = element_text(size = 24)
      )
dev.off()
# cor <- cmb_filter |> group_by(DHS) |> summarize(cor = cor(rank_obs,rank_pred, method = "spearman"))                  

cor_data <- do.call(rbind, cor_list)
boxplot(cor_data$cor)

cor_list_pro = list()  
for(i in 1:5){
      test = read_csv(paste0("data/dat_discovery/wgCERES-gRNAs-k562-discovery-screen-pro_baseMean125-binary-", i, "-test-forrank.csv"))
      # pred_xgb <- read_csv(paste0("/proj/milovelab/mu/dukeproj/data/dat_discovery/result/wgCERES-gRNAs-k562-discovery-screen-pro_baseMean125-binary-",i,"-test-top-seqanno.csv"),col_names=F)
      pred_xgb <- read_csv(paste0("/proj/milovelab/mu/dukeproj/data/dat_discovery/promoter/result/wgCERES-gRNAs-k562-discovery-screen-pro_baseMean125-binary-",i,"-test-top-seqanno.csv"),col_names=F)
      pred_cnn <- read_csv(paste0("/proj/yunligrp/users/tianyou/gRNA/result/binary_fivefold/gRNA_binary-pro-BCE-seq-topannot-fold",i,"-Mar11.csv"),col_names=T)
      cmb <- data.frame(test,XGBoost = pred_xgb$X1, CNN = pred_cnn$predict, fold=i, row.names = NULL, check.rows = FALSE, check.names = TRUE, stringsAsFactors = default.stringsAsFactors())
      # cmb_filter <- cmb  |> filter(nsig1==1)  |> 
      #                   select(DHS, padj, pred) |> 
      #                   group_by(DHS) |> mutate(rank_obs = dense_rank(padj),rank_pred = dense_rank(pred))
      cmb_filter <- cmb  |> filter(nsig1==1)  |> 
                        select(DHS, log2FoldChange, padj, XGBoost,CNN, fold) |> 
                        group_by(DHS) |> mutate(
                              top1_XGBoost = ifelse(XGBoost == max(XGBoost),1,0) |> as.factor(),
                              top1_CNN = ifelse(CNN == max(CNN),1,0) |> as.factor(),
                              top2_XGBoost = ifelse(XGBoost %in% sort(XGBoost, na.last = TRUE, decreasing = TRUE)[1:2],1,0) |> as.factor(),
                              top2_CNN = ifelse(CNN %in% sort(CNN, na.last = TRUE, decreasing = TRUE)[1:2],1,0) |> as.factor()) 
      cor_list_pro[[i]] = cmb_filter
}
cmb_filter <- do.call(rbind, cor_list_pro) 
wilcox.test(cmb_filter$padj[cmb_filter$top1_XGBoost==1],cmb_filter$padj[cmb_filter$top1_XGBoost==0],alternative = "less")
summary(cmb_filter$padj[cmb_filter$top1_XGBoost==1])
#     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.000000 0.001899 0.657808 0.508284 0.905512 0.999560 
summary(cmb_filter$padj[cmb_filter$top1_XGBoost==0])
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  0.0000  0.5749  0.8349  0.6956  0.9436  1.0000 
wilcox.test(cmb_filter$padj[cmb_filter$top2_XGBoost==1],cmb_filter$padj[cmb_filter$top2_XGBoost==0],alternative = "less")
summary(cmb_filter$padj[cmb_filter$top2_XGBoost==1])
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.00000 0.01401 0.72916 0.55704 0.91746 0.99986 
summary(cmb_filter$padj[cmb_filter$top2_XGBoost==0])
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  0.0000  0.5984  0.8411  0.7075  0.9453  1.0000 
wilcox.test(cmb_filter$padj[cmb_filter$top1_CNN==1],cmb_filter$padj[cmb_filter$top1_CNN==0],alternative = "less")
summary(cmb_filter$padj[cmb_filter$top1_CNN==1])
#     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.000000 0.001873 0.662709 0.507113 0.902329 0.999870 
summary(cmb_filter$padj[cmb_filter$top1_CNN==0])
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  0.0000  0.5761  0.8349  0.6957  0.9437  1.0000 
wilcox.test(cmb_filter$padj[cmb_filter$top2_CNN==1],cmb_filter$padj[cmb_filter$top2_CNN==0],alternative = "less")
summary(cmb_filter$padj[cmb_filter$top2_CNN==1])
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.00000 0.01522 0.73211 0.55884 0.91867 0.99994 
summary(cmb_filter$padj[cmb_filter$top2_CNN==0])
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  0.0000  0.5974  0.8404  0.7071  0.9452  1.0000 
wilcox.test(abs(cmb_filter$log2FoldChange[cmb_filter$top1_XGBoost==1]),abs(cmb_filter$log2FoldChange[cmb_filter$top1_XGBoost==0]),alternative = "greater")
summary(abs(cmb_filter$log2FoldChange[cmb_filter$top1_XGBoost==1]))
#     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.000198 0.020052 0.058351 0.365309 0.519787 4.278335 
summary(abs(cmb_filter$log2FoldChange[cmb_filter$top1_XGBoost==0]))
#     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.000004 0.013041 0.030441 0.159827 0.073895 3.790058
wilcox.test(abs(cmb_filter$log2FoldChange[cmb_filter$top2_XGBoost==1]),abs(cmb_filter$log2FoldChange[cmb_filter$top2_XGBoost==0]),alternative = "greater")
wilcox.test(abs(cmb_filter$log2FoldChange[cmb_filter$top1_CNN==1]),abs(cmb_filter$log2FoldChange[cmb_filter$top1_CNN==0]),alternative = "greater")
summary(abs(cmb_filter$log2FoldChange[cmb_filter$top1_CNN==1]))
#     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.000041 0.020675 0.060346 0.356805 0.521759 4.278335 
summary(abs(cmb_filter$log2FoldChange[cmb_filter$top1_CNN==0]))
#     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.000004 0.013031 0.030400 0.160828 0.073660 3.790058 
wilcox.test(abs(cmb_filter$log2FoldChange[cmb_filter$top2_CNN==1]),abs(cmb_filter$log2FoldChange[cmb_filter$top2_CNN==0]),alternative = "greater")
    
p1 <- ggplot(cmb_filter)+
      geom_density(aes(padj,fill=top1_XGBoost,col=top1_XGBoost),alpha=0.8) +
      labs(x = "adjusted p-value", fill = "Top 1 (XGBoost)", col = "Top 1 (XGBoost)")+
      scale_fill_manual(values = wes_palette("Moonrise3"), labels = c("Non-top 1", "Top 1"))+
      scale_colour_manual(values = wes_palette("Moonrise3"), labels = c("Non-top 1", "Top 1"))+
      theme_classic()+
      theme(panel.background = element_rect(fill = "white", colour = "grey50"),
            legend.position = c(.05, .99),
            text = element_text(size = 20),
            legend.text = element_text(size = 24),
            legend.justification = c("left", "top"),
            legend.box.just = "right",
            legend.margin = margin(6, 6, 6, 6))
p3 <- ggplot(cmb_filter)+
      geom_density(aes(padj,fill=top1_CNN,col=top1_CNN),alpha=0.8) +
      labs(x = "adjusted p-value", fill = "Top 1 (CNN)", col = "Top 1 (CNN)")+
      scale_fill_manual(values = wes_palette("GrandBudapest2"), labels = c("Non-top 1", "Top 1"))+
      scale_colour_manual(values = wes_palette("GrandBudapest2"), labels = c("Non-top 1", "Top 1"))+
      theme_classic()+
      theme(panel.background = element_rect(fill = "white", colour = "grey50"),
            legend.position = c(.05, .99),
            text = element_text(size = 20),
            legend.text = element_text(size = 24),
            legend.justification = c("left", "top"),
            legend.box.just = "right",
            legend.margin = margin(6, 6, 6, 6))
p2 <- ggplot(cmb_filter)+
      geom_density(aes(log2FoldChange,fill=top1_XGBoost,col=top1_XGBoost),alpha=0.8) +
      labs(x = "log2FC", fill = "Top 1 (XGBoost)", col = "Top 1 (XGBoost)")+
      scale_fill_manual(values = wes_palette("Moonrise3"), labels = c("Non-top 1", "Top 1"))+
      scale_colour_manual(values = wes_palette("Moonrise3"), labels = c("Non-top 1", "Top 1"))+
      theme_classic()+
      theme(panel.background = element_rect(fill = "white", colour = "grey50"),
            legend.position = c(.05, .99),
            text = element_text(size = 20),
            legend.text = element_text(size = 24),
            legend.justification = c("left", "top"),
            legend.box.just = "right",
            legend.margin = margin(6, 6, 6, 6))
p4 <- ggplot(cmb_filter)+
      geom_density(aes(log2FoldChange,fill=top1_CNN,col=top1_CNN),alpha=0.8) +
      labs(x = "log2FC",  fill = "Top 1 (CNN)", col = "Top 1 (CNN)")+
      scale_fill_manual(values = wes_palette("GrandBudapest2"), labels = c("Non-top 1", "Top 1"))+
      scale_colour_manual(values = wes_palette("GrandBudapest2"), labels = c("Non-top 1", "Top 1"))+
      theme_classic()+
      theme(panel.background = element_rect(fill = "white", colour = "grey50"),
            legend.position = c(.05, .99),
            text = element_text(size = 20),
            legend.text = element_text(size = 24),
            legend.justification = c("left", "top"),
            legend.box.just = "right",
            legend.margin = margin(6, 6, 6, 6))

jpeg(file="/proj/milovelab/mu/dukeproj/data/dat_discovery/result/plots/top1_comparison_pro.jpg",width = 10, height = 10,units = "in",res=450)
(p3 +p1) /(p4+p2) & theme(text = element_text(size = 28),
      plot.tag = element_text(size = 24)
      )
dev.off()
(p1 + p2) / (p3 +p4) + plot_annotation(tag_levels = 'A') &
theme(text = element_text(size = 24),
      plot.tag = element_text(size = 16)
      )

p5 <- ggplot(cmb_filter)+
      geom_density(aes(padj,fill=top2_XGBoost,col=top2_XGBoost),alpha=0.8) +
      labs(x = "adjusted p-value", fill = "Top 2 (XGBoost)", col = "Top 2 (XGBoost)")+
      scale_fill_manual(values = wes_palette("Moonrise3"), labels = c("Non-top 2", "Top 2"))+
      scale_colour_manual(values = wes_palette("Moonrise3"), labels = c("Non-top 2", "Top 2"))+
      theme_classic()+
      theme(panel.background = element_rect(fill = "white", colour = "grey50"),
            legend.position = c(.05, .99),
            text = element_text(size = 20),
            legend.text = element_text(size = 24),
            legend.justification = c("left", "top"),
            legend.box.just = "right",
            legend.margin = margin(6, 6, 6, 6))
p6 <- ggplot(cmb_filter)+
      geom_density(aes(padj,fill=top2_CNN,col=top2_CNN),alpha=0.8) +
      labs(x = "adjusted p-value", fill = "Top 2 (CNN)", col = "Top 2 (CNN)")+
      scale_fill_manual(values = wes_palette("GrandBudapest2"), labels = c("Non-top 2", "Top 2"))+
      scale_colour_manual(values = wes_palette("GrandBudapest2"), labels = c("Non-top 2", "Top 2"))+
      theme_classic()+
      theme(panel.background = element_rect(fill = "white", colour = "grey50"),
            legend.position = c(.05, .99),
            text = element_text(size = 20),
            legend.text = element_text(size = 24),
            legend.justification = c("left", "top"),
            legend.box.just = "right",
            legend.margin = margin(6, 6, 6, 6))
p7 <- ggplot(cmb_filter)+
      geom_density(aes(log2FoldChange,fill=top2_XGBoost,col=top2_XGBoost),alpha=0.8) +
      labs(x = "log2FC", fill = "Top 2 (XGBoost)", col = "Top 2 (XGBoost)")+
      scale_fill_manual(values = wes_palette("Moonrise3"), labels = c("Non-top 2", "Top 2"))+
      scale_colour_manual(values = wes_palette("Moonrise3"), labels = c("Non-top 2", "Top 2"))+
      theme_classic()+
      theme(panel.background = element_rect(fill = "white", colour = "grey50"),
            legend.position = c(.05, .99),
            text = element_text(size = 20),
            legend.text = element_text(size = 24),
            legend.justification = c("left", "top"),
            legend.box.just = "right",
            legend.margin = margin(6, 6, 6, 6))
p8 <- ggplot(cmb_filter)+
      geom_density(aes(log2FoldChange,fill=top2_CNN,col=top2_CNN),alpha=0.8) +
      labs(x = "log2FC",  fill = "Top 2 (CNN)", col = "Top 2 (CNN)")+
      scale_fill_manual(values = wes_palette("GrandBudapest2"), labels = c("Non-top 2", "Top 2"))+
      scale_colour_manual(values = wes_palette("GrandBudapest2"), labels = c("Non-top 2", "Top 2"))+
      theme_classic()+
      theme(panel.background = element_rect(fill = "white", colour = "grey50"),
            legend.position = c(.05, .99),
            text = element_text(size = 20),
            legend.text = element_text(size = 24),
            legend.justification = c("left", "top"),
            legend.box.just = "right",
            legend.margin = margin(6, 6, 6, 6))
jpeg(file="/proj/milovelab/mu/dukeproj/data/dat_discovery/result/plots/top2_comparison_pro.jpg",width = 10, height = 10,units = "in",res=450)
(p6 +p5) /(p8+p7) & theme(text = element_text(size = 28),
      plot.tag = element_text(size = 24)
      )
dev.off()