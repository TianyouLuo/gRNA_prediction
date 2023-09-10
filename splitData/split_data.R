library(tidyverse)
require(rtracklayer)
# extraCols_narrowPeak <- c(singnalValue = "numeric", pValue = "numeric",
#                           qValue = "numeric", peak = "integer")
# gr_narrowPeak <- import("C:/Users/wancenmu/OneDrive - University of North Carolina at Chapel Hill/Lab/Yun project/duke_result/K562-results-julien/gRNA_prediction/data_new/Hichip/k562_H3K27ac_encode_replicated_peaks_narrowpeaks_ENCFF044JNJ.bed", format = "BED",
#                         extraCols = extraCols_narrowPeak)
hic <- read_tsv("/proj/milovelab/mu/dukeproj/K562-results-julien/gRNA_prediction/K562_combined_r1r2r3.10k.2.peaks.bedpe")
library(GenomicRanges)
library(InteractionSet)
# Converting data.frame to GInteraction
convertToGI <- function(df){
  row.regions <- GRanges(df$chr1, IRanges(df$start1,df$end1))# interaction start
  col.regions <- GRanges(df$chr2, IRanges(df$start2,df$end2))# interaction end
  gi <- GInteractions(row.regions, col.regions)
  mcols(gi) <- df[,7:14] # Interaction frequencies
  return(gi)
}

hic.gi <- convertToGI(hic)
hic.gi$log10fdr <- -log10(hic.gi$fdr)
hic.gi$log10fdr[is.infinite(hic.gi$log10fdr)] <- max(hic.gi$log10fdr[is.finite(hic.gi$log10fdr)])
library(readr)
library(plyranges)
setwd("/proj/milovelab/mu/dukeproj/data/dat_April_resplit/")
## add H3K27ac to unique table1
dat<- read_table2("/proj/milovelab/mu/dukeproj/K562-supp-tables/screen_hg19_uniqe.txt")

dat1 = read_csv("/proj/milovelab/mu/dukeproj/K562-supp-tables/discovery_screen_k562_sgrna_deseq2_results_hg19_new.csv", guess_max = 10000)
dat1 = dat1[!duplicated(dat1$protospacer),]
dat<-dat %>% left_join(dat1[,7:ncol(dat1)],by="protospacer")

# dat = read_csv("./data_new/wgCERES-gRNAs-k562-discovery-screen-pro_rawp0.05.csv")
# colnames(dat)
gr <- makeGRangesFromDataFrame(dat, keep.extra.columns=TRUE,strand.field = "strand.y",
                                start.field="chromStart",end.field="chromEnd")
gr <- makeGRangesFromDataFrame(dat, keep.extra.columns=TRUE,strand.field = "strand",
                                start.field="start",end.field="end")

olap <- findOverlaps(gr,hic.gi) %>% as.data.frame()
olap$log10fdr <- hic.gi$log10fdr[olap$subjectHits]
sumlogfdr <- olap %>% group_by(queryHits) %>% summarise(sumlog10fdr=sum(log10fdr))

countolap <- countOverlaps(gr,hic.gi)
ref = read_tsv("/proj/milovelab/mu/dukeproj/K562-results-julien/gRNA_prediction/gene_info.gencode.v28lift37.txt")
refgr <- makeGRangesFromDataFrame(ref, keep.extra.columns=TRUE,
                               start.field="Promoter_Start",end.field="Promoter_End")

out <- linkOverlaps(hic.gi, gr, refgr) %>% as.data.frame()
out <- unique(out[,c(1,2)])
out$log10fdr <- hic.gi$log10fdr[out$query]
prom <- out %>% dplyr::count(subject1)
promfdr <- out %>% group_by(subject1) %>% summarise(sumlog10fdr=sum(log10fdr))


dat$promnumber <- rep(0,nrow(dat))
dat$promnumber[prom$subject1]<-prom$n
dat$promlog10fdr <- rep(0,nrow(dat))
dat$promlog10fdr[promfdr$subject1]<-promfdr$sumlog10fdr
dat$olapnumber <- countolap
dat$sumlog10fdr <- rep(0,nrow(dat))
dat$sumlog10fdr[sumlogfdr$queryHits]<-sumlogfdr$sumlog10fdr

## deltagb and deltagh
protospacer <- dat$protospacer
strand = dat$strand.y
position2 <- matrix(NA,ncol=19,nrow=length(protospacer)) %>% as.data.frame
for(i in 1:19){
  position2[,i] <- as.factor(substr(protospacer,i,i+1))
}
colnames(position2) <- paste("position",seq_len(19),sep = "_")

delta <- c(-1,-2.1,-1.8,-0.9,-0.9,-2.1,-1.7,-0.9,-1.3,-2.7,-2.9,-1.1,-0.6,-1.5,-1.6,-0.2)
names(delta)<-c("TT","GT","CT","AT","TG","GG","CG","AG","TC","GC","CC","AC","TA","GA","CA","AA")
weight= c(1.80, 1.96, 1.90, 2.13, 1.38, 1.46, 1.00, 1.39, 1.51, 1.98, 1.88, 1.72, 2.02, 1.93, 2.08, 1.94, 2.15, 2.04, 2.25)
deltagh<- sapply(1:length(protospacer), function(i){sum(weight * delta[position2[i,] %>% unname() %>% unlist()] )})
dat$deltagh <-deltagh

rna_rna <- c(-0.93,-2.24,-2.08,-1.1,-2.11,-3.26,-2.36,-2.08,-2.35,-3.42,-3.26,-2.24,-1.33,-2.35,-2.11,-0.93)
dna_dna <- c(-1,-1.44,-1.28,-0.88,-1.45,-1.84,-2.17,-1.28,-1.3,-2.24,-1.84,-1.44,-0.58,-1.3,-1.45,-1)
protospacer2 <- ifelse(strand=="+",protospacer,chartr('ATGC', 'TACG', protospacer))
protospacer2 <- ifelse(strand=="+",protospacer2,reverse(protospacer2))

position3 <- matrix(NA,ncol=19,nrow=length(protospacer2)) %>% as.data.frame
for(i in 1:19){
  position3[,i] <- as.factor(substr(protospacer2,i,i+1))
}
colnames(position3) <- paste("position",seq_len(19),sep = "_")
deltagu<- sapply(1:length(protospacer2), function(i){sum(rna_rna[position3[i,] %>% unname() %>% unlist()] )})

deltago<- sapply(1:length(protospacer), function(i){sum(dna_dna[position2[i,] %>% unname() %>% unlist()] )})
deltagb<- deltagh-deltago-deltagu
dat$deltagb <-deltagb

write.table(dat, file = "/proj/milovelab/mu/dukeproj/K562-supp-tables/screen_hg19_uniqe.txt", sep = "\t",quote=FALSE,
            row.names = FALSE, col.names = TRUE)

dat <- read_table(file = "/proj/milovelab/mu/dukeproj/K562-supp-tables/screen_hg19_uniqe.txt")
## read H3K27ac data and change it to GRanges
H3K27ac_cpm1kb<-read.table("/proj/milovelab/mu/dukeproj/K562-supp-tables/H3K27ac_cpm1kb.csv",header = T)
H3K27ac_cpm1kb_gr <- makeGRangesFromDataFrame(H3K27ac_cpm1kb, keep.extra.columns=TRUE)

## read DNAse data
DNAse_cpm1kb<-read.table("/proj/milovelab/mu/dukeproj/K562-supp-tables/DNase_cpm1kb.csv",header = T)
DNAse_cpm1kb_gr <- makeGRangesFromDataFrame(DNAse_cpm1kb, keep.extra.columns=TRUE)

## read ATAC-seq data
ATAC_cpm1kb<-read.table("/proj/milovelab/mu/dukeproj/K562-supp-tables/ATAC_cpm1kb.csv",header = T)
ATAC_cpm1kb_gr <- makeGRangesFromDataFrame(ATAC_cpm1kb, keep.extra.columns=TRUE)

## read H3Kme4 data
H3Kme4_cpm1kb<-read.table("/proj/milovelab/mu/dukeproj/K562-supp-tables/H3Kme4_cpm1kb.csv",header = T)
H3Kme4_cpm1kb_gr <- makeGRangesFromDataFrame(H3Kme4_cpm1kb, keep.extra.columns=TRUE)

## read TF_GATA2 data
TF_GATA2_cpm1kb<-read.table("/proj/milovelab/mu/dukeproj/K562-supp-tables/TF_GATA2_cpm1kb.csv",header = T)
TF_GATA2_cpm1kb_gr <- makeGRangesFromDataFrame(TF_GATA2_cpm1kb, keep.extra.columns=TRUE)

## read TF_TAL1 data
TF_TAL1_cpm1kb<-read.table("/proj/milovelab/mu/dukeproj/K562-supp-tables/TF_TAL1_cpm1kb.csv",header = T)
TF_TAL1_cpm1kb_gr <- makeGRangesFromDataFrame(TF_TAL1_cpm1kb, keep.extra.columns=TRUE)

## read TF_GATA2 data
TF_MYC_cpm1kb<-read.table("/proj/milovelab/mu/dukeproj/K562-supp-tables/TF_MYC_cpm1kb.csv",header = T)
TF_MYC_cpm1kb_gr <- makeGRangesFromDataFrame(TF_MYC_cpm1kb, keep.extra.columns=TRUE)

library(rtracklayer)
path = system.file(package="liftOver", "extdata", "hg38ToHg19.over.chain")
ch = import.chain(path)
H3K27ac_cpm1kb_gr = liftOver(H3K27ac_cpm1kb_gr, ch)
H3K27ac_cpm1kb_gr = unlist(H3K27ac_cpm1kb_gr)
genome(H3K27ac_cpm1kb_gr) = "hg19"

DNAse_cpm1kb_gr = liftOver(DNAse_cpm1kb_gr, ch)
DNAse_cpm1kb_gr = unlist(DNAse_cpm1kb_gr)
genome(DNAse_cpm1kb_gr) = "hg19"

ATAC_cpm1kb_gr = liftOver(ATAC_cpm1kb_gr, ch)
ATAC_cpm1kb_gr = unlist(ATAC_cpm1kb_gr)
genome(ATAC_cpm1kb_gr) = "hg19"

H3Kme4_cpm1kb_gr = liftOver(H3Kme4_cpm1kb_gr, ch)
H3Kme4_cpm1kb_gr = unlist(H3Kme4_cpm1kb_gr)
genome(H3Kme4_cpm1kb_gr) = "hg19"

TF_GATA2_cpm1kb_gr = liftOver(TF_GATA2_cpm1kb_gr, ch)
TF_GATA2_cpm1kb_gr = unlist(TF_GATA2_cpm1kb_gr)
genome(TF_GATA2_cpm1kb_gr) = "hg19"

TF_TAL1_cpm1kb_gr = liftOver(TF_TAL1_cpm1kb_gr, ch)
TF_TAL1_cpm1kb_gr = unlist(TF_TAL1_cpm1kb_gr)
genome(TF_TAL1_cpm1kb_gr) = "hg19"

TF_MYC_cpm1kb_gr = liftOver(TF_MYC_cpm1kb_gr, ch)
TF_MYC_cpm1kb_gr = unlist(TF_MYC_cpm1kb_gr)
genome(TF_MYC_cpm1kb_gr) = "hg19"

# table5 <- read_table2(file="/proj/milovelab/mu/dukeproj/K562-supp-tables/supplementary_table_5_DHS_summary_results.csv")
table5 <- read_table2("/proj/milovelab/mu/dukeproj/K562-supp-tables/supplementary_table_5_DHS_summary_results.txt.gz")
table5_gr <- makeGRangesFromDataFrame(table5, keep.extra.columns=TRUE,
                                      start.field="chromStart",end.field="chromEnd")

olap <- findOverlaps(table5_gr,H3K27ac_cpm1kb_gr) %>% as.data.frame()
olap$H3k27ac_CPM_1Kb_new <- H3K27ac_cpm1kb_gr$cpm1kb[olap$subjectHits]
H3k27ac_CPM_1Kb_new <- olap %>% group_by(queryHits) %>% summarise(H3k27ac_CPM_1Kb_new=mean(H3k27ac_CPM_1Kb_new))
table5$H3k27ac_CPM_1Kb_new <-0
table5$H3k27ac_CPM_1Kb_new[H3k27ac_CPM_1Kb_new$queryHits]<-H3k27ac_CPM_1Kb_new$H3k27ac_CPM_1Kb_new

olap <- findOverlaps(table5_gr,DNAse_cpm1kb_gr) %>% as.data.frame()
olap$DNAse_CPM_1Kb_new <- DNAse_cpm1kb_gr$cpm1kb[olap$subjectHits]
DNAse_CPM_1Kb_new <- olap %>% group_by(queryHits) %>% summarise(DNAse_CPM_1Kb_new=mean(DNAse_CPM_1Kb_new))
table5$DNAse_CPM_1Kb_new <-0
table5$DNAse_CPM_1Kb_new[DNAse_CPM_1Kb_new$queryHits]<-DNAse_CPM_1Kb_new$DNAse_CPM_1Kb_new

olap <- findOverlaps(table5_gr,ATAC_cpm1kb_gr) %>% as.data.frame()
olap$ATAC_CPM_1Kb_new <- ATAC_cpm1kb_gr$cpm1kb[olap$subjectHits]
ATAC_CPM_1Kb_new <- olap %>% group_by(queryHits) %>% summarise(ATAC_CPM_1Kb_new=mean(ATAC_CPM_1Kb_new))
table5$ATAC_CPM_1Kb_new <-0
table5$ATAC_CPM_1Kb_new[ATAC_CPM_1Kb_new$queryHits]<-ATAC_CPM_1Kb_new$ATAC_CPM_1Kb_new

olap <- findOverlaps(table5_gr,H3Kme4_cpm1kb_gr) %>% as.data.frame()
olap$H3Kme4_CPM_1Kb_new <- H3Kme4_cpm1kb_gr$cpm1kb[olap$subjectHits]
H3Kme4_CPM_1Kb_new <- olap %>% group_by(queryHits) %>% summarise(H3Kme4_CPM_1Kb_new=mean(H3Kme4_CPM_1Kb_new))
table5$H3Kme4_CPM_1Kb_new <-0
table5$H3Kme4_CPM_1Kb_new[H3Kme4_CPM_1Kb_new$queryHits]<-H3Kme4_CPM_1Kb_new$H3Kme4_CPM_1Kb_new

olap <- findOverlaps(table5_gr,TF_GATA2_cpm1kb_gr) %>% as.data.frame()
olap$TF_GATA2_CPM_1Kb_new <- TF_GATA2_cpm1kb_gr$cpm1kb[olap$subjectHits]
TF_GATA2_CPM_1Kb_new <- olap %>% group_by(queryHits) %>% summarise(TF_GATA2_CPM_1Kb_new=mean(TF_GATA2_CPM_1Kb_new))
table5$TF_GATA2_CPM_1Kb_new <-0
table5$TF_GATA2_CPM_1Kb_new[TF_GATA2_CPM_1Kb_new$queryHits]<-TF_GATA2_CPM_1Kb_new$TF_GATA2_CPM_1Kb_new

olap <- findOverlaps(table5_gr,TF_TAL1_cpm1kb_gr) %>% as.data.frame()
olap$TF_TAL1_CPM_1Kb_new <- TF_TAL1_cpm1kb_gr$cpm1kb[olap$subjectHits]
TF_TAL1_CPM_1Kb_new <- olap %>% group_by(queryHits) %>% summarise(TF_TAL1_CPM_1Kb_new=mean(TF_TAL1_CPM_1Kb_new))
table5$TF_TAL1_CPM_1Kb_new <-0
table5$TF_TAL1_CPM_1Kb_new[TF_TAL1_CPM_1Kb_new$queryHits]<-TF_TAL1_CPM_1Kb_new$TF_TAL1_CPM_1Kb_new

olap <- findOverlaps(table5_gr,TF_MYC_cpm1kb_gr) %>% as.data.frame()
olap$TF_MYC_CPM_1Kb_new <- TF_MYC_cpm1kb_gr$cpm1kb[olap$subjectHits]
TF_MYC_CPM_1Kb_new <- olap %>% group_by(queryHits) %>% summarise(TF_MYC_CPM_1Kb_new=mean(TF_MYC_CPM_1Kb_new))
table5$TF_MYC_CPM_1Kb_new <-0
table5$TF_MYC_CPM_1Kb_new[TF_MYC_CPM_1Kb_new$queryHits]<-TF_MYC_CPM_1Kb_new$TF_MYC_CPM_1Kb_new

# write.table(table5,file="/proj/milovelab/mu/dukeproj/K562-supp-tables/supplementary_table_5_DHS_summary_results.csv",row.names = F)
# table5 <- read_csv("/proj/milovelab/mu/dukeproj/K562-supp-tables/supplementary_table_5_DHS_summary_results.csv")

table1 <- read_table2("/proj/milovelab/mu/dukeproj/K562-results-julien/wgCERES-gRNAs.annotated.txt.gz")

table5_select <- table5[,c("name","annotation_wg","H3k27ac_CPM_1Kb_new", "DNAse_CPM_1Kb_new", "ATAC_CPM_1Kb_new",
                           "H3Kme4_CPM_1Kb_new","TF_GATA2_CPM_1Kb_new", "TF_TAL1_CPM_1Kb_new","TF_MYC_CPM_1Kb_new")]
# install.packages("picante")
library(picante)
#Insert the name of your dataset in the code below
cor <- cor.table(table5_select[,3:9], cor.method="pearson")

table5_select$annotation_wg_orig = table5_select$annotation_wg
table5_select$annotation_wg <- ifelse(table5_select$annotation_wg == "Promoter", "Promoter","enhancer")

table1_select <- table1[,-c(1:6,8)]

data <- dat %>% left_join(table1_select,by="gRNAid")
data <- data %>% left_join(table5_select, by=c("DHS"="name"))

# write_csv(data, "./data_Apirl_resplit/wgCERES-gRNAs-k562-discovery-screen-full.csv")
write_csv(data, "./data/dat_discovery/wgCERES-gRNAs-k562-discovery-screen-full_mean.csv")
data <- read_csv("wgCERES-gRNAs-k562-discovery-screen-full.csv")



a<-data[,c(29:30,72:78)]
a[is.na(a)] <- 0
cor <- cor.table(a, cor.method="pearson")

hist(log(data$H3k27ac_CPM_1Kb_new[which(data$annotation_wg=="Promoter")]+1),breaks=40,xlab = "log(H3K27ac_cpm_per1kb)",main = "promoter")
hist(log(data$H3k27ac_CPM_1Kb_new[which(data$annotation_wg=="enhancer")]+1),xlab = "log(H3K27ac_cpm_per1kb)",main = "enhancer")
hist(log(data$H3k27ac_CPM_1Kb_new[which(data$annotation_wg_orig=="Intron")]+1),xlab = "log(H3K27ac_cpm_per1kb)",main = "enhancer")

hist(log(data$DNase_CPM_per_1kbp[which(data$annotation_wg=="Promoter")]+1),breaks=40,xlab = "log(DNAse_cpm_per1kb)",main = "promoter")
hist(log(data$DNase_CPM_per_1kbp[which(data$annotation_wg=="enhancer")]+1),breaks=40,xlab = "log(DNAse_cpm_per1kb)",main = "enhancer")

hist(log(data$ATAC_CPM_1Kb_new[which(data$annotation_wg=="Promoter")]+1),breaks=40,xlab = "log(ATAC_cpm_per1kb)",main = "promoter")
hist(log(data$ATAC_CPM_1Kb_new[which(data$annotation_wg=="enhancer")]+1),breaks=40,xlab = "log(ATAC_cpm_per1kb)",main = "enhancer")

hist(log(data$H3Kme4_CPM_1Kb_new[which(data$annotation_wg=="Promoter")]+1),breaks=40,xlab = "log(H3Kme4_cpm_per1kb)",main = "promoter")
hist(log(data$H3Kme4_CPM_1Kb_new[which(data$annotation_wg=="enhancer")]+1),breaks=40,xlab = "log(H3Kme4_cpm_per1kb)",main = "enhancer")

hist(log(data$TF_GATA2_CPM_1Kb_new[which(data$annotation_wg=="Promoter")]+1),breaks=40,xlab = "log(TF_GATA2_cpm_per1kb)",main = "promoter")
hist(log(data$TF_GATA2_CPM_1Kb_new[which(data$annotation_wg=="enhancer")]+1),breaks=40,xlab = "log(TF_GATA2_cpm_per1kb)",main = "enhancer")

hist(log(data$TF_TAL1_CPM_1Kb_new[which(data$annotation_wg=="Promoter")]+1),breaks=40,xlab = "log(TF_TAL1_cpm_per1kb)",main = "promoter")
hist(log(data$TF_TAL1_CPM_1Kb_new[which(data$annotation_wg=="enhancer")]+1),breaks=40,xlab = "log(TF_TAL1_cpm_per1kb)",main = "enhancer")

hist(log(data$TF_MYC_CPM_1Kb_new[which(data$annotation_wg=="Promoter")]+1),breaks=40,xlab = "log(TF_MYC_cpm_per1kb)",main = "promoter")
hist(log(data$TF_MYC_CPM_1Kb_new[which(data$annotation_wg=="enhancer")]+1),breaks=40,xlab = "log(TF_MYC_cpm_per1kb)",main = "enhancer")

plot(x=log(data$H3k27ac_CPM_1Kb_new[which(data$annotation_wg=="Promoter")]+1),
     y=data$log2FoldChange[which(data$annotation_wg=="Promoter")],
     xlab = "log(H3K27ac_cpm_per1kb)",main = "promoter",ylab="logFC")
plot(x=log(data$H3k27ac_CPM_1Kb_new[which(data$annotation_wg=="enhancer")]+1),
     y=data$log2FoldChange[which(data$annotation_wg=="enhancer")],
     xlab = "log(H3K27ac_cpm_per1kb)",main = "enhancer",ylab="logFC")

plot(x=log(data$DNAse_CPM_1Kb_new[which(data$annotation_wg=="Promoter")]+1),
     y=data$log2FoldChange[which(data$annotation_wg=="Promoter")],
     xlab = "log(DNAse_cpm_per1kb)",main = "promoter",ylab="logFC")
plot(x=log(data$DNAse_CPM_1Kb_new[which(data$annotation_wg=="enhancer")]+1),
     y=data$log2FoldChange[which(data$annotation_wg=="enhancer")],
     xlab = "log(DNAse_cpm_per1kb)",main = "enhancer",ylab="logFC")

plot(x=log(data$ATAC_CPM_1Kb_new[which(data$annotation_wg=="Promoter")]+1),
     y=data$log2FoldChange[which(data$annotation_wg=="Promoter")],
     xlab = "log(ATAC_cpm_per1kb)",main = "promoter",ylab="logFC")
plot(x=log(data$ATAC_CPM_1Kb_new[which(data$annotation_wg=="enhancer")]+1),
     y=data$log2FoldChange[which(data$annotation_wg=="enhancer")],
     xlab = "log(ATAC_cpm_per1kb)",main = "enhancer",ylab="logFC")

plot(x=log(data$TF_GATA2_CPM_1Kb_new[which(data$annotation_wg=="Promoter")]+1),
     y=data$log2FoldChange[which(data$annotation_wg=="Promoter")],
     xlab = "log(TF_GATA2_cpm_per1kb)",main = "promoter",ylab="logFC")
plot(x=log(data$TF_TAL1_CPM_1Kb_new[which(data$annotation_wg=="enhancer")]+1),
     y=data$log2FoldChange[which(data$annotation_wg=="enhancer")],
     xlab = "log(TF_TAL1_cpm_per1kb)",main = "enhancer",ylab="logFC")

plot(x=log(data$TF_MYC_CPM_1Kb_new[which(data$annotation_wg=="Promoter")]+1),
     y=data$log2FoldChange[which(data$annotation_wg=="Promoter")],
     xlab = "log(TF_GATA2_cpm_per1kb)",main = "promoter",ylab="logFC")
plot(x=log(data$TF_MYC_CPM_1Kb_new[which(data$annotation_wg=="enhancer")]+1),
     y=data$log2FoldChange[which(data$annotation_wg=="enhancer")],
     xlab = "log(TF_TAL1_cpm_per1kb)",main = "enhancer",ylab="logFC")
# protospacer <- dat[,"protospacer"]
# a <- unique(sapply(protospacer, nchar))
# position <- matrix(NA,ncol=a,nrow=nrow(dat)) %>% as.data.frame
# for(i in 1:a){
#   position[,i] <- as.factor(substr(protospacer$protospacer,i,i))
# }
# colnames(position) <- paste("position",seq_len(20),sep = "_")
# position2 <- matrix(NA,ncol=19,nrow=nrow(dat)) %>% as.data.frame
# for(i in 1:19){
#   position2[,i] <- as.factor(substr(protospacer$protospacer,i,i+1))
# }
# colnames(position2) <- paste("position",seq_len(19),sep = "_")
# 
# library(caret)
# f <- paste("~", paste(colnames(position), collapse="+"))
# dummy <- dummyVars(f, data=position,sep = "_")
# newdata <- data.frame(predict(dummy, newdata = position))
# 
# f <- paste("~", paste(colnames(position2), collapse="+"))
# dummy <- dummyVars(f, data=position2,sep = "_")
# newdata2 <- data.frame(predict(dummy, newdata = position2))
# 
# dinucleotide<- levels(position2$position_1)
# count <- matrix(NA,ncol=16,nrow=nrow(dat)) %>% as.data.frame
# for(i in 1:16){
#   count[,i] <- str_count(protospacer$protospacer,dinucleotide[i])
# }
# colnames(count) <- paste0(dinucleotide,"count")

## add converse 
data.bed <- data[,c("chrom","chromStart","chromEnd","gRNAid")]
write.table(data.bed,file = "wgCERES-gRNAs-k562-discovery-screen-full.bed",row.names = FALSE,col.names = FALSE, sep="\t", quote=FALSE)
converse <- rtracklayer::import("./data/dat_April_resplit/overlap.bed", format="bed")
converse2 <- converse %>% group_by(name) %>% reduce_ranges()
write_bed(converse2,file = "overlap_reduce.bed")
converse2 <- import("./data/dat_April_resplit/overlap_reduce.bed", format="bed")

gr <- makeGRangesFromDataFrame(data.bed, keep.extra.columns=TRUE,strand.field = "strand.y",
                               start.field="chromStart",end.field="chromEnd")
# converse2$length <- width(converse2)
converse2$name <- as.numeric(converse2$name)
left_rng <- join_overlap_left(gr, converse2) %>% 
  # mutate(score.ind = name*length) %>% 
  group_by(gRNAid) %>% 
  summarise(length.total=n(), 
            score.total = sum(name),
            score.max = max(name),
            score.median = median(name)) %>% 
  as.data.frame() %>% 
  mutate(score.mean = score.total/length.total)
left_rng[is.na(left_rng)]<-0
data2 <- data %>% 
  left_join(left_rng %>% select(gRNAid,score.mean,score.max,score.median)) 
# overlap.index<- findOverlaps(gr,converse) %>% as.data.frame() %>% 
#   modify_if(is.numeric,as.character)
# converse.df <- converse %>% as.data.frame() %>% 
#   rownames_to_column(var = "subjectHits") %>% 
#   right_join(overlap.index)
# converse.df$name <- as.numeric(converse.df$name)
# data.score <- converse.df %>% 
#   mutate(score.ind= width*name) %>% 
#   group_by(queryHits) %>% 
#   summarise(sumwidth = sum(width),score = sum(score.ind)/sumwidth)

# data2 <- data.bed %>% 
#   rownames_to_column(var = "queryHits") %>% 
#   right_join(data.score) %>% 
#   select(gRNAid,score) %>% 
#   right_join(data)

data2 = data2 %>%
  mutate(strand.y = ifelse(strand.y == "+", 1, 
                         ifelse(strand.y == "-", -1, NA)),
         ploidyZhou = as.numeric(str_sub(ploidyZhou, 3, 3)),  ## imputation with mode
         LossHetZhou = as.numeric(LossHetZhou),
         SV_Zhou = as.numeric(SV_Zhou),
         SNV_Zhou = as.numeric(SNV_Zhou),
         Conserved = as.numeric(Conserved),
         Gquad_same_strand = as.numeric(Gquad_same_strand),
         Gquad_other_strand = as.numeric(Gquad_other_strand),
         HartEssential = as.numeric(HartEssential),
         vc_sqrt_sum = as.numeric(as.character(data$vc_sqrt_sum)))
names(data2)[names(data2) == 'H3Kme4_CPM_1Kb_new'] <- "H3K4me3_CPM_1Kb_new"
names(data2)[names(data2) == 'strand.y'] <- "strand"
# write_csv(data2, "wgCERES-gRNAs-k562-discovery-screen-full-add-converse.csv")
write_csv(data2, "./data/dat_discovery/wgCERES-gRNAs-k562-discovery-screen-full-add-converse.csv")
data2 <- read_csv("./data/dat_discovery/wgCERES-gRNAs-k562-discovery-screen-full-add-converse.csv")

## select promoter
dat_pro = data2 %>% filter(annotation_wg=="Promoter")

ggplot(dat_pro %>% filter(baseMean<1000),aes(x= baseMean, y=log10(pvalue)))+
  geom_point()+
  geom_hline(yintercept=log10(0.05), linetype="dashed", color = "red")

## logFC prediction
# library(caret)
# set.seed(2021)
dat_padj_sel = dat_pro %>%
  # filter(padj <= 0.2)
  filter(padj <= 0.1)
#### transform extreme value
data_new = dat_padj_sel %>%
  # rename(strand = `strand.y`) %>%
  select(protospacer:deltagb, olapnumber:padj, deltagh:Gcount, 
         ploidyZhou:Gquad_n_overlap_other_strand,
         distance:vc_sqrt_sum, 
         annotation_wg:TF_MYC_CPM_1Kb_new,score.mean:score.median)
library(DescTools)
for (i in c(34,46,48:54,9:10)){
  print(colnames(data_new)[i])
  data_new[,i] = Winsorize(data_new[,i], probs = c(0,0.95), na.rm = T)
}


neg_thresh <- quantile(data_new$log2FoldChange[data_new$log2FoldChange<0],0.15)
data_new$class <- ifelse(data_new$log2FoldChange<=neg_thresh,0,ifelse(data_new$log2FoldChange>0,2,1))

## choose chr21, chrX, chr10, chr14, chr18 chr21 as test data approximate equal to 20.3%
# testid <- data_new$chrom %in% c("chr2","chr6","chr10","chr14","chr18","chr21")
# train = data_new[!testid,]
# test = data_new[testid,]
# write_csv(train, "wgCERES-gRNAs-k562-discovery-screen-pro_logFC0.2-train.csv")
# write_csv(test, "wgCERES-gRNAs-k562-discovery-screen-pro_logFC0.2-test.csv")

# write_csv(train, "wgCERES-gRNAs-k562-discovery-screen-pro_logFC0.1-train.csv")
# write_csv(test, "wgCERES-gRNAs-k562-discovery-screen-pro_logFC0.1-test.csv")

##### binary outcome prediction: filter out padj between 0.05 and 0.2
dat_sel = dat_pro %>%
  filter(baseMean >125)

# ggplot(dat_sel,aes(x=log2FoldChange))+
#   geom_density()+
#   facet_wrap(~significant)

write_csv(dat_sel, "wgCERES-gRNAs-k562-discovery-screen-pro_baseMean125.csv")

dat_sel <- read_csv("./data/dat_April_resplit/wgCERES-gRNAs-k562-discovery-screen-pro_baseMean125.csv")
dat_sel_binary = dat_sel %>%
  filter(padj > 0.2 | padj <= 0.05) %>% 
  mutate(significance = ifelse(padj<=0.05,1,0))
data_new = dat_sel_binary %>%
  # rename(strand = strand.y) %>%
  select(protospacer:deltagb, olapnumber:padj, deltagh:Gcount, 
         ploidyZhou:Gquad_n_overlap_other_strand,
         distance:vc_sqrt_sum, 
         annotation_wg:TF_MYC_CPM_1Kb_new,significance,score.mean:score.median)

# library(TxDb.Hsapiens.UCSC.hg19.knownGene)
# # Load the GenomicRanges package for working with genomic intervals
# library(GenomicRanges)

# # Extract protein-coding exons from the TxDb object
# txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
# columns(txdb)
# keytypes(txdb)
# gene_ids <- keys(txdb)

# # Ensure that gene_ids is a character vector
# gene_ids <- as.character(gene_ids)

# # Use the select function to extract the protein coding gene regions
# gene_regions <- select(txdb, keys = gene_ids, keytype = "GENEID", 
#                         columns = c("EXONID","TXSTART","TXEND","TXSTRAND"), overlap = "any")
# makeGRangesFromDataFrame(TF_GATA2_cpm1kb, keep.extra.columns=TRUE)
# # Collapse exons into protein-coding gene regions
# gr_protein_coding <- reduce(gene_regions)                 


for (i in c(34,46,48:54,9:10)){
  print(colnames(data_new)[i])
  data_new[,i] = Winsorize(data_new[,i], probs = c(0,0.95), na.rm = T)
}


prop.table(table(data_new$chrom,data_new$significance),1)
# set.seed(2021)
# trainid <- createDataPartition(dat_sel_binary$significant, p=0.8, list=FALSE)
# trainid = createDataPartition(dat_sel_binary$direction_sig, p=0.8, list=FALSE)
freq <- table(data_new$chrom)/nrow(data_new)*100
freq2 <- sample(freq)
cumsum(freq2) # split the fold based on cumsum
split <- list("a" = c("chr21","chrX","chr10","chr4","chr16","chr17"), "b" = c("chr22","chr14","chr1"),
        "c" =c("chr7","chr2","chr11","chr13"),"d"= c("chr19","chr18", "chr5", "chr6", "chr3"),
        "e" = c("chr15", "chr8", "chr20", "chr9", "chr12"))

for(i in 1:5){
      testid <- data_new$chrom %in% split[[i]]
      train = data_new[!testid,]
      test = data_new[testid,]
      write_csv(train, paste0("./data/dat_discovery/wgCERES-gRNAs-k562-discovery-screen-pro_baseMean125-binary-", i, "-train2.csv"))
      write_csv(test, paste0("./data/dat_discovery/wgCERES-gRNAs-k562-discovery-screen-pro_baseMean125-binary-", i, "-test2.csv"))
      # write_csv(train, paste0("../dat_fivefold/wgCERES-gRNAs-k562-discovery-screen-pro_baseMean125-binary-", i, "-train.csv"))
      # write_csv(test, paste0("../dat_fivefold/wgCERES-gRNAs-k562-discovery-screen-pro_baseMean125-binary-", i, "-test.csv"))
}

for(i in 1:5){
      test = read_csv(paste0("../dat_discovery/wgCERES-gRNAs-k562-discovery-screen-pro_baseMean125-binary-", i, "-test.csv"))
      # create covariate for selecting DHS that has significant results
      test_rank <- test |> group_by(DHS)  |> summarize(nsig = sum(significance)) |> 
              mutate(DHS_select = ifelse(nsig>=3,1,0)) |> 
              ungroup()
      test2 <- test  |> left_join(test_rank)
      write_csv(test2, paste0("../dat_discovery/wgCERES-gRNAs-k562-discovery-screen-pro_baseMean125-binary-", i, "-test-forrank.csv"))
}
## select enhancer
dat_enh = data2 %>% filter(annotation_wg=="enhancer")
# write_csv(dat_nonsig, "./data_new/wgCERES-gRNAs-k562-discovery-screen-enh_rawp-nonsignificant.csv")


dat_padj_sel = dat_enh %>%
  filter(padj <= 0.2)
  # filter(padj <= 0.1)
#### transform extreme value
data_new = dat_padj_sel %>%
  # rename(strand = strand.y) %>%
  select(protospacer:deltagb, promnumber:padj, deltagh:Gcount, 
         ploidyZhou:Gquad_n_overlap_other_strand,
         distance:vc_sqrt_sum, 
         annotation_wg:TF_MYC_CPM_1Kb_new,score.mean:score.median)

for (i in c(36,48,50:56,9:12)){
  print(colnames(data_new)[i])
  data_new[,i] = Winsorize(data_new[,i], probs = c(0,0.95), na.rm = T)
}

table(data_new$chrom)/nrow(data_new)*100
neg_thresh <- quantile(data_new$log2FoldChange[data_new$log2FoldChange<0],0.15)
data_new$class <- ifelse(data_new$log2FoldChange<=neg_thresh,0,ifelse(data_new$log2FoldChange>0,2,1))
## choose chr2, chr6, chr10, chr14, chr18, chr22 as test data
testid <- data_new$chrom %in% c("chr2","chr6","chr10","chr14","chr18","chr21")
train = data_new[!testid,]
test = data_new[testid,]

write_csv(train, "wgCERES-gRNAs-k562-discovery-screen-enh_logFC0.2-train.csv")
write_csv(test, "wgCERES-gRNAs-k562-discovery-screen-enh_logFC0.2-test.csv")

write_csv(train, "wgCERES-gRNAs-k562-discovery-screen-enh_logFC0.1-train.csv")
write_csv(test, "wgCERES-gRNAs-k562-discovery-screen-enh_logFC0.1-test.csv")

##### binary outcome prediction: filter out padj between 0.05 and 0.2
dat_sel = dat_enh %>%
  filter(baseMean>125)

write_csv(dat_sel, "wgCERES-gRNAs-k562-discovery-screen-enh_baseMean125.csv")
dat_sel <- read_csv("wgCERES-gRNAs-k562-discovery-screen-enh_baseMean125.csv")

dat_sel_binary = dat_sel %>%
  filter(padj > 0.2 | padj <= 0.05) %>% 
  mutate(significance = ifelse(padj<=0.05,1,0))
data_new = dat_sel_binary %>%
  # rename(strand = strand.y) %>%
  select(protospacer:deltagb, promnumber:padj, deltagh:Gcount, 
         ploidyZhou:Gquad_n_overlap_other_strand,
         distance:vc_sqrt_sum, 
         annotation_wg:TF_MYC_CPM_1Kb_new,significance,score.mean:score.median)

for (i in c(36,48,50:56,9:12)){
  print(colnames(data_new)[i])
  data_new[,i] = Winsorize(data_new[,i], probs = c(0,0.95), na.rm = T)
}
prop.table(table(data_new$chrom,data_new$significance),1)

for(i in 1:5){
      testid <- data_new$chrom %in% split[[i]]
      train = data_new[!testid,]
      test = data_new[testid,]

      write_csv(train, paste0("../dat_fivefold/wgCERES-gRNAs-k562-discovery-screen-enh_baseMean125-binary-", i, "-train.csv"))
      write_csv(test, paste0("../dat_fivefold/wgCERES-gRNAs-k562-discovery-screen-enh_baseMean125-binary-", i, "-test.csv"))
}

for(i in 1:5){
      test = read_csv(paste0("../dat_discovery/wgCERES-gRNAs-k562-discovery-screen-enh_baseMean125-binary-", i, "-test.csv"))
      # create covariate for selecting DHS that has significant results
      test_rank <- test |> group_by(DHS)  |> summarize(nsig = sum(significance)) |> 
              mutate(DHS_select = ifelse(nsig>=3,1,0)) |> 
              ungroup()
      test2 <- test  |> left_join(test_rank)
      write_csv(test2, paste0("../dat_discovery/wgCERES-gRNAs-k562-discovery-screen-enh_baseMean125-binary-", i, "-test-forrank.csv"))
}
