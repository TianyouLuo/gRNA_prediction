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
  mcols(gi) <- df[,7] # Interaction frequencies
  return(gi)
}

hic.gi <- convertToGI(hic)
library(vroom)
library(plyranges)
library(data.table)
dat <- vroom("/proj/milovelab/mu/dukeproj/data/mhc/k562/k562.latentvar.MAST_expect_cells.grna_de.markers.aggr.txt.gz")
# dat2 <- vroom("/proj/milovelab/mu/dukeproj/data/mhc/k562.expect_cells.grna.de.markers.MAST.annotatedfull.final.update20220117.LRB.tsv")
dat2 <- vroom("/proj/milovelab/mu/dukeproj/data/mhc/k562/k562.expect_cells.grna.de.markers.MAST.annotatedfull.final.update20220117.LRB.tsv")

# dat3 <- dat %>% left_join(dat2[,7:47], by = c("grna" = "protospacer"))
colnames(dat)
unique(dat$grna) %>% length()
dat3 <- dat %>% group_by(grna) %>% slice_min(p_val,n=1,with_ties = FALSE) %>% 
  left_join(dat2[,c(7:20,40:47)], by = c("grna" = "protospacer", "gene_symbol"))
dat3 <- dat3[!duplicated(dat3), ] %>% filter(type == "targeting")
# dat2 <- setDT(dat)[ , .SD[which.min(p_val)], by = protospacer]   # Min of groups
gr <- makeGRangesFromDataFrame(dat3, keep.extra.columns=TRUE,
                               seqnames.field="grna.chr",
                               strand.field = "grna.strand",
                               start.field="grna.start",end.field="grna.end")

olap <- findOverlaps(gr,hic.gi) %>% as.data.frame()
# olap$log10fdr <- hic.gi$log10fdr[olap$subjectHits]
# sumlogfdr <- olap %>% group_by(queryHits) %>% summarise(sumlog10fdr=sum(log10fdr))

countolap <- countOverlaps(gr,hic.gi)
ref = read_tsv("/proj/milovelab/mu/dukeproj/K562-results-julien/gRNA_prediction/gene_info.gencode.v28lift37.txt")
refgr <- makeGRangesFromDataFrame(ref, keep.extra.columns=TRUE,
                                  start.field="Promoter_Start",end.field="Promoter_End")

out <- linkOverlaps(hic.gi, gr, refgr) %>% as.data.frame()
out <- unique(out[,c(1,2)])
# out$log10fdr <- hic.gi$log10fdr[out$query]
prom <- out %>% dplyr::count(subject1)
# promfdr <- out %>% group_by(subject1) %>% summarise(sumlog10fdr=sum(log10fdr))

dat<-dat3
dat$promnumber <- rep(0,nrow(dat))
dat$promnumber[prom$subject1]<-prom$n
# dat$promlog10fdr <- rep(0,nrow(dat))
# dat$promlog10fdr[promfdr$subject1]<-promfdr$sumlog10fdr
dat$olapnumber <- countolap
# dat$sumlog10fdr <- rep(0,nrow(dat))
# dat$sumlog10fdr[sumlogfdr$queryHits]<-sumlogfdr$sumlog10fdr

## 581 dhs
unique(dat$dhs) %>% length()
dhs_gr <- makeGRangesFromDataFrame(dat, keep.extra.columns=TRUE,
                                   seqnames.field="dhs.chr",
                                   strand.field = "grna.strand",
                                   start.field="dhs.start",end.field="dhs.end")
## read H3K27ac data and change it to GRanges

H3K27ac_cpm1kb<-read.table("/proj/milovelab/mu/dukeproj/K562-supp-tables/H3K27ac_cpm1kb.csv",header = T)
H3K27ac_cpm1kb_gr <- makeGRangesFromDataFrame(H3K27ac_cpm1kb, keep.extra.columns=TRUE)

## read ATAC-seq data
ATAC_cpm1kb<-read.table("/proj/milovelab/mu/dukeproj/K562-supp-tables/ATAC_cpm1kb.csv",header = T)
ATAC_cpm1kb_gr <- makeGRangesFromDataFrame(ATAC_cpm1kb, keep.extra.columns=TRUE)

## read H3K4me3 data
H3K4me3_cpm1kb<-read.table("/proj/milovelab/mu/dukeproj/K562-supp-tables/H3Kme4_cpm1kb.csv",header = T)
H3K4me3_cpm1kb_gr <- makeGRangesFromDataFrame(H3K4me3_cpm1kb, keep.extra.columns=TRUE)

## overlap with dhs level
olap <- findOverlaps(dhs_gr,H3K27ac_cpm1kb_gr) %>% as.data.frame()
olap$H3k27ac_CPM_1Kb_new <- H3K27ac_cpm1kb_gr$cpm1kb[olap$subjectHits]
H3k27ac_CPM_1Kb_new <- olap %>% group_by(queryHits) %>% summarise(H3k27ac_CPM_1Kb_new=mean(H3k27ac_CPM_1Kb_new))
dat$H3k27ac_CPM_1Kb_new <-0
dat$H3k27ac_CPM_1Kb_new[H3k27ac_CPM_1Kb_new$queryHits]<-H3k27ac_CPM_1Kb_new$H3k27ac_CPM_1Kb_new

olap <- findOverlaps(dhs_gr,ATAC_cpm1kb_gr) %>% as.data.frame()
olap$ATAC_CPM_1Kb_new <- ATAC_cpm1kb_gr$cpm1kb[olap$subjectHits]
ATAC_CPM_1Kb_new <- olap %>% group_by(queryHits) %>% summarise(ATAC_CPM_1Kb_new=mean(ATAC_CPM_1Kb_new))
dat$ATAC_CPM_1Kb_new <-0
dat$ATAC_CPM_1Kb_new[ATAC_CPM_1Kb_new$queryHits]<-ATAC_CPM_1Kb_new$ATAC_CPM_1Kb_new

olap <- findOverlaps(dhs_gr,H3K4me3_cpm1kb_gr) %>% as.data.frame()
olap$H3K4me3_CPM_1Kb_new <- H3K4me3_cpm1kb_gr$cpm1kb[olap$subjectHits]
H3K4me3_CPM_1Kb_new <- olap %>% group_by(queryHits) %>% summarise(H3K4me3_CPM_1Kb_new=mean(H3K4me3_CPM_1Kb_new))
dat$H3K4me3_CPM_1Kb_new <-0
dat$H3K4me3_CPM_1Kb_new[H3K4me3_CPM_1Kb_new$queryHits]<-H3K4me3_CPM_1Kb_new$H3K4me3_CPM_1Kb_new

## deltagb and deltagh
protospacer <- dat$grna
strand = dat$grna.strand
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
write.table(dat,file="/proj/milovelab/mu/dukeproj/data/mhc/k562/k562_onegene_full.csv",row.names = F)

dat <- read.table(file="/proj/milovelab/mu/dukeproj/data/mhc/k562/k562_onegene_full.csv", header = T)
dat %>%
  select(grna.chr, grna.start, grna.end, grna, dhs, grna.strand) %>%
  write_tsv("/proj/milovelab/mu/dukeproj/data/mhc/k562/k562_coordinates.tsv", col_names = F)


library(data.table)
library(tidyverse)

dat <- read.table(file="/proj/milovelab/mu/dukeproj/data/mhc/k562/k562_onegene_full.csv", header = T)
dir2 = "/proj/milovelab/mu/dukeproj/data/mhc/k562"
OGEE = fread("/proj/milovelab/mu/dukeproj/data/mhc/k562/k562_OGEE.txt")
OGEE_downstream = fread("/proj/milovelab/mu/dukeproj/data/mhc/k562/k562_OGEE_downstream.txt")
dat = dat %>% 
  left_join(OGEE %>% select(grna, distance, OGEE_prop_Essential), by = "grna") %>% 
  left_join(OGEE_downstream %>% select(grna, distance, OGEE_prop_Essential) %>% 
              dplyr::rename(distance_downstream = distance, OGEE_prop_Essential_downstream = OGEE_prop_Essential), by = "grna")
write_tsv(dat, file.path(dir2, "k562_full.tsv"))

dat_sel = dat %>%
  filter(pval_fdr_corrected <= 0.05 | pval_fdr_corrected >= 0.2) %>%
  mutate(significant = ifelse(pval_fdr_corrected <= 0.05, 1, 0))
table(dat_sel$significant)
write_csv(dat_sel, file.path(dir2, "k562-rawp0.05.csv"))

## randomly sample dhs to split train and test data
## Binary classfication
set.seed(2050)
dhs <-dat_sel$dhs %>% unique() %>% sort()
sample = sample(1:length(dhs), length(dhs))
for (i in 1:5){
  s1 = ceiling(length(dhs)/5) * (i-1) + 1
  s2 = min(ceiling(length(dhs)/5) * i, length(dhs))
  dhs_fold = dhs[sample[s1:s2]]
  test = dat_sel %>% filter(dhs %in% dhs_fold)
  train = dat_sel %>% filter(!(dhs %in% dhs_fold))
  print(paste0("Fold ", i, " # of gRNA: ", dim(test)[1]))
  write_csv(train, file.path(dir2, paste0("k562-pfdr0.05-pfdr0.2-binary-fold",i,"-train.csv")))
  write_csv(test, file.path(dir2, paste0("k562-pfdr0.05-pfdr0.2-binary-fold",i,"-test.csv")))
}


