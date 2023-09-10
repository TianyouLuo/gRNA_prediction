library(vroom)
library(plyranges)
library(data.table)
# have gRNA location information
# dat_old <- vroom("/proj/milovelab/mu/dukeproj/data/scOct4th/grna.de.markers.nFeature_RNA.MAST.txt")
dat <- vroom("/proj/milovelab/mu/dukeproj/data/scOct4th/grna.de.markers.nFeature_RNA.with_grna_coords.txt.gz")
dat2 <- vroom("/pine/scr/t/i/tianyou/Patrick/single-cell/single-cell-data.csv") # for DHS location information
colnames(dat)
unique(dat$grna) %>% length()  
library(tidyverse, help, pos = 2, lib.loc = NULL)
# dat3 <- dat %>% group_by(grna) %>% slice_min(p_val,n=1,with_ties = FALSE) %>% 
#   left_join(dat2[,c(7:18)], by = c("grna", "gene_symbol"))
dat2 <- dat2[,c(7,14:17)] %>% as.data.table() %>% unique()
dat3 <- dat %>% group_by(grna) %>% slice_min(p_val,n=1,with_ties = FALSE) %>% inner_join(dat2, by = c("grna"))
dat3 <- dat3  |> rename(grna.chrom = chrom, grna.start = start, grna.end = end)
dhs_gr <- makeGRangesFromDataFrame(dat3, keep.extra.columns=TRUE, seqnames.field = "dhs_chrom",
                                      start.field="dhs_start",end.field="dhs_end")

## read H3K27ac data and change it to GRanges
H3K27ac_cpm1kb<-read.table("/proj/milovelab/mu/dukeproj/K562-supp-tables/H3K27ac_cpm1kb.csv",header = T)
H3K27ac_cpm1kb_gr <- makeGRangesFromDataFrame(H3K27ac_cpm1kb, keep.extra.columns=TRUE)

## read DNAse data
DNAse_cpm1kb<-read.table("/proj/milovelab/mu/dukeproj/K562-supp-tables/DNase_cpm1kb.csv",header = T)
DNAse_cpm1kb_gr <- makeGRangesFromDataFrame(DNAse_cpm1kb, keep.extra.columns=TRUE)

## read ATAC-seq data
ATAC_cpm1kb<-read.table("/proj/milovelab/mu/dukeproj/K562-supp-tables/ATAC_cpm1kb.csv",header = T)
ATAC_cpm1kb_gr <- makeGRangesFromDataFrame(ATAC_cpm1kb, keep.extra.columns=TRUE)

## read H3K4me3 data
H3K4me3_cpm1kb<-read.table("/proj/milovelab/mu/dukeproj/K562-supp-tables/H3Kme4_cpm1kb.csv",header = T)
H3K4me3_cpm1kb_gr <- makeGRangesFromDataFrame(H3K4me3_cpm1kb, keep.extra.columns=TRUE)

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

H3K4me3_cpm1kb_gr = liftOver(H3K4me3_cpm1kb_gr, ch)
H3K4me3_cpm1kb_gr = unlist(H3K4me3_cpm1kb_gr)
genome(H3K4me3_cpm1kb_gr) = "hg19"

TF_GATA2_cpm1kb_gr = liftOver(TF_GATA2_cpm1kb_gr, ch)
TF_GATA2_cpm1kb_gr = unlist(TF_GATA2_cpm1kb_gr)
genome(TF_GATA2_cpm1kb_gr) = "hg19"

TF_TAL1_cpm1kb_gr = liftOver(TF_TAL1_cpm1kb_gr, ch)
TF_TAL1_cpm1kb_gr = unlist(TF_TAL1_cpm1kb_gr)
genome(TF_TAL1_cpm1kb_gr) = "hg19"

TF_MYC_cpm1kb_gr = liftOver(TF_MYC_cpm1kb_gr, ch)
TF_MYC_cpm1kb_gr = unlist(TF_MYC_cpm1kb_gr)
genome(TF_MYC_cpm1kb_gr) = "hg19"

olap <- findOverlaps(dhs_gr,H3K27ac_cpm1kb_gr) %>% as.data.frame()
olap$H3k27ac_CPM_1Kb_new <- H3K27ac_cpm1kb_gr$cpm1kb[olap$subjectHits]
H3k27ac_CPM_1Kb_new <- olap %>% group_by(queryHits) %>% summarise(H3k27ac_CPM_1Kb_new=mean(H3k27ac_CPM_1Kb_new))
dat3$H3k27ac_CPM_1Kb_new <-0
dat3$H3k27ac_CPM_1Kb_new[H3k27ac_CPM_1Kb_new$queryHits]<-H3k27ac_CPM_1Kb_new$H3k27ac_CPM_1Kb_new

olap <- findOverlaps(dhs_gr,DNAse_cpm1kb_gr) %>% as.data.frame()
olap$DNAse_CPM_1Kb_new <- DNAse_cpm1kb_gr$cpm1kb[olap$subjectHits]
DNAse_CPM_1Kb_new <- olap %>% group_by(queryHits) %>% summarise(DNAse_CPM_1Kb_new=mean(DNAse_CPM_1Kb_new))
dat3$DNAse_CPM_1Kb_new <-0
dat3$DNAse_CPM_1Kb_new[DNAse_CPM_1Kb_new$queryHits]<-DNAse_CPM_1Kb_new$DNAse_CPM_1Kb_new

olap <- findOverlaps(dhs_gr,ATAC_cpm1kb_gr) %>% as.data.frame()
olap$ATAC_CPM_1Kb_new <- ATAC_cpm1kb_gr$cpm1kb[olap$subjectHits]
ATAC_CPM_1Kb_new <- olap %>% group_by(queryHits) %>% summarise(ATAC_CPM_1Kb_new=mean(ATAC_CPM_1Kb_new))
dat3$ATAC_CPM_1Kb_new <-0
dat3$ATAC_CPM_1Kb_new[ATAC_CPM_1Kb_new$queryHits]<-ATAC_CPM_1Kb_new$ATAC_CPM_1Kb_new

olap <- findOverlaps(dhs_gr,H3K4me3_cpm1kb_gr) %>% as.data.frame()
olap$H3K4me3_CPM_1Kb_new <- H3K4me3_cpm1kb_gr$cpm1kb[olap$subjectHits]
H3K4me3_CPM_1Kb_new <- olap %>% group_by(queryHits) %>% summarise(H3K4me3_CPM_1Kb_new=mean(H3K4me3_CPM_1Kb_new))
dat3$H3K4me3_CPM_1Kb_new <-0
dat3$H3K4me3_CPM_1Kb_new[H3K4me3_CPM_1Kb_new$queryHits]<-H3K4me3_CPM_1Kb_new$H3K4me3_CPM_1Kb_new

olap <- findOverlaps(dhs_gr,TF_GATA2_cpm1kb_gr) %>% as.data.frame()
olap$TF_GATA2_CPM_1Kb_new <- TF_GATA2_cpm1kb_gr$cpm1kb[olap$subjectHits]
TF_GATA2_CPM_1Kb_new <- olap %>% group_by(queryHits) %>% summarise(TF_GATA2_CPM_1Kb_new=mean(TF_GATA2_CPM_1Kb_new))
dat3$TF_GATA2_CPM_1Kb_new <-0
dat3$TF_GATA2_CPM_1Kb_new[TF_GATA2_CPM_1Kb_new$queryHits]<-TF_GATA2_CPM_1Kb_new$TF_GATA2_CPM_1Kb_new

olap <- findOverlaps(dhs_gr,TF_TAL1_cpm1kb_gr) %>% as.data.frame()
olap$TF_TAL1_CPM_1Kb_new <- TF_TAL1_cpm1kb_gr$cpm1kb[olap$subjectHits]
TF_TAL1_CPM_1Kb_new <- olap %>% group_by(queryHits) %>% summarise(TF_TAL1_CPM_1Kb_new=mean(TF_TAL1_CPM_1Kb_new))
dat3$TF_TAL1_CPM_1Kb_new <-0
dat3$TF_TAL1_CPM_1Kb_new[TF_TAL1_CPM_1Kb_new$queryHits]<-TF_TAL1_CPM_1Kb_new$TF_TAL1_CPM_1Kb_new

olap <- findOverlaps(dhs_gr,TF_MYC_cpm1kb_gr) %>% as.data.frame()
olap$TF_MYC_CPM_1Kb_new <- TF_MYC_cpm1kb_gr$cpm1kb[olap$subjectHits]
TF_MYC_CPM_1Kb_new <- olap %>% group_by(queryHits) %>% summarise(TF_MYC_CPM_1Kb_new=mean(TF_MYC_CPM_1Kb_new))
dat3$TF_MYC_CPM_1Kb_new <-0
dat3$TF_MYC_CPM_1Kb_new[TF_MYC_CPM_1Kb_new$queryHits]<-TF_MYC_CPM_1Kb_new$TF_MYC_CPM_1Kb_new

## deltagb and deltagh
dat3$protospacer<- str_split(dat3$grna,"_",simplify = T)[,2]
protospacer <- dat3$protospacer
strand = dat3$strand
position2 <- matrix(NA,ncol=19,nrow=length(protospacer)) %>% as.data.frame
for(i in 1:19){
  position2[,i] <- as.factor(substr(protospacer,i,i+1))
}
colnames(position2) <- paste("position",seq_len(19),sep = "_")

delta <- c(-1,-2.1,-1.8,-0.9,-0.9,-2.1,-1.7,-0.9,-1.3,-2.7,-2.9,-1.1,-0.6,-1.5,-1.6,-0.2)
names(delta)<-c("TT","GT","CT","AT","TG","GG","CG","AG","TC","GC","CC","AC","TA","GA","CA","AA")
weight= c(1.80, 1.96, 1.90, 2.13, 1.38, 1.46, 1.00, 1.39, 1.51, 1.98, 1.88, 1.72, 2.02, 1.93, 2.08, 1.94, 2.15, 2.04, 2.25)
deltagh<- sapply(1:length(protospacer), function(i){sum(weight * delta[position2[i,] %>% unname() %>% unlist()] )})
dat3$deltagh <-deltagh

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
dat3$deltagb <-deltagb

write.table(dat3,file="/proj/milovelab/mu/dukeproj/data/scOct4th/sc_full.csv",row.names = F)

dat <- read.table(file="/proj/milovelab/mu/dukeproj/data/scOct4th/sc_full.csv", header = T)
# dat %>%
#   select(grna.chr, grna.start, grna.end, grna, dhs, grna.strand) %>%
#   write_tsv("/proj/milovelab/mu/dukeproj/data/scOct4th/sc_coordinates.tsv", col_names = F)


library(data.table)
library(tidyverse)
library(DescTools)
dir2 = "/proj/milovelab/mu/dukeproj/data/scOct4th/"
OGEE = fread("/proj/yunligrp/users/tianyou/gRNA/OGEE/sc_OGEE.txt")
dat = dat %>% 
  left_join(OGEE %>% select(grna, distance, OGEE_prop_Essential), by = "grna") 


data_new = dat %>%
  filter(pval_fdr_corrected > 0.2 | pval_fdr_corrected <= 0.05) %>% 
  mutate(significance = ifelse(pval_fdr_corrected<=0.05,1,0))

# for (i in c(17:23)){
#   print(colnames(data_new)[i])
#   data_new[,i] = Winsorize(data_new[,i], probs = c(0,0.95), na.rm = T)
# }
write_tsv(data_new, file.path(dir2, "sc_full_filtered.tsv"))

prop.table(table(data_new$chrom,data_new$significance),1)
# set.seed(2021)
# trainid <- createDataPartition(dat_sel_binary$significant, p=0.8, list=FALSE)
# trainid = createDataPartition(dat_sel_binary$direction_sig, p=0.8, list=FALSE)
freq <- table(data_new$chrom)/nrow(data_new)*100
freq2 <- sample(freq)
cumsum(freq2) # split the fold based on cumsum
split <- list("a" = c("chr12","chr6","chr21","chr15","chr5","chr22"), "b" = c("chr16","chr10","chr2","chr4","chr18"),
        "c" =c("chr1","chr7","chrX"),"d"= c("chr8","chr9", "chr20", "chr14", "chr3"),
        "e" = c("chr13", "chr19", "chr17", "chr11"))

for(i in 1:5){
      testid <- data_new$grna.chrom %in% split[[i]]
      train = data_new[!testid,]
      test = data_new[testid,]

      write_csv(train, paste0("data/scOct4th/sc_fdr_filtered-binary-", i, "-train.csv"))
      write_csv(test, paste0("data/scOct4th/sc_fdr_filtered-binary-", i, "-test.csv"))
}



#######################################################################################

library(Hmisc)
library("viridis")  
a <- data2 %>% mutate(total = cut2(gene_count_wo_grna, g=10,digits = 2), target = cut2(gene_count_w_grna, g=10)) %>% 
  group_by(total, target) %>% 
  dplyr::summarize(frac_sig = mean(pval_fdr_corrected < 0.05), number = n()) %>% 
  ggplot(aes(total, target, fill=frac_sig)) + 
  theme_grey(base_size=7)+ theme(axis.ticks=element_line(),
                                 panel.grid.major = element_blank(),
                                 # panel.grid.minor = element_blank(),
                                 # plot.background=element_blank(), 
                                 # panel.border=element_blank(),
                                 plot.title = element_text(color="black",hjust = 0),
                                 axis.text.x=element_text(angle=90, size = 10,hjust = 0.5, vjust = 1),
                                 axis.text.y=element_text(size = 10))+
  geom_tile()+
  scale_fill_viridis(option="C")+
  geom_point(aes(size=number))

ggplot(data2%>% filter(total<3000),aes(x= total, y=log10(pval_fdr_corrected)))+
  geom_point()+
  geom_hline(yintercept=log10(0.05), linetype="dashed", color = "red")

# For association analysis
## move grna at promoter region
neg <-cut2(sc$avg_logFC[which(sc$avg_logFC<0)],g=3,onlycuts = T,digits = 1)
pos <-cut2(sc$avg_logFC[which(sc$avg_logFC>=0)],g=3,onlycuts = T)

sc<- sc %>% mutate(effect = cut2(avg_logFC, c(neg[2:3],0,pos[-1])), 
                   pvalue = cut2(pval_fdr_corrected, c(0.05,0.1,0.2,0.5,1)),
                   prom = ifelse(distance_to_tss_of_linked_gene<=500,1,0)) 

sc_noprom <-sc %>% filter(prom==0) %>% 
  unite(pvalue_effect_class, pvalue, effect, sep = "_",remove = F)
## categorized grna-gene pairs
sc %>% group_by(effect, pvalue) %>% 
  dplyr::summarize(number = n(),prom_number = sum(prom)) %>% 
  ggplot(aes(effect, pvalue, fill=prom_number)) + 
  theme_grey(base_size=7)+ theme(axis.ticks=element_line(),
                                 panel.grid.major = element_blank(),
                                 axis.title = element_text(size = 10),
                                 # panel.grid.minor = element_blank(),
                                 # plot.background=element_blank(), 
                                 # panel.border=element_blank(),
                                 plot.title = element_text(color="black",hjust = 0),
                                 axis.text.x=element_text(angle=45, size = 10,hjust = 0.5, vjust = 0.5),
                                 axis.text.y=element_text(size = 10))+
  geom_tile()+
  scale_fill_viridis(option="C")+
  geom_point(aes(size=number),color="gray")

## overlap with hichip data
hic <- read_tsv("/proj/milovelab/mu/dukeproj/K562-results-julien/gRNA_prediction/K562_combined_r1r2r3.10k.2.peaks.bedpe")
library(GenomicRanges)
library(InteractionSet)
# Converting data.frame to GInteraction
convertToGI <- function(df,metadata){
  row.regions <- GRanges(df$chr1, IRanges(df$start1,df$end1))# interaction start
  col.regions <- GRanges(df$chr2, IRanges(df$start2,df$end2))# interaction end
  gi <- GInteractions(row.regions, col.regions)
  mcols(gi) <- df[,metadata] # Interaction frequencies
  return(gi)
}

convertToGI2 <- function(df,metadata){
  row.regions <- GRanges(df$dhs_chrom, IRanges(df$dhs_start,df$dhs_end))# interaction start
  col.regions <- GRanges(df$chrom, IRanges(df$start,df$end),strand = df$strand)# interaction gene
  col.regions <- promoters(col.regions,downstream = 500, upstream = 1500)
  gi <- GInteractions(row.regions, col.regions)
  mcols(gi) <- df[,metadata] # Interaction frequencies
  return(gi)
}
hic.gi <- convertToGI(hic,metadata=7:14)


# mutate(group_indices = group_indices(., pvalue, prom))
sc.gi <- convertToGI2(sc_noprom,metadata=c(1:9,18:22))

# olap <- findOverlaps(hic.gi, sc.gi)
# olap

sc_noprom$overlap_number <- countOverlaps(sc.gi,hic.gi)
sc_noprom2 <- sc_noprom %>% mutate(overlap = ifelse(overlap_number>0,1,0)) %>% 
  group_by(effect, pvalue) %>% 
  dplyr::summarize(number = n(),overlap_total = sum(overlap_number),overlap_number = sum(overlap)) %>% 
  dplyr::mutate(count_overlap = overlap_total / number, rate_overlap = overlap_number / number) 
plt_rate <- sc_noprom2 %>% ggplot(aes(effect, pvalue, fill=rate_overlap)) + 
  theme_grey(base_size=7)+ theme(axis.ticks=element_line(),
                                 panel.grid.major = element_blank(),
                                 axis.title = element_text(size = 10),
                                 # panel.grid.minor = element_blank(),
                                 # plot.background=element_blank(), 
                                 # panel.border=element_blank(),
                                 plot.title = element_text(color="black",hjust = 0),
                                 axis.text.x=element_text(angle=45, size = 10,hjust = 0.5, vjust = 0.5),
                                 axis.text.y=element_text(size = 10))+
  geom_tile()+
  scale_fill_viridis(option="C")+
  geom_point(aes(size=number),color="gray")
plt_rate
plt_count <- sc_noprom2 %>% ggplot(aes(effect, pvalue, fill=count_overlap)) + 
  theme_grey(base_size=7)+ theme(axis.ticks=element_line(),
                                 panel.grid.major = element_blank(),
                                 axis.title = element_text(size = 10),
                                 # panel.grid.minor = element_blank(),
                                 # plot.background=element_blank(), 
                                 # panel.border=element_blank(),
                                 plot.title = element_text(color="black",hjust = 0),
                                 axis.text.x=element_text(angle=45, size = 10,hjust = 0.5, vjust = 0.5),
                                 axis.text.y=element_text(size = 10))+
  geom_tile()+
  scale_fill_viridis(option="C")+
  geom_point(aes(size=number),color="gray")
plt_count
