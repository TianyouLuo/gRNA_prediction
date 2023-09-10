library(vroom)
library(GenomicRanges)
library(stringr)
library(tidyverse)
library(viridis)
library(ggpubr)
# Alex <- vroom("/pine/scr/t/i/tianyou/Patrick/single-cell/single-cell-data.csv")
Alex <- vroom("/proj/milovelab/mu/dukeproj/data/scOct4th/grna.de.markers.nFeature_RNA.with_grna_coords.txt.gz")
Alex$gRNAid = str_split(Alex$grna,"_",simplify = T)[,1]
Alex2 <-Alex %>% group_by(gRNAid) %>% filter(p_val== min(p_val)) 
Alex2$log10pfdr <- -log10(Alex2$pval_fdr_corrected)
Alex2$log10padj <- -log10(Alex2$p_val_adj)

enh_nosig <- vroom("/proj/milovelab/mu/dukeproj/data/dat_April_resplit/wgCERES-gRNAs-k562-discovery-screen-enh_baseMean125_nonsig.csv")
hist(enh_nosig$log2FoldChange,breaks = 50)

enh_sig <- vroom("/proj/milovelab/mu/dukeproj/data/dat_April_resplit/wgCERES-gRNAs-k562-discovery-screen-enh_baseMean125_sig.csv")
hist(enh_sig$log2FoldChange,breaks = 50)

pro_nosig <- vroom("/proj/milovelab/mu/dukeproj/data/dat_April_resplit/wgCERES-gRNAs-k562-discovery-screen-pro_baseMean125_nonsig.csv")
hist(pro_nosig$log2FoldChange,breaks = 50)
nosig_Alex_pro <- Alex2 %>% inner_join(pro_nosig,by="gRNAid")
hist(nosig_Alex_pro$avg_logFC,breaks = 50)
hist(abs(nosig_Alex_pro$avg_logFC),breaks = 50)
quantile(abs(nosig_Alex_pro$avg_logFC))

pro_sig <- vroom("/proj/milovelab/mu/dukeproj/data/dat_April_resplit/wgCERES-gRNAs-k562-discovery-screen-pro_baseMean125_sig.csv")
hist(pro_sig$log2FoldChange,breaks = 50)
pro_sig2 <-pro_sig %>% filter(padj<0.05)
sig_Alex_pro <- Alex2 %>% inner_join(pro_sig,by="gRNAid")
hist(sig_Alex_pro$avg_logFC,breaks = 50)
hist(abs(sig_Alex_pro$avg_logFC),breaks = 50)
quantile(abs(sig_Alex_pro$avg_logFC))

sig <- rbind(enh_sig,pro_sig)
nosig <- rbind(enh_nosig,pro_nosig)

## Draw venn diagram
library(data.table)
library(tidyverse)
dir2 = "/proj/milovelab/mu/dukeproj/data/mhc/ipsc"
ipsc <- read_csv(file.path(dir2, "ipsc-rawp0.05.csv"))
nrow(ipsc) # 9570
set.seed(2050)
dhs_ipsc <-ipsc$dhs %>% unique() %>% sort() #580

dir2 = "/proj/milovelab/mu/dukeproj/data/mhc/k562"
k562 <- read_csv(file.path(dir2, "k562-rawp0.05.csv"))
nrow(k562) # 9550
set.seed(2050)
dhs_k562 <-k562$dhs %>% unique() %>% sort() 
length(unique(dhs_k562)) # 581

dir2 = "/proj/milovelab/mu/dukeproj/data/mhc/npc"
npc <- read_csv(file.path(dir2, "npc-rawp0.05.csv"))
nrow(npc) # 9336
set.seed(2050)
dhs_npc <-npc$dhs %>% unique() %>% sort() 
length(unique(dhs_npc)) # 581

### 3 cell types overlap
mhc_all <- npc %>% inner_join(k562,by="grna")  |> inner_join(ipsc,by="grna")
nrow(mhc_all) # 6846
length(unique(mhc_all$dhs)) 

### k562 and npc
k562_npc <- npc %>% inner_join(k562,by="grna")  |> anti_join(ipsc,by="grna")
nrow(k562_npc) # 7973
length(unique(k562_npc$dhs.x)) # 581

### k562 and ipsc
k562_ipsc <- k562 %>% inner_join(ipsc,by="grna")  |> anti_join(npc,by="grna")
nrow(k562_ipsc) # 8165
length(unique(k562_ipsc$dhs.x)) # 580

### ipsc and npc
npc_ipsc <- ipsc %>% inner_join(npc,by="grna")  |> anti_join(k562,by="grna")
nrow(npc_ipsc) # 8165
length(unique(k562_ipsc$dhs.x)) # 580

### ipsc alone
ipsc_alone <- ipsc %>% anti_join(npc,by="grna")  |> anti_join(k562,by="grna")
nrow(ipsc_alone) # 8165

### npc alone
npc_alone <- npc %>% anti_join(ipsc,by="grna")  |> anti_join(k562,by="grna")
nrow(npc_alone) # 204

### k562 alone
k562_alone <- k562 %>% anti_join(ipsc,by="grna")  |> anti_join(npc,by="grna")
nrow(k562_alone) # 258

table5 <- read_table2("/proj/milovelab/mu/dukeproj/K562-supp-tables/supplementary_table_5_DHS_summary_results.txt.gz")
pro <- read_csv("data/dat_April_resplit/wgCERES-gRNAs-k562-discovery-screen-pro_baseMean125.csv")
dim(pro)
pro_binary = pro %>%
  filter(padj > 0.2 | padj <= 0.05) %>% 
  mutate(significance = ifelse(padj<=0.05,1,0))
dim(pro_binary) #280080
length(unique(pro_binary$DHS))
head(pro_binary$chromEnd)
pro_binary_dhs = pro_binary  |> select("DHS")  |> inner_join(table5, by=c("DHS"="name"))
dim(pro_binary_dhs)
length(unique(pro_binary$DHS)) # 29678

pro_binary_dhs$dhs = paste(pro_binary_dhs$chrom, ":", pro_binary_dhs$chromStart, "-", pro_binary_dhs$chromEnd, sep = "")
eng <- read_csv("data/dat_April_resplit/wgCERES-gRNAs-k562-discovery-screen-enh_baseMean125.csv")
enh_binary = eng %>%
  filter(padj > 0.2 | padj <= 0.05) %>% 
  mutate(significance = ifelse(padj<=0.05,1,0))
dim(enh_binary) #733551
length(unique(enh_binary$DHS)) #79800
enh_binary_dhs = enh_binary  |> select("DHS")  |> inner_join(table5, by=c("DHS"="name"))
enh_binary_dhs$dhs = paste(enh_binary_dhs$chrom, ":", enh_binary_dhs$chromStart, "-", enh_binary_dhs$chromEnd, sep = "")
pro_binary$region = "promoter"
enh_binary$region = "enhancer"

### promoter vs enhancer 
overlap_values <- intersect(unique(enh_binary$DHS), unique(pro_binary$DHS))
### sc data
sc <- read.table(file="/proj/milovelab/mu/dukeproj/data/scOct4th/sc_full.csv", header = T)
sc_binary = sc %>%
  filter(pval_fdr_corrected > 0.2 | pval_fdr_corrected <= 0.05) %>% 
  mutate(significance = ifelse(pval_fdr_corrected<=0.05,1,0))
dim(sc_binary) # 2597
head(sc_binary) # 2597
sc_binary$gRNAid = str_split(sc_binary$grna,"_",simplify = T)[,1]
sc_binary$dhs = paste(sc_binary$dhs_chrom, ":", sc_binary$dhs_start, "-", sc_binary$dhs_end, sep = "")
head(sc_binary$gRNAid)
pro_sc <- sc_binary %>% inner_join(pro_binary,by="gRNAid")
dim(pro_sc) #153
length(unique(pro_sc$gRNAid)) #153
enh_sc <- sc_binary %>% inner_join(enh_binary,by="gRNAid")
dim(enh_sc) #1072
length(unique(enh_sc$gRNAid)) #1072
sc_only <- sc_binary %>% anti_join(enh_binary,by="gRNAid")  |>  anti_join(pro_binary,by="gRNAid")
dim(sc_only) # 1372

### Read in MHC region
mhc <- read.table(file="/proj/milovelab/mu/dukeproj/data/mhc/ipsc/ipsc_onegene_full.csv", header = T)
mhc_binary = mhc %>%
  filter(pval_fdr_corrected <= 0.05 | pval_fdr_corrected >= 0.2) %>%
  mutate(significant = ifelse(pval_fdr_corrected <= 0.05, 1, 0))
dim(mhc_binary) #9570
head(mhc_binary)
mhc_sc <- mhc_binary %>% inner_join(sc_binary,by="dhs")
dim(mhc_sc)
mhc_sc <- mhc_binary %>% inner_join(enh_binary_dhs,by="dhs")
dim(mhc_sc)
mhc_dhs <-mhc_binary$dhs %>% unique() %>% sort()
length(mhc_dhs)
### read in abundance
count <- read_csv("/proj/yunligrp/users/tianyou/gRNA/data/Maria_WT_count/Plasmid_WT_count_simplified_20220812.csv")
count_dhs <- count |> group_by(DHS) |> summarise(n= n())
count_dhs |> dim() # total 373
dhs_enh_count <- intersect(unique(enh_binary$DHS), count_dhs$DHS)
length(dhs_enh_count) # 350
dhs_pro_count <- intersect(unique(pro_binary$DHS), count_dhs$DHS)
length(dhs_pro_count) # 22

dhs_sc_count <- intersect(unique(sc_binary$dhs_id), count_dhs$DHS)
length(dhs_sc_count) # 164
dhs_pro_sc_count <- intersect(unique(pro_sc$dhs_id), count_dhs$DHS)
length(dhs_pro_sc_count) # 1

dhs_enh_sc_count <- intersect(unique(enh_sc$dhs_id), count_dhs$DHS)
length(dhs_enh_sc_count) # 88

sc_only_count <- count %>% inner_join(sc_binary,by="gRNAid") |> anti_join(enh_binary,by="gRNAid")  |> anti_join(pro_binary, by = "gRNAid")
dhs_sc_only_count <- intersect(unique(sc_only_count$DHS),unique(count_dhs$DHS))
length(dhs_sc_only_count)
dhs_sc_count <- intersect(unique(sc_binary$dhs_id), count_dhs$DHS)
length(dhs_sc_count) # 164

grnaid <- str_split(count$gRNA_label,"__",simplify = TRUE)[,2]
count$gRNAid = paste(count$DHS,grnaid,sep = ".")
head(count$gRNAid)
dim(count) #9875
pro_count <- count %>% inner_join(pro_binary,by="gRNAid")
dim(pro_count) #70
pro_count$DHS.x |> unique() |> length() #21
enh_count <- count %>% inner_join(enh_binary,by="gRNAid")
dim(enh_count) #2450
enh_count$DHS.x |> unique() |> length() #345
count_only <- count %>% anti_join(enh_binary,by="gRNAid")  |>  anti_join(pro_binary,by="gRNAid") |>  anti_join(sc_binary,by="gRNAid")
dim(count_only) # 7316
count_only$DHS |> unique() |> length() #298

sc_count <- count %>% inner_join(sc_binary,by="gRNAid")
dim(sc_count) # 105
sc_pro_count <- count %>% inner_join(sc_binary,by="gRNAid") %>% inner_join(pro_binary,by="gRNAid")
dim(sc_pro_count) # 105
sc_enh_count <- count %>% inner_join(sc_binary,by="gRNAid") %>% inner_join(enh_binary,by="gRNAid")
dim(sc_enh_count) # 66
sc_enh_count$DHS.x |> unique() |> length() #66

enh_nosc_count <- count %>% inner_join(enh_binary,by="gRNAid") %>% anti_join(sc_binary,by="gRNAid")
dim(enh_nosc_count) # 2384
enh_nosc_count$DHS.x |> unique() |> length() # 343

## enrichment of discovery and single cell
nosig_Alex <- Alex2 %>% inner_join(nosig,by="gRNAid")
sig_Alex <- Alex2 %>% inner_join(sig,by="gRNAid")

library(ggplot2)
library(ggsci)
p1 <- ggplot(data.frame(value = c(nosig_Alex$avg_logFC,sig_Alex$avg_logFC),
                  status = c(rep("nosig",length(abs(nosig_Alex$avg_logFC))),rep("sig",length(abs(sig_Alex$avg_logFC)))))
) +
  geom_density(aes(value,fill=status,color=status),alpha=0.4)+
  labs(x = "avg_logFC")+
  scale_fill_viridis(discrete = TRUE, alpha=0.3) +
  coord_cartesian(xlim = c(0,1), expand = TRUE) +
  # theme_bw()+
  theme(
    panel.background = element_rect(color=NA, fill = "white"),
    panel.border = element_rect(),
    # panel.grid.major = element_line(size = 0.5, linetype = "solid", colour = "black"),
    # axis.line = element_line(size = 1, linetype = "solid", colour = "black"),
    text = element_text(size = 20),
    legend.text = element_text(size = 24),
    legend.background = element_rect(color = NA, fill = NA),
    legend.position = c(0.95,0.95),
    legend.justification = c("right", "top"))
  # scale_color_igv()

p2 <- ggplot(data.frame(value = c(nosig_Alex$log10pfdr,sig_Alex$log10pfdr),
                  status = c(rep("nosig",length(abs(nosig_Alex$log10pfdr))),rep("sig",length(abs(sig_Alex$log10pfdr)))))
) +
  geom_density(aes(value,fill=status,color=status),alpha=0.4)+
  labs(x = "-log10(pvalue)")+
  scale_fill_viridis(discrete = TRUE, alpha=0.3) +
  coord_cartesian(xlim = c(0,1), expand = TRUE) +
  # theme_bw()+
  theme(
    panel.background = element_rect(color=NA, fill = "white"),
    panel.border = element_rect(),
    # panel.grid.major = element_line(size = 0.5, linetype = "solid", colour = "black"),
    # axis.line = element_line(size = 1, linetype = "solid", colour = "black"),
    text = element_text(size = 20),
    legend.text = element_text(size = 24),
    legend.background = element_rect(color = NA, fill = NA),
    legend.position = c(0.95,0.95),
    legend.justification = c("right", "top"))

dat <- data.frame(value = c(nosig_Alex$pval_fdr_corrected,sig_Alex$pval_fdr_corrected),
                  Groups = c(rep("insignificant",length(abs(nosig_Alex$pval_fdr_corrected))),rep("significant",length(abs(sig_Alex$pval_fdr_corrected)))))


dat |> filter(Groups == "significant") |> mutate(FDR0.1 = ifelse(value <= 0.15,1,0)) |> 
    group_by(FDR0.1) |> summarise(counts = n())
dat |> filter(Groups == "insignificant") |> mutate(FDR0.1 = ifelse(value <= 0.15,1,0)) |> 
    group_by(FDR0.1) |> summarise(counts = n())

p3 <- ggplot(data.frame(value = c(nosig_Alex$pval_fdr_corrected,sig_Alex$pval_fdr_corrected),
                  Groups = c(rep("insignificant",length(abs(nosig_Alex$pval_fdr_corrected))),rep("significant",length(abs(sig_Alex$pval_fdr_corrected)))))
) +
  geom_density(aes(value,fill=Groups,col=Groups),alpha=0.4) +
  labs(x = "adjusted p-value")+
  scale_fill_viridis(discrete = TRUE, alpha=0.3) +
  scale_color_viridis(discrete = TRUE) +
  coord_cartesian(xlim = c(0,1), expand = TRUE) +
  theme_bw()+
  theme(
    panel.background = element_rect(fill = "white"),
    panel.border = element_rect(fill = NA, color = "black", size = 2),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    # panel.grid.major = element_line(size = 0.5, linetype = "solid", colour = "black"),
    # axis.line = element_line(size = 1, linetype = "solid", colour = "black"),
    text = element_text(size = 20),
    legend.text = element_text(size = 24),
    legend.background = element_rect(color = NA, fill = NA),
    legend.position = c(0.75,0.95),
    legend.justification = c("right", "top"))

### look at effect # of genes
nosig_Alex <- Alex %>% inner_join(nosig,by=c("gRNAid"))%>% mutate(sig = pval_fdr_corrected<0.2)
nosig_number_gene<- nosig_Alex %>% group_by(gRNAid) %>% summarise(n=sum(sig))
nosigsummary<-nosig_number_gene %>%summarise(cut1off = n>0) %>%  group_by(cut1off) %>% tally()

sig_Alex <- Alex %>% inner_join(sig,by=c("gRNAid"))%>% mutate(sig = pval_fdr_corrected<0.2)
sig_number_gene<- sig_Alex %>% group_by(gRNAid) %>%  summarise(n=sum(sig))

sigsummary<-sig_number_gene %>%summarise(cut1off = n>0) %>%  group_by(cut1off) %>% tally()
data = data.frame(value = c(nosig_number_gene$n,sig_number_gene$n),
                  Groups = c(rep("insignificant",nrow(nosig_number_gene)),rep("significant",nrow(sig_number_gene))))
p4 <- ggplot( data)+
  geom_histogram(aes(value,fill=Groups,col=Groups),position = "stack", alpha=0.4)+
  labs(x = "# of significant genes")+
  scale_fill_viridis(discrete = TRUE, alpha=0.3) +
  scale_color_viridis(discrete = TRUE) +
  theme_bw()+
  theme(
    panel.background = element_rect(fill = "white"),
    panel.border = element_rect(fill = NA, color = "black", size = 2),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    # panel.grid.major = element_line(size = 0.5, linetype = "solid", colour = "black"),
    # axis.line = element_line(size = 0, linetype = "solid", colour = "black"),
    axis.line = element_blank(),
    text = element_text(size = 20),
    legend.text = element_text(size = 24),
    legend.background = element_rect(color = NA, fill = NA),
    legend.position = c(0.75,0.95),
    legend.justification = c("right", "top"))

library(patchwork)
pdf(file="plots/enrichment1.pdf",height = 185/25.4, width = 250/25.4)
p3 + p4
dev.off()

ggsave("plots/enrichment1.jpg", plot = p3 + p4, width = 20, height = 10, dpi = 450)

table(data)
data1 <-matrix(c(413,586,144,295),nrow = 2) #0.2
chisq.test(data1) # p-value = 0.002677
fisher.test(data1)
## one or more significant gene
data1 <-matrix(c(502,769,55,112),nrow = 2) #0.05
chisq.test(data1) # 0.12
fisher.test(data1)

data1 <-matrix(c(543,840,14,41),nrow = 2) #0.05
chisq.test(data1) # p-value = 0.002677
fisher.test(data1)
# larger than 3 gene 
data1 <-matrix(c(555,860,2,21),nrow = 2) #0.05
chisq.test(data1)
fisher.test(data1) # p-value = 0.001996

data1 <-matrix(c(484,728,73,153),nrow = 2) #0.05
chisq.test(data1) # p-value = 0.00267

## look at mediation effect



## read validation screen
valid <- vroom("/proj/milovelab/mu/dukeproj/K562-supp-tables/supplementary_table_6_validation_screen_k562_sgrna_deseq2_results_hg19.csv.gz")
valid <- vroom("/proj/milovelab/mu/dukeproj/K562-supp-tables/validationv2_k562_sgrna_deseq2_results_hg19.csv")
valid$guideID = paste(valid$DHS,str_split(valid$gRNA_label,"_",simplify = T)[,4],sep=".")

## new valid
nosig_valid <- valid %>% inner_join(nosig,by=c("guideID"="gRNAid")) %>% filter(pvalue.y>0.05) ## 462 nosig in discovery

nosig_Alex <- Alex2 %>% inner_join(nosig_valid,by=c("gRNAid"="guideID"))
hist(nosig_Alex$avg_logFC,breaks = 50)
hist(abs(nosig_Alex$avg_logFC),breaks = 50)
quantile(abs(nosig_Alex$avg_logFC))

## 2909 gRNAs are sig in both valid and discover
sig_valid <- valid %>% inner_join(sig,by=c("guideID"="gRNAid")) %>% filter(pvalue.y<0.05& pvalue.x<0.05 & sign(log2FoldChange.x) == sign(log2FoldChange.y))
sig_Alex <- Alex2 %>% inner_join(sig_valid,by=c("gRNAid"="guideID"))
hist(sig_Alex$avg_logFC,breaks = 50)
hist(abs(sig_Alex$avg_logFC),breaks = 50)
quantile(abs(sig_Alex$avg_logFC))

library(ggplot2)
library(ggsci)
p1 <- ggplot(data.frame(value = c(nosig_Alex$avg_logFC,sig_Alex$avg_logFC),
                  status = c(rep("nosig",length(abs(nosig_Alex$avg_logFC))),rep("sig",length(abs(sig_Alex$avg_logFC)))))
) +
  geom_density(aes(value,fill=status,color=status),alpha=0.4)+
  labs(x = "avg_logFC")+
  theme_bw()+
  scale_fill_igv()+
  scale_color_igv()

p2 <- ggplot(data.frame(value = c(nosig_Alex$log10pfdr,sig_Alex$log10pfdr),
                  status = c(rep("nosig",length(abs(nosig_Alex$log10pfdr))),rep("sig",length(abs(sig_Alex$log10pfdr)))))
) +
  geom_density(aes(value,fill=status,color=status),alpha=0.4)+
  labs(x = "-log10(pvalue)")+
  theme_bw()+
  scale_fill_igv()+
  scale_color_igv()

p3 <- ggplot(data.frame(value = c(nosig_Alex$pval_fdr_corrected,sig_Alex$pval_fdr_corrected),
                  status = c(rep("nosig",length(abs(nosig_Alex$pval_fdr_corrected))),rep("sig",length(abs(sig_Alex$pval_fdr_corrected)))))
) +
  geom_density(aes(value,fill=status,col=status),alpha=0.4) +
  labs(x = "p_fdr_corrected")+
  theme_bw()+
  scale_fill_igv()+
  scale_color_igv()

library(patchwork)
p1+p2+p3+p4
### look at effect # of genes
nosig_Alex <- Alex %>% inner_join(nosig_valid,by=c("gRNAid"="guideID"))%>% mutate(sig = pval_fdr_corrected<0.2)
nosig_number_gene<- nosig_Alex %>% group_by(gRNAid) %>% summarise(n=sum(sig))
nosigsummary<-nosig_number_gene %>%summarise(cut1off = n>0) %>%  group_by(cut1off) %>% tally()

sig_Alex <- Alex %>% inner_join(sig_valid,by=c("gRNAid"="guideID"))%>% mutate(sig = pval_fdr_corrected<0.2)
sig_number_gene<- sig_Alex %>% group_by(gRNAid) %>%  summarise(n=sum(sig))

sigsummary<-sig_number_gene %>%summarise(cut1off = n>0) %>%  group_by(cut1off) %>% tally()
data = data.frame(value = c(nosig_number_gene$n,sig_number_gene$n),
                  status = c(rep("nosig",nrow(nosig_number_gene)),rep("sig",nrow(sig_number_gene))))
p4 <- ggplot( data)+
  geom_histogram(aes(value,fill=status,col=status),alpha=0.4)+
  labs(x = "# of significant genes")+
  theme_bw()+
  scale_fill_igv()+
  scale_color_igv()

data <-matrix(c(335,713,127,350),nrow = 2)
chisq.test(data)

data <-matrix(c(335,713,127,350),nrow = 2)
chisq.test(data)

ggplot(data, aes(x = value, fill = status)) + 
  geom_histogram(position = "identity",alpha=0.4)


### compare valid and discovery
data = data.frame(valid_LFC = c(nosig_valid_enh$log2FoldChange.x,sig_valid_enh$log2FoldChange.x),
                  discover_LFC = c(nosig_valid_enh$log2FoldChange.y,sig_valid_enh$log2FoldChange.y),
                  valid_padj=c(nosig_valid_enh$padj.x,sig_valid_enh$padj.x),
                  discover_padj=c(nosig_valid_enh$padj.y,sig_valid_enh$padj.y),
                  status = c(rep("nosig",length(abs(nosig_valid_enh$log2FoldChange.x))),rep("sig",length(abs(sig_valid_enh$log2FoldChange.x)))))
plot(x=data$discover_LFC,y=data$valid_LFC,col=c(rep("gray",length(abs(nosig_valid_enh$log2FoldChange.x))),rep("red",length(abs(sig_valid_enh$log2FoldChange.x)))))
identify(data$discover_LFC, data$valid_LFC, labels=row.names(data)) # identify points

ggplot(data) +geom_point(aes(x=discover_LFC,y=valid_LFC,col=status))

nosig_valid_pro <- valid %>% inner_join(pro_nosig,by=c("guideID"="gRNAid"))
sig_valid_pro <- valid %>% inner_join(pro_sig,by=c("guideID"="gRNAid"))
data2 = data.frame(valid_LFC = c(nosig_valid_pro$log2FoldChange.x,sig_valid_pro$log2FoldChange.x),
                  discover_LFC = c(nosig_valid_pro$log2FoldChange.y,sig_valid_pro$log2FoldChange.y),
                  valid_padj=c(nosig_valid_pro$padj.x,sig_valid_pro$padj.x),
                  discover_padj=c(nosig_valid_pro$padj.y,sig_valid_pro$padj.y),
                  status = c(rep("nosig",length(abs(nosig_valid_pro$log2FoldChange.x))),rep("sig",length(abs(sig_valid_pro$log2FoldChange.x)))))
plot(x=data2$discover_LFC,y=data2$valid_LFC,col=c(rep("gray",length(abs(nosig_valid_pro$log2FoldChange.x))),rep("red",length(abs(sig_valid_pro$log2FoldChange.x)))))
identify(data2$discover_LFC, data2$valid_LFC, labels=row.names(data2)) # identify points

