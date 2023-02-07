library(tidyverse)
library(ggrepel)
library(ggpubr)
library(here)
library(ggrastr)

#read in data file
dataDir = "diff";
dataFile = "diff.significant.glm.ATACseq.txt";
data = read.delim(file=here(dataDir, dataFile), header=T, sep="\t");

#Output directory
output_dir <- "volcano_plots"

#Finding genes with significant expression levels to include
gene_cut = 1.5

rppm <- select(data, peak, ends_with(".rppm"))
rppm <- rppm[,c(1, grep("CM|EM|EMRA", names(rppm)))]
rppm$allmean <- rowMeans(rppm[,-1])
gene_df<- subset(rppm, allmean>=gene_cut)
gene_check <- as.character(gene_df$peak)

#Subsetting differential data
#data2 <- data[data$peak%in%gene_check,c(1,2,grep("(CD8|CD4)_(EM|CM|EMRA)_blood_unstim.v.\\1_(EM|CM|EMRA)_blood_stim", names(data)), 1860)]
#data2 <- data[data$peak%in%gene_check,]
data2 <- data

#Comparisons to Naive
##########################################################################################
# CD4_CM <- select(data2, peak, starts_with("CD4_CM_blood_unstim.v.CD4_CM_")) %>% select(peak, ends_with(".logFC"), ends_with(".fdr"), ends_with(".sig"))
# names(CD4_CM) <- c("gene", "logFC", "fdr", "sig")
# CD4_EM <- select(data2, peak, starts_with("CD4_EM_blood_unstim.v.CD4_EM_")) %>% select(peak, ends_with(".logFC"), ends_with(".fdr"), ends_with(".sig"))
# names(CD4_EM) <- c("gene", "logFC", "fdr", "sig")
# CD8_CM <- select(data2, peak, starts_with("CD8_Nav_blood_unstim.v.CD8_CM_"), starts_with("k.5.clus")) %>% select(peak, ends_with(".logFC"), ends_with(".fdr"), ends_with(".sig"), starts_with("k.5.clus"))
# names(CD8_CM)[1:4] <- c("gene", "logFC", "fdr", "sig")
# CD8_CM$k.5.cluster <- factor(CD8_CM$k.5.cluster)
# CD8_EM <- select(data2, peak, starts_with("CD8_Nav_blood_unstim.v.CD8_EM_"), starts_with("k.5.clus")) %>% select(peak, ends_with(".logFC"), ends_with(".fdr"), ends_with(".sig"), starts_with("k.5.clus"))
# names(CD8_EM)[1:4] <- c("gene", "logFC", "fdr", "sig")
# CD8_EM$k.5.cluster <- factor(CD8_EM$k.5.cluster)
# CD8_EMRA <- select(data2, peak, starts_with("CD8_Nav_blood_unstim.v.CD8_EMRA_"), starts_with("k.5.clus")) %>% select(peak, ends_with(".logFC"), ends_with(".fdr"), ends_with(".sig"), starts_with("k.5.clus"))
# names(CD8_EMRA)[1:4] <- c("gene", "logFC", "fdr", "sig")
# CD8_EMRA$k.5.cluster <- factor(CD8_EMRA$k.5.cluster)

###########################################################################################
#Intra EM EMRA plot
CD8_EM_EMRA <- select(data2, peak, starts_with("CD8_EMRA_blood_unstim.v.CD8_EM_")) %>% 
  select(peak, ends_with(".logFC"), ends_with(".fdr"), ends_with(".sig")) %>%
  drop_na()
names(CD8_EM_EMRA)[1:4] <- c("peak", "logFC", "fdr", "sig")
CD8_EM_EMRA <- CD8_EM_EMRA %>% mutate(neglog10fdr=-log(fdr, base=10))
CD8_EM_EMRAsig <- subset(CD8_EM_EMRA, sig==T)

#Intra CD8 EM CM plot
CD8_EM_CM <- select(data2, peak, starts_with("CD8_EM_blood_unstim.v.CD8_CM_"), starts_with("k.5.clus")) %>% 
  select(peak, ends_with(".logFC"), ends_with(".fdr"), ends_with(".sig")) %>%
  drop_na()
names(CD8_EM_CM)[1:4] <- c("peak", "logFC", "fdr", "sig")
CD8_EM_CM <- CD8_EM_CM %>% mutate(neglog10fdr=-log(fdr, base=10))
CD8_EM_CMsig <- subset(CD8_EM_CM, sig==T)

CD4_EM_CM <- select(data2, peak, starts_with("CD4_EM_blood_unstim.v.CD4_CM_"), starts_with("k.5.clus")) %>% 
  select(peak, ends_with(".logFC"), ends_with(".fdr"), ends_with(".sig")) %>%
  drop_na()
names(CD4_EM_CM)[1:4] <- c("peak", "logFC", "fdr", "sig")
CD4_EM_CM <- CD4_EM_CM %>% mutate(neglog10fdr=-log(fdr, base=10))
CD4_EM_CMsig <- subset(CD4_EM_CM, sig==T)

categorize <- function(x){
  for (i in 1:length(x$peak)){
    if(x$sig[i]==FALSE){
      x$sig_cat[i] = "Not significant"
    } else if (x$sig[i]==TRUE & x$logFC[i] >=0){
      x$sig_cat[i] = "EM" 
    } else {
      x$sig_cat[i] = "CM"
    }
  }
  return(x)
}

categorize_emra <- function(x){
  for (i in 1:length(x$peak)){
    if(x$sig[i]==FALSE){
      x$sig_cat[i] = "Not significant"
    } else if (x$sig[i]==TRUE & x$logFC[i] >=0){
      x$sig_cat[i] = "EMRA" 
    } else {
      x$sig_cat[i] = "EM"
    }
  }
  return(x)
}

CD8_EM_CM <- categorize(CD8_EM_CM)
CD4_EM_CM <- categorize(CD4_EM_CM)
CD8_EM_EMRA <- categorize_emra(CD8_EM_EMRA)

CD8_sig_counts <- table(CD8_EM_CMsig$logFC>=0)
CD4_sig_counts <- table(CD4_EM_CMsig$logFC>=0)
CD8_emra_sigcounts <- table(CD8_EM_EMRAsig$logFC>=0)
# lst <- list(CD8_CM=CD8_CM,
#             CD8_EM=CD8_EM,
#             CD8_EMRA=CD8_EMRA)

# rm(CD8_CM, CD8_EM, CD8_EMRA)
# gc()

#Adding new variables
# lst <- map(lst, function(x) mutate(x, neglogFC=-logFC, 
#                                    neglog10fdr=-log(fdr, base=10), 
#                                    rank=sign(neglogFC)*neglog10fdr
#                                    )
#           )

#Plotting cutoffs for gene labels
#first entry = graph 1, etc
# logFC_cut=c(2, 1, 7, 3, 3)
# logp_cut=c(2, 1.3, 10, 3, 3)
# logFC_topcut=c(1, 1, 3, 1, 1)
# logp_topcut=c(4, 1.3, 20, 10, 10)
# label_size=3

x=10
y=35

theme_set(theme_light())

p_CD8_EM_EMRA <- ggplot(CD8_EM_EMRA, aes(x=logFC, y=neglog10fdr))
p_CD8_EM_CM <- ggplot(CD8_EM_CM, aes(x=logFC, y=neglog10fdr))
p_CD4_EM_CM <- ggplot(CD4_EM_CM, aes(x=logFC, y=neglog10fdr))

p_CD8_EM_EMRA + rasterize(geom_point(aes(color=sig_cat), size=0.5, alpha=0.7)) +
  #geom_text_repel(aes(label=if_else((abs(logFC)>=logFC_cut[1]&abs(neglog10fdr)>=logp_cut[1])|((abs(logFC)>=logFC_topcut[1]&abs(neglog10fdr)>=logp_topcut[1])), as.character(gene), "")),min.segment.length = 0,size=label_size) +
  geom_hline(yintercept = 1.3, linetype="dashed") +
  geom_vline(xintercept = 1, linetype="dashed") +
  geom_vline(xintercept = -1, linetype="dashed") +
  xlab("log2(FC) (EMRA/EM)") +
  ylab("-log10(FDR)") +
  scale_color_manual(values=c("Not significant"="grey", "EM"="#E55100FF", "EMRA"="#870D4EFF")) +
  scale_y_continuous(limits=c(0, 35)) +
  scale_x_continuous(limits=c(-8, 8)) +
  labs(color="") +
  annotate("text", x = -5, y = 25, label = as.character(CD8_emra_sigcounts[[1]]), color="#E55100FF") +
  annotate("text", x = 5, y = 25, label = as.character(CD8_emra_sigcounts[[2]]), color="#870D4EFF")
ggsave(filename=here(output_dir, "CD8_EMRA_EM_volcano_raster.pdf"), height = 3, width=5)

p_CD8_EM_CM + rasterize(geom_point(aes(color=sig_cat), size=0.5, alpha=0.7)) +
  #geom_text_repel(aes(label=if_else((gene%in%gene_check & abs(logFC)>=logFC_cut[2]&abs(neglog10fdr)>=logp_cut[2])|(gene%in%gene_check & (abs(logFC)>=logFC_topcut[2]&abs(neglog10fdr)>=logp_topcut[2])), as.character(gene), "")),min.segment.length = 0,size=label_size) +
  geom_hline(yintercept = 1.3, linetype="dashed") +
  geom_vline(xintercept = 1, linetype="dashed") +
  geom_vline(xintercept = -1, linetype="dashed") +
  xlab("log2(FC) (EM/CM)") +
  ylab("-log10(FDR)") +
  scale_color_manual(values=c("Not significant"="grey", "EM"="#E55100FF", "CM"="#0C46A0FF")) +
  scale_y_continuous(limits=c(0, 35)) +
  scale_x_continuous(limits=c(-8, 8)) +
  labs(color="") +
  annotate("text", x = -5, y = 25, label = as.character(CD8_sig_counts[[1]]), color="#0C46A0FF") +
  annotate("text", x = 5, y = 25, label = as.character(CD8_sig_counts[[2]]), color="#E55100FF")
ggsave(filename=here(output_dir, "CD8_EMCM_volcano_stimunstim_raster.pdf"), height = 3, width=5)

p_CD4_EM_CM + rasterize(geom_point(aes(color=sig_cat), size=0.5)) +
  #geom_text_repel(aes(label=if_else((gene%in%gene_check & abs(logFC)>=logFC_cut[2]&abs(neglog10fdr)>=logp_cut[2])|(gene%in%gene_check & (abs(logFC)>=logFC_topcut[2]&abs(neglog10fdr)>=logp_topcut[2])), as.character(gene), "")),min.segment.length = 0,size=label_size) +
  geom_hline(yintercept = 1.3, linetype="dashed") +
  geom_vline(xintercept = 1, linetype="dashed") +
  geom_vline(xintercept = -1, linetype="dashed") +
  xlab("log2(FC) (EM/CM)") +
  ylab("-log10(FDR)") +
  scale_color_manual(values=c("Not significant"="grey", "EM"="#FFB74CFF", "CM"="#2096F2FF")) +
  scale_y_continuous(limits=c(0, 35)) +
  scale_x_continuous(limits=c(-8, 8)) +
  labs(color="") +
  annotate("text", x = -5, y = 25, label = as.character(CD4_sig_counts[[1]]), color="#2096F2FF") +
  annotate("text", x = 5, y = 25, label = as.character(CD4_sig_counts[[2]]), color="#FFB74CFF")
ggsave(filename=here(output_dir, "CD4_EMCM_volcano_stimunstim_raster.pdf"), height = 3, width=5)

#####################################################################################
#Interactive




#############################################################################
#CM
lst_CM_ggbase <- list()
lst_CM_ggplot <- list()
lst_CM_rnk <- list()

p_CD8_CM <- ggplot(lst$CD8_CM, aes(x=neglogFC, y=neglog10fdr))

for (i in 1:length(levels(lst$CD8_CM$k.5.cluster))){
  sub_tmp <- subset(lst$CD8_CM, k.5.cluster==i)
  lst_CM_ggbase[[i]] <- ggplot(sub_tmp, aes(x=neglogFC, y=neglog10fdr))
  lst_CM_rnk[[i]] <- select(sub_tmp, gene, rank)
  lst_CM_rnk[[i]]$gene <- toupper(lst_CM_rnk[[i]]$gene)
  names(lst_CM_rnk[[i]]) <- c("GENE", "rnk")
}

CM_Vol <- p_CD8_CM + geom_point(aes(color=sig)) +
    #geom_text_repel(aes(label=if_else((abs(neglogFC)>=logFC_cut[3]&abs(neglog10fdr)>=logp_cut[3])|((abs(neglogFC)>=logFC_topcut[3]&abs(neglog10fdr)>=logp_topcut[3])), as.character(gene), "")),min.segment.length = 0,size=label_size) +
    geom_hline(yintercept = 1.3, linetype="dashed") +
    geom_vline(xintercept = 1, linetype="dashed") +
    geom_vline(xintercept = -1, linetype="dashed") +
    xlab("log2(FC) (Mem/Nav)") +
    ylab("-log10(FDR)") +
    scale_color_manual(values=c("FALSE"="grey", "TRUE"="#0C46A0FF")) +
    facet_wrap(vars(k.5.cluster))

for (c in 1:length(lst_CM_ggbase)){
  
  lst_CM_ggplot[[c]] <- lst_CM_ggbase[[c]] + geom_point(aes(color=sig)) +
    geom_text_repel(aes(label=if_else((gene%in%gene_check & abs(neglogFC)>=logFC_cut[c]&abs(neglog10fdr)>=logp_cut[c])|(gene%in%gene_check & (abs(neglogFC)>=logFC_topcut[c]&abs(neglog10fdr)>=logp_topcut[c])), as.character(gene), "")),min.segment.length = 0,size=label_size) +
    geom_hline(yintercept = 1.3, linetype="dashed") +
    geom_vline(xintercept = 1, linetype="dashed") +
    geom_vline(xintercept = -1, linetype="dashed") +
    xlab("log2(FC) (Mem/Nav)") +
    ylab("-log10(FDR)") +
    scale_color_manual(values=c("FALSE"="grey", "TRUE"="#0C46A0FF")) +
    ggtitle(paste0("CD8 CM, Cluster ", c)) +
    scale_x_continuous(limits=c(-x,x)) +
    scale_y_continuous(limits=c(0,y))
  
}
# ggsave(filename=paste(dataDir, "CD8_CM_volcano_stimunstim.pdf",sep=""), height = 7, width=9)

######################################################################################
#EM

lst_EM_ggbase <- list()
lst_EM_ggplot <- list()
lst_EM_rnk <- list()

p_CD8_EM <- ggplot(lst$CD8_EM, aes(x=neglogFC, y=neglog10fdr))

for (i in 1:length(levels(lst$CD8_EM$k.5.cluster))){
  sub_tmp <- subset(lst$CD8_EM, k.5.cluster==i)
  lst_EM_ggbase[[i]] <- ggplot(sub_tmp, aes(x=neglogFC, y=neglog10fdr))
  lst_EM_rnk[[i]] <- select(sub_tmp, gene, rank)
  lst_EM_rnk[[i]]$gene <- toupper(lst_EM_rnk[[i]]$gene)
  names(lst_EM_rnk[[i]]) <- c("GENE", "rnk")
}

EM_Vol <- p_CD8_EM + geom_point(aes(color=sig)) +
    #geom_text_repel(aes(label=if_else((abs(neglogFC)>=logFC_cut[3]&abs(neglog10fdr)>=logp_cut[3])|((abs(neglogFC)>=logFC_topcut[3]&abs(neglog10fdr)>=logp_topcut[3])), as.character(gene), "")),min.segment.length = 0,size=label_size) +
    geom_hline(yintercept = 1.3, linetype="dashed") +
    geom_vline(xintercept = 1, linetype="dashed") +
    geom_vline(xintercept = -1, linetype="dashed") +
    xlab("log2(FC) (Mem/Naive)") +
    ylab("-log10(FDR)") +
    scale_color_manual(values=c("FALSE"="grey", "TRUE"="#E55100FF"))  +
    facet_wrap(vars(k.5.cluster))
# ggsave(filename=paste(dataDir, "CD8_EM_volcano_stimunstim.pdf",sep=""), height = 7, width=9)

for (c in 1:length(lst_EM_ggbase)){
  
  lst_EM_ggplot[[c]] <- lst_EM_ggbase[[c]] + geom_point(aes(color=sig)) +
    geom_text_repel(aes(label=if_else((gene%in%gene_check & abs(neglogFC)>=logFC_cut[c]&abs(neglog10fdr)>=logp_cut[c])|(gene%in%gene_check &(abs(neglogFC)>=logFC_topcut[c]&abs(neglog10fdr)>=logp_topcut[c])), as.character(gene), "")),min.segment.length = 0,size=label_size) +
    geom_hline(yintercept = 1.3, linetype="dashed") +
    geom_vline(xintercept = 1, linetype="dashed") +
    geom_vline(xintercept = -1, linetype="dashed") +
    xlab("log2(FC) (Mem/Nav)") +
    ylab("-log10(FDR)") +
    scale_color_manual(values=c("FALSE"="grey", "TRUE"="#E55100FF")) +
    ggtitle(paste0("CD8 EM, Cluster ", c))  +
    scale_x_continuous(limits=c(-x,x)) +
    scale_y_continuous(limits=c(0,y))
  
}

#########################################################################################
#EMRA

lst_EMRA_ggbase <- list()
lst_EMRA_ggplot <- list()
lst_EMRA_rnk <- list()

p_CD8_EMRA <- ggplot(lst$CD8_EMRA, aes(x=neglogFC, y=neglog10fdr))

for (i in 1:length(levels(lst$CD8_EMRA$k.5.cluster))){
  sub_tmp <- subset(lst$CD8_EMRA, k.5.cluster==i)
  lst_EMRA_ggbase[[i]] <- ggplot(sub_tmp, aes(x=neglogFC, y=neglog10fdr))
  lst_EMRA_rnk[[i]] <- select(sub_tmp, gene, rank)
  lst_EMRA_rnk[[i]]$gene <- toupper(lst_EMRA_rnk[[i]]$gene)
  names(lst_EMRA_rnk[[i]]) <- c("GENE", "rnk")
}

EMRA_Vol <- p_CD8_EMRA + geom_point(aes(color=sig)) +
  #geom_text_repel(aes(label=if_else((abs(neglogFC)>=logFC_cut&abs(neglog10fdr)>=logp_cut)|((abs(neglogFC)>=logFC_topcut&abs(neglog10fdr)>=logp_topcut)), as.character(gene), "")),min.segment.length = 0,size=label_size) +
  geom_hline(yintercept = 1.3, linetype="dashed") +
  geom_vline(xintercept = 1, linetype="dashed") +
  geom_vline(xintercept = -1, linetype="dashed") +
  xlab("log2(FC) (Mem/Naive)") +
  ylab("-log10(FDR)") +
  scale_color_manual(values=c("FALSE"="grey", "TRUE"="#870D4EFF")) +
  facet_wrap(vars(k.5.cluster))
# ggsave(filename=paste(dataDir, "CD8_EMRA_volcano_stimunstim.pdf",sep=""), height = 7, width=9)

for (c in 1:length(lst_EMRA_ggbase)){

lst_EMRA_ggplot[[c]] <- lst_EMRA_ggbase[[c]] + geom_point(aes(color=sig)) +
  geom_text_repel(aes(label=if_else((gene%in%gene_check & abs(neglogFC)>=logFC_cut[c]&abs(neglog10fdr)>=logp_cut[c]) |(gene%in%gene_check & (abs(neglogFC)>=logFC_topcut[c]&abs(neglog10fdr)>=logp_topcut[c])), as.character(gene), "")),min.segment.length = 0,size=label_size) +
  geom_hline(yintercept = 1.3, linetype="dashed") +
  geom_vline(xintercept = 1, linetype="dashed") +
  geom_vline(xintercept = -1, linetype="dashed") +
  xlab("log2(FC) (Mem/Nav)") +
  ylab("-log10(FDR)") +
  scale_color_manual(values=c("FALSE"="grey", "TRUE"="#870D4EFF")) +
  ggtitle(paste0("CD8 EMRA, Cluster ", c))  +
  scale_x_continuous(limits=c(-x,x)) +
  scale_y_continuous(limits=c(0,y))

}

c=1
ggarrange(lst_CM_ggplot[[1]], lst_EM_ggplot[[1]], lst_EMRA_ggplot[[1]], ncol=3, common.legend = T)
ggsave(filename = here(output_dir, "clus1_composite_vol.pdf"), height=5, width=15)
c=2
ggarrange(lst_CM_ggplot[[2]], lst_EM_ggplot[[2]], lst_EMRA_ggplot[[2]], ncol=3, common.legend = T)
ggsave(filename = here(output_dir, "clus2_composite_vol.pdf"), height=5, width=15)
c=3
ggarrange(lst_CM_ggplot[[3]], lst_EM_ggplot[[3]], lst_EMRA_ggplot[[3]], ncol=3, common.legend = T)
ggsave(filename = here(output_dir, "clus3_composite_vol.pdf"), height=5, width=15)
c=4
ggarrange(lst_CM_ggplot[[4]], lst_EM_ggplot[[4]], lst_EMRA_ggplot[[4]], ncol=3, common.legend = T)
ggsave(filename = here(output_dir, "clus4_composite_vol.pdf"), height=5, width=15)
c=5
ggarrange(lst_CM_ggplot[[5]], lst_EM_ggplot[[5]], lst_EMRA_ggplot[[5]], ncol=3, common.legend = T)
ggsave(filename = here(output_dir, "clus5_composite_vol.pdf"), height=5, width=15)

#####################################################################################
#Saving rnk files for GSEA analysis of Mem/Nav comparisons
for (i in 1:5){
  write.table(lst_CM_rnk[[i]], file=here(output_dir, paste0("CD8_CM_v_Nav_clus", i, ".rnk")), row.names=F, col.names=F, sep="\t", quote=F)
  write.table(lst_EM_rnk[[i]], file=here(output_dir, paste0("CD8_EM_v_Nav_clus", i, ".rnk")), row.names=F, col.names=F, sep="\t", quote=F)
  write.table(lst_EMRA_rnk[[i]], file=here(output_dir, paste0("CD8_EMRA_v_Nav_clus", i, ".rnk")), row.names=F, col.names=F, sep="\t", quote=F)
}




