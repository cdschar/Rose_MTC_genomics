#Script for comparing logFC of pagernk scores vs logFC of Expression via scatter plot
#Edited for stim 03Mar2021
#jrose

setwd("~/Documents/Emory/Boss/data/ATAC/PageRank/pgrank_stim")

library(tidyverse)
library(here)
library(ggrepel)

pgrnk_input_dir <- "input"
diff_input_dir <- "input"
output_dir <- "output"

#Read in PageRank data
d <- read.table(here(pgrnk_input_dir,"TFRanks.tsv"),header=T,row.names=1)
#Filter
x1 = apply(d, 1, max)>1e-4
d<-d[x1,]

#Calculate folds for each score over naive
d1<-data.frame(d$CD4_CM/d$CD4_N,
               d$CD4_EM/d$CD4_N, 
               d$CD4_EMRA/d$CD4_N, 
               d$CD8_CM/d$CD8_N, 
               d$CD8_EM/d$CD8_N, 
               d$CD8_EMRA/d$CD8_N,
               d$CD4_CMStim/d$CD4_NStim,
               d$CD4_EMStim/d$CD4_NStim, 
               d$CD4_EMRAStim/d$CD4_NStim, 
               d$CD8_CMStim/d$CD8_NStim, 
               d$CD8_EMStim/d$CD8_NStim, 
               d$CD8_EMRAStim/d$CD8_NStim,
               d$CD4_CMStim/d$CD4_N,
               d$CD4_EMStim/d$CD4_N, 
               d$CD4_EMRAStim/d$CD4_N, 
               d$CD8_CMStim/d$CD8_N, 
               d$CD8_EMStim/d$CD8_N, 
               d$CD8_EMRAStim/d$CD8_N)

#Log2 transformation of folds
d1 <- apply(d1, c(1,2), function(x) log(x,base=2)) %>% data.frame()
d$gene <- rownames(d)
d1$gene <- rownames(d)


#Read in expression data
RNA_diff <- read.delim(here(diff_input_dir,"diff.significant.glm.Human_Tcell_RNA.StimUnstimRev.txt"), header=T, sep="\t")

# rpkm <- cbind(RNA_diff[,c("SYMBOL")], RNA_diff[,grep("(CD8)|(CD4)_((CM)|(EM)|(EMRA))Stim.rpkm", colnames(RNA_diff))])
# colnames(rpkm_CD8_EMRA )[1] <- "SYMBOL"
# mean_rpkm_CD8_EMRA  <- pivot_longer(rpkm_CD8_EMRA , 
#                                     cols=ends_with('.rpkm'),
#                                     names_to=c("Sample", "Celltype", "Subtype"),
#                                     names_pattern= "RNA_(.*)_(.*)_(.*).rpkm",
#                                     values_to="rpkm"
# ) %>%
#   group_by(SYMBOL, Subtype) %>%
#   summarise(mean_rpkm=mean(rpkm)) %>%
#   pivot_wider(names_from=Subtype, values_from=mean_rpkm) %>%
#   select(SYMBOL, CMStim, EMStim, EMRAStim)

#Assigning column numbers for groups of interest
SYM <- 2
col1 <- 209  #CD4_CM_unstim sig column    #Col # for base R
col2 <- 134  #CD4_CM_stim sig column

col3 <- 439  #CD4_EM_unstim sig  column
col4 <- 384  #CD4_EM_stim sig column

col5 <- 629 #CD8_CM_unstim sig column
col6 <- 594 #CD8_CM_stim sig column

col7 <- 699 #CD8_EM_unstim sig column
col8 <- 684 #CD8_EM_stim sig column

col9 <- 674 #CD8_EMRA_unstim sig column
col10 <- 649 #CD8_EMRA_stim sig column

sig_cols <- c(SYM, col1, col2, col3, col4, col5, col6, col7, col8, col9, col10)
sig_vars <- names(RNA_diff)[sig_cols]  #Variable name for dplyr

fc_cols <- c(2, sig_cols[-1] - 3)
fc_vars <- names(RNA_diff)[fc_cols]

fdr_cols <- c(2, sig_cols[-1] - 1)
fdr_vars <- names(RNA_diff)[fdr_cols]

#rpkm gene expression columns
col_CD4CM <- grep("CD4_CM.rpkm", names(RNA_diff))
names(col_CD4CM) <- names(RNA_diff)[col_CD4CM]
RNA_diff$CD4CM_mean_rpkm <- rowMeans(RNA_diff[,col_CD4CM])
col_stim_CD4CM <- grep("CD4_CMStim.rpkm", names(RNA_diff))
names(col_stim_CD4CM) <- names(RNA_diff)[col_stim_CD4CM]
RNA_diff$CD4CM_meanstim_rpkm <- rowMeans(RNA_diff[,col_stim_CD4CM])

col_CD4EM <- grep("CD4_EM.rpkm", names(RNA_diff))
names(col_CD4EM) <- names(RNA_diff)[col_CD4EM]
RNA_diff$CD4EM_mean_rpkm <- rowMeans(RNA_diff[,col_CD4EM])
col_stim_CD4EM <- grep("CD4_EMStim.rpkm", names(RNA_diff))
names(col_stim_CD4EM) <- names(RNA_diff)[col_stim_CD4EM]
RNA_diff$CD4EM_meanstim_rpkm <- rowMeans(RNA_diff[,col_stim_CD4EM])

col_CD8CM <- grep("CD8_CM.rpkm", names(RNA_diff))
names(col_CD8CM) <- names(RNA_diff)[col_CD8CM]
RNA_diff$CD8CM_mean_rpkm <- rowMeans(RNA_diff[,col_CD8CM])
col_stim_CD8CM <- grep("CD8_CMStim.rpkm", names(RNA_diff))
names(col_stim_CD8CM) <- names(RNA_diff)[col_stim_CD8CM]
RNA_diff$CD8CM_meanstim_rpkm <- rowMeans(RNA_diff[,col_stim_CD8CM])

col_CD8EM <- grep("CD8_EM.rpkm", names(RNA_diff))
names(col_CD8EM) <- names(RNA_diff)[col_CD8EM]
RNA_diff$CD8EM_mean_rpkm <- rowMeans(RNA_diff[,col_CD8EM])
col_stim_CD8EM <- grep("CD8_EMStim.rpkm", names(RNA_diff))
names(col_stim_CD8EM) <- names(RNA_diff)[col_stim_CD8EM]
RNA_diff$CD8EM_meanstim_rpkm <- rowMeans(RNA_diff[,col_stim_CD8EM])

col_CD8EMRA <- grep("CD8_EMRA.rpkm", names(RNA_diff))
names(col_CD8EMRA) <- names(RNA_diff)[col_CD8EMRA]
RNA_diff$CD8EMRA_mean_rpkm <- rowMeans(RNA_diff[,col_CD8EMRA])
col_stim_CD8EMRA <- grep("CD8_EMRAStim.rpkm", names(RNA_diff))
names(col_stim_CD8EMRA) <- names(RNA_diff)[col_stim_CD8EMRA]
RNA_diff$CD8EMRA_meanstim_rpkm <- rowMeans(RNA_diff[,col_stim_CD8EMRA])

#r_fc_cols <- c(r_fc_cols, 246:255)

#List of genes to plot in heatmap
genes_df <- read.csv(here(diff_input_dir,"topgenes.csv"), header=T, colClasses = 'character', na.strings="")

genes_df <- select(genes_df,ATAC_cat_Tfs) %>% pivot_longer(cols=everything(),names_to="path", values_to = "gene") %>% drop_na(gene)
genes_df <- genes_df[!duplicated(genes_df$gene),]
TFs <- genes_df$gene %>% unlist() %>% factor(ordered = T)

#Subsetting and putting bothdatasets together
RNA_TF <- RNA_diff[RNA_diff$SYMBOL %in% d1$gene,] %>% select(SYMBOL,ends_with(".logFC"), ends_with("mean_rpkm"), ends_with("meanstim_rpkm"))
RNA_TF_sub <- RNA_TF[RNA_TF$SYMBOL %in% TFs,]
pgRnk_RNA <- inner_join(d1,RNA_TF_sub, by = c("gene" = "SYMBOL"))

#Plotting

#Thresholds for RNA (exludes below thresold) and combined PgRnk + RNAlogFC values (labels above)
RNA <- 0
combined <- 1.5

theme_set(theme_gray()+ theme(axis.line = element_line(size=0.5),
                              panel.background = element_rect(fill=NA,size=rel(20)),
                              panel.grid.minor = element_line(colour = NA),
                              panel.grid.major = element_line(color=NA),
                              axis.text = element_text(size=10), axis.title = element_text(size=12)))


pCD8_EMRA <- subset(pgRnk_RNA, CD8EMRA_mean_rpkm >=RNA) %>%
  ggplot(aes(x=CD8_EMRA_blood_unstim.v.CD8_Nav_blood_unstim.logFC, y=d.CD8_EMRA.d.CD8_N))
pCD8_EMRA + geom_point() +
  xlab("Expression log2FC (Mem Unstim/Naive Unstim)") +
  ylab("PageRank log2FC (Mem Unstim/Naive Unstim)") +
  geom_hline(yintercept = 0, linetype="dashed") +
  geom_vline(xintercept = 0, linetype="dashed") +
  geom_text_repel(aes(label=ifelse(abs(CD8_EMRA_blood_unstim.v.CD8_Nav_blood_unstim.logFC +d.CD8_EMRA.d.CD8_N >=combined), as.character(gene),"")), hjust=1, vjust=0, size=3)
ggsave(filename=here(output_dir,"CD8_EMRA_pgnrkscatter.pdf"), height=7, width=7)

pCD8_EM <- subset(pgRnk_RNA, CD8EM_mean_rpkm >=RNA) %>%
  ggplot(aes(x=CD8_EM_blood_unstim.v.CD8_Nav_blood_unstim.logFC, y=d.CD8_EM.d.CD8_N))
pCD8_EM + geom_point() +
  xlab("Expression log2FC (Mem Unstim/Naive Unstim)") +
  ylab("PageRank log2FC (Mem Unstim/Naive Unstim)") +
  geom_hline(yintercept = 0, linetype="dashed") +
  geom_vline(xintercept = 0, linetype="dashed") +
  geom_text_repel(aes(label=ifelse(abs(CD8_EM_blood_unstim.v.CD8_Nav_blood_unstim.logFC +d.CD8_EM.d.CD8_N >=combined), as.character(gene),"")), hjust=1, vjust=0, size=3)
ggsave(filename=here(output_dir,"CD8_EM_pgnrkscatter.pdf"), height=7, width=7)

pCD8_CM <- subset(pgRnk_RNA, CD8CM_mean_rpkm >=RNA) %>%
  ggplot(aes(x=CD8_CM_blood_unstim.v.CD8_Nav_blood_unstim.logFC, y=d.CD8_CM.d.CD8_N))
pCD8_CM + geom_point() +
  xlab("Expression log2FC (Mem Unstim/Naive Unstim)") +
  ylab("PageRank log2FC (Mem Unstim/Naive Unstim)") +
  geom_hline(yintercept = 0, linetype="dashed") +
  geom_vline(xintercept = 0, linetype="dashed") +
  geom_text_repel(aes(label=ifelse(abs(CD8_CM_blood_unstim.v.CD8_Nav_blood_unstim.logFC +d.CD8_CM.d.CD8_N >=combined), as.character(gene),"")), hjust=1, vjust=0, size=3)
ggsave(filename=here(output_dir,"CD8_CM_pgnrkscatter.pdf"), height=7, width=7)

pCD4_EM <- subset(pgRnk_RNA, CD4EM_mean_rpkm >=RNA) %>%
  ggplot(aes(x=CD4_EM_blood_unstim.v.CD4_Nav_blood_unstim.logFC, y=d.CD4_EM.d.CD4_N))
pCD4_EM + geom_point() +
  xlab("Expression log2FC (Mem Unstim/Naive Unstim)") +
  ylab("PageRank log2FC (Mem Unstim/Naive Unstim)") +
  geom_hline(yintercept = 0, linetype="dashed") +
  geom_vline(xintercept = 0, linetype="dashed") +
  geom_text_repel(aes(label=ifelse(abs(CD4_EM_blood_unstim.v.CD4_Nav_blood_unstim.logFC +d.CD4_EM.d.CD4_N >=combined), as.character(gene),"")), hjust=1, vjust=0, size=3)
ggsave(filename=here(output_dir,"CD4_EM_pgnrkscatter.pdf"), height=7, width=7)

pCD4_CM <- subset(pgRnk_RNA, CD4CM_mean_rpkm >=RNA) %>% 
  ggplot(aes(x=CD4_CM_blood_unstim.v.CD4_Nav_blood_unstim.logFC, y=d.CD4_CM.d.CD4_N))
pCD4_CM + geom_point() +
  xlab("Expression log2FC (Mem Unstim/Naive Unstim)") +
  ylab("PageRank log2FC (Mem Unstim/Naive Unstim)") +
  geom_hline(yintercept = 0, linetype="dashed") +
  geom_vline(xintercept = 0, linetype="dashed") +
  geom_text_repel(aes(label=ifelse(abs(CD4_CM_blood_unstim.v.CD4_Nav_blood_unstim.logFC +d.CD4_CM.d.CD4_N >=combined), as.character(gene),"")), hjust=1, vjust=0, size=3)
ggsave(filename=here(output_dir,"CD4_CM_pgnrkscatter.pdf"), height=7, width=7)

pCD8_EMStim <- subset(pgRnk_RNA, CD8EM_meanstim_rpkm >=RNA) %>%
  ggplot(aes(x=CD8_EM_blood_stim.v.CD8_Nav_blood_stim.logFC, y=d.CD8_EMStim.d.CD8_NStim))
pCD8_EMStim + geom_point() +
  xlab("Expression log2FC (Mem Stim/Naive Stim)") +
  ylab("PageRank log2FC (Mem Stim/Naive Stim)") +
  geom_hline(yintercept = 0, linetype="dashed") +
  geom_vline(xintercept = 0, linetype="dashed") +
  scale_x_continuous(limits = c(-4,4)) +
  scale_y_continuous(limits = c(-4,4)) +
  #geom_text_repel(aes(label=gene), hjust=1, vjust=0, size=3)
  geom_text_repel(aes(label=ifelse(abs(CD8_EM_blood_stim.v.CD8_Nav_blood_stim.logFC +d.CD8_EMStim.d.CD8_NStim) >=combined | gene=="TBX21", as.character(gene),"")), hjust=1, vjust=0, size=3, min.segment.length = 0)
ggsave(filename=here(output_dir,"CD8_EMstim_pgnrkscatter.pdf"), height=4, width=4)

pCD8_CMStim <- subset(pgRnk_RNA, CD8CM_meanstim_rpkm >=RNA) %>%
  ggplot(aes(x=CD8_CM_blood_stim.v.CD8_Nav_blood_stim.logFC, y=d.CD8_CMStim.d.CD8_NStim))
pCD8_CMStim + geom_point() +
  xlab("Expression log2FC (Mem Stim/Naive Stim)") +
  ylab("PageRank log2FC (Mem Stim/Naive Stim)") +
  geom_hline(yintercept = 0, linetype="dashed") +
  geom_vline(xintercept = 0, linetype="dashed") +
  geom_text_repel(aes(label=ifelse(abs(CD8_CM_blood_stim.v.CD8_Nav_blood_stim.logFC +d.CD8_CMStim.d.CD8_NStim) >=combined, as.character(gene),"")), hjust=1, vjust=0, size=3)
ggsave(filename=here(output_dir,"CD8_CMstim_pgnrkscatter.pdf"), height=7, width=7)

pCD4_EMStim <- subset(pgRnk_RNA, CD4EM_meanstim_rpkm >=RNA) %>%
  ggplot(aes(x=CD4_EM_blood_stim.v.CD4_Nav_blood_stim.logFC, y=d.CD4_EMStim.d.CD4_NStim))
pCD4_EMStim + geom_point() +
  xlab("Expression log2FC (Mem Stim/Naive Stim)") +
  ylab("PageRank log2FC (Mem Stim/Naive Stim)") +
  geom_hline(yintercept = 0, linetype="dashed") +
  geom_vline(xintercept = 0, linetype="dashed") +
  scale_x_continuous(limits = c(-4,4)) +
  scale_y_continuous(limits = c(-4,4)) +
  geom_text_repel(aes(label=ifelse(abs(CD4_EM_blood_stim.v.CD4_Nav_blood_stim.logFC +d.CD4_EMStim.d.CD4_NStim) >=combined, as.character(gene),"")), hjust=1, vjust=0, size=3)
ggsave(filename=here(output_dir,"CD4_EMstim_pgnrkscatter.pdf"), height=7, width=7)

pCD4_CMStim <- subset(pgRnk_RNA, CD4CM_meanstim_rpkm >=RNA) %>%
  ggplot(aes(x=CD4_CM_blood_stim.v.CD4_Nav_blood_stim.logFC, y=d.CD4_CMStim.d.CD4_NStim))
pCD4_CMStim + geom_point() +
  xlab("Expression log2FC (Mem Stim/Naive Stim)") +
  ylab("PageRank log2FC (Mem Stim/Naive Stim)") +
  geom_hline(yintercept = 0, linetype="dashed") +
  geom_vline(xintercept = 0, linetype="dashed") +
  geom_text_repel(aes(label=ifelse(abs(CD4_CM_blood_stim.v.CD4_Nav_blood_stim.logFC +d.CD4_CMStim.d.CD4_NStim) >=combined, as.character(gene),"")), hjust=1, vjust=0, size=3)
ggsave(filename=here(output_dir,"CD4_CMstim_pgnrkscatter.pdf"), height=7, width=7)