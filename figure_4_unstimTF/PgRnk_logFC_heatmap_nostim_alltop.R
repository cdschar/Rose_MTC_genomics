#Script for quick heatmap of PageRank scores of known list of genes
setwd("~/Documents/Emory/Boss/data/ATAC")

library(tidyverse)
library(here)
library(pheatmap)
library(ComplexHeatmap)
library(circlize)
library(viridis)
library(RColorBrewer)

input_dir <- "PageRank/Run1"
output_dir <- "PageRank/Run1/heatmap"

load(file=here(input_dir, "scatter_logFCNav", "topOrdered_TmemSubsets_PgRnk_list.rda"))
bottom_TFs_df <- map_dfr(bottom, "gene")
TFs_df <- map_dfr(top, "gene")
#TFs_df <- TFs_df[!grepl("Stim", names(TFs_df))]
#write.csv(TFs_df, here(input_dir, "Top_TFs.csv"))

expr_genes <- read.csv("/home/jim/Documents/Emory/Boss/data/RNA/expressedgenes.csv", sep=",", header=F)
expr_genes <- drop_na(expr_genes, V1) %>% select(V2) %>% unlist() %>% as.character()

d <- read.table(here(input_dir,"TFRanks.tsv"),header=T,row.names=1)
#Filter out low values
x1 = apply(d, 1, max)>1e-4
d<-d[x1,]

#Calculate folds for each score over naive
d1<-data.frame(d$CD4_CM/d$CD4_N,
               d$CD4_EM/d$CD4_N, 
               d$CD4_EMRA/d$CD4_N, 
               d$CD8_CM/d$CD8_N, 
               d$CD8_EM/d$CD8_N, 
               d$CD8_EMRA/d$CD8_N,
               # d$CD4_CMStim/d$CD4_NStim,
               # d$CD4_EMStim/d$CD4_NStim, 
               # d$CD4_EMRAStim/d$CD4_NStim, 
               # d$CD8_CMStim/d$CD8_NStim, 
               # d$CD8_EMStim/d$CD8_NStim, 
               # d$CD8_EMRAStim/d$CD8_NStim,
               d$CD4_CMStim/d$CD4_N,
               d$CD4_EMStim/d$CD4_N, 
               d$CD4_EMRAStim/d$CD4_N, 
               d$CD8_CMStim/d$CD8_N, 
               d$CD8_EMStim/d$CD8_N, 
               d$CD8_EMRAStim/d$CD8_N,
               d$CD4_N/d$CD4_N,
               d$CD8_N/d$CD8_N,
               d$CD4_NStim/d$CD4_N,
               d$CD8_NStim/d$CD8_N)

#Log2 transformation of folds
d1 <- apply(d1, c(1,2), function(x) log(x,base=2)) %>% data.frame()
rownames(d1) <- rownames(d)
names(d1) <- gsub(".d.CD[48]_N", "", names(d1))
names(d1) <- gsub("d.", "", names(d1))

# TFs <- c("HIF1A",
#           "MYBL2",
#           "XBP1",
#           "CEBPD",
#           "ASCL2",
#           "MSC",
#           "ARID5B",
#           "GSX2",
#           "SMAD3",
#           "IRF6",
#           "TCF7",
#           "LEF1",
#           "AHR", 
#           "EOMES",
#          "PRDM1",
#          "TBX21")

#TFs <- as.character(unique(unlist(TFs_df[1:10,])))
TFs_dfalltop <- TFs_df[,14:18]

exclude <- function(x){
        index <- x%in%expr_genes
        return(x[index])
}
TFs_checked <- apply(TFs_dfalltop, c(2), exclude)

#TFs <- as.character(unique(unlist(TFs_dfalltop[1:10,])))
#TFs <- unique(unlist(lapply(TFs_checked, "[", c(1:10))))
#TFs <- TFs[!grepl("MSC", TFs)]
TFs <- c("CREB3","SMAD1","ARID3B","BCL6", "XBP1", "MSC", "AHR", "MEOX1", "TBX21", "FOSL2", "EOMES", "CEBPD", "RORC", "FOXP3", "BHLHE40", "PRDM1", "LEF1", "TCF7", "ASCL2", "SMAD3", "TERF2", "ARID5A", "EGR4", "HIF1A", "PPARD", "STAT5B", "GMEB1")

index <- rownames(d1)%in%TFs
d1<- d1[index,]
d1 <- d1[,!grepl("CD4_EMRA", colnames(d1))]

d1_CD4 <- d1[,grep("CD4", colnames(d1))]
d1_CD8 <- d1[,grep("CD8", colnames(d1))]
d1_nostim <- d1[,!grepl("Stim", colnames(d1))]
d1_CD8_nostim <- d1_CD8[,!grepl("Stim", colnames(d1_CD8))]

# d_log <- apply(d, c(1,2), function(x) log(x,base=10))
# d_log_nostim <- d_log[,!grepl("Stim", colnames(d_log))]
# d_log_nostim_CD8 <- d_log_nostim[,grep("CD8", colnames(d_log_nostim))]
d_nostim_noNav <- d1_nostim[,!grepl("_N", colnames(d1_nostim))]

# pdf(file=here(output_dir, "top10_PgRnkHeat_raw.pdf"), useDingbats = F, height=7, width=7)
# pheatmap(d, fontsize_row = 7)
# dev.off()
# 
# pdf(file=here(output_dir, "top10_PgRnkHeat_rawzscore.pdf"), useDingbats=F, height=7, width=7)
# pheatmap(d, scale="row", fontsize_row = 7)
# dev.off()
# 
# pdf(file=here(output_dir, "top10_CD4_PgRnkHeat_rawzscore.pdf"), useDingbats=F, height=7, width=7)
# pheatmap(d_CD4, scale="row", fontsize_row = 7)
# dev.off()
# 
# pdf(file=here(output_dir, "top10_CD8_PgRnkHeat_rawzscore.pdf"), useDingbats=F, height=7, width=7)
# pheatmap(d_CD8, scale="row", fontsize_row = 7)
# dev.off()

# pdf(file=here(output_dir, "top10_unstim_PgRnkHeat_logFC_nav_noMSC.pdf"), useDingbats=F, height=7, width=7)
# pheatmap(d1_nostim,
#          #scale= "row",
#          fontsize_row=7
#          )
# dev.off()
# 
# pdf(file=here(output_dir, "top10_unstim_PgRnkHeat_logFC_nav_zscore_noMSC.pdf"), useDingbats=F, height=7, width=7)
# pheatmap(d1_nostim,
#          scale= "row",
#          fontsize_row=7
# )
# dev.off()
# 
# pdf(file=here(output_dir, "top10_CD8unstim_PgRnkHeat_logFC_nav_noMSC.pdf"), useDingbats=F, height=7, width=7)
# pheatmap(d1_CD8_nostim,
#          #scale= "row",
#          fontsize_row=7
# )
# dev.off()
# 
# pdf(file=here(output_dir, "top10_CD8unstim_PgRnkHeat_logFC_nav_zscore_noMSC.pdf"), useDingbats=F, height=7, width=7)
# pheatmap(d1_CD8_nostim,
#          scale= "row",
#          fontsize_row=7
# )
# dev.off()
# 
# pdf(file=here(output_dir, "top10_CD8unstim_PgRnkHeat_logFC_nav_noNav.pdf"), useDingbats=F, height=4, width=2.5)
# pheatmap(d_nostim_noNav,
#          scale= "row",
#          fontsize_row=7,
#          treeheight_col = 25,
#          treeheight_row = 25
#         )
# dev.off()
# pheatmap(d_CD4, scale="row")
# pheatmap(d_CD8, scale="row")
# pheatmap(d_CD8_nostim, scale="row", fontsize_row = 7)
# pheatmap(d_nostim)
# pheatmap(d_nostim, scale= "row", fontsize_row=7)
# pheatmap(d_log)
# pheatmap(d_log_nostim)
# pheatmap(d_log_nostim_CD8)

#Colors
custom_pal <- c("#90A4ADFF","#5F7D8BFF", "#2096F2FF", "#0C46A0FF", "#FFB74CFF", "#E55100FF", "#870D4EFF")
names(custom_pal) <- c("CD4_Nav", "CD8_Nav", "CD4_CM", "CD8_CM", "CD4_EM", "CD8_EM", "CD8_EMRA")

#ColorRamp for viridis
col = colorRamp2(c(-8, 0, 8), c("navy", "grey98", "darkgoldenrod2"), space = "sRGB")
col2 = colorRamp2(c(8,0,-8), brewer.pal(3, "PiYG"), space="RGB")

pdf(file=here(output_dir, "top10_CD8unstim_PgRnkcurrated_logFC_nav_noNav2.pdf"), useDingbats=F, height=5, width=3.7)
Heatmap(d_nostim_noNav, 
        row_names_gp = gpar(fontsize=7),
        cluster_columns = F,
        column_dend_height = unit(7, "mm"),
        col = col,
        top_annotation = HeatmapAnnotation(cellSubset=colnames(d_nostim_noNav),
                                           col=list(cellSubset=custom_pal)
                                            ),
        show_column_names = F
        )
dev.off()
