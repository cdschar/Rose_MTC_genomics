#Heatmap for pre-specified list of genes
#v3 uses compelx heatmap
#20Mar20
#jrose

library(here)
library(tidyverse)
library(ComplexHeatmap)
# library(Rmpfr)
library(circlize)
library(viridis)

input_dir <- "heatmap"
output_dir <- "heatmap/all"

#Plotting variables
filename = "HIF1A_enz_hyp_ad"
TF <- "HIF1A"
row_font = 12
show_row = T
show_col = F
clus_row = T
clus_col = T
custom_pal <- c("CD4_Nav"="#90A4ADFF","CD8_Nav"="#455964FF", "CD4_CM"="#2096F2FF", "CD8_CM"="#0C46A0FF", "CD4_EM"="#FFB74CFF", "CD8_EM"="#E55100FF", "CD8_EMRA"="#870D4EFF")

#List of genes to plot in heatmap
genes_df <- read.csv(here("topgenes.csv"), header=T, colClasses = 'character', na.strings="")

genes_df <- select(genes_df, HIF1A_enzyme,HIF1A_hypoxia,HIF1A_adhesion) %>% pivot_longer(cols=everything(),names_to="path", values_to = "gene") %>% drop_na(gene)
genes_df <- genes_df[!duplicated(genes_df$gene),]
genes <- genes_df$gene %>% unlist() %>% factor(ordered = T)

#genes <- c("MYC","TCF7","LEF1", "INFG","TNF","CCR7", "IL6R", "EOMES", "PRDM1", "PRF1", "MYBL2", "TOX", "TBX21", "SMAD3", "IL2RB", "MYBL1", "KLRG1", "CCR4", "S1PR5", "ARID5B", "ARID5A", "ITGB1", "FAS", "GZMK", "ANXA2", "PYHIN1")
#genes <- c("IL2RA", "IL2RB")

#Load in RNA data sample manifest
fqDir = "/home/jim/Documents/Emory/Boss/data/"
fqFile = "RNAseq.sample.manifest.txt";
files = read.table(paste0(fqDir, fqFile), sep = "\t", header = T, as.is = T);
files = files[grepl("EM|CM|Nav", files$group),]
files = files[files$stim=="unstim",]
files = files[!grepl("CD4_EMRA", files$group),]

#read in data file
dataDir = "/home/jim/Documents/Emory/Boss/data/RNA/naive_v_mem_overlap/";
dataFile = "diff.significant.glm.Human_Tcell_RNA.rev.txt";
dataFile2 = "diff.glm.Human_Tcell_RNA.rev.txt";
data = read.delim(file=paste0(dataDir, dataFile), header=T, sep="\t");

#use genes from geneset only
data.heat = data[data$SYMBOL %in% genes,]
rownames <- factor(data.heat$SYMBOL, 
                   levels = genes, 
                   #ordered = T
                  )

#Getting average expression data by subtype for MA plot
col_CD8Nav <- grep("CD8_N.rpkm", names(data.heat))
data.heat$CD8_Nav_mean_rpkm <- rowMeans(data.heat[,col_CD8Nav])

col_CD8CM <- grep("CD8_CM.rpkm", names(data.heat))
data.heat$CD8_CM_mean_rpkm <- rowMeans(data.heat[,col_CD8CM])

col_CD8EM <- grep("CD8_EM.rpkm", names(data.heat))
data.heat$CD8_EM_mean_rpkm <- rowMeans(data.heat[,col_CD8EM])

col_CD8EMRA <- grep("CD8_EMRA.rpkm", names(data.heat))
data.heat$CD8_EMRA_mean_rpkm <- rowMeans(data.heat[,col_CD8EMRA])

col_CD4Nav <- grep("CD4_N.rpkm", names(data.heat))
data.heat$CD4_Nav_mean_rpkm <- rowMeans(data.heat[,col_CD4Nav])

col_CD4CM <- grep("CD4_CM.rpkm", names(data.heat))
data.heat$CD4_CM_mean_rpkm <- rowMeans(data.heat[,col_CD4CM])

col_CD4EM <- grep("CD4_EM.rpkm", names(data.heat))
data.heat$CD4_EM_mean_rpkm <- rowMeans(data.heat[,col_CD4EM])

#Matrix of averages
data.sub <- select(data.heat, CD8_Nav_mean_rpkm, CD8_CM_mean_rpkm, CD8_EM_mean_rpkm, CD8_EMRA_mean_rpkm, CD4_Nav_mean_rpkm, CD4_CM_mean_rpkm, CD4_EM_mean_rpkm)
data.ave = t(apply(data.sub, 1, scale))
rownames(data.ave) <- rownames
colnames(data.ave) <- gsub("_mean_rpkm", "", colnames(data.sub))
data.ave = data.ave[order(match(rownames,genes)),]


#select data columns & removing suffix
data.heat = data.heat[, grep(".rpkm", names(data.heat))]
data.heat = data.heat[ ,grep(paste(files$sample, collapse = "|"), names(data.heat))]
data.heat = data.heat[, !grepl("Stim", names(data.heat))]
data.heat = data.heat[, !grepl("CD4_EMRA", names(data.heat))]
colnames <- colnames(data.heat)
data.heat = t(apply(data.heat, 1, scale))
rownames(data.heat) <- rownames
colnames(data.heat) <- gsub(".rpkm", "", colnames)
data.heat = data.heat[order(match(rownames,genes)),]


#Annotation File
ph_anno <- select(files, cellType, cellSubset) %>% mutate(MemSub=paste(cellType, cellSubset, sep="_")) %>% select(MemSub)
rownames(ph_anno) <- files$sample
ph_anno <- ph_anno[match(colnames(data.heat), rownames(ph_anno)),]

row_anno <- data.frame(path=genes_df$path[match(rownames(data.heat), genes_df$gene)])


#ColorRamp for viridis
col = colorRamp2(c(-3.5, 0, 3.5), c("#440154FF", "#29AF7FFF", "#FDE725FF"), space = "sRGB")

data2 = read.delim(file=paste0(dataDir, dataFile2), header=T, sep="\t");
rpkm <- data2[,grepl(paste(files$sample, collapse = "|"),colnames(data)) & !grepl("Stim", colnames(data))]
rpkm <- rpkm[,!grepl("CD4_EMRA", colnames(rpkm))]
#rownames(rpkm) <- data$SYMBOL
rpkm$SYMBOL <- data2$SYMBOL
mean_rpkm <- pivot_longer(rpkm,
                          cols=ends_with('.rpkm'),
                          names_to=c("Sample", "Celltype", "Subtype"),
                          names_pattern= "RNA_(.*)_(.*)_(.*).rpkm",
                          values_to="rpkm"
) %>%
  mutate(CellSubtype=paste(Celltype, Subtype, sep="_")) %>%
  group_by(SYMBOL, CellSubtype) %>%
  summarise(mean_rpkm=mean(rpkm)) %>%
  pivot_wider(names_from=CellSubtype, values_from=mean_rpkm) %>%
  ungroup()

colnames(mean_rpkm) <- gsub("_N", "_Nav", colnames(mean_rpkm))

#Colors
CD4.nav.col = "#90A4ADFF"
CD4.cm.col =  "#2096F2FF"
CD4.em.col = "#FFB74CFF"
CD8.nav.col = "#5F7D8BFF"
CD8.cm.col =  "#0C46A0FF"
CD8.em.col = "#E55100FF"
CD8.emra.col ="#870D4EFF"
cellcolors = c("CD4_Nav"=CD4.nav.col,
               "CD4_CM"=CD4.cm.col,
               "CD4_EM"=CD4.em.col,
               "CD8_Nav"=CD8.nav.col,
               "CD8_CM"=CD8.cm.col,
               "CD8_EM"=CD8.em.col,
               "CD8_EMRA"=CD8.emra.col)
#Stolen from here: https://gist.github.com/hrbrmstr/42c21f512fa4c4a81aa4
# # CreateAdjacencyMatrix <- function(x) {
#   
#   s <- gsub("\\.", "", x)
#   m <- matrix(0, 10, 10)
#   
#   for (i in 1:(nchar(s)-1)) {
#     m[as.numeric(substr(s, i, i))+1,
#       as.numeric(substr(s, i+1, i+1))+1] <-
#       m[as.numeric(substr(s, i, i))+1,
#         as.numeric(substr(s, i+1, i+1))+1]+1
#     
#   }
#   
#   rownames(m) = 0:9
#   colnames(m) = 0:9
#   
#   m
#   
# }
# 
# m1 = CreateAdjacencyMatrix(formatMpfr(Const("pi",2000)))
# col = colorRamp2(quantile(m1, seq(0, 1, by = 0.25)), viridis(5))

#Plot heatmaps
pdf(file=here(output_dir, paste(filename,"_all_Complexheatmap.pdf", sep="")), useDingbats = F)
  Heatmap(data.heat,
            show_row_names = show_row,
            show_column_names = show_col,
            row_names_gp = gpar(fontsize = row_font),
            show_heatmap_legend = T,
            cluster_rows = clus_row,
            cluster_columns = clus_col,
            row_split = row_anno,
            border = T,
            #col = col,
            top_annotation = HeatmapAnnotation(cellType=ph_anno,
                                                col=list(cellType=custom_pal)
                                               )
    )
dev.off()

pdf(file=here(output_dir, paste(filename,"_mean_all_Complexheatmap.pdf", sep="")), useDingbats = F)

  Heatmap(data.ave,
          show_row_names = show_row,
          show_column_names = show_col,
          row_names_gp = gpar(fontsize = row_font),
          show_heatmap_legend = T,
          cluster_rows = clus_row,
          cluster_columns = clus_col,
          row_split = row_anno,
          border = T,
          top_annotation = HeatmapAnnotation(TF_expression=anno_barplot(mean_rpkm[grep(paste0("^",TF, "$"),mean_rpkm$SYMBOL),colnames(data.ave)] %>% unlist(),
                                                                        gp=gpar(fill=cellcolors[colnames(data.ave)])),
                                             cellType= c('CD8_Nav', 'CD8_CM', "CD8_EM", 'CD8_EMRA', 'CD4_Nav', 'CD4_CM', 'CD4_EM'),
                                             col=list(cellType=c("CD8_Nav"="#455964FF", "CD8_CM"="#0C46A0FF","CD8_EM"="#E55100FF","CD8_EMRA"="#870D4EFF", "CD4_EM"="#FFB74CFF", "CD4_CM"="#2096F2FF", "CD4_Nav"="#90A4ADFF"))
          )
          #col= col
  )
dev.off()
