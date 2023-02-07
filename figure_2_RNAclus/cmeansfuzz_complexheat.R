#Script for visualizing results (txt file with rpkm and cluster vlaues) of k-means clustering script
#jrose
#Edited 06Aug2020

library(here)
library(ComplexHeatmap)
library(tidyverse)
library(som)
library(ggpubr)
library(cowplot)

#read in the manifest file
fqDir = "/home/jim/Documents/Emory/Boss/data/";
fqFile = "RNAseq.sample.manifest.txt";
files = read.table(paste0(fqDir, fqFile), sep = "\t", header = T, as.is = T);
files = files[grepl("EM|CM|Nav", files$group),]
files = files[!grepl("CD4", files$group),]
#files = files[!grepl("CD4_EMRA", files$group),]
files = files[which(files$include),]
files = files[files$stim=="unstim",]
files = mutate(files, cellSubtype = paste(cellType, cellSubset, stim, sep="_"))

#Reorder based on prespecified order
smpl_order = c("CD8_Nav_unstim", "CD8_CM_unstim", "CD8_EM_unstim", "CD8_EMRA_unstim")
files = files[order(match(files$cellSubtype, smpl_order)),]

#Read in dataset with k-means cluster values
input_dir <- "kmeans"
input_file <- "CD8_unstim_labfilter_RNA_fuzzmeans_clusterdata_v2.txt"
output_dir <- "kmeans/k5/fuzz"
dataDiff <- read.table(file=here(input_dir, input_file), header=T)
data <- drop_na(dataDiff)

#Read in list of DEGs and other genes of interest
genes_df <- read.csv(here("topgenes.csv"), header=T, colClasses = 'character', na.strings="")
#Subset to only DEGs
data <- data[data$SYMBOL%in%genes_df$DEG,]

#Genes to highlight on heatmap
genes_highlight <- c("CEBPD","SATB1", "DNMT3A", "IGF1R","TCF7","RORA","RORC", "LDHB","PFKFB2","BACH1","BACH2", "TCF3","GABPA","NRF1","PDK1","JUNB", "JUND", "AHR", "HIF1A", "TP73", "PDCD1", "OGDH", "PGAM1", "MKI67", "DUSP4",  "GAPDH", "IFNG", "PDP1", "IL2RA","GZMB", "LEF1", "ICOS", "PRF1", "CD28", "CTLA4", "CCR7",  "IL7R", "PRDM1", "KLRG1", "SELL", "CD27", "TBX21", "EOMES", "BATF")

#Select k
k = 5
subset = "CD8_unstim_labfilter07Oct"

#Arragne dataset by selected cluster in accending order
data <- arrange(data, !!(paste("k", k,"cluster", sep=".")))
#Select clusters to use
clusters <- select(data, !!(paste("k", k,"cluster", sep=".")))
#cluster_col_index <- grep(paste("k", k,"cluster", sep="."), names(data))
#Order data by cluster
#data <- data[order(data[,cluster_col_index], decreasing = F),]

#Parameters for plots
show_row = F
show_col = F
row_font = 8
show_row_title = F

#Colors
nav.col = "#5F7D8BFF"
cm.col =  "#0C46A0FF"
em.col = "#E55100FF"
emra.col ="#870D4EFF"
cellcolors = c("CD8_Nav_unstim"=nav.col,
               "CD8_CM_unstim"=cm.col,
               "CD8_EM_unstim"=em.col,
               "CD8_EMRA_unstim"=emra.col)


#Subset data to matrix of RNA counts
rpkm <- data[,grep(".rpkm", names(data))]
#Normalize by row
rpkm <- som::normalize(rpkm)
rownames(rpkm) <- data$SYMBOL
colnames(rpkm) <- colnames(data)[grep(".rpkm", names(data))]

#Greyscale for clusters on plot
grayscale <- function(x, k){
  gray((k-x)/k)
}
clus <- unique(clusters)[order(unique(clusters)),]
clustercolors <- grayscale(clus, k=k)
names(clustercolors) <- clus
#pie(rep(1, length(clustercolors)),col=clustercolors)

#GGplot Averages for dot plot
###############################################################################

ave_data <- cbind(rpkm, data[,c(1,2,grep(".clus", names(data)))])
rpkm_tidy <- pivot_longer(ave_data, cols=ends_with(".rpkm"), names_to = "sample", values_to = "rpkm")
rpkm_tidy$sample <- gsub(".rpkm", "", rpkm_tidy$sample)
rpkm_tidy <- left_join(rpkm_tidy, select(files, sample, group, cellType, tissue, stim, cellSubset, cellSubtype), by="sample")
names(rpkm_tidy)[grep(paste("k", k,"cluster", sep="."), names(rpkm_tidy))] <- "main.cluster"
means <- rpkm_tidy %>% group_by(main.cluster, cellSubtype) %>% summarise(mean=mean(rpkm), sd=sd(rpkm))
means$cellSubtype <- factor(means$cellSubtype, levels=smpl_order, ordered=T)

#Creating matrix for means heatmap
gene_means <- rpkm_tidy %>% group_by(SYMBOL, cellSubtype, main.cluster) %>% 
  summarise(mean=mean(rpkm)) %>% 
  pivot_wider(names_from = cellSubtype, values_from=mean) %>%
  select(SYMBOL, main.cluster, smpl_order) %>%
  arrange(main.cluster)
gene_means <- gene_means[gene_means$SYMBOL%in%genes_df$DEG,]

means_rpkm <- as.matrix(gene_means[,c(-1, -2)])
rownames(means_rpkm) <- gene_means$SYMBOL


#Plotting 
###############################################################################
#line chart
theme_set(theme_light()+ theme(axis.line = element_line(size=0.5),
                              axis.text = element_blank(),
                              axis.title = element_blank(),
                              legend.text = element_text(size=7),
                              legend.title = element_text(size=8)
                              )
         )
p <- ggplot(means, aes(x=cellSubtype, y=mean))
pline <- p + geom_point(aes(color=as.factor(cellSubtype), size=5)) +
 geom_line(aes(group=as.factor(main.cluster))) +
 geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd, group=cellSubtype)) +
 xlab("Cell Subtype") +
 ylab("Average Normalized RPKM") +
 facet_wrap(vars(main.cluster), ncol=1) +
 scale_color_manual(name="Cell Subtype", values=cellcolors) +
 guides(size=FALSE)

pline2 <- p + geom_point(aes(color=as.factor(cellSubtype), size=5)) +
  geom_line(aes(group=as.factor(main.cluster))) +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd, group=cellSubtype)) +
  xlab("Cell Subtype") +
  ylab("Average Normalized RPKM") +
  facet_wrap(vars(main.cluster), ncol=1) +
  scale_color_manual(name="Cell Subtype", values=cellcolors) +
  guides(size=FALSE, color=FALSE)

pline
ggsave(filename=here(output_dir,paste(subset, "_linePlot_gg_fuzz_k", k, ".pdf", sep="")), width=3, height=6.8)
dev.off()


lbl_index <- vector()
for (i in 1:length(genes_highlight)){
  if (genes_highlight[i]%in%rownames(rpkm)) {
    lbl_index[i] <- grep(genes_highlight[i], rownames(rpkm))
  }
}

#Heatmaps
pdf(file=here(output_dir, paste(subset, "_heatmap_gg_fuzz_k", k, ".pdf", sep="")), useDingbats = F, height=7)
  Heatmap(rpkm,
          show_row_names = show_row,
          show_column_names = show_col,
          row_names_gp = gpar(fontsize = row_font),
          show_heatmap_legend = T,
          cluster_rows = F,
          cluster_columns = F,
          row_split = clusters,
          use_raster = T,
          raster_device = "png",
          #column_split = cols[1:colnum],
          #row_title_rot = 0,
          #row_title = show_row_title
          #border = T,
          #col = col,
          top_annotation = HeatmapAnnotation(cellSubset=as.vector(files$cellSubtype),
                                             col=list(cellSubset=cellcolors)
                                             ),
          left_annotation = rowAnnotation(cluster=clusters[,1],
                                          col=list(cluster=clustercolors)
                                          ),
          right_annotation = rowAnnotation(gene = anno_mark(at = lbl_index, labels = genes_highlight))
  )
dev.off()

#Heatmap of means
pdf(file=here(output_dir, paste(subset, "_heatmap_means_gg_fuzz_k", k, ".pdf", sep="")), useDingbats = F, height=7)
  Heatmap(means_rpkm,
          show_row_names = show_row,
          show_column_names = show_col,
          row_names_gp = gpar(fontsize = row_font),
          show_heatmap_legend = T,
          cluster_rows = F,
          cluster_columns = F,
          row_split = gene_means$main.cluster,
          #column_split = cols[1:colnum],
          #row_title_rot = 0,
          #row_title = show_row_title
          border = T,
          #col = col,
          top_annotation = HeatmapAnnotation(cellSubset=smpl_order,
                                             col=list(cellSubset=cellcolors)
          ),
          left_annotation = rowAnnotation(cluster=gene_means$main.cluster,
                                          col=list(cluster=clustercolors)
          ),
          #right_annotation = rowAnnotation(foo=anno)
  )
dev.off()
  
#Plotting individuals cluster heatmaps & membership plot
#####################################################################################

#Creating tidy membership dataframe
data_mem <- cbind(rpkm, data[,grep(paste("k", k,"clus.", sep="."), names(data))])
names(data_mem) <- gsub(paste("k.", k, ".", sep=""), "", names(data_mem))
data_mem$gene <- factor(rownames(data_mem))

colnum = dim(files)[1]

theme_set(theme_light()+ theme(axis.line = element_line(size=0.5),
                               axis.text = element_text(size=8),
                               axis.title = element_blank(),
                               legend.text = element_text(size=7),
                               legend.title = element_text(size=8)
                                )
        )


for (h in 1:k){
  
  ###
  #Membership plot
  data_mem_tmp <- subset(data_mem, cluster==h) %>% pivot_longer(1:colnum, names_to="sample")
  names(data_mem_tmp)[grep(paste("clus", h,sep="."), names(data_mem_tmp))] <- "membership"
  
  #Add in cellSubtype from files dataframe
  data_mem_tmp$sample <- gsub(".rpkm", "", data_mem_tmp$sample)
  data_mem_tmp <- left_join(data_mem_tmp,  select(files, sample, cellSubtype), by="sample")
  data_mem_tmp$cellSubtype <- factor(data_mem_tmp$cellSubtype, levels=smpl_order, ordered=T)
  
  data_mem_tmp <- data_mem_tmp[order(data_mem_tmp$cellSubtype),]
  
  data_mem_tmp_mean <- data_mem_tmp %>% group_by(gene, cellSubtype) %>% summarise(mean=mean(value), membership = mean(membership),gene=gene, cluster=cluster)
  #order the dataframe by score
  #data_mem_tmp <- data_mem_tmp[order(data_mem_tmp$membership),]
  #set the order by setting the factors using forcats
  levels(data_mem_tmp$gene) = fct_inorder(data_mem$gene)
  
  #Plot 
  ggplot(data_mem_tmp_mean, aes(x=cellSubtype, y=mean)) +
    geom_line(aes(color=membership, group=gene, alpha=membership)) +
    scale_color_gradientn(colors=c('blue1','gold')) +
    stat_summary(aes(group=1), fun = "mean", geom="line")
  ggsave(filename=here(output_dir,paste(subset, "_member_gg_fuzz_clus", h, ".pdf", sep="")))
  
  ###
  
  #Heatmap
  #Subsets out genes with membership lower than 0.1
  sub_rpkmtidy <- subset(rpkm_tidy, main.cluster==h & !!paste("k",k,"clus",h,"mem", sep=".")>=0.7)
  sub_genes <- unique(sub_rpkmtidy$SYMBOL)
  sub_rpkm <- rpkm[rownames(rpkm)%in%sub_genes,]
  
  print({
  pdf(file=here(output_dir, paste(subset, "_heatmap_clusteronly_fuzz_k", h, ".pdf", sep="")), useDingbats = F, height=4)
    Heatmap(sub_rpkm,
            show_row_names = T,
            show_column_names = T,
            row_names_gp = gpar(fontsize = 3),
            column_names_gp = gpar(fontsize=4),
            show_heatmap_legend = T,
            cluster_rows = T,
            cluster_columns = F, 
            top_annotation = HeatmapAnnotation(cellSubset=as.vector(files$cellSubtype),
                                               col=list(cellSubset=cellcolors)
            ),
    )
  })
  dev.off()
  
  ###
  #Membership plot with colors based on single clusters at a time
  data_mem_tidy <- pivot_longer(data_mem, 1:colnum, names_to="sample")
  names(data_mem_tidy)[grep(paste("clus", h,sep="."), names(data_mem_tidy))] <- "membership"
  
  #Add in cellSubtype from files dataframe
  data_mem_tidy$sample <- gsub(".rpkm", "", data_mem_tidy$sample)
  data_mem_tidy <- left_join(data_mem_tidy,  select(files, sample, cellSubtype), by="sample")
  data_mem_tidy$cellSubtype <- factor(data_mem_tidy$cellSubtype, levels=smpl_order, ordered=T)
  
  data_mem_tidy <- data_mem_tidy[order(data_mem_tidy$cellSubtype),]
  
  data_mem_tidy_mean <- data_mem_tidy %>% group_by(cluster, gene, cellSubtype) %>% summarise(mean=mean(value), membership = mean(membership),gene=gene, cluster=cluster)
  #order the dataframe by score
  #data_mem_tidy <- data_mem_tidy[order(data_mem_tidy$membership),]
  #set the order by setting the factors using forcats
  levels(data_mem_tidy$gene) = fct_inorder(data_mem$gene)
  
  #Plot 
  ggplot(data_mem_tidy_mean, aes(x=cellSubtype, y=mean)) +
    geom_line(aes(color=membership, group=gene, alpha=membership)) +
    scale_color_gradientn(colors=c('blue1','gold')) +
    stat_summary(aes(group=1), fun = "mean", geom="line") + 
    facet_wrap(vars(cluster)) +
    theme(axis.text = element_text(angle=90, size=6))
  ggsave(filename=here(output_dir,paste(subset, "_member_gg_fuzz_facetbyclus", h, ".pdf", sep="")))
}

#Faceted membership plot with correct memberships for each facet
data_mem_tidy = NULL
data_mem_tidy <- pivot_longer(data_mem, 1:colnum, names_to="sample") %>% pivot_longer(cols=starts_with("clus."), names_to="mem_clus", values_to="membership")

data_mem_tidy$sample <- gsub(".rpkm", "", data_mem_tidy$sample)
data_mem_tidy <- left_join(data_mem_tidy,  select(files, sample, cellSubtype), by="sample")
data_mem_tidy$cellSubtype <- factor(data_mem_tidy$cellSubtype, levels=smpl_order, ordered=T)

data_mem_tidy <- data_mem_tidy[order(data_mem_tidy$cellSubtype),]

data_mem_tidy_mean <- data_mem_tidy %>% group_by(cluster, gene, mem_clus,cellSubtype) %>% summarise(mean=mean(value), membership = mean(membership),gene=gene, cluster=cluster)
levels(data_mem_tidy$gene) = fct_inorder(data_mem$gene)

matchr <- function(x, y){
  return(grepl(x, y))
}
index = NULL
for (i in 1:length(data_mem_tidy_mean$cluster)){
  index[i] <- matchr(data_mem_tidy_mean$cluster[i], data_mem_tidy_mean$mem_clus[i])  
}
data_mem_tidy_mean  <- data_mem_tidy_mean[index,]

#Plot 
pmem<-ggplot(data_mem_tidy_mean, aes(x=cellSubtype, y=mean)) +
  geom_line(aes(color=membership, group=gene, alpha=membership)) +
  scale_color_gradientn(colors=c('blue1','gold')) +
  stat_summary(aes(group=1), fun = "mean", geom="line") + 
  facet_wrap(vars(cluster), ncol=1) #+
  #theme(axis.text = element_text(angle=90, size=6))
pmem
ggsave(filename=here(output_dir,paste(subset, "_member_gg_fuzz_facet.pdf", sep="")), width=3, height=6.8)

#Combined plot
ggarrange(pline2+ theme(axis.text = element_blank()), 
          pmem+ theme(axis.text = element_blank()),
          widths = c(1,2)
          ) 
ggsave(filename=here(output_dir,paste(subset, "line_fuzz_facet.pdf", sep="")), width=6, height=6.8)

#####################################################################################
#Saving diffFile with appended cluster information of selected K

cluster_final <- select(dataDiff, starts_with(paste("k.",k,".", sep="")))
data_save <- cbind(dataDiff[,!grepl("clus", names(dataDiff))], cluster_final)
write.table(data_save, file=here(output_dir, paste0(subset,"_diff_fcm_k.",k,".txt")), sep="\t", row.names=F, quote=F)
