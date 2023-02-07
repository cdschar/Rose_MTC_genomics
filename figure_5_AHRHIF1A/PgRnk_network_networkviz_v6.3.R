#Script for creating a network viz figure of TF gene networks
#v3.1 Trying to add PgRnk score as size of nodes compared to v2.0
#v4.0 Either showing only TF from network or using shapes to distinguish
#v5 Exlcudes non-TF network genes and incorporates more network connections
#v6 Compares CM to EM in terms of logFC
  #v6.3 adding direction of pagerank to plot

setwd("~/Documents/Emory/Boss/data/ATAC")

library(tidyverse)
library(scales)
library(network)
library(ggnetwork)
library(here)
library(RColorBrewer)

input_dir <- "PageRank/Run1/Network"
RNA_input_dir <- "/home/jim/Documents/Emory/Boss/data/ATAC/PageRank/Run1/Network/network_viz/diff.significant.glm.Human_Tcell_RNA.rev.txt"
output_dir <- "PageRank/Run1/Network/network_viz"
PgRnk_input_dir <- "PageRank/Run1"

cutoff=2
#^Min rpkm required across all samples in a group to be included in plot
diff_cut = 1.5
#^Min logFC requirement (abs) to be included in plot
min_plt = -4
#^min for logFC color scale
max_plt = 4
#^max

#TFs <- c("MSC", "BATF", "GATA3", "GFI1","LEF1", "TCF7", "TBX21", "EOMES", "AHR", "HIF1A", "ARID5B", "IRF4", "IRF8")
#TFs <- c("MSC", "LEF1", "TCF7", "IRF8", "TBX21", "EOMES", "FOXO1", "BATF", "IRF4", "ID2", "STAT4", "PRDM1", "NFATC1", "BCL6", "STAT3", "NFATC2", "ZEB1", "ZEB2")
#TFs <- c("MSC", "AHR", "HIF1A", "BATF", "LEF1", "TBX21", "IRF8", "TCF7", "EOMES")
#TFs <- c("EGR4","MSC", "CEBPD","FOSL2","TBX21","EOMES","ASCL2","AHR","HIF1A","BHLHE40", "BATF", "LEF1", "TCF7", "RUNX3", "BCL6", "STAT5B", "SMAD3")
#TFs <- c("MSC", "TCF7", "LEF1", "TBX21", "HIF1A", "AHR", "BATF","FOSL2")

#TFs <-c("AHR","EGR4","KLF15","DMBX1","SMAD3","HOXA11","HES7","GSC2","HES2","GLI3","CEBPD","ZIC2","HIF1A","GFI1B","SMAD1","NKX2-4","ARID3C","MNT","ZSCAN16","EOMES","HOXB7","HES5","HOXB8","ZNF333","E2F5","AC226150.2","HEY1","ONECUT1","HOXC8","ZHX1","MIXL1","GSX2","NR1H4","HSF4","GLIS1","HEYL","TOPORS","LBX1","ATOH7","POU5F1","ARID5B","SMAD2","POU1F1","TWIST2","ONECUT2","HIC2","NR1I3","CREB3L4","ONECUT3")
#^Generated from top 50 CD8 EM/CD8 CM PgRnk LogFC

#TF <- "MSC"
TFs <- c("AHR")

fqDir = "/home/jim/Documents/Emory/Boss/data/"
fqFile = "RNAseq.sample.manifest.txt";
files = read.table(paste0(fqDir, fqFile), sep = "\t", header = T, as.is = T);
files = files[grepl("EM|CM|Nav", files$group),]
#files = files[!grepl("CD8", files$group),]
files = files[!grepl("CD4_EMRA", files$group),]
files = files[which(files$include),]
files = files[files$stim!="stim",]
files = mutate(files, cellSubtype = paste(cellType, cellSubset, stim, sep="_"))

#TF <- TFtxt[i,1] %>% as.character()


#Create TF specfic output sub dir
# if(!dir.exists(here(output_dir, TF))){
#   dir.create(here(output_dir, TF))
# }
output_dir <- "PageRank/Run1/Network/network_viz/PgRnk_TF_Networks"

lst <- list(
  CD8_Nav = read.delim(file=here(input_dir, "CD8_N_ATAC_network.tsv"), header=F, col.names=c("Target", "TFs")) %>% mutate(subtype="CD8_Nav"),
  CD8_EM = read.delim(file=here(input_dir, "CD8_EM_ATAC_network.tsv"), header=F, col.names=c("Target", "TFs")) %>% mutate(subtype="CD8_EM"),
  CD8_CM = read.delim(file=here(input_dir, "CD8_CM_ATAC_network.tsv"), header=F, col.names=c("Target", "TFs")) %>% mutate(subtype="CD8_CM"),
  CD8_EMRA = read.delim(file=here(input_dir, "CD8_EMRA_ATAC_network.tsv"), header=F, col.names=c("Target", "TFs")) %>% mutate(subtype="CD8_EMRA"),
  CD4_Nav = read.delim(file=here(input_dir, "CD4_N_ATAC_network.tsv"), header=F, col.names=c("Target", "TFs")) %>% mutate(subtype="CD4_Nav"),
  CD4_EM = read.delim(file=here(input_dir, "CD4_EM_ATAC_network.tsv"), header=F, col.names=c("Target", "TFs")) %>% mutate(subtype="CD8_EM"),
  CD4_CM = read.delim(file=here(input_dir, "CD4_CM_ATAC_network.tsv"), header=F, col.names=c("Target", "TFs")) %>% mutate(subtype="CD8_CM"),
  CD8_NavStim = read.delim(file=here(input_dir, "CD8_NStim_ATAC_network.tsv"), header=F, col.names=c("Target", "TFs")) %>% mutate(subtype="CD8_NavStim"),
  CD8_EMStim = read.delim(file=here(input_dir, "CD8_EMStim_ATAC_network.tsv"), header=F, col.names=c("Target", "TFs")) %>% mutate(subtype="CD8_EMStim"),
  CD8_CMStim = read.delim(file=here(input_dir, "CD8_CMStim_ATAC_network.tsv"), header=F, col.names=c("Target", "TFs")) %>% mutate(subtype="CD8_CMStim"),
  CD8_EMRAStim = read.delim(file=here(input_dir, "CD8_EMRAStim_ATAC_network.tsv"), header=F, col.names=c("Target", "TFs")) %>% mutate(subtype="CD8_EMRAStim"),
  CD4_NavStim = read.delim(file=here(input_dir, "CD4_NStim_ATAC_network.tsv"), header=F, col.names=c("Target", "TFs")) %>% mutate(subtype="CD4_NavStim"),
  CD4_EMStim = read.delim(file=here(input_dir, "CD4_EMStim_ATAC_network.tsv"), header=F, col.names=c("Target", "TFs")) %>% mutate(subtype="CD4_EMStim"),
  CD4_CMStim = read.delim(file=here(input_dir, "CD4_CMStim_ATAC_network.tsv"), header=F, col.names=c("Target", "TFs")) %>% mutate(subtype="CD4_CMStim")
)


#lst <- lst[!grepl("Stim", names(lst))]

#Function for extracting list of targets for each subtype of selected transcription factors (TFs)
##############################################################################
extract_multiTF_targets <- function(x){
  output_df <- data.frame()
  for (z in 1:length(TFs)){
    TF <- TFs[z]
    targets <- x$Target[grep(TFs[z],x$TFs)]
    tmp_df <- data.frame(network = targets)
    tmp_df <- tmp_df %>% mutate(TF=TF) %>% arrange(network)
    tmp_df <- tmp_df[!grepl(TF, tmp_df$network),]
    output_df <- bind_rows(output_df, tmp_df) 
  }
  return(output_df)
}
TF_lst <- map(lst, extract_multiTF_targets)

#Getting vector of all unique targets for this TF (regardless of subtype)
##############################################################################
All_targets <- vector()
for (i in 1:length(TF_lst)){
  All_targets <- c(All_targets, as.character(TF_lst[[i]]$network))
}
All_targets <- unique(All_targets)


#Filtering target list by DEG and rpkm cutoff
############################################################################
RNA_diff <- read.delim(RNA_input_dir, header=T, sep="\t") %>% drop_na(SYMBOL)

rpkm <- RNA_diff[,grepl(paste(files$sample, collapse = "|"),colnames(RNA_diff)) & !grepl("Stim", colnames(RNA_diff))]
rpkm <- rpkm[,!grepl("CD4_EMRA", colnames(rpkm))]
rownames(rpkm) <- RNA_diff$SYMBOL

#Removing genes with low expression
source("/home/jim/Code/Boss/bistools/bisTools.R");
expr_test = apply(rpkm[,grepl(paste(files$sample, collapse = "|"), colnames(rpkm))], 1, function(x) filter_detected(x, files$group, cutoff) )
print('Genes before rpkm filter')
dim(rpkm)
rpkm <- rpkm[expr_test,]
print('Genes after rpkm filter')
dim(rpkm)

All_targets <- All_targets[All_targets %in% RNA_diff$SYMBOL[RNA_diff$SYMBOL %in% rownames(rpkm)]]

#Creating edges
# edges_df <- data.frame(network=All_targets) %>% mutate(TF=TF) %>% arrange(network)
# edges_df <- edges_df[!grepl(TF, edges_df$network),]
# ^Already done in extract function in this script!

#Getting logFC of network genes (vs Nav)
RNA_network <- rbind(RNA_diff[RNA_diff$SYMBOL %in% All_targets, ], RNA_diff[RNA_diff$SYMBOL %in% TFs,])

#PageRank Scores
d <- read.table(here(PgRnk_input_dir,"TFRanks.tsv"),header=T,row.names=1)

RNA_network <- RNA_network[RNA_network$SYMBOL %in% rownames(d) | RNA_network$SYMBOL %in% TFs,]

#Single Subset Network LogFC values
logFC_lst <- list()
logFC_lst$CD8_EM <- RNA_network %>% 
                select(SYMBOL, CD8_CM_blood_unstim.v.CD8_EM_blood_unstim.logFC) %>%
                subset(abs(CD8_CM_blood_unstim.v.CD8_EM_blood_unstim.logFC)>diff_cut | SYMBOL %in% TFs) %>%
                arrange(SYMBOL) %>%
                dplyr::rename(logFC = ends_with('.logFC')) %>%
                mutate(logFC = -logFC)
logFC_lst$CD8_CM <- RNA_network %>% 
                select(SYMBOL, CD8_CM_blood_unstim.v.CD8_EM_blood_unstim.logFC) %>%
                subset(abs(CD8_CM_blood_unstim.v.CD8_EM_blood_unstim.logFC)>diff_cut | SYMBOL %in% TFs) %>%
                arrange(SYMBOL) %>%
                dplyr::rename(logFC = ends_with('.logFC'))
logFC_lst$CD8_EMRA <- RNA_network %>% 
                select(SYMBOL, CD8_CM_blood_unstim.v.CD8_EMRA_blood_unstim.logFC) %>%
                subset(abs(CD8_CM_blood_unstim.v.CD8_EMRA_blood_unstim.logFC)>diff_cut | SYMBOL %in% TFs) %>%
                arrange(SYMBOL) %>%
                dplyr::rename(logFC = ends_with('.logFC')) %>%
                mutate(logFC = -logFC)
logFC_lst$CD4_EM <- RNA_network %>% 
                select(SYMBOL, CD4_CM_blood_unstim.v.CD4_EM_blood_unstim.logFC) %>%
                subset(abs(CD4_CM_blood_unstim.v.CD4_EM_blood_unstim.logFC)>diff_cut | SYMBOL %in% TFs) %>%
                arrange(SYMBOL) %>%
                dplyr::rename(logFC = ends_with('.logFC'))%>%
                mutate(logFC = -logFC)
logFC_lst$CD4_CM <- RNA_network %>% 
                select(SYMBOL, CD4_CM_blood_unstim.v.CD4_EM_blood_unstim.logFC) %>%
                subset(abs(CD4_CM_blood_unstim.v.CD4_EM_blood_unstim.logFC)>diff_cut | SYMBOL %in% TFs) %>%
                arrange(SYMBOL) %>%
                dplyr::rename(logFC = ends_with('.logFC'))

logFC_passed_genes <- map(logFC_lst, "SYMBOL") %>% unlist() %>% unique() %>% as.character()

#Constant color scale
palette <- colorRampPalette(rev(brewer.pal(11, name="RdBu")))
color_scale <- scale_color_gradientn(colors=palette(50), limits=c(min_plt, max_plt), aesthetics = c("color", "fill"))
#size_scale <- scale_size_continuous(limits=c(0, max_plt))

layout = "kamadakawai"
#^Network layout algorithm


rescale <- function(nchar, low, high){
  min_d <- min(nchar)
  max_d <- max(nchar)
  rscl <- ((high-low)*(nchar-min_d))/(max_d-min_d)+low
  return(rscl)
}

for (i in 1:length(logFC_lst)) {
  #Creating node attributes
  subset <- names(logFC_lst)[[i]]
  ####logFC
  logFC <- as.matrix(logFC_lst[[i]]$logFC)
  rownames(logFC) <- logFC_lst[[i]]$SYMBOL
  logFC <- as.matrix(logFC[!duplicated(rownames(logFC)),1])
  
  for (z in 1:length(TFs)){
    if (length(grep(TFs[z], rownames(logFC)))==0){
      names_save <- rownames(logFC)
      logFC <- c(logFC, 0) %>% as.matrix()
      rownames(logFC) <- c(names_save, TFs[z])
    }
  }
  #logFC_edges <- logFC[rownames(logFC)!=TF,] %>% as.matrix()
  edges_df_tmp <- TF_lst[[subset]][TF_lst[[subset]]$network %in% rownames(logFC),]
  #edges_df_tmp$logFC <- logFC_edges[order(match(rownames(logFC_edges), edges_df_tmp$network)),]
  
  ####PgRnk
  if(grepl("CD4_CM", names(logFC_lst)[[i]])){
      denom=d[,"CD4_EM"]
  } else if (grepl("CD4_EM", names(logFC_lst)[[i]])) {
        denom=d[,"CD4_CM"]
  } else if (grepl("CD8_CM", names(logFC_lst)[[i]])) {
    denom=d[,"CD8_EM"]
  } else if (grepl("CD8_EM", names(logFC_lst)[[i]])) {
    denom=d[,"CD8_CM"]
  } else if (grepl("CD8_EMRA", names(logFC_lst)[[i]])) {
    denom=d[,"CD8_CM"]
  } else {
    break
  }
  
  PgRnk <- as.matrix(log(d[,names(logFC_lst)[[i]]] / denom, base=2)
                     )
  #PgRnk <- abs(PgRnk)
  #^Setting shape min at 1
  rownames(PgRnk) <- rownames(d)
  PgRnk <- PgRnk[rownames(PgRnk)%in%edges_df_tmp$network | rownames(PgRnk)%in%edges_df_tmp$TF,] %>%
            as.matrix()
  PgRnk_nas <- as.matrix(rep(1, length(which(!edges_df_tmp$network%in%rownames(PgRnk)))))
  rownames(PgRnk_nas) <- edges_df_tmp$network[!edges_df_tmp$network%in%rownames(PgRnk)]
  PgRnk_full <- as.matrix(c(PgRnk, PgRnk_nas))
  rownames(PgRnk_full) <- c(rownames(PgRnk), rownames(PgRnk_nas))
  
  ####Is TF?
  # TF_att <- rep(FALSE, length(PgRnk_full)) %>% as.matrix()
  # rownames(TF_att) <- rownames(PgRnk_full)
  # TF_att[rownames(TF_att)%in%rownames(PgRnk),1] = TRUE
  
  #Node list
  nodes <- data.frame(label=unique(c(as.character(edges_df_tmp$network), as.character(edges_df_tmp$TF)))) %>% rowid_to_column("id")
  nodes$label <- as.character(nodes$label)
  
  #Creating network object
  if (length(edges_df_tmp$network)>2){
      ############There is an error here!! Directions are reversed???######################
      edges <- edges_df_tmp %>% left_join(nodes, by=c("TF"="label")) %>% dplyr::rename(from=id)
      edges <- edges %>% left_join(nodes, by=c("network"="label")) %>% dplyr::rename(to=id) %>%
          select(from, to)
      tf.net <- network(edges, vertext.attr=nodes, matrix.type="edgelist", ignore.eval=FALSE, directed = T, multiple = F)
  } else {
    print(paste("Not enough regulated DEG for", names(logFC_lst)[[i]], "to create network"))
    next
  }
  
  vert_name_index <- network.vertex.names(tf.net) 
  network.vertex.names(tf.net) <- nodes$label[order(match(nodes$id, vert_name_index))]
  tf.net %v% "logFC_vert" <-  logFC[network.vertex.names(tf.net),]
  tf.net %v% "PgRnk" <-  PgRnk_full[network.vertex.names(tf.net),]
  #tf.net %v% "TF" <- TF_att[network.vertex.names(tf.net),]
  
  ELSE <- TRUE
  
  network_plot <- ggnetwork(tf.net, 
                            layout = layout, 
                            weights="PgRnk", 
                            arrow.gap=0.03
                            ) %>% 
                  mutate(PgRnk_Sig = case_when(
                                      PgRnk >=0 ~ "Pos",
                                      ELSE ~ "Neg"
                                              )
                        )
  
  ggplot(data = network_plot,
         aes(x, y, xend = xend, yend = yend)) +
    geom_edges(color="grey55", arrow = arrow(length = unit(4, "pt"))) +
    #geom_edges(aes(color=logFC, linetype=TF), size=1) +
    geom_nodes(aes(shape=PgRnk_Sig,size=rescale(abs(PgRnk), 1, 5)*1.1), color="black") +
    geom_nodes(aes(shape=PgRnk_Sig,size=rescale(abs(PgRnk), 1, 5),
                   color=logFC_vert,
                   fill = logFC_vert
                   #shape=TF
                   )
                ) +
    geom_nodelabel(aes(label = vertex.names),
                  size = 3, vjust = -0.5) +
    color_scale + 
    #size_scale +
    #scale_colour_distiller(type="div", palette="RdBu", breaks=c(-8, -5, 0, 5, 8)) +
    scale_shape_manual(values=c('Pos'=17, 'Neg'=25)) +
    #scale_linetype_manual(values=c('TRUE'='solid', 'FALSE'='dashed')) +
    labs(color="RNA log2FC", size="PgRnk log2FC") +
    guides(fill="none") +
    xlim(c(-0.05, 1.05)) +
    theme_blank() +
    theme(legend.position = "bottom")
  ggsave(filename=here(output_dir, paste0(subset,"_PgRnk_network_v6.3.pdf")), height=6, width=6.5)
}


#Heatmaps

# library(ComplexHeatmap)
# rpkm$SYMBOL <- rownames(rpkm)
# mean_rpkm <- pivot_longer(rpkm, 
#                           cols=ends_with('.rpkm'),
#                           names_to=c("Sample", "Celltype", "Subtype"),
#                           names_pattern= "RNA_(.*)_(.*)_(.*).rpkm",
#                           values_to="rpkm"
#                           ) %>%
#   mutate(CellSubtype=paste(Celltype, Subtype, sep="_")) %>%
#   group_by(SYMBOL, CellSubtype) %>%
#   summarise(mean_rpkm=mean(rpkm)) %>%
#   pivot_wider(names_from=CellSubtype, values_from=mean_rpkm) %>%
#   ungroup()
# 
# rpkm_z <- apply(select(mean_rpkm, -SYMBOL), c(1), scale)
# rpkm_z <- t(rpkm_z)
# rownames(rpkm_z) <- mean_rpkm$SYMBOL
# colnames(rpkm_z) <- colnames(mean_rpkm)[-c(1)]
# rpkm_z_sub <- rpkm_z[rownames(rpkm_z) %in% logFC_passed_genes,]
# 
# genes_highlight <- c(
#   "AHR", 
#   #"DUSP5",
#   "ITGB1",
#   "GATA3",
#   "CCR4",
#   "LEF1",
#   "MYC",
#   "BATF",
#   "CCR7",
#   "TCF7",
#   "ITGB2",
#   "TRABD2A",
#   #"ARID5B",
#   "GZMK",
#   "HNRNPLL",
#   "PECAM1",
#   #"IL17RA",
#   "BCL2",
#   "FOXO1",
#   "VCL", 
#   "GFI1",
#   "NFATC2",
#   "PTPRJ",
#   "IL12RB1",
#   "APOBEC3G"
# )
# 
# lbl_index <- vector()
# for (i in 1:length(genes_highlight)){
#   if (genes_highlight[i]%in%rownames(rpkm_z_sub)) {
#     lbl_index[i] <- grep(genes_highlight[i], rownames(rpkm_z_sub))
#   }
# }
# 
# pdf(file=here(output_dir,TF,paste0(TF,"_heatmap.pdf")), height=6, width=6)
# Heatmap(rpkm_z_sub, 
#         show_row_names = F,
#         show_column_names = T,
#         row_names_gp = gpar(fontsize = 3),
#         right_annotation = rowAnnotation(gene = anno_mark(at = lbl_index, labels = genes_highlight, labels_gp=gpar(fontsize=10))
#                                          ),
#         top_annotation = HeatmapAnnotation(TF_expression=anno_barplot(mean_rpkm[grep(paste0("^",TF, "$"),mean_rpkm$SYMBOL),-c(1)] %>% unlist())
#                                           )
#         )
# dev.off()

