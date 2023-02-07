#Heatmap for leading edge genes
#20Mar20
#jrose

setwd("~/Documents/Emory/Boss/data/RNA")

library(here)
library(tidyverse)
library(ComplexHeatmap)

input_dir <- "GSEA/PageRank_Networks"
output_dir <- "GSEA/Analysis/PageRank_Networks"

#List of directories for xls files
comparison_dirs <- c("CD8_CM_EM_PgRnk_hits.GseaPreranked.1634131408035")

#Specifying geneset used in GSEA
geneset <- "HIF1A"
TF <- "HIF1A"

#Creating list of core-enriched genes for heatmap
core_enrch <- list()

for (i in 1:length(comparison_dirs)){
  tmp <- read.table(here(input_dir, comparison_dirs[i], paste(geneset, ".xls", sep="")), sep="\t", header=T) %>%
      subset(CORE.ENRICHMENT=="Yes") %>%
      mutate(Comp=comparison_dirs[i])
  core_enrch[[i]] <- tmp %>% select(PROBE, Comp)
}

enriched <- data.frame()
for (i in 1:length(core_enrch)){
  enriched <- bind_rows(enriched, core_enrch[[i]])
}
genes <- unique(enriched$PROBE)

#Load in RNA data sample manifest
fqDir = "/home/jim/Documents/Emory/Boss/data/"
fqFile = "RNAseq.sample.manifest.txt";
files = read.table(paste0(fqDir, fqFile), sep = "\t", header = T, as.is = T);
files = files[grepl("EM|CM|Nav", files$group),]
files = files[files$stim=="unstim",]
files = files[!grepl("CD4_EMRA", files$group),]

#read in data file
dataDir = "/home/jim/Documents/Emory/Boss/data/";
dataFile = "Tcell_RNA.combined.RPKM.csv";
data = read.table(file=paste0(dataDir, dataFile), header=T, sep=",", comment.char="", quote="");

#Extract rpkm and average by group
rpkm <- data[,grepl(paste(files$sample, collapse = "|"),colnames(data)) & 
             !grepl("Stim", colnames(data))
             ]
rpkm <- rpkm[,!grepl("CD4_EMRA", colnames(rpkm))]
rpkm$SYMBOL <- data$SYMBOL

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

rpkm_z <- apply(select(mean_rpkm, -SYMBOL), c(1), scale)
rpkm_z <- t(rpkm_z)
rownames(rpkm_z) <- mean_rpkm$SYMBOL
colnames(rpkm_z) <- colnames(mean_rpkm)[-c(1)]
rpkm_z_sub <- rpkm_z[rownames(rpkm_z) %in% genes,]

#Annotation File
# ph_anno <- select(files, cellSubset, cellType)
# rownames(ph_anno) <- files$sample

#heatmap
pdf(file=here(output_dir, paste(geneset, "leading_edge_all_heatmap.pdf", sep="_")), useDingbats = F)
Heatmap(rpkm_z_sub,
        show_row_names = T,
        show_column_names = T,
        row_names_gp = gpar(fontsize = 6),
        # right_annotation = rowAnnotation(gene = anno_mark(at = lbl_index, labels = genes_highlight, labels_gp=gpar(fontsize=10))
        #                                  ),
        top_annotation = HeatmapAnnotation(TF_expression=anno_barplot(mean_rpkm[grep(paste0("^",TF, "$"),mean_rpkm$SYMBOL),-c(1)] %>% unlist())
                                          )
        )
dev.off()
