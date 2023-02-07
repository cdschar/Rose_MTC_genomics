
setwd("~/Documents/Emory/Boss/data/RNA")

library(tidyverse)
library(here)
library(ggpubr)

#Setting directories
# cts_input_dir <- "coverage"
man_input_dir <- "/home/jim/Documents/Emory/Boss/data"
diff_input_dir <- "naive_v_mem_overlap"
output_dir <- "barplots/fig2"

#Loading data and sample manifest
RNA_diff <- read.delim(here(diff_input_dir, "diff.glm.Human_Tcell_RNA.rev.txt"), header=T, sep="\t")
samples <- read.table(paste(man_input_dir, "RNAseq.sample.manifest.txt", sep="/"), sep = "\t", header = T, as.is = T)
samples <- samples[grepl("CM|EM|EMRA|Nav", samples$group),]
samples <- samples[!grepl("_stim", samples$group),]
samples <- samples[!grepl("CD4_EMRA", samples$group),]

#Selecting just rpkm columns
rpkm <- RNA_diff[, c(1, 2, match(paste0(samples$sample, ".rpkm"), colnames(RNA_diff)))]

test <- pivot_longer(rpkm, cols=ends_with(".rpkm"), names_to = "sample", values_to = "rpkm")
test$sample <- gsub(".rpkm", "", test$sample)
test2 <- left_join(test, select(samples, sample, group, cellType, tissue, stim, cellSubset), by="sample")

# gene=c("SELL", 
#        "CMKLR1",
#        "PECAM1",
#        "ICAM1",
#        "ITGB1",
#        "ITGA1",
#        "CX3CR1",
#        "XCL2",
#        "XCL1",
#        "VCAM1",
#        "ITGAL",
#        "SELP",
#        "CXCR3",
#        "CCL5",
#        "CXCR4",
#        "IL6R",
#        "CCR2",
#        "CCR4",
#        "CCR5",
#        "CCL5",
#        "CCR6",
#        "CD69",
#        "CCR7",
#        "S1PR1",
#        "S1PR5",
#        "S1PR2"
#        )

#gene=c("ITGAE", "CCR8", "CCR9", "CXCL9", "CCL3", "CCL4", "CXCR5")
#gene=c("FCGR3A", "FCGR3B", "SH2D1B", "KLRF1", "KLRC1", "KLRC2", "KLRC3", "LYN", "VAV3", "MYBL1", "KLRB1", "MYD88", "NKG7", "NCAM1", "ITGB2", "KIR2DL3", "KIR3DL1", "KIR3DL2", "SH3BP2", "IL2RA", "IL7R", "CD28", "MYC", "IL6R")
#gene=c("CCDC167","UST","ANTXR2","CARD17","PHLDA1","HNRNPLL","ATP2B1","IFNG−AS1","NABP1","CENPM","WDR86−AS1","CD40LG","CCR4","AP3M2","AHR","LRIG1","COTL1","CD82","TRADD","TRAT1","LINC00426","PDE4B","LOC100130476","FAM45A","TP53INP1","ELOVL5","GPR171","GALM","CYB561","GLB1","PTGER2","SLAMF1","REEP3","MDFIC","RTKN2","CXCR3","LIMS1","LINC00239","RILPL2","IGFBP4","TNFAIP3","NINJ2","SAMSN1","ANXA1","KLF6","SIAH2","GPRIN3","TYMP","LPAR6","SUMO4","DPP4","DUSP16","MTM1","FAAH2","TNFRSF13C","INPP4B","ICOS","CD28","TIAM1","GPR15","TTC39C−AS1","CLDND1","TIFA","CCR8","F5","CTLA4","ZNF165","ACVR1","EPHA4","GATA3","CAPG","CD69","PI16","NPDC1","GPR183","ZC3H12D","MAP3K1","CDC14A")
gene=c("IL2RA", "IL2RB")
gene=c("IFNG", "CTLA4", "ICOS", "CD28", "AHR", "IL2RA", "LDHB", "S1PR5", "PDP1", "ELOVL5")

#Names of groups from sample manifest
grpOrder = c(
  "CD4_Nav_blood_unstim", 
  # "CD4_Nav_blood_stim", 
  "CD4_CM_blood_unstim", 
  # "CD4_CM_blood_stim", 
  "CD4_EM_blood_unstim",
  # "CD4_EM_blood_stim",
  "CD8_Nav_blood_unstim",
  #"CD8_Nav_blood_stim",
  "CD8_CM_blood_unstim",
  #"CD8_CM_blood_stim",
  "CD8_EM_blood_unstim",
  #"CD8_EM_blood_stim",
  "CD8_EMRA_blood_unstim"
  #"CD8_EMRA_blood_stim"
)


#Comparisons to t-test
my_comparisons <- list( c("CD8_Nav", "CD8_CM"),
                        c("CD8_CM", "CD8_EM"),
                        c("CD8_EM", "CD8_EMRA"),
                        c("CD8_Nav", "CD8_EM"),
                        c("CD8_Nav", "CD8_EMRA"),
                       c("CD4_Nav", "CD4_CM"),
                       c("CD4_CM", "CD4_EM"),
                       c("CD4_Nav", "CD4_EM")
)

#Color pallete (matched to group names)
custom_pal <- c(
                "CD4_Nav"="#90A4ADFF",
                "CD4_CM"="#2096F2FF", 
                "CD4_EM"="#FFB74CFF", 
                "CD8_Nav"="#455964FF", 
                "CD8_CM"="#0C46A0FF", 
                "CD8_EM"="#E55100FF", 
                "CD8_EMRA"="#870D4EFF")
#pie(rep(1,length(custom_pal)),col=custom_pal, labels=names(custom_pal))

for (i in 1:length(gene)){  

  plot_data <- subset(test2, SYMBOL==gene[i])
  
  if(length(plot_data$SYMBOL)==0){
    print(paste("No data for ", gene[i], sep=""))
    next
  }
  
  #Clean up group names
  plot_data$group <- gsub("_blood_unstim", "", plot_data$group)
  grpOrder <- gsub("_blood_unstim", "", grpOrder)
  grpOrder <- gsub("_blood_stim", "", grpOrder)
  plot_data$group <- factor(plot_data$group, levels=grpOrder, ordered=T)
  
  p <- ggbarplot(plot_data, x = "group", y = "rpkm",
                 fill = "group", 
                 palette = custom_pal,
                 x.text.angle= 90,
                 add = c("mean_sd","jitter")
                )
  #p
  
  #p +  stat_compare_means(label = "p.signif", method = "t.test",
  #                        ref.group = "CD8_Nav_blood_unstim")  
  ggpar(p+ stat_compare_means(comparisons=my_comparisons, label="p.signif") + 
          ggtitle(gene[i]), legend="right", legend.title = "Subtype")
  ggsave(here(output_dir, paste(gene[i], ".pdf", sep="")), width=4.5, height=5)
  ggsave(here(output_dir, paste(gene[i], ".jpg", sep="")), width=4.5, height=5)
}   
