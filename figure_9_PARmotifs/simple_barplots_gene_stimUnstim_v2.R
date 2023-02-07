
setwd("~/Documents/Emory/Boss/data/RNA")

library(tidyverse)
library(here)
library(ggpubr)

#Setting directories
# cts_input_dir <- "coverage"
man_input_dir <- "/home/jim/Documents/Emory/Boss/data"
diff_input_dir <- "naive_v_mem_overlap"
output_dir <- "barplots/fig5"

#Loading data and sample manifest
RNA_diff <- read.delim(here(diff_input_dir, "diff.glm.Human_Tcell_RNA.rev.txt"), header=T, sep="\t")
samples <- read.table(paste(man_input_dir, "RNAseq.sample.manifest.txt", sep="/"), sep = "\t", header = T, as.is = T)
samples <- samples[grepl("CM|EM|Nav", samples$group),]
#samples <- samples[!grepl("_stim", samples$group),]
samples <- samples[!grepl("EMRA", samples$group),]

#Selecting just rpkm columns
rpkm <- RNA_diff[, c(1, 2, match(paste0(samples$sample, ".rpkm"), colnames(RNA_diff)))]

test <- pivot_longer(rpkm, cols=ends_with(".rpkm"), names_to = "sample", values_to = "rpkm")
test$sample <- gsub(".rpkm", "", test$sample)
test2 <- left_join(test, select(samples, sample, group, cellType, tissue, stim, cellSubset), by="sample")
test2 <- mutate(test2, cellSubtype = paste(cellType, cellSubset, sep="_"))

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
#        "S1PR2",
#        )

# gene=c(
#       #"S1PR1", "S1PR2", "S1PR5", "CD69",
#       # "ITGAL", "ICAM1", "ITGA4", "VCAM1",
#        #"XCL2", "XCL1",
#        "SELP", "ITGB1", "CCR6", "CCR4", "CCR2",
#        "SELL", "GPR15",
#        )

#gene=c("NKG7")
#gene=c("CCDC167","UST","ANTXR2","CARD17","PHLDA1","HNRNPLL","ATP2B1","IFNG−AS1","NABP1","CENPM","WDR86−AS1","CD40LG","CCR4","AP3M2","AHR","LRIG1","COTL1","CD82","TRADD","TRAT1","LINC00426","PDE4B","LOC100130476","FAM45A","TP53INP1","ELOVL5","GPR171","GALM","CYB561","GLB1","PTGER2","SLAMF1","REEP3","MDFIC","RTKN2","CXCR3","LIMS1","LINC00239","RILPL2","IGFBP4","TNFAIP3","NINJ2","SAMSN1","ANXA1","KLF6","SIAH2","GPRIN3","TYMP","LPAR6","SUMO4","DPP4","DUSP16","MTM1","FAAH2","TNFRSF13C","INPP4B","ICOS","CD28","TIAM1","GPR15","TTC39C−AS1","CLDND1","TIFA","CCR8","F5","CTLA4","ZNF165","ACVR1","EPHA4","GATA3","CAPG","CD69","PI16","NPDC1","GPR183","ZC3H12D","MAP3K1","CDC14A","PERP","GBP5","RTKN2","PASK","KCNN4","DPYSL2","PDP1","CDCA7","TMEM64","UNG","RNF214","CDCA4","CD84","OPTN","FAM45A","CLDND1","TP53INP1","NABP1","STAM","CKS2","PYHIN1","STX11","ELOVL5","CD2","IL2RA","ELOVL1","CARD17","TMEM14A","NDC80","LOC100419583","MIEN1","AQP3","GPR15","ICOS","CTLA4","F5","FBXO5","MIR21","RNF19A","SYT11","GPRIN3","ITGB1","PHTF2","SLC9A9","MYO5A","ATXN1","WEE1","SELENOS","CCR4","LRIG1","CERKL","DUSP16","CCDC167","CASP8","FAM160B1","MT1F","C16orf87","CENPM","SMC4","NINJ2","ZC3H12D","HMGB2","TICAM1","GPR25","LIMS1","LIMS3","C1orf216","ZC2HC1A","TIFA")
#gene=c("IL1A","IL1B","IL1RL1","IL2","IL3","IL4","IL5","IL6","IL7","IL8","IL9","IL10","IL11","IL12","IL13","IL14","IL15","IL16","IL17A","IL17B","IL17C","IL17D","IL17E","IL17F","IL17E","IL18","IL19","IL20","IL21","IL22","IL23","IL24","IL25","IL26","IL27","IL28","IL29","IL30","IL31","IL32","IL33","IL34","IL35","IL36","IL37","IL38","IFNG", "TNF", "TGFB1")
#gene=c("NCOA2", "AHR", "HIF1A", "MSC", "MSC-ASI1")
#gene=c("AHR", "ASB2", "ICOS", "CCL20", "CTLA4", "CD80", "IRF8", "RORA", "CLEC2B", "P2RY8", "CTSW", "HCST", "PELI1", "CD28", "UCP2", "NCOA2")
#gene=c("PSCM4","HNRNPU","EIF4G2","C1QBP","GNL3","HDAC2","CCT7","IMPDH2","DHX15","HNRNPC","CCT5","G3BP1","ABCE1","SRSF1","VDAC1","CCT4","DEK","CCT2","XPOT","PRPS2","TFDP1")
#gene=c("LGALS1","LGALS3","IFNG","FOXP3","EOMES","PRDM1","LAG3","CTLA4","IL23R","APOBEC3D","CCR4","CD70","PRF1","CCR5", "TOX", "TRAT1", "ANXA1", "P2RY8", "S100A4", "APOBEC3G", "APOBEC3H")
#gene=c("IL7R")
#gene=c("JUNB","FOS","TCF7","LEF1","TCF3","PRDM1","STAT3", "GATA1", "GATA2", "GATA3")
#gene=c("TBX21","BACH1","NFE2L2","BACH2","IRF8","IRF2","STAT5A","STAT5B","RORC", "EOMES")
gene=c("IL22", "IL4", "IL1A", "ADGRG1", "PLEKHG3", "CCL5", "EOMES", "LGALS3", "CXCR6", "MSC", "PHLDA2", "GZMH", "NKG7", "CXCR3", "TBX21", "MNDA")
#gene=c("IRF4","IRF3","BATF","NFATC1","NFATC2","NFATC3","NFKB1","NFKB2")

gene = unique(gene)

#Names of groups from sample manifest
grpOrder = c(
  "CD4_Nav_blood_unstim", 
   "CD4_Nav_blood_stim", 
  "CD4_CM_blood_unstim", 
   "CD4_CM_blood_stim", 
  "CD4_EM_blood_unstim",
   "CD4_EM_blood_stim",
  "CD8_Nav_blood_unstim",
  "CD8_Nav_blood_stim",
  "CD8_CM_blood_unstim",
  "CD8_CM_blood_stim",
  "CD8_EM_blood_unstim",
  "CD8_EM_blood_stim",
  "CD8_EMRA_blood_unstim",
  "CD8_EMRA_blood_stim"
)


#Comparisons to t-test
my_comparisons <- list( c("CD8_Nav_unstim", "CD8_Nav_stim"),
                        c("CD8_CM_unstim", "CD8_CM_stim"),
                        c("CD8_EM_unstim", "CD8_EM_stim"),
                        c("CD8_EMRA_unstim", "CD8_EMRA_stim"),
                        c("CD4_Nav_unstim", "CD4_Nav_stim"),
                        c("CD4_CM_unstim", "CD4_CM_stim"),
                        c("CD4_EM_unstim", "CD4_EM_stim")
)

#Color pallete (matched to group names)
custom_pal <- c("CD4_Nav"="#90A4ADFF",
                "CD4_CM"="#2096F2FF", 
                "CD4_EM"="#FFB74CFF", 
                "CD8_Nav"="#455964FF", 
                "CD8_CM"="#0C46A0FF", 
                "CD8_EM"="#E55100FF", 
                "CD8_EMRA"="#870D4EFF",
                "stim" ="grey30",
                "unstim"="white")
#pie(rep(1,length(custom_pal)),col=custom_pal, labels=names(custom_pal))

for (i in 1:length(gene)){  

  plot_data <- subset(test2, SYMBOL==gene[i])
  
  #Clean up group names
  plot_data$group <- gsub("_blood", "", plot_data$group)
  grpOrder <- gsub("_blood", "", grpOrder)
  plot_data$group <- factor(plot_data$group, levels=grpOrder, ordered=T)
  
  p <- ggbarplot(plot_data, x = "group", y = "rpkm",
                 fill = "stim",
                 color = "cellSubtype",
                 palette = custom_pal,
                 x.text.angle= 90,
                 add = c("mean_sd","jitter")
                ) +
              theme(legend.key = element_rect(fill = "white")) +
              guides(color=guide_legend(override.aes = list(fill='white')))
  #p
  
  #p +  stat_compare_means(label = "p.signif", method = "t.test",
  #                        ref.group = "CD8_Nav_blood_unstim")  
  ggpar(p + stat_compare_means(comparisons=my_comparisons, label="p.signif") +
          ggtitle(gene[i]), legend="right", legend.title = "Subtype")
  ggsave(here(output_dir, paste(gene[i], "stimunstim.pdf", sep="")), width=7, height=5)
  ggsave(here(output_dir, paste(gene[i], "stimunstim.jpg", sep="")), width=7, height=5)
}   
