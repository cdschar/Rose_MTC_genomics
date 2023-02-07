setwd("~/Documents/Emory/Boss/data/RNA")

library(tidyverse)
library(here)
library(ggpubr)

#Setting directories
# cts_input_dir <- "coverage"
man_input_dir <- "/home/jim/Documents/Emory/Boss/data"
diff_input_dir <- "naive_v_mem_overlap"
output_dir <- "barplots/fig6"

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

gene <- c("CXCR3", "IL4", "MSC", "TBX21","EOMES","AHR","HIF1A")

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

my_comparisons <- list( c("CD8_Nav_unstim", "CD8_Nav_stim"),
                        c("CD8_CM_unstim", "CD8_CM_stim"),
                        c("CD8_EM_unstim", "CD8_EM_stim"),
                        #c("CD8_EMRA_unstim", "CD8_EMRA_stim"),
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