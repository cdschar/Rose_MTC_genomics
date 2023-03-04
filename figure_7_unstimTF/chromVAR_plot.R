#Refining plots of chromVAR results
#30Oct20
#jrose
setwd("~/Documents/Emory/Boss/data/ATAC")

library("chromVAR")
library("tidyverse")
library("ggrepel")
library("here")
library("ggpubr")

output_dir <- "chromVAR"

load(file=here("chromVAR/Tmem_unstim_dev.rda"))

var <- variability %>% mutate(rank = rank(-variability),neglog10p = -log(p_value_adj, base=10)) %>%
  arrange(rank)

labels = c("EOMES", "TBX21", "FOS::JUN", "JUNB", "BATF::JUN", "FOSL1", "FOSL2", "TCF7L2")
labels2 = c("LEF1", "EOMES", "TBX21", "FOS::JUN", "JUNB", "BATF::JUN", "FOSL1", "FOSL2", "TCF7L2", "TCF4", "RUNX3", "RUNX2", "RORA", "NFYB", "FOS","CTCF", "IRF8")

ggplot(var, aes(x=rank, y=variability)) + geom_point() +
  geom_text_repel(aes(label=ifelse(variability>8, as.character(name), "")), size=3) +
  #geom_text_repel(aes(label=ifelse(name%in%labels, as.character(name), "")), nudge_y = 0.1,nudge_x=2, min.segment.length = 0) +
  scale_x_continuous(limits=c(0,50)) +
  scale_y_continuous(limits=c(0,25)) + 
  theme_classic()
ggsave(here(output_dir, "Tmem_unstim_varplot2.pdf"), width=3.5, height=3.5)

set.seed(307)
tsne_results <- deviationsTsne(dev, threshold = 1.5, perplexity = 9)
dev$celltype <- gsub("_blood_unstim", "", dev$celltype)
tsne_plots <- plotDeviationsTsne(dev, tsne_results, 
                                 annotation_name = labels2, 
                                 sample_column = "celltype", 
                                 shiny = FALSE)
filenames = c("All", labels2)
for (i in 1:length(tsne_plots)){
  set.seed(307)
  pdf(file=here(output_dir, paste(filenames[i], ".pdf", sep="")), height=3.5, width=3.5)
  print(tsne_plots[[i]])
  dev.off()
}

custom_pal <- c("#90A4ADFF","#5F7D8BFF", "#2096F2FF", "#0C46A0FF", "#FFB74CFF", "#E55100FF", "#870D4EFF")
names(custom_pal) <- c("CD4_Nav", "CD8_Nav", "CD4_CM", "CD8_CM", "CD4_EM", "CD8_EM", "CD8_EMRA")


pdf(file=here(output_dir, paste("All", ".pdf", sep="")), height=3.5, width=4.5)
tsne_plots[[1]] + scale_color_manual(values=custom_pal)
dev.off()
