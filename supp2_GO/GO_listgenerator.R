# Geneset analysis of shared EM or CM between CD4 and CD8 respectively

library(tidyverse)
library(here)

load("~/Documents/Emory/Boss/data/RNA/naive_v_mem_overlap/EMCM_unique_shared_DEGs_Naive_v_Mem.rda")
output_dir <- "naive_v_mem_overlap/CD4CD8_comp"

for (i in 1:length(names(CM_naive_Mem_DEGs))){
  outgenelst <- CM_naive_Mem_DEGs[[i]]$SYMBOL
  write.csv(outgenelst, file=here(output_dir, paste(names(CM_naive_Mem_DEGs)[i], "genelist.csv", sep="_")))
}

for (i in 1:length(names(EM_naive_Mem_DEGs))){
  outgenelst <- EM_naive_Mem_DEGs[[i]]$SYMBOL
  write.csv(outgenelst, file=here(output_dir, paste(names(EM_naive_Mem_DEGs)[i], "genelist.csv", sep="_")))
}

#Panther GO analysis run via browswer using the txt files exported

table(CM_naive_Mem_DEGs$CM_shared$SYMBOL%in%EM_naive_Mem_DEGs$EM_shared$SYMBOL)
