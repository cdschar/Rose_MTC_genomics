setwd("~/Documents/Emory/Boss/data/RNA")

library(tidyverse)
library(here)
library(ggpubr)

#Setting directories
# cts_input_dir <- "coverage"
man_input_dir <- "/home/jim/Documents/Emory/Boss/data"
diff_input_dir <- "StimUnstim"
output_dir <- "barplots/fig5"

#Loading data and sample manifest
RNA_diff <- read.delim(here(diff_input_dir, "diff.significant.glm.Human_Tcell_RNA.StimUnstimRev.txt"), header=T, sep="\t")

RNA_diff %>% select(SYMBOL,CD4_Nav_blood_stim.v.CD4_Nav_blood_unstim.sig,CD4_CM_blood_stim.v.CD4_CM_blood_unstim.sig,CD4_EM_blood_stim.v.CD4_EM_blood_unstim.sig,CD8_Nav_blood_stim.v.CD8_Nav_blood_unstim.sig,CD8_CM_blood_stim.v.CD8_CM_blood_unstim.sig,CD8_EM_blood_stim.v.CD8_EM_blood_unstim.sig,CD8_EMRA_blood_stim.v.CD8_EMRA_blood_unstim.sig) %>%
  subset(SYMBOL=="AHR"|SYMBOL=="HIF1A"|SYMBOL=="EOMES"|SYMBOL=="TBX21"|SYMBOL=="MSC")
