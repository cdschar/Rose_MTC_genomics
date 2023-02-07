# R function for subsetting DEG file by comparison and finding shared or unique DEGs
# of naive vs meory T cell subsets
#jrose
#12Feb20
#Edited: 23Sep20

library(tidyverse)
library(openxlsx)
library(here)

# Loading differentially expressed genes
diff <- read.delim(here("/analysis/unstim/diff_unstim/diff.significant.glm.Human_Tcell_RNA.NAremoved.txt"), header=T, sep="\t")

CM_CD4only <- subset(diff, CD4_Nav_blood_unstim.v.CD4_CM_blood_unstim.sig =="TRUE" &
                       CD8_Nav_blood_unstim.v.CD8_CM_blood_unstim.sig =="FALSE") %>%
  mutate(overlap="CD4_CM")
CM_CD8only <- subset(diff, CD4_Nav_blood_unstim.v.CD4_CM_blood_unstim.sig =="FALSE" &
                       CD8_Nav_blood_unstim.v.CD8_CM_blood_unstim.sig =="TRUE") %>%
  mutate(overlap="CD8_CM")
CM_shared <- subset(diff, CD4_Nav_blood_unstim.v.CD4_CM_blood_unstim.sig =="TRUE" &
                      CD8_Nav_blood_unstim.v.CD8_CM_blood_unstim.sig =="TRUE") %>%
  mutate(overlap="CD4CD8_CM")

EM_CD4only <- subset(diff, CD4_Nav_blood_unstim.v.CD4_EM_blood_unstim.sig =="TRUE" &
                       CD8_Nav_blood_unstim.v.CD8_EM_blood_unstim.sig =="FALSE") %>%
  mutate(overlap="CD4_EM")
EM_CD8only <- subset(diff, CD4_Nav_blood_unstim.v.CD4_EM_blood_unstim.sig =="FALSE" &
                       CD8_Nav_blood_unstim.v.CD8_EM_blood_unstim.sig =="TRUE") %>%
  mutate(overlap="CD8_EM")
EM_shared <- subset(diff, CD4_Nav_blood_unstim.v.CD4_EM_blood_unstim.sig =="TRUE" &
                      CD8_Nav_blood_unstim.v.CD8_EM_blood_unstim.sig =="TRUE") %>%
  mutate(overlap="CD4CD8_EM")


#Make lists
CM_naive_Mem_DEGs <- list(CM_CD4only <- CM_CD4only,
                          CM_CD8only <- CM_CD8only,
                          CM_shared <- CM_shared)
names(CM_naive_Mem_DEGs) <- c("CM_CD4only", "CM_CD8only", "CM_shared")

EM_naive_Mem_DEGs <- list(EM_CD4only <- EM_CD4only,
                          EM_CD8only <- EM_CD8only,
                          EM_shared <- EM_shared)
names(EM_naive_Mem_DEGs) <- c("EM_CD4only", "EM_CD8only", "EM_shared")

#Save
save(CM_naive_Mem_DEGs,EM_naive_Mem_DEGs, file=here("/analysis/unstim/diff_unstim/", "EMCM_unique_shared_DEGs_Naive_v_Mem.rda"))

