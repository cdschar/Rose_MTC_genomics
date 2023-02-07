#Trims the output of ATAC_peak_cat.R
#v1.1 includes adding Memory category
#18Dec20
#jrose

library(tidyverse)
library(readr)
library(here)

input_dir <- "analysis/Stim_Unstim/diff"
output_dir <- "analysis/Stim_Unstim/diff"

diff <- read_tsv(here(input_dir, "diff.glm.stim_unstim_cat.ATACseq_v1_1.txt"), col_names=TRUE)

#Need to load manifest file for this to work
fqDir = "/Volumes/ESB/Boss_Projects/Human_Tcell_subsets/ATAC_combined/pipeline/";
fqFile = "ATACseq.sample.manifest.withFRiP.txt";
files = read.table(paste0(fqDir, fqFile), sep = "\t", header = T, as.is = T);
files = files[grepl("EM|CM|Nav", files$group),]
#files = files[!grepl("CD4", files$group),]

remove <- paste(files$sample[!files$include], ".rppm", sep="")
files = files[which(files$include),]

#files = files[files$stim=="unstim",]
files = mutate(files, cellSubtype = paste(cellType, cellSubset, stim, sep="_"))
grps_consv <- unique(files$group)

#Induced only looks at stim samples
files_induc <- files[files$stim=="stim",]
grps_induc <- unique(files_induc$group)

#Naive only looks at Naive unstim
files_nav <- files[grepl("Nav_unstim", files$cellSubtype),]
grps_nav <- unique(files_nav$group)

#Primed only looks at mem unstim (also used for Memory group)
files_prime <- files[grepl("((CM)|(EM)|(EMRA))_unstim", files$cellSubtype),]
grps_prime <- unique(files_prime$group)

#Remove samples not with include in manifest from diff file
diff <- diff[,!colnames(diff)%in%remove]

#############################################################################################
#Break up the diff dataframe into lists of each peak type with separate entries for each memory sample

consv <- list()
induc <- list()
naive <- list()
primed <- list()

vars <- names(diff)[grep("peak_cat", names(diff))]

for (i in 1:length(vars)){
  consv[[i]] <- diff[diff[[vars[i]]]=="Conserved" | diff[[vars[i]]]=="Other",]
  #^NOTE: Consv also contains other group peaks!
  induc[[i]] <- diff[diff[[vars[i]]]=="Induced",]
  naive[[i]] <- diff[diff[[vars[i]]]=="Naive",]
  primed[[i]] <- diff[diff[[vars[i]]]=="Primed",]
}

names(consv) <- gsub("_peak_cat", "", vars)
names(induc) <- gsub("_peak_cat", "", vars)
names(naive) <- gsub("_peak_cat", "", vars)
names(primed) <- gsub("_peak_cat", "", vars)

print("Before trimming")
print("#############################################################################################")
str(consv, max.level = 1)
str(induc, max.level = 1)
str(naive, max.level = 1)
str(primed, max.level = 1)
print("#############################################################################################")
rm(diff)
gc()

#############################################################################################

topcut <- 3
trimmer <- function(data, mani, cut) {
  #filter_detected funtion to filter for detected genes based on a cutoff and groups from bistools
  filter_detected = function(x, group, cutoff){
    
    one_gene = data.frame(group=group, rpm=as.vector(t(x)))
    one_gene = aggregate(rpm~group, data=one_gene, FUN=function(x) sum(x>cutoff))
    
    if( any( one_gene$rpm >= table(group) ) )
    { return(TRUE) }
    else
    { return(FALSE) }
  }
  
  if (dim(data)[1] > 0){
    access_test = apply(data[,grepl(paste(paste0(mani$sample, ".rppm"), collapse = "|"), colnames(data))], 1, function(x) filter_detected(x, mani$group, cut) )
    #access_test <- apply(data,c(1), function(x)any(which(x>=3)))
    data <- data[access_test,]
  } else {
    print("Empty dataframe")
  }
  
  return(data)
}

consv_trim <- map(consv, function(x) trimmer(x, files, topcut))
print("Consv finished")
timestamp()
induc_trim <- map(consv, function(x) trimmer(x, files_induc, topcut))
print("Induc finished")
timestamp()
naive_trim <- map(consv, function(x) trimmer(x, files_nav, topcut))
print("Naive finished")
timestamp()
primed_trim <- map(consv, function(x) trimmer(x, files_prime, topcut))
print("Primed finished")
timestamp()
mem_trim <- map(mem, function(x) trimmer(x, files_prime, topcut))
print("Memory finished")


print("After trimming")
print("#############################################################################################")
str(consv_trim, max.level = 1)
str(induc_trim, max.level = 1)
str(naive_trim, max.level = 1)
str(primed_trim, max.level = 1)
print("#############################################################################################")

#############################################################################################
#Reconstruct dataframes of all peak types and save

#vars2 <- gsub("_peak_cat", "", vars)

diff_trim_CD4_CM <- bind_rows(consv_trim[["CD4_CM"]], induc_trim[["CD4_CM"]], naive_trim[["CD4_CM"]], primed_trim[["CD4_CM"]])
diff_trim_CD4_EM <- bind_rows(consv_trim[["CD4_EM"]], induc_trim[["CD4_EM"]], naive_trim[["CD4_EM"]], primed_trim[["CD4_EM"]])
diff_trim_CD4_EMRA <- bind_rows(consv_trim[["CD4_EMRA"]], induc_trim[["CD4_EMRA"]], naive_trim[["CD4_EMRA"]], primed_trim[["CD4_EMRA"]])
diff_trim_CD8_CM <- bind_rows(consv_trim[["CD8_CM"]], induc_trim[["CD8_CM"]], naive_trim[["CD8_CM"]], primed_trim[["CD8_CM"]])
diff_trim_CD8_EM <- bind_rows(consv_trim[["CD8_EM"]], induc_trim[["CD8_EM"]], naive_trim[["CD8_EM"]], primed_trim[["CD8_EM"]])
diff_trim_CD8_EMRA <- bind_rows(consv_trim[["CD8_EMRA"]], induc_trim[["CD8_EMRA"]], naive_trim[["CD8_EMRA"]], primed_trim[["CD8_EMRA"]])

write.table(diff_trim_CD4_CM,file=here(output_dir, "diff.glm.CD4_CM_cat_trim.ATACseq.txt"), sep="\t", quote=F, row.names=F)
write.table(diff_trim_CD4_CM,file=here(output_dir, "diff.glm.CD4_EM_cat_trim.ATACseq.txt"), sep="\t", quote=F, row.names=F)
write.table(diff_trim_CD4_EMRA,file=here(output_dir, "diff.glm.CD4_EMRA_cat_trim.ATACseq.txt"), sep="\t", quote=F, row.names=F)
write.table(diff_trim_CD8_CM,file=here(output_dir, "diff.glm.CD8_CM_cat_trim.ATACseq.txt"), sep="\t", quote=F, row.names=F)
write.table(diff_trim_CD8_EM,file=here(output_dir, "diff.glm.CD8_EM_cat_trim.ATACseq.txt"), sep="\t", quote=F, row.names=F)
write.table(diff_trim_CD8_EMRA,file=here(output_dir, "diff.glm.CD8_EMRA_cat_trim.ATACseq.txt"), sep="\t", quote=F, row.names=F)

sessionInfo()