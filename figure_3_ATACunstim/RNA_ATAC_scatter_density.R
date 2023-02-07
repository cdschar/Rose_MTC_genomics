# Script for plotting RNA and ATAC scatter plot with genes
#jrose
#15Mar21

###############################################################################################
#Setup and data directories (change as needed)

library(here)
library(tidyverse)
library(scales)
library(ggrepel)
library(awtools)
library(ggrastr)

#ATAC_input_dir <- "analysis/unstim/diff"
ATAC_input_dir <- "/home/jim/Documents/Emory/Boss/data/ATAC/diff"
ATAC_file <- "diff.significant.glm.ATACseq_rev.txt"

#RNA_input_dir <- "/Volumes/ESB/Boss_Projects/Human_Tcell_subsets/RNA_combined/analysis/unstim/diff_unstim"
RNA_input_dir <- "/home/jim/Documents/Emory/Boss/data/RNA/naive_v_mem_overlap"
RNA_file <- "diff.significant.glm.Human_Tcell_RNA.rev.txt"

#output_dir <- "analysis/unstim/diff/concord"
output_dir <- "tests"

bistools_dir <- "/home/jim/Code/Boss/bistools/bisTools.R"

#fqDir = "/Volumes/ESB/Boss_Projects/Human_Tcell_subsets/RNA_combined/pipeline/";
fqDir = "/home/jim/Documents/Emory/Boss/data/"
fqFile = "RNAseq.sample.manifest.txt"

###############################################################################################
#Loading data:
ATAC_diff <- read.table(paste(ATAC_input_dir,ATAC_file, sep="/"), sep="\t", header=T, quote= "", comment.char="", check.names=F)
RNA_diff <- read.delim(paste(RNA_input_dir, RNA_file, sep="/"), header=T, sep="\t")

#read in the manifest file & subset as desirec

files = read.table(paste0(fqDir, fqFile), sep = "\t", header = T, as.is = T);
files = files[grepl("EM|CM|Nav", files$group),]
files = files[!grepl("CD4", files$group),]
files = files[which(files$include),]
files = files[files$stim=="unstim",]
files = mutate(files, cellSubtype = paste(cellType, cellSubset, stim, sep="_"))

###############################################################################################
#Data wrangling, filtering by group-wise expression

#Creating table of rpkm values from samples of interest
rpkm <- RNA_diff[,paste(files$sample, ".rpkm", sep="")]
rownames(rpkm) <- RNA_diff$ENTREZID

#Removing genes with low expression
source(bistools_dir);
cutoff=3
#^Change this cutoff based on how you want to filter DEG expression

expr_test = apply(rpkm[,grepl(paste(files$sample, collapse = "|"), colnames(rpkm))], 1, function(x) filter_detected(x, files$group, cutoff) )
print('DEGs before filter')
dim(rpkm)
rpkm <- rpkm[expr_test,]
print('DEGs after filter')
dim(rpkm)
RNA_diff <- RNA_diff[RNA_diff$ENTREZID%in%rownames(rpkm),]

###############################################################################################
#Assigning column numbers for groups of interest in RNA data

comps <- c("CD4_CM_blood_unstim.v.CD4_Nav_blood_unstim",
           "CD4_EM_blood_unstim.v.CD4_Nav_blood_unstim"#,
           #"CD8_CM_blood_unstim.v.CD8_Nav_blood_unstim",
           #"CD8_EM_blood_unstim.v.CD8_Nav_blood_unstim",
           #"CD8_EMRA_blood_unstim.v.CD8_Nav_blood_unstim"
           )
#^Change the entries in this vector to your comparisons of interest
names(comps) <- c("CD4_CM", "CD4_EM"#,
                  #"CD8_EMRA"
                  )
#^Change entries in this vector to the labels you want in your filenames

r_sig_cols <- vector()
a_sig_cols <- vector()

for (i in 1:length(comps)){
  r_sig_cols[i] <- grep(paste0(comps[i],".sig"), names(RNA_diff))
  a_sig_cols[i] <- grep(paste0(comps[i],".sig"), names(ATAC_diff))
}

r_ent_ID <- grep("ENTREZID", names(RNA_diff))
SYM <- grep("SYMBOL", names(RNA_diff))
a_ent_ID <- grep("Entrez.ID", names(ATAC_diff))

# r_sig_cols <- c(r_ent_ID, SYM, r_sig_cols)
# a_sig_cols <- c(a_ent_ID, a_sig_cols)

r_sig_vars <- names(RNA_diff)[r_sig_cols]  #Variable name for dplyr
a_sig_vars <- names(ATAC_diff)[a_sig_cols]  #Variable name for dplyr

r_fc_cols <- r_sig_cols - 3
r_fc_vars <- names(RNA_diff)[r_fc_cols]
a_fc_cols <- a_sig_cols - 3
a_fc_vars <- names(ATAC_diff)[a_fc_cols]


###############################################################################################
#For each memory subset select out relevant logFC and sig columns, join together RNA and ATAC datasets
#Here is am also totaling up the signifcance T/F columns and adding a variable called diff for plotting later

folds <- list()

for (i in 1:length(comps)){
  folds[[i]] <- select(RNA_diff, SYMBOL, ENTREZID, !!r_fc_vars[i], !!r_sig_vars[i]) %>%
    left_join(y=select(ATAC_diff, Entrez.ID, !!a_fc_vars[i], !!a_fc_vars[i], !!a_sig_vars[i]),
              by=c("ENTREZID" = "Entrez.ID"),
              suffix = c(".RNA", ".ATAC")) %>%
    drop_na(!!paste(a_sig_vars[i], ".ATAC", sep="")) %>%
    select(SYMBOL, !!paste(r_fc_vars[i], ".RNA", sep=""), !!paste(a_fc_vars[i], ".ATAC", sep=""), !!paste(r_sig_vars[i], ".RNA", sep=""), !!paste(a_sig_vars[i], ".ATAC", sep=""))
  names(folds[[i]]) <- c("SYMBOL", "RNA_logFC", "ATAC_logFC", "RNA_sig", "ATAC_sig")
  folds[[i]]$sigtotal <- select(folds[[i]], RNA_sig, ATAC_sig) %>% rowSums() %>% factor()
  folds[[i]] <- mutate(folds[[i]], diff=abs(RNA_logFC-ATAC_logFC))
}

#Function for labeling number of significant DEGs/DARs between RNA and ATAC
sig_class <- function(X){
  X$sig <- vector(mode="character", length=length(X$SYMBOL))
  
  for (i in 1:length(X$SYMBOL)){
    if (X$sigtotal[i]==1 & X[i,4]==TRUE){
      X$sig[i]=names(X)[4]
    } else if (X$sigtotal[i]==1 & X[i,5]==TRUE){
      X$sig[i]=names(X)[5]
    } else if (X$sigtotal[i]==2) {
      X$sig[i]="Both"
    } else if (X$sigtotal[i]==0) {
      X$sig[i]="None"
    } else {
      X$sig[i]=NA
    }
  }
  X$sig <- as.factor(X$sig)
  return(X)
}
folds <- map(folds, sig_class)
folds_2sig <- map(folds, function(x) subset(x, sigtotal==2))

###############################################################################################
#Plotting

#Color pallete
base_color_pal <- mpalette[c(4,5,8,3)]
names(base_color_pal) <- levels(folds[[1]]$sig)

# Concordance plots
theme_set(theme_gray()+ theme(axis.line = element_line(size=0.5),
                              panel.background = element_rect(fill=NA,size=rel(20)),
                              panel.grid.minor = element_line(colour = NA),
                              axis.text = element_text(size=16),
                              axis.title = element_text(size=18)))

#Difference between logFCs of states which are labeled
diff_lvl = 4
diff_lvl2 = 2.5

for (i in 1:length(comps)){
  ggplot(folds[[i]], aes(x=RNA_logFC, y=ATAC_logFC, label=SYMBOL)) + 
    rasterize(geom_point(color="grey", alpha= 0.4)) +
    scale_y_continuous(limits=c(-9, 9)) +
    scale_x_continuous(limits=c(-9, 9)) +
    geom_density_2d() +
    #geom_text(aes(label=ifelse(SYMBOL%in%genes & (abs(RNA_logFC)>=1 & abs(ATAC_logFC)>=1) , as.character(SYMBOL),"")), hjust=1, vjust=0, color="cornflower blue") +
    #guides(color=F) +
    geom_hline(yintercept=0) +
    geom_hline(yintercept = -1, linetype="dashed") +
    geom_hline(yintercept = 1, linetype="dashed") +
    geom_vline(xintercept=0) +
    geom_vline(xintercept=1, linetype="dashed") +
    geom_vline(xintercept = -1, linetype="dashed") + 
    annotate("text", x=4, y=5, 
             label=paste("r = ",
                         as.character(round(cor(folds[[i]]$RNA_logFC, folds[[i]]$ATAC_logFC), 2)),
                         sep=""),
             size=8
            )
  ggsave(paste(output_dir, paste0(names(comps)[i],"_RNA_ATAC_concord_labfilter_density_cor.pdf"), sep="/"),height=8, width=8)
  ggsave(paste(output_dir, paste0(names(comps)[i],"_RNA_ATAC_concord_labfilter_density_cor.jpg"), sep="/"), height=8, width=8)
}

###############################################################################################
#End of Script
print("########################################################################")
print("")
print("")

sessionInfo()
