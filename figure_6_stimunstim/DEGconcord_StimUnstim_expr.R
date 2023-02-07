#Concordance plots of Naive vs Memory Subsets
#13Mar20
#jrose

setwd("~/Documents/Emory/Boss/data/RNA")

library(tidyverse)
library(here)
library(awtools)
library(ggrepel)

#Loading differentially expressed genes

# input_dir <- "analysis/unstim/diff_unstim"
# output_dir <- "analysis/unstim/diff_unstim/StimUnstim"
input_dir <- "StimUnstim"
output_dir <- "StimUnstim"

diff <- read.delim(here(input_dir, "diff.significant.glm.Human_Tcell_RNA.StimUnstimRev.txt"), header=T, sep="\t")

#Assigning column numbers for groups of interest
SYM <- 2
col1 <- 209  #CD4_CM_unstim sig column    #Col # for base R
col2 <- 134  #CD4_CM_stim sig column

col3 <- 439  #CD4_EM_unstim sig  column
col4 <- 384  #CD4_EM_stim sig column

col5 <- 629 #CD8_CM_unstim sig column
col6 <- 594 #CD8_CM_stim sig column

col7 <- 699 #CD8_EM_unstim sig column
col8 <- 684 #CD8_EM_stim sig column

col9 <- 674 #CD8_EMRA_unstim sig column
col10 <- 649 #CD8_EMRA_stim sig column

sig_cols <- c(SYM, col1, col2, col3, col4, col5, col6, col7, col8, col9, col10)
sig_vars <- names(diff)[sig_cols]  #Variable name for dplyr

fc_cols <- c(2, sig_cols[-1] - 3)
fc_vars <- names(diff)[fc_cols]

fdr_cols <- c(2, sig_cols[-1] - 1)
fdr_vars <- names(diff)[fdr_cols]

#rpkm gene expression columns
col_CD4CM <- grep("CD4_CM.rpkm", names(diff))
names(col_CD4CM) <- names(diff)[col_CD4CM]
CD4CM_mean_rpkm <- rowMeans(diff[,col_CD4CM])
col_stim_CD4CM <- grep("CD4_CMStim.rpkm", names(diff))
names(col_stim_CD4CM) <- names(diff)[col_stim_CD4CM]
CD4CM_meanstim_rpkm <- rowMeans(diff[,col_stim_CD4CM])

col_CD4EM <- grep("CD4_EM.rpkm", names(diff))
names(col_CD4EM) <- names(diff)[col_CD4EM]
CD4EM_mean_rpkm <- rowMeans(diff[,col_CD4EM])
col_stim_CD4EM <- grep("CD4_EMStim.rpkm", names(diff))
names(col_stim_CD4EM) <- names(diff)[col_stim_CD4EM]
CD4EM_meanstim_rpkm <- rowMeans(diff[,col_stim_CD4EM])

col_CD8CM <- grep("CD8_CM.rpkm", names(diff))
names(col_CD8CM) <- names(diff)[col_CD8CM]
CD8CM_mean_rpkm <- rowMeans(diff[,col_CD8CM])
col_stim_CD8CM <- grep("CD8_CMStim.rpkm", names(diff))
names(col_stim_CD8CM) <- names(diff)[col_stim_CD8CM]
CD8CM_meanstim_rpkm <- rowMeans(diff[,col_stim_CD8CM])

col_CD8EM <- grep("CD8_EM.rpkm", names(diff))
names(col_CD8EM) <- names(diff)[col_CD8EM]
CD8EM_mean_rpkm <- rowMeans(diff[,col_CD8EM])
col_stim_CD8EM <- grep("CD8_EMStim.rpkm", names(diff))
names(col_stim_CD8EM) <- names(diff)[col_stim_CD8EM]
CD8EM_meanstim_rpkm <- rowMeans(diff[,col_stim_CD8EM])

col_CD8EMRA <- grep("CD8_EMRA.rpkm", names(diff))
names(col_CD8EMRA) <- names(diff)[col_CD8EMRA]
CD8EMRA_mean_rpkm <- rowMeans(diff[,col_CD8EMRA])
col_stim_CD8EMRA <- grep("CD8_EMRAStim.rpkm", names(diff))
names(col_stim_CD8EMRA) <- names(diff)[col_stim_CD8EMRA]
CD8EMRA_meanstim_rpkm <- rowMeans(diff[,col_stim_CD8EMRA])

#Subsetting just columns for specific comparisons
expr_cutoff = 3

CD4_CM <- select(diff,SYMBOL, !!fc_vars[2], !!fc_vars[3], !!sig_vars[2], !!sig_vars[3])
names(CD4_CM) <- c("SYMBOL", "Unstim_logFC", "Stim_logFC", "Unstim_sig", "Stim_sig")
CD4_CM$sigtotal <- select(CD4_CM, Unstim_sig, Stim_sig) %>% rowSums() %>% factor()
CD4_CM <- mutate(CD4_CM, diff=abs(Unstim_logFC-Stim_logFC), meanCM=CD4CM_mean_rpkm, meanCMstim=CD4CM_meanstim_rpkm) %>%
  subset(meanCM >=expr_cutoff | meanCMstim >=expr_cutoff)

CD4_EM<- select(diff,SYMBOL, !!fc_vars[4], !!fc_vars[5], !!sig_vars[4], !!sig_vars[5])
names(CD4_EM) <- c("SYMBOL", "Unstim_logFC", "Stim_logFC", "Unstim_sig", "Stim_sig")
CD4_EM$sigtotal <- select(CD4_EM, Unstim_sig, Stim_sig) %>% rowSums() %>% factor()
CD4_EM<- mutate(CD4_EM, diff=abs(Unstim_logFC-Stim_logFC), meanEM=CD4EM_mean_rpkm, meanEMstim=CD4EM_meanstim_rpkm) %>%
  subset(meanEM >=expr_cutoff | meanEMstim >=expr_cutoff)

CD8_CM <- select(diff,SYMBOL, !!fc_vars[6], !!fc_vars[7], !!sig_vars[6], !!sig_vars[7])
names(CD8_CM) <- c("SYMBOL", "Unstim_logFC", "Stim_logFC", "Unstim_sig", "Stim_sig")
CD8_CM$sigtotal <- select(CD8_CM, Unstim_sig, Stim_sig) %>% rowSums() %>% factor()
CD8_CM <- mutate(CD8_CM, diff=abs(Unstim_logFC-Stim_logFC), meanCM=CD8CM_mean_rpkm, meanCMstim=CD8CM_meanstim_rpkm) %>%
  subset(meanCM >=expr_cutoff | meanCMstim >=expr_cutoff)

CD8_EM <- select(diff,SYMBOL, !!fc_vars[8], !!fc_vars[9], !!sig_vars[8], !!sig_vars[9] )
names(CD8_EM) <- c("SYMBOL", "Unstim_logFC", "Stim_logFC", "Unstim_sig", "Stim_sig")
CD8_EM$sigtotal <- select(CD8_EM, Unstim_sig, Stim_sig) %>% rowSums() %>% factor()
CD8_EM <- mutate(CD8_EM, diff=abs(Unstim_logFC-Stim_logFC), meanEM=CD8EM_mean_rpkm, meanEMstim=CD8EM_meanstim_rpkm) %>%
  subset(meanEM >=expr_cutoff | meanEMstim >=expr_cutoff)

CD8_EMRA <- select(diff,SYMBOL, !!fc_vars[10], !!fc_vars[11], !!sig_vars[10], !!sig_vars[11])
names(CD8_EMRA) <- c("SYMBOL", "Unstim_logFC", "Stim_logFC", "Unstim_sig", "Stim_sig")
CD8_EMRA$sigtotal <- select(CD8_EMRA, Unstim_sig, Stim_sig) %>% rowSums() %>% factor()
CD8_EMRA <- mutate(CD8_EMRA, diff=abs(Unstim_logFC-Stim_logFC), meanEMRA=CD8EMRA_mean_rpkm, meanEMRAstim=CD8EMRA_meanstim_rpkm) %>%
  subset(meanEMRA >=expr_cutoff | meanEMRAstim >=expr_cutoff)

#Concatinating into list object
folds <- list(CD4_CM=CD4_CM,
              CD4_EM=CD4_EM,
              CD8_EM=CD8_EM,
              CD8_CM=CD8_CM,
              CD8_EMRA=CD8_EMRA)

#Classifying by significance in each comparison in the sig variable

##Function definition
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

#Mapping through folds list
folds <- map(folds, sig_class)

#Saving folds list
save(folds, file=here(output_dir, "foldchanges.rda"))

#Color pallete
base_color_pal <- mpalette[c(5,8,1,9)]
names(base_color_pal) <- levels(folds$CD4_CM$sig)

# Concordance plots

diff_cut <- 10.5

theme_set(theme_gray()+ theme(axis.line = element_line(size=0.5),
                              panel.background = element_rect(fill=NA,size=rel(20)),
                              panel.grid.minor = element_line(colour = NA),
                              axis.text = element_text(size=10), axis.title = element_text(size=12)))


p_CD4_CM <- ggplot(subset(folds$CD4_CM, sig !="None"), aes(x=Unstim_logFC, y=Stim_logFC, label=SYMBOL))
p_CD4_CM + geom_point(aes(color=sig), size=1) +
  geom_abline(slope=1, intercept=0, ) +
  geom_abline(slope=1, intercept=-2, linetype="dashed") +
  geom_abline(slope=1, intercept=2, linetype="dashed") +
  scale_y_continuous(limits=c(-9, 9)) +
  scale_x_continuous(limits=c(-9, 9)) +
  scale_color_manual(name="FDR Sig", values=base_color_pal) +
  geom_text(aes(label=ifelse(diff>diff_cut & sigtotal==1, as.character(SYMBOL),"")), hjust=1, vjust=0)
ggsave(here(output_dir, "CD4_CM_StimUnstim_concord.pdf"), width=4.5, height=3.5)


p_CD4_EM <- ggplot(subset(folds$CD4_EM, sig !="None"), aes(x=Unstim_logFC, y=Stim_logFC, label=SYMBOL))
p_CD4_EM + geom_point(aes(color=sig), size=1) +
  geom_abline(slope=1, intercept=0, ) +
  geom_abline(slope=1, intercept=-2, linetype="dashed") +
  geom_abline(slope=1, intercept=2, linetype="dashed") +
  scale_y_continuous(limits=c(-9, 9)) +
  scale_x_continuous(limits=c(-9, 9)) +
  scale_color_manual(name="FDR Sig", values=base_color_pal) +
  geom_text(aes(label=ifelse(diff>diff_cut & sigtotal==1, as.character(SYMBOL),"")), hjust=1, vjust=0)
ggsave(here(output_dir, "CD4_EM_StimUnStim_concord.pdf"), width=4.5, height=3.5)


p_CD8_CM <- ggplot(subset(folds$CD8_CM, sig !="None"), aes(x=Unstim_logFC, y=Stim_logFC, label=SYMBOL))
p_CD8_CM + geom_point(aes(color=sig), size=1) +
  geom_abline(slope=1, intercept=0, ) +
  geom_abline(slope=1, intercept=-2, linetype="dashed") +
  geom_abline(slope=1, intercept=2, linetype="dashed") +
  scale_y_continuous(limits=c(-9, 9)) +
  scale_x_continuous(limits=c(-9, 9)) +
  scale_color_manual(name="FDR Sig", values=base_color_pal) +
  geom_text(aes(label=ifelse(diff>diff_cut & sigtotal==1, as.character(SYMBOL),"")), hjust=1, vjust=0)
ggsave(here(output_dir, "CD8_CM_StimUnStim_concord.pdf"), width=4.5, height=3.5)


p_CD8_EM <- ggplot(subset(folds$CD8_EM, sig !="None"), aes(x=Unstim_logFC, y=Stim_logFC, label=SYMBOL))
p_CD8_EM + geom_point(aes(color=sig), size=1) +
  geom_abline(slope=1, intercept=0, ) +
  geom_abline(slope=1, intercept=-2, linetype="dashed") +
  geom_abline(slope=1, intercept=2, linetype="dashed") +
  scale_y_continuous(limits=c(-9, 9)) +
  scale_x_continuous(limits=c(-9, 9)) +
  scale_color_manual(name="FDR Sig", values=base_color_pal) +
  geom_text(aes(label=ifelse(diff>diff_cut & sigtotal==1, as.character(SYMBOL),"")), hjust=1, vjust=0)
ggsave(here(output_dir, "CD8_EM_StimUnstim_concord.pdf"), width=4.5, height=3.5)

p_CD8_EMRA <- ggplot(subset(folds$CD8_EMRA, sig !="None"), aes(x=Unstim_logFC, y=Stim_logFC, label=SYMBOL))
p_CD8_EMRA + geom_point(aes(color=sig), size=1) +
  geom_abline(slope=1, intercept=0, ) +
  geom_abline(slope=1, intercept=-2, linetype="dashed") +
  geom_abline(slope=1, intercept=2, linetype="dashed") +
  scale_y_continuous(limits=c(-9, 9)) +
  scale_x_continuous(limits=c(-9, 9)) +
  scale_color_manual(name="FDR Sig", values=base_color_pal) +
  geom_text(aes(label=ifelse(diff>diff_cut & sigtotal==1, as.character(SYMBOL),"")), hjust=1, vjust=0)
ggsave(here(output_dir, "CD8_EMRA_StimUnstim_concord.pdf"), width=4.5, height=3.5)


#Tally up the poised and resting DEGs
sig_full <- function(x){
  mutate(x, 
         sig_full = case_when(
              sig == "Unstim_sig" & Unstim_logFC > 0 ~ "Resting_up",
              sig == "Unstim_sig" & Unstim_logFC < 0 ~ "Resting_dn",
              sig == "Stim_sig" & Stim_logFC > 0 ~ "Poised_up",
              sig == "Stim_sig" & Stim_logFC < 0 ~ "Poised_dn",
              sig == "Both" & Unstim_logFC > 0 & Stim_logFC >0 ~ "Both_up_up",
              sig == "Both" & Unstim_logFC < 0 & Stim_logFC >0 ~ "Both_dn_up",
              sig == "Both" & Unstim_logFC < 0 & Stim_logFC <0 ~ "Both_dn_dn",
              sig == "Both" & Unstim_logFC > 0 & Stim_logFC <0 ~ "Both_up_dn"
                                )
  )
}

folds <- map(folds, sig_full)

sig_counts <- bind_rows(table(folds$CD4_CM$sig_full), 
      table(folds$CD4_EM$sig_full),
      table(folds$CD8_CM$sig_full),
      table(folds$CD8_EM$sig_full),
      table(folds$CD8_EMRA$sig_full)
      ) %>% data.frame()
sig_counts$Poised_dn <- -sig_counts$Poised_dn
sig_counts$Resting_dn <- -sig_counts$Resting_dn
sig_counts$Both_dn_dn <- -sig_counts$Both_dn_dn
sig_counts$subset <- c("CD4_CM", "CD4_EM", "CD8_CM", "CD8_EM", "CD8_EMRA")

sig_counts_long <- pivot_longer(sig_counts, !subset, names_to="DEG_state_dir", values_to="count") %>% 
  replace_na(list(count=0)) %>%
  mutate(DEG_state = case_when(
         DEG_state_dir == "Poised_dn" | DEG_state_dir == "Poised_up" ~ "Poised",
         DEG_state_dir == "Resting_dn" | DEG_state_dir == "Resting_up" ~ "Resting",
         DEG_state_dir == "Both_dn_dn" | DEG_state_dir =="Both_up_up" | DEG_state_dir =="Both_dn_up" | DEG_state_dir =="Both_up_dn" ~ "Both"
          ) 
  ) #%>%
#replace_na(list(DEG_state="Both"))
sig_counts_long$count <- as.numeric(sig_counts_long$count)
sig_counts_long$DEG_state <- factor(sig_counts_long$DEG_state, levels=c("Resting", "Poised", "Both"), ordered=T)

# ggplot(subset(sig_counts_long, DEG_state != "NA"), aes(x=subset, y=count)) + geom_col(aes(fill=DEG_state_dir)) +
#   facet_wrap(vars(DEG_state)) +
#   theme(axis.text.x=element_text(angle=90)) +
#   scale_fill_manual(values=c("Resting_dn"="cornflower blue", "Resting_up"="firebrick", "Poised_dn"="light blue", "Poised_up"="orange"))
# ggsave(here(output_dir, "resting_poised_barplot.pdf"), width=4.5, height=3.5)
# 
# ggplot(subset(sig_counts_long, DEG_state=="Resting"), aes(x=subset, y=count)) + geom_col(aes(fill=DEG_state_dir)) +
#   theme(axis.text.x=element_text(angle=90)) +
#   scale_fill_manual(values=c("Resting_dn"="cornflower blue", "Resting_up"="firebrick", "Poised_dn"="light blue", "Poised_up"="orange"))
# ggsave(here(output_dir, "resting_barplot.pdf"), width=4.5, height=3.5)
# 
# ggplot(subset(sig_counts_long, DEG_state=="Poised"), aes(x=subset, y=count)) + geom_col(aes(fill=DEG_state_dir)) +
#   theme(axis.text.x=element_text(angle=90)) +
#   scale_fill_manual(values=c("Resting_dn"="cornflower blue", "Resting_up"="firebrick", "Poised_dn"="light blue", "Poised_up"="orange"))
# ggsave(here(output_dir, "poised_barplot.pdf"), width=4.5, height=3.5)

# ggplot(sig_counts_long, aes(x=subset, y=count)) + geom_col(aes(fill=DEG_state_dir)) +
#   theme(axis.text.x=element_text(angle=90)) +
#   scale_fill_manual(values=c("Both" = "purple", "Resting_dn"="cornflower blue", "Resting_up"="firebrick", "Poised_dn"="light blue", "Poised_up"="orange"))
# ggsave(here(output_dir, "stacked_barplot.pdf"), width=4.5, height=3.5)

ggplot(subset(sig_counts_long, DEG_state !="Both"), aes(x=subset, y=count)) + geom_col(aes(fill=DEG_state_dir)) +
  theme(axis.text.x=element_text(angle=90)) +
  scale_fill_manual(values=c("Both" = "purple", "Resting_dn"="cornflower blue", "Resting_up"="firebrick", "Poised_dn"="light blue", "Poised_up"="orange"))
ggsave(here(output_dir, "stacked_noboth_barplot.pdf"), width=4.5, height=3.5)

ggplot(subset(sig_counts_long, DEG_state == "Both"), aes(x=subset, y=count)) + geom_col(aes(fill=DEG_state_dir)) +
  theme(axis.text.x=element_text(angle=90))
ggsave(here(output_dir, "stacked_justboth_barplot.pdf"), width=4.5, height=3.5)

ggplot(subset(sig_counts_long, DEG_state_dir !="Both_up_dn" & DEG_state_dir !="Both_dn_up"), aes(x=subset, y=count)) + geom_col(aes(fill=DEG_state_dir)) +
  theme(axis.text.x=element_text(angle=90)) +
  scale_fill_manual(values=c("Both_up_up" = "purple", "Both_dn_dn"="violet","Resting_dn"="cornflower blue", "Resting_up"="firebrick", "Poised_dn"="light blue", "Poised_up"="orange"))
ggsave(here(output_dir, "stacked_both_barplot.pdf"), width=4.5, height=3.5)

ggplot(subset(sig_counts_long, DEG_state !="Resting" & DEG_state_dir !="Both_up_dn" & DEG_state_dir !="Both_dn_up"), aes(x=subset, y=count)) + geom_col(aes(fill=DEG_state_dir)) +
  theme(axis.text.x=element_text(angle=90)) +
  scale_fill_manual(values=c("Both_up_up" = "purple", "Both_dn_dn"="violet","Resting_dn"="cornflower blue", "Resting_up"="firebrick", "Poised_dn"="light blue", "Poised_up"="orange"))
ggsave(here(output_dir, "stacked_norest_barplot.pdf"), width=4.5, height=3.5)

#Volcano of poised, resting, both

small_diff <- select(diff, SYMBOL, ENTREZID, !!fc_vars[-1], !!fdr_vars[-1], !!sig_vars[-1])

join_diff <- function(x){
  left_join(x, small_diff, by="SYMBOL") %>%
    mutate(CD4_CM_blood_unstim.v.CD4_Nav_blood_unstim.neglog10p=-log(CD4_CM_blood_unstim.v.CD4_Nav_blood_unstim.fdr, base=10),
           CD4_CM_blood_stim.v.CD4_Nav_blood_stim.neglog10p=-log(CD4_CM_blood_stim.v.CD4_Nav_blood_stim.fdr, base=10),
           CD4_EM_blood_unstim.v.CD4_Nav_blood_unstim.neglog10p=-log(CD4_EM_blood_unstim.v.CD4_Nav_blood_unstim.fdr, base=10),
           CD4_EM_blood_stim.v.CD4_Nav_blood_stim.neglog10p=-log(CD4_EM_blood_stim.v.CD4_Nav_blood_stim.fdr, base=10),
           CD8_CM_blood_unstim.v.CD8_Nav_blood_unstim.neglog10p=-log(CD8_CM_blood_unstim.v.CD8_Nav_blood_unstim.fdr, base=10),
           CD8_CM_blood_stim.v.CD8_Nav_blood_stim.neglog10p=-log(CD8_CM_blood_stim.v.CD8_Nav_blood_stim.fdr, base=10),
           CD8_EM_blood_unstim.v.CD8_Nav_blood_unstim.neglog10p=-log(CD8_EM_blood_unstim.v.CD8_Nav_blood_unstim.fdr, base=10),
           CD8_EM_blood_stim.v.CD8_Nav_blood_stim.neglog10p=-log(CD8_EM_blood_stim.v.CD8_Nav_blood_stim.fdr, base=10),
           CD8_EMRA_blood_unstim.v.CD8_Nav_blood_unstim.neglog10p=-log(CD8_EMRA_blood_unstim.v.CD8_Nav_blood_unstim.fdr, base=10),
           CD8_EMRA_blood_stim.v.CD8_Nav_blood_stim.neglog10p=-log(CD8_EMRA_blood_stim.v.CD8_Nav_blood_stim.fdr, base=10)
           ) %>%
    replace_na(list(sig_full="None"))
}
full_folds <- map(folds, join_diff)

unstim_comps <- c("CD4_CM_blood_unstim.v.CD4_Nav_blood_unstim", 
                  "CD4_EM_blood_unstim.v.CD4_Nav_blood_unstim",
                  "CD8_CM_blood_unstim.v.CD8_Nav_blood_unstim",
                  "CD8_EM_blood_unstim.v.CD8_Nav_blood_unstim",
                  "CD8_EMRA_blood_unstim.v.CD8_Nav_blood_unstim")
stim_comps <- gsub("unstim", "stim", unstim_comps)

#Volcano plots of poised genes
for (i in 1:length(stim_comps)){
  short_comp <- paste(strsplit(unstim_comps[i], "_")[[1]][1], strsplit(unstim_comps[i], "_")[[1]][2], sep="_")
  varlogFC <- paste(stim_comps[i],"logFC", sep=".")
  varneglog10p <- paste(stim_comps[i],"neglog10p", sep=".")
  
  tmp <- subset(full_folds[[short_comp]], 
                sig_full=="Poised_up"| sig_full =="Poised_dn" | sig_full=="None") %>%
    rename(logFC=matches(varlogFC), neglog10p=matches(varneglog10p))
  
  print(
    ggplot(tmp, aes(x=logFC, y=neglog10p)) +
      geom_point(aes(color=sig_full)) +
      geom_hline(yintercept = 1.3, linetype="dashed") +
      geom_vline(xintercept = 1, linetype="dashed") +
      geom_vline(xintercept = -1, linetype="dashed") +
      scale_color_manual(values=c("None"="grey", "Poised_dn"="light blue", "Poised_up"="orange")) +
      annotate("text", x = -5, y = 25, label = as.character(abs(sig_counts[grep(short_comp, sig_counts$subset),"Poised_dn"])), color="light blue") +
      annotate("text", x = 5, y = 25, label = as.character(sig_counts[grep(short_comp, sig_counts$subset),"Poised_up"]), color="orange") #+
      #geom_text_repel(aes(label=ifelse(logFC >4 | neglog10p > 10, as.character(SYMBOL), "" )),size=2)
  )
  ggsave(filename=here(output_dir, "volcano", paste0(stim_comps[i], "_volcano.pdf")), height=4, width=5)
}

#Volcano plots of poised genes with gene labels
for (i in 1:length(stim_comps)){
  short_comp <- paste(strsplit(unstim_comps[i], "_")[[1]][1], strsplit(unstim_comps[i], "_")[[1]][2], sep="_")
  varlogFC <- paste(stim_comps[i],"logFC", sep=".")
  varneglog10p <- paste(stim_comps[i],"neglog10p", sep=".")
  
  tmp <- subset(full_folds[[short_comp]], 
                sig_full=="Poised_up"| sig_full =="Poised_dn" | sig_full=="None") %>%
        rename(logFC=matches(varlogFC), neglog10p=matches(varneglog10p))
  
  print(
    ggplot(tmp, aes(x=logFC, y=neglog10p)) +
    geom_point(aes(color=sig_full)) +
    geom_hline(yintercept = 1.3, linetype="dashed") +
    geom_vline(xintercept = 1, linetype="dashed") +
    geom_vline(xintercept = -1, linetype="dashed") +
    scale_color_manual(values=c("None"="grey", "Poised_dn"="light blue", "Poised_up"="orange")) +
    annotate("text", x = -5, y = 25, label = as.character(abs(sig_counts[grep(short_comp, sig_counts$subset),"Poised_dn"])), color="light blue") +
    annotate("text", x = 5, y = 25, label = as.character(sig_counts[grep(short_comp, sig_counts$subset),"Poised_up"]), color="orange") +
    geom_text_repel(aes(label=ifelse(logFC >4 | neglog10p > 10, as.character(SYMBOL), "" )),size=2)
  )
  ggsave(filename=here(output_dir, "volcano", paste0(stim_comps[i], "_volcanolabel.pdf")), height=4, width=5)
}

#Volcano plots of both category genes
for (i in 1:length(stim_comps)){
  short_comp <- paste(strsplit(unstim_comps[i], "_")[[1]][1], strsplit(unstim_comps[i], "_")[[1]][2], sep="_")
  varlogFC <- paste(stim_comps[i],"logFC", sep=".")
  varneglog10p <- paste(stim_comps[i],"neglog10p", sep=".")
  
  tmp <- subset(full_folds[[short_comp]], 
                sig_full=="Both_up_up" | sig_full=="Both_dn_dn" | sig_full=="None") %>%
    rename(logFC=matches(varlogFC), neglog10p=matches(varneglog10p))
  
  print(
    ggplot(tmp, aes(x=logFC, y=neglog10p)) +
      geom_point(aes(color=sig)) +
      geom_hline(yintercept = 1.3, linetype="dashed") +
      geom_vline(xintercept = 1, linetype="dashed") +
      geom_vline(xintercept = -1, linetype="dashed") +
      scale_color_manual(values=c("None"="grey", "Poised_dn"="light blue", "Poised_up"="orange", "Both"="purple")) +
      #annotate("text", x = -5, y = 25, label = as.character(abs(sig_counts[grep(short_comp, sig_counts$subset),"Poised_dn"])), color="light blue") +
      #annotate("text", x = 5, y = 25, label = as.character(sig_counts[grep(short_comp, sig_counts$subset),"Poised_up"]), color="orange") #+
      geom_text_repel(aes(label=ifelse(logFC >4 | neglog10p > 10, as.character(SYMBOL), "" )),size=2)
  )
  ggsave(filename=here(output_dir, "volcano", paste0(stim_comps[i], "both_volcano.pdf")), height=4, width=5)
}


#Recolored concordance

p_CD4_CM_colors <- ggplot(subset(full_folds$CD4_CM, sig !="None"), aes(x=Unstim_logFC, y=Stim_logFC, label=SYMBOL))
p_CD4_CM_colors + geom_point(aes(color=sig_full), size=1) +
  geom_abline(slope=1, intercept=0, ) +
  geom_abline(slope=1, intercept=-2, linetype="dashed") +
  geom_abline(slope=1, intercept=2, linetype="dashed") +
  scale_y_continuous(limits=c(-9, 9)) +
  scale_x_continuous(limits=c(-9, 9)) +
  scale_color_manual(name="DEG Category", values=c("None"="grey", "Poised_dn"="light blue", "Poised_up"="orange", "Both_dn_dn"="violet","Both_up_up"="purple", "Both_dn_up"="purple","Both_up_dn"="purple","Resting_dn"="cornflower blue", "Resting_up"="firebrick"))
  #geom_text(aes(label=ifelse(diff>diff_cut & sigtotal==1, as.character(SYMBOL),"")), hjust=1, vjust=0)
ggsave(here(output_dir, "CD4_CM_StimUnstim_Colorconcord.pdf"), width=4.5, height=3.5)

diff_cut <- 4

p_CD4_EM_colors <- ggplot(subset(full_folds$CD4_EM, sig !="None"), aes(x=Unstim_logFC, y=Stim_logFC, label=SYMBOL))
p_CD4_EM_colors + geom_point(aes(color=sig_full), size=1) +
  #geom_abline(slope=1, intercept=0, ) +
  #geom_abline(slope=1, intercept=-2, linetype="dashed") +
  #geom_abline(slope=1, intercept=2, linetype="dashed") +
  geom_hline(yintercept = 0, linetype="dashed") +
  geom_vline(xintercept = 0, linetype="dashed") +
  scale_y_continuous(limits=c(-9, 9)) +
  scale_x_continuous(limits=c(-9, 9)) +
  scale_color_manual(name="DEG Category", values=c("None"="grey", "Poised_dn"="light blue", "Poised_up"="orange", "Both_dn_dn"="violet","Both_up_up"="purple", "Both_dn_up"="purple","Both_up_dn"="purple","Resting_dn"="cornflower blue", "Resting_up"="firebrick"))+
  geom_text(aes(label=ifelse(diff>diff_cut & sigtotal==1, as.character(SYMBOL),"")), hjust=1, vjust=0)
ggsave(here(output_dir, "CD4_EM_StimUnstim_Colorconcord.pdf"), width=4.5, height=3.5)

p_CD8_CM_colors <- ggplot(subset(full_folds$CD8_CM, sig !="None"), aes(x=Unstim_logFC, y=Stim_logFC, label=SYMBOL))
p_CD8_CM_colors + geom_point(aes(color=sig_full), size=1) +
  geom_abline(slope=1, intercept=0, ) +
  geom_abline(slope=1, intercept=-2, linetype="dashed") +
  geom_abline(slope=1, intercept=2, linetype="dashed") +
  scale_y_continuous(limits=c(-9, 9)) +
  scale_x_continuous(limits=c(-9, 9)) +
  scale_color_manual(name="DEG Category", values=c("None"="grey", "Poised_dn"="light blue", "Poised_up"="orange", "Both_dn_dn"="violet","Both_up_up"="purple", "Both_dn_up"="purple","Both_up_dn"="purple","Resting_dn"="cornflower blue", "Resting_up"="firebrick"))
#geom_text(aes(label=ifelse(diff>diff_cut & sigtotal==1, as.character(SYMBOL),"")), hjust=1, vjust=0)
ggsave(here(output_dir, "CD8_CM_StimUnstim_Colorconcord.pdf"), width=4.5, height=3.5)

p_CD8_EM_colors <- ggplot(subset(full_folds$CD8_EM, sig !="None"), aes(x=Unstim_logFC, y=Stim_logFC, label=SYMBOL))
p_CD8_EM_colors + geom_point(aes(color=sig_full), size=1) +
  geom_abline(slope=1, intercept=0, ) +
  geom_abline(slope=1, intercept=-2, linetype="dashed") +
  geom_abline(slope=1, intercept=2, linetype="dashed") +
  scale_y_continuous(limits=c(-9, 9)) +
  scale_x_continuous(limits=c(-9, 9)) +
  scale_color_manual(name="DEG Category", values=c("None"="grey", "Poised_dn"="light blue", "Poised_up"="orange", "Both_dn_dn"="violet","Both_up_up"="purple", "Both_dn_up"="purple","Both_up_dn"="purple","Resting_dn"="cornflower blue", "Resting_up"="firebrick"))
#geom_text(aes(label=ifelse(diff>diff_cut & sigtotal==1, as.character(SYMBOL),"")), hjust=1, vjust=0)
ggsave(here(output_dir, "CD8_EM_StimUnstim_Colorconcord.pdf"), width=4.5, height=3.5)

p_CD8_EMRA_colors <- ggplot(subset(full_folds$CD8_EMRA, sig !="None"), aes(x=Unstim_logFC, y=Stim_logFC, label=SYMBOL))
p_CD8_EMRA_colors + geom_point(aes(color=sig_full), size=1) +
  geom_abline(slope=1, intercept=0, ) +
  geom_abline(slope=1, intercept=-2, linetype="dashed") +
  geom_abline(slope=1, intercept=2, linetype="dashed") +
  scale_y_continuous(limits=c(-9, 9)) +
  scale_x_continuous(limits=c(-9, 9)) +
  scale_color_manual(name="DEG Category", values=c("None"="grey", "Poised_dn"="light blue", "Poised_up"="orange", "Both_dn_dn"="violet","Both_up_up"="purple", "Both_dn_up"="purple","Both_up_dn"="purple","Resting_dn"="cornflower blue", "Resting_up"="firebrick"))
#geom_text(aes(label=ifelse(diff>diff_cut & sigtotal==1, as.character(SYMBOL),"")), hjust=1, vjust=0)
ggsave(here(output_dir, "CD8_EMRA_StimUnstim_Colorconcord.pdf"), width=4.5, height=3.5)

#Extra

ggplot(diff, aes(x=CD8_CM_blood_stim.v.CD8_EM_blood_stim.logFC, y=-log(CD8_CM_blood_stim.v.CD8_EM_blood_stim.fdr, base=10))) +
  geom_point(aes(color=CD8_CM_blood_stim.v.CD8_EM_blood_stim.sig)) +
  scale_x_continuous(limits=c(-10,10)) +
  geom_text_repel(aes(label=ifelse(CD8_CM_blood_stim.v.CD8_EM_blood_stim.sig=="TRUE", as.character(SYMBOL), "")))

#Saving important objects

save(full_folds, sig_counts_long, file=here(output_dir, "stim_unstim_all_full_folds.rda"))


CD4_EM_plotdata <- subset(full_folds$CD4_EM, sig !="None")
CD4_EM_plotdata$sig_full <- factor(CD4_EM_plotdata$sig_full)
CD4_EM_plotdata$sig_full <- recode(CD4_EM_plotdata$sig_full, Both_dn_dn="Mem_Augmented_dndn",
                                   Both_up_dn="Mem_Augmented_updn",
                                   Both_up_up="Mem_Augmented_upup",
                                   Poised_dn="Mem_Induced_dn",
                                   Poised_up="Mem_Induced_up",
                                   Resting_dn="Mem_Expressed_dn",
                                   Resting_up="Mem_Expressed_up")

genes <- c("IL22", "IL31", "IL4", "IL1A", "ADGRG1", "PLEKHG3", "CCL5", "EOMES", "LGALS3", "CXCR6", "MSC", "PHLDA2", "GZMH", "NKG7", "LMNA", "CXCR3", "TBX21")

ggplot(CD4_EM_plotdata, aes(x=Unstim_logFC, y=Stim_logFC, label=SYMBOL)) +
  geom_point(aes(color=sig_full), size=1) +
  geom_hline(yintercept=0, linetype="dashed") + 
  geom_vline(xintercept=0, linetype="dashed") + 
  theme_light() +
  theme(axis.title =element_text(size=8)) +
  scale_y_continuous(limits=c(-10, 10)) +
  scale_x_continuous(limits=c(-10, 10)) +
  xlab("RNA log2FC (Mem Unstim/Naive Unstim)") +
  ylab("RNA log2FC (Mem Stim/Naive Stim)") +
  scale_color_manual(name="DEG Category", values=c("None"="grey", "Mem_Induced_dn"="light blue", "Mem_Induced_up"="orange", "Mem_Augmented_dndn"="violet","Mem_Augmented_upup"="purple", "Mem_Augmented_dnup"="purple","Mem_Augmented_updn"="purple","Mem_Expressed_dn"="cornflower blue", "Mem_Expressed_up"="firebrick")) +
  #geom_text_repel(aes(label=ifelse(diff>diff_cut, as.character(SYMBOL),"")), hjust=1, vjust=0) +
  #geom_text_repel(aes(label=ifelse(abs(Unstim_logFC)>2.5 & abs(Stim_logFC)>2.5, as.character(SYMBOL),"")), hjust=1, vjust=0)
  geom_text_repel(aes(label=ifelse(SYMBOL %in% genes, as.character(SYMBOL), "")))
ggsave(here(output_dir, "CD4_EM_StimUnstim_Colorconcord2.pdf"), width=9.5, height=7.5)
