#Heatmap of NES scores for Hallmark genesets of both Stim and Unstim
#20Nov20
#jrose

library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(here)

unstim_dir <- "GSEA/Hallmark_24Feb20"
stim_dir <- "GSEA/Hallmark_Stim"

#Loading unstim GSEA data and wrangling out NES
load(here(unstim_dir, "NES_data.rda"))
dat_unstim <- dat %>% mutate(condition="Unstim")
rm(dat)

unstim_NES_mat <- dat_unstim %>% select(NAME, analysis, NES) %>%
    pivot_wider(names_from=analysis, values_from=NES)

#Loading stim GSEA data and wranlging out NES
load(here(stim_dir, "NES_Stim_hallmark_data.rda"))
dat_stim <- dat %>% mutate(condition="Stim")
rm(dat)

stim_NES_mat <- dat_stim %>% select(NAME, analysis, NES) %>%
  pivot_wider(names_from=analysis, values_from=NES)

#Joining datasets
NES_df <- inner_join(unstim_NES_mat, stim_NES_mat, by="NAME") %>%
  select(NAME, 
         CD4_CM_Nav, CD4_EM_Nav, CD8_CM_Nav, CD8_EM_Nav, CD8_EMRA_Nav,
         CD4_CM_Nav_stim, CD4_EM_Nav_stim, CD8_CM_Nav_stim, CD8_EM_Nav_stim, CD8_EMRA_Nav_stim)
NES_mat <- as.matrix(NES_df[,-1])
rownames(NES_mat) <- gsub("HALLMARK_", "", NES_df$NAME)
  
#Heatmap

ph_anno <- data.frame(condition=c(rep("Unstim", 5), rep("Stim", 5)),
                      subtype=c(rep(c("CD4_CM", "CD4_EM", "CD8_CM", "CD8_EM", "CD8_EMRA"),2)))
rownames(ph_anno) <- colnames(NES_mat)

col_fun <- colorRamp2(breaks=c(min(NES_mat),0, max(NES_mat)), colors = c("lightblue", "black","khaki1"), space = "sRGB")
custom_pal <- c("#90A4ADFF","#5F7D8BFF", "#2096F2FF", "#0C46A0FF", "#FFB74CFF", "#E55100FF", "#870D4EFF")
names(custom_pal) <- c("CD4_Nav", "CD8_Nav", "CD4_CM", "CD8_CM", "CD4_EM", "CD8_EM", "CD8_EMRA")

pdf(file=here(stim_dir, "NES_heat_StimUnstim.pdf"), useDingbats = F, height=8, width=7)
Heatmap(NES_mat,
        cluster_columns = F,
        col=col_fun,
        row_names_gp = gpar(fontsize = 8),
        show_column_names = F,
        column_split = rev(ph_anno$condition),
        column_title = NULL,
        top_annotation = HeatmapAnnotation(Stim=ph_anno$condition,
                                           Subtype=ph_anno$subtype,
                                           col=list(Stim=c("Unstim"="grey", "Stim"="black"),
                                                    Subtype=custom_pal)
                                          )
        )
dev.off()


#Highlight only certain pathways
path_highlight = c("IL2_STAT5_SIGNALING", "COMPLEMENT", "INTERFERON_GAMMA_RESPONSE",
                   "ALLOGRAFT_REJECTION", "CHOLESTEROL_HOMEOSTASIS", "INTERFERON_ALPHA_RESPONSE",
                   "INFLAMMATORY_RESPONSE", "FATTY_ACID_METABOLISM", "PI3K_AKT_MTOR_SIGNALING",
                   "PEROXISOME", "TGF_BETA_SIGNALING", "PROTEIN_SECRETION", "HYPOXIA", "TNFA_SIGNALING_VIA_NFKB",
                   "UNFOLDED_PROTEIN_RESPONSE", "MTORC1_SIGNALING", "COAGULATION", "KRAS_SIGNALING_DN",
                   "GLYCOLYSIS", "G2M_CHECKPOINT", "E2F_TARGETS", "OXIDATIVE_PHOSPHORYLATION",
                   "NOTCH_SIGNALING", "MYC_TARGETS_V1", "WNT_BETA_CATENIN_SIGNALING")

lbl_index <- vector()
for (i in 1:length(path_highlight)){
  if (path_highlight[i]%in%rownames(NES_mat)) {
    lbl_index[i] <- grep(path_highlight[i], rownames(NES_mat))
  }
}


pdf(file=here(stim_dir, "NES_heat_StimUnstim_highlight.pdf"), useDingbats = F, height=4, width=5)
Heatmap(NES_mat,
        cluster_columns = F,
        col=col_fun,
        row_names_gp = gpar(fontsize = 8),
        show_column_names = F,
        show_row_names = F,
        column_split = rev(ph_anno$condition),
        column_title = NULL,
        top_annotation = HeatmapAnnotation(Stim=ph_anno$condition,
                                           Subtype=ph_anno$subtype,
                                           col=list(Stim=c("Unstim"="grey", "Stim"="black"),
                                                    Subtype=custom_pal)
        ),
        right_annotation = rowAnnotation(gene = anno_mark(at = lbl_index, 
                                                          labels = path_highlight,
                                                          labels_gp = gpar(fontsize=7)))
)
dev.off()