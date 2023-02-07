
library(tidyverse)
library(ggpubr)
library(fmsb)
library(here)

load(here("kmeans/GO/CD4_unstim_labfilter_c5allv7_data.rda"))
outdir="kmeans/GO"

df <- data.frame()

for (i in 1:length(df_clus_lst)){
  df <- bind_rows(df, df_clus_lst[[i]])
}

plots <- list()
for (i in 1:length(df_clus_polish)){
  plots[[i]] <- ggplot(df_clus_polish[[i]] %>% top_n(10, neglog10P), aes(x=neglog10P, y=geneset)) +
    geom_col(fill="cornflowerblue") +
    geom_text(aes(label=geneset), hjust=0.25, size=2) +
    scale_x_continuous(limits=c(0,15)) +
    theme_light() +
    theme(axis.text.y = element_blank()) +
    geom_vline(xintercept=1.3, linetype="dashed")
}

plots[[1]]
ggsave(filename=here(outdir, "CD4_unstim_labfilter_cluster1fisher.pdf"), height=4, width=4)
plots[[2]]
ggsave(filename=here(outdir, "CD4_unstim_labfilter_cluster2fisher.pdf"), height=4, width=4)
plots[[3]]
ggsave(filename=here(outdir, "CD4_unstim_labfilter_cluster3fisher.pdf"), height=4, width=4)
plots[[4]]
ggsave(filename=here(outdir, "CD4_unstim_labfilter_cluster4fisher.pdf"), height=4, width=4)
plots[[5]]
ggsave(filename=here(outdir, "CD4_unstim_labfilter_cluster5fisher.pdf"), height=4, width=4)



genesets=c("GO_METHYLATION","GO_LYMPHOCYTE_COSTIMULATION","GO_T_CELL_DIFFERENTIATION","GO_REGULATION_OF_T_CELL_DIFFERENTIATION_IN_THYMUS","GO_POSITIVE_REGULATION_OF_T_HELPER_1_TYPE_IMMUNE_RESPONSE","GO_T_HELPER_1_CELL_DIFFERENTIATION","GO_T_HELPER_2_CELL_DIFFERENTIATION","GO_T_HELPER_17_TYPE_IMMUNE_RESPONSE","GO_REGULATION_OF_T_CELL_DIFFERENTIATION","GO_TYPE_2_IMMUNE_RESPONSE","GO_T_CELL_MEDIATED_IMMUNITY","GO_T_CELL_MEDIATED_CYTOTOXICITY","GO_REGULATORY_T_CELL_DIFFERENTIATION","GO_ACTIVATION_OF_INNATE_IMMUNE_RESPONSE","GO_NATURAL_KILLER_CELL_MEDIATED_IMMUNITY","GO_POSITIVE_REGULATION_OF_CHEMOTAXIS","GO_LYMPHOCYTE_CHEMOTAXIS","GO_POSITIVE_REGULATION_OF_LYMPHOCYTE_CHEMOTAXIS","GO_T_CELL_MIGRATION","GO_HETEROTYPIC_CELL_CELL_ADHESION","GO_ANTIGEN_PROCESSING_AND_PRESENTATION_OF_PEPTIDE_ANTIGEN","GO_T_CELL_APOPTOTIC_PROCESS","GO_REGULATION_OF_T_CELL_APOPTOTIC_PROCESS","GO_POSITIVE_REGULATION_OF_INTERFERON_GAMMA_PRODUCTION","GO_REGULATION_OF_CYTOKINE_PRODUCTION_INVOLVED_IN_IMMUNE_RESPONSE","GO_CD4_POSITIVE_ALPHA_BETA_T_CELL_CYTOKINE_PRODUCTION","GO_POSITIVE_REGULATION_OF_LEUKOCYTE_MIGRATION","GO_NEGATIVE_REGULATION_OF_T_CELL_MEDIATED_IMMUNITY","GO_INTERLEUKIN_2_PRODUCTION","GO_INTRACELLULAR_LIPID_TRANSPORT","GO_OXIDATIVE_PHOSPHORYLATION","GO_ATP_SYNTHESIS_COUPLED_ELECTRON_TRANSPORT","GO_MITOCHONDRIAL_ELECTRON_TRANSPORT_NADH_TO_UBIQUINONE","GO_CARBOHYDRATE_CATABOLIC_PROCESS","GO_SULFUR_COMPOUND_CATABOLIC_PROCESS","GO_ACETYL_COA_BIOSYNTHETIC_PROCESS","GO_UNSATURATED_FATTY_ACID_BIOSYNTHETIC_PROCESS","GO_LONG_CHAIN_FATTY_ACYL_COA_BIOSYNTHETIC_PROCESS","GO_PURINERGIC_NUCLEOTIDE_RECEPTOR_SIGNALING_PATHWAY","GO_AMINE_CATABOLIC_PROCESS","GO_PROTEIN_KINASE_A_SIGNALING","GO_PROTEIN_KINASE_C_SIGNALING","GO_REGULATION_OF_PROTEIN_TYROSINE_KINASE_ACTIVITY","GO_ACTIVATION_OF_MAPK_ACTIVITY","GO_POSITIVE_REGULATION_OF_MAP_KINASE_ACTIVITY","GO_MYD88_DEPENDENT_TOLL_LIKE_RECEPTOR_SIGNALING_PATHWAY","GO_ENDOPLASMIC_RETICULUM_CALCIUM_ION_HOMEOSTASIS","GO_REGULATION_OF_VOLTAGE_GATED_CALCIUM_CHANNEL_ACTIVITY","GO_LYMPHOCYTE_ACTIVATION_INVOLVED_IN_IMMUNE_RESPONSE","GO_ALPHA_BETA_T_CELL_ACTIVATION","GO_REGULATION_OF_T_CELL_ACTIVATION","GO_T_CELL_ACTIVATION","GO_REGULATION_OF_T_CELL_RECEPTOR_SIGNALING_PATHWAY","GO_T_CELL_RECEPTOR_SIGNALING_PATHWAY","GO_NEGATIVE_REGULATION_OF_INFLAMMATORY_RESPONSE","GO_NEGATIVE_REGULATION_OF_CD4_POSITIVE_ALPHA_BETA_T_CELL_ACTIVATION","GO_CELL_AGING","GO_CYTOLYSIS","GO_REGULATION_OF_CELL_KILLING","GO_CELL_KILLING","GO_LEUKOCYTE_MEDIATED_CYTOTOXICITY","GO_REGULATION_OF_INNATE_IMMUNE_RESPONSE", "GO_IMMUNE_EFFECTOR_PROCESS", "GO_SECRETORY_GRANULE","GO_MEMBRANE_LIPID_BIOSYNTHETIC_PROCESS", "GO_COFACTOR_BIOSYNTHETIC_PROCESS", "GO_CYTOKINE_SECRETION", "GO_POSITIVE_REGULATION_OF_CELL_GROWTH")
#genesets=c("GO_LYMPHOCYTE_COSTIMULATION", "GO_CYTOLYSIS", "GO_INTERLEUKIN_2_PRODUCTION","GO_METHYLATION","GO_POSITIVE_REGULATION_OF_MAP_KINASE_ACTIVITY","GO_REGULATION_OF_T_CELL_DIFFERENTIATION_IN_THYMUS","GO_CELL_KILLING","GO_REGULATION_OF_CYTOKINE_PRODUCTION_INVOLVED_IN_IMMUNE_RESPONSE","GO_T_CELL_DIFFERENTIATION","GO_T_HELPER_2_CELL_DIFFERENTIATION","GO_T_CELL_MEDIATED_CYTOTOXICITY","GO_ANTIGEN_PROCESSING_AND_PRESENTATION_OF_PEPTIDE_ANTIGEN","GO_LEUKOCYTE_MEDIATED_CYTOTOXICITY","GO_PROTEIN_KINASE_C_SIGNALING","GO_T_CELL_APOPTOTIC_PROCESS","GO_ACTIVATION_OF_MAPK_ACTIVITY","GO_REGULATION_OF_PROTEIN_TYROSINE_KINASE_ACTIVITY","GO_T_CELL_ACTIVATION","GO_T_CELL_RECEPTOR_SIGNALING_PATHWAY","GO_REGULATION_OF_T_CELL_ACTIVATION")

df_currated <- subset(df, geneset%in%genesets)%>% mutate(neglog10P=-log(pvalue, base=10))
df_currated_radar <- subset(df, geneset%in%genesets) %>% mutate(neglog10P=-log(pvalue, base=10)) %>%
  select(k.5.cluster,neglog10P, geneset) %>%
  pivot_wider(names_from=k.5.cluster, values_from=neglog10P)

df_currated_sig <- subset(df_currated, neglog10P>=1.3) %>% arrange(desc(neglog10P))
df_currated$geneset <- factor(df_currated$geneset, levels=unique(df_currated_sig$geneset), ordered=T)
df_currated_sig <- arrange(df_currated_sig, neglog10P)
df_currated_sig$geneset <- factor(df_currated_sig$geneset, levels=unique(df_currated_sig$geneset), ordered=T)

df_currated$neglog10P[df_currated$neglog10P<0.1] <- 0.1

ggplot(df_currated_sig, aes(x=geneset, y=neglog10P)) +
  geom_col(aes(group=k.5.cluster, fill=as.factor(k.5.cluster)),position=position_dodge(width = 1)) +
  geom_hline(yintercept=1.3, linetype="dashed") +
  #scale_y_continuous(trans="reverse") +
  coord_flip() +
  theme_light() + theme(axis.text.x=element_text(hjust=1, vjust=1)) +
  guides(fill=guide_legend(title="Module"))
ggsave(filename=here(outdir, "CD4_unstim_labfilter_cmeanscolorfisherplot.pdf"), height=7, width=10)


ggplot(df_currated, aes(x=geneset, y=neglog10P)) +
  geom_col(aes(group=k.5.cluster, fill=as.factor(k.5.cluster)),position=position_dodge2(width = 1, padding=0.3)) +
  geom_hline(yintercept=1.3, linetype="dashed") +
  #scale_y_continuous(trans="reverse") +
  #coord_flip() +
  theme_light() + theme(axis.text.x=element_text(hjust=1, vjust=1, angle=45))
ggsave(filename=here(outdir, "CD4_unstim_labfilter_widefisherplot.pdf"), height=7, width=15)
