#Script for comparing logFC of pagernk scores via scatter plot
setwd("~/Documents/Emory/Boss/data/ATAC")

library(tidyverse)
library(ggrepel)
library(here)

input_dir <- "PageRank/Run1"
output_dir <- "PageRank/Run1/scatter_logFCNav"

d <- read.table(here(input_dir,"TFRanks.tsv"),header=T,row.names=1)

#Filter
x1 = apply(d, 1, max)>1e-4
d<-d[x1,]

#Calculate folds for each score over naive
d1<-data.frame(d$CD4_CM/d$CD4_N,
               d$CD4_EM/d$CD4_N, 
               d$CD4_EMRA/d$CD4_N, 
               d$CD8_CM/d$CD8_N, 
               d$CD8_EM/d$CD8_N, 
               d$CD8_EMRA/d$CD8_N,
               d$CD4_CMStim/d$CD4_NStim,
               d$CD4_EMStim/d$CD4_NStim, 
               d$CD4_EMRAStim/d$CD4_NStim, 
               d$CD8_CMStim/d$CD8_NStim, 
               d$CD8_EMStim/d$CD8_NStim, 
               d$CD8_EMRAStim/d$CD8_NStim,
               d$CD4_CMStim/d$CD4_N,
               d$CD4_EMStim/d$CD4_N, 
               d$CD4_EMRAStim/d$CD4_N, 
               d$CD8_CMStim/d$CD8_N, 
               d$CD8_EMStim/d$CD8_N, 
               d$CD8_EMRAStim/d$CD8_N)

#Comparisons
# CD8_EM_EMRA
# CD8_CM_EMRA
# CD8_EM_CM
# CD4_EM_CM
# CD8_EMRAStim_EMStim
# CD8_EMRAStim_CMStim
# CD8_EMStim_CMStim
# CD4_EMStim_CMStim
# CD8_EMStim_EM
# CD8_CMStim_CM
# CD8_EMRAStim_EMRA
# CD4_EMStim_EM
# CD4_CMStim_CM

#Log2 transformation of folds
d1 <- apply(d1, c(1,2), function(x) log(x,base=2)) %>% data.frame()
d$gene <- rownames(d)
d1$gene <- rownames(d)

#Combining into large df of both logFC and raw scores
d_full <- left_join(d1, d, by="gene") %>% select(gene, everything())

#Adding differences between folds (for labeling purposes) and mean scores per group

d_full <- mutate(d_full, CD8_EMRA_EM_diff=abs(d.CD8_EMRA.d.CD8_N-d.CD8_EM.d.CD8_N),
             CD8_EMRA_EM_mean = log((CD8_EMRA+CD8_EM)/2, base=10),
             CD8_EMRA_CM_diff=abs(d.CD8_EMRA.d.CD8_N-d.CD8_CM.d.CD8_N),
             CD8_EMRA_CM_mean = log((CD8_EMRA+CD8_CM)/2, base=10),
             CD8_EM_CM_diff = abs(d.CD8_EM.d.CD8_N-d.CD8_CM.d.CD8_N),
             CD8_EM_CM_mean = log((CD8_EM+CD8_CM)/2, base=10),
             # CD4_EMRA_EM_diff=abs(d.CD4_EMRA.d.CD4_N-d.CD4_EM.d.CD4_N),
             # CD4_EMRA_EM_mean = mean(CD4_EMRA, CD4_EM),
             # CD4_EMRA_CM_diff=abs(d.CD4_EMRA.d.CD4_N-d.CD4_CM.d.CD4_N),
             # CD4_EMRA_CM_mean = mean(CD4_EMRA, CD4_CM),
             CD4_EM_CM_diff = abs(d.CD4_EM.d.CD4_N-d.CD4_CM.d.CD4_N),
             CD4_EM_CM_mean = log((CD4_EM+CD4_CM)/2, base=10),
             CD8_EMRAStim_EMStim_diff=abs(d.CD8_EMRAStim.d.CD8_NStim-d.CD8_EMStim.d.CD8_NStim),
             CD8_EMRAStim_EMStim_mean = log((CD8_EMRAStim+CD8_EMStim)/2, base=10),
             CD8_EMRAStim_CMStim_diff=abs(d.CD8_EMRAStim.d.CD8_NStim-d.CD8_CMStim.d.CD8_NStim),
             CD8_EMRAStim_CMStim_mean = log((CD8_EMRAStim+CD8_CMStim)/2, base=10),
             CD8_EMStim_CMStim_diff = abs(d.CD8_EMStim.d.CD8_NStim-d.CD8_CMStim.d.CD8_NStim),
             CD8_EMStim_CMStim_mean = log((CD8_EMStim+CD8_CMStim)/2, base=10),
             # CD4_EMRAStim_EMStim_diff=abs(d.CD4_EMRAStim.d.CD4_NStim-d.CD4_EMStim.d.CD4_NStim),
             # CD4_EMRAStim_EMStim_mean = mean(CD4_EMRAStim, CD4_EMStim),
             # CD4_EMRAStim_CMStim_diff=abs(d.CD4_EMRAStim.d.CD4_NStim-d.CD4_CMStim.d.CD4_NStim),
             # CD4_EMRAStim_CMStim_mean = mean(CD4_EMRAStim, CD4_CMStim),
             CD4_EMStim_CMStim_diff = abs(d.CD4_EMStim.d.CD4_NStim-d.CD4_CMStim.d.CD4_NStim),
             CD4_EMStim_CMStim_mean = log((CD4_EMStim+CD4_CMStim)/2, base=10),
             CD8_EMStim_EM_diff = abs(d.CD8_EMStim.d.CD8_N-d.CD8_EM.d.CD8_N),
             CD8_EMStim_EM_mean = log((CD8_EMStim+CD8_EM)/2, base=10),
             CD8_CMStim_CM_diff = abs(d.CD8_CMStim.d.CD8_N-d.CD8_CM.d.CD8_N),
             CD8_CMStim_CM_mean = log((CD8_CMStim+CD8_CM)/2, base=10),
             CD8_EMRAStim_EMRA_diff = abs(d.CD8_EMRAStim.d.CD8_N-d.CD8_EMRA.d.CD8_N),
             CD8_EMRAStim_EMRA_mean = log((CD8_EMRAStim+CD8_EMRA)/2, base=10),
             CD4_EMStim_EM_diff = abs(d.CD4_EMStim.d.CD4_N-d.CD4_EM.d.CD4_N),
             CD4_EMStim_EM_mean = log((CD4_EMStim+CD4_EM)/2, base=10),
             CD4_CMStim_CM_diff = abs(d.CD4_CMStim.d.CD4_N-d.CD4_CM.d.CD4_N),
             CD4_CMStim_CM_mean = log((CD4_CMStim+CD4_CM)/2, base=10),
             )

top <- list()
top$CD8_EM_EMRA <- arrange(d_full, desc(CD8_EMRA_EM_mean)) %>% top_n(50, wt=CD8_EMRA_EM_diff)
top$CD8_CM_EMRA <- arrange(d_full, desc(CD8_EMRA_CM_mean)) %>% top_n(50, wt=CD8_EMRA_CM_diff)
top$CD8_EM_CM <- arrange(d_full, desc(CD8_EM_CM_mean)) %>% top_n(50, wt=CD8_EM_CM_diff)
top$CD4_EM_CM <- arrange(d_full, desc(CD4_EM_CM_mean)) %>% top_n(50, wt=CD4_EM_CM_diff)
top$CD8_EMRAStim_EMStim <- arrange(d_full, desc(CD8_EMRAStim_EMStim_mean)) %>% top_n(50, wt=CD8_EMRAStim_EMStim_diff)
top$CD8_EMRAStim_CMStim <- arrange(d_full, desc(CD8_EMRAStim_CMStim_mean)) %>% top_n(50, wt=CD8_EMRAStim_CMStim_diff)
top$CD8_EMStim_CMStim <- arrange(d_full, desc(CD8_EMStim_CMStim_mean)) %>% top_n(50, wt=CD8_EMStim_CMStim_diff)
top$CD4_EMStim_CMStim <- arrange(d_full, desc(CD4_EMStim_CMStim_mean)) %>% top_n(50, wt=CD4_EMStim_CMStim_diff)
top$CD8_EMStim_EM <- arrange(d_full, desc(CD8_EMStim_EM_mean)) %>% top_n(50, wt=CD8_EMStim_EM_diff)
top$CD8_CMStim_CM <- arrange(d_full, desc(CD8_CMStim_CM_mean)) %>% top_n(50, wt=CD8_CMStim_CM_diff)
top$CD8_EMRAStim_EMRA <- arrange(d_full, desc(CD8_EMRAStim_EMRA_mean)) %>% top_n(50, wt=CD8_EMRAStim_EMRA_diff)
top$CD4_EMStim_EM <- arrange(d_full, desc(CD4_EMStim_EM_mean)) %>% top_n(50, wt=CD4_EMStim_EM_diff)
top$CD4_CMStim_CM <- arrange(d_full, desc(CD4_CMStim_CM_mean)) %>% top_n(50, wt=CD4_CMStim_CM_diff)

save(top, file=here(output_dir,"topOrdered_TmemSubsets_PgRnk_list.rda"))

#Plots

theme_set(theme_gray()+ theme(axis.line = element_line(size=0.5),
                              panel.background = element_rect(fill=NA,size=rel(20)),
                              panel.grid.minor = element_line(colour = NA),
                              axis.text = element_text(size=16), axis.title = element_text(size=13)))

p_CD8_EM_EMRA <- ggplot(d_full, aes(x=d.CD8_EM.d.CD8_N, y=d.CD8_EMRA.d.CD8_N))
# p_CD8_EM_EMRA + geom_point(aes(size=CD8_EMRA_EM_mean, color=CD8_EMRA_EM_mean)) +
#   geom_abline(slope=1, intercept=0, ) +
#   geom_text(aes(label=ifelse(CD8_EMRA_EM_diff>1.5, as.character(gene),"")), hjust=1, vjust=0, size=3) +
#   xlab("CD8 EM vs Naive") +
#   ylab("CD8 EMRA vs Naive") +
#   scale_x_continuous(limits=c(-5,5)) +
#   scale_y_continuous(limits=c(-5,5)) +
#   guides(size=F, color=F)
p_CD8_EM_EMRA + geom_point(size=1) +
  geom_abline(slope=1, intercept=0, ) +
  geom_abline(slope=1, intercept=1.5, linetype="dashed") +
  geom_abline(slope=1, intercept=-1.5, linetype="dashed") +
  #geom_text(aes(label=ifelse(gene %in% top$CD8_EM_CM$gene, as.character(gene),"")), hjust=1, vjust=0, size=3) +
  geom_text_repel(aes(label=ifelse(CD8_EMRA_EM_diff>=1.5, as.character(gene),"")), 
                  #hjust=1.1, 
                  #vjust=0, 
                  size=3) +
  xlab("PageRank log2FC (CD8 EM vs Naive)") +
  ylab("PageRank log2FC (CD8 EMRA vs Naive)")+
  scale_x_continuous(limits=c(-5,5)) +
  scale_y_continuous(limits=c(-5,5)) +
  guides(size=F, color=F)

ggsave(filename=here(output_dir,"CD8_EM_EMRA_pgrnklofFC_pub.pdf"), height=4, width=4)

#ggsave(filename=here(output_dir,"CD8_EM_EMRA_pgrnklofFC_2.pdf"), height=7, width=7)

p_CD8_CM_EMRA <- ggplot(d_full, aes(x=d.CD8_CM.d.CD8_N, y=d.CD8_EMRA.d.CD8_N))
# p_CD8_CM_EMRA + geom_point(aes(size=CD8_EMRA_CM_mean, color=CD8_EMRA_CM_mean)) +
#   geom_abline(slope=1, intercept=0) +
#   geom_text(aes(label=ifelse(CD8_EMRA_CM_diff>1.5, as.character(gene),"")), hjust=1, vjust=0, size=3) +
#   xlab("CD8 CM vs Naive") +
#   ylab("CD8 EMRA vs Naive") +
#   scale_x_continuous(limits=c(-5,5)) +
#   scale_y_continuous(limits=c(-5,5)) +
#   guides(size=F, color=F)

p_CD8_CM_EMRA + geom_point(aes(size=CD8_EMRA_CM_mean, color=CD8_EMRA_CM_mean)) +
  geom_abline(slope=1, intercept=0) +
  geom_text(aes(label=ifelse(gene %in% top$CD8_CM_EMRA$gene, as.character(gene),"")), hjust=1, vjust=0, size=3) +
  xlab("CD8 CM vs Naive") +
  ylab("CD8 EMRA vs Naive") +
  scale_x_continuous(limits=c(-5,5)) +
  scale_y_continuous(limits=c(-5,5)) +
  guides(size=F, color=F)

ggsave(filename=here(output_dir,"CD8_CM_EMRA_pgrnklofFC_2.pdf"), height=7, width=7)

p_CD8_EM_CM <- ggplot(d_full, aes(x=d.CD8_EM.d.CD8_N, y=d.CD8_CM.d.CD8_N))
p_CD8_EM_CM + geom_point(size=1) +
  geom_abline(slope=1, intercept=0, ) +
  geom_abline(slope=1, intercept=1.5, linetype="dashed") +
  geom_abline(slope=1, intercept=-1.5, linetype="dashed") +
  #geom_text(aes(label=ifelse(gene %in% top$CD8_EM_CM$gene, as.character(gene),"")), hjust=1, vjust=0, size=3) +
  geom_text_repel(aes(label=ifelse(CD8_EM_CM_diff>=1.5, as.character(gene),"")), 
                  #hjust=1.1, 
                  #vjust=0, 
                  size=3) +
  xlab("PageRank log2FC (CD8 EM vs Naive)") +
  ylab("PageRank log2FC (CD8 CM vs Naive)")+
  scale_x_continuous(limits=c(-5,5)) +
  scale_y_continuous(limits=c(-5,5)) +
  guides(size=F, color=F)

ggsave(filename=here(output_dir,"CD8_EM_CM_pgrnklofFC_pub.pdf"), height=4, width=4)

#CD4
p_CD4_EM_CM <- ggplot(d_full, aes(x=d.CD4_EM.d.CD4_N, y=d.CD4_CM.d.CD4_N))
p_CD4_EM_CM + geom_point(size=1) +
  geom_abline(slope=1, intercept=0, ) +
  geom_abline(slope=1, intercept=1.5, linetype="dashed") +
  geom_abline(slope=1, intercept=-1.5, linetype="dashed") +
  #geom_text(aes(label=ifelse(gene %in% top$CD8_EM_CM$gene, as.character(gene),"")), hjust=1, vjust=0, size=3) +
  geom_text_repel(aes(label=ifelse(CD4_EM_CM_diff>=1.5, as.character(gene),"")), 
                  #hjust=1.1, 
                  #vjust=0, 
                  size=3) +
  xlab("PageRank log2FC (CD4 EM vs Naive)") +
  ylab("PageRank log2FC (CD4 CM vs Naive)")+
  scale_x_continuous(limits=c(-5,5)) +
  scale_y_continuous(limits=c(-5,5)) +
  guides(size=F, color=F)

ggsave(filename=here(output_dir,"CD4_CM_EM_pgrnklofFC_pub.pdf"), height=4, width=4)

#Stimulated
p_CD8_EMStim_EMRAStim <- ggplot(d_full, aes(x=d.CD8_EMStim.d.CD8_NStim, y=d.CD8_EMRAStim.d.CD8_NStim))
p_CD8_EMStim_EMRAStim + geom_point(aes(size=CD8_EMRAStim_EMStim_mean, color=CD8_EMRAStim_EMStim_mean)) +
  geom_abline(slope=1, intercept=0, ) +
  geom_text(aes(label=ifelse(gene %in% top$CD8_EMRAStim_EMStim$gene, as.character(gene),"")), hjust=1, vjust=0, size=3) +
  xlab("CD8 EMStim vs NaiveStim") +
  ylab("CD8 EMRAStim vs NaiveStim") +
  scale_x_continuous(limits=c(-5,5)) +
  scale_y_continuous(limits=c(-5,5)) +
  guides(size=F, color=F)
ggsave(filename=here(output_dir,"CD8_EMStim_EMRAStim_pgrnklofFC_2.pdf"), height=7, width=7)

p_CD8_CMStim_EMRAStim <- ggplot(d_full, aes(x=d.CD8_CMStim.d.CD8_NStim, y=d.CD8_EMRAStim.d.CD8_NStim))
p_CD8_CMStim_EMRAStim + geom_point(aes(size=CD8_EMRAStim_CMStim_mean, color=CD8_EMRAStim_CMStim_mean)) +
  geom_abline(slope=1, intercept=0, ) +
  geom_text(aes(label=ifelse(gene %in% top$CD8_EMRAStim_CMStim$gene, as.character(gene),"")), hjust=1, vjust=0, size=3) +
  xlab("CD8 CMStim vs NaiveStim") +
  ylab("CD8 EMRAStim vs NaiveStim") +
  scale_x_continuous(limits=c(-5,5)) +
  scale_y_continuous(limits=c(-5,5)) +
  guides(size=F, color=F)
ggsave(filename=here(output_dir,"CD8_CMStim_EMRAStim_pgrnklofFC_2.pdf"), height=7, width=7)

p_CD8_EMStim_CMStim <- ggplot(d_full, aes(x=d.CD8_EMStim.d.CD8_NStim, y=d.CD8_CMStim.d.CD8_NStim))
p_CD8_EMStim_CMStim + geom_point(aes(size=CD8_EMStim_CMStim_mean, color=CD8_EMStim_CMStim_mean)) +
  geom_abline(slope=1, intercept=0, ) +
  geom_text(aes(label=ifelse(gene %in% top$CD8_EMStim_CMStim$gene, as.character(gene),"")), hjust=1, vjust=0, size=3) +
  xlab("CD8 EMStim vs NaiveStim") +
  ylab("CD8 CMStim vs NaiveStim") +
  scale_x_continuous(limits=c(-5,5)) +
  scale_y_continuous(limits=c(-5,5)) +
  guides(size=F, color=F)
ggsave(filename=here(output_dir,"CD8_EMStim_CMStim_pgrnklofFC_2.pdf"), height=7, width=7)

p_CD4_EMStim_CMStim <- ggplot(d_full, aes(x=d.CD4_EMStim.d.CD4_NStim, y=d.CD4_CMStim.d.CD4_NStim))
p_CD4_EMStim_CMStim + geom_point(aes(size=CD4_EMStim_CMStim_mean, color=CD4_EMStim_CMStim_mean)) +
  geom_abline(slope=1, intercept=0, ) +
  geom_text(aes(label=ifelse(gene %in% top$CD4_EMStim_CMStim$gene, as.character(gene),"")), hjust=1, vjust=0, size=3) +
  xlab("CD4 EMStim vs NaiveStim") +
  ylab("CD4 CMStim vs NaiveStim") +
  scale_x_continuous(limits=c(-5,5)) +
  scale_y_continuous(limits=c(-5,5)) +
  guides(size=F, color=F)

ggsave(filename=here(output_dir,"CD4_EMStim_CMStim_pgrnklofFC_2.pdf"), height=7, width=7)

#Stim vs Unstim

stim_TF <- c("ELK1",
             "FOS",
             "JUN",
             "NFATC1",
             "NFATC2",
             "NFATC3",
             "NFATC4",
             "NFKB1",
             "NFKBIA",
             "RELA")

p_CD8_EMStim_EM <- ggplot(d_full, aes(x=d.CD8_EMStim.d.CD8_N, y=d.CD8_EM.d.CD8_N))
p_CD8_EMStim_EM + geom_point(aes(size=CD8_EMStim_EM_mean, color=CD8_EMStim_EM_mean)) +
  geom_abline(slope=1, intercept=0, ) +
  geom_text(aes(label=ifelse(gene %in% top$CD8_EMStim_EM$gene, as.character(gene),"")), hjust=1, vjust=0, size=3) +
  xlab("CD8 EMStim vs Naive") +
  ylab("CD8 EM vs Naive") +
  scale_x_continuous(limits=c(-5,5)) +
  scale_y_continuous(limits=c(-5,5)) +
  guides(size=F, color=F)
ggsave(filename=here(output_dir,"CD8_EMStim_EM_pgrnklofFC_2.pdf"), height=7, width=7)

p_CD8_EMStim_EM + geom_point(aes(size=CD8_EMStim_EM_mean, color=CD8_EMStim_EM_mean)) +
  geom_abline(slope=1, intercept=0, ) +
  geom_text(aes(label=ifelse(d_full$gene %in% stim_TF==T, as.character(gene),"")), hjust=1, vjust=0, size=3) +
  xlab("CD8 EMStim vs Naive") +
  ylab("CD8 EM vs Naive") +
  scale_x_continuous(limits=c(-5,5)) +
  scale_y_continuous(limits=c(-5,5)) +
  guides(size=F, color=F)
ggsave(filename=here(output_dir,"CD8_EMStim_EM_pgrnklofFC_stimTFs.pdf"), height=7, width=7)



p_CD8_CMStim_CM <- ggplot(d_full, aes(x=d.CD8_CMStim.d.CD8_N, y=d.CD8_CM.d.CD8_N))
p_CD8_CMStim_CM + geom_point(aes(size=CD8_CMStim_CM_mean, color=CD8_CMStim_CM_mean)) +
  geom_abline(slope=1, intercept=0, ) +
  geom_text(aes(label=ifelse(gene %in% top$CD8_CMStim_CM$gene, as.character(gene),"")), hjust=1, vjust=0, size=3) +
  xlab("CD8 CMStim vs Naive") +
  ylab("CD8 CM vs Naive") +
  scale_x_continuous(limits=c(-5,5)) +
  scale_y_continuous(limits=c(-5,5)) +
  guides(size=F, color=F)
ggsave(filename=here(output_dir,"CD8_CMStim_CM_pgrnklofFC_2.pdf"), height=7, width=7)

p_CD8_CMStim_CM + geom_point(aes(size=CD8_CMStim_CM_mean, color=CD8_CMStim_CM_mean)) +
  geom_abline(slope=1, intercept=0, ) +
  geom_text(aes(label=ifelse(d_full$gene %in% stim_TF==T, as.character(gene),"")), hjust=1, vjust=0, size=3) +
  xlab("CD8 CMStim vs Naive") +
  ylab("CD8 CM vs Naive") +
  scale_x_continuous(limits=c(-5,5)) +
  scale_y_continuous(limits=c(-5,5)) +
  guides(size=F, color=F)
ggsave(filename=here(output_dir,"CD8_CMStim_CM_pgrnklofFC_stimTFs.pdf"), height=7, width=7)

p_CD8_EMRAStim_EMRA <- ggplot(d_full, aes(x=d.CD8_EMRAStim.d.CD8_N, y=d.CD8_EMRA.d.CD8_N))
p_CD8_EMRAStim_EMRA + geom_point(aes(size=CD8_EMRAStim_EMRA_mean, color=CD8_EMRAStim_EMRA_mean)) +
  geom_abline(slope=1, intercept=0, ) +
  geom_text(aes(label=ifelse(gene %in% top$CD8_EMRAStim_EMRA$gene, as.character(gene),"")), hjust=1, vjust=0, size=3) +
  xlab("CD8 EMRAStim vs Naive") +
  ylab("CD8 EMRA vs Naive") +
  scale_x_continuous(limits=c(-5,5)) +
  scale_y_continuous(limits=c(-5,5)) +
  guides(size=F, color=F)
ggsave(filename=here(output_dir,"CD8_EMRAStim_EMRA_pgrnklofFC_2.pdf"), height=7, width=7)

p_CD8_EMRAStim_EMRA + geom_point(aes(size=CD8_EMRAStim_EMRA_mean, color=CD8_EMRAStim_EMRA_mean)) +
  geom_abline(slope=1, intercept=0, ) +
  geom_text(aes(label=ifelse(d_full$gene %in% stim_TF==T, as.character(gene),"")), hjust=1, vjust=0, size=3) +
  xlab("CD8 EMRAStim vs Naive") +
  ylab("CD8 EMRA vs Naive") +
  scale_x_continuous(limits=c(-5,5)) +
  scale_y_continuous(limits=c(-5,5)) +
  guides(size=F, color=F)
ggsave(filename=here(output_dir,"CD8_EMRAStim_EMRA_pgrnklofFC_stimTFs.pdf"), height=7, width=7)

p_CD4_EMStim_EM <- ggplot(d_full, aes(x=d.CD4_EMStim.d.CD4_N, y=d.CD4_EM.d.CD4_N))
p_CD4_EMStim_EM + geom_point(aes(size=CD4_EMStim_EM_mean, color=CD4_EMStim_EM_mean)) +
  geom_abline(slope=1, intercept=0, ) +
  geom_text(aes(label=ifelse(gene %in% top$CD4_EMStim_EM$gene, as.character(gene),"")), hjust=1, vjust=0, size=3) +
  xlab("CD4 EMStim vs Naive") +
  ylab("CD4 EM vs Naive") +
  scale_x_continuous(limits=c(-5,5)) +
  scale_y_continuous(limits=c(-5,5)) +
  guides(size=F, color=F)
ggsave(filename=here(output_dir,"CD4_EMStim_EM_pgrnklofFC_2.pdf"), height=7, width=7)

p_CD4_EMStim_EM + geom_point(aes(size=CD4_EMStim_EM_mean, color=CD4_EMStim_EM_mean)) +
  geom_abline(slope=1, intercept=0, ) +
  geom_text(aes(label=ifelse(d_full$gene %in% stim_TF==T, as.character(gene),"")), hjust=1, vjust=0, size=3) +
  xlab("CD4 EMStim vs Naive") +
  ylab("CD4 EM vs Naive") +
  scale_x_continuous(limits=c(-5,5)) +
  scale_y_continuous(limits=c(-5,5)) +
  guides(size=F, color=F)
ggsave(filename=here(output_dir,"CD4_EMStim_EM_pgrnklofFC_stimTFs.pdf"), height=7, width=7)

p_CD4_CMStim_CM <- ggplot(d_full, aes(x=d.CD4_CMStim.d.CD4_N, y=d.CD4_CM.d.CD4_N))
p_CD4_CMStim_CM + geom_point(aes(size=CD4_CMStim_CM_mean, color=CD4_CMStim_CM_mean)) +
  geom_abline(slope=1, intercept=0, ) +
  geom_text(aes(label=ifelse(gene %in% top$CD4_CMStim_CM$gene, as.character(gene),"")), hjust=1, vjust=0, size=3) +
  xlab("CD4 CMStim vs Naive") +
  ylab("CD4 CM vs Naive") +
  scale_x_continuous(limits=c(-5,5)) +
  scale_y_continuous(limits=c(-5,5)) +
  guides(size=F, color=F)
ggsave(filename=here(output_dir,"CD4_CMStim_CM_pgrnklofFC_2.pdf"), height=7, width=7)

p_CD4_CMStim_CM + geom_point(aes(size=CD4_CMStim_CM_mean, color=CD4_CMStim_CM_mean)) +
  geom_abline(slope=1, intercept=0, ) +
  geom_text(aes(label=ifelse(d_full$gene %in% stim_TF==T, as.character(gene),"")), hjust=1, vjust=0, size=3) +
  xlab("CD4 CMStim vs Naive") +
  ylab("CD4 CM vs Naive") +
  scale_x_continuous(limits=c(-5,5)) +
  scale_y_continuous(limits=c(-5,5)) +
  guides(size=F, color=F)
ggsave(filename=here(output_dir,"CD4_CMStim_CM_pgrnklofFC_stimTFs.pdf"), height=7, width=7)