
library(tidyverse)
library(ggpubr)
library(here)
library(triwise)

input_dir <- "naive_v_mem_overlap"
output_dir <- "naive_v_mem_overlap/triwise"

#Loading differentially expressed genes
diff <- read.delim(here(input_dir, "diff.significant.glm.Human_Tcell_RNA.rev.txt"), header=T, sep="\t")

#################################################################################################
#CD8 EMRA, EM, CM plot

#Creating matrix of mean rpkm counts
rpkm_CD8_EMRA <- cbind(diff[,c("SYMBOL")], diff[,grep("CD8_((CM)|(EM)|(EMRA)).rpkm", colnames(diff))])
colnames(rpkm_CD8_EMRA )[1] <- "SYMBOL"
mean_rpkm_CD8_EMRA  <- pivot_longer(rpkm_CD8_EMRA , 
                          cols=ends_with('.rpkm'),
                          names_to=c("Sample", "Celltype", "Subtype"),
                          names_pattern= "RNA_(.*)_(.*)_(.*).rpkm",
                          values_to="rpkm"
                          ) %>%
  group_by(SYMBOL, Subtype) %>%
  summarise(mean_rpkm=mean(rpkm)) %>%
  pivot_wider(names_from=Subtype, values_from=mean_rpkm)

#Creating matrix of sig (T/F)
sig_CD8 <- cbind(diff[c("SYMBOL")], diff[,grep("CD8_((CM_)|(EM_)|(EMRA_))blood_unstim.v.CD8_Nav_blood_unstim.sig", colnames(diff))]) %>% 
  drop_na(SYMBOL) %>% 
  drop_na(ends_with("sig"))
colnames(sig_CD8) <- gsub("_blood_unstim.v.CD8_Nav_blood_unstim.sig", "sig", colnames(sig_CD8))
sig_CD8$bin <- rowSums(sig_CD8[,-1])
sig_CD8_sig <- subset(sig_CD8, bin!=0)

sig_CD8_omni <- cbind(diff[c("SYMBOL")], diff[,grep("CD8_((CM_)|(EM_)|(EMRA_))blood_unstim.v.CD8_((CM_)|(EM_)|(EMRA_)|(Nav_))blood_unstim.sig", colnames(diff))]) %>% 
  drop_na(SYMBOL)  %>% 
  drop_na(ends_with("sig"))
#colnames(sig_CD8) <- gsub("_blood_unstim.v.CD8_((CM_)|(EM_)|(EMRA_)|(Nav_))blood_unstim.sig", "sig", colnames(sig_CD8))
sig_CD8_omni$bin <- rowSums(sig_CD8_omni[,-1])
sig_CD8_omnisig <- subset(sig_CD8_omni, bin!=0)

#Various combos
sig_CD8_EMEMRA <- subset(sig_CD8, CD8_CMsig==F & CD8_EMsig==T & CD8_EMRAsig==T)
sig_CD8_all <- subset(sig_CD8, CD8_CMsig==T & CD8_EMsig==T & CD8_EMRAsig==T )
sig_CD8_CMEM <- subset(sig_CD8, CD8_CMsig==T & CD8_EMsig==T & CD8_EMRAsig==F )
sig_CD8_EM <- subset(sig_CD8, CD8_CMsig==F & CD8_EMsig==T & CD8_EMRAsig==F )
sig_CD8_EMRA <- subset(sig_CD8, CD8_CMsig==F & CD8_EMsig==F & CD8_EMRAsig==T )
sig_CD8_CM <- subset(sig_CD8, CD8_CMsig==T & CD8_EMsig==F & CD8_EMRAsig==F )

#Matrix for triwise
mean_rpkm_mtx_CD8_EMRA<- as.matrix(mean_rpkm_CD8_EMRA[,-1])
rownames(mean_rpkm_mtx_CD8_EMRA) <- mean_rpkm_CD8_EMRA$SYMBOL
mean_rpkm_mtx_CD8_EMRA2 <- mean_rpkm_mtx_CD8_EMRA[rowMeans(mean_rpkm_mtx_CD8_EMRA)>3,]
mean_rpkm_mtx_CD8_EMRAsub <- mean_rpkm_mtx_CD8_EMRA[rownames(mean_rpkm_mtx_CD8_EMRA)%in%sig_CD8_sig$SYMBOL,]
#mean_rpkm_mtx <- mean_rpkm_mtx[order(match(rownames(mean_rpkm_mtx), sig_CD8_EM$SYMBOL)),]

#Log2 transformations
mean_rpkm_mtx_CD8_EMRA <- apply(mean_rpkm_mtx_CD8_EMRA, c(1,2), function(x) log(x, base=2))
mean_rpkm_mtx_CD8_EMRA2 <- apply(mean_rpkm_mtx_CD8_EMRA2, c(1,2), function(x) log(x, base=2))
mean_rpkm_mtx_CD8_EMRAsub <- apply(mean_rpkm_mtx_CD8_EMRAsub, c(1,2), function(x) log(x, base=2))

#Barycoordinate transformation
barycoords_CD8_EMRA = transformBarycentric(mean_rpkm_mtx_CD8_EMRA) %>% drop_na()
barycoords_CD8_EMRAsub = transformBarycentric(mean_rpkm_mtx_CD8_EMRAsub) %>% drop_na()
barycoords_CD8_EMRA2 = transformBarycentric(mean_rpkm_mtx_CD8_EMRA2) %>% drop_na()

plotDotplot(barycoords_CD8_EMRAsub,
            #rmax=12,
            Gdiffexp=sig_CD8_sig$SYMBOL,
            #colorby="diffexp", 
            sizevalues = stats::setNames(c(0.5, 2), c(FALSE, TRUE)),
            alphavalues = stats::setNames(c(0.3, 0.8), c(FALSE, TRUE)),
            list(diffexpall="#000000", nodiffexpall="#AAAAAA", nodiffexpgset="#FFAAAA", diffexpgset="#FF0000")
            )

plotRoseplot(barycoords_CD8_EMRA,relative = T, rmax=1)

pCD8_EMRA_bary  <- plotDotplot(barycoords_CD8_EMRA2, showlabels = F) + guides(color=F, size=F, alpha=F)
pCD8_EMRA_rose <- plotRoseplot(barycoords_CD8_EMRA2, relative = T)

#################################################################################################
#CD8 Nav, EM, CM plot

#Creating matrix of mean rpkm counts
rpkm_CD8 <- cbind(diff[,c("SYMBOL")], diff[,grep("CD8_((CM)|(EM)|(N)).rpkm", colnames(diff))])
colnames(rpkm_CD8 )[1] <- "SYMBOL"
mean_rpkm_CD8  <- pivot_longer(rpkm_CD8 , 
                                    cols=ends_with('.rpkm'),
                                    names_to=c("Sample", "Celltype", "Subtype"),
                                    names_pattern= "RNA_(.*)_(.*)_(.*).rpkm",
                                    values_to="rpkm"
) %>%
  group_by(SYMBOL, Subtype) %>%
  summarise(mean_rpkm=mean(rpkm)) %>%
  pivot_wider(names_from=Subtype, values_from=mean_rpkm)

# #Creating matrix of sig (T/F)
# sig_CD8 <- cbind(diff[c("SYMBOL")], diff[,grep("CD8_((CM_)|(EM_)|(EMRA_))blood_unstim.v.CD8_Nav_blood_unstim.sig", colnames(diff))]) %>% drop_na(SYMBOL)
# colnames(sig_CD8) <- gsub("_blood_unstim.v.CD8_Nav_blood_unstim.sig", "sig", colnames(sig_CD8))
# sig_CD8$bin <- rowSums(sig_CD8[,-1])
# sig_CD8_sig <- subset(sig_CD8, bin!=0)
# #Various combos
# sig_CD8_EMEMRA <- subset(sig_CD8, CD8_CMsig==F & CD8_EMsig==T & CD8_EMRAsig==T)
# sig_CD8_all <- subset(sig_CD8, CD8_CMsig==T & CD8_EMsig==T & CD8_EMRAsig==T )
# sig_CD8_CMEM <- subset(sig_CD8, CD8_CMsig==T & CD8_EMsig==T & CD8_EMRAsig==F )
# sig_CD8_EM <- subset(sig_CD8, CD8_CMsig==F & CD8_EMsig==T & CD8_EMRAsig==F )
# sig_CD8_EMRA <- subset(sig_CD8, CD8_CMsig==F & CD8_EMsig==F & CD8_EMRAsig==T )
# sig_CD8_CM <- subset(sig_CD8, CD8_CMsig==T & CD8_EMsig==F & CD8_EMRAsig==F )

sig_CD8_noEMRA <- cbind(diff[c("SYMBOL")], diff[,grep("CD8_((CM_)|(EM_))blood_unstim.v.CD8_Nav_blood_unstim.sig", colnames(diff))]) %>% 
  drop_na(SYMBOL)  %>% 
  drop_na(ends_with("sig"))
colnames(sig_CD8_noEMRA) <- gsub("_blood_unstim.v.CD8_Nav_blood_unstim.sig", "sig", colnames(sig_CD8_noEMRA))
sig_CD8_noEMRA$bin <- rowSums(sig_CD8_noEMRA[,-1])
sig_CD8_noEMRA_sig <- subset(sig_CD8_noEMRA, bin!=0)

#Matrix for triwise
mean_rpkm_mtx_CD8<- as.matrix(mean_rpkm_CD8[,-1])
rownames(mean_rpkm_mtx_CD8) <- mean_rpkm_CD8$SYMBOL
mean_rpkm_mtx_CD8sub <- mean_rpkm_mtx_CD8[rownames(mean_rpkm_mtx_CD8)%in%sig_CD8_noEMRA_sig$SYMBOL,]
#mean_rpkm_mtx <- mean_rpkm_mtx[order(match(rownames(mean_rpkm_mtx), sig_CD8_EM$SYMBOL)),]

#Log2 transform
mean_rpkm_mtx_CD8 <- apply(mean_rpkm_mtx_CD8, c(1,2), function(x) log(x, base=2))
mean_rpkm_mtx_CD8sub <- apply(mean_rpkm_mtx_CD8sub, c(1,2), function(x) log(x, base=2))

barycoords_CD8 = transformBarycentric(mean_rpkm_mtx_CD8) %>% drop_na()
barycoords_CD8sub = transformBarycentric(mean_rpkm_mtx_CD8sub) %>% drop_na()
 
plotDotplot(barycoords_CD8,
            rmax=8,
            Gdiffexp=sig_CD8_noEMRA_sig$SYMBOL,
            #colorby="diffexp", 
            list(diffexpall="#000000", nodiffexpall="#AAAAAA", nodiffexpgset="#FFAAAA", diffexpgset="#FF0000")
) #+ geom_label(aes(label=ifelse(gene=="LEF1", as.character(gene), "")))

plotRoseplot(barycoords_CD8_EMRA2)

pCD8_bary  <- plotDotplot(barycoords_CD8sub, rmax=8, showlabels = F) + guides(color=F, size=F, alpha=F)
pCD8_rose <- plotRoseplot(barycoords_CD8, Gdiffexp=sig_CD8_noEMRA_sig$SYMBOL, relative = T)


#################################################################################################
#CD4
rpkm_CD4 <- cbind(diff[,c("SYMBOL")], diff[,grep("CD4_((CM)|(EM)|(N)).rpkm", colnames(diff))])
colnames(rpkm_CD4)[1] <- "SYMBOL"
mean_rpkm_CD4 <- pivot_longer(rpkm_CD4, 
                          cols=ends_with('.rpkm'),
                          names_to=c("Sample", "Celltype", "Subtype"),
                          names_pattern= "RNA_(.*)_(.*)_(.*).rpkm",
                          values_to="rpkm"
) %>%
  group_by(SYMBOL, Subtype) %>%
  summarise(mean_rpkm=mean(rpkm)) %>%
  pivot_wider(names_from=Subtype, values_from=mean_rpkm)

#Creating matrix of sig (T/F)
sig_CD4 <- cbind(diff[c("SYMBOL")], diff[,grep("CD4_((CM_)|(EM_))blood_unstim.v.CD4_Nav_blood_unstim.sig", colnames(diff))]) %>% 
  drop_na(SYMBOL)  %>% 
  drop_na(ends_with("sig"))
colnames(sig_CD4) <- gsub("_blood_unstim.v.CD4_Nav_blood_unstim.sig", "sig", colnames(sig_CD4))
sig_CD4$bin <- rowSums(sig_CD4[,-1])
sig_CD4_sig <- subset(sig_CD4, bin!=0)

#Matrix for triwise
mean_rpkm_mtx_CD4 <- as.matrix(mean_rpkm_CD4[,-1])
rownames(mean_rpkm_mtx_CD4) <- mean_rpkm_CD4$SYMBOL
mean_rpkm_mtx_CD4sub <- mean_rpkm_mtx_CD4[rownames(mean_rpkm_mtx_CD4)%in%sig_CD4_sig$SYMBOL,]

#Log2 transform
mean_rpkm_mtx_CD4 <- apply(mean_rpkm_mtx_CD4, c(1,2), function(x) log(x, base=2))
mean_rpkm_mtx_CD4sub <- apply(mean_rpkm_mtx_CD4sub, c(1,2), function(x) log(x, base=2))

barycoords_CD4 = transformBarycentric(mean_rpkm_mtx_CD4) %>% drop_na()
barycoords_CD4sub = transformBarycentric(mean_rpkm_mtx_CD4sub) %>% drop_na()

plotDotplot(barycoords_CD4,
            rmax=8,
            Gdiffexp=sig_CD4_sig$SYMBOL,
            #colorby="diffexp", 
            list(diffexpall="#000000", nodiffexpall="#AAAAAA", nodiffexpgset="#FFAAAA", diffexpgset="#FF0000")
)
plotRoseplot(barycoords_CD4,
             rmax=1,
             #Gdiffexp=sig_CD4_sig$SYMBOL
              )  
pCD4_bary <- plotDotplot(barycoords_CD4sub, rmax=8, showlabels = F) + guides(color=F, size=F, alpha=F)
pCD4_rose <- plotRoseplot(barycoords_CD4sub, relative = T)

#################################################################################################
#CD8 Nav, EM, EMRA plot

#Creating matrix of mean rpkm counts
rpkm_CD8_noCM <- cbind(diff[,c("SYMBOL")], diff[,grep("CD8_((EM)|(EMRA)|(N)).rpkm", colnames(diff))])
colnames(rpkm_CD8_noCM)[1] <- "SYMBOL"
mean_rpkm_CD8_noCM  <- pivot_longer(rpkm_CD8_noCM, 
                               cols=ends_with('.rpkm'),
                               names_to=c("Sample", "Celltype", "Subtype"),
                               names_pattern= "RNA_(.*)_(.*)_(.*).rpkm",
                               values_to="rpkm"
) %>%
  group_by(SYMBOL, Subtype) %>%
  summarise(mean_rpkm=mean(rpkm)) %>%
  pivot_wider(names_from=Subtype, values_from=mean_rpkm)

sig_CD8_noCM <- cbind(diff[c("SYMBOL")], diff[,grep("CD8_((EM_)|(EMRA_))blood_unstim.v.CD8_Nav_blood_unstim.sig", colnames(diff))]) %>% 
  drop_na(SYMBOL) %>% 
  drop_na(ends_with("sig"))
colnames(sig_CD8_noCM) <- gsub("_blood_unstim.v.CD8_Nav_blood_unstim.sig", "sig", colnames(sig_CD8_noCM))
sig_CD8_noCM$bin <- rowSums(sig_CD8_noCM[,-1])
sig_CD8_noCM_sig <- subset(sig_CD8_noCM, bin!=0)

#Matrix for triwise
mean_rpkm_mtx_CD8_noCM<- as.matrix(mean_rpkm_CD8_noCM[,-1])
rownames(mean_rpkm_mtx_CD8_noCM) <- mean_rpkm_CD8_noCM$SYMBOL
mean_rpkm_mtx_CD8_noCMsub <- mean_rpkm_mtx_CD8_noCM[rownames(mean_rpkm_mtx_CD8_noCM)%in%sig_CD8_noCM_sig$SYMBOL,]
#mean_rpkm_mtx <- mean_rpkm_mtx[order(match(rownames(mean_rpkm_mtx), sig_CD8_EM$SYMBOL)),]

#Log2 trans
mean_rpkm_mtx_CD8_noCM <- apply(mean_rpkm_mtx_CD8_noCM, c(1,2), function(x) log(x, base=2))
mean_rpkm_mtx_CD8_noCMsub <- apply(mean_rpkm_mtx_CD8_noCMsub, c(1,2), function(x) log(x, base=2))

barycoords_CD8_noCM = transformBarycentric(mean_rpkm_mtx_CD8_noCM) %>% drop_na()
barycoords_CD8_noCMsub = transformBarycentric(mean_rpkm_mtx_CD8_noCMsub) %>% drop_na()

plotDotplot(barycoords_CD8_noCM,
            rmax=8,
            Gdiffexp=sig_CD8_noCM_sig$SYMBOL,
            #colorby="diffexp", 
            list(diffexpall="#000000", nodiffexpall="#AAAAAA", nodiffexpgset="#FFAAAA", diffexpgset="#FF0000")
)

plotRoseplot(barycoords_CD8_noCMsub, relative = F)

pCD8_noCM_bary  <- plotDotplot(barycoords_CD8_noCMsub, rmax=8, showlabels = F) + guides(color=F, size=F, alpha=F)
pCD8_noCM_rose <- plotRoseplot(barycoords_CD8_noCMsub, relative = T)




#################################################################################################
#Plotting and saving

pCD8_bary
ggsave(filename=here(output_dir,"CD8_triwise_dot.pdf"), width=3.5, height=3)

pCD8_rose
ggsave(filename=here(output_dir,"CD8_triwise_rose.pdf"), width=3.5, height=3)

pCD4_bary
ggsave(filename=here(output_dir,"CD4_triwise_dot.pdf"), width=3.5, height=3)

pCD4_rose
ggsave(filename=here(output_dir,"CD4_triwise_rose.pdf"), width=3.5, height=3)

ggarrange(pCD8_bary, pCD4_bary)
ggsave(filename=here(output_dir,"CD8_CD4_triwise.pdf"), width=7, height=3)

pCD8_EMRA_bary
ggsave(filename=here(output_dir,"CD8_EMRA_triwise_dot.pdf"), width=5.5, height=5)

pCD8_EMRA_rose
ggsave(filename=here(output_dir,"CD8_EMRA_triwise_rose.pdf"), width=5.5, height=5)

ggarrange(pCD8_bary, pCD8_noCM_bary, pCD8_EMRA_bary, nrow=1)

ggarrange(pCD8_bary, pCD4_bary,pCD8_EMRA_bary, pCD8_rose, pCD4_rose,pCD8_EMRA_rose,
         labels=c("CD8 Nav", "CD4 Nav", "CD8 EMRA", "", "", ""))
ggsave(filename=here(output_dir,"CD8_CD4_EMRA_triwise.pdf"), width=10.5, height=9)

ggarrange(pCD8_bary, pCD8_noCM_bary, pCD8_EMRA_bary, pCD8_rose, pCD8_noCM_rose, pCD8_EMRA_rose)

#################################################################################################
#Plots by genesets

genes_df <- read.csv(here("topgenes.csv"), header=T)
genes_df <- genes_df[,-1]


#barycoords_CD8_EMRA_geneset <- barycoords_CD8_EMRA[rownames(barycoords_CD8_EMRA)%in%geneset,]

dotplot = interactiveDotplot(mean_rpkm_mtx_CD8_EMRA2, 
                             Gdiffexp=sig_CD8_omnisig$SYMBOL,
                            #Goi=geneset,
                            #rmax=12
                            )
dotplot

geneset_name <- "Migration"
geneset=select(genes_df, !!geneset_name) %>% unlist() %>% unique()%>%as.character()

plotDotplot(barycoords_CD8_EMRA2,
            #rmax=12,
            Gdiffexp=sig_CD8_omnisig$SYMBOL,
            Goi=list(geneset=geneset),
            #colorby="diffexp", 
            sizevalues = stats::setNames(c(0.5, 2), c(FALSE, TRUE)),
            alphavalues = stats::setNames(c(0.3, 0.8), c(FALSE, TRUE)),
            showlabels = F,
            colorvalues=list(diffall="black", diffgeneset="red", NAgeneset ="black", nodiffall="grey", nodiffgeneset="pink")
) + guides(color=F, size=F, alpha=F) + ggtitle(label=geneset_name)
ggsave(filename=here(output_dir, paste0("CD8_unstim_",geneset_name, "_triwise.pdf")), height=4, width=4)

plotDotplot(barycoords_CD8sub,
            #rmax=12,
            Gdiffexp=sig_CD8_sig$SYMBOL,
            Goi=list(geneset=geneset),
            #colorby="diffexp", 
            sizevalues = stats::setNames(c(0.5, 2), c(FALSE, TRUE)),
            alphavalues = stats::setNames(c(0.3, 0.8), c(FALSE, TRUE)),
            showlabels = F,
            colorvalues=list(diffall="black", diffgeneset="red", NAgeneset ="black", nodiffall="grey", nodiffgeneset="pink")
)

plotDotplot(barycoords_CD8_noCMsub,
            #rmax=12,
            Gdiffexp=sig_CD8_sig$SYMBOL,
            Goi=list(geneset=geneset),
            #colorby="diffexp", 
            sizevalues = stats::setNames(c(0.5, 2), c(FALSE, TRUE)),
            alphavalues = stats::setNames(c(0.3, 0.8), c(FALSE, TRUE)),
            showlabels = F,
            colorvalues=list(diffall="black", diffgeneset="red", NAgeneset ="black", nodiffall="grey", nodiffgeneset="pink")
)
