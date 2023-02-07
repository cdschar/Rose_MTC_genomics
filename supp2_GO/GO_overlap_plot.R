library(here)
library(tidyverse)
library(ggrepel)

input_dir <- 'naive_v_mem_overlap/CD4CD8_comp/GO_output'
output_dir <- 'naive_v_mem_overlap/CD4CD8_comp/GO_output'
files <- c('CM_CD4only_GO.txt', 'CM_CD8only_GO.txt', 'CMshared_GO.txt')

######################################################################################
#Read in list of currated pathways
pathways <- read.delim(file=here(input_dir, 'pathways_CM.txt'), sep="\t", header=F, 
                       colClasses=c("factor", "character"),
                       col.names=c("Pathway", "Access")
                       )

######################################################################################
#Read in GO output files into list and fix issues
GOlst <- list()
for (i in 1:length(files)){
  GOlst[[i]] <- read.delim(file=here(input_dir,files[i]), header = T)
  names(GOlst[[i]]) <- gsub("Client.Text.Box.Input..", "", names(GOlst[[i]]))
  GOlst[[i]]$fold.Enrichment. <- as.numeric(GOlst[[i]]$fold.Enrichment.)
  GOlst[[i]] <- separate(GOlst[[i]], GO.biological.process.complete, c('Pathway', 'Access'), '\\(GO:')
  GOlst[[i]]$Access <- gsub(")", "", GOlst[[i]]$Access)
}

names(GOlst) <- gsub('_GO.txt', '', files)
for (i in 1:length(GOlst)){
  GOlst[[i]] <- mutate(GOlst[[i]], comp=names(GOlst)[i])
}

######################################################################################
#Adding log trans of FDR and fold enrichment
addvars <- function(x){
  return(mutate(x, neglog10FDR=-log(FDR., base=10), log2Fenrich=log(fold.Enrichment., base=2)))
}

GOlst <- map(GOlst, addvars)

######################################################################################
#Ploting top pathways

library(plotly)

#ggplot(subset(GOlst[[1]], Homo.sapiens...REFLIST..20851.<=50), aes(y=neglog10FDR, x=fold.Enrichment.)) + geom_point() + geom_text_repel(aes(label=if_else(fold.Enrichment.>100 & neglog10FDR >=0.75, as.character(Pathway), "")))
#ggplot(subset(GOlst[[2]], Homo.sapiens...REFLIST..20851.<=50), aes(y=neglog10FDR, x=fold.Enrichment.)) + geom_point() + geom_text_repel(aes(label=if_else(fold.Enrichment.>1 & neglog10FDR >=1.3, as.character(Pathway), "")))
#ggplot(subset(GOlst[[3]], Homo.sapiens...REFLIST..20851.<=50), aes(y=neglog10FDR, x=fold.Enrichment.)) + geom_point() + geom_text_repel(aes(label=if_else(fold.Enrichment.<300 & neglog10FDR >=2, as.character(Pathway), "")))
    #+ geom_text_repel(aes(label=if_else(log2Fenrich>7 & neglog10FDR >=6, as.character(Pathway), "")))

plot_ly(data=subset(GOlst[[1]], Homo.sapiens...REFLIST..20851.<=100), 
        x=~fold.Enrichment., 
        y= ~neglog10FDR,
        text = ~Pathway)

plot_ly(data=subset(GOlst[[2]], Homo.sapiens...REFLIST..20851.<=100), 
        x=~fold.Enrichment., 
        y= ~neglog10FDR,
        text = ~Pathway)

plot_ly(data=subset(GOlst[[3]], Homo.sapiens...REFLIST..20851.<=100), 
        x=~fold.Enrichment., 
        y= ~neglog10FDR,
        text = ~Pathway)

plot_ly(data=subset(GOlst[[3]], Homo.sapiens...REFLIST..20851.<=100), 
        x=~Homo.sapiens...REFLIST..20851., 
        y= ~neglog10FDR,
        text = ~Pathway)

######################################################################################
#Subsetting to just currated pathawys and combining to one df

sub_path <- function(x){
  y <- x[x$Access%in%pathways$Access,]
  return(y)
}

GOlst_sub <- map(GOlst, sub_path)

GOdf <- data.frame()
for (i in 1:length(GOlst_sub)){
  GOdf <- bind_rows(GOdf, GOlst_sub[[i]])
}

######################################################################################
#Plot currated pathways

# p1 <- ggplot(GOdf, aes(x=comp, y=Pathway)) +
#     geom_point(aes(color=neglog10FDR, size=log2Fenrich)) +
#     theme_light()
# p1

#GOdf <- GOdf %>% arrange(neglog10FDR)
GOdf$Pathway <- substr(GOdf$Pathway, 1, nchar(GOdf$Pathway)-1)
GOdf$Pathway <- factor(GOdf$Pathway, levels=pathways$Pathway, ordered = T)

p2 <- ggplot(GOdf, aes(x=Pathway, y=neglog10FDR, group=comp)) +
    geom_col(aes(fill=comp),position=position_dodge(width = 0.8)) +
    coord_flip() +
    scale_x_discrete(limits = rev(levels(GOdf$Pathway)))+
    theme_light() +
    ylab("FDR (-log10)") +
    labs(fill="Overlap") +
    geom_hline(yintercept=1.3, linetype="dashed")
p2
ggsave(here(output_dir, "CM_GO_overlap.pdf"), width=7, height=4)

######################################################################################
#Biggest difference in pvalue

# GObigdf <- data.frame()
# for (i in 1:length(GOlst)){
#   GObigdf <- bind_rows(GObigdf, GOlst[[i]])
# }

