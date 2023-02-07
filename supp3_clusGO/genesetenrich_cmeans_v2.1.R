#Script for performing fisher exact tests on membership in gene cluster and genesets
#v2 is the same as v1, but can handle larger genesets
#v2.1 update added gene names to the dataframe for reference 
#jrose
#10Sep2020

library(tidyverse)
library(qusage)
library(magrittr)
library(here)

geneset_name <- "Tcell_curated"
input_name <- "CD8_unstim_labfilter"

#Loading results of clustering
input_dir <- "analysis/kmeans"
input_file <- "CD8_unstim_labfilter_RNA_fuzzmeans_clusterdata_v2.txt"
output_dir <- 'analysis/kmeans/enrich'
dataDiff <- read.table(file=here(input_dir, input_file), header=T)

#loading geneset as list object
geneset_dir <- 'analysis/unstim/GSEA/genesets'
#genesets <- read.gmt(here(geneset_dir, 'hallmark.gmt'))
genesets <- read.gmt(here(geneset_dir, 'Tcell_curated.gmt'))

#Creating contingency tables by cluster and geneset
#Custom definied function to: '
  # - Group by cluster identity
  # - Create tables of inclusion (yes/no) for each cluster
conting <- function(x){  
  dataDiff %>% group_by(k.5.cluster) %>% 
  summarise(no=table(SYMBOL%in%x)[1], yes=table(SYMBOL%in%x)[2], genes=list(x[x%in%SYMBOL])) %>%
  replace_na(list(yes=0, no=0))
}
#Mapping conting function through geneset list
suppressMessages(
  tables <- map(genesets, conting)
) #group_by is throwing a nuicanse warning


#Function creating 2x2 contingency tables and performing fisher test
#Meant to be mapped through tables list above on contingency tables
fisher_testing <- function(x){
  tests <- list()
  
  for (i in 1:dim(x)[1]){
    tmp <- matrix(c(x[i,"no"], x[i, "yes"], sum(x[-i,"no"]), sum(x[-i,"yes"])), nrow=2, ncol=2, byrow=T)
    tmp <- apply(tmp, c(1,2), as.integer)
    colnames(tmp) <- c("No", "Yes")
    rownames(tmp) <- c("In Clus", "Out Clus")
    
    tests[[i]] <- fisher.test(tmp)
  }
  
  #Extracting pvalues and estimates from tests
  y <- tests %>% {
    tibble(
      pvalue = map_dbl(., "p.value"),
      #conf.int = map_dbl(., "conf.int"),
      estimate = map_dbl(., "estimate")
          )
  }
  
  #Combining test statistics with tables table (ha!)
  return(cbind(x, y))
      
}
tables_tests <- map(tables, fisher_testing)

for (n in names(tables_tests)){
  tables_tests[[n]]['geneset'] = n
}

#Extract list to dataframe
df <- map_dfr(tables_tests, extract, c("k.5.cluster", "pvalue", "no", "yes", "genes", "estimate", "null.value", "geneset"))

#Remake list of dataframes by cluster
df_clus_lst <- list()
for (c in unique(df$k.5.cluster)){
  df_clus_lst[[c]] = subset(df, k.5.cluster==c)
}

#Polish up each dataframe to pull out the most significant pathways in each  
polish <- function(x){
  tmp <- x %>% subset(pvalue<=0.05) %>% 
    mutate(neglog10P = -log(pvalue, base=10)) %>%
    arrange(neglog10P) %>%
    top_n(50, neglog10P)
  tmp$geneset <- factor(tmp$geneset, levels=unique(tmp$geneset), ordered=T)
  return(tmp)
}
df_clus_polish <- map(df_clus_lst, polish)

plots <- list()

for (i in 1:length(df_clus_polish)){
  plots[[i]] <- ggplot(df_clus_polish[[i]], aes(x=neglog10P, y=geneset)) +
            geom_col() +
            theme(axis.text.y = element_text(size=5))
  ggsave(here(output_dir, paste0(input_name,"_", geneset_name,"_cluster", i, "fisher.pdf")))
}

save(df_clus_polish, df_clus_lst, file=here(output_dir, paste0(input_name,"_", geneset_name, "_data.rda")))

# ggplot(df2, aes(x=neglog10P, y=geneset)) +
#   geom_col() +
#   facet_wrap(vars(k.5.cluster))
# ggsave(here(output_dir, "CD8_unstim_labfilter_Hallmark_fisher.pdf"))