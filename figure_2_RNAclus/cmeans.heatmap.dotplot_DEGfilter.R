#Function looks at optimal Fuzzy c means clustering groups
#jrose 13Jul2020
#Unstimulated CD8+ samples only

setwd("/Volumes/ESB/Boss_Projects/Human_Tcell_subsets/RNA_combined/analysis/kmeans/")

library("data.table")
library("GenomicRanges")
library("som")
library('tidyverse')
library('factoextra')
library('e1071')


##########################################################
#Read in data and wrangle

#read in the manifest file
fqDir = "/Volumes/ESB/Boss_Projects/Human_Tcell_subsets/RNA_combined/pipeline/";
fqFile = "RNAseq.sample.manifest.txt";
files = read.table(paste0(fqDir, fqFile), sep = "\t", header = T, as.is = T);
files = files[grepl("EM|CM|Nav", files$group),]
files = files[!grepl("CD4", files$group),]
files = files[which(files$include),]
files = files[files$stim=="unstim",]
files = mutate(files, cellSubtype = paste(cellType, cellSubset, stim, sep="_"))

#Prefix for output files
subset = "CD8_unstim_DEGfilter"

#Reorder based on prespecified group order
smpl_order = c("CD8_Nav_unstim", "CD8_CM_unstim", "CD8_EM_unstim", "CD8_EMRA_unstim")
files = files[order(match(files$cellSubtype, smpl_order)),]

#read in diff file
dataDir = "/Volumes/ESB/Boss_Projects/Human_Tcell_subsets/RNA_combined/analysis/diff/";
dataFile = "diff.significant.glm.Human_Tcell_RNA.txt";
diffRpm = read.table(file=paste0(dataDir, dataFile), header=T, sep="\t", comment.char="", quote="");

#Subset just significant in desired group comparisons
grps <- unique(files$group)
smpl_comp <- apply(expand.grid(grps, grps), 1, paste, collapse=".v.")
smpl_comp_sig <- paste(smpl_comp, ".sig", sep="")
smpl_comp <- c(paste(smpl_comp, ".logFC", sep=""),
               paste(smpl_comp, ".logCPM", sep=""),
               paste(smpl_comp, ".PValue", sep=""),
               paste(smpl_comp, ".fdr", sep=""),
               paste(smpl_comp, ".sig", sep="")
                )
keep <-diffRpm[apply(diffRpm[,names(diffRpm)%in%smpl_comp_sig],1,function(x) any(x==TRUE)),'SYMBOL']
diffRpm = diffRpm[diffRpm$SYMBOL%in%keep,]

#log transform data
#data = log(diffRpm[,paste(files$sample, ".rpkm", sep="")]+1, 2)
#range(data);

#Matrix of counts
data <- diffRpm[,paste(files$sample, ".rpkm", sep="")]

#Dataframe for saving
data_save <- select(diffRpm, ENTREZID, SYMBOL, sigAny) %>% cbind(data, diffRpm[,names(diffRpm)%in%smpl_comp])

#Exclude genes with lower average expression than threshold
#expr_test <- rowMeans(data)
#data <- data[expr_test>=5,]
#data_save <- data_save[expr_test>=5,]

#Normalize rowwise
data[1:length(files$sample)] = som::normalize(data)

#Remove NAs
rownames(data) <- diffRpm$ENTREZID
data <- drop_na(data)
data_save <- data_save[data_save$ENTREZID%in%rownames(data),]

##########################################################
#Determine fuzzifiern parameter (follows this paper: https://www.ncbi.nlm.nih.gov/pubmed/20880957)
mestimate<- function(df){
  N <-  dim(df)[[1]]
  D <- dim(df)[[2]]
  m.sj <- 1 + (1418/N + 22.05)*D^(-2) + (12.33/N +0.243)*D^(-0.0406*log(N) - 0.1134)
  return(m.sj)
}

m <- mestimate(data)
m

##########################################################
#Determine number of k clusters (must look at output and change krange variable later!)

#helper function for the within sum of squared error
sumsqr <- function(x, clusters){
  sumsqr <- function(x) sum(scale(x, scale = FALSE)^2)
  wss <- sapply(split(as.data.frame(x), clusters), sumsqr)
  return(wss)
}

#get the wss for repeated clustering
iterate_fcm_WSS <- function(df,m){
  totss <- numeric()
  for (i in 2:20){
    FCMresults <- cmeans(df,centers=i,m=m)
    totss[i] <- sum(sumsqr(df,FCMresults$cluster))
  }
  return(totss)
}
wss_2to20 <- iterate_fcm_WSS(data,m)

cairo_pdf(paste0(subset,"_fuzzmeans.optimization_nolog.pdf"))
plot(1:20, wss_2to20[1:20], type="b", xlab="Number of Clusters", ylab="wss")
dev.off()

##########################################################
#Actual clustering

#plot heatmaps from ideal k range
source("/Volumes/GRAID/seqTools/R.Libraries/heatmap3/heatmap3.R")

#Heatmap parameters
height = 1.5;
width = 1.5;
res = 1200;
scale = NA; #scale for heatmap
norm = NA;

#Colors
nav.col = "#5F7D8BFF"
cm.col =  "#0C46A0FF"
em.col = "#E55100FF"
emra.col ="#870D4EFF"
cols = c(rep(nav.col, 4),rep(cm.col, 4),rep(em.col, 4),rep(emra.col, 4));
cellcolors = c("CD8_Nav_unstim"=nav.col,
               "CD8_CM_unstim"=cm.col,
               "CD8_EM_unstim"=em.col,
               "CD8_EMRA_unstim"=emra.col)

rampColsDiv = c(rgb(31, 73, 125, maxColorValue = 255), rgb(1,1,1), rgb(149, 55, 53, maxColorValue = 255));
#rampCols = c(rgb(1,1,1), rgb(149, 55, 53, maxColorValue = 255));
#rampCols = c("white","orange","dark red");

iter.max = 100;
krange = 5:8
fms = NULL
fms_centroids = NULL
fms_centroids_cor = NULL
fms_centroids_long = NULL
fms_centroids_sum = NULL

colnum = dim(data)[2]

theme_set(theme_light()+ theme(axis.line = element_line(size=0.5),
                               axis.text = element_blank(),
                               axis.title = element_blank(),
                               legend.text = element_text(size=7),
                               legend.title = element_text(size=8)
                               )
          )

for (k in krange) {
  
  print(paste("K", k, Sys.time()))
  set.seed(1)
  fms[[k]] = cmeans(as.matrix(data), iter.max = iter.max, centers = k, m=m)
  
  #Looking at centroids
  fms_centroids[[k]] <- data.frame(fms[[k]]$centers)
  fms_centroids[[k]]$cluster <- row.names(fms_centroids[[k]])
  #How correlated are each of the centroids?
  fms_centroids_cor[[k]] <- cor(t(fms_centroids[[k]][1:colnum]))
  print(fms_centroids_cor[[k]])
  
  #Tyding up for average plot
  fms_centroids_long[[k]] <- pivot_longer(fms_centroids[[k]], 1:colnum, names_to="sample")
  fms_centroids_long[[k]]$sample <- gsub(".rpkm", "", fms_centroids_long[[k]]$sample)
  fms_centroids_long[[k]] <- left_join(fms_centroids_long[[k]],  select(files, sample, cellSubtype), by="sample")
  fms_centroids_long[[k]]$cellSubtype <- factor(fms_centroids_long[[k]]$cellSubtype, levels=smpl_order, ordered=T)
  
  #Means and std dev by group
  fms_centroids_sum[[k]] <- fms_centroids_long[[k]] %>% group_by(cluster, cellSubtype) %>% summarise(mean=mean(value), sd=sd(value))
  
  #Plot centroids by cluster  
  ggplot(fms_centroids_sum[[k]], aes(x=cellSubtype,y=mean, group=cluster)) + 
    geom_line() +
    geom_point(aes(color=cellSubtype, size=5)) +
    geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd, group=cellSubtype)) +
    xlab("Subset") +
    ylab("Expression") +
    labs(title= "Cluster Expression by Subset",color = "Cluster") +
    scale_color_manual(name="Cell Subtype", values=cellcolors) +
    facet_wrap(vars(cluster), ncol=1) +
    guides(size=FALSE)
  ggsave(filename=paste("plots/", "dotplot.fuzzy.", k, ".", subset, ".pdf", sep=""), width=3, height=6.8)
  
  #annotate cluster to output file
  data_save = eval(parse(text = paste0("cbind(data_save, k.", k, ".cluster = fms[[", k, "]]$cluster)")))
  
  #Annotate membership to output file
  for (i in 1:k){
    data_save = eval(parse(text = paste0("cbind(data_save, k.", k, ".clus.", i, ".mem = fms[[", k, "]]$membership[,",i,"])")))
  }
  
  #Heatmap
  #order data according k clustering
   dOrd = data[order(fms[[k]]$cluster, decreasing = T),]
   annotCol = gray((k - fms[[k]]$cluster[order(fms[[k]]$cluster, decreasing = T)]) / max(fms[[k]]$cluster))
   
   #plot heatmap
   scale = c(-1.5,1.5);
   ramp = colorRampPalette(rampColsDiv);   #Default colorRamp palette
   heatFile = paste0("plots/heatmap.fuzzy.k.", k, ".", subset, ".pdf")
   heatmap3(as.matrix(dOrd[1:colnum]), outFile = heatFile, rowClust = F, colClust = F, annotCol = cols, annotRow = annotCol, col = ramp(100), heatDim = c(2, 2), annotWid = 0.25, clusterWid = 0.4, key = F, keyWid = 0.3, flatten = T, hcMethod = "average", norm = norm, scale = scale, res = res, labDim = c(0.1, 0.1), height = height, width = width, outMai = c(0.1, 0.1, 0.1, 0.1), flattenRowClust = T, flattenColClust = T, flattenKey = T, labCol=colnames(dOrd[1:colnum]))
  
}

#Save data
write.table(data_save, file = paste0(subset, "_RNA_fuzzmeans_clusterdata_v2.txt"), sep="\t", row.names=F, quote=F)


