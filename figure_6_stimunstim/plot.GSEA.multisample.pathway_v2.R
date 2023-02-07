#Edited by jrose 20Mar20

#set working directory
homeDir = "/home/jim/Documents/Emory/Boss/data/RNA/GSEA/Stim_unstim/Hallmark/";
#setwd(homeDir);

#read in collapsed PreRanked Gene list used for GSEA (resides in the "edb" folder in the output directory)(
geneListDir = "/home/jim/Documents/Emory/Boss/data/RNA/GSEA/Stim_rnk_files/";
geneListFile = c("CD8_CM_blood_stim.v.CD8_Nav_blood_stim.rnk",
                 "CD8_EM_blood_stim.v.CD8_Nav_blood_stim.rnk",
                 "CD8_EMRA_blood_stim.v.CD8_Nav_blood_stim.rnk",
                 "CD4_CM_blood_stim.v.CD4_Nav_blood_stim.rnk",
                 "CD4_EM_blood_stim.v.CD4_CM_blood_stim.rnk") 
geneList <- list()
for (i in 1:length(geneListFile)){
geneList[[i]] = read.table(paste0(geneListDir, geneListFile[i]), sep = "\t")
geneList[[i]] = geneList[[i]][order(geneList[[i]][, 2], decreasing = T), ]
}

#Name of output file
output="TNFA_NFKB"
#plot title
main <- "TNFA Signaling via NFKB"

#Name of geneset
geneset <- "HALLMARK_TNFA_SIGNALING_VIA_NFKB"

#set colors for each group if desired
upCol = NA
dwnCol = NA

#set group names if desired
upPheno = NA
dwnPheno = NA

###
#point to .xls file for the gene set to plot
# dataDir = "/home/jim/Documents/Emory/Boss/data/RNA/GSEA/Stim_Tcurrated/";
# dataFile = c(paste("CD8_CM_Nav_stim.GseaPreranked.1605117807129/", geneset,".xls", sep=""),
#              paste("CD8_EM_Nav_stim.GseaPreranked.1605117791900/", geneset,".xls", sep=""),
#              paste("CD8_EMRA_Nav_stim.GseaPreranked.1605117778503/", geneset,".xls", sep=""),
#              paste("CD4_CM_Nav_stim.GseaPreranked.1605117911753/", geneset,".xls", sep=""),
#              paste("CD4_EM_Nav_stim.GseaPreranked.1605117816182/", geneset,".xls", sep="")
#              )
# dataDir = "/home/jim/Documents/Emory/Boss/data/RNA/GSEA/Stim_Analysis/";
# dataFile = c(paste("CD8_CM_Nav_stim.GseaPreranked.1605116177644/", geneset,".xls", sep=""),
#                           paste("CD8_EM_Nav_stim.GseaPreranked.1605116209182/", geneset,".xls", sep=""),
#                           paste("CD8_EMRA_Nav_stim.GseaPreranked.1605116218788/", geneset,".xls", sep=""),
#                           paste("CD4_CM_Nav_stim.GseaPreranked.1605116127445/", geneset,".xls", sep=""),
#                           paste("CD4_EM_Nav_stim.GseaPreranked.1605116139107/", geneset,".xls", sep="")
#                           )

dataDir = "/home/jim/Documents/Emory/Boss/data/RNA/GSEA/Stim_unstim/Hallmark/";
dataFile = c(paste("CD8_CM_stim_unstim.GseaPreranked.1605640668505/", geneset,".xls", sep=""),
                          paste("CD8_EM_stim_unstim.GseaPreranked.1605640676504/", geneset,".xls", sep=""),
                          paste("CD8_EMRA_stim_unstim.GseaPreranked.1605640682201/", geneset,".xls", sep=""),
                          paste("CD4_CM_stim_unstim.GseaPreranked.1605640636229/", geneset,".xls", sep=""),
                          paste("CD4_EM_stim_unstim.GseaPreranked.1605640658266/", geneset,".xls", sep="")
                          )


data <- list()
for (i in 1:length(dataFile)){
data[[i]] = read.table(paste0(dataDir, dataFile[i]), sep = "\t", header = T, comment.char="", quote="")
}

#set columns names to parse file
esCol = "RUNNING.ES"
rankCol = "RANK.IN.GENE.LIST"

#set axis range based on length of PreRanked gene list length
xlim = c(0, dim(geneList[[1]])[1])

#Edit ylim depending on whether they are positively (0, 1) or negatively (-1, 0) enriched
ylim = c(-1, 1)

# #set color ramp for ranks ranked gene list plot
# if(is.na(upCol)) 
# 	{upCol = rgb(31, 73, 125, maxColorValue = 255)}
# if(is.na(dwnCol)) {dwnCol = rgb(192, 80, 77, maxColorValue = 255)}
# 
# rampCols = c(upCol, rgb(0.5, 0.5,0.5), rgb(1,1,1), rgb(0.5, 0.5, 0.5), dwnCol)    # 5 colors scale
# #rampCols = c(rgb(0.5, 0.5,0.5))      # all gray color scale
# ramp = colorRampPalette(rampCols)

#List of colors for samples
smp_col <- c("#0C46A0FF","#E55100FF", "#870D4EFF", "#2096F2FF","#FFB74CFF")

#ploting parameters
height = 5; 
width  = 5;

#set the sample names for the plot by simplifying the file name
names <- vector()
for (i in 1:length(dataFile)){
names[i] = tolower(gsub("_", " ", gsub("\\..*", "", dataFile[[i]])))
}

#set the output file name or output to screen
#X11(height = height, width = width)
pdf(paste0(homeDir, output, ".pdf"), height = height, width = width); 

#set plot layout for 3 sections
layout(matrix(1:(length(names)+1), nrow = (length(names)+1)), widths = 1, heights = c(2, rep(0.2, length(names))))

#section 1: enrichment plot
#section 1 layout options
par(mai = c(0, 0.5, 0.25, 0.25), mgp = c(1.2, 0.4, 0), bty = "n", xaxs = "i", yaxs = "i")
#plot axis
plot(NA, xlim = xlim, ylim = ylim, xaxt = "n", xlab = NA, ylab = "Enrichment Score", main=main, cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.1, font.lab = 2);
lines(c(0, dim(geneList[[1]])[1]), c(0,0), lwd = 1, col = rgb(0.5, 0.5, 0.5))
#plot the green enrichment line
for (i in 1:length(geneList)){
lines(c(0, data[[i]][, rankCol], dim(geneList[[i]])[1]), c(0, data[[i]][, esCol], 0), lwd = 3, col = smp_col[i])
}
#legend
legend(x="bottomleft",legend=names, col=smp_col, lty=1, cex=1.2)

#section 2: tick section showing where all the genes are in your gene set
#section 2 layout options
for (i in 1:length(data)){
  par(mai = c(0, 0.5, 0, 0.25), mgp = c(0.2, 0, 0), bty = "o", xaxs = "i", yaxs = "i", las = 0.9)
  #plot box
  plot(NA, xlim = xlim, ylim = c(0, 1), ylab = NA, xlab = NA, yaxt = "n", cex.axis = 0.8, font.lab = 2, xaxt = "n");
  #loop through and plot a tick for each gene
  for (j in 1:dim(data[[i]])[1]) lines(rep(data[[i]][j, rankCol], 2), c(0, 1), lwd = 0.2, col = smp_col[i])
}
#section 3: barplot section showing the distribution of rank values
#section 3 layout options
#par(mai = c(0.25, 0.5, 0, 0.25), mgp = c(0.2, 1, 0), bty = "o", xaxs = "i", yaxt = "n", las = 1)
#loop through and plot a tick for each gene
#barplot(geneList[,2], ylab = NA, xlab = "Preranked Genes", cex.lab = 0.8, cex.axis = 0.8, font.lab = 2, border=ramp(length(geneList[,2])))
#med = median(which(geneList[,2]==0))
#lines(c(med,med), range(geneList[,2]), lwd = 0.5, col = rgb(0, 0, 0), lty=2)
#add phenotype labels
#if(!is.na(dwnPheno)) {text(dim(geneList)[1]-(dim(geneList)[1] * 0.001), geneList[dim(geneList)[1],2] * 0.5, dwnPheno, col = dwnCol)}
#if(!is.na(upPheno)) {text(1+(dim(geneList)[1] * 0.2), geneList[1,2] * 0.5, upPheno, col = upCol)}

dev.off()

