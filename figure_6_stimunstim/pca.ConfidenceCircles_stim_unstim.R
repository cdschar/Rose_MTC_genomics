#pca plot with 99% confidence circles
#Stim samples from Memory subsets only
#CD4 Emra group removed due to bad sorting gates

library("som")
library("data.table")
library("ggsci");
library("tidyverse");
library("vegan");

#set working directory as coverage subfolder of the analysis folder
homeDir = "/Volumes/ESB/Boss_Projects/Human_Tcell_subsets/RNA_combined/analysis/pca/";
setwd(homeDir);

#read in the manifest file
filesDir = "/Volumes/ESB/Boss_Projects/Human_Tcell_subsets/RNA_combined/pipeline/";
filesFile = "RNAseq.sample.manifest.txt";
files = read.table(paste0(filesDir, filesFile), sep = "\t", header = TRUE, as.is = T);
files = files[grepl("EM|CM|Nav", files$group),]
files = files[!grepl("CD4_EMRA", files$group),]
#files = files[files$stim=="stim",]

#Shortening the group labels
files$group <- gsub("_blood_unstim", "", files$group)
files$group <- gsub("_blood_stim", "_Stim", files$group)

#setting color scheme
pal <- c(rep("#2096F2FF",2), rep("#FFB74CFF",2), rep("#90A4ADFF",2), rep("#0C46A0FF",2),rep("#E55100FF",2),rep("#870D4EFF",2), rep("#455964FF",2));

#colors for each group, should be in same order as manifest
#unique(files$group) #check order of groups
grpCols = pal[unique(as.numeric(factor(files$group)))]

#colors for each sample, should be in same order as the manifest
cols = pal[as.numeric(factor(files$group))]

#read in peak file
dataDir = "/Volumes/ESB/Boss_Projects/Human_Tcell_subsets/RNA_combined/analysis/diff/";
dataFile = "diff.significant.glm.Human_Tcell_RNA.txt";
data = read.table(file=paste0(dataDir, dataFile), header=T, sep="\t", comment.char="", quote="");

#select data to use
data.pca = data[ ,grep(paste(files$sample, collapse = "|"), names(data))]

#use sig only peaks/genes if needed
data.sig = data[which(data$sigAny == "TRUE"), ];
data.pca = data.sig[ ,grep(paste(files$sample, collapse = "|"), names(data.sig))]

#remove extensions for headers
colnames(data.pca) = gsub(".rppm|.rpm|.fpkm|.rpkm", "", names(data.pca))
#make the same order as the manifest
data.pca = data.pca[,as.vector(files$sample)]

#z-score normalize data by row
range(data.pca);
data.norm <- som::normalize(data.pca, byrow=TRUE);
range(data.norm);
#remove NaN values if necessary
data.norm = na.omit(data.norm)

#perform PCA analysis
pca = rda(data.norm, scale = F)

#extract variance
pca.var = eigenvals(pca)/sum(eigenvals(pca))

#select components to plot
p1 = 1; #x-axis
p2 = 2; #y-axis

#axis labels with variance
xlab = paste0("PC", p1, " ", signif(pca.var[p1] * 100, digits = 3), "%")
ylab = paste0("PC", p2, " ", signif(pca.var[p2] * 100, digits = 3), "%")

#access ordination values for select pc
pcas = scores(pca, choices = c(p1, p2))

#set plot range
xlim = c(-5, 4)
ylim = c(-5, 3.5)

#plot
cairo_pdf(paste0(gsub(".txt", "", dataFile), ".PCAconfidence.", dim(data)[1], ".items.pc", p1, ".v.pc", p2, "stim_unstim", ".pdf"), height = 5, width = 5)
par(mai = c(0.8, 0.8, 0.5, 0.5), mgp = c(2, 0.8, 0), family = "Arial")

#pca plot
plot(pcas$species[, 1], pcas$species[, 2], pch = c(17,19)[as.numeric(as.factor(files$stim))],
	xlab = xlab, ylab = ylab,
	cex.main = 1, cex.lab = 0.8, font.lab = 2, cex.axis = 0.6, col = cols,
	main = paste("PCA of", dim(data.norm)[1], "DEGs"), cex = 1, ylim = ylim, xlim = xlim)
#add in 99% confidence intervals
ordiellipse(pca, factor(files$group, level=unique(files$group)), conf = 0.99, kind = "se", display = "species", choices = c(p1, p2), col = grpCols)
legend("bottomright", legend = unique(files$group), text.col = grpCols, col = grpCols, pch =rep(c(17,19), 7), bty = "n", cex=0.6)
#text(pcas$species[, 1], pcas$species[, 2], labels = files$group, cex = 0.3, col = cols)
dev.off();
