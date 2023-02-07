#pca plot with 99% confidence circles

library("som")
library("data.table")
#install.packages("vegan")
library("vegan");
library("preprocessCore")

#set working directory as coverage subfolder of the analysis folder
homeDir = "/Volumes/ESB/Boss_Projects/Human_Tcell_subsets/RNA_combined/analysis/pca/";
setwd(homeDir);

#read in the manifest file
filesDir = "/Volumes/ESB/Boss_Projects/Human_Tcell_subsets/RNA_combined/pipeline/";
filesFile = "RNAseq.sample.manifest.txt";
files = read.table(paste0(filesDir, filesFile), sep = "\t", header = TRUE, as.is = T);

#remove samples if necessary
files = files[which(files$cellType != "Monocyte"), ]
files = files[which(files$stim == "unstim"), ]

#sample color scheme
#grpOrder = c("CD8_bulk_blood_unstim", "CD8_bulk_blood_stim", "CD8_Nav_blood_unstim", "CD8_Nav_blood_stim", "CD8_EM_blood_unstim", "CD8_EM_blood_stim",
#            "CD8_EMRA_blood_unstim", "CD8_EMRA_blood_stim", "CD8_CM_blood_unstim", "CD8_CM_blood_stim", "CD8_bulk_tonsil",
#             "CD4_bulk_blood_unstim", "CD4_bulk_blood_stim", "CD4_Nav_blood_unstim", "CD4_Nav_blood_stim", "CD4_EM_blood_unstim", "CD4_EM_blood_stim", 
#             "CD4_EMRA_blood_unstim", "CD4_EMRA_blood_stim", "CD4_CM_blood_unstim", "CD4_CM_blood_stim", "CD4_TFH_blood_unstim", "CD4_TFH_blood_stim",
#             "CD4_nonTFH_tonsil", "CD4_TFH_tonsil", "Monocyte_blood_unstim", "Monocyte_blood_stim")

#set colors
cd4.bl.un = "orange"
cd4.bl.stim = "red"
cd8.bl.un = "blue"
cd8.bl.stim = "dark blue"
ton.col = "dark green"
mono.col = "purple"

grpcols = function(colnms) {if(grepl("CD8", colnms)) cd8.bl.un else if (grepl("CD4", colnms)) cd4.bl.un
  else if (grepl("Monocyte", colnms)) mono.col}

cols = unlist(lapply(files$cellType, grpcols))

grpCols = unlist(lapply(unique(files$group), grpcols))

#read in peak file
dataDir = "/Volumes/ESB/Boss_Projects/Human_Tcell_subsets/RNA_combined/coverage/";
dataFile = "Tcell_RNA.combined.RPKM.csv";
data = read.table(file=paste0(dataDir, dataFile), header=T, sep=",", comment.char="", quote="");

#select data to use
data.pca = data[ ,grep(paste(files$sample, collapse = "|"), names(data))]

#quantile normalize data
norm = normalize.quantiles(as.matrix(data.pca))
norm = data.frame(norm)
colnames(norm) = gsub(".rppm|.rpm|.fpkm|.rpkm", "", names(data.pca))

#remove extensions for headers
#colnames(data.pca) = gsub(".rppm|.rpm|.fpkm|.rpkm", "", names(data.pca))
#make the same order as the manifest
data.pca = norm[,as.vector(files$sample)]

#z-score normalize data by row
range(data.pca);
data.norm <- som::normalize(data.pca, byrow=TRUE);
range(data.norm);
#remove NaN values if necessary
#data.norm = na.omit(data.norm) 

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
xlim = c(-3, 3)
ylim = c(-3, 5)

#plot
cairo_pdf(paste0("All.nonStim.cellTypes.items.pc", p1, ".v.pc", p2, ".pdf"), height = 6, width = 6)
par(mai = c(0.8, 0.8, 0.5, 0.5), mgp = c(2, 0.8, 0), family = "Arial")

#pca plot
plot(pcas$species[, 1], pcas$species[, 2], pch = 19, xlab = xlab, ylab = ylab,
	cex.main = 1, cex.lab = 0.8, font.lab = 2, cex.axis = 0.6, col = cols,
	main = paste("PCA", dim(data.norm)[1]), cex = 0.3)
#add in 99% confidence intervals
#ordiellipse(pca, factor(files$group, level=unique(files$group)), conf = 0.99, kind = "se", display = "species", choices = c(p1, p2), col = grpCols)
#legend("topleft", legend = unique(files$group), text.col = grpCols, bty = "n")
text(pcas$species[, 1], pcas$species[, 2], labels = files$group, cex = 0.3, col = cols)
dev.off();


#run tsne on data
perplexity = 18 #this can be changed depending on the number of samples
tsne <- Rtsne(t(data.norm), dims = 2, perplexity = perplexity, verbose=TRUE, max_iter = 500)

cairo_pdf(file = "tSNE.zscore.pdf", height = 6, width = 6)
par(mai = c(0.8, 0.8, 0.5, 0.5), mgp = c(2, 0.8, 0), family = "Arial")
plot(tsne$Y, t='n', ylab = "tSNE2", xlab = "tSNE1")
points(tsne$Y, col = cols, cex = .3, pch = 19)
text(tsne$Y, labels = files$group, cex = 0.3, col = cols)
#ordiellipse(tsne.z$Y, factor(files$group, level=unique(files$group)), conf = 0.99, kind = "se", display = "species", choices = c(1, 2), col = grpCols)
#legend("bottomleft", legend = unique(files$group), text.col = grpCols, bty = "n")
dev.off();
