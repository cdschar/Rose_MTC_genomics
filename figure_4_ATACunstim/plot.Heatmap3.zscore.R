#HeatMap plotting script
library("som")

setwd("/Volumes/ESB/Boss_Projects/Human_Tcell_subsets/ATAC_combined/analysis/unstim/heatmap/")

#read in heatmap3 function
source("/Volumes/GRAID/seqTools/R.Libraries/heatmap3/heatmap3.R")

#set file name
experimentName = "T_cell_ATAC"

#read in the manifest file
fqDir = "/Volumes/ESB/Boss_Projects/Human_Tcell_subsets/ATAC_combined/pipeline/"
fqFile = "ATACseq.sample.manifest.txt";
files = read.table(paste0(fqDir, fqFile), sep = "\t", header = T, as.is = T);
files = files[grepl("EM|CM|Nav", files$group),]
files = files[files$stim=="unstim",]
files = files[!grepl("CD4_EMRA", files$group),]

#sample color scheme
g1.col = rgb(200,10,10, maxColorValue = 255);
g2.col = rgb(10,10,200, maxColorValue = 255);

#set up column colors to match number of samples in each group
#Colcols = c(rep(g1.col, 3), rep(g2.col, 3));
Colcols = gsub("CD4_CM_blood_unstim", "#2096F2FF", files$group)
#Colcols = gsub("CD4_EMRA_blood_unstim", "red", Colcols)
Colcols = gsub("CD4_EM_blood_unstim", "#FFB74CFF", Colcols)
Colcols = gsub("CD4_Nav_blood_unstim", "#90A4ADFF", Colcols)
Colcols = gsub("CD8_CM_blood_unstim", "#0C46A0FF", Colcols)
Colcols = gsub("CD8_EMRA_blood_unstim", "#870D4EFF", Colcols)
Colcols = gsub("CD8_EM_blood_unstim", "#E55100FF", Colcols)
Colcols = gsub("CD8_Nav_blood_unstim", "#5F7D8BFF", Colcols)

#set up row colors if desired
Rowcols = NA

#read in data file
dataDir = "/Volumes/ESB/Boss_Projects/Human_Tcell_subsets/ATAC_combined/analysis/unstim/diff/";
dataFile = "diff.significant.glm.ATACseq.txt";
data = read.table(file=paste0(dataDir, dataFile), header=T, sep="\t", comment.char="", quote="");

#use sig data only
#data.sig = data[which(data$sigAny == "TRUE"), ];

#select data columns
data.heat = data[ ,grep(paste(files$sample, collapse = "|"), names(data))]

#order data to be same as manifest
colnames(data.heat) = gsub(".rppm|.rpm|.fpkm|.rpkm", "", names(data.heat))
data.heat = data.heat[,as.vector(files$sample)]

#z-score normalize data by row
range(data.heat);
data.norm <- som::normalize(data.heat, byrow=TRUE);
range(data.norm);
#remove NaN values if necessary
#data.norm = na.omit(data.norm)

#heatmap size, resolution, flatten options
height = 2.5 #pdf size in inches
width = 2.5 #pdf size in inches
flatten = T #overall control for heatmap, key, and dendrograms
res = 600 #resolution of heatmap

#controls for internal margins
labDim = c(0.7, 0.4) #change the size of the label area for rows or columns
outMai = c(0.2, 0.1, 0.1, 0.1)
heatDim = c(3, 3)
annotSp = c(0.025)

#heatmap key options
key = TRUE #controls if key is plotted
keyWid = 0.5 #controls size of key
flattenKey = flatten #determines if key is flattened

#controls for clustering and dendrogram display/size
clustMeth = "average" #options: ward.D, ward.D2, single, complete, average, mcquitty, median ,or centroid
rowClust = TRUE #determines if rows are clustered
colClust = TRUE #determines if columns are clustered
annotWid = 0.25
clusterWid = 0.4
flattenColClust = flatten #flatten column dendrogram options
flattenRowClust = flatten #flatten row dendrogram options
colDen.lwd = 1 #line width for column dendrogram
rowDen.lwd = 0.1 #line width for row dendrogram

#heatmap colors (blue-white-red)
rampCols = c(rgb(31, 73, 125, maxColorValue = 255), rgb(1,1,1), rgb(192, 80, 77, maxColorValue = 255));
ramp = colorRampPalette(rampCols);	#Default colorRamp palette
breaks = 100 #number of bins in color scale
scale = c(-2, 2) #cap Z-score range of data

#label options
colLab = gsub("[[:punct:]]", " ", names(data.heat))
rowLab = NA
labCex = 0.6 #size of row and column lable font

#heatmap output file
heatFile = paste0(experimentName, ".heatmap3.pdf")

#Generate Heatmap
heatmap3(as.matrix(data.norm), outFile = heatFile, rowClust = T, colClust = T,
	col = ramp(breaks), annotCol = Colcols, annotRow = Rowcols, heatDim = heatDim, annotWid = annotWid, clusterWid = clusterWid,
	key = key, keyWid = keyWid, flatten = flatten, res = res, labDim = labDim, height = height, width = width,
	outMai = outMai, flattenColClust = flattenColClust, flattenRowClust = flattenRowClust, clustMeth = clustMeth,
	scale = scale, colLab = colLab, rowLab = rowLab, labCex = labCex, colDen.lwd = colDen.lwd, rowDen.lwd = rowDen.lwd)
