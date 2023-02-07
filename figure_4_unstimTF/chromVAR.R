#Script for calculating deviation and variability of ATAC data at TFBS
#jrose
#24Aug2020

library("here")
library("rtracklayer")
library("GenomicRanges")
library("SummarizedExperiment")
library("chromVAR")
library("BSgenome.Hsapiens.UCSC.hg38")
library("BiocParallel")
library("motifmatchr")
#library("tictoc")

#Directory
home_dir <- "analysis/chromVAR"

#Parallelization
register(MulticoreParam(8, progressbar = TRUE))

#Loading in peaks bed file and resizing to fixed width 250bp
org_peaks <- import.bed(here("analysis/unstim/diff/peaks_nohead.bed"))
peaks <- resize(org_peaks, width = 250, fix='center')

#Loading sample manifest
files = read.table(here('pipeline', 'ATACseq.sample.manifest.txt'), sep = "\t", header = T, as.is = T);
files = files[grepl("EM|CM|Nav", files$group),]
files = files[files$stim=="unstim",]
files = files[!grepl("CD4_EMRA", files$group),]

#Vector of bam file from ESB server & groups for column data
bamfiles=paste0(files$dir, files$bamFile)
groups=files$group

#Generating counts object
#tic()
data_obj <- getCounts(bamfiles, 
                         peaks, 
                         paired =  TRUE, 
                         by_rg = F, 
                         format = "bam", 
                         colData = DataFrame(celltype = groups)
                          )
#toc()

data_obj <- addGCBias(data_obj, 
                            genome = BSgenome.Hsapiens.UCSC.hg38)

save(data_obj, file=here(home_dir, "Tmem_unstim_sumexpfiles_gc.rda"))
#The steps above (creating count SumExp object) take the most time, if doing downstream analysis it's faster just to load SummExp obj instead of remaking it
#load(here(home_dir, "Tmem_unstim_sumexpfiles_gc.rda"))

#Since this is bulk I am skipping the filterSamples step which is used to exclude samples with low coverage
#This needs to be done with scATACseq though
# counts_filtered <- filterSamples(data_obj, min_depth = 1500, 
#                                  min_in_peaks = 0.15, shiny = T)
# 
# filtering_plot <- filterSamplesPlot(data_obj, min_depth = 1500, 
#                                     min_in_peaks = 0.15, use_plotly = T)
# filtering_plot

#I will filter out peaks with low counts though
counts_filtered <- filterPeaks(data_obj, non_overlapping = TRUE)

#Finding motifs from JASPER database
motifs <- getJasparMotifs()
#Matching to counts data
motif_ix <- matchMotifs(motifs, counts_filtered, 
                        genome = BSgenome.Hsapiens.UCSC.hg38)

#Finding the deviation and variability of TFBS
dev <- computeDeviations(object = counts_filtered, annotations = motif_ix)
variability <- computeVariability(dev)


#Saving deviation object
save(dev, variability, file=here(home_dir, "Tmem_unstim_dev.rda"))

#Plot variability graph
pdf(file=here(home_dir, "Tmem_unstim_variabilityplot.pdf"), useDingbats = F)
plotVariability(variability, use_plotly = F, n=20) 
dev.off()

#Tsne of results
tsne_results <- deviationsTsne(dev, threshold = 1.5, perplexity = 9)

tsne_plots <- plotDeviationsTsne(dev, tsne_results, 
                                 #annotation_name = "TBX21", 
                                 sample_column = "celltype", 
                                 shiny = FALSE)

pdf(file=here(home_dir, "Tmem_unstim_chromVARtsneplot.pdf"), useDingbats = F)
tsne_plots[[1]]
dev.off()