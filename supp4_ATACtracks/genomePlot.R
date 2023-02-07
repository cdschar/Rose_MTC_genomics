# Import plotting function
source("/home/boss_lab/Apps/bitbucket/genomePlot/Plot.region.lib.R") # Magic

projectDir = "/home/boss_lab/Projects/Boss_Projects/Human_Tcell_subsets/ATAC_combined/" # Replace, must have '/' at the end

home.dir = paste0(projectDir, "analysis/genomePlots/jrose_genomePlots/")
setwd(home.dir)
outputDir = paste0(home.dir, "plots/") # Replace, if necessary

# Output format, options: "png", "pdf", "svg", "screen"
output = "pdf" # Replace

# Set genome, options: "Hsapiens.UCSC.hg38", "Mmusculus.UCSC.mm10", "Hsapiens.UCSC.hg19", "Mmusculus.UCSC.mm9"
org = "Hsapiens.UCSC.hg38" # Replace

# Define tracks to plot
track.dir = home.dir # Replace, if necessary
track.file = "tracks_AHRHIF.txt" # Replace, will be the output of create_tracks.R
trackFile = paste0(track.dir, track.file)

# Define regions to plot
regions.dir = home.dir # Replace, if necessary
regions.file = "regions.txt" # Replace, if necessary
regions = read.table(paste0(regions.dir, regions.file), header = TRUE, sep = "\t", as.is = TRUE)

# Limit regions to plot
regions.plot = regions[regions$plot, ]

#annotate a genomic sequence
plot.pattern = FALSE
pattern = "CTGCAG"  #this could be CG, a restriction site, or any other sequence
pattern.height = 1 #controls the height of the pattern ticks
pattern.col = c(rgb(0, 0.5, 0)) #sets color for pattern, default dark green
pattern.lwd = 1

#plot options
height = 3 #set height in inches
fixed.width = c(TRUE) #TRUE uses the width option below. FALSE scales to bp.length option below
width = 5
bp.length = c(0.001) #bp length for FALSE fixed.width setting
outer.margin = c(0.4, 0.8, 0.25, 0.25)
track.margin = c(0.02, 0, 0.02, 0)
mgp = c(3, 1, 0)

#gene TSS direction arrows
arrow.length = 0.05
arrow.lwd = 1
arrow.height = 0.8
arrowhead.size = 0.1

#axis options
axis.lwd = c(0.5)
axis.cex = 0.5 #set master font size for labels

#x-axis options
xlab.size = 1
xaxis.size = 1
xlab.line = 1
xaxis.win.height = 0.05
xlab.inc = c(5000)
xlab.inc.on = c(FALSE)

#y-axis options
ylab.inc = c(0.5)
ylab.inc.on = c(FALSE)

#legend
legend.horiz = TRUE
legend.loc = c("topleft")
legend.size = 1

#draw 10kb scale bar (TRUE) or axis bar across bottom (FALSE)
scale = TRUE

for (region in 1:dim(regions.plot)[1]) {
  
  #set output file
  outFile = paste(outputDir, regions.plot$name[region], ".", output, sep = "");
  
  #call plot function
  plotRegion(trackFile, chr = as.character(regions.plot$chr[region]), start = as.numeric(regions.plot$start[region]), end = as.numeric(regions.plot$end[region]), 
             output = output, output.file = outFile, org = org, plot.fwd = regions.plot$plot.fwd[region], width = width, plot.pattern = plot.pattern, pattern = pattern,  arrow.length = arrow.length, 
             xlab.size = xlab.size, xaxis.size = xaxis.size, xlab.line = xlab.line, xaxis.win.height = xaxis.win.height, outer.margin = outer.margin, track.margin = track.margin, legend.horiz = legend.horiz,
             scale.inc.on = scale, pattern.lwd = pattern.lwd);
  
}

#capture R session info
sessionDir = "sessionInfo/"
if (!file.exists(sessionDir)) dir.create(sessionDir);
#plots = paste0(regions.plot$name, collapse = ".")
capture.output(sessionInfo(), file = paste0(sessionDir, "Rsession.Info.", gsub("\\D", "", Sys.time()), ".txt"))
