#Creating simple bar chart of up and down regulated DEGs from Naive vs Memory cells
#jrose
#14Feb20

library(tidyverse)
library(here)

# Loading data

input_dir <- "analysis/unstim/diff_unstim"
output_dir <- "analysis/unstim/diff_unstim/naive_v_mem_overlap"

counts <- read.delim(here(input_dir, "deg_counts.txt"), header=T, sep="\t")
counts_sub <- counts[grep("^CD(.)_Nav_blood_unstim.v.CD\\1_", counts$Comp),]
#^Using a regex expression to find only Nav vs Mem comparison of the same lineage (CD4/CD8)

counts_sub <- counts_sub[!grepl("CD4_Nav_blood_unstim.v.CD4_EMRA_blood_unstim.sig", counts_sub$Comp),]

#Adding additional data to df
counts_sub <- counts_sub %>% mutate(neg_neg = -neg)
subtype <- strsplit(as.character(counts_sub$Comp), "_")
for (i in 1:length(subtype)){
  counts_sub$subtype[i] <- paste(subtype[[i]][1],subtype[[i]][5], sep="_")
}

#Plotting...yes I'm "plotting" something...

theme_set(theme_gray()+ theme(axis.line = element_line(size=0.5),panel.background =
                                element_rect(fill=NA,size=rel(20)),
				panel.grid.minor = element_line(colour = NA),
                              	axis.text = element_text(size=12, angle=45, hjust=1),
				axis.title = element_text(size=12),
				legend.text = element_text(size=12)
				)
	 )

#Colors (found from awtools:mpalette)
upcol <- "#017a4a"
downcol <- "#ff363c"

p <- ggplot(counts_sub, aes(x=subtype))
p + geom_bar(aes(y=pos, fill="up-regulated vs Naive"), stat='identity') + 
  geom_bar(aes(y=neg_neg, fill="down-regulated vs Naive"), stat='identity') + 
  scale_fill_manual("",values = c("up-regulated vs Naive" = upcol, "down-regulated vs Naive" = downcol),
                    limits=c("up-regulated vs Naive", "down-regulated vs Naive")) +
  ylab("Number of Differentially Expressed Genes (FDR<0.05)") + xlab('Memory Subtype') +
  scale_y_continuous(limits=c(-1100,1100))

ggsave(here(output_dir, "DEGbarchart_noCD4emra.pdf"))
