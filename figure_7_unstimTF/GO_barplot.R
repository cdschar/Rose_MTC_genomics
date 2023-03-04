library(here)
library(tidyverse)
library(RColorBrewer)

input_dir <- 'PageRank/Run1/Network/heatmaps'
output_dir <- 'PageRank/Run1/Network/heatmaps'
file <- 'MSC_GOanalysis.txt'

subset <- "MSC"

GOdf <- read.delim(file=here(input_dir,file), header = T)
names(GOdf) <- gsub("Client.Text.Box.Input..", "", names(GOdf))
head(GOdf)

GOdf$fold.Enrichment. <- as.numeric(GOdf$fold.Enrichment.)

data <- GOdf %>% mutate(neglog10FDR=-log(FDR., base=10), log2Fenrich=log(fold.Enrichment., base=2))
data <- arrange(data, desc(log2Fenrich))

p_tmp <- ggplot(data, aes(x=log2Fenrich, y=neglog10FDR)) +
     geom_point()
p_tmp

data_sub <- subset(data, log2Fenrich>4 & neglog10FDR>2)
#data_sub <- data
data_sub <- top_n(data_sub, 10, log2Fenrich) %>% arrange(log2Fenrich)
data_sub$GO.biological.process.complete <- factor(data_sub$GO.biological.process.complete, levels=data_sub$GO.biological.process.complete, ordered = T)


p <- ggplot(data=data_sub, aes(x=GO.biological.process.complete, y=neglog10FDR)) +
  geom_col(aes(fill=log2Fenrich)) +
  xlab("Pathway") +
  ylab("-log10(Pvalue)") +
  scale_fill_gradient(low="white", high="firebrick", limits=c(0,10)) +
  scale_y_continuous(limits=c(0,8)) +
  coord_flip() + 
  theme(axis.text = element_text(size=5))
p
#ggsave(here(output_dir, paste0(subset,"_GO_barplot_top20_clus4.pdf")), height=7, width=7)


data_2 <- top_n(data, 10, neglog10FDR) %>% arrange(neglog10FDR)
data_2$GO.biological.process.complete <- factor(data_2$GO.biological.process.complete, levels=data_2$GO.biological.process.complete, ordered = T)
p2 <- ggplot(data=data_2, aes(x=GO.biological.process.complete, y=neglog10FDR)) +
  geom_col(aes(fill=log2Fenrich)) +
  xlab("Pathway") +
  ylab("-log10(Pvalue)") +
  scale_fill_gradient(low="white", high="firebrick", limits=c(0,10)) +
  #scale_y_continuous(limits=c(0,10)) +
  coord_flip() +
  guides(fill=guide_legend(title="Fold Enrichment (log2)")) + 
  theme(axis.text.y=element_text(size=15),axis.text.x=element_text(size=10), legend.text = element_text(size=10))
p2
#ggsave(here(output_dir, paste0(subset,"_GO_barplot_topFDR_clus5.pdf")), height=7, width=10)

theme_set(theme_gray()+ theme(axis.line = element_line(size=0.5),
                              panel.background = element_rect(fill=NA,size=rel(20)),
                              panel.grid.minor = element_line(colour = NA),
                              axis.text = element_text(size=10), axis.title = element_text(size=12)))



data_3 <- subset(data, `94.`<=25) %>% top_n(15, neglog10FDR) %>% arrange(neglog10FDR)
data_3$GO.biological.process.complete <- factor(data_3$GO.biological.process.complete, levels=data_3$GO.biological.process.complete, ordered = T)
p3 <- ggplot(data=data_3, aes(x=GO.biological.process.complete, y=neglog10FDR)) +
  geom_col(aes(fill=log2Fenrich)) +
  xlab("Pathway") +
  ylab("-log10(Pvalue)") +
  scale_fill_gradient(low="white", high="firebrick", limits=c(0,10)) +
  #scale_y_continuous(limits=c(0,10)) +
  coord_flip() +
  guides(fill=guide_legend(title="Fold Enrichment (log2)")) + 
  theme(axis.text.y=element_text(size=15),axis.text.x=element_text(size=10), legend.text = element_text(size=10))
p3
ggsave(here(output_dir, paste0(subset,"_GO_barplot_topFDR.pdf")), height=5, width=10)

ELSE <- TRUE
data_4 <- subset(data, `94.`<=25 & neglog10FDR >=1.3)  %>% top_n(100, neglog10FDR) %>% 
  mutate(category = case_when(
                          grepl("activation", GO.biological.process.complete, ignore.case = T) ~ "Cell Activation",
                          grepl("adhesion", GO.biological.process.complete, ignore.case = T) ~ "Cell Adhesion", 
                          grepl("differentiation|selection|hemopoiesis", GO.biological.process.complete, ignore.case = T) ~ "T Cell Differentiation",
                          grepl("locomotion|migration|motility|movement", GO.biological.process.complete, ignore.case = T) ~ "Migration",
                          grepl("cytokine", GO.biological.process.complete, ignore.case = T) ~ "Cytokine Response/Production",
                          grepl("death|apopto", GO.biological.process.complete, ignore.case = T) ~ "Cell Migration",
                          grepl("immune system process|immune response", GO.biological.process.complete, ignore.case = T) ~ "Immune Process",
                          ELSE ~ "Other"
                        )
  
        ) %>% 
        group_by(category) %>%
        summarise(number=n(),
                  ave_foldenrich = mean(fold.Enrichment.),
                  ave_neglog10FDR=mean(neglog10FDR),
                  sd_neglog10FDR=sd(neglog10FDR)
                  ) %>%
        arrange(ave_neglog10FDR)
data_4$category <- factor(data_4$category, levels=data_4$category, ordered = T)
#data_4$number <- factor(data_4$number)

data_5 <- subset(data, `94.`<=25 & neglog10FDR >=1.3)  %>% top_n(100, neglog10FDR) %>% 
  mutate(category = case_when(
    grepl("activation", GO.biological.process.complete, ignore.case = T) ~ "Cell Activation",
    grepl("adhesion", GO.biological.process.complete, ignore.case = T) ~ "Cell Adhesion", 
    grepl("differentiation|selection|hemopoiesis", GO.biological.process.complete, ignore.case = T) ~ "T Cell Differentiation",
    grepl("locomotion|migration|motility|movement", GO.biological.process.complete, ignore.case = T) ~ "Migration",
    grepl("cytokine", GO.biological.process.complete, ignore.case = T) ~ "Cytokine Response/Production",
    grepl("death|apopto", GO.biological.process.complete, ignore.case = T) ~ "Cell Migration",
    grepl("immune system process|immune response", GO.biological.process.complete, ignore.case = T) ~ "Immune Process",
    ELSE ~ "Other"
          )
    )

data_5$category <- factor(data_5$category, levels=levels(data_4$category), ordered = T)
ggplot(data_5, aes(x=category,y=neglog10FDR)) + geom_jitter() + coord_flip()

colors <- brewer.pal(9,'YlOrRd')

ggplot(data_4, aes(x=category,y=ave_neglog10FDR)) + geom_col(aes(fill=number)) +
  coord_flip() + 
  scale_fill_stepsn(colors=colors, labels=c("<10", "10-20", ">30")) + 
  geom_jitter(data=data_5, aes(x=category, y=neglog10FDR)) +
  geom_errorbar(aes(ymin=ave_neglog10FDR-sd_neglog10FDR, ymax=ave_neglog10FDR+sd_neglog10FDR)) +
  guides(fill=guide_legend(title="# Pathways Enriched")) + 
  theme(axis.text.y=element_text(size=15),
        axis.text.x=element_text(size=12), 
        legend.text = element_text(size=10),
        axis.title.y=element_blank()
        )
ggsave(here(output_dir, paste0(subset,"_GO_barplot_overlaps.pdf")), height=5, width=7)
