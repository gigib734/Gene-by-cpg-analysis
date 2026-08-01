methylation_data <- read.csv("C:/Users/gbloc/OneDrive/Desktop/PCDHG_Gila/Coding/Methylation & Transcription Correlation/Meta-analysis results DMPs blood PCDH genes - Meta-analysis results DMPs blood PCDH genes (1).csv")
transcription_data <- read.csv("C:/Users/gbloc/OneDrive/Desktop/PCDHG_Gila/Coding/Methylation & Transcription Correlation/Meta-analysis results DEGs blood PCDH genes - Meta-analysis results DEGs blood PCDH genes.csv.csv")
library(dplyr)
library(tidyverse)
library(readr)
library(ggplot2)

#effect size plot
library(tidyverse)

#first need to correct for multiple testing
transcription_data$FDR=p.adjust(transcription_data$P.value, method="fdr", n=length(transcription_data$P.value))
mRNA <- transcription_data%>%
  #filter(FDR<0.05)%>%
  dplyr::select("MarkerName" ,
                `Effect`)%>%
  dplyr::rename(Gene="MarkerName")%>%
  dplyr::filter(Gene %in% c("PCDHGA1", "PCDHGA2", "PCDHGA3", "PCDHGA4", "PCDHGB1"))


DNAm <- methylation_data %>%
  filter(FDR<0.05)%>%
  dplyr::rename(Gene = `Annotated_genes`)%>%
  separate_rows(Gene,
                sep=";")%>%
  dplyr::select(Gene,
                `Effect_size`)%>%
  dplyr::filter(Gene %in% c("PCDHGA1", "PCDHGA2", "PCDHGA3", "PCDHGA4", "PCDHGB1"))


merged <- left_join(DNAm,
                    mRNA)%>%
  drop_na()



setwd("//ad.monash.edu/home/User050/mjac0029/Desktop/Papers in progress/Methylation Multi OMICS skeletal muscle VO2 and exercise/Datasets VOmax/")
library(ggrepel)
tiff('Effect size DNAm vs mRNA for common genes.tiff',
     width =100,
     height = 80,
     units = 'mm',
     res=300)
ggplot(data = merged,
       aes(x = `Effect_size`,
           y = `Effect`,
           label = Gene)) +
  xlab("% DNAm change per chronological age")+
  ylab("mRNA log2FC per unit of chronological age")+
  geom_point(colour="#1AEDD8")+
  geom_label_repel(data = merged,
                   size = 2,
                   box.padding = unit(0.45, "lines"),
                   point.padding = unit(0.45, "lines"),
                   max.overlaps = 30)+
  geom_hline(yintercept = 0,
             lty="dashed")+
  geom_vline(xintercept = 0,
             lty="dashed")+
  theme_classic()
dev.off()

#based on chrom state

DNAm <- methylation_data %>%
  filter(FDR<0.05)%>%
  dplyr::rename(Gene = `Annotated_genes`)%>%
  separate_rows(Gene,
                sep=";")%>%
  dplyr::select(Gene,
                `Effect_size`,
                chr_state,
                CGI_position)%>%
  dplyr::filter(Gene %in% c("PCDHGA1", "PCDHGA2", "PCDHGA3", "PCDHGA4", "PCDHGB1"))%>%
  dplyr::filter(chr_state %in% c("Active TSS","Flanking bivalent TSS/Enh","Flanking active TSS","Bivalent/poised TSS"))

merged <- left_join(DNAm,
                    mRNA)%>%
  drop_na()



setwd("//ad.monash.edu/home/User050/mjac0029/Desktop/Papers in progress/Methylation Multi OMICS skeletal muscle VO2 and exercise/Datasets VOmax/")
library(ggrepel)
tiff('Effect size DNAm vs mRNA for common genes.tiff',
     width =100,
     height = 80,
     units = 'mm',
     res=300)
ggplot(data = merged,
       aes(x = `Effect_size`,
           y = `Effect`,
           label = Gene)) +
  xlab("% DNAm change per chronological age")+
  ylab("mRNA log2FC per unit of chronological age")+
  geom_point(colour="#1AEDD8")+
  geom_label_repel(data = merged,
                   size = 2,
                   box.padding = unit(0.45, "lines"),
                   point.padding = unit(0.45, "lines"),
                   max.overlaps = 30)+
  geom_hline(yintercept = 0,
             lty="dashed")+
  geom_vline(xintercept = 0,
             lty="dashed")+
  theme_classic()
dev.off()