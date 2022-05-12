#import libraries
library(tidyverse)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
#import norm. counts
dna_cnts <- read_csv('dnadmg_deseq_norm_counts.csv')%>% column_to_rownames("...1")
ahr_cnts <- read_csv('ahr_deseq_norm_counts.csv') %>% column_to_rownames("...1")
car_cnts <- read_csv('carpxr_deseq_norm_counts.csv') %>% column_to_rownames("...1")

#plot heatmaps-
pheatmap(as.matrix(dna_cnts), show_rownames = FALSE, scale = 'row', color=colorRampPalette(c("red","white","blue"))(100))
pheatmap(as.matrix(ahr_cnts), show_rownames = FALSE, scale = 'row', color=colorRampPalette(c("red","white","blue"))(100))
pheatmap(as.matrix(car_cnts), show_rownames = FALSE, scale = 'row', color=colorRampPalette(c("red","white","blue"))(100))