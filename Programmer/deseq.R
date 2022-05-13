#import libraries
library(tidyverse)
library(ggplot2)
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
library(DESeq2)
library(apeglm)
#import sample counts and metadata
ctrl <- read_csv('control_counts.csv') %>% select(Geneid, SRR1178004, SRR1178006, SRR1178013, SRR1178050, SRR1178061, SRR1178063)
trtmt <- read_csv('counts.csv')
meta <- read_csv('toxgroup_3_rna_info.csv')
#create combined counts tibble and filter out rows genes with 0 counts
tot_counts <- trtmt %>% inner_join(ctrl, by = c("gene_id"="Geneid")) %>% select(-1)
tot_counts <- subset(tot_counts, rowSums(tot_counts==0)==0) %>% column_to_rownames("gene_id")
#make subset counts for the different mode of actions-
dnadmg_cnts <- tot_counts %>% select(SRR1177981, SRR1177982, SRR1177983, SRR1178004, SRR1178006, SRR1178013)
dnadmg_meta <- meta %>% filter(vehicle=="SALINE_100_%")
ahr_cnts <- tot_counts %>% select(SRR1178008, SRR1178009, SRR1178010, SRR1178050, SRR1178061, SRR1178063)
ahr_meta <- meta %>% filter(vehicle=="CORN_OIL_100_%") %>% filter(mode_of_action != "CAR/PXR")
carpxr_cnts <- tot_counts %>% select(SRR1178014, SRR1178021, SRR1178047, SRR1178050, SRR1178061, SRR1178063)
carpxr_meta <- meta %>% filter(vehicle=="CORN_OIL_100_%") %>% filter(mode_of_action != "AhR")
#make deseq objects
dnadmg_dds <- DESeqDataSetFromMatrix(countData = dnadmg_cnts, colData = dnadmg_meta, design = ~mode_of_action)
dnadmg_dds$mode_of_action <- relevel(dnadmg_dds$mode_of_action, ref="Control")
ahr_dds <- DESeqDataSetFromMatrix(countData = ahr_cnts, colData = ahr_meta, design = ~mode_of_action)
ahr_dds$mode_of_action <- relevel(ahr_dds$mode_of_action, ref="Control")
carpxr_dds <- DESeqDataSetFromMatrix(countData = carpxr_cnts, colData = carpxr_meta, design = ~mode_of_action)
carpxr_dds$mode_of_action <- relevel(carpxr_dds$mode_of_action, ref="Control")
#peforming deseq-
dnadmg_dds <- DESeq(dnadmg_dds)
ahr_dds <- DESeq(ahr_dds)
carpxr_dds <- DESeq(carpxr_dds)
#extract results
dnadmg_res <- results(dnadmg_dds, contrast = c("mode_of_action", "DNA_Damage", "Control"))
dnadmg_res <- lfcShrink(dnadmg_dds, coef = 2)
ahr_res <- results(ahr_dds, contrast = c("mode_of_action", "AhR", "Control"))
ahr_res <- lfcShrink(ahr_dds, coef = 2)
carpxr_res <- results(carpxr_dds, contrast = c("mode_of_action", "CAR/PXR", "Control"))
carpxr_res <- lfcShrink(carpxr_dds, coef = 2)
#export results-
write.csv(dnadmg_res, "dna_damage_results.csv")
write.csv(ahr_res, "ahr_results.csv")
write.csv(carpxr_res, "carpxr_results.csv")
#export normalized counts
write.csv(counts(dnadmg_dds,normalized=TRUE),'dnadmg_deseq_norm_counts.csv')
write.csv(counts(ahr_dds,normalized=TRUE),'ahr_deseq_norm_counts.csv')
write.csv(counts(carpxr_dds,normalized=TRUE),'carpxr_deseq_norm_counts.csv')
# convert to tibbles-
dnadmg_tib <- as_tibble(dnadmg_res, rownames= NA) %>% rownames_to_column("genes")
ahr_tib <- as_tibble(ahr_res, rownames= NA) %>% rownames_to_column("genes")
carpxr_tib <- as_tibble(carpxr_res, rownames= NA) %>% rownames_to_column("genes")
#filter to only significant genes
dnadmg_tib <- filter(dnadmg_tib, padj <=0.05)
ahr_tib <- filter(ahr_tib, padj <=0.05)
carpxr_tib <- filter(carpxr_tib, padj <= 0.05)
#plot histograms
ggplot(dnadmg_tib, aes(log2FoldChange))+geom_histogram(fill='cadetblue', color='black')+theme_bw()+annotate("text", x=-8, y=16, label="A")
ggplot(ahr_tib, aes(log2FoldChange))+geom_histogram(fill='forestgreen', color='black')+theme_bw()+annotate("text", x=-10, y=400, label="B")
ggplot(carpxr_tib, aes(log2FoldChange))+geom_histogram(fill='coral', color='black')+theme_bw()+annotate("text", x=-9, y=1000, label="C")
#plot scatter plots
ggplot(dnadmg_tib, aes(log2FoldChange, -log10(padj)))+geom_point()+theme_bw()+annotate("text", x=-8, y=60, label="A")
ggplot(ahr_tib, aes(log2FoldChange, -log10(padj)))+geom_point()+theme_bw()+annotate("text", x=-10, y=60, label="B")
ggplot(carpxr_tib, aes(log2FoldChange, -log10(padj)))+geom_point()+theme_bw()+annotate("text", x=-9, y=150, label="C")
