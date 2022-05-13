#import libraries
library(tidyverse)
library(ggplot2)
#import the first file
df <- read.delim("SRR1177981.txt", skip=1)
#create character vector of remaining file names to loop over
files <- c("SRR1177982.txt", "SRR1177983.txt", "SRR1178008.txt", "SRR1178009.txt", "SRR1178010.txt", "SRR1178014.txt", "SRR1178021.txt", "SRR1178047.txt")
#initialize counts df by selecting required columns
counts_tib <- select(df, Geneid, last_col()) %>% rename(gene_id = Geneid, SRR1177981 = colnames(df)[7])
#loop over files, joining using gene ID
for(f in files){
  temp <- read.delim(f, skip = 1) %>% select(Geneid, last_col())
  #temp <- rename(temp, f = colnames(temp)[7])
  counts_tib <- counts_tib %>%
    full_join(temp, by = c("gene_id"="Geneid"))
}
#rename columns to just sample names
names(counts_tib) <- sub("X.projectnb.bf528.users.tinman_2022.project_3.star_results.", "", names(counts_tib))
names(counts_tib) <- sub("Aligned.sortedByCoord.out.bam", "", names(counts_tib))
#write out csv file
write.csv(counts_tib, "counts.csv")
#create plotting df by converting data to long form
plot_tib <- pivot_longer(counts_tib, cols=starts_with("SRR"), names_to = "Sample", values_to = "Counts")
#plot box plot
ggplot(plot_tib, aes(Sample, Counts))+
  geom_boxplot()+
  scale_y_log10()+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90))
