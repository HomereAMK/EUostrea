# Gene search on the gff file ~



# Cleans the environment ~ 
rm(list=ls())

# Loads required packages ~
pacman::p_load(tidyverse, extrafont, lemon, data.table, MetBrewer, cowplot, GenomicRanges, IRanges, RcppCNPy, ggrepel, seqinr, GenomicFeatures)


# Only run this once
gff_Genome_MOR <- read_tsv("/Users/homere/Desktop/IGV/gff_oedulis.tsv", skip = 1, col_names = F) %>%
  filter(X3=="gene")

# Specify colnames
colnames(gff_Genome_MOR) <- c("seqid","source", "feature", "start", "end",  "score", "strand", "phase", "attributes")

# Filter rows to only include scaffold4 and start between 1 and 10**7 invReg
InvReg4 <- gff_Genome_MOR %>%
  filter(seqid == "scaffold4" & start >= 1 & end <= 10**7)

# Filter rows to only include scaffold5 and start between 1 and 20000000 invReg
InvReg5 <- gff_Genome_MOR %>%
  filter(seqid == "scaffold5" & start >= 1 & end <= 20000000)

# Filter rows to only include scaffold8 and start between 32000000 and 58502465 invReg
InvReg8 <- gff_Genome_MOR %>%
  filter(seqid == "scaffold8" & start >= 32000000 & end <= 58502465)

