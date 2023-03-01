## The BEGINNING ~~~~~
##
#  FastaFromAlleleCounts Mtgenome phylogeography | First written by Nicolas Lou with later modifications by Hom?re J. Alves Monteiro

# Cleans the environment  
rm(list=ls())
library(tidyverse)
library(cowplot)
library(ape)
library(knitr)
library(pegas)
library(RColorBrewer)

source("~/Desktop/Scripts/Flat_oysters/04_local_R/00_scripts/NIC_mtgenome_functions.R")

# Wrangling the BAM list
bam_list_file <- "~/Desktop/Scripts/Data/MtGenome_EUostrea/Mt_HapNetwork_10feb23_15pops_bamlist.txt"
bams <- read.table(bam_list_file)[, 1]
bams <- gsub(".bam", "", bams)
bams <- gsub(".+/", "", bams)

# Defining the population variable
pop <- ifelse(grepl("^Lurida", bams), substr(bams, 1, 6), substr(bams, 1, 4))
# Defining the individual identifier
ind <- ifelse(grepl("^Lurida", bams), substr(bams, 1, 9), substr(bams, 1, 7))



# Analyze depth count across all invididuals
depth_count <- read_tsv("~/Desktop/Scripts/Data/MtGenome_EUostrea/Mt_HapNetwork_10feb23.depth_counts.counts.gz") %>%
  dplyr::select(-144) %>%
  t() %>%
  as_tibble() %>%
  bind_cols(ind=ind, population=pop, .) %>%
  pivot_longer(3:16356, names_to = "position", values_to="depth") %>%
  mutate(position=as.numeric(substring(position, 2)))

# Average depth per individual
set.seed(1)
ind_average_depth <- group_by(depth_count, ind, population) %>%
  summarise(average_depth=mean(depth))
ggplot(ind_average_depth, mapping=aes(x=population, y=average_depth, group=population)) +
  geom_boxplot(outlier.alpha = 0) +
  geom_jitter() +
  theme_cowplot() +
  coord_flip()

## samples with the lowest depth
group_by(depth_count, ind, population) %>%
  summarise(average_depth=mean(depth, na.rm = TRUE)) %>%
  ungroup() %>%
  arrange(average_depth) %>%
  head(n=20) %>%
  kable()

#Average depth per position
group_by(depth_count, position) %>%
  summarise(average_depth=mean(depth, na.rm = TRUE)) %>%
  ggplot(aes(x=position, y=average_depth)) +
  geom_line() +
  theme_cowplot()


#Generate fasta files from allele counts
# minimum depth 4, minimum major allele frequency 0.75
concensus <- convert_count_to_concensus(x="~/Desktop/Scripts/Data/MtGenome_EUostrea/Mt_HapNetwork_10feb23.allele_counts.counts.gz", min_depth=4, min_maf=0.8)
ind_label <- ind
concensus_to_fasta(concensus, ind_label, "~/Desktop/Scripts/Data/MtGenome_EUostrea/Mt_HapNetwork_10feb23_mindepth4_minmaf80.fasta")
concensus_wide <- as_tibble(t(concensus))
colnames(concensus_wide) <- ind_label
concensus_wide <- mutate(concensus_wide, position=1:nrow(concensus_wide))
concensus_long <- pivot_longer(concensus_wide, cols = 1:length(ind_label), names_to = "sample_id", values_to = "allele")
write_tsv(concensus_long, "~/Desktop/Scripts/Data/MtGenome_EUostrea/Mt_HapNetwork_10feb23_mindepth4_minmaf80.tsv")


#Convert fasta file to nexus for PopArt (relaxed filter)
fasta <- read.dna("~/Desktop/Scripts/Data/MtGenome_EUostrea/Mt_HapNetwork_10feb23_mindepth4_minmaf80.fasta", format="fasta")
concensus_long <- read_tsv("~/Desktop/Scripts/Data/MtGenome_EUostrea/Mt_HapNetwork_10feb23_mindepth4_minmaf80.tsv")
## Filter out invariant sites
snp <- concensus_long %>%
  filter(allele!="N") %>%
  group_by(position) %>%
  filter(length(unique(allele))>1) %>%
  .$position %>%
  unique()
## Filter out singleton mutations
non_singleton_snp <- concensus_long %>%
  filter(position %in% snp, allele!="N") %>%
  count(position, allele) %>%
  group_by(position) %>%
  slice_max(order_by=n, n=2, with_ties=FALSE) %>%
  slice_min(order_by=n, n=1, with_ties=FALSE) %>%
  filter(n>1) %>%
  .$position
## Get a list of sites with average depth lower than 30
low_depth_site <- group_by(depth_count, position) %>%
  summarise(average_depth=mean(depth, na.rm = TRUE)) %>%
  filter(average_depth < 30) %>%
  .$position
## Exclude sites with low depth from the non-singleton SNP list
site_to_include <- setdiff(non_singleton_snp, low_depth_site)
## Get a list of sample with at least 2 Ns at the selected loci
samples_to_exclude <- concensus_long %>%
  filter(position %in% site_to_include, allele=="N") %>%
  group_by(sample_id) %>%
  summarise(position_missing=n()) %>%
  arrange(desc(position_missing)) %>%
  filter(position_missing > 1) %>%
  .$sample_id
## Subset the sample table
sample_table<- data.frame(bams, ind_label, pop)
network_sample_table <- sample_table %>%
  filter(!ind_label %in% samples_to_exclude)
## Subset the fasta file 
network_fasta <- fasta[(! sample_table$ind_label %in% samples_to_exclude), site_to_include]

# Convert to nexus format
write.nexus.data(network_fasta, file="~/Desktop/Scripts/Data/MtGenome_EUostrea/Mt_HapNetwork_10feb23_mindepth4_minmaf80.nex", format="dna")

#Generate trait tables for popart

sample_table %>%
  dplyr::filter(!(ind_label %in% samples_to_exclude)) %>%
  dplyr::select(ind_label, pop) %>%
  mutate(temp=1) %>%
  spread(key = pop, value = temp) %>%
  mutate_all(~replace(., is.na(.), 0)) %>%
  write_csv("~/Desktop/Scripts/Data/MtGenome_EUostrea/Mt_HapNetwork_10feb23_mindepth4_minmaf80_traitPpopart.csv")






concensus_to_fasta(concensus, ind_label, 
                   "~/Desktop/Scripts/Flat_oysters/04_local_R/03_results/bam_list_realigned_mtgenome_sorted_filtered_minq20_mindepth2_minmaf75.fasta")
tail(concensus)

#### minimum depth 2, minimum major allele frequency 0.75 #### 
concensus <- convert_count_to_concensus(x="~/Desktop/", min_depth=2, min_maf=0.75)
summary(concensus)


#### minimum depth 4, minimum major allele frequency 0.75 ####
concensus <- convert_count_to_concensus(x="~/Desktop/Scripts/Flat_oysters/04_local_R/02_data/bam_list_realigned_mtgenome_sorted_filtered_minq20.allele.counts.gz", min_depth=4, min_maf=0.75)
concensus_to_fasta(concensus, ind_label, "~/Desktop/Scripts/Flat_oysters/04_local_R/03_results/bam_list_realigned_mtgenome_sorted_filtered_minq20_mindepth4_minmaf75.fasta")

#Convert fasta file to nexus for PopArt (relaxed filter)
fasta <- read.dna("~/Desktop/Scripts/Flat_oysters/04_local_R/03_results/bam_list_realigned_mtgenome_sorted_filtered_minq20_mindepth2_minmaf75.fasta", format="fasta")
# get depth count
depth_count <- read_tsv("~/Desktop/Scripts/Flat_oysters/04_local_R/02_data/bam_list_realigned_mtgenome_sorted_filtered_minq20.depth.counts.gz") %>%
  dplyr::select(-X481) %>%
  t() %>%
  as_tibble() %>%
  bind_cols(ind=ind_label, population=pop, .) %>%
  gather(key = position, value="depth", 3:16356) %>%
  mutate(position=as.numeric(substring(position, 2)))

#### exclude samples depth>2 ####
samples_to_exclude <- filter(depth_count, depth>2) %>%
  count(ind) %>%
  mutate(position_missing = 16356-n) %>%
  arrange(desc(position_missing)) %>%
  filter(position_missing > 0.3*16356) %>%
  .$ind
fasta <- subset(fasta, subset = !(attr(fasta,"dimnames")[[1]] %in% samples_to_exclude))
fasta
write.nexus.data(fasta, file="~/Desktop/Scripts/Flat_oysters/04_local_R/03_results/bam_list_realigned_mtgenome_sorted_filtered_minq20_mindepth2_minmaf75.nex", format="dna")


#### exclude samples depth>10 ####
samples_to_exclude2 <- filter(depth_count, depth>10) %>%
  count(ind) %>%
  mutate(position_missing = 16356-n) %>%
  arrange(desc(position_missing)) %>%
  filter(position_missing > 0.3*16356) %>%
  .$ind
fasta <- subset(fasta, subset = !(attr(fasta,"dimnames")[[1]] %in% samples_to_exclude2))
fasta
write.nexus.data(fasta, file="~/Desktop/Scripts/Flat_oysters/04_local_R/03_results/bam_list_realigned_mtgenome_sorted_filtered_minq20_mindepth10_minmaf75.nex", format="dna")



###### Average depth per individual ####
set.seed(1)
dfyes <- depth_count %>%
  filter(! str_detect(population, c("USAM|MORL")))
ind_average_depth <- group_by(dfyes, ind, population) %>%
  summarise(average_depth=mean(depth))
ggplot(ind_average_depth, mapping=aes(x=population, y=average_depth, group=population)) +
  geom_boxplot() +
  geom_jitter() +
  theme_cowplot() +
  coord_flip()
last_plot()
ggsave(filename = "Desktop/Scripts/Flat_oysters/04_local_R/03_results/Average_depthMt_pop.pdf") #change path


#### Average depth per position ####
group_by(depth_count, position) %>%
  summarise(average_depth=mean(depth)) %>%
  ggplot(aes(x=position, y=average_depth)) +
  geom_line() +
  theme_cowplot()
last_plot()
ggsave(filename = "Desktop/Scripts/Flat_oysters/04_local_R/03_results/Average_depthMt_perPosition.pdf") #change path

sample_table<- data.frame(bams, ind_label, pop)


####traits table for popart####
sample_table %>%
  dplyr::filter(!(ind_label %in% samples_to_exclude)) %>%
  dplyr::select(ind_label, pop) %>%
  mutate(temp=1) %>%
  spread(key = pop, value = temp) %>%
  mutate_all(~replace(., is.na(.), 0)) %>%
  write_csv("~/Desktop/Scripts/Flat_oysters/04_local_R/03_results/bam_list_realigned_mtgenome_sorted_filtered_population_trait_table_for_popart.csv")

####traits table for popart mindepth>10####
sample_table %>%
  dplyr::filter(!(ind_label %in% samples_to_exclude2)) %>%
  dplyr::select(ind_label, pop) %>%
  mutate(temp=1) %>%
  spread(key = pop, value = temp) %>%
  mutate_all(~replace(., is.na(.), 0)) %>%
  write_csv("~/Desktop/Scripts/Flat_oysters/04_local_R/03_results/bam_list_realigned_mtgenome_sorted_filtered_population_trait_table_for_popart_mindepth10.csv")
