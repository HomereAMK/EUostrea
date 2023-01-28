rm(list=ls(all.names = TRUE))

library(tidyverse)
basedir="~/Desktop"  

chr<- c("scaffold1:", "scaffold2:", "scaffold3:", "scaffold4:", "scaffold5:",
        "scaffold6:", "scaffold7:", "scaffold8:", "scaffold9:", "scaffold10:")

pruned_position <- read_lines(paste0(basedir, "/Scripts/Data/LDpruning/EUostrea/Dec22.LD.AllCHRs.min_weight0.5.pruneIN")) %>% 
  str_remove("scaffold10:") %>% str_remove("scaffold1:") %>% str_remove("scaffold2:") %>% 
  str_remove("scaffold3:") %>% str_remove("scaffold4:") %>% str_remove("scaffold5:") %>%
  str_remove("scaffold6:") %>% str_remove("scaffold7:") %>% str_remove("scaffold8:") %>%
  str_remove("scaffold9:") %>% 
  as.integer()
pruned_position


pruned_snp_list <- read_tsv(paste0(basedir, "/Scripts/")) %>%
  dplyr::select(1:4) %>%
  filter(position %in% pruned_position)

write_tsv(pruned_snp_list, paste0(basedir, "/Scripts/"), col_names = F)


######  LD pruned SNPs with global MAF >= 0.05 ########
sample_table <- read_tsv("", col_names = T)
pops <- unique(sample_table$group)
global_maf <- read_tsv("~/Desktop/Scripts/Data/LDpruning/EUostrea/Jan23_NICdepthfilters_het_angsd.mafs.gz", skip=1, col_names = FALSE)
ld_pruned_snp_list <- pruned_snp_list %>%
  dplyr::select(1:2) %>%
  left_join(global_maf, by = c("X1", "X2")) %>%
  filter(X6>=0.05) 

mafs_filtered <- mutate(mafs, 
                        keep = knownEM >= 0.01,
                        row_number = ifelse(keep, cumsum(keep), -1),
                        #keep = ifelse(row_number%%3!=0, F, keep),
                        keep = ifelse(!str_detect(chromo, "scaffold"), F, keep)) %>%
  dplyr::select(-row_number)
## Number of SNPs left
filter(mafs_filtered, keep==T) %>% nrow()

mafs_filtered %>%
  filter(keep==T) %>%
  ggplot(aes(x=position, y=chromo, fill=chromo)) +
  ggridges::geom_density_ridges(alpha=0.5) +
  scale_fill_viridis_d() +
  theme_cowplot() +
  theme(legend.position = "none")


##### 
