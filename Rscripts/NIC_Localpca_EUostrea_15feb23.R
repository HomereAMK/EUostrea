### The BEGINNING ~~~~~
##
# ~ Plots Local_pca | First written by Nicolas Lou with later modifications by Homère Alves Monteiro

#### Cleans the environment #### 
rm(list=ls())
pacman::p_load(optparse, tidyverse, plyr, RColorBrewer, extrafont, ggforce, ggstar, RcppCNPy)

#### Assemble local_pca input #### 
lg_list <- read_lines("~/Desktop/Scripts/Data/Local_PCA/FullDataset/lg_list.txt")
i=1
for (lg in lg_list){
  pca_summary_temp <- read_tsv(paste0("~/Desktop/Scripts/Data/LocalPCA_EUostrea/pca_summary14feb23_1000snp_", lg, "_2pc.tsv"), col_names = F)
  snp_position_temp <- read_tsv(paste0("~/Desktop/Scripts/Data/LocalPCA_EUostrea/snp_position14feb23_1000snp_", lg, "_2pc.tsv"), col_names = F)
  if (i == 1) {
    pca_summary <- pca_summary_temp
    snp_position <- snp_position_temp
  } else {
    pca_summary <- bind_rows(pca_summary, pca_summary_temp)
    snp_position <- bind_rows(snp_position, snp_position_temp)
  }
  i <- i+1
}

write_tsv(pca_summary, "~/Desktop/Scripts/Data/LocalPCA_EUostrea/pca_summary14feb23_1000snp_2pc.tsv", col_names = F)
write_tsv(snp_position, "~/Desktop/Scripts/Data/LocalPCA_EUostrea/snp_position14feb23_1000snp_2pc.tsv", col_names = F)

#### Run pc_dist and assemble the output #### 
#install.packages("data.table")
#devtools::install_github("petrelharp/local_pca/lostruct")
library(lostruct)
# Read the input
pca_summary <- read_tsv("~/Desktop/Scripts/Data/LocalPCA_EUostrea/pca_summary14feb23_1000snp_2pc.tsv", col_names = F) 
# Run pc_dist with pca_summary
pca_summary <- as.matrix(pca_summary)
attr(pca_summary, 'npc') <- 2
dist <- pc_dist(pca_summary) #Is taking a long ass time.
write_tsv(as.data.frame(dist), '~/Desktop/Scripts/Data/LocalPCA_EUostrea/run_pc_dist_1000snp_2pc.14feb23.tsv', col_names = F)

#### Analysis with all windows ####
## Read the SNP position file
snp_position <- read_tsv("~/Desktop/Scripts/Data/LocalPCA_EUostrea/snp_position14feb23_1000snp_2pc.tsv", col_names = F) %>%
  transmute(lg = sub("\\_.*", "", X1), 
            start = as.integer(sub(".*\\_", "", X1)), 
            stop = as.integer(sub(".*\\_", "", X2)), 
            center = (start+stop)/2)
## Read in the window distance matrix
dist <- read_tsv("~/Desktop/Scripts/Data/LocalPCA_EUostrea/run_pc_dist_1000snp_2pc.14feb23.tsv", col_names = F) %>%
  as.matrix()
## Do MDS with the distance matrix
mds <- cmdscale(as.dist(dist), k=10, eig = T)
proportion_variance <- round(mds$eig/sum(mds$eig)*100,2)
mds <- mds$points %>%
   as.data.frame() %>%
   set_names(paste0("dist_", 1:10)) %>%
   bind_cols(snp_position, .)
#reorder scaffold
# Reorders Population ~
mds$lg <- factor(mds$lg, ordered = T,
                            levels = c("scaffold1", "scaffold2", "scaffold3", "scaffold4", "scaffold5", "scaffold6", "scaffold7",
                                       "scaffold8", "scaffold9", "scaffold10"))
#pdf("~/Desktop/Scripts/Flat_oysters/04_local_R/03_results/Local_pca/DistWinMDS1-10_10k_localpca.feb22.pdf", width = 12, height = 8)
Local_pca <- mds %>%
  gather(key = axis, value = value, 5:14) %>%
  mutate(axis=fct_relevel(axis, paste0("dist_", 1:10))) %>%
  ggplot(aes(x=center/10^6, y=value, color=value)) +
  geom_point(size=0.5) +
  facet_grid(axis~lg, scales="free_x", space="free_x") +
  theme_bw() +
  xlab("position (Mbp)") +
  ylab("principal coordinate value") +
  theme(panel.spacing = unit(0.1, "lines"),
        axis.title.x=element_text(),
        legend.position="none",
        text = element_text(size=10),
        axis.text = element_text(size=6),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  theme(legend.key = element_blank()) +
  theme(legend.title=element_blank()) +
  theme(axis.title.x = element_text(size = 10, color="#000000", face="bold", margin = margin(t = 20, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size = 10, color="#000000", face="bold", margin = margin(t = 0, r = 20, b = 0, l = 0))) +
  theme(legend.text=element_text(size=11)) +
  theme(panel.background = element_rect(fill = '#FAFAFA')) +
  theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank()) +
  theme(axis.line = element_line(colour = "#000000", size = 0.3)) +
  theme(panel.border = element_blank())

dev.off()
ggsave(plot=Local_pca, file = "~/Desktop/Scripts/EUostrea/Figures/Local_PCA/DistWinMDS1-10_1000snps_2pc_localpca.feb23.pdf", device = "pdf",  scale = 1.1, width = 12, height = 8, dpi = 600)

#mds
mds1vs2<-ggplot(mds, aes(x=dist_1, y=dist_2, color=lg,  size=1, alpha=0.8)) +
  geom_point() +
  #scale_shape_manual(values = c(rep(c(15,16,17,18),7), 15, 16)) +
  scale_color_manual(values = MetBrewer::met.brewer("Redon", n = 10, type = "discrete")) +
  scale_size(guide = "none") +
  scale_alpha(guide = "none") +
  xlab(paste0("PCo1 (", proportion_variance[1], "%)")) +
  ylab(paste0("PCo2 (", proportion_variance[2], "%)")) +
  labs(color = "Scaffolds") +  # Renaming the color scale in the legend title
  theme(legend.key = element_blank()) +
  theme(legend.title=element_blank()) +
  theme(axis.title.x = element_text(size = 10, color="#000000", face="bold", margin = margin(t = 20, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size = 10, color="#000000", face="bold", margin = margin(t = 0, r = 20, b = 0, l = 0))) +
  theme(legend.text=element_text(size=11)) +
  theme(panel.background = element_rect(fill = '#FAFAFA')) +
  theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank()) +
  theme(axis.line = element_line(colour = "#000000", size = 0.3)) +
  theme(panel.border = element_blank()) +
  cowplot::theme_cowplot()
ggsave(mds1vs2, file = "~/Desktop/Scripts/EUostrea/Figures/Local_PCA/MDS1vs2--DistWin_1000snps_2pc_localpca.feb24.pdf", device = "pdf",  scale = 1.1, width = 12, height = 8, dpi = 600)

mds3vs4<-ggplot(mds, aes(x=dist_3, y=dist_4, color=lg, shape=lg)) +
  geom_point() +
  scale_shape_manual(values = c(rep(c(15,16,17,18),7), 15, 16)) +
  xlab(paste0("PCo3 (", proportion_variance[3], "%)")) +
  ylab(paste0("PCo4 (", proportion_variance[4], "%)")) +
  cowplot::theme_cowplot()
ggsave(mds3vs4, file = "~/Desktop/Scripts/Flat_oysters/04_local_R/03_results/Local_pca/MDS3vs4--DistWin_10k_localpca.feb22.pdf", device = "pdf",  scale = 1.1, width = 12, height = 8, dpi = 600)

cowplot::theme_cowplot()
mds1vs3<-ggplot(mds, aes(x=dist_1, y=dist_3, color=lg, shape=lg)) +
  geom_point() +
  scale_shape_manual(values = c(rep(c(15,16,17,18),7), 15, 16)) +
  xlab(paste0("PCo1 (", proportion_variance[1], "%)")) +
  ylab(paste0("PCo3 (", proportion_variance[3], "%)")) +
  cowplot::theme_cowplot()
ggsave(mds1vs3, file = "~/Desktop/Scripts/Flat_oysters/04_local_R/03_results/Local_pca/MDS1vs3--DistWin_10k_localpca.feb22.pdf", device = "pdf",  scale = 1.1, width = 12, height = 8, dpi = 600)
ggplot(mds, aes(x=dist_5, y=dist_6, color=lg, shape=lg)) +
  geom_point() +
  scale_shape_manual(values = c(rep(c(15,16,17,18),7), 15, 16)) +
  xlab(paste0("PCo5 (", proportion_variance[5], "%)")) +
  ylab(paste0("PCo6 (", proportion_variance[6], "%)")) +
  cowplot::theme_cowplot()
ggplot(mds, aes(x=dist_7, y=dist_8, color=lg, shape=lg)) +
  geom_point() +
  scale_shape_manual(values = c(rep(c(15,16,17,18),7), 15, 16)) +
  xlab(paste0("PCo7 (", proportion_variance[7], "%)")) +
  ylab(paste0("PCo8 (", proportion_variance[8], "%)")) +
  cowplot::theme_cowplot()



## Find inversion start and stop positions
LG04_inversion_start <- filter(mds, dist_4 > 0.5, lg=="scaffold4") %>%
  .$start %>%
  min()
LG04_inversion_stop <- filter(mds, dist_4 > 0.5, lg=="scaffold4") %>%
  .$stop %>%
  max()
LG05_inversion_start <- filter(mds, dist_1 > 0.5, lg=="scaffold5") %>%
  .$start %>%
  min()
LG05_inversion_stop <- filter(mds, dist_1 > 0.5, lg=="scaffold5") %>%
  .$stop %>%
  max()
LG08_inversion_start <- filter(mds, dist_3 < -0.4, lg=="scaffold8") %>%
  .$start %>%
  min()
LG08_inversion_stop <- filter(mds, dist_3 < -0.4, lg=="scaffold8") %>%
  .$stop %>%
  max()
## Filter out inversions in pca_summary 
pca_summary <- read_tsv('~/Desktop/Scripts/Data/Local_PCA/FullDataset/DatasetIpca_summary_10000snp_2pc.21feb22.tsv', col_names = F) 
combined_filtered <- bind_cols(snp_position, pca_summary) %>%
  filter(!(lg=="scaffold4" & start> LG04_inversion_start & start < LG04_inversion_stop)) %>%
  filter(!(lg=="scaffold5" & start> LG05_inversion_start & start < LG05_inversion_stop)) %>%
  filter(!(lg=="scaffold8" & start> LG08_inversion_start & start < LG08_inversion_stop))
pca_summary_filtered <-  combined_filtered[,-(1:4)] %>% 
  as.matrix()
attr(pca_summary_filtered, 'npc') <- 2
snp_position_filtered <-  combined_filtered[,(1:4)]
write_tsv(snp_position_filtered, "~/Desktop/Scripts/Data/Local_PCA/FullDataset/DatasetIsnp_position_filtered_10000snp_2pc.24feb22.tsv", col_names = F)

## Run pc_dist with pca_summary_filtered
dist_filtered <- pc_dist(pca_summary_filtered)
write_tsv(as.data.frame(dist_filtered), "~/Desktop/Scripts/Data/Local_PCA/FullDataset/DatasetIwindow_dist_filtered_10000snp_2pc.24feb22.tsv", col_names = F)

###### Analysis with inversions filtered out ########
## Read the SNP position file
snp_position_filtered <- read_tsv("~/Desktop/Scripts/Data/Local_PCA/FullDataset/DatasetIsnp_position_filtered_10000snp_2pc.24feb22.tsv", col_names = F) %>%
  transmute(lg = X1, 
            start = X2, 
            stop = X3, 
            center = X4)
## Read in the window distance matrix
dist_filtered<- read_tsv("~/Desktop/Scripts/Data/Local_PCA/FullDataset/DatasetIwindow_dist_filtered_10000snp_2pc.24feb22.tsv", col_names = F) %>%
  as.matrix()
## Do MDS with the distance matrix
mds_filtered <- cmdscale(as.dist(dist_filtered), k=10, eig = T)
proportion_variance_filtered <- round(mds_filtered$eig/sum(mds_filtered$eig)*100,2)
mds_filtered <- mds_filtered$points %>%
  as.data.frame() %>%
  set_names(paste0("dist_", 1:10)) %>%
  bind_cols(snp_position_filtered, .)
mds_filtered$lg <- factor(mds_filtered$lg, ordered = T,
                 levels = c("scaffold1", "scaffold2", "scaffold3", "scaffold4", "scaffold5", "scaffold6", "scaffold7",
                            "scaffold8", "scaffold9", "scaffold10"))
pdf("~/Desktop/Scripts/Flat_oysters/04_local_R/03_results/Local_pca/DistWinMDS1-10_10k_localpca--FilteredInversions.feb22.pdf", width = 12, height = 8)
mds_filtered %>%
  gather(key = axis, value = value, 5:14) %>%
  mutate(axis=fct_relevel(axis, paste0("dist_", 1:10))) %>%
  ggplot(aes(x=center/10^6, y=value, color=axis)) +
  geom_point(size=0.5) +
  facet_grid(axis~lg, scales="free_x", space="free_x") +
  theme_bw() +
  xlab("position (Mbp)") +
  ylab("principal coordinate value") +
  theme(panel.spacing = unit(0.1, "lines"),
        axis.title.x=element_text(),
        legend.position="none",
        text = element_text(size=10),
        axis.text = element_text(size=6),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) 
dev.off()

mdsFIL1vs2<-ggplot(mds_filtered, aes(x=dist_1, y=dist_2, color=lg, shape=lg)) +
  geom_point() +
  scale_shape_manual(values = c(rep(c(15,16,17,18),7), 15, 16)) +
  xlab(paste0("PCo1 (", proportion_variance_filtered[1], "%)")) +
  ylab(paste0("PCo2 (", proportion_variance_filtered[2], "%)")) +
  cowplot::theme_cowplot()
ggsave(mdsFIL1vs2, file = "~/Desktop/Scripts/Flat_oysters/04_local_R/03_results/Local_pca/mds1vs2--DistWin_10k_localpca--FilteredInversions.feb22.pdf", device = "pdf",  scale = 1.1, width = 12, height = 8, dpi = 600)
mdsFIL3vs4<-ggplot(mds_filtered, aes(x=dist_3, y=dist_4, color=lg, shape=lg)) +
  geom_point() +
  scale_shape_manual(values = c(rep(c(15,16,17,18),7), 15, 16)) +
  xlab(paste0("PCo3 (", proportion_variance_filtered[3], "%)")) +
  ylab(paste0("PCo4 (", proportion_variance_filtered[4], "%)")) +
  cowplot::theme_cowplot()
ggsave(mdsFIL3vs4, file = "~/Desktop/Scripts/Flat_oysters/04_local_R/03_results/Local_pca/mds1vs2--DistWin_10k_localpca--FilteredInversions.feb22.pdf", device = "pdf",  scale = 1.1, width = 12, height = 8, dpi = 600)


###### Check the PCA pattern with local PCA outliers including the inversions ######
#Subset the list of all beagle files to only include outlier windows
snp_position <- read_tsv("~/Desktop/Scripts/Data/Local_PCA/Dataset_Isnp_position_10000snp_2pc.tsv", col_names = c("start", "stop"))
write_tsv(mds, "~/Desktop/Scripts/Data/Local_PCA/trialMDS_datasetI_3chr.tsv", col_names = T)
beagle_list <- read_tsv("~/Desktop/Scripts/Data/Local_PCA/Trial//DatasetI_beagle_list.txt", col_names = F)
mds_table_filtered_joined <- left_join(snp_position, mds_table_filtered, by=c("start", "stop"))

## dist1 outliers
dist1_outlier_index <- which(mds$dist_1 < -0.15)
beagle_list_dist1_outlier <- beagle_list[dist1_outlier_index,]
write_tsv(beagle_list_dist1_outlier, "~/Desktop/Scripts/Data/Local_PCA/DatasetI_trial_beagle_list_dist1_outlier.txt", col_names = F)
## dist2 outliers
dist2_outlier_index <- which(mds$dist_2 > -0.25)
beagle_list_dist2_outlier <- beagle_list[dist2_outlier_index,]
write_tsv(beagle_list_dist2_outlier, "~/Desktop/Scripts/Data/Local_PCA/DatasetI_trialbeagle_list_dist2_outlier.txt", col_names = F)




