### The BEGINNING ~~~~~
##
# ~ Plots Population structure results | First written by Nicolas R, Lou  with later modifications by  Homère J. Alves Monteiro


# Cleans the environment ~ 
rm(list=ls())

# Sets working directory ~
setwd("~/Desktop/Scripts/")
#install.packages("RcppCNPy")
# Loads required packages ~
pacman::p_load(vegan, tidyverse, RcppCNPy, pheatmap, extrafont, ggforce, ggrepel, ggstar, np, reticulate, cowplot)


## Plotting CovMat and IbsMat from ANGSD
cov_mat <- as.matrix(read.table("~/Desktop/Scripts/Data/PopStruct_EUostrea/28jan23_prunedLDminweight0.5_PopStruct.covMat")) 
annot <- read.table("../Scripts/EUostrea/01_infofiles/bamlist_EUostrea.annot", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
source("~/Desktop/Scripts/Flat_oysters/04_local_R/00_scripts/individual_pca_functions_hjam_dec22.R")
# Reorders Population ~
annot$V2 <- factor(annot$V2, ordered = T,
                   levels = c("MOLU", "ZECE", "CRES",
                              "ORIS","CORS", "PONT",  "RIAE",
                              "MORL",
                              "USAM",
                              "TOLL", "COLN", "BARR",
                              "TRAL", "CLEW",
                              "RYAN",
                              "GREV", "WADD",
                              "NISS","LOGS","VENO", "HALS", "THIS",
                              "KALV", "HYPP",
                              "LANG", "BUNN", "DOLV", "HAUG", "HAFR",
                              "INNE","VAGS", "AGAB", "OSTR"))
#Plot genome-wide PCA with the covMat matrix
Pca_axis1_2<-PCA(cov_mat, annot$V1, annot$V2, 1, 2, show.ellipse = FALSE, show.label = FALSE)
Pca_axis2_3<-PCA(cov_mat, annot$V1, annot$V2, 2, 3, show.ellipse = FALSE, show.label = FALSE)

#ggsave(, file = "~/Desktop/Scripts/EUostrea/Figures/PopulationStructure/.pdf", device = cairo_pdf, scale = 1.1, width = 12, height = 8, dpi = 300)
#dev.off()

#Plot genome-wide PCoA with the ibsMat matrix
ibs_mat <- read_tsv("~/Desktop/Scripts/Data/PopStruct_EUostrea/28jan23_prunedLDminweight0.5_PopStruct.ibsMat", col_names = F) %>% 
  dplyr::select(1:nrow(.)) %>%
  as.matrix()
PCoA(ibs_mat, annot$V1, annot$V2, 153, 1, 2, show.ellipse = F)


## Plotting CovMat and IbsMat from Pcangsd
#LD pruned SNP list
genome_cov <- read_delim("Data/PopStruct_EUostrea/30jan23_prunedLDminweight0.5_PopStruct_pcangsd.cov", delim = " ", col_names = FALSE) %>%
  as.matrix()
PCA(genome_cov, annot$V1, annot$V2, 1, 2, show.ellipse = F)

pca_table_ld_pruned <- pca_table
pca_table_joined <- pca_table_ld_pruned %>%
  dplyr::select(1:10) %>%
  left_join(annot, by=c("individual"="V1", "population"="V2"))
pca_table_joined %>%
  ggplot(aes(x=PC1, y=PC2)) +
  geom_point(data=dplyr::select(pca_table_joined, -population), size = 0.1, color="grey") +
  geom_point(size=1, mapping = aes(color=population)) +
  facet_wrap(~population) +
  theme_cowplot() +
  theme(axis.text = element_blank())



## Plot selection scan result with PCAngsd (e not specified) 
# .npy and .sites files from PCangsd
genome_selection <- read_tsv("~/Desktop/Scripts/Data/PopStruct_EUostrea/30jan23_prunedLDminweight0.5_PopStruct_pcangsd.selection.tsv", col_names = F)
genome_selection_sites <- read_tsv("Data/PopStruct_EUostrea/LDprunedlist_rightmafs_AllCHRs.min_weight0.5_23jan23", col_names = c("lg", "pos", "major", "minor")) %>%
  bind_cols(read_table("Data/PopStruct_EUostrea/30jan23_prunedLDminweight0.5_PopStruct_pcangsd.sites", col_names = "keep")) %>%
  filter(keep==1) %>%
  mutate(chi_squared=genome_selection[,1]) %>%
  dplyr::select(lg, pos, chi_squared) %>%
  mutate(neg_log_p_value=-log(1-pchisq(chi_squared$X1, df=1))) 


#reorder chr
genome_selection_sites$lg <- factor(genome_selection_sites$lg , ordered = T,
                   levels = c("scaffold1", "scaffold2", "scaffold3", "scaffold4", "scaffold5",
                              "scaffold6", "scaffold7", "scaffold8", "scaffold9", "scaffold10"))

## Extract sites exceeding  p-value cut off of 0.05 after Bonferroni correction and export them
genome_selection_sites_outliers <- filter(genome_selection_sites, neg_log_p_value > -log(0.05/dim(genome_selection)[1]))
## Plot with one axis
genome_selection_plot <- genome_selection_sites %>%
  ggplot(aes(x=pos/10^6, y=neg_log_p_value)) +
  geom_point(size=0.2, alpha=0.3) +
  geom_smooth(color="purple", se =F) +
  geom_hline(yintercept = -log(0.05/dim(genome_selection)[1]), linetype = "dashed") + # This line marks a p-value cut off of 0.05 after Bonferroni correction
  theme_cowplot() +
  scale_x_continuous(breaks=seq(0, 150, 20)) +
  xlab("position (Mbp)") +
  ylab("-log(p)") +
  facet_grid(~lg, scales="free_x", space="free_x") +
  theme(panel.spacing = unit(0.1, "lines"),
        axis.title.x=element_text(),
        legend.position="none",
        text = element_text(size=10),
        axis.text = element_text(size=6))+
  theme(legend.key = element_blank()) +
  theme(legend.title=element_blank()) +
  theme(axis.title.x = element_text(size = 10, color="#000000", face="bold", margin = margin(t = 20, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size = 10, color="#000000", face="bold", margin = margin(t = 0, r = 20, b = 0, l = 0))) +
  theme(legend.text=element_text(size=11)) +
  theme(panel.background = element_rect(fill = '#FAFAFA')) +
  theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank()) +
  theme(axis.line = element_line(colour = "#000000", size = 0.3)) +
  theme(panel.border = element_blank())

ggsave(filename  = "EUostrea/Figures/Genome_scan/2feb23_PCAngsd_LDpruned_genomescan.pdf", 
       plot=genome_selection_plot, width = 40, height = 10, units = "cm", pointsize = 20, dpi = 300)



## plot with 3 axis
genome_selection_sites <- read_tsv("Data/PopStruct_EUostrea/LDprunedlist_rightmafs_AllCHRs.min_weight0.5_23jan23", col_names = c("lg", "pos", "major", "minor")) %>%
  bind_cols(read_table("Data/PopStruct_EUostrea/30jan23_prunedLDminweight0.5_PopStruct_pcangsd.sites", col_names = "keep")) %>%
  filter(keep==1) %>%
  mutate(pc1=genome_selection[,1],
         pc2=genome_selection[,2],
         pc3=genome_selection[,3]) %>%
  pivot_longer(names_to = "pc", values_to = "chi_squared", cols = 6:8) %>%
  dplyr::select(lg, pos, pc, chi_squared) %>%
  mutate(chi_squared = ifelse(pc == "pc1", chi_squared$X1, 
                              ifelse(pc == "pc2", chi_squared$X2,
                                     ifelse(pc == "pc3", chi_squared$X3, NA)))) %>%
  filter(!is.na(chi_squared)) %>%
  mutate(neg_log_p_value=-log(1-pchisq(chi_squared, df=1))) 

#reorder chr
genome_selection_sites$lg <- factor(genome_selection_sites$lg , ordered = T,
                                    levels = c("scaffold1", "scaffold2", "scaffold3", "scaffold4", "scaffold5",
                                               "scaffold6", "scaffold7", "scaffold8", "scaffold9", "scaffold10"))


## Plot
genome_selection_plot <- genome_selection_sites %>%
  ggplot(aes(x=pos/10^6, y=neg_log_p_value)) +
  geom_point(size=0.2, alpha=0.5) +
  geom_smooth(color="purple", se =F) +
  geom_hline(yintercept = -log(0.05/dim(genome_selection)[1]), linetype = "dashed") + # This line marks a p-value cut off of 0.05 after Bonferroni correction
  theme_cowplot() +
  scale_x_continuous(breaks=seq(0, 100, 10)) +
  coord_cartesian(ylim=c(0, 40), expand = F) +
  xlab("position (Mbp)") +
  ylab("-log(p)") +
  facet_grid(pc~lg, scales="free_x", space="free_x") +
  theme(panel.spacing = unit(0.1, "lines"),
        axis.title.x=element_text(),
        legend.position="none",
        text = element_text(size=10),
        axis.text = element_text(size=6))

ggsave(filename  = "EUostrea/Figures/Genome_scan/2feb23_PCAngsd_LDpruned_genomescan_3pcs.pdf", 
       plot=genome_selection_plot, width = 40, height = 10, units = "cm", pointsize = 20, dpi = 300)

  