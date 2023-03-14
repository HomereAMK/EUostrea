### The BEGINNING ~~~~~
##
# ~ Plots Pcangsd results from --selection flag 
# ~ | First written by Nicolas R, Lou  with later modifications by  Homère J. Alves Monteiro

# Cleans the environment ~ 
rm(list=ls())

# Sets working directory ~
setwd("~/Desktop/Scripts/")

# Loads required packages ~
pacman::p_load(vegan, tidyverse, RcppCNPy, pheatmap, extrafont, ggforce, ggrepel, ggstar, np, reticulate, cowplot)
pacman::p_load(tidyverse, extrafont, lemon, data.table, MetBrewer, cowplot, GenomicRanges, IRanges, RcppCNPy, ggrepel, seqinr, GenomicFeatures)


## plot with 10 axis 1feb23_GlobalVariants
genome_selection <- read_tsv("~/Desktop/Scripts/Data/GenomeScan_EUostrea/1feb23_GlobalVariants_selectionscan_pcangsd_e10.selection.tsv", col_names = F)
genome_selection_sites <- read_tsv("Data/GenomeScan_EUostrea/1feb23_GlobalVariants_selectionscan__minmaf0.05_List.txt", 
                                   col_names = c("lg", "pos", "major", "minor")) %>%
  bind_cols(read_table("Data/GenomeScan_EUostrea/1feb23_GlobalVariants_selectionscan__minmaf0.05_pcangsd_e10.sites", col_names = "keep")) %>%
  filter(keep==1) %>%
  mutate(PC1=genome_selection[,1],
         PC2=genome_selection[,2],
         PC3=genome_selection[,3],
         PC4=genome_selection[,4],
         PC5=genome_selection[,5],
         PC6=genome_selection[,6],
         PC7=genome_selection[,7],
         PC8=genome_selection[,8],
         PC9=genome_selection[,9],
         PC10=genome_selection[,10]) %>%
  pivot_longer(names_to = "pc", values_to = "chi_squared", cols = 6:15) %>%
  dplyr::select(lg, pos, pc, chi_squared) %>%
  mutate(chi_squared = ifelse(pc == "PC1", chi_squared$X1, 
                              ifelse(pc == "PC2", chi_squared$X2,
                                     ifelse(pc == "PC3", chi_squared$X3, 
                                            ifelse(pc == "PC4", chi_squared$X4, 
                                                   ifelse(pc == "PC5", chi_squared$X5,
                                                          ifelse(pc == "PC6", chi_squared$X6,
                                                                 ifelse(pc == "PC7", chi_squared$X7,
                                                                        ifelse(pc == "PC8", chi_squared$X8, 
                                                                               ifelse(pc == "PC9", chi_squared$X9,
                                                                                      ifelse(pc == "PC10", chi_squared$X10,NA))))))))))) %>%
  filter(!is.na(chi_squared)) %>%
  mutate(neg_log_p_value=-log(1-pchisq(chi_squared, df=1))) 

#reorder chr
genome_selection_sites$lg <- factor(genome_selection_sites$lg , ordered = T,
                                    levels = c("scaffold1", "scaffold2", "scaffold3", "scaffold4", "scaffold5",
                                             "scaffold6", "scaffold7", "scaffold8", "scaffold9", "scaffold10"))
genome_selection_sites$pc <- factor(genome_selection_sites$pc , ordered = T,
                                    levels = c("PC1", "PC2", "PC3", "PC4", "PC5",
                                               "PC6", "PC7", "PC8", "PC9", "PC10"))


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


#ggsave(filename  = "EUostrea/Figures/Genome_scan/4feb23_PCAngsd_Global_VC_genomescan_10pcs.pdf", 
 #      plot=genome_selection_plot, width = 30, height = 40, units = "cm", pointsize = 20, dpi = 250)
ggsave(filename  = "EUostrea/Figures/Genome_scan/4feb23_PCAngsd_Global_VC_genomescan_10pcs.png", 
       plot=genome_selection_plot, width = 40, height = 50, units = "cm", pointsize = 20, dpi = 250)
dev.off()

##### Outliers #####
## Extract sites exceeding  p-value cut off of 0.05 after Bonferroni correction and export them
genome_selection_sites_outliers <- filter(genome_selection_sites, neg_log_p_value > -log(0.05/dim(genome_selection)[1]))
head(genome_selection_sites_outliers)
max(genome_selection_sites_outliers$neg_log_p_value)


genome_selection_plot_red <- genome_selection_sites_outliers %>%
  ggplot(aes(x=pos/10^6, y=neg_log_p_value, color="red")) +
  geom_point(size=0.2, alpha=0.5) +
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

  ggsave(filename  = "EUostrea/Figures/Genome_scan/13mar23_PCAngsd_Global_VC_genomescan_10pcs_Outliers.png", 
         plot=genome_selection_plot_red, width = 40, height = 50, units = "cm", pointsize = 20, dpi = 250)
  dev.off()
  


## For PC9, scaffold2 Extract sites exceeding  p-value cut off of 0.05 after Bonferroni correction
genome_selection_sites_outliers_PC9 <- filter(genome_selection_sites_outliers, 
                                              (lg == "scaffold2" | lg == "scaffold10") & pc == "PC9")

## For PC2, scaffold2 Extract sites exceeding  p-value cut off of 0.05 after Bonferroni correction
genome_selection_sites_outliers_PC2_sca2 <- filter(genome_selection_sites_outliers, 
                                              lg == "scaffold2"  & pc == "PC2")
## For PC2, scaffold2 Extract sites exceeding  p-value cut off of 0.05 after Bonferroni correction
genome_selection_sites_outliers_PC2_sca1 <- filter(genome_selection_sites_outliers, 
                                                   lg == "scaffold1"  & pc == "PC2")

# For PC2, scaffold3 Extract sites exceeding  p-value cut off of 0.05 after Bonferroni correction
genome_selection_sites_outliers_PC2_sca3 <- filter(genome_selection_sites_outliers, 
                                                   lg == "scaffold3"  & pc == "PC2")
# Only run this once
gff_Genome_MOR <- read_tsv("/Users/homere/Desktop/IGV/gff_oedulis.tsv", skip = 1, col_names = F) %>%
  filter(X3=="gene")

# Specify colnames
colnames(gff_Genome_MOR) <- c("seqid","source", "feature", "start", "end",  "score", "strand", "phase", "attributes")

## Extract outlier positions genome_selection_sites_outliers_PC2_sca2
genome_selection_sites_outliers_PC2_sca2 <- genome_selection_sites_outliers_PC2_sca2 %>% 
  mutate(pos_bin = cut(pos, breaks = (1:40)*10000000)) %>% 
  group_by(lg, pos_bin) %>%
  summarise(min_pos = min(pos), max_pos = max(pos)) %>%
  ungroup() %>%
  mutate(length = max_pos-min_pos, lg_pos_bin = paste0(lg, pos_bin))
genome_selection_sites_outliers_PC2_sca2

## Extract outlier positions genome_selection_sites_outliers_PC2_sca1 neglogpvalue >23
genome_selection_sites_outliers_PC2_sca1_trans <-  filter(genome_selection_sites_outliers_PC2_sca1,neg_log_p_value>23)  %>% 
  mutate(pos_bin = cut(pos, breaks = (1:40)*10000000)) %>% 
  group_by(lg, pos_bin) %>%
  summarise(min_pos = min(pos), max_pos = max(pos)) %>%
  ungroup() %>%
  mutate(length = max_pos-min_pos, lg_pos_bin = paste0(lg, pos_bin))
genome_selection_sites_outliers_PC2_sca1_trans

## Extract outlier positions genome_selection_sites_outliers_PC2_sca3 neglogpvalue >29
genome_selection_sites_outliers_PC2_sca3_trans <-  filter(genome_selection_sites_outliers_PC2_sca3,neg_log_p_value>29)  %>% 
  mutate(pos_bin = cut(pos, breaks = (1:40)*10000000)) %>% 
  group_by(lg, pos_bin) %>%
  summarise(min_pos = min(pos), max_pos = max(pos)) %>%
  ungroup() %>%
  mutate(length = max_pos-min_pos, lg_pos_bin = paste0(lg, pos_bin))
genome_selection_sites_outliers_PC2_sca3_trans

## Functions ~
#Define a function to find genes close to outlier SNPs
find_nearby_genes <- function(outliers, annotation, distance){
  annotation_range <- with(annotation, GRanges(seqid, IRanges(start, end), strand))
  for (i in seq_len(dim(outliers)[1])){
    outlier <- outliers[i,]
    outlier_lg_pos_bin <- outlier$lg_pos_bin
    outlier_range <- with(outlier, GRanges(lg, IRanges(min_pos - distance, max_pos + distance)))
    overlapping_genes_index <- countOverlaps(annotation_range, outlier_range)>0
    overlapping_genes <- annotation[overlapping_genes_index, ] %>%
      mutate(lg_pos_bin = outlier_lg_pos_bin)
    if (i==1){
      overlapping_genes_final <- overlapping_genes
    } else {
      overlapping_genes_final <- bind_rows(overlapping_genes_final, overlapping_genes)
    }
  }
  return(overlapping_genes_final)
}
#Define a function to plot outlier position and nearby genes
plot_nearby_genes <- function(outlier_sites, nearby_genes){
  ggplot(nearby_genes, aes(ymin=start, ymax=end, y=(start+end)/2, x=0)) +
    geom_linerange(size=1, position=position_jitter(1, 0, 1)) +
    geom_label(aes(label=attributes, x=-0.2), size=5, position=position_jitter(1, 0, 1)) +
    geom_linerange(data=outlier_sites, aes(ymin=min_pos, ymax=max_pos, x=1, y=(min_pos+max_pos)/2), color="red", size=5) +
    facet_wrap(~lg_pos_bin, ncol = 1,scales="free") +
    coord_flip() +
    theme_cowplot()
}

plot_nearby_genes <- function(outlier_sites, nearby_genes){
  lg_labeller <- function(lg_pos_bin){
    return(str_extract(lg_pos_bin, "[^(]+"))
  }
  mutate(nearby_genes, attributes=ifelse(nchar(attributes>10), str_extract(attributes, "[^ ]+"), attributes), attributes=ifelse(attributes=="unknown", NA, attributes)) %>%
    ggplot() +
    geom_segment(aes(x=start, xend=end, y=strand, yend=strand), size=1, arrow = arrow(length = unit(0.2,"cm"))) +
    geom_segment(aes(x=end, xend=start, y=strand, yend=strand), size=1, arrow = arrow(length = unit(0.2,"cm"))) +
    geom_text(aes(label=attributes, x=(start+end)/2, y=strand), size=4, position = position_nudge(y=0.5)) +
    geom_segment(data=outlier_sites, aes(x=min_pos, xend=max_pos, y="outlier", yend="outlier"), color="red", size=2, arrow = arrow(length = unit(0.2,"cm"))) +
    geom_segment(data=outlier_sites, aes(x=max_pos, xend=min_pos, y="outlier", yend="outlier"), color="red", size=2, arrow = arrow(length = unit(0.2,"cm"))) +
    facet_wrap(~lg_pos_bin, ncol = 1, scales="free_x", strip.position="right", labeller = labeller(lg_pos_bin = lg_labeller)) +
    xlab("position") +
    theme_cowplot()
}


selection_scan_outlier_nearby_genes_100kb <- find_nearby_genes(genome_selection_sites_outliers_PC2_sca2, gff_Genome_MOR, 1*10^5)
set.seed(1)
Outliers_Sca2_pc2 <- plot_nearby_genes(genome_selection_sites_outliers_PC2_sca2, selection_scan_outlier_nearby_genes_100kb)
ggsave(filename  = "EUostrea/Figures/Genome_scan/Outliers_Sca2_pc2.png", 
       plot=Outliers_Sca2_pc2, width = 40, height = 50, units = "cm", pointsize = 20, dpi = 250)
dev.off()

# Plot Sca1 PC2 top 3 outliers
selection_scan_outlier_nearby_genes_100kb <- find_nearby_genes(genome_selection_sites_outliers_PC2_sca1_trans, gff_Genome_MOR, 1*10^5)
set.seed(1)
PlotOutliers_Sca1_pc2_top3 <- plot_nearby_genes(genome_selection_sites_outliers_PC2_sca1_trans, selection_scan_outlier_nearby_genes_100kb)
ggsave(filename  = "EUostrea/Figures/Genome_scan/Outliers_Sca1_pc2_top3.png", 
       plot=PlotOutliers_Sca1_pc2_top3, width = 40, height = 50, units = "cm", pointsize = 20, dpi = 250)
dev.off()

# Plot Sca3 PC2 top 3 outliers
selection_scan_outlier_nearby_genes_100kb <- find_nearby_genes(genome_selection_sites_outliers_PC2_sca3_trans, gff_Genome_MOR, 1*10^5)
set.seed(1)
PlotOutliers_Sca3_pc2_top3 <- plot_nearby_genes(genome_selection_sites_outliers_PC2_sca3_trans, selection_scan_outlier_nearby_genes_100kb)
ggsave(filename  = "EUostrea/Figures/Genome_scan/tOutliers_Sca1_pc2_top1.png", 
       plot=PlotOutliers_Sca3_pc2_top3, width = 40, height = 50, units = "cm", pointsize = 20, dpi = 250)
dev.off()










# Group the data by "lg" and calculate the min and max positions
genome_selection_sites_outliers_PC9_mm <- genome_selection_sites_outliers_PC9 %>%
  group_by(lg) %>%
  summarize(min_pos = min(pos), max_pos = max(pos))

# Filter the result for "scaffold2" and "scaffold10"
result_filtered <- genome_selection_sites_outliers_PC9_mm %>%
  filter(lg %in% c("scaffold2", "scaffold10"))

# Plot resulting PCA of the Outliers SNPS on the PC9 of the PCAngsd --selection genome scan analysis
cov_mat <- as.matrix(read.table("~/Desktop/Scripts/Data/GenomeScan_EUostrea/Outliers/Outliers_PC9_sca2_10_10feb23.covMat")) 
annot <- read.table("../Scripts/EUostrea/01_infofiles/bamlist_EUostrea.annot", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
# Loads functions
source("~/Desktop/Scripts/Flat_oysters/04_local_R/00_scripts/individual_pca_functions_hjam_feb23.R")
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
Outliers_PC9_sca2_10_pcangsd<-PCA(cov_mat, annot$V1, annot$V2, 1, 2, show.ellipse = F, show.label = FALSE)
ggsave(filename  = "EUostrea/Figures/Genome_scan/Outliers/PCA_SNPOutliers_PC9_sca2_10_10feb23.png", 
       plot=Outliers_PC9_sca2_10_pcangsd, width = 60, height = 40, units = "cm", pointsize = 20, dpi = 250)
dev.off()

#
##
### The END ~~~~~


