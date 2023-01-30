### The BEGINNING ~~~~~
##
# ~ Plots Population structure results | First written by Nicolas R, Lou  with later modifications by  Homère J. Alves Monteiro


# Cleans the environment ~ 
rm(list=ls())

# Sets working directory ~
setwd("~/Desktop/Scripts/")

# Loads required packages ~
pacman::p_load(vegan, tidyverse, RcppCNPy, pheatmap, extrafont, ggforce, ggrepel, ggstar, RcppCNPy)


## Plotting CovMat and IbsMat from ANGSD
cov_mat <- as.matrix(read.table("~/Desktop/Scripts/Data/covMat")) 
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
PCA(cov_mat, annot$V1, annot$V2, 1, 2, show.ellipse = FALSE, show.label = FALSE)
PCA(cov_mat, annot$V1, annot$V2, 2, 3, show.ellipse = FALSE, show.label = FALSE)
DAPC(n = 50, x_axis=1, y_axis = 2)

#ggsave(, file = "~/Desktop/Scripts/EUostrea/Figures/PopulationStructure/.pdf", device = cairo_pdf, scale = 1.1, width = 12, height = 8, dpi = 300)
#dev.off()

#Plot genome-wide PCoA with the ibsMat matrix
ibs_mat <- read_tsv("~/.ibsMat", col_names = F) %>% 
  dplyr::select(1:nrow(.)) %>%
  as.matrix()
PCoA(ibs_mat, annot$V1, annot$V2, 153, 1, 2, show.ellipse = F)


## Plotting CovMat and IbsMat from PCangsd
#LD pruned SNP list
genome_cov <- read_delim(".cov", delim = " ", col_names = FALSE) %>%
  as.matrix()
PCA(genome_cov, annot$V1, annot$V2, 1, 2, show.ellipse = F)

pca_table_ld_pruned <- pca_table
pca_table_joined <- pca_table_ld_pruned %>%
  dplyr::select(1:10) %>%
  left_join(sample_table, by=c("individual"="sample_id", "population"="population"))
pca_table_joined %>%
  ggplot(aes(x=PC1, y=PC2)) +
  geom_point(data=dplyr::select(pca_table_joined, -population), size = 0.1, color="grey") +
  geom_point(size=1, mapping = aes(color=island_shore)) +
  facet_wrap(~population) +
  theme_cowplot() +
  theme(axis.text = element_blank())


## Plot selection scan result with PCAngsd (e not specified)
# .npy and .sites files from PCangsd
genome_selection <- npyLoad(".selection.npy")
genome_selection_sites <- read_tsv("_mindp368_maxdp928_minind167_minq20_downsampled_unlinked.txt", col_names = c("lg", "pos", "major", "minor")) %>%
  bind_cols(read_table(".sites", col_names = "keep")) %>%
  filter(keep==1) %>%
  mutate(chi_squared=genome_selection[,1]) %>%
  filter(substr(lg, 1, 2) == 'LG') %>%
  dplyr::select(lg, pos, chi_squared) %>%
  mutate(neg_log_p_value=-log(1-pchisq(chi_squared, df=1))) 
## Extract sites exeeding  p-value cut off of 0.05 after Bonferroni correction and export them
genome_selection_sites_outliers <- filter(genome_selection_sites, neg_log_p_value > -log(0.05/dim(genome_selection)[1]))
## Plot
genome_selection_plot <- genome_selection_sites %>%
  ggplot(aes(x=pos/10^6, y=neg_log_p_value)) +
  geom_point(size=0.2, alpha=0.5) +
  geom_smooth(color="blue", se =F) +
  geom_hline(yintercept = -log(0.05/dim(genome_selection)[1]), linetype = "dashed") + # This line marks a p-value cut off of 0.05 after Bonferroni correction
  theme_cowplot() +
  scale_x_continuous(breaks=seq(0, 40, 5)) +
  xlab("position (Mbp)") +
  ylab("-log(p)") +
  facet_grid(~lg, scales="free_x", space="free_x") +
  theme(panel.spacing = unit(0.1, "lines"),
        axis.title.x=element_text(),
        legend.position="none",
        text = element_text(size=10),
        axis.text = element_text(size=6))
ggsave(filename  = "", plot=genome_selection_plot, width = 40, height = 10, units = "cm", pointsize = 20)
