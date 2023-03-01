### The BEGINNING ~~~~~
##
# ~ Plots --HAPfreq_PCA Inversion | First written by Nicolas Lou with later modifications by Homère J. Alves Monteiro
rm(list = ls(all = TRUE))
# Loads functions
source("~/Desktop/Scripts/Flat_oysters/04_local_R/00_scripts/individual_pca_functions_hjam_dec22.R")
# Loads required packages ~
pacman::p_load(tidyverse, cowplot, knitr, scales, MetBrewer)

# Reorders Population ~
annot <- read.table("~/Desktop/Scripts/EUostrea/01_infofiles/bamlist_EUostrea.annot", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
annot$V2 <- factor(annot$V2, ordered = T,
                   levels = c("MOLU", "ZECE", "CRES","ORIS","CORS", "PONT",  "RIAE", "MORL",
                              "USAM", "TOLL", "COLN", "BARR","TRAL", "CLEW", "RYAN",
                              "GREV", "WADD", "NISS","LOGS","VENO", "HALS", "THIS",
                              "KALV", "HYPP","LANG", "BUNN", "DOLV", "HAUG", "HAFR",
                              "INNE","VAGS", "AGAB", "OSTR"))


# Inv Reg scaffold4 
cov_mat <- as.matrix(read.table("~/Desktop/Scripts/Data/InvReg_EUostrea/7feb23_scaffold4_InvReg_pcangsd.cov")) 
PCA(cov_mat, annot$V1, annot$V2, 1, 2, show.ellipse = F)
lg04_pca_table <- pca_table
pca_table %>%
  ggplot(aes(x=PC1, fill=population)) +
  geom_histogram(color="black") +
  theme_cowplot()
# Inv Reg scaffold5
cov_mat <- as.matrix(read.table("~/Desktop/Scripts/Data/InvReg_EUostrea/7feb23_scaffold5_InvReg_pcangsd.cov")) 
PCA(cov_mat, annot$V1, annot$V2, 1, 2, show.ellipse = F)
lg05_pca_table <- pca_table
pca_table %>%
  ggplot(aes(x=PC1, fill=population)) +
  geom_histogram(color="black") +
  theme_cowplot()

# Inv Reg scaffold5
cov_mat <- as.matrix(read.table("~/Desktop/Scripts/Data/InvReg_EUostrea/7feb23_scaffold8_InvReg_pcangsd.cov")) 
PCA(cov_mat, annot$V1, annot$V2, 1, 2, show.ellipse = F)
lg08_pca_table <- pca_table
pca_table %>%
  ggplot(aes(x=PC1, fill=population)) +
  geom_histogram(color="black") +
  theme_cowplot()


# Population haplotype frequencies histogram
inversion_pca_table <- mutate(lg04_pca_table, lg="Reg04", PC1=PC1, min_cutoff=-0.03, max_cutoff=0.03) %>%
  bind_rows(mutate(lg05_pca_table, lg="Reg05", PC1=PC1,  min_cutoff=-0.03, max_cutoff=0.03)) %>%
  bind_rows(mutate(lg08_pca_table, lg="Reg08", PC1=PC1,  min_cutoff=0, max_cutoff=0.06)) %>%
  mutate(genotype=ifelse(PC1<min_cutoff, "homo_1", "het"),
         genotype=ifelse(PC1>max_cutoff, "homo_2", genotype)) %>%
  mutate(genotype=factor(genotype, levels = c("homo_1", "het", "homo_2"))) %>%
  mutate(minor_allele_count = 3-as.numeric(genotype)) %>%
  dplyr::select(individual, population, lg, PC1, min_cutoff, max_cutoff, genotype,minor_allele_count)
Gen_freq_histo <- inversion_pca_table %>%
  ggplot(aes(x=PC1, fill=population)) +
  scale_fill_manual(values =c( "#A02353", "#A02353", "#A02353",
                               "#AD5B35",
                               "#ad7358",
                               "#CC480C",  "#CC480C",
                               "#969696",
                               "#000000",
                               "#D38C89", "#D38C89", "#D38C89",
                               "#C89AD1", "#C89AD1",
                               "#7210A0",
                               "#91BD96", "#91BD96",
                               "#02630C","#02630C","#02630C", "#02630C", "#02630C",
                               "#45D1F7", "#45D1F7",
                               "#588cad", "#588cad", "#588cad", "#588cad", "#588cad",
                               "#240377", "#240377", "#240377", "#240377" ))+
  geom_histogram(color="black") +
  geom_vline(aes(xintercept = min_cutoff)) +
  geom_vline(aes(xintercept = max_cutoff)) +
  facet_wrap(~lg, scales = "free") +
  theme_cowplot() +
  theme(panel.background = element_rect(fill = "#ffffff"),
        panel.border = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.position = "right",
        legend.title = element_text(color = "#000000", size = 13),
        legend.text = element_text(size = 11),
        axis.title.x = element_text(size = 18, face = "bold", margin = margin(t = 20, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size = 18, face = "bold", margin = margin(t = 0, r = 20, b = 0, l = 0)),
        axis.text = element_text(color = "#000000", size = 13),
        axis.ticks = element_line(color = "#000000", size = 0.3),
        axis.line = element_line(colour = "#000000", size = 0.3)) +
  guides(fill = guide_legend(title = "Populations", title.theme = element_text(size = 15, face = "bold"),
                             label.theme = element_text(size = 14)))
ggsave(Gen_freq_histo, file = "~/Desktop/Scripts/EUostrea/Figures/InvReg/PopulationHaplotypesFrequencies_InvReg040508_histo.pdf", device = cairo_pdf, scale = 1.1, width = 12, height = 8, dpi = 300)
dev.off()

# Genotype frequencies
inversion_pca_table %>%
  count(population, lg, genotype) %>%
  ggplot(aes(x=population, y=n)) +
  geom_col(aes(fill=genotype), color="black") +
  scale_fill_manual(values = c("#5a0c30","#973fc9","#91ccff")) +
  facet_wrap(~lg, scales = "free") +
  theme_cowplot()

# Save dataframe for the Map pie chart 
inversion_pca_table %>%
  group_by(population, lg) %>%
  summarise(prop_homo_1 = mean(genotype == "homo_1"),
            prop_homo_2 = mean(genotype == "homo_2"),
            prop_het = mean(genotype == "het"), .groups = 'drop') %>% 
  write.table(., "~/Desktop/Scripts/Data/InvReg_EUostrea/PropHomo1_het_homo2_InvReg_27feb23.tsv", sep="\t", row.names=FALSE)

#
Genfreq1 <- inversion_pca_table %>%
  count(population, population, lg, genotype) %>%
  ggplot(aes(x=population, y=n)) +
  geom_col(aes(fill=genotype), color="black", position="fill") +
  scale_fill_manual(values = c("#5a0c30","#973fc9","#91ccff")) +
  facet_wrap(~lg, scales = "free") +
  theme_cowplot()+
  coord_flip()+
  theme(panel.background = element_rect(fill = "#ffffff"),
        panel.grid.major.x = element_line(colour = "#ededed", linetype = "dashed", size = .00005),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(),
        axis.line = element_line(colour = "#000000", size = .3),
        axis.title.x = element_blank(),
        axis.text.x = element_text(colour = "#000000", size = 10),
        axis.text.y = element_text(colour = c( "#A02353", "#A02353", "#A02353", "#AD5B35", "#ad7358",
                                               "#CC480C",  "#CC480C",
                                               "#969696",
                                               "#000000",
                                               "#D38C89", "#D38C89", "#D38C89",
                                               "#C89AD1", "#C89AD1",
                                               "#7210A0",
                                               "#91BD96", "#91BD96",
                                               "#02630C","#02630C","#02630C", "#02630C", "#02630C",
                                               "#45D1F7", "#45D1F7",
                                               "#588cad", "#588cad", "#588cad", "#588cad", "#588cad",
                                               "#240377", "#240377", "#240377", "#240377" ), size = 12, face = "bold", angle = -45, vjust = 1, hjust = 1),
        axis.ticks.x = element_line(colour = "#000000", size = .3),
        axis.ticks.y = element_line(colour = "#000000", size = .3),
        strip.text = element_text(colour = "#000000", size = 12, face = "bold"),
        legend.position = "top",
        legend.key = element_rect(fill = NA),
        legend.background =element_blank()) +
  guides(fill = guide_legend(title = "Haplotypes for putative inversion regions", title.theme = element_text(size = 12, face = "bold"),
                             label = TRUE,
                             label.theme = element_text(size = 14),
                             override.aes = list(size = 5, alpha = .9)), colour = "none")
ggsave(Genfreq1, file = "~/Desktop/Scripts/EUostrea/Figures/InvReg/3GenotypesFrequencies_InvReg040508_histo.pdf", device = cairo_pdf, scale = 1.1, width = 12, height = 8, dpi = 300)
dev.off()

inversion_allele_frequency <- inversion_pca_table %>%
  count(population, lg, genotype) %>%
  pivot_wider(names_from = genotype, values_from = n, values_fill=0) %>%
  mutate(minor=2*homo_1+het, major=2*homo_2+het) %>%
  dplyr::select(-c(homo_1, het, homo_2)) %>%
  pivot_longer(cols = c(minor, major), names_to = "genotype", values_to = "n") %>%
  group_by(lg, population) %>%
  mutate(sum_n = sum(n), p = n/sum_n) %>%
  ungroup()
Genfreq2 <- inversion_allele_frequency %>% 
  ggplot(aes(x=population, y=n)) +
  geom_col(aes(fill=genotype), color="black", position="fill") +
  scale_fill_manual(values = c("#91ccff", "#5a0c30")) +
  facet_wrap(~lg, scales = "free") +
  theme_cowplot()+
  coord_flip()+
  theme(panel.background = element_rect(fill = "#ffffff"),
        panel.grid.major.x = element_line(colour = "#ededed", linetype = "dashed", size = .00005),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(),
        axis.line = element_line(colour = "#000000", size = .3),
        axis.title.x = element_blank(),
        axis.text.x = element_text(colour = "#000000", size = 10),
        axis.text.y = element_text(colour = c( "#A02353", "#A02353", "#A02353", "#AD5B35", "#ad7358",
                                                       "#CC480C",  "#CC480C",
                                                       "#969696",
                                                       "#000000",
                                                       "#D38C89", "#D38C89", "#D38C89",
                                                       "#C89AD1", "#C89AD1",
                                                       "#7210A0",
                                                       "#91BD96", "#91BD96",
                                                       "#02630C","#02630C","#02630C", "#02630C", "#02630C",
                                                       "#45D1F7", "#45D1F7",
                                                       "#588cad", "#588cad", "#588cad", "#588cad", "#588cad",
                                                       "#240377", "#240377", "#240377", "#240377" ), size = 12, face = "bold", angle = -45, vjust = 1, hjust = 1),
        axis.ticks.x = element_line(colour = "#000000", size = .3),
        axis.ticks.y = element_line(colour = "#000000", size = .3),
        strip.text = element_text(colour = "#000000", size = 12, face = "bold"),
        legend.position = "top",
        legend.key = element_rect(fill = NA),
        legend.background =element_blank()) +
  guides(fill = guide_legend(title = "Alleles for putative inversion regions", title.theme = element_text(size = 12, face = "bold"),
                             label = TRUE,
                             label.theme = element_text(size = 14),
                             override.aes = list(size = 5, alpha = .9)), colour = "none")
ggsave(Genfreq2, file = "~/Desktop/Scripts/EUostrea/Figures/InvReg/MajorMinorGenotypesFrequencies_InvReg040508_histo.pdf", device = cairo_pdf, scale = 1.1, width = 12, height = 8, dpi = 300)
dev.off()


# Heterozygosity at inverted regions
inversion_h <- inversion_pca_table %>%
  count(population, lg, genotype, minor_allele_count) %>%
  group_by(population, lg) %>%
  mutate(genotype_frequency=n/sum(n), 
         het=ifelse(minor_allele_count==1, 1, 0),
         p=minor_allele_count*genotype_frequency/2) %>%
  group_by(population, lg) %>%
  summarise(observed=sum(genotype_frequency*het), expected=2*sum(p)*(1-sum(p))) %>%
  ungroup() %>%
  pivot_longer(cols = c(observed, expected), names_to = "type", values_to = "heterozygosity")
inversion_heterozygosity_plot <- inversion_h %>%
  ggplot(aes(x = population, y = heterozygosity, color = type, linetype = type, alpha = 0.3)) +
  geom_point(size = 3) +
  scale_color_manual(values = MetBrewer::met.brewer("Archambault", n = 2, type = "discrete")) +
  facet_wrap(~ lg) +  
  coord_flip() +
  labs(title = "Heterozygosity values at putative inversion-regions by population",
       color = "Type",
       x = "Population",
       y = "Heterozygosity") +
  theme_cowplot() +
  theme(panel.background = element_rect(fill = "#ffffff"),
        panel.grid.major.x = element_line(colour = "#ededed", linetype = "dashed", size = .00005),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(),
        axis.line = element_line(colour = "#000000", size = .3),
        axis.text.x = element_text(colour = "#000000", size = 10),
        axis.text.y = element_text(colour = c("#A02353", "#A02353", "#A02353", "#AD5B35", "#ad7358",
                                              "#CC480C", "#CC480C",
                                              "#969696",
                                              "#000000",
                                              "#D38C89", "#D38C89", "#D38C89",
                                              "#C89AD1", "#C89AD1",
                                              "#7210A0",
                                              "#91BD96", "#91BD96",
                                              "#02630C", "#02630C", "#02630C", "#02630C", "#02630C",
                                              "#45D1F7", "#45D1F7",
                                              "#588cad", "#588cad", "#588cad", "#588cad", "#588cad",
                                              "#240377", "#240377", "#240377", "#240377"), size = 12, face = "bold", angle = -45, vjust = 1, hjust = 1),
        axis.ticks = element_line(colour = "#000000", size = .3),
        strip.text = element_text(colour = "#000000", size = 12, face = "bold"),
        legend.text = element_text(colour = "#000000", size = 12, face = "bold"),
        legend.position = "top",
        legend.background = element_blank())
ggsave(inversion_heterozygosity_plot, file = "~/Desktop/Scripts/EUostrea/Figures/InvReg/HetatInvReg_InvReg040508_histo.pdf", device = cairo_pdf, scale = 1.1, width = 12, height = 8, dpi = 300)
dev.off()






#
##
### The END ~~~~~

























# 
# 
# #### Scaffold4 InvReg ####
# #Loads data
# #cov_mat <- read_tsv("~/Desktop/Scripts/Data/PCA_inv/8sept22_subset/8sept22_subset_scaffold4_inv_pcangsd.cov", col_names = FALSE) %>% 
# #  dplyr::select(1:nrow(.)) %>%
# #  as.matrix()
# 
# cov_mat <- as.matrix(read.table("~/Desktop/Scripts/Data/InvReg_EUostrea/7feb23_scaffold4_InvReg_pcangsd.cov")) 
# annot <- read.table("../Scripts/EUostrea/01_infofiles/bamlist_EUostrea.annot", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
# 
# # Reorders Population ~
# annot$V2 <- factor(annot$V2, ordered = T,
#                    levels = c("MOLU", "ZECE", "CRES",
#                               "ORIS","CORS", "PONT",  "RIAE",
#                               "MORL",
#                               "USAM",
#                               "TOLL", "COLN", "BARR",
#                               "TRAL", "CLEW",
#                               "RYAN",
#                               "GREV", "WADD",
#                               "NISS","LOGS","VENO", "HALS", "THIS",
#                               "KALV", "HYPP",
#                               "LANG", "BUNN", "DOLV", "HAUG", "HAFR",
#                               "INNE","VAGS", "AGAB", "OSTR"))
# #Plot genome-wide PCA with the covMat matrix
# PCA(cov_mat, annot$V1, annot$V2, 1, 2, show.ellipse = FALSE, show.label = TRUE)
# PCA(cov_mat, annot$V1, annot$V2, 2, 3, show.ellipse = FALSE, show.label = FALSE)
# 
# str(unique(pca_table$population))
# 
# #Plot genome-wide PCoA with the ibsMat matrix
# ibs_mat <- read_tsv("~/Desktop/Scripts/Data/InvReg_EUostrea/6feb23_scaffold4_InvReg.ibsMat", col_names = F) %>% 
#   dplyr::select(1:nrow(.)) %>%
#   as.matrix()
# PCoA(ibs_mat, annot$V1, annot$V2, 153, 1, 2, show.ellipse = F)
# 
# # Reorders Population ~
# pca_table$population <- factor(pca_table$population, ordered = T,
#                                levels = c("MOLU", "ZECE", "CRES",
#                                           "ORIS","CORS", "PONT",  "RIAE",
#                                           "MORL",
#                                           "USAM",
#                                           "TOLL", "COLN", "BARR",
#                                           "TRAL", "CLEW",
#                                           "RYAN",
#                                           "GREV", "WADD",
#                                           "NISS","LOGS","VENO", "HALS", "THIS",
#                                           "KALV", "HYPP",
#                                           "LANG", "BUNN", "DOLV", "HAUG", "HAFR",
#                                           "INNE","VAGS", "AGAB", "OSTR")) 
# 
# #histogram 
# pca_table %>%
#   ggplot(aes(x=PC1, fill=population)) +
#   scale_fill_manual(values =c( "#A02353", "#A02353", "#A02353",
#                                  "#AD5B35",
#                                  "#ad7358",
#                                  "#CC480C",  "#CC480C",
#                                  "#969696",
#                                  "#000000",
#                                  "#D38C89", "#D38C89", "#D38C89",
#                                  "#C89AD1", "#C89AD1",
#                                  "#7210A0",
#                                  "#91BD96", "#91BD96",
#                                  "#02630C","#02630C","#02630C", "#02630C", "#02630C",
#                                  "#45D1F7", "#45D1F7",
#                                  "#588cad", "#588cad", "#588cad", "#588cad", "#588cad",
#                                  "#240377", "#240377", "#240377", "#240377" ))+
#   geom_histogram(color="black") +
#   theme_cowplot()
# 
# pca_table %>%
#   ggplot(aes(x=PC1, fill=population, group=individual)) +
#   scale_fill_manual(values =c( "#A02353", "#A02353", "#A02353",
#                                "#AD5B35",
#                                "#ad7358",
#                                "#CC480C",  "#CC480C",
#                                "#969696",
#                                "#000000",
#                                "#D38C89", "#D38C89", "#D38C89",
#                                "#C89AD1", "#C89AD1",
#                                "#7210A0",
#                                "#91BD96", "#91BD96",
#                                "#02630C","#02630C","#02630C", "#02630C", "#02630C",
#                                "#45D1F7", "#45D1F7",
#                                "#588cad", "#588cad", "#588cad", "#588cad", "#588cad",
#                                "#240377", "#240377", "#240377", "#240377" ))+
#   
#   geom_histogram(color="black") +
#   facet_wrap(~population) +
#   theme_cowplot()
# 
# Gen_freq_histo <- pca_table %>%
#   ggplot(aes(x=PC1, fill=population, group=individual)) +
#   scale_fill_manual(values =c("#A02353", "#A02353", "#A02353",
#                               "#AD5B35",
#                               "#ad7358",
#                               "#CC480C",  "#CC480C",
#                               "#969696",
#                               "#000000",
#                               "#D38C89", "#D38C89", "#D38C89",
#                               "#C89AD1", "#C89AD1",
#                               "#7210A0",
#                               "#91BD96", "#91BD96",
#                               "#02630C","#02630C","#02630C", "#02630C", "#02630C",
#                               "#45D1F7", "#45D1F7",
#                               "#588cad", "#588cad", "#588cad", "#588cad", "#588cad",
#                               "#240377", "#240377", "#240377", "#240377" ))+
#   geom_histogram(color="black") +
#   facet_wrap(~population) +
#   theme_cowplot()
# ggsave(Gen_freq_histo, file = "~/Desktop/Scripts/EUostrea/Figures/InvReg/scaffold4/Scaffold4inversion_genotypes_frequencies_histo.pdf", device = cairo_pdf, scale = 1.1, width = 12, height = 8, dpi = 300)
# dev.off()
# 
# #### Plot PCAngsd selection scan on InvReg #### 
# genome_selection <- read_tsv("~/Desktop/Scripts/Data/InvReg_EUostrea/1feb23_GlobalVariants_InvRegscaffold8_pcangsd.selection.tsv", col_names = F)
# genome_selection_sites <- read_table("Data/InvReg_EUostrea/scaffold8_InvReg_awk.txt", 
#                                    col_names = c("lg", "pos", "major", "minor")) %>%
#   bind_cols(read_table("Data/InvReg_EUostrea/7feb23_scaffold8_InvReg_pcangsd.sites", col_names = "keep")) %>%
#   filter(keep==1) %>%
#   mutate(PC1=genome_selection[,1],
#          PC2=genome_selection[,2],
#          PC3=genome_selection[,3],
#          PC4=genome_selection[,4],
#          PC5=genome_selection[,5],
#          PC6=genome_selection[,6],
#          PC7=genome_selection[,7],
#          PC8=genome_selection[,8],
#          PC9=genome_selection[,9],
#          PC10=genome_selection[,10]) %>%
#   pivot_longer(names_to = "pc", values_to = "chi_squared", cols = 6:15) %>%
#   dplyr::select(lg, pos, pc, chi_squared) %>%
#   mutate(chi_squared = ifelse(pc == "PC1", chi_squared$X1, 
#                               ifelse(pc == "PC2", chi_squared$X2,
#                                      ifelse(pc == "PC3", chi_squared$X3, 
#                                             ifelse(pc == "PC4", chi_squared$X4, 
#                                                    ifelse(pc == "PC5", chi_squared$X5,
#                                                           ifelse(pc == "PC6", chi_squared$X6,
#                                                                  ifelse(pc == "PC7", chi_squared$X7,
#                                                                         ifelse(pc == "PC8", chi_squared$X8, 
#                                                                                ifelse(pc == "PC9", chi_squared$X9,
#                                                                                       ifelse(pc == "PC10", chi_squared$X10,NA))))))))))) %>%
#   filter(!is.na(chi_squared)) %>%
#   mutate(neg_log_p_value=-log(1-pchisq(chi_squared, df=1))) 
# 
# #reorder chr
# genome_selection_sites$lg <- factor(genome_selection_sites$lg , ordered = T,
#                                     levels = c("scaffold1", "scaffold2", "scaffold3", "scaffold4", "scaffold5",
#                                                "scaffold6", "scaffold7", "scaffold8", "scaffold9", "scaffold10"))
# genome_selection_sites$pc <- factor(genome_selection_sites$pc , ordered = T,
#                                     levels = c("PC1", "PC2", "PC3", "PC4", "PC5",
#                                                "PC6", "PC7", "PC8", "PC9", "PC10"))
# 
# genome_selection_plot <- genome_selection_sites %>%
#   ggplot(aes(x=pos/10^6, y=neg_log_p_value)) +
#   geom_point(size=0.2, alpha=0.5) +
#   geom_smooth(color="purple", se =F) +
#   geom_hline(yintercept = -log(0.05/dim(genome_selection)[1]), linetype = "dashed") + # This line marks a p-value cut off of 0.05 after Bonferroni correction
#   theme_cowplot() +
#   scale_x_continuous(breaks=seq(0, 100, 10)) +
#   coord_cartesian(ylim=c(0, 40), expand = F) +
#   xlab("position (Mbp)") +
#   ylab("-log(p)") +
#   facet_grid(pc~lg, scales="free_x", space="free_x") +
#   theme(panel.spacing = unit(0.1, "lines"),
#         axis.title.x=element_text(),
#         legend.position="none",
#         text = element_text(size=10),
#         axis.text = element_text(size=6))+
#   theme(legend.key = element_blank()) +
#   theme(legend.title=element_blank()) +
#   theme(axis.title.x = element_text(size = 10, color="#000000", face="bold", margin = margin(t = 20, r = 0, b = 0, l = 0)),
#         axis.title.y = element_text(size = 10, color="#000000", face="bold", margin = margin(t = 0, r = 20, b = 0, l = 0))) +
#   theme(legend.text=element_text(size=11)) +
#   theme(panel.background = element_rect(fill = '#FAFAFA')) +
#   theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank()) +
#   theme(axis.line = element_line(colour = "#000000", size = 0.3)) +
#   theme(panel.border = element_blank())
# 
# ggsave(genome_selection_plot, file = "~/Desktop/Scripts/EUostrea/Figures/InvReg/scaffold8/9feb23_PCAngsdSelection_InvReg_scaffold8.png", device = cairo_pdf, scale = 1.1, width = 12, height = 8, dpi = 300)
# dev.off()
# 
# ## Extract sites exceeding  p-value cut off of 0.05 after Bonferroni correction and export them
# genome_selection_sites_outliers_InvRegSca8 <- filter(genome_selection_sites, neg_log_p_value > -log(0.05/dim(genome_selection)[1]))
# write_tsv(genome_selection_sites_outliers_InvRegSca4, file = "~/Desktop/Scripts/EUostrea/Figures/InvReg/scaffold8/9feb23_PCAngsdSelection_InvReg_scaffold8_SNPSpvalueaboveBonf.tsv")
# 
# #### Genotype frequencies ####
# 
# inversion_pca_table <- mutate(pca_table, lg="scaffold4", PC1=PC1, min_cutoff=-0.03, max_cutoff=0.03) %>% #change min_cutoff=-0.05, max_cutoff=0.01
#   mutate(genotype=ifelse(PC1<min_cutoff, "homo_1", "het"),
#          genotype=ifelse(PC1>max_cutoff, "homo_2", genotype)) %>%
#   mutate(genotype=factor(genotype, levels = c("homo_1", "het", "homo_2"))) %>%
#   mutate(minor_allele_count = 3-as.numeric(genotype)) %>%
#   dplyr::select(individual, population, lg, PC1, min_cutoff, max_cutoff, genotype,minor_allele_count)
# cut_off <- inversion_pca_table %>%
#   ggplot(aes(x=PC1, fill=population)) +
#   geom_histogram(color="black") +
#   geom_vline(aes(xintercept = min_cutoff)) +
#   geom_vline(aes(xintercept = max_cutoff)) +
#   facet_wrap(~lg, scales = "free") +
#   theme_cowplot()
# ggsave(cut_off, file = "~/Desktop/Scripts/EUostrea/Figures/InvReg/scaffold4/9feb23_CutoffPCA_InvRegSaffold4.pdf", device = cairo_pdf, scale = 1.1, width = 12, height = 8, dpi = 300)
# dev.off()
# inversion_pca_table %>%
#   count(population, population, lg, genotype) %>%
#   ggplot(aes(x=population, y=n)) +
#   geom_col(aes(fill=genotype), color="black") +
#   scale_fill_viridis_d(begin = 0.5, option = "A") +
#   facet_wrap(~lg, scales = "free") 
# Gen_freq_prop <- inversion_pca_table %>%
#   count(population, population, lg, genotype) %>%
#   ggplot(aes(x=population, y=n)) +
#   geom_col(aes(fill=genotype), color="black", position="fill") +
#   scale_fill_viridis_d(begin = 0.5, option = "A") +
#   facet_wrap(~lg, scales = "free")
# ggsave(Gen_freq_prop, file = "~/Desktop/Scripts/Flat_oysters/04_local_R/03_results/PCA_inv/8sept22_subset/Scaffold8inversion_proportion_genotypefrequencies.pdf", device = cairo_pdf, scale = 1.1, width = 12, height = 8, dpi = 300)
# dev.off()
# inversion_allele_frequency <- inversion_pca_table %>%
#   count(population, lg, genotype) %>%
#   pivot_wider(names_from = genotype, values_from = n, values_fill=0) %>%
#   mutate(minor=2*homo_1+het, major=2*homo_2+het) %>%
#   dplyr::select(-c(homo_1, het, homo_2)) %>%
#   pivot_longer(cols = c(minor, major), names_to = "genotype", values_to = "n") %>%
#   group_by(lg, population) %>%
#   mutate(sum_n = sum(n), p = n/sum_n) %>%
#   ungroup() 
# major_minor <- inversion_allele_frequency %>% 
#   ggplot(aes(x=population, y=n)) +
#   geom_col(aes(fill=genotype), color="black", position="fill") +
#   scale_fill_manual(values = MetBrewer::met.brewer("Nattier", n = 2, type = "discrete")) +
#   facet_wrap(~lg, scales = "free") 
# ggsave(major_minor, file = "~/Desktop/Scripts/Flat_oysters/04_local_R/03_results/PCA_inv/8sept22_subset/Scaffold8inversion_proportion_majorminorgenotype.pdf", device = cairo_pdf, scale = 1.1, width = 12, height = 8, dpi = 300)
# dev.off()
# inversion_allele_frequency_plot <- inversion_allele_frequency %>%
#   filter(genotype=="minor") %>%
#   ggplot(aes(x=population, y=p, color=population, shape=lg)) +
#   geom_point(size=3) + 
#   scale_color_manual(values =c( "#000000",
#                                "#A02353",
#                                "#AD5B35",
#                                "#CC480C",
#                                "#969696",
#                                "#D38C89",
#                                "#C89AD1","#C89AD1",
#                                "#7210A0" , "#7210A0",
#                                "#91BD96","#91BD96",
#                                "#02630C", "#02630C", 
#                                "#45D1F7",
#                                "#588CAD","#588CAD",
#                                "#240377","#240377","#240377","#240377" ))+  
#   facet_wrap(~"Allele frequencies at scaffold8 inversion region \n pos:40000018-58502426") +
#   theme_cowplot() +
#   theme(legend.position = "none")+
#   theme(title = element_text(size = 10, color="#000000", face="bold"),
#         axis.title.x = element_text(size = 10, color="#000000", face="bold", margin = margin(t = 20, r = 0, b = 0, l = 0)),
#         axis.title.y = element_text(size = 10, color="#000000", face="bold", margin = margin(t = 0, r = 20, b = 0, l = 0))) +
#   theme(panel.background = element_rect(fill = '#FAFAFA')) +
#   theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank()) +
#   theme(axis.line = element_line(colour = "#000000", size = 0.3)) +
#   theme(panel.border = element_blank())
# ggsave(inversion_allele_frequency_plot, file = "~/Desktop/Scripts/Flat_oysters/04_local_R/03_results/PCA_inv/8sept22_subset/Scaffold8inversion_allele_frequency_plot.pdf", device = cairo_pdf, scale = 1.1, width = 12, height = 8, dpi = 300)
# dev.off()
# ### Heterozygosity at inverted regions ###
# inversion_h <- inversion_pca_table %>%
#   count(population, lg, genotype, minor_allele_count) %>%
#   group_by(population, lg) %>%
#   mutate(genotype_frequency=n/sum(n), 
#          het=ifelse(minor_allele_count==1, 1, 0),
#          p=minor_allele_count*genotype_frequency/2) %>%
#   group_by(population, lg) %>%
#   summarise(observed=sum(genotype_frequency*het), expected=2*sum(p)*(1-sum(p))) %>%
#   ungroup() %>%
#   pivot_longer(cols = c(observed, expected), names_to = "type", values_to = "heterozygosity")
# inversion_h$population <- factor(inversion_h$population, ordered = T,
#                                levels = c("USAM",
#                                           "ZECE", 
#                                           "ORIS", 
#                                           "PONT",
#                                           "MORL",
#                                           "COLN",
#                                           "TRAL","CLEW",
#                                           "RYAN", "NELL",
#                                           "WADD","GREV",
#                                           "NISS", "HALS",
#                                           "HYPP",
#                                           "HAFR","LANG",
#                                           "AGAB", "INNE", "OSTR", "VAGS")) 
# 
# 
# 
# inversion_heterozygosity_plot <- inversion_h %>%
#   ggplot(aes(x=population, y=heterozygosity, color=type, linetype=type, alpha=0.3)) +
#   geom_point(size=3) +
#   scale_color_manual(values = MetBrewer::met.brewer("Isfahan1", n = 2, type = "discrete")) +
#   facet_wrap(~"Heterozygosity at inversions scaffold8 region") +
#   theme(panel.background = element_rect(fill = "#ffffff"),
#         panel.grid.major.x = element_line(colour = "#ededed", linetype = "dashed", size = .00005),
#         panel.grid.major.y = element_blank(),
#         panel.grid.minor = element_blank(), 
#         panel.border = element_blank(),
#         axis.line = element_line(colour = "#000000", size = .3),
#         axis.title.x = element_blank(),
#         axis.title.y = element_text(colour = "#000000", size = 16, face = "bold"),
#         axis.text.x = element_text(colour = "#000000", size = 12, face = "bold", angle = 45, vjust = 1, hjust = 1),
#         axis.text.y = element_text(colour = "#000000", size = 12),
#         axis.ticks.x = element_line(colour = "#000000", size = .3),
#         axis.ticks.y = element_line(colour = "#000000", size = .3),
#         strip.background.y = element_rect(colour = "#000000", fill = "#d6d6d6", size = 0.3),
#         strip.text = element_text(colour = "#000000", size = 12, face = "bold"),
#         legend.position = "top",
#         legend.margin = margin(t = 0, b = 0, r = 0, l = 0),
#         legend.box.margin = margin(t = 10, b = 20, r = 0, l = 0),
#         legend.key = element_rect(fill = NA),
#         legend.background =element_blank())
# inversion_heterozygosity_plot
# write_tsv(inversion_pca_table, "~/Desktop/Scripts/Data/PCA_inv/8sept22_subset/scaffold8inversion_genotype.tsv")
# ggsave(inversion_heterozygosity_plot, file = "~/Desktop/Scripts/Flat_oysters/04_local_R/03_results/PCA_inv/8sept22_subset/Scaffold8inversion_inversion_heterozygosity_plot.pdf", device = cairo_pdf, scale = 1.1, width = 12, height = 8, dpi = 300)
# dev.off()
# 
# 
# 
# #### Scaffold5 Inv without Lurida - 8septsubset ####
# #Loads data
# rm(list = ls(all = TRUE))
# # Loads functions
# source("~/Desktop/Scripts/Flat_oysters/04_local_R/00_scripts/individual_pca_functions.R")
# # Loads required packages ~
# cov_mat <- as.matrix(read.table("~/Desktop/Scripts/Data/PCA_inv/8sept22_subset/8sept22_subset_scaffold5_inv_pcangsd.cov")) 
# annot <- read.table("~/Desktop/Scripts/Data/PCA_inv/8sept22_subset/Bamlist_forInv_no_LURI_8sept22.annot", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
# #Plot genome-wide PCA with the covMat matrix
# ## Reorders Population ~
# PCA(cov_mat, annot$V1, annot$V2, 1, 2, show.ellipse = F)
# PCA(cov_mat, annot$V1, annot$V2, 3, 4, show.ellipse = F)
# #Plot genome-wide PCoA with the ibsMat matrix
# ibs_mat <- read_tsv("~/Desktop/Scripts/Data/PCA_inv/8sept22_subset/8sept22_subset_scaffold5_inv.ibsMat", col_names = F) %>% 
#   dplyr::select(1:nrow(.)) %>%
#   as.matrix()
# PCoA(ibs_mat, annot$V1, annot$V2, 153, 1, 2, show.ellipse = F)
# 
# # Reorders Population ~
# pca_table$population <- factor(pca_table$population, ordered = T,
#                                levels = c("USAM",
#                                           "ZECE", 
#                                           "ORIS", 
#                                           "PONT",
#                                           "MORL",
#                                           "COLN",
#                                           "TRAL","CLEW",
#                                           "RYAN", "NELL",
#                                           "WADD","GREV",
#                                           "NISS", "HALS",
#                                           "HYPP",
#                                           "HAFR","LANG",
#                                           "AGAB", "INNE", "OSTR", "VAGS")) 
# 
# #histogram 
# pca_table %>%
#   ggplot(aes(x=PC1, fill=population)) +
#   scale_fill_manual(values =c(   "#000000",
#                                  "#A02353",
#                                  "#AD5B35",
#                                  "#CC480C",
#                                  "#969696",
#                                  "#D38C89",
#                                  "#C89AD1","#C89AD1",
#                                  "#7210A0" , "#7210A0",
#                                  "#91BD96","#91BD96",
#                                  "#02630C", "#02630C", 
#                                  "#45D1F7",
#                                  "#588CAD","#588CAD",
#                                  "#240377","#240377","#240377","#240377"))+
#   geom_histogram(color="black") +
#   theme_cowplot()
# Gen_freq_histo <- pca_table %>%
#   ggplot(aes(x=PC1, fill=population, group=individual)) +
#   scale_fill_manual(values =c( "#000000",
#                                "#A02353",
#                                "#AD5B35",
#                                "#CC480C",
#                                "#969696",
#                                "#D38C89",
#                                "#C89AD1","#C89AD1",
#                                "#7210A0" , "#7210A0",
#                                "#91BD96","#91BD96",
#                                "#02630C", "#02630C", 
#                                "#45D1F7",
#                                "#588CAD","#588CAD",
#                                "#240377","#240377","#240377","#240377" ))+
#   geom_histogram(color="black") +
#   facet_wrap(~population) +
#   theme_cowplot()
# ggsave(Gen_freq_histo, file = "~/Desktop/Scripts/Flat_oysters/04_local_R/03_results/PCA_inv/8sept22_subset/Scaffold5inversion_genotypes_frequencies_histo.pdf", device = cairo_pdf, scale = 1.1, width = 12, height = 8, dpi = 300)
# dev.off()
# 
# genome_selection_sites <- read_tsv("~/Desktop/Scripts/Data/PCA_inv/8sept22_subset/8sept22_subset_scaffold5_inv.mafs.gz") %>%
#   dplyr::select(1:4)
# min(genome_selection_sites$position)
# 
# ### Genotype frequencies ###
# inversion_pca_table <- mutate(pca_table, lg="scaffold5 inversion region pos:199-24998832", PC1=PC1, min_cutoff=-0.03, max_cutoff=0.05) %>%
#   mutate(genotype=ifelse(PC1<min_cutoff, "homo_1", "het"),
#          genotype=ifelse(PC1>max_cutoff, "homo_2", genotype)) %>%
#   mutate(genotype=factor(genotype, levels = c("homo_1", "het", "homo_2"))) %>%
#   mutate(minor_allele_count = 3-as.numeric(genotype)) %>%
#   dplyr::select(individual, population, lg, PC1, min_cutoff, max_cutoff, genotype,minor_allele_count)
# cut_off <- inversion_pca_table %>%
#   ggplot(aes(x=PC1, fill=population)) +
#   geom_histogram(color="black") +
#   geom_vline(aes(xintercept = min_cutoff)) +
#   geom_vline(aes(xintercept = max_cutoff)) +
#   facet_wrap(~lg, scales = "free") +
#   theme_cowplot()
# cut_off
# ggsave(cut_off, file = "~/Desktop/Scripts/Flat_oysters/04_local_R/03_results/PCA_inv/8sept22_subset/Scaffold5inversion_cutoff_genotypefrequencies.pdf", device = cairo_pdf, scale = 1.1, width = 12, height = 8, dpi = 300)
# dev.off()
# inversion_pca_table %>%
#   count(population, population, lg, genotype) %>%
#   ggplot(aes(x=population, y=n)) +
#   geom_col(aes(fill=genotype), color="black") +
#   scale_fill_viridis_d(begin = 0.5, option = "A") +
#   facet_wrap(~lg, scales = "free") 
# Gen_freq_prop <- inversion_pca_table %>%
#   count(population, population, lg, genotype) %>%
#   ggplot(aes(x=population, y=n)) +
#   geom_col(aes(fill=genotype), color="black", position="fill") +
#   scale_fill_viridis_d(begin = 0.5, option = "A") +
#   facet_wrap(~lg, scales = "free")
# #ggsave(Gen_freq_prop, file = "~/Desktop/Scripts/Flat_oysters/04_local_R/03_results/PCA_inv/8sept22_subset/Scaffold5inversion_proportion_genotypefrequencies.pdf", device = cairo_pdf, scale = 1.1, width = 12, height = 8, dpi = 300)
# #dev.off()
# inversion_allele_frequency <- inversion_pca_table %>%
#   count(population, lg, genotype) %>%
#   pivot_wider(names_from = genotype, values_from = n, values_fill=0) %>%
#   mutate(minor=2*homo_1+het, major=2*homo_2+het) %>%
#   dplyr::select(-c(homo_1, het, homo_2)) %>%
#   pivot_longer(cols = c(minor, major), names_to = "genotype", values_to = "n") %>%
#   group_by(lg, population) %>%
#   mutate(sum_n = sum(n), p = n/sum_n) %>%
#   ungroup() 
# major_minor <- inversion_allele_frequency %>% 
#   ggplot(aes(x=population, y=n)) +
#   geom_col(aes(fill=genotype), color="black", position="fill") +
#   scale_fill_manual(values = MetBrewer::met.brewer("Nattier", n = 2, type = "discrete")) +
#   facet_wrap(~lg, scales = "free") 
# ggsave(major_minor, file = "~/Desktop/Scripts/Flat_oysters/04_local_R/03_results/PCA_inv/8sept22_subset/Scaffold5inversion_proportion_majorminorgenotype.pdf", device = cairo_pdf, scale = 1.1, width = 12, height = 8, dpi = 300)
# dev.off()
# 
# inversion_allele_frequency_plot <- inversion_allele_frequency %>%
#   filter(genotype=="minor") %>%
#   ggplot(aes(x=population, y=p, color=population, shape=lg)) +
#   geom_point(size=3) + 
#   scale_color_manual(values =c( "#000000",
#                                 "#A02353",
#                                 "#AD5B35",
#                                 "#CC480C",
#                                 "#969696",
#                                 "#D38C89",
#                                 "#C89AD1","#C89AD1",
#                                 "#7210A0" , "#7210A0",
#                                 "#91BD96","#91BD96",
#                                 "#02630C", "#02630C", 
#                                 "#45D1F7",
#                                 "#588CAD","#588CAD",
#                                 "#240377","#240377","#240377","#240377" ))+  
#   facet_wrap(~"Allele frequencies at scaffold5 inversion region \n pos:199-24998832") +
#   theme_cowplot() +
#   theme(legend.position = "none")+
#   theme(title = element_text(size = 10, color="#000000", face="bold"),
#         axis.title.x = element_text(size = 10, color="#000000", face="bold", margin = margin(t = 20, r = 0, b = 0, l = 0)),
#         axis.title.y = element_text(size = 10, color="#000000", face="bold", margin = margin(t = 0, r = 20, b = 0, l = 0))) +
#   theme(panel.background = element_rect(fill = '#FAFAFA')) +
#   theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank()) +
#   theme(axis.line = element_line(colour = "#000000", size = 0.3)) +
#   theme(panel.border = element_blank())
# ggsave(inversion_allele_frequency_plot, file = "~/Desktop/Scripts/Flat_oysters/04_local_R/03_results/PCA_inv/8sept22_subset/Scaffold5inversion_allele_frequency_plot.pdf", device = cairo_pdf, scale = 1.1, width = 12, height = 8, dpi = 300)
# dev.off()
# ### Heterozygosity at inverted regions ###
# inversion_h <- inversion_pca_table %>%
#   count(population, lg, genotype, minor_allele_count) %>%
#   group_by(population, lg) %>%
#   mutate(genotype_frequency=n/sum(n), 
#          het=ifelse(minor_allele_count==1, 1, 0),
#          p=minor_allele_count*genotype_frequency/2) %>%
#   group_by(population, lg) %>%
#   summarise(observed=sum(genotype_frequency*het), expected=2*sum(p)*(1-sum(p))) %>%
#   ungroup() %>%
#   pivot_longer(cols = c(observed, expected), names_to = "type", values_to = "heterozygosity")
# inversion_h$population <- factor(inversion_h$population, ordered = T,
#                                  levels = c("USAM",
#                                             "ZECE", 
#                                             "ORIS", 
#                                             "PONT",
#                                             "MORL",
#                                             "COLN",
#                                             "TRAL","CLEW",
#                                             "RYAN", "NELL",
#                                             "WADD","GREV",
#                                             "NISS", "HALS",
#                                             "HYPP",
#                                             "HAFR","LANG",
#                                             "AGAB", "INNE", "OSTR", "VAGS")) 
# 
# 
# 
# inversion_heterozygosity_plot <- inversion_h %>%
#   ggplot(aes(x=population, y=heterozygosity, color=type, linetype=type, alpha=0.3)) +
#   geom_point(size=3) +
#   scale_color_manual(values = MetBrewer::met.brewer("Isfahan2", n = 2, type = "discrete")) +
#   facet_wrap(~"Heterozygosity at inversions scaffold5 region") +
#   theme(panel.background = element_rect(fill = "#ffffff"),
#         panel.grid.major.x = element_line(colour = "#ededed", linetype = "dashed", size = .00005),
#         panel.grid.major.y = element_blank(),
#         panel.grid.minor = element_blank(), 
#         panel.border = element_blank(),
#         axis.line = element_line(colour = "#000000", size = .3),
#         axis.title.x = element_blank(),
#         axis.title.y = element_text(colour = "#000000", size = 16, face = "bold"),
#         axis.text.x = element_text(colour = "#000000", size = 12, face = "bold", angle = 45, vjust = 1, hjust = 1),
#         axis.text.y = element_text(colour = "#000000", size = 12),
#         axis.ticks.x = element_line(colour = "#000000", size = .3),
#         axis.ticks.y = element_line(colour = "#000000", size = .3),
#         strip.background.y = element_rect(colour = "#000000", fill = "#d6d6d6", size = 0.3),
#         strip.text = element_text(colour = "#000000", size = 12, face = "bold"),
#         legend.position = "top",
#         legend.margin = margin(t = 0, b = 0, r = 0, l = 0),
#         legend.box.margin = margin(t = 10, b = 20, r = 0, l = 0),
#         legend.key = element_rect(fill = NA),
#         legend.background =element_blank())
# inversion_heterozygosity_plot
# write_tsv(inversion_pca_table, "~/Desktop/Scripts/Data/PCA_inv/8sept22_subset/scaffold5inversion_genotype.tsv")
# ggsave(inversion_heterozygosity_plot, file = "~/Desktop/Scripts/Flat_oysters/04_local_R/03_results/PCA_inv/8sept22_subset/Scaffold5inversion_inversion_heterozygosity_plot.pdf", device = cairo_pdf, scale = 1.1, width = 12, height = 8, dpi = 300)
# dev.off()
# 
# 
# 
# #### Scaffold4 Inv without Lurida - 8septsubset ####
# #Loads data
# rm(list = ls(all = TRUE))
# # Loads functions
# source("~/Desktop/Scripts/Flat_oysters/04_local_R/00_scripts/individual_pca_functions.R")
# # Loads required packages ~
# cov_mat <- as.matrix(read.table("~/Desktop/Scripts/Data/PCA_inv/8sept22_subset/8sept22_subset_scaffold4_inv_pcangsd.cov")) 
# annot <- read.table("~/Desktop/Scripts/Data/PCA_inv/8sept22_subset/Bamlist_forInv_no_LURI_8sept22.annot", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
# #Plot genome-wide PCA with the covMat matrix
# ## Reorders Population ~
# PCA(cov_mat, annot$V1, annot$V2, 1, 2, show.ellipse = F)
# PCA(cov_mat, annot$V1, annot$V2, 3, 4, show.ellipse = F)
# #Plot genome-wide PCoA with the ibsMat matrix
# ibs_mat <- read_tsv("~/Desktop/Scripts/Data/PCA_inv/8sept22_subset/8sept22_subset_scaffold4_inv.ibsMat", col_names = F) %>% 
#   dplyr::select(1:nrow(.)) %>%
#   as.matrix()
# PCoA(ibs_mat, annot$V1, annot$V2, 153, 1, 2, show.ellipse = F)
# 
# # Reorders Population ~
# pca_table$population <- factor(pca_table$population, ordered = T,
#                                levels = c("USAM",
#                                           "ZECE", 
#                                           "ORIS", 
#                                           "PONT",
#                                           "MORL",
#                                           "COLN",
#                                           "TRAL","CLEW",
#                                           "RYAN", "NELL",
#                                           "WADD","GREV",
#                                           "NISS", "HALS",
#                                           "HYPP",
#                                           "HAFR","LANG",
#                                           "AGAB", "INNE", "OSTR", "VAGS")) 
# 
# #histogram 
# pca_table %>%
#   ggplot(aes(x=PC1, fill=population)) +
#   scale_fill_manual(values =c(   "#000000",
#                                  "#A02353",
#                                  "#AD5B35",
#                                  "#CC480C",
#                                  "#969696",
#                                  "#D38C89",
#                                  "#C89AD1","#C89AD1",
#                                  "#7210A0" , "#7210A0",
#                                  "#91BD96","#91BD96",
#                                  "#02630C", "#02630C", 
#                                  "#45D1F7",
#                                  "#588CAD","#588CAD",
#                                  "#240377","#240377","#240377","#240377"))+
#   geom_histogram(color="black") +
#   theme_cowplot()
# Gen_freq_histo <- pca_table %>%
#   ggplot(aes(x=PC1, fill=population, group=individual)) +
#   scale_fill_manual(values =c( "#000000",
#                                "#A02353",
#                                "#AD5B35",
#                                "#CC480C",
#                                "#969696",
#                                "#D38C89",
#                                "#C89AD1","#C89AD1",
#                                "#7210A0" , "#7210A0",
#                                "#91BD96","#91BD96",
#                                "#02630C", "#02630C", 
#                                "#45D1F7",
#                                "#588CAD","#588CAD",
#                                "#240377","#240377","#240377","#240377" ))+
#   geom_histogram(color="black") +
#   facet_wrap(~population) +
#   theme_cowplot()
# ggsave(Gen_freq_histo, file = "~/Desktop/Scripts/Flat_oysters/04_local_R/03_results/PCA_inv/8sept22_subset/Scaffold4inversion_genotypes_frequencies_histo.pdf", device = cairo_pdf, scale = 1.1, width = 12, height = 8, dpi = 300)
# dev.off()
# 
# genome_selection_sites <- read_tsv("~/Desktop/Scripts/Data/PCA_inv/8sept22_subset/8sept22_subset_scaffold4_inv.mafs.gz") %>%
#   dplyr::select(1:4)
# length(genome_selection_sites$position)
# 
# ### Genotype frequencies ###
# inversion_pca_table <- mutate(pca_table, lg="scaffold4 inversion region pos:1736-12999529", PC1=PC1, min_cutoff=-0.03, max_cutoff=0.05) %>%
#   mutate(genotype=ifelse(PC1<min_cutoff, "homo_1", "het"),
#          genotype=ifelse(PC1>max_cutoff, "homo_2", genotype)) %>%
#   mutate(genotype=factor(genotype, levels = c("homo_1", "het", "homo_2"))) %>%
#   mutate(minor_allele_count = 3-as.numeric(genotype)) %>%
#   dplyr::select(individual, population, lg, PC1, min_cutoff, max_cutoff, genotype,minor_allele_count)
# cut_off <- inversion_pca_table %>%
#   ggplot(aes(x=PC1, fill=population)) +
#   geom_histogram(color="black") +
#   geom_vline(aes(xintercept = min_cutoff)) +
#   geom_vline(aes(xintercept = max_cutoff)) +
#   facet_wrap(~lg, scales = "free") +
#   theme_cowplot()
# cut_off
# ggsave(cut_off, file = "~/Desktop/Scripts/Flat_oysters/04_local_R/03_results/PCA_inv/8sept22_subset/Scaffold4inversion_cutoff_genotypefrequencies.pdf", device = cairo_pdf, scale = 1.1, width = 12, height = 8, dpi = 300)
# dev.off()
# inversion_pca_table %>%
#   count(population, population, lg, genotype) %>%
#   ggplot(aes(x=population, y=n)) + 
#   geom_col(aes(fill=genotype), color="black") +
#   scale_fill_viridis_d(begin = 0.5, option = "B") +
#   facet_wrap(~lg, scales = "free") 
# Gen_freq_prop <- inversion_pca_table %>%
#   count(population, population, lg, genotype) %>%
#   ggplot(aes(x=population, y=n)) +
#   geom_col(aes(fill=genotype), color="black", position="fill") +
#   scale_fill_viridis_d(begin = 0.5, option = "D") +
#   facet_wrap(~lg, scales = "free")
# ggsave(Gen_freq_prop, file = "~/Desktop/Scripts/Flat_oysters/04_local_R/03_results/PCA_inv/8sept22_subset/Scaffold4inversion_proportion_genotypefrequencies.pdf", device = cairo_pdf, scale = 1.1, width = 12, height = 8, dpi = 300)
# dev.off()
# inversion_allele_frequency <- inversion_pca_table %>%
#   count(population, lg, genotype) %>%
#   pivot_wider(names_from = genotype, values_from = n, values_fill=0) %>%
#   mutate(minor=2*homo_1+het, major=2*homo_2+het) %>%
#   dplyr::select(-c(homo_1, het, homo_2)) %>%
#   pivot_longer(cols = c(minor, major), names_to = "genotype", values_to = "n") %>%
#   group_by(lg, population) %>%
#   mutate(sum_n = sum(n), p = n/sum_n) %>%
#   ungroup() 
# major_minor <- inversion_allele_frequency %>% 
#   ggplot(aes(x=population, y=n)) +
#   geom_col(aes(fill=genotype), color="black", position="fill") +
#   scale_fill_manual(values = MetBrewer::met.brewer("Nattier", n = 2, type = "discrete")) +
#   facet_wrap(~lg, scales = "free") 
# ggsave(major_minor, file = "~/Desktop/Scripts/Flat_oysters/04_local_R/03_results/PCA_inv/8sept22_subset/Scaffold4inversion_proportion_majorminorgenotype.pdf", device = cairo_pdf, scale = 1.1, width = 12, height = 8, dpi = 300)
# dev.off()
# 
# inversion_allele_frequency_plot <- inversion_allele_frequency %>%
#   filter(genotype=="minor") %>%
#   ggplot(aes(x=population, y=p, color=population, shape=lg)) +
#   geom_point(size=3) + 
#   scale_color_manual(values =c( "#000000",
#                                 "#A02353",
#                                 "#AD5B35",
#                                 "#CC480C",
#                                 "#969696",
#                                 "#D38C89",
#                                 "#C89AD1","#C89AD1",
#                                 "#7210A0" , "#7210A0",
#                                 "#91BD96","#91BD96",
#                                 "#02630C", "#02630C", 
#                                 "#45D1F7",
#                                 "#588CAD","#588CAD",
#                                 "#240377","#240377","#240377","#240377" ))+  
#   facet_wrap(~"Allele frequencies at scaffold4 inversion region pos:1736-12999529") +
#   theme_cowplot() +
#   theme(legend.position = "none")+
#   theme(title = element_text(size = 10, color="#000000", face="bold"),
#         axis.title.x = element_text(size = 10, color="#000000", face="bold", margin = margin(t = 20, r = 0, b = 0, l = 0)),
#         axis.title.y = element_text(size = 10, color="#000000", face="bold", margin = margin(t = 0, r = 20, b = 0, l = 0))) +
#   theme(panel.background = element_rect(fill = '#FAFAFA')) +
#   theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank()) +
#   theme(axis.line = element_line(colour = "#000000", size = 0.3)) +
#   theme(panel.border = element_blank())
# ggsave(inversion_allele_frequency_plot, file = "~/Desktop/Scripts/Flat_oysters/04_local_R/03_results/PCA_inv/8sept22_subset/Scaffold4inversion_allele_frequency_plot.pdf", device = cairo_pdf, scale = 1.1, width = 12, height = 8, dpi = 300)
# dev.off()
#  ### Heterozygosity at inverted regions ###
# inversion_h <- inversion_pca_table %>%
#   count(population, lg, genotype, minor_allele_count) %>%
#   group_by(population, lg) %>%
#   mutate(genotype_frequency=n/sum(n), 
#          het=ifelse(minor_allele_count==1, 1, 0),
#          p=minor_allele_count*genotype_frequency/2) %>%
#   group_by(population, lg) %>%
#   summarise(observed=sum(genotype_frequency*het), expected=2*sum(p)*(1-sum(p))) %>%
#   ungroup() %>%
#   pivot_longer(cols = c(observed, expected), names_to = "type", values_to = "heterozygosity")
# inversion_h$population <- factor(inversion_h$population, ordered = T,
#                                  levels = c("USAM",
#                                             "ZECE", 
#                                             "ORIS", 
#                                             "PONT",
#                                             "MORL",
#                                             "COLN",
#                                             "TRAL","CLEW",
#                                             "RYAN", "NELL",
#                                             "WADD","GREV",
#                                             "NISS", "HALS",
#                                             "HYPP",
#                                             "HAFR","LANG",
#                                             "AGAB", "INNE", "OSTR", "VAGS")) 
# 
# 
# 
# inversion_heterozygosity_plot <- inversion_h %>%
#   ggplot(aes(x=population, y=heterozygosity, color=type, linetype=type, alpha=0.3, shape=lg)) +
#   geom_point(size=3) +
#   scale_color_manual(values = MetBrewer::met.brewer("Hokusai2", n = 2, type = "discrete")) +
#   facet_wrap(~"Heterozygosity at inversions scaffold4 region") +
#   theme(panel.background = element_rect(fill = "#ffffff"),
#         panel.grid.major.x = element_line(colour = "#ededed", linetype = "dashed", size = .00005),
#         panel.grid.major.y = element_blank(),
#         panel.grid.minor = element_blank(), 
#         panel.border = element_blank(),
#         axis.line = element_line(colour = "#000000", size = .3),
#         axis.title.x = element_blank(),
#         axis.title.y = element_text(colour = "#000000", size = 16, face = "bold"),
#         axis.text.x = element_text(colour = "#000000", size = 12, face = "bold", angle = 45, vjust = 1, hjust = 1),
#         axis.text.y = element_text(colour = "#000000", size = 12),
#         axis.ticks.x = element_line(colour = "#000000", size = .3),
#         axis.ticks.y = element_line(colour = "#000000", size = .3),
#         strip.background.y = element_rect(colour = "#000000", fill = "#d6d6d6", size = 0.3),
#         strip.text = element_text(colour = "#000000", size = 12, face = "bold"),
#         legend.position = "top",
#         legend.margin = margin(t = 0, b = 0, r = 0, l = 0),
#         legend.box.margin = margin(t = 10, b = 20, r = 0, l = 0),
#         legend.key = element_rect(fill = NA),
#         legend.background =element_blank())
# inversion_heterozygosity_plot
# write_tsv(inversion_pca_table, "~/Desktop/Scripts/Data/PCA_inv/8sept22_subset/scaffold4inversion_genotype.tsv")
# ggsave(inversion_heterozygosity_plot, file = "~/Desktop/Scripts/Flat_oysters/04_local_R/03_results/PCA_inv/8sept22_subset/Scaffold4inversion_inversion_heterozygosity_plot.pdf", device = cairo_pdf, scale = 1.1, width = 12, height = 8, dpi = 300)
# dev.off()
