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


## plot with 10 axis
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

## Extract sites exceeding  p-value cut off of 0.05 after Bonferroni correction and export them
genome_selection_sites_outliers <- filter(genome_selection_sites, neg_log_p_value > -log(0.05/dim(genome_selection)[1]))
head(genome_selection_sites_outliers)

## For PC9, scaffold2 Extract sites exceeding  p-value cut off of 0.05 after Bonferroni correction
genome_selection_sites_outliers_PC9 <- filter(genome_selection_sites_outliers, 
                                                   (lg == "scaffold2" | lg == "scaffold10") & pc == "PC9")
# Group the data by "lg" and calculate the min and max positions
genome_selection_sites_outliers_PC9_mm <- genome_selection_sites_outliers_PC9 %>%
  group_by(lg) %>%
  summarize(min_pos = min(pos), max_pos = max(pos))

# Filter the result for "scaffold2" and "scaffold10"
result_filtered <- genome_selection_sites_outliers_PC9_mm %>%
  filter(lg %in% c("scaffold2", "scaffold10"))

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



