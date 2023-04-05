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
pca1 <- PCA(cov_mat, annot$V1, annot$V2, 1, 2, show.ellipse = T, show.label = T)
pca2 <- PCA(cov_mat, annot$V1, annot$V2, 2, 3, show.ellipse = FALSE, show.label = FALSE)
pca3 <- PCA(cov_mat, annot$V1, annot$V2, 1, 3, show.ellipse = FALSE, show.label = FALSE)

ggsave(pca3, file = "~/Desktop/Scripts/EUostrea/Figures/PopulationStructure/PCA/15mar23_prunedLDminweight0.5covMat1vs3.pdf", device = cairo_pdf, scale = 1.1, width = 12, height = 8, dpi = 300)
#dev.off()

#Plot genome-wide PCoA with the ibsMat matrix
ibs_mat <- read_tsv("~/Desktop/Scripts/Data/PopStruct_EUostrea/28jan23_prunedLDminweight0.5_PopStruct.ibsMat", col_names = F) %>% 
  dplyr::select(1:nrow(.)) %>%
  as.matrix()
PCoA(ibs_mat, annot$V1, annot$V2, 153, 1, 4, show.ellipse = F)


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



##### Plot individual labels for PC4 and PC9 outliers of the WADD population #####
# Runs PCA ~
PCA <- eigen(cov_mat)
# Merges the first 3 PCs with annot ~
PCA_Annot <- as.data.frame(cbind(annot, PCA$vectors[, c(1:10)]))
colnames(PCA_Annot) <- c("Sample_ID", "population", "PCA_1", "PCA_2", "PCA_3", "PCA_4", "PCA_5", "PCA_6", "PCA_7", "PCA_8", "PCA_9", "PCA_10")
# Gets Eigenvalues of each Eigenvectors ~
PCA_Eigenval_Sum <- sum(PCA$values)
varPC1 <-(PCA$values[1]/PCA_Eigenval_Sum)*100
varPC2 <-(PCA$values[2]/PCA_Eigenval_Sum)*100
varPC3 <-(PCA$values[3]/PCA_Eigenval_Sum)*100
varPC4 <-(PCA$values[4]/PCA_Eigenval_Sum)*100
varPC5 <-(PCA$values[5]/PCA_Eigenval_Sum)*100
varPC6 <-(PCA$values[6]/PCA_Eigenval_Sum)*100
varPC7 <-(PCA$values[7]/PCA_Eigenval_Sum)*100
varPC8 <-(PCA$values[8]/PCA_Eigenval_Sum)*100
varPC9 <-(PCA$values[9]/PCA_Eigenval_Sum)*100
varPC10 <-(PCA$values[10]/PCA_Eigenval_Sum)*100

WADD_PC4 <- ggplot(PCA_Annot, aes(PCA_1, PCA_4, label = ifelse(population == "WADD", Sample_ID, ""))) +
  geom_point() +
  geom_label_repel(aes(fill = population), max.overlaps = Inf) +
  scale_x_continuous(paste("PC1, ", round(varPC1, digits = 3), "% variation", sep = ""),
                     expand = c(0.015, 0.015)) +
  scale_y_continuous(paste("PC4, ", round(varPC4, digits = 3), "% variation", sep = ""),
                     expand = c(0.015, 0.015)) +
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
  guides(fill = guide_legend(title = "WADD population on PC4", title.theme = element_text(size = 15, face = "bold"),
                             label.theme = element_text(size = 14)))
ggsave(filename  = "EUostrea/Figures/PopulationStructure/WADD_PCA/WADD_pc4_PCA_LDpruned_15feb23.pdf", 
       plot=WADD_PC4, width = 60, height = 40, units = "cm", pointsize = 20, dpi = 300)

WADD_PC9 <- ggplot(PCA_Annot, aes(PCA_1, PCA_9, label = ifelse(population == "WADD", Sample_ID, ""))) +
  geom_point() +
  geom_label_repel(aes(fill = population), max.overlaps = Inf) +
  scale_x_continuous(paste("PC1, ", round(varPC1, digits = 3), "% variation", sep = ""),
                     expand = c(0.015, 0.015)) +
  scale_y_continuous(paste("PC9, ", round(varPC9, digits = 3), "% variation", sep = ""),
                     expand = c(0.015, 0.015)) +
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
  guides(fill = guide_legend(title = "WADD population on PC9", title.theme = element_text(size = 15, face = "bold"),
                             label.theme = element_text(size = 14)))
ggsave(filename  = "EUostrea/Figures/PopulationStructure/WADD_PCA/WADD_pc9_PCA_LDpruned_15feb23.pdf", 
       plot=WADD_PC9, width = 60, height = 40, units = "cm", pointsize = 20, dpi = 300)



# Missing Data ~
MissingData <- read.table("Data/PopStruct_EUostrea/30jan23_prunedLDminweight0.5_PopStruct.GL-MissingData.txt", sep = "\t", header = FALSE)

colnames(MissingData) <- c("Sample_ID", "Sample_ID2", "NumberMissing", "PercentageMissing")

# Runs PCA ~
PCA <- eigen(cov_mat)
# Merges the first 3 PCs with annot ~
PCA_Annot <- as.data.frame(cbind(annot, PCA$vectors[, c(1:3)]))
colnames(PCA_Annot) <- c("Sample_ID", "Population", "PCA_1", "PCA_2", "PCA_3")

# Binds the 2 DFs based on common columns ~
fulldf <- merge(PCA_Annot, MissingData, by = "Sample_ID")

# PercentageMissing as Numeric ~
fulldf$PercentageMissing <- as.numeric(as.character(fulldf$PercentageMissing))


# Expands MissingData by adding MissingCategory ~
fulldf$MissingCategory <- ifelse(fulldf$PercentageMissing <= 10, "< 10%",
                                 ifelse(fulldf$PercentageMissing <= 15, "< 15%",
                                        ifelse(fulldf$PercentageMissing <= 25, "< 25%",
                                               ifelse(fulldf$PercentageMissing <= 40, "< 40%",
                                                      ifelse(fulldf$PercentageMissing <= 50, "< 50%",
                                                             ifelse(fulldf$PercentageMissing <= 60, "< 60%", "> 60%"))))))


# Gets Eigenvalues of each Eigenvectors ~
PCA_Eigenval_Sum <- sum(PCA$values)
varPC1 <-(PCA$values[1]/PCA_Eigenval_Sum)*100
varPC2 <-(PCA$values[2]/PCA_Eigenval_Sum)*100
varPC3 <-(PCA$values[3]/PCA_Eigenval_Sum)*100

#### Creates PCA plot Missing Data ~ #### 
PCA_MissData1_3 <-ggplot(fulldf, aes_string(x = "PCA_1", y = "PCA_3", fill = "MissingCategory")) +
  geom_point(alpha = .9, size = 2.75, shape = 21, colour = "#000000") +
  scale_fill_manual(values = c("#fef0d9", "#fdd49e", "#fdbb84", "#fc8d59", "#e34a33", "#b30000", "#000000")) +
  scale_x_continuous(paste("PC1, ", round(varPC1,digits =3), "% variation", sep=""),
                     #breaks = c(0.99, 1, 1.01),
                     #labels = c("0.99", "1", "1.01"),
                     #limits = c(-0.25, 0.15),
                     expand = c(.015, .015)) +
  scale_y_continuous(paste("PC3, ", round(varPC3,digits =3), "% variation", sep=""),
                     #breaks = c(-0.05, -0.025, 0, 0.025, 0.05), 
                     #labels = c("-0.05", "-0.025", "0", "0.025", "0.05"), 
                     #limits = c(-0.0525, 0.0525),
                     expand = c(.015, .015)) +
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
  guides(fill = guide_legend(title = "Missingness", title.theme = element_text(size = 15, face = "bold"),
                             label.theme = element_text(size = 14)))
last_plot()

# Saves plot ~
ggsave(PCA_MissData, file = "~/Desktop/Scripts/EUostrea/Figures/PopulationStructure/PCA/21mar23_prunedLDminweight0.5covMat1vs2--MissingCategory.pdf", device = cairo_pdf, scale = 1.1, width = 12, height = 8, dpi = 600)
ggsave(PCA_MissData1_3, file = "~/Desktop/Scripts/EUostrea/Figures/PopulationStructure/PCA/21mar23_prunedLDminweight0.5covMat1vs3--MissingCategory.pdf", device = cairo_pdf, scale = 1.1, width = 12, height = 8, dpi = 600)







##### Mtgenome #####
cov_mat <- as.matrix(read.table("~/Desktop/Scripts/Data/MtGenome_EUostrea/MtOXgenome_HapNetwork_24mar23_minInd210.covMat")) 
bamlist <- read.table("~/Desktop/Scripts/Data/MtGenome_EUostrea/Mtbamlist_24mar23.txt", as.is = TRUE)
bams <- bamlist[, 1]
bams <- gsub(".bam", "", bams)
bams <- gsub(".+/", "", bams)
# Defining the population variable
population <- ifelse(grepl("^Lurida", bams), substr(bams, 1, 6), substr(bams, 1, 4))
# Defining the individual identifier
ind <- ifelse(grepl("^Lurida", bams), substr(bams, 1, 9), substr(bams, 1, 7))
# Combining population and ind in the same table
Mt_annot <- data.frame(population = population, ind = ind)
count(unique(Mt_annot$population))
source("~/Desktop/Scripts/Flat_oysters/04_local_R/00_scripts/individual_pca_functions_hjam_dec22_MTDNA.R")
# Reorders Population ~
Mt_annot$population <- factor(Mt_annot$population, ordered = T,
                                 levels = c("MOLU", "ZECE", "CRES",
                                            "ORIS","CORS", "PONT",  "RIAE",
                                            "MORL",
                                            "TOLL", "COLN", "BARR",
                                            "TRAL", "CLEW",
                                            "RYAN",
                                            "GREV", "WADD",
                                            "NISS","LOGS","VENO", "HALS", "THIS",
                                            "KALV", "HYPP",
                                            "LANG", "BUNN", "DOLV", "HAUG", "HAFR",
                                            "INNE","VAGS", "AGAB", "OSTR"))

#Plot genome-wide PCA with the covMat matrix
pca1 <- PCA(cov_mat, Mt_annot$ind, Mt_annot$population, 1, 2, show.ellipse = T, show.label = F)
ggsave(pca1, file = "~/Desktop/Scripts/EUostrea/Figures/Mtgenome/24mar23_MTDNA_401snps_PCA1vs2.pdf", device = cairo_pdf, scale = 1.1, width = 12, height = 8, dpi = 600)
























##### Scandinavia #####
## Plotting CovMat and IbsMat from ANGSD
cov_mat <- as.matrix(read.table("~/Desktop/Scripts/Data/PopStruct_EUostrea/SCAND/1mar23_prunedLDminweight0.5snps_SCAND.covMat")) 
bamlist <- read.table("~/Desktop/Scripts/EUostrea/01_infofiles/bamlist_EUostrea_Scandinavia.txt", as.is = TRUE)
bams <- bamlist[, 1]
bams <- gsub(".bam", "", bams)
bams <- gsub(".+/", "", bams)
# Defining the population variable
population <- ifelse(grepl("^Lurida", bams), substr(bams, 1, 6), substr(bams, 1, 4))
# Defining the individual identifier
ind <- ifelse(grepl("^Lurida", bams), substr(bams, 1, 9), substr(bams, 1, 7))
# Combining population and ind in the same table
Scand_annot <- data.frame(population = population, ind = ind)
count(unique(Scand_annot$population))
source("~/Desktop/Scripts/Flat_oysters/04_local_R/00_scripts/individual_pca_functions_hjam_mar23_SCAND.R")
# Reorders Population ~
Scand_annot$population <- factor(Scand_annot$population, ordered = T,
                   levels = c("GREV", "WADD",
                              "NISS","LOGS","VENO", "HALS", "THIS",
                              "KALV", "HYPP",
                              "LANG", "BUNN", "DOLV", "HAUG", "HAFR",
                              "INNE","VAGS", "AGAB", "OSTR"))

#Plot genome-wide PCA with the covMat matrix
pca1 <- PCA(cov_mat, Scand_annot$ind, Scand_annot$population, 1, 2, show.ellipse = T, show.label = T)
ggsave(pca1, file = "~/Desktop/Scripts/EUostrea/Figures/PopulationStructure/PCA_SCAND/pca1vs2_scand_24mar23.pdf", device = cairo_pdf, width = 16, height = 8, dpi = 600)
pca2 <- PCA(cov_mat, Scand_annot$ind, Scand_annot$population, 1, 3, show.ellipse = T, show.label = T)
ggsave(pca2, file = "~/Desktop/Scripts/EUostrea/Figures/PopulationStructure/PCA_SCAND/pca1vs3_scand_24mar23.pdf", device = cairo_pdf, width = 16, height = 8, dpi = 600)
pca3 <- PCA(cov_mat, Scand_annot$ind, Scand_annot$population, 1, 4, show.ellipse = FALSE, show.label = FALSE)




















#
##
### The END ~~~~~

