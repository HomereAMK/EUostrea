### The BEGINNING ~~~~~
##
# ~ Plots Sanity-check/ dendograms | First written by Misha Matz  with later modifications by  Homère J. Alves Monteiro


# Cleans the environment ~ 
rm(list=ls())

# Sets working directory ~
setwd("~/Desktop/Scripts/")

# Loads required packages ~
pacman::p_load(vegan, tidyverse, RcppCNPy, pheatmap, extrafont, ggforce, ggrepel, ggstar, RcppCNPy)

# Bams file
bams=read.table("Data/PCA/SanityCheck_chr1_EUostrea/bamlist_EUostrea_replicates_corr.txt")[,1] # list of bam files
bams0=bams
# removing leading / trailing filename bits from sample names
bams=sub(".bam","",bams)
bams=sub(".+/","",bams)
length(bams)
bams

#------ PCAngsd magic
# population designations (first letter of sample names)
pops=substr(bams,0,4)
pops
indivi=substr(bams,0,30)
indivi
# reading pcangsd covariance matrix, converting to correlation-based distances
pcc = as.matrix(read.table("Data/PCA/SanityCheck_chr1_EUostrea/Corr_SC_chr1_Tyler_minMapQ20_minInd0.25_setMinDepthInd1_setMinDepth7_rmTriallelic0.05minMaf0.05__setMaxDepth20000_SNPpval1e-6_minMaf0.05_nov22.covMat"))
dimnames(pcc)=list(indivi,indivi)
pccd=1-cov2cor(pcc)

# reading the ibs matrix
ma = as.matrix(read.table("Data/PCA/SanityCheck_chr1_EUostrea/Corr_SC_chr1_Tyler_minMapQ20_minInd0.25_setMinDepthInd1_setMinDepth7_rmTriallelic0.05minMaf0.05__setMaxDepth20000_SNPpval1e-6_minMaf0.05_nov22.ibsMat"))
dimnames(ma)=list(indivi,indivi)

# heatmaps of distances
pheatmap(ma)

# hierarchical clustering trees
pdf(file = "~/Desktop/Scripts/EUostrea/Figures/SanityCheck/Corr_DendrocovMat_SC_chr1_Tyler_minMapQ20_minInd0.25_setMinDepthInd1_setMinDepth7_rmTriallelic0.05minMaf0.05__setMaxDepth20000_SNPpval1e-6_minMaf0.05_nov22.pdf",   # The directory you want to save the file in
    width = 40, # The width of the plot in inches
    height = 10 )
hc=hclust(as.dist(ma),"ave")
plot(hc,cex=0.3) # clustering of samples by IBS 
dev.off()
pdf(file = "~/Desktop/Scripts/EUostrea/Figures/SanityCheck/Corr_DendroIbsMat_SC_chr1_Tyler_minMapQ20_minInd0.25_setMinDepthInd1_setMinDepth7_rmTriallelic0.05minMaf0.05__setMaxDepth20000_SNPpval1e-6_minMaf0.05_nov22.pdf",   # The directory you want to save the file in
    width = 40, # The width of the plot in inches
    height = 10 )
hc=hclust(as.dist(pccd),"ave")
plot(hc,cex=0.3) # clustering of samples by SNP correlation (from pcangsd) 
dev.off()


# Plotting PCA, PCoA, MDS
# PCoA
ordi=capscale(ma~1)
#ordi=capscale(pccd~1)
# eigenvectors
plot(ordi$CA$eig) 
# how many "interesting" ordination axes do we see? => 
cov_mat <- as.matrix(read.table("~/Desktop/Scripts/Data/PCA/SanityCheck_chr1_EUostrea/Corr_SC_chr1_Tyler_minMapQ20_minInd0.25_setMinDepthInd1_setMinDepth7_rmTriallelic0.05minMaf0.05__setMaxDepth20000_SNPpval1e-6_minMaf0.05_nov22.covMat")) 
#annot <- read.table("Data/PCA/SanityCheck_chr1_EUostrea/bamlist_EUostrea_replicates.annot", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
#MissingData <- read.table("Data/PCA/SanityCheck_chr1_EUostrea/SC_chr1_Tyler_minMapQ20_minInd0.25_setMinDepthInd1_setMinDepth7_rmTriallelic0.05minMaf0.05__setMaxDepth20000_SNPpval1e-6_minMaf0.05_nov22.GL-MissingData.txt", sep = "\t", header = FALSE)
colnames(MissingData) <- c("Sample_ID", "NumberMissing", "PercentageMissing")
source("~/Desktop/Scripts/Flat_oysters/04_local_R/00_scripts/individual_pca_functions_hjam.R")
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
                               "KALV", "HYPP", "HAVS",
                              "LANG", "BUNN", "DOLV", "HAUG", "HAFR",
                              "INNE","VAGS", "AGAB", "OSTR"))
#Plot genome-wide PCA with the covMat matrix
PCA(cov_mat, annot$V1, annot$V2, 1, 2, show.ellipse = F)
ggsave(SC_chr1_pca, file = "~/Desktop/Scripts/EUostrea/Figures/SanityCheck/Corr_PCA_SC_chr1_Tyler_minMapQ20_minInd0.25_setMinDepthInd1_setMinDepth7_rmTriallelic0.05minMaf0.05__setMaxDepth20000_SNPpval1e-6_minMaf0.05_nov22.pdf", device = cairo_pdf, scale = 1.1, width = 12, height = 8, dpi = 300)
dev.off()

# Missing data wrangling ~
PCA <- eigen(cov_mat)
PCA_Annot <- as.data.frame(cbind(annot, PCA$vectors[, c(1:3)]))
colnames(PCA_Annot) <- c("Sample_ID", "Population", "PCA_1", "PCA_2", "PCA_3")
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

# Creates PCA plot Missing Data ~ 
SC_chr1__MissData <-ggplot(fulldf, aes_string(x = "PCA_1", y = "PCA_2", fill = "MissingCategory")) +
  geom_point(alpha = .9, size = 2.75, shape = 21, colour = "#000000") +
  scale_fill_manual(values = c("#fef0d9", "#fdd49e", "#fdbb84", "#fc8d59", "#e34a33", "#b30000", "#7210A0")) +
  scale_x_continuous(paste("PC1, ", round(varPC1,digits =3), "% variation", sep=""),
                     #breaks = c(0.99, 1, 1.01),
                     #labels = c("0.99", "1", "1.01"),
                     #limits = c(-0.25, 0.15),
                     expand = c(.015, .015)) +
  scale_y_continuous(paste("PC2, ", round(varPC2,digits =3), "% variation", sep=""),
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
  guides(fill = guide_legend(title = "Population", title.theme = element_text(size = 15, face = "bold"),
                             label.theme = element_text(size = 14)))
ggsave(SC_chr1__MissData, file = "~/Desktop/Scripts/EUostrea/Figures/SanityCheck/MissingPCA_SC_chr1_Tyler_minMapQ20_minInd0.25_setMinDepthInd1_setMinDepth7_rmTriallelic0.05minMaf0.05__setMaxDepth20000_SNPpval1e-6_minMaf0.05_nov22.pdf", device = cairo_pdf, scale = 1.1, width = 12, height = 8, dpi = 300)
ggsave(SC_chr1__MissData, file = "~/Desktop/Scripts/EUostrea/Figures/SanityCheck/MissingPCA_SC_chr1_Tyler_minMapQ20_minInd0.25_setMinDepthInd1_setMinDepth7_rmTriallelic0.05minMaf0.05__setMaxDepth20000_SNPpval1e-6_minMaf0.05_nov22.png",  dpi = 300)
dev.off()


# Plot with labels ~
fulldf2 <- fulldf %>% filter(Population %in% c("GREV", "WADD", 
                                               "NISS","LOGS","VENO", "HALS", "THIS",
                                               "KALV", "HYPP",
                                               "LANG", "BUNN", "DOLV", "HAUG", "HAFR",  
                                               "INNE","VAGS", "AGAB", "OSTR"))


identify  <- ggplot(fulldf2, aes(PCA_1, PCA_2, label = Sample_ID)) +
  geom_point() +
  geom_label_repel(aes(fill = Population),max.overlaps=Inf, size=0.8)+
  scale_fill_manual(values = c("#91BD96", "#91BD96",
                               "#02630C","#02630C","#02630C", "#02630C", "#02630C",
                               "#45D1F7", "#45D1F7",
                               "#588cad", "#588cad", "#588cad", "#588cad", "#588cad",
                               "#240377", "#240377", "#240377", "#240377"))+
  scale_x_continuous(paste("PC1, ", round(varPC1,digits =3), "% variation", sep=""),
                     #breaks = c(0.99, 1, 1.01),
                     #labels = c("0.99", "1", "1.01"),
                     #limits = c(-0.25, 0.15),
                     expand = c(.015, .015)) +
  scale_y_continuous(paste("PC2, ", round(varPC2,digits =3), "% variation", sep=""),
                     #breaks = c(-0.05, -0.025, 0, 0.025, 0.05), 
                     #labels = c("-0.05", "-0.025", "0", "0.025", "0.05"), 
                     #limits = c(-0.0525, 0.0525),
                     expand = c(.015, .015)) +
  theme(legend.position="none")
ggsave(identify, file = "~/Desktop/Scripts/EUostrea/Figures/SanityCheck/SubsetLabelPCA_SC_chr1_Tyler_minMapQ20_minInd0.25_setMinDepthInd1_setMinDepth7_rmTriallelic0.05minMaf0.05__setMaxDepth20000_SNPpval1e-6_minMaf0.05_nov22.pdf", device = cairo_pdf, scale = 1.1, width = 12, height = 8, dpi = 600)
ev.off()



#### Plotting 32.000 snps ####
cov_mat <- as.matrix(read.table("~/Desktop/Scripts/Data/PCA/EUostrea/angsd0.937_htslib1.16_minMapQ20minQ20_minInd145.25_setMinDepthInd5_setMinDepth600setMaxDepth1200/angsd0.937_htslib1.16_minMapQ20minQ20_minInd145.25_setMinDepthInd5_setMinDepth600setMaxDepth1200.covMat")) 
annot <- read.table("../Scripts/EUostrea/01_infofiles/bamlist_EUostrea.annot", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
source("~/Desktop/Scripts/Flat_oysters/04_local_R/00_scripts/individual_pca_functions_hjam.R")
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
                               "KALV", "HYPP", "HAVS",
                              "LANG", "BUNN", "DOLV", "HAUG", "HAFR",
                              "INNE","VAGS", "AGAB", "OSTR"))
#Plot genome-wide PCA with the covMat matrix
PCA(cov_mat, annot$V1, annot$V2, 1, 2, show.ellipse = F, show.label = F)
ggsave(SC_chr1_pca, file = "~/Desktop/Scripts/EUostrea/Figures/SanityCheck/Corr_PCA_SC_chr1_Tyler_minMapQ20_minInd0.25_setMinDepthInd1_setMinDepth7_rmTriallelic0.05minMaf0.05__setMaxDepth20000_SNPpval1e-6_minMaf0.05_nov22.pdf", device = cairo_pdf, scale = 1.1, width = 12, height = 8, dpi = 300)
dev.off()

#### Plotting 10M snps ####
cov_mat <- as.matrix(read.table("~/Desktop/Scripts/Data/PCA/EUostrea/A940_minMapQ20minQ20_NOMININD_setMinDepthInd1_setMinDepth100setMaxDepth1400/A940_minMapQ20minQ20_NOMININD_setMinDepthInd1_setMinDepth100setMaxDepth1400.covMat")) 
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
                              "KALV", "HYPP", "HAVS",
                              "LANG", "BUNN", "DOLV", "HAUG", "HAFR",
                              "INNE","VAGS", "AGAB", "OSTR"))
#Plot genome-wide PCA with the covMat matrix
PCA(cov_mat, annot$V1, annot$V2, 1, 2, show.ellipse = F, show.label = F)
ggsave(SC_chr1_pca, file = "~/Desktop/Scripts/EUostrea/Figures/SanityCheck/Corr_PCA_SC_chr1_Tyler_minMapQ20_minInd0.25_setMinDepthInd1_setMinDepth7_rmTriallelic0.05minMaf0.05__setMaxDepth20000_SNPpval1e-6_minMaf0.05_nov22.pdf", device = cairo_pdf, scale = 1.1, width = 12, height = 8, dpi = 300)
dev.off()

#Plot genome-wide PCoA with the ibsMat matrix
ibs_mat <- read_tsv("~/Desktop/Scripts/Data/PCA/EUostrea/A940_minMapQ20minQ20_NOMININD_setMinDepthInd1_setMinDepth100setMaxDepth1400/A940_minMapQ20minQ20_NOMININD_setMinDepthInd1_setMinDepth100setMaxDepth1400.ibsMat", col_names = F) %>% 
  dplyr::select(1:nrow(.)) %>%
  as.matrix()
PCoA(ibs_mat, annot$V1, annot$V2, 153, 1, 2, show.ellipse = F)
