### The BEGINNING ~~~~~
##
# ~ Plots ngsAdmix EUostrea | First written by Jose Samaniego with later modifications by George Pacheco and later later modifications by Homère J. Alves Monteiro.

# Cleans the environment ~ 
rm(list=ls())


# Sets working directory ~
#setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd("~/Desktop/Scripts/Data/NgsAdmix_EUostrea/SCAND/")

# Loads required packages ~
pacman::p_load(tidyverse, MetBrewer, scales, optparse, plyr, RColorBrewer, extrafont, gtable, grid, mdthemes, ggtext, glue)


# Creates colour palette ~

nb.cols <- 10
MyColours <- colorRampPalette(brewer.pal(9, "Paired"))(nb.cols)


# Loads the data ~
samples <- read.table("~/Desktop/Scripts/Data/NgsAdmix_EUostrea/prunedLDminweight0.5_SNPs_13mar23-AllSamples--popfile.txt", stringsAsFactors = FALSE, sep = "\t")


# Reads the annotation file ~
ids <- read.table("~/Desktop/Scripts/EUostrea/01_infofiles/bamlist_EUostrea.annot", stringsAsFactors = FALSE, sep = "\t", header = FALSE)


# Adds column ids names ~
colnames(ids) <- c("Sample_ID", "Population")



# Ask Sama ~
fulldf <- data.frame()


# Ask Sama 2 ~
x <- list(c(9,5,3,6,4,8,2,10,7,1),
  c(2,3,5,9,4,7,6,8,1),
  c(5,7,6,2,1,4,8,3),
  c(4,1,5,6,2,3,7),
  c(6,5,3,1,2,4),
  c(3,1,2,5,4),
  c(3,1,2,4),
  c(1,3,2),
  c(1,2))

# Defines samples' IDs ~
sampleid = "Sample_ID"


# Ask Sama 3 ~
for (j in 1:length(samples[,1])){
  data <- read.table(samples[j,1])[, x[[j]]]
  for (i in 1:dim(data)[2]) { 
    temp <- data.frame(Ancestry = data[, i])
    temp$K <- as.factor(rep(i, times = length(temp$Ancestry)))
    temp[sampleid] <- as.factor(ids[sampleid][,1])
    temp$K_Value <- as.factor(rep(paste("K = ", dim(data)[2], sep = ""), times = length(temp$Ancestry)))
    temp <- merge(ids, temp)
    fulldf <- rbind(fulldf, temp)}}


# Specifies xbalebels colours ~
#fulldf$xlabelscolours <- ifelse(fulldf$PercentageMissing > 50, "red", "black")



# Reorders Population ~
fulldf$Population <- factor(fulldf$Population, ordered = T,
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



#fulldf$xlabelscolours_Ordered <- fulldf$xlabelscolours[order(fulldf$Population)]


# Defines the target to be plotted ~
target = "Population"

# Creates the plots ~
ngsAdmix <-
 ggplot(fulldf, aes(x = Sample_ID, y = Ancestry, fill = K)) +
  geom_bar(stat = "identity", width = .85) +
   facet_grid(K_Value ~ get(target), space = "free_x", scales = "free_x") +
   scale_x_discrete(expand = c(0, 0)) + 
   scale_y_continuous(expand = c(0, 0), breaks = NULL) +
   scale_fill_manual(values = MetBrewer::met.brewer("Redon", n = 10, type = "discrete")) +
  theme(panel.background = element_rect(fill = "#ffffff"),
        panel.border = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = margin(t = 0.005, b = 0.005, r = .4, l = .4, unit = "cm"),
        axis.line = element_blank(),
        axis.text = element_text(colour = "#000000", size = 10, face = "bold"),
        axis.text.x.bottom = element_text(colour = "#000000", face = "bold", angle = 90, vjust = .5, hjust = .5, size = 1),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        strip.background = element_rect(colour = c("#A02353", "#A02353", "#AD5B35", "#ad7358",
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
                                                   "#240377", "#240377", "#240377", "#240377"), fill = "#FAFAFA", size = .05),
        strip.text.x = element_text(colour = "black", face = "bold", size = 10, angle = 90, margin = margin(.75, 0, .75, 0, "cm")),
        legend.position = "right",
        legend.key = element_rect(fill = NA),
        legend.background = element_blank(),
        legend.title = element_text(colour = "#000000", size = 10, face = "bold"),
        legend.text = element_text(colour = "#000000", size = 10, face="bold"))


ngsAdmix

# Saves the final plot ~
ggsave(ngsAdmix, file = "~/Desktop/Scripts/EUostrea/Figures/PopulationStructure/Admixture/prunedLDminweight0.5_13mar23--ngsAdmix_OdilonRedon.pdf", device = cairo_pdf, width = 16, height = 8, dpi = 600)











#### SCAND admixture #####
# Cleans the environment ~ 
rm(list=ls())

# Creates colour palette ~
setwd("~/Desktop/Scripts/Data/NgsAdmix_EUostrea/SCAND/")

# Loads the data ~
samples <- read.table("~/Desktop/Scripts/Data/NgsAdmix_EUostrea/SCAND/qopt_popfile_scand.txt", stringsAsFactors = FALSE, sep = "\t")

# Reads the annotation file ~
ids <- read.table("~/Desktop/Scripts/EUostrea/01_infofiles/bamlist_EUostrea_Scandinavia.txt", stringsAsFactors = FALSE, sep = "\t", header = FALSE)
bams <- ids[, 1]
bams <- gsub(".bam", "", bams)
bams <- gsub(".+/", "", bams)
# Defining the population variable
population <- ifelse(grepl("^Lurida", bams), substr(bams, 1, 6), substr(bams, 1, 4))
# Defining the individual identifier
ind <- ifelse(grepl("^Lurida", bams), substr(bams, 1, 9), substr(bams, 1, 7))
# Combining population and ind in the same table
ids <- data.frame(population = population, ind = ind)
# Adds column ids names ~
colnames(ids) <- c( "Population","Sample_ID")
(unique(ids$Population))

# Ask Sama ~
fulldf <- data.frame()
# Ask Sama 2 ~
x <- list(c(4,3,7,2,6,8,5,1),
  c(6,2,3,7,1,4,5),
  c(2,5,1,6,3,4),
  c(3,4,5,2,1),
  c(4,3,2,1),
  c(3,1,2),
  c(1,2))

# Defines samples' IDs ~
sampleid = "Sample_ID"
# Ask Sama 3 ~
for (j in 1:length(samples[,1])){
  data <- read.table(samples[j,1])[, x[[j]]]
  for (i in 1:dim(data)[2]) { 
    temp <- data.frame(Ancestry = data[, i])
    temp$K <- as.factor(rep(i, times = length(temp$Ancestry)))
    temp[sampleid] <- as.factor(ids[sampleid][,1])
    temp$K_Value <- as.factor(rep(paste("K = ", dim(data)[2], sep = ""), times = length(temp$Ancestry)))
    temp <- merge(ids, temp)
    fulldf <- rbind(fulldf, temp)}}

# Reorders Population ~
fulldf$Population <- factor(fulldf$Population, ordered = T,
                            levels = c("GREV", "WADD",
                                       "NISS","LOGS","VENO", "HALS", "THIS",
                                       "KALV", "HYPP",
                                       "LANG", "BUNN", "DOLV", "HAUG", "HAFR",
                                       "INNE","VAGS", "AGAB", "OSTR"))



#fulldf$xlabelscolours_Ordered <- fulldf$xlabelscolours[order(fulldf$Population)]


# Defines the target to be plotted ~
target = "Population"

# Creates the plots ~
ngsAdmix <-
  ggplot(fulldf, aes(x = Sample_ID, y = Ancestry, fill = K)) +
  geom_bar(stat = "identity", width = .85) +
  facet_grid(K_Value ~ get(target), space = "free_x", scales = "free_x") +
  scale_x_discrete(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0), breaks = NULL) +
  scale_fill_manual(values = MetBrewer::met.brewer("Redon", n = 8, type = "discrete")) +
  theme(panel.background = element_rect(fill = "#ffffff"),
        panel.border = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = margin(t = 0.005, b = 0.005, r = .4, l = .4, unit = "cm"),
        axis.line = element_blank(),
        axis.text = element_text(colour = "#000000", size = 10, face = "bold"),
        axis.text.x.bottom = element_text(colour = "#000000", face = "bold", angle = 90, vjust = .5, hjust = .5, size = 1),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        strip.background = element_rect(colour = c(
                                                   "#91BD96", "#91BD96",
                                                   "#02630C", "#025b0b", "#02530a", "#024b09", "#024308",
                                                   "#45D1F7", "#45D1F7",
                                                   "#588cad", "#4d7f96", "#416c7e", "#35586c", "#294458",
                                                   "#240377", "#1e026d", "#1a0163", "#160059"), fill = "#FAFAFA", size = .05),
        strip.text.x = element_text(colour = "black", face = "bold", size = 10, angle = 90, margin = margin(.75, 0, .75, 0, "cm")),
        legend.position = "right",
        legend.key = element_rect(fill = NA),
        legend.background = element_blank(),
        legend.title = element_text(colour = "#000000", size = 10, face = "bold"),
        legend.text = element_text(colour = "#000000", size = 10, face="bold"))

ngsAdmix
ggsave(ngsAdmix, file = "~/Desktop/Scripts/EUostrea/Figures/PopulationStructure/Admixture_SCAND/K2_scand_24mar23.pdf", device = cairo_pdf, width = 16, height = 8, dpi = 600)





















 ## evalAdmix test ##
#Assessing model fit for SCAND K=2
source("~/Desktop/Scripts/EUostrea/Rscripts/visFuns.R") # import some funcitons to help in visualization
bamlist <- read.table("~/Desktop/Scripts/EUostrea/01_infofiles/bamlist_EUostrea_Scandinavia.txt", as.is = TRUE)
bams <- bamlist[, 1]
bams <- gsub(".bam", "", bams)
bams <- gsub(".+/", "", bams)
# Defining the population variable
population <- ifelse(grepl("^Lurida", bams), substr(bams, 1, 6), substr(bams, 1, 4))
# Defining the individual identifier
ind <- ifelse(grepl("^Lurida", bams), substr(bams, 1, 9), substr(bams, 1, 7))

# Combining population and ind in the same table
Scand_data <- data.frame(population = population, ind = ind)

q<-read.table("~/Desktop/Scripts/Data/NgsAdmix_EUostrea/SCAND/1mar23_prunedLDminweight0.5snps_SCAND_NGSadmix2kiter.2.qopt")
r <- as.matrix(read.table("~/Desktop/Scripts/Data/NgsAdmix_EUostrea/SCAND/1mar23_prunedLDminweight0.5snps_SCAND_NGSadmix2kiter.2.corres"))
str(pop)
pop2 <- as.table(pop)
ord<-orderInds(Scand_data = Scand_data[,1], q=q) # sort indivduals by population and within population by admixture proportion
plotCorRes(r, Scand_data=Scand_data[,2], ord=ord, max_z = 0.2)



# # Adds grob ~
# ngsAdmix_G <- ggplotGrob(ngsAdmix)
# ngsAdmix_G <- gtable_add_rows(ngsAdmix_G, unit(1.25, "cm"), pos = 5)
# 
# 
# # Adds top strips ~
# ngsAdmix_G <- gtable_add_grob(ngsAdmix_G, list(rectGrob(gp = gpar(col = "#000000", fill = "#44AA99", alpha = .9, size = .75, lwd = 0.25)),
#                textGrob("Remote Localities Within Natural Range", gp = gpar(cex = 1.4, fontface = 'bold', col = "black"))),
#                t = 6, l = 4, b = 6, r = 17, name = c("a", "b"))
# ngsAdmix_G <- gtable_add_grob(ngsAdmix_G, list(rectGrob(gp = gpar(col = "#000000", fill = "#F0E442", alpha = .9, size = .75, lwd = 0.25)),
#                textGrob("Urban Localities Within Natural Range", gp = gpar(cex = 1.4, fontface = 'bold', col = "black"))),
#                t = 6, l = 19, b = 6, r = 49, name = c("a", "b"))
# ngsAdmix_G <- gtable_add_grob(ngsAdmix_G, list(rectGrob(gp = gpar(col = "#000000", fill = "#E69F00", alpha = .9, size = .5, lwd = 0.25)),
#                textGrob("Localities Outside Natural Range", gp = gpar(cex = 1.4, fontface = 'bold', col = "black"))),
#                t = 6, l = 51, b = 6, r = 75, name = c("a", "b"))
# ngsAdmix_G <- gtable_add_grob(ngsAdmix_G, list(rectGrob(gp = gpar(col = "#000000", fill = "#56B4E9", alpha = .9, size = .5, lwd = 0.25)),
#                textGrob("Captives", gp = gpar(cex = 1.4, fontface = 'bold', col = "black"))),
#                t = 6, l = 77, b = 6, r = 81, name = c("a", "b"))
# 
# 
# # Controls separation ~
# ngsAdmix_G <- gtable_add_rows(ngsAdmix_G, unit(2/10, "line"), 6)
# 
# 
# # Creates the final plot ~
# grid.newpage()
# grid.draw(ngsAdmix_G)
# 
# 
# # Saves the final plot ~
# ggsave(ngsAdmix_G, file = "FPG--ngsAdmix_RColours.pdf", device = cairo_pdf, width = 40, height = 15, dpi = 600)
# 
# 
# 
# 
# #### Admixture Scandinavia 
# # Cleans the environment ~ 
# rm(list=ls())
# 
# 
# # Sets working directory ~
# #setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
# setwd("~/Desktop/Scripts/Data/NGSadmix/Scandinavia_globalsnpList200k/")
# 
# # Creates colour palette ~
# nb.cols <- 9
# #MyColours <- colorRampPalette(brewer.pal(9, "Paired"))(nb.cols)
# MyColours <- colorRampPalette(MetBrewer::met.brewer("Renoir", n = 9, type = "discrete"))
# 
# 
# # Loads the data ~
# samples <- read.table("~/Desktop/Scripts/Data/NGSadmix/Scandinavia_globalsnpList200k/Admixture_scandinavia_11oct22_globalsnpList200k--poptext_file.txt", stringsAsFactors = FALSE, sep = "\t")
# 
# 
# # Reads the annotation file ~
# ids <- read.table("~/Desktop/Scripts/Data/NGSadmix/Scandinavia_globalsnpList200k/Bam_list_Scandinavia_22sept22.annot", stringsAsFactors = FALSE, sep = "\t", header = FALSE)
# 
# 
# # Adds column ids names ~
# colnames(ids) <- c("Sample_ID", "Population")
# 
# # Ask Sama ~
# fulldf <- data.frame()
# 
# x <- list(c(9,3,4,7,8,2,1,6,5),
#   c(2,5,6,1,4,3,8,7),
#   c(3,2,1,5,6,7,4),
#   c(3,5,1,4,6,2),
#   c(5,1,2,3,4),
#   c(2,3,4,1),
#   c(1,3,2),
#   c(2,1))
# 
# 
# # Defines samples' IDs ~
# sampleid = "Sample_ID"
# 
# 
# # Ask Sama 3 ~
# for (j in 1:length(samples[,1])){
#   data <- read.table(samples[j,1])[, x[[j]]]
#   for (i in 1:dim(data)[2]) { 
#     temp <- data.frame(Ancestry = data[, i])
#     temp$K <- as.factor(rep(i, times = length(temp$Ancestry)))
#     temp[sampleid] <- as.factor(ids[sampleid][,1])
#     temp$K_Value <- as.factor(rep(paste("K = ", dim(data)[2], sep = ""), times = length(temp$Ancestry)))
#     temp <- merge(ids, temp)
#     fulldf <- rbind(fulldf, temp)}}
# 
# # Reorders Population ~
# fulldf$Population <- factor(fulldf$Population, ordered = T,
#                             levels = c(
#                                        "GREV", "WADD", 
#                                        "FURI", "NISS","LOGS","VENO", "HALS", "THIS",
#                                        "HAVS", "KALV",  "HFJO", "RAMS", "ORNE", "HYPP", "LILL", "SVALL", "SYDK", "VADE", "BOVA",
#                                        "LANG", "BUNN", "DOLV", "HAUG", "HAFR",  
#                                        "INNE","VAGS", "AGAB", "OSTR"))
# 
# # Defines the target to be plotted ~
# target = "Population"
# 
# fulldf <- fulldf %>% filter(Population != "HVAD")
# #Remove KALV_11 duplicate
# fulldf <- fulldf %>% filter(Sample_ID != "KALV_11")
# 
# # Creates the plots ~
# ngsAdmix <-
#   ggplot(fulldf, aes(x = Sample_ID, y = Ancestry, fill = K)) +
#   geom_bar(stat = "identity", width = .85) +
#   facet_grid(K_Value ~ get(target), space = "free_x", scales = "free_x") +
#   scale_x_discrete(expand = c(0, 0)) + 
#   scale_y_continuous(expand = c(0, 0), breaks = NULL) +
#   scale_fill_manual(values = MetBrewer::met.brewer("Renoir", n = 9, type = "discrete")) +
#   theme(
#     panel.background = element_rect(fill = "#ffffff"),
#     panel.grid.minor.x = element_blank(),
#     panel.grid.major = element_blank(),
#     panel.spacing = unit(.2, "lines"),
#     plot.title = element_blank(),
#     axis.title = element_blank(),
#     axis.text.x.bottom = element_text(colour = "#000000", face = "bold", angle = 90, vjust = .5, hjust = .5, size = 1),
#     axis.text.y = element_blank(),
#     axis.ticks = element_blank(),
#     strip.background = element_rect(colour = "#000000", fill = "#FAFAFA", size = .05),
#     strip.text.x = element_text(colour = "#000000", face = "bold", size = 6, margin = margin(.75, 0, .75, 0, "cm")),
#     strip.text.y = element_text(colour = "#000000", face = "bold", size = 9, angle = 90, margin = margin(0, .1, 0, .1, "cm")),
#     legend.position = "none")
# # Saves the final plot ~
# ggsave(ngsAdmix, file = "~/Desktop/Scripts/Flat_oysters/04_local_R/03_results/Admixture/Admixture_scandinavia_11oct22_globalsnpList200k/Admixture_scandinavia_11oct22_globalsnpList200k.pdf", device = cairo_pdf, width = 16, height = 8, dpi = 300)


#
##
### The END ~~~~~



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20),
#c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19),
#c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18),
#c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17),
#c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16),
#c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15),
#c(1,2,3,4,5,6,7,8,9,10,11,12,13,14),
#c(1,2,3,4,5,6,7,8,9,10,11,12,13),
#c(1,2,3,4,5,6,7,8,9,10,11,12),
#c(1,2,3,4,5,6,7,8,9,10,11),
#c(1,2,3,4,5,6,7,8,9,10),
#c(1,2,3,4,5,6,7,8,9),
#c(1,2,3,4,5,6,7,8),
#c(1,2,3,4,5,6,7),
#c(1,2,3,4,5,6),
#c(1,2,3,4,5),
#c(1,2,3,4),
#c(1,2,3),
#c(1,2))
