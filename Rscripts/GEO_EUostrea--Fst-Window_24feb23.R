### The BEGINNING ~~~~~
##
# ~ Plots Windowed Fst | By George Pacheco and Homere J. Alves Monteiro


# Cleans the environment ~ 
rm(list=ls())
# 
# Loads required packages ~
pacman::p_load(tidyverse, extrafont, lemon, ggridges, pheatmap, RColorBrewer, gtools)

# List all the files in the working directory
setwd("~/Desktop/Scripts/Data/Fst_EUostrea/Fst_SL15kb/")

# List all the files in the working directory
filenames <- list.files()

# Loop to load each file into a separate data frame and add a "Pops" column
for (filename in filenames) {
  # Extract the population name from the filename
  pop_name <- unlist(strsplit(filename, "_"))[5]
  
  # Load the file into a data frame
  df <- read.table(filename, header = FALSE)
  
  # Add colnames
  colnames(df) <- c("CHR", "SNP", "gPoint", "END", "NumberOfSites", "Fst")
  
  # Extract the gPoint_c variable
  df$gPoint_c <- seq(15000, by = 15000, length.out = nrow(df))
  
  # Add a "Pops" column to the data frame
  df$Pops <- factor(paste(pop_name))
  
  # CHR state
  df <- df[order(as.numeric(substr(df$CHR, 9, nchar(df$CHR)))), ] #super important
  df_2 <- as.data.frame(unique(df$CHR)); colnames(df_2) <- c("CHR")
  df_2$CHR_IDs <- seq.int(nrow(df_2))
  df_3 <- merge(df, df_2, by = "CHR")
  df_3 <- df_3 %>% arrange(CHR_IDs)
  df_3$CHR_State <- ifelse(df_3$CHR_IDs %% 2 == 0, "Even", "Odd")
  df_3 <- df_3[order(as.numeric(substr(df_3$CHR, 9, nchar(df_3$CHR)))), ] #super important
  
  # Assign the data frame to a variable with the appropriate name
  assign(pop_name, df_3)
  
}

# Column names 
Fst_Window_ColumnNames <- colnames(ZECE.AGAB)


# Merges DFs ~
fulldf_wadd <- rbind(ZECE.WADD, ORIS.WADD, RIAE.WADD, BARR.WADD, TRAL.WADD, RYAN.WADD, MORL.WADD, GREV.WADD, WADD.THIS, WADD.HALS, WADD.KALV, WADD.DOLV, WADD.LANG, WADD.AGAB, WADD.OSTR)
fulldf_niss <- rbind(MOLU.NISS, ORIS.NISS, RIAE.NISS, MORL.NISS, BARR.NISS, TRAL.NISS, RYAN.NISS, WADD.NISS, NISS.HALS, NISS.HYPP, NISS.BUNN, NISS.DOLV, NISS.BUNN, NISS.HAUG, NISS.INNE, NISS.OSTR)

# Reorders Species ~
y_strip_labels <- c("scaffold1" = "CHR 01", "scaffold2" = "CHR 02", "scaffold3" = "CHR 03", "scaffold4" = "CHR 04",
                    "scaffold5" = "CHR 05", "scaffold6" = "CHR 06", "scaffold7" = "CHR 07", "scaffold8" = "CHR 08",
                    "scaffold9" = "CHR 09", "scaffold10" = "CHR 10")

#### Creates Manhattan panel geom_line ~ ####
# geom_line
ggplot() +
  geom_line(data = AGAB.OSTR, aes(x = gPoint, y = Fst, colour = Pops), linetype = 1, size = .2) +
  facet_rep_grid(CHR ~. , scales = "free", labeller = labeller(CHR = y_strip_labels)) +
  scale_x_continuous("Genomic Position",
                     breaks = c( 15000000, 30000000, 40000000, 50000000, 60000000, 70000000, 80000000, 90000000, 100000000, 110000000, 115000000), 
                     labels = c( "15Mb", "30Mb", "40Mb", "50Mb", "60Mb", "70Mb", "80Mb","90Mb","100Mb", "110Mb", "115Mb"),
                     limits = c(0,115000000 ),
                     expand = c(0.01, 0.01)) +
  scale_y_continuous("Fst Across Chromosomes",
                     breaks = c(.25, .5, .75), 
                     labels = c(".25", ".50", ".75"),
                     limits = c(0, 1),
                     expand = c(0.01, 0.01)) +
  scale_colour_manual(values = c("#083359", "#BF820F",  "#1BBC9B")) +
  theme(panel.background = element_rect(fill = "#ffffff"),
        panel.border = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "#000000", size = .3),
        axis.title.x = element_text(size = 16, face = "bold", color = "#000000", margin = margin(t = 30, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size = 16, face = "bold", color = "#000000", margin = margin(t = 0, r = 30, b = 0, l = 0)),
        axis.text = element_text(colour = "#000000", size = 10),
        axis.ticks = element_line(color = "#000000", size = .3),
        strip.background.y = element_rect(colour = "#000000", fill = "#d6d6d6", size = 0.3),
        strip.text = element_text(colour = "#000000", size = 10, face = "bold"),
        legend.position = "top",
        legend.margin = margin(t = 0, b = 0, r = 0, l = 0),
        legend.box.margin = margin(t = 30, b = 25, r = 0, l = 0),
        legend.key = element_rect(fill = NA),
        legend.background =element_blank()) +
  guides(colour = guide_legend(title = "Fst Comparisons:", title.theme = element_text(size = 16, face = "bold"),
                             label.theme = element_text(size = 16), override.aes = list(size = 1.4)),
         fill = "none")

# Saves Manhattan plot ~
ggsave(Fst_Window, file = "~/Desktop/Scripts/Flat_oysters/04_local_R/03_results/FST/Local_Limfjorden/Niss_Løgst_ThiS_Hals_SLWin15K_Jun22list.pdf", device = cairo_pdf, scale = 1, width = 26, height = 20, dpi = 600)
dev.off()

# geom_area
ggplot(data = fulldf, aes(x = gPoint, y = Fst, fill = Pops)) +
  geom_area(alpha = 0.5) +
  facet_rep_grid(CHR ~. , scales = "free", labeller = labeller(CHR = y_strip_labels)) +
  scale_x_continuous("Genomic Position",
                     breaks = c( 15000000, 30000000, 40000000, 50000000, 60000000, 70000000, 80000000, 90000000, 100000000, 110000000, 115000000), 
                     labels = c( "15Mb", "30Mb", "40Mb", "50Mb", "60Mb", "70Mb", "80Mb","90Mb","100Mb", "110Mb", "115Mb"),
                     limits = c(0,115000000 ),
                     expand = c(0.01, 0.01)) +
  scale_y_continuous("Fst Across Chromosomes",
                     breaks = c(.25, .5, .75), 
                     labels = c(".25", ".50", ".75"),
                     limits = c(0, 1),
                     expand = c(0.01, 0.01)) +
  scale_fill_manual(values = c("#083359", "#BF820F",  "#1BBC9B")) +
  theme(panel.background = element_rect(fill = "#ffffff"),
        panel.border = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "#000000", size = .3),
        axis.title.x = element_text(size = 16, face = "bold", color = "#000000", margin = margin(t = 30, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size = 16, face = "bold", color = "#000000", margin = margin(t = 0, r = 30, b = 0, l = 0)),
        axis.text = element_text(colour = "#000000", size = 10),
        axis.ticks = element_line(color = "#000000", size = .3),
        strip.background.y = element_rect(colour = "#000000", fill = "#d6d6d6", size = 0.3),
        strip.text = element_text(colour = "#000000", size = 10, face = "bold"),
        legend.position = "top",
        legend.margin = margin(t = 0, b = 0, r = 0, l = 0),
        legend.box.margin = margin(t = 30, b = 25, r = 0, l = 0),
        legend.key = element_rect(fill = NA),
        legend.background =element_blank()) +
  guides(colour = guide_legend(title = "Fst Comparisons:", title.theme = element_text(size = 16, face = "bold"),
                               label.theme = element_text(size = 16), override.aes = list(size = 1.4)),
         fill = "none")

# geom_area sca1
ggplot(data = scaffold1_df, aes(x = gPoint, y = Fst, fill = Pops)) +
  geom_area(alpha = 0.8) +
  scale_x_continuous("Genomic Position",
                     breaks = c( 15000000, 30000000, 40000000, 50000000, 60000000, 70000000, 80000000, 90000000, 100000000, 110000000, 115000000), 
                     labels = c( "15Mb", "30Mb", "40Mb", "50Mb", "60Mb", "70Mb", "80Mb","90Mb","100Mb", "110Mb", "115Mb"),
                     limits = c(0,115000000 ),
                     expand = c(0.01, 0.01)) +
  scale_y_continuous("Fst Across Chromosomes",
                     breaks = c(.25, .5, .75), 
                     labels = c(".25", ".50", ".75"),
                     limits = c(0, 1),
                     expand = c(0.01, 0.01)) +
  scale_fill_manual(values = c("#083359", "#BF820F",  "#1BBC9B")) 
  

#### Creates Manhattan panel geom_point ~ 



#
##
### The END ~~~~~






####

# Wadd full comparison ~
fulldf_wadd$Pops <- factor(fulldf_wadd$Pops, ordered = TRUE,
                           levels = c("ZECE.WADD", "ORIS.WADD",  "RIAE.WADD", "MORL.WADD", "BARR.WADD", "TRAL.WADD", "RYAN.WADD", "GREV.WADD", "WADD.THIS", 
                                      "WADD.HALS", "WADD.KALV", "WADD.DOLV", "WADD.LANG", "WADD.AGAB", "WADD.OSTR"))

fulldf_wadd$CHR_State <- as.character(factor(fulldf_wadd$CHR_State, ordered = TRUE,
                                             levels = c("Odd", "Even")))

fulldf_wadd$CHR <- factor(fulldf_wadd$CHR, ordered = TRUE,
                          levels = c("scaffold1",
                                     "scaffold2",
                                     "scaffold3",
                                     "scaffold4",
                                     "scaffold5",
                                     "scaffold6",
                                     "scaffold7",
                                     "scaffold8",
                                     "scaffold9",
                                     "scaffold10"))
fullcompWadd <- ggplot() +
  geom_point(data = fulldf_wadd,
             aes(x = gPoint_c, y = Fst,  colour = CHR_State, fill= CHR), shape = 21, size = 0.6, alpha = 0.6) +
  facet_rep_grid(Pops~. , scales = "free_x") +
  scale_x_continuous("Chromosomes",
                     expand = c(.005, .5)) +
  scale_y_continuous("Fst (15Kb Sliding Windows)",
                     breaks = c(.30, .60, .90), 
                     labels = c(".30", ".60", ".90"),
                     limits = c(0, .99),
                     expand = c(0.01, 0.01)) +
  scale_colour_brewer() +
  scale_fill_manual(values = c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00", "#ffff33", "#a65628", "#f781bf", "#999999", "#8dd3c7")) +
  theme(panel.background = element_rect(fill = "#ffffff"),
        panel.border = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "#000000", size = .3),
        axis.title.x = element_text(size = 8, face = "bold", color = "#000000", margin = margin(t = 30, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size = 8, face = "bold", color = "#000000", margin = margin(t = 0, r = 30, b = 0, l = 0)),
        axis.text.x = element_text(colour = "#000000", size = 4),
        axis.text.y = element_text(colour = "#000000", size = 8, face = "bold"),
        axis.ticks.x = element_line(color = "#000000", size = .3),
        axis.ticks.y = element_line(color = "#000000", size = .3),
        strip.background.y = element_rect(colour = "#000000", fill = "#FFFFFF", size = .5),
        strip.text = element_text(colour = "#000000", size = 7, face = "bold"),
        legend.position = "top",
        legend.margin = margin(t = 0, b = 0, r = 0, l = 0),
        legend.box.margin = margin(t = 30, b = 25, r = 0, l = 0),
        legend.key = element_rect(fill = NA),
        legend.background = element_blank()) +
  guides(colour = "none", fill = "none")
fullcompWadd
ggsave(filename  = "~/Desktop/Scripts/EUostrea/Figures/Fst/15kbSLwin/Wadd_fullcomp_15kbwin_27feb23.png", 
       plot=fullcompWadd, width = 40, height = 50, units = "cm", pointsize = 20, dpi = 250)
dev.off()

# Niss full comparison ~
fulldf_niss$Pops <- factor(fulldf_niss$Pops, ordered = TRUE,
                           levels = c("MOLU.NISS", "ORIS.NISS", "RIAE.NISS", "MORL.NISS", "BARR.NISS", "TRAL.NISS", 
                                      "RYAN.NISS", "WADD.NISS", "NISS.HALS", "NISS.HYPP", "NISS.BUNN", "NISS.DOLV", "NISS.HAUG", "NISS.INNE", "NISS.OSTR"))
fullcompNiss <- ggplot() +
  geom_point(data = fulldf_niss,
             aes(x = gPoint_c, y = Fst,  colour = CHR_State, fill= CHR), shape = 21, size = 0.6, alpha = 0.6) +
  facet_rep_grid(Pops~. , scales = "free_x") +
  scale_x_continuous("Chromosomes",
                     expand = c(.005, .5)) +
  scale_y_continuous("Fst (15Kb Sliding Windows)",
                     breaks = c(.30, .60, .90), 
                     labels = c(".30", ".60", ".90"),
                     limits = c(0, .99),
                     expand = c(0.01, 0.01)) +
  scale_colour_brewer() +
  scale_fill_manual(values = c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00", "#ffff33", "#a65628", "#f781bf", "#999999", "#8dd3c7")) +
  theme(panel.background = element_rect(fill = "#ffffff"),
        panel.border = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "#000000", size = .3),
        axis.title.x = element_text(size = 8, face = "bold", color = "#000000", margin = margin(t = 30, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size = 8, face = "bold", color = "#000000", margin = margin(t = 0, r = 30, b = 0, l = 0)),
        axis.text.x = element_text(colour = "#000000", size = 4),
        axis.text.y = element_text(colour = "#000000", size = 8, face = "bold"),
        axis.ticks.x = element_line(color = "#000000", size = .3),
        axis.ticks.y = element_line(color = "#000000", size = .3),
        strip.background.y = element_rect(colour = "#000000", fill = "#FFFFFF", size = .5),
        strip.text = element_text(colour = "#000000", size = 7, face = "bold"),
        legend.position = "top",
        legend.margin = margin(t = 0, b = 0, r = 0, l = 0),
        legend.box.margin = margin(t = 30, b = 25, r = 0, l = 0),
        legend.key = element_rect(fill = NA),
        legend.background = element_blank()) +
  guides(colour = "none", fill = "none")
ggsave(filename  = "~/Desktop/Scripts/EUostrea/Figures/Fst/15kbSLwin/Niss_fullcomp_15kbwin_27feb23.png", 
       plot=fullcompNiss, width = 40, height = 50, units = "cm", pointsize = 20, dpi = 250)
dev.off()

##
### The END ~~~~~