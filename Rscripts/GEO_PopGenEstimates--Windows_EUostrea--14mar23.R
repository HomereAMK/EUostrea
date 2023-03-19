### The BEGINNING ~~~~~
##
# ~ Plots BSG_Turbot--PopGenEstimates | designed George Pacheco modified by Homère J. Alves Monteiro


# Cleans the environment ~ 
rm(list=ls())


# List all the files in the working directory
setwd("~/Desktop/Scripts/Data/PopGenEstimates_Windows_EUostrea/")


# Loads required packages ~
pacman::p_load(tidyverse, lemon, extrafont)

# List all the files in the working directory
filenames <- list.files()

# Create an empty data frame to store merged data
merged_data <- data.frame()

# Load each file into a separate data frame and add a "Pops" column
for (filename in filenames) {
  # Extract the population name from the filename
  pop_name <- gsub(".tsv", "", unlist(strsplit(filename, "_"))[4])
  
  # Load the file into a data frame
  df <- read.table(filename, header = TRUE)
  
  # Add a "Pops" column to the data frame
  df$Pops <- factor(paste(pop_name))
  
  # Order CHR
  df$CHR <- factor(df$CHR, ordered = TRUE,
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
  
  
  # Append the data frame to the merged data
  merged_data <- bind_rows(merged_data, df)
  
  # Assign the data frame to a variable with the appropriate name
  assign(pop_name, df)
  
}



y_strip_labels <- c("scaffold1" = "CHR 01", "scaffold2" = "CHR 02", "scaffold3" = "CHR 03", "scaffold4" = "CHR 04",
                    "scaffold5" = "CHR 05", "scaffold6" = "CHR 06", "scaffold7" = "CHR 07", "scaffold8" = "CHR 08",
                    "scaffold9" = "CHR 09", "scaffold10" = "CHR 10")


# Gets column names ~
PopGenEstimates_Windows_ColumnNames <- colnames(BARR)
# Merges DFs ~
fulldf <- full_join(NISS, OSTR, by = PopGenEstimates_Windows_ColumnNames)

fulldf_InvReg <- filter(fulldf, 
                        (CHR == "scaffold4" | CHR == "scaffold5" | CHR == "scaffold8"))

# Creates Tp plot for USAM NISS OSTR ~
π_NISS_OSTR_InvReg  <- ggplot() +
  geom_line(data = fulldf_InvReg, aes(x = gPoint, y = Tp, colour = Pops), linetype = 1, size = 1, alpha=0.4) +
  facet_rep_grid(CHR ~. , scales = "free", labeller = labeller(CHR = y_strip_labels)) +
  geom_hline(yintercept = 0, linetype = "dotted", size = .2, color = "#FF6545") +
  scale_x_continuous("Genomic Position",
                     breaks = c( 15000000, 30000000, 40000000, 50000000, 60000000, 70000000, 80000000, 90000000, 100000000, 110000000, 115000000), 
                     labels = c( "15Mb", "30Mb", "40Mb", "50Mb", "60Mb", "70Mb", "80Mb","90Mb","100Mb", "110Mb", "115Mb"),
                     limits = c(0,115000000 ),
                     expand = c(0.01, 0.01)) +
  scale_y_continuous("Nucleotide Diversity Across Chromosomes (20k windows)",
                     breaks = c(.005, .015, 0.025), 
                     labels = c(".005", ".015", ".025" ),
                     limits = c(0, .025),
                     expand = c(0.0025, 0.0025)) +
  scale_colour_manual(values = c( "#02630C", "#240377")) +
  theme(panel.background = element_rect(fill = "#ffffff"),
        panel.border = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "#000000", size = .3),
        axis.title.x = element_text(size = 20, face = "bold", color = "#000000", margin = margin(t = 30, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size = 20, face = "bold", color = "#000000", margin = margin(t = 0, r = 30, b = 0, l = 0)),
        axis.text = element_text(colour = "#000000", size = 15, face = "bold"),
        axis.ticks = element_line(color = "#000000", size = .3),
        strip.background.y = element_rect(colour = "#000000", fill = "#d6d6d6", size = .3),
        strip.text = element_text(colour = "#000000", size = 11.5, face = "bold"),
        legend.position = "top",
        legend.margin = margin(t = 0, b = 0, r = 0, l = 0),
        legend.box.margin = margin(t = 30, b = 25, r = 0, l = 0),
        legend.key = element_rect(fill = NA),
        legend.background =element_blank()) +
  guides(colour = guide_legend(title = "Populations:", title.theme = element_text(size = 21, face = "bold"),
                               label.theme = element_text(size = 19), override.aes = list(size = 1.4)),
         fill = "none")

# Saves plot ~
ggsave(π_NISS_OSTR_InvReg, file = "~/Desktop/Scripts/EUostrea/Figures/PopGenEstimates--Windows/π_NISS_OSTR_CHR4_5_8--14mar23.pdf", 
       device = cairo_pdf, scale = 1, width = 26, height = 30, dpi = 600)
dev.off()
  