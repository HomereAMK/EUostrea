### The BEGINNING ~~~~~
##
# ~ Plots Relatedness Genome-wide | By Homère J. Alves Monteiro

# Cleans the environment ~ 
rm(list=ls())

# Loads required packages ~
pacman::p_load(latex2exp,scales, extrafont, dplyr, grid, lubridate, cowplot, egg, tidyverse, viridis, lemon, stringr, reshape)

# List all the files in the working directory
setwd("~/Desktop/Scripts/Data/Relatedness_EUostrea/")

# List all the files in the working directory
filenames <- list.files()

# Create an empty data frame to store merged data
merged_data <- data.frame()


# Load each file into a separate data frame and add a "Pops" column
for (filename in filenames) {
  # Extract the population name from the filename
  pop_name <- gsub(".res", "", unlist(strsplit(filename, "_"))[5])
  
  # Load the file into a data frame
  df <- read.table(filename, header = TRUE)
  
  # Add a "Pops" column to the data frame
  df$Pops <- factor(paste(pop_name))
  
  # Append the data frame to the merged data
  merged_data <- bind_rows(merged_data, df)
  
  # Assign the data frame to a variable with the appropriate name
  assign(pop_name, df)
  
}

# Expands MissingData by adding MissingCategory ~
merged_data$rab_c <- ifelse(merged_data$rab <= 0.01, "0-0.01",
                       ifelse(merged_data$rab <= 0.03, "0.01-0.03",
                              ifelse(merged_data$rab <= 0.06, "0.03-0.06",
                                     ifelse(merged_data$rab <= 0.12, "0.06-0.12",
                                            ifelse(merged_data$rab <= 1, "0.12-1")))))


# First summarise and transform your data ~
merged_data_trans <- merged_data %>% 
  group_by(Pops, rab_c) %>% 
  summarise(count = n()) %>% 
  mutate(perc = count/sum(count))

# Reorders Population ~
merged_data_trans$Pops <- factor(merged_data_trans$Pops, ordered = T,
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

# Plot the per pop relatedness 
rxy_stacked_Plot <-
  ggplot(merged_data_trans, aes(x = factor(Pops), y = perc*100, fill = factor(rab_c))) +
  geom_bar(stat="identity", width = 0.7) +
  coord_flip()+
  scale_fill_grey() +
  labs(x = "Populations", y = "Proportion of pairwise relatedness", fill = "Relatedness class") +
  theme(panel.background = element_rect(fill = "#ffffff"),
        plot.margin = margin(t = 0.005, b = 0.005, r = .2, l = .2, unit = "cm"),
        axis.title = element_text(colour = "#000000", size = 12, face = "bold"),
        axis.text.x = element_text(colour = "#000000", size = 10, face = "bold", angle = 45, vjust = 1, hjust = 1),
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
                                              "#240377", "#240377", "#240377", "#240377"), size = 10, face = "bold"),
        axis.ticks = element_line(color = "#000000", size = 1),
        legend.position = "right",
        legend.margin = margin(t = 0, b = 0, r = 0, l = 0),
        legend.box.margin = margin(t = 15, b = 25, r = 0, l = 0),
        legend.key = element_rect(fill = NA),
        legend.background = element_blank(),
        legend.title = element_text(colour = "#000000", size = 12, face = "bold"),
        legend.text = element_text(colour = "#000000", size = 10, face="bold"))+
  guides(fill = guide_legend(TeX(r'(Hedrick & Lacy's relatedness)'), title.theme = element_text(size = 7, face = "bold"),
                             label.theme = element_text(size = 7)))

# Saves plot ~
ggsave(rxy_stacked_Plot, file = "~/Desktop/Scripts/EUostrea/Figures/Relatedness/Rxy_stacked_Plot_allpop_27feb23.pdf", device = cairo_pdf, scale = 1.1, width = 8, height = 10, dpi = 600)
dev.off()

# Sibship in the WADD population
str(WADD)
head(WADD)
WADD %>% write.table(., "~/Desktop/Scripts/Data/Relatedness_EUostrea/NGSrelate_output_WADD_2mar23.tsv", 
                     sep="\t", row.names=FALSE)
# Filter out rows where a and b are the same
WADD_filtered <- WADD %>% filter(a != b)
# Create scatterplot
# Split ab variable into a and b variables
WADD$a <- as.numeric(gsub("^(\\d+)_.*", "\\1", WADD$ab))
WADD$b <- as.numeric(gsub("^\\d+_(\\d+)$", "\\1", WADD$ab))

# Sort a and b variables in ascending order
WADD$ab_sorted <- apply(WADD[, c("a", "b")], 1, function(x) paste(sort(x), collapse = "_"))

# Remove rows with duplicate ab_sorted values
#WADD_unique <- unique(WADD, by = "ab_sorted")
WADD$a <- ifelse(WADD$a == 0, 1, WADD$a + 1)
WADD$b <- ifelse(WADD$b == 0, 1, WADD$b + 1)

# Create a subset of WADD with only the highest rab values
WADD_high_rab <- WADD[order(WADD$rab, decreasing = TRUE)[1:10],]

# Create a plot with dots and labels
ggplot(WADD, aes(x = nSites, y = rab, color=rab)) +
  geom_point(size = 3) +
  scale_color_gradient(low = "blue", high = "red", name = "") +
  geom_text(data = WADD_high_rab, aes(label = ab_sorted), size = 3, vjust = -0.5, color = "black") +
  theme(panel.background = element_rect(fill = "#ffffff"),
        panel.border = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = margin(t = 0.005, b = 0.005, r = .4, l = .4, unit = "cm"),
        axis.line = element_blank(),
        axis.text = element_text(colour = "#000000", size = 10, face = "bold"),
        axis.text.x = element_text( size = 10, face = "bold", angle = 45, vjust = 1, hjust = 1),
        legend.position = "right",
        legend.key = element_rect(fill = NA),
        legend.background = element_blank(),
        legend.title = element_text(colour = "#000000", size = 10, face = "bold"),
        legend.text = element_text(colour = "#000000", size = 10, face="bold"))
str(WADD$ab)


#
##
### The END ~~~~~


