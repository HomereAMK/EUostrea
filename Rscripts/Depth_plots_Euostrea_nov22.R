### The BEGINNING ~~~~~
##
# ~ Plots depth for Low-cov libraries | Written by  Homère J. Alves Monteiro

# Cleans the environment ~ 
rm(list=ls())

# Sets working directory ~
setwd("~/Desktop/Scripts/")

# Loads required packages ~
pacman::p_load( tidyverse,cowplot,knitr,RColorBrewer,MetBrewer)

# Load the dataset ~
df<-read_csv("EUostrea/02_data/Depth/tmp_EUostreadepth_nov22.csv")

# Create a Pop column
pop<-substr(df$bamfile,0,4) #name only 4characters
df <- cbind(df, pop)

# Plot
tmp_MeanDepth <-ggplot(df, aes(x=pop, y=mean_depth, color=pop)) +
  geom_jitter(aes(color=pop), height=0) +
  geom_boxplot(alpha=0.5, outlier.colour=NA) +
  #scale_colour_gradient(colors = met.brewer("")) +
  coord_flip() +
  theme_cowplot()+
  theme(
    panel.background = element_rect(fill = "#ffffff"),
    panel.grid.minor.x = element_blank(),
    panel.grid.major = element_blank(),
    panel.spacing = unit(.2, "lines"),
    plot.title = element_blank(),
    strip.background = element_rect(colour = "#000000", fill = "#FAFAFA", size = .05),
    strip.text.x = element_text(colour = "#000000", face = "bold", size = 6, margin = margin(.75, 0, .75, 0, "cm")),
    strip.text.y = element_text(colour = "#000000", face = "bold", size = 9, angle = 90, margin = margin(0, .1, 0, .1, "cm")),
    legend.position = "none")
last_plot()
# Saves the final plot ~
ggsave(tmp_MeanDepth, file = "~/Desktop/Scripts/EUostrea/Figures/Depth/tmp_MeanDepth_nov22.pdf", device = cairo_pdf, width = 16, height = 8, dpi = 300)

tmp_Gcov <-ggplot(df, aes(x=pop, y=proportion_of_reference_covered , color=pop)) +
  geom_jitter(aes(color=pop), height=0) +
  geom_boxplot(alpha=0.5, outlier.colour=NA) +
  #scale_colour_gradient(colors = met.brewer("")) +
  coord_flip() +
  theme_cowplot()+
  theme(
    panel.background = element_rect(fill = "#ffffff"),
    panel.grid.minor.x = element_blank(),
    panel.grid.major = element_blank(),
    panel.spacing = unit(.2, "lines"),
    plot.title = element_blank(),
    strip.background = element_rect(colour = "#000000", fill = "#FAFAFA", size = .05),
    strip.text.x = element_text(colour = "#000000", face = "bold", size = 6, margin = margin(.75, 0, .75, 0, "cm")),
    strip.text.y = element_text(colour = "#000000", face = "bold", size = 9, angle = 90, margin = margin(0, .1, 0, .1, "cm")),
    legend.position = "none")
last_plot()
# Saves the final plot ~
ggsave(tmp_MeanDepth, file = "~/Desktop/Scripts/EUostrea/Figures/Depth/tmp_Gcov_nov22.pdf", device = cairo_pdf, width = 16, height = 8, dpi = 300)
