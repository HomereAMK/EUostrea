# install and load necessary packages
# Loads required packages ~
setwd(dir = "~/Desktop/Scripts/Map_ala_N.Lou/00_scripts/")

# Cleans the environment ~ 
rm(list=ls())
pacman::p_load(devtools, tidyverse, viridis,ggrepel, knitr, rgdal, marmap, ggmap, sdmpredictors, ggcorrplot, raster, gdistance, ade4, cowplot, ggplot2, maps, scatterpie)
# install the rgdal, raster, and gdistance packages
# install.packages("rgdal")
# install.packages("raster")
# install.packages("gdistance")

#Load the populations
pop_info <- read_csv("~/Desktop/Scripts/Shucking/01_infofiles/Master_list_pop.csv")%>%
  mutate(lat=lat, long=lon) %>%
  mutate(long=round(lon, 2), lat=round(lat,2)) %>%
  filter(Tag != "HVAD")%>%
  filter(Tag != "NELL")%>%
  filter(Tag != "FURI")%>%
  filter(Tag != "KLEV")%>%
  filter(Tag != "SYDK")%>%
  filter(Tag != "VADE")%>%
  filter(Tag != "GASO")%>%
  filter(Tag != "HFJO")%>%
  filter(Tag != "LILL")%>%
  filter(Tag != "BOVA")%>%
  filter(Tag != "ORNE")%>%
  filter(Tag != "Havs")%>%
  filter(Tag != "HAVS")%>%
  filter(Tag != "SVAL")%>%
  filter(Tag != "RAMS")%>%
  filter(Tag != "USAM")

# Generate random genotype frequencies but generate random genotype frequencies while ensuring that the total frequency of all genotypes adds up to 1 for each population
# populations <- populations %>%
#   group_by(Tag) %>%
#   mutate(
#     frequency = runif(1),
#     genotype_1_frequency = frequency,
#     genotype_2_frequency = 1 - frequency,
#     genotype_3_frequency = 0
#   ) %>%
#   ungroup() %>%
#   dplyr::select(-frequency)


# Load the full haplotypes for the Inversions Regions dataset
prop <- read_tsv("~/Desktop/Scripts/Data/InvReg_EUostrea/PropHomo1_het_homo2_InvReg_27feb23.tsv")

# Merge pop_info and prop dataset
head(pop_info)
head(prop)
prop <- prop %>% dplyr::rename( Tag =population )
info_prop <- merge(pop_info, prop, by = "Tag")

# Plot for INvReg Reg04
# InvRegPieplot_sim <- readOGR(dsn = "ne_10m_admin_0_countries.shp", layer = "ne_10m_admin_0_countries") %>%
#   fortify() %>%
#   filter(lat>35, lat<73) %>%
#   ggplot(aes(x =long, y = lat)) +
#   geom_polygon(aes(group = group), fill="grey90", color="black", size=0.1)+
#   coord_map(projection = "azequidistant", xlim = c(-15, 20), ylim = c(36, 71))+
#   geom_scatterpie(data=info_prop %>%
#                     filter(lg == "Reg04"), 
#                     aes(x=long, y=lat, group=Locality),
#                     cols=c("prop_homo_1","prop_homo_2", "prop_het"), color="black", fill=(values = c("#5a0c30","#973fc9","#91ccff")))+
#   theme_bw() 
# head(info_prop)

# Reorders Population ~
info_prop$Tag <- factor(info_prop$Tag, ordered = T,
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

# Plot the pie charts corresponding to InvReg04 haplotype frequencies EUROPE
map_pie_Reg05 <-  readOGR(dsn = "ne_10m_admin_0_countries.shp", layer = "ne_10m_admin_0_countries") %>%
  fortify() %>%
  filter(lat>35, lat<73) %>%ggplot(aes(x =long, y = lat)) +
    geom_polygon(aes(group = group), fill="grey90", color="black", size=0.1)+
    coord_map(projection = "azequidistant", xlim = c(-13, 19), ylim = c(36, 71))+
    geom_scatterpie(data=info_prop %>%
                      filter(lg == "Reg05"), 
                    aes(x=long, y=lat, group=Locality),
                    cols=c("prop_homo_1","prop_homo_2", "prop_het"), color="black", alpha=.7, size=2, pie_scale=1.2)+
  scale_fill_brewer(palette = "Accent", 
                    labels = c("Homo1", "Homo2", "Het"), 
                    name = "Haplotype") +
  geom_label_repel(data=info_prop %>%
                       filter(!(Country %in% c("Norway", "Sweden", "Denmark"))) %>%
                       filter(lg == "Reg08"), aes(label=Tag), size=3, 
                     seed = 10, min.segment.length = 0, force = 0, segment.curvature = 1, segment.angle = 3,
                     nudge_x = 1, nudge_y = 0.8, max.overlaps = 2) +
    labs(title = "Haplotypes distribution and frequencies of the SVs at scaffold05",
         x = "Longitude",
         y = "Latitude") +
    theme(panel.spacing = unit(0.1, "lines"),
          axis.title.x=element_text(),
          legend.position="right",
          text = element_text(size=10),
          axis.text = element_text(size=6))+
    theme(legend.key = element_blank()) +
    theme(legend.title=element_blank()) +
    theme(axis.title.x = element_text(size = 10, color="#000000", face="bold", margin = margin(t = 20, r = 0, b = 0, l = 0)),
          axis.title.y = element_text(size = 10, color="#000000", face="bold", margin = margin(t = 0, r = 20, b = 0, l = 0)),
          title = element_text(size = 10, color="#000000", face="bold")) +
    theme(legend.text=element_text(size=11)) +
    theme(panel.background = element_rect(fill = '#FAFAFA')) +
    theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank()) +
    theme(axis.line = element_line(colour = "#000000", size = 0.3)) +
    theme(panel.border = element_blank()) +
    guides(fill = guide_legend(title = "Haplotype", title.theme = element_text(size = 15, face = "bold"),
                             label.theme = element_text(size = 14)))

map_pie_Reg04
ggsave(map_pie_Reg05, file = "~/Desktop/Scripts/EUostrea/Figures/InvReg/scaffold5/28mar23MapPie_InvReg05_populations.pdf", device = cairo_pdf, scale = 1.1, width = 25, height = 24, dpi = 300)
  
  

# Plot the pie charts corresponding to InvReg05 haplotype frequencies SCANDINAVIA
scan_reg08 <-  readOGR(dsn = "ne_10m_admin_0_countries.shp", layer = "ne_10m_admin_0_countries") %>%
    fortify() %>%
    filter(lat>35, lat<73) %>%ggplot(aes(x =long, y = lat)) +
    geom_polygon(aes(group = group), fill="grey90", color="black", size=0.1)+
    coord_map(projection = "azequidistant", xlim = c(0, 17), ylim = c(50, 65))+
  geom_scatterpie(data = info_prop %>%
                    filter(Country %in% c("Norway", "Netherlands", "Sweden", "Denmark")) %>%
                    filter(lg == "Reg08"),
                  aes(x = long, y = lat, group = Locality),
                  cols = c("prop_homo_1", "prop_homo_2", "prop_het"), color = "black", alpha = .7, size = 2, pie_scale = 2) +
  scale_fill_brewer(palette = "Accent", 
                    labels = c("Homo1", "Homo2", "Het"), 
                    name = "Haplotype") +
    geom_label_repel(data=info_prop %>%
                       filter(Country %in% c("Norway", "Netherlands", "Sweden", "Denmark")) %>%
                       filter(lg == "Reg04"), aes(label=Tag), size=3, 
                     seed = 10, min.segment.length = 0, force = 0, segment.curvature = 0, segment.angle = 0,
                     nudge_x = -1, nudge_y = 0.2, max.overlaps = 1) +
    #scale_colour_manual(values =c("#91BD96", "#91BD96",
     #                              "#02630C","#02630C","#02630C", "#02630C", "#02630C",
        #                           "#45D1F7", "#45D1F7",
        #                           "#588cad", "#588cad", "#588cad", "#588cad", "#588cad",
              #                     "#240377", "#240377", "#240377", "#240377" ))+
    labs(title = "Haplotypes distribution and frequencies of the SVs at scaffold08",
         x = "Longitude",
         y = "Latitude"
         ) +
    theme(panel.spacing = unit(0.1, "lines"),
          axis.title.x=element_text(),
          legend.position="right",
          text = element_text(size=10),
          axis.text = element_text(size=6))+
    theme(legend.key = element_blank()) +
    theme(legend.title=element_blank()) +
    theme(axis.title.x = element_text(size = 10, color="#000000", face="bold", margin = margin(t = 20, r = 0, b = 0, l = 0)),
          axis.title.y = element_text(size = 10, color="#000000", face="bold", margin = margin(t = 0, r = 20, b = 0, l = 0)),
          title = element_text(size = 10, color="#000000", face="bold")) +
    theme(legend.text=element_text(size=11)) +
    theme(panel.background = element_rect(fill = '#FAFAFA')) +
    theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank()) +
    theme(axis.line = element_line(colour = "#000000", size = 0.3)) +
    theme(panel.border = element_blank()) +
    guides(fill = guide_legend(title = "Haplotype", title.theme = element_text(size = 15, face = "bold"),
                             label.theme = element_text(size = 14)))
scan_reg08
ggsave(scan_reg08 , file = "~/Desktop/Scripts/EUostrea/Figures/InvReg/scaffold8/28mar23--SCANDMapPie_InvReg08_populations.pdf", device = cairo_pdf, scale = 1.1, width = 25, height = 24, dpi = 300)




## HWE test











# Gen_freq <- info_prop %>% dplyr::select(Tag, lg, prop_homo_1, prop_het, prop_homo_2)
# head(Gen_freq)
# # Function to perform the Hardy-Weinberg equilibrium test
# hardy_weinberg_test <- function(df) {
# 
#   # Calculate allele frequencies
#   p <- prop_homo_1 * 0.5 + prop_het * 0.5
#   q <- prop_homo_2 * 0.5 + prop_het * 0.5
#   
#   # Calculate expected genotype frequencies under Hardy-Weinberg equilibrium
#   exp_homo1 <- p^2
#   exp_het <- 2 * p * q
#   exp_homo2 <- q^2
#   
#   # Calculate the chi-squared statistic
#   observed <- c(prop_homo_1, prop_het, prop_homo_2)
#   expected <- c(exp_homo_1, exp_het, exp_homo_2)
#   chi_sq <- sum((observed - expected)^2 / expected)
#   
#   # Perform the chi-squared test
#   df <- length(observed) - 2 - 1
#   p_value <- pchisq(chi_sq, df, lower.tail = FALSE)
#   
#   return(p_value)
# }
# 
# 
# hardy_weinberg_test <- function(df, homo1_col, het_col, homo2_col) {
#   
#   # Add a small constant to genotype frequencies
#   eps <- 0.00000001
#   df[[homo1_col]] <- df[[homo1_col]] + eps
#   df[[het_col]] <- df[[het_col]] + eps
#   df[[homo2_col]] <- df[[homo2_col]] + eps
#   
#   # Calculate allele frequencies
#   p <- df[[homo1_col]] * 0.5 + df[[het_col]] * 0.5
#   q <- df[[homo2_col]] * 0.5 + df[[het_col]] * 0.5
#   
#   # Calculate expected genotype frequencies under Hardy-Weinberg equilibrium
#   exp_homo1 <- p^2
#   exp_het <- 2 * p * q
#   exp_homo2 <- q^2
#   
#   # Calculate the chi-squared statistic
#   observed <- c(df[[homo1_col]], df[[het_col]], df[[homo2_col]])
#   expected <- c(exp_homo1, exp_het, exp_homo2)
#   chi_sq <- sum((observed - expected)^2 / expected)
#   
#   # Perform the chi-squared test
#   df <- length(observed) - 2 - 1
#   p_value <- chisq.test(observed, p = expected)
#   
#   return(p_value)
# }
# 
# 
# 
# result <- Gen_freq %>%
#   group_by(Tag, lg) %>%
#   mutate(p_value = hardy_weinberg_test(., "prop_homo_1", "prop_het", "prop_homo_2"))
# result
# 
# 
# 
# 
# 
# # Apply the function to each group (Tag, lg) and store the p-values
# result <- Gen_freq %>%
#   group_by(Tag, lg) %>%
#   dplyr::summarise(p_value = hardy_weinberg_test(cur_data()))
# 
# # Print the result
# print(result)
# 
# # Interpret the results
# alpha <- 0.05
# result <- result %>%
#   mutate(significance = ifelse(p_value < alpha, "Significant deviation", "No significant deviation"))
# 
# # Print the interpreted results
# print(result)

#
##
### The END ~~~~~
