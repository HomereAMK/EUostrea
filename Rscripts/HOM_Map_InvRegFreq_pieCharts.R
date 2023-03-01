# install and load necessary packages
# Loads required packages ~
setwd(dir = "~/Desktop/Scripts/Map_ala_N.Lou/00_scripts/")

# Cleans the environment ~ 
rm(list=ls())
pacman::p_load(devtools, tidyverse, ggrepel, knitr, rgdal, marmap, ggmap, sdmpredictors, ggcorrplot, raster, gdistance, ade4, cowplot, ggplot2, maps, scatterpie)
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
info_prop <- merge(pop_info, prop, by = "Tag")

# Plot for INvReg Reg04
InvRegPieplot_sim <- readOGR(dsn = "ne_10m_admin_0_countries.shp", layer = "ne_10m_admin_0_countries") %>%
  fortify() %>%
  filter(lat>35, lat<73) %>%
  ggplot(aes(x =long, y = lat)) +
  geom_polygon(aes(group = group), fill="grey90", color="black", size=0.1)+
  coord_map(projection = "azequidistant", xlim = c(-15, 20), ylim = c(36, 71))+
  geom_scatterpie(data=info_prop %>%
                    filter(lg == "Reg04"), 
                    aes(x=long, y=lat, group=Locality),
                    cols=c("prop_homo_1","prop_homo_2", "prop_het"), color="black", fill=(values = c("#5a0c30","#973fc9","#91ccff")))+
  theme_bw() 
head(info_prop)

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
  readOGR(dsn = "ne_10m_admin_0_countries.shp", layer = "ne_10m_admin_0_countries") %>%
  fortify() %>%
  filter(lat>35, lat<73) %>%ggplot(aes(x =long, y = lat)) +
    geom_polygon(aes(group = group), fill="grey90", color="black", size=0.1)+
    coord_map(projection = "azequidistant", xlim = c(-13, 19), ylim = c(36, 71))+
    geom_scatterpie(data=info_prop %>%
                      filter(lg == "Reg08"), 
                    aes(x=long, y=lat, group=Locality),
                    cols=c("prop_homo_1","prop_homo_2", "prop_het"), color="black", alpha=.7, size=2, pie_scale=1.2)+
    scale_fill_manual(values=c("#5a0c30","#91ccff","#973fc9"))+
    geom_label_repel(data=info_prop %>%
                       filter(!(Country %in% c("Norway", "Sweden", "Denmark"))) %>%
                       filter(lg == "Reg08"), aes(color=Tag,label=Tag), size=3, 
                     seed = 10, min.segment.length = 10, force = 1, segment.curvature = 1, segment.angle = 3,
                     nudge_x = 1.5, nudge_y = 0.8, max.overlaps = Inf) +
    scale_colour_manual(values =c( "#A02353", "#A02353", "#A02353",
                                   "#AD5B35",
                                   "#ad7358",
                                   "#CC480C",  "#CC480C",
                                   "#969696",
                                   "#D38C89", "#D38C89", "#D38C89",
                                   "#C89AD1", "#C89AD1",
                                   "#7210A0",
                                   "#91BD96", "#91BD96",
                                   "#02630C","#02630C","#02630C", "#02630C", "#02630C",
                                   "#45D1F7", "#45D1F7",
                                   "#588cad", "#588cad", "#588cad", "#588cad", "#588cad",
                                   "#240377", "#240377", "#240377", "#240377" ))+
    labs(title = "Haplotypes distribution and frequencies of the putative inversion region at scaffold08",
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
    guides(colour = FALSE)
  
  
ggsave(, file = "~/Desktop/Scripts/EUostrea/Figures/InvReg/HetatInvReg_InvReg040508_histo.pdf", device = cairo_pdf, scale = 1.1, width = 12, height = 8, dpi = 300)
dev.off()
  
# Plot the pie charts corresponding to InvReg04 haplotype frequencies SCANDINAVIA
  readOGR(dsn = "ne_10m_admin_0_countries.shp", layer = "ne_10m_admin_0_countries") %>%
    fortify() %>%
    filter(lat>35, lat<73) %>%ggplot(aes(x =long, y = lat)) +
    geom_polygon(aes(group = group), fill="grey90", color="black", size=0.1)+
    coord_map(projection = "azequidistant", xlim = c(0, 17), ylim = c(50, 65))+
    geom_scatterpie(data=info_prop %>%
                      filter(Country %in% c("Norway", "Netherlands", "Sweden", "Denmark")) %>%
                      filter(lg == "Reg04"),
                      aes(x=long, y=lat, group=Locality),
                    cols=c("prop_homo_1","prop_homo_2", "prop_het"), color="black", alpha=.7, size=2, pie_scale=2.2)+
    scale_fill_manual(values=c("#5a0c30","#91ccff","#973fc9"))+
    geom_label_repel(data=info_prop %>%
                       filter(Country %in% c("Norway", "Netherlands", "Sweden", "Denmark")) %>%
                       filter(lg == "Reg04"), aes(color=Tag,label=Tag), size=2, 
                     seed = 10, min.segment.length = 0, force = 0, segment.curvature = 0, segment.angle = 0,
                     nudge_x = 0.6, nudge_y = 0.2, max.overlaps = Inf) +
    scale_colour_manual(values =c("#91BD96", "#91BD96",
                                   "#02630C","#02630C","#02630C", "#02630C", "#02630C",
                                   "#45D1F7", "#45D1F7",
                                   "#588cad", "#588cad", "#588cad", "#588cad", "#588cad",
                                   "#240377", "#240377", "#240377", "#240377" ))+
    labs(title = "Haplotypes distribution and frequencies of the putative inversion region at scaffold04",
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
    guides(colour = FALSE)
  
  

#
##
### The END ~~~~~
