# install and load necessary packages
# Loads required packages ~
setwd(dir = "~/Desktop/Scripts/Map_ala_N.Lou/00_scripts/")

# Cleans the environment ~ 
rm(list=ls())
library(ggplot2)
library(ggmap)
library(ggrepel)
pacman::p_load(devtools, tidyverse, ggrepel, knitr, rgdal, marmap, sdmpredictors, ggcorrplot, raster, gdistance, ade4, cowplot, ggplot2, maps)
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

# Generate random genotype frequencies for each population
populations <- pop_info
populations <- populations %>%
  group_by(Locality) %>%
  mutate(
    genotype_1_frequency = runif(1),
    genotype_2_frequency = runif(1),
    genotype_3_frequency = runif(1)
  )


readOGR(dsn = "ne_10m_admin_0_countries.shp",
              layer = "ne_10m_admin_0_countries") %>%
  fortify() %>%
  filter(lat>35, lat<73) %>%
  ggplot(aes(x =long, y = lat)) +
  geom_polygon(aes(group = group), fill="grey90", color="black", size=0.1)+
  coord_map(projection = "azequidistant", xlim = c(-15, 20), ylim = c(36, 71))+
  geom_scatterpie(aes(x=long, y=lat, group=Locality), data=populations,
                  cols=c("genotype_1_frequency","genotype_2_frequency", "genotype_3_frequency"), color=NA)+
  theme_bw() +
  theme(axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank())+
  guides(fill = guide_legend(title = "Populations", title.theme = element_text(size = 10.5, face = "bold", family = "Helvetica"),
                             label.theme = element_text(size = 8, face = "italic", family = "Helvetica"),
                             override.aes = list(size = 3, alpha = 0.9)), colour = "none")

