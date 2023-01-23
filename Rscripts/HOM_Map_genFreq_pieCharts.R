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


# read in the shapefile
readOGR(dsn = "ne_10m_admin_0_countries.shp", layer = "ne_10m_admin_0_countries") %>%
# convert the shapefile to a dataframe
fortify() %>%
# filter the dataframe to only include latitudes between 35 and 73
filter(lat>35, lat<73) %>%
# create a ggplot object with the x axis as longitude and the y axis as latitude
ggplot(aes(x =long, y = lat)) +
# add a polygon layer to the ggplot object, filling the polygons with grey90, and colouring the edges black
geom_polygon(aes(group = group), fill="grey90", color="black", size=0.1)+
# set the projection to azimuthal equidistant, and set the x and y limits
coord_map(projection = "azequidistant", xlim = c(-15, 20), ylim = c(36, 71))+
# add a scatterpie layer to the ggplot object, using the populations dataframe, and colouring the pies according to the genotype frequencies
geom_scatterpie(aes(x=long, y=lat, group=Locality), data=populations, cols=c("genotype_1_frequency","genotype_2_frequency", "genotype_3_frequency"), color=NA)+
# set the theme to black and white
theme_bw() +
# remove the axis titles, ticks and text
theme(axis.title = element_blank(), axis.ticks = element_blank(), axis.text = element_blank())+
# add a legend, with the title "Populations"
guides(fill = guide_legend(title = "Populations", title.theme = element_text(size = 10.5, face = "bold", family = "Helvetica"), label.theme = element_text(size = 8, face = "italic", family = "Helvetica"), override.aes = list(size = 3, alpha = 0.9)), colour = "none")

