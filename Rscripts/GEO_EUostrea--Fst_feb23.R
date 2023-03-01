### The BEGINNING ~~~~
##
# ~ Plots BSG_Turbot--Fst | By George Pacheco.


# Cleans the environment ~ 
rm(list=ls())




# Loads required packages ~
pacman::p_load(pheatmap, tidyverse, reshape2)


# Loads Fst table ~
Fst <- read.table("~/Desktop/Scripts/Data/Fst_EUostrea/Fst_matrix/10feb23--mindInd0.25_globalList_Unfolded_Fst.tsv", sep = "\t", header = FALSE, stringsAsFactors = FALSE)


# Adds column names ~
colnames(Fst) <- c("Pop1", "Pop2", "NumberOfSites", "Unweighted", "Weighted")


# Melts datasets ~
Fst_Pops = union(Fst$Pop1, Fst$Pop2)
n = length(Fst_Pops)


# Creates Fst-Sites matrix ~
FstSites <- matrix(0, nrow = n, ncol = n, dimnames = list(Fst_Pops, Fst_Pops))
for (i in 1:nrow(Fst)) {
  FstSites[Fst[i, "Pop1"], Fst[i, "Pop2"]] = Fst[i, "NumberOfSites"]
  FstSites[Fst[i, "Pop2"], Fst[i, "Pop1"]] = Fst[i, "Weighted"]}


# Writes Fst-Sites matrix ~
write.table(FstSites, "~/Desktop/Scripts/Data/Fst_EUostrea/Fst_matrix/10feb23--mindInd0.25_globalList_Unfolded_Fst-Sites.txt", sep = "\t", row.names = TRUE, col.names = TRUE)


# Creates Fst-Fst matrix ~
Fst_Fst <- matrix(0, nrow = n, ncol = n, dimnames = list(Fst_Pops, Fst_Pops))
for (i in 1:nrow(Fst)) {
  Fst_Fst[Fst[i, "Pop1"], Fst[i, "Pop2"]] = Fst[i, "Weighted"]
  Fst_Fst[Fst[i, "Pop2"], Fst[i, "Pop1"]] = Fst[i, "Weighted"]}


# Writes Fst-Fst matrix ~
write.table(Fst_Fst, "~/Desktop/Scripts/Data/Fst_EUostrea/Fst_matrix/10feb23--mindInd0.25_globalList_Unfolded_Fst-Fst.txt", sep = "\t", row.names = TRUE, col.names = TRUE)


# Gets Fst midpoint ~
Pops <- which(upper.tri(Fst_Fst), arr.ind = TRUE)
Fst_df <- data.frame(Site1 = dimnames(Fst_Fst)[[2]][Pops[, 2]],
                     Site2 = dimnames(Fst_Fst)[[1]][Pops[, 1]],
                     Value = Fst_Fst[Pops] %>% round(digits = 6))
Fstmiddle = max(Fst_df$Value) / 2


# Corrects x-axis labels ~
#levels(Fst$Pop1 <- sub("NorthSea", "North Sea", Fst$Pop1))
#levels(Fst$Pop2 <- sub("BalticSea", "Baltic Sea", Fst$Pop2))


# Gets Fst label ~
Fst.label = expression(italic("F")[ST])

# Reorder the pop ~
Fst$Pop2 <- factor(Fst$Pop2, ordered = T,
                   levels = c("MOLU", "ZECE", "CRES","ORIS","CORS", "PONT",  "RIAE", "MORL",
                              "USAM", "TOLL", "COLN", "BARR","TRAL", "CLEW", "RYAN",
                              "GREV", "WADD", "NISS","LOGS","VENO", "HALS", "THIS",
                              "KALV", "HYPP","LANG", "BUNN", "DOLV", "HAUG", "HAFR",
                              "INNE","VAGS", "AGAB", "OSTR"))

# Creates plot ~
plotFstmat <- ggplot(data = Fst, aes(x = Pop1, y = Pop2, fill = Weighted)) +
  geom_tile(color = "#ffffff", lwd = 1.5, linetype = 1, width = 1, height = 1) +
  coord_fixed() +
  geom_text(aes(label = round(Weighted,digits =3)), color = "#ffffff", size = 2) +
  scale_fill_gradientn(colors = c("#fde0dd", "#fa9fb5", "#c51b8a"), 
                       name = "Weighted Fst", 
                       breaks = c(0,  0.05, 0.10,  0.15, 0.2, 0.25),
                       labels = c("0", "0.05", "0.10", "0.15", "0.2", "0.25"),
                       limits = c(0, .25)) +
  scale_x_discrete(expand = c(0,0))+
  scale_y_discrete(expand = c(0,0), position = "right")+
  labs(x = "Population 1", y = "Population 2", 
       title = "Weighted Fst Values between Populations", 
       fill = "Weighted Fst") +
  theme(panel.background = element_rect(fill = "#ffffff"),
        panel.border = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = margin(t = 0.005, b = 0.005, r = .4, l = .4, unit = "cm"),
        axis.line = element_blank(),
        axis.text = element_text(colour = "#000000", size = 10, face = "bold"),
        axis.text.x = element_text(colour = c("#A02353", "#A02353", "#A02353", "#AD5B35", "#ad7358",
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
                                              "#240377", "#240377", "#240377", "#240377"), size = 10, face = "bold", angle = 45, vjust = 1, hjust = 1),
        axis.text.y = element_text(colour = c( "#A02353", "#A02353", "#AD5B35", "#ad7358",
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
        axis.ticks = element_line(color = "#000000", size = .3),
        legend.position = "right",
        legend.key = element_rect(fill = NA),
        legend.background = element_blank(),
        legend.title = element_text(colour = "#000000", size = 10, face = "bold"),
        legend.text = element_text(colour = "#000000", size = 10, face="bold"))


# Saves plot ~
ggsave(plotFstmat, file = "~/Desktop/Scripts/EUostrea/Figures/Fst/FstMatrix/FstValuePlot_23feb23.pdf",
       device = cairo_pdf, scale = 1, width = 12, height = 12, dpi = 600)
dev.off()
  

#
##
### The END ~~~~~
### The END ~~~~~