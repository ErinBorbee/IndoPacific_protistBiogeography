#load packages
library(vegan)
library(ggplot2)
library(ggpubr)

#set working directory
setwd("~/Desktop/JPhyc_sub/")

#import data
asv_table <- read.csv("~/Desktop/DataAnalysis/threshold95/SAR/SAR-percent-table-95.csv", row.names = 1)
asv_table <- asv_table[,-c(1:8)]
asv_transformed <- as.data.frame(t(asv_table))
otu_table <- read.csv("~/Desktop/DataAnalysis/threshold95/SAR/SAR_OTU_percent.csv", row.names = 1)
otu_transformed <- as.data.frame(t(otu_table))
otu_filtered <- subset(otu_transformed, rownames(otu_transformed) %in% rownames(metadata))

metadata <- read.csv("~/Desktop/DataAnalysis/Master_V9_MappingFile.csv", row.names = 1)
names(metadata)

numeric_metadata <- metadata[,c(12,17,21:29,32:33,35:37,39:42)]

metadata_filtered <- as.data.frame(subset(numeric_metadata, rownames(numeric_metadata) %in% rownames(asv_transformed)))

#run dbRDA
dbRDA <- dbrda(asv_transformed ~ npp_mean + wave_mean + avg_complexity + cover + detritAbund + carniAbund + 
                 herbAbund + invertAbund + omniAbund + planktAbund + coraliAbund + hardCoral_percent + 
                 availSubstrate_percent + rubble_percent + algae_percent + softCoral_percent + 
                 abundance + target_biomass + nontarget_biomass + total_biomass, metadata_filtered, 
               distance = "bray", na.action = na.exclude)
summary(dbRDA)
anova(dbRDA, permutations = 999)
anova(dbRDA, by = "term", permutations = 999)
anova(dbRDA, by = "axis", permutations = 999)

dbRDA <- dbrda(asv_transformed ~ planktAbund, 
               metadata_filtered, distance = "bray", na.action = na.exclude)
summary(dbRDA)
anova(dbRDA, permutations = 999)
anova(dbRDA, by = "term", permutations = 999)
anova(dbRDA, by = "axis", permutations = 999)

#subset metadata and asvs by depth and filter size
metadata_M04 <- subset(metadata, subset == "M04")
metadata_M12 <- subset(metadata, subset == "M12")
metadata_S04 <- subset(metadata, subset == "S04")
metadata_S12 <- subset(metadata, subset == "S12")

asvs_M04 <- as.data.frame(subset(asv_transformed, rownames(asv_transformed) %in% rownames(metadata_M04)))
asvs_M12 <- as.data.frame(subset(asv_transformed, rownames(asv_transformed) %in% rownames(metadata_M12)))
asvs_S04 <- as.data.frame(subset(asv_transformed, rownames(asv_transformed) %in% rownames(metadata_S04)))
asvs_S12 <- as.data.frame(subset(asv_transformed, rownames(asv_transformed) %in% rownames(metadata_S12)))

metadata_S04 <- subset(metadata_S04, rownames(metadata_S04) %in% rownames(asvs_S04))

#Run dbRDA by subset
dbRDA_M04 <- dbrda(asvs_M04 ~ wave_mean + cover + rubble_percent, 
                   metadata_M04, distance = "bray", na.action = na.exclude)
summary(dbRDA_M04)
anova(dbRDA_M04, permutations = 999)
anova(dbRDA_M04, by = "term", permutations = 999)
anova(dbRDA_M04, by = "axis", permutations = 999)
vif.cca(dbRDA_M04)

dbRDA_M12 <- dbrda(asvs_M12 ~ npp_mean + detritAbund, 
                   metadata_M12, distance = "bray", na.action = na.exclude)
summary(dbRDA_M12)
anova(dbRDA_M12, permutations = 999)
anova(dbRDA_M12, by = "term", permutations = 999)
anova(dbRDA_M12, by = "axis", permutations = 999)
vif.cca(dbRDA_M12)

dbRDA_S04 <- dbrda(asvs_S04 ~ wave_mean + detritAbund + coraliAbund, 
                   metadata_S04, distance = "bray", na.action = na.exclude)
summary(dbRDA_S04)
anova(dbRDA_S04, permutations = 999)
anova(dbRDA_S04, by = "term", permutations = 999)
anova(dbRDA_S04, by = "axis", permutations = 999)
vif.cca(dbRDA_S04)

dbRDA_S12 <- dbrda(asvs_S12 ~ npp_mean + cover + invertAbund + 
                 rubble_percent, metadata_S12, 
               distance = "bray", na.action = na.exclude)
summary(dbRDA_S12)
anova(dbRDA_S12, permutations = 999)
anova(dbRDA_S12, by = "term", permutations = 999)
anova(dbRDA_S12, by = "axis", permutations = 999)
vif.cca(dbRDA_S12)

#extract coordinates for each subset for plotting
M04_coords <- as.data.frame(scores(dbRDA_M04, display = "sites"))
M12_coords <- as.data.frame(scores(dbRDA_M12, display = "sites"))
S04_coords <- as.data.frame(scores(dbRDA_S04, display = "sites"))
S12_coords <- as.data.frame(scores(dbRDA_S12, display = "sites"))

M04_arrows <- as.data.frame(scores(dbRDA_M04, display = "bp"))
M12_arrows <- as.data.frame(scores(dbRDA_M12, display = "bp"))
S04_arrows <- as.data.frame(scores(dbRDA_S04, display = "bp"))
S12_arrows <- as.data.frame(scores(dbRDA_S12, display = "bp"))

#merge coordinates with metadata
M04_merged <- merge(M04_coords, metadata_M04, by = "row.names")
M12_merged <- merge(M12_coords, metadata_M12, by = "row.names")
S04_merged <- merge(S04_coords, metadata_S04, by = "row.names")
S12_merged <- merge(S12_coords, metadata_S12, by = "row.names")

#plot dbRDAs
M04_plot <- ggplot() + 
  geom_point(data = M04_merged, aes(dbRDA1, dbRDA2, fill = region), pch = 21, size = 3, alpha = 0.8) + 
  stat_ellipse(data = M04_merged, aes(dbRDA1, dbRDA2, color = region)) +
  geom_segment(data = M04_arrows, aes(x = 0, y = 0, xend = dbRDA1, yend = dbRDA2), arrow = arrow(length = unit(0.5, "cm"))) +
  geom_label(data = M04_arrows, aes(x = dbRDA1, y = dbRDA2, label = rownames(M04_arrows)), fill = "white", color = "black") +
  theme_bw() + theme(panel.grid = element_blank()) + 
  scale_fill_manual(values = c('Lombok' = "#8e86c5",
                               'Wakatobi' = "#007cb5",
                               'Misool' = "#a6c54b",
                               'Waigeo' = "#f73e42")) + 
  scale_color_manual(values = c('Lombok' = "#8e86c5",
                                'Wakatobi' = "#007cb5",
                                'Misool' = "#a6c54b",
                                'Waigeo' = "#f73e42")) + 
  xlab("dbRDA1 [3.844%]") + ylab("dbRDA2 [1.915%]")
M04_plot 

M12_plot <- ggplot() + 
  geom_point(data = M12_merged, aes(dbRDA1, dbRDA2, fill = region), pch = 21, size = 3, alpha = 0.8) + 
  stat_ellipse(data = M12_merged, aes(dbRDA1, dbRDA2, color = region)) +
  geom_segment(data = M12_arrows, aes(x = 0, y = 0, xend = dbRDA1, yend = dbRDA2), arrow = arrow(length = unit(0.5, "cm"))) +
  geom_label(data = M12_arrows, aes(x = dbRDA1, y = dbRDA2, label = rownames(M12_arrows)), fill = "white", color = "black") +
  theme_bw() + theme(panel.grid = element_blank()) + 
  scale_fill_manual(values = c('Lombok' = "#8e86c5",
                               'Wakatobi' = "#007cb5",
                               'Misool' = "#a6c54b",
                               'Waigeo' = "#f73e42")) + 
  scale_color_manual(values = c('Lombok' = "#8e86c5",
                                'Wakatobi' = "#007cb5",
                                'Misool' = "#a6c54b",
                                'Waigeo' = "#f73e42")) + 
  xlab("dbRDA1 [2.854%]") + ylab("dbRDA2 [1.834%]")
M12_plot

S04_plot <- ggplot() + 
  geom_point(data = S04_merged, aes(dbRDA1, dbRDA2, fill = region), pch = 21, size = 3, alpha = 0.8) + 
  stat_ellipse(data = S04_merged, aes(dbRDA1, dbRDA2, color = region)) +
  geom_segment(data = S04_arrows, aes(x = 0, y = 0, xend = dbRDA1, yend = dbRDA2), arrow = arrow(length = unit(0.5, "cm"))) +
  geom_label(data = S04_arrows, aes(x = dbRDA1, y = dbRDA2, label = rownames(S04_arrows)), fill = "white", color = "black") +
  theme_bw() + theme(panel.grid = element_blank()) + 
  scale_fill_manual(values = c('Lombok' = "#8e86c5",
                               'Wakatobi' = "#007cb5",
                               'Misool' = "#a6c54b",
                               'Waigeo' = "#f73e42")) + 
  scale_color_manual(values = c('Lombok' = "#8e86c5",
                                'Wakatobi' = "#007cb5",
                                'Misool' = "#a6c54b",
                                'Waigeo' = "#f73e42")) + 
  xlab("dbRDA1 [3.10%]") + ylab("dbRDA2 [2.181%]")
S04_plot

S12_plot <- ggplot() + 
  geom_point(data = S12_merged, aes(dbRDA1, dbRDA2, fill = region), pch = 21, size = 3, alpha = 0.8) + 
  stat_ellipse(data = S12_merged, aes(dbRDA1, dbRDA2, color = region)) +
  geom_segment(data = S12_arrows, aes(x = 0, y = 0, xend = dbRDA1, yend = dbRDA2), arrow = arrow(length = unit(0.5, "cm"))) +
  geom_label(data = S12_arrows, aes(x = dbRDA1, y = dbRDA2, label = rownames(S12_arrows)), fill = "white", color = "black") +
  theme_bw() + theme(panel.grid = element_blank()) + 
  scale_fill_manual(values = c('Lombok' = "#8e86c5",
                               'Wakatobi' = "#007cb5",
                               'Misool' = "#a6c54b",
                               'Waigeo' = "#f73e42")) + 
  scale_color_manual(values = c('Lombok' = "#8e86c5",
                                'Wakatobi' = "#007cb5",
                                'Misool' = "#a6c54b",
                                'Waigeo' = "#f73e42")) + 
  xlab("dbRDA1 [2.853%]") + ylab("dbRDA2 [1.983%]")
S12_plot

#arrange plots into single figure
final_plot <- ggarrange(M04_plot, M12_plot, S04_plot, S12_plot, nrow = 2, ncol = 2,
                        common.legend = TRUE, align = "hv", legend = "bottom")
final_plot

#export plot and make final edits in illustrator to adjust label positions
ggsave("dbRDA_fig.pdf", final_plot, height = 8, width = 8)

#Run same analysis but with OTUs
#run dbRDA
dbRDA_OTU <- dbrda(otu_filtered ~ npp_mean + wave_mean + avg_complexity + cover + detritAbund + carniAbund + 
                 herbAbund + invertAbund + omniAbund + planktAbund + coraliAbund + hardCoral_percent + 
                 availSubstrate_percent + rubble_percent + algae_percent + softCoral_percent + 
                 abundance + target_biomass + nontarget_biomass + total_biomass, metadata_filtered, 
               distance = "bray", na.action = na.exclude)
summary(dbRDA_OTU)
anova(dbRDA_OTU, permutations = 999)
anova(dbRDA_OTU, by = "term", permutations = 999)
anova(dbRDA_OTU, by = "axis", permutations = 999)

dbRDA_OTU <- dbrda(asv_transformed ~ planktAbund, 
               metadata_filtered, distance = "bray", na.action = na.exclude)
summary(dbRDA_OTU)
anova(dbRDA_OTU, permutations = 999)
anova(dbRDA_OTU, by = "term", permutations = 999)
anova(dbRDA_OTU, by = "axis", permutations = 999)

#subset metadata and asvs by depth and filter size
metadata_M04 <- subset(metadata, subset == "M04")
metadata_M12 <- subset(metadata, subset == "M12")
metadata_S04 <- subset(metadata, subset == "S04")
metadata_S12 <- subset(metadata, subset == "S12")

otus_M04 <- as.data.frame(subset(otu_filtered, rownames(otu_filtered) %in% rownames(metadata_M04)))
otus_M12 <- as.data.frame(subset(otu_filtered, rownames(otu_filtered) %in% rownames(metadata_M12)))
otus_S04 <- as.data.frame(subset(otu_filtered, rownames(otu_filtered) %in% rownames(metadata_S04)))
otus_S12 <- as.data.frame(subset(otu_filtered, rownames(otu_filtered) %in% rownames(metadata_S12)))

metadata_S04 <- subset(metadata_S04, rownames(metadata_S04) %in% rownames(asvs_S04))

#Run dbRDA by subset
dbRDA_OTU_M04 <- dbrda(otus_M04 ~ rubble_percent, 
                     metadata_M04, distance = "bray", na.action = na.exclude)
summary(dbRDA_OTU_M04)
anova(dbRDA_OTU_M04, permutations = 999)
anova(dbRDA_OTU_M04, by = "term", permutations = 999)
anova(dbRDA_OTU_M04, by = "axis", permutations = 999)
vif.cca(dbRDA_OTU_M04)

dbRDA_OTU_M12 <- dbrda(otus_M12 ~ detritAbund, 
                       metadata_M12, distance = "bray", na.action = na.exclude)
summary(dbRDA_OTU_M12)
anova(dbRDA_OTU_M12, permutations = 999)
anova(dbRDA_OTU_M12, by = "term", permutations = 999)
anova(dbRDA_OTU_M12, by = "axis", permutations = 999)
vif.cca(dbRDA_OTU_M12)

dbRDA_OTU_S04 <- dbrda(otus_S04 ~ detritAbund, 
                     metadata_S04, distance = "bray", na.action = na.exclude)
summary(dbRDA_OTU_S04)
anova(dbRDA_OTU_S04, permutations = 999)
anova(dbRDA_OTU_S04, by = "term", permutations = 999)
anova(dbRDA_OTU_S04, by = "axis", permutations = 999)
vif.cca(dbRDA_OTU_S04)

dbRDA_OTU_S12 <- dbrda(otus_S12 ~ detritAbund + invertAbund,
                       metadata_S12, distance = "bray", na.action = na.exclude)
summary(dbRDA_OTU_S12)
anova(dbRDA_OTU_S12, permutations = 999)
anova(dbRDA_OTU_S12, by = "term", permutations = 999)
anova(dbRDA_OTU_S12, by = "axis", permutations = 999)
vif.cca(dbRDA_OTU_S12)
