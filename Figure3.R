#Load packages
library(ggplot2)
library(ggpubr)

#Set working directory
setwd("~/Desktop/DataAnalysis/threshold95/SAR/PCoA")
#Import data
data_M04 <- read.csv("PCoA_coords_M04.csv")
data_M12 <- read.csv("PCoA_coords_M12.csv")
data_S04 <- read.csv("PCoA_coords_S04.csv")
data_S12 <- read.csv("PCoA_coords_S12.csv")

data_M04$location <- factor(data_M04$location, levels = c("South Lombok",
                                                          "North Lombok",
                                                          "Wakatobi",
                                                          "Misool",
                                                          "Waigeo"))
data_M12$location <- factor(data_M12$location, levels = c("South Lombok",
                                                          "North Lombok",
                                                          "Wakatobi",
                                                          "Misool",
                                                          "Waigeo"))
data_S04$location <- factor(data_S04$location, levels = c("South Lombok",
                                                          "North Lombok",
                                                          "Wakatobi",
                                                          "Misool",
                                                          "Waigeo"))
data_S12$location <- factor(data_S12$location, levels = c("South Lombok",
                                                          "North Lombok",
                                                          "Wakatobi",
                                                          "Misool",
                                                          "Waigeo"))

data_M04$region <- factor(data_M04$region, levels = c("Lombok",
                                                      "Wakatobi",
                                                      "Misool",
                                                      "Waigeo"))
data_M12$region <- factor(data_M12$region, levels = c("Lombok",
                                                      "Wakatobi",
                                                      "Misool",
                                                      "Waigeo"))
data_S04$region <- factor(data_S04$region, levels = c("Lombok",
                                                      "Wakatobi",
                                                      "Misool",
                                                      "Waigeo"))
data_S12$region <- factor(data_S12$region, levels = c("Lombok",
                                                      "Wakatobi",
                                                      "Misool",
                                                      "Waigeo"))

#Define colors
colors <- c('South Lombok' = "#464586",
            'North Lombok' = "#8e86c5",
            'Wakatobi' = "#007cb5",
            'Misool' = "#a6c54b",
            'Waigeo' = "#f73e42")

colors <- c('Lombok' = "#8e86c5",
            'Wakatobi' = "#007cb5",
            'Misool' = "#a6c54b",
            'Waigeo' = "#f73e42")

#Plot PCoAs
M04_plot <- ggplot(data_M04, aes(Axis10.888827, Axis10.273335, fill = location)) +
                  geom_point(size = 3, alpha = 0.8, pch = 21, color = "black") + stat_ellipse(aes(color = location)) + theme_bw() +
                  theme(strip.background = element_rect(color = "black", fill = "black"), strip.text = element_text(color = "white", face = "bold", size = 12)) + 
                  theme(panel.grid = element_blank()) + scale_fill_manual(values = colors) + scale_color_manual(values = colors) +
                  xlab("PC1 [10.89%]") + ylab("PC2 [10.27%]")
M04_plot

M12_plot <- ggplot(data_M12, aes(Axis13.778467, Axis6.232042, fill = location)) +
  geom_point(size = 3, alpha = 0.8, pch = 21, color = "black") + stat_ellipse(aes(color = location)) + theme_bw() +
  theme(strip.background = element_rect(color = "black", fill = "black"), strip.text = element_text(color = "white", face = "bold", size = 12)) + 
  theme(panel.grid = element_blank()) + scale_fill_manual(values = colors) + scale_color_manual(values = colors) +
  xlab("PC1 [13.78%]") + ylab("PC2 [6.23%]")
M12_plot

S04_plot <- ggplot(data_S04, aes(Axis7.757232, Axis6.925851, fill = location)) +
  geom_point(size = 3, alpha = 0.8, pch = 21, color = "black") + stat_ellipse(aes(color = location)) + theme_bw() +
  theme(strip.background = element_rect(color = "black", fill = "black"), strip.text = element_text(color = "white", face = "bold", size = 12)) + 
  theme(panel.grid = element_blank()) + scale_fill_manual(values = colors) + scale_color_manual(values = colors) +
  xlab("PC1 [7.76%]") + ylab("PC2 [6.93%]")
S04_plot

S12_plot <- ggplot(data_S12, aes(Axis6.92244, Axis4.881823, fill = location)) +
  geom_point(size = 3, alpha = 0.8, pch = 21, color = "black") + stat_ellipse(aes(color = location)) + theme_bw() +
  theme(strip.background = element_rect(color = "black", fill = "black"), strip.text = element_text(color = "white", face = "bold", size = 12)) + 
  theme(panel.grid = element_blank()) + scale_fill_manual(values = colors) + scale_color_manual(values = colors) +
  xlab("PC1 [6.92%]") + ylab("PC2 [4.88%]")
S12_plot

plot <- ggarrange(M04_plot, M12_plot, S04_plot, S12_plot, 
                  ncol = 2, nrow = 2, common.legend = TRUE, 
                  legend = "right", align = "hv")
plot

ggsave("PCoA_SAR_95_all_by_location.png", plot, width = 10, height = 8)

#Change working directory
setwd("~/Desktop/DataAnalysis/threshold95/SAR/ANOSIM/distances")
#Import distance data
distances_M04 <- read.csv("M04-distances-region.csv")
distances_M12 <- read.csv("M12-distances-region.csv")
distances_S04 <- read.csv("S04-distances-region.csv")
distances_S12 <- read.csv("S12-distances-region.csv")

#Reorder location columns
distances_M04$Group1 <- factor(distances_M04$Group1, levels = c("Lombok",
                                                                "Wakatobi",
                                                                "Misool",
                                                                "Waigeo"))
distances_M12$Group1 <- factor(distances_M12$Group1, levels = c("Lombok",
                                                                "Wakatobi",
                                                                "Misool",
                                                                "Waigeo"))
distances_S04$Group1 <- factor(distances_S04$Group1, levels = c("Lombok",
                                                                "Wakatobi",
                                                                "Misool",
                                                                "Waigeo"))
distances_S12$Group1 <- factor(distances_S12$Group1, levels = c("Lombok",
                                                                "Wakatobi",
                                                                "Misool",
                                                                "Waigeo"))

distances_M04$Group2 <- factor(distances_M04$Group2, levels = c("Lombok",
                                                                "Wakatobi",
                                                                "Misool",
                                                                "Waigeo"))
distances_M12$Group2 <- factor(distances_M12$Group2, levels = c("Lombok",
                                                                "Wakatobi",
                                                                "Misool",
                                                                "Waigeo"))
distances_S04$Group2 <- factor(distances_S04$Group2, levels = c("Lombok",
                                                                "Wakatobi",
                                                                "Misool",
                                                                "Waigeo"))
distances_S12$Group2 <- factor(distances_S12$Group2, levels = c("Lombok",
                                                                "Wakatobi",
                                                                "Misool",
                                                                "Waigeo"))

distances_M04$Group1 <- factor(distances_M04$Group1, levels = c("South Lombok",
                                                                "North Lombok",
                                                                "Wakatobi",
                                                                "Misool",
                                                                "Waigeo"))
distances_M12$Group1 <- factor(distances_M12$Group1, levels = c("South Lombok",
                                                                "North Lombok",
                                                                "Wakatobi",
                                                                "Misool",
                                                                "Waigeo"))
distances_S04$Group1 <- factor(distances_S04$Group1, levels = c("South Lombok",
                                                                "North Lombok",
                                                                "Wakatobi",
                                                                "Misool",
                                                                "Waigeo"))
distances_S12$Group1 <- factor(distances_S12$Group1, levels = c("South Lombok",
                                                                "North Lombok",
                                                                "Wakatobi",
                                                                "Misool",
                                                                "Waigeo"))

distances_M04$Group2 <- factor(distances_M04$Group2, levels = c("South Lombok",
                                                                "North Lombok",
                                                                "Wakatobi",
                                                                "Misool",
                                                                "Waigeo"))
distances_M12$Group2 <- factor(distances_M12$Group2, levels = c("South Lombok",
                                                                "North Lombok",
                                                                "Wakatobi",
                                                                "Misool",
                                                                "Waigeo"))
distances_S04$Group2 <- factor(distances_S04$Group2, levels = c("South Lombok",
                                                                "North Lombok",
                                                                "Wakatobi",
                                                                "Misool",
                                                                "Waigeo"))
distances_S12$Group2 <- factor(distances_S12$Group2, levels = c("South Lombok",
                                                                "North Lombok",
                                                                "Wakatobi",
                                                                "Misool",
                                                                "Waigeo"))

#Define colors
colors <- c('South Lombok' = "#464586",
            'North Lombok' = "#8e86c5",
            'Wakatobi' = "#007cb5",
            'Misool' = "#a6c54b",
            'Waigeo' = "#f73e42")

colors <- c('Lombok' = "#8e86c5",
            'Wakatobi' = "#007cb5",
            'Misool' = "#a6c54b",
            'Waigeo' = "#f73e42")

#plot boxplots
M04_distance_plot <- ggplot(distances_M04, aes(Group1, Distance, color = Group2, fill = Group2)) +
                            geom_boxplot(position = "dodge2", color = "black", alpha = 0.8) +
                            theme_bw() + theme(panel.grid = element_blank()) + theme(axis.title.x = element_blank()) +
                            theme(strip.background = element_rect(color = "black", fill = "black"), strip.text = element_text(color = "white", face = "bold", size = 12)) +
                            scale_color_manual(values = colors) + scale_fill_manual(values = colors, na.value = "white") + ylab("Bray-Curtis Distance")
M04_distance_plot

M12_distance_plot <- ggplot(distances_M12, aes(Group1, Distance, fill = Group2, color = Group2, alpha = 0.5)) +
  geom_boxplot(position = "dodge2", color = "black", alpha = 0.8) +
  theme_bw() + theme(panel.grid = element_blank()) + theme(axis.title.x = element_blank()) +
  theme(strip.background = element_rect(color = "black", fill = "black"), strip.text = element_text(color = "white", face = "bold", size = 12)) +
  scale_color_manual(values = colors) + scale_fill_manual(values = colors) + ylab("Bray-Curtis Distance")
M12_distance_plot

S04_distance_plot <- ggplot(distances_S04, aes(Group1, Distance, fill = Group2, color = Group2, alpha = 0.5)) +
  geom_boxplot(position = "dodge2", color = "black", alpha = 0.8) +
  theme_bw() + theme(panel.grid = element_blank()) + theme(axis.title.x = element_blank()) +
  theme(strip.background = element_rect(color = "black", fill = "black"), strip.text = element_text(color = "white", face = "bold", size = 12)) +
  scale_color_manual(values = colors) + scale_fill_manual(values = colors) + ylab("Bray-Curtis Distance")
S04_distance_plot

S12_distance_plot <- ggplot(distances_S12, aes(Group1, Distance, fill = Group2, color = Group2, alpha = 0.5)) +
  geom_boxplot(position = "dodge2", color = "black", alpha = 0.8) +
  theme_bw() + theme(panel.grid = element_blank()) + theme(axis.title.x = element_blank()) +
  theme(strip.background = element_rect(color = "black", fill = "black"), strip.text = element_text(color = "white", face = "bold", size = 12)) +
  scale_color_manual(values = colors) + scale_fill_manual(values = colors) + ylab("Bray-Curtis Distance")
S12_distance_plot

plot <- ggarrange(M04_plot, M04_distance_plot, M12_plot, M12_distance_plot,
                  S04_plot, S04_distance_plot, S12_plot, S12_distance_plot,
                  ncol = 2, nrow = 4, widths = c(1,2), common.legend = TRUE,
                  legend = "bottom", align = "hv")
plot

ggsave("BrayCurtis_combined_location_region.png", plot, height = 10, width = 10)


