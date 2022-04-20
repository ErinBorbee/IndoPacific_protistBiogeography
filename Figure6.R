#Load packages
library(vegan)
library(ggplot2)
library(reshape2)
library(ExcelFunctionsR)

#Set workings directory
setwd("~/Desktop/DataAnalysis/threshold95/richness")

#import data
data <- read.csv("richness.csv", row.names = 1)
dataM04 <- subset(data, subset == "M04")
dataM12 <- subset(data, subset == "M12")
dataS04 <- subset(data, subset == "S04")
dataS12 <- subset(data, subset == "S12")

data_long <- read.csv("richness_long.csv") 
  
data_long <- data[,-c(10:14)]
data_long <- data_long[,-c(8:12)]
data_long <- melt(data)
write.csv(data_long,"richnessLong.csv")

data_long <- read.csv("richnessLong.csv")
data_long <- subset(data_long, type == "OTUrichness")

plot <- ggplot(data_long, aes(management_type,value)) + geom_boxplot() + facet_grid(region~variable, scales = "free")
plot

data$region <- factor(data$region, levels = c("Lombok","Wakatobi","Misool","Waigeo"))
data_long$region <- factor(data_long$region, levels = c("Lombok","Wakatobi","Misool","Waigeo"))
data_long$variable <- factor(data_long$variable, levels = c("All","Dinoflagellata","Bacillariophyta"))

#Create vectors for each richness measure
SAR_ASV_richness <- data$SAR_ASVrichness
SAR_OTU_richness <- data$SAR_OTUrichness

diatom_ASV_richness <- data$diatom_ASVrichness
diatom_OTU_richness <- data$diatom_OTUrichness

dinoAll_ASV_richness <- data$dinoAll_ASVrichness
dinoAll_OTU_richness <- data$dinoAll_OTUrichness

dinophyceae_ASV_richness <- data$dinophyceae_ASVrichness
dinophyceae_OTU_richness <- data$dinophyceae_OTUrichness

#Plot histograms and test for normality in each richness variable
#If data are not normal, transform and test for normality again
#Normal data = p > 0.05 
hist(SAR_ASV_richness, main = "Richness", breaks = 10)
shapiro.test(SAR_ASV_richness)

SAR_ASV_richness_trans <- (SAR_ASV_richness)^0.75
hist(SAR_ASV_richness_trans, main = "Richness", breaks = 10)
shapiro.test(SAR_ASV_richness_trans)
#All: Normal after sqrt transformation
#  M: Normal after sqrt transformation
#M12: Normal after sqrt transformation
#M04: Normal after ^1.25 transformation
#  S: Normal after ^0.75 transformation
#S12: Not normal after transformation
#S04: Normal without transformation

hist(SAR_OTU_richness, main = "Richness", breaks = 10)
shapiro.test(SAR_OTU_richness)

SAR_OTU_richness_trans <- (SAR_OTU_richness)^0.8
hist(SAR_OTU_richness_trans, main = "Richness", breaks = 10)
shapiro.test(SAR_OTU_richness_trans)
#All: Not normal after transformation
#  M: Not normal after transformation
#M12: Normal after sqrt transformation
#M04: Normal after ^1.25 transformation
#  S: Not normal after transformation
#S12: Not normal after transformation
#S04: Normal without transformation

hist(diatom_ASV_richness, main = "Richness", breaks = 10)
shapiro.test(diatom_ASV_richness)

diatom_ASV_richness_trans <- (diatom_ASV_richness)^0.5
hist(diatom_ASV_richness_trans, main = "Richness", breaks = 10)
shapiro.test(diatom_ASV_richness_trans)
#All: Not normal after transformation
#  M: Not normal after transformation
#M12: Not normal after transformation
#M04: Normal after sqrt transformation
#  S: Normal without transformation
#S12: Normal without transformation
#S04: Normal after sqrt transformation

hist(diatom_OTU_richness, main = "Richness", breaks = 10)
shapiro.test(diatom_OTU_richness)

diatom_OTU_richness_trans <- (diatom_OTU_richness)^0.5
hist(diatom_OTU_richness_trans, main = "Richness", breaks = 10)
shapiro.test(diatom_OTU_richness_trans)
#All: Normal after sqrt transformation
#  M: Normal after ^0.75 transformation
#M12: Normal after sqrt transformation
#M04: Normal after ^0.9 transformation
#  S: Normal after sqrt transformation
#S12: Normal after sqrt transformation
#S04: Normal without transformation

hist(dinoAll_ASV_richness, main = "Richness", breaks = 10)
shapiro.test(dinoAll_ASV_richness)

dinoAll_ASV_richness_trans <- (dinoAll_ASV_richness)^0.5
hist(dinoAll_ASV_richness_trans, main = "Richness", breaks = 10)
shapiro.test(dinoAll_ASV_richness_trans)
#All: Not normal after transformation
#  M: Not normal after transformation
#M12: Not normal after transformation
#M04: Normal without transformation
#  S: Normal after sqrt transformation
#S12: Normal after sqrt transformation
#S04: Normal after sqrt transformation

hist(dinoAll_OTU_richness, main = "Richness", breaks = 10)
shapiro.test(dinoAll_OTU_richness)

dinoAll_OTU_richness_trans <- (dinoAll_ASV_richness_trans)^0.75
hist(dinoAll_OTU_richness_trans, main = "Richness", breaks = 10)
shapiro.test(dinoAll_OTU_richness_trans)
#All: Not normal after transformation
#  M: Not normal after transformation
#M12: Normal after sqrt transformation
#M04: Normal after sqrt transformation
#  S: Normal after ^0.75 transformation
#S12: Normal after ^0.75 transformation
#S04: Normal without transformation

hist(dinophyceae_ASV_richness, main = "Richness", breaks = 10)
shapiro.test(dinophyceae_ASV_richness)

dinophyceae_ASV_richness_trans <- sqrt(dinophyceae_ASV_richness)
hist(dinophyceae_ASV_richness_trans, main = "Richness", breaks = 10)
shapiro.test(dinophyceae_ASV_richness_trans)
#All: Normal after sqrt transformation
#  M: Normal after sqrt transformation
#M12: Normal after sqrt transformation
#M04: Normal without transformation
#  S: Normal after sqrt transformation
#S12: Normal after sqrt transformation
#S04: Normal after sqrt transformation

hist(dinophyceae_OTU_richness, main = "Richness", breaks = 10)
shapiro.test(dinophyceae_OTU_richness)

dinophyceae_OTU_richness_trans <- (dinophyceae_OTU_richness)^0.95
hist(dinophyceae_OTU_richness_trans, main = "Richness", breaks = 10)
shapiro.test(dinophyceae_OTU_richness_trans)
#All: Not normal after transformation
#  M: Normal after sqrt transformation
#M12: Normal after sqrt transformation
#M04: Normal without transformation
#  S: Not normal after transformation
#S12: Not normal after transformation
#S04: Normal without transformation

#Since most variable were not normal after transformation we use nonparametric stats
#Parametric = ANOVA & Tukey
#Nonparametric = Kruskal-Wallis & Wilcoxon

#For comparison of richness by designated metadata variable we use Kruskal-Wallis
#If Kruskal-Wallis is significant (p < 0.05) then we use Wilcoxon for pairwise comparisons
kruskal.test(SAR_ASVrichness ~ region, data = data)
#All: p-value is significant (p < 2.2e-16)
#  M: p-value is significant (p = 6.013e-08)
#M12: p-value is significant (p = 4.057e-06)
#M04: p-value is significant (p = 0.0007241)
#  S: p-value significant (p = 7.54e-13)
#S12: p-value significant (p = 1.54e-07)
#S04: p-value significant (p = 1.08e-05)
pairwise.wilcox.test(dataS12$SAR_ASVrichness, dataS12$region, p.adjust.method = "BH")
#All: Waigeo is only region with significantly different ASV richness for SAR
#  M: WGO differs from all other regions
#M12: WGO + WK, WGO + MIS, and LBK + WK are significantly different
#M04: WGO differs from all other regions
#  S: WGO differs from all other regions
#S12: WGO differs from all other regions
#S04: WGO differs from all other regions

kruskal.test(SAR_OTUrichness ~ region, data = data)
#All: p-value is significant (p < 2.2e-16)
#  M: p-value is significant (p = 6.209e-08)
#M12: p-value is significant (p = 5.514e-07)
#M04: p-value is significant (p = 0.0007414)
#  S: p-value is significant (p = 3.478e-13)
#S12: p-value is significant (p = 1.064e-07)
#S12: p-value is significant (p = 1.277e-06)
pairwise.wilcox.test(dataS12$SAR_OTUrichness, dataS12$region, p.adjust.method = "BH")
#All: Waigeo is only region with significantly different OTU richness for SAR
#  M: WGO differs from all other regions
#M12: WGO is different from all regions, LBK + WK are also different
#M04: WGO differs from all other regions
#  S: WGO differs from all other regions
#S12: WGO differs from all other regions
#S04: WGO differs from all other regions and also WK + LBK differ

kruskal.test(diatom_ASVrichness ~ region, data = data)
#All: p-value not significant (p = 0.3266)
#  M: p-value significant (p = 0.02277)
#M12: p-value not significant (p = 0.4415)
#M04: p-value significant (p = 0.02118)
#  S: p-value significant (p = 0.0253)
#S12: p-value significant (p = 0.008323)
#S04: p-value significant (p = 0.01591)
pairwise.wilcox.test(dataS12$diatom_ASVrichness, dataS12$region, p.adjust.method = "BH")
#All: No significant differences between any pairs
#  M: WGO differs from MIS and LBK
#M12: No significant differences between any pairs
#M04: No significant differences between any pairs
#  S: No significant differences between any pairs
#S12: WK differs from all other regions
#S12: WGO differs from all other regions except WK

kruskal.test(diatom_OTUrichness ~ region, data = data )
#All: p-value not significant (p = 0.1012)
#  M: p-value not significant (p = 0.454)
#M12: p-value significant (p = 0.008303)
#M04: p-value significant (p = 0.04889)
#  S: p-value not significant (p = 0.08192)
#S12: p-value significant (p = 0.001446)
#S04: p-value not significant (p = 0.4841)
pairwise.wilcox.test(dataS12$diatom_OTUrichness, dataS12$region, p.adjust.method = "BH")
#All: No signicant differences between any pairs
#  M: No signicant differences between any pairs
#M12: LBK + MIS and WK + WGO are the same
#M04: No significant differences between any pairs
#  S: No significant differences between any pairs
#S12: WGO differs from all other regions
#S04: No significant differences between any pairs

kruskal.test(dinoAll_ASVrichness ~ region, data = data)
#All: p-value significant (p < 2.2e-16)
#  M: p-value significant (p = 1.93e-09)
#M12: p-value significant (p = 3.031e-09)
#M04: p-value significant (p = 0.0004574)
#  S: p-value significant (p = 1.2e-07)
#S12: p-value significant (p = 4.894e-05)
#S04: p-value significant (p = 0.001349)
pairwise.wilcox.test(dataS12$dinoAll_ASVrichness, dataS12$region, p.adjust.method = "BH")
#All: WK + LBK are only pair that are the same
#  M: WGO differs from all regions
#M12: WGO differs from all and also LBK + MIS differ
#M04: WGO differs from all but WK
#  S: WGO differs from all and also LBK + MIS differ
#S12: WGO differs from all regions
#S04: WGO differs from all but LBK

kruskal.test(dinoAll_OTUrichness ~ region, data = data)
#All: p-value significant (p < 2.2e-16)
#M12: p-value significant (p = 6.756e-10)
#M12: p-value significant (p = 1.597e-08)
#M04: p-value significant (p = 9.547e-05)
#  S: p-value significant (p = 9.016e-10)
#S12: p-value significant (p = 2.01e-05)
#S04: p-value significant (p = 3.509e-05)
pairwise.wilcox.test(dataS12$dinoAll_OTUrichness, dataS12$region, p.adjust.method = "BH")
#All: Waigeo significantly differs from all regions and LBK and MIS are also different
#  M: WGO differs from all regions
#M12: WGO is different from all and also LBK + MIS differ
#M04: WGO differs from all but WK
#  S: WGO differs from all and LBK + MIS also differ
#S12: WGO differs from all regions
#S04: All pairs differ except WGO + LBK and MIS + WK

kruskal.test(dinophyceae_ASVrichness ~ region, data = data)
#All: p-value significant (p < 2.2e-16)
#  M: p-value significant (p = 2.541e-11)
#M12: p-value significant (p = 8.2e-08)
#M04: p-value significant (p = 0.0002627)
#  S: p-value significant (p = 1.869e-07)
#S12: p-value significant (p = 5.934e-05)
#S04: p-value significant (p = 0.007605)
pairwise.wilcox.test(data$dinophyceae_ASVrichness, data$region, p.adjust.method = "BH")
#All: All pairs significantly different except for LBK and WK
#  M: WGO and MIS differ from each other and all other regions
#M12: LBK + WK is the only non-significant pair
#M04: MIS differs from WGO and WK
#  S: WGO differs from all and LBK + MIS also differ
#S12: WGO differs from all and LBK + MIS also differ
#S04: WGO differs from all but LBK

kruskal.test(dinophyceae_OTUrichness ~ region, data = data)
#All: p-value significant (p < 2.2e-16)
#  M: p-value significant (p = 4.581e-13)
#M12: p-value significant (p = 7.797e-09)
#M04: p-value significant (p = 7.34e-05)
#  S: p-value significant (p = 4.561e-09)
#S12: p-value significant (p = 5.074e-05)
#S04: p-value significant (p = 6.606e-05)
pairwise.wilcox.test(data$dinophyceae_OTUrichness, data$region, p.adjust.method = "BH")
#All: All pairs significantly different
#  M: WGO and MIS differ from each other and all other regions
#M12: LBK + WK is the only non-significant pair
#M04: WGO differs from all
#  S: LBK and WGO differ from each other and all other regions
#S12: WGO differs from all and LBK + MIS also differ
#S04: WGO differs from all regions

kruskal.test(ciliate_ASVrichness ~ region, data = data)
#All: p-value significant (p < 2.2e-16)
#  M: p-value significant (p = 2.541e-11)
#M12: p-value significant (p = 8.2e-08)
#M04: p-value significant (p = 0.0002627)
#  S: p-value significant (p = 1.869e-07)
#S12: p-value significant (p = 5.934e-05)
#S04: p-value significant (p = 0.007605)
pairwise.wilcox.test(dataS12$ciliate_ASVrichness, dataS12$region, p.adjust.method = "BH")
#All: All pairs significantly different except for LBK and WK
#  M: WGO and MIS differ from each other and all other regions
#M12: LBK + WK is the only non-significant pair
#M04: MIS differs from WGO and WK
#  S: WGO differs from all and LBK + MIS also differ
#S12: WGO differs from all and LBK + MIS also differ
#S04: WGO differs from all but LBK

kruskal.test(ciliate_OTUrichness ~ region, data = data)
#All: p-value significant (p < 2.2e-16)
#  M: p-value significant (p = 4.581e-13)
#M12: p-value significant (p = 7.797e-09)
#M04: p-value significant (p = 7.34e-05)
#  S: p-value significant (p = 4.561e-09)
#S12: p-value significant (p = 5.074e-05)
#S04: p-value significant (p = 6.606e-05)
pairwise.wilcox.test(dataM04$ciliate_OTUrichness, dataM04$region, p.adjust.method = "BH")
#All: All pairs significantly different
#  M: WGO and MIS differ from each other and all other regions
#M12: LBK + WK is the only non-significant pair
#M04: WGO differs from all
#  S: LBK and WGO differ from each other and all other regions
#S12: WGO differs from all and LBK + MIS also differ
#S04: WGO differs from all regions

colors <- c('high' = "#f73e42",
            'low' = "#007cb5",
            'none' = "white")

SAR_plot <- ggplot(data, aes(region, SAR_OTUrichness, color = region, fill = region, alpha = 0.5)) +
              geom_boxplot() + facet_grid(depth~filter) + theme_bw() + theme(panel.grid = element_blank()) +
              theme(strip.background = element_rect(fill = "grey20", color = "grey20", size = 1), strip.text = element_text(color = "white", face = "bold")) +
              theme(panel.border = element_rect(color = "grey20", size = 1)) + scale_fill_manual(values = colors) + scale_color_manual(values = colors) +
              theme(legend.position = "none") + theme(axis.title.x = element_blank()) +
              ylab("OTU Richness\n")
SAR_plot
ggsave("SAR_OTUrichness_all_by_location.png", SAR_plot, width = 9, height = 8)

diatom_plot <- ggplot(data, aes(region, diatom_OTUrichness, color = region, fill = region, alpha = 0.5)) +
  geom_boxplot() + facet_grid(depth~filter) + theme_bw() + theme(panel.grid = element_blank()) +
  theme(strip.background = element_rect(fill = "grey20"), strip.text = element_text(color = "white", face = "bold")) +
  scale_fill_manual(values = colors) + scale_color_manual(values = colors)
diatom_plot
ggsave("diatom_OTUrichness_all_by_location.png", diatom_plot, width = 9, height = 8)

dinoAll_plot <- ggplot(data, aes(region, dinoAll_OTUrichness, color = region, fill = region, alpha = 0.5)) +
  geom_boxplot() + facet_grid(depth~filter) + theme_bw() + theme(panel.grid = element_blank()) +
  theme(strip.background = element_rect(fill = "grey20"), strip.text = element_text(color = "white", face = "bold")) +
  scale_fill_manual(values = colors) + scale_color_manual(values = colors)
dinoAll_plot
ggsave("dinoAll_OTUrichness_all_by_location.png", dinoAll_plot, width = 9, height = 8)

dinophyceae_plot <- ggplot(data, aes(region, dinophyceae_OTUrichness, color = region, fill = region, alpha = 0.5)) +
  geom_boxplot() + facet_grid(depth~filter) + theme_bw() + theme(panel.grid = element_blank()) +
  theme(strip.background = element_rect(fill = "grey20"), strip.text = element_text(color = "white", face = "bold")) +
  scale_fill_manual(values = colors) + scale_color_manual(values = colors)
dinophyceae_plot
ggsave("dinophyceae_OTUrichness_all_by_location.png", dinophyceae_plot, width = 9, height = 8)

plot <- ggplot(data_long, aes(region, value, fill = color)) +
              geom_boxplot(color = "black", alpha = 0.8) + facet_grid(variable~subset, scales = "free") + theme_bw() + theme(panel.grid = element_blank()) +
              theme(strip.background = element_rect(fill = "grey20", color = "grey20", size = 1), strip.text = element_text(color = "white", face = "bold")) +
              theme(panel.border = element_rect(color = "grey20", size = 1)) + scale_fill_manual(values = colors) + scale_color_manual(values = colors) +
              theme(legend.position = "none") + theme(axis.title.x = element_blank()) +
              ylab("OTU Richness\n")
plot

ggsave("richness_boxplot_by_subset_all.png", plot, width = 10, height = 7)
