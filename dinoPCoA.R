#load packages
library(vegan)
library(ape)
library(reshape2)
library(phyloseq)
library(ggplot2)
library(dplyr)

#set working directory
setwd("~/Desktop/DataAnalysis/threshold95/dino_all/")

#import data for PCoA graph
data <- read.csv("../SAR/SAR-percent-table-95.csv", row.names = 1, na.strings = c("","NA"))
dino_table <- subset(data, L2 == "Dinoflagellata")
dino_table <- dino_table[-c(1:8)]
dino_table_trans <- as.data.frame(t(dino_table))

#format data for phyloseq
otu <- otu_table(dino_table, taxa_are_rows = TRUE)
sample_data <- read.csv("~/Desktop/DataAnalysis/Master_V9_MappingFile.csv", row.names = 1)
tax_table <- read.csv("~/Desktop/DataAnalysis/pr2-taxonomy.csv", row.names = 1)
taxmat <- as.matrix(tax_table)
tree <- read.tree("~/Desktop/DataAnalysis/tree.nwk")
meta <- sample_data(sample_data)
tax <- tax_table(taxmat)

#construct phyloseq objects and filter objects to different subsets
physeq <- phyloseq(otu,meta,tax,tree)
physeq04 <- subset_samples(physeq, filter == "0.4um")
physeq12 <- subset_samples(physeq, filter == "12um")

physeqM <- subset_samples(physeq, depth == "M")
physeqM04 <- subset_samples(physeqM, filter == "0.4um")
physeqM12 <- subset_samples(physeqM, filter == "12um")

physeqS <- subset_samples(physeq, depth == "S")
physeqS04 <- subset_samples(physeqS, filter == "0.4um")
physeqS12 <- subset_samples(physeqS, filter == "12um")

#define region colors for PCoA graph
colors <- c('Lombok' = "#8e86c5",
            'Wakatobi' = "#007cb5",
            'Misool' = "#a6c54b",
            'Waigeo' = "#f73e42")

#run PCoA, extract scores, and plot ordination
pcoa <- ordinate(physeq04, "PCoA", "bray")
pcoa_coords <- as.data.frame(pcoa[["vectors"]])
pcoa_coords <- pcoa_coords[1:2]
pcoa_coords$sampleID <- rownames(pcoa_coords)
pcoa_coords_merged <- merge(pcoa_coords, sample_data, by = "row.names")
pcoa_plot <- ggplot(pcoa_coords_merged, aes(Axis.1, Axis.2)) + 
  stat_ellipse(aes(color = region)) + geom_point(size = 4, alpha = 0.85, 
                                                 aes(shape = interaction(subset, region), 
                                                     color = interaction(subset, region), 
                                                     fill = interaction(subset, region)), stroke = 0.75) + 
  theme_bw() + theme(panel.grid = element_blank()) +
  scale_fill_manual(values = c('M04.Lombok' = "#8e86c5", 'M04.Misool' = "#a6c54b", 'M04.Waigeo' = "#f73e42", 'M04.Wakatobi' = "#007cb5",
                               'S04.Lombok' = "white", 'S04.Misool' = "white", 'S04.Waigeo' = "white", 'S04.Wakatobi' = "white")) + 
  scale_shape_manual(values = c('M04.Lombok' = 21, 'M04.Misool' = 22, 'M04.Waigeo' = 23, 'M04.Wakatobi' = 24,
                                'S04.Lombok' = 21, 'S04.Misool' = 22, 'S04.Waigeo' = 23, 'S04.Wakatobi' = 24)) + 
  xlab("PC1 (12.1%)") + ylab("PC2 (7.8%)") + 
  scale_color_manual(values = c('M04.Lombok' = "black", 'M04.Misool' = "black", 'M04.Waigeo' = "black", 'M04.Wakatobi' = "black",
                                'S04.Lombok' = "#8e86c5", 'S04.Misool' = "#a6c54b", 'S04.Waigeo' = "#f73e42", 'S04.Wakatobi' = "#007cb5",
                                'Lombok' = "#8e86c5",'Wakatobi' = "#007cb5",'Misool' = "#a6c54b",'Waigeo' = "#f73e42"))
pcoa_plot

#import and format data for barplots (rows = ASVs, columns = sampleID)
data <- read.csv("~/Desktop/DataAnalysis/threshold95/SAR/SAR_OTU_percent_L3.csv")
#subset to only dinoflagellates
dino_data <- subset(data, L3 %in% c("Dinoflagellata","Dinophyceae","Ellobiophyceae","Noctilucophyceae","Syndiniales","Dinophyta_X"))
#transform dataset to long format
dino_data_long <- melt(dino_data)
#rename variable column to sampleID
names(dino_data_long)[2] <- "sampleID"
#import metadata
metadata <- read.csv("~/Desktop/DataAnalysis/Master_V9_mappingFile.csv")
#merge percent and metadata
merged_table <- merge(dino_data_long, metadata, by = "sampleID")
#sort region column into geographic order for plotting
merged_table$region <- factor(merged_table$region, levels = c("Lombok","Wakatobi","Misool","Waigeo"))

#average by region, subset, and taxa
mean_data <- aggregate(value~L3+region+subset, data = merged_table, FUN = mean)

#subset to only 0.4um water samples
mean_04_data <- subset(mean_data, subset %in% c("M04","S04"))

#define taxa colors
taxa_colors <- c('Dinoflagellata' = "black",
                 'Dinophyceae' = "#8dba5d",
                 'Ellobiophyceae' = "#35698f",
                 'Noctilucophyceae' = "#fcac31",
                 'Syndiniales' = "#85262a")

#plot barplot of dinoflagellates out of total SAR community
SARplot <- ggplot(mean_04_data, aes(region, value, fill = L3)) + geom_bar(stat = "identity") + facet_grid(~subset, space = "free", scales = "free") +
  ylim(0,1) + theme_bw() + theme(panel.grid = element_blank()) + theme(strip.background = element_rect(fill = "black"), strip.text = element_text(colour = "white", face = "bold")) +
  scale_fill_manual(values = taxa_colors) + theme(axis.text.x = element_blank()) + theme(axis.title.x = element_blank()) + ylab("Mean relative abundance") 
SARplot

#plot dinoflagellate community barplot
dinoplot <- ggplot(mean_04_data, aes(region, value, fill = L3)) + geom_bar(stat = "identity", position = "fill") + facet_grid(~subset, space = "free", scales = "free") +
  theme_bw() + theme(panel.grid = element_blank()) + theme(strip.background = element_rect(fill = "black"), strip.text = element_text(colour = "white", face = "bold")) +
  scale_fill_manual(values = taxa_colors) + theme(axis.text.x = element_text(angle = 90)) + theme(axis.title.x = element_blank()) + ylab("Mean relative abundance")
dinoplot


#arrange plots into one big plot with three panels
comboPlot <- ggarrange(pcoa_plot, ggarrange(SARplot, dinoplot, nrow = 2, ncol = 1, common.legend = TRUE,
                                            legend = "right", labels = c("b","c"), align = "hv"), nrow = 1, ncol = 2, 
                       legend = "none", labels = c("a", ""), align = "v")
comboPlot

ggsave("~/Desktop/dinoflagellateFig.pdf", comboPlot, width = 10, height = 6)
