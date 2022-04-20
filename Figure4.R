#Load packages
library(phyloseq)
library(ape)
library(ggplot2)
library(vegan)
library(ggpubr)

#Set working directory
setwd("~/Desktop/Dissertation/Chapter2/")

#Import data
otu_table <- read.csv("~/Desktop/DataAnalysis/threshold95/SAR/SAR-percent-table-95.csv", row.names = 1)
otu_table <- otu_table[,-c(1:8)]
sample_data <- read.csv("~/Desktop/DataAnalysis/Master_V9_MappingFile.csv", row.names = 1)
tax_table <- read.csv("~/Desktop/DataAnalysis/pr2-taxonomy.csv", row.names = 1)
taxmat <- as.matrix(tax_table)
tree <- read.tree("~/Desktop/DataAnalysis/tree.nwk")

#Order the location column of metadata file
sample_data$location <- factor(sample_data$location, levels = c("South Lombok","North Lombok","Wakatobi","Misool","Waigeo"))

#Convert data into phyloseq format
otu <- otu_table(otu_table, taxa_are_rows = TRUE)
meta <- sample_data(sample_data)
tax <- tax_table(taxmat)

#Generate phyloseq object
physeq <- phyloseq(otu,meta,tax,tree)

#Subset phyloseq object
physeq_LBK <- subset_samples(physeq, region == "Lombok")
physeq_WK <- subset_samples(physeq, region == "Wakatobi")
physeq_MIS <- subset_samples(physeq, region == "Misool")
physeq_WGO <- subset_samples(physeq, region == "Waigeo")

physeq_LBK_M04 <- subset_samples(physeq_LBK, subset == "M04")
physeq_LBK_M12 <- subset_samples(physeq_LBK, subset == "M12")
physeq_LBK_S04 <- subset_samples(physeq_LBK, subset == "S04")
physeq_LBK_S12 <- subset_samples(physeq_LBK, subset == "S12")

physeq_WK_M04 <- subset_samples(physeq_WK, subset == "M04")
physeq_WK_M12 <- subset_samples(physeq_WK, subset == "M12")
physeq_WK_S04 <- subset_samples(physeq_WK, subset == "S04")
physeq_WK_S12 <- subset_samples(physeq_WK, subset == "S12")

physeq_MIS_M04 <- subset_samples(physeq_MIS, subset == "M04")
physeq_MIS_M12 <- subset_samples(physeq_MIS, subset == "M12")
physeq_MIS_S04 <- subset_samples(physeq_MIS, subset == "S04")
physeq_MIS_S12 <- subset_samples(physeq_MIS, subset == "S12")

physeq_WGO_M04 <- subset_samples(physeq_WGO, subset == "M04")
physeq_WGO_M12 <- subset_samples(physeq_WGO, subset == "M12")
physeq_WGO_S04 <- subset_samples(physeq_WGO, subset == "S04")
physeq_WGO_S12 <- subset_samples(physeq_WGO, subset == "S12")

#Run and plot PCoA ordination
PCoA_LBK_M04 <- ordinate(physeq_LBK_M04, "PCoA","bray")

plot_ordination(physeq_LBK_M04, PCoA_LBK_M04, type = "sites", color = "restriction") +
  geom_point(size = 4, aes(shape = restriction)) + theme_bw() + scale_shape_manual(values = c('fishing prohibited' = 21,
                                                                                              'fishing gear restriction' = 22,
                                                                                              'unrestricted fishing' = 24)) +
  facet_grid(region~subset, scales = "free") + stat_ellipse(mapping = aes())

#Extract coordinates from ordination
PCoA_LBK_M04_scores <- as.data.frame(PCoA_LBK_M04[["vectors"]])
PCoA_LBK_M04_scores <- PCoA_LBK_M04_scores[1:2]
names(PCoA_LBK_M04_scores)[1] <- "PC1"
names(PCoA_LBK_M04_scores)[2] <- "PC2"

#Merge coordinates with metadata
PCoA_LBK_M04_scores <- merge(sample_data, PCoA_LBK_M04_scores, by = "row.names", all = FALSE, no.dups = TRUE)

#Plot formatted ordination
PCoA_LBK_M04_plot <- ggplot(PCoA_LBK_M04_scores, aes(PC1, PC2, fill = restriction, shape = restriction)) + 
  geom_point(size = 4, alpha = 0.75, color = "black") +
  stat_ellipse(aes(color = restriction)) + theme_bw() + theme(panel.grid = element_blank()) +
  scale_shape_manual(values = c('fishing prohibited' = 21,
                                'fishing gear restriction' = 22,
                                'unrestricted fishing' = 24)) +
  scale_fill_manual(values = c('fishing prohibited' = "#f73e42",
                               'fishing gear restriction' = "#8e86c5",
                               'unrestricted fishing' = "#007cb5")) +
  scale_color_manual(values = c('fishing prohibited' = "#f73e42",
                               'fishing gear restriction' = "#8e86c5",
                               'unrestricted fishing' = "#007cb5")) +
  xlab("PC1 [25.8%]") + ylab("PC2 [21.2%]") +
  facet_grid(region ~ subset, scales = "free") + 
  theme(strip.text.x = element_text(face = "bold", color = "white"), 
        strip.background.x = element_rect(fill = "black"), 
        strip.background.y = element_blank(), 
        strip.text.y = element_blank()) 
PCoA_LBK_M04_plot

#Run and plot PCoA ordination
PCoA_LBK_M12 <- ordinate(physeq_LBK_M12, "PCoA","bray")

plot_ordination(physeq_LBK_M12, PCoA_LBK_M12, type = "sites", color = "restriction") +
  geom_point(size = 4, aes(shape = restriction)) + theme_bw() + scale_shape_manual(values = c('fishing prohibited' = 21,
                                                                                              'fishing gear restriction' = 22,
                                                                                              'unrestricted fishing' = 24)) +
  facet_grid(region~subset, scales = "free") + stat_ellipse(mapping = aes())

#Extract coordinates from ordination
PCoA_LBK_M12_scores <- as.data.frame(PCoA_LBK_M12[["vectors"]])
PCoA_LBK_M12_scores <- PCoA_LBK_M12_scores[1:2]
names(PCoA_LBK_M12_scores)[1] <- "PC1"
names(PCoA_LBK_M12_scores)[2] <- "PC2"

#Merge coordinates with metadata
PCoA_LBK_M12_scores <- merge(sample_data, PCoA_LBK_M12_scores, by = "row.names", all = FALSE, no.dups = TRUE)

#Plot formatted ordination
PCoA_LBK_M12_plot <- ggplot(PCoA_LBK_M12_scores, aes(PC1, PC2, fill = restriction, shape = restriction)) + 
  geom_point(size = 4, alpha = 0.75, color = "black") +
  stat_ellipse(aes(color = restriction)) + theme_bw() + theme(panel.grid = element_blank()) +
  scale_shape_manual(values = c('fishing prohibited' = 21,
                                'fishing gear restriction' = 22,
                                'unrestricted fishing' = 24)) +
  scale_fill_manual(values = c('fishing prohibited' = "#f73e42",
                               'fishing gear restriction' = "#8e86c5",
                               'unrestricted fishing' = "#007cb5")) +
  scale_color_manual(values = c('fishing prohibited' = "#f73e42",
                                'fishing gear restriction' = "#8e86c5",
                                'unrestricted fishing' = "#007cb5")) +
  xlab("PC1 [18.2%]") + ylab("PC2 [15.9%]") +
  facet_grid(region ~ subset, scales = "free") + 
  theme(strip.text.x = element_text(face = "bold", color = "white"), 
        strip.background.x = element_rect(fill = "black"), 
        strip.background.y = element_blank(), 
        strip.text.y = element_blank()) 
PCoA_LBK_M12_plot

#Run and plot PCoA ordination
PCoA_LBK_S04 <- ordinate(physeq_LBK_S04, "PCoA","bray")

plot_ordination(physeq_LBK_S04, PCoA_LBK_S04, type = "sites", color = "restriction") +
  geom_point(size = 4, aes(shape = restriction)) + theme_bw() + scale_shape_manual(values = c('fishing prohibited' = 21,
                                                                                              'fishing gear restriction' = 22,
                                                                                              'unrestricted fishing' = 24)) +
  facet_grid(region~subset, scales = "free") + stat_ellipse(mapping = aes())

#Extract coordinates from ordination
PCoA_LBK_S04_scores <- as.data.frame(PCoA_LBK_S04[["vectors"]])
PCoA_LBK_S04_scores <- PCoA_LBK_S04_scores[1:2]
names(PCoA_LBK_S04_scores)[1] <- "PC1"
names(PCoA_LBK_S04_scores)[2] <- "PC2"

#Merge coordinates with metadata
PCoA_LBK_S04_scores <- merge(sample_data, PCoA_LBK_S04_scores, by = "row.names", all = FALSE, no.dups = TRUE)

#Plot formatted ordination
PCoA_LBK_S04_plot <- ggplot(PCoA_LBK_S04_scores, aes(PC1, PC2, fill = restriction, shape = restriction)) + 
  geom_point(size = 4, alpha = 0.75, color = "black") +
  stat_ellipse(aes(color = restriction)) + theme_bw() + theme(panel.grid = element_blank()) +
  scale_shape_manual(values = c('fishing prohibited' = 21,
                                'fishing gear restriction' = 22,
                                'unrestricted fishing' = 24)) +
  scale_fill_manual(values = c('fishing prohibited' = "#f73e42",
                               'fishing gear restriction' = "#8e86c5",
                               'unrestricted fishing' = "#007cb5")) +
  scale_color_manual(values = c('fishing prohibited' = "#f73e42",
                                'fishing gear restriction' = "#8e86c5",
                                'unrestricted fishing' = "#007cb5")) +
  xlab("PC1 [14.9%]") + ylab("PC2 [8.7%]") +
  facet_grid(region ~ subset, scales = "free") + 
  theme(strip.text.x = element_text(face = "bold", color = "white"), 
        strip.background.x = element_rect(fill = "black"), 
        strip.background.y = element_blank(), 
        strip.text.y = element_blank()) 
PCoA_LBK_S04_plot

#Run and plot PCoA ordination
PCoA_LBK_S12 <- ordinate(physeq_LBK_S12, "PCoA","bray")

plot_ordination(physeq_LBK_S12, PCoA_LBK_S12, type = "sites", color = "restriction") +
  geom_point(size = 4, aes(shape = restriction)) + theme_bw() + scale_shape_manual(values = c('fishing prohibited' = 21,
                                                                                              'fishing gear restriction' = 22,
                                                                                              'unrestricted fishing' = 24)) +
  facet_grid(region~subset, scales = "free") + stat_ellipse(mapping = aes())

#Extract coordinates from ordination
PCoA_LBK_S12_scores <- as.data.frame(PCoA_LBK_S12[["vectors"]])
PCoA_LBK_S12_scores <- PCoA_LBK_S12_scores[1:2]
names(PCoA_LBK_S12_scores)[1] <- "PC1"
names(PCoA_LBK_S12_scores)[2] <- "PC2"

#Merge coordinates with metadata
PCoA_LBK_S12_scores <- merge(sample_data, PCoA_LBK_S12_scores, by = "row.names", all = FALSE, no.dups = TRUE)

#Plot formatted ordination
PCoA_LBK_S12_plot <- ggplot(PCoA_LBK_S12_scores, aes(PC1, PC2, fill = restriction, shape = restriction)) + 
  geom_point(size = 4, alpha = 0.75, color = "black") +
  stat_ellipse(aes(color = restriction)) + theme_bw() + theme(panel.grid = element_blank()) +
  scale_shape_manual(values = c('fishing prohibited' = 21,
                                'fishing gear restriction' = 22,
                                'unrestricted fishing' = 24)) +
  scale_fill_manual(values = c('fishing prohibited' = "#f73e42",
                               'fishing gear restriction' = "#8e86c5",
                               'unrestricted fishing' = "#007cb5")) +
  scale_color_manual(values = c('fishing prohibited' = "#f73e42",
                                'fishing gear restriction' = "#8e86c5",
                                'unrestricted fishing' = "#007cb5")) +
  xlab("PC1 [14.5%]") + ylab("PC2 [8.8%]") +
  facet_grid(region ~ subset, scales = "free") + 
  theme(strip.text = element_text(face = "bold", color = "white"), 
        strip.background = element_rect(fill = "black")) 
PCoA_LBK_S12_plot

#Run and plot PCoA ordination
PCoA_MIS_M04 <- ordinate(physeq_MIS_M04, "PCoA","bray")

plot_ordination(physeq_MIS_M04, PCoA_MIS_M04, type = "sites", color = "restriction") +
  geom_point(size = 4, aes(shape = restriction)) + theme_bw() + scale_shape_manual(values = c('fishing prohibited' = 21,
                                                                                              'fishing gear restriction' = 22,
                                                                                              'unrestricted fishing' = 24)) +
  facet_grid(region~subset, scales = "free") + stat_ellipse(mapping = aes())

#Extract coordinates from ordination
PCoA_MIS_M04_scores <- as.data.frame(PCoA_MIS_M04[["vectors"]])
PCoA_MIS_M04_scores <- PCoA_MIS_M04_scores[1:2]
names(PCoA_MIS_M04_scores)[1] <- "PC1"
names(PCoA_MIS_M04_scores)[2] <- "PC2"

#Merge coordinates with metadata
PCoA_MIS_M04_scores <- merge(sample_data, PCoA_MIS_M04_scores, by = "row.names", all = FALSE, no.dups = TRUE)

#Plot formatted ordination
PCoA_MIS_M04_plot <- ggplot(PCoA_MIS_M04_scores, aes(PC1, PC2, fill = restriction, shape = restriction)) + 
  geom_point(size = 4, alpha = 0.75, color = "black") +
  stat_ellipse(aes(color = restriction)) + theme_bw() + theme(panel.grid = element_blank()) +
  scale_shape_manual(values = c('fishing prohibited' = 21,
                                'fishing gear restriction' = 22,
                                'unrestricted fishing' = 24)) +
  scale_fill_manual(values = c('fishing prohibited' = "#f73e42",
                               'fishing gear restriction' = "#8e86c5",
                               'unrestricted fishing' = "#007cb5")) +
  scale_color_manual(values = c('fishing prohibited' = "#f73e42",
                                'fishing gear restriction' = "#8e86c5",
                                'unrestricted fishing' = "#007cb5")) +
  xlab("PC1 [24.1%]") + ylab("PC2 [13.1%]") +
  facet_grid(region ~ subset, scales = "free") + 
  theme(strip.text.x = element_blank(), 
        strip.background.x = element_blank(), 
        strip.background.y = element_blank(), 
        strip.text.y = element_blank()) 
PCoA_MIS_M04_plot

#Run and plot PCoA ordination
PCoA_MIS_M12 <- ordinate(physeq_MIS_M12, "PCoA","bray")

plot_ordination(physeq_MIS_M12, PCoA_MIS_M12, type = "sites", color = "restriction") +
  geom_point(size = 4, aes(shape = restriction)) + theme_bw() + scale_shape_manual(values = c('fishing prohibited' = 21,
                                                                                              'fishing gear restriction' = 22,
                                                                                              'unrestricted fishing' = 24)) +
  facet_grid(region~subset, scales = "free") + stat_ellipse(mapping = aes())

#Extract coordinates from ordination
PCoA_MIS_M12_scores <- as.data.frame(PCoA_MIS_M12[["vectors"]])
PCoA_MIS_M12_scores <- PCoA_MIS_M12_scores[1:2]
names(PCoA_MIS_M12_scores)[1] <- "PC1"
names(PCoA_MIS_M12_scores)[2] <- "PC2"

#Merge coordinates with metadata
PCoA_MIS_M12_scores <- merge(sample_data, PCoA_MIS_M12_scores, by = "row.names", all = FALSE, no.dups = TRUE)

#Plot formatted ordination
PCoA_MIS_M12_plot <- ggplot(PCoA_MIS_M12_scores, aes(PC1, PC2, fill = restriction, shape = restriction)) + 
  geom_point(size = 4, alpha = 0.75, color = "black") +
  stat_ellipse(aes(color = restriction)) + theme_bw() + theme(panel.grid = element_blank()) +
  scale_shape_manual(values = c('fishing prohibited' = 21,
                                'fishing gear restriction' = 22,
                                'unrestricted fishing' = 24)) +
  scale_fill_manual(values = c('fishing prohibited' = "#f73e42",
                               'fishing gear restriction' = "#8e86c5",
                               'unrestricted fishing' = "#007cb5")) +
  scale_color_manual(values = c('fishing prohibited' = "#f73e42",
                                'fishing gear restriction' = "#8e86c5",
                                'unrestricted fishing' = "#007cb5")) +
  xlab("PC1 [26.9%]") + ylab("PC2 [9.4%]") +
  facet_grid(region ~ subset, scales = "free") + 
  theme(strip.text.x = element_blank(), 
        strip.background.x = element_blank(), 
        strip.background.y = element_blank(), 
        strip.text.y = element_blank()) 
PCoA_MIS_M12_plot

#Run and plot PCoA ordination
PCoA_MIS_S04 <- ordinate(physeq_MIS_S04, "PCoA","bray")

plot_ordination(physeq_MIS_S04, PCoA_MIS_S04, type = "sites", color = "restriction") +
  geom_point(size = 4, aes(shape = restriction)) + theme_bw() + scale_shape_manual(values = c('fishing prohibited' = 21,
                                                                                              'fishing gear restriction' = 22,
                                                                                              'unrestricted fishing' = 24)) +
  facet_grid(region~subset, scales = "free") + stat_ellipse(mapping = aes())

#Extract coordinates from ordination
PCoA_MIS_S04_scores <- as.data.frame(PCoA_MIS_S04[["vectors"]])
PCoA_MIS_S04_scores <- PCoA_MIS_S04_scores[1:2]
names(PCoA_MIS_S04_scores)[1] <- "PC1"
names(PCoA_MIS_S04_scores)[2] <- "PC2"

#Merge coordinates with metadata
PCoA_MIS_S04_scores <- merge(sample_data, PCoA_MIS_S04_scores, by = "row.names", all = FALSE, no.dups = TRUE)

#Plot formatted ordination
PCoA_MIS_S04_plot <- ggplot(PCoA_MIS_S04_scores, aes(PC1, PC2, fill = restriction, shape = restriction)) + 
  geom_point(size = 4, alpha = 0.75, color = "black") +
  stat_ellipse(aes(color = restriction)) + theme_bw() + theme(panel.grid = element_blank()) +
  scale_shape_manual(values = c('fishing prohibited' = 21,
                                'fishing gear restriction' = 22,
                                'unrestricted fishing' = 24)) +
  scale_fill_manual(values = c('fishing prohibited' = "#f73e42",
                               'fishing gear restriction' = "#8e86c5",
                               'unrestricted fishing' = "#007cb5")) +
  scale_color_manual(values = c('fishing prohibited' = "#f73e42",
                                'fishing gear restriction' = "#8e86c5",
                                'unrestricted fishing' = "#007cb5")) +
  xlab("PC1 [12.9%]") + ylab("PC2 [8.5%]") +
  facet_grid(region ~ subset, scales = "free") + 
  theme(strip.text.x = element_blank(), 
        strip.background.x = element_blank(), 
        strip.background.y = element_blank(), 
        strip.text.y = element_blank()) 
PCoA_MIS_S04_plot

#Run and plot PCoA ordination
PCoA_MIS_S12 <- ordinate(physeq_MIS_S12, "PCoA","bray")

plot_ordination(physeq_MIS_S12, PCoA_MIS_S12, type = "sites", color = "restriction") +
  geom_point(size = 4, aes(shape = restriction)) + theme_bw() + scale_shape_manual(values = c('fishing prohibited' = 21,
                                                                                              'fishing gear restriction' = 22,
                                                                                              'unrestricted fishing' = 24)) +
  facet_grid(region~subset, scales = "free") + stat_ellipse(mapping = aes())

#Extract coordinates from ordination
PCoA_MIS_S12_scores <- as.data.frame(PCoA_MIS_S12[["vectors"]])
PCoA_MIS_S12_scores <- PCoA_MIS_S12_scores[1:2]
names(PCoA_MIS_S12_scores)[1] <- "PC1"
names(PCoA_MIS_S12_scores)[2] <- "PC2"

#Merge coordinates with metadata
PCoA_MIS_S12_scores <- merge(sample_data, PCoA_MIS_S12_scores, by = "row.names", all = FALSE, no.dups = TRUE)

#Plot formatted ordination
PCoA_MIS_S12_plot <- ggplot(PCoA_MIS_S12_scores, aes(PC1, PC2, fill = restriction, shape = restriction)) + 
  geom_point(size = 4, alpha = 0.75, color = "black") +
  stat_ellipse(aes(color = restriction)) + theme_bw() + theme(panel.grid = element_blank()) +
  scale_shape_manual(values = c('fishing prohibited' = 21,
                                'fishing gear restriction' = 22,
                                'unrestricted fishing' = 24)) +
  scale_fill_manual(values = c('fishing prohibited' = "#f73e42",
                               'fishing gear restriction' = "#8e86c5",
                               'unrestricted fishing' = "#007cb5")) +
  scale_color_manual(values = c('fishing prohibited' = "#f73e42",
                                'fishing gear restriction' = "#8e86c5",
                                'unrestricted fishing' = "#007cb5")) +
  xlab("PC1 [9.6%]") + ylab("PC2 [8.6%]") +
  facet_grid(region ~ subset, scales = "free") + 
  theme(strip.text.y = element_text(face = "bold", color = "white"), 
        strip.background.y = element_rect(fill = "black"),
        strip.text.x = element_blank(),
        strip.background.x = element_blank()) 
PCoA_MIS_S12_plot

#Run and plot PCoA ordination
PCoA_WGO_M04 <- ordinate(physeq_WGO_M04, "PCoA","bray")

plot_ordination(physeq_WGO_M04, PCoA_WGO_M04, type = "sites", color = "restriction") +
  geom_point(size = 4, aes(shape = restriction)) + theme_bw() + scale_shape_manual(values = c('fishing prohibited' = 21,
                                                                                              'fishing gear restriction' = 22,
                                                                                              'unrestricted fishing' = 24)) +
  facet_grid(region~subset, scales = "free") + stat_ellipse(mapping = aes())

#Extract coordinates from ordination
PCoA_WGO_M04_scores <- as.data.frame(PCoA_WGO_M04[["vectors"]])
PCoA_WGO_M04_scores <- PCoA_WGO_M04_scores[1:2]
names(PCoA_WGO_M04_scores)[1] <- "PC1"
names(PCoA_WGO_M04_scores)[2] <- "PC2"

#Merge coordinates with metadata
PCoA_WGO_M04_scores <- merge(sample_data, PCoA_WGO_M04_scores, by = "row.names", all = FALSE, no.dups = TRUE)

#Plot formatted ordination
PCoA_WGO_M04_plot <- ggplot(PCoA_WGO_M04_scores, aes(PC1, PC2, fill = restriction, shape = restriction)) + 
  geom_point(size = 4, alpha = 0.75, color = "black") +
  stat_ellipse(aes(color = restriction)) + theme_bw() + theme(panel.grid = element_blank()) +
  scale_shape_manual(values = c('fishing prohibited' = 21,
                                'fishing gear restriction' = 22,
                                'unrestricted fishing' = 24)) +
  scale_fill_manual(values = c('fishing prohibited' = "#f73e42",
                               'fishing gear restriction' = "#8e86c5",
                               'unrestricted fishing' = "#007cb5")) +
  scale_color_manual(values = c('fishing prohibited' = "#f73e42",
                                'fishing gear restriction' = "#8e86c5",
                                'unrestricted fishing' = "#007cb5")) +
  xlab("PC1 [14.6%]") + ylab("PC2 [9.3%]") +
  facet_grid(region ~ subset, scales = "free") + 
  theme(strip.text.x = element_blank(), 
        strip.background.x = element_blank(), 
        strip.background.y = element_blank(), 
        strip.text.y = element_blank()) 
PCoA_WGO_M04_plot

#Run and plot PCoA ordination
PCoA_WGO_M12 <- ordinate(physeq_WGO_M12, "PCoA","bray")

plot_ordination(physeq_WGO_M12, PCoA_WGO_M12, type = "sites", color = "restriction") +
  geom_point(size = 4, aes(shape = restriction)) + theme_bw() + scale_shape_manual(values = c('fishing prohibited' = 21,
                                                                                              'fishing gear restriction' = 22,
                                                                                              'unrestricted fishing' = 24)) +
  facet_grid(region~subset, scales = "free") + stat_ellipse(mapping = aes())

#Extract coordinates from ordination
PCoA_WGO_M12_scores <- as.data.frame(PCoA_WGO_M12[["vectors"]])
PCoA_WGO_M12_scores <- PCoA_WGO_M12_scores[1:2]
names(PCoA_WGO_M12_scores)[1] <- "PC1"
names(PCoA_WGO_M12_scores)[2] <- "PC2"

#Merge coordinates with metadata
PCoA_WGO_M12_scores <- merge(sample_data, PCoA_WGO_M12_scores, by = "row.names", all = FALSE, no.dups = TRUE)

#Plot formatted ordination
PCoA_WGO_M12_plot <- ggplot(PCoA_WGO_M12_scores, aes(PC1, PC2, fill = restriction, shape = restriction)) + 
  geom_point(size = 4, alpha = 0.75, color = "black") +
  stat_ellipse(aes(color = restriction)) + theme_bw() + theme(panel.grid = element_blank()) +
  scale_shape_manual(values = c('fishing prohibited' = 21,
                                'fishing gear restriction' = 22,
                                'unrestricted fishing' = 24)) +
  scale_fill_manual(values = c('fishing prohibited' = "#f73e42",
                               'fishing gear restriction' = "#8e86c5",
                               'unrestricted fishing' = "#007cb5")) +
  scale_color_manual(values = c('fishing prohibited' = "#f73e42",
                                'fishing gear restriction' = "#8e86c5",
                                'unrestricted fishing' = "#007cb5")) +
  xlab("PC1 [13.9%]") + ylab("PC2 [10.2%]") +
  facet_grid(region ~ subset, scales = "free") + 
  theme(strip.text.x = element_blank(), 
        strip.background.x = element_blank(), 
        strip.background.y = element_blank(), 
        strip.text.y = element_blank()) 
PCoA_WGO_M12_plot

#Run and plot PCoA ordination
PCoA_WGO_S04 <- ordinate(physeq_WGO_S04, "PCoA","bray")

plot_ordination(physeq_WGO_S04, PCoA_WGO_S04, type = "sites", color = "restriction") +
  geom_point(size = 4, aes(shape = restriction)) + theme_bw() + scale_shape_manual(values = c('fishing prohibited' = 21,
                                                                                              'fishing gear restriction' = 22,
                                                                                              'unrestricted fishing' = 24)) +
  facet_grid(region~subset, scales = "free") + stat_ellipse(mapping = aes())

#Extract coordinates from ordination
PCoA_WGO_S04_scores <- as.data.frame(PCoA_WGO_S04[["vectors"]])
PCoA_WGO_S04_scores <- PCoA_WGO_S04_scores[1:2]
names(PCoA_WGO_S04_scores)[1] <- "PC1"
names(PCoA_WGO_S04_scores)[2] <- "PC2"

#Merge coordinates with metadata
PCoA_WGO_S04_scores <- merge(sample_data, PCoA_WGO_S04_scores, by = "row.names", all = FALSE, no.dups = TRUE)

#Plot formatted ordination
PCoA_WGO_S04_plot <- ggplot(PCoA_WGO_S04_scores, aes(PC1, PC2, fill = restriction, shape = restriction)) + 
  geom_point(size = 4, alpha = 0.75, color = "black") +
  stat_ellipse(aes(color = restriction)) + theme_bw() + theme(panel.grid = element_blank()) +
  scale_shape_manual(values = c('fishing prohibited' = 21,
                                'fishing gear restriction' = 22,
                                'unrestricted fishing' = 24)) +
  scale_fill_manual(values = c('fishing prohibited' = "#f73e42",
                               'fishing gear restriction' = "#8e86c5",
                               'unrestricted fishing' = "#007cb5")) +
  scale_color_manual(values = c('fishing prohibited' = "#f73e42",
                                'fishing gear restriction' = "#8e86c5",
                                'unrestricted fishing' = "#007cb5")) +
  xlab("PC1 [11.8%]") + ylab("PC2 [9.4%]") +
  facet_grid(region ~ subset, scales = "free") + 
  theme(strip.text.x = element_blank(), 
        strip.background.x = element_blank(), 
        strip.background.y = element_blank(), 
        strip.text.y = element_blank()) 
PCoA_WGO_S04_plot

#Run and plot PCoA ordination
PCoA_WGO_S12 <- ordinate(physeq_WGO_S12, "PCoA","bray")

plot_ordination(physeq_WGO_S12, PCoA_WGO_S12, type = "sites", color = "restriction") +
  geom_point(size = 4, aes(shape = restriction)) + theme_bw() + scale_shape_manual(values = c('fishing prohibited' = 21,
                                                                                              'fishing gear restriction' = 22,
                                                                                              'unrestricted fishing' = 24)) +
  facet_grid(region~subset, scales = "free") + stat_ellipse(mapping = aes())

#Extract coordinates from ordination
PCoA_WGO_S12_scores <- as.data.frame(PCoA_WGO_S12[["vectors"]])
PCoA_WGO_S12_scores <- PCoA_WGO_S12_scores[1:2]
names(PCoA_WGO_S12_scores)[1] <- "PC1"
names(PCoA_WGO_S12_scores)[2] <- "PC2"

#Merge coordinates with metadata
PCoA_WGO_S12_scores <- merge(sample_data, PCoA_WGO_S12_scores, by = "row.names", all = FALSE, no.dups = TRUE)

#Plot formatted ordination
PCoA_WGO_S12_plot <- ggplot(PCoA_WGO_S12_scores, aes(PC1, PC2, fill = restriction, shape = restriction)) + 
  geom_point(size = 4, alpha = 0.75, color = "black") +
  stat_ellipse(aes(color = restriction)) + theme_bw() + theme(panel.grid = element_blank()) +
  scale_shape_manual(values = c('fishing prohibited' = 21,
                                'fishing gear restriction' = 22,
                                'unrestricted fishing' = 24)) +
  scale_fill_manual(values = c('fishing prohibited' = "#f73e42",
                               'fishing gear restriction' = "#8e86c5",
                               'unrestricted fishing' = "#007cb5")) +
  scale_color_manual(values = c('fishing prohibited' = "#f73e42",
                                'fishing gear restriction' = "#8e86c5",
                                'unrestricted fishing' = "#007cb5")) +
  xlab("PC1 [10.6%]") + ylab("PC2 [8.9%]") +
  facet_grid(region ~ subset, scales = "free") + 
  theme(strip.text.y = element_text(face = "bold", color = "white"), 
        strip.background.y = element_rect(fill = "black"),
        strip.text.x = element_blank(),
        strip.background.x = element_blank()) 
PCoA_WGO_S12_plot

#Combine plots into a single grid and save plot
plot <- ggarrange(PCoA_LBK_M04_plot, PCoA_LBK_M12_plot, PCoA_LBK_S04_plot, PCoA_LBK_S12_plot,
                  PCoA_MIS_M04_plot, PCoA_MIS_S04_plot, PCoA_MIS_S04_plot, PCoA_MIS_S12_plot,
                  PCoA_WGO_M04_plot, PCoA_WGO_M12_plot, PCoA_WGO_S04_plot, PCoA_WGO_S12_plot,
                  nrow = 3, ncol = 4, align = "hv", legend = "bottom", common.legend = TRUE)
plot
ggsave("PCoA_by_restriction.pdf", plot, height = 8, width = 12)
