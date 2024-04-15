#load packages
library(ggplot2)
library(ggpubr)
library(ggvenn)
library(reshape2)

#set working directory
setwd("~/Desktop/JPhyc_sub/")

#import data
asv_table <- read.csv("~/Desktop/DataAnalysis/threshold95/SAR/SAR-percent-table-95.csv", row.names = 1)
asv_table <- asv_table[,-c(1:8)]
asv_transformed <- as.data.frame(t(asv_table))

otu_table <- read.csv("~/Desktop/DataAnalysis/threshold95/SAR/SAR_OTU_percent.csv", row.names = 1)
otu_transformed <- as.data.frame(t(otu_table))

metadata <- read.csv("~/Desktop/DataAnalysis/Master_V9_MappingFile.csv", row.names = 1)

#create separate metadata files for each region
LBK_metadata <- subset(metadata, region == "Lombok")
WK_metadata <- subset(metadata, region == "Wakatobi")
MIS_metadata <- subset(metadata, region == "Misool")
WGO_metadata <- subset(metadata, region == "Waigeo")

#create vectors for each region for ASVs present from ASV table
LBK_table <- subset(asv_transformed, rownames(asv_transformed) %in% rownames(LBK_metadata))
LBK_table <- as.data.frame(t(LBK_table))
LBK_table$sums <- rowSums(LBK_table)
LBK_filtered <- subset(LBK_table, sums > 0)
LBK_ASVs <- rownames(LBK_filtered)

WK_table <- subset(asv_transformed, rownames(asv_transformed) %in% rownames(WK_metadata))
WK_table <- as.data.frame(t(WK_table))
WK_table$sums <- rowSums(WK_table)
WK_filtered <- subset(WK_table, sums > 0)
WK_ASVs <- rownames(WK_filtered)

MIS_table <- subset(asv_transformed, rownames(asv_transformed) %in% rownames(MIS_metadata))
MIS_table <- as.data.frame(t(MIS_table))
MIS_table$sums <- rowSums(MIS_table)
MIS_filtered <- subset(MIS_table, sums > 0)
MIS_ASVs <- rownames(MIS_filtered)

WGO_table <- subset(asv_transformed, rownames(asv_transformed) %in% rownames(WGO_metadata))
WGO_table <- as.data.frame(t(WGO_table))
WGO_table$sums <- rowSums(WGO_table)
WGO_filtered <- subset(WGO_table, sums > 0)
WGO_ASVs <- rownames(WGO_filtered)

#build venn diagram
venn_asv <- list(Lombok = LBK_ASVs,
             Wakatobi = WK_ASVs,
             Misool = MIS_ASVs,
             Waigeo = WGO_ASVs)

venn_asvPlot <- ggvenn(venn_asv, c("Lombok","Wakatobi","Misool","Waigeo"), fill_color = c("#8e86c5",
                                                                                          "#007cb5",
                                                                                          "#a6c54b",
                                                                                          "#f73e42")) 
venn_asvPlot

#Run same analysis on OTUs
#create vectors for each region for ASVs present from ASV table
LBK_otu_table <- subset(otu_transformed, rownames(otu_transformed) %in% rownames(LBK_metadata))
LBK_otu_table <- as.data.frame(t(LBK_otu_table))
LBK_otu_table$sums <- rowSums(LBK_otu_table)
LBK_otu_filtered <- subset(LBK_otu_table, sums > 0)
LBK_OTUs <- rownames(LBK_otu_filtered)

WK_otu_table <- subset(otu_transformed, rownames(otu_transformed) %in% rownames(WK_metadata))
WK_otu_table <- as.data.frame(t(WK_otu_table))
WK_otu_table$sums <- rowSums(WK_otu_table)
WK_otu_filtered <- subset(WK_otu_table, sums > 0)
WK_OTUs <- rownames(WK_otu_filtered)

MIS_otu_table <- subset(otu_transformed, rownames(otu_transformed) %in% rownames(MIS_metadata))
MIS_otu_table <- as.data.frame(t(MIS_otu_table))
MIS_otu_table$sums <- rowSums(MIS_otu_table)
MIS_otu_filtered <- subset(MIS_otu_table, sums > 0)
MIS_OTUs <- rownames(MIS_otu_filtered)

WGO_otu_table <- subset(otu_transformed, rownames(otu_transformed) %in% rownames(WGO_metadata))
WGO_otu_table <- as.data.frame(t(WGO_otu_table))
WGO_otu_table$sums <- rowSums(WGO_otu_table)
WGO_otu_filtered <- subset(WGO_otu_table, sums > 0)
WGO_OTUs <- rownames(WGO_otu_filtered)

#build venn diagram
venn_otu <- list(Lombok = LBK_OTUs,
                 Wakatobi = WK_OTUs,
                 Misool = MIS_OTUs,
                 Waigeo = WGO_OTUs)

venn_otuPlot <- ggvenn(venn_otu, c("Lombok","Wakatobi","Misool","Waigeo"), fill_color = c("#8e86c5",
                                                                                          "#007cb5",
                                                                                          "#a6c54b",
                                                                                          "#f73e42"))
venn_otuPlot

#format dataframe of venn output to plot barplot of results
venn_data <- data.frame(row.names = c("LBK","WK","MIS","WGO","LWM","LWW","LMW",
                                      "WMW","LWk","LM","LWg","WM","WW","MW","All"))
venn_data$regions <- c("LBK","WK","MIS","WGO","LWM","LWW","LMW",
                       "WMW","LWk","LM","LWg","WM","WW","MW","All")
venn_data$venn <- c(2562,1385,2534,2021,220,99,208,125,238,394,146,186,96,213,681)
venn_data$total <- 11108 - venn_data$venn
venn_long <- melt(venn_data)
venn_long$variable <- factor(venn_long$variable, levels = c("total", "venn"))
venn_long$regions <- factor(venn_long$regions, levels = c(c("All","LBK","WK","MIS","WGO","LMW","LWM","LWW",
                                                            "WMW","MW","LM","LWg","LWk","WM","WW")))

#plot barplot
barplot <- ggplot(venn_data, aes(reorder(rownames(venn_data), -venn), venn)) + geom_bar(stat = "identity", fill = "grey20") +
  theme_bw() + theme(panel.grid = element_blank()) + theme(axis.title.x = element_blank(), axis.text.x = element_blank()) +
  ylab("Number of unique ASVs per \nregion or region groupings") + scale_y_continuous(expand = c(0,0), limits = c(0,2600))
barplot

stacked_barplot <- ggplot(venn_long, aes(regions, value, fill = variable)) + 
  geom_bar(stat = "identity", position = "fill") + theme_bw() + theme(panel.grid = element_blank()) +
  ylab("Percent ASVs in dataset unique to each \nregion/region grouping") + scale_y_continuous(expand = c(0,0)) +
  scale_fill_manual(values = c("grey75","grey20")) + theme(axis.title.x = element_blank(), axis.text.x = element_blank()) +
  theme(legend.position = "none")
stacked_barplot

# build data frame to plot region bubble plot from
region_data <- data.frame(row.names = c("LBK","MIS","WGO","WK","All","LM","LWk","LWM","MW","LMW",
                                        "WM","LWg","WMW","LWW","WW"))
region_data$regions <- c("LBK","MIS","WGO","WK","All","LM","LWk","LWM","MW","LMW",
                         "WM","LWg","WMW","LWW","WW")
region_data$Lombok <-   c(1,0,0,0,1,1,1,1,0,1,0,1,0,1,0)
region_data$Wakatobi <- c(0,0,0,1,1,0,1,1,0,0,1,0,1,1,1)
region_data$Misool <-   c(0,1,0,0,1,1,0,1,1,1,1,0,1,0,0)
region_data$Waigeo <-   c(0,0,1,0,1,0,0,0,1,1,0,1,1,1,1)

region_data_long <- melt(region_data)
region_data_long <- subset(region_data_long, value > 0)
region_data_long$variable <- factor(region_data_long$variable, levels = c("Waigeo", "Misool","Wakatobi","Lombok"))
region_data_long$regions <- factor(region_data_long$regions, levels = c("All","LBK","WK","MIS","WGO","LMW","LWM","LWW",
                                                                          "WMW","MW","LM","LWg","LWk","WM","WW"))

region_plot <- ggplot(region_data_long, aes(regions, variable)) + geom_point(size = 3) + 
  theme_bw() + theme(axis.title = element_blank(), axis.text.x = element_blank(), 
                     axis.ticks.x = element_blank(), axis.text.y = element_text(size = 12, face = "bold"))
region_plot

#filtering ASV table to build data frame with percent of ASVs in overall dataset
asv_counts <- read.csv("~/Desktop/DataAnalysis/threshold95/SAR/SAR-count-table-95.csv", row.names = 1)

LBK_only <- subset(asv_counts, rownames(asv_counts) %in% attributes(list_by_regions)$intersections$Lombok)
WK_only <- subset(asv_counts, rownames(asv_counts) %in% attributes(list_by_regions)$intersections$Wakatobi)
MIS_only <- subset(asv_counts, rownames(asv_counts) %in% attributes(list_by_regions)$intersections$Misool)
WGO_only <- subset(asv_counts, rownames(asv_counts) %in% attributes(list_by_regions)$intersections$Waigeo)

shared <- subset(asv_counts, rownames(asv_counts) %in% attributes(list_by_regions)$intersections$`Lombok:Wakatobi:Misool:Waigeo`)
LWM <- subset(asv_counts, rownames(asv_counts) %in% attributes(list_by_regions)$intersections$`Lombok:Wakatobi:Misool`)
LMW <- subset(asv_counts, rownames(asv_counts) %in% attributes(list_by_regions)$intersections$`Lombok:Misool:Waigeo`)
WMW <- subset(asv_counts, rownames(asv_counts) %in% attributes(list_by_regions)$intersections$`Wakatobi:Misool:Waigeo`)
LWW <- subset(asv_counts, rownames(asv_counts) %in% attributes(list_by_regions)$intersections$`Lombok:Wakatobi:Waigeo`)

LWk <- subset(asv_counts, rownames(asv_counts) %in% attributes(list_by_regions)$intersections$`Lombok:Wakatobi`)
LM <- subset(asv_counts, rownames(asv_counts) %in% attributes(list_by_regions)$intersections$`Lombok:Misool`)
LWg <- subset(asv_counts, rownames(asv_counts) %in% attributes(list_by_regions)$intersections$`Lombok:Waigeo`)
WM <- subset(asv_counts, rownames(asv_counts) %in% attributes(list_by_regions)$intersections$`Wakatobi:Misool`)
WW <- subset(asv_counts, rownames(asv_counts) %in% attributes(list_by_regions)$intersections$`Wakatobi:Waigeo`)
MW <- subset(asv_counts, rownames(asv_counts) %in% attributes(list_by_regions)$intersections$`Misool:Waigeo`)

sum(asv_counts)
sum(LBK_only)
sum(WK_only)
sum(MIS_only)
sum(WGO_only)

sum(shared)
sum(LWM)
sum(LMW)
sum(WMW)
sum(LWW)

sum(LWk)
sum(LM)
sum(LWg)
sum(WM)
sum(WW)
sum(MW)

#assemble data frame with counts for each region grouping
count_by_region <- data.frame(row.names = c("LBK","MIS","WGO","WK","All","LM","LWk","LWM","MW","LMW",
                                            "WM","LWg","WMW","LWW","WW"))
count_by_region$regions <- c("LBK","MIS","WGO","WK","All","LM","LWk","LWM","MW","LMW",
                                            "WM","LWg","WMW","LWW","WW")
count_by_region$count <- c(101650,95870,210259,47638,2561608,51703,32883,97987,60328,167930,
                           21253,34029,78536,94879,15036)
count_by_region$total <- 3671589 - count_by_region$count
count_by_region_long <- melt(count_by_region)

count_by_region_long$variable <- factor(count_by_region_long$variable, levels = c("total", "count"))
count_by_region_long$regions <- factor(count_by_region_long$regions, levels = c("All","LBK","WK","MIS","WGO","LMW","LWM","LWW",
                                                                                "WMW","MW","LM","LWg","LWk","WM","WW"))

#plot percent reads out of whole dataset
read_count_barplot <- ggplot(count_by_region_long, aes(regions, value, fill = variable)) +
  geom_bar(stat = "identity", position = "fill") + scale_y_continuous(expand = c(0,0)) + 
  theme_bw() + theme(panel.grid = element_blank()) + scale_fill_manual(values = c("grey75","grey20")) +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank()) + theme(legend.position = "none") +
  ylab("Percent reads in total dataset unique to \nregion/region groupings")
read_count_barplot

#read ASV table with taxonomy
asv_table_with_tax <- read.csv("~/Desktop/DataAnalysis/threshold95/SAR/SAR-percent-table-95.csv", row.names = 1)
asv_table_with_tax <- asv_table_with_tax[,1:8]
asv_counts_with_tax <- merge(asv_table_with_tax, asv_counts, by = "row.names")
asv_counts_with_tax$sum <- rowSums(asv_counts)
asv_counts_with_tax <- asv_counts_with_tax[,c(1:9,315)]

phyla <- asv_table_with_tax[,3:4]
tax_nodups <- phyla %>% distinct(L2, L3, .keep_all = TRUE)

shared <- subset(asv_counts_with_tax, Row.names %in% rownames(shared))
shared_count <- shared %>% count(L1,L2,L3)

LBK <- subset(asv_counts_with_tax, Row.names %in% rownames(LBK_only))
LBK_count <- LBK_treemap %>% count(L1,L2,L3)

WK <- subset(asv_counts_with_tax, Row.names %in% rownames(WK_only))
WK_count <- WK_treemap %>% count(L1,L2,L3)

MIS <- subset(asv_counts_with_tax, Row.names %in% rownames(MIS_only))
MIS_count <- MIS_treemap %>% count(L1,L2,L3)

WGO <- subset(asv_counts_with_tax, Row.names %in% rownames(WGO_only))
WGO_count <- WGO_treemap %>% count(L1,L2,L3)

#stacked barplot of asv counts unique or shared amongst regions
shared_count$region <- "shared"
LBK_count$region <- "LBK"
WK_count$region <- "WK"
MIS_count$region <- "MIS"
WGO_count$region <- "WGO"

merged_data <- rbind(shared_count,LBK_count,WK_count,MIS_count,WGO_count)


#write csv to fix color categories then reimport file
write.csv(merged_data, "merged_data.csv")
merged_data <- read.csv("merged_data.csv", row.names = 1)
unique(merged_data$color)
merged_data$region <- factor(merged_data$region, levels = c("WGO", "MIS", "WK", "LBK", "shared"))
merged_data$color <- factor(merged_data$color, levels = c("Alveolata undefined",
                                                          "Apicomplexa",
                                                          "Ciliophora",
                                                          "Dinoflagellata undefined",
                                                          "Dinophyceae",
                                                          "Ellobiophyceae",
                                                          "Noctilucophyceae",
                                                          "Syndiniales",
                                                          "Cercozoa",
                                                          "Foraminifera",
                                                          "Radiolaria",
                                                          "Stramenopiles undefined",
                                                          "Ochrophyta undefined",
                                                          "Bacillariophyta",
                                                          "Bolidophyceae",
                                                          "Chrysophyceae",
                                                          "Dictyochophyceae",
                                                          "MOCH",
                                                          "Pelagophyceae",
                                                          "Phaeophyceae",
                                                          "Raphidophyceae",
                                                          "Opalozoa",
                                                          "Pseudofungi",
                                                          "Sagenista",
                                                          "Perkinsea",
                                                          "Eustigmatophyceae",
                                                          "Pinguiophyceae",
                                                          "Synurophyceae"))

stacked_barplot_regions <- ggplot(merged_data, aes(region, n, fill = color)) + 
  geom_bar(stat = "identity", position = "fill") + coord_flip() + theme_bw() + theme(panel.grid = element_blank()) + 
  theme(legend.position = "bottom", legend.spacing = unit(0.01, "cm"), legend.key.size = unit(0.5, "cm")) + scale_fill_manual(values = c('Alveolata undefined' = "#c83b4e",
                                                                                                                                          'Apicomplexa' = "#5b081a",
                                                                                                                                          'Ciliophora' = "#af1c3f",
                                                                                                                                          'Dinoflagellata undefined' = "#ff829c",
                                                                                                                                          'Dinophyceae' = "#fd3f65",
                                                                                                                                          'Ellobiophyceae' = "#9d002f",
                                                                                                                                          'Noctilucophyceae' = "#720000",
                                                                                                                                          'Syndiniales' = "#e1204d",
                                                                                                                                          'Cercozoa' = "#68992e",
                                                                                                                                          'Foraminifera' = "#adc57d",
                                                                                                                                          'Radiolaria' = "#2e6100",
                                                                                                                                          'Stramenopiles undefined' = "#008494",
                                                                                                                                          'Ochrophyta undefined' = "#006781",
                                                                                                                                          'Bacillariophyta' = "#56b7d0",
                                                                                                                                          'Bolidophyceae' = "#46819c",
                                                                                                                                          'Chrysophyceae' = "#9be1f5",
                                                                                                                                          'Dictyochophyceae' = "#006aa2",
                                                                                                                                          'MOCH' = "#005181",
                                                                                                                                          'Pelagophyceae' = "#65affe",
                                                                                                                                          'Phaeophyceae' = "#356fbe",
                                                                                                                                          'Raphidophyceae' = "#00326c",
                                                                                                                                          'Opalozoa' = "#003346",
                                                                                                                                          'Pseudofungi' = "#6a83c5",
                                                                                                                                          'Sagenista' = "#a4bfff",
                                                                                                                                          'Perkinsea' = "#646bb1",
                                                                                                                                          'Eustigmatophyceae' = "#0f316b",
                                                                                                                                          'Pinguiophyceae' = "#8374b9",
                                                                                                                                          'Synurophyceae' = "#9a99e0"
                                                                                                                                          )) +
  scale_y_continuous(expand = c(0,0)) + ylab("Relative abundance of ASV counts shared by all \nor unique to each region by taxa") +
  theme(axis.title.y = element_blank())
stacked_barplot_regions
#combine barplot and region plot
combo_plot <- ggarrange(ggarrange(venn_asvPlot, stacked_barplot_regions, align = "hv", nrow = 2, ncol = 1,
                                  heights = c(1,1), labels = c("a","b")), ggarrange(read_count_barplot, stacked_barplot, region_plot, 
                                            nrow = 3, ncol = 1, heights = c(3,3,2), align = "hv"),
                        align = "hv", nrow = 1, ncol = 2, labels = c("","c"))
combo_plot
ggsave("venn_figure.pdf", combo_plot, width = 12, height = 8)



