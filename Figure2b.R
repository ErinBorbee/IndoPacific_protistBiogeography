#load packages
library(reshape2)
library(ggplot2)

#Import data
data <- read.csv("SAR_OTU_percent_L2.csv")
data_long <- melt(data)
write.csv(data_long,"SAR_OTU_L2_long.csv")

#Import data after adding and collapsing by region
data <- read.csv("SAR_OTU_L2_long_by_subset.csv")

data$region <- factor(data$region, levels = c("Lombok", "Wakatobi", "Misool", "Waigeo"))

data$taxa <- factor(data$taxa, levels = c("Dinoflagellata",
                                      "Ciliophora",
                                      "Apicomplexa",
                                      "Bacillariophyta",
                                      "Chrysophyceae",
                                      "Pelagophyceae",
                                      "MAST",
                                      "Radiolaria",
                                      "Cercozoa",
                                      "Foraminifera",
                                      "Other"))

#Define colors
colors <- c('Dinoflagellata' = "#fd3f65",
            'Ciliophora' = "#af1c3f",
            'Apicomplexa' = "#5b081a",
            'Bacillariophyta' = "#56b7d0",
            'Chrysophyceae' = "#2f8db1",
            'Pelagophyceae' = "#215a8c",
            'MAST' = "#00325c",
            'Radiolaria' = "#adc57d",
            'Cercozoa' = "#70983a",
            'Foraminifera' = "#364619",
            'Other' = "grey50")

#plot stacked barplot
plot <- ggplot(data, aes(region, percent, fill = taxa)) + 
        geom_bar(stat = "identity", position = "stack") + 
        theme_bw() + theme(panel.grid = element_blank()) + 
        theme(panel.border = element_rect(color = "grey20", size = 1)) +
        scale_fill_manual(values = colors) + 
        theme(axis.ticks.x = element_blank(), axis.title.x = element_blank(),
              axis.text.x = element_text(face = "bold", size = 12)) +
        theme(axis.ticks.length.y = unit(0.1, "in"),
              axis.text.y = element_text(size = 12),
              axis.title.y = element_text(size = 12)) + 
        ylab("Relative Abundance \n") +
        theme(legend.text = element_text(size = 12),
              legend.title = element_text(size = 12), 
              legend.position = "bottom") + 
        facet_grid(~subset) + 
        theme(strip.background = element_rect(fill = "grey20", colour = "grey20", size = 1),
              strip.text = element_text(face = "bold", color = "white", size = 12))
plot
ggsave("SAR_barplot_by_subset.png", width = 15, height = 5)
