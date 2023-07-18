#load packages
library(ggplot2)
library(scatterpie)

#import data
data <- read.csv("SAR_piechart.csv")

data$taxa <- factor(data$taxa, levels = c("Other", 
                                          "Foraminifera",
                                          "Cercozoa",
                                          "Radiolaria",
                                          "Apicomplexa",
                                          "Ciliophora",
                                          "Dinoflagellata",
                                          "MAST",
                                          "Pelagophyceae",
                                          "Chrysophyceae",
                                          "Bacillariophyta"))

data$region <- factor(data$region, levels = c("Lombok", "Wakatobi", "Misool","Waigeo"))

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

#plot piechart
ggplot(data, aes(x = "", y = percent, fill = taxa)) + 
  geom_bar(stat = "identity", width = 1, color = "white", size = 0.25) + 
  geom_text(aes(y = ypos, label = label), color = "white", size = 5) + coord_polar("y", start = 0) +
  theme_bw() + theme(panel.grid = element_blank()) + 
  theme(panel.border = element_rect(color = "grey20", size = 1)) +
  scale_fill_manual(values = colors) + 
  theme(axis.ticks.x = element_blank(), axis.title.x = element_blank(),
        axis.text.x = element_blank()) +
  theme(axis.ticks.length.y = unit(0.1, "in"),
        axis.text.y = element_blank(),
        axis.title.y = element_blank()) + 
  theme(legend.position = "none") + 
  facet_grid(~region) + 
  theme(strip.background = element_rect(fill = "grey20", colour = "grey20", size = 1),
        strip.text = element_text(face = "bold", color = "white", size = 12))
ggsave("SAR_piecharts.png", width = 10, height = 3)
