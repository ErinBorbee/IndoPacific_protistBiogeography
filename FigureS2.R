library(ggplot2)
library(ggpubr)

M04_data <- read.csv("M04-observed_otus.csv")
M12_data <- read.csv("M12-observed_otus.csv")
S04_data <- read.csv("S04-observed_otus.csv")
S12_data <- read.csv("S12-observed_otus.csv")

M04_data <- melt(M04_data)
M12_data <- melt(M12_data)
S04_data <- melt(S04_data)
S12_data <- melt(S12_data)

write.csv(M04_data, "M04_observed_otus_long.csv")
write.csv(M12_data, "M12_observed_otus_long.csv")
write.csv(S04_data, "S04_observed_otus_long.csv")
write.csv(S12_data, "S12_observed_otus_long.csv")

M04_data_long <- read.csv("M04_observed_otus_long.csv")
M12_data_long <- read.csv("M12_observed_otus_long.csv")
S04_data_long <- read.csv("S04_observed_otus_long.csv")
S12_data_long <- read.csv("S12_observed_otus_long.csv")

colors <- c('Lombok' = "#8e86c5",
            'Wakatobi' = "#007cb5",
            'Misool' = "#a6c54b",
            'Waigeo' = "#f73e42")

M04_plot <- ggplot(M04_data_long, aes(depth, value, color = region)) + geom_smooth(se = FALSE) +
  theme_bw() + theme(panel.grid = element_blank()) + scale_x_continuous(expand = c(0,0), limits = c(0,2500)) + scale_y_continuous(expand = c(0,0), limits = c(0,250)) +
  scale_color_manual(values = colors) + ylab("Observed OTUs\n") + theme(axis.title = element_blank())
M04_plot

M12_plot <- ggplot(M12_data_long, aes(depth, value, color = region)) + geom_smooth(se = FALSE) +
  theme_bw() + theme(panel.grid = element_blank()) + scale_x_continuous(expand = c(0,0), limits = c(0,2500)) + scale_y_continuous(expand = c(0,0), limits = c(0,250)) +
  scale_color_manual(values = colors) + ylab("Observed OTUs\n") + theme(axis.title = element_blank())
M12_plot

S04_plot <- ggplot(S04_data_long, aes(depth, value, color = region)) + geom_smooth(se = FALSE) +
  theme_bw() + theme(panel.grid = element_blank()) + scale_x_continuous(expand = c(0,0), limits = c(0,2500)) + scale_y_continuous(expand = c(0,0), limits = c(0,250)) +
  scale_color_manual(values = colors) + ylab("Observed OTUs\n") + theme(axis.title = element_blank())
S04_plot

S12_plot <- ggplot(S12_data_long, aes(depth, value, color = region)) + geom_smooth(se = FALSE) +
  theme_bw() + theme(panel.grid = element_blank()) + scale_x_continuous(expand = c(0,0), limits = c(0,2500)) + scale_y_continuous(expand = c(0,0), limits = c(0,250)) +
  scale_color_manual(values = colors) + ylab("Observed OTUs\n") + theme(axis.title = element_blank())
S12_plot

plot <- ggarrange(M04_plot,M12_plot,S04_plot,S12_plot, ncol = 1, nrow = 4, common.legend = TRUE,
                  align = "hv", legend = "bottom")
plot

plot <- annotate_figure(plot, left = "Observed OTUs",)
plot

ggsave("rarefaction_byregion.png", plot, width = 8, height = 10) 
