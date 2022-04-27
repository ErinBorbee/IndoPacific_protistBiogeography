#load packages
library(vegan)
library(EnvStats)
library(ggplot2)
library(reshape2)

#Set working directory
setwd("~/Desktop/Dissertation")

#Import data
data <- read.csv("~/Desktop/DataAnalysis/UVC/fish_abundance_by_site.csv")

#detect outliers and filter from each column
#carnivores
rosnerTest(data$carnivore, k = 5)
data$carnivore <- replace(data$carnivore, data$carnivore>37, NA)
#herbivores
rosnerTest(data$herbivore, k = 8)
data$herbivore <- replace(data$herbivore, data$herbivore>41, NA)
#benthic invertivores
rosnerTest(data$benthic.invertivore, k = 6)
data$benthic.invertivore <- replace(data$benthic.invertivore, data$benthic.invertivore>67, NA)
#omnivores
rosnerTest(data$omnivore, k = 2)
data$omnivore <- replace(data$omnivore, data$omnivore>43, NA)
#planktivores
rosnerTest(data$planktivore, k = 2)
data$planktivore <- replace(data$planktivore, data$planktivore>200, NA)
#corallivores
rosnerTest(data$corallivore, k = 6)
data$corallivore <- replace(data$corallivore, data$corallivore>600, NA)
#fish biomass
rosnerTest(data$fishBiomass, k = 5)
data$fishBiomass <- replace(data$fishBiomass, data$fishBiomass>100000, NA)
#target biomass
rosnerTest(data$target_biomass, k = 10)
data$target_biomass <- replace(target_biomass, data$target_biomass>52000, NA)
#nontarget biomass
rosnerTest(data$nontarget_biomass, k = 10)
data$nontarget_biomass <- replace(nontarget_biomass, data$nontarget_biomass>48000, NA)

#Test for normality in each column and transform to normal when needed
#carnivores
hist(data$carnivore, main = "carnivore", breaks = 10)
shapiro.test(data$carnivore)
#Normal without transformation

#herbivores
hist(data$herbivore, main = "herbivore", breaks = 10)
shapiro.test(data$herbivore)

data$herbivore_trans <- sqrt(data$herbivore)
hist(data$herbivore_trans, main = "herbivore transformed", breaks = 10)
shapiro.test(data$herbivore_trans)

#benthic invertivores
hist(data$benthic.invertivore, main = "benthic invertivore", breaks = 10)
shapiro.test(data$benthic.invertivore)

data$benthic.invertivore_trans <- 1/sqrt(data$benthic.invertivore)
hist(data$benthic.invertivore_trans, main = "benthic invertivore transformed", breaks = 10)
shapiro.test(data$benthic.invertivore_trans)

#omnivores
hist(data$omnivore, main = "omnivore", breaks = 10)
shapiro.test(data$omnivore)
#Normal without transformation

#planktivores
hist(data$planktivore, main = "planktivore", breaks = 10)
shapiro.test(data$planktivore)

data$planktivore_trans <- sqrt(data$planktivore)
hist(data$planktivore_trans, main = "planktivore transformed", breaks = 10)
shapiro.test(data$planktivore_trans)

#corallivore
hist(data$corallivore, main = "corallivore", breaks = 10)
shapiro.test(data$corallivore)

data$corallivore_trans <- sqrt(data$corallivore)
hist(data$corallivore_trans, main = "corallivore transformed", breaks = 10)
shapiro.test(data$corallivore_trans)

#fishBiomass
hist(data$fishBiomass, main = "fishBiomass", breaks = 10)
shapiro.test(data$fishBiomass)

data$fishBiomass_trans <- log10(data$fishBiomass)
hist(data$fishBiomass_trans, main = "fishBiomass transformed", breaks = 10)
shapiro.test(data$fishBiomass_trans)

#target biomass
hist(data$target_biomass, main = "target biomass", breaks = 10)
shapiro.test(data$target_biomass)

data$target_biomass_trans <- log10(data$target_biomass)
hist(data$target_biomass_trans, main = "target biomass transformed", breaks = 10)
shapiro.test(data$target_biomass_trans)

#nontarget biomass
hist(data$nontarget_biomass, main = "nontarget biomass", breaks = 10)
shapiro.test(data$nontarget_biomass)

data$nontarget_biomass_trans <- log10(data$nontarget_biomass)
hist(data$nontarget_biomass_trans, main = "nontarget biomass transformed", breaks = 10)
shapiro.test(data$nontarget_biomass_trans)

#Test for significant differences between individual columns and location
#Both ANOVA/Tukey and Kruskal-Wallis/Wilcoxon used for consistency throughout manuscript
#carnivore
carnivore_anova <- aov(carnivore ~ region, data = data)
summary(carnivore_anova)
TukeyHSD(carnivore_anova)

kruskal.test(carnivore ~ region, data = data)
pairwise.wilcox.test(data$carnivore, data$region, p.adjust.method = "BH")

#herbivore
herbivore_anova <- aov(herbivore_trans ~ region, data = data)
summary(herbivore_anova)
TukeyHSD(herbivore_anova)

kruskal.test(herbivore ~ region, data = data)
pairwise.wilcox.test(data$herbivore, data$region, p.adjust.method = "BH")

#benthic invertivore
benthic.invertivore_anova <- aov(benthic.invertivore_trans ~ region, data = data)
summary(benthic.invertivore_anova)
TukeyHSD(benthic.invertivore_anova)

kruskal.test(benthic.invertivore ~ region, data = data)
pairwise.wilcox.test(data$benthic.invertivore, data$region, p.adjust.method = "BH")

#omnivore
omnivore_anova <- aov(omnivore ~ region, data = data)
summary(omnivore_anova)
TukeyHSD(omnivore_anova)

kruskal.test(omnivore ~ region, data = data)
pairwise.wilcox.test(data$omnivore, data$region, p.adjust.method = "BH")

#planktivore
planktivore_anova <- aov(planktivore_trans ~ region, data = data)
summary(planktivore_anova)
TukeyHSD(planktivore_anova)

kruskal.test(planktivore ~ region, data = data)
pairwise.wilcox.test(data$planktivore, data$region, p.adjust.method = "BH")

#corallivore
corallivore_anova <- aov(corallivore_trans ~ region, data = data)
summary(corallivore_anova)
TukeyHSD(corallivore_anova)

kruskal.test(corallivore ~ region, data = data)
pairwise.wilcox.test(data$corallivore, data$region, p.adjust.method = "BH")

#fishBiomass
fishBiomass_anova <- aov(fishBiomass_trans ~ region, data = data)
summary(fishBiomass_anova)
TukeyHSD(fishBiomass_anova)

kruskal.test(fishBiomass ~ region, data = data)
pairwise.wilcox.test(data$fishBiomass, data$region, p.adjust.method = "BH")

#target_biomass
target_biomass_anova <- aov(target_biomass_trans ~ region, data = data)
summary(target_biomass_anova)
TukeyHSD(target_biomass_anova)

kruskal.test(target_biomass ~ region, data = data)
pairwise.wilcox.test(data$target_biomass, data$region, p.adjust.method = "BH")

#nontarget_biomass
nontarget_biomass_anova <- aov(nontarget_biomass_trans ~ region, data = data)
summary(nontarget_biomass_anova)
TukeyHSD(nontarget_biomass_anova)

kruskal.test(nontarget_biomass ~ region, data = data)
pairwise.wilcox.test(data$nontarget_biomass, data$region, p.adjust.method = "BH")

#plot boxplots
names(data)
data <- data[-c(4,7,11,13)]
data_long <- melt(data)
data_long$region <- factor(data_long$region, levels = c("Lombok",
                                                        "Wakatobi",
                                                        "Misool",
                                                        "Waigeo"))
#Export csv, manually enter color column for significance and import
write.csv(data_long, "data_long.csv")
data_long <- read.csv("data_long.csv")
data_long$region <- factor(data_long$region, levels = c("Lombok",
                                                        "Wakatobi",
                                                        "Misool",
                                                        "Waigeo"))

#define colors
colors <- c('high' = "#f73e42",
            'low' = "#007cb5",
            'none' = "white")

plot <- ggplot(data_long, aes(region, value, fill = color)) + geom_boxplot(alpha = 0.75) + 
  facet_wrap(~variable, scales = "free") + scale_fill_manual(values = colors) + theme_bw() +
  theme(panel.grid = element_blank()) + theme(strip.background = element_rect(fill = "black")) +
  theme(strip.text = element_text(face = "bold", color = "white")) + theme(legend.position = "none") +
  ylab("abundance") + theme(axis.title.x = element_blank())
plot
ggsave("Figure7.pdf", plot, height = 6, width = 7)


