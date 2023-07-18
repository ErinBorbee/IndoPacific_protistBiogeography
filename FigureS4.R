#load packages
library(vegan)
library(EnvStats)
library(reshape2)
library(ggplot2)

setwd("~/Desktop/Dissertation/Chapter1/")

#import data
data <- read.csv("~/Desktop/DataAnalysis/UVC/fish_abundance_by_site.csv", row.names = 1)

#Detect and filter outliers
rosnerTest(data$carnivore, k = 2)
data$carnivore <- replace(data$carnivore, data$carnivore>37, NA)

rosnerTest(data$herbivore, k = 8)
data$herbivore <- replace(data$herbivore, data$herbivore>41, NA)

rosnerTest(data$benthic.invertivore, k = 6)
data$benthic.invertivore <- replace(data$benthic.invertivore, data$benthic.invertivore>67, NA)

rosnerTest(data$omnivore, k = 2)
data$omnivore <- replace(data$omnivore, data$omnivore>43, NA)

rosnerTest(data$planktivore, k = 2)
data$planktivore <- replace(data$planktivore, data$planktivore>200, NA)

rosnerTest(data$corallivore, k = 6)
data$corallivore <- replace(data$corallivore, data$corallivore>600, NA)

rosnerTest(data$protist, k = 5)

#test each column of dataset for normality (p > 0.05 = normal)
#carnivores
hist(data$carnivore, main = "carnivore", breaks = 10)
shapiro.test(data$carnivore)
#Normal without transformation

#Repeat above two steps for each column
#herbivores
hist(data$herbivore, main = "herbivore", breaks = 10)
shapiro.test(data$herbivore)

#transform non-normal data and re-test for normality
data$herbivore_trans <- sqrt(data$herbivore)
hist(data$herbivore_trans, main = "herbivore transformed", breaks = 10)
shapiro.test(data$herbivore_trans)

#benthic invertivores
hist(data$benthic.invertivore, main = "benthic invertivore", breaks = 10)
shapiro.test(data$benthic.invertivore)

data$benthic.invertivore_trans <- log10(data$benthic.invertivore)
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

#corallivores
hist(data$corallivore, main = "corallivore", breaks = 10)
shapiro.test(data$corallivore)

data$corallivore_trans <- sqrt(data$corallivore)
hist(data$corallivore_trans, main = "corallivore transformed", breaks = 10)
shapiro.test(data$corallivore_trans)

#protists
hist(data$protist_water_97, main = "protist OTU richness", breaks = 10)
shapiro.test(data$protist_97)

data$protist_trans <- sqrt(data$protist)
hist(data$protist_trans, main = "protist OTU richness transformed", breaks = 10)
shapiro.test(data$protist_trans)

#Test for significance across fishing restrictions using ANOVA (parametric with transformed data)
#Kruskal-Wallis is also used as non-parametric on non-transformed data for consistency with
#other statistical comparisons in the manuscript
#If ANOVA is significant use Tukey to run pairwise comparisons (when using Kruskal-Wallis,
#use Wilcoxon for pairwise comparisons)
#carnivores
carnivore_anova <- aov(carnivore ~ restriction, data = data)
summary(carnivore_anova)
TukeyHSD(carnivore_anova)

kruskal.test(carnivore ~ restriction, data = data)
pairwise.wilcox.test(data$carnivore, data$restriction, p.adjust.method = "BH")

#herbivores
herbivore_anova <- aov(herbivore_trans ~ restriction, data = data)
summary(herbivore_anova)
TukeyHSD(herbivore_anova)

kruskal.test(herbivore ~ restriction, data = data)
pairwise.wilcox.test(data$herbivore, data$restriction, p.adjust.method = "BH")

#benthic invertivores
benthic.invertivore_anova <- aov(benthic.invertivore_trans ~ restriction, data = data)
summary(benthic.invertivore_anova)
TukeyHSD(benthic.invertivore_anova)

kruskal.test(benthic.invertivore ~ restriction, data = data)
pairwise.wilcox.test(data$benthic.invertivore, data$restriction, p.adjust.method = "BH")

#omnivores
omnivore_anova <- aov(omnivore ~ restriction, data = data)
summary(omnivore_anova)
TukeyHSD(omnivore_anova)

kruskal.test(omnivore ~ restriction, data = data)
pairwise.wilcox.test(data$omnivore, data$restriction, p.adjust.method = "BH")

#planktivores
planktivore_anova <- aov(planktivore_trans ~ restriction, data = data)
summary(planktivore_anova)
TukeyHSD(planktivore_anova)

kruskal.test(planktivore ~ restriction, data = data)
pairwise.wilcox.test(data$planktivore, data$restriction, p.adjust.method = "BH")

#corallivore
corallivore_anova <- aov(corallivore_trans ~ restriction, data = data)
summary(corallivore_anova)
TukeyHSD(corallivore_anova)

kruskal.test(corallivore ~ restriction, data = data)
pairwise.wilcox.test(data$corallivore, data$restriction, p.adjust.method = "BH")

#protists
protist_anova <- aov(protist ~ restriction, data = data)
summary(protist_anova)
TukeyHSD(protist_anova)

kruskal.test(protist ~ restriction, data = data)
pairwise.wilcox.test(data$protist, data$restriction, p.adjust.method = "BH")

#subset data to individual regions to identify patterns within sampling region
LBK_data <- subset(data, region == "Lombok")
WK_data <- subset(data, region == "Wakatobi")
MIS_data <- subset(data, region == "Misool")
WGO_data <- subset(data, region == "Waigeo")

#Rerun stats from above on individual regions
#test each column of LBK_data for normality (p > 0.05 = normal)
#carnivores
hist(LBK_data$carnivore, main = "carnivore", breaks = 10)
shapiro.test(LBK_data$carnivore)

#transform non-normal LBK_data and re-test for normality
LBK_data$carnivore_trans <- (LBK_data$carnivore)^1.1
hist(LBK_data$carnivore_trans, main = "carnivore transformed", breaks = 10)
shapiro.test(LBK_data$carnivore_trans)

#Repeat above two steps for each column
#herbivores
hist(LBK_data$herbivore, main = "herbivore", breaks = 10)
shapiro.test(LBK_data$herbivore)

LBK_data$herbivore_trans <- sqrt(LBK_data$herbivore)
hist(LBK_data$herbivore_trans, main = "herbivore transformed", breaks = 10)
shapiro.test(LBK_data$herbivore_trans)

#benthic invertivores
hist(LBK_data$benthic.invertivore, main = "benthic invertivore", breaks = 10)
shapiro.test(LBK_data$benthic.invertivore)

LBK_data$benthic.invertivore_trans <- sqrt(LBK_data$benthic.invertivore)
hist(LBK_data$benthic.invertivore_trans, main = "benthic invertivore transformed", breaks = 10)
shapiro.test(LBK_data$benthic.invertivore_trans)

#omnivores
hist(LBK_data$omnivore, main = "omnivore", breaks = 10)
shapiro.test(LBK_data$omnivore)
#Normal without transformation

#planktivores
hist(LBK_data$planktivore, main = "planktivore", breaks = 10)
shapiro.test(LBK_data$planktivore)
#Normal without transformation

#corallivores
hist(LBK_data$coralivore, main = "corallivore", breaks = 10)
shapiro.test(LBK_data$coralivore)

LBK_data$corallivore_trans <- log10(LBK_data$coralivore)
hist(LBK_data$corallivore_trans, main = "corallivore transformed", breaks = 10)
shapiro.test(LBK_data$corallivore_trans)

#protists
hist(LBK_data$protist, main = "protist OTU richness", breaks = 10)
shapiro.test(LBK_data$protist)
#Normal without transformation

#Test for significance across fishing restrictions using ANOVA (parametric with transformed LBK_data)
#Kruskal-Wallis is also used as non-parametric on non-transformed LBK_data for consistency with
#other statistical comparisons in the manuscript
#If ANOVA is significant use Tukey to run pairwise comparisons (when using Kruskal-Wallis,
#use Wilcoxon for pairwise comparisons)
#carnivores
carnivore_anova <- aov(carnivore_trans ~ restriction, data = LBK_data)
summary(carnivore_anova)
TukeyHSD(carnivore_anova)

kruskal.test(carnivore ~ restriction, data = LBK_data)
pairwise.wilcox.test(LBK_data$carnivore, LBK_data$restriction, p.adjust.method = "BH")

#herbivores
herbivore_anova <- aov(herbivore_trans ~ restriction, data = LBK_data)
summary(herbivore_anova)
TukeyHSD(herbivore_anova)

kruskal.test(herbivore ~ restriction, data = LBK_data)
pairwise.wilcox.test(LBK_data$herbivore, LBK_data$restriction, p.adjust.method = "BH")

#benthic invertivores
benthic.invertivore_anova <- aov(benthic.invertivore_trans ~ restriction, data = LBK_data)
summary(benthic.invertivore_anova)
TukeyHSD(benthic.invertivore_anova)

kruskal.test(benthic.invertivore ~ restriction, data = LBK_data)
pairwise.wilcox.test(LBK_data$benthic.invertivore, LBK_data$restriction, p.adjust.method = "BH")

#omnivores
omnivore_anova <- aov(omnivore_trans ~ restriction, LBK_data = LBK_data)
summary(omnivore_anova)
TukeyHSD(omnivore_anova)

kruskal.test(omnivore ~ restriction, data = LBK_data)
pairwise.wilcox.test(LBK_data$omnivore, LBK_data$restriction, p.adjust.method = "BH")

#planktivores
planktivore_anova <- aov(planktivore ~ restriction, data = LBK_data)
summary(planktivore_anova)
TukeyHSD(planktivore_anova)

kruskal.test(planktivore ~ restriction, data = LBK_data)
pairwise.wilcox.test(LBK_data$planktivore, LBK_data$restriction, p.adjust.method = "BH")

#corallivore
corallivore_anova <- aov(corallivore_trans ~ restriction, data = LBK_data)
summary(corallivore_anova)
TukeyHSD(corallivore_anova)

kruskal.test(coralivore ~ restriction, data = LBK_data)
pairwise.wilcox.test(LBK_data$coralivore, LBK_data$restriction, p.adjust.method = "BH")

#protists
protist_anova <- aov(protist ~ restriction, data = LBK_data)
summary(protist_anova)
TukeyHSD(protist_anova)

kruskal.test(protist ~ restriction, data = LBK_data)
pairwise.wilcox.test(LBK_data$protist, LBK_data$restriction, p.adjust.method = "BH")

#Rerun stats from above on individual regions
#test each column of WK_data for normality (p > 0.05 = normal)
#carnivores
hist(WK_data$carnivore, main = "carnivore", breaks = 10)
shapiro.test(WK_data$carnivore)
#Normal without transformation

#herbivores
hist(WK_data$herbivore, main = "herbivore", breaks = 10)
shapiro.test(WK_data$herbivore)

#transform non-normal WK_data and re-test for normality
WK_data$herbivore_trans <- sqrt(WK_data$herbivore)
hist(WK_data$herbivore_trans, main = "herbivore transformed", breaks = 10)
shapiro.test(WK_data$herbivore_trans)

#benthic invertivores
hist(WK_data$benthic.invertivore, main = "benthic invertivore", breaks = 10)
shapiro.test(WK_data$benthic.invertivore)

WK_data$benthic.invertivore_trans <- log10(WK_data$benthic.invertivore)
hist(WK_data$benthic.invertivore_trans, main = "benthic invertivore transformed", breaks = 10)
shapiro.test(WK_data$benthic.invertivore_trans)

#omnivores
hist(WK_data$omnivore, main = "omnivore", breaks = 10)
shapiro.test(WK_data$omnivore)

WK_data$omnivore_trans <- log10(WK_data$omnivore)
hist(WK_data$omnivore_trans, main = "omnivore transformed", breaks = 10)
shapiro.test(WK_data$omnivore_trans)

#planktivores
hist(WK_data$planktivore, main = "planktivore", breaks = 10)
shapiro.test(WK_data$planktivore)
#Normal without transformation

#corallivores
hist(WK_data$coralivore, main = "corallivore", breaks = 10)
shapiro.test(WK_data$coralivore)
#Normal without transformation

#protists
hist(WK_data$protist, main = "protist OTU richness", breaks = 10)
shapiro.test(WK_data$protist)
#Normal without transformation

#Test for significance across fishing restrictions using ANOVA (parametric with transformed WK_data)
#Kruskal-Wallis is also used as non-parametric on non-transformed WK_data for consistency with
#other statistical comparisons in the manuscript
#If ANOVA is significant use Tukey to run pairwise comparisons (when using Kruskal-Wallis,
#use Wilcoxon for pairwise comparisons)
#carnivores
carnivore_anova <- aov(carnivore_trans ~ restriction, data = WK_data)
summary(carnivore_anova)
TukeyHSD(carnivore_anova)

kruskal.test(carnivore ~ restriction, data = WK_data)
pairwise.wilcox.test(WK_data$carnivore, WK_data$restriction, p.adjust.method = "BH")

#herbivores
herbivore_anova <- aov(herbivore_trans ~ restriction, data = WK_data)
summary(herbivore_anova)
TukeyHSD(herbivore_anova)

kruskal.test(herbivore ~ restriction, data = WK_data)
pairwise.wilcox.test(WK_data$herbivore, WK_data$restriction, p.adjust.method = "BH")

#benthic invertivores
benthic.invertivore_anova <- aov(benthic.invertivore_trans ~ restriction, data = WK_data)
summary(benthic.invertivore_anova)
TukeyHSD(benthic.invertivore_anova)

kruskal.test(benthic.invertivore ~ restriction, data = WK_data)
pairwise.wilcox.test(WK_data$benthic.invertivore, WK_data$restriction, p.adjust.method = "BH")

#omnivores
omnivore_anova <- aov(omnivore_trans ~ restriction, WK_data = WK_data)
summary(omnivore_anova)
TukeyHSD(omnivore_anova)

kruskal.test(omnivore ~ restriction, data = WK_data)
pairwise.wilcox.test(WK_data$omnivore, WK_data$restriction, p.adjust.method = "BH")

#planktivores
planktivore_anova <- aov(planktivore ~ restriction, data = WK_data)
summary(planktivore_anova)
TukeyHSD(planktivore_anova)

kruskal.test(planktivore ~ restriction, data = WK_data)
pairwise.wilcox.test(WK_data$planktivore, WK_data$restriction, p.adjust.method = "BH")

#corallivore
corallivore_anova <- aov(corallivore_trans ~ restriction, data = WK_data)
summary(corallivore_anova)
TukeyHSD(corallivore_anova)

kruskal.test(coralivore ~ restriction, data = WK_data)
pairwise.wilcox.test(WK_data$coralivore, WK_data$restriction, p.adjust.method = "BH")

#protists
protist_anova <- aov(protist ~ restriction, data = WK_data)
summary(protist_anova)
TukeyHSD(protist_anova)

kruskal.test(protist ~ restriction, data = WK_data)
pairwise.wilcox.test(WK_data$protist, WK_data$restriction, p.adjust.method = "BH")

#Rerun stats from above on individual regions
#test each column of MIS_data for normality (p > 0.05 = normal)
#carnivores
hist(MIS_data$carnivore, main = "carnivore", breaks = 10)
shapiro.test(MIS_data$carnivore)
#Normal without transformation

#herbivores
hist(MIS_data$herbivore, main = "herbivore", breaks = 10)
shapiro.test(MIS_data$herbivore)

#transform non-normal MIS_data and re-test for normality
MIS_data$herbivore_trans <- 1/sqrt(MIS_data$herbivore)
hist(MIS_data$herbivore_trans, main = "herbivore transformed", breaks = 10)
shapiro.test(MIS_data$herbivore_trans)

#benthic invertivores
hist(MIS_data$benthic.invertivore, main = "benthic invertivore", breaks = 10)
shapiro.test(MIS_data$benthic.invertivore)

MIS_data$benthic.invertivore_trans <- log10(MIS_data$benthic.invertivore)
hist(MIS_data$benthic.invertivore_trans, main = "benthic invertivore transformed", breaks = 10)
shapiro.test(MIS_data$benthic.invertivore_trans)

#omnivores
hist(MIS_data$omnivore, main = "omnivore", breaks = 10)
shapiro.test(MIS_data$omnivore)
#Normal without transformation

#planktivores
hist(MIS_data$planktivore, main = "planktivore", breaks = 10)
shapiro.test(MIS_data$planktivore)
#Normal without transformation

#corallivores
hist(MIS_data$coralivore, main = "corallivore", breaks = 10)
shapiro.test(MIS_data$coralivore)

MIS_data$corallivore_trans <- sqrt(MIS_data$coralivore)
hist(MIS_data$corallivore_trans, main = "corallivore transformed", breaks = 10)
shapiro.test(MIS_data$corallivore_trans)

#protists
hist(MIS_data$protist, main = "protist OTU richness", breaks = 10)
shapiro.test(MIS_data$protist)
#Normal without transformation

#Test for significance across fishing restrictions using ANOVA (parametric with transformed MIS_data)
#Kruskal-Wallis is also used as non-parametric on non-transformed MIS_data for consistency with
#other statistical comparisons in the manuscript
#If ANOVA is significant use Tukey to run pairwise comparisons (when using Kruskal-Wallis,
#use Wilcoxon for pairwise comparisons)
#carnivores
carnivore_anova <- aov(carnivore_trans ~ restriction, data = MIS_data)
summary(carnivore_anova)
TukeyHSD(carnivore_anova)

kruskal.test(carnivore ~ restriction, data = MIS_data)
pairwise.wilcox.test(MIS_data$carnivore, MIS_data$restriction, p.adjust.method = "BH")

#herbivores
herbivore_anova <- aov(herbivore_trans ~ restriction, data = MIS_data)
summary(herbivore_anova)
TukeyHSD(herbivore_anova)

kruskal.test(herbivore ~ restriction, data = MIS_data)
pairwise.wilcox.test(MIS_data$herbivore, MIS_data$restriction, p.adjust.method = "BH")

#benthic invertivores
benthic.invertivore_anova <- aov(benthic.invertivore_trans ~ restriction, data = MIS_data)
summary(benthic.invertivore_anova)
TukeyHSD(benthic.invertivore_anova)

kruskal.test(benthic.invertivore ~ restriction, data = MIS_data)
pairwise.wilcox.test(MIS_data$benthic.invertivore, MIS_data$restriction, p.adjust.method = "BH")

#omnivores
omnivore_anova <- aov(omnivore_trans ~ restriction, MIS_data = MIS_data)
summary(omnivore_anova)
TukeyHSD(omnivore_anova)

kruskal.test(omnivore ~ restriction, data = MIS_data)
pairwise.wilcox.test(MIS_data$omnivore, MIS_data$restriction, p.adjust.method = "BH")

#planktivores
planktivore_anova <- aov(planktivore ~ restriction, data = MIS_data)
summary(planktivore_anova)
TukeyHSD(planktivore_anova)

kruskal.test(planktivore ~ restriction, data = MIS_data)
pairwise.wilcox.test(MIS_data$planktivore, MIS_data$restriction, p.adjust.method = "BH")

#corallivore
corallivore_anova <- aov(corallivore_trans ~ restriction, data = MIS_data)
summary(corallivore_anova)
TukeyHSD(corallivore_anova)

kruskal.test(coralivore ~ restriction, data = MIS_data)
pairwise.wilcox.test(MIS_data$coralivore, MIS_data$restriction, p.adjust.method = "BH")

#protists
protist_anova <- aov(protist ~ restriction, data = MIS_data)
summary(protist_anova)
TukeyHSD(protist_anova)

kruskal.test(protist ~ restriction, data = MIS_data)
pairwise.wilcox.test(MIS_data$protist, MIS_data$restriction, p.adjust.method = "BH")

#Rerun stats from above on individual regions
#test each column of MIS_data for normality (p > 0.05 = normal)
#carnivores
hist(MIS_data$carnivore, main = "carnivore", breaks = 10)
shapiro.test(MIS_data$carnivore)
#Normal without transformation

#herbivores
hist(MIS_data$herbivore, main = "herbivore", breaks = 10)
shapiro.test(MIS_data$herbivore)

#transform non-normal MIS_data and re-test for normality
MIS_data$herbivore_trans <- sqrt(MIS_data$herbivore)
hist(MIS_data$herbivore_trans, main = "herbivore transformed", breaks = 10)
shapiro.test(MIS_data$herbivore_trans)

#benthic invertivores
hist(MIS_data$benthic.invertivore, main = "benthic invertivore", breaks = 10)
shapiro.test(MIS_data$benthic.invertivore)

MIS_data$benthic.invertivore_trans <- log10(MIS_data$benthic.invertivore)
hist(MIS_data$benthic.invertivore_trans, main = "benthic invertivore transformed", breaks = 10)
shapiro.test(MIS_data$benthic.invertivore_trans)

#omnivores
hist(MIS_data$omnivore, main = "omnivore", breaks = 10)
shapiro.test(MIS_data$omnivore)

MIS_data$omnivore_trans <- log10(MIS_data$omnivore)
hist(MIS_data$omnivore_trans, main = "omnivore transformed", breaks = 10)
shapiro.test(MIS_data$omnivore_trans)

#planktivores
hist(MIS_data$planktivore, main = "planktivore", breaks = 10)
shapiro.test(MIS_data$planktivore)
#Normal without transformation

#corallivores
hist(MIS_data$coralivore, main = "corallivore", breaks = 10)
shapiro.test(MIS_data$coralivore)
#Normal without transformation

#protists
hist(MIS_data$protist, main = "protist OTU richness", breaks = 10)
shapiro.test(MIS_data$protist)
#Normal without transformation

#Test for significance across fishing restrictions using ANOVA (parametric with transformed MIS_data)
#Kruskal-Wallis is also used as non-parametric on non-transformed MIS_data for consistency with
#other statistical comparisons in the manuscript
#If ANOVA is significant use Tukey to run pairwise comparisons (when using Kruskal-Wallis,
#use Wilcoxon for pairwise comparisons)
#carnivores
carnivore_anova <- aov(carnivore ~ restriction, data = MIS_data)
summary(carnivore_anova)
TukeyHSD(carnivore_anova)

kruskal.test(carnivore ~ restriction, data = MIS_data)
pairwise.wilcox.test(MIS_data$carnivore, MIS_data$restriction, p.adjust.method = "BH")

#herbivores
herbivore_anova <- aov(herbivore_trans ~ restriction, data = MIS_data)
summary(herbivore_anova)
TukeyHSD(herbivore_anova)

kruskal.test(herbivore ~ restriction, data = MIS_data)
pairwise.wilcox.test(MIS_data$herbivore, MIS_data$restriction, p.adjust.method = "BH")

#benthic invertivores
benthic.invertivore_anova <- aov(benthic.invertivore_trans ~ restriction, data = MIS_data)
summary(benthic.invertivore_anova)
TukeyHSD(benthic.invertivore_anova)

kruskal.test(benthic.invertivore ~ restriction, data = MIS_data)
pairwise.wilcox.test(MIS_data$benthic.invertivore, MIS_data$restriction, p.adjust.method = "BH")

#omnivores
omnivore_anova <- aov(omnivore_trans ~ restriction, MIS_data = MIS_data)
summary(omnivore_anova)
TukeyHSD(omnivore_anova)

kruskal.test(omnivore ~ restriction, data = MIS_data)
pairwise.wilcox.test(MIS_data$omnivore, MIS_data$restriction, p.adjust.method = "BH")

#planktivores
planktivore_anova <- aov(planktivore ~ restriction, data = MIS_data)
summary(planktivore_anova)
TukeyHSD(planktivore_anova)

kruskal.test(planktivore ~ restriction, data = MIS_data)
pairwise.wilcox.test(MIS_data$planktivore, MIS_data$restriction, p.adjust.method = "BH")

#corallivore
corallivore_anova <- aov(coralivore ~ restriction, data = MIS_data)
summary(corallivore_anova)
TukeyHSD(corallivore_anova)

kruskal.test(coralivore ~ restriction, data = MIS_data)
pairwise.wilcox.test(MIS_data$coralivore, MIS_data$restriction, p.adjust.method = "BH")

#protists
protist_anova <- aov(protist ~ restriction, data = MIS_data)
summary(protist_anova)
TukeyHSD(protist_anova)

kruskal.test(protist ~ restriction, data = MIS_data)
pairwise.wilcox.test(MIS_data$protist, MIS_data$restriction, p.adjust.method = "BH")

#Rerun stats from above on individual regions
#test each column of MIS_data for normality (p > 0.05 = normal)
#carnivores
hist(MIS_data$carnivore, main = "carnivore", breaks = 10)
shapiro.test(MIS_data$carnivore)
#Normal without transformation

#herbivores
hist(MIS_data$herbivore, main = "herbivore", breaks = 10)
shapiro.test(MIS_data$herbivore)

#transform non-normal MIS_data and re-test for normality
MIS_data$herbivore_trans <- 1/sqrt(MIS_data$herbivore)
hist(MIS_data$herbivore_trans, main = "herbivore transformed", breaks = 10)
shapiro.test(MIS_data$herbivore_trans)

#benthic invertivores
hist(MIS_data$benthic.invertivore, main = "benthic invertivore", breaks = 10)
shapiro.test(MIS_data$benthic.invertivore)

MIS_data$benthic.invertivore_trans <- log10(MIS_data$benthic.invertivore)
hist(MIS_data$benthic.invertivore_trans, main = "benthic invertivore transformed", breaks = 10)
shapiro.test(MIS_data$benthic.invertivore_trans)

#omnivores
hist(MIS_data$omnivore, main = "omnivore", breaks = 10)
shapiro.test(MIS_data$omnivore)
#Normal without transformation

#planktivores
hist(MIS_data$planktivore, main = "planktivore", breaks = 10)
shapiro.test(MIS_data$planktivore)
#Normal without transformation

#corallivores
hist(MIS_data$coralivore, main = "corallivore", breaks = 10)
shapiro.test(MIS_data$coralivore)

MIS_data$corallivore_trans <- sqrt(MIS_data$coralivore)
hist(MIS_data$corallivore_trans, main = "corallivore transformed", breaks = 10)
shapiro.test(MIS_data$corallivore_trans)

#protists
hist(MIS_data$protist, main = "protist OTU richness", breaks = 10)
shapiro.test(MIS_data$protist)
#Normal without transformation

#Test for significance across fishing restrictions using ANOVA (parametric with transformed MIS_data)
#Kruskal-Wallis is also used as non-parametric on non-transformed MIS_data for consistency with
#other statistical comparisons in the manuscript
#If ANOVA is significant use Tukey to run pairwise comparisons (when using Kruskal-Wallis,
#use Wilcoxon for pairwise comparisons)
#carnivores
carnivore_anova <- aov(carnivore_trans ~ restriction, data = MIS_data)
summary(carnivore_anova)
TukeyHSD(carnivore_anova)

kruskal.test(carnivore ~ restriction, data = MIS_data)
pairwise.wilcox.test(MIS_data$carnivore, MIS_data$restriction, p.adjust.method = "BH")

#herbivores
herbivore_anova <- aov(herbivore_trans ~ restriction, data = MIS_data)
summary(herbivore_anova)
TukeyHSD(herbivore_anova)

kruskal.test(herbivore ~ restriction, data = MIS_data)
pairwise.wilcox.test(MIS_data$herbivore, MIS_data$restriction, p.adjust.method = "BH")

#benthic invertivores
benthic.invertivore_anova <- aov(benthic.invertivore_trans ~ restriction, data = MIS_data)
summary(benthic.invertivore_anova)
TukeyHSD(benthic.invertivore_anova)

kruskal.test(benthic.invertivore ~ restriction, data = MIS_data)
pairwise.wilcox.test(MIS_data$benthic.invertivore, MIS_data$restriction, p.adjust.method = "BH")

#omnivores
omnivore_anova <- aov(omnivore_trans ~ restriction, MIS_data = MIS_data)
summary(omnivore_anova)
TukeyHSD(omnivore_anova)

kruskal.test(omnivore ~ restriction, data = MIS_data)
pairwise.wilcox.test(MIS_data$omnivore, MIS_data$restriction, p.adjust.method = "BH")

#planktivores
planktivore_anova <- aov(planktivore ~ restriction, data = MIS_data)
summary(planktivore_anova)
TukeyHSD(planktivore_anova)

kruskal.test(planktivore ~ restriction, data = MIS_data)
pairwise.wilcox.test(MIS_data$planktivore, MIS_data$restriction, p.adjust.method = "BH")

#corallivore
corallivore_anova <- aov(corallivore_trans ~ restriction, data = MIS_data)
summary(corallivore_anova)
TukeyHSD(corallivore_anova)

kruskal.test(coralivore ~ restriction, data = MIS_data)
pairwise.wilcox.test(MIS_data$coralivore, MIS_data$restriction, p.adjust.method = "BH")

#protists
protist_anova <- aov(protist ~ restriction, data = MIS_data)
summary(protist_anova)
TukeyHSD(protist_anova)

kruskal.test(protist ~ restriction, data = MIS_data)
pairwise.wilcox.test(MIS_data$protist, MIS_data$restriction, p.adjust.method = "BH")

#Rerun stats from above on individual regions
#test each column of WGO_data for normality (p > 0.05 = normal)
#carnivores
hist(WGO_data$carnivore, main = "carnivore", breaks = 10)
shapiro.test(WGO_data$carnivore)
#Normal without transformation

#herbivores
hist(WGO_data$herbivore, main = "herbivore", breaks = 10)
shapiro.test(WGO_data$herbivore)

#transform non-normal WGO_data and re-test for normality
WGO_data$herbivore_trans <- 1/sqrt(WGO_data$herbivore)
hist(WGO_data$herbivore_trans, main = "herbivore transformed", breaks = 10)
shapiro.test(WGO_data$herbivore_trans)

#benthic invertivores
hist(WGO_data$benthic.invertivore, main = "benthic invertivore", breaks = 10)
shapiro.test(WGO_data$benthic.invertivore)

WGO_data$benthic.invertivore_trans <- log10(WGO_data$benthic.invertivore)
hist(WGO_data$benthic.invertivore_trans, main = "benthic invertivore transformed", breaks = 10)
shapiro.test(WGO_data$benthic.invertivore_trans)

#omnivores
hist(WGO_data$omnivore, main = "omnivore", breaks = 10)
shapiro.test(WGO_data$omnivore)
#Normal without transformation

#planktivores
hist(WGO_data$planktivore, main = "planktivore", breaks = 10)
shapiro.test(WGO_data$planktivore)
#Normal without transformation

#corallivores
hist(WGO_data$coralivore, main = "corallivore", breaks = 10)
shapiro.test(WGO_data$coralivore)

WGO_data$corallivore_trans <- sqrt(WGO_data$coralivore)
hist(WGO_data$corallivore_trans, main = "corallivore transformed", breaks = 10)
shapiro.test(WGO_data$corallivore_trans)

#protists
hist(WGO_data$protist, main = "protist OTU richness", breaks = 10)
shapiro.test(WGO_data$protist)
#Normal without transformation

#Test for significance across fishing restrictions using ANOVA (parametric with transformed WGO_data)
#Kruskal-Wallis is also used as non-parametric on non-transformed WGO_data for consistency with
#other statistical comparisons in the manuscript
#If ANOVA is significant use Tukey to run pairwise comparisons (when using Kruskal-Wallis,
#use Wilcoxon for pairwise comparisons)
#carnivores
carnivore_anova <- aov(carnivore_trans ~ restriction, data = WGO_data)
summary(carnivore_anova)
TukeyHSD(carnivore_anova)

kruskal.test(carnivore ~ restriction, data = WGO_data)
pairwise.wilcox.test(WGO_data$carnivore, WGO_data$restriction, p.adjust.method = "BH")

#herbivores
herbivore_anova <- aov(herbivore_trans ~ restriction, data = WGO_data)
summary(herbivore_anova)
TukeyHSD(herbivore_anova)

kruskal.test(herbivore ~ restriction, data = WGO_data)
pairwise.wilcox.test(WGO_data$herbivore, WGO_data$restriction, p.adjust.method = "BH")

#benthic invertivores
benthic.invertivore_anova <- aov(benthic.invertivore_trans ~ restriction, data = WGO_data)
summary(benthic.invertivore_anova)
TukeyHSD(benthic.invertivore_anova)

kruskal.test(benthic.invertivore ~ restriction, data = WGO_data)
pairwise.wilcox.test(WGO_data$benthic.invertivore, WGO_data$restriction, p.adjust.method = "BH")

#omnivores
omnivore_anova <- aov(omnivore ~ restriction, data = WGO_data)
summary(omnivore_anova)
TukeyHSD(omnivore_anova)

kruskal.test(omnivore ~ restriction, data = WGO_data)
pairwise.wilcox.test(WGO_data$omnivore, WGO_data$restriction, p.adjust.method = "BH")

#planktivores
planktivore_anova <- aov(planktivore ~ restriction, data = WGO_data)
summary(planktivore_anova)
TukeyHSD(planktivore_anova)

kruskal.test(planktivore ~ restriction, data = WGO_data)
pairwise.wilcox.test(WGO_data$planktivore, WGO_data$restriction, p.adjust.method = "BH")

#corallivore
corallivore_anova <- aov(corallivore_trans ~ restriction, data = WGO_data)
summary(corallivore_anova)
TukeyHSD(corallivore_anova)

kruskal.test(coralivore ~ restriction, data = WGO_data)
pairwise.wilcox.test(WGO_data$coralivore, WGO_data$restriction, p.adjust.method = "BH")

#protists
protist_anova <- aov(protist ~ restriction, data = WGO_data)
summary(protist_anova)
TukeyHSD(protist_anova)

kruskal.test(protist ~ restriction, data = WGO_data)
pairwise.wilcox.test(WGO_data$protist, WGO_data$restriction, p.adjust.method = "BH")

#transform data frame to long format for plotting purposes
data <- data[,c(1:2,4:9,11)]
data <- data[,-c(5:6)]
data$restriction <- factor(data$restriction, levels = c("unrestricted fishing",
                                                        "fishing gear restriction",
                                                        "fishing prohibited"))
data_long <- melt(data)

#Define colors
colors <- c('fishing prohibited' = "#f73e42",
            'fishing gear restriction' = "#8e86c5",
            'unrestricted fishing' = "#007cb5")
#plot boxplots for abundance/richness by fishing gear restriction
plot <- ggplot(data_long, aes(restriction, value, fill = restriction)) + geom_boxplot() +
  facet_wrap(~variable, scales = "free", nrow = 3, ncol = 2) + theme_bw() +
  theme(strip.background = element_rect(fill = "black"), strip.text = element_text(face = "bold", color = "white")) +
  theme(panel.grid = element_blank()) + scale_fill_manual(values = colors)
plot

ggsave("../Figure5.pdf", plot, height = 8, width = 7.5)
