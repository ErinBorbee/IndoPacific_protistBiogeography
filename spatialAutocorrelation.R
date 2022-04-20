library(vegan)
library(ade4)
library(ggplot2)

setwd("~/Desktop/DataAnalysis/threshold95/SAR/")

meta <- read.csv("~/Desktop/DataAnalysis/Master_V9_MappingFile.csv")
table <- read.csv("SAR-percent-table-95.csv", row.names = 1)
table <- table[,-c(1:8)]
table <- as.data.frame(t(table))

meta <- subset(meta, meta$sampleID %in% rownames(table))
metaM04 <- subset(meta, subset == "M04")
metaM12 <- subset(meta, subset == "M12")
metaS04 <- subset(meta, subset == "S04")
metaS12 <- subset(meta, subset == "S12")

metaM04_LBK <- subset(metaM04, region == "Lombok")
metaM04_WK <- subset(metaM04, region == "Wakatobi")
metaM04_MIS <- subset(metaM04, region == "Misool")
metaM04_WGO <- subset(metaM04, region == "Waigeo")

metaM12_LBK <- subset(metaM12, region == "Lombok")
metaM12_WK <- subset(metaM12, region == "Wakatobi")
metaM12_MIS <- subset(metaM12, region == "Misool")
metaM12_WGO <- subset(metaM12, region == "Waigeo")

metaS04_LBK <- subset(metaS04, region == "Lombok")
metaS04_WK <- subset(metaS04, region == "Wakatobi")
metaS04_MIS <- subset(metaS04, region == "Misool")
metaS04_WGO <- subset(metaS04, region == "Waigeo")

metaS12_LBK <- subset(metaS12, region == "Lombok")
metaS12_WK <- subset(metaS12, region == "Wakatobi")
metaS12_MIS <- subset(metaS12, region == "Misool")
metaS12_WGO <- subset(metaS12, region == "Waigeo")

tableM04 <- subset(table, rownames(table) %in% metaM04$sampleID)
tableM12 <- subset(table, rownames(table) %in% metaM12$sampleID)
tableS04 <- subset(table, rownames(table) %in% metaS04$sampleID)
tableS12 <- subset(table, rownames(table) %in% metaS12$sampleID)

tableM04_LBK <- subset(tableM04, rownames(tableM04) %in% metaM04_LBK$sampleID)
tableM12_LBK <- subset(tableM12, rownames(tableM12) %in% metaM12_LBK$sampleID)
tableS04_LBK <- subset(tableS04, rownames(tableS04) %in% metaS04_LBK$sampleID)
tableS12_LBK <- subset(tableS12, rownames(tableS12) %in% metaS12_LBK$sampleID)

tableM04_WK <- subset(tableM04, rownames(tableM04) %in% metaM04_WK$sampleID)
tableM12_WK <- subset(tableM12, rownames(tableM12) %in% metaM12_WK$sampleID)
tableS04_WK <- subset(tableS04, rownames(tableS04) %in% metaS04_WK$sampleID)
tableS12_WK <- subset(tableS12, rownames(tableS12) %in% metaS12_WK$sampleID)

tableM04_MIS <- subset(tableM04, rownames(tableM04) %in% metaM04_MIS$sampleID)
tableM12_MIS <- subset(tableM12, rownames(tableM12) %in% metaM12_MIS$sampleID)
tableS04_MIS <- subset(tableS04, rownames(tableS04) %in% metaS04_MIS$sampleID)
tableS12_MIS <- subset(tableS12, rownames(tableS12) %in% metaS12_MIS$sampleID)

tableM04_WGO <- subset(tableM04, rownames(tableM04) %in% metaM04_WGO$sampleID)
tableM12_WGO <- subset(tableM12, rownames(tableM12) %in% metaM12_WGO$sampleID)
tableS04_WGO <- subset(tableS04, rownames(tableS04) %in% metaS04_WGO$sampleID)
tableS12_WGO <- subset(tableS12, rownames(tableS12) %in% metaS12_WGO$sampleID)

#Spatial autocorrelation in full dataset
siteDist <- vegdist(cbind(meta$long,meta$lat), method = "bray")
sampleDist <- vegdist(table, method = "bray")

mantel.rtest(siteDist, sampleDist, nrepet = 9999)

#Spatial autocorrelation in M04 samples
siteDist_M04 <- vegdist(cbind(metaM04$long,metaM04$lat), method = "bray")
sampleDist_M04 <- vegdist(tableM04, method = "bray")

mantel.rtest(siteDist_M04, sampleDist_M04, nrepet = 9999)

#Spatial autocorrelation in M12 samples
siteDist_M12 <- vegdist(cbind(metaM12$long,metaM12$lat), method = "bray")
sampleDist_M12 <- vegdist(tableM12, method = "bray")

mantel.rtest(siteDist_M12, sampleDist_M12, nrepet = 9999)

#Spatial autocorrelation in S04 samples
siteDist_S04 <- vegdist(cbind(metaS04$long,metaS04$lat), method = "bray")
sampleDist_S04 <- vegdist(tableS04, method = "bray")

mantel.rtest(siteDist_S04, sampleDist_S04, nrepet = 9999)

#Spatial autocorrelation in S12 samples
siteDist_S12 <- vegdist(cbind(metaS12$long,metaS12$lat), method = "bray")
sampleDist_S12 <- vegdist(tableS12, method = "bray")

mantel.rtest(siteDist_S12, sampleDist_S12, nrepet = 9999)

#Spatial autocorrelation in LBK M04 samples
siteDist_M04_LBK <- vegdist(cbind(metaM04_LBK$long,metaM04_LBK$lat), method = "bray")
sampleDist_M04_LBK <- vegdist(tableM04_LBK, method = "bray")

mantel.rtest(siteDist_M04_LBK, sampleDist_M04_LBK, nrepet = 9999)

#Spatial autocorrelation in LBK M12 samples
siteDist_M12_LBK <- vegdist(cbind(metaM12_LBK$long,metaM12_LBK$lat), method = "bray")
sampleDist_M12_LBK <- vegdist(tableM12_LBK, method = "bray")

mantel.rtest(siteDist_M12_LBK, sampleDist_M12_LBK, nrepet = 9999)

#Spatial autocorrelation in LBK S04 samples
siteDist_S04_LBK <- vegdist(cbind(metaS04_LBK$long,metaS04_LBK$lat), method = "bray")
sampleDist_S04_LBK <- vegdist(tableS04_LBK, method = "bray")

mantel.rtest(siteDist_S04_LBK, sampleDist_S04_LBK, nrepet = 9999)

#Spatial autocorrelation in LBK S12 samples
siteDist_S12_LBK <- vegdist(cbind(metaS12_LBK$long,metaS12_LBK$lat), method = "bray")
sampleDist_S12_LBK <- vegdist(tableS12_LBK, method = "bray")

mantel.rtest(siteDist_S12_LBK, sampleDist_S12_LBK, nrepet = 9999)

#Spatial autocorrelation in WK M04 samples
siteDist_M04_WK <- vegdist(cbind(metaM04_WK$long,metaM04_WK$lat), method = "bray")
sampleDist_M04_WK <- vegdist(tableM04_WK, method = "bray")

mantel.rtest(siteDist_M04_WK, sampleDist_M04_WK, nrepet = 9999)

#Spatial autocorrelation in WK M12 samples
siteDist_M12_WK <- vegdist(cbind(metaM12_WK$long,metaM12_WK$lat), method = "bray")
sampleDist_M12_WK <- vegdist(tableM12_WK, method = "bray")

mantel.rtest(siteDist_M12_WK, sampleDist_M12_WK, nrepet = 9999)

#Spatial autocorrelation in WK S04 samples
siteDist_S04_WK <- vegdist(cbind(metaS04_WK$long,metaS04_WK$lat), method = "bray")
sampleDist_S04_WK <- vegdist(tableS04_WK, method = "bray")

mantel.rtest(siteDist_S04_WK, sampleDist_S04_WK, nrepet = 9999)

#Spatial autocorrelation in WK S12 samples
siteDist_S12_WK <- vegdist(cbind(metaS12_WK$long,metaS12_WK$lat), method = "bray")
sampleDist_S12_WK <- vegdist(tableS12_WK, method = "bray")

mantel.rtest(siteDist_S12_WK, sampleDist_S12_WK, nrepet = 9999)

#Spatial autocorrelation in MIS M04 samples
siteDist_M04_MIS <- vegdist(cbind(metaM04_MIS$long,metaM04_MIS$lat), method = "bray")
sampleDist_M04_MIS <- vegdist(tableM04_MIS, method = "bray")

mantel.rtest(siteDist_M04_MIS, sampleDist_M04_MIS, nrepet = 9999)

#Spatial autocorrelation in MIS M12 samples
siteDist_M12_MIS <- vegdist(cbind(metaM12_MIS$long,metaM12_MIS$lat), method = "bray")
sampleDist_M12_MIS <- vegdist(tableM12_MIS, method = "bray")

mantel.rtest(siteDist_M12_MIS, sampleDist_M12_MIS, nrepet = 9999)

#Spatial autocorrelation in MIS S04 samples
siteDist_S04_MIS <- vegdist(cbind(metaS04_MIS$long,metaS04_MIS$lat), method = "bray")
sampleDist_S04_MIS <- vegdist(tableS04_MIS, method = "bray")

mantel.rtest(siteDist_S04_MIS, sampleDist_S04_MIS, nrepet = 9999)

#Spatial autocorrelation in MIS S12 samples
siteDist_S12_MIS <- vegdist(cbind(metaS12_MIS$long,metaS12_MIS$lat), method = "bray")
sampleDist_S12_MIS <- vegdist(tableS12_MIS, method = "bray")

mantel.rtest(siteDist_S12_MIS, sampleDist_S12_MIS, nrepet = 9999)

#Spatial autocorrelation in WGO M04 samples
siteDist_M04_WGO <- vegdist(cbind(metaM04_WGO$long,metaM04_WGO$lat), method = "bray")
sampleDist_M04_WGO <- vegdist(tableM04_WGO, method = "bray")

mantel.rtest(siteDist_M04_WGO, sampleDist_M04_WGO, nrepet = 9999)

#Spatial autocorrelation in WGO M12 samples
siteDist_M12_WGO <- vegdist(cbind(metaM12_WGO$long,metaM12_WGO$lat), method = "bray")
sampleDist_M12_WGO <- vegdist(tableM12_WGO, method = "bray")

mantel.rtest(siteDist_M12_WGO, sampleDist_M12_WGO, nrepet = 9999)

#Spatial autocorrelation in WGO S04 samples
siteDist_S04_WGO <- vegdist(cbind(metaS04_WGO$long,metaS04_WGO$lat), method = "bray")
sampleDist_S04_WGO <- vegdist(tableS04_WGO, method = "bray")

mantel.rtest(siteDist_S04_WGO, sampleDist_S04_WGO, nrepet = 9999)

#Spatial autocorrelation in WGO S12 samples
siteDist_S12_WGO <- vegdist(cbind(metaS12_WGO$long,metaS12_WGO$lat), method = "bray")
sampleDist_S12_WGO <- vegdist(tableS12_WGO, method = "bray")

mantel.rtest(siteDist_S12_WGO, sampleDist_S12_WGO, nrepet = 9999)
