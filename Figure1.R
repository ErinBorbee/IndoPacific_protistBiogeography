#Load packages
library(maptools)
library(raster)
library(maps)
library(GISTools)
library(prettymapr)

#Set working directory
setwd("~/Desktop/Dissertation/")

#Load coordinate data and remove duplicate site rows
meta <- read.csv("~/Desktop/DataAnalysis/Master_V9_MappingFile.csv")
meta <- meta[!duplicated(meta$siteID),]
#Subset sites by management gear restriction
noTake_meta <- subset(meta, restriction == "fishing prohibited")
restricted_meta <- subset(meta, restriction == "fishing gear restriction")
openAccess_meta <- subset(meta, restriction == "unrestricted fishing")

#Load country shapefile data
Indo <- getData('GADM', country ='IDN', level = 0)

#Plot Lombok map with site points (each management with different shape)
plot(Indo, col = "grey57", border = "white", bg = "white", axes=T, xlim = c(115.5,117), ylim = c(-9, -8))
points(noTake_meta$long, noTake_meta$lat, pch = 21, bg = "black", col = "white", cex = 3, lwd = 1)
points(restricted_meta$long, restricted_meta$lat, pch = 22, bg = "black", col = "white", cex = 3, lwd = 1)
points(openAccess_meta$long, openAccess_meta$lat, pch = 24, bg = "black", col = "white", cex = 3, lwd = 1)
LBK_map <- recordPlot()

#Add scalebar and save
plot.new()
png("~/Desktop/LBK_Map_with_sites.png", height = 4000, width = 4000, res = 300)
LBK_map

addscalebar(plotunit = NULL, plotepsg = NULL, widthhint = 0.25,
            unitcategory = "metric", htin = 0.1, padin = c(2, 2),
            style = "bar", bar.cols = c("black", "white"), lwd = 0.5,
            linecol = "grey20", tick.cex = 0.5, labelpadin = 0.1, label.cex = 1,
            label.col = "grey20", pos = "bottomleft")
dev.off()

#Plot Wakatobi map with site points
plot(Indo, col = "grey57", border = NA, bg = "white", axes=T, xlim = c(123,124.5), ylim = c(-6.3, -5.1))
points(noTake_meta$long, noTake_meta$lat, pch = 21, bg = "black", col = "white", cex = 15, lwd = 5)
points(restricted_meta$long, restricted_meta$lat, pch = 22, bg = "black", col = "white", cex = 15, lwd = 5)
points(openAccess_meta$long, openAccess_meta$lat, pch = 24, bg = "black", col = "white", cex = 15, lwd = 5)
WK_map <- recordPlot()

#Add scalebar and save
plot.new()
png("WK_Map_with_sites.png", height = 4000, width = 4000)
WK_map

addscalebar(plotunit = NULL, plotepsg = NULL, widthhint = 0.25,
            unitcategory = "metric", htin = 0.75, padin = c(2, 2),
            style = "bar", bar.cols = c("grey20", "white"), lwd = 6,
            linecol = "grey20", tick.cex = 0.7, labelpadin = 1, label.cex = 8,
            label.col = "grey20", pos = "bottomleft")
dev.off()

#Plot Waigeo map with site points
plot(Indo, col = "grey57", border = NA, bg = "white", axes=T, xlim = c(130,131.5), ylim = c(-1, 0.25))
points(noTake_meta$long, noTake_meta$lat, pch = 21, bg = "black", col = "white", cex = 15, lwd = 5)
points(restricted_meta$long, restricted_meta$lat, pch = 22, bg = "black", col = "white", cex = 15, lwd = 5)
points(openAccess_meta$long, openAccess_meta$lat, pch = 24, bg = "black", col = "white", cex = 15, lwd = 5)
WGO_map <- recordPlot()

#Add scalebar and save
plot.new()
png("WGO_Map_with_sites.png", height = 4000, width = 4000)
WGO_map

addscalebar(plotunit = NULL, plotepsg = NULL, widthhint = 0.25,
            unitcategory = "metric", htin = 0.75, padin = c(2, 2),
            style = "bar", bar.cols = c("grey20", "white"), lwd = 6,
            linecol = "grey20", tick.cex = 0.7, labelpadin = 1, label.cex = 8,
            label.col = "grey20", pos = "bottomleft")
dev.off()

#Plot Misool map and site points
plot(Indo, col = "grey57", border = NA, bg = "white", axes=T, xlim = c(129.5,131.5), ylim = c(-2.3, -1.4))
points(noTake_meta$long, noTake_meta$lat, pch = 21, bg = "black", col = "white", cex = 15, lwd = 5)
points(restricted_meta$long, restricted_meta$lat, pch = 22, bg = "black", col = "white", cex = 15, lwd = 5)
points(openAccess_meta$long, openAccess_meta$lat, pch = 24, bg = "black", col = "white", cex = 15, lwd = 5)
MIS_map <- recordPlot()

#Add scale bar to map and save
plot.new()
png("MIS_Map_with_sites.png", height = 4000, width = 4000)
MIS_map

addscalebar(plotunit = NULL, plotepsg = NULL, widthhint = 0.25,
            unitcategory = "metric", htin = 0.75, padin = c(2, 2),
            style = "bar", bar.cols = c("grey20", "white"), lwd = 6,
            linecol = "grey20", tick.cex = 0.7, labelpadin = 1, label.cex = 8,
            label.col = "grey20", pos = "bottomleft")
dev.off()
