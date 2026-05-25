# post beast 
install.packages("devtools")
library(devtools)
install_github("sdellicour/seraphim/unix_OS")
library(seraphim)
install.packages("diagram")
library(diagram)

#création du raster 

setwd("/Users/salomeembise/Desktop/épidémio/TP/Binome_03/génomique/ProjetBinome3-Arbovirus-Risque/doc pour beauti ") 

areastud = read.table("EBOV_OK_loc.tsv", header=T, sep="\t")
xmin = min(areastud$longitude)
xmax = max(areastud$longitude)
ymin = min(areastud$latitude)
ymax = max(areastud$latitude)


margin = 0.1 # ajout d'une marge
ext = extent(xmin - margin, xmax + margin, ymin - margin, ymax + margin)

#raster
study_area = raster(ext, resolution=c(0.5, 0.5))
values(study_area) = 1
writeRaster(study_area, "Study_area.asc", format="ascii")

setwd("/Users/salomeembise/Desktop/épidémio/TP/Binome_03/génomique/ProjetBinome3-Arbovirus-Risque")

library(seraphim)

#Extraction d'informations spatio-temporelles intégrées dans l'arbre MCC 
#renommer location rate dans .tree pour etre reconnu par la fn mccTreeExtractions
tree_txt=readLines("EBOV_MCC.tree")
tree_txt=gsub("Location.rate1", "location2", tree_txt)
tree_txt=gsub("Location.rate2", "location1", tree_txt)
writeLines(tree_txt, "EBOV_MCC_OK.tree")

#suite du protocole
mcc_tre = readAnnotatedNexus("EBOV_MCC_OK.tree")
mostRecentSamplingDatum = 2025.00
mcc_tab = mccTreeExtractions(mcc_tre, mostRecentSamplingDatum)

# Prmblm cas où startYear > endYear donc:
for (i in 1:nrow(mcc_tab)) {
  if (mcc_tab[i, "startYear"] > mcc_tab[i, "endYear"]) {
    # Inverse TOUT : années ET coordonnées
    mcc_tab[i, c("startYear", "endYear")] = mcc_tab[i, c("endYear", "startYear")]
    mcc_tab[i, c("startLon", "endLon")] = mcc_tab[i, c("endLon", "startLon")]
    mcc_tab[i, c("startLat", "endLat")] = mcc_tab[i, c("endLat", "startLat")]
    mcc_tab[i, c("startNodeL", "endNodeL")] = mcc_tab[i, c("endNodeL", "startNodeL")]
  }
}

write.csv(mcc_tab, "EBOV_MCC.csv", row.names=F, quote=F)
mcc_tab = read.csv("EBOV_MCC.csv", head=T)


#Extraction d'informations spatio-temporelles intégrées dans les arbres postérieurs
localTreesDirectory = "Extracted_trees" 
allTrees = scan(file="EBOV_OK.trees", what="", sep="\n", quiet=T) 
burnIn = 0
randomSampling = FALSE 
nberOfTreesToSample = 10000  
coordinateAttributeName = "Location" 
treeExtractions(localTreesDirectory, allTrees, burnIn, randomSampling, nberOfTreesToSample, mostRecentSamplingDatum, coordinateAttributeName)

#Estimation de la région HPD pour chaque tranche de temps
library(seraphim)

nberOfExtractionFiles = nberOfTreesToSample
prob = 0.80
startDatum = (2024 + (36/365))
precision = 1/12

polygons = suppressWarnings(spreadGraphic2(localTreesDirectory,
                                           nberOfExtractionFiles,
                                           prob,
                                           startDatum,
                                           precision))

#correction des coordonnées : inversion longitude et latitude
polygons = lapply(polygons, function(poly) {
  for (p in seq_along(poly@polygons)) {
    for (r in seq_along(poly@polygons[[p]]@Polygons)) {
      coords = poly@polygons[[p]]@Polygons[[r]]@coords
      slot(slot(slot(poly, "polygons")[[p]], "Polygons")[[r]], "coords") = coords[, c(2, 1)]
    }
  }
  return(poly)
})

saveRDS(polygons, "polygons.rds")

# Forcer la sauvegarde en écrasant
saveRDS(polygons_corriges, file = "polygons.rds", compress = FALSE)

# Vérification immédiate en relisant le fichier
polygons_relu = readRDS("polygons.rds")
head(polygons_relu[[1]]@polygons[[1]]@Polygons[[1]]@coords)


#Définir les différentes échelles de couleurs à utiliser
colour_scale = colorRampPalette(brewer.pal(11,"RdYlGn"))(141)[21:121]
minYear = 2024; maxYear = mostRecentSamplingDatum 
endYears_indices = (((mcc_tab[,"endYear"]-minYear)/(maxYear-minYear))*100)+1
endYears_colours = colour_scale[endYears_indices]
polygons_colours = rep(NA, length(polygons))
for (i in 1:length(polygons)) {
  date = as.numeric(names(polygons[[i]]))
  polygon_index = round((((date-minYear)/(maxYear-minYear))*100)+1)
  polygons_colours[i] = paste0(colour_scale[polygon_index],"40")
  }
saveRDS(polygons_colours, "polygons_colours.rds")

#Co-tracer les régions HPD et l'arbre MCC
mcc_tab = read.csv("EBOV_MCC.csv", head=T)
polygons = readRDS("polygons.rds")
head(polygons[[1]]@polygons[[1]]@Polygons[[1]]@coords)
polygons_colours = readRDS("polygons_colours.rds")

setwd("/Users/salomeembise/Desktop/épidémio/TP/Binome_03/génomique/ProjetBinome3-Arbovirus-Risque")

library(seraphim)
library(diagram)
library(raster)
library(sf)
library(RColorBrewer)


guinea = st_read("gadm41_GIN_shp/gadm41_GIN_0.shp")
liberia = st_read("gadm41_LBR_shp/gadm41_LBR_0.shp")
sierra_leone = st_read("gadm41_SLE_shp/gadm41_SLE_0.shp")
cotediv= st_read("gadm41_CIV_shp/gadm41_CIV_0.shp")
mali= st_read("gadm41_MLI_shp/gadm41_MLI_0.shp")
senegal= st_read("gadm41_SEN_shp/gadm41_SEN_0.shp")
guineebis= st_read("gadm41_GNB_shp/gadm41_GNB_0.shp")
gambie= st_read("gadm41_GMB_shp/gadm41_GMB_0.shp")
burk= st_read("gadm41_BFA_shp/gadm41_BFA_0.shp")

#combine
borders = rbind(guinea, liberia, sierra_leone,cotediv,mali,senegal,guineebis,gambie, burk)

#raster
template_raster = raster("Study_area.asc")

#bordures raster
borders = st_crop(borders, st_as_sfc(st_bbox(extent(template_raster))))
borders = st_simplify(borders, dTolerance = 0.01)

####VISUALISATION
dev.off()
dev.new(width=6, height=6.3)
par(mar=c(0,0,0,0), oma=c(1.2,3.5,1,0), mgp=c(0,0.4,0), lwd=0.2, bty="o")

#fond blanc
plot(template_raster, col="grey80", box=F, axes=F, colNA="grey90", legend=F)# Lis chaque shapefile


plot(borders, add=T, lwd=0.8, border="gray30", col="white")

for (i in 1:length(polygons)) {
  plot(polygons[[i]], axes=F, col=polygons_colours[i], add=T, border=NA)
}

for (i in 1:dim(mcc_tab)[1]) {
  curvedarrow(cbind(mcc_tab[i,"startLon"], mcc_tab[i,"startLat"]),
              cbind(mcc_tab[i,"endLon"], mcc_tab[i,"endLat"]), 
              arr.length=0, arr.width=0, lwd=0.2, lty=1, 
              lcol="gray10", arr.col=NA, arr.pos=F, curve=0.1, 
              dr=NA, endhead=F)
}

for (i in 1:dim(mcc_tab)[1]) {
  if (i == 1) {
    points(mcc_tab[i,"startLon"], mcc_tab[i,"startLat"], pch=16,
           col=colour_scale[1], cex=0.8)
    points(mcc_tab[i,"startLon"], mcc_tab[i,"startLat"], pch=1,
           col="gray10", cex=0.8)
  }
  points(mcc_tab[i,"endLon"], mcc_tab[i,"endLat"], pch=16,
         col=endYears_colours[i], cex=0.8)
  points(mcc_tab[i,"endLon"], mcc_tab[i,"endLat"], pch=1,
         col="gray10", cex=0.8)
}

#cadre
rect(xmin(template_raster), ymin(template_raster), 
     xmax(template_raster), ymax(template_raster), 
     xpd=T, lwd=0.2)

#axe X 
axis(1, 
     at = c(ceiling(xmin(template_raster)), floor(xmax(template_raster))),
     labels = c(paste(abs(ceiling(xmin(template_raster))), "°W"),
                paste(abs(floor(xmax(template_raster))), "°W")),
     pos=ymin(template_raster), 
     mgp=c(0,0.4,0), cex.axis=0.7, lwd=0, lwd.tick=0.2,
     padj=-0.8, tck=-0.01, col.axis="gray30")

#axe Y
axis(2, 
     at = c(ceiling(ymin(template_raster)), floor(ymax(template_raster))),
     labels = c(paste(ceiling(ymin(template_raster)), "°N"),
                paste(floor(ymax(template_raster)), "°N")),
     pos=xmin(template_raster), 
     mgp=c(0,0.6,0), cex.axis=0.7, lwd=0, lwd.tick=0.2,
     padj=1, tck=-0.01, col.axis="gray30")

#légende
rast = raster(matrix(nrow=1, ncol=2))
rast[1] = min(mcc_tab[,"startYear"])
rast[2] = max(mcc_tab[,"endYear"])

plot(rast, legend.only=TRUE, add=TRUE, col=colour_scale,
     legend.width=0.25,
     legend.shrink=0.18,
     smallplot=c(0.049, 0.50, 0.105, 0.125),
     horizontal=TRUE,
     legend.args=list(text="", cex=0.7, line=0.1, col="gray30"),
     axis.args=list(cex.axis=0.7, lwd=0, lwd.tick=0.2, tck=-0.4,
                    col.axis="gray30", line=-0.2, mgp=c(0, -0.1, 0),
                    at=seq(min(mcc_tab[,"startYear"]), max(mcc_tab[,"endYear"]), 0.2)))


setwd("/Users/salomeembise/Desktop/épidémio/TP/Binome_03/génomique/ProjetBinome3-Arbovirus-Risque")

localTreesDirectory = "Extracted_trees" 
nberOfExtractionFiles = 1000
timeSlices = 100
onlyTipBranches = FALSE
showingPlots = T
outputName = "EBOV"
nberOfCores = 7
slidingWindow = 0.1

spreadStatistics(localTreesDirectory, nberOfExtractionFiles, timeSlices,
                 onlyTipBranches, showingPlots, outputName, nberOfCores, 
                 slidingWindow)

#Ne(t)
log_data = read.table("EBOV_OK.log", header=TRUE)

skyline_moy=rev(colMeans(log_data[, grep("skygrid.logPopSize", colnames(log_data))]))

#taille efficace: Ne = exp(logPopSize)
Ne = exp(skyline_moy)

#graphe Ne(t)
time_points= seq(2024 + (26/365), 2025, length.out=47)

plot(time_points, Ne, type="l", lwd=2, 
     xlab="t(année)", ylab="Ne",
     main="Ne(t)")

#graphe R0(t)
#param du calcul
cut_off=0.95
num_intervals=47
delta_t=cut_off / num_intervals  # Durée d'un intervalle 
D =21 / 365 # Durée EBOLA

#diff entre chaque log
diff_logs=diff(skyline_moy)

#r et R0
r=diff_logs / delta_t
r0= 1 + (r * D)

# plot
time_points <- seq(2024 + (26/365), 2025, length.out = length(r0))

plot(time_points, r0, 
     type = "o",    
     pch = 16,           
     cex = 0.6,            
     col = "black",      
     lwd = 1,       
     xlab = "t(année)", 
     ylab = "R0 ",
     main = "R0(t)")
abline(h = 1, col = "red", lty = 2, lwd = 1.5)


