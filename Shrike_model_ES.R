###################################################
#
# Title: Shrike_model_ES.R
# Purpose: run a Species Distribution Model for the red-backed shrike in Catalonia, Spain, to quantify the effect of grassland-based Agri-Environmental Schemes (AES)
#          on the shrike's occurrence, and to test whether this effect is moderated by landscape complexity. 
# Author: Stephanie Roilo, Dresden University of Technology (TUD)
#         Last updated on the 3rd of August 2023
#
###################################################

# set the language to EN
Sys.setenv(LANGUAGE="en")
# load packages
library(sf)
library(mapview)
library(dplyr)
library(psych)   # to plot nice pairs.plots 
library(MuMIn)
library(visreg)
library(DHARMa)  # to simulate residuals and check GLM assumptions
library(ncf)  # to check for spatial autocorrelation in the model residuals
library(data.table)
library(ggplot2)
library(interflex) # to check whether interaction terms in GLM are trustworthy (i.e. enough common support and linear interactive effect)
library(gridExtra)

# DATA PREPARATION ------------------------------
# load Catalan bird data - this grid contains the proportions of cover calculated from the mosaiced raster "Rast_IACS_ES2019"
# with rasterized land cover information from the IACS/LPIS data of 2019 on top of which the Small Woody Features were overlaid
gridES = st_read("C:/Users/sroilo/Desktop/BESTMAP documents/Data/Across_CS/ES/rasterES/Grid_ES_20221215.shp")
# the first columns are older data, let's get rid of those
gridES = gridES[,c(1,2, 11,12,14,16,20:45)]
names(gridES) <- c("CODI","TOTAL_C","Tmax","Precip","Elevation", "GrM", "NaN", "Arable", "orchard", "greenhouse", "fruit",
                   "citrus", "nuts", "olive", "vineyard", "unproductive", "forest",  "pasture with trees", "pasture with shrubs",
                   "pastureland", "buildings", "urban area", "vials",  "water sources", "olive and fruits", "vineyard and olive",
                   "vineyard and fruits", "nuts and vineyard",  "nuts and olive", "security zone",  "SWF", "geometry")
# compute grouped land-cover classes
gridES$Grass = gridES$pastureland
gridES$Wooded = rowSums(st_drop_geometry(gridES[,c("forest", "pasture with trees", "pasture with shrubs")]), na.rm=T)
gridES$PermCult = rowSums(st_drop_geometry(gridES[,c("orchard","fruit", "citrus", "nuts", "olive", "vineyard","olive and fruits", "vineyard and olive",
                                                     "vineyard and fruits", "nuts and vineyard",  "nuts and olive")]), na.rm=T)
gridES$Urban = rowSums(st_drop_geometry(gridES[,c("greenhouse","buildings", "urban area", "vials","security zone")]), na.rm=T)
# compute land-cover diversity on the following land cover classes: 
grid = st_drop_geometry(gridES)
for ( i in seq_along(grid$CODI)) {
  prop = as.vector(unlist(grid[i,c("NaN", "Arable", "PermCult", "Urban", "forest","pasture with trees","pasture with shrubs", "pastureland","water sources","unproductive", "SWF")]))
  prop = prop[which(prop>0)]
  gridES$LCdiv[i] = ifelse(length(prop)>0, -sum(sapply(X = prop, FUN = function(x) {x * log(x)}, simplify = "vector")), 0)
}
mapview(gridES, zcol="LCdiv")
head(gridES)
# get rid of some columns
gridES = gridES[,c(1:6,8,31,33:37)]

# now, update temperature and precipitation to cover only May to July
r = rast(crs="epsg:3035",extent=ext(gridES), resolution=c(1000,1000))
tasmax = rast(c("X:/Stephanie/Data/Europe/Climate/CHELSA/Climatologies/CHELSA_tasmax_05_1981-2010_V.2.1.tif",
                "X:/Stephanie/Data/Europe/Climate/CHELSA/Climatologies/CHELSA_tasmax_06_1981-2010_V.2.1.tif",
                "X:/Stephanie/Data/Europe/Climate/CHELSA/Climatologies/CHELSA_tasmax_07_1981-2010_V.2.1.tif")) %>% project(r) %>% app(max)
precip = rast(c("X:/Stephanie/Data/Europe/Climate/CHELSA/Climatologies/CHELSA_pr_05_1981-2010_V.2.1.tif",
                "X:/Stephanie/Data/Europe/Climate/CHELSA/Climatologies/CHELSA_pr_06_1981-2010_V.2.1.tif",
                "X:/Stephanie/Data/Europe/Climate/CHELSA/Climatologies/CHELSA_pr_07_1981-2010_V.2.1.tif")) %>% project(r) %>% app(sum)
gridES$Tmax = extract(tasmax, vect(gridES), fun="mean")[,2]
gridES$Precip = extract(precip, vect(gridES), fun="mean")[,2]
# save to file
st_write(gridES, "C:/Users/sroilo/Desktop/BESTMAP documents/Data/Across_CS/ES/Grid_vars_env_20221216.shp")
# TOTAL_C is the sum of all land cover types; remove the ones significantly lower than 1 (i.e. for which land-cover information is incomplete)
gridES = gridES[gridES$TOTAL_C>0.99,]   # 859 entries removed
# merge bird data from the Global Biodiversity Information Facility (GBIF) for 2019
ESbirds = st_read("C:/Users/sroilo/Desktop/BESTMAP documents/Data/ES_Catalonia/Biodiversity/GBIF/Birds_ES_ICO.shp")
ESbirds = ESbirds[ESbirds$year==2019,]
ESbirds = st_join(ESbirds[,c("species", "eventDt", "year")], gridES[, c(1,3:13)])
mapview(ESbirds, zcol="Grass") # some NA, remove them
ESbirds = na.omit(ESbirds)
#save to file
ESbirds = cbind(ESbirds, st_coordinates(ESbirds))
st_write(ESbirds, "C:/Users/sroilo/Desktop/BESTMAP documents/Data/Across_CS/datasets_birds/ESbirds_2019_20221216.shp")

### MODELLING: data filtering ------------------------
# load the gridded dataset holding bird observations and environmental information
ESbirds = st_read("C:/Users/sroilo/Desktop/BESTMAP documents/Data/Across_CS/datasets_birds/ESbirds_2019_20221216.shp")
# Filter out grid cells with less than 1 ha of grassland cover, as this is the usual territory size of the red-backed shrike
# As Catalonia is topographically very heterogeneous, we also filter out grid cells that are very low of high in elevation, based on 
# the shrike's suitable elevation range from the book 'Atles dels Ocells nidificants de Catalunya'
grbir = ESbirds[which(ESbirds$Grass>=0.01 & ESbirds$Elevation>200 & ESbirds$Elevation<2000),]
# check spatial distribution of presence-absences 
# grbir$presence = ifelse(grbir$species=="Lanius collurio", 1, 0)
# plot(grbir[,"presence"], pch=20, pal=hcl.colors(2, "Blue-Yellow", rev=T))
# drop the geometry, clean up the dataframe
grbir = st_drop_geometry(grbir)
grbir = cbind(grbir[,16:17], grbir[,1:7], grbir[,c(9,11,12,13,10,15,8)])
names(grbir)[7:16] <- c("TMAX", "PRECIP", "ELEVATION", "ARABLE","GRASS", "FOREST","ORCHARD", "SHRUBS", "LANDDIV", "AES")
## make a dataframe to model Lanius collurio
lces = grbir
lces$presence = ifelse(lces$species=="Lanius collurio", 1, 0)
# get rid of duplicates in the presence and absence points
lces = lces[order(lces$presence, decreasing=T),]
lces = lces[-which(duplicated(lces[,c("X","Y")])),]  ; sum(lces$presence) #2236 points, 233 presences
# AES was calculated from shapefiles, while GRASS from raster data. Make sure that there are no major discrepancies
nrow(lces[which(lces$AES > lces$GRASS),]) # 2 grid cells in which AES > GRASS; correct 
lces$AES = ifelse(lces$AES > lces$GRASS, lces$GRASS, lces$AES)

# check whether there is enough common support to investigate interactions between AES and SHRUBS and LANDDIV
# see interflex vignette here: https://yiqingxu.org/packages/interflex/RGuide.html 
interflex(estimator = "raw",                  
          Y = "presence", #outcome variable
          D = "AES", # treatment variable
          X = "SHRUBS", # moderating variable
          data = lces, 
          treat.type="continuous",
          method="logit",
          ncols=3,theme.bw = TRUE, span=30)  
interflex(estimator = "raw",
          Y = "presence", #outcome variable
          D = "AES", # treatment variable
          X = "LANDDIV", # moderating variable
          data = lces, 
          treat.type="continuous",
          method="logit",
          ncols=3, theme.bw = TRUE, span=30)   

### TRIM the dataset to reach common support ------------------------------
# use three different thresholds, in a sensitivity analysis
quantile(lces$AES, probs = c(0.9, 0.95, 0.99))           # 90% -> 0.06500 | 95% -> 0.13150 | 99% -> 0.28265 
lces = lces[which(lces$AES <= 0.06500),]  # grid-cells left:         2016 |       2124     |         2213 
# repeat the common support check above via the interflex function and move to the next cutoff value for the trimming

# produce scatter plots to check for correlation among variables
lceso = lces[order(lces$presence, decreasing=F),]     
pairs.panels(lceso[,c("presence", "AES", "SHRUBS", "LANDDIV", "GRASS", "ARABLE", "FOREST", "ORCHARD", "TMAX", "PRECIP", "ELEVATION")], 
             bg=c("blue","yellow")[as.factor(lceso$presence)], pch=21,
             method = "spearman", hist.col = "darkorange", density = TRUE, ellipses = FALSE )
# normalise variables in advance of the modelling:
dat = cbind(lces[,1:6], lces[,17], lces[,7:16])
names(dat)[7] <-"presence"
dat[,8:17] = scale(dat[,8:17])
# check for correlation and create a matrix which returns TRUE if the absolute value of correlation coefficient is lower than 0.7
cmat <- abs(cor(dat[,8:17], method="spearman")) <= 0.7
cmat[!lower.tri(cmat, diag=F)] <- NA 

### MODEL SELECTION ------------
fit <- glm(presence ~ AES + ARABLE + GRASS + SHRUBS + LANDDIV + TMAX + PRECIP + ELEVATION + FOREST + ORCHARD + AES:LANDDIV + AES:SHRUBS, 
           data = dat, family="binomial", na.action="na.fail")
# dredge the full model
dr1 = dredge(fit, subset=cmat, extra = c("R^2", "adjR^2"))  
write.table(dr1, sep=";", dec=",", row.names=F, "C:/Users/sroilo/Desktop/BESTMAP documents/Data/Across_CS/ES/July_trials/LC_trim0132_20230714.csv")
# get variable importance of highly supported models (deltaAIC<2)
drsub = dr1[which(dr1$delta<2),]
imp = sw(drsub)
data.frame(imp)
write.table(data.frame(imp),  sep=";", dec=",", row.names=T, "C:/Users/sroilo/Desktop/BESTMAP documents/Data/Across_CS/ES/July_trials/LC_trim0132_weights_20230714.csv")

## HEAVIEST TRIMMING (90% quantile):
# fit best model:
best1 = glm(presence ~ AES + ELEVATION + FOREST + GRASS + LANDDIV + ORCHARD + SHRUBS + AES:SHRUBS,
            data = dat, family="binomial", na.action="na.fail")
summary(best1)  # Orchard and LCdiv non-significant
# diagnostic plots with DHARMa package
simulateResiduals(best1, plot=T)   
# check for spatial autocorrelation in the residuals
max.dist = max(dist(dat[,c("X","Y")])) * 2/3
scor = ncf::spline.correlog(x=best1$data$X, y = best1$data$Y, z=best1$residuals, resamp=100, xmax = max.dist)
plot(scor)  

# refit the model with the original (non-scaled) data to produce conditional plots
best1 = glm(presence ~ AES + ELEVATION + FOREST + GRASS + LANDDIV + ORCHARD + SHRUBS + AES:SHRUBS,
            data = lces, family="binomial", na.action="na.fail")
# plot conditional plots for each predictor and interaction term
p1 = visreg(best1, "AES", scale="response", gg=TRUE, line=list(col="black"),
       fill=list(fill="lightblue", alpha=0.5),
       points = list(pch=NA)) + theme_bw() + ylim(0,0.5) +
       ylab("Occurrence probability")
p2 = visreg(best1, "ELEVATION", scale="response", gg=TRUE, line=list(col="black"),
       fill=list(fill="lightblue", alpha=0.5),
       points = list(pch=NA)) + theme_bw() + ylim(0,0.5) +
       ylab("Occurrence probability")
p3 = visreg(best1, "FOREST", scale="response", gg=TRUE, line=list(col="black"),
            fill=list(fill="lightblue", alpha=0.5),
            points = list(pch=NA)) + theme_bw() + ylim(0,0.5) +
            ylab("Occurrence probability")
p4 = visreg(best1, "GRASS", scale="response", gg=TRUE, line=list(col="black"),
            fill=list(fill="lightblue", alpha=0.5),
            points = list(pch=NA)) + theme_bw() + ylim(0,0.5) +
            ylab("Occurrence probability")
p5 = visreg(best1, "LANDDIV", scale="response", gg=TRUE, line=list(col="black"),
            fill=list(fill="lightblue", alpha=0.5),
            points = list(pch=NA)) + theme_bw() + ylim(0,0.5) +
            ylab("Occurrence probability")
p6 = visreg(best1, "ORCHARD", scale="response", gg=TRUE, line=list(col="black"),
            fill=list(fill="lightblue", alpha=0.5),
            points = list(pch=NA)) + theme_bw() + ylim(0,0.5) +
            ylab("Occurrence probability")
p7 = visreg(best1, "SHRUBS", scale="response", gg=TRUE, line=list(col="black"),
            fill=list(fill="lightblue", alpha=0.5),
            points = list(pch=NA)) + theme_bw() + ylim(0,0.5) +
            ylab("Occurrence probability")
p8 = visreg(best1, "AES", scale="response", by="SHRUBS", gg=TRUE, line=list(col="black"),
            fill=list(fill="lightblue", alpha=0.5),
            points = list(pch=NA)) + theme_bw() + ylim(0,0.5) +
            ylab("Occurrence probability")
allp = grid.arrange(p1, p2, p3, p4, p5, p6, p7, p8, nrow=3, ncol=3)  

# plot also intrinsic interactive effect of LANDDIV on AES
p8 = visreg(best1, "AES", scale="response", by="LANDDIV", gg=TRUE, line=list(col="black"),
            fill=list(fill="lightblue", alpha=0.5),
            points = list(pch=NA)) + theme_bw() + ylim(0,0.5) +
  ylab("Occurrence probability")


## EFFECT MAPS -------------------------------------------------------------
# make an effect map to identify where it would be most effective to convert grassland to AES
# load the whole Catalan grid
gridES = st_read("C:/Users/sroilo/Desktop/BESTMAP documents/Data/Across_CS/ES/Grid_vars_env_20221216.shp")
# filter to only retain the grid cells for which land-cover data is complete 
gridES = gridES[gridES$TOTAL_C>0.99,] 
# filter out grid cells with less than 1 ha of grassland and outside the elevation limits of the shrike
grbir = gridES[which(gridES$Grass>=0.01 & gridES$Elevation>200 & gridES$Elevation<2000),]
names(grbir)[3:13] = c("TMAX", "PRECIP","ELEVATION", "AES", "ARABLE","SHRUBS", "GRASS","FOREST","ORCHARD", "URBAN","LANDDIV")

# Model prediction on the data of AES adoption of 2019
#first, correct for eventual discrepancies in AES > GRASS data
grbir$AES = ifelse(grbir$AES > grbir$GRASS, grbir$GRASS, grbir$AES)
# project the model on the current scenario (as of 2019) of AES adoption
grbir$pred_true = predict.glm(best1, newdata=grbir, type="response")
# set AES equal to grassland -> project model in a scenario in which all grassland is converted to AES
grbir$oldAES = grbir$AES # save old AES information in a different column
grbir$AES = grbir$GRASS
grbir$predAES = predict.glm(best1, newdata=grbir, type="response")
# calculate the arithmetic difference between the two projections and plot the effect map 
grbir$diff = grbir$predAES - grbir$pred_true
summary(grbir$diff)  # check which thresholds to use as breaks for the different colour classes
mapview(grbir, zcol="diff", col.regions = hcl.colors(8, palette="Red-Green"), at = c(-0.7, -0.3, -0.1, 0, 0.1, 0.3, 1))
mapview(grbir, zcol="SHRUBS")

# make some proper maps with tmap
library(tmap)
# load Catalan boundary
catalonia = st_read("C:/Users/sroilo/Desktop/BESTMAP documents/Data/ES_Catalonia/Catalonia boundary/Catalonia_EPSG3035.shp")
# save without titles which override the maps
tmap1 = tm_shape(catalonia) + tm_fill() +
  tm_shape(grbir) + tm_polygons(col = "pred_true", palette="viridis", border.alpha = 0.5, title = expression("")) + 
  tm_scale_bar(breaks=c(0,0.2,0.4,0.6,0.8,1))
tmap2 = tm_shape(catalonia) + tm_fill() +
  tm_shape(grbir) + tm_polygons(col = "predAES", palette="viridis", border.alpha = 0.5, title = expression("")) + 
  tm_scale_bar(breaks=c(0,0.2,0.4,0.6,0.8,1))
tmap3 = tm_shape(catalonia) + tm_fill() +
  tm_shape(grbir) + tm_polygons(col = "diff", palette="RdYlGn", midpoint=0,
                                breaks = c(-0.7, -0.4, 0, 0.4, 0.7, 1),
                                border.alpha = 0.5, title = expression(""))
all = tmap_arrange(tmap1, tmap2, tmap3)
tmap_save(all, width= 4, height=12, filename="C:/Users/sroilo/Desktop/BESTMAP documents/Data/Across_CS/ES/effect_maps_noTitle_20230715_vert2.jpeg")
# save in square format
tmap_save(all, width= 8, height=8, filename="C:/Users/sroilo/Desktop/BESTMAP documents/Data/Across_CS/ES/effect_maps_noTitle_20230715_sq2.jpeg")
# save only change map
tmap_save(all, width= 5, height=5, filename="C:/Users/sroilo/Desktop/BESTMAP documents/Data/Across_CS/ES/Change_map_noTitle_202307152.jpeg")

#save the shapefile to file 
st_write(grbir, "C:/Users/sroilo/Desktop/BESTMAP documents/Data/Across_CS/ES/July_trials/LC_predictions_AES_20230715.shp")


rm(list=ls())
setwd("C:/Users/sroilo/Documents")
