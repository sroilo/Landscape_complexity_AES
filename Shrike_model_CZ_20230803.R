###################################################
#
# Title: Shrike_model_CZ_20230803.R
# Purpose: run a Species Distribution Model for the red-backed shrike in South Moravia, Czech Republic, to quantify the effect of grassland-based Agri-Environmental Schemes (AES)
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
library(units)
library(terra)
library(psych)
library(MatchIt)
library(MuMIn)
library(visreg)
library(DHARMa)
library(exactextractr)
library(data.table)
library(ncf)
library(interflex)
library(gridExtra)

### DATA PREPARATION -----------------
# rasterize the IACS data of 2019 over the S2GLC map, to make sure grassland and arable land cover are accurate
# prepare a 1km grid for CZ
bound = st_read("C:/Users/sroilo/Desktop/BESTMAP documents/Data/CZ_SouthMoravia/South Moravia boundary/SouthMoravia_EPSG3035.shp")
grid = st_make_grid(bound, cellsize=1000, square=T, what="polygons") %>%  st_as_sf()
grid = grid[which(st_intersects(grid, bound, sparse=F) == T),]
grid$FID = c(1:nrow(grid))
# load small woody feature map to be burnt in the land cover map
swf = rast("X:/Stephanie/Data/CZ_SouthMoravia/Land cover/SWF_CZ.tif")
# load LC map, disaggregate it to 5 m resolution and crop to grid extent
lcmap = rast("U:/Stephanie/Land cover maps/Europe/S2GLC/S2GLC_Europe_2017_v1.2.tif") %>% resample(swf, method="near")
# load IACS data
in9 = st_read("C:/Users/sroilo/Desktop/BESTMAP documents/Data/CZ_SouthMoravia/IACS data/edited/CZ_fields2019.shp")
# only rasterize the fields that are cropland or grassland fields
# KULTURA field: 2 = arable, 7 = permanent grassland, 11 = grassland on arable land, 10 = fallow land
in9 = in9[which(in9$KULTURA %in% c(2,7,11,10)),]
# rename the column according to the S2GLC classes before rasterizing
in9$class = ifelse(in9$KULTURA %in% c(2,11,10), 73, 102)
in9r = rasterize(vect(in9), y=lcmap, field="class", update=T)
# burn in the SWF
inrm = mask(x=in9r, mask=swf, maskvalues=1, updatevalue=5, filename="C:/Users/sroilo/Desktop/BESTMAP documents/Data/Across_CS/CZ/CZraster/LCmap_iacs2019_SWF.tif", overwrite=T)
# extract the % cover of each land use type from the raster, for each polygon in the grid
# code taken from here: https://tmieno2.github.io/R-as-GIS-for-Economists/demo3.html 
cdl_extracted <- exact_extract(inrm, grid) %>%
  lapply(., function(x) data.table(x)[, .N, by = value]) %>%
  #--- combine the list of data.tables into one data.table ---#
  rbindlist(idcol = TRUE) %>%
  #--- find the share of each land use type ---#
  .[, share := N / sum(N), by = .id] %>%
  .[, N := NULL] 
#from long to wide format
cdl_extracted <- cdl_extracted %>% dcast(.id ~ value, value.var = "share")
head(cdl_extracted)
# substitute the NAs with zeros
cdl_extracted[is.na(cdl_extracted)] <-0
# now add the info on the different land cover classes to the grid
grid$SWF = cdl_extracted$`5`
grid$Arable = cdl_extracted$`73`
grid$PermCult = cdl_extracted$`75`
grid$Grass = cdl_extracted$`102`
grid$Wooded = rowSums(cdl_extracted[,c("82", "83","103", "105")])  # sum of broadleaf and conifer forest, moors and heathland
#grid$Unstb = rowSums(cdl_extracted[,c("62",  "82", "83", "162")]) # sum of artificial, forests and water bodies
# compute land-cover diversity
for ( i in seq_along(cdl_extracted$.id)) {
  prop = as.vector(unlist(cdl_extracted[i,2:11]))
  prop = prop[which(prop>0)]
  grid$LCdiv[i] = ifelse(length(prop)>0, -sum(sapply(X = prop, FUN = function(x) {x * log(x)}, simplify = "vector")), 0)
}
mapview(grid, zcol="LCdiv")
rm(swf, i , prop, cdl_extracted, lcmap, in9, in9r)
# plot and check the data
plot(grid[,3:8], borders=0)
# now, extract temperature, precipitation, elevation
r = rast(crs="epsg:3035", ext(vect(grid)), resolution=c(1000,1000))
tasmax = rast(c("X:/Stephanie/Data/Europe/Climate/CHELSA/Climatologies/CHELSA_tasmax_05_1981-2010_V.2.1.tif",
                "X:/Stephanie/Data/Europe/Climate/CHELSA/Climatologies/CHELSA_tasmax_06_1981-2010_V.2.1.tif",
                "X:/Stephanie/Data/Europe/Climate/CHELSA/Climatologies/CHELSA_tasmax_07_1981-2010_V.2.1.tif")) %>% project(r) %>% app(max)
precip = rast(c("X:/Stephanie/Data/Europe/Climate/CHELSA/Climatologies/CHELSA_pr_05_1981-2010_V.2.1.tif",
                "X:/Stephanie/Data/Europe/Climate/CHELSA/Climatologies/CHELSA_pr_06_1981-2010_V.2.1.tif",
                "X:/Stephanie/Data/Europe/Climate/CHELSA/Climatologies/CHELSA_pr_07_1981-2010_V.2.1.tif")) %>% project(r) %>% app(sum)
grid$Tmax = terra::extract(tasmax, vect(grid), fun="mean")[,2]
grid$Precip = terra::extract(precip, vect(grid), fun="mean")[,2]
mapview(grid, zcol="Tmax")
dem = rast("X:/Stephanie/Data/DE_Mulde/Topography/eu_dem_v11_E40N30/eu_dem_v11/eu_dem_v11_E40N20.TIF")
grid$Elevation = terra::extract(dem, vect(grid), fun="mean", na.rm=T)[,2]
mapview(grid, zcol="Elevation")

# prepare the bird dataset
# read bird dataset (observation points of bird species)
birds = st_read("C:/Users/sroilo/Desktop/BESTMAP documents/Data/CZ_SouthMoravia/Biodiversity/Birds_Butterflies_12112020/edited/CZ_2010-2020_Farmland_birds_NDOP.shp") 
# translate names into English and get rid of some more columns
birds = birds[,c("DRUH","DATUM_OD","geometry")]
birds$DATUM_OD = as.Date(birds$DATUM_OD, format = "%d.%m.%Y")
### January 17th, 2023, add the extra data from Tomas, CZ-specific list of species was taken from here: https://www.tandfonline.com/doi/full/10.1080/00063657.2015.1048423 
extrab = st_read("C:/Users/sroilo/Desktop/BESTMAP documents/Data/CZ_SouthMoravia/Biodiversity/Birds_addition_11012023/Birds_addition_11012023.shp") %>% st_transform(3035)
extrab = extrab[,c("DRUH","DATUM_OD","geometry")]
# remove species with imprecise taxonomy
extrab = extrab[-which(extrab$DRUH=="Corvus corone/cornix"),]
extrab$DATUM_OD = as.Date(as.character(extrab$DATUM_OD), format = "%Y%m%d")
# merge the two datasets together
birds = bind_rows(birds, extrab)
# append the coordinates to the df
birds = cbind(st_coordinates(birds), birds)
names(birds) = c("X", "Y", "species", "date", "geometry" )
# get rid of entries without geometry
birds = na.omit(birds)
# order by decreasing date 
birds = birds[order(birds$date, decreasing=T),]
# subselect only 2016-2019 data
birds$year = strftime(birds$date, format = "%Y") # get year
birds = birds[which(birds$year %in% c(2016,2017,2018,2019)),]
# get summary of included species
data.frame(table(birds$species))
# prepare dataset for modelling of Lanius collurio
birds$presence = ifelse(birds$species=="Lanius collurio",1,0)

# rasterize presence data of Lanius collurio per 1km square, divided by year:
lc6 = rasterize(x=vect(birds[birds$year==2016,]), y=r, field="presence", fun="max")
lc7 = rasterize(x=vect(birds[birds$year==2017,]), y=r, field="presence", fun="max")
lc8 = rasterize(x=vect(birds[birds$year==2018,]), y=r, field="presence", fun="max")
lc9 = rasterize(x=vect(birds[birds$year==2019,]), y=r, field="presence", fun="max")
# extract values for each year, check if any place was a presence on one year and absence on another
grid$LC2016 = terra::extract(lc6, vect(grid), method="simple")[,2]
# check if everything went fine
plot(grid[,"LC2016"]); plot(lc6)   # YES!
# do the same for the other years
grid$LC2017 = terra::extract(lc7, vect(grid), method="simple")[,2]
grid$LC2018 = terra::extract(lc8, vect(grid), method="simple")[,2]
grid$LC2019 = terra::extract(lc9, vect(grid), method="simple")[,2]
# now, compute the cross-year presences and absences:
# find always most recent observation
grid$LCmr = ifelse(is.na(grid$LC2019), grid$LC2018, grid$LC2019)
grid$LCmr = ifelse(is.na(grid$LCmr), grid$LC2017, grid$LCmr)
grid$LCmr = ifelse(is.na(grid$LCmr), grid$LC2016, grid$LCmr)
# plot everything 
plot(grid[,c(12:16)], border="lightgrey")
# now, calculate the proportion of GrM in each 1km square in the 4 years
for (y in c(2016, 2017, 2018, 2019)) {
  # load IACS layers
  indat = st_read(paste0("C:/Users/sroilo/Desktop/BESTMAP documents/Data/CZ_SouthMoravia/IACS data/edited/CZ_fields", y, ".shp"))
  # retain only fields with some GrM
  indat = indat[which(indat$GrM>0),]
  #calculate the ratio of GrM area on the total field area (adjust for errors in the data where GrM>LP_area)
  indat$GrMratio = ifelse(indat$GrM/indat$LP_area_ha > 1, 1, indat$GrM/indat$LP_area_ha)
  # compute intersection between grid and IACS data
  int = st_intersection(x = grid,  y = st_make_valid(indat))
  int$area = set_units(st_area(int), NULL)
  # start looping through the observations and compute metrics from the IACS data
  for (i in seq_along(grid$FID)) {
    subset = st_drop_geometry(int[int$FID == grid$FID[i],])
    grid[i, paste0("GrM_", y)] = ifelse((nrow(subset) == 0), 0, sum(subset$GrMratio*subset$area)/1000000)
  }
}
# find the GrM proportion corresponding to the most recent bird observation
grid$GrMmr = NA
grid$GrMmr = ifelse(is.na(grid$LC2016), grid$GrMmr, grid$GrM_2016)
grid$GrMmr = ifelse(is.na(grid$LC2017), grid$GrMmr, grid$GrM_2017)
grid$GrMmr = ifelse(is.na(grid$LC2018), grid$GrMmr, grid$GrM_2018)
grid$GrMmr = ifelse(is.na(grid$LC2019), grid$GrMmr, grid$GrM_2019)

# add the coordinates of the centroids of each grid cell to be used in the modelling to control for spatial autocorrelation
cents = st_centroid(grid)
grid <- cbind(st_coordinates(cents), grid)
# check how many grid cells were monitored each year, and for several years:
sum(!is.na(grid$LC2016))  # 516
sum(!is.na(grid$LC2017))  # 476 
sum(!is.na(grid$LC2018))  # 437 
sum(!is.na(grid$LC2019))  # 376 
sum(!is.na(grid$LCmr))    # 1008
multip = rowSums(is.na(grid[,13:16]))
table(multip)# 1253 grid cells NEVER monitored, 527 monitored only once, 261 monitored twice, 124 monitored 3 times, 96 monitored 4 times
plot(grid[,13:17], border="lightgrey")

# save everything to file
st_write(grid, "C:/Users/sroilo/Desktop/BESTMAP documents/Data/Across_CS/CZ/Grid_allVars_20230120.shp")
rm ( indat, int, inrm, lc6, lc7, lc8, lc9, precip, tasmax, subset, r, multip, y, i, cents, birds, extrab, dem)


### MODELLING ---------------------
grid = st_read("C:/Users/sroilo/Desktop/BESTMAP documents/Data/Across_CS/CZ/Grid_allVars_20230120.shp")
# use as presences (and absences) of Lanius collurio the most recent sightings
grid$presence = grid$LCmr
# filter to only retain the monitored grid cells
grid = grid[-which(is.na(grid$presence)),]  # 126 presences
grid = cbind(st_drop_geometry(grid[,1:3]), grid$presence, grid[, c(4:12, 22)])
# rename some columns, set GrM name as AES
names(grid)[4:14] <- c("presence", "SHRUBS","ARABLE","ORCHARD","GRASS","FOREST", "LANDDIV", "TMAX","PRECIP","ELEVATION","AES")
# filter to only retain grid cells with min. 1 ha of grassland
lces = grid[which(grid$GRASS>=0.01),]  #877 grid cells left
# AES was calculated from shapefiles, while GRASS from raster data. Make sure that there are no major discrepancies
nrow(lces[which(lces$AES > lces$GRASS),]) # 4 grid cells in which AES > GRASS; correct 
lces$AES = ifelse(lces$AES > lces$GRASS, lces$GRASS, lces$AES)

# check whether there is enough common support to investigate interactions between AES and SWF and LCdiv
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
quantile(lces$AES, probs = c(0.9, 0.95, 0.99))           # 90% -> 0.2552287 | 95% -> 0.4261586 | 99% -> 0.6324632
lces = lces[which(lces$AES <= 0.2552287),]  # grid-cells left:     789      |            833   |       868  
# repeat the common support check above via the interflex function and move to the next cutoff value for the trimming

# plot scatterplots
lceso = lces[order(lces$presence, decreasing=F),]
pairs.panels(lceso[,c("presence", "AES", "SHRUBS", "LANDDIV", "GRASS", "ARABLE", "FOREST", "ORCHARD", "TMAX", "PRECIP", "ELEVATION")],
             bg=c("blue","yellow")[as.factor(lceso$presence)], pch=21,
             method = "spearman", hist.col = "darkorange", density = TRUE, ellipses = FALSE )

# normalise variables in advance of the modelling:
dat = cbind(lces[,1:4], scale(lces[,5:14]) )
# create correlation matrix to exclude highly correlated variables from the same model
cmat <- abs(cor(dat[,5:14], method="spearman")) <= 0.7
cmat[!lower.tri(cmat, diag=F)] <- NA  

### MODEL SELECTION ----------------
fit = glm(presence ~ AES + ARABLE + GRASS + SHRUBS + LANDDIV +  FOREST + ORCHARD + AES:LANDDIV + AES:SHRUBS + TMAX + PRECIP + ELEVATION, 
          data = dat, family="binomial", na.action="na.fail") 
dr1 = dredge(fit, subset=cmat, extra = c("R^2", "adjR^2")) 
write.table(dr1, sep=";", dec=",", row.names=F,  "C:/Users/sroilo/Desktop/BESTMAP documents/Data/Across_CS/CZ/July_trials/LC_trim0426_20230715.csv")
# get variable importance of highly supported models (deltaAIC<2)
drsub = dr1[which(dr1$delta<2),]
imp = sw(drsub)
data.frame(imp)
write.table(data.frame(imp),  sep=";", dec=",", row.names=T, "C:/Users/sroilo/Desktop/BESTMAP documents/Data/Across_CS/CZ/July_trials/LC_trim0426_20230715_weights.csv")

# fit best model overall:
best1 = glm(presence ~ AES + ARABLE + FOREST + GRASS + ORCHARD + SHRUBS + TMAX,
            data = dat, family="binomial", na.action="na.fail")
summary(best1)
# diagnostic plots with DHARMa package
simulateResiduals(best1, plot=T)  # all good here
# check for spatial autocorrelation in the residuals
max.dist = max(dist(dat[,c("X","Y")])) * 2/3
scor = ncf::spline.correlog(x=best1$data$X, y = best1$data$Y, z=best1$residuals, resamp=100, xmax = max.dist) 
plot(scor)  # no spatial autocorrelation

# refit the model with the original (non-scaled) data to produce conditional plots
best1 = glm(presence ~ AES + ARABLE + FOREST + GRASS + ORCHARD + SHRUBS + TMAX,
            data = lces, family="binomial", na.action="na.fail")
# plot conditional plots for each predictor and interaction term
p1 = visreg(best1, "AES", scale="response", gg=TRUE, line=list(col="black"),
            fill=list(fill="lightblue", alpha=0.5),
            points = list(pch=NA)) + theme_bw() + ylim(0,1) +
  ylab("Occurrence probability")
p2 = visreg(best1, "ARABLE", scale="response", gg=TRUE, line=list(col="black"),
            fill=list(fill="lightblue", alpha=0.5),
            points = list(pch=NA)) + theme_bw() + ylim(0,1) +
  ylab("Occurrence probability")
p3 = visreg(best1, "FOREST", scale="response", gg=TRUE, line=list(col="black"),
            fill=list(fill="lightblue", alpha=0.5),
            points = list(pch=NA)) + theme_bw() + ylim(0,1) +
  ylab("Occurrence probability")
p4 = visreg(best1, "GRASS", scale="response", gg=TRUE, line=list(col="black"),
            fill=list(fill="lightblue", alpha=0.5),
            points = list(pch=NA)) + theme_bw() + ylim(0,1) +
  ylab("Occurrence probability")
p5 = visreg(best1, "ORCHARD", scale="response", gg=TRUE, line=list(col="black"),
            fill=list(fill="lightblue", alpha=0.5),
            points = list(pch=NA)) + theme_bw() + ylim(0,1) +
  ylab("Occurrence probability")
p6 = visreg(best1, "SHRUBS", scale="response", gg=TRUE, line=list(col="black"),
            fill=list(fill="lightblue", alpha=0.5),
            points = list(pch=NA)) + theme_bw() + ylim(0,1) +
  ylab("Occurrence probability")
p7 = visreg(best1, "TMAX", scale="response", gg=TRUE, line=list(col="black"),
            fill=list(fill="lightblue", alpha=0.5),
            points = list(pch=NA)) + theme_bw() + ylim(0,1) +
  ylab("Occurrence probability")
p8 = visreg(best1, "AES", scale="response", by="SHRUBS", gg=TRUE, line=list(col="black"),
            fill=list(fill="lightblue", alpha=0.5),
            points = list(pch=NA)) + theme_bw() + ylim(0,1) +
  ylab("Occurrence probability")
allp = grid.arrange(p1, p2, p3, p4, p5, p6, p7, p8, nrow=3, ncol=3)  

### EFFECT MAP  --------------------------------------------
# make an effect map to identify where it would be most effective to convert grassland to AES
grid = st_read("C:/Users/sroilo/Desktop/BESTMAP documents/Data/Across_CS/CZ/Grid_allVars_20230120.shp") 
grid= grid[grid$Grass>=0.01,]
# set the AES variable as GrM_2019 (correct for eventual discrepancies with AES > GRASS)
grid$AES = ifelse(grid$GrM_2019 > grid$Grass, grid$Grass, grid$GrM_2019)
# rename variables
names(grid)[4:12] = c("SHRUBS", "ARABLE", "ORCHARD", "GRASS", "FOREST", "LANDDIV", "TMAX", "PRECIP", "ELEVATION")
# prediction on actual data
grid$pred_true = predict.glm(best1, newdata=grid, type="response")
# set AES equal to grassland -> if all grassland were extensive
grid$oldAES = grid$AES
grid$AES = grid$GRASS
grid$predAES = predict.glm(best1, newdata=grid, type="response")
# calculate difference and plot effect map 
grid$diff = grid$predAES - grid$pred_true
summary(grid$diff)

# make some proper maps
library(tmap)
# load South Moravian boundary
southmor = st_read("C:/Users/sroilo/Desktop/BESTMAP documents/Data/CZ_SouthMoravia/South Moravia boundary/SouthMoravia_EPSG3035.shp")
tmap1 = tm_shape(southmor) + tm_fill() +  
  tm_shape(grid) + tm_polygons(col = "pred_true", palette="viridis", border.alpha = 0.5, title = expression("")) + 
  tm_scale_bar(breaks=c(0,0.2,0.4,0.6,0.8,1))
tmap2 = tm_shape(southmor) + tm_fill() +  
  tm_shape(grid) + tm_polygons(col = "predAES", palette="viridis", border.alpha = 0.5, title = expression("")) + 
  tm_scale_bar(breaks=c(0,0.2,0.4,0.6,0.8,1))
tmap3 = tm_shape(southmor) + tm_fill() +  
  tm_shape(grid) + tm_polygons(col = "diff", palette="RdYlGn", midpoint=0, #breaks = c(-0.3, 0, 0.1, 0.3, 0.7),
                               border.alpha = 0.5, title = expression(""), legend.position="none")
all = tmap_arrange(tmap1, tmap2, tmap3)
tmap_save(all, width= 4, height=12, filename="C:/Users/sroilo/Desktop/BESTMAP documents/Data/Across_CS/CZ/effect_maps_bestM_20230715_vert.jpeg")

#save the shapefile to file 
st_write(grid, "C:/Users/sroilo/Desktop/BESTMAP documents/Data/Across_CS/CZ/LC_predictions_AES_20230715_bestM.shp")

rm(list=ls())
setwd("C:/Users/sroilo/Documents")
