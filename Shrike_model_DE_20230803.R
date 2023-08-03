###################################################
#
# Title: Shrike_model_DE_20230803.R
# Purpose: run a Species Distribution Model for the red-backed shrike in the Mulde River Basin, Germany, to quantify the effect of grassland-based Agri-Environmental Schemes (AES)
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
library(MuMIn)
library(visreg)
library(DHARMa)
library(exactextractr)
library(data.table)
library(ncf)
library(mgcv)
library(interflex)
library(gridExtra)

### DATA PREPARATION ------------------------------------------------
# load Mulde 1km grid
grid = st_read("C:/Users/sroilo/Desktop/BESTMAP documents/Data/DE_Mulde/Mulde/Mulde_1km_grid.shp") %>% st_transform(3035)
# load the SWF map to be overlaid on top of the land cover map
swf = rast("C:/Users/sroilo/Desktop/BESTMAP documents/Data/DE_Mulde/Land cover/Small Woody Features 2015_DE030/swf_2015_Mulde_large.tif")
# load LC map
lcmap = rast("U:/Stephanie/Land cover maps/Europe/S2GLC/S2GLC_Europe_2017_v1.2.tif") %>% resample(swf, method="near")
# load invekos data
in9 = st_read("C:/Users/sroilo/Desktop/BESTMAP documents/Data/DE_Mulde/INVEKOS_data/edited/InVeKoS_2019_Sc.shp")
# load file with crop-codes and corresponding crop names
codes <- read.table("C:/Users/sroilo/Desktop/BESTMAP documents/Data/DE_Mulde/INVEKOS_data/NC_LIST_SN_2020_sc.txt", header = T, sep="\t", stringsAsFactors = FALSE)
codes[1,1] <- "050"
# match crop groups
in9$BNK <- codes$BNK[match(in9$HKCODE, codes$NC)]
# only rasterize the fields that are cropland or grassland fields
in9 = in9[which(in9$BNK %in% c("AL", "DGL")),]
# rename the column according to the S2GLC classes before rasterizing
in9$class = ifelse(in9$BNK=="AL", 73, 102)
in9r = rasterize(vect(in9), y=lcmap, field="class", update=T)
# burn in the SWF
inrm = mask(x=lcmap, mask=swf, maskvalues=1, updatevalue=5, filename="C:/Users/sroilo/Desktop/BESTMAP documents/Data/Across_CS/DE/DEraster/LCmap_iacs2019_SWF.tif", overwrite=T)
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
grid$Grass = cdl_extracted$`102`
#grid$Forest = rowSums(cdl_extracted[,c("82", "83")], na.rm=T) 
grid$Wooded = rowSums(cdl_extracted[,c("82", "83","103", "105", "106")], na.rm=T) # this includes forests, moors, heathland, marshes and peatbogs
grid$Urban = cdl_extracted$`62` 
grid$PermCult = cdl_extracted$`75` # vineyards
# compute land cover diversity
for ( i in seq_along(cdl_extracted$.id)) {
  prop = as.vector(unlist(cdl_extracted[i,2:14]))
  prop = prop[which(prop>0)]
  grid$LCdiv[i] = ifelse(length(prop)>0, -sum(sapply(X = prop, FUN = function(x) {x * log(x)}, simplify = "vector")), 0)
}
mapview(grid, zcol="LCdiv")

# now, extract temperature, precipitation, elevation
r = rast(crs="epsg:3035",extent=ext(grid), resolution=c(1000,1000))
tasmax = rast(c("X:/Stephanie/Data/Europe/Climate/CHELSA/Climatologies/CHELSA_tasmax_05_1981-2010_V.2.1.tif",
                "X:/Stephanie/Data/Europe/Climate/CHELSA/Climatologies/CHELSA_tasmax_06_1981-2010_V.2.1.tif",
                "X:/Stephanie/Data/Europe/Climate/CHELSA/Climatologies/CHELSA_tasmax_07_1981-2010_V.2.1.tif")) %>% project(r) %>% app(max)
precip = rast(c("X:/Stephanie/Data/Europe/Climate/CHELSA/Climatologies/CHELSA_pr_05_1981-2010_V.2.1.tif",
                "X:/Stephanie/Data/Europe/Climate/CHELSA/Climatologies/CHELSA_pr_06_1981-2010_V.2.1.tif",
                "X:/Stephanie/Data/Europe/Climate/CHELSA/Climatologies/CHELSA_pr_07_1981-2010_V.2.1.tif")) %>% project(r) %>% app(sum)
grid$Tmax = terra::extract(tasmax, vect(grid), fun="mean")[,2]
grid$Precip = terra::extract(precip, vect(grid), fun="mean")[,2]
mapview(grid, zcol="Tmax")
dem = rast("X:/Stephanie/Data/DE_Mulde/Topography/eu_dem_v11_E40N30/EU-DEM_v11_Mulde.tif")
grid$Elevation = terra::extract(dem, vect(grid), fun="mean")[,2]
rm(precip, tasmax, dem, i , prop, cdl_extracted)
# plot and check the data
plot(grid[,3:11], borders=0)

# prepare the bird dataset
# read bird dataset (observation points of bird species)
birds = st_read("C:/Users/sroilo/Desktop/BESTMAP documents/Data/DE_Mulde/Biodiversity/ZenA Anfrage BESTMAP/by groups/Birds_2000-2020.shp") 
# translate names into English and get rid of some more columns
birds = birds[,c(2,3,10,23,25,26,27,28,29,31,32,34,35,37)]
names(birds) = c("X", "Y", "species", "uncertainty", "date", "month", "year", "source", "behaviour", "ind.count", "unit", "reproduction",    
                 "quality", "geometry" )
# order by decreasing date 
birds$date = as.Date(birds$date, format = "%d.%m.%Y")
birds = birds[order(birds$date, decreasing=T),]
# subselect only 2016-2019 data
birds = birds[which(birds$year %in% c(2016, 2017, 2018, 2019)),]
# select only species from Busch et al. or Gamero et al. 
blist = read.csv("C:/Users/sroilo/Desktop/BESTMAP documents/Data/Across_CS/DE/Bird_selection_list.csv", sep=";")
birds = birds[which(birds$species %in% blist$Species),]
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
#check if everything went fine
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
plot(grid[,c(3:11)], border=0)
plot(grid[,c(12:16)], border=0)
# now, calculate the proportion of GrM in each 1km square in the 4 years
# load invekos layers
in6 = st_read("C:/Users/sroilo/Desktop/BESTMAP documents/Data/DE_Mulde/INVEKOS_data/edited/InVeKoS_2016_Sc.shp")
in7 = st_read("C:/Users/sroilo/Desktop/BESTMAP documents/Data/DE_Mulde/INVEKOS_data/edited/InVeKoS_2017_Sc.shp")
in8 = st_read("C:/Users/sroilo/Desktop/BESTMAP documents/Data/DE_Mulde/INVEKOS_data/edited/InVeKoS_2018_Sc_NEU.shp")
in9 = st_read("C:/Users/sroilo/Desktop/BESTMAP documents/Data/DE_Mulde/INVEKOS_data/edited/InVeKoS_2019_Sc.shp")
# filter invekos layers to only retain the GrM fields
in9 = in9[which(grepl(pattern="GL1|GL2|GL4|GL5", x = in9$AUM_INFO1)| grepl(pattern="GL1|GL2|GL4|GL5", x = in9$AUM_INFO2)),]
in8 = in8[which(grepl(pattern="GL1|GL2|GL4|GL5", x = in8$AUM_INFO1)| grepl(pattern="GL1|GL2|GL4|GL5", x = in8$AUM_INFO2)),]
in7 = in7[which(grepl(pattern="GL1|GL2|GL4|GL5", x = in7$AUM_INFO1)| grepl(pattern="GL1|GL2|GL4|GL5", x = in7$AUM_INFO2)),]
in6 = in6[which(grepl(pattern="GL1|GL2|GL4|GL5", x = in6$AUM_INFO1)| grepl(pattern="GL1|GL2|GL4|GL5", x = in6$AUM_INFO2)),]
# year 2019
int = st_intersection(x = grid,  y = st_make_valid(in9))
int$area = set_units(st_area(int), NULL)
# start looping through the observations and compute metrics from the invekos data
for (i in seq_along(grid$FID)) {
    subset = st_drop_geometry(int[int$FID == grid$FID[i],])
    grid$GrM_2019[i] = ifelse((nrow(subset) == 0), 0, sum(subset$area)/1000000)
}
# year 2018
int = st_intersection(x = grid,  y = st_make_valid(in8))
int$area = set_units(st_area(int), NULL)
# start looping through the observations and compute metrics from the invekos data
for (i in seq_along(grid$FID)) {
  subset = st_drop_geometry(int[int$FID == grid$FID[i],])
  grid$GrM_2018[i] = ifelse((nrow(subset) == 0), 0, sum(subset$area)/1000000)
}
# year 2017
int = st_intersection(x = grid,  y = st_make_valid(in7))
int$area = set_units(st_area(int), NULL)
# start looping through the observations and compute metrics from the invekos data
for (i in seq_along(grid$FID)) {
  subset = st_drop_geometry(int[int$FID == grid$FID[i],])
  grid$GrM_2017[i] = ifelse((nrow(subset) == 0), 0, sum(subset$area)/1000000)
}
# year 2016
int = st_intersection(x = grid,  y = st_make_valid(in6))
int$area = set_units(st_area(int), NULL)
# start looping through the observations and compute metrics from the invekos data
for (i in seq_along(grid$FID)) {
  subset = st_drop_geometry(int[int$FID == grid$FID[i],])
  grid$GrM_2016[i] = ifelse((nrow(subset) == 0), 0, sum(subset$area)/1000000)
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
sum(!is.na(grid$LC2016))  #327
sum(!is.na(grid$LC2017))  #696
sum(!is.na(grid$LC2018))  #477
sum(!is.na(grid$LC2019))  #519
sum(!is.na(grid$LCmr))   # 1301
multip = rowSums(is.na(grid[,14:17]))
table(multip)# 4854 grid cells NEVER monitored, 834 monitored only once, 278 monitored twice, 132 monitored 3 times, 60 monitored 4 times

#save to new file
st_write(grid, "C:/Users/sroilo/Desktop/BESTMAP documents/Data/Across_CS/DE/DEGrid_allVars_20230119.shp")

##### MODELLING -------------------------------------
grid = st_read("C:/Users/sroilo/Desktop/BESTMAP documents/Data/Across_CS/DE/DEGrid_allVars_20230119.shp")
# use as presences (and absences) of Lanius collurio the most recent sightings
grid$presence = grid$LCmr
# filter to only retain the monitored grid cells
grid = grid[-which(is.na(grid$presence)),]  #1304 rows, 425 presences
grid = grid[,c(1:7,10:13, 23:25)]
names(grid)  <- c("X", "Y", "FID","SHRUBS","ARABLE","GRASS","FOREST","LANDDIV","TMAX","PRECIP","ELEVATION", "AES", "geometry", "presence")
# filter to only retain grid cells with min. 1ha grassland
lces = st_drop_geometry(grid[which(grid$GRASS>=0.01),])  #1237 rows
# AES was calculated from shapefiles, while GRASS from raster data. Make sure that there are no major discrepancies
nrow(lces[which(lces$AES > lces$GRASS),]) # 18 grid cells in which AES > GRASS; correct 
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

### TRIM the dataset to reach common support  ------------------------------
# use three different thresholds, like a sensitivity analysis
quantile(lces$AES, probs = c(0.9, 0.95, 0.99))           # 90% -> 0.1734424 | 95% -> 0.2561157 | 99% -> 0.4692268 
lces = lces[which(lces$AES <= 0.1734424),]  # landscapes left:        1113      |           1175   |      1224       |  full dataset: 1237
# repeat the common support check above via the interflex function and move to the next cutoff value for the trimming

# plot scatterplots
lceso = lces[order(lces$presence, decreasing=F),]
pairs.panels(lceso[,c("presence", "AES", "SHRUBS", "LANDDIV", "GRASS", "ARABLE", "FOREST", "TMAX", "PRECIP", "ELEVATION")],
             bg=c("blue","yellow")[as.factor(lceso$presence)], pch=21,
             method = "spearman", hist.col = "darkorange", density = TRUE, ellipses = FALSE )

# normalise variables in advance of the modelling:
dat = lces
dat = cbind(dat[,1:3], dat[,13], dat[,c(4:12)])
names(dat)[4] <-"presence"
dat[,c(5:13)] = scale(dat[,c(5:13)])
# create correlation matrix to exclude highly correlated variables from the same model
cmat <- abs(cor(dat[,5:13], method="spearman")) <= 0.7
cmat[!lower.tri(cmat, diag=F)] <- NA  

### MODEL SELECTION -----------------
fit = glm(presence ~ AES + ARABLE + GRASS + SHRUBS + LANDDIV + FOREST + AES:LANDDIV + AES:SHRUBS + TMAX + PRECIP + ELEVATION,
          data = dat, family="binomial", na.action="na.fail") 
dr1 = dredge(fit, subset=cmat, extra = c("adjR^2")) 
write.table(dr1, sep=";", dec=",", row.names=F,  "C:/Users/sroilo/Desktop/BESTMAP documents/Data/Across_CS/DE/July_trials/LC_trim0173_20230715.csv")
# get variable importance of highly supported models (deltaAIC<2)
drsub = dr1[which(dr1$delta<2),]
imp = sw(drsub)
data.frame(imp)
#write.table(data.frame(imp), sep=";", dec=",", row.names=T,  "C:/Users/sroilo/Desktop/BESTMAP documents/Data/Across_CS/DE/July_trials/LC_trim0173_20230715_weights.csv")
# fit best model with interaction term:
best1 = glm(presence ~ AES + FOREST + LANDDIV + PRECIP + SHRUBS, 
            data = dat, family="binomial", na.action="na.fail")
# diagnostic plots with DHARMa package
simulateResiduals(best1, plot=T)  # not great plots here...
# check for spatial autocorrelation in the residuals
max.dist = max(dist(dat[,c("X","Y")])) * 2/3
scor = ncf::spline.correlog(x=best1$data$X, y = best1$data$Y, z=best1$residuals, resamp=100, xmax = max.dist) # high SAC
plot(scor)  # significant correlation in the residuals!
# repeat using GEneralised Additive Models (GAMs), adding a spatial spline
library(mgcv)
fit = gam(presence ~ AES + ARABLE + GRASS + SHRUBS + LANDDIV + FOREST + AES:LANDDIV + AES:SHRUBS + TMAX + PRECIP + ELEVATION + s(X, Y, bs="gp", k= 150, m=2), 
          data = dat, family="binomial", na.action="na.fail") 
dr1 = dredge(fit, subset=cmat, extra = c("adjR^2")) 
write.table(dr1, sep=";", dec=",", row.names=F,  "C:/Users/sroilo/Desktop/BESTMAP documents/Data/Across_CS/DE/July_trials/LC_trim0256_20230715_GAMk150.csv")
# get variable importance of highly supported models (deltaAIC<2)
drsub = dr1[which(dr1$delta<2),]
imp = sw(drsub)
data.frame(imp)
write.table(data.frame(imp),  sep=";", dec=",", row.names=T, "C:/Users/sroilo/Desktop/BESTMAP documents/Data/Across_CS/DE/July_trials/LC_trim0256_20230715_GAMk150_weights.csv")
# fit best model with interaction term:
best1 = gam(presence ~ AES + PRECIP + SHRUBS + s(X, Y, bs="gp", k= 150, m=2), 
            data = dat, family="binomial", na.action="na.fail")
# diagnostic plots with DHARMa package
summary(best1) 
simulateResiduals(best1, plot=T)  # all good
gam.check(best1)
# check for spatial autocorrelation in the residuals
max.dist = max(dist(dat[,c("X","Y")])) * 2/3
scor = ncf::spline.correlog(x=dat$X, y = dat$Y, z=best1$residuals, resamp=100, xmax = max.dist) 
plot(scor)  
visreg(best1, scale="response")
plot(best1)  

# plot the conditional plots
best1 = gam(presence ~ AES + PRECIP + SHRUBS + s(X, Y, bs="gp", k= 150, m=2), 
            data = lces, family="binomial", na.action="na.fail")
p1 = visreg(best1, "AES", scale="response", gg=TRUE, line=list(col="black"),
            fill=list(fill="lightblue", alpha=0.5),
            points = list(pch=NA)) + theme_bw() + ylim(0,1) +
  ylab("Occurrence probability")
p2 = visreg(best1, "PRECIP", scale="response", gg=TRUE, line=list(col="black"),
            fill=list(fill="lightblue", alpha=0.5),
            points = list(pch=NA)) + theme_bw() + ylim(0,1) +
  ylab("Occurrence probability")
p3 = visreg(best1, "SHRUBS", scale="response", gg=TRUE, line=list(col="black"),
            fill=list(fill="lightblue", alpha=0.5),
            points = list(pch=NA)) + theme_bw() + ylim(0,1) +
  ylab("Occurrence probability")
p4 = visreg(best1, "X", scale="response", gg=TRUE, line=list(col="black"),
            fill=list(fill="lightblue", alpha=0.5),
            points = list(pch=NA)) + theme_bw() + ylim(0,1) + scale_x_continuous(breaks=c(4500000, 4535000, 4570000)) +
  ylab("Occurrence probability")
p5 = visreg(best1, "Y", scale="response", gg=TRUE, line=list(col="black"),
            fill=list(fill="lightblue", alpha=0.5),
            points = list(pch=NA)) + theme_bw() + ylim(0,1) + scale_x_continuous(breaks=c(3030000,3100000, 3170000)) +
  ylab("Occurrence probability")
p8 = visreg(best1, "AES", scale="response", by="SHRUBS", gg=TRUE, line=list(col="black"),
            fill=list(fill="lightblue", alpha=0.5),
            points = list(pch=NA)) + theme_bw() + ylim(0,1) +
  ylab("Occurrence probability")
allp = grid.arrange(p1, p2, p3, p4, p5, p8, nrow=3, ncol=3)  

### EFFECT MAPS ---------------
# make an effect map to identify where it would be most effective to convert grassland to AES
grid = st_read("C:/Users/sroilo/Desktop/BESTMAP documents/Data/Across_CS/DE/DEGrid_allVars_20230119.shp")
# set the AES variable as GrM_2019 (correct for eventual discrepancies with AES > GRASS)
grid$AES = ifelse(grid$GrM_2019 > grid$Grass, grid$Grass, grid$GrM_2019)
grid = grid[,c(1:7,10:13, 25)]
names(grid) <- c("X", "Y", "FID","SHRUBS","ARABLE","GRASS","FOREST","LANDDIV","TMAX","PRECIP","ELEVATION", "AES", "geometry")
# filter to only retain grid cells with min. 1ha grassland
grid = grid[which(grid$GRASS>=0.01),]
# prediction on actual data
grid$pred_true = predict.gam(best1, newdata=grid, type="response")
grid$oldAES = grid$AES
# set AES equal to grassland -> prediction in a scenario in which all grassland is converted to AES
grid$AES = grid$GRASS
grid$predAES = predict.gam(best1, newdata=grid, type="response")
# calculate arithmetic difference between the scenarios and plot effect map 
grid$diff = grid$predAES - grid$pred_true
summary(grid$diff)

# make some proper maps with tmap
library(tmap)
# load Mulde boundary
mulde = st_read("C:/Users/sroilo/Desktop/BESTMAP documents/Data/DE_Mulde/Mulde/Mulde_EPSG3035.shp")
tmap1 = tm_shape(mulde) + tm_fill() +  
  tm_shape(grid) + tm_polygons(col = "pred_true", palette="viridis", border.alpha = 0.5, title = expression("")) + 
  tm_scale_bar(breaks=c(0,0.2,0.4,0.6,0.8,1))
tmap2 = tm_shape(mulde) + tm_fill() +  
  tm_shape(grid) + tm_polygons(col = "predAES", palette="viridis", border.alpha = 0.5, title = expression("")) + 
  tm_scale_bar(breaks=c(0,0.2,0.4,0.6,0.8,1))
tmap3 = tm_shape(mulde) + tm_fill() +  
  tm_shape(grid) + tm_polygons(col = "diff", palette="RdYlGn", midpoint=0, 
                               border.alpha = 0.5, title = expression(""))
all = tmap_arrange(tmap1, tmap2, tmap3)
tmap_save(all, width= 4, height=12, filename="C:/Users/sroilo/Desktop/BESTMAP documents/Data/Across_CS/DE/effect_maps_bestM_20230715.jpeg")
#save the results to file 
st_write(grid, "C:/Users/sroilo/Desktop/BESTMAP documents/Data/Across_CS/DE/LC_predictions_AES_20230715.shp")

### CROSS-REGIONAL comparisons -------------------------
# 2023.07.07 - read and put together all predictors' relative importance measures from different regions
# load results from heaviest-trimmed datasets
wes = read.table("C:/Users/sroilo/Desktop/BESTMAP documents/Data/Across_CS/ES/July_trials/LC_trim0065_weights_20230714.csv", header=T, sep=";", dec=",")
wde = read.table("C:/Users/sroilo/Desktop/BESTMAP documents/Data/Across_CS/DE/July_trials/LC_trim0173_20230715_GAMk150_weights.csv", header=T, sep=";", dec=",")
wcz = read.table("C:/Users/sroilo/Desktop/BESTMAP documents/Data/Across_CS/CZ/July_trials/LC_trim0255_20230715_weights.csv", header=T, sep=";", dec=",")
wcz$Row.names = row.names(wcz)
allw = merge(wes, wde, by="row.names", all=T)
allw = merge(allw, wcz, by="Row.names", all=T)
names(allw) = c("Variable", "ES", "DE", "CZ")
write.table(allw, "C:/Users/sroilo/Desktop/BESTMAP documents/Data/Across_CS/All_weights_90trim_20230715.csv", sep=";", dec=",", row.names=F)

# load results from mid-trimmed datasets
wes = read.table("C:/Users/sroilo/Desktop/BESTMAP documents/Data/Across_CS/ES/July_trials/LC_trim0132_weights_20230714.csv", header=T, sep=";", dec=",")
wde = read.table("C:/Users/sroilo/Desktop/BESTMAP documents/Data/Across_CS/DE/July_trials/LC_trim0173_20230715_GAMk150_weights.csv", header=T, sep=";", dec=",")
wcz = read.table("C:/Users/sroilo/Desktop/BESTMAP documents/Data/Across_CS/CZ/July_trials/LC_trim0426_20230715_weights.csv", header=T, sep=";", dec=",")
wcz$Row.names = row.names(wcz)
allw = merge(wes, wde, by="row.names", all=T)
allw = merge(allw, wcz, by="Row.names", all=T)
names(allw) = c("Variable", "ES", "DE", "CZ")
write.table(allw, "C:/Users/sroilo/Desktop/BESTMAP documents/Data/Across_CS/All_weights_mediumTrim_20230715.csv", sep=";", dec=",", row.names=F)

# load results from light-trimmed datasets
wes = read.table("C:/Users/sroilo/Desktop/BESTMAP documents/Data/Across_CS/ES/July_trials/LC_trim0283_weights_20230714.csv", header=T, sep=";", dec=",")
wde = read.table("C:/Users/sroilo/Desktop/BESTMAP documents/Data/Across_CS/DE/July_trials/LC_trim0469_20230715_GAMk150_weights.csv", header=T, sep=";", dec=",")
wcz = read.table("C:/Users/sroilo/Desktop/BESTMAP documents/Data/Across_CS/CZ/July_trials/LC_trim0632_20230715_weights.csv", header=T, sep=";", dec=",")
wcz$Row.names = row.names(wcz)
allw = merge(wes, wde, by="row.names", all=T)
allw = merge(allw, wcz, by="Row.names", all=T)
names(allw) = c("Variable", "ES", "DE", "CZ")
write.table(allw, "C:/Users/sroilo/Desktop/BESTMAP documents/Data/Across_CS/All_weights_lightTrim_20230707.csv", sep=";", dec=",", row.names=F)



rm(list=ls())
setwd("C:/Users/sroilo/Documents")




