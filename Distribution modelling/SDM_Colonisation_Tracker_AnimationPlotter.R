################################ PLOT ANIMATION OF COLONISATION ###############################

##### Load libraries
library("geosphere")
library("animation")
library("AUC")
library("boot")
library("dbscan")
library("dismo")
library("dplyr")
library("ecospat")
#library("ENMTools")
library("FNN")
library("gam")
library("ggplot2")
library("gridExtra")
library("maxnet")
library("PresenceAbsence")
library("randomForest")
library("raster")
library("rasterVis")
library("RColorBrewer")
#library("rgdal")
#library("rgeos")

### Set working directory and import data
setwd("/Users/luqman/Desktop/SDM/Dsylvestris/Data")
#setwd("/cluster/project/gdc/people/lhirzi/SDM")

##### Define run parameters
### The number of time points define the number of time steps from the LGM to the preset. Note that PMIP LGM data refers to 21,000 years ago.
no.time.points <- 211
### max_dispersal defines the maximum distance (radius) from a cluster that new points can emerge. SDM presences beyond this max distance are ignored. I.e. this adds dispersal ability to the SDM model.
## Be aware that the effect of this parameter is time dependent, since it defines the max_dispersal distance per time step.
max_dispersal <- 0.15
### Load results from SDM_Colonisation_Tracker.R (cluster.timeseries.AllLineages.raster)
#cluster.timeseries.AllLineages.raster <- stack("TIMESERIES.DF.ANIMATION.TIMEINTERVALS211.DISPERSAL0.15.3LineagesSplit.2020-02-14_11_54_25.grd")
cluster.timeseries.AllLineages.raster <- stack("TIMESERIES.DF.ANIMATION.TIMEINTERVALS211.DISPERSAL0.15.3LineagesSplit.BAC.1.grd")

### Load the DEM background raster
GLACIER_LANDSEA_DEM_TIMESERIES <-stack("Chelsa_TraCE_Glacier_LandSea_DEM_TimeSeries.grd")
e_alps <- extent(3, 23, 37.5, 48.5)
GLACIER_LANDSEA_DEM_TIMESERIES <- crop(GLACIER_LANDSEA_DEM_TIMESERIES, e_alps)

## Define the time ranges, to show in the plot
data.range.LGM_today = seq((no.time.points-1)/10, 0, -0.1)

### Define plotting parameters. Here we have three probabilistic projections, but we want to colour them differently, i.e. we want to have three separate colour gradients for the three different lineages.
## We do this by custom defining the z-breaks and colours, according to the min and max values of the rasters.
## For the time-series, we'll want to have a fixed z-value legend, so you may want to define this as custom
min_zValue <- c(min(minValue(GLACIER_LANDSEA_DEM_TIMESERIES)[which(minValue(GLACIER_LANDSEA_DEM_TIMESERIES) >-32768)]), 5000.4, 5001.4, 5002.4)
#min_zValue <- c(-758, 5000.4, 5001.4, 5002.4)
max_zValue <- c(max(maxValue(GLACIER_LANDSEA_DEM_TIMESERIES)), 5000.9, 5001.9, 5002.9)
#max_zValue <- c(4627, 5000.9, 5001.9, 5002.9)
breakpoints_DEM <- c(min_zValue[1], min_zValue[1]+(max_zValue[1] - min_zValue[1])/9*1, min_zValue[1]+(max_zValue[1] - min_zValue[1])/9*2, min_zValue[1]+(max_zValue[1] - min_zValue[1])/9*3, min_zValue[1]+(max_zValue[1] - min_zValue[1])/9*4, min_zValue[1]+(max_zValue[1] - min_zValue[1])/9*5, min_zValue[1]+(max_zValue[1] - min_zValue[1])/9*6, min_zValue[1]+(max_zValue[1] - min_zValue[1])/9*7, min_zValue[1]+(max_zValue[1] - min_zValue[1])/9*8, min_zValue[1]+(max_zValue[1] - min_zValue[1])/9*9, 5000)
breakpoints_Apennines <- c(min_zValue[2], min_zValue[2]+(max_zValue[2] - min_zValue[2])/9*1, min_zValue[2]+(max_zValue[2] - min_zValue[2])/9*2, min_zValue[2]+(max_zValue[2] - min_zValue[2])/9*3, min_zValue[2]+(max_zValue[2] - min_zValue[2])/9*4, min_zValue[2]+(max_zValue[2] - min_zValue[2])/9*5, min_zValue[2]+(max_zValue[2] - min_zValue[2])/9*6, min_zValue[2]+(max_zValue[2] - min_zValue[2])/9*7, min_zValue[2]+(max_zValue[2] - min_zValue[2])/9*8, min_zValue[2]+(max_zValue[2] - min_zValue[2])/9*9, 5001)
breakpoints_Balkans <- c(min_zValue[3], min_zValue[3]+(max_zValue[3] - min_zValue[3])/9*1, min_zValue[3]+(max_zValue[3] - min_zValue[3])/9*2, min_zValue[3]+(max_zValue[3] - min_zValue[3])/9*3, min_zValue[3]+(max_zValue[3] - min_zValue[3])/9*4, min_zValue[3]+(max_zValue[3] - min_zValue[3])/9*5, min_zValue[3]+(max_zValue[3] - min_zValue[3])/9*6, min_zValue[3]+(max_zValue[3] - min_zValue[3])/9*7, min_zValue[3]+(max_zValue[3] - min_zValue[3])/9*8, min_zValue[3]+(max_zValue[3] - min_zValue[3])/9*9, 5002)
breakpoints_CentralAlps <- c(min_zValue[4], min_zValue[4]+(max_zValue[4] - min_zValue[4])/9*1, min_zValue[4]+(max_zValue[4] - min_zValue[4])/9*2, min_zValue[4]+(max_zValue[4] - min_zValue[4])/9*3, min_zValue[4]+(max_zValue[4] - min_zValue[4])/9*4, min_zValue[4]+(max_zValue[4] - min_zValue[4])/9*5, min_zValue[4]+(max_zValue[4] - min_zValue[4])/9*6, min_zValue[4]+(max_zValue[4] - min_zValue[4])/9*7, min_zValue[4]+(max_zValue[4] - min_zValue[4])/9*8, min_zValue[4]+(max_zValue[4] - min_zValue[4])/9*9, 5003)
breakpoints <- c(breakpoints_DEM, breakpoints_Apennines, breakpoints_Balkans, breakpoints_CentralAlps)
custom_colors <- c(brewer.pal(9, "Greys"), "grey95", "grey95", brewer.pal(9, "Blues"), "grey95", "grey95", brewer.pal(9, "Greens"), "grey95", "grey95", brewer.pal(9, "YlOrBr"), "grey95")

## To output this sequence of plots as an animation:
# You may want to plot the sea level change (as a transparent single colour layer) in your final plot!
time_stamp <- gsub(":","_",gsub(" ", "_", Sys.time()))
saveGIF({
  for (i in seq(1,no.time.points,1)) {
    ### And finally we plot. Here rasterVis levelplot looks nicer than the normal raster plot. For final plot, adjust maxpixels = 1e7
    ### Add additional legend for predicted probabilities (e.g. see: http://r-sig-geo.2731867.n2.nabble.com/Multiple-legends-with-levelplot-and-spplot-td7587300.html)
    #plot(cluster.timeseries.AllLineages.raster[[i]], breaks=breakpoints_pooled, col=custom_colors,  main="Predicted LGM Main Refugia - Pooled lineage SDM", legend.args = list(text = "Predicted Probability", side = 2, font = 4, cex = 0.9, line = 0.5), legend.width = 2, legend.mar = 5, legend.shrink = 0.9)
    #maskedCluster.lineages.rasters.transformed.final.plot <- levelplot(cluster.timeseries.AllLineages.raster[[i]], col.regions = custom_colors, at=breakpoints,  maxpixels = 1e7, margin = FALSE, main=list(label=paste("Predicted colonisation trajectory from LGM refugia (Lineage-Specific Ensemble Models, DBSCAN,  dispersal parameter = ", max_dispersal, ") - ", round(data.range.LGM_today[i], 1), "K years before present", sep = ""), cex = 2.5), legend=list(top=list(fun=grid::textGrob("Elevation", y=1, x=1.06))))
    maskedCluster.lineages.rasters.transformed.final.plot <- levelplot(cluster.timeseries.AllLineages.raster[[i]], col.regions = custom_colors, at=breakpoints,  maxpixels = 1e7, margin = FALSE, main=list(label=paste0("Reconstruction of colonisation (lineage-specific SDM w/ dispersal kernel (", max_dispersal, ") & competitive exclusion) - ", round(data.range.LGM_today[i], 1), "Kya"), cex = 3), legend=list(top=list(fun=grid::textGrob("Elevation", y=1, x=1.06))))
    plot(maskedCluster.lineages.rasters.transformed.final.plot)
    ani.pause()
  }
}, ani.height = dim(cluster.timeseries.AllLineages.raster)[1], ani.width = dim(cluster.timeseries.AllLineages.raster)[2], interval = 0.03, nmax = no.time.points, movie.name = paste0("TIMESERIES.ANIMATION.TIMEINTERVALS",no.time.points,".DISPERSAL",max_dispersal,".3LineagesSplit.",time_stamp,".gif"))



## In case you'd like to plot specific time slices e.g. the first and last time point:
# Crop rasters (if desired)
cluster.timeseries.AllLineages.raster.cropped <- crop(cluster.timeseries.AllLineages.raster, extent(3.5,22,39.25,48))
# Redefine colour palette
colfunc_yellows <- colorRampPalette(c("#FCF4EB","#F5A216","#683B00"))
colfunc_yellows(9)
custom_colors <- c(brewer.pal(9, "Greys"), "grey95", "grey95", colfunc_yellows(9), "grey95", "grey95", brewer.pal(9, "Greens"), "grey95", "grey95", brewer.pal(9, "Reds"), "grey95")
# Plot specific time slices
i <- 1
maskedCluster.lineages.rasters.transformed.final.plot <- levelplot(cluster.timeseries.AllLineages.raster.cropped[[i]], col.regions = custom_colors, at=breakpoints,  maxpixels = 1e7, margin = FALSE, legend=list(top=list(fun=grid::textGrob("Elevation", y=1, x=1.06))))
par(mar = c(2, 2, 2, 2), mfrow = c(1, 1))
plot(maskedCluster.lineages.rasters.transformed.final.plot)
