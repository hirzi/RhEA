
################################ TRACK LINEAGE-SPECIFIC EXPANSION/COLONISATION FROM REFUGIA BY CLUSTER ASSIGNMENT OF TIME SERIES DATA ###############################

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
## This is in units of latitude (at 44N) degrees, where 1 degree or unit is equivalent to 80km dispersal per 100 years (if using 100 year time intervals). E.g. 0.1 is equivalent to 80m/year. 0.1 or 0.15 seem to be appropriate values.
max_dispersal <- 0.15
### Response variable type
responce_variable_type <- "binary"
### Define CRS projection
prj_longlat <- "+init=epsg:4326"
# Choose which maxent version to run. "maxent" uses (sources) the Maxent Java package (and hence requires maxent installation for external .jar file). As an alternative (same authors, same methodology), you can now use the R package Package ‘maxnet’ (independent R package, with some slightly different implementations).
maxent_maxnet <- "maxnet"
# To change competitive/niche exclusion order, change lines: 172-175 and 384-387

##### Import environmental rasters
### Current
ENV_RASTERS <-stack("ENV_RASTERS_TODAY_EUROPE_CHELSA_LARGE.grd") 
### LGM - PMIP3
climate_PMIP3_model_selector <- 2
clim_PMIP3_model_list.names <- c("ENSEMBLE", "ENSEMBLE_NCAR_MIROC_MRI_MPI", "NCAR-CCSM4", "MIROC-ESM", "MRI-CGCM3", "MPI-ESM-P", "CESS-FGOALS-g2", "IPSL-CM5A-LR", "CNRM-CM5")
#clim_PMIP3_model_list.data <- c("ENV_RASTERS_LGM_EUROPE_PMIP_ENSEMBLE_ALL.rds", "ENV_RASTERS_LGM_EUROPE_PMIP_ENSEMBLE_NCAR_MIROC_MRI_MPI.rds", "ENV_RASTERS_LGM_EUROPE_PMIP3_NCAR_CCSM4.rds", "ENV_RASTERS_LGM_EUROPE_PMIP_MIROC_ESM.rds", "ENV_RASTERS_LGM_EUROPE_PMIP_MRI-CGCM3.rds", "ENV_RASTERS_LGM_EUROPE_PMIP_MPI_ESM_P.rds", "ENV_RASTERS_LGM_EUROPE_PMIP_CESS_FGOALS_g2.rds", "ENV_RASTERS_LGM_EUROPE_PMIP_IPSL_CM5A_LR.rds", "ENV_RASTERS_LGM_EUROPE_PMIP_CNRM_CM5.rds")
clim_PMIP3_model_list.data <- c("ENV_RASTERS_LGM_EUROPE_PMIP_ENSEMBLE_ALL_LARGE.grd", "ENV_RASTERS_LGM_EUROPE_PMIP_ENSEMBLE_NCAR_MIROC_MRI_MPI_LARGE.grd", "ENV_RASTERS_LGM_EUROPE_PMIP3_NCAR_CCSM4_LARGE.grd", "ENV_RASTERS_LGM_EUROPE_PMIP3_MIROC_ESM_LARGE.grd", "ENV_RASTERS_LGM_EUROPE_PMIP3_MRI-CGCM3_LARGE.grd", "ENV_RASTERS_LGM_EUROPE_PMIP3_MPI_ESM_P_LARGE.grd", "ENV_RASTERS_LGM_EUROPE_PMIP3_CESS_FGOALS_g2_LARGE.grd", "ENV_RASTERS_LGM_EUROPE_PMIP3_IPSL_CM5A_LR_LARGE.grd", "ENV_RASTERS_LGM_EUROPE_PMIP3_CNRM_CM5_LARGE.grd")
climate_model_PMIP3 <- clim_PMIP3_model_list.data[climate_PMIP3_model_selector]
#climate_model_PMIP3 <- "ENV_RASTERS_EUROPE_CHELSA_TraCE_XLARGE_-192.grd"
ENV_LGM_RASTERS <- stack(climate_model_PMIP3)
### And crop to desired extent
e_alps <- extent(3, 23, 37.5, 48.5)
ENV_RASTERS <- crop(ENV_RASTERS, e_alps)
ENV_LGM_RASTERS <- crop(ENV_LGM_RASTERS, e_alps)
### Define the DEM background raster for the intial time point
## Fixed DEM (Chelsa PaleoDEM)
# base_map_dem_init <- ENV_LGM_RASTERS$ELEV.GMTED
## Chelsa TraCE glacier plus landsea time series DEM
GLACIER_LANDSEA_DEM_TIMESERIES <-stack("Chelsa_TraCE_Glacier_LandSea_DEM_TimeSeries.grd")
GLACIER_LANDSEA_DEM_TIMESERIES <- crop(GLACIER_LANDSEA_DEM_TIMESERIES, e_alps)
base_map_dem_init <- GLACIER_LANDSEA_DEM_TIMESERIES[[1]]
base_map_dem_init[base_map_dem_init == -32768] <- NA # assign NA value

##### Initialise empty list for dataframes and progress bar
cluster.timeseries.Apennines.raster <- stack()
cluster.timeseries.Apennines.df <- list()
cluster.timeseries.Balkans.raster <- stack()
cluster.timeseries.Balkans.df <- list()
cluster.timeseries.CentralAlps.raster <- stack()
cluster.timeseries.CentralAlps.df <- list()
cluster.timeseries.AllLineages.raster <- stack()


################################ INITIAL (LGM) TIME POINT, T=0  ################################

##### We first want to extract the relevant clusters from the DBScan results, and fill in the initial time point in our time series rasters
### Make empty stack to append results (across the three lineages)
lineage_cluster.intial.rasters <- stack()
### Start count for for loop.
count_initial <- 1

### Loop across the lineages
for (lineage in c("Apennines", "Balkans", "CentralAlps")) {
  
  ## Load saved outputs of the SDM probabilistic projection (raster) and DBScan clustering results (dataframe), from running SDM_Dsyl_X.R script.
  if (lineage == "Apennines") {
    pred.bin.LGM.clusters.df <- readRDS("pred.bin.LGM.clusters.df.Apennines.selectVar.optimalThreshold.2020-02-13_16_11_12.rds")
    pred.prob.LGM <- readRDS("pred.prob.LGM.Apennines.selectVar.optimalThreshold.2020-02-13_16_11_12.rds")
  } else if (lineage == "Balkans") {
    pred.bin.LGM.clusters.df <- readRDS("pred.bin.LGM.clusters.df.Balkans.selectVar.optimalThreshold.2020-02-13_15_22_23.rds")
    pred.prob.LGM <- readRDS("pred.prob.LGM.Balkans.selectVar.optimalThreshold.2020-02-13_15_22_23.rds")
  } else if (lineage == "CentralAlps") {
    pred.bin.LGM.clusters.df <- readRDS("pred.bin.LGM.clusters.df.CentralAlps.selectVar.optimalThreshold.2020-02-13_14_51_11.rds")
    pred.prob.LGM <- readRDS("pred.prob.LGM.CentralAlps.selectVar.optimalThreshold.2020-02-13_14_51_11.rds")
  }
  
  ## Transform dataframe to spatial points dataframe, and apply CRS projection
  pred.bin.LGM.clusters <- pred.bin.LGM.clusters.df
  coordinates(pred.bin.LGM.clusters) <- ~x+y
  proj4string(pred.bin.LGM.clusters) = CRS(prj_longlat)
  
  ## Extract relevant cluster. Here we do so automatically, rather than via manual selection, via defining focal geographic points.
  ## First we define origins around which we'll propose selection of cluster.
  if (lineage == "Apennines") {
    cluster_origin_x <- 11.5
    cluster_origin_y <- 44
  } else if (lineage == "Balkans") {
    cluster_origin_x <- 20
    cluster_origin_y <- 42.5
  } else if (lineage == "CentralAlps") {
    cluster_origin_x <- 10.5
    cluster_origin_y <- 45.75
  }
  cluster_origin_xy <- data.frame(cluster_origin_x, cluster_origin_y)
  colnames(cluster_origin_xy) <- c("x", "y")
  
  ## We calculate the distances between the defined origins (focal point) and the data points
  dist_from_origin <- distm(cluster_origin_xy, pred.bin.LGM.clusters.df[,c("x","y")], fun = distGeo)
  
  ## Extract which cluster is closest to defined origin point
  closest_to_origin_x_idx <- which(dist_from_origin == min(dist_from_origin))
  cluster_origin_xy$cluster <- pred.bin.LGM.clusters.df$cluster[max.col(-dist_from_origin)]
  closest_to_origin_x <- pred.bin.LGM.clusters.df$cluster[closest_to_origin_x_idx]
  pred.bin.LGM.single_cluster.df <- pred.bin.LGM.clusters.df[pred.bin.LGM.clusters.df$cluster == closest_to_origin_x,]
  
  ## Convert to spatial points dataframe and apply CRS projection.
  pred.bin.LGM.single_cluster <- pred.bin.LGM.single_cluster.df
  coordinates(pred.bin.LGM.single_cluster) <- ~x+y
  proj4string(pred.bin.LGM.single_cluster) = CRS(prj_longlat)
  ## Convert SPDF to raster
  r <- raster(nrows=nrow(ENV_RASTERS$ISOTHERMALITY), ncols=ncol(ENV_RASTERS$ISOTHERMALITY), ext=e_alps)
  #proj4string(r) = CRS(prj_longlat)
  pred.bin.LGM.single_cluster.raster <- rasterize(pred.bin.LGM.single_cluster, r, background=0, field = 1)
  
  ## Append lineage specific results to list (inital time point)
  if (lineage == "Apennines") {
    cluster.timeseries.Apennines.raster <- stack(cluster.timeseries.Apennines.raster, pred.bin.LGM.single_cluster.raster)
    cluster.timeseries.Apennines.df[[1]] <- pred.bin.LGM.single_cluster.df
  } else if (lineage == "Balkans") {
    cluster.timeseries.Balkans.raster <- stack(cluster.timeseries.Balkans.raster, pred.bin.LGM.single_cluster.raster)
    cluster.timeseries.Balkans.df[[1]] <- pred.bin.LGM.single_cluster.df
  } else if (lineage == "CentralAlps") {
    cluster.timeseries.CentralAlps.raster <- stack(cluster.timeseries.CentralAlps.raster, pred.bin.LGM.single_cluster.raster)
    cluster.timeseries.CentralAlps.df[[1]] <- pred.bin.LGM.single_cluster.df
  }
  
  ## We now mask the probabilistic projection by this single cluster spdf.
  pred.prob.LGM.masked_Cluster <- mask(pred.prob.LGM, pred.bin.LGM.single_cluster)
  
  ## Finally we append the results to the list and stack.
  lineage_cluster.intial.rasters <- stack(lineage_cluster.intial.rasters, pred.prob.LGM.masked_Cluster)
  names(lineage_cluster.intial.rasters)[count_initial] <- lineage
  ## And update the count. 
  count_initial <- count_initial + 1
  
}

### Let's now differentiate the raster z values by adding 1 and 2 to layers 2 and 3. We do this to differentiate the lineage colours for plotting later. We further add 5000 to be able to later superimpose and distinguish it from the base DEM layer.
lineage_cluster.intial.rasters.transformed <- lineage_cluster.intial.rasters
lineage_cluster.intial.rasters.transformed$Apennines <- lineage_cluster.intial.rasters.transformed$Apennines + 5000
lineage_cluster.intial.rasters.transformed$Balkans <- lineage_cluster.intial.rasters.transformed$Balkans + 5000 + 1
lineage_cluster.intial.rasters.transformed$CentralAlps <- lineage_cluster.intial.rasters.transformed$CentralAlps + 5000 + 2

### Add/superimpose the 3 lineage rasters and base map together, and add CRS projection.
## Use the raster function cover to combine the 3 masked rasters in the desired order (Balkans > Appenines > CentralAlps), add base map together, and add CRS projection.
## Note that the order or sequence of "covering" dictates the order of competitive exclusion among the lineages such that e.g. here Balkans > Appenines > CentralAlps or Balkans > CentralAlps > Appenines, where > = outcompetes. As an alternative, consider a order-free approach.
lineage_cluster.intial.rasters.transformed_temp1 <- cover(lineage_cluster.intial.rasters.transformed$Balkans, lineage_cluster.intial.rasters.transformed$Apennines)
#lineage_cluster.intial.rasters.transformed_temp1 <- cover(lineage_cluster.intial.rasters.transformed$Balkans, lineage_cluster.intial.rasters.transformed$CentralAlps)
lineage_cluster.intial.rasters.transformed_temp2 <- cover(lineage_cluster.intial.rasters.transformed_temp1, lineage_cluster.intial.rasters.transformed$CentralAlps)
#lineage_cluster.intial.rasters.transformed_temp2 <- cover(lineage_cluster.intial.rasters.transformed_temp1, lineage_cluster.intial.rasters.transformed$Apennines)
lineage_cluster.intial.rasters.transformed.final <- cover(lineage_cluster.intial.rasters.transformed_temp2, base_map_dem_init)
proj4string(lineage_cluster.intial.rasters.transformed.final) = CRS(prj_longlat)

### And append results
cluster.timeseries.AllLineages.raster <- stack(cluster.timeseries.AllLineages.raster, lineage_cluster.intial.rasters.transformed.final)


################################ SUBSEQUENT TIME POINTS, T > 0 ################################

##### Now we perform lineage-specific SDMs by applying the estimated SDM models on the time series data. We then add a dispersal kernel (via knnxdist) and niche exclusion (raster:cover).
### Let's add a progress bar (updated per time step)
pb = txtProgressBar(min = 0, max = no.time.points, initial = 0)

### Loop across time steps
i <- 1
for (i in seq(1,no.time.points-1,1)) {
  
  ## Update environmental rasters. Here, we do a naive (linear) interpolation. Recall that the first time point has already been added, so we start from 1 rather than 0.
  ISOTHERMALITY <- ENV_LGM_RASTERS$ISOTHERMALITY + (i * ((ENV_RASTERS$ISOTHERMALITY - ENV_LGM_RASTERS$ISOTHERMALITY) / (no.time.points-1) ))
  TEMP_SEASONALITY <- ENV_LGM_RASTERS$TEMP_SEASONALITY + (i * ((ENV_RASTERS$TEMP_SEASONALITY - ENV_LGM_RASTERS$TEMP_SEASONALITY) / (no.time.points-1) ))
  TEMP_MAX_WARMEST_MONTH <- ENV_LGM_RASTERS$TEMP_MAX_WARMEST_MONTH + (i * ((ENV_RASTERS$TEMP_MAX_WARMEST_MONTH - ENV_LGM_RASTERS$TEMP_MAX_WARMEST_MONTH) / (no.time.points-1) ))
  TEMP_MEAN_WETTEST_QUARTER <- ENV_LGM_RASTERS$TEMP_MEAN_WETTEST_QUARTER + (i * ((ENV_RASTERS$TEMP_MEAN_WETTEST_QUARTER - ENV_LGM_RASTERS$TEMP_MEAN_WETTEST_QUARTER) / (no.time.points-1) ))
  TEMP_MEAN_DRIEST_QUARTER <- ENV_LGM_RASTERS$TEMP_MEAN_DRIEST_QUARTER + (i * ((ENV_RASTERS$TEMP_MEAN_DRIEST_QUARTER - ENV_LGM_RASTERS$TEMP_MEAN_DRIEST_QUARTER) / (no.time.points-1) ))
  PREC_SEASONALITY <- ENV_LGM_RASTERS$PREC_SEASONALITY + (i * ((ENV_RASTERS$PREC_SEASONALITY - ENV_LGM_RASTERS$PREC_SEASONALITY) / (no.time.points-1) ))
  PREC_WARMEST_QUARTER <- ENV_LGM_RASTERS$PREC_WARMEST_QUARTER + (i * ((ENV_RASTERS$PREC_WARMEST_QUARTER - ENV_LGM_RASTERS$PREC_WARMEST_QUARTER) / (no.time.points-1) ))
  PREC_COLDEST_QUARTER <- ENV_LGM_RASTERS$PREC_COLDEST_QUARTER + (i * ((ENV_RASTERS$PREC_COLDEST_QUARTER - ENV_LGM_RASTERS$PREC_COLDEST_QUARTER) / (no.time.points-1) ))
  ENV_RASTERS_TIMESERIES <- stack(ISOTHERMALITY, TEMP_SEASONALITY, TEMP_MAX_WARMEST_MONTH, TEMP_MEAN_WETTEST_QUARTER, TEMP_MEAN_DRIEST_QUARTER, PREC_SEASONALITY, PREC_WARMEST_QUARTER, PREC_COLDEST_QUARTER, ENV_RASTERS$PH_5cm, ENV_RASTERS$SLOPE)
  names(ENV_RASTERS_TIMESERIES) <- c("ISOTHERMALITY","TEMP_SEASONALITY","TEMP_MAX_WARMEST_MONTH","TEMP_MEAN_WETTEST_QUARTER","TEMP_MEAN_DRIEST_QUARTER","PREC_SEASONALITY","PREC_WARMEST_QUARTER","PREC_COLDEST_QUARTER","PH_5cm","SLOPE")
  
  ##### Load a DEM background raster
  ## Fixed DEM (Chelsa PaleoDEM)
  # ELEV.GMTED <- ENV_LGM_RASTERS$ELEV.GMTED + (i * ((ENV_RASTERS$ELEV.GMTED - ENV_LGM_RASTERS$ELEV.GMTED) / (no.time.points-1) )) # for base DEM map
  # base_map_dem <- ELEV.GMTED
  ## Chelsa TraCE glacier plus landsea time series DEM
  base_map_dem <- GLACIER_LANDSEA_DEM_TIMESERIES[[i+1]] # this indexing "i+1" is only applicable when no.time.points = 211
  base_map_dem[base_map_dem == -32768] <- NA # assign NA value
  
  ## Make the LGM environmental projection dataframe
  Projection.timeseries<- as.data.frame(rasterToPoints(ENV_RASTERS_TIMESERIES))
  Projection.timeseries <- na.omit(Projection.timeseries)
  Projection.timeseries.raster <- rasterFromXYZ(Projection.timeseries) 
  
  ### First we'll need to project the lineage-specific models (with dispersal kernel) for each lineage separately. We do this in the following for loop. 
  ## Define for loop variables (lineages)
  #maskedCluster.lineages.list <- list()
  maskedCluster.lineages.rasters <- stack()
  count <- 1
  ## Loop across lineages
  for (lineage in c("Apennines", "Balkans", "CentralAlps")) {
    
    if (lineage == "Apennines") {
      ## The below RData file should load: glm.multi.step, gam.step, rf.prob, maxent.prob, thres.glm.multi.step, thres.gam.step, thres.rf, thres.maxent, thres.mean, auc_glm, auc_gam, auc_maxent, auc_rf, auc_mean (see: http://www.sthda.com/english/wiki/saving-data-into-r-data-format-rds-and-rdata)
      load("colonisation_tracker_metaData_Apennines_2020-02-13_16_11_12.RData")
      cluster.timeseries.raster <- cluster.timeseries.Apennines.raster
      cluster.timeseries.df <- cluster.timeseries.Apennines.df
      dist_correction_factor <- 111/80
    } else if (lineage == "Balkans") {
      load("colonisation_tracker_metaData_Balkans_2020-02-13_15_22_23.RData")
      cluster.timeseries.raster <- cluster.timeseries.Balkans.raster
      cluster.timeseries.df <- cluster.timeseries.Balkans.df
      dist_correction_factor <- 111/82
    } else if (lineage == "CentralAlps") {
      load("colonisation_tracker_metaData_CentralAlps_2020-02-13_14_51_11.RData")
      cluster.timeseries.raster <- cluster.timeseries.CentralAlps.raster
      cluster.timeseries.df <- cluster.timeseries.CentralAlps.df
      dist_correction_factor <- 111/78
    }
    
    ## Predict onto projection raster
    ## Sometimes an error pops up here for no real reason. Close R and rerun if this is the case.
    pred.glm.prob.timeseries.currentstep <- predict(object = Projection.timeseries.raster, model = glm.multi.step, type = "response")
    pred.gam.prob.timeseries.currentstep <- predict(object = Projection.timeseries.raster, model = gam.step, type = "response")
    pred.rf.prob.timeseries.currentstep <- predict(object = Projection.timeseries.raster, model = rf.prob, type = "response")
    if (responce_variable_type == "binary") {
      if (maxent_maxnet == "maxent") {
        pred.maxent.prob.timeseries.currentstep <- predict(object = Projection.timeseries.raster, model = maxent.prob, type = "response")
      } else if (maxent_maxnet == "maxnet") {
        pred.maxent.prob.timeseries.currentstep <- predict(object = Projection.timeseries.raster, model = maxent.prob, type = "cloglog")  # cloglog (default) and logistic give very similar results
      }
      pred.maxent.bin.timeseries.currentstep <- (pred.maxent.prob.timeseries.currentstep > thres.maxent)
      pred.glm.bin.timeseries.currentstep <- pred.glm.prob.timeseries.currentstep > thres.glm.multi.step[4, 2]
      pred.gam.bin.timeseries.currentstep <- pred.gam.prob.timeseries.currentstep > thres.gam.step[4, 2]
      pred.rf.bin.timeseries.currentstep <- (pred.rf.prob.timeseries.currentstep > thres.rf)
    }
    ## Then for the ensemble model
    if (responce_variable_type == "binary") {
      pred.prob.timeseries.currentstep <-(auc_glm*pred.glm.prob.timeseries.currentstep + auc_gam*pred.gam.prob.timeseries.currentstep + auc_maxent*pred.maxent.prob.timeseries.currentstep + auc_rf*pred.rf.prob.timeseries.currentstep) / (auc_mean*4)
      pred.bin.timeseries.currentstep <- pred.prob.timeseries.currentstep > thres.mean
    } else if (responce_variable_type == "frequency") {
      pred.prob.timeseries.currentstep <-(pred.glm.prob.timeseries.currentstep + pred.gam.prob.timeseries.currentstep + pred.rf.prob.timeseries.currentstep) / (3)
    }
    
    ## Calculate difference in predictions of t and t+1
    pred.bin.timeseries.laststep <- cluster.timeseries.raster[[i]]
    pred.bin.timeseries.laststep[cluster.timeseries.raster[[i]] > 0] <- 2
    turnover.bin <- pred.bin.timeseries.currentstep + pred.bin.timeseries.laststep
    #freq(turnover.bin)
    ## We generate a turnover/range change table
    #table(getValues(pred.bin.timeseries.currentstep),getValues(pred.bin.timeseries.laststep))
    
    ## Calculate positive turnover (i.e. colonisation points) and negative turnover (i.e. contraction points)
    positive.turnover <- turnover.bin
    positive.turnover[positive.turnover != 1] <- NA
    positive.turnover.df <- as.data.frame(rasterToPoints(positive.turnover))
    positive.turnover.df <- positive.turnover.df[,c(1,2)]
    negative.turnover <- turnover.bin
    negative.turnover[negative.turnover != 2] <- NA
    negative.turnover.df <- as.data.frame(rasterToPoints(negative.turnover))
    negative.turnover.df <- negative.turnover.df[,c(1,2)]
    
    ## Perform k-Nearest Neighbour Classification, conditional upon maximum dispersal distance beyond which classification (assignment to cluster) is not performed and presence is ignored (NA)
    ## Remove unassigned individuals (annotated as cluster 0) from the training data set.
    pred.bin.clusters.timeseries.laststep.train <- cluster.timeseries.df[[i]]
    pred.bin.clusters.timeseries.laststep.train[pred.bin.clusters.timeseries.laststep.train == 0] <- NA
    pred.bin.clusters.timeseries.laststep.train <- na.omit(pred.bin.clusters.timeseries.laststep.train)
    pred.bin.clusters.timeseries.laststep.train.coords <- pred.bin.clusters.timeseries.laststep.train[,c(1,2)]
    pred.bin.clusters.timeseries.laststep.train.labels <- pred.bin.clusters.timeseries.laststep.train[,c(3)]
    ## Condition on maximum dispersal distance.
    # Note 1: we make the approximation of Euclidean space, which is fair for small geographic scales where the curvature is minimal.
    # Note 2: however, we correct for the fact that at 44 degrees latitude and 10 degrees longitude, 1 degree latitude is approx 11/8 times longer than 1 degree longitude (111km vs 80km; see https://www.nhc.noaa.gov/gccalc.shtml) 
    positive.turnover.df.adj <- positive.turnover.df
    positive.turnover.df.adj$x <- positive.turnover.df.adj$x * (dist_correction_factor)
    pred.bin.clusters.timeseries.laststep.train.coords.adj <- pred.bin.clusters.timeseries.laststep.train.coords
    pred.bin.clusters.timeseries.laststep.train.coords.adj$x <- pred.bin.clusters.timeseries.laststep.train.coords.adj$x * (dist_correction_factor)
    knn_dist.currentstep <- knnx.dist(data = pred.bin.clusters.timeseries.laststep.train.coords.adj, query = positive.turnover.df.adj, k=1)
    positive.turnover.train <- positive.turnover.df
    positive.turnover.train["NN_dist"] <- knn_dist.currentstep
    positive.turnover.train <- positive.turnover.train[positive.turnover.train$NN_dist < max_dispersal, ]
    positive.turnover.train <- positive.turnover.train[,c(1,2)]
    ## # kNN classification. Since we're dealing with each lineage specifically, this is no longer needed. DEPRECIATED
    # knn_timeseries.currentstep <- knn(train = pred.bin.clusters.timeseries.laststep.train.coords, test = positive.turnover.train, cl = pred.bin.clusters.timeseries.laststep.train.labels, k=3)
    # positive.turnover.cluster.temp <- positive.turnover.train
    # positive.turnover.cluster.temp["cluster"] <- as.integer(knn_timeseries.currentstep)
    # positive.turnover.unassigned <- anti_join(positive.turnover.df, positive.turnover.train, by = c("x", "y"))
    # if (nrow(positive.turnover.unassigned) > 0) {
    #   positive.turnover.unassigned["cluster"] <- as.integer(0)
    #   positive.turnover.cluster <- rbind(positive.turnover.cluster.temp, positive.turnover.unassigned)
    # } else if (nrow(positive.turnover.unassigned) == 0) {
    #   positive.turnover.cluster <- positive.turnover.cluster.temp
    # }
    
    ## Make updated cluster-assigned dataframe for current step. First we add the positive turnover to the presences of the last time step, then we subtract the negative turnover; to obtain all cluster-assigned presences of the current time step. the 
    positive.turnover.cluster <- positive.turnover.train
    positive.turnover.cluster["cluster"] <- cluster.timeseries.df[[i]][3][1,]
    pred.bin.clusters.timeseries.currentstep.temp <- rbind(cluster.timeseries.df[[i]], positive.turnover.cluster)
    pred.bin.clusters.timeseries.currentstep <- anti_join(pred.bin.clusters.timeseries.currentstep.temp,negative.turnover.df, by = c("x", "y"))
    
    ## Convert DF to SpatialPointsDF
    pred.bin.clusters.timeseries.currentstep.sp <- pred.bin.clusters.timeseries.currentstep
    coordinates(pred.bin.clusters.timeseries.currentstep.sp) <- ~x+y
    proj4string(pred.bin.clusters.timeseries.currentstep.sp) = CRS(prj_longlat)
    ## Convert SPDF to raster
    #r <- raster(nrows=nrow(ENV_RASTERS$ISOTHERMALITY), ncols=ncol(ENV_RASTERS$ISOTHERMALITY), ext=e_alps)
    #proj4string(r) = CRS(prj_longlat)
    pred.bin.clusters.timeseries.currentstep.raster <- rasterize(pred.bin.clusters.timeseries.currentstep.sp, r, background=0, field = 1)
    
    ## Mask the probabilistic projection by this single cluster spdf.
    pred.prob.maskedCluster.timeseries.currentstep.raster <- mask(pred.prob.timeseries.currentstep, pred.bin.clusters.timeseries.currentstep.sp)
    ## Convert masked raster to dataframe - ONLY FOR OPTION 1 BELOW!
    #pred.prob.maskedCluster.timeseries.currentstep.df <- as.data.frame(rasterToPoints(pred.prob.maskedCluster.timeseries.currentstep.raster))
    
    ## Assign result to loop variable
    #maskedCluster.lineages.list[[count]] <- pred.prob.maskedCluster.timeseries.currentstep.df # OPTION 1
    #names(maskedCluster.lineages.list) <- lineage
    maskedCluster.lineages.rasters <- stack(maskedCluster.lineages.rasters, pred.prob.maskedCluster.timeseries.currentstep.raster) # OPTION 2
    names(maskedCluster.lineages.rasters)[count] <- lineage
    
    ## Assign lists
    if (lineage == "Apennines") {
      cluster.timeseries.Apennines.raster <- stack(cluster.timeseries.Apennines.raster, pred.bin.clusters.timeseries.currentstep.raster)
      cluster.timeseries.Apennines.df[[i+1]] <- pred.bin.clusters.timeseries.currentstep
    } else if (lineage == "Balkans") {
      cluster.timeseries.Balkans.raster <- stack(cluster.timeseries.Balkans.raster, pred.bin.clusters.timeseries.currentstep.raster)
      cluster.timeseries.Balkans.df[[i+1]] <- pred.bin.clusters.timeseries.currentstep
    } else if (lineage == "CentralAlps") {
      cluster.timeseries.CentralAlps.raster <- stack(cluster.timeseries.CentralAlps.raster, pred.bin.clusters.timeseries.currentstep.raster)
      cluster.timeseries.CentralAlps.df[[i+1]] <- pred.bin.clusters.timeseries.currentstep
    }
    
    ## Update count
    count = count + 1
    
    ## Remove (reset) variables
    rm(glm.multi.step, gam.step, rf.prob, maxent.prob, thres.glm.multi.step, thres.gam.step, thres.rf, thres.maxent, thres.mean, auc_glm, auc_gam, auc_maxent, auc_rf, auc_mean)
    ## Remove (reset) variables
    rm(dist_correction_factor, cluster.timeseries.raster, cluster.timeseries.df)
    
  }
  
  ##### Now we combine the lineage-specific results into one raster, which we'll append to the time series results.
  
  ### To do this we have two options.
  ### OPTION 1 (UNFINISHED)
  ## Combine the results of the separate lineages together
  ## cbind the three  pred.prob.LGM.maskedCluster.timeseries.currentstep.df dataframes (first adding new column to indicate lineage), and identify duplicate points.
  ## if a point is duplicate, select only one (among the duplicates) based on competitive exclusion order (Apennines -> Balkans -> CentralAlps)
  ## Remove lineage column (entries are now distinguishable exclusively by z values)  and convert combined 3 column dataframe to raster and add CRS projection
  
  ### OPTION 2 (WORKING,PREFERRED)
  ## Differentiate the raster z (probability) values by adding 1 and 2 to lineages 2 and 3. We do this to differentiate the lineage colours for plotting later. We further add 5000 to be able to later superimpose and distinguish it from the base DEM layer.
  maskedCluster.lineages.rasters.transformed <- maskedCluster.lineages.rasters
  maskedCluster.lineages.rasters.transformed$Apennines <- maskedCluster.lineages.rasters.transformed$Apennines + 5000
  maskedCluster.lineages.rasters.transformed$Balkans <- maskedCluster.lineages.rasters.transformed$Balkans + 5000 + 1
  maskedCluster.lineages.rasters.transformed$CentralAlps <- maskedCluster.lineages.rasters.transformed$CentralAlps + 5000 + 2
  
  ### Combine the results of the separate lineages together
  ## Use the raster function cover to combine the 3 masked rasters in the desired order (Balkans > Appenines > CentralAlps), add base map together, and add CRS projection.
  ## Note that the order or sequence of "covering" dictates the order of competitive exclusion among the lineages such that e.g. here Balkans > Appenines > CentralAlps or Balkans > CentralAlps > Appenines, where > = outcompetes. As an alternative, consider a order-free approach
  maskedCluster.lineages.rasters.transformed_temp1 <- cover(maskedCluster.lineages.rasters.transformed$Balkans, maskedCluster.lineages.rasters.transformed$Apennines)
  #maskedCluster.lineages.rasters.transformed_temp1 <- cover(maskedCluster.lineages.rasters.transformed$Balkans, maskedCluster.lineages.rasters.transformed$CentralAlps)
  maskedCluster.lineages.rasters.transformed_temp2 <- cover(maskedCluster.lineages.rasters.transformed_temp1, maskedCluster.lineages.rasters.transformed$CentralAlps)
  #maskedCluster.lineages.rasters.transformed_temp2 <- cover(maskedCluster.lineages.rasters.transformed_temp1, maskedCluster.lineages.rasters.transformed$Apennines)
  maskedCluster.lineages.rasters.transformed.final <- cover(maskedCluster.lineages.rasters.transformed_temp2, base_map_dem)
  proj4string(maskedCluster.lineages.rasters.transformed.final) = CRS(prj_longlat)
  
  ### And finally we append the results to the result lists and stacks.
  cluster.timeseries.AllLineages.raster <- stack(cluster.timeseries.AllLineages.raster, maskedCluster.lineages.rasters.transformed.final)
  
  ### Update progress bar
  setTxtProgressBar(pb,i)
}

# Output results
time_stamp <- gsub(":","_",gsub(" ", "_", Sys.time()))
#saveRDS(cluster.timeseries.AllLineages.df, paste0("TIMESERIES.DF.ANIMATION.",no.time.points,".DISPERSAL",max_dispersal,"3LineagesSplit.rds"))
writeRaster(cluster.timeseries.AllLineages.raster, paste0("TIMESERIES.DF.ANIMATION.TIMEINTERVALS",no.time.points,".DISPERSAL",max_dispersal,".3LineagesSplit.",time_stamp))

################################ PLOT ANIMATION OF COLONISATION ###############################

# ## Define the time ranges, to show in the plot
# data.range.LGM_today = seq((no.time.points-1)/10, 0, -0.1)
# 
# ### Define plotting parameters. Here we have three probabilistic projections, but we want to colour them differently, i.e. we want to have three separate colour gradients for the three different lineages.
# ## We do this by custom defining the z-breaks and colours, according to the min and max values of the rasters.
# ## For the time-series, we'll want to have a fixed z-value legend, so you may want to define this as custom
# min_zValue <- c(min(minValue(GLACIER_LANDSEA_DEM_TIMESERIES)[which(minValue(GLACIER_LANDSEA_DEM_TIMESERIES) >-32768)]), 5000.4, 5001.4, 5002.4)
# #min_zValue <- c(-758, 5000.4, 5001.4, 5002.4)
# max_zValue <- c(max(maxValue(GLACIER_LANDSEA_DEM_TIMESERIES)), 5000.9, 5001.9, 5002.9)
# #max_zValue <- c(4627, 5000.9, 5001.9, 5002.9)
# breakpoints_DEM <- c(min_zValue[1], min_zValue[1]+(max_zValue[1] - min_zValue[1])/9*1, min_zValue[1]+(max_zValue[1] - min_zValue[1])/9*2, min_zValue[1]+(max_zValue[1] - min_zValue[1])/9*3, min_zValue[1]+(max_zValue[1] - min_zValue[1])/9*4, min_zValue[1]+(max_zValue[1] - min_zValue[1])/9*5, min_zValue[1]+(max_zValue[1] - min_zValue[1])/9*6, min_zValue[1]+(max_zValue[1] - min_zValue[1])/9*7, min_zValue[1]+(max_zValue[1] - min_zValue[1])/9*8, min_zValue[1]+(max_zValue[1] - min_zValue[1])/9*9, 5000)
# breakpoints_Apennines <- c(min_zValue[2], min_zValue[2]+(max_zValue[2] - min_zValue[2])/9*1, min_zValue[2]+(max_zValue[2] - min_zValue[2])/9*2, min_zValue[2]+(max_zValue[2] - min_zValue[2])/9*3, min_zValue[2]+(max_zValue[2] - min_zValue[2])/9*4, min_zValue[2]+(max_zValue[2] - min_zValue[2])/9*5, min_zValue[2]+(max_zValue[2] - min_zValue[2])/9*6, min_zValue[2]+(max_zValue[2] - min_zValue[2])/9*7, min_zValue[2]+(max_zValue[2] - min_zValue[2])/9*8, min_zValue[2]+(max_zValue[2] - min_zValue[2])/9*9, 5001)
# breakpoints_Balkans <- c(min_zValue[3], min_zValue[3]+(max_zValue[3] - min_zValue[3])/9*1, min_zValue[3]+(max_zValue[3] - min_zValue[3])/9*2, min_zValue[3]+(max_zValue[3] - min_zValue[3])/9*3, min_zValue[3]+(max_zValue[3] - min_zValue[3])/9*4, min_zValue[3]+(max_zValue[3] - min_zValue[3])/9*5, min_zValue[3]+(max_zValue[3] - min_zValue[3])/9*6, min_zValue[3]+(max_zValue[3] - min_zValue[3])/9*7, min_zValue[3]+(max_zValue[3] - min_zValue[3])/9*8, min_zValue[3]+(max_zValue[3] - min_zValue[3])/9*9, 5002)
# breakpoints_CentralAlps <- c(min_zValue[4], min_zValue[4]+(max_zValue[4] - min_zValue[4])/9*1, min_zValue[4]+(max_zValue[4] - min_zValue[4])/9*2, min_zValue[4]+(max_zValue[4] - min_zValue[4])/9*3, min_zValue[4]+(max_zValue[4] - min_zValue[4])/9*4, min_zValue[4]+(max_zValue[4] - min_zValue[4])/9*5, min_zValue[4]+(max_zValue[4] - min_zValue[4])/9*6, min_zValue[4]+(max_zValue[4] - min_zValue[4])/9*7, min_zValue[4]+(max_zValue[4] - min_zValue[4])/9*8, min_zValue[4]+(max_zValue[4] - min_zValue[4])/9*9, 5003)
# breakpoints <- c(breakpoints_DEM, breakpoints_Apennines, breakpoints_Balkans, breakpoints_CentralAlps)
# custom_colors <- c(brewer.pal(9, "Greys"), "grey95", "grey95", brewer.pal(9, "Blues"), "grey95", "grey95", brewer.pal(9, "Greens"), "grey95", "grey95", brewer.pal(9, "YlOrBr"), "grey95")
# 
# ## To output this sequence of plots as an animation:
# # You may want to plot the sea level change (as a transparent single colour layer) in your final plot!
# time_stamp <- gsub(":","_",gsub(" ", "_", Sys.time()))
# saveGIF({
#   for (i in seq(1,no.time.points,1)) {
#     ### And finally we plot. Here rasterVis levelplot looks nicer than the normal raster plot. For final plot, adjust maxpixels = 1e7
#     ### Add additional legend for predicted probabilities (e.g. see: http://r-sig-geo.2731867.n2.nabble.com/Multiple-legends-with-levelplot-and-spplot-td7587300.html)
#     #plot(cluster.timeseries.AllLineages.raster[[i]], breaks=breakpoints_pooled, col=custom_colors,  main="Predicted LGM Main Refugia - Pooled lineage SDM", legend.args = list(text = "Predicted Probability", side = 2, font = 4, cex = 0.9, line = 0.5), legend.width = 2, legend.mar = 5, legend.shrink = 0.9)
#     maskedCluster.lineages.rasters.transformed.final.plot <- levelplot(cluster.timeseries.AllLineages.raster[[i]], col.regions = custom_colors, at=breakpoints,  maxpixels = 1e7, margin = FALSE, main=list(label=paste("Predicted colonisation trajectory from LGM refugia (Lineage-Specific Ensemble Models, DBSCAN,  dispersal parameter = ", max_dispersal, ") - ", round(data.range.LGM_today[i], 1), "K years before present", sep = ""), cex = 2.5), legend=list(top=list(fun=grid::textGrob("Elevation", y=1, x=1.06))))
#     plot(maskedCluster.lineages.rasters.transformed.final.plot)
#     ani.pause()
#   }
# }, ani.height = dim(cluster.timeseries.AllLineages.raster)[1], ani.width = dim(cluster.timeseries.AllLineages.raster)[2], interval = 0.03, nmax = no.time.points, movie.name = paste0("TIMESERIES.ANIMATION.TIMEINTERVALS",no.time.points,".DISPERSAL",max_dispersal,".3LineagesSplit.",time_stamp,".gif"))

