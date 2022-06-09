################################ SDM SCRIPT WITH CONTEMPORARY PREDICTION, HINDCASTING & FORECASTING ###############################
## Notes/revisions ###
# 1. Due to legacy reasons, casted variable names are annotated with 'LGM', regardless of whether it refers to LGM (hindcasted) variables or not
# 2. GLM appears inadequate to model allele frequencies, thus consider excluding this method when running for AF.
# 3. When forecasting allele frequencies to future climate, GLM (and to a lesser but still significant extent) GAM appear inadequate to model allele frequencies, thus consider excluding these methods when forecasting for AF.

############################### INITIALISE ENVIRONMENT ###############################
## Load packages
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
#library("maxnet") # loaded later depending on usage
library("PresenceAbsence")
library("randomForest")
library("raster")
library("rasterVis")
library("RColorBrewer")
library("rgdal")
library("rgeos")
library("extendedForest")

### Set working directory and import data
setwd("/Users/luqman/Desktop/SDM/Dsylvestris/Data")

############################### DEFINE TYPE OF ANALYSIS ###############################

# Define species: "Dsyl" or "Dcar"
species <- "Dsyl"
#species <- "Dcar"

# Define whether to hindcast (to LGM) or forecast (to future)
#cast <- "forecast"
cast <- "hindcast"

# Define whether response variable is "binary" (i.e. Presence-Absence) or "frequency" (i.e. continuous between 0 and 1)
responce_variable_type <- "binary"
#responce_variable_type <- "frequency"
### When running "frequency", you may want to replot the rasterVis plots (lines 694 & 797); for some reason, they don't plot automatically when sourced.

# Define whether absence data is presence in input file. If absent (absence_data==FALSE), pseudo-absences will be generated. 
absence_data <- FALSE 
#absence_data <- TRUE
if (responce_variable_type == "frequency") {
  absence_data <- TRUE
}

# Define occurrence data file. This file must contain either 1) presence-only data, 2) presence-absence data or 3) (allele) frequency data. In case of (allele) frequency data, the total number of alleles is required (for weighting); see https://stats.stackexchange.com/questions/233366/how-to-fit-a-mixed-model-with-response-variable-between-0-and-1
# The file must also contain geograpic coordinates (longitude-latitude)
if (species == "Dsyl" && responce_variable_type == "binary") {
  occurrence_file <- "./Occurrence_database_3lineages_wFrenchLongicaulis_revised.csv" # New Dsyl species occurrence data incl. new occurrence data from Italy and the Balkan
  #occurrence_file <- "./Occurrence_database.csv" # Old Dsyl species occurrence data
} else if (species == "Dsyl" && responce_variable_type == "frequency") {
  occurrence_file <- "./CEN_occurrence_revised.csv" # Dsyl CEN data
} else if (species == "Dcar" && responce_variable_type == "binary") {
  occurrence_file <- "/Users/luqman/Desktop/SDM/Dcarthusianorum/Data/Occurrence_database.csv" # New Dsyl species occurrence data incl. new occurrence data from Italy and the Balkan
}

# Define parameters specific for frequency type response data
if (responce_variable_type == "frequency") {
  # Remove variable assignments from previous runs
  #rm(list=setdiff(ls(), c("species", "pred.bin.current.species", "pred.bin.LGM.species", "responce_variable_type", "occurrence_file")))
  # Import species SDM for masking
  pred.bin.current.species <- readRDS("/Users/luqman/Desktop/SDM/Dsylvestris/Data/SDM_CentralAlps_forecast_defaultAbsenceRegime/pred.bin.current.CentralAlps.selectVar.optimalThreshold.2021-04-16_16_32_57.rds")
  if (cast == "forecast") {
    pred.bin.LGM.species <- readRDS("/Users/luqman/Desktop/SDM/Dsylvestris/Data/SDM_CentralAlps_forecast_defaultAbsenceRegime/pred.bin.LGM.CentralAlps.selectVar.optimalThreshold.2021-04-16_16_32_57.rds")
  } else if (cast == "hindcast") {
    pred.bin.LGM.species <- readRDS("/Users/luqman/Desktop/SDM/Dsylvestris/Data/CentralAlps_hindcast_goodExample/pred.bin.LGM.clusters.df.CentralAlps.selectVar.optimalThreshold.2021-05-04_11_43_09.rds")
  }
  # Define the column header of the response variable.
  response_variable <- "FreqHigh"
  # Import a mask for the Central Alpine populations
  CentralAlps_current_mask <- readRDS("CENTRALALPS.CURRENT.RASTER.MASK.210.DISPERSAL0.333.rds")
}

# Define the lineage: "All", "Apennines", "Balkans", "CentralAlps"
if (species == "Dsyl") {
  #lineage <- "All"
  #lineage <- "Apennines"
  #lineage <- "Balkans"
  lineage <- "CentralAlps"
} else if (species == "Dcar") {
  #lineage <- "All"
  lineage <- "CentralAlps"
}

# For Simone's SDM masked_enviro_variable maps, we want to build an SDM mask based on the Central Alps lineage, as that is what his manuscript addresses. Since we don't know the genetic structure of D.carthusianorum, and hence can't do lineage-informed SDMs, we rather constrain the SDM (presence/absence) data and projection to a smaller/contrained extent within the Central Alpine region.
# For you SDM runs, define as FALSE! Only define as TRUE when making Simone's SDM masks.
if (species == "Dcar") {
  if (lineage == "CentralAlps") {
    #constrain_CentralAlps_forSimonesPlots <- FALSE
    constrain_CentralAlps_forSimonesPlots <- TRUE
  } else {
    constrain_CentralAlps_forSimonesPlots <- FALSE
  }
} else if (species == "Dsyl") {
  constrain_CentralAlps_forSimonesPlots <- FALSE
}
 
  
# Define a list of uncorrelated and important predictor variables
if (species == "Dsyl" && responce_variable_type == "binary") {
  # Apennine and Balkan lineage are inferred to have exactly the same set of variables matching our criteria (by GLM D-squared and Pearson's correlation > 0.7). For the central alps, it's very similar, and indeed can be made identical by substituting highly correlated variables (specifically PREC_DRIEST_QUARTER with PREC_COLDEST_QUARTER). The same applies for ALL lineages, substituting PREC_WETTEST_QUARTER with PREC_WARMEST_QUARTER)
  #pred.var.selection <- c("TEMP_MEAN_ANNUAL","TEMP_ANNUAL_RANGE","PREC_ANNUAL","PH_5cm","SLOPE")
  #pred.var.selection <- c("TEMP_MEAN_ANNUAL","ISOTHERMALITY","TEMP_SEASONALITY","PREC_SEASONALITY","PH_5cm","SLOPE") # This may be more appropriate for CentralAlps, as the excluded variables explain only minimal deviance compared to this retained set.
  #pred.var.selection <- c("TEMP_MEAN_ANNUAL","ISOTHERMALITY","TEMP_SEASONALITY","PREC_SEASONALITY","SNOW_COVER_DURATION","PH_15cm","SLOPE") # This may be more appropriate for CentralAlps, as the excluded variables explain only minimal deviance compared to this retained set.
  pred.var.selection <- c("ISOTHERMALITY","TEMP_SEASONALITY","TEMP_MAX_WARMEST_MONTH","TEMP_MEAN_WETTEST_QUARTER","TEMP_MEAN_DRIEST_QUARTER","PREC_SEASONALITY","PREC_WARMEST_QUARTER","PREC_COLDEST_QUARTER","PH_5cm","SLOPE")
  #pred.var.selection <- c("ISOTHERMALITY","TEMP_SEASONALITY","TEMP_MAX_WARMEST_MONTH","TEMP_MEAN_WETTEST_QUARTER","TEMP_MEAN_DRIEST_QUARTER","PREC_SEASONALITY","PREC_WARMEST_QUARTER","PREC_COLDEST_QUARTER","SNOW_COVER_DURATION","ELEV.GMTED","PH_15cm","SLOPE")
} else if (species == "Dsyl" && responce_variable_type == "frequency") {
  # For CEN gene. 
  # For present-day exlusive projection, we can include longitude and latitude to (partially) control for spatial autocorrelation. In this case, this is justified in that there is strong correlation and high r^2 between longitude and latitide and PC1 and PC2 (from the genomic PCA).
  # When hind-casting to the LGM however, we do not expect (LGM) structure to reflect that of the present-day. So we exclude.
  #pred.var.selection <- c("ISOTHERMALITY","TEMP_SEASONALITY","TEMP_MAX_WARMEST_MONTH","TEMP_MEAN_WETTEST_QUARTER","TEMP_MEAN_DRIEST_QUARTER","PREC_SEASONALITY","PREC_WARMEST_QUARTER","PREC_COLDEST_QUARTER","PH_5cm","SLOPE","LONGITUDE","LATITUDE")
  pred.var.selection <- c("ISOTHERMALITY","TEMP_SEASONALITY","TEMP_MAX_WARMEST_MONTH","TEMP_MEAN_WETTEST_QUARTER","TEMP_MEAN_DRIEST_QUARTER","PREC_SEASONALITY","PREC_WARMEST_QUARTER","PREC_COLDEST_QUARTER","PH_5cm","SLOPE")
} else if (species == "Dcar" && responce_variable_type == "binary") {
  # For species.
  #pred.var.selection <- c("TEMP_MEAN_ANNUAL","TEMP_SEASONALITY","PREC_SEASONALITY","PREC_WARMEST_QUARTER","PH_5cm","SLOPE")
  #pred.var.selection <- c("TEMP_MEAN_ANNUAL","TEMP_ANNUAL_RANGE","TEMP_MEAN_WETTEST_QUARTER","TEMP_MEAN_DRIEST_QUARTER","PREC_ANNUAL","PH_5cm","SLOPE")
  pred.var.selection <- c("ISOTHERMALITY","TEMP_SEASONALITY","TEMP_MAX_WARMEST_MONTH","TEMP_MEAN_WETTEST_QUARTER","TEMP_MEAN_DRIEST_QUARTER","PREC_SEASONALITY","PREC_WARMEST_QUARTER","PREC_COLDEST_QUARTER","PH_5cm","SLOPE")
}

# Define projection extent: "europe", "alps"
if (constrain_CentralAlps_forSimonesPlots == FALSE) {
  projection_extent <- "europe"
  #projection_extent <- "alps"
} else if (constrain_CentralAlps_forSimonesPlots == TRUE) {
  projection_extent <- "alps"
}

if (cast == "hindcast") {
  cast.var.name <- "LGM"
} else if (cast == "forecast") {
  cast.var.name <- "RCP45.2061-2080"
}

if (cast == "hindcast") {
  # Define LGM (PMIP3) climate model. Note, TEMP_ANNUAL_RANGE, ISOTHERMALITY, TEMP_SEASONALITY & TEMP_ANNUAL_RANGE are highly variable between the models, ESPECIALLY between CESS-FGOALS-g2, IPSL-CM5A-LR, CNRM-CM5 against the others. The other 4 models are quite consistent with each other, while the aforementioned 3 models are highly variable both within and between these models for these variables.
  # The large variation in these variables among (these models) means that taking the ensemble of all results in weird results (specifically, the SDM is driven by the variation in the TEMP_ANNUAL_RANGE predictor, which is HIGHLY variable between the models (esp 3 models above). Use Check_models.R to visualise these differences. 
  # Thus, taking the ensemble of the more consistent (stable) models is preferred. 
  # Add weighted-ensemble model!
  climate_PMIP3_model_selector <- 2
  clim_PMIP3_model_list.names <- c("ENSEMBLE", "ENSEMBLE_NCAR_MIROC_MRI_MPI", "NCAR-CCSM4", "MIROC-ESM", "MRI-CGCM3", "MPI-ESM-P", "CESS-FGOALS-g2", "IPSL-CM5A-LR", "CNRM-CM5")
  #clim_PMIP3_model_list.data <- c("ENV_RASTERS_LGM_EUROPE_PMIP_ENSEMBLE_ALL.rds", "ENV_RASTERS_LGM_EUROPE_PMIP_ENSEMBLE_NCAR_MIROC_MRI_MPI.rds", "ENV_RASTERS_LGM_EUROPE_PMIP3_NCAR_CCSM4.rds", "ENV_RASTERS_LGM_EUROPE_PMIP_MIROC_ESM.rds", "ENV_RASTERS_LGM_EUROPE_PMIP_MRI-CGCM3.rds", "ENV_RASTERS_LGM_EUROPE_PMIP_MPI_ESM_P.rds", "ENV_RASTERS_LGM_EUROPE_PMIP_CESS_FGOALS_g2.rds", "ENV_RASTERS_LGM_EUROPE_PMIP_IPSL_CM5A_LR.rds", "ENV_RASTERS_LGM_EUROPE_PMIP_CNRM_CM5.rds")
  clim_PMIP3_model_list.data <- c("ENV_RASTERS_LGM_EUROPE_PMIP_ENSEMBLE_ALL_LARGE.grd", "ENV_RASTERS_LGM_EUROPE_PMIP_ENSEMBLE_NCAR_MIROC_MRI_MPI_LARGE.grd", "ENV_RASTERS_LGM_EUROPE_PMIP3_NCAR_CCSM4_LARGE.grd", "ENV_RASTERS_LGM_EUROPE_PMIP3_MIROC_ESM_LARGE.grd", "ENV_RASTERS_LGM_EUROPE_PMIP3_MRI-CGCM3_LARGE.grd", "ENV_RASTERS_LGM_EUROPE_PMIP3_MPI_ESM_P_LARGE.grd", "ENV_RASTERS_LGM_EUROPE_PMIP3_CESS_FGOALS_g2_LARGE.grd", "ENV_RASTERS_LGM_EUROPE_PMIP3_IPSL_CM5A_LR_LARGE.grd", "ENV_RASTERS_LGM_EUROPE_PMIP3_CNRM_CM5_LARGE.grd")
  climate_model_PMIP3 <- clim_PMIP3_model_list.data[climate_PMIP3_model_selector]
  #climate_model_PMIP3 <- "ENV_RASTERS_EUROPE_CHELSA_TraCE_XLARGE_-192.grd"
  #climate_model_PMIP3 <- "ENV_RASTERS_EUROPE_CHELSA_TraCE_XLARGE_2021_identicalTo2020VersionExceptUpdatedDEM-192.grd"
} else if (cast == "forecast") {
  # Define future (CMIP5) climate model. 
  # Here selection follows the recommendation of Knutti & Caldwell, 2015. I.e. here we select >5 models that are the most distant or more specifically whoe the lowest amount of interdependence between each other. In their Fig 4, this is simply selecting models from right to left on the x-axis
  # Here we include the RCP45 scenario, for the time interval 2061-2080.
  climate_CMIP5_model_selector <- 1
  clim_CMIP5_model_list.names <- c("ENSEMBLE", "CESM1_BGC", "MPI_ESM_MR", "MIROC5", "CMCC_CM", "CESM1_CAM5", "IPSL_CM5A_MR", "ACCESS1_0")
  clim_CMIP5_model_list.data <- c("ENV_RASTERS_CMIP5_EUROPE_ENSEMBLE_ALL.rds", "ENV_RASTERS_CMIP5_EUROPE_CESM1-BGC_RCP45_2061-2080.rds", "ENV_RASTERS_CMIP5_EUROPE_MPI-ESM-MR_RCP45_2061-2080.rds", "ENV_RASTERS_CMIP5_EUROPE_MIROC5_RCP45_2061-2080.rds","ENV_RASTERS_CMIP5_EUROPE_CMCC-CM_RCP45_2061-2080.rds", "ENV_RASTERS_CMIP5_EUROPE_CESM1-CAM5_RCP45_2061-2080.rds", "ENV_RASTERS_CMIP5_EUROPE_IPSL-CM5A-MR_RCP45_2061-2080.rds", "ENV_RASTERS_CMIP5_EUROPE_ACCESS1-0_RCP45_2061-2080.rds")
  climate_model_CMIP5 <- clim_CMIP5_model_list.data[climate_CMIP5_model_selector]
}

## Define the CRS (coordinate reference system). We define both longitude-latitude and a UTM (Zone 32) CRS.
prj_longlat <- "+init=epsg:4326"
prj_utm32 <- "+init=epsg:32632"

# Choose which maxent version to run. "maxent" uses (sources) the Maxent Java package (and hence requires maxent installation for external .jar file). As an alternative (same authors, same methodology), you can now use the R package Package ‘maxnet’ (independent R package, with some slightly different implementations).
# While maxent generally gives reliable, consistent results with the other methods, maxnet may give strange, errorneous results (see: https://www.researchgate.net/post/Is_now_a_good_moment_for_an_optimal_transition_to_the_new_version_of_Maxent_and_its_R_package_maxnet)
# Generally we'll want to use maxent, however consider maxnet should you choose to run on the server. Note however, maxnet sometimes encounters an error (see: https://github.com/bobmuscarella/ENMeval/issues/62) which seems to commonly pop-up with Dcar. Hence for Dcar, we use the maxent algorithm instead.
#maxent_maxnet <- "maxnet"
maxent_maxnet <- "maxent" 

# Define number of SDM runs. Given stochastisity in sampling (pseudo-observed) data, different runs will produce different results. To account for this stochasticity, we can run multple iterations of the SDMs (while perturbing initial conditions) to produce an 'average' result.
num_iterations <- 1
#num_iterations <- 10 # Final settings

# Choose whether to plot intermediate results
if (num_iterations <= 1) {
  plot_all <- TRUE
} else if (num_iterations > 1) {
  plot_all <- FALSE
}

############################### MAIN SDM FUNCTION ###############################

# Initialise empty raster stack and list to store results
pred.prob.LGM.stack <- stack()
if (responce_variable_type == "binary") {
  thres.stack <- c()
  auc.stack <- c()
}

for (run in seq(1,num_iterations)) {
  
  # Update progress of loop
  print(paste("Generating SDM iteration:", run, sep=" "))
  
  ############################### IMPORT AND PREPARE THE OCCURRENCE DATA ###############################
  ## Load the occurrence data:
  Occurrences <- read.csv(occurrence_file, h=T)
  
  if (responce_variable_type == "frequency") {
    # Filter occurrence data to remove outliers (as identified in the model diagnoses below)
    Occurrences <- Occurrences[-c(65), ]
    row.names(Occurrences) <- seq(1,nrow(Occurrences))
  }
  
  # Import stack of environmental rasters
  ENV_RASTERS <-stack("ENV_RASTERS_TODAY_EUROPE_CHELSA_LARGE.grd") 
  #ENV_RASTERS <-stack("ENV_RASTERS_EUROPE_CHELSA_TraCE_XLARGE_20.grd")
  #ENV_RASTERS <-stack("ENV_RASTERS_EUROPE_CHELSA_TraCE_XLARGE_2021_identicalTo2020VersionExceptUpdatedDEM20.grd")
  
  ## Prepare the data
  if (responce_variable_type == "binary") {
    
    # Split the data into according to species.
    if (species == "Dsyl") {
      if (lineage == "All") {
        # We first need to balance the data, as we have many more occurrences in the Central Alps compared to the Apennines and Balkans
        Occurrences.CBNA <- Occurrences[Occurrences$Source == "CBNA",] 
        Occurrences.others <- Occurrences[Occurrences$Source != "CBNA",] 
        Occurrences.CBNA <- Occurrences.CBNA[sample(nrow(Occurrences.CBNA), length(Occurrences.others[,1])/5), ]
        Occurrences <- rbind(Occurrences.others, Occurrences.CBNA)
        Occurrences.CentralAlps <- Occurrences[Occurrences$Species == "Dianthus sylvestris ssp sylvestris",] 
        Occurrences.otherLineages <- Occurrences[Occurrences$Species != "Dianthus sylvestris ssp sylvestris",] 
        Occurrences.CentralAlps <- Occurrences.CentralAlps[sample(nrow(Occurrences.CentralAlps), length(Occurrences.otherLineages[,1])/1.5), ]
        Occurrences <- rbind(Occurrences.otherLineages, Occurrences.CentralAlps)
        rownames(Occurrences) <- NULL
      } else if (lineage == "CentralAlps") {
        Occurrences <- Occurrences[Occurrences$Species == "Dianthus sylvestris ssp sylvestris",]
        # If necessary (e.g. with the CBNA and CBNMed databases), balance the data
        Occurrences.CBNA <- Occurrences[Occurrences$Source == "CBNA",] 
        Occurrences.others <- Occurrences[Occurrences$Source != "CBNA",] 
        Occurrences.CBNA <- Occurrences.CBNA[sample(nrow(Occurrences.CBNA), length(Occurrences.others[,1])/5), ]
        Occurrences <- rbind(Occurrences.others, Occurrences.CBNA)
        if (constrain_CentralAlps_forSimonesPlots == TRUE) {
          Occurrences<-Occurrences[Occurrences$x <= 11,]
          Occurrences<-Occurrences[Occurrences$x >= 6,]
          Occurrences<-Occurrences[Occurrences$y <= 48,]
          Occurrences<-Occurrences[Occurrences$y >= 45,]
        }
        rownames(Occurrences) <- NULL
      } else if (lineage == "Apennines") {
        Occurrences<-Occurrences[Occurrences$Species == "Dianthus sylvestris ssp longicaulis",]
        Occurrences.CBNA <- Occurrences[Occurrences$Source == "CBNA",] 
        Occurrences.others <- Occurrences[Occurrences$Source != "CBNA",] 
        Occurrences.CBNA <- Occurrences.CBNA[sample(nrow(Occurrences.CBNA), length(Occurrences.others[,1])/5), ]
        Occurrences <- rbind(Occurrences.others, Occurrences.CBNA)
      } else if (lineage == "Balkans") {
        Occurrences<-Occurrences[Occurrences$Species == "Dianthus sylvestris ssp agg.",]
        #Occurrences<-Occurrences[Occurrences$y <= 45.8,] # should we want to exclude the Julian Alps cluster/populations. We do this because I believe they have a distinct niche compared to the rest of the Balkan sub-lineages (with a niche closer to the Central Alp populations), and because they also occupy a distinct grouping within the Balkan lineage.
      }
    } else if (species == "Dcar") {
      if (lineage == "All") {
        Occurrences<-Occurrences[Occurrences$x <= 23,]
        Occurrences<-Occurrences[Occurrences$x >= 3,]
        Occurrences<-Occurrences[Occurrences$y <= 48.5,]
        Occurrences<-Occurrences[Occurrences$y >= 37.5,]
      } else if (lineage == "CentralAlps" && constrain_CentralAlps_forSimonesPlots == FALSE) {
        Occurrences<-Occurrences[Occurrences$x <= 15,]
        Occurrences<-Occurrences[Occurrences$x >= 4,]
        Occurrences<-Occurrences[Occurrences$y <= 48,]
        Occurrences<-Occurrences[Occurrences$y >= 43,]
      } else if (lineage == "CentralAlps" && constrain_CentralAlps_forSimonesPlots == TRUE) {
        Occurrences<-Occurrences[Occurrences$x <= 11,]
        Occurrences<-Occurrences[Occurrences$x >= 6,]
        Occurrences<-Occurrences[Occurrences$y <= 48,]
        Occurrences<-Occurrences[Occurrences$y >= 45,]
      }
    }

    # Revise dataframe by excluding unnecessary columns and retaining only x-y coordinates and presence/absence. Remove duplicate entries.
    Occurrences<-Occurrences[,c(3,4,8)]
    Occurrences<-unique(Occurrences)
    
  } else if (responce_variable_type == "frequency") {
    Occurrences<-Occurrences[,c(4,3,2,6,7,8,9,10,11,12)]
    Occurrences<-unique(Occurrences)
    names(Occurrences)[names(Occurrences) == response_variable] <- "Presence"
    Occurrences.desagg <- Occurrences
  }
  
  ############################### Generating pseudo-absence data ###############################
  if (absence_data == FALSE) {
    ## Since we have only presence data, we'll need to generate some pseudo-absence data. We do this by randomly sample coordinates within the geographic extent we think they can inhabit.
    ## First, we'll need a raster map to draw random samples from. Here we'll use a 30 arcsec DEM raster
    #ELEV.GMTED.init <- raster("./GMTED2010N30E000_300/30n000e_20101117_gmted_mea300.tif")
    ELEV.GMTED.init <- ENV_RASTERS$ELEV.GMTED
    ## Crop the map to the extent of interest
    if (lineage == "All") {
      e <- extent(3, 23, 37.5, 48.5)
      e_absence <- extent(4.7, 21.8, 37.7, 47.8)
      e_absence_cropping_polygon_x_coords_1 <- c(7, 10, 10, 7, 7)
      e_absence_cropping_polygon_y_coords_1 <- c(43.5, 43.5, 38.5, 38.5, 43.5)
      e_absence_cropping_polygon_xy_matrix_1 <- cbind(e_absence_cropping_polygon_x_coords_1,e_absence_cropping_polygon_y_coords_1)
      e_absence_cropping_polygon_xy_matrix_list <- list("Poly1" = e_absence_cropping_polygon_xy_matrix_1)
    } else if (lineage == "CentralAlps" && species == "Dsyl" && constrain_CentralAlps_forSimonesPlots == FALSE) {
      e <- extent(4, 15, 43, 48)
      #e <- extent(3, 17.5, 42, 48.5)
      e_absence <- extent(4.7, 13.4, 43.6, 47.6)
      #e_absence <- extent(3.3, 17.2, 42.2, 48.3)
      e_absence_cropping_polygon_x_coords_1 <- c(4.7, 13.4, 13.4, 4.7, 4.7) # this is redundant but avoids adding unnecessary code later
      e_absence_cropping_polygon_y_coords_1 <- c(47.6, 47.6, 43.6, 43.6, 47.6) # this is redundant but avoids adding unnecessary code later
      e_absence_cropping_polygon_x_coords_2 <- c(4.7, 9.9, 9.9, 4.7, 4.7)
      e_absence_cropping_polygon_y_coords_2 <- c(43.1, 43.1, 38.8, 38.8, 43.1)
      e_absence_cropping_polygon_xy_matrix_1 <- cbind(e_absence_cropping_polygon_x_coords_1,e_absence_cropping_polygon_y_coords_1)
      e_absence_cropping_polygon_xy_matrix_2 <- cbind(e_absence_cropping_polygon_x_coords_2,e_absence_cropping_polygon_y_coords_2)
      e_absence_cropping_polygon_xy_matrix_list <- list("Poly1" = e_absence_cropping_polygon_xy_matrix_1, "Poly2" = e_absence_cropping_polygon_xy_matrix_2)
    } else if (lineage == "CentralAlps" && species == "Dsyl" && constrain_CentralAlps_forSimonesPlots == TRUE) {
      e <- extent(6, 11, 45, 48)
      e_absence <- extent(6.1, 10.9, 45.1, 47.9)
      e_absence_cropping_polygon_x_coords_1 <- c(6.1, 10.9, 10.9, 6.1, 6.1) # this is redundant but avoids adding unnecessary code later
      e_absence_cropping_polygon_y_coords_1 <- c(47.9, 47.9, 45.1, 45.1, 47.9) # this is redundant but avoids adding unnecessary code later
      e_absence_cropping_polygon_xy_matrix_1 <- cbind(e_absence_cropping_polygon_x_coords_1,e_absence_cropping_polygon_y_coords_1)
      e_absence_cropping_polygon_xy_matrix_list <- list("Poly1" = e_absence_cropping_polygon_xy_matrix_1)
    } else if (lineage == "CentralAlps" && species == "Dcar" && constrain_CentralAlps_forSimonesPlots == FALSE) {
      e <- extent(4, 15, 43, 48)
      #e <- extent(3, 21, 38, 48)
      e_absence <- extent(4.7, 13.4, 43.6, 47.6)
      e_absence_cropping_polygon_x_coords_1 <- c(4.7, 13.4, 13.4, 4.7, 4.7) # this is redundant but avoids adding unnecessary code later
      e_absence_cropping_polygon_y_coords_1 <- c(47.6, 47.6, 43.6, 43.6, 47.6) # this is redundant but avoids adding unnecessary code later
      e_absence_cropping_polygon_xy_matrix_1 <- cbind(e_absence_cropping_polygon_x_coords_1,e_absence_cropping_polygon_y_coords_1)
      e_absence_cropping_polygon_xy_matrix_list <- list("Poly1" = e_absence_cropping_polygon_xy_matrix_1)
    } else if (lineage == "CentralAlps" && species == "Dcar" && constrain_CentralAlps_forSimonesPlots == TRUE) {
      e <- extent(6, 11, 45, 48)
      e_absence <- extent(6.1, 10.9, 45.1, 47.9)
      e_absence_cropping_polygon_x_coords_1 <- c(6.1, 10.9, 10.9, 6.1, 6.1) # this is redundant but avoids adding unnecessary code later
      e_absence_cropping_polygon_y_coords_1 <- c(47.9, 47.9, 45.1, 45.1, 47.9) # this is redundant but avoids adding unnecessary code later
      e_absence_cropping_polygon_xy_matrix_1 <- cbind(e_absence_cropping_polygon_x_coords_1,e_absence_cropping_polygon_y_coords_1)
      e_absence_cropping_polygon_xy_matrix_list <- list("Poly1" = e_absence_cropping_polygon_xy_matrix_1)
    } else if (lineage == "Apennines") {
      e <- extent(4, 18.8, 37.6, 47.5)
      e_absence <- extent(4.7, 18.7, 37.9, 47)
      #e_absence_cropping_polygon_x_coords_1 <- c(7.8, 13.3, 13.3, 14.6, 14.6, 17.8, 17.8, 18.7, 18.7, 15.5, 15.5, 10, 10, 7.8, 7.8) # polygon of Apennine peninsula
      e_absence_cropping_polygon_x_coords_1 <- c(9.5, 21.8, 21.8, 18.6, 18.6, 14.9, 14.9, 13, 13, 9.5, 9.5)
      #e_absence_cropping_polygon_y_coords_1 <- c(45, 45, 44, 44, 42.5, 42.5, 40.8, 40.8, 37.7, 37.7, 39.8, 39.8, 43.7, 43.7, 45) # polygon of Apennine peninsula
      e_absence_cropping_polygon_y_coords_1 <- c(47, 47.3, 39.3, 39.3, 42.2, 42.2, 44, 44, 44.4, 44.4, 47)
      e_absence_cropping_polygon_x_coords_2 <- c(4.7, 9.9, 9.9, 15.5, 15.5, 4.7, 4.7)
      e_absence_cropping_polygon_y_coords_2 <- c(43.1, 43.1, 38.8, 38.8, 37.9, 37.9, 43.1)
      e_absence_cropping_polygon_xy_matrix_1 <- cbind(e_absence_cropping_polygon_x_coords_1,e_absence_cropping_polygon_y_coords_1)
      e_absence_cropping_polygon_xy_matrix_2 <- cbind(e_absence_cropping_polygon_x_coords_2,e_absence_cropping_polygon_y_coords_2)
      e_absence_cropping_polygon_xy_matrix_list <- list("Poly1" = e_absence_cropping_polygon_xy_matrix_1, "Poly2" = e_absence_cropping_polygon_xy_matrix_2)
    } else if (lineage == "Balkans") {
      #e <- extent(11.8, 22.2, 39, 47.9)
      e <- extent(9.8, 23, 37.8, 48.2)
      #e_absence <- extent(12, 21.8, 39.3, 47.8)
      e_absence <- extent(10, 23, 38, 48)
      #e_absence_cropping_polygon_x_coords_1 <- c(12, 21.8, 21.8, 18.6, 18.6, 14.9, 14.9, 13, 13, 12, 12) # pseudo-absence range restricted to presence range
      e_absence_cropping_polygon_x_coords_1 <- c(9.5, 23, 23, 18.6, 18.6, 14.9, 14.9, 13, 13, 9.5, 9.5) # expanded pseudo-absence range
      #e_absence_cropping_polygon_y_coords_1 <- c(47.8, 47.8, 39.3, 39.3, 42.2, 42.2, 44, 44, 44.4, 44.4, 47.8) # pseudo-absence range restricted to presence range
      e_absence_cropping_polygon_y_coords_1 <- c(48, 48, 38, 38, 42.2, 42.2, 44, 44, 45.5, 45.5, 48) # expanded pseudo-absence range
      e_absence_cropping_polygon_x_coords_2 <- c(17.8, 23, 23, 21.3, 21.3, 19.6, 19.6, 17.8, 17.8)
      e_absence_cropping_polygon_y_coords_2 <- c(48, 48, 45.1, 45.1, 46, 46, 47, 47, 48)
      e_absence_cropping_polygon_xy_matrix_1 <- cbind(e_absence_cropping_polygon_x_coords_1,e_absence_cropping_polygon_y_coords_1)
      e_absence_cropping_polygon_xy_matrix_2 <- cbind(e_absence_cropping_polygon_x_coords_2,e_absence_cropping_polygon_y_coords_2)
      #e_absence_cropping_polygon_xy_matrix_list <- list("Poly1" = e_absence_cropping_polygon_xy_matrix_1)
      e_absence_cropping_polygon_xy_matrix_list <- list("Poly1" = e_absence_cropping_polygon_xy_matrix_1, "Poly2" = e_absence_cropping_polygon_xy_matrix_2)
    }
    ELEV.GMTED.alps.small <- crop(ELEV.GMTED.init, e_absence)
    ELEV.GMTED.alps <- crop(ELEV.GMTED.init, e)
    ## We need to limit the random sampling of coordinates to habitable zones (e.g. excluding the sea). We can achieve this by masking all values below e.g. 1m as "NA" using the mask funtion
    ELEV.GMTED.alps.small[ELEV.GMTED.alps.small < 50] <- NA
    ## We further want to contrain the pseudo-absences to lie within certain areas
    e_absence_cropping_polygon <- lapply(e_absence_cropping_polygon_xy_matrix_list, Polygon)
    e_absence_cropping_polygons <- lapply(seq_along(e_absence_cropping_polygon), function(i) Polygons(list(e_absence_cropping_polygon[[i]]), ID = names(e_absence_cropping_polygon)[i]))
    e_absence_cropping_SP <- SpatialPolygons(e_absence_cropping_polygons, proj4string=CRS(prj_longlat))
    if (lineage == "All" || lineage == "Apennines") {
      ELEV.GMTED.alps.small <- raster::mask(ELEV.GMTED.alps.small, e_absence_cropping_SP[1,], inverse=TRUE)
    } else {
      ELEV.GMTED.alps.small <- raster::mask(ELEV.GMTED.alps.small, e_absence_cropping_SP[1,])
    }
    if (species == "Dsyl" && constrain_CentralAlps_forSimonesPlots == FALSE) {
      if (lineage == "Balkans" || lineage == "Apennines" || lineage == "CentralAlps") {
        ELEV.GMTED.alps.small <- raster::mask(ELEV.GMTED.alps.small, e_absence_cropping_SP[2,], inverse=TRUE)
      }
    }
    
    ## To randomly sample coordinates, we use the sampleRandom function in the Raster package.
    ## Initially, we'll sample more absences than presences (here 5x). Finally, we'll want to have approximately equal numner of absences to presences, but we'll remove some later (with the pseudoabsence function)
    num_presences<-as.numeric(length(Occurrences[,1]))
    ## We then randomly sample (without replacement) the same number of coordinates as presences.
    if (constrain_CentralAlps_forSimonesPlots == FALSE && cast == "hindcast") {
      absences<-as.data.frame(sampleRandom(ELEV.GMTED.alps.small, size=200*num_presences, na.rm = TRUE, cells=FALSE, sp=FALSE, xy=TRUE))
    } else if (constrain_CentralAlps_forSimonesPlots == TRUE || cast == "forecast") {
      absences<-as.data.frame(sampleRandom(ELEV.GMTED.alps.small, size=10*num_presences, na.rm = TRUE, cells=FALSE, sp=FALSE, xy=TRUE))
    }
    ######## POTENTIAL ERROR: in some case, an error will pop up here. Either try running again or change/reduce the size argument.
    ## Note that this random sampling samples much more from the lowlands than the mountain peaks, given that there is much more land area in the lowlands than near the mountain peaks. 
    ## This will lead however to pseudo-absence points more often that not falling in the lowlands rather than the mountain peaks, which leads to an unbounded (lower tempearture side) in the response curves, and biased representation in LGM projection (where high elevation glacial ice sheets predominate). 
    #hist(absences$ELEV.GMTED)
    ## Want way to solve this would be to random sample not simply by area, but also by elevation, e.g. give equal sampling across the distribution of elevation. Because in the end, what we're interested in is sampling across environmental space moreso than geographic space (and especially past/future space, proxied by elevation, if we plan to hind/forecast)
    if (constrain_CentralAlps_forSimonesPlots == FALSE && cast == "hindcast") {
      if (lineage == "CentralAlps"  && species == "Dsyl") {
        absences_0_1000 <- absences[match(sample(absences[absences$ELEV.GMTED < 1000,][,3], num_presences*1.25), absences$ELEV.GMTED),]
        absences_1000_2000 <- absences[match(sample(absences[(absences$ELEV.GMTED > 1000) & (absences$ELEV.GMTED < 2000), ][,3], num_presences*1.25), absences$ELEV.GMTED),]
        absences_2000_3000 <- absences[match(sample(absences[(absences$ELEV.GMTED > 2000) & (absences$ELEV.GMTED < 3000), ][,3], num_presences), absences$ELEV.GMTED),]
        absences_3000_above <- absences[match(sample(absences[(absences$ELEV.GMTED > 3000), ][,3], num_presences), absences$ELEV.GMTED),]
        absences <- rbind(absences_0_1000, absences_1000_2000, absences_2000_3000, absences_3000_above)
      } else if (lineage == "CentralAlps"  && species == "Dcar") {
        absences_0_1000 <- absences[match(sample(absences[absences$ELEV.GMTED < 1000,][,3], num_presences), absences$ELEV.GMTED),]
        absences_1000_2000 <- absences[match(sample(absences[(absences$ELEV.GMTED > 1000) & (absences$ELEV.GMTED < 2000), ][,3], num_presences), absences$ELEV.GMTED),]
        absences_2000_3000 <- absences[match(sample(absences[(absences$ELEV.GMTED > 2000) & (absences$ELEV.GMTED < 3000), ][,3], num_presences*0.75), absences$ELEV.GMTED),]
        absences_3000_above <- absences[match(sample(absences[(absences$ELEV.GMTED > 3000), ][,3], num_presences*0.75), absences$ELEV.GMTED),]
        absences <- rbind(absences_0_1000, absences_1000_2000, absences_2000_3000, absences_3000_above)
      } else {
        absences_0_1000 <- absences[match(sample(absences[absences$ELEV.GMTED < 1000,][,3], num_presences*1.5), absences$ELEV.GMTED),]
        absences_1000_2000 <- absences[match(sample(absences[(absences$ELEV.GMTED > 1000) & (absences$ELEV.GMTED < 2000), ][,3], num_presences*1.5), absences$ELEV.GMTED),]
        absences_2000_above <- absences[match(sample(absences[(absences$ELEV.GMTED > 2000), ][,3], num_presences), absences$ELEV.GMTED),]
        absences <- rbind(absences_0_1000, absences_1000_2000, absences_2000_above)
      }
    }
  
    absences<-absences[,c(1,2,3)]
    colnames(absences) <- c("x", "y", "Presence")
    absences$Presence [absences$Presence > 0] <- 0
    
    ## We desaggregate the occurrence (and absence) data, to avoid biases in imbalanced sampling. Recall, 1 degree is approx 111,139m. Let's take a minimum distance of 5.5km (0.05 degrees) for the pooled species and 3.33 for the separated lineages (smaller for the separated lineages because we don't have so many occurrence points and we want to coserve them). Remember, it does NOT make so much sense to take a minimum distance less than the resolution of the predictor variables (<900m).
    if (species == "Dsyl") {
      if (lineage == "All" || lineage == "CentralAlps") {
        min_desagg_dist <- 0.05
      } else {
        min_desagg_dist <- 0.03
      }
    } else if (species == "Dcar") {
      min_desagg_dist <- 0.08
    }
    if (plot_all == TRUE) {
      Occurrences.desagg.temp <- ecospat.occ.desaggregation(xy=Occurrences, min.dist=min_desagg_dist, by=NULL)
    } else if (plot_all == FALSE) {
      body(ecospat.occ.desaggregation)[[8]] <- substitute(result <- result) # this line edits the ecospat.occ.desaggregation function to withhold the printing of results
      Occurrences.desagg.temp <- ecospat.occ.desaggregation(xy=Occurrences, min.dist=min_desagg_dist, by=NULL)
    }
    num_presences.desagg<-as.numeric(length(Occurrences.desagg.temp[,1]))
    
    # Now we will make a selection of the pseudo-absence points. We do this conditional on a minimum distance to a presence point (below which we cannot draw an absence point). For this we will use the 'ecospat.rand.pseudoabsences' function, which randomly subsamples a absence/pseudoabsence list (doesn't upsample/interpolate/extrapolate) given a set of arguments.
    absence.presence.desagg <- as.data.frame(t(cbind(t(Occurrences.desagg.temp), t(absences))))
    presences.desagg <- absence.presence.desagg[absence.presence.desagg[,3]!=0,1:2]
    # We return slightly more than num_presences.desagg (e.g. 1.3x) because some coordinates may by sampled twice (duplicates)
    min_absence_desagg_dist <- 0.02
    pseudo.absences <- ecospat.rand.pseudoabsences (nbabsences=1.3*num_presences.desagg, glob=absence.presence.desagg , colxyglob=1:2, colvar = NULL, presence=presences.desagg, colxypresence=1:2, mindist=min_absence_desagg_dist)
    pseudo.absences$Presence <- 0
    pseudo.absences <- unique(pseudo.absences)
    
    # Plot the presence and pseudoabsence data
    if (plot_all == TRUE) {
      par(mfrow=c(1,2), mar = c(1, 1, 1, 1))
      plot(ELEV.GMTED.alps, main="Presences")
      points(Occurrences.desagg.temp)
      plot(ELEV.GMTED.alps, main="Pseudo-absences")
      points(pseudo.absences)
    }
    
    # We then append the absence data to the presence data, remove any duplicate entries/coordinates if there exists any, rename the rows
    Occurrences.desagg <- as.data.frame(t(cbind(t(Occurrences.desagg.temp), t(pseudo.absences))))
    Occurrences.desagg <- Occurrences.desagg[!duplicated(Occurrences.desagg[,1:2]),]
    row.names(Occurrences.desagg) <- 1:nrow(Occurrences.desagg)
    
    # Remove temporary files
    rm(ELEV.GMTED.init, ELEV.GMTED.alps, absence.presence.desagg, Occurrences.desagg.temp, e)
    rm(e_absence, ELEV.GMTED.alps.small, Occurrences.CBNA, Occurrences.others, Occurrences.CentralAlps, Occurrences.otherLineages)
  }
  
  ############################### PREPARE THE ENVIRONMENTAL DATA ###############################
  
  # Extract variable raster stack and import spatial polygons layer for country borders
  if (projection_extent == "alps") {
    if (lineage == "All") {
      e_alps <- extent(3, 23, 37.5, 48.5)
    } else if (lineage == "CentralAlps" && constrain_CentralAlps_forSimonesPlots == FALSE) {
      e_alps <- extent(4, 15, 43, 48)
    } else if (lineage == "CentralAlps" && constrain_CentralAlps_forSimonesPlots == TRUE) {
      e_alps <- extent(6, 11, 45, 48)
    } else if (lineage == "Apennines") {
      e_alps <- extent(7.7, 18.8, 37.6, 45.2)
    } else if (lineage == "Balkans") {
      e_alps <- extent(11.8, 22.2, 39, 47.9)
    }
    ENV_RASTERS <- crop(ENV_RASTERS, e_alps)
    borders <- readRDS("alps.rds")
  } else if (projection_extent == "europe") {
    if (cast == "forecast" || responce_variable_type == "frequency") {
      e_europe <- extent(3, 21, 38, 48)
    } else {
      e_europe <- extent(3, 23, 37.5, 48.5)
    }
    ENV_RASTERS <- crop(ENV_RASTERS, e_europe)
    borders <- readRDS("europe.rds")
  }
  #plot(ENV_RASTERS)
  
  # Extract values from raster stack, and convert to dataframe. Don't forget to convert non-numeric data types to factors
  Climate <- as.data.frame(extract(ENV_RASTERS, Occurrences.desagg[,c("x","y")]))
  Climate[,26] <- as.factor(Climate[,26])
  Climate[,27] <- as.factor(Climate[,27])
  
  ## Append climatic data to occurrences data and discard incomplete observations (make sure to check your other fields/columns don't have N/As, e.g. in elevation, source etc.). This will represent the main DATASET for this SDM analysis.
  Occurrences.Climate <- cbind(Occurrences.desagg, Climate)
  Occurrences.Climate <- na.omit(Occurrences.Climate)
  ## Check that the variables are coded in the right datatype, e.g. num for continous variables and (unordered factor) for categorical variables.
  #str(Occurrences.Climate)
  
  ############################### SELECT THE PREDICTOR (ENVIRONMENTAL) VARIABLES ###############################
  
  ## We check the predictors in our data set, which we extracted from the grid-stack. We check for correlation between variables; we first load the pairs.hist.cor function. We can exclude the categorical (unordered factor) variables here.
  if (plot_all == TRUE) {
    # Open a pdf file
    #pdf(paste0("/Users/luqman/Desktop/Corr_",lineage,"plot.pdf"), width = 30, height = 30)
    load("pairs.hist.cor.rda")
    if (responce_variable_type == "binary") {
      pairs.hist.cor(Occurrences.Climate[, 4:28], cor.method = "pearson")
    } else if (responce_variable_type == "frequency") { 
      pairs.hist.cor(Occurrences.Climate[, 11:30], cor.method = "pearson")
    }
    # Close the pdf file
    #dev.off()
    # Visualise correlation between variables as a hierarchical cluster tree
    if (responce_variable_type == "binary") {
      round(cor(Occurrences.Climate[, 4:28]),2)
      varCor <- cor(Occurrences.Climate[, 4:28])
      allDistNew <- abs(as.dist(cor(Occurrences.Climate[, 4:28])))
    } else if (responce_variable_type == "frequency") { 
      round(cor(Occurrences.Climate[, 11:30]),2)
      varCor <- cor(Occurrences.Climate[, 11:30])
      allDistNew <- abs(as.dist(cor(Occurrences.Climate[, 11:30])))
    }
    allClusNew <- hclust(1 - allDistNew)
    par(mar = c(1, 1, 1, 1), mfrow = c(1, 1)) 
    plot(allClusNew, hang=-1)
  }

  ## Next, we write a simple function to test the predictive power of each variable alone (utilising GLM). Multiply by 100 to get percentage rather than fraction.
  if (responce_variable_type == "binary") {
    cont.pred.var <- sapply(names(Occurrences.Climate)[4:28], function(x) ecospat.adj.D2.glm(glm(Occurrences.Climate$Presence ~ Occurrences.Climate[, x]+I(Occurrences.Climate[, x]^2), family = "binomial")))
    cat.pred.var <- sapply(names(Occurrences.Climate)[29:30], function(x) ecospat.adj.D2.glm(glm(Occurrences.Climate$Presence ~ Occurrences.Climate[, x], family = "binomial")))
    int.pred.var <- ecospat.adj.D2.glm(glm(Occurrences.Climate$Presence ~ Occurrences.Climate$SLOPE * Occurrences.Climate$ASPECT, family = "binomial", maxit = 100))
  } else if (responce_variable_type == "frequency") { cont.pred.var <- sapply(names(Occurrences.Climate)[11:35], function(x) ecospat.adj.D2.glm(glm(Occurrences.Climate$Presence ~ Occurrences.Climate[, x]+I(Occurrences.Climate[, x]^2), family = "binomial", weights = Occurrences.Climate$Nr_Total)))
  cat.pred.var <- sapply(names(Occurrences.Climate)[36:37], function(x) ecospat.adj.D2.glm(glm(Occurrences.Climate$Presence ~ Occurrences.Climate[, x], family = "binomial", weights = Occurrences.Climate$Nr_Total)))
  int.pred.var <- ecospat.adj.D2.glm(glm(Occurrences.Climate$Presence ~ Occurrences.Climate$SLOPE * Occurrences.Climate$ASPECT, family = "binomial", maxit = 100, weights = Occurrences.Climate$Nr_Total))
  }
  if (plot_all == TRUE) {
    cont.pred.var.df <- as.data.frame(cont.pred.var)
    cat.pred.var.df <- as.data.frame(cat.pred.var)
    colnames(cont.pred.var.df) <- "D2"
    colnames(cat.pred.var.df) <- "D2"
    GLM_uniVar_deviance <- rbind(cont.pred.var.df,cat.pred.var.df)
    print(GLM_uniVar_deviance)
    print(paste("SLOPE x ASPECT : ", round(int.pred.var, 7), sep = " "))
  }
  
  ## We can also test the importance of each predictor e.g. in a random forest model:
  rforest.allVars.form <- as.formula(paste("Presence", paste(names(ENV_RASTERS), collapse = " + "), sep = " ~ "))
  rf.prob.allVars <- randomForest(rforest.allVars.form, data=Occurrences.Climate)
  if (plot_all == TRUE) {
    rf.prob.allVars.importance <- importance(rf.prob.allVars)
    print(rf.prob.allVars.importance)
    randomForest::varImpPlot(rf.prob.allVars)
  }
  
  ## Based on the correlation and predictive power of the predictor variables, we make a preliminary selection of variables.
  Climate <- Climate[,pred.var.selection]
  
  ## Append climatic data to occurrences data and discard incomplete observations (make sure to check your other fields/columns don't have N/As, e.g. in elevation, source etc.). This will represent the main DATASET for this SDM analysis.
  Occurrences.Climate <- cbind(Occurrences.desagg, Climate)
  Occurrences.Climate <- na.omit(Occurrences.Climate)
  
  # In case you run against a single predictor variable
  #colnames(Occurrences.Climate)[4] <- "TEMP_MEAN_ANNUAL"
  #colnames(Occurrences.Climate)[4] <- "PREC_SEASONALITY"
  
  ## Let's test the predictive power of each variable alone, again.
  if (responce_variable_type == "binary") {
    first.pred.var.idx <- 4
  } else if (responce_variable_type == "frequency") {
    #first.pred.var.idx <- match("TEMP_MEAN_ANNUAL", colnames(Occurrences.Climate))
    first.pred.var.idx <- match(pred.var.selection[1], colnames(Occurrences.Climate))
  }
  last.pred.var.idx <- length(Occurrences.Climate)
  if (responce_variable_type == "binary") {
    sapply(names(Occurrences.Climate)[first.pred.var.idx:last.pred.var.idx], function(x) ecospat.adj.D2.glm(glm(Occurrences.Climate$Presence ~ Occurrences.Climate[, x]+I(Occurrences.Climate[, x]^2), family = "binomial"))) * 100
  } else if (responce_variable_type == "frequency"){ sapply(names(Occurrences.Climate)[first.pred.var.idx:last.pred.var.idx], function(x) ecospat.adj.D2.glm(glm(Occurrences.Climate$Presence ~ Occurrences.Climate[, x]+I(Occurrences.Climate[, x]^2), family = "binomial", weights = Occurrences.Climate$Nr_Total))) * 100
  }
  
  ## To allow GLM and GAM to work with allele frequencies, we define the response variable to be the allele frequencies AND we weight the frequencies by the total number(occurrence) of alleles.
  # Recall: if the response variable is a fraction p=m/n of two integers and all ns are known, then one can use standard logistic regression, aka binomial GLM, given that the weights are supplied as such:
  # glm(p ~ a+b+c, myData, family="binomial", weights=n), assuming that n is a vector of N values for each data point. See: https://stats.stackexchange.com/questions/233366/how-to-fit-a-mixed-model-with-response-variable-between-0-and-1
  if (responce_variable_type == "binary") {
    GLM.weights <- NULL
  } else if (responce_variable_type == "frequency") {
    GLM.weights <- Occurrences.Climate$Nr_Total # This is the header for the total number of haplotypes in the data table
  }
  
  ############################### Split the data into training and evaluation datasets, for later cross-validation and model evaluation (for Random Forest and Maxent) ###############################
  
  # Let's make training and testing datasets, for model evaluation
  if (responce_variable_type == "binary") {
    pseudo.absences.temp <- pseudo.absences[,c(1,2)]
    row.names(pseudo.absences.temp) <- 1:nrow(pseudo.absences.temp)
    pseudo.absences.group <- kfold(pseudo.absences.temp, 5)
    row.names(presences.desagg) <- 1:nrow(presences.desagg)
    presences.group <- kfold(presences.desagg, 5)
    presence_train <- presences.desagg[presences.group != 1, ]
    presence_test <- presences.desagg[presences.group == 1, ]
    absence_train <- pseudo.absences.temp[presences.group != 1, ]
    absence_test <- pseudo.absences.temp[presences.group == 1, ]
    presence_train$Presence <- 1
    absence_train$Presence <- 0
  } #else if (responce_variable_type == "frequency") {
    # presence.col.index <- match("Presence", colnames(Occurrences.desagg))
    # absences.temp <- Occurrences.desagg[Occurrences.desagg$Presence == 0, ]
    # absences.temp <- absences.temp[,c(1,2,presence.col.index)]
    # row.names(absences.temp) <- 1:nrow(absences.temp)
    # absences.group <- kfold(absences.temp, 5)
    # presences.temp <- Occurrences.desagg[Occurrences.desagg$Presence > 0, ]
    # presences.temp <- presences.temp[,c(1,2,presence.col.index)]
    # row.names(presences.temp) <- 1:nrow(presences.temp)
    # presences.group <- kfold(presences.temp, 5)
    # presence_train <- presences.temp[presences.group != 1, ]
    # presence_test <- presences.temp[presences.group == 1, ]
    # absence_train <- absences.temp[presences.group != 1, ]
    # absence_test <- absences.temp[presences.group == 1, ]
  #}
  
  if (responce_variable_type == "binary") {
    # Extract predictor variables from the ENV raster for the testing datasets
    CLIMATE.presence.test <- as.data.frame( extract(ENV_RASTERS, presence_test[,c("x","y")]) )
    CLIMATE.abence.test <- as.data.frame( extract(ENV_RASTERS, absence_test[,c("x","y")]) )
    CLIMATE.presence.test <- CLIMATE.presence.test[,pred.var.selection]
    CLIMATE.abence.test <- CLIMATE.abence.test[,pred.var.selection]
    CLIMATE.presence.test <- na.omit(CLIMATE.presence.test)
    CLIMATE.abence.test <- na.omit(CLIMATE.abence.test)
    
    # Extract predictor variables from the ENV raster for the training datasets
    Occurrences.train <- as.data.frame(t(cbind(t(presence_train), t(absence_train))))
    Climate.train <- as.data.frame(extract(ENV_RASTERS, Occurrences.train[,c("x","y")]))
    Climate.train <- Climate.train[,pred.var.selection]
    Occurrences.Climate.train <- cbind(Occurrences.train, Climate.train)
    Occurrences.Climate.train <- na.omit(Occurrences.Climate.train)
  }
  
  ############################### Define the population coordinates for plotting ###############################
  
  coordinates.pops <- Occurrences.Climate[,c("x","y","Presence")]
  if (responce_variable_type == "binary") {
    coordinates.presence <- coordinates.pops[coordinates.pops$Presence == 1, ]
  } else if (responce_variable_type == "frequency") {
    coordinates.presence <- coordinates.pops[coordinates.pops$Presence > 0, ]
  }
  coordinates.absence <- coordinates.pops[coordinates.pops$Presence == 0, ]
  coordinates(coordinates.pops) <- ~x+y
  coordinates(coordinates.presence) <- ~x+y
  if (nrow(coordinates.absence) > 0) {
    coordinates(coordinates.absence) <- ~x+y
    proj4string(coordinates.absence) = CRS(prj_longlat)
  }
  proj4string(coordinates.pops) = CRS(prj_longlat)
  proj4string(coordinates.presence) = CRS(prj_longlat)
  
  ############################### FIT THE DATA TO A GLM MODEL (UNIVARIATE) ###############################
  
  ## Let's create a vector of strings with all predictor names to aid in building the GLM formula:
  #predictors_all <- c("TEMP_MEAN_ANNUAL", "TEMP_DIURNAL_RANGE", "ISOTHERMALITY", "TEMP_SEASONALITY", "TEMP_MAX_WARMEST_MONTH", "TEMP_MIN_COLDEST_MONTH", "TEMP_ANNUAL_RANGE", "TEMP_MEAN_WETTEST_QUARTER", "TEMP_MEAN_DRIEST_QUARTER", "TEMP_MEAN_WARMEST_QUARTER", "TEMP_MEAN_COLDEST_QUARTER", "PREC_ANNUAL", "PREC_WETTEST_MONTH", "PREC_DRIEST_MONTH", "PREC_SEASONALITY", "PREC_WETTEST_QUARTER", "PREC_DRIEST_QUARTER", "PREC_WARMEST_QUARTER", "PREC_COLDEST_QUARTER", "ELEV.GMTED", "PH_5cm", "SLOPE", "ASPECT", "SOIL")
  #predictors_con_sel <- c("TEMP_MEAN_ANNUAL","TEMP_ANNUAL_RANGE","PREC_ANNUAL","PH_5cm","SLOPE")
  #predictors_cat <- c("ASPECT", "SOIL")
  predictors_con_sel <- pred.var.selection
  
  ## Fit univariate Poisson regression models for each predictor and use this model to get an estimate of presence/absence for the values of the predictors calculated previously and plot response curves:
  if (plot_all == TRUE) {
    par(mar = c(2, 1, 2, 1), mfrow = c(3, 3)) 
    for (var in predictors_con_sel) {
      # choose predictor
      vari <- Occurrences.Climate[, var]
      mi <- glm(as.formula(paste("Presence ~ poly(", var, ",2)")), family = binomial,  maxit = 100, weights = GLM.weights, data = Occurrences.Climate)
      data_plot <-data.frame(cbind(Occurrences.Climate[, var],predict(mi,type = "response")))
      sort1 <-na.omit(data_plot[order(data_plot[, 1], decreasing = FALSE),])
      data_plot <-data.frame(cbind(sort1[, 1], sort1[, 2]))
      plot(data_plot[, 1], data_plot[, 2], xlab = "", ylab = "Probability of occurrence", main = var, frame.plot = F,type = "l", ylim =c(0, 1))
      # Plot the presences in red and absences in blue
      if (responce_variable_type == "binary") {
        points(Occurrences.Climate[, var][Occurrences.Climate$Presence == 1], Occurrences.Climate$Presence[Occurrences.Climate$Presence ==1], col = "red")
        points(Occurrences.Climate[, var][Occurrences.Climate$Presence == 0], Occurrences.Climate$Presence[Occurrences.Climate$Presence ==0], col = "blue")
      } else if (responce_variable_type == "frequency") {
        points(Occurrences.Climate[, var][Occurrences.Climate$Presence > 0], Occurrences.Climate$Presence[Occurrences.Climate$Presence > 0], col = "red")
        points(Occurrences.Climate[, var][Occurrences.Climate$Presence == 0], Occurrences.Climate$Presence[Occurrences.Climate$Presence == 0], col = "blue")
      }
    } 
  }
  
  ############################### FIT THE DATA TO A GLM MODEL (MULTIVARIATE) ###############################
  
  ## We fit a Binomial regression model with all the predictors that are available:
  # E.g.:    glm.multi <- glm(Presence ~  poly(TEMP,2) + poly(TEMPRANGE,2) + poly(PREC,2) + poly(ISOTHERM,2), data=Occurrences.Climate, family=binomial, maxit = 100)
  
  # We can write this more generally by first generating the formula we will use
  continuous_var_form <- as.formula(paste("Presence", paste(paste("poly(", predictors_con_sel, ",2)", sep = ""), collapse = " + "), sep = " ~ "))
  #continuous_var <- paste("Presence", paste(paste("poly(", predictors_con_sel, ",2)", sep = ""), collapse = " + "), sep = " ~ ")
  #categorical_var <- paste(" + ", paste(paste(predictors_cat_sel, sep = ""), collapse = " + "), sep = "")
  #interaction_var <- paste(" + ", "SLOPE * ASPECT")
  #all_vars_form <- as.formula (paste(continuous_var, categorical_var, interaction_var, sep = ""))
  
  ## If including categorical data. Before we run the GLM, check how the contrast is set. Our categorical data (aspect and soil) are unordered, hence we have two options: contr.treatment and contr.sum
  #options("contrasts")
  #options(contrasts=c(unordered='contr.sum',ordered='contr.poly'))
  #options(contrasts=c('contr.treatment','contr.poly'))
  
  ## We then run the GLM and explore the results
  glm.multi <- glm(continuous_var_form, data = Occurrences.Climate, family = binomial, maxit = 100, weights = GLM.weights)
  if (plot_all == TRUE) {
    summary(glm.multi)
    ecospat.adj.D2.glm(glm.multi)
  }
  
  # Next, we get rid of the non-significant variables, to make the model "parsimonious"
  glm.multi.step <- step(glm.multi, directions = "both", trace = 0)
  if (plot_all == TRUE) {
    summary(glm.multi.step)
    ecospat.adj.D2.glm(glm.multi.step)
  }
  
  ################ Model diagnostics ################
  
  ## First we obtain Cook’s distance and leverage for all observations using glm.diag function:
  #glm diagnostics- in particular for cook and leverage plot diags <- glm.diag(glm.multi)
  diags <- glm.diag(glm.multi.step)
  cook <- diags$cook
  levg <- diags$h
  stdr.levg <- levg/(1-levg)
  ## Next we set cook and leverage thresholds:
  n <- nrow(Occurrences.Climate)
  p <- glm.multi.step$df.null-glm.multi.step$df.residual+1
  cook.thrsh <- 8/(n-2*p)
  levg.thrsh <- 2*p/(n-2*p)
  ## Next we obtain linear predictors βT xi and Pearson residuals ri for all observations i:
  xx.fit <- predict(glm.multi.step, type="response") 
  yy.fit <- residuals(glm.multi.step, type="pearson")
  ## Now we create smoother of linear predictor vs. Pearson residuals:
  ls.fit <- loess.smooth(xx.fit, yy.fit, family="gaussian")
  if (plot_all == TRUE) {
    ## Finally we create Tukey-Anscombe and standardized levarage vs. Cookś distance plots:
    # Tukey-Anscombe plot with pearson residuals
    par(mar = c(5, 5, 5, 5), mfrow = c(1, 2)) 
    # Plot linear predictor vs. pearson residuals
    plot(xx.fit, yy.fit, xlab="linear predictor", ylab="Pearson residuals", main="Tukey-Anscombe plot")
    # Add smoother for real Pearson residuals (possibly having structure) vs. linear predictor lines(ls.fit$x, ls.fit$y, col="red")
    abline(h=0, lty=3)
    # Cook distance and leverage plot
    plot(stdr.levg, cook, xlab="Standardized leverage", ylab="Cook statistic", main="Cook statistic - leverage plot") 
    abline(h=cook.thrsh, v=levg.thrsh, lty=2)
  }
  
  # Identify influential outliers:
  indx.out <- which(stdr.levg > levg.thrsh & cook > cook.thrsh) 
  length(indx.out)/nrow(Occurrences.Climate)*100 # this is the percentage of points which are considered outliers
  outs <- Occurrences.Climate[indx.out,]
  
  if (plot_all == TRUE) {
    # Plot outliers:
    par(mfrow=c(1,1))
    par(mar=c(3,3,3,3))
    plot(ENV_RASTERS$TEMP_MEAN_ANNUAL, col = rev(heat.colors(20)), box = F, axes = F, main="GLM outlier populations") 
    #plot(ch.sp, add=T)
    if (responce_variable_type == "binary") {
      points(Occurrences.Climate$x, Occurrences.Climate$y, col=c("blue","green4")[Occurrences.Climate$Presence+1], pch=3, cex=0.5) 
      points(outs$x, outs$y, cex=3, col="black")
      legend(7.5e5, 9e4, legend=c("absent","present", "outlier"), pch=c(3,3,1),col=c("blue","green4","black"))
    } else if (responce_variable_type == "frequency") {
      points(Occurrences.Climate$x, Occurrences.Climate$y, col="blue", pch=3, cex=0.5) 
      points(outs$x, outs$y, cex=3, col="black")
      legend(7.5e5, 9e4, legend=c("population", "outlier"), pch=c(3,3,1),col=c("blue","black"))
    }
    ## Consider removing outliers, if outliers seem erroneous.
  }
  
  ################ Model evaluation, optimisation and projection ################
  
  if (responce_variable_type == "binary") {
    ## First create a dataframe containing the observed presence/absence and the predicted probability of occurrence for the response variable at each sampled location:
    data_validation.glm <- data.frame(cbind(Occurrences.Climate$Presence, predict(glm.multi.step, type = "response")))
    colnames(data_validation.glm) <- c("Observed", "Projected")
    
    if (plot_all == TRUE) {
      ## Load glm.model.eval.r file with necessary functions for model evaluation. Compute the table of accuracy statistics with an arbitrary threshold of 0.5:
      source("./glm.model.eval.r")
      meva.table(data_validation.glm$Projected, data_validation.glm$Observed, 0.5)
    }
    
    ## Now, we want to test, how good the model actually is. For this, we first generate a function called 'cv.model()'. Don't worry, if you don't understand what it does. Be aware, that you need to have the 'dismo' package installed for this function
    cv.model <- function(model, K, data = model$data){
      require(dismo) # to ensure the dismo package is attached
      ks <- kfold(model$data, k = K, by = model$data[,as.character(formula(model)[2])])
      cvpreds <- data.frame(row = row.names(data), cvpred = numeric(length = nrow(data)))
      for(i in 1:K){
        train <- data[ks!=i,]
        test <- data[ks==i,]
        modtmp <- update(model, data = train)
        cvpreds[which(ks==i),2] <- predict(modtmp, newdata = test, type = 'response')
      }
      cvpreds
    }
    
    ## Now we apply the crossvalidation function to our step-wise optimized model with a 5-fold splitting into training and testing
    xval.glm.multi.step <- cv.model(glm.multi.step, K = 5)
    
    if (plot_all == TRUE) {
      ## And we plot the fitted values from the step-wise and the cross-validated models
      par(mfrow=c(1,2), mar = c(5, 5, 5, 5))
      plot(glm.multi.step$fitted.values, xval.glm.multi.step$cvpred, xlab = 'fitted values from stepwise optimize model', ylab = 'predicted values form cross-validation', main="Model evaluation - GLM")
      abline(0, 1, lwd = 3, col = "red")
    }
    
    ## Now we want to explore what threshold gives best predictions. For this we use the PresenceAbsence package. To prepare, we generate three columns: ID, observed, predicted!
    # We do so for the stepwise optimized model, and for the xval predicted model test.
    glm.multi.step.test <- data.frame(ID = 1:nrow(Occurrences.Climate), observed = Occurrences.Climate$Presence,
                                      predicted = glm.multi.step$fitted)
    xval.glm.multi.step.test <- data.frame(ID = 1:nrow(Occurrences.Climate), observed = Occurrences.Climate$Presence,
                                           predicted = xval.glm.multi.step$cvpred)
    
    ## Evaluate the optimal thresholds for the stepwise optimized model. We've changed the threshhold from 1001 to the default of 0.5.
    thres.glm.multi.step <- optimal.thresholds(glm.multi.step.test, opt.methods = 1:9)
    #thres.glm.multi.step
    # ... and for the crossvalidated model
    thres.glm.multi.xval <- optimal.thresholds(xval.glm.multi.step.test, opt.methods = 1:9)
    #thres.glm.multi.xval
    
    if (plot_all == TRUE) {
      ## Print max Kappa for both, the stepwise optimized and the crossvalidated model output
      Kappa(cmx(glm.multi.step.test, threshold = thres.glm.multi.step[4,2]))
      Kappa(cmx(xval.glm.multi.step.test, threshold = thres.glm.multi.xval[4,2]))
      
      ## Calculate threshold-independent ROC plot and Area Under (Receiver Operator) Curve (AUC)
      str(roc(glm.multi.step.test$predicted, as.factor(glm.multi.step.test$observed)))
    }
    
    ## Now we calculate the ROC Curves (AUC) from the stepwise and the x-validated models
    roc.glm.multi.step <- roc(glm.multi.step.test$predicted, as.factor(glm.multi.step.test$observed))
    roc.glm.multi.xval <- roc(xval.glm.multi.step.test$predicted, as.factor(xval.glm.multi.step.test$observed))
    
    if (plot_all == TRUE) {
      ## Now let’s compare two GLMs in terms of their ROC-curves. . .
      plot(roc.glm.multi.step, main="ROC curve - GLM", col = "grey20", lwd = 5, lty = 1)
      plot(roc.glm.multi.xval, col = "grey70", lwd = 3, lty = 2, add = TRUE)
      legend("bottomright", fill = c("grey20", "grey70"), legend = c("Stepwise optimised model", "Cross-validated model"), xpd = TRUE,bty = "n")
      
      ## ... and let's print the AUC values for the two models:
      AUC::auc(roc.glm.multi.step)
      AUC::auc(roc.glm.multi.xval)
    }
  }
  
  ## Projection of the species distribution model over the Alps
  # Old version (may encounter errors in large raster datasets) - if errors encountered, see revised code below!
  # We first need the environmental rasters for which we'll project the model on. 
  Projection<- as.data.frame(rasterToPoints(ENV_RASTERS))
  Projection <- Projection[,c("x","y",pred.var.selection)]
  Projection <- na.omit(Projection)
  # Don't forget to convert non-numeric data types (categorical variable) to factors, 
  #Projection0[,10] <- as.factor(Projection0[,10])
  # We convert the projected environmental dataframe to a raster
  Projection.raster <- rasterFromXYZ(Projection)
  
  # # New, revised version (error-free, faster for large datasets) - written for Carex
  # # We first need the environmental rasters for which we'll project the model on.
  # # We must make sure that there are no NA values in the raster layers, otherwise the predict function will throw errors.
  # # We do this by synchronising NA among the layers of the raster stack. We achieve this by recording the location of all cells with an NA value in at least one of the stack's layers. mask()then allows us to apply that mask to every layer in the stack.
  # # Don't forget to convert non-numeric data types (categorical variable) to factors (if applicable)
  # ENV_RASTERS.SEL <- subset(ENV_RASTERS, pred.var.selection)
  # Projection.raster <- mask(ENV_RASTERS.SEL, calc(ENV_RASTERS.SEL,fun = sum))
  # #rm(ENV_RASTERS.SEL, ENV_RASTERS)
  
  ## We have fitted the model from a few hundred points. Let's now project the model to the whole Alps. For that, we need all predictor information per cell, which is what we have generated above. For this, we use the 'predict' function.
  pred.glm.multi.prob <- predict(object = Projection.raster, model = glm.multi.step,  type = 'response')
  
  if (responce_variable_type == "binary") {
    ## If we cut probabilities with the optimal threshold, then we can plot a presence/absence map
    pred.glm.multi.bin <- pred.glm.multi.prob > thres.glm.multi.xval[4,2]
  }
  
  ############################### FIT THE DATA TO A GAM MODEL (MULTIVARIATE) ###############################
  
  ## Now we repeat the same procedure for GAM; for this we load the library 'gam':
  # How many degrees of freedom should we allow?
  #continuous_var.gam <- paste("Presence", paste("s(", predictors_con_sel, ", df = 2)", collapse = " + "), sep = " ~ ")
  #categorical_var.gam <- paste(" + ", paste(paste(predictors_cat_sel, sep = ""), collapse = " + "), sep = "")
  #form.gam <- as.formula (paste(continuous_var.gam, categorical_var.gam, sep = ""))
  form.gam <- as.formula(paste("Presence", paste("s(", predictors_con_sel, ", df = 2)", collapse = " + "), sep = " ~ "))
  gam.full <- gam(form.gam, data = Occurrences.Climate, family = "binomial", weights = GLM.weights)
  if (plot_all == TRUE) {
    ecospat.adj.D2.glm(gam.full)
    summary(gam.full)
  }
  
  ## Now we apply stepwise variable selection, which, for GAMs, requires to set up a so called scope. The scope defines the models to be considered and can be produced with ghe gam-builtin function gam.scope:
  scope.gam <- gam.scope(Occurrences.Climate[, predictors_con_sel], arg = "df = 2")
  gam.step <- step.Gam(gam.full, scope = scope.gam, direction = "both", trace = 0)
  if (plot_all == TRUE) {
    ecospat.adj.D2.glm(gam.step)
    summary(gam.step)
  }
  
  if (responce_variable_type == "binary") {
    ## Model evaluation - compute the table of accuracy statistics with an arbitrary threshold of 0.5:
    data_validation.gam <- data.frame(cbind(Occurrences.Climate$Presence, predict(gam.step, type = "response")))
    colnames(data_validation.gam) <- c("Observed", "Projected")
    if (plot_all == TRUE) {
      meva.table(data_validation.gam$Projected, data_validation.gam$Observed, 0.5)
    }
    
    ## And we test the model by means of a cross-validation.
    xval.gam.step <- cv.model(gam.step, K = 5)
    if (plot_all == TRUE) {
      par(mfcol = c(1, 2))
      plot(gam.step$fitted.values, xval.gam.step$cvpred, xlab = 'fitted values from stepwise optimize model', ylab = 'predicted values form cross-validation', main="Model evaluation - GAM")
      abline(0, 1, lwd = 3, col = "red")
    }
    
    ## Now we want to explore what threshold gives best predictions. For this we use the PresenceAbsence package. To prepare, we generate three columns: ID, observed, predicted!
    # We do so for the stepwise optimized model, and for the xval predicted model test.
    gam.step.test <- data.frame(ID = 1:nrow(Occurrences.Climate), observed = Occurrences.Climate$Presence,
                                predicted = gam.step$fitted)
    xval.gam.step.test <- data.frame(ID = 1:nrow(Occurrences.Climate), observed = Occurrences.Climate$Presence,
                                     predicted = xval.gam.step$cvpred)
    
    ## Evaluate the optimal thresholds for the stepwise optimized model. We've changed the threshhold from 1001 to the default of 0.5.
    thres.gam.step <- optimal.thresholds(gam.step.test, opt.methods = 1:9)
    #thres.gam.step
    # ... and for the crossvalidated model
    thres.gam.xval <- optimal.thresholds(xval.gam.step.test, opt.methods = 1:9)
    #thres.gam.xval
    
    if (plot_all == TRUE) {
      ## Print max Kappa for both, the stepwise optimized and the crossvalidated model output
      Kappa(cmx(gam.step.test, threshold = thres.gam.step[4,2]))
      Kappa(cmx(xval.gam.step.test, threshold = thres.gam.xval[4,2]))
      
      ## Calculate threshold-independent ROC plot and Area Under (Receiver Operator) Curve (AUC)
      str(roc(gam.step.test$predicted, as.factor(gam.step.test$observed)))
    }
    
    ## Now we calculate the ROC Curves (AUC) from the stepwise and the x-validated models
    roc.gam.step <- roc(gam.step.test$predicted, as.factor(gam.step.test$observed))
    roc.gam.xval <- roc(xval.gam.step.test$predicted, as.factor(xval.gam.step.test$observed))
    
    if (plot_all == TRUE) {
      ## Now let’s compare two GLMs in terms of their ROC-curves. . .
      plot(roc.gam.step, main="ROC curve - GAM", col = "grey20", lwd = 5, lty = 1)
      plot(roc.gam.xval, col = "grey70", lwd = 3, lty = 2, add = TRUE)
      legend("bottomright", fill = c("grey20", "grey70"), legend = c("Stepwise optimised model", "Cross-validated model"), xpd = TRUE,bty = "n")
      
      ## ... and let's print the AUC values for the two models:
      AUC::auc(roc.gam.step)
      AUC::auc(roc.gam.xval)
    }
  }
  
  ## Finally, we generate the spatial map projection of the optimized and tested GAM model
  pred.gam.prob <- predict(object = Projection.raster, model = gam.step,type = "response")
  
  if (responce_variable_type == "binary") {
    pred.gam.bin <- pred.gam.prob > thres.gam.step[4,2]
  }
  
  ############################### FIT THE DATA TO A MAXIMUM ENTROPY MODEL ###############################
  
  # NOTE 1 - This implementation currently uses (sources) the Maxent Java package. As an alternative (same authors, same methodology), you can now use the R package Package ‘maxnet’. 
  # NOTE 2 - Maxent only works with presence data!
  if (responce_variable_type == "binary") {
    
    # Choose between running maxent (require maxent installation and external .jar file) or maxnet (independent R package, with some slightly different implementations)
    if (maxent_maxnet == "maxent") {
      
      # First we need to call the maxent.jar file
      jar <- paste(system.file(package="dismo"), "/java/maxent.jar", sep='')
      
      # Fit a maximum entropy model to the data. In case you have a categorical variable, don't forget to e.g. add the argument: factors='SOIL'
      maxent.prob <- maxent(x=Occurrences.Climate[,c(first.pred.var.idx:last.pred.var.idx)], p=Occurrences.Climate$Presence)
      
      if (plot_all == TRUE) {
        # Plot variable contribution
        par(mfrow=c(1,1), mar = c(5, 5, 5, 5))
        plot(maxent.prob, main = "Variable contribution - Maxent")
        # Plot univariate response curves
        response(maxent.prob, main = "Response curve - Maxent")
      }
      
      # Evaluate model. Note that there seems to be a bug with evaluating the maxnet model.
      maxent.train <- maxent(x=Occurrences.Climate.train[,c(first.pred.var.idx:last.pred.var.idx)], p=Occurrences.Climate.train$Presence)
      eval.maxent <- evaluate(CLIMATE.presence.test, CLIMATE.abence.test, maxent.train)
      if (plot_all == TRUE) {
        eval.maxent
      }
      
      # Predict and project to climate rasters
      pred.maxent.prob <- predict(object = Projection.raster, model = maxent.prob,  type = 'response')
      
      # Set threshold for binary projection
      thres.maxent <- threshold(eval.maxent, 'kappa')
      pred.maxent.bin <- (pred.maxent.prob > thres.maxent)
      
    } else if (maxent_maxnet == "maxnet") {
      
      # Load maxnet library
      library("maxnet")
      
      # Fit a maximum entropy model to the data. In case you have a categorical variable, don't forget to e.g. add the argument: factors='SOIL'
      maxent.prob <- maxnet(data=Occurrences.Climate[,c(first.pred.var.idx:last.pred.var.idx)], p=Occurrences.Climate$Presence)
      
      if (plot_all == TRUE) {
        # Plot univariate response curves
        par(mfrow=c(1,1), mar = c(5, 5, 5, 5))
        plot(maxent.prob, main = "Variable contribution - Maxent")
      }
      
      # Evaluate model. Note that there seems to be a bug with evaluating the maxnet model.
      maxent.train <- maxnet(data=Occurrences.Climate.train[,c(first.pred.var.idx:last.pred.var.idx)], p=Occurrences.Climate.train$Presence)
      eval.maxent <- evaluate(CLIMATE.presence.test, CLIMATE.abence.test, maxent.train)
      if (plot_all == TRUE) {
        eval.maxent # this gives weird results for max TPR+TNR
      }
      
      # Predict and project to climate rasters
      pred.maxent.prob <- predict(object = Projection.raster, model = maxent.prob,  type = 'cloglog') # cloglog (default) and logistic give very similar results
      
      # Set threshold for binary projection
      thres.maxent <- threshold(eval.maxent, 'prevalence') # because of erroneous output for max TPR+TNR, we use another statistic instead to estimate threshold
      if (thres.maxent < 0.25 || thres.maxent > 0.75) {  # the eval function has a bug with the new maxnet implementation. In case it return an unrealistic value, we default the threshold to 0.5
        thres.maxent <- 0.5
      }
      pred.maxent.bin <- (pred.maxent.prob > thres.maxent)
    }
  }
  
  ############################### FIT THE DATA TO A RANDOM FOREST MODEL ###############################
  
  # Fit a random forest model to the data. Here we use the default values for mtry and ntrees, which seem reasonable.
  rforest.form <- as.formula(paste("Presence", paste(predictors_con_sel, collapse = " + "), sep = " ~ "))
  #rforest.form <- as.formula(paste("Presence", paste(names(ENV_RASTERS), collapse = " + "), sep = " ~ "))
  rf.prob <- randomForest(rforest.form, data=Occurrences.Climate)
  rf.prob.importance <- importance(rf.prob)
  
  if (plot_all == TRUE) {
    # Plot random forest error vs ntrees
    par(mfrow = c(1, 1), mar = c(5, 5, 5, 5))
    plot(rf.prob, main = "Tree optimisation - Random Forest")
    randomForest::varImpPlot(rf.prob)
    #print(importance(rf.prob))
  }
  
  # Cross-validation for feature selection
  eval.rf.cv <- rfcv(trainx = Occurrences.Climate[,c(first.pred.var.idx:last.pred.var.idx)], trainy = Occurrences.Climate$Presence, cv.fold = 5)
  if (plot_all == TRUE) {
    head(eval.rf.cv, n = 2L)
  }
  # Evaluate model
  if (responce_variable_type == "binary") {
    rf.train <- randomForest(rforest.form, data=Occurrences.Climate.train)
    eval.rf <- evaluate(CLIMATE.presence.test, CLIMATE.abence.test, rf.train)
    eval.rf
  }
  
  # Tune randomforest model by optimising mtry parameter
  #fgl.res <- tuneRF(Occurrences.Climate[,4:8], Occurrences.Climate[,3], stepFactor=3)
  
  # Predict and project to climate rasters
  pred.rf.prob <- predict(object = Projection.raster, model = rf.prob,  type = 'response')
  
  if (responce_variable_type == "binary") {
    # Set threshold for binary projection
    thres.rf <- threshold(eval.rf, 'kappa')
    pred.rf.bin <- (pred.rf.prob > thres.rf)
  }
  
  ############################### Plot all models (CONTEMPORARY CLIMATE) ###############################
  
  # Make a raster stack of of all model projections
  if (responce_variable_type == "binary") {
    pred.allModels.prob <- stack(pred.glm.multi.prob, pred.gam.prob, pred.maxent.prob, pred.rf.prob)
    pred.allModels.bin <- stack(pred.glm.multi.bin, pred.gam.bin, pred.maxent.prob > thres.maxent, pred.rf.prob > thres.rf)
    names(pred.allModels.prob) <- c("GLM.Today","GAM.Today","MAXENT.Today","Random Forest.Today")
    names(pred.allModels.bin) <- c("GLM.Today)","GAM.Today","MAXENT.Today","RandomForest.Today")
    if (plot_all == TRUE) {
      # And plot
      plot(pred.allModels.prob, col = brewer.pal(9, "YlOrBr"), legend.args = list(text = "Predicted Probability", side = 2, font = 4, cex = 0.9, line = 0.5), legend.width = 2, legend.mar = 5, zlim = c(0, 1), legend.shrink = 0.9)
      plot(pred.allModels.bin, col = c("grey90", "green4"), legend.width = 2, legend.mar=5, legend=FALSE)
      legend("bottomleft", fill = c("green4", "grey90"), legend = c("Presence", "Absence"), xpd = TRUE,bty = "n")
    }
  } else if (responce_variable_type == "frequency") {
    #pred.allModels.prob <- stack(pred.glm.multi.prob, pred.gam.prob, pred.rf.prob)
    #pred.allModels.prob <- stack(pred.gam.prob, pred.rf.prob)
    pred.allModels.prob <- stack(pred.rf.prob)
    #names(pred.allModels.prob) <- c("GLM.Today","GAM.Today","Random Forest.Today")
    #names(pred.allModels.prob) <- c("GAM.Today","Random Forest.Today")
    names(pred.allModels.prob) <- c("Random Forest.Today")
    if (plot_all == TRUE) {
      plot(pred.allModels.prob, col = brewer.pal(9, "YlOrBr"), legend.args = list(text = "Predicted Probability", side = 2, font = 4, cex = 0.9, line = 0.5), legend.width = 2, legend.mar = 5, zlim = c(0, 1), legend.shrink = 0.9)
    }
  }
  
  ## Generate an ensemble model
  #  Consider if you want to weight all models equally, or if you want to weight the models e.g. by AUC. Here we weight the model contribution by their AUC values.
  if (responce_variable_type == "binary") {
    auc_glm <- AUC::auc(roc.glm.multi.step)
    auc_gam <- AUC::auc(roc(gam.step.test$predicted, as.factor(gam.step.test$observed)))
    auc_maxent <- eval.maxent@auc
    auc_rf <- eval.rf@auc
    auc_mean <- (auc_glm + auc_gam + auc_maxent + auc_rf) / 4
    thres.mean <- (auc_glm*thres.glm.multi.step[4, 2] + auc_gam*thres.gam.step[4, 2] + auc_maxent*thres.maxent + auc_rf*thres.rf) / (auc_mean*4)
  }
  # Take the weighted average of the predictions
  if (responce_variable_type == "binary") {
    pred.prob.current <-(auc_glm*pred.glm.multi.prob + auc_gam*pred.gam.prob + auc_maxent*pred.maxent.prob + auc_rf*pred.rf.prob) / (auc_mean*4)
    par(mfrow=c(1,2), mar = c(3, 3, 3, 3))
  } else if (responce_variable_type == "frequency") {
    #pred.prob.current <- (pred.glm.multi.prob + pred.gam.prob + pred.rf.prob) / 3
    #pred.prob.current <- (pred.gam.prob + pred.rf.prob) / 2
    pred.prob.current <- pred.rf.prob
    par(mfrow=c(1,1), mar = c(3, 3, 3, 3))
  }
  # Plot the ensemble model
  if (plot_all == TRUE) {
    plot(pred.prob.current, col = brewer.pal(9, "YlOrBr"), main = "Continuous Spatial Projection (Ensemble Model - Today)",
         legend.args = list(text = "Predicted Probability", side = 2, font = 4, cex = 0.9, line = 0.5), legend.width = 2, legend.mar = 5, zlim = c(0, 1), legend.shrink = 0.9)
    #plot(borders, add=TRUE,lwd=0.2,lty=1)
    points(coordinates.presence, pch = 16, cex = 0.6)
  }
  if (responce_variable_type == "binary") {
    pred.bin.current <- pred.prob.current > thres.mean
    if (plot_all == TRUE) {
      plot(pred.bin.current, col = c("grey90", "green4"), main="Binary Spatial Projection (Ensemble Model - Today)", legend.width = 2, legend.mar=5, legend=FALSE)
      #plot(borders, add=TRUE,lwd=0.2,lty=1)
      legend("bottomleft", fill = c("green4", "grey90"), legend = c("Presence", "Absence"), xpd = TRUE,bty = "n")
      points(coordinates.presence, pch = 16, cex = 0.6)
    }
  }
  
  # To mask allele distribution models, we need the species distribution model. Hence run this first.
  if (responce_variable_type == "binary") {
    pred.bin.current.species <- pred.bin.current
  }
  
  if (responce_variable_type == "frequency")  {
    if (exists("pred.bin.current.species") == TRUE) {
      # Constrain the allele projected distribution to where the species occurs (the binary species projection)
      # To mask projected CEN raster by the project species SDM raster, we use the mask function. We mask as such because the occurrence of the alleles is contingent on the occurence of the species.
      # E.g. if we defined the projected CEN raster as pred.bin.current and the species raster as pred.bin.current.species, we can mask as follows (for the binary):
      #pred.bin.current.masked<-mask(pred.bin.current, pred.bin.current.species, maskvalue=0, updatevalue = 0)
      # And similarly for the continious prediction (whole species mask)
      pred.prob.current.masked<-mask(pred.prob.current, pred.bin.current.species, maskvalue=0, updatevalue = 0)
      if (plot_all == TRUE) {
        # We can then plot the masked prediction raster (using whole species distribution mask)
        proj4string(pred.prob.current.masked) = CRS(prj_longlat)
        zeroCol <-"grey90"
        #reds <- brewer.pal('YlOrRd', n = 9)
        myPalette<-brewer.pal(11,"RdYlBu")
        my.at=seq(0, 1, by=0.1)
        #my.brks=seq(0, 1, by=0.1)
        my.brks=seq(0.05, 1.05, by=0.1)
        myLabels <- c("NA (species not present)", "0.9 - low allele", "0.8 - low allele", "0.7 - low allele", "0.6 - low allele", "0.5 - high/low allele", "0.6 - high allele", "0.7 - high allele", "0.8 - high allele", "0.9 - high allele", "1.0 - high allele")
        #myLabels <- c("NA (species not present)", "0.1/0.9 - high/low allele", "0.2/0.8 - high/low allele", "0.3/0.7 - high/low allele", "0.4/0.6 - high/low allele", "0.5/0.5 - high/low allele", "0.6/0.4 - high/low allele", "0.7/0.3 - high/low allele", "0.8/0.2 - high/low allele", "0.9/0.1 - high/low allele", "1.0/0.0 - high/low allele")
        #myColorkey <- list(at=my.brks, labels=list(at=my.brks, labels=my.at))
        myColorkey <- list(at=my.brks, labels=list(at=my.at, labels=myLabels))
        #myTheme <- rasterTheme(region = c(zeroCol, reds))
        myTheme <- rasterTheme(region = c(zeroCol, myPalette))
        #e_cen_plot <- extent(5, 14, 43, 48)
        e_cen_plot <- extent(5, 12, 43.5, 48)
        pred.prob.current.masked <- crop(pred.prob.current.masked, e_cen_plot)
        pred.prob.current.masked.plot <- levelplot(pred.prob.current.masked, par.settings = myTheme, at=my.at, colorkey=myColorkey, maxpixels = 1e7, margin = FALSE, main = "Continuous Spatial Projection (Masked Ensemble Model - Today)")
        pred.prob.current.masked.plot
      }
    }
    # For Central alpine mask
    pred.prob.current.masked.CentralAlps <-mask(pred.prob.current, CentralAlps_current_mask, maskvalue=0, updatevalue = 0)
    if (plot_all == TRUE) {
      # We can then plot the masked prediction raster (using Central alpine distribution mask)
      proj4string(pred.prob.current.masked.CentralAlps) = CRS(prj_longlat)
      pred.prob.current.masked.CentralAlps <- crop(pred.prob.current.masked.CentralAlps, e_cen_plot)
      pred.prob.current.masked.CentralAlps.plot <- levelplot(pred.prob.current.masked.CentralAlps, par.settings = myTheme, at=my.at, colorkey=myColorkey, maxpixels = 1e7, margin = FALSE, main = "Continuous Spatial Projection (Masked Ensemble Model - Today)")
      pred.prob.current.masked.CentralAlps.plot
    }
    #grid.arrange(pred.prob.current.masked.plot, pred.prob.current.masked.CentralAlps.plot, ncol=2)
  }
  
  if (responce_variable_type == "binary") {
    if (plot_all == TRUE) {
      # Print AUC of all models
      print(paste("AUC (GLM): ", auc_glm, sep = " "))
      print(paste("AUC (GAM): ", auc_gam, sep = " "))
      print(paste("AUC (MAXENT): ", auc_maxent, sep = " "))
      print(paste("AUC (RF): ", auc_rf, sep = " "))
      print(paste("AUC (ENSEMBLE): ", auc_mean, sep = " "))
    }
  }
  
  time_stamp <- gsub(":","_",gsub(" ", "_", Sys.time()))
  saveRDS(pred.prob.current, paste0("pred.prob.current.",lineage,".selectVar.optimalThreshold.",time_stamp,".rds"))
  if (responce_variable_type == "binary") {
    saveRDS(pred.bin.current, paste0("pred.bin.current.",lineage,".selectVar.optimalThreshold.",time_stamp,".rds"))
  }
  
  ############################### HINDCAST THE MODELS TO THE LGM / FORECAST THE MODELS TO THE FUTURE ###############################

  ## We load the climactic rasters for the LGM (these tiff files are already defined with CRS = WGS84). These are the 19 Bioclim variables from Chelsa (http://chelsa-climate.org/), and in 2 byte integers (integer) precision (higher precision floating point also available but this takes much more space)

  # Define whether to hindcast or to forecast
  if (cast == "hindcast") {
    #ENV_LGM_RASTERS <- readRDS(climate_model_PMIP3)
    ENV_LGM_RASTERS <- stack(climate_model_PMIP3)
  } else if (cast == "forecast") {
    ENV_LGM_RASTERS <- readRDS(climate_model_CMIP5)
  }

  if (projection_extent == "alps") {
    ENV_LGM_RASTERS <- crop(ENV_LGM_RASTERS, e_alps)
  } else if (projection_extent == "europe") {
    ENV_LGM_RASTERS <- crop(ENV_LGM_RASTERS, e_europe)
  }
  #plot(ENV_LGM_RASTERS)

  # Make the LGM environmental projection dataframe
  Projection.LGM<- as.data.frame(rasterToPoints(ENV_LGM_RASTERS))
  Projection.LGM <- Projection.LGM[,c("x","y",pred.var.selection)]
  Projection.LGM <- na.omit(Projection.LGM)
  Projection.LGM.raster <- rasterFromXYZ(Projection.LGM)

  ## Now let’s see how the three models project the past distribution of Dianthus sylvestris
  # For this, we just create predicted distributions for the past using the fitted models and past environmental data.
  pred.glm.prob.LGM <- predict(object = Projection.LGM.raster, model = glm.multi.step, type = "response")
  pred.gam.prob.LGM <- predict(object = Projection.LGM.raster, model = gam.step, type = "response")
  pred.rf.prob.LGM <- predict(object = Projection.LGM.raster, model = rf.prob, type = "response")

  if (responce_variable_type == "binary") {
    if (maxent_maxnet == "maxent") {
      pred.maxent.prob.LGM <- predict(object = Projection.LGM.raster, model = maxent.prob, type = "response")
    } else if (maxent_maxnet == "maxnet") {
      pred.maxent.prob.LGM <- predict(object = Projection.LGM.raster, model = maxent.prob, type = "cloglog") # cloglog (default) and logistic give very similar results
    }
    pred.maxent.bin.LGM <- (pred.maxent.prob.LGM > thres.maxent)
    pred.glm.bin.LGM <- pred.glm.prob.LGM > thres.glm.multi.step[4, 2]
    pred.gam.bin.LGM <- pred.gam.prob.LGM > thres.gam.step[4, 2]
    pred.rf.bin.LGM <- (pred.rf.prob.LGM > thres.rf)
  }

  # Plot all models for the LGM
  if (responce_variable_type == "binary") {
    pred.allModels.prob.LGM <- stack(pred.glm.prob.LGM, pred.gam.prob.LGM, pred.maxent.prob.LGM, pred.rf.prob.LGM)
    pred.allModels.bin.LGM <- stack(pred.glm.bin.LGM, pred.gam.bin.LGM, pred.maxent.bin.LGM, pred.rf.bin.LGM)
    names(pred.allModels.prob.LGM) <- c(paste0("GLM.",cast.var.name),paste0("GAM.",cast.var.name),paste0("MAXENT.",cast.var.name),paste0("RandomForest.",cast.var.name))
    names(pred.allModels.bin.LGM) <- c(paste0("GLM.",cast.var.name),paste0("GAM.",cast.var.name),paste0("MAXENT.",cast.var.name),paste0("RandomForest.",cast.var.name))
    if (plot_all == TRUE) {
      plot(pred.allModels.prob.LGM, col = brewer.pal(9, "YlOrBr"), legend.args = list(text = "Predicted Probability", side = 2, font = 4, cex = 0.9, line = 0.5), legend.width = 2, legend.mar = 5, zlim = c(0, 1), legend.shrink = 0.9)
      plot(pred.allModels.bin.LGM, col = c("grey90", "green4"), legend.width = 2, legend.mar=5, legend=FALSE)
      legend("bottomleft", fill = c("green4", "grey90"), legend = c("Presence", "Absence"), xpd = TRUE,bty = "n")
    }
  } else if (responce_variable_type == "frequency") {
    #pred.allModels.prob.LGM <- stack(pred.glm.prob.LGM, pred.gam.prob.LGM, pred.rf.prob.LGM)
    #pred.allModels.prob.LGM <- stack(pred.gam.prob.LGM, pred.rf.prob.LGM)
    pred.allModels.prob.LGM <- stack(pred.rf.prob.LGM)
    #names(pred.allModels.prob.LGM) <- c("GLM.LGM","GAM.LGM","Random Forest.LGM")
    #names(pred.allModels.prob.LGM) <- c("GAM.LGM","Random Forest.LGM")
    names(pred.allModels.prob.LGM) <- c("Random Forest.LGM")
    if (plot_all == TRUE) {
      plot(pred.allModels.prob.LGM, col = brewer.pal(9, "YlOrBr"), legend.args = list(text = "Predicted Probability", side = 2, font = 4, cex = 0.9, line = 0.5), legend.width = 2, legend.mar = 5, zlim = c(0, 1), legend.shrink = 0.9)
    }
  }

  # Generate an ensemble model for the LGM
  if (responce_variable_type == "binary") {
    pred.prob.LGM <-(auc_glm*pred.glm.prob.LGM + auc_gam*pred.gam.prob.LGM + auc_maxent*pred.maxent.prob.LGM + auc_rf*pred.rf.prob.LGM) / (auc_mean*4)
    par(mfrow=c(1,2), mar = c(3, 3, 3, 3))
  } else if (responce_variable_type == "frequency") {
    #pred.prob.LGM <- (pred.glm.prob.LGM + pred.gam.prob.LGM + pred.rf.prob.LGM) / 3
    #pred.prob.LGM <- (pred.gam.prob.LGM + pred.rf.prob.LGM) / 2
    pred.prob.LGM <- pred.rf.prob.LGM
    par(mfrow=c(1,1), mar = c(3, 3, 3, 3))
  }
  # Plot the ensemble model - for LGM
  if (plot_all == TRUE) {
    plot(pred.prob.LGM, col = brewer.pal(9, "YlOrBr"), main = paste0("Continuous Spatial Projection (Ensemble Model - ", cast.var.name, ")"),
         legend.args = list(text = "Predicted Probability", side = 2, font = 4, cex = 0.9, line = 0.5),
         legend.width = 2, legend.mar = 5, zlim = c(0, 1),
         legend.shrink = 0.9)
  }
  #plot(borders, add=TRUE,lwd=0.2,lty=1)
  if (responce_variable_type == "binary") {
    # if (lineage == "CentralAlps") {
    #   thres.mean <- 0.4 # we may want to lower/relax the threshold if model is overfitted and hindcast to restricted
    # }
    pred.bin.LGM <- pred.prob.LGM > thres.mean
    if (plot_all == TRUE) {
      plot(pred.bin.LGM, col = c("grey90", "green4"), main=paste0("Binary Spatial Projection (Ensemble Model - ", cast.var.name, ")"),
           legend.width = 2, legend.mar=5, legend=FALSE)
      #plot(borders, add=TRUE,lwd=0.2,lty=1)
      legend("bottomleft", fill = c("green4", "grey90"),
             legend = c("Presence", "Absence"), xpd = TRUE,bty = "n")
    }
  }

  # Append results (pred.prob.LGM and treshold values) to stack
  pred.prob.LGM.stack <- stack(pred.prob.LGM.stack, pred.prob.LGM)
  if (responce_variable_type == "binary")  {
    thres.stack[run] <- thres.mean
    auc.stack[run] <- auc_mean
  }
}

############################### END OF MAIN SDM FUNCTION ###############################

############################### CALCULATE AND PLOT AVERAGE (CONSENSUS) SDM OVER N RUNS ###############################

# Find the average across iterations
if (num_iterations > 1) {
  if (responce_variable_type == "binary") {
    pred.prob.LGM <- weighted.mean(pred.prob.LGM.stack, auc.stack, na.rm=TRUE)
    thres.mean <- mean(thres.stack)
    par(mfrow=c(1,2), mar = c(3, 3, 3, 3))
  } else if (responce_variable_type == "frequency") {
    pred.prob.LGM <- mean(pred.prob.LGM.stack, na.rm=TRUE)
    par(mfrow=c(1,1), mar = c(3, 3, 3, 3))
  }
}

# Set plotting parameters
if (responce_variable_type == "binary") {
  par(mfrow=c(1,2), mar = c(3, 3, 3, 3))
} else if (responce_variable_type == "frequency") {
  par(mfrow=c(1,1), mar = c(3, 3, 3, 3))
}

# Plot the ensemble model - for LGM
plot(pred.prob.LGM, col = brewer.pal(9, "YlOrBr"), main = paste("Continuous Spatial Projection (Ensemble Model -", cast.var.name, ",", num_iterations, "runs)", sep=" "),
     legend.args = list(text = "Predicted Probability", side = 2, font = 4, cex = 0.9, line = 0.5),
     legend.width = 2, legend.mar = 5, zlim = c(0, 1),
     legend.shrink = 0.9)
#plot(borders, add=TRUE,lwd=0.2,lty=1)
if (responce_variable_type == "binary") {
  pred.bin.LGM <- pred.prob.LGM > thres.mean
  plot(pred.bin.LGM, col = c("grey90", "green4"), main=paste("Binary Spatial Projection (Ensemble Model -", cast.var.name, ",", num_iterations, "runs)", sep=" "),
       legend.width = 2, legend.mar=5, legend=FALSE)
  #plot(borders, add=TRUE,lwd=0.2,lty=1)
  legend("bottomleft", fill = c("green4", "grey90"),
         legend = c("Presence", "Absence"), xpd = TRUE,bty = "n")
}

# To mask allele distribution models, we need the species distribution model. Hence run this first.
if (responce_variable_type == "binary") {
  pred.bin.LGM.species <- pred.bin.LGM
}

if (responce_variable_type == "frequency")  {
  if (exists("pred.bin.LGM.species") == TRUE) {
    # Constrain the allele projected distribution to where the species occurs (the binary species projection)
    if (cast == "hindcast") {
      # Convert to spatial object
      pred.bin.LGM.species <- pred.bin.LGM.species[pred.bin.LGM.species$cluster == 1,]
      coordinates(pred.bin.LGM.species) <- ~x+y
      proj4string(pred.bin.LGM.species) = CRS(prj_longlat)
    }
    # To mask projected CEN raster by the project species SDM raster, we use the mask function. We mask as such because the occurrence of the alleles is contingent on the occurence of the species.
    # E.g. if we defined the projected CEN raster as pred.bin.LGM and the species raster as pred.bin.LGM.species, we can mask as follows (for the binary):
    #pred.bin.LGM.masked<-mask(pred.bin.LGM, pred.bin.LGM.species, maskvalue=0, updatevalue = 0)
    # And similarly for the continious prediction
    if (cast == "hindcast") {
      pred.prob.LGM.masked<-mask(pred.prob.LGM, pred.bin.LGM.species, updatevalue = 0)
    } else if (cast == "forecast") {
      pred.prob.LGM.masked<-mask(pred.prob.LGM, pred.bin.LGM.species, maskvalue=0, updatevalue = 0)
    }
    #pred.prob.LGM.masked<-mask(pred.prob.LGM, pred.bin.current.species, maskvalue=0, updatevalue = 0)

    # We can then plot the masked prediction rasters
    proj4string(pred.prob.LGM.masked) = CRS(prj_longlat)
    zeroCol <-"grey90"
    #reds <- brewer.pal('YlOrRd', n = 9)
    myPalette<-brewer.pal(11,"RdYlBu")
    my.at=seq(0, 1, by=0.1)
    #my.brks=seq(0, 1, by=0.1)
    my.brks=seq(0.05, 1.05, by=0.1)
    myLabels <- c("NA (species not present)", "0.9 - low allele", "0.8 - low allele", "0.7 - low allele", "0.6 - low allele", "0.5 - high/low allele", "0.6 - high allele", "0.7 - high allele", "0.8 - high allele", "0.9 - high allele", "1.0 - high allele")
    #myLabels <- c("NA (species not present)", "0.1/0.9 - high/low allele", "0.2/0.8 - high/low allele", "0.3/0.7 - high/low allele", "0.4/0.6 - high/low allele", "0.5/0.5 - high/low allele", "0.6/0.4 - high/low allele", "0.7/0.3 - high/low allele", "0.8/0.2 - high/low allele", "0.9/0.1 - high/low allele", "1.0/0.0 - high/low allele")
    #myColorkey <- list(at=my.brks, labels=list(at=my.brks, labels=my.at))
    myColorkey <- list(at=my.brks, labels=list(at=my.at, labels=myLabels))
    #myTheme <- rasterTheme(region = c(zeroCol, reds))
    myTheme <- rasterTheme(region = c(zeroCol, myPalette))
    pred.prob.LGM.masked <- crop(pred.prob.LGM.masked, e_cen_plot)
    pred.prob.LGM.masked.plot <- levelplot(pred.prob.LGM.masked, par.settings = myTheme, at=my.at, colorkey=myColorkey, maxpixels = 1e7, margin = FALSE, main = paste0("Continuous Spatial Projection (Masked Ensemble Model - ", cast.var.name, ")"))
    pred.prob.LGM.masked.plot
  }
  pred.prob.LGM.masked.CentralAlps <-mask(pred.prob.LGM, CentralAlps_current_mask, maskvalue=0, updatevalue = 0)
  if (plot_all == TRUE) {
    # We can then plot the masked prediction raster (using Central alpine distribution mask)
    proj4string(pred.prob.LGM.masked.CentralAlps) = CRS(prj_longlat)
    pred.prob.LGM.masked.CentralAlps <- crop(pred.prob.LGM.masked.CentralAlps, e_cen_plot)
    pred.prob.LGM.masked.CentralAlps.plot <- levelplot(pred.prob.LGM.masked.CentralAlps, par.settings = myTheme, at=my.at, colorkey=myColorkey, maxpixels = 1e7, margin = FALSE, main = paste0("Continuous Spatial Projection (Masked Ensemble Model - ", cast.var.name, ")"))
    pred.prob.LGM.masked.CentralAlps.plot
  }
}

saveRDS(pred.prob.LGM, paste0("pred.prob.LGM.",lineage,".selectVar.optimalThreshold.",time_stamp,".rds"))
saveRDS(pred.bin.LGM, paste0("pred.bin.LGM.",lineage,".selectVar.optimalThreshold.",time_stamp,".rds"))
save(glm.multi.step, gam.step, rf.prob, maxent.prob, thres.glm.multi.step, thres.gam.step, thres.rf, thres.maxent, thres.mean, auc_glm, auc_gam, auc_maxent, auc_rf, auc_mean, file = paste0("colonisation_tracker_metaData_",lineage,"_",time_stamp,".RData"))

############################### CLUSTER ANALYSIS FOR REFUGIA ###############################

if (responce_variable_type == "binary" & cast == "hindcast") {
  # Perform cluster analysis on predicted binary projection LGM, to define refugia.
  # First convert the raster to a dataframe. Then select only presence points
  pred.bin.LGM.df <- as.data.frame(rasterToPoints(pred.bin.LGM))
  pred.bin.LGM.clusters.df <- pred.bin.LGM.df[pred.bin.LGM.df$layer == 1, ]
  pred.bin.LGM.clusters.df <- pred.bin.LGM.clusters.df[,c(1,2)]

  # Perform density-based spatial clustering (DBSCAN). eps defines the size of the epsilon neighborhood and minPts defines number of minimum points in the eps region (for core points).
  #DBSCAN <- dbscan(pred.bin.LGM.clusters.df, eps = 0.5, minPts = nrow(pred.bin.LGM.clusters.df)/100)
  DBSCAN <- dbscan(pred.bin.LGM.clusters.df, eps = 0.5, minPts = nrow(pred.bin.LGM.clusters.df)/30)

  # Convert to SpatialPointsDataframe
  pred.bin.LGM.clusters.df["cluster"] <- DBSCAN$cluster
  pred.bin.LGM.clusters <- pred.bin.LGM.clusters.df
  coordinates(pred.bin.LGM.clusters) <- ~x+y
  proj4string(pred.bin.LGM.clusters) = CRS(prj_longlat)
  proj4string(pred.bin.LGM) = CRS(prj_longlat)

  # And plot
  par(mfrow=c(1,1), mar = c(3, 3, 3, 3))
  plot(pred.bin.LGM, col = c("grey90", "grey90"), main="Predicted LGM Refugial Clusters (DBSCAN - Ensemble Model)", legend = FALSE)
  points(pred.bin.LGM.clusters, col = pred.bin.LGM.clusters$cluster, cex = 0.01)
  legend("bottomleft", col=unique(pred.bin.LGM.clusters$cluster), legend = unique(pred.bin.LGM.clusters$cluster), xpd = TRUE,bty = "n", pch=20, title = "Clusters")
}

saveRDS(pred.bin.LGM.clusters.df, paste0("pred.bin.LGM.clusters.df.",lineage,".selectVar.optimalThreshold.",time_stamp,".rds"))
#saveRDS(pred.bin.LGM, "pred.bin.LGM.plottingTemplate.rds")

############################### INVESTIGATE TEMPORAL TURNOVER ###############################

# Finally, we can visualise the difference in the contemporary and hindcasted projections
# We look at the turnover, i.e. the different kinds of change (gain, loss, no change) that we can expect. E.g. here for GAM model
# For continuous difference
if (responce_variable_type == "binary") {
  current.subtract.LGM <- pred.prob.current - pred.prob.LGM
  LGM.subtract.current <- pred.prob.LGM - pred.prob.current
} else if (responce_variable_type == "frequency") {
  current.subtract.LGM <- pred.prob.current.masked - pred.prob.LGM.masked
  LGM.subtract.current <- pred.prob.LGM.masked - pred.prob.current.masked
}

current.subtract.LGM[current.subtract.LGM < 0] <- NA
LGM.subtract.current[current.subtract.LGM < 0] <- NA
LGM.subtract.current.masked <- mask(LGM.subtract.current, current.subtract.LGM, inverse=TRUE)
current.subtract.LGM.masked <- mask(current.subtract.LGM, LGM.subtract.current)
if (responce_variable_type == "binary") {
  # For binary difference
  pred.bin.LGM.rcl <- pred.bin.LGM
  pred.bin.LGM.rcl[pred.bin.LGM > 0] <- 2
  turnover.bin <- pred.bin.current + pred.bin.LGM.rcl
  freq(turnover.bin)
  # We generate a turnover/range change table
  table(getValues(pred.bin.current),getValues(pred.bin.LGM.rcl))
}

# Plot the continous and binary turnover in suitable habitat between current and future for the ensemble model
# Make one raster all positive values and the other all negative values, so that we can combine the two into one plot (raster)
LGM.subtract.current.masked.negative <- 0 - LGM.subtract.current.masked
# Use the 'cover' function to overlay the two rasters together
current.continuous.turnover <- cover(LGM.subtract.current.masked.negative, current.subtract.LGM.masked)
proj4string(current.continuous.turnover) = CRS(prj_longlat)
# And plot
if (cast == "hindcast") {
  current.continuous.turnover.plot <- levelplot(current.continuous.turnover,par.settings=RdBuTheme(), at=seq(-max(abs(cellStats(current.continuous.turnover, range))), max(abs(cellStats(current.continuous.turnover, range))), len=100),  maxpixels = 1e7, margin = FALSE, main = expression("Expected Temporal Turnover (Continous)"), legend=list(top=list(fun=grid::textGrob("Turnover", y=1, x=1.06))))
} else if (cast == "forecast") {
  current.continuous.turnover.plot <- levelplot(current.continuous.turnover,par.settings=BuRdTheme(), at=seq(-max(abs(cellStats(current.continuous.turnover, range))), max(abs(cellStats(current.continuous.turnover, range))), len=100),  maxpixels = 1e7, margin = FALSE, main = expression("Expected Temporal Turnover (Continous)"), legend=list(top=list(fun=grid::textGrob("Turnover", y=1, x=1.06))))
}

if (responce_variable_type == "binary") {
  # A raster that contains categorical data can be defined with the ratify function
  turnover.bin.cat <- ratify(turnover.bin)
  proj4string(turnover.bin.cat) = CRS(prj_longlat)
  # The levels are stored in the “Raster Attribute Table” (RAT) that can be manipulated with the levels function:
  turnover.bin.cat.ID <- levels(turnover.bin.cat)[[1]]
  if (cast == "hindcast") {
    turnover.bin.cat.ID$State <- c('Absence', 'Colonisation', 'Contraction', 'Presence')
  } else if (cast == "forecast") {
    turnover.bin.cat.ID$State <- c('Absence', 'Contraction', 'Colonisation', 'Presence')
  }
  levels(turnover.bin.cat) <- turnover.bin.cat.ID
  # And plot
  if (cast == "hindcast") {
    current.bin.turnover.plot <- levelplot(turnover.bin.cat, att = "State", col.regions=c("#E5E4E2", "#3182bd", "#de2d26", "#85BB65"), maxpixels = 1e7, margin = FALSE, main = expression("Expected Temporal Turnover (Binary)"))
  } else if (cast == "forecast") {
    current.bin.turnover.plot <- levelplot(turnover.bin.cat, att = "State", col.regions=c("#E5E4E2", "#de2d26", "#3182bd", "#85BB65"), maxpixels = 1e7, margin = FALSE, main = expression("Expected Temporal Turnover (Binary)"))
  }
}

if (responce_variable_type == "binary") {
  grid.arrange(current.continuous.turnover.plot, current.bin.turnover.plot, ncol=2)
} else if (responce_variable_type == "frequency") {
  current.continuous.turnover.plot
}

############################### PLOT IN ENVIRONMENTAL SPACE ###############################

# Calculate PCA of predictor variables. We do this because for visualisation, as we are limited to 2 dimensions.

# Variable selection - NO NEED FOR THIS, because PCA already converts a set of observations of possibly correlated variables into a set of values of linearly uncorrelated variables, using orthogonal transformations.
#pred.var.selection <- c("TEMP_MEAN_ANNUAL","TEMP_ANNUAL_RANGE","PREC_ANNUAL","PH_5cm","SLOPE")
#ENV_LGM_RASTERS.selection <- subset(ENV_LGM_RASTERS, pred.var.selection)
#ENV_RASTERS.selection <- subset(ENV_RASTERS, pred.var.selection)

# Select number of principal components to plot
num_PCs <- 2

# First calculate for current climate raster. This code is taken from ENMtools 'raster.pca' function (type 'raster.pca' on the command prompt to see how this function is defined), modified to use the prcomp function rather than the princomp function.
env.val <- getValues(ENV_RASTERS)
env.val.keepers <- which(complete.cases(env.val))
env.val.nas <- which(!complete.cases(env.val))
env.pca <- prcomp(env.val[env.val.keepers, ], scale. = TRUE, center = TRUE)
ENV_PCA_TODAY <- ENV_RASTERS[[1:num_PCs]]
ENV_PCA_TODAY[env.val.nas] <- NA
for (i in 1:num_PCs) {
  ENV_PCA_TODAY[[i]][env.val.keepers] <- env.pca$x[, i]
}
names(ENV_PCA_TODAY) <- paste0("PC", 1:num_PCs)
ENV_PCA_TODAY <- setMinMax(ENV_PCA_TODAY)

# Then calculate the PCA of the LGM environment raster. Note that here we project this onto PCA of the current climate (i.e. that this LGM PCA will be transformed, scaled and centered based on the current climate PCA), for apples-to-apples comparison
env.val.LGM <- getValues(ENV_LGM_RASTERS)
env.val.keepers.LGM <- which(complete.cases(env.val.LGM))
env.val.nas.LGM <- which(!complete.cases(env.val.LGM))
# Here we transform the data by applying the scaling and centering of the current climate PCA. This code was modified from Simon's mapPCA.r function.
env.val.LGM.tranformed <- scale(env.val.LGM[env.val.keepers.LGM, ][, rownames(env.pca$rotation)], scale = FALSE, center = env.pca$center)
env.val.LGM.tranformed <- scale(env.val.LGM.tranformed, scale = env.pca$scale, center = FALSE)
# Then map this (LGM) data (transformed as the current climate data) in the PCA space of the current climate data
env.LGM.pca <-  as.matrix(env.val.LGM.tranformed) %*% as.matrix(env.pca$rotation)
# The build the raster
ENV_PCA_LGM <- ENV_LGM_RASTERS[[1:num_PCs]]
ENV_PCA_LGM[env.val.nas.LGM] <- NA
for (i in 1:num_PCs) {
  ENV_PCA_LGM[[i]][env.val.keepers.LGM] <- env.LGM.pca[, i]
}
names(ENV_PCA_LGM) <- paste0("PC", 1:num_PCs)
ENV_PCA_LGM <- setMinMax(ENV_PCA_LGM)

# Let's also get the eigenvalues for these PCAs
PCA.eigenvalues.today <- env.pca$sdev^2/sum(env.pca$sdev^2)*100
PCA.eigenvalues.today <- PCA.eigenvalues.today[seq(1,2)]
env.LGM.pca.df <- as.data.frame(env.LGM.pca)
env.LGM.pca.sd <- sapply(env.LGM.pca.df, sd)
PCA.eigenvalues.LGM <- env.LGM.pca.sd^2/sum(env.LGM.pca.sd^2)*100
PCA.eigenvalues.LGM <- PCA.eigenvalues.LGM[seq(1,2)]

# Extract predictor PCA components at occurrence sites
Occurrences.Climate.absences <- Occurrences.Climate[Occurrences.Climate$Presence == 0,]
if (responce_variable_type == "binary") {
  Occurrences.Climate.presences <- Occurrences.Climate[Occurrences.Climate$Presence == 1,]
} else if (responce_variable_type == "frequency") {
  Occurrences.Climate.presences <- Occurrences.Climate[Occurrences.Climate$Presence > 0,]
}
Data.ENV_PCA.presences <- na.omit(as.data.frame(extract(ENV_PCA_TODAY, Occurrences.Climate.presences[,c("x","y")])))
Data.ENV_PCA.absences <- na.omit(as.data.frame(extract(ENV_PCA_TODAY, Occurrences.Climate.absences[,c("x","y")])))

# Extract predictor PCA components at all sites (aka convert full raster to points)
PCA1_ENV_TODAY <- as.data.frame(rasterToPoints(ENV_PCA_TODAY[[1]]))
PCA2_ENV_TODAY <- as.data.frame(rasterToPoints(ENV_PCA_TODAY[[2]]))
PCA1_ENV_LGM <- as.data.frame(rasterToPoints(ENV_PCA_LGM[[1]]))
PCA2_ENV_LGM <- as.data.frame(rasterToPoints(ENV_PCA_LGM[[2]]))
# Combine the PC1 and PC2 dataframes
ENV_PCA1PCA2_TODAY <- cbind(PCA1_ENV_TODAY, PCA2_ENV_TODAY[3])
ENV_PCA1PCA2_LGM <- cbind(PCA1_ENV_LGM, PCA2_ENV_LGM[3])
# We round the numeric vector because for merging dataframes later, we want to avoid floating number artifacts (see: https://stackoverflow.com/questions/9508518/why-are-these-numbers-not-equal) and (https://stackoverflow.com/questions/40183163/merging-with-all-equal)
ENV_PCA1PCA2_TODAY$x=round(ENV_PCA1PCA2_TODAY$x,6)
ENV_PCA1PCA2_TODAY$y=round(ENV_PCA1PCA2_TODAY$y,6)
ENV_PCA1PCA2_LGM$x=round(ENV_PCA1PCA2_LGM$x,6)
ENV_PCA1PCA2_LGM$y=round(ENV_PCA1PCA2_LGM$y,6)

# We want to colour the above PCA in environmental space by the suitabilty of each cell. We get this from the predicted SDM model. We again round the coordinates to allow merging dataframes later, we want to avoid floating number artifacts
pred.prob.current.df <- as.data.frame(rasterToPoints(pred.prob.current))
pred.prob.LGM.df <- as.data.frame(rasterToPoints(pred.prob.LGM))
pred.prob.current.df$x=round(pred.prob.current.df$x,6)
pred.prob.current.df$y=round(pred.prob.current.df$y,6)
pred.prob.LGM.df$x=round(pred.prob.LGM.df$x,6)
pred.prob.LGM.df$y=round(pred.prob.LGM.df$y,6)

# Finally we merge the PCA components to the predicted probability of occurrence (e.g. suitability)
ENV_PCA1PCA2_Pred_TODAY <- merge(ENV_PCA1PCA2_TODAY, pred.prob.current.df, by.x=c("x", "y"), by.y=c("x", "y"))
names(ENV_PCA1PCA2_Pred_TODAY)[names(ENV_PCA1PCA2_Pred_TODAY) == "layer"] <- "presence_prob"
ENV_PCA1PCA2_Pred_LGM <- merge(ENV_PCA1PCA2_LGM, pred.prob.LGM.df, by.x=c("x", "y"), by.y=c("x", "y"))
names(ENV_PCA1PCA2_Pred_LGM)[names(ENV_PCA1PCA2_Pred_LGM) == "layer"] <- "presence_prob"

# Take a random subset of dataframe for plotting
ENV_PCA1PCA2_Pred_TODAY <- ENV_PCA1PCA2_Pred_TODAY[sample(nrow(ENV_PCA1PCA2_Pred_TODAY), 300000), ]
ENV_PCA1PCA2_Pred_LGM <- ENV_PCA1PCA2_Pred_LGM[sample(nrow(ENV_PCA1PCA2_Pred_LGM), 300000), ]

# To produce biplot (add variable vectors to final PCA), we extract variable vectors from output of prcomp
# Code for biplot variable vectors taken from source code of: getAnywhere(biplot.prcomp) and getAnywhere(biplot.default)
# For reference, see: https://stat.ethz.ch/R-manual/R-devel/library/stats/html/biplot.princomp.html
# Define scale (scale has to be between 0 and 1)
scale = 1
choices = 1L:2L
scores <- env.pca$x
lam <- env.pca$sdev[choices]
n_scores <- NROW(scores)
lam <- lam * sqrt(n_scores)
if (scale != 0) {
  lam <- lam^scale
} else {
  lam <- 1
}

# Apply scaling and produce the observations and variables matrices for biplot.default
observations_matrix <- t(t(scores[, choices])/lam)
variables_matrix <- t(t(env.pca$rotation[, choices]) * lam)

# Plot the presence points and the suitability of all cells in environmental space (PC1 vs PC2)

title.today <- paste("Projected (Today, Ensemble) model in native (Today) environment space - PC1"," (",signif(PCA.eigenvalues.today[1], digits=3),"%)"," / PC2"," (",signif(PCA.eigenvalues.today[2], digits=3),"%)",sep="",collapse="")
title.LGM <- paste("Projected (LGM, Ensemble) model in projected (Today) environment space - PC1"," (",signif(PCA.eigenvalues.LGM[1], digits=3),"%)"," / PC2"," (",signif(PCA.eigenvalues.LGM[2], digits=3),"%)",sep="",collapse="")
# In case warning message pops up stating removal of rows, reduce the arrow_scale to a lower value (error is because the arrows and labels lie outside plotting range)
arrow_scale <- 60

ENV.space.today <- ggplot() +
  geom_point(data = ENV_PCA1PCA2_Pred_TODAY, aes(x = PC1, y = PC2, col = presence_prob)) +
  geom_point(data = Data.ENV_PCA.presences, aes(x = PC1, y = PC2, fill = "green4"), shape = 21) +
  geom_point(data = Data.ENV_PCA.absences, aes(x = PC1, y = PC2, fill = "red"), shape = 21) +
  scale_colour_distiller(palette = "YlOrBr", direction = 1, name = "Presence probabilty") + scale_fill_manual(values = c("green4", "red"), labels = c("Presence", "Absence"), name = "Observed data") +
  geom_segment(aes(x = 0, y = 0, xend = PC1/arrow_scale, yend = PC2/arrow_scale), data = as.data.frame(variables_matrix), arrow=arrow(length=unit(0.5,"cm"))) +
  geom_text(data=as.data.frame(variables_matrix), mapping=aes(x=PC1/arrow_scale, y=PC2/arrow_scale, label=rownames(variables_matrix)), size=4) +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_hline(yintercept=0) + geom_vline(xintercept = 0) +
  xlim(-25, 25) + ylim(-11,6) + xlab('PC1') + ylab('PC2') + ggtitle(title.today)

ENV.space.LGM <- ggplot() +
  geom_point(data = ENV_PCA1PCA2_Pred_LGM, aes(x = PC1, y = PC2, col = presence_prob)) +
  scale_colour_distiller(palette = "YlOrBr", direction = 1, name = "Presence probabilty") +
  geom_segment(aes(x = 0, y = 0, xend = PC1/arrow_scale, yend = PC2/arrow_scale), data = as.data.frame(variables_matrix), arrow=arrow(length=unit(0.5,"cm"))) +
  geom_text(data=as.data.frame(variables_matrix), mapping=aes(x=PC1/arrow_scale, y=PC2/arrow_scale, label=rownames(variables_matrix)), size=4) +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_hline(yintercept=0) + geom_vline(xintercept = 0) +
  xlim(-25, 25) + ylim(-11,6) + xlab('PC1') + ylab('PC2') + ggtitle(title.LGM)

grid.arrange(ENV.space.today, ENV.space.LGM, ncol=2)

############################### FOR INTERACTIVE PLOTTING OF MAPS USING LEAFLET ###############################

# library(leaflet)
#
# # Define colour palette (to see available colur palettes, see: http://leaflet-extras.github.io/leaflet-providers/preview/index.html)
# pal <- colorNumeric(palette = "YlOrBr", values(pred.prob.current),na.color = "transparent")
#
# # Define CRS
# proj4string(pred.prob.current) = CRS(prj_longlat)
#
# # Generate interactive leaflet map (see: https://rstudio.github.io/leaflet/raster.html)
# m <- leaflet() %>% addTiles() %>%
#   addProviderTiles(providers$Esri.WorldPhysical) %>%
#   addRasterImage(pred.prob.current, colors = brewer.pal(9, "YlOrBr"), opacity = 0.8) %>%
#   addLegend(pal = pal, values = values(pred.prob.current), title = "Probability of presence") %>%
#   fitBounds(3, 38, 21, 48)
# #setMaxBounds(3, 38, 21, 48)
#
# # View plot
# m
