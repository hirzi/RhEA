##### This script performs gradient forest modelling to investigate composition of whole-genome adaptive genetic variation across space and time (current and LGM) #####

# Note 1: You'll probably want to mask to the rasters to their predicted SDM distributions. 
# Note 2: Consider if you want to remove SNPs that are polymorphic in <5 populations (to ensure robust regression).
# Note 3: Consider if it makes sense to control for neutral structure by running GF analysis first on a putative neutral set of SNPs, and then removing the effect of this (e.g. rotated, rescaled, and finding the difference) from the GF adaptive set.
# Note 4: Think about the correct interpretation of the multi-dimensional Euclidean distance of transformed environmental predictors between present and LGM - e.g. is this the difference in aggregate allele frequencies? the difference in composition of SNP allele frequencies?

# Load libraries
library(raster)
library(gradientForest)
library(gtools)
library(FNN)
library(concaveman)
library(gtools)

# Define run parameters
working_dir <- "/Users/luqman/Documents/Hirzi/ETHPHD/PopulationGenetics_wholeGenome/Variant calling for Central Alps/"
#working_dir <- "/cluster/work/gdc/shared/p461/secondNGSdataset/Population_Genetics/GF_GDM/"
setwd(working_dir)
# Define gradientForest run parameters
num_SNPs_perRun <- 10000
number_total_SNPs <- 390262
# Define/explore results based on differing LGM refugial extents ("big" or "small")
LGM_mask_extent <- "small"
# Define hull method for masks
hull_method <- "concave"
# Define plotting parameters
plot_figures <- "PARTIAL" # ALL, PARTIAL, NONE
# Default (r=PC1, g=PC2, b=PC3) or alternative (r=PC1 + PC2, g=-PC2, b=PC3 + PC2 - PC1)
col_palette <- "default"
## Import population coordinate data
#pop_geodata <- read.csv(paste0(working_dir,"popInfo_w14indsCorrectHe_NEW.csv"))
pop_geodata <- read.csv(paste0(working_dir,"popInfo_w14indsCorrectHe_NEW_RevisedR2ENV.csv"))
pop_geodata_GF <- pop_geodata
# Define whether object is of class gradientForest or combinedGradientForest, i.e. whether you want to run gradientForest (FALSE) or whether you want to use an imported a combinedGradientForest object (TRUE)
combined_gradientForest <- TRUE
# Import elevation base map
#ELEV.GMTED <- raster("/Users/luqman/Desktop/SDM/Dsylvestris/Data/GMTED2010N30E000_075/30n000e_20101117_gmted_mea075.tif")
#ELEV.GMTED <- raster("/Users/luqman/Desktop/SDM/Dsylvestris/Data/GMTED2010N30E000_300/30n000e_20101117_gmted_mea300.tif")
ELEV_RASTERS <-stack("/Users/luqman/Desktop/SDM/Dsylvestris/Data/ENV_RASTERS_TODAY_EUROPE_CHELSA_LARGE.grd") 
ELEV_RASTERS_LGM <- stack("/Users/luqman/Desktop/SDM/Dsylvestris/Data/ENV_RASTERS_LGM_EUROPE_PMIP_ENSEMBLE_NCAR_MIROC_MRI_MPI_LARGE.grd")
ELEV.alps <- ELEV_RASTERS$ELEV.GMTED
ELEV.LGM.alps <- ELEV_RASTERS_LGM$ELEV.GMTED
prj_longlat <- "+init=epsg:4326"
crs(ELEV.alps) <- crs(prj_longlat)
crs(ELEV.LGM.alps) <- crs(prj_longlat)
#e_alps <- extent(5,13,43.5,47.5)
e_alps <- extent(4.75,13.3,43.75,47.9)
ELEV.alps <- crop(ELEV.alps, e_alps)
ELEV.LGM.alps <- crop(ELEV.LGM.alps, e_alps)

# If running gradientForest:
if (combined_gradientForest == FALSE) { 
  ## Import allele frequency data
  allele_data <- read.table(paste0(working_dir,"AllPops_innerMerged.minQ30minDP1204minAC14miss10_moreStringent.prim.SNP.exons_minDepth7.AF.AFcolumnsOnly"))
  # Transpose
  allele_data <- as.data.frame(t(allele_data))
  # Convert first row to header
  names(allele_data) <- as.character(unlist(allele_data[1,]))
  allele_data <- allele_data[-1,]
  colnames(allele_data)[1] <- "pop"
  # Rename pops (remove "AF_" prefix)
  allele_data$pop <- gsub('AF_', '', allele_data$pop)
  
  ## Merge population coordinate and SNP allele frequency data
  allele_data <- merge(pop_geodata, allele_data, by.x="pop", by.y="pop", all.x=TRUE)
  # Convert first column to rownames
  allele_data_temp <- allele_data[,-1]
  rownames(allele_data_temp) <- allele_data[,1]
  allele_data <- allele_data_temp
  rm(allele_data_temp)
  # Convert all classes to numeric
  allele_data[] <- lapply(allele_data, function(x) {
    if(is.factor(x)) as.numeric(as.character(x)) else x
  })
  #sapply(allele_data, class)
}

## Import current environmental data
CLIM.alps <-stack(paste0(working_dir,"ENV_RASTERS_TODAY_EUROPE_CHELSA_LARGE.grd"))
# Choose predictor variables. 
# Here we take uncorrelated variables (< 0.7) assessed at a larger extent of occurrences (SDM occurrences) and evaluated using pairs.hist.cor(Occurrences.Climate[, 4:28], cor.method = "pearson"), and round(cor(Occurrences.Climate[, 4:28]),2); varCor <- cor(Occurrences.Climate[, 4:28]); allDistNew <- abs(as.dist(cor(Occurrences.Climate[, 4:28]))); allClusNew <- hclust(1 - allDistNew); plot(allClusNew, hang=-1)
# To choose among the correlated variables, we refer to the importance of each variable evaluated from full (all 27 predictor variables) GDM and GF models (preliminary analyses). In case there is disagreement, we choose the optimum (most informative) compromise. For more details, see Predictor variable selection.xlsx
# We also assess the variance inflation factors (VIF) of the variables (conditional on VIF < 10).
#pred.var.selection <- c("TEMP_DIURNAL_RANGE", "TEMP_SEASONALITY", "TEMP_MIN_COLDEST_MONTH", "TEMP_MEAN_WETTEST_QUARTER", "TEMP_MEAN_DRIEST_QUARTER", "PREC_SEASONALITY", "PREC_WARMEST_QUARTER", "PREC_COLDEST_QUARTER", "PH_5cm", "SLOPE", "SOLAR_RAD", "LONGITUDE", "LATITUDE")
pred.var.selection <- c("TEMP_DIURNAL_RANGE", "TEMP_SEASONALITY", "TEMP_MIN_COLDEST_MONTH", "TEMP_MEAN_WETTEST_QUARTER", "TEMP_MEAN_DRIEST_QUARTER", "PREC_SEASONALITY", "PREC_WARMEST_QUARTER", "PREC_COLDEST_QUARTER", "PH_5cm", "SLOPE", "LONGITUDE", "LATITUDE")
CLIM.alps <- subset(CLIM.alps, pred.var.selection)
# And crop to desired geographic extent
e <- extent(4, 15, 43, 48)
CLIM.alps <- crop(CLIM.alps, e)
CLIM.alps.gf <- CLIM.alps
# Similarly, import LGM environmental data
CLIM.alps.LGM <-stack(paste0(working_dir,"ENV_RASTERS_LGM_EUROPE_PMIP_ENSEMBLE_NCAR_MIROC_MRI_MPI_LARGE.grd"))
CLIM.alps.LGM <- subset(CLIM.alps.LGM, pred.var.selection)
CLIM.alps.LGM <- crop(CLIM.alps.LGM, e)

# Import masks for contraining predictions to SDM extent.
# In the end, we'll start from the output of the SDM_Colonisation_Tracker, specifically we require the first and last time steps. We import this.
# To explore the dependence on LGM geographic extent, we choose between preliminary large and small LGM inferred refugial extents.
if (LGM_mask_extent == "big") {
  current.mask <- readRDS("/Users/luqman/Desktop/SDM/Dsylvestris/Data/currentMask.CentralAlps.dummy.rds")
  current.mask.original <- current.mask
  current.mask.original <- crop(current.mask.original, e)
  current.mask[current.mask > 0] <- NA
  LGM.mask <- readRDS("/Users/luqman/Desktop/SDM/Dsylvestris/Data/LGMMask.CentralAlps.dummy3.rds")
  LGM.mask.original <- LGM.mask
  LGM.mask.original[LGM.mask.original == 2] <- 1
  LGM.mask.original <- crop(LGM.mask.original, e)
  LGM.mask[LGM.mask > 0] <- NA
} else if (LGM_mask_extent == "small") {
  # Import time series SDMs (out of SDM_Colonisation_Tracker)
  cluster.timeSeries.rasters <- stack("/Users/luqman/Desktop/SDM/Dsylvestris/SDM_Euler/TIMESERIES.DF.ANIMATION.TIMEINTERVALS211.DISPERSAL0.2.3LineagesSplit.BAC.grd")
  # Extract the first (LGM) and last (current) time points.
  current.SDMtimeSeries.CentralAlps.rast <- cluster.timeSeries.rasters[[211]]
  current.SDMtimeSeries.CentralAlps.rast[current.SDMtimeSeries.CentralAlps.rast < 5002] <- 0
  current.mask.original <- current.SDMtimeSeries.CentralAlps.rast
  current.mask.original[current.mask.original >= 5002] <- 1
  # Additional correction
  current.mask.original.temp <- crop(current.mask.original, extent(15.7,16.8,43.8,44.8))
  current.mask.original.temp[current.mask.original.temp > 0] <- 0
  current.mask.original.temp1 <- crop(current.mask.original, extent(3,4.95,42,48))
  current.mask.original.temp1[current.mask.original.temp1 > 0] <- 0
  current.mask.original.temp2 <- crop(current.mask.original, extent(7,12,47.77,48.5))
  current.mask.original.temp2[current.mask.original.temp2 > 0] <- 0
  current.mask.original.temp3 <- crop(current.mask.original, extent(5,7,47.4,48.5))
  current.mask.original.temp3[current.mask.original.temp3 > 0] <- 0
  current.mask.original.temp4 <- crop(current.mask.original, extent(13.2,15,46,48))
  current.mask.original.temp4[current.mask.original.temp4 > 0] <- 0
  current.mask.original.temp5 <- crop(current.mask.original, extent(7.9,9,43.9,45))
  current.mask.original.temp5[current.mask.original.temp5 > 0] <- 0
  current.mask.original <- merge(current.mask.original.temp, current.mask.original)
  current.mask.original <- merge(current.mask.original.temp1, current.mask.original)
  current.mask.original <- merge(current.mask.original.temp2, current.mask.original)
  current.mask.original <- merge(current.mask.original.temp3, current.mask.original)
  current.mask.original <- merge(current.mask.original.temp4, current.mask.original)
  current.mask.original <- merge(current.mask.original.temp5, current.mask.original)
  current.mask.original <- crop(current.mask.original, e)
  # Extract Central Alpine lineage. See SDM_Adaptation_Prediction_LineageSpecific.R for cluster extraction based on raster value.
  current.SDMtimeSeries.CentralAlps.rast[current.SDMtimeSeries.CentralAlps.rast >= 5002] <- NA
  #plot(current.SDMtimeSeries.CentralAlps.rast)
  current.mask <- current.SDMtimeSeries.CentralAlps.rast
  # And similarly for the LGM
  # Import LGM mask clustered and filtered using DBSCAN to retain only eastern and dense/contiguous refugial clusters, given results from psi analysis and to remove sparse, disconnected presence points.
  #pred.bin.LGM.mask.clusters.raster <- readRDS("/Users/luqman/Desktop/SDM/Dsylvestris/Data/pred.bin.LGM.mask.clusters.raster.2021-08-28_14_21_31.rds") # larger
  pred.bin.LGM.mask.clusters.raster <- readRDS("/Users/luqman/Desktop/SDM/Dsylvestris/Data/pred.bin.LGM.mask.clusters.raster.2021-08-28_12_08_53.rds") # smaller, eastern clusters
  LGM.SDMtimeSeries.CentralAlps.rast.original <- cluster.timeSeries.rasters[[1]]
  LGM.SDMtimeSeries.CentralAlps.rast <- LGM.SDMtimeSeries.CentralAlps.rast.original
  # For genomic glacial offset mask
  LGM.SDMtimeSeries.CentralAlps.rast[LGM.SDMtimeSeries.CentralAlps.rast < 5002] <- 0
  LGM.mask.original <- LGM.SDMtimeSeries.CentralAlps.rast
  LGM.mask.original[LGM.mask.original >= 5002] <- 1
  LGM.mask.original <- crop(LGM.mask.original, e)
  LGM.SDMtimeSeries.CentralAlps.rast[LGM.SDMtimeSeries.CentralAlps.rast >= 5002] <- NA
  #plot(LGM.SDMtimeSeries.CentralAlps.rast)
  LGM.mask <- LGM.SDMtimeSeries.CentralAlps.rast
  # For plotting mask
  LGM.SDMtimeSeries.CentralAlps.rast.plot <- LGM.SDMtimeSeries.CentralAlps.rast.original
  LGM.SDMtimeSeries.CentralAlps.rast.plot[LGM.SDMtimeSeries.CentralAlps.rast.plot < 9999] <- 0
  LGM.SDMtimeSeries.CentralAlps.rast.cropped <- crop(LGM.SDMtimeSeries.CentralAlps.rast.plot, e)
  pred.bin.LGM.mask.clusters.raster.extended <- extend(pred.bin.LGM.mask.clusters.raster, e, value=NA)
  LGM.mask.original.plot <- mask(LGM.SDMtimeSeries.CentralAlps.rast.cropped, pred.bin.LGM.mask.clusters.raster.extended, maskvalue=1, updatevalue = 1)
  pred.bin.LGM.mask.clusters.raster.extended2 <- extend(pred.bin.LGM.mask.clusters.raster, extent(LGM.SDMtimeSeries.CentralAlps.rast.plot), value=NA)
  LGM.mask.plot <- mask(LGM.SDMtimeSeries.CentralAlps.rast.plot, pred.bin.LGM.mask.clusters.raster.extended2, maskvalue=1, updatevalue = NA)
}

# Should you wish to regenerate dbscaned LGM mask
dbscan_LGM_mask <- FALSE

if (combined_gradientForest == FALSE) {
  
  ## Combine environmental data and allele frequency data
  CLIM.extract <- data.frame(extract(CLIM.alps, allele_data[,c("x","y")]))
  CLIM_allele <- data.frame(CLIM.extract, allele_data) # let's use lat, lon as a proxy for PCA1, PCA2 (as correlation between these are high)
  preds <- c(names(CLIM.extract))
  loci <- names(allele_data[,-c(1,2)])

# To speed up analysis, we run gradientForest in parallel subsets (e.g. on the cluster)
  idx_start <- ((run_number-1)*num_SNPs_perRun) + 1
  idx_end <- run_number*num_SNPs_perRun
  if (run_number < ceiling(number_total_SNPs/num_SNPs_perRun)) {
    loci <- loci[idx_start:idx_end]
  } else if (run_number == ceiling(number_total_SNPs/num_SNPs_perRun)) {
    loci <- loci[idx_start:number_total_SNPs]
  }

  ## Run gradientForest
  gfVars <- gradientForest(CLIM_allele,
                           predictor.vars=preds,
                           response.vars=loci,
                           ntree=500,
                           maxLevel=floor(log2(0.368*nrow(CLIM_allele)/2)),
                           trace=T,
                           corr.threshold=0.50) # or 0.7
  # Often a warning will prompt. Type 'warnings()' to see them. Warning refers to whether you want to use regression on columns with fewer than n unique values.
  
  # Save output
  #time_stamp <- gsub(":","_",gsub(" ", "_", Sys.time()))"
  #saveRDS(gfVars, paste0("gfModel_selectedPredVars.split", run_number,".rds"))
  #gfVars <- readRDS(paste0("gfModel_allPredVars.",time_stamp,".rds"))
}


# ## Combine parallel (subset) runs of gradientForest into one gradientForest model (object).
# # This function is very resource heavy (>100GB memory for this job) so run on Euler!
#  
# # List gradientForest parallel runs
# #files <- list.files(path = working_dir, pattern = "^gfModel_allPredVars.split.*\\.10000SNPsets.rds$")
# files <- list.files(path = working_dir, pattern = "^gfModel_selectedPredVars.split.*\\.rds$")
# files <- mixedsort(files)
# 
# # Read results into R
# gfModel_subset_list <- lapply(files, readRDS)
# 
# # Combine parallel (subset) runs into an aggregate gradientForest model (object)
# gfModel_combined <- combinedGradientForest(gfModel_subset_list[[1]],gfModel_subset_list[[2]],gfModel_subset_list[[3]],gfModel_subset_list[[4]],gfModel_subset_list[[5]],gfModel_subset_list[[6]],gfModel_subset_list[[7]],gfModel_subset_list[[8]],gfModel_subset_list[[9]],gfModel_subset_list[[10]],gfModel_subset_list[[11]],gfModel_subset_list[[12]],gfModel_subset_list[[13]],gfModel_subset_list[[14]],gfModel_subset_list[[15]],gfModel_subset_list[[16]],gfModel_subset_list[[17]],gfModel_subset_list[[18]],gfModel_subset_list[[19]],gfModel_subset_list[[20]],gfModel_subset_list[[21]],gfModel_subset_list[[22]],gfModel_subset_list[[23]],gfModel_subset_list[[24]],gfModel_subset_list[[25]],gfModel_subset_list[[26]],gfModel_subset_list[[27]],gfModel_subset_list[[28]],gfModel_subset_list[[29]],gfModel_subset_list[[30]],gfModel_subset_list[[31]],gfModel_subset_list[[32]],gfModel_subset_list[[33]],gfModel_subset_list[[34]],gfModel_subset_list[[35]],gfModel_subset_list[[36]],gfModel_subset_list[[37]],gfModel_subset_list[[38]],gfModel_subset_list[[39]],gfModel_subset_list[[40]])
# 
# # Save object
# #saveRDS(gfModel_combined, "gfModel_allPredVars_combined_allSplits.rds")
# saveRDS(gfModel_combined, "gfModel_selectedPredVars_combined_allSplits.rds")

#if (combined_gradientForest == TRUE) { 
# Choose default options (method=2, standardize="before"). Choose between cor=0.5,0.7 - results appear very similar.
gfVars <- readRDS(paste0(working_dir,"gfModel_cor0.5_selectedPredVars_woSolarRad_combined_method2standardizeBefore_allSplits.rds"))
#gfVars <- readRDS(paste0(working_dir,"gfModel_cor0.5_selectedPredVars_combined_method2standardizeBefore_allSplits.rds"))
#gfVars <- readRDS(paste0(working_dir,"gfModel_cor0.5_selectedPredVars_combined_method2standardizeAfter_allSplits.rds"))
#gfVars <- readRDS(paste0(working_dir,"gfModel_cor0.7_selectedPredVars_combined_method2standardizeBefore_allSplits.rds"))
#gfVars <- readRDS(paste0(working_dir,"gfModel_cor0.7_selectedPredVars_combined_method2standardizeAfter_allSplits.rds"))
#gfVars <- readRDS(paste0(working_dir,"gfModel_cor0.7_allPredVars_combined_method2standardizeBefore_allSplits.rds"))

# Print results
# Predictor importance (recall that the accuracy important is how well each variable is at predicting the response)
vars <- names(importance(gfVars))
if (plot_figures == "ALL") {
  importance(gfVars)
  plot(gfVars, plot.type="O")
}

# Make df of variable importance
var_importance <- gfVars$imp.rsq
var_importance <- var_importance[,-1]
var_importance.df <- as.data.frame(rowSums(var_importance))
var_importance.df[,2] <- 0
var_importance.sorted.df <- var_importance.df[rev(order(var_importance.df[1])),]
predVar_colSum <- colSums(var_importance.sorted.df)[1]
var_importance.sorted.df[,2] <- var_importance.sorted.df[,1]/predVar_colSum
colnames(var_importance.sorted.df) <- c("importance", "normalized_importance")

if (plot_figures == "ALL") {
  # Measure of model fit (R2) by SNP (gradientForest) or SNP set (combinedGradientForest)
  if (combined_gradientForest == FALSE) {
    plot(gfVars, plot.type="P", show.names=T, horizontal=F)
  } else if (combined_gradientForest == TRUE) {
    #plot(gfVars, plot.type="Performance") # default plot
    # Extract total-R2 and ENV(only)-R2 from gfvars. The former can be found directly at gfVars$rsq, however the latter needs to be calculated.
    r2_allpreds <- as.data.frame(t(gfVars$imp.rsq))[-c(1),]
    r2_allpreds[] <- lapply(r2_allpreds, function(x) {
      if(is.factor(x)) as.numeric(as.character(x)) else x
    })
    #sapply(r2_allpreds, class)
    # Calculate total-R2 and ENV(only)-R2
    r2_allpreds[,ncol(r2_allpreds)+1] <- rowSums(r2_allpreds[,1:10]) # without SolarRad
    #r2_allpreds[,ncol(r2_allpreds)+1] <- rowSums(r2_allpreds[,1:11])
    colnames(r2_allpreds)[ncol(r2_allpreds)] <- "r2_ENV"
    r2_allpreds[,ncol(r2_allpreds)+1] <- rowSums(r2_allpreds[,1:12])
    #r2_allpreds[,ncol(r2_allpreds)+1] <- rowSums(r2_allpreds[,1:13]) # without SolarRad              
    colnames(r2_allpreds)[ncol(r2_allpreds)] <- "r2_total"
    # Format data for boxplot
    #r2_factors <- split(gfVars$rsq,factor(rep(names(gfVars$nspec),gfVars$nspec)))
    #r2_factors <- split(r2_allpreds$r2_total,factor(rep(names(gfVars$nspec),gfVars$nspec))) # for total R2, identical to above
    r2_factors <- split(r2_allpreds$r2_ENV,factor(rep(names(gfVars$nspec),gfVars$nspec))) # for env-only R2 (i.e. excluding latlon)
    #r2_factors.sorted <- r2_factors[mixedsort(names(r2_factors))] # for horizontal
    r2_factors.sorted <- r2_factors[rev(mixedsort(names(r2_factors)))] # for vertical
    #boxplot(r2_factors.sorted,plot=TRUE, col=rgb(0.5,0.65,0.85,0.85), pars=list(outpch = 21, outcol="grey10", outbg = "brown1", outcex=1), las=2, ylab="R2", xlab="gradientForest runs") # nicer plot
    boxplot(r2_factors.sorted,plot=TRUE, col=rgb(0.5,0.65,0.85,0.85), pars=list(outpch = 21, outcol="grey10", outbg = "brown1", outcex=1.65), las=2, xlab="R2", ylab="gradientForest runs", horizontal = TRUE, ylim = c(0, 0.9)) # nicer plot
  }
  
  # Plot split densities of predictor variables
  if (combined_gradientForest == FALSE) {
    plot(gfVars, plot.type="S")
  } else if (combined_gradientForest == TRUE) {
    plot(gfVars, plot.type="Predictor.Density")
  }

  # Cumulative importance plot
  if (combined_gradientForest == FALSE) {
    # Species-level  
    plot(gfVars, plot.type="C", show.overall=F, leg.nspecies=5, cex.lab=1, common.scale=T, par.args=list(mgp=c(1.5,0.5,0), mar=c(2.5,1,0.1,0.5), omi=c(0,0.3,0,0)))
    # Community-level
    plot(gfVars, plot.type="C", show.species=F, leg.nspecies=5, cex.lab=1, common.scale=T, par.args=list(mgp=c(1.5,0.5,0), mar=c(2.5,1,0.1,0.5), omi=c(0,0.3,0,0)))
  } else if (combined_gradientForest == TRUE) {
    # Cumulative importance plot - species-level  
    plot(gfVars, plot.type="C", cex.lab=1, common.scale=T, par.args=list(mgp=c(1.5,0.5,0), mar=c(2.5,1,0.1,0.5), omi=c(0,0.3,0,0)))
  }
}

## Project model in space
CLIM.alps.df <- as.data.frame(CLIM.alps, xy=TRUE, na.rm=TRUE)
CLIM.alps.df.GF <- data.frame(CLIM.alps.df[,c("x","y")], predict(gfVars,CLIM.alps.df[,vars]))

# The multi-dimensional biological space can most effectively be represented by taking the principle components of the transformed grid and presenting the first two dimensions in a biplot.
# To remove neutral structure (here represented in proxy via longitude and latitude) from the visualisation, we calculate and plot the PCs excluding contributions from latitude and longitude.
vars.woLatLon <- vars [! vars %in% c("LATITUDE", "LONGITUDE")]
vars.woLatLon.gf <- vars.woLatLon
# Consider if you want to scale the PCA, or if you'd like to retain/preserve differences in the magnitude of genetic importance among the environmental variables (Fitzpatrick & Keller 2015, pg.6). Let's choose the latter (without scaling).
#PCs <- prcomp(CLIM.alps.df.GF[, vars.woLatLon], scale. = TRUE, center = TRUE)
#PCs <- prcomp(CLIM.alps.df.GF[, vars.woLatLon], scale. = FALSE, center = TRUE)
PCs <- prcomp(CLIM.alps.df.GF[, vars.woLatLon], scale. = FALSE, center = TRUE, rank. = 8) # To keep consistent with GDM analysis, and allow straightforward Procrustes superpositioning and generation of ensemble model, we only calculate the first 8 PCs here (since 3 variables non-important and removed in GDM). This is perfectly fine because we anyway only visualise the first 3 PCs.
PCs.full.gf <- PCs # Output for ensemble model generation
PCs.gf <- PCs$x # Output for ensemble model generation

# Let's convert this to a raster
PCs.df <- as.data.frame(PCs$x)
PCs.df <- cbind(CLIM.alps.df.GF[,c("x","y")], PCs.df)
PCs.rast <- rasterFromXYZ(PCs.df)
crs(PCs.rast) <- crs(CLIM.alps)

## Hindcast model to LGM
# Project model in LGM space
CLIM.alps.LGM.df <- as.data.frame(CLIM.alps.LGM, xy=TRUE, na.rm=TRUE)
CLIM.alps.LGM.df.GF <- data.frame(CLIM.alps.LGM.df[,c("x","y")], predict(gfVars,CLIM.alps.LGM.df[,vars]))
# Then calculate the PCA of the LGM biological space. We do so in the PCA space of the current climate (i.e. that this LGM PCA will be transformed, scaled and centered based on the current climate PCA), for apples-to-apples comparison
CLIM.alps.LGM.df.GF.tranformed <- scale(CLIM.alps.LGM.df.GF[, rownames(PCs$rotation)], scale = FALSE, center = PCs$center)
# Consider if you want to scale the PCA, or if you'd like to retain/preserve differences in the magnitude of genetic importance among the environmental variables (Fitzpatrick & Keller 2015, pg.6). Let's choose the latter (without scaling).
#CLIM.alps.LGM.df.GF.tranformed <- scale(CLIM.alps.LGM.df.GF.tranformed, scale = PCs$scale, center = FALSE)
PCs.LGM <-  as.matrix(CLIM.alps.LGM.df.GF.tranformed) %*% as.matrix(PCs$rotation)
PCs.LGM.gf <- PCs.LGM # Output for ensemble model generation

# Let's convert this to a raster
PCs.LGM.df <- as.data.frame(PCs.LGM)
PCs.LGM.df <- cbind(CLIM.alps.LGM.df.GF[,c("x","y")], PCs.LGM.df)
PCs.LGM.rast <- rasterFromXYZ(PCs.LGM.df)
crs(PCs.LGM.rast) <- crs(CLIM.alps.LGM)

# Define colour palette
if (col_palette == "alternative") {
  # For current
  a1 <- PCs.rast$PC1
  a2 <- PCs.rast$PC2
  a3 <- PCs.rast$PC3
  r <- a1 + a2
  g <- -a2
  b <- a3 + a2 - a1
  # For LGM
  a1.lgm <- PCs.LGM.rast$PC1
  a2.lgm <- PCs.LGM.rast$PC2
  a3.lgm <- PCs.LGM.rast$PC3
  r.lgm <- a1.lgm + a2.lgm
  g.lgm <- -a2.lgm
  b.lgm <- a3.lgm + a2.lgm - a1.lgm
  # Scale rasters by common scale
  # Current
  r.scaled <- (r-min(r@data@min, r.lgm@data@min)) / (max(r@data@max, r.lgm@data@max)-min(r@data@min, r.lgm@data@min))*255
  g.scaled <- (g-min(g@data@min, g.lgm@data@min)) / (max(g@data@max, g.lgm@data@max)-min(g@data@min, g.lgm@data@min))*255
  b.scaled <- (b-min(b@data@min, b.lgm@data@min)) / (max(b@data@max, b.lgm@data@max)-min(b@data@min, b.lgm@data@min))*255
  pcaRast.rgb <- stack(r.scaled, g.scaled, b.scaled)
  # LGM
  r.lgm.scaled <- (r.lgm-min(r@data@min, r.lgm@data@min)) / (max(r@data@max, r.lgm@data@max)-min(r@data@min, r.lgm@data@min))*255
  g.lgm.scaled <- (g.lgm-min(g@data@min, g.lgm@data@min)) / (max(g@data@max, g.lgm@data@max)-min(g@data@min, g.lgm@data@min))*255
  b.lgm.scaled <- (b.lgm-min(b@data@min, b.lgm@data@min)) / (max(b@data@max, b.lgm@data@max)-min(b@data@min, b.lgm@data@min))*255
  pcaLGMRast.rgb <- stack(r.lgm.scaled, g.lgm.scaled, b.lgm.scaled)
} else if (col_palette == "default") {
  # Scale rasters by common scale
  # Current
  r.scaled <- (PCs.rast[[1]]-min(PCs.rast[[1]]@data@min, PCs.LGM.rast[[1]]@data@min)) / (max(PCs.rast[[1]]@data@max, PCs.LGM.rast[[1]]@data@max)-min(PCs.rast[[1]]@data@min, PCs.LGM.rast[[1]]@data@min))*255
  g.scaled <- (PCs.rast[[2]]-min(PCs.rast[[2]]@data@min, PCs.LGM.rast[[2]]@data@min)) / (max(PCs.rast[[2]]@data@max, PCs.LGM.rast[[2]]@data@max)-min(PCs.rast[[2]]@data@min, PCs.LGM.rast[[2]]@data@min))*255
  b.scaled <- (PCs.rast[[3]]-min(PCs.rast[[3]]@data@min, PCs.LGM.rast[[3]]@data@min)) / (max(PCs.rast[[3]]@data@max, PCs.LGM.rast[[3]]@data@max)-min(PCs.rast[[3]]@data@min, PCs.LGM.rast[[3]]@data@min))*255
  pcaRast.rgb <- stack(r.scaled, g.scaled, b.scaled)
  # LGM
  r.lgm.scaled <- (PCs.LGM.rast[[1]]-min(PCs.rast[[1]]@data@min, PCs.LGM.rast[[1]]@data@min)) / (max(PCs.rast[[1]]@data@max, PCs.LGM.rast[[1]]@data@max)-min(PCs.rast[[1]]@data@min, PCs.LGM.rast[[1]]@data@min))*255
  g.lgm.scaled <- (PCs.LGM.rast[[2]]-min(PCs.rast[[2]]@data@min, PCs.LGM.rast[[2]]@data@min)) / (max(PCs.rast[[2]]@data@max, PCs.LGM.rast[[2]]@data@max)-min(PCs.rast[[2]]@data@min, PCs.LGM.rast[[2]]@data@min))*255
  b.lgm.scaled <- (PCs.LGM.rast[[3]]-min(PCs.rast[[3]]@data@min, PCs.LGM.rast[[3]]@data@min)) / (max(PCs.rast[[3]]@data@max, PCs.LGM.rast[[3]]@data@max)-min(PCs.rast[[3]]@data@min, PCs.LGM.rast[[3]]@data@min))*255
  pcaLGMRast.rgb <- stack(r.lgm.scaled, g.lgm.scaled, b.lgm.scaled)
}

## Prepare basemap
ELEV.alps[ELEV.alps < 1] <- NA
ELEV.LGM.alps[ELEV.LGM.alps < 1] <- NA
color.palette=colorRampPalette(c('gray92','gray55'))
#color.palette2=colorRampPalette(c('gray98','gray15'))

## Map in geographic space
crs(pcaRast.rgb) <- crs(CLIM.alps)
crs(pcaLGMRast.rgb) <- crs(CLIM.alps.LGM)
# Crop
pcaRast.rgb <-  crop(pcaRast.rgb, e_alps)
pcaLGMRast.rgb <-  crop(pcaLGMRast.rgb, e_alps)
# The following map plots predicted PC scores in geographic coordinates and represents continuous changes in inferred compositional patterns associated with the predictors.
if (plot_figures == "ALL" || plot_figures == "PARTIAL" ) {
  #par(mar = c(2, 2, 2, 2), mfrow = c(1, 2))
  par(mfrow=c(1,1), mar = c(4, 4, 3, 3))
  #par(mar = c(2, 2, 2, 5), mfrow = c(1, 1))
  plotRGB(pcaRast.rgb, r=1, g=2, b=3)
  #points(pop_geodata[,c(2,3)], pch=19, cex=2)
  plotRGB(pcaLGMRast.rgb, r=1, g=2, b=3)
}

## Should you wish to constrain the prediction to SDM predictions, we can apply SDM masks:
# For current time
current.mask <- crop(current.mask, extent(pcaRast.rgb))
pcaRast.masked <-mask(pcaRast.rgb, current.mask, maskvalue=0, updatevalue = NA)
pcaRast.masked.cropped <- crop(pcaRast.masked, e_alps)
#plotRGB(pcaRast.masked, r=1, g=2, b=3)
par(mfrow=c(1,1), mar = c(4, 4, 3, 3))
plot(ELEV.alps, col = color.palette(100), legend = FALSE)
plotRGB(pcaRast.masked.cropped, r=1, g=2, b=3, bgalpha = 0, add=TRUE)
#plot(current.mask, legend=FALSE, add=TRUE, col = "grey90")
#points(pop_geodata[,c(2,3)], pch=19, cex=2)
# For LGM
LGM.mask.plotSDM <- crop(LGM.mask.plot, extent(pcaLGMRast.rgb))
pcaLGMRast.masked <-mask(pcaLGMRast.rgb, LGM.mask.plotSDM, maskvalue=0, updatevalue = NA)
pcaLGMRast.masked.cropped <- crop(pcaLGMRast.masked, e_alps)
plot(ELEV.LGM.alps, col = color.palette(100), legend = FALSE)
#plotRGB(pcaLGMRast.masked, r=1, g=2, b=3)
plotRGB(pcaLGMRast.masked.cropped, r=1, g=2, b=3, bgalpha = 0, add=TRUE)
#plot(LGM.mask, legend=FALSE, add=TRUE, col = "grey90")

## We want to build a hull around our predicted presence points, to include the environment in areas around our predicted presences (and to not be constrained/overfitted to our SDM prediction)
# Acquire coordinates of presence points
# For current
current.mask.values <- getValues(current.mask.original)
current.mask.presencePoints.cellNumbers <- which(current.mask.values == 1)
current.mask.presencePoints.coords.temp <- xyFromCell(current.mask.original, current.mask.presencePoints.cellNumbers)
# For LGM
LGM.mask.values <- getValues(LGM.mask.original.plot)
LGM.mask.presencePoints.cellNumbers <- which(LGM.mask.values == 1)
LGM.mask.presencePoints.coords.temp <- xyFromCell(LGM.mask.original.plot, LGM.mask.presencePoints.cellNumbers)
# For LGM - updated. Here, we retain only eastern and denser/contiguous refugial clusters, given results from psi analysis and to remove sparse, disconnected presence points. We do this via dbscan.
# We've created and saved this prior. To recreate this clustered mask raster, see below:
# if(dbscan_LGM_mask == TRUE) {
#   # First convert the raster to a dataframe. Then select only presence points
#   pred.bin.LGM.mask.df <- as.data.frame(rasterToPoints(LGM.mask.original))
#   pred.bin.LGM.mask.clusters.df <- pred.bin.LGM.mask.df[pred.bin.LGM.mask.df$layer == 1, ]
#   pred.bin.LGM.mask.clusters.df <- pred.bin.LGM.mask.clusters.df[,c(1,2)]
#   # Perform density-based spatial clustering (DBSCAN). eps defines the size of the epsilon neighborhood and minPts defines number of minimum points in the eps region (for core points).
#   #DBSCAN.lgm.mask <- dbscan(pred.bin.LGM.mask.clusters.df, eps = 0.25, minPts = nrow(pred.bin.LGM.mask.clusters.df)/40)
#   DBSCAN.lgm.mask <- dbscan(pred.bin.LGM.mask.clusters.df, eps = 0.15, minPts = nrow(pred.bin.LGM.mask.clusters.df)/40)
#   # Format dataframe
#   pred.bin.LGM.mask.clusters.df["cluster"] <- DBSCAN.lgm.mask$cluster
#   # Filter according to identified clusters (here, we retain only the eastern ones, i.e. >0.8 psi TDOA probability)
#   pred.bin.LGM.mask.clusters <- pred.bin.LGM.mask.clusters.df[(pred.bin.LGM.mask.clusters.df$cluster == 4) | (pred.bin.LGM.mask.clusters.df$cluster == 6) | (pred.bin.LGM.mask.clusters.df$cluster == 7) | (pred.bin.LGM.mask.clusters.df$cluster == 8), ]
#   pred.bin.LGM.mask.clusters <- pred.bin.LGM.mask.clusters.df
#   # Convert to SpatialPointsDataframe
#   coordinates(pred.bin.LGM.mask.clusters) <- ~x+y
#   proj4string(pred.bin.LGM.mask.clusters) = CRS(prj_longlat)
#   # Specify as gridded
#   gridded(pred.bin.LGM.mask.clusters) <- TRUE
#   # Convert to raster
#   pred.bin.LGM.mask.clusters.raster <- raster(pred.bin.LGM.mask.clusters)
#   plot(pred.bin.LGM.mask.clusters.raster, col = brewer.pal(length(table(DBSCAN.lgm.mask$cluster)), "Paired"))
#   # Convert raster presence values to 1
#   pred.bin.LGM.mask.clusters.raster[(pred.bin.LGM.mask.clusters.raster == 4) | (pred.bin.LGM.mask.clusters.raster == 6) | (pred.bin.LGM.mask.clusters.raster == 7) | (pred.bin.LGM.mask.clusters.raster == 8)] <- 1
#   #plot(pred.bin.LGM.mask.clusters.raster)
# }

if (hull_method == "convex") {
  # Convex hull
  current.mask.convexHull <- chull(current.mask.presencePoints.coords.temp)
  LGM.mask.convexHull <- chull(LGM.mask.presencePoints.coords.temp)
  current.mask.presencePoints.coords <- current.mask.presencePoints.coords.temp[c(current.mask.convexHull, current.mask.convexHull[1]), ]  # closed polygon
  LGM.mask.presencePoints.coords <- LGM.mask.presencePoints.coords.temp[c(LGM.mask.convexHull, LGM.mask.convexHull[1]), ]  # closed polygon
} else if (hull_method == "concave") {
  # Concave hull
  current.mask.concaveHull <- concaveman(current.mask.presencePoints.coords.temp, 1)
  LGM.mask.concaveHull <- concaveman(LGM.mask.presencePoints.coords.temp, 2)
  current.mask.presencePoints.coords <- current.mask.concaveHull
  LGM.mask.presencePoints.coords <- LGM.mask.concaveHull
}

# Plot to check hull
if (plot_figures == "ALL") {
  par(mar = c(2, 2, 2, 2), mfrow = c(1, 2))
  plot(current.mask.presencePoints.coords.temp, pch=19)
  lines(current.mask.presencePoints.coords, col="red")
  plot(LGM.mask.presencePoints.coords.temp, pch=19)
  lines(LGM.mask.presencePoints.coords, col="red")
}

# Convert to SpatialPolygon
current.mask.polygon <- SpatialPolygons(list(Polygons(list(Polygon(current.mask.presencePoints.coords)), ID=1)), proj4string = crs(CLIM.alps))
LGM.mask.polygon <- SpatialPolygons(list(Polygons(list(Polygon(LGM.mask.presencePoints.coords)), ID=1)), proj4string = crs(CLIM.alps))

# Add buffer, to include surrounding areas (here width is in units of degrees)
#current.mask.polygon.wBuffer <- current.mask.polygon
#LGM.mask.polygon.wBuffer <- LGM.mask.polygon
current.mask.polygon.wBuffer <- buffer(current.mask.polygon, width = 0.025)
LGM.mask.polygon.wBuffer <- buffer(LGM.mask.polygon, width = 0.025)

# Plot maps with SDM hulls
par(mfrow=c(1,1), mar = c(4, 4, 3, 3))
plotRGB(pcaRast.rgb, r=1, g=2, b=3, maxpixels=1e7)
plot(current.mask.polygon.wBuffer, add = TRUE, lwd=10)
plotRGB(pcaLGMRast.rgb, r=1, g=2, b=3, maxpixels=1e7)
plot(LGM.mask.polygon.wBuffer, add=TRUE, lwd=10)

# Alternatively, we can mask the maps by the SDM hulls
# Plot (current)
pcaRast.polygonMasked <- mask(pcaRast.rgb, current.mask.polygon.wBuffer, updatevalue = NA)
plot(ELEV.alps, col = color.palette(100), legend = FALSE)
plotRGB(pcaRast.polygonMasked, r=1, g=2, b=3, bgalpha = 0, add=TRUE)
#points(pop_geodata[,c(2,3)], pch=19, cex=2)
# LGM
pcaLGMRast.polygonMasked <- mask(pcaLGMRast.rgb, LGM.mask.polygon.wBuffer, updatevalue = NA)
plot(ELEV.LGM.alps, col = color.palette(100), legend = FALSE)
plotRGB(pcaLGMRast.polygonMasked, r=1, g=2, b=3, bgalpha = 0, add=TRUE)

## Predicting biological change (difference in allele frequencies) through time - naive strategy (aka agnostic to LGM refugia)
# Calculate difference in allele frequencies (multi-dimensional Euclidean distance of transformed environmental predictors) between present and LGM at xy coordinates
# First, check that coordinates are ordered the same between the dataframes
all((CLIM.alps.df.GF[,c(1,2)] == CLIM.alps.LGM.df.GF[,c(1,2)]))
# We calculate difference in allele frequencie for both the full model and excluding latitude and longitude, the former to represent the expected AF shift from all contributing factors and the latter to represent the expected AF shift from only environmental variables.
CLIM.alps.GF.dist.naive <- as.data.frame((rowSums((CLIM.alps.df.GF[,-c(1,2)] - CLIM.alps.LGM.df.GF[,-c(1,2)])^2))^0.5)
# However note that in the naive strategy, w/ and wo/ latlon are identical because the lat and lon are the same in present and past times. Thus we can skip applying it here and this only needs to be applied in the refugia-explicit strategy.
#CLIM.alps.GF.dist.naive.woLatlon <- as.data.frame((rowSums((CLIM.alps.df.GF[,-c(1,2,3,4)] - CLIM.alps.LGM.df.GF[,-c(1,2,3,4)])^2))^0.5)
# Make dataframe of x,y,dist
CLIM.alps.GF.dist.coords.naive <- data.frame(CLIM.alps.df.GF[,c(1,2)], CLIM.alps.GF.dist.naive)
#CLIM.alps.GF.dist.coords.naive.woLatlon <- data.frame(CLIM.alps.df.GF[,c(1,2)], CLIM.alps.GF.dist.naive.woLatlon)
# Convert to raster
CLIM.alps.GF.dist.naive.rast <- rasterFromXYZ(CLIM.alps.GF.dist.coords.naive)
#CLIM.alps.GF.dist.naive.rast.woLatlon <- rasterFromXYZ(CLIM.alps.GF.dist.coords.naive.woLatlon)
crs(CLIM.alps.GF.dist.naive.rast) <- crs(CLIM.alps)
#crs(CLIM.alps.GF.dist.naive.rast.woLatlon) <- crs(CLIM.alps)

# Define plotting parameters
col_pal <- colorRampPalette(c('dodgerblue', 'red', 'gray10'))
color_levels=40
col_breaks <- seq(0,maxValue(CLIM.alps.GF.dist.naive.rast),length.out=color_levels+1)

# And plot
if (plot_figures == "PARTIAL" || plot_figures == "ALL") {
  par(mar = c(2, 2, 2, 5), mfrow = c(1, 1))
  plot(CLIM.alps.GF.dist.naive.rast, col=col_pal(n=color_levels), breaks=col_breaks, axes=FALSE, legend=F)
  plot(CLIM.alps.GF.dist.naive.rast, zlim=c(round(minValue(CLIM.alps.GF.dist.naive.rast),3), round(maxValue(CLIM.alps.GF.dist.naive.rast),3)), col=col_pal(n=color_levels), legend.only=T)
  #plot(CLIM.alps.GF.dist.naive.rast.woLatlon, col=col_pal(n=color_levels), breaks=col_breaks, axes=FALSE, legend=F)
  #plot(CLIM.alps.GF.dist.naive.rast.woLatlon, zlim=c(round(minValue(CLIM.alps.GF.dist.naive.rast.woLatlon),3), round(maxValue(CLIM.alps.GF.dist.naive.rast.woLatlon),3)), col=col_pal(n=color_levels), legend.only=T)
}

# Plot difference. Recall that in this naive strategy, w/ and wo/ latlon are identical because the lat and lon are the same in present and past times.
#plot(CLIM.alps.GF.dist.naive.rast - CLIM.alps.GF.dist.naive.rast.woLatlon, col=col_pal(n=color_levels), breaks=col_breaks, axes=FALSE, legend=F)

## Predicting biological change (difference in allele frequencies) through time - considering inferred LGM refugia
# Here, we calculate the difference in allele frequencies (multi-dimensional Euclidean distance of transformed environmental predictors) between all cells (CLIM.alps) and the closest cell within the SDM inferred LGM extent (pcaLGMRast.masked.offset).
# Note 1: To allow for a fast KNN strategy, we make the approximation of Euclidean space, which is fair for small geographic scales where the curvature is minimal.
# Note 2: however, we correct for the fact that at 44 degrees latitude and 10 degrees longitude, 1 degree latitude is approx 11/8 times longer than 1 degree longitude (111km vs 80km; see https://www.nhc.noaa.gov/gccalc.shtml) 
LGM.mask <- crop(LGM.mask, extent(pcaLGMRast.rgb))
pcaLGMRast.masked.offset <-mask(pcaLGMRast.rgb, LGM.mask, maskvalue=0, updatevalue = NA)
dist_correction_factor <- 111/78
current_xy <- na.omit(CLIM.alps.df.GF[,1:2])
LGM_SDM_xy <- na.omit(as.data.frame(rasterToPoints(pcaLGMRast.masked.offset[[1]])[,1:2]))
current_xy.adj <- current_xy
current_xy.adj$x <- current_xy.adj$x * dist_correction_factor
LGM_SDM_xy.adj <- LGM_SDM_xy
LGM_SDM_xy.adj$x <- LGM_SDM_xy.adj$x * dist_correction_factor
# We find the KNN nearest-neighbour index
knn_index <- knnx.index(data = LGM_SDM_xy.adj, query = current_xy.adj, k=1)
# And produce a dataframe which holds the coordinates and environmental variables for these closest points in LGM refugia (the row order of which corresponds to the present-time dataframe; i.e. matching rows in CLIM.alps.LGM.df.GF.closestPoints and CLIM.alps.df.GF hold closest pairs)
CLIM.alps.LGM.df.GF.rast <- rasterFromXYZ(CLIM.alps.LGM.df.GF)
closestPoints_inLGMrefugia <- LGM_SDM_xy[knn_index,]
CLIM.alps.LGM.df.GF.closestPoints <- data.frame(closestPoints_inLGMrefugia, extract(CLIM.alps.LGM.df.GF.rast, y=closestPoints_inLGMrefugia[,c("x","y")]))
# Calculate multi-dimensional Euclidean distance of transformed environmental predictors between all cells of CLIM.alps.df.GF and their respective closest cell within the SDM inferred LGM extent (CLIM.alps.LGM.df.GF.closestPoints))
# We do so for both the full model and excluding latitude and longitude, the former to represent the expected AF shift from all contributing factors and the latter to represent the expected AF shift from only environmental variables.
CLIM.alps.GF.dist <- as.data.frame((rowSums((CLIM.alps.df.GF[,-c(1,2)] - CLIM.alps.LGM.df.GF.closestPoints[,-c(1,2)])^2))^0.5)
CLIM.alps.GF.dist.woLatlon <- as.data.frame((rowSums((CLIM.alps.df.GF[,-c(1,2,3,4)] - CLIM.alps.LGM.df.GF.closestPoints[,-c(1,2,3,4)])^2))^0.5)
CLIM.alps.GF.dist.coords <- data.frame(CLIM.alps.df.GF[,c(1,2)], CLIM.alps.GF.dist)
CLIM.alps.GF.dist.coords.woLatlon <- data.frame(CLIM.alps.df.GF[,c(1,2)], CLIM.alps.GF.dist.woLatlon)
CLIM.alps.GF.dist.rast <- rasterFromXYZ(CLIM.alps.GF.dist.coords)
CLIM.alps.GF.dist.rast.woLatlon <- rasterFromXYZ(CLIM.alps.GF.dist.coords.woLatlon)
crs(CLIM.alps.GF.dist.rast) <- crs(CLIM.alps)
crs(CLIM.alps.GF.dist.rast.woLatlon) <- crs(CLIM.alps)
# Crop
CLIM.alps.GF.dist.rast <-  crop(CLIM.alps.GF.dist.rast, e_alps)
# Define plotting parameters
#col_breaks_1 <- seq(0,maxValue(CLIM.alps.GF.dist.rast),length.out=color_levels+1)
#col_breaks_1 <- seq(0,0.04358562,length.out=color_levels+1)
col_breaks_1 <- seq(0,0.04,length.out=color_levels+1)

# And plot
if (plot_figures == "PARTIAL" || plot_figures == "ALL") {
  par(mar = c(2, 2, 2, 5), mfrow = c(1, 1))
  plot(CLIM.alps.GF.dist.rast, col=col_pal(n=color_levels), breaks=col_breaks_1, axes=TRUE, legend=F, maxpixels=1e7)
  #plot(CLIM.alps.GF.dist.rast, zlim=c(round(minValue(CLIM.alps.GF.dist.rast),3), round(maxValue(CLIM.alps.GF.dist.rast),3)), col=col_pal(n=color_levels), legend.only=T)
  plot(CLIM.alps.GF.dist.rast, zlim=c(0,0.04), col=col_pal(n=color_levels), legend.only=T)
  #plot(CLIM.alps.GF.dist.rast, zlim=c(0,0.04358562), col=col_pal(n=color_levels), legend.only=T)
  plot(current.mask.polygon.wBuffer, add = TRUE, lwd=10)
  plot(CLIM.alps.GF.dist.rast.woLatlon, col=col_pal(n=color_levels), breaks=col_breaks_1, axes=FALSE, legend=F)
  plot(CLIM.alps.GF.dist.rast.woLatlon, zlim=c(round(minValue(CLIM.alps.GF.dist.rast.woLatlon),3), round(maxValue(CLIM.alps.GF.dist.rast.woLatlon),3)), col=col_pal(n=color_levels), legend.only=T)
}
# Let's mask by SDM and plot
CLIM.alps.GF.dist.rast.masked <- mask(CLIM.alps.GF.dist.rast, current.mask, maskvalue=0, updatevalue = NA)
CLIM.alps.GF.dist.rast.masked.cropped <-  crop(CLIM.alps.GF.dist.rast.masked, e_alps)
par(mar = c(4, 4, 3, 5), mfrow = c(1, 1))
plot(ELEV.alps, col = color.palette(20), legend = FALSE)
plot(CLIM.alps.GF.dist.rast.masked.cropped, col=col_pal(n=color_levels), breaks=col_breaks_1, axes=FALSE, legend=F, add=TRUE)
plot(CLIM.alps.GF.dist.rast, zlim=c(round(minValue(CLIM.alps.GF.dist.rast),3), round(maxValue(CLIM.alps.GF.dist.rast),3)), col=col_pal(n=color_levels), legend.only=T)

# Plot difference
#plot(CLIM.alps.GF.dist.rast - CLIM.alps.GF.dist.rast.woLatlon, col=col_pal(n=color_levels), breaks=col_breaks, axes=FALSE, legend=F)

# ## Biplot of biological space
# # Define colour palette
# if (col_palette == "alternative") {
#   a1.biplot <- PCs$x[, 1]
#   a2.biplot <- PCs$x[, 2]
#   a3.biplot <- PCs$x[, 3]
#   r.biplot <- a1.biplot + a2.biplot
#   g.biplot <- -a2.biplot
#   b.biplot <- a3.biplot + a2.biplot - a1.biplot
#   # Scale between 0-255. Use same standardisation scale as above. Recall that PCs.rast[[1]]@data@min, for all i.
#   r.biplot <- (r.biplot - min(r@data@min, r.lgm@data@min)) / (max(r@data@max, r.lgm@data@max)-min(r@data@min, r.lgm@data@min))*255
#   g.biplot <- (g.biplot - min(g@data@min, g.lgm@data@min)) / (max(g@data@max, g.lgm@data@max)-min(g@data@min, g.lgm@data@min))*255
#   b.biplot <- (b.biplot - min(b@data@min, b.lgm@data@min)) / (max(b@data@max, b.lgm@data@max)-min(b@data@min, b.lgm@data@min))*255
# } else if (col_palette == "default") {
#   r.biplot <- PCs$x[, 1]
#   g.biplot <- PCs$x[, 2]
#   b.biplot <- PCs$x[, 3]
#   # Scale between 0-255. Use same standardisation scale as above. Recall that PCs.rast[[1]]@data@min == min(PCs$x[, i]), for all i.
#   r.biplot <- (r.biplot - min(PCs.rast[[1]]@data@min, PCs.LGM.rast[[1]]@data@min)) / (max(PCs.rast[[1]]@data@max, PCs.LGM.rast[[1]]@data@max)-min(PCs.rast[[1]]@data@min, PCs.LGM.rast[[1]]@data@min))*255
#   g.biplot <- (g.biplot - min(PCs.rast[[2]]@data@min, PCs.LGM.rast[[2]]@data@min)) / (max(PCs.rast[[2]]@data@max, PCs.LGM.rast[[2]]@data@max)-min(PCs.rast[[2]]@data@min, PCs.LGM.rast[[2]]@data@min))*255
#   b.biplot <- (b.biplot - min(PCs.rast[[3]]@data@min, PCs.LGM.rast[[3]]@data@min)) / (max(PCs.rast[[3]]@data@max, PCs.LGM.rast[[3]]@data@max)-min(PCs.rast[[3]]@data@min, PCs.LGM.rast[[3]]@data@min))*255
# }
# 
# # The environmental variables may be shown as vectors, e.g. here limited to the most important predictors — in this example, variables to show as vectors are selected.
# nvs <- dim(PCs$rotation)[1]
# lv <- length(vars.woLatLon)
# vind <- rownames(PCs$rotation) %in% vars.woLatLon
# scal <- 40
# xrng <- range(PCs$x[, 1], PCs$rotation[, 1]/scal) * 1.1
# yrng <- range(PCs$x[, 2], PCs$rotation[, 2]/scal) * 1.1
# # And plot. Different coordinate positions in the biplot represent differing compositions, as associated with the predictors.
# par(mar = c(2, 2, 2, 2), mfrow = c(1, 1))
# plot((PCs$x[, 1:2]), xlim = xrng, ylim = yrng, pch = ".", cex = 4, col = rgb(r.biplot, g.biplot, b.biplot, max = 255), asp = 1, xaxt="n", yaxt="n")
# points(PCs$rotation[!vind, 1:2]/scal, pch = "+")
# arrows(rep(0, lv), rep(0, lv), PCs$rotation[vars.woLatLon, 1]/scal, PCs$rotation[vars.woLatLon, 2]/scal, length = 0.1, lwd = 2.5)
# jit <- 0.0015
# text(PCs$rotation[vars.woLatLon, 1]/scal + jit * sign(PCs$rotation[vars.woLatLon, 1]), PCs$rotation[vars.woLatLon, 2]/scal + jit * sign(PCs$rotation[vars.woLatLon, 2]), labels = vars.woLatLon)
# # Draw axes through origin
# axis(1, pos=0)
# axis(2, pos=0)
# # Add the location of populations in biological space (circles)
# # First transform the population environmental predictors, which are available from gf$X.
# Trns_site <- predict(gfVars)
# PCsites <- predict(PCs, Trns_site[, vars.woLatLon])
# PCsites <- unique(PCsites) # since we combined n gradientForest objects via combinedGradientForest, we have n*pops; i.e. duplicate rows.
# points(PCsites[, 1:2])
# # Add the weight mean location of species (crosses).
# if (combined_gradientForest == FALSE) {
#   SpsWtd <- sweep(gfVars$Y, 2, apply(gfVars$Y, 2, min), "-")
#   SpsWtdPCs <- (t(SpsWtd) %*% (PCsites[, 1:2]))/colSums(SpsWtd)
#   points(SpsWtdPCs, col = "red", pch = "+")
#   #If required, the frequency of any given SNP may be plotted on the biplot. For example the first SNP from gf$Y:
#   #sp <- colnames(SpsWtd)[1]
#   #points(PCsites[, 1:2], col = "blue", cex = SpsWtd[, sp]/2)
# }
