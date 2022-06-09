##### This script calculates (gradientForest) R2-weighted neutrality statistics (e.g. Tajima's D, Fay & Wu's H, Zeng's E)

## Load libraries
library(gradientForest)

## Define functions
# Define function to extract scaffold-position names
extract_names <- function(x) {  
  strsplit(x, "\\.")[[1]][2]
}

# Define function to calculate the (n-1)th harmonic number
An <- function(n){
  sum_h <- 0
  for (i in seq(1,n-1)) {sum_h <- sum_h + 1/i}
  return(sum_h)
}

# Define function to calculate the (n-1)th 2nd-order harmonic number
An2 <- function(n){
  sum_h <- 0
  for (i in seq(1,n-1)) {sum_h <- sum_h + (1/i)**2}
  return(sum_h)
}

# Define function to calculate Tajima's D
# For Tajima's D (Source: Hamilton 2009. Population Genetics; https://en.wikipedia.org/wiki/Tajima%27s_D, https://ocw.mit.edu/courses/health-sciences-and-technology/hst-508-quantitative-genomics-fall-2005/study-materials/tajimad1.pdf):
calc_tajimaD <- function(n, thetaW, thetaPi) {
  # ANGSD per site thetas are given in natural log scale. So first, convert these thetas to original scale
  thetaW <- exp(thetaW)
  thetaPi <-  exp(thetaPi)
  # Then, let's define some identities to below to make things simpler 
  b1 = ((n+1) / (3*(n-1)))
  b2 = ((2*((n**2)+n+3)) / (9 * n * (n-1)))
  c1 = (b1 - (1/An(n)))
  c2 = (b2 - ((n+2) / (An(n)*n)) + (An2(n) / An(n)**2))
  e1 = (c1 / An(n))
  e2 = (c2 / ((An(n)**2) + An2(n)))
  # Then, calculate Tajima's D
  S = thetaW * An(n)                  # number of segregating sites. See: https://en.wikipedia.org/wiki/Watterson_estimator
  Var_d = ((e1*S) + (e2*S*(S-1)))     # this if the variance of Tajima's d
  SD_d = (Var_d**0.5)                 # this is the standard deviation of Tajima's d
  tajima_d = thetaPi - thetaW         # Tajima's d
  tajimaD = tajima_d / SD_d           # Tajima's D (normalised Tajima's d)
  # Return output
  return(tajimaD)
}

# Define function to calculate Fay & Wu's H.
# For Fay & Wu's H. (see Fay & Wu 2000, Zeng et al.2006, Ferretti et al. 2017 for formula)
# This seems to be off from the weighted scaffold mean by a factor of 2. Likely a missing coefficient somewhere. As a correction, we multiply our final estimate by 1/2 at the end.
calc_fayWuH <- function(n, thetaW, thetaPi, thetaL) {
  # ANGSD per site thetas are given in natural log scale. So first, convert these thetas to original scale
  thetaW <- exp(thetaW)
  thetaPi <-  exp(thetaPi)
  thetaL <-  exp(thetaL)
  r1 = (n-2)/(6*(n-1))
  r2 = ( ((18*(n^2))*((3*n)+2)*An2(n+1)) - ((88*(n^3)) + (9*(n^2)) - (13*n) + 6) ) / (9*n*((n-1)^2))
  S = thetaW * An(n)      # number of segregating sites. See: https://en.wikipedia.org/wiki/Watterson_estimator
  thet2 = (S*(S-1)) / (((An(n))^2) + An2(n))
  #thet2 = theta_W_gW^2
  Var_H = (r1*thetaW) + (r2*thet2)
  #Var_H = (r1*(thet2^0.5)) + (r2*thet2)
  SD_H = (Var_H**0.5)
  FayWu_H = (thetaPi - thetaL) / SD_H # Fay & Wu's normalised H
  FayWu_H = 0.5*FayWu_H # Correction (factor of 1/2)
  # Return output
  return(FayWu_H)
}

# Define function to calculate Zeng's E.
# For Zeng's E. (see Zeng et al.2006, Ferretti et al. 2017 for formula)
# This seems to be off from the weighted scaffold mean by a factor of 2. Likely a missing coefficient somewhere. As a correction, we multiply our final estimate by 1/2 at the end.
calc_zengE <- function(n, thetaW, thetaL) {
  # ANGSD per site thetas are given in natural log scale. So first, convert these thetas to original scale
  thetaW <- exp(thetaW)
  thetaL <-  exp(thetaL)
  r1 = (n/(2*(n-1))) - (1/An(n))
  r2 = (An2(n)/((An(n))^2)) + (2*((n/(n-1))^2)*An2(n)) - (((2*((n*An2(n)) - n + 1))/((n-1)*An(n))) - (((3*n)+1)/(n-1)))
  S = thetaW * An(n)      # number of segregating sites. See: https://en.wikipedia.org/wiki/Watterson_estimator & Zeng et al.2006
  thet2 = (S*(S-1)) / (((An(n))^2) + An2(n))
  #thet2 = theta_W_gW^2
  Var_E = (r1*thetaW) + (r2*thet2)
  #Var_E = (r1*(thet2^0.5)) + (r2*thet2)
  SD_E = (Var_E**0.5)
  Zeng_E = (thetaL - thetaW) / SD_E # Zeng's normalised E
  Zeng_E = 0.5*Zeng_E # Correction (factor of 1/2)
  # Return output
  return(Zeng_E)
}

## Define input variables

# Read arguments from command line
args <- commandArgs(trailingOnly = TRUE)
# Assign argument to population name
pop <- as.character(args[1])
#pop <- "Nevache"

# Define directories
#out_dir <- "/Users/luqman/Desktop/"
out_dir <- "/cluster/work/gdc/shared/p461/secondNGSdataset/Population_Genetics/GF_GDM/14indPerPop_unfolded_ANC_Lusitanus_Sylvestris_GL2_exonRegions_GF_R2weighted_thetas_woSolarRad_revised/"
#gf_dir <- "/Users/luqman/Documents/Hirzi/ETHPHD/PopulationGenetics_wholeGenome/Variant calling for Central Alps/"
gf_dir <- "/cluster/work/gdc/shared/p461/secondNGSdataset/Population_Genetics/GF_GDM/"
#thetas_dir <- "/Users/luqman/Desktop/"
thetas_dir <- "/cluster/work/gdc/shared/p461/secondNGSdataset/Population_Genetics/FST_GDM/14indPerPop_unfolded_ANC_Lusitanus_Sylvestris_SAF_files_ANGSD_GL2_exonRegions/"

# Define the number of samples
n <- 28

# Select gradientForest results
# Choose default options (method=2, standardize="before"). Choose between cor=0.5,0.7 - results appear very similar.
gfVars <- readRDS(paste0(gf_dir,"gfModel_cor0.5_selectedPredVars_woSolarRad_combined_method2standardizeBefore_allSplits.rds"))
#gfVars <- readRDS(paste0(gf_dir,"gfModel_cor0.5_selectedPredVars_combined_method2standardizeBefore_allSplits.rds"))
#gfVars <- readRDS(paste0(gf_dir,"gfModel_cor0.5_selectedPredVars_combined_method2standardizeAfter_allSplits.rds"))
#gfVars <- readRDS(paste0(gf_dir,"gfModel_cor0.7_selectedPredVars_combined_method2standardizeBefore_allSplits.rds"))
#gfVars <- readRDS(paste0(gf_dir,"gfModel_cor0.7_selectedPredVars_combined_method2standardizeAfter_allSplits.rds"))
#gfVars <- readRDS(paste0(gf_dir,"gfModel_cor0.7_allPredVars_combined_method2standardizeBefore_allSplits.rds"))

## Format gradientForest results and ANGSD per-site thetas dataframes
r2_allpreds_raw <- as.data.frame(t(gfVars$imp.rsq))[-c(1),]
r2_allpreds_raw[] <- lapply(r2_allpreds_raw, function(x) {
  if(is.factor(x)) as.numeric(as.character(x)) else x
})
#sapply(r2_allpreds_raw, class)
r2_allpreds_raw[,ncol(r2_allpreds_raw)+1] <- rowSums(r2_allpreds_raw[,1:10])
colnames(r2_allpreds_raw)[ncol(r2_allpreds_raw)] <- "r2_ENV"
r2_allpreds_raw[,ncol(r2_allpreds_raw)+1] <- rowSums(r2_allpreds_raw[,1:12])                
colnames(r2_allpreds_raw)[ncol(r2_allpreds_raw)] <- "r2_total"
r2_allpreds_raw[,ncol(r2_allpreds_raw)+1] <- as.data.frame(sapply(rownames(r2_allpreds_raw),extract_names, simplify = TRUE))
colnames(r2_allpreds_raw)[ncol(r2_allpreds_raw)] <- "position"
thetaPerSite_table <- read.table(paste0(thetas_dir, pop, ".thetas.table"), header = TRUE, comment.char="!")
#thetaPerSite_table <- read.table(paste0(thetas_dir, pop, ".thetas.table"), header = TRUE, comment.char="!")
thetaPerSite_table[,ncol(thetaPerSite_table)+1] <- paste0(thetaPerSite_table$X.Chromo,"_pos",thetaPerSite_table$Pos)
colnames(thetaPerSite_table)[ncol(thetaPerSite_table)] <- "position"

## Need to check that position system is the same (between gf and angsd). I.e. 0-based or 1-based.
# -> ANGSD is 1-indexed (http://www.popgen.dk/angsd/index.php/Sites). So is VCF (https://www.internationalgenome.org/wiki/Analysis/Variant%20Call%20Format/vcf-variant-call-format-version-41, http://samtools.github.io/hts-specs/VCFv4.1.pdf, http://samtools.github.io/hts-specs/VCFv4.2.pdf. My guess is that popStats then is also 1-indexed (as it's derived from vcf)
# So both are likely 1-indexed and can continue.

## Merge dataframes (inner join) by scaff-position
rsq_perSite <- r2_allpreds_raw[,c(13:15)]
thetaR2_PerSite_df <- merge(thetaPerSite_table, rsq_perSite, by = "position")
colnames(thetaR2_PerSite_df)[2] <- "scaffold"

## Calculate per site neutrality statitics
thetaR2_PerSite_df[,ncol(thetaR2_PerSite_df)+1] <- calc_tajimaD(n, thetaW=thetaR2_PerSite_df$Watterson, thetaPi=thetaR2_PerSite_df$Pairwise)
thetaR2_PerSite_df[,ncol(thetaR2_PerSite_df)+1] <- calc_fayWuH(n, thetaW=thetaR2_PerSite_df$Watterson, thetaPi=thetaR2_PerSite_df$Pairwise, thetaL=thetaR2_PerSite_df$thetaL)
thetaR2_PerSite_df[,ncol(thetaR2_PerSite_df)+1] <- calc_zengE(n, thetaW=thetaR2_PerSite_df$Watterson, thetaL=thetaR2_PerSite_df$thetaL)
colnames(thetaR2_PerSite_df)[ncol(thetaR2_PerSite_df)-2] <- "TajD"
colnames(thetaR2_PerSite_df)[ncol(thetaR2_PerSite_df)-1] <- "FayWuH"
colnames(thetaR2_PerSite_df)[ncol(thetaR2_PerSite_df)] <- "ZengE"

## Calculate R2-weighted neutrality statistics
thetaR2_PerSite_df[,ncol(thetaR2_PerSite_df)+1] <- thetaR2_PerSite_df$r2_total * thetaR2_PerSite_df$TajD
thetaR2_PerSite_df[,ncol(thetaR2_PerSite_df)+1] <- thetaR2_PerSite_df$r2_total * thetaR2_PerSite_df$FayWuH
thetaR2_PerSite_df[,ncol(thetaR2_PerSite_df)+1] <- thetaR2_PerSite_df$r2_total * thetaR2_PerSite_df$ZengE
colnames(thetaR2_PerSite_df)[ncol(thetaR2_PerSite_df)-2] <- "TajDxR2"
colnames(thetaR2_PerSite_df)[ncol(thetaR2_PerSite_df)-1] <- "FayWuHxR2"
colnames(thetaR2_PerSite_df)[ncol(thetaR2_PerSite_df)] <- "ZengExR2"
#TajD_weightedR2 <- sum(thetaR2_PerSite_df$TajDxR2)/sum(thetaR2_PerSite_df$r2)
#FayWuH_weightedR2 <- sum(thetaR2_PerSite_df$FayWuHxR2)/sum
thetaR2_PerSite_df[,ncol(thetaR2_PerSite_df)+1] <- thetaR2_PerSite_df$r2_ENV * thetaR2_PerSite_df$TajD
thetaR2_PerSite_df[,ncol(thetaR2_PerSite_df)+1] <- thetaR2_PerSite_df$r2_ENV * thetaR2_PerSite_df$FayWuH
thetaR2_PerSite_df[,ncol(thetaR2_PerSite_df)+1] <- thetaR2_PerSite_df$r2_ENV * thetaR2_PerSite_df$ZengE
colnames(thetaR2_PerSite_df)[ncol(thetaR2_PerSite_df)-2] <- "TajDxR2ENV"
colnames(thetaR2_PerSite_df)[ncol(thetaR2_PerSite_df)-1] <- "FayWuHxR2ENV"
colnames(thetaR2_PerSite_df)[ncol(thetaR2_PerSite_df)] <- "ZengExR2ENV"

## Write out table
write.table(thetaR2_PerSite_df, paste0(out_dir, pop,".thetas.R2weighted"), col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
