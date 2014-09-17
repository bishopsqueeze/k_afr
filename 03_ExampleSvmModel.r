##------------------------------------------------------------------
## Notes on the data:
##------------------------------------------------------------------
## Data fields
##
## SOC, pH, Ca, P, Sand are the five target variables for predictions.
##
## The data have been monotonously transformed from the original
## measurements and thus include negative values.
##
##  PIDN: unique soil sample identifier
##  SOC: Soil organic carbon
##  pH: pH values
##  Ca: Mehlich-3 extractable Calcium
##  P: Mehlich-3 extractable Phosphorus
##  Sand: Sand content
##
##  m7497.96 - m599.76: There are 3,578 mid-infrared absorbance measurements.
##      For example, the "m7497.96" column is the absorbance at wavenumber
##      7497.96 cm-1. We suggest you to remove spectra CO2 bands which are
##      in the region m2379.76 to m2352.76, but you do not have to.
##
##  Depth: Depth of the soil sample (2 categories: "Topsoil", "Subsoil")
##
##  We have also included some potential spatial predictors from remote
##  sensing data sources. Short variable descriptions are provided below
##  and additional descriptions can be found at AfSIS data. The data
##  have been mean centered and scaled.
##
##  BSA:    average long-term Black Sky Albedo measurements from MODIS
##          satellite images (BSAN = near-infrared, BSAS = shortwave, BSAV = visible)
##  CTI:    compound topographic index calculated from Shuttle Radar
##          Topography Mission elevation data
##  ELEV:   Shuttle Radar Topography Mission elevation data
##  EVI:    average long-term Enhanced Vegetation Index from MODIS satellite images.
##  LST:    average long-term Land Surface Temperatures from MODIS satellite images
##          (LSTD = day time temperature, LSTN = night time temperature)
##  Ref:    average long-term Reflectance measurements from MODIS satellite images
##          (Ref1 = blue, Ref2 = red, Ref3 = near-infrared, Ref7 = mid-infrared)
##  Reli:   topographic Relief calculated from Shuttle Radar Topography mission elevation data
##  TMAP & TMFI: average long-term Tropical Rainfall Monitoring Mission data
##          (TMAP = mean annual precipitation, TMFI = modified Fournier index)
##------------------------------------------------------------------

##------------------------------------------------------------------
## Load libraries
##------------------------------------------------------------------
library(data.table)
library(caret)
library(foreach)
library(doMC)

##------------------------------------------------------------------
## register cores
##------------------------------------------------------------------
registerDoMC(4)

##------------------------------------------------------------------
## Clear the workspace
##------------------------------------------------------------------
rm(list=ls())

##------------------------------------------------------------------
## Set the working directory
##------------------------------------------------------------------
setwd("/Users/alexstephens/Development/kaggle/africa/data/proc")

##------------------------------------------------------------------
## Load raw data
##------------------------------------------------------------------
#load("01_AfricaRawTest.Rdata")
load("01_AfricaRawTrain.Rdata")

##------------------------------------------------------------------
## Define target variables
##------------------------------------------------------------------
target.vars <- c("ca","p","ph","soc","sand")


##******************************************************************
## Simple beat-the-benchmark model
##******************************************************************

##------------------------------------------------------------------
## Step 1: Identify features in the data
##------------------------------------------------------------------

## identify known problem areas
spectra.co2_bands   <- c("m2379.76", "m2377.83", "m2375.9", "m2373.97", "m2372.04",
"m2370.11", "m2368.18", "m2366.26","m2364.33", "m2362.4",
"m2360.47", "m2358.54", "m2356.61", "m2354.68", "m2352.76")

## define indices
spec.idx            <- 2:3579
target.idx          <- which(names(train.raw) %in% target.vars)
nonspec.idx         <- setdiff(2:3600, union(spec.idx, target.idx))
spectra.co2_index   <- which(colnames(train.raw) %in% spectra.co2_bands)

## create a binary variables for depth
train.raw$depth     <- ifelse(train.raw$depth == "Topsoil", 1, 0)


##------------------------------------------------------------------
## Step 2:  Use caret to find the optimal tuning params for each fit
##------------------------------------------------------------------

ctrl        <- trainControl(method = "cv", savePred=TRUE)
svm.grid    <- expand.grid(C = 2^(seq(-2,10,1)))
tmp.res     <- list()

## loop over each target variable
for (i in 1:length(target.vars)) {
    
    ## load all raw variables
    tmp.data            <- train.raw[, c(target.idx[i], spec.idx, nonspec.idx)]
    
    ## use all variables in the file
    tmp.formula         <- as.formula(paste0(target.vars[i], " ~ ."))
    
    ## save the results
    tmp.res[[i]]        <- train(tmp.formula, data=tmp.data, method = "svmLinear", verbose=FALSE, trControl = ctrl, tuneGrid=svm.grid)
}


## minmum RMSE (all variables)
## [1] = 0.335 +- 0.113, C = 1
## [2] = 0.830 +- 0.355, C = 16
## [3] = 0.326 +- 0.058, C = 8
## [4] = 0.303 +- 0.072, C = 1
## [5] = 0.333 +- 0.049, C = 1

## use the propectr package to pre-process some of the spectra


##------------------------------------------------------------------
## Save results
##------------------------------------------------------------------
##save.image("02_AfricaProcessedData.Rdata")