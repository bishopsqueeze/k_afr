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
##  Reli: topographic Relief calculated from Shuttle Radar Topography mission elevation data
##  TMAP & TMFI: average long-term Tropical Rainfall Monitoring Mission data
##          (TMAP = mean annual precipitation, TMFI = modified Fournier index)
##------------------------------------------------------------------

##------------------------------------------------------------------
## Load Libraries
##------------------------------------------------------------------
library(forecast)

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
load("01_AfricaRawTest.Rdata")
load("01_AfricaRawTrain.Rdata")

##------------------------------------------------------------------
## Append NA(s) for the targets to the test data
##------------------------------------------------------------------
test.raw[, c("ca","p","ph","soc","sand")]   <- NA

##------------------------------------------------------------------
## Combine the two dataset
##------------------------------------------------------------------
comb.raw    <- rbind(train.raw, test.raw)


##******************************************************************
## Process the spectra
##******************************************************************

##------------------------------------------------------------------
## Step 1: Pre-process the data
##------------------------------------------------------------------

## extract raw spectra
spectra.raw <- as.matrix(comb.raw[ , c(2:3579)])

## compute a moving average
spectra.ma3 <- t(apply(spectra.raw, 1, ma, 3))


## identify known problem areas
spectra.co2_bands   <- c("m2379.76", "m2377.83", "m2375.9", "m2373.97", "m2372.04",
"m2370.11", "m2368.18", "m2366.26","m2364.33", "m2362.4",
"m2360.47", "m2358.54", "m2356.61", "m2354.68", "m2352.76")

## identify problem columns
spectra.co2_index   <- which(colnames(spectra.raw) %in% spectra.co2_bands)

## extract wavenumbers and nanometers
spectra.wn  <- as.numeric(gsub("m","",colnames(spectra.raw)))  ## spectral wave numbers
spectra.nm  <- 1 / (spectra.wn * (1.0e-07))                 ## in nano meters


##------------------------------------------------------------------
## Step 2:  Compute derivatives of the spectra
##------------------------------------------------------------------

## components used to compute the finite differences
fx              <- spectra.raw
fx_plus_h       <- cbind(NA, spectra.raw)[, -(dim(spectra.raw)[2]+1)]
fx_minus_h      <- cbind(spectra.raw[, -1], NA)

## compute the first/second derivatives
spectra.der.1   <- fx_plus_h - fx
spectra.der.2   <- fx_plus_h - 2*fx + fx_minus_h


##------------------------------------------------------------------
## Step 3:  Mean-subtract the data chunks
##------------------------------------------------------------------

## mean-subtract the spectra
#spectra.raw.mean    <- apply(spectra.ma3, 2, mean)
#spectra.raw.msub    <- as.matrix(spectra.ma3 - spectra.raw.mean)

## mean-subtract the first derivative
#spectra.der.1.mean    <- apply(spectra.der.1, 2, mean)
#spectra.der.1.msub    <- as.matrix(spectra.der.1 - spectra.der.1.mean)

## mean-subtract the second derivative
#spectra.der.2.mean    <- apply(spectra.der.2, 2, mean)
#spectra.der.2.msub    <- as.matrix(spectra.der.2 - spectra.der.2.mean)


##------------------------------------------------------------------
## Step 4:  Zero-out the co2 band & and edge effects
##------------------------------------------------------------------

## co2 bands
#spectra.raw.msub[, spectra.co2_index]       <- 0
#spectra.der.1.msub[, spectra.co2_index]     <- 0
#spectra.der.2.msub[, spectra.co2_index]     <- 0

## NAs
spectra.der.1.msub[is.na(spectra.der.1.msub)] <- 0
spectra.der.2.msub[is.na(spectra.der.2.msub)] <- 0


##------------------------------------------------------------------
## STEP 5: Perform the SVD
##------------------------------------------------------------------
##  - Column of U are eigenvectors of AAT
##  - Columns of V are eigenvectors of ATA
##  - S is a diagonal matrix containing the square rroots of eigenvalues from U or V
##------------------------------------------------------------------

## raw spectra
spectra.raw.svd         <- svd(t(spectra.raw))
spectra.pct_var         <- cumsum(spectra.raw.svd$d / sum(spectra.raw.svd$d))

## 1st derivatives
#spectra.der.1.svd       <- svd(t(spectra.der.1.msub))
#spectra.der.1.pct_var   <- cumsum(spectra.der.1.svd$d / sum(spectra.der.1.svd$d))

## 2nd derivatives
#spectra.der.2.svd       <- svd(t(spectra.der.2.msub))
#spectra.der.2.pct_var   <- cumsum(spectra.der.2.svd$d / sum(spectra.der.2.svd$d))








##------------------------------------------------------------------
## Save results
##------------------------------------------------------------------
##save.image("02_AfricaProcessedData.Rdata")