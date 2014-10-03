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
load("01_AfricaRawData.Rdata")

##------------------------------------------------------------------
## Target Vars
##------------------------------------------------------------------
target.vars     <- c("ca", "p", "ph", "soc", "sand")    ## raw
target.bc.vars  <- paste0(target.vars, ".bc")           ## box-cox x-formed


##******************************************************************
## Modify base features
##******************************************************************

##------------------------------------------------------------------
## Clone the raw data
##------------------------------------------------------------------
train.tmp  <- train.raw
test.tmp   <- test.raw

##------------------------------------------------------------------
## Compute BoxCox lambdas; add 2 to ensure positive values
##------------------------------------------------------------------
target.bc.lambda    <- vector(, length=length(target.vars))
for (i in 1:length(target.vars)) {
    tmp.var                                 <- target.vars[i]
    target.bc.lambda[i]                     <- BoxCox.lambda((train.tmp[, tmp.var]+2))
    train.tmp[, paste0(tmp.var, ".bc")]     <- BoxCox((train.tmp[, tmp.var]+2), target.bc.lambda[i])
}
names(target.bc.lambda) <- target.vars

##------------------------------------------------------------------
## Append NA(s) for the targets to the test data
##------------------------------------------------------------------
test.tmp[, target.vars]     <- NA
test.tmp[, target.bc.vars]  <- NA

##------------------------------------------------------------------
## Append source tags
##------------------------------------------------------------------
train.tmp$id    <- 1
test.tmp$id     <- 0

##------------------------------------------------------------------
## Create binary flags for the depth "factor"
##  Topsoil = 1
##  Subsoil = 0
##------------------------------------------------------------------
train.tmp$depth     <- ifelse(train.tmp$depth == "Topsoil", 1, 0)
test.tmp$depth      <- ifelse(test.tmp$depth == "Topsoil", 1, 0)

##------------------------------------------------------------------
## Combine the two datasets into a single file
##------------------------------------------------------------------
comb.tmp           <- rbind(train.tmp, test.tmp)


##******************************************************************
## Process the spectra
##******************************************************************

##------------------------------------------------------------------
## Step 1: Extract spectra & basic feautres
##------------------------------------------------------------------

## define spectra index
spectra.index   <- 2:3579

## extract raw spectra
spectra <- as.matrix(comb.tmp[ , spectra.index])

## identify known problem areas
spectra.co2_bands   <-
c("m2379.76", "m2377.83", "m2375.9",  "m2373.97", "m2372.04",
  "m2370.11", "m2368.18", "m2366.26", "m2364.33", "m2362.4",
  "m2360.47", "m2358.54", "m2356.61", "m2354.68", "m2352.76")

## identify problem columns
spectra.co2_index   <- which(colnames(spectra) %in% spectra.co2_bands)

## extract wavenumbers and wavelengths (in nanometers)
spectra.wn      <- as.numeric(gsub("m","",colnames(spectra)))   ## spectral wave numbers
spectra.nm      <- 1 / (spectra.wn * (1.0e-07))                 ## in nanometers


##------------------------------------------------------------------
## Step 2: Compute derivatives of the spectra
##------------------------------------------------------------------

## components used to compute the finite differences
fx              <- spectra
fx_plus_h       <- cbind(NA, spectra)[, -(ncol(spectra)+1)]
fx_minus_h      <- cbind(spectra[, -1], NA)

## compute the first/second derivatives
spectra.der.1   <- fx_plus_h - fx
spectra.der.2   <- fx_plus_h - 2*fx + fx_minus_h

## append column names
colnames(spectra.der.1) <- gsub("m", "d1_", colnames(spectra))
colnames(spectra.der.2) <- gsub("m", "d2_", colnames(spectra))


##------------------------------------------------------------------
## Step 3:  Zero-out the co2 band & and edge effects
##------------------------------------------------------------------

## NA(s)
spectra.der.1[is.na(spectra.der.1)] <- 0
spectra.der.2[is.na(spectra.der.2)] <- 0


##------------------------------------------------------------------
## Step 4: Append the derivatives to the original data
##------------------------------------------------------------------
all.proc    <- data.frame(
                    comb.tmp[, -spectra.index],
                    spectra[, -spectra.co2_index],
                    spectra.der.1[, -spectra.co2_index],
                    spectra.der.2[, -spectra.co2_index])

train.proc  <- all.proc[(all.proc$id == 1), ]
test.proc   <- all.proc[(all.proc$id == 0), ]



##******************************************************************
## Save results
##******************************************************************
save(
    all.proc,
    train.proc,
    test.proc,
    spectra.co2_bands,
    spectra.wn,
    spectra.nm,
    spectra.index,
    target.vars,
    target.bc.vars,
    target.bc.lambda,
    file="02_AfricaProcessedData.Rdata")






