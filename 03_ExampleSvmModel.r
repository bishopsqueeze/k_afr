##------------------------------------------------------------------
## Notes on the data:
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
## Load processed data
##------------------------------------------------------------------
load("02_AfricaProcessedData.Rdata")

##------------------------------------------------------------------
## Define target variables
##------------------------------------------------------------------
target.vars     <- c("ca.bc","p.bc","ph.bc","soc.bc","sand.bc")
nonspec.vars    <- c("bsan","bsas","bsav","cti","elev","evi","lstd","lstn","ref1","ref2","ref3","ref7","reli","tmap","tmfi","depth")
spec.vars       <- colnames(train.proc)[grep("^m", colnames(train.proc))]
der1.vars       <- colnames(train.proc)[grep("^d1", colnames(train.proc))]
der2.vars       <- colnames(train.proc)[grep("^d2", colnames(train.proc))]
id.vars         <- c("pidn", "id")

##******************************************************************
## Simple beat-the-benchmark model
##******************************************************************

## define indices
target.idx      <- which(names(train.proc) %in% target.vars)
predictor.idx   <- which(!(names(train.proc) %in% c(target.vars, id.vars)))


##------------------------------------------------------------------
## Preprocess chunks of the data
##------------------------------------------------------------------

nonspec.x       <- train.proc[, c(nonspec.vars)]
preProcValues   <- preProcess(nonspec.x, method = c("center", "scale"))
nonspec.sc      <- predict(preProcValues, nonspec.x)



nonspec.zsd     <- which( apply(nonspec.x, 2, sd) == 0 )
if (length(nonspec.zsd) > 0) {
    nonspec.x <- nonspec.x[, -nonspec.zsd]
}

spec.x       <- train.proc[, c(spec.vars)]
spec.zsd     <- which( apply(spec.x, 2, sd) == 0 )
if (length(spec.zsd) > 0) {
    spec.x <- spec.x[, -spec.zsd]
}
preProcValues   <- preProcess(spec.x, method = c("center", "scale"))
spec.sc      <- predict(preProcValues, spec.x)



der1.x       <- train.proc[, c(der1.vars)]
der1.zsd     <- which( apply(der1.x, 2, sd) == 0 )
if (length(der1.zsd) > 0) {
    der1.x <- der1.x[, -der1.zsd]
}

der2.x       <- train.proc[, c(der2.vars)]
der2.zsd     <- which( apply(der2.x, 2, sd) == 0 )
if (length(der2.zsd) > 0) {
    der2.x <- der2.x[, -der2.zsd]
}



##------------------------------------------------------------------
## Step 1:  Use caret to find the optimal tuning params for each fit
##------------------------------------------------------------------

ctrl        <- trainControl(method = "cv", savePred=FALSE)
svm.grid    <- expand.grid(C = 2^(seq(-1,6,1)))
#svm.grid    <- expand.grid(C = 1)
gbm.grid    <- expand.grid(interaction.depth = c(1, 5, 9, 13),
                        n.trees = 150,
                        shrinkage = 0.1)
tmp.res     <- list()

## loop over each target variable
#for (i in 1:length(target.vars)) {
for (i in 2:2) {
    
    ## load all raw variables
    tmp.x                   <- data.frame(spec.x)
    tmp.y                   <- train.proc[, c(target.vars[i])]
    tmp.data                <- data.frame(tmp.y, tmp.x)
    colnames(tmp.data)[1]   <- target.vars[i]
    
    ## use all variables in the file
    tmp.formula         <- as.formula(paste0(target.vars[i], " ~ ."))
    
    ## save the results
    tmp.res[[i]]        <- train(tmp.formula, data=tmp.data, method = "svmLinear", verbose=FALSE, trControl = ctrl, tuneGrid=svm.grid)
    #tmp.res[[i]]        <- train(tmp.formula, data=tmp.data, method = "svmRadial", verbose=FALSE, trControl = ctrl)
    #tmp.res[[i]]        <- train(tmp.formula, data=tmp.data, method = "gbm", verbose=FALSE, trControl = ctrl, tuneGrid=gbm.grid)

}




## minmum RMSE (all variables; svmLinear)
## [1] = 0.335 +- 0.113, C = 1
## [2] = 0.830 +- 0.355, C = 16
## [3] = 0.326 +- 0.058, C = 8
## [4] = 0.303 +- 0.072, C = 1
## [5] = 0.333 +- 0.049, C = 1

## minmum RMSE (only spectra; svmLinear)
## [1] = 0.323 +- 0.090, C = 1
## [2] = 0.xxx +- 0.xxx, C = 1
## [3] = 0.xxx +- 0.xxx, C = 1
## [4] = 0.xxx +- 0.xxx, C = 1
## [5] = 0.xxx +- 0.xxx, C = 1


## minmum RMSE (only spectra; centered/scaled; svmLinear)
## [1] = 0.333 +- 0.151, C = 1
## [2] = 0.xxx +- 0.xxx, C = 1
## [3] = 0.xxx +- 0.xxx, C = 1
## [4] = 0.xxx +- 0.xxx, C = 1
## [5] = 0.xxx +- 0.xxx, C = 1


## minmum RMSE (only der1; svmLinear)
## [1] = 0.450 +- 0.100, C = 1
## [2] = 0.xxx +- 0.xxx, C = 1
## [3] = 0.xxx +- 0.xxx, C = 1
## [4] = 0.xxx +- 0.xxx, C = 1
## [5] = 0.xxx +- 0.xxx, C = 1

## minmum RMSE (only der2; svmLinear)
## [1] = 0.491 +- 0.142, C = 1
## [2] = 0.xxx +- 0.xxx, C = 1
## [3] = 0.xxx +- 0.xxx, C = 1
## [4] = 0.xxx +- 0.xxx, C = 1
## [5] = 0.xxx +- 0.xxx, C = 1

## minmum RMSE (spectra + der1; svmLinear)
## [1] = 0.461 +- 0.102, C = 1
## [2] = 0.xxx +- 0.xxx, C = 1
## [3] = 0.xxx +- 0.xxx, C = 1
## [4] = 0.xxx +- 0.xxx, C = 1
## [5] = 0.xxx +- 0.xxx, C = 1



## use the propectr package to pre-process some of the spectra


##------------------------------------------------------------------
## Save results
##------------------------------------------------------------------
##save.image("02_AfricaProcessedData.Rdata")