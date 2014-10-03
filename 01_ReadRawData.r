##------------------------------------------------------------------
## Load the raw data files from the Kaggle website into Rdata files
## http://www.kaggle.com/c/afsis-soil-properties
##------------------------------------------------------------------

##------------------------------------------------------------------
## Load libraries
##------------------------------------------------------------------
library(data.table)

##------------------------------------------------------------------
## Clear the workspace
##------------------------------------------------------------------
rm(list=ls())

##------------------------------------------------------------------
## Set the working directory
##------------------------------------------------------------------
setwd("/Users/alexstephens/Development/kaggle/africa/data/raw")

##------------------------------------------------------------------
## Read the test and training datasets
##------------------------------------------------------------------

## read the raw data
test.raw    <- read.csv("test.csv", header=TRUE)
train.raw   <- read.csv("training.csv", header=TRUE)

## lowercase the headers
colnames(test.raw)  <- tolower(colnames(test.raw))
colnames(train.raw) <- tolower(colnames(train.raw))

##------------------------------------------------------------------
## Save results as a data.table
##------------------------------------------------------------------

## change to a new data directory
setwd("/Users/alexstephens/Development/kaggle/africa/data/proc")

## save results to separate files given the size
save(train.raw, test.raw,  file="01_AfricaRawData.Rdata")

