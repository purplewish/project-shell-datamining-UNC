## Header file

rm(list=ls())
options(scipen=999)
#options(java.parameters = "-Xmx10000m")
options(java.parameters = "-Xmx80g")

library(randomForest)
library(ggplot2)
library(reshape2)
library(dplyr)
library(bartMachine)

# path to project folder
#repo_path = "C:/Apps/projects/DataMiningUNC"
repo_path = "Z:/project/DataMiningUNC"

