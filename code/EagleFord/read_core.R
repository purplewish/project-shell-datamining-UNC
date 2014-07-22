#################################################################################################
## Read in data
setwd("C:/Apps/projects/Corrosion/")

d <- read.csv("C:/Apps/projects/DataMiningUNC/kaggle/02_Input data/Eagle Ford/", header=T)

a <- read.csv("./data/dataAdd.csv", header=T)
a <- a[a[,1] <= max(d[,1]), ] # T<=300
T  <- a[, 1]
pw <- a[, 2]
pm <- a[, 3]
w  <- a[, 4]