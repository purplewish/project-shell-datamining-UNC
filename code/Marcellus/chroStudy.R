######################################################################################################################
# Chronological Study
######################################################################################################################

#---------------------------------------------------------------------------------------------------------------------
### Setup
#---------------------------------------------------------------------------------------------------------------------
# code path
setwd("Z:/GitHup/project-shell-datamining-UNC/code/Marcellus")

source("header.R")
source("plotFuns.R")
source("loadData.R")

# Results directory
setwd(file.path(repo_path, "Code/results/marcellus"))



#---------------------------------------------------------------------------------------------------------------------
### Chronological Study
#---------------------------------------------------------------------------------------------------------------------
#@@ Heatmap of prediction

# generate grid point for the production prediction
d <- a %>% select(Longitude, Latitude)
grid <- grid4HeatmapProd(d)
write.csv(grid, file="grid.csv", row.names=FALSE)
grid <- read.csv("grid.csv", as.is=T)

# True production
d <- a %>% select(Longitude, Latitude, Norm.Lat.12.Month.Gas) %>% rename(Production=Norm.Lat.12.Month.Gas)

# plotWellProd(d)
distmat <- rdist(grid, d[,1:2])
nnTarget <- idw(d[,3], distmat, 2, 100)
grid.pred <- cbind(grid, Production=nnTarget)

plotHeatmapProd(grid.pred, d[,1:2], long.range=c(-82, -75.5), lat.range=c(37.6, 42), time="Through Q1 2015")



q <- quantile(d$Production, probs = c(0.25, 0.75, 0.8, 0.9, 0.95, 0.975, 0.99, 0.995, 0.998))
top.q <- d[d$Production>=q[2],] 
bot.q <- d[d$Production<=q[1],] 

plotWellProd(top.q) 
plotWellProd(bot.q) 

