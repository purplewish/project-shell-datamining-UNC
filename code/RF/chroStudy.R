######################################################################################################################
# Chronological Study
######################################################################################################################

#---------------------------------------------------------------------------------------------------------------------
### Setup
#---------------------------------------------------------------------------------------------------------------------
# code path
setwd("Z:/GitHup/project-shell-datamining-UNC/code/RF")

source("header.R")
source("runRF.R")
source("plotFuns.R")
source("loadChronIHS1st12monData.R")

# Results directory
setwd(file.path(repo_path, "Code/RF/results"))



#---------------------------------------------------------------------------------------------------------------------
### Chronological Study
#---------------------------------------------------------------------------------------------------------------------
#@@ First 12 month IHS & Core Lab data by different cutoff dates
# RF using Loc + Cov + completion par (m=8)
# Spearman correlation; Top Q1 / Bottom Q1

#@@ RF model
# 1. X w/ loc
# 2. X w/o loc
spearman.cor <- NULL
q1.top.rec <- NULL
q1.bot.rec <- NULL

for (i in 1:nrow(cut.info)){
  
  cut <- as.POSIXlt(as.Date(cut.info$cut.off.date[i]))
  cut$year <- cut$year-1
  cut <- as.Date(cut)
  
  # data of a particular cutoff date
  d <- chro.dat %>% filter(pct.IHS==cut.info$percentage.IHS.used[i])      
  
  # Traing RF model and predict wells after cutoff date
  #set.seed(777)
  set.seed(222) # only for m=6
  sol <- runRFReg4(d, cutoff=cut, model=formula.reg.chro.noloc, m=6, no.tree=1000, ntrace=500)
  #sol <- runRFReg4(d, cutoff=cut, model=formula.reg.chro.geo, m=10, no.tree=1000, ntrace=500)

  x <- sol[[2]] %>% select(-UWI) %>% rename(Target=Norm.Lat.12.Month.Liquid, RF=Pred)

  spearman.cor <- c(spearman.cor, cor(x[,1], x[,2], method="spearman") )  # spearman correlation
  
  q1.top.rec0 <- qRecCurv(x) * 100
  q1.top.rec <- c(q1.top.rec, q1.top.rec0[ceiling(nrow(x) * 0.25),2])  # top quartile
  
  q1.bot.rec0 <- qRecCurvBottom(x) * 100
  q1.bot.rec <- c(q1.bot.rec, q1.bot.rec0[ceiling(nrow(x) * 0.25),2])  # bottorm quartile
}

spearman.cor  <- data.frame(cutoff=cut.info$cut.off.date, pct.IHS=cut.info$percentage.IHS.used, pct.Core=cut.info$precentage.CoreLabs.used, spearman.cor=spearman.cor)
q1.top.rec.rf <- data.frame(cutoff=cut.info$cut.off.date, pct.IHS=cut.info$percentage.IHS.used, pct.Core=cut.info$precentage.CoreLabs.used, q1.top.rf=q1.top.rec)
q1.bot.rec.rf <- data.frame(cutoff=cut.info$cut.off.date, pct.IHS=cut.info$percentage.IHS.used, pct.Core=cut.info$precentage.CoreLabs.used, q1.bot.rf=q1.bot.rec)

# saveRDS(spearman.cor,  "first12mon_q1_spearman_cor_rf_reg_diff_cutoff_noloc_m6.rds")  
# saveRDS(q1.top.rec.rf, "first12mon_q1_top_rec_rf_reg_diff_cutoff_noloc_m6.rds")  
# saveRDS(q1.bot.rec.rf, "first12mon_q1_bot_rec_rf_reg_diff_cutoff_noloc_m6.rds")


#@@ Kriging model
spearman.cor <- NULL
q1.top.rec <- NULL
q1.bot.rec <- NULL

for (i in 1:nrow(cut.info)){
  cut <- as.POSIXlt(as.Date(cut.info$cut.off.date[i]))
  cut$year <- cut$year-1
  cut <- as.Date(cut)
  
  # data of a particular cutoff date
  d <- chro.dat %>% filter(pct.IHS==cut.info$percentage.IHS.used[i])   
  d <- d[d$Date.Production.Start>cut, ]  # out of sample data (test data)
  
  x <- d %>% select(Norm.Lat.12.Month.Liquid, Kriged.Production) %>% rename(Target=Norm.Lat.12.Month.Liquid, kriged=Kriged.Production)
  
  spearman.cor <- c(spearman.cor, cor(x[,1], x[,2], method="spearman") )  # spearman correlation
  
  q1.top.rec0 <- qRecCurv(x) * 100
  q1.top.rec <- c(q1.top.rec, q1.top.rec0[ceiling(nrow(x) * 0.25),2])  # top quartile
  
  q1.bot.rec0 <- qRecCurvBottom(x) * 100
  q1.bot.rec <- c(q1.bot.rec, q1.bot.rec0[ceiling(nrow(x) * 0.25),2])  # bottorm quartile
}

spearman.cor  <- data.frame(cutoff=cut.info$cut.off.date, pct.IHS=cut.info$percentage.IHS.used, pct.Core=cut.info$precentage.CoreLabs.used, spearman.cor.kriged=spearman.cor)
q1.top.rec.kriged <- data.frame(cutoff=cut.info$cut.off.date, pct.IHS=cut.info$percentage.IHS.used, pct.Core=cut.info$precentage.CoreLabs.used, q1.top.kriged=q1.top.rec)
q1.bot.rec.kriged <- data.frame(cutoff=cut.info$cut.off.date, pct.IHS=cut.info$percentage.IHS.used, pct.Core=cut.info$precentage.CoreLabs.used, q1.bot.kriged=q1.bot.rec)



#@@ Chronological comparison plots RF(noloc,m=6) vs. Kriged
spearman.cor.rf <- readRDS("first12mon_q1_spearman_cor_rf_reg_diff_cutoff_noloc_m6.rds")  
q1.top.rf <- readRDS("first12mon_q1_top_rec_rf_reg_diff_cutoff_noloc_m6.rds")  
q1.bot.rf <- readRDS("first12mon_q1_bot_rec_rf_reg_diff_cutoff_noloc_m6.rds")

spearman.cor.kriged <- readRDS("first12mon_q1_spearman_cor_kriged_diff_cutoff.rds")  
q1.top.kriged <- readRDS("first12mon_q1_top_rec_kriged_diff_cutoff.rds")
q1.bot.kriged <- readRDS("first12mon_q1_bot_rec_kriged_diff_cutoff.rds")


# spearman correlation
sp <- cbind(spearman.cor.rf, spearman.cor.kriged=spearman.cor.kriged[,4])
sp2 <- sp[c(1,4,7,10,13),]
sp3 <- sp[c(1,3,5,7,9,11,13),]

# top quartile
q1.top <- cbind(q1.top.rf, q1.top.kriged=q1.top.kriged[,4])
q1.top2 <- q1.top[c(1,4,7,10,13),]
q1.top3 <- q1.top[c(1,3,5,7,9,11,13),]

# bot quartile
q1.bot <- cbind(q1.bot.rf, q1.bot.kriged=q1.bot.kriged[,4])
q1.bot2 <- q1.bot[c(1,4,7,10,13),]
q1.bot3 <- q1.bot[c(1,3,5,7,9,11,13),]


#@@spearman correlation plot y(0.47, 0.7)
dat <- rbind(data.frame(cut=sp2[,1], sp=sp2[,4], method="Random Forest"), 
             data.frame(cut=sp2[,1], sp=sp2[,5], method="Kriged")
            )
plotMLine(dat,"", "Spearman Correlation", NULL, ylim=c(0.4,0.7))

#@@top q1 recovery rate y(50, 66)
dat <- rbind(data.frame(cut=q1.top2[,1], sp=q1.top2[,4], method="Random Forest"), 
             data.frame(cut=q1.top2[,1], sp=q1.top2[,5], method="Kriged")
)
plotMLine(dat,"", "Top Quartile Recovery Rate", NULL, c(45, 66))

#@@bot q1 recovery rate y(34, 60)
dat <- rbind(data.frame(cut=q1.bot2[,1], sp=q1.bot2[,4], method="Random Forest"), 
             data.frame(cut=q1.bot2[,1], sp=q1.bot2[,5], method="Kriged")
)
plotMLine(dat,"", "Bottom Quartile Recovery Rate", NULL, c(34,60))



#@@ Heatmap of prediction

# generate grid point for the production prediction
# d <- chro.dat %>% filter(pct.IHS==cut.info$percentage.IHS.used[1]) %>% select(Latitude, Longitude)
# grid <- grid4HeatmapProd(d)
# write.csv(grid, file = "grid.csv",row.names=FALSE)

grid <- read.csv("grid.csv", as.is=T)

pred <- NULL
for (i in c(1,4,7,10,13)){
  
  cut <- as.POSIXlt(as.Date(cut.info$cut.off.date[i]))
  cut$year <- cut$year-1
  cut <- as.Date(cut)
  
  # Training data of a particular cutoff date
  d <- chro.dat %>% filter(pct.IHS==cut.info$percentage.IHS.used[i])
  train <- d[d$Date.Production.Start<=cut, ]
  test  <- dplyr::setdiff(d, train)
  
  # Traing RF model
  set.seed(222) # only for m=6
  rf.model <- randomForest(formula.reg.chro.noloc, data=train, mtry=6, do.trace=500, ntree=1000)  
  
  test.pred <- cbind(test[,c(3,4)], Production=predict(rf.model, newdata=test))  # Uwi, Target, Pred
  train.obs <- train %>% select(Latitude, Longitude, Norm.Lat.12.Month.Liquid) %>% rename(Production=Norm.Lat.12.Month.Liquid)
  both <- rbind(test.pred, train.obs)
    
  distmat <- rdist(grid, both[,2:1])
  
  nnTarget <- idw(both[,3], distmat, 2, 40)  # nearest neighbor est
    
# grid prediction
#   d <- grid.dat %>% filter(pct.IHS==cut.info$percentage.IHS.used[i]) %>% select(-pct.IHS, -pct.Core)  
#   grid.pred <- cbind(grid, Production=predict(rf.model, newdata=d))
 
  grid.pred <- cbind(grid, Production=nnTarget)
  pred <- c(pred, list(grid.pred))
}

j <- 0
for (i in c(1,4,7,10,13)){
  j <- j+1
  plotHeatmapProd(pred[[j]], long.range=c(-100.6, -97), lat.range=c(27.6, 30.3), time=cut.info$cut.off.date[i])
}

# True production
d <- chro.dat %>% filter(pct.IHS==cut.info$percentage.IHS.used[i])
d <- d %>% select( Longitude, Latitude, Norm.Lat.12.Month.Liquid) %>% rename(Production=Norm.Lat.12.Month.Liquid)

#plotWellProd(d)
distmat <- rdist(grid, d[,1:2])
nnTarget <- idw(d[,3], distmat, 2, 40)
grid.pred <- cbind(grid, Production=nnTarget)

plotHeatmapProd(grid.pred, long.range=c(-100.6, -97), lat.range=c(27.6, 30.3), time="True Production")




#@@ Heatmap of prediction based soly on geology pars

grid <- read.csv("grid.csv", as.is=T)

pred <- NULL
for (i in c(1,4,7,10,13)){
  
  cut <- as.POSIXlt(as.Date(cut.info$cut.off.date[i]))
  cut$year <- cut$year-1
  cut <- as.Date(cut)
  
  # Training data of a particular cutoff date
  d <- chro.dat %>% filter(pct.IHS==cut.info$percentage.IHS.used[i])
  train <- d[d$Date.Production.Start<=cut, ]
  test  <- dplyr::setdiff(d, train)
  
  # Traing RF model
  set.seed(222) # only for m=6
  rf.model <- randomForest(formula.reg.chro.geo, data=train, mtry=6, do.trace=500, ntree=1000)  
  
  test.pred <- cbind(test[,c(3,4)], Production=predict(rf.model, newdata=test))  # Uwi, Target, Pred
  train.obs <- train %>% select(Latitude, Longitude, Norm.Lat.12.Month.Liquid) %>% rename(Production=Norm.Lat.12.Month.Liquid)
  both <- rbind(test.pred, train.obs)
  
  distmat <- rdist(grid, both[,2:1])
  
  nnTarget <- idw(both[,3], distmat, 2, 40)  # nearest neighbor est
  
  # grid prediction
  #   d <- grid.dat %>% filter(pct.IHS==cut.info$percentage.IHS.used[i]) %>% select(-pct.IHS, -pct.Core)  
  #   grid.pred <- cbind(grid, Production=predict(rf.model, newdata=d))
  
  grid.pred <- cbind(grid, Production=nnTarget)
  pred <- c(pred, list(grid.pred))
}

j <- 0
for (i in c(1,4,7,10,13)){
  j <- j+1
  plotHeatmapProd(pred[[j]], long.range=c(-100.6, -97), lat.range=c(27.6, 30.3), time=cut.info$cut.off.date[i])
}






