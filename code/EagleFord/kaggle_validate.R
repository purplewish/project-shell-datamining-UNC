

options(scipen=999)  # remove scientific notation
require("sqldf")

# Set data directory -----------------------------------------------------------------------------------
path.old <- getwd()
path <- "C:/Apps/projects/DataMiningUNC/kaggle/Final/Documentation and Input Files September 22 2014/Documentation and Input Files/"  # data path
setwd(path)



# EF Oil ----------------------------------------------------------------------------------------------
oil.in <- read.csv("EagleFordOilInput.csv", header=T, as.is=T)
oil.out <- read.csv("EagleFordOilOutput.csv", header=T, as.is=T)
oil.out <- cbind(oil.out, mean=rowMeans(oil.out[,4:6]), diff=0)
oil.all <- cbind(oil.in, oil.out[,4:7])
#length(unique(oil.all$Uwi))  # tot=5368 wells =5222 unique + 146 duplicates
oil.all.uniq <- oil.all[!duplicated(oil.all[,1]),]  # rm duplicates
num.col <- ncol(oil.all)
write.table(oil.all.uniq[,c(1, num.col-5, num.col)], file="C:/Apps/projects/DataMiningUNC/data/analysis/EagleFordKaggleOilUniq.csv", row.names=F, sep=",")

# Find duplicate Uwi
dup <- oil.all[duplicated(oil.all[,1]),]
write.table(dup[,1], file="C:/Apps/projects/DataMiningUNC/data/analysis/EagleFordKaggleOilDuplicateUwi.csv", row.names=F, col.names="Duplicate.Uwi")


#oil_out[oil_out[,8]==max(oil_out[,8]),]
#oil_out[oil_out[,8]==min(oil_out[,8]),]
#oil_all <- sqldf("Select * From oil_in a Left Join oil_out b on (a.Uwi = b.Uwi)")



# Analysis ----------------------------------------------------------------------------------------------
# Kaggle's model results
path <- "C:/Apps/projects/DataMiningUNC/data/analysis/"  # analysis path
setwd(path)

true <- read.csv("EagleFordRawOilUniq.csv", header=T, as.is=T)  #3475
pred <- read.csv("EagleFordKaggleOilUniq.csv", header=T, as.is=T)  # 5222

# Find pred Uwi in raw data (true production)
pred_true <- sqldf("Select * from pred a, true b where a.Uwi = b.API")
pred.true <- sqldf("Select Uwi, mean, 1.0*oil_12m/Producer_EstimatedLength_Joined from pred_true")
pred.true <- cbind(pred.true, diff=pred.true[,3]-pred.true[,2])  # diff min=-20.8 max=464

top.q.num <- round(nrow(pred.true)*0.25)  #780

pred.top.q <- pred.true[order(pred.true[,2], decreasing=T),]
pred.top.q <- pred.top.q[1:top.q.num, 1:2]
pred.top.uwi <- sort(pred.top.q[,1])  # pred top 25% uwi

true.top.q <- pred.true[order(pred.true[,3], decreasing=T),]
true.top.q <- true.top.q[1:top.q.num, c(1,3)]
true.top.uwi <- sort(true.top.q[,1])  # true top 25% uwi

overlap.true.pred <- intersect(pred.top.uwi, true.top.uwi)  #398


# Rule-based approach
names(oil.in)
oil.toc <- oil.in$Core.Toc.Kriged
summary(oil.toc)

oil.mat<-oil.in$Core.RoCalculated.Kriged
summary(oil.mat)

oil.por<-oil.in$Core.GriTotalPorosity.Kriged
summary(oil.por)



summary(gas.in)






