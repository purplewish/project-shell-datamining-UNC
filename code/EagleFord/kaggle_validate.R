

options(scipen=999)  # remove scientific notation
require("sqldf")
require("ggplot2")
require('dplyr')

# Set data directory -----------------------------------------------------------------------------------
path.old <- getwd()
path <- "C:/Apps/projects/DataMiningUNC/kaggle/Final/Documentation and Input Files September 22 2014/Documentation and Input Files/"  # data path
setwd(path)


# EF oil data -------------------------------------------------------------------------------------------
oil.in <- read.csv("EagleFordOilInput.csv", header=T, as.is=T)
oil.out <- read.csv("EagleFordOilOutput.csv", header=T, as.is=T)
oil.out <- cbind(oil.out, mean=rowMeans(oil.out[,4:6]))
oil.all <- cbind(oil.in, oil.out[,4:7])
length(unique(oil.all$Uwi))  # tot=5368 = 5222 unique + 146 duplicates
oil.all.uniq <- oil.all[!duplicated(oil.all[,1]),]  # rm duplicates
num.col <- ncol(oil.all)
##write.table(oil.all.uniq[,c(1, num.col-5, (num.col-3):num.col)], file="C:/Apps/projects/DataMiningUNC/data/analysis/EagleFordKaggleOilUniq.csv", row.names=F, sep=",")

# Find duplicate Uwi
# dup <- oil.all[duplicated(oil.all[,1]),]
##write.table(dup[,1], file="C:/Apps/projects/DataMiningUNC/data/analysis/EagleFordKaggleOilDuplicateUwi.csv", row.names=F, col.names="Duplicate.Uwi")


# Analysis ----------------------------------------------------------------------------------------------
# Kaggle's model results
#path <- "C:/Apps/projects/DataMiningUNC/data/analysis/"  # analysis path
path <- "Z:/project/DataMiningUNC/Data/analysis"
setwd(path)

true <- read.csv("EagleFordRawOilUniq.csv", header=T, as.is=T)  #3494 Note: use aggregate datatab (Prod_IP_Cum_Norm)
#true.all <- read.csv("EagleFordRawOilGasUniq.csv", header=T, as.is=T)  #5814
pred <- read.csv("EagleFordKaggleOilUniq.csv", header=T, as.is=T)  # 5222

# Find pred Uwi does not have 12month raw production data
#true_all <- true.all
#pred.notin.true.all <- sqldf("Select * from pred a where a.Uwi not in (select API from true_all)")  # 193 Note: these ID maybe avaialbe in monthly prod tab
##write.table(pred.notin.true.all, file="./EagleFord.PredUwi.Notin.RawProd.csv", row.names=F, sep=",")


# Find pred Uwi in raw data (true production)
pred_true <- sqldf("Select * from pred a, true b where a.Uwi = b.API")  #3137 UWI has both pred and true value
#pred_true <- cbind(ID=seq(nrow(pred_true)), pred_true, oil_12m_norm=pred_true$oil_12m/pred_true$Producer_EstimatedLength_Joined)

#ggplot(pred_true, aes(x=ID, y=oil_12m_norm)) + geom_line()
#ggplot(pred_true, aes(x=ID, y=oil_12m)) + geom_bar(stat="identity")

#pred_true_ord <- pred_true[order(pred_true$oil_12m_norm),]
#pred_true_ord$ID <- seq(nrow(pred_true))
#ggplot(pred_true_ord, aes(x=ID, y=oil_12m_norm)) + geom_line()
#ggplot(pred_true_ord[2900:3000,], aes(x=ID, y=oil_12m_norm)) + geom_ribbon(aes(ymin=Production_P90,  ymax=Production_P10), alpha=0.3) + geom_line()
#ggplot(pred_true_ord, aes(x=ID, y=oil_12m)) + geom_bar(stat="identity")


pred.true <- sqldf("Select Uwi, mean, 1.0*oil_12m/Producer_EstimatedLength_Joined, Surface_Latitude,  Surface_Longitude from pred_true")  # 3137
names(pred.true)[2] <- "pred_norm"
names(pred.true)[3] <- "true_prod_norm"
names(pred.true)[4] <- "Latitude"
names(pred.true)[5] <- "Longitude"

#pred.true <- cbind(pred.true, diff=pred.true[,3]-pred.true[,2])  # diff min=-20.8 max=464

#pred.notin.true <- sqldf("Select * from pred a where a.Uwi not in (select API from true)")  # 2100
#write.table(pred.notin.true, file="./EagleFord.PredUwi.Notin.RawOilProd.csv", row.names=F, sep=",")

top.q.num <- round(nrow(pred.true)*0.25)  #784

pred.top.q <- pred.true[order(pred.true[,2], decreasing=T),]
pred.top.q <- pred.top.q[1:top.q.num, c(-2,-3)]

pred.top.uwi <- sort(pred.top.q[,1])  # pred top 25% uwi

true.top.q <- pred.true[order(pred.true[,3], decreasing=T),]
true.top.q <- true.top.q[1:top.q.num, c(-2,-3)]
true.top.uwi <- sort(true.top.q[,1])  # true top 25% uwi

overlap.true.pred <- intersect(pred.top.uwi, true.top.uwi)  #400 ~51% overlap

pred.top.q <- cbind(pred.top.q, type=ifelse(pred.top.q$Uwi %in% overlap.true.pred, "Overlap", "Pred") )
#false.pos <- pred.top.q[pred.top.q[,4]=="Pred",]  
#dim(false.pos)  #384
true.top.q <- cbind(true.top.q, type=ifelse(true.top.q$Uwi %in% overlap.true.pred, "Overlap", "True") )


all <- rbind(pred.top.q, true.top.q)  #1568
all.nodup <- all[!duplicated(all[,1]),]  #1168

top.q.plot <- ggplot(all.nodup, aes(x=Longitude, y=Latitude,  color=type)) + geom_point(size=4)
p<-top.q.plot + theme(
                   #legend.title=element_blank(),
                   #legend.text = element_text(size = 20),
                   #legend.position=c(.1, .9),
                   legend.position="none",
                   #panel.background = element_blank(),
                   axis.title.x = element_text(size=28),
                   axis.title.y = element_text(size=28)
                   )
p + ylim(26,30.5)

          
          
ggplot(pred.top.q, aes(x=Surface_Latitude, y=Surface_Longitude)) + geom_point()
ggplot(true.top.q, aes(x=Surface_Latitude, y=Surface_Longitude)) + geom_point()



# Rule-based approach-----------------------------------------------------------------
path <- "C:/Apps/projects/DataMiningUNC/kaggle/Final/RulesBasedApproach/RulesBasedApproach/"
setwd(path)

rule.oil.in <- read.csv("Rules features as of Jan 2014.csv", header=T, as.is=T)
a <- rule.oil.in
a[is.na(a)] <- 999999

a$TOC.S <- 100
a[a[,2]<=8, 13] <- 50
a[a[,2]<=6, 13] <- 0
a[a[,2]<=4, 13] <- -50
a[a[,2]<=3, 13] <- -100
a[a[,2]==999999, 13] <- -100

a$RoMeasured.S    <- 100
a[a[,3]<=1, 14]   <- 50
a[a[,3]<=0.9, 14] <- 0
a[a[,3]<=0.8, 14] <- -50
a[a[,3]<=0.75, 14]<- -100
a[a[,3]==999999, 14] <- -100

a$oilApiGravity.S<- 100
a[a[,4]<=40, 15] <- 50
a[a[,4]<=35, 15] <- 0
a[a[,4]<=30, 15] <- -50
a[a[,4]<=25, 15] <- -100
a[a[,4]==999999, 15] <- -100

a$GriSaturationsSo.S<- 100
a[(a[,5]/100)<=0.7, 16] <- 50
a[(a[,5]/100)<=0.6, 16] <- 0
a[(a[,5]/100)<=0.5, 16] <- -50
a[(a[,5]/100)<=0.4, 16] <- -100
a[a[,5]==999999, 16] <- -100

a$ln_formationThickness.S<- 100
a[(exp(a[,6]))<=150, 17] <- 50
a[(exp(a[,6]))<=100, 17] <- 0
a[(exp(a[,6]))<=50, 17] <- -50
a[(exp(a[,6]))<=15, 17] <- -100
a[a[,6]==999999, 17] <- -100

a$GriTotalPorosity.S<- 100
a[a[,7]<=10,18] <- 50
a[a[,7]<=8, 18] <- 0
a[a[,7]<=6, 18] <- -50
a[a[,7]<=4, 18] <- -100
a[a[,7]==999999, 18] <- -100

a$XrdQuartz.S    <- 100
a[a[,8]<=70,19]  <- 50
a[a[,8]<=60, 19] <- 0
a[a[,8]<=50, 19] <- -50
a[a[,8]<=30, 19] <- -100
a[a[,8]==999999, 19] <- -100

a$mudPressure.S    <- 100
a[a[,9]<=0.7, 20]   <- 50
a[a[,9]<=0.5, 20]  <- 0
a[a[,9]<=0.475,20]<- -50
a[a[,9]<=0.45, 20] <- -100
a[a[,9]==999999, 20] <- -100


rule.s <- a[, 13:20]
a$Score <- rowSums(rule.s)
#diff.s <- a[,10] - a[,21]
#a[1:2, ]

c <- a[,c(1, 11, 12, 21)]
perf <- NULL
pred.n <- NULL  # # of overlap between pred and true
rule.n <- NULL  # # of overlap between rule and true

path <- "C:/Apps/projects/DataMiningUNC/data/analysis/"  # analysis path
setwd(path)
#write.table(b, file="./EagleFordOilRule.csv", row.names=F, sep=",")

n.top <- round(nrow(a) * 0.25)  # 658

top.q.true <- c[order(c[,3], decreasing=T),]
top.q.true <- top.q.true[1:n.top, ]
true.id <- as.data.frame(top.q.true[,1])
names(true.id)[1] <- "Top.Q.True"
true.id <- sort(true.id[,1])
#write.table(true.id, file="./EagleFordOilTrueQ1.csv", row.names=F, sep=",")  # true Q1

top.q.pred <- c[order(c[,2], decreasing=T),]
top.q.pred <- top.q.pred[1:n.top, ]
pred.id <- as.data.frame(top.q.pred[,1])
names(pred.id)[1] <- "Top.Q.Pred"
pred.id <- sort(pred.id[,1])
#write.table(pred.id, file="./EagleFordOilPredQ1.csv", row.names=F, sep=",")  # pred Q1


for(i in c(1:1000))
{

  b <- c[sample(1:nrow(a)), ]

  top.q.rule <- b[order(b[,4], decreasing=T),]
  top.q.rule <- top.q.rule[1:n.top, ]
  rule.id <- as.data.frame(top.q.rule[,1])
  names(rule.id)[1] <- "Top.Q.Rule"
  #write.table(rule.id, file="./EagleFordOilRuleQ1.csv", row.names=F, sep=",")  # rule based Q1


  rule.id <- sort(rule.id[,1])


  # overlap between rule and true
  overlap.true.rule <- intersect(true.id, rule.id)  #177/658 ~27% overlap
  overlap.true.pred <- intersect(true.id, pred.id)  #366/658 ~56% overlap  x2.07 times better than rule based results

  perf <- c(perf, length(overlap.true.pred)/length(overlap.true.rule) )
  pred.n <- c(pred.n, length(overlap.true.pred))
  rule.n <- c(rule.n, length(overlap.true.rule))

}

mean(perf)  # 1.63 Kaggle's model/rule
mean(pred.n)
mean(rule.n)

# ------------------------------------------------------------------------------------------
## Plot of Kaggle's prediction vs. true Q1
path <- "C:/Apps/projects/DataMiningUNC/kaggle/Final/RulesBasedApproach/RulesBasedApproach/"
setwd(path)

rule.oil.in <- read.csv("Rules features as of Jan 2014.csv", header=T, as.is=T)
kag <- rule.oil.in[, c(1,11,12)]  # kaggle's prediction and true value 2631

# Load locations of wells
path <- "C:/Apps/projects/DataMiningUNC/data/analysis/"  # analysis path
setwd(path)
loc <- read.csv("EagleFordProdLoc.csv", header=T, as.is=T)

pred.true <- sqldf("Select * from kag a, loc b where a.Uwi = b.API")  
b <- pred.true[,-4]

n.q1 <- round(nrow(kag) * 0.25)  # 658

q1.true <- b[order(b[,3], decreasing=T),]
q1.true <- q1.true[1:n.q1, ]
true.id <- as.data.frame(q1.true[,1])
names(true.id)[1] <- "Q1.True"
#write.table(true.id, file="./EagleFordOilTrueQ1.csv", row.names=F, sep=",")  # true Q1


q1.pred <- b[order(b[,2], decreasing=T),]
q1.pred <- q1.pred[1:n.q1, ]
pred.id <- as.data.frame(q1.pred[,1])
names(pred.id)[1] <- "Q1.Pred"
#write.table(pred.id, file="./EagleFordOilKaggleQ1.csv", row.names=F, sep=",")  # pred Q1

overlap.true.pred <- intersect(true.id[,1], pred.id[,1])  # 366 overlap 366/658=0.556

q1.true <- cbind(q1.true, type=ifelse(q1.true$Uwi %in% overlap.true.pred, "Overlap", "True") )
q1.pred <- cbind(q1.pred, type=ifelse(q1.pred$Uwi %in% overlap.true.pred, "Overlap", "Pred") )

all <- rbind(q1.true, q1.pred)  # 1316
all.nodup <- all[!duplicated(all[,1]),]  # 950

q1.plot <- ggplot(all.nodup, aes(x=Longitude, y=Latitude, color=type)) + geom_point(size=4)
q1.plot + theme(
  #legend.title=element_blank(),
  #legend.text = element_text(size = 20),
  #legend.position=c(.1, .9),
  legend.position="none",
  #panel.background = element_blank(),
  axis.title.x = element_text(size=28),
  axis.title.y = element_text(size=28)
)










