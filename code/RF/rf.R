
# code path
setwd("Z:/GitHup/project-shell-datamining-UNC/code/RF")
source("header.R")
source("loadData.R")
source("runRF.R")
source("plotFuns.R")

# Results directory
setwd(file.path(repo_path, "Code/RF/results"))



#================================================================================================================================
# RF Classification Approach ###
#================================================================================================================================

## RF model on all data with selected pars
set.seed(777)
rf <- runRF(dat=all, train.pct=1, model=formula.class2, m=5, no.tree=1000, nrep=1)
sol <- rf[[1]]  # Pred accuracy on test data
rf.mod <- rf[[2]]  # Last rf model obj
# saveRDS(rf.mod, "rfMod_sweetspot.rds")
# rf.mod <- readRDS("rfMod_sweetspot.rds")

## Partial plots
#par(mfrow=c(2,2))
#partialPlot(rf.mod, all, Core.Tmax.Kriged, xlab="Tmax", main="")
#partialPlot(rf.mod, all, Producer.DepthTrueVertical.Joined, xlab="True Vertical Depth", main="")
#partialPlot(rf.mod, all, Core.S2.Kriged, xlab="S2", main="")
#partialPlot(rf.mod, all, Core.S3.Kriged, xlab="S3", main="")

## Plot sweet-spots
target <- select(all, Uwi, Target.Q4, Latitude, Longitude)
a <- rf.mod$predicted
b <- as.numeric(names(a))
c <- data.frame(b, a)
d <- arrange(c, b)  # order predicted value from RF model
dat <- cbind(target, Target.Q4.Pred=d[,2])  # Uwi Target.Q4 Latitude Longitude pred

plotSweetspot(dat)


#-------------------------------------------------------------------------------------------------------------------------
## RF model 5-fold CV with ALL vars 
set.seed(777)
rf <- runRFCV(dat=all, model=formula.class2, m=5, no.tree=1000, k=5)
sol <- rf[[1]]  # CV pred accuracy 
pred <- rf[[2]]  # CV pred results

## Plot sweetspots
# saveRDS(sol, "sweetspot_5CV_accy.rds")
# sol <- readRDS("sweetspot_5CV_accy.rds")
# saveRDS(pred, "sweetspot_5CV_pred.rds")
# pred <- readRDS("sweetspot_5CV_pred.rds")
plotSweetspot(pred)

#------------------------------------------

## RF model 5-fold CV with 4 selected vars
formula.class.top10overlap.rmRo <- readRDS("formula_4impvar_p.rds")
set.seed(777)
rf <- runRFCV2(dat=all, model=formula.class.top10overlap.rmRo, no.tree=1000, k=5)  # use default setting m
sol <- rf[[1]]  # CV pred accuracy 
pred <- rf[[2]]  # CV pred results

## Plot sweetspots
# saveRDS(sol, "sweetspot_5CV_accy_4vars.rds")
# sol <- readRDS("sweetspot_5CV_accy_4vars.rds")
# saveRDS(pred, "sweetspot_5CV_pred_4vars.rds")
# pred <- readRDS("sweetspot_5CV_pred_4vars.rds")
plotSweetspot(pred)


#-------------------------------------------------------------------------------------------------------------------------
## RF model on vars importance (avg importance score with rep runs)
rep <- 50
train.pct <- 1
mda <- NULL; mdg <-NULL;

set.seed(777)
for (i in 1:rep){
  rf <- runRF(dat=all, train.pct=train.pct, model=formula.class2, m=5, no.tree=1000)
  rf.mod <- rf[[2]]   # Last rf model obj
  
  Predictor <- rownames(importance(rf.mod))
  #Predictor <- gsub("(Core.|.Kriged|.Joined)",'',rownames(importance(rf.mod)))
  dat <- data.frame(importance(rf.mod))
  mda <- cbind(mda, dat[,3]) # mean decrease in accuracy
  mdg <- cbind(mdg, dat[,4]) # mean decrease in gini
  
  print(paste0("rep=",i), quote=F)
}

dat <- data.frame(Predictor, mda=rowMeans(mda), mdg=rowMeans(mdg))  # average nrep runs importance scores
imp.accy <- dat %>% select(Predictor, mda) %>% arrange(desc(mda))   # mean decrease in accuracy
imp.gini <- dat %>% select(Predictor, mdg) %>% arrange(desc(mdg))   # mean decrease in gini
# saveRDS(imp.accy, "imp_vars_accy_50rep.rds")
# saveRDS(imp.gini, "imp_vars_gini_50rep.rds")
# imp.accy <- readRDS("imp_vars_accy_50rep.rds")
# imp.gini <- readRDS("imp_vars_gini_50rep.rds")

## Plot var importance
imp.accy$Predictor <- gsub("(Core.|.Kriged|.Joined)",'',imp.accy$Predictor)
imp.gini$Predictor <- gsub("(Core.|.Kriged|.Joined)",'',imp.gini$Predictor)
plotRFVarImp3(imp.accy)
plotRFVarImp3(imp.gini)

#----------------------------------------------------
# RF model: variables importance at a fix training % 
set.seed(777)
rf <- runRF(dat=all, train.pct=0.3, model=formula.class2, m=5, no.tree=1000, nrep=1)
rf.mod <- rf[[2]]

plotRFVarImp(rf.mod)
plotRFVarImp2(rf.mod)


#-------------------------------------------------------------------------------------------------------------------------
## RF model 5-fold CV using default mtry setting with top n important vars

imp.accy <- readRDS("imp_vars_accy_50rep.rds")
imp.gini <- readRDS("imp_vars_gini_50rep.rds")
sols <- NULL

set.seed(777)
for (top.n in nrow(imp.accy):1){
  sel.vars.accy <- imp.accy$Predictor[1:top.n]  
  sel.vars.gini <- imp.gini$Predictor[1:top.n]  
  
  formula.class.imp.accy <- formula(paste("Target.Q4~", paste(sel.vars.accy,collapse="+")))
  formula.class.imp.gini <- formula(paste("Target.Q4~", paste(sel.vars.gini,collapse="+"))) # class:Q4 ~Q4, topQ vs. ~topQ
  
  rf.accy <- runRFCV2(dat=all, model=formula.class.imp.accy, no.tree=1000, k=5)  # default mtry setting
  rf.gini <- runRFCV2(dat=all, model=formula.class.imp.gini, no.tree=1000, k=5)
  
  sol.accy <- data.frame(method="accy", topn=top.n, rf.accy[[1]])
  sol.gini <- data.frame(method="gini", topn=top.n, rf.gini[[1]])
  sols <- rbind(sols, sol.accy, sol.gini)
  
  print(paste0("top.n=",top.n), quote=F)
}

# saveRDS(sols, "topn_impvars_5CV_default_mtry.rds")
sols <- readRDS("topn_impvars_5CV_default_mtry.rds")


accy <- sols %>% filter(method=="accy") %>% select(topn, Accy, TP, TN) %>% arrange(topn)
gini <- sols %>% filter(method=="gini") %>% select(topn, Accy, TP, TN) %>% arrange(topn)

plotRFAcc(accy,"Top K Important Variables", "Test Classification Accuracy", xtick=seq(1,31,1))
plotRFAcc(gini,"Top K Important Variables", "Test Classification Accuracy", xtick=seq(1,31,1))


#----------------------------------------------------
## RF model on selected important vars with 5 fold CV
# top n vars overlapped between gini and accy measure
imp.accy <- readRDS("imp_vars_accy_50rep.rds")
imp.gini <- readRDS("imp_vars_gini_50rep.rds")

#-----------------------------
# top 10 and 6 overlapped vars
top.n <-10
sel.accy <- imp.accy$Predictor[1:top.n]
sel.gini <- imp.gini$Predictor[1:top.n]
overlap <- intersect(sel.accy,sel.gini)
#saveRDS(overlap, "top10_impvar_overlap.rds")
#saveRDS(overlap, "top6_impvar_overlap.rds")
top10.overlap <- readRDS("top10_impvar_overlap.rds")
top6.overlap <- readRDS("top6_impvar_overlap.rds")

# top 10
formula.class.top10overlap <- formula(paste("Target.Q4~", paste(top10.overlap,collapse="+"))) # class:Q4 ~Q4, topQ vs. ~topQ
set.seed(777)
rf.top10.ov <- runRFCV2(dat=all, model=formula.class.top10overlap, no.tree=1000, k=5)  # default mtry setting
sol.10ov <- rf.top10.ov[[1]]
# saveRDS(sol.10ov, "top10_impvar_overlap_pred.rds")
sol.10ov <- readRDS("top10_impvar_overlap_pred.rds")

# top 6 
formula.class.top6overlap <- formula(paste("Target.Q4~", paste(top6.overlap,collapse="+"))) # class:Q4 ~Q4, topQ vs. ~topQ
set.seed(777)
rf.top6.ov <- runRFCV2(dat=all, model=formula.class.top6overlap, no.tree=1000, k=5)  # default mtry setting
sol.6ov <-rf.top6.ov[[1]]
# saveRDS(sol.6ov, "top6_impvar_overlap_pred.rds")
sol.6ov <- readRDS("top6_impvar_overlap_pred.rds")


#-----------------------------------------
# rm high corr variables of top 10 overlap
corr.10ov <- readRDS("top10_overlap_corr.rds")
top10.overlap <- readRDS("top10_impvar_overlap.rds")
top10.overlap.rmRo <- top10.overlap[-1]  # rm RoCalculated
top10.overlap.rmRoS3 <- top10.overlap.rmRo [-4]  # rm RoCalculated
top10.overlap.rmRoS3S2 <- top10.overlap.rmRoS3 [-4]  # rm RoCalculated

# rm Ro
formula.class.top10overlap.rmRo <- formula(paste("Target.Q4~", paste(top10.overlap.rmRo,collapse="+"))) # class:Q4 ~Q4, topQ vs. ~topQ
# saveRDS(formula.class.top10overlap.rmRo, "formula_4impvar_p.rds")
set.seed(777)
rf.top10.ov.rmRo <- runRFCV2(dat=all, model=formula.class.top10overlap.rmRo, no.tree=1000, k=5)  # default mtry setting
sol.10ov.rmRo <-rf.top10.ov.rmRo[[1]]


# rm Ro S3
formula.class.top10overlap.rmRoS3 <- formula(paste("Target.Q4~", paste(top10.overlap.rmRoS3,collapse="+"))) # class:Q4 ~Q4, topQ vs. ~topQ
set.seed(777)
rf.top10.ov.rmRoS3 <- runRFCV2(dat=all, model=formula.class.top10overlap.rmRoS3, no.tree=1000, k=5)  # default mtry setting
sol.10ov.rmRoS3 <-rf.top10.ov.rmRoS3[[1]]


#-------------------------------------------------------------------------------------------------------------------------
## RF model on selected important vars with different training pct
nrep=1;
top.n <- 3  # top.n important vars 
train.pct.seq <- seq(0.2,0.8,0.1); 

sols <- NULL
set.seed(777)
for (train.pct in train.pct.seq){
  
  rf <- runRF(dat=all, train.pct=train.pct, model=formula.class2, m=5, no.tree=500, nrep=nrep)
  sol <- data.frame(rf[[1]], method="all")  # Pred accuracy on test data
  rf.mod <- rf[[2]]   # Last rf model obj

  dat <- data.frame(rownames(importance(rf.mod)),round(importance(rf.mod),2))
  names(dat)[c(1,ncol(dat)-1,ncol(dat))] <- c("Predictor","mda","mdg")
  rownames(dat) <- NULL
  pred.accy <- dat %>% select(Predictor, mda) %>% arrange(desc(mda))  # mean decrease in accuracy
  pred.gini <- dat %>% select(Predictor, mdg) %>% arrange(desc(mdg))  # mean decrease in gini

  sel.vars.gini <- pred.gini$Predictor[1:top.n]  
  sel.vars.accy <- pred.accy$Predictor[1:top.n]  
  
  formula.class.imp.gini <- formula(paste("Target.Q4~", paste(sel.vars.gini,collapse="+"))) # class:Q4 ~Q4, topQ vs. ~topQ
  formula.class.imp.accy <- formula(paste("Target.Q4~", paste(sel.vars.accy,collapse="+")))
  
  rf <- runRF(dat=all, train.pct=train.pct, model=formula.class.imp.gini, m=5, no.tree=500, nrep=nrep)
  sol.gini <- data.frame(rf[[1]], method="gini.sel.vars")  # Pred accuracy on test data
  rf.mod.gini <- rf[[2]]   # Last rf model obj
  
  rf <- runRF(dat=all, train.pct=train.pct, model=formula.class.imp.accy, m=5, no.tree=500, nrep=nrep)
  sol.accy <- data.frame(rf[[1]], method="accy.sel.vars")  # Pred accuracy on test data
  rf.mod.accy <- rf[[2]]   # Last rf model obj
  
  sols <- rbind(sols, sol, sol.gini, sol.accy)
}

# saveRDS(sols, "top10_imp_vars.rds")
# sols.top10vars <- readRDS("top10_imp_vars.rds")

# saveRDS(sols, "top5_imp_vars.rds")
#sols.top5vars <- readRDS("top5_imp_vars.rds")
# 
# saveRDS(sols, "top3_imp_vars.rds")
sols.top3vars <- readRDS("top3_imp_vars.rds")

# top10.diff <- setdiff(sel.vars.accy, sel.vars.gini)  # n.top=10
# saveRDS(top10.diff, "top10_imp_var_diff_accy_gini.rds")
# top10.diff.accy.gini <- readRDS("top10_imp_var_diff_accy_gini.rds")

# top10.diff <- setdiff(sel.vars.gini, sel.vars.accy)  # n.top=10
# saveRDS(top10.diff, "top10_imp_var_diff_gini_accy.rds")
# top10.diff.gini.accy <- readRDS("top10_imp_var_diff_gini_accy.rds")

# saveRDS(sel.vars.accy, "top10_accy.rds")
# saveRDS(sel.vars.gini, "top10_gini.rds")
# top10.accy <- readRDS("top10_accy.rds")
# top10.gini <- readRDS("top10_gini.rds")

# overlap <- intersect(sel.vars.accy,sel.vars.gini)
# saveRDS(overlap, "top10_overlap.rds")
#top10.overlap <- readRDS("top10_overlap.rds")  # 5 overlaped
# accy.add <- setdiff(sel.vars.accy, top10.overlap) # overlap Xs rank: 1,2,3,4,8
# gini.add <- setdiff(sel.vars.gini, top10.overlap) # overlap Xs rank: 1,2,6,7,10

## overlapped vars of top 10 important vars(gini and accy); totally 5 vars are overlapped
top10.overlap <- readRDS("top10_overlap.rds")  # 5 overlaped vars
formula.class.top10overlap <- formula(paste("Target.Q4~", paste(top10.overlap,collapse="+"))) # class:Q4 ~Q4, topQ vs. ~topQ

nrep=1;
train.pct.seq <- seq(0.1,0.9,0.1);

sols <- NULL
set.seed(777)
for (train.pct in train.pct.seq){
  
  rf <- runRF(dat=all, train.pct=train.pct, model=formula.class2, m=5, no.tree=1500, nrep=nrep)
  sol <- data.frame(rf[[1]], method="all")  # Pred accuracy on test data
      
  rf <- runRF(dat=all, train.pct=train.pct, model=formula.class.top10overlap, m=5, no.tree=1500, nrep=nrep)
  sol.overlap <- data.frame(rf[[1]], method="top10overlap")  # Pred accuracy on test data
    
  sols <- rbind(sols, sol, sol.overlap)
}

saveRDS(sols, "top10_overlap_vars.rds")
sols.top10overlap <- readRDS("top10_overlap_vars.rds")


#-------------------------------------------------------------------------------------------------------------------------
## Parameter selection
# RF model: effect of number of trees (select no.tree=1000)
num.tree <- 3000

set.seed(123)
rf <- runRF(dat=all, train.pct=1, model=formula.class2, m=3, no.tree=num.tree, nrep=1)
sol <- rf[[1]]  # Pred accuracy on test data
rf.mod <- rf[[2]]  # Last rf model obj

# saveRDS(rf.mod, "rf_ntrees_m3.rds")  # m3 m5 m10
# rf.mod <- readRDS("rf_ntrees_m3.rds")
plotRFOOBErr(rf.mod)  # use oob error, no need to do cross validation


#--------------------------------------
# RF model: effect of mtry (select m=5)
m.seq <- c(3, 5, 10)  # mtry=sqrt(31)=5

set.seed(789)
sol.all <- NULL
for(m in m.seq){
  rf <- runRF(dat=all, train.pct=1, model=formula.class2, m=m, no.tree=1000, nrep=1)
  rf.mod <- rf[[2]]
  sol.all <- rbind(sol.all, c(m=m, rf.mod$err.rate[1000,]))
}
sol.all
#saveRDS(sol.all, "mtry_train_100pct.rds")
#sol.all <- readRDS("mtry_train_100pct.rds")


#-------------------------------------------------------------------------------------------------------------------------------
## RF model: effect of different train % 
# Using all data without CV
train.pct.seq <- seq(0.1,0.9,0.1)

set.seed(151)
sol.all <- NULL
for(train.pct in train.pct.seq){
  rf <- runRF(dat=all, train.pct=train.pct, model=formula.class2, m=5, no.tree=1000, nrep=10)
  sol.all <- rbind(sol.all, rf[[1]])
}
sol.all
# write.csv(sol.all, "./train_perc_errrate_avg.csv", row.names=F)
# saveRDS(sol.all, "train_pct_errrate_avg.rds")
# sol.all <- readRDS("train_pct_errrate_avg.rds")
sol <- sol.all[,2:5]
plotRFAcc(sol)

#-------------------------------------------------------------
# RF model: effect of different train % using CV
K <- c(10, 5, 3, 2, 3, 5, 10)  # K-fold CV
Rev <- c(T, T, T, T, F, F, F)  # Reverse CV or not
Train.Perc <- c(0.1, 0.2, 0.333, 0.5, 0.66, 0.8, 0.9)  # Training pct <-> K-fold CV

set.seed(999)
sol.all <- NULL
for(i in 1:length(K)){
  sol <- runRFCV(dat=all, model=formula.class2, m=5, no.tree=1000, K[i], Rev[i])
  sol.all <- rbind(sol.all, sol[[1]])
}
sol.all
#saveRDS(sol.all, "train_pct_errrate_cv.rds")
sol.all <- readRDS("train_pct_errrate_cv.rds")
sol <- cbind(Train.Perc=Train.Perc, sol.all[,3:5])
plotRFAcc(sol, "Percentage of Training Data", "Test Classification Accuracy", xtick=seq(0.1,0.9,0.1))


#-----------------------------------------------------------------
# RF model: effect of different train % using Cross Validation
# Top 10 overlap gini and accy without RoCalculated (only 4 vars)
formula.class.top10overlap.rmRo <- readRDS("formula_4impvar_p.rds")

K <- c(10, 5, 3, 2, 3, 5, 10)  # K-fold CV
Rev <- c(T, T, T, T, F, F, F)  # Reverse CV or not
Train.Perc <- c(0.1, 0.2, 0.333, 0.5, 0.66, 0.8, 0.9)  # Training pct <-> K-fold CV

set.seed(777)
sol.all <- NULL
for(i in 1:length(K)){
  sol <- runRFCV2(dat=all, model=formula.class.top10overlap.rmRo, no.tree=1000, K[i], Rev[i])  # use default setting m
  sol.all <- rbind(sol.all, sol[[1]])
}
sol.all
#saveRDS(sol.all, "train_pct_errrate_cv_4vars.rds")
sol.all <- readRDS("train_pct_errrate_cv_4vars.rds")
sol <- cbind(Train.Perc=Train.Perc, sol.all[,3:5])
plotRFAcc(sol, "Percentage of Training Data", "Test Classification Accuracy", xtick=seq(0.1,0.9,0.1))



#-------------------------------------------------------------------------------------------------------------------------------------
# Prediction based on timeline (cutoff at time point)

# Cutoff=q1, test well > cutoff+365
q <- summary(all$Date.Production.Start) # check quantiles of production starting date
train <- filter(all, Date.Production.Start < q[2]) # early 25% wells
test <- filter(all, Date.Production.Start > q[2]+365) # ~ 10% wells 

set.seed(777)
rf <- runRF2(train=train, test=test, model=formula.class2, m=5, no.tree=1500)
sol <- rf[[1]]  # Pred accuracy on test data
rf.mod <- rf[[2]]  # rf model obj
plotRFOOBErr(rf.mod)


#-----
# Cutoff=q1, test well = rest
q <- summary(all$Date.Production.Start) # check quantiles of production starting date
train <- filter(all, Date.Production.Start < q[2]) # early 25% wells
test <- setdiff(all, train) # rest wells

set.seed(777)
rf <- runRF2(train=train, test=test, model=formula.class2, m=5, no.tree=500)
sol <- rf[[1]]  # Pred accuracy on test data
rf.mod <- rf[[2]]  # rf model obj
plotRFOOBErr(rf.mod)


# randome sample 25%
train <- sample_frac(all, 0.25, replace=F)
test <- setdiff(all, train)

set.seed(777)
rf <- runRF2(train=train, test=test, model=formula.class2, m=5, no.tree=500)
sol <- rf[[1]]  # Pred accuracy on test data
rf.mod <- rf[[2]]  # rf model obj
plotRFOOBErr(rf.mod)


#-----
# Cutoff=5% 10% 15% 20% ... test well = rest
q <- quantile(as.numeric(all$Date.Production.Start),  probs=c(5, 10, 15, 20, 25, 30, 35, 40, 50)/100)
class(q)="Date"

#5% Accy=0.773
train <- filter(all, Date.Production.Start < q[1]) # early 5% wells
test <- filter(all, Date.Production.Start > q[1]+365) # ~ 10% wells 

set.seed(777)
rf <- runRF2(train=train, test=test, model=formula.class2, m=5, no.tree=500)
sol <- rf[[1]]  # Pred accuracy on test data
rf.mod <- rf[[2]]  # rf model obj
plotRFOOBErr(rf.mod)

#-----
#10% Accy=0.802
train <- filter(all, Date.Production.Start < q[2]) # early 5% wells
test <- filter(all, Date.Production.Start > q[2]+365) # ~ 10% wells 

set.seed(777)
rf <- runRF2(train=train, test=test, model=formula.class2, m=5, no.tree=500)
sol <- rf[[1]]  # Pred accuracy on test data
rf.mod <- rf[[2]]  # rf model obj
plotRFOOBErr(rf.mod)

#-----
#15% Accy=0.817
train <- filter(all, Date.Production.Start < q[3]) # early 5% wells
test <- filter(all, Date.Production.Start > q[3]+365) # ~ 10% wells 

set.seed(777)
rf <- runRF2(train=train, test=test, model=formula.class2, m=5, no.tree=500)
sol <- rf[[1]]  # Pred accuracy on test data
rf.mod <- rf[[2]]  # rf model obj
plotRFOOBErr(rf.mod)

#-----
#20% Accy=0.811
train <- filter(all, Date.Production.Start < q[4]) # early 5% wells
test <- filter(all, Date.Production.Start > q[4]+365) # ~ 10% wells 

set.seed(777)
rf <- runRF2(train=train, test=test, model=formula.class2, m=5, no.tree=500)
sol <- rf[[1]]  # Pred accuracy on test data
rf.mod <- rf[[2]]  # rf model obj
plotRFOOBErr(rf.mod)

#-----
#25% Accy=0.816
train <- filter(all, Date.Production.Start < q[5]) # early 5% wells
test <- filter(all, Date.Production.Start > q[5]+365) # ~ 10% wells 

set.seed(777)
rf <- runRF2(train=train, test=test, model=formula.class2, m=5, no.tree=500)
sol <- rf[[1]]  # Pred accuracy on test data
rf.mod <- rf[[2]]  # rf model obj
plotRFOOBErr(rf.mod)

#-----
#30%
train <- filter(all, Date.Production.Start < q[6]) # early 5% wells
test <- filter(all, Date.Production.Start > q[6]+365) # ~ 10% wells 

set.seed(777)
rf <- runRF2(train=train, test=test, model=formula.class2, m=5, no.tree=500)
sol <- rf[[1]]  # Pred accuracy on test data
rf.mod <- rf[[2]]  # rf model obj
plotRFOOBErr(rf.mod)

#-----
#35%
train <- filter(all, Date.Production.Start < q[7]) # early 5% wells
test <- filter(all, Date.Production.Start > q[7]+365) # ~ 10% wells 

set.seed(777)
rf <- runRF2(train=train, test=test, model=formula.class2, m=5, no.tree=500)
sol <- rf[[1]]  # Pred accuracy on test data
rf.mod <- rf[[2]]  # rf model obj
plotRFOOBErr(rf.mod)

#-------------------
# 5%~25% results
accy<-data.frame(train.pct=seq(5,25,5), accy=c(0.77,0.80,0.82,0.81,0.82))
plotLine(accy, "Percentage of Training Data", "Test Classification Accuracy")


#-------------------------------------------------------------------------------------------------------------------------------
# Correlation study
library(ellipse)

X <- select(all, -Uwi, -Target, -Target.Q, -Target.Q4, -Latitude, -Longitude, -Date.Production.Start)
names(X) <- gsub("(Core.|.Kriged|.Joined)",'',names(X))  # strip unneccsary characters
plotCorr(X)

# top 10 overlap importance vars between accuracy and gini measurement
X.ov10 <- select(X, Core.RoCalculated.Kriged,Producer.DepthTrueVertical.Joined, Core.S2.Kriged, Core.Tmax.Kriged, Core.S3.Kriged)
names(X.ov10) <- gsub("(Core.|.Kriged|.Joined)",'',names(X.ov10))  # strip unneccsary characters
plotCorr(X.ov10)
corr <- as.matrix(cor(X.ov10))
# saveRDS(corr, "top10_overlap_corr.rds")
# corr <- readRDS("top10_overlap_corr.rds")

# h.corr <- which(corr>0.8&corr<1, arr.in=TRUE)
# 
# X.none.overlap <- select(all, match(accy.add, names(all)), match(gini.add, names(all)))
# corr.none.overlap <- as.matrix(cor(X.none.overlap))
# colnames(corr.none.overlap)=NULL
# rownames(corr.none.overlap)=NULL
# plotcorr(corr.none.overlap, numbers=T)



#================================================================================================================================
# Data overview
#================================================================================================================================

# Production well plot
dat <- select(all, Longitude, Latitude, Production=Target)
plotWellProd(dat)


# Production well + Core plot
core.loc <- core.loc %>% select(Longitude, Latitude) %>% mutate(ID="core")
prod.loc <- dat %>% select(Longitude, Latitude) %>% mutate(ID="prod")
core.prod <- rbind(core.loc, prod.loc)
plotCoreProd(core.prod)

