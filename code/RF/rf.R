
# code path
setwd("Z:/Mingqi.Wu/GitHup/project-shell-datamining-UNC/code/RF")
source("header.R")
source("loadData.R")
source("runRF.R")
source("plotFuns.R")

# Results directory
setwd(file.path(repo_path, "Code/RF/results"))

#-------------------------------------------------------------------------------------------------------------------------
## RF model on one set of pars 
set.seed(777)
rf <- runRF(dat=all, train.pct=1, model=formula.class2, m=5, no.tree=500, nrep=1)
sol <- rf[[1]]  # Pred accuracy on test data
rf.mod <- rf[[2]]  # Last rf model obj

## Plot sweet-spots
target <- select(all, Uwi, Target.Q4, Latitude, Longitude)
plotSweetspot(rf.mod, target)


#-------------------------------------------------------------------------------------------------------------------------
## RF model on selected important vars
nrep=50;
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
#sols.top10vars <- readRDS("top10_imp_vars.rds")

# saveRDS(sols, "top5_imp_vars.rds")
# sols.top5vars <- readRDS("top5_imp_vars.rds")
# 
# saveRDS(sols, "top3_imp_vars.rds")
# sols.top3vars <- readRDS("top3_imp_vars.rds")

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
# top10.overlap <- readRDS("top10_overlap.rds")  # 5 overlaped
# accy.add <- setdiff(sel.vars.accy, top10.overlap) # overlap Xs rank: 1,2,3,4,8
# gini.add <- setdiff(sel.vars.gini, top10.overlap) # overlap Xs rank: 1,2,6,7,10

## overlapped vars of top 10 important vars(gini and accy); totally 5 vars are overlapped
top10.overlap <- readRDS("top10_overlap.rds")  # 5 overlaped vars
formula.class.top10overlap <- formula(paste("Target.Q4~", paste(top10.overlap,collapse="+"))) # class:Q4 ~Q4, topQ vs. ~topQ

nrep=50;
train.pct.seq <- seq(0.2,0.8,0.1);

sols <- NULL
set.seed(777)
for (train.pct in train.pct.seq){
  
  rf <- runRF(dat=all, train.pct=train.pct, model=formula.class2, m=5, no.tree=500, nrep=nrep)
  sol <- data.frame(rf[[1]], method="all")  # Pred accuracy on test data
      
  rf <- runRF(dat=all, train.pct=train.pct, model=formula.class.top10overlap, m=5, no.tree=500, nrep=nrep)
  sol.overlap <- data.frame(rf[[1]], method="top10overlap")  # Pred accuracy on test data
    
  sols <- rbind(sols, sol, sol.overlap)
}

saveRDS(sols, "top10_overlap_vars.rds")
sols.top10overlap <- readRDS("top10_overlap_vars.rds")


#-------------------------------------------------------------------------------------------------------------------------
## RF model: effect of number of trees (select no.tree=500)
num.tree <- 3000

set.seed(123)
rf <- runRF(dat=all, train.pct=0.75, model=formula.class2, m=5, no.tree=num.tree, nrep=1)
sol <- rf[[1]]  # Pred accuracy on test data
rf.mod <- rf[[2]]  # Last rf model obj

plotRFOOBErr(rf.mod)


#----------------------------------------------------------------------------------------------------------------------------
# RF model: effect of mtry (select m=5)
m.seq <- c(3, 5, 6, 10)  # mtry=sqrt(31)=5

set.seed(789)
sol.all <- NULL
for(m in m.seq){
  rf <- runRF(dat=all, train.pct=0.75, model=formula.class2, m=m, no.tree=500, nrep=1)
  sol.all <- rbind(sol.all, rf[[1]])
}
sol.all
# write.csv(sol.all, "./mtry.csv", row.names=F)


#-------------------------------------------------------------------------------------------------------------------------------
# RF model: effect of different train % 
train.pct.seq <- seq(0.1,0.9,0.1)

set.seed(151)
sol.all <- NULL
for(train.pct in train.pct.seq){
  rf <- runRF(dat=all, train.pct=train.pct, model=formula.class2, m=5, no.tree=500, nrep=5)
  sol.all <- rbind(sol.all, rf[[1]])
}
sol.all
# write.csv(sol.all, "./train_perc_errrate_avg.csv", row.names=F)

plotRFAcc(sol.all)


#-------------------------------------------------------------------------------------------------------------------------------------
# RF model: variables importance 

set.seed(777)
rf <- runRF(dat=all, train.pct=0.75, model=formula.class2, m=5, no.tree=500, nrep=1)
rf.mod <- rf[[2]]

plotRFVarImp(rf.mod)
plotRFVarImp2(rf.mod)


#-------------------------------------------------------------------------------------------------------------------------------------
# Prediction based on timeline (cutoff at time point)

# Cutoff=q1, test well > cutoff+365
q <- summary(all$Date.Production.Start) # check quantiles of production starting date
train <- filter(all, Date.Production.Start < q[2]) # early 25% wells
test <- filter(all, Date.Production.Start > q[2]+365) # ~ 10% wells 

set.seed(777)
rf <- runRF2(train=train, test=test, model=formula.class2, m=5, no.tree=500)
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
corr <- as.matrix(cor(X))
h.corr <- which(corr>0.8&corr<1, arr.in=TRUE)

X.none.overlap <- select(all, match(accy.add, names(all)), match(gini.add, names(all)))
corr.none.overlap <- as.matrix(cor(X.none.overlap))
colnames(corr.none.overlap)=NULL
rownames(corr.none.overlap)=NULL
plotcorr(corr.none.overlap, numbers=T)
