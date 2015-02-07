
setwd("C:/Apps/projects/GitHub/project-shell-datamining-UNC/code/RF")
source("header.R")
source("loadData.R")
source("runRF.R")
source("plotFuns.R")

# Results directory
setwd("C:/Apps/projects/DataMiningUNC/Code/RF/results")

#-------------------------------------------------------------------------------------------------------------------------
## RF model on one set of pars 
set.seed(777)
rf <- runRF(dat=all, train.pct=0.25, model=formula.class2, m=5, no.tree=500, nrep=1)
sol <- rf[[1]]  # Pred accuracy on test data
rf.mod <- rf[[2]]  # Last rf model obj


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
# Prediction based on time cutoff

# Cutoff=q1
date.q <- summary(all$Date.Production.Start) # check quantiles of production starting date
q1 <- filter(all, Date.Production.Start < date.q[2])
q1.test <- filter(all, Date.Production.Start > date.q[2]+365)

set.seed(777)
rf <- runRF2(train=q1, test=q1.test, model=formula.class2, m=5, no.tree=500)
sol <- rf[[1]]  # Pred accuracy on test data
rf.mod <- rf[[2]]  # rf model obj

plotRFOOBErr(rf.mod)



