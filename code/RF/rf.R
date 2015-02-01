
setwd("C:/Apps/projects/GitHub/project-shell-datamining-UNC/code/RF")
source("header.R")
source("loadData.R")
source("runRF.R")
source("RFVarImp.R")

# Results directory
setwd("C:/Apps/projects/DataMiningUNC/Code/RF/results")

#-------------------------------------------------------------------------------------------------------------------------
## RF model on one set of pars 
set.seed(777)
rf <- runRF(dat=all, train.pct=0.75, model=formula.class2, m=5, no.tree=500, nrep=1)
sol <- rf[[1]]  # Pred accuracy on test data
rf.mod <- rf[[2]]  # Last rf model obj

#-------------------------------------------------------------------------------------------------------------------------
## RF model: effect of number of trees (select no.tree=500)
set.seed(123)
rf <- runRF(dat=all, train.pct=0.75, model=formula.class2, m=5, no.tree=3000, nrep=1)
sol <- rf[[1]]  # Pred accuracy on test data
rf.mod <- rf[[2]]  # Last rf model obj

# OOB error rate at differnt number of trees
oob.err <- as.data.frame(cbind(No.Trees=1:3000, rf.mod$err.rate))
names(oob.err)[2] <- "Classification error rate"
names(oob.err)[3] <- "~Q1 error rate"
names(oob.err)[4] <- "Q1 error rate"
dat <- melt(oob.err, id="No.Trees")

# Plot error rate vs. number of trees
ggplot(data=dat, aes(x=No.Trees, y=value, colour=variable)) + #+ xlim(0,1000) 
  geom_line(size=1.1) + geom_point(size=1) + 
  xlab("Number of trees") + ylab("Error rate") +
  theme(
    axis.title.x = element_text(size=24),
    axis.title.y = element_text(size=24),
    axis.text.x = element_text(colour="grey20",size=15),
    axis.text.y = element_text(colour="grey20",size=15),
    legend.title=element_blank(),
    legend.text = element_text(size = 20),
    legend.justification=c(1,0), legend.position=c(1,0.9),
    legend.background = element_rect(fill="gray90", size=.5, linetype="dotted")
  )

#----------------------------------------------------------------------------------------------------------------------------
# RF model: effect of mtry (select m=5)
set.seed(789)
m.seq <- c(3, 5, 6, 10)  # mtry=sqrt(31)=5
sol.all <- NULL
for(m in m.seq){
  rf <- runRF(dat=all, train.pct=0.75, model=formula.class2, m=m, no.tree=500, nrep=1)
  sol.all <- rbind(sol.all, rf[[1]])
}
sol.all
# write.csv(sol.all, "./mtry.csv", row.names=F)

#-------------------------------------------------------------------------------------------------------------------------------
# RF model: @different train % 
set.seed(151)
sol.all <- NULL
train.pct.seq <- seq(0.1,0.9,0.1)
for(train.pct in train.pct.seq){
  rf <- runRF(dat=all, train.pct=train.pct, model=formula.class2, m=5, no.tree=500, nrep=50)
  sol.all <- rbind(sol.all, rf[[1]])
}
sol.all
# write.csv(sol.all, "./train_perc_errrate_avg.csv", row.names=F)

# Plot accuracy at different training % (Avg 50 runs results)
err <- as.data.frame(sol.all[,seq(3,6)])
names(err)[2] <- "Classification Accuracy"
names(err)[3] <- "True Positive Rate"
names(err)[4] <- "True Negative Rate"
dat <- melt(err, id="Train.Perc")
dat.accy <- dat[1:9,]

ggplot(data=dat, aes(x=Train.Perc, y=value, colour=variable)) + 
  geom_line(size=1.1) + geom_point(size=4) + 
  ylim(0.5,1) + xlim(0.1,0.9) + scale_x_continuous(breaks=seq(0.1,0.9,0.1)) +
  xlab("Percentage of Training Data") + ylab("Test Classification Accuracy") + 
  theme(
    axis.title.x = element_text(size=24),
    axis.title.y = element_text(size=24),
    axis.text.x = element_text(colour="grey20",size=15),
    axis.text.y = element_text(colour="grey20",size=15),
    legend.title=element_blank(),
    legend.text = element_text(size = 20),
    legend.justification=c(1,0), legend.position=c(1,0),
    legend.background = element_rect(fill="gray90", size=.5, linetype="dotted")
  )


# Plot accuracy at training % = 70%  
dat.train.70pct <- read.csv("./train_perc_0.7_errrate.csv", header=T)
dat2 <- dat.train.70pct[,-2]
names(dat2)[2] <- "Classification Accuracy"
names(dat2)[3] <- "True Positive Rate"
names(dat2)[4] <- "Ture Negative Rate"
dat2 <- melt(dat2, id="Index")

ggplot(data=dat2, aes(x=Index, y=value, colour=variable)) + 
  geom_line(size=1.1) + geom_point(size=4) + 
  ylim(0.3,1) + xlim(1,30) + scale_x_continuous(breaks=seq(1,30,1)) +
  xlab("Runs for 70% Training Dataset") + ylab("Test Classification Accuracy") + 
  theme(
    axis.title.x = element_text(size=24),
    axis.title.y = element_text(size=24),
    axis.text.x = element_text(colour="grey20",size=15),
    axis.text.y = element_text(colour="grey20",size=15),
    legend.title=element_blank(),
    legend.text = element_text(size = 20),
    legend.justification=c(1,0), legend.position=c(1,0),
    legend.background = element_rect(fill="gray90", size=.5, linetype="dotted")
  )


#-------------------------------------------------------------------------------------------------------------------------------------
# RF model: variables importance 

set.seed(777)
rf <- runRF(dat=all, train.pct=0.75, model=formula.class2, m=5, no.tree=500, nrep=1)
rf.mod <- rf[[2]]

plotRFVarImp(rf.mod)


#pdf(file = file, width=11, height=8.5)
#dev.off()

# varImpPlot(rf.model, n.var=10)
# varImp1 <- sort(rf.model$importance[,5], dec=T)
# plot(varImp1, type="h")

# size of trees in an ensemble
# hist(treesize(rf.model))


