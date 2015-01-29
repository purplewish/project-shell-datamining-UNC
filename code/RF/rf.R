
rm(list=ls())
options(scipen=999)
options(java.parameters = "-Xmx10000m")

library(randomForest)
library(ggplot2)
library(reshape2)
library(dplyr)


#---------------------------------------------------------------------------------------------------------------------------------
## Prepare Dataset for RF model

# Load predictors Xs
setwd("C:/Apps/projects/DataMiningUNC/Kaggle/Final/Documentation and Input Files September 22 2014/Documentation and Input Files")

x <- read.csv("EagleFordOilInput.csv")
x <- distinct(x, Uwi)  # rm duplicate records (5222 x 35)
x <- select(x, -Latitude, -Longitude, -Producer.EstimatedLength.Joined)
x.vars <- names(x)
x.vars <- x.vars[-1]  # rm Uwi

# Load response Y
setwd("C:/Apps/projects/DataMiningUNC/Kaggle/Final/RulesBasedApproach Oct 8/RulesBasedApproach Oct 8")

y <- read.csv("Rules features using recent and Jan 2014 data.csv")
y <- select(y, Uwi, Target)  # 2631 x 2

# Merge X and Y
all <- left_join(x, y, by="Uwi")
all <- filter(all, !is.na(Target))  # rm missing production record

# Classify quantiles based on production
all$Target.Q <- cut(all$Target,breaks=quantile(all$Target),labels=paste0("Q",1:4),include.lowest=T)  # class: Q1 Q2 Q3 Q4
all <- mutate(all, Target.Q4=factor(Target.Q=="Q4"))  # class: Q4 ~Q4 TRUE FALSE

# Classification formula
formula.class1 <- formula(paste("Target.Q~", paste(x.vars,collapse="+")))  # class:Q1 Q2 Q3 Q4
formula.class2 <- formula(paste("Target.Q4~", paste(x.vars,collapse="+"))) # class:Q4 ~Q4

# Results directory
setwd("C:/Apps/projects/DataMiningUNC/Code/RF/results")


#-------------------------------------------------------------------------------------------------------------------------
## RF model: one set of pars 

# Pars
n <- 1  # reps for same set of pars
m <- 5  # mtry=sqrt(31)=5
no.tree <- 500
train.perc <- 0.75
pred.df <- data.frame(Index=1:n,Accy=NA,TP=NA,TN=NA)  # Accy: classification accuracy; TP/TN: True positive/True negative

set.seed(777)
sol <- NULL
for(i in 1:n){
  
  print(paste0("i=",i, " Train%=", train.perc), quote=F)
  
  # Split data into train/test set
  train <- sample_frac(all, train.perc, replace=F)
  test  <- dplyr::setdiff(all, train)
  
  #######################################################################################################
  rf.model <- randomForest(formula.class2, data=train, importance=T, mtry=m, do.trace=500, ntree=no.tree)  
  #######################################################################################################
  
  # Predict test set and calculate error rate
  test.pred <- cbind( test[, c(1,35)], Target.Q4.Pred=predict(rf.model, newdata=test) )  # Uwi, Target.Q4, Target.Q4.Pred
  pred.df$Accy[i] <- sum(test.pred[,2]==test.pred[,3])/nrow(test.pred)  # Classification accuracy
  pred.df$TP[i] <- sum(test.pred[,2]==test.pred[,3] & test.pred[,2]=="TRUE")/sum(test.pred[,2]=="TRUE")  # TP
  pred.df$TN[i] <- sum(test.pred[,2]==test.pred[,3] & test.pred[,2]=="FALSE")/sum(test.pred[,2]=="FALSE")  # TN
  
}

# pred.df
# Avg accuracy on test data at fixed training %
sol <- rbind(sol, c(m=m, n.Tree=no.tree, Train.Perc=train.perc, Accy=mean(pred.df[,2]), TP=mean(pred.df[,3]), TN=mean(pred.df[,4])) )
sol
# names(rf.model)


#-------------------------------------------------------------------------------------------------------------------------
## RF model: effect of number of trees (set no.tree=500)

# Pars
n <- 1  # reps for same set of pars
m <- 5  # mtry=sqrt(31)=5
no.tree <- 3000
train.perc <- 0.75
pred.df <- data.frame(Index=1:n,Accy=NA,TP=NA,TN=NA)  # Accy: classification accuracy; TP/Tn: True positive/True negative

set.seed(123)
sol <- NULL
for(i in 1:n){
  
  print(paste0("i=",i, " Train%=", train.perc), quote=F)
  
  # Split data into train/test set
  train <- sample_frac(all, train.perc, replace=F)
  test  <- dplyr::setdiff(all, train)
  
  #######################################################################################################
  rf.model <- randomForest(formula.class2, data=train, importance=T, mtry=m, do.trace=500, ntree=no.tree)  
  #######################################################################################################
  
  # Predict test set and calculate error rate
  test.pred <- cbind( test[, c(1,35)], Target.Q4.Pred=predict(rf.model, newdata=test) )  # Uwi, Target.Q4, Target.Q4.Pred
  pred.df$Accy[i] <- sum(test.pred[,2]==test.pred[,3])/nrow(test.pred)  # Classification accuracy
  pred.df$TP[i] <- sum(test.pred[,2]==test.pred[,3] & test.pred[,2]=="TRUE")/sum(test.pred[,2]=="TRUE")  # TP
  pred.df$TN[i] <- sum(test.pred[,2]==test.pred[,3] & test.pred[,2]=="FALSE")/sum(test.pred[,2]=="FALSE")  # TN
  
}

# pred.df
# Avg accuracy on test data at fixed training %
sol <- rbind(sol, c(m=m, n.Tree=no.tree, Train.Perc=train.perc, Accy=mean(pred.df[,2]), TP=mean(pred.df[,3]), TN=mean(pred.df[,4])) )  
sol

# OOB error rate at differnt number of trees
oob.err <- as.data.frame(cbind(No.Trees=1:no.tree, rf.model$err.rate))
names(oob.err)[2] <- "Classification error rate"
names(oob.err)[3] <- "~Q1 error rate"
names(oob.err)[4] <- "Q1 error rate"
dat <- melt(oob.err, id="No.Trees")

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
# RF model: effect of mtry (set m=5)

# Pars
n <- 1  # reps for same set of pars
m.seq <- c(3, 5, 6, 10)  # mtry=sqrt(31)=5
no.tree <- 500
train.perc <- 0.75
pred.df <- data.frame(Index=1:n,Accy=NA,TP=NA,TN=NA)  # Accy: classification accuracy; TP/TN: True positive/True negative

set.seed(456)
sol.all <- NULL
for(m in m.seq){
  
  sol <- NULL
  for(i in 1:n){
    
    print(paste0("i=",i, " Train%=", train.perc), quote=F)
    
    # Split data into train/test set
    train <- sample_frac(all, train.perc, replace=F)
    test  <- dplyr::setdiff(all, train)
    
    #######################################################################################################
    rf.model <- randomForest(formula.class2, data=train, importance=T, mtry=m, do.trace=500, ntree=no.tree)  
    #######################################################################################################
    
    # Predict test set and calculate error rate
    test.pred <- cbind( test[, c(1,35)], Target.Q4.Pred=predict(rf.model, newdata=test) )  # Uwi, Target.Q4, Target.Q4.Pred
    pred.df$Accy[i] <- sum(test.pred[,2]==test.pred[,3])/nrow(test.pred)  # Classification accuracy
    pred.df$TP[i] <- sum(test.pred[,2]==test.pred[,3] & test.pred[,2]=="TRUE")/sum(test.pred[,2]=="TRUE")  # TP
    pred.df$TN[i] <- sum(test.pred[,2]==test.pred[,3] & test.pred[,2]=="FALSE")/sum(test.pred[,2]=="FALSE")  # TN
    
  }
  sol <- rbind(sol, c(m=m, n.Tree=no.tree, Train.Perc=train.perc, Accy=mean(pred.df[,2]), TP=mean(pred.df[,3]), TN=mean(pred.df[,4])) )  # Avg accuracy on test data at fixed training %
  
  sol.all <- rbind(sol.all, sol)
}

sol.all
# write.csv(sol.all, "./mtry.csv", row.names=F)


#-------------------------------------------------------------------------------------------------------------------------------
# RF model: @different train % 

# Pars
n <- 50   # reps for same set of pars
m <- 5    # mtry=sqrt(31)=5
no.tree <- 500
train.perc.seq <- seq(0.1,0.9,0.1)
pred.df <- data.frame(Index=1:n,Train.Perc=NA,Accy=NA,TP=NA,TN=NA)  # Accy: classification accuracy; TP/Tn: True positive/True negative

set.seed(151)
sol.all <- NULL
for(train.perc in train.perc.seq){
  
  sol <- NULL
  for(i in 1:n){
    
    print(paste0("i=",i, " Train%=", train.perc), quote=F)
    
    # Split data into train/test set
    train <- sample_frac(all, train.perc, replace=F)
    test  <- dplyr::setdiff(all, train)
    
    #######################################################################################################
    rf.model <- randomForest(formula.class2, data=train, importance=T, mtry=m, do.trace=500, ntree=no.tree)  
    #######################################################################################################
    
    # Predict test set and calculate error rate
    test.pred <- cbind( test[, c(1,35)], Target.Q4.Pred=predict(rf.model, newdata=test) )  # Uwi, Target.Q4, Target.Q4.Pred
    
    pred.df$Accy[i] <- sum(test.pred[,2]==test.pred[,3])/nrow(test.pred)  # Classification accuracy
    pred.df$TP[i] <- sum(test.pred[,2]==test.pred[,3] & test.pred[,2]=="TRUE")/sum(test.pred[,2]=="TRUE")  # TP
    pred.df$TN[i] <- sum(test.pred[,2]==test.pred[,3] & test.pred[,2]=="FALSE")/sum(test.pred[,2]=="FALSE")  # TN
    pred.df$Train.Perc[i] <- train.perc
    
  }
  # write.table(pred.df, file="./train_perc_errrate.csv", row.names=F, append=T, sep=",")
  # Avg accuracy on test data at fixed training %
  sol <- rbind(sol, c(m=m, n.Tree=no.tree, Train.Perc=train.perc, Accy=mean(pred.df[,3]), TP=mean(pred.df[,4]), TN=mean(pred.df[,5])) )  
  sol.all <- rbind(sol.all, sol)
  
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

# Pars
m <- 5  # mtry=sqrt(31)=5
no.tree <- 500
train.perc <- 1

# Split data into train/test set
set.seed(777)
train <- sample_frac(all, train.perc, replace=F)
test  <- dplyr::setdiff(all, train)
  
#######################################################################################################
rf.model <- randomForest(formula.class2, data=train, importance=T, mtry=m, do.trace=500, ntree=no.tree)  
#######################################################################################################

# Importance data
d <- data.frame(rownames(importance(rf.model)),round(importance(rf.model),2))
# measure 1:mean decrease in accuracy  2:mean decrease in gini
names(d)[c(1,ncol(d)-1,ncol(d))] <- c("Predictor","mda","mdg")  
rownames(d) <- NULL

pred.acc  <- select(d, Predictor, mda)
pred.gini <- select(d, Predictor, mdg)

# Var importance plot function
importancePlot <- function(d,ylb,fontsize){
  fontsize <- as.numeric(fontsize)
  d <- d[order(d[,2],decreasing=T),]
  d$Predictor <- factor(as.character(d$Predictor),levels=rev(as.character(d$Predictor)))
  rownames(d) <- NULL
  abs.min <- abs(min(d[,2]))
  g1 <- ggplot(data=d,aes_string(x="Predictor",y=ylb,group="Predictor",colour="Predictor",fill="Predictor")) + geom_bar(stat="identity") + theme_grey(base_size=fontsize)
  #g1 <- ggplot(data=d,aes_string(x="Predictor",y=ylb,group="Predictor",color="black",fill="black")) + geom_bar(stat="identity") + theme_grey(base_size=fontsize)
  if(ylb=="mda") g1 <- g1 + labs(y="Mean decrease in accuracy") else if(ylb=="mdg") g1 <- g1 + labs(y="Mean decrease in Gini")
  g1 <- g1 + theme(axis.text.x = element_text(angle=90,hjust=1,vjust=0.4)) + geom_hline(yintercept=abs.min,linetype="dashed",colour="black") + coord_flip()
  print(g1)
}

importancePlot(d=pred.acc, ylb="mda", 20)
importancePlot(d=pred.gini, ylb="mdg", 20)

#pdf(file = file, width=11, height=8.5)
#dev.off()


# varImpPlot(rf.model, n.var=10)
# varImp1 <- sort(rf.model$importance[,5], dec=T)
# plot(varImp1, type="h")

# size of trees in an ensemble
# hist(treesize(rf.model))
