
rm(list=ls())
options(scipen=999)
options(java.parameters = "-Xmx10000m")
  
library(randomForest)
library(sqldf)
library(ggplot2)
library(reshape2)

setwd("C:/Users/Mingqi.Wu/Desktop/rf")


# Load predictors
x <- read.csv("EagleFordOilInput.csv")
x <- x[!duplicated(x$Uwi),]  # 5222 x 35
x <- subset(x, select=-c(Latitude, Longitude, Producer.EstimatedLength.Joined))
x.vars <- names(x)
x.vars <- gsub('.','_', x.vars, fixed=TRUE)
x.vars <- x.vars[-1]  # rm Uwi


# Load response
y <- read.csv("Rules features using recent and Jan 2014 data.csv")
y <- y[,c("Uwi", "Target")]  # 2631 x 2


# Merge X and Y
all <- sqldf("Select * from x a Left Join y b on (a.Uwi=b.Uwi)")
all <- all[!is.na(all$Target),]  # rm missing production record
all <- all[,-33]  # rm duplicate Uwi


# Classify quantiles based on production
all$Target.Q <- cut(all$Target,breaks=quantile(all$Target),labels=paste0("Q",1:4),include.lowest=T)  # class: Q1 Q2 Q3 Q4
all$Target.Q4 <- factor(all$Target.Q=="Q4")  # class: Q4 ~Q4 TRUE FALSE


# Classification formula
#formula.class1 <- formula(paste("Target.Q~",paste(x.vars,collapse="+")))  # class:Q1 Q2 Q3 Q4
formula.class2 <- formula(paste("Target.Q4~",paste(x.vars,collapse="+")))  # class:Q4 ~Q4


#--------------------------------------------------------------------------------------------------------
# RF model: one set of pars 

# Pars
n <- 3  # reps for same set of pars
m <- 5  # mtry=sqrt(31)=5
no.tree <- 500
train.perc <- 0.75
pred.df <- data.frame(Index=1:n,Accy=NA,TP=NA,TN=NA)  # Accy: classification accuracy; TP/Tn: True positive/True negative

set.seed(777)
sol <- NULL
for(i in 1:n){
  
  print(paste0("i=",i, " Train%=", train.perc), quote=F)
      
  train.samp <- sort(sample(nrow(all),train.perc*nrow(all)))  # training data index
  test.samp <- setdiff(1:nrow(all),train.samp)  # test data index
  
  rf.model <- randomForest(formula.class2,data=all[train.samp,],importance=T, mtry=m, do.trace=500, ntree=no.tree)  
    
  test.pred <- cbind( all[test.samp, c(1,35)], Target.Q4.Pred=predict(rf.model,newdata=all[test.samp,]) )  # Uwi, Target.Q4, Target.Q4.Pred
  pred.df$Accy[i] <- sum(test.pred[,2]==test.pred[,3])/nrow(test.pred)  # Classification accuracy
  pred.df$TP[i] <- sum(test.pred[,2]==test.pred[,3] & test.pred[,2]=="TRUE")/sum(test.pred[,2]=="TRUE")  # TP
  pred.df$TN[i] <- sum(test.pred[,2]==test.pred[,3] & test.pred[,2]=="FALSE")/sum(test.pred[,2]=="FALSE")  # TN
  
}

pred.df
sol <- rbind(sol, c(m=m, n.Tree=no.tree, Train.Perc=train.perc, Accy=mean(pred.df[,2]), TP=mean(pred.df[,3]), TN=mean(pred.df[,4])) )  # Avg accuracy on test data at fixed training %
sol
# names(rf.model)

#--------------------------------------------------------------------------------------------------------
# RF model: effect of number of trees (set no.tree=1000)
library("reshape2")

# Pars
n <- 2  # reps for same set of pars
m <- 5  # mtry=sqrt(31)=5
no.tree <- 3000
train.perc <- 0.75
pred.df <- data.frame(Index=1:n,Accy=NA,TP=NA,TN=NA)  # Accy: classification accuracy; TP/Tn: True positive/True negative

set.seed(123)
sol <- NULL
for(i in 1:n){
  
  print(paste0("i=",i, " Train%=", train.perc), quote=F)
  
  train.samp <- sort(sample(nrow(all),train.perc*nrow(all)))  # training data index
  test.samp <- setdiff(1:nrow(all),train.samp)  # test data index
  
  rf.model <- randomForest(formula.class2,data=all[train.samp,],importance=T, mtry=m, do.trace=500, ntree=no.tree)  
  
  test.pred <- cbind( all[test.samp, c(1,35)], Target.Q4.Pred=predict(rf.model,newdata=all[test.samp,]) )  # Uwi, Target.Q4, Target.Q4.Pred
  pred.df$Accy[i] <- sum(test.pred[,2]==test.pred[,3])/nrow(test.pred)  # Classification accuracy
  pred.df$TP[i] <- sum(test.pred[,2]==test.pred[,3] & test.pred[,2]=="TRUE")/sum(test.pred[,2]=="TRUE")  # TP
  pred.df$TN[i] <- sum(test.pred[,2]==test.pred[,3] & test.pred[,2]=="FALSE")/sum(test.pred[,2]=="FALSE")  # TN
  
}

pred.df
sol <- rbind(sol, c(m=m, n.Tree=no.tree, Train.Perc=train.perc, Accy=mean(pred.df[,2]), TP=mean(pred.df[,3]), TN=mean(pred.df[,4])) )  # Avg accuracy on test data at fixed training %
sol

oob.err <- as.data.frame(cbind(No.Trees=1:no.tree, rf.model$err.rate))
names(oob.err)[2] <- "Classification error rate"
names(oob.err)[3] <- "NotQ1 error rate"
names(oob.err)[4] <- "Q1 error rate"
dat = melt(oob.err, id="No.Trees")

ggplot(data=dat, aes(x=No.Trees, y=value, colour=variable))+ geom_line() #+ xlim(0,1000)


#--------------------------------------------------------------------------------------------------------
# RF model: effect of mtry (set m=10)
library("reshape2")

# Pars
n <- 3  # reps for same set of pars
m.seq <- c(3, 5, 10)  # mtry=sqrt(31)=5
no.tree <- 1000
train.perc <- 0.75
pred.df <- data.frame(Index=1:n,Accy=NA,TP=NA,TN=NA)  # Accy: classification accuracy; TP/Tn: True positive/True negative

set.seed(456)
sol.all <- NULL
for(m in m.seq){

  sol <- NULL
  for(i in 1:n){
    
    print(paste0("i=",i, " Train%=", train.perc), quote=F)
    
    train.samp <- sort(sample(nrow(all),train.perc*nrow(all)))  # training data index
    test.samp <- setdiff(1:nrow(all),train.samp)  # test data index
    
    rf.model <- randomForest(formula.class2,data=all[train.samp,],importance=T, mtry=m, do.trace=500, ntree=no.tree)  
    
    test.pred <- cbind( all[test.samp, c(1,35)], Target.Q4.Pred=predict(rf.model,newdata=all[test.samp,]) )  # Uwi, Target.Q4, Target.Q4.Pred
    pred.df$Accy[i] <- sum(test.pred[,2]==test.pred[,3])/nrow(test.pred)  # Classification accuracy
    pred.df$TP[i] <- sum(test.pred[,2]==test.pred[,3] & test.pred[,2]=="TRUE")/sum(test.pred[,2]=="TRUE")  # TP
    pred.df$TN[i] <- sum(test.pred[,2]==test.pred[,3] & test.pred[,2]=="FALSE")/sum(test.pred[,2]=="FALSE")  # TN
    
  }
  sol <- rbind(sol, c(m=m, n.Tree=no.tree, Train.Perc=train.perc, Accy=mean(pred.df[,2]), TP=mean(pred.df[,3]), TN=mean(pred.df[,4])) )  # Avg accuracy on test data at fixed training %

  sol.all <- rbind(sol.all, sol)
}

sol.all
#write.csv(sol.all, "./mtry.csv", row.names=F)


#--------------------------------------------------------------------------------------------------------
# RF model: @different train % 

# Pars
n <- 30   # reps for same set of pars
m <- 10  # mtry=sqrt(31)=5
no.tree <- 1000
train.perc.seq <- seq(0.1,0.9,0.1)
pred.df <- data.frame(Index=1:n,Train.Perc=NA,Accy=NA,TP=NA,TN=NA)  # Accy: classification accuracy; TP/Tn: True positive/True negative

set.seed(151)
sol.all <- NULL
for(train.perc in train.perc.seq){
  
  sol <- NULL
  for(i in 1:n){
    
    print(paste0("i=",i, " Train%=", train.perc), quote=F)
    
    train.samp <- sort(sample(nrow(all),train.perc*nrow(all)))  # training data index
    test.samp <- setdiff(1:nrow(all),train.samp)  # test data index
    
    rf.model <- randomForest(formula.class2,data=all[train.samp,],importance=T, mtry=m, do.trace=1000, ntree=no.tree)  
    
    test.pred <- cbind( all[test.samp, c(1,35)], Target.Q4.Pred=predict(rf.model,newdata=all[test.samp,]) )  # Uwi, Target.Q4, Target.Q4.Pred
    pred.df$Accy[i] <- sum(test.pred[,2]==test.pred[,3])/nrow(test.pred)  # Classification accuracy
    pred.df$TP[i] <- sum(test.pred[,2]==test.pred[,3] & test.pred[,2]=="TRUE")/sum(test.pred[,2]=="TRUE")  # TP
    pred.df$TN[i] <- sum(test.pred[,2]==test.pred[,3] & test.pred[,2]=="FALSE")/sum(test.pred[,2]=="FALSE")  # TN
    pred.df$Train.Perc[i] <- train.perc
    
  }
 # write.table(pred.df, file="./train_perc_errrate.csv", row.names=F, append=T, sep=",")
  sol <- rbind(sol, c(m=m, n.Tree=no.tree, Train.Perc=train.perc, Accy=mean(pred.df[,3]), TP=mean(pred.df[,4]), TN=mean(pred.df[,5])) )  # Avg accuracy on test data at fixed training %

  sol.all <- rbind(sol.all, sol)

}

sol.all
# write.csv(sol.all, "./train_perc_errrate_avg.csv", row.names=F)

# Plot accuracy at different training % Avg results
err <- as.data.frame(sol.all[,seq(3,6)])
names(err)[2] <- "Classification Accuracy"
names(err)[3] <- "True Positive Rate"
names(err)[4] <- "Ture Negative Rate"
dat <- melt(err, id="Train.Perc")
dat.accy <- dat[1:9,]

ggplot(data=dat, aes(x=Train.Perc, y=value, colour=variable)) + 
  geom_line(size=1.1) + geom_point(size=4) + 
  ylim(0.5,1) + xlim(0.1,0.9) + scale_x_continuous(breaks=seq(0.1,0.9,0.1)) +
  xlab("Percentage of Training Data") + ylab("Test Classification Accuracy") + theme(
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
  xlab("Runs for 70% Training Dataset") + ylab("Test Classification Accuracy") + theme(
    axis.title.x = element_text(size=24),
    axis.title.y = element_text(size=24),
    axis.text.x = element_text(colour="grey20",size=15),
    axis.text.y = element_text(colour="grey20",size=15),
    legend.title=element_blank(),
    legend.text = element_text(size = 20),
    legend.justification=c(1,0), legend.position=c(1,0),
    legend.background = element_rect(fill="gray90", size=.5, linetype="dotted")
  )













n=10
sol <- nULL
#Perc.TrainSeq <- seq(0.1,0.95,0.05)
Perc.TrainSeq <- c(0.1, 0.75)
set.seed(777)
for (Perc.Train in Perc.TrainSeq){
  
  pred.df <- data.frame(Index=1:n,PredQ4=nA)
  
  for(i in 1:n){

    print(i)
    train.samp <- sort(sample(nrow(all),Perc.Train*nrow(all)))  # training data index
    test.samp <- setdiff(1:nrow(all),train.samp)  # test data index
    
    rf.Q <- randomForest(form.Class,data=all[train.samp,],importance=T, mtry=10, do.trace=100)  # mtry=sqrt(31)
    # print(rf.Q)  # trees=500, # no. of vars tried at each split=5
    
    pred.df$PredQ4[i] <- sum(all[test.samp,]$Target.Q=="Q4" & predict(rf.Q,newdata=all[test.samp,])=="Q4")/sum(all[test.samp,]$Target.Q=="Q4")  # accuracy on test data
  
  }
  
  sol <- rbind(sol, c(TrainPerc=Perc.Train, Accuracy=mean(pred.df[,2])) )  # Avg accuracy on test data at fixed training %

}

sol
# write.csv(sol, "./trainPerc_Q4Accuracy.csv", row.names=F)

# size of trees in an ensemble
hist(treesize(rf.Q))

varImpPlot(rf.Q)

varImp1 <- sort(rf.Q$importance[,5], dec=T)
plot(varImp1, type="h")


