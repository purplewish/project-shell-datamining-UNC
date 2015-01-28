
rm(list=ls())

setwd("C:/Users/Mingqi.Wu/Desktop/rf")

options( java.parameters = "-Xmx10000m" )

library(randomForest)
library(e1071)
library(bartMachine)

XDat<-read.csv("EagleFordOilInput.csv")
XDat<-XDat[!duplicated(XDat$Uwi),]
names(XDat)[names(XDat)=="Uwi"]<-"Uwi.x"
XModVars<-setdiff(names(XDat),c("Uwi.x","Latitude","Longitude","Producer.EstimatedLength.Joined"))

YDat<-read.csv("Rules features using recent and Jan 2014 data.csv")
YDat<-YDat[,c("Uwi","Rules.Prediction","Kaggle.Prediction","Target")]
names(YDat)[names(YDat)=="Uwi"]<-"Uwi.y"

AllDat<-cbind(XDat,YDat[match(as.character(XDat$Uwi.x),as.character(YDat$Uwi.y)),])
AllDat<-AllDat[!is.na(AllDat$Uwi.y),]
AllDat$Target.Quartile<-cut(AllDat$Target,breaks=quantile(AllDat$Target),labels=paste0("Q",1:4),include.lowest = T)
AllDat$Target.Quartile.Q4<- factor(AllDat$Target.Quartile=="Q4")

form.Cat<-formula(paste("Target.Quartile~",paste(XModVars,collapse="+")))

N=30
Percent.Train<-1
my.df<-data.frame(Index=1:N,PredQ4=NA)

for(i in 1:N){
print(i)
Train.Samp<-sort(sample(nrow(AllDat),Percent.Train*nrow(AllDat)))  # training data index
Test.Samp<-setdiff(1:nrow(AllDat),Train.Samp)  # test data index

(rf.Quartile<-randomForest(form.Cat,data=AllDat[Train.Samp,],importance=T))

#my.df$PredQ4[i]<-sum(AllDat[Test.Samp,]$Target.Quartile=="Q4" & predict( rf.Quartile,newdata=AllDat[Test.Samp,])=="Q4")/sum(AllDat[Test.Samp,]$Target.Quartile=="Q4")

my.df$PredQ4[i]<-sum(AllDat[1:nrow,]$Target.Quartile=="Q4" & predict( rf.Quartile,newdata=AllDat[Test.Samp,])=="Q4")/sum(AllDat[Test.Samp,]$Target.Quartile=="Q4")

}

summary(my.df$PredQ4)

hist(my.df$PredQ4)



form.Cat<-formula(paste("Target.Quartile~",paste(XModVars,collapse="+")))
form.Cont<-formula(paste("Target~",paste(XModVars,collapse="+")))

form.Cat<-formula(paste("Target.Quartile~",paste(XModVars,collapse="+")))
tr<-randomForest(form.Cat,data=AllDat,importance=T)
imp<-importance(tr)
impvar <- rownames(imp)[order(imp[, "MeanDecreaseGini"], decreasing=TRUE)]
tr

#impvar<-sample(impvar)
impvar<-rev(impvar)

out.df<-data.frame(Index=1:length(impvar),var=impvar,ErrorRate=NA)
for(i in 1:length(impvar)){
  print
  form.Cat1<-formula(paste("Target.Quartile~",paste(impvar[1:i],collapse="+")))
  tr1<-randomForest(form.Cat1,data=AllDat,importance=T,ntree=100)
  #out.df$ErrorRate[i]<-as.numeric(tr1$err.rate[nrow(tr1$err.rate),"OOB"])
  out.df$ErrorRate[i]<-sum(AllDat$Target.Quartile=="Q4" & predict(tr1)=="Q4")/sum(AllDat$Target.Quartile=="Q4")
  
  
}


plot(ErrorRate~Index,out.df)



