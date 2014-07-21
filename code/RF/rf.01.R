
library(randomForest)
library(MASS)

####################
## RF for regression
####################

data(Boston)
y <- Boston[,14]
x <- Boston[,1:13]

rf <- randomForest(x,y)
print(rf)

# Out-of-bag est
plot( predict(rf), y)
abline(c(0,1),col=2)

# Training data as new data
plot( predict(rf,newdata=x), y)
abline(c(0,1),col=2)

# mtry is the only par to be tuned
(rf.2v <- randomForest(x,y,mtry=2)) # 2 vars @ each split
(rf.4v <- randomForest(x,y,mtry=4)) # 4 vars @ each split
(rf.8v <- randomForest(x,y,mtry=8)) # 8 vars @ each split
(rf.10v <- randomForest(x,y,mtry=10)) # 10 vars @ each split


(rf.var.imp <- randomForest(x,y,importance=TRUE)) # assess importance of predictors
varImpPlot(rf.var.imp)


########################
## RF for classification
########################
data(fgl)
X <- fgl[,-10]
Y <- fgl[,10]

rf.cls <- randomForest(X, Y, ntree=500, proximity=TRUE, oob.prox=TRUE)

MDSplot(rf.cls, Y)




