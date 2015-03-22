
library(dplyr)
library(randomForest)

###############################################################################
runRF <- function(dat, train.pct, model, m, no.tree, nrep=1, ntrace=500){
# Run RF model n times with given pars based on full dataset
#
# Args:
#   dat: data set for rf, will be split to train and test set
#   train.pct: training data percentage 
#   model: formula for RF model
#   m: Number of variables randomly sampled as candidates at each split; mtry
#   no.tree: Number of trees
#   nrep: Number of replication to run RF model under the same set of pars
#   ntrace: running output is printed for every ntrace trees.
#
# Returns:
#	1. Prediction accuracy on test data; 2. Last RF model object
	sol <- NULL
	# Prediction data frame (Accy: classification accuracy; TP/TN: True+/True-)
	pred.df <- data.frame(Index=1:nrep,Accy=NA,TP=NA,TN=NA) 

	for(i in 1:nrep){
		print(paste0("i=",i, " Train%=", train.pct), quote=F)

		# Split data into train/test set
		train <- sample_frac(dat, train.pct, replace=F)
		test  <- dplyr::setdiff(dat, train)

		#################################################################################################
		rf.model <- randomForest(model, data=train, importance=T, mtry=m, do.trace=ntrace, ntree=no.tree)  
		#################################################################################################

		# Predict test dataset and calculate error rate
		test.pred <- cbind(test[,c(1,35)],Target.Q4.Pred=predict(rf.model,newdata=test))  # Uwi, Target.Q4, Target.Q4.Pred
		pred.df$Accy[i] <- sum(test.pred[,2]==test.pred[,3])/nrow(test.pred)  # Classification accuracy
		pred.df$TP[i] <- sum(test.pred[,2]==test.pred[,3] & test.pred[,2]=="TRUE")/sum(test.pred[,2]=="TRUE")  # TP
		pred.df$TN[i] <- sum(test.pred[,2]==test.pred[,3] & test.pred[,2]=="FALSE")/sum(test.pred[,2]=="FALSE")  # TN
	}

	# Avg accuracy over nrep on test dataset at fixed training %
	sol <- rbind(sol, c(m=m, n.Tree=no.tree, Train.Perc=train.pct, Accy=mean(pred.df[,2]), TP=mean(pred.df[,3]), TN=mean(pred.df[,4])) )
	
	return(list(sol, rf.model)) 
}


#######################################################################################################################################
runRF2 <- function(train, test, model, m, no.tree, ntrace=500){
# Run RF model based on training and testing datesets
#
# Args:
#   train: training data 
#   test: testing data 
#   model: formula for RF model
#   m: Number of variables randomly sampled as candidates at each split; mtry
#   no.tree: Number of trees
#   ntrace: running output is printed for every ntrace trees.
#
# Returns:
# 1. Prediction accuracy on test data; 2. RF model object
  
#################################################################################################
rf.model <- randomForest(model, data=train, importance=T, mtry=m, do.trace=ntrace, ntree=no.tree)  
#################################################################################################
    
# Predict test dataset and calculate error rate (Accy: classification accuracy; TP/TN: True+/True-)
test.pred <- cbind(test[,c(1,35)],Target.Q4.Pred=predict(rf.model,newdata=test))  # Uwi, Target.Q4, Target.Q4.Pred
pred.Accy <- sum(test.pred[,2]==test.pred[,3])/nrow(test.pred)  # Classification accuracy
pred.TP   <- sum(test.pred[,2]==test.pred[,3] & test.pred[,2]=="TRUE")/sum(test.pred[,2]=="TRUE")  # TP
pred.TN   <- sum(test.pred[,2]==test.pred[,3] & test.pred[,2]=="FALSE")/sum(test.pred[,2]=="FALSE")  # TN
sol <- c(m=m, n.Tree=no.tree, Accy=pred.Accy, TP=pred.TP, TN=pred.TN)

return(list(sol, rf.model)) 
}


####################################################################################################
runRFCV <- function(dat, model, m, no.tree, k, rev=FALSE, ntrace=500){
  # Run RF model with Cross Validation
  #
  # Args:
  #   dat: data set for rf, will be split to train and test set
  #   model: formula for RF model
  #   m: Number of variables randomly sampled as candidates at each split; mtry
  #   no.tree: Number of trees
  #   k:   K-fold CV 
  #   rev: reverse CV or not
  #   ntrace: running output is printed for every ntrace trees.
  #
  # Returns:
  #   1. Prediction accuracy based on CV  2. Prediction results for each fold
  
  folds <- cvFolds(nrow(dat), K=k)
  accy <-0; tp <-0; tn <-0;
  pred <- NULL; sol <- NULL;
  
  for(i in 1:k){  
    # Split data into train/test set
    if(rev==TRUE){
      train <- dat[folds$subsets[folds$which==i],]
      test  <- dplyr::setdiff(dat, train)
    } else {
      test  <- dat[folds$subsets[folds$which==i],]
      train <- dplyr::setdiff(dat, test)
    }
    
    #################################################################################################
    rf.model <- randomForest(model, data=train, importance=T, mtry=m, do.trace=ntrace, ntree=no.tree)  
    #################################################################################################
    
    # Predict test dataset and calculate error rate
    test.pred <- cbind(test[,c(1,35:37)],Target.Q4.Pred=predict(rf.model,newdata=test))  # Uwi, Target.Q4, Latitude, Longitude, Target.Q4.Pred
    accy <- accy + sum(test.pred[,2]==test.pred[,5])/nrow(test.pred)  # Accuracy
    tp   <- tp   + sum(test.pred[,2]==test.pred[,5] & test.pred[,2]=="TRUE")/sum(test.pred[,2]=="TRUE")  # TP
    tn   <- tn   + sum(test.pred[,2]==test.pred[,5] & test.pred[,2]=="FALSE")/sum(test.pred[,2]=="FALSE")  # TN  
  
    pred <- rbind(pred, test.pred)  # save prediction results for fold i
    
    print(paste0("K=", k, " i=", i, " Rev=", rev), quote=F)
  }
  # CV results
  sol <- data.frame(K=k, Rev=rev, Accy=accy/k, TP=tp/k, TN=tn/k, m=m, n.Tree=no.tree) 
  
  return(list(sol, pred))
  
}


####################################################################################################
runRFCV2 <- function(dat, model, no.tree, k, rev=FALSE, ntrace=500){
  # Run RF model with Cross Validation using default setting mtry
  #
  # Args:
  #   dat: data set for rf, will be split to train and test set
  #   model: formula for RF model
  #   no.tree: Number of trees
  #   k:   K-fold CV 
  #   rev: reverse CV or not
  #   ntrace: running output is printed for every ntrace trees.
  #
  # Returns:
  #   1. Prediction accuracy based on CV  2. Prediction results for each fold
  
  folds <- cvFolds(nrow(dat), K=k)
  accy <-0; tp <-0; tn <-0;
  pred <- NULL; sol <- NULL;
  
  for(i in 1:k){  
    # Split data into train/test set
    if(rev==TRUE){
      train <- dat[folds$subsets[folds$which==i],]
      test  <- dplyr::setdiff(dat, train)
    } else {
      test  <- dat[folds$subsets[folds$which==i],]
      train <- dplyr::setdiff(dat, test)
    }
    
    #################################################################################################
    rf.model <- randomForest(model, data=train, importance=T, do.trace=ntrace, ntree=no.tree)  
    #################################################################################################
    
    # Predict test dataset and calculate error rate
    test.pred <- cbind(test[,c(1,35:37)],Target.Q4.Pred=predict(rf.model,newdata=test))  # Uwi, Target.Q4, Latitude, Longitude, Target.Q4.Pred
    accy <- accy + sum(test.pred[,2]==test.pred[,5])/nrow(test.pred)  # Accuracy
    tp   <- tp   + sum(test.pred[,2]==test.pred[,5] & test.pred[,2]=="TRUE")/sum(test.pred[,2]=="TRUE")  # TP
    tn   <- tn   + sum(test.pred[,2]==test.pred[,5] & test.pred[,2]=="FALSE")/sum(test.pred[,2]=="FALSE")  # TN  
    
    pred <- rbind(pred, test.pred)  # save prediction results for fold i
    
    print(paste0("K=", k, " i=", i, " Rev=", rev), quote=F)
  }
  # CV results
  sol <- data.frame(K=k, Rev=rev, Accy=accy/k, TP=tp/k, TN=tn/k, n.Tree=no.tree) 
  
  return(list(sol, pred))
  
}
