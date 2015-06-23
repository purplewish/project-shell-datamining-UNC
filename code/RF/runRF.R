
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



###############################################################################
runRFReg <- function(dat, train.pct, model, m, no.tree, ntrace=500){
# Run RF model on regression with given pars 
#
# Args:
#   dat: data set for rf, will be split to train and test set
#   train.pct: training data percentage 
#   model: formula for RF model
#   m: Number of variables randomly sampled as candidates at each split; mtry
#   no.tree: Number of trees
#   ntrace: running output is printed for every ntrace trees.
#
# Returns:
#	1. Prediction mse on test data; 2. Last RF model object
	sol <- NULL
	
	print(paste0("Train%=", train.pct), quote=F)

	# Split data into train/test set
	train <- sample_frac(dat, train.pct, replace=F)
	test  <- dplyr::setdiff(dat, train)

	#################################################################################################
	rf.model <- randomForest(model, data=train, importance=T, mtry=m, do.trace=ntrace, ntree=no.tree)  
	#################################################################################################

	# Predict test dataset and calculate mse
	test.pred <- cbind(test[,c(1,33)],Pred=predict(rf.model,newdata=test))  # Uwi, Target, Pred
	mse <- sum((test.pred[,2]-test.pred[,3])^2)/nrow(test.pred)  # mean square of residuals

	# Avg accuracy over nrep on test dataset at fixed training %
	sol <- rbind(sol, c(m=m, n.Tree=no.tree, Train.Perc=train.pct, mse=mse))
	
	return(list(sol, rf.model)) 
}




###############################################################################
runRFReg2 <- function(dat, cutoff, model, m, no.tree, ntrace=500){
  # Run RF model on regression with given pars 
  #
  # Args:
  #   dat: data set for rf, will be split to train and test set
  #   cutoff: cutoff date, wells before cutoff date will be used for training
  #   model: formula for RF model
  #   m: Number of variables randomly sampled as candidates at each split; mtry
  #   no.tree: Number of trees
  #   ntrace: running output is printed for every ntrace trees.
  #
  # Returns:
  #  1. Prediction mse on test data; 2.RF model object
  sol <- NULL
  pred <- NULL
  
  # Split data into train/test set based on cutoff date
  train <- dat[dat$Date.Production.Start<=cutoff, ]
  test  <- dplyr::setdiff(dat, train)
  
  print(paste0("Cutoff date=", cutoff, " Training%=", nrow(train)/nrow(dat)), quote=F)
  #################################################################################################
  rf.model <- randomForest(model, data=train, mtry=m, do.trace=ntrace, ntree=no.tree)  
  #################################################################################################
  
  # Predict test dataset and calculate mse
  test.pred <- cbind(test[,c(1,33)],Pred=predict(rf.model, newdata=test))  # Uwi, Target, Pred
  mse <- sum((test.pred[,2]-test.pred[,3])^2)/nrow(test.pred)  # mean square of residuals
  
  # Prediction results (combine observed with predict)
  pred <- rbind(cbind(train[,c(1,33)], Pred=train[,33]), test.pred)
  
  # Accuracy on test dataset at given cutoff date
  sol <- data.frame(cutoff=cutoff, train.pct=nrow(train)/nrow(dat), mse=mse, rmse=sqrt(mse), m=m, n.Tree=no.tree)
  
  return(list(sol, pred)) 
}



###############################################################################
runRFReg3 <- function(dat, cutoff, model, m, no.tree, ntrace=500){
  # Run RF model on regression with given pars 
  #
  # Args:
  #   dat: data set for rf, will be split to train and test set
  #   cutoff: cutoff date, wells before cutoff date will be used for training
  #   model: formula for RF model
  #   m: Number of variables randomly sampled as candidates at each split; mtry
  #   no.tree: Number of trees
  #   ntrace: running output is printed for every ntrace trees.
  #
  # Returns:
  #  1. Prediction mse on test data; 2.RF model object
  sol <- NULL
  pred <- NULL
  
  # Split data into train/test set based on cutoff date
  train <- dat[dat$Date.Production.Start<=cutoff, ]
  test  <- dplyr::setdiff(dat, train)
  
  print(paste0("Cutoff date=", cutoff, " Training%=", nrow(train)/nrow(dat)), quote=F)
  #################################################################################################
  rf.model <- randomForest(model, data=train, mtry=m, do.trace=ntrace, ntree=no.tree)  
  #################################################################################################
  
  # Predict test dataset and calculate mse
  test.pred <- cbind(test[,c(1,35)],Pred=predict(rf.model, newdata=test))  # Uwi, Target, Pred
  mse <- sum((test.pred[,2]-test.pred[,3])^2)/nrow(test.pred)  # mean square of residuals
  
  # Prediction results (combine observed with predict)
  pred <- rbind(cbind(train[,c(1,35)], Pred=train[,35]), test.pred)
  
  # Accuracy on test dataset at given cutoff date
  sol <- data.frame(cutoff=cutoff, train.pct=nrow(train)/nrow(dat), mse=mse, rmse=sqrt(mse), m=m, n.Tree=no.tree)
  
  return(list(sol, pred)) 
}



###############################################################################
runRFReg4 <- function(dat, cutoff, model, m, no.tree, ntrace=500){
  # Run RF model on regression with given pars 
  #
  # Args:
  #   dat: data set for rf, will be split to train and test set
  #   cutoff: cutoff date, wells before cutoff date will be used for training
  #   model: formula for RF model
  #   m: Number of variables randomly sampled as candidates at each split; mtry
  #   no.tree: Number of trees
  #   ntrace: running output is printed for every ntrace trees.
  #
  # Returns:
  #  1. Prediction mse on test data; 2.RF model object
  sol <- NULL
  pred <- NULL
  
  # Split data into train/test set based on cutoff date
  train <- dat[dat$Date.Production.Start<=cutoff, ]
  test  <- dplyr::setdiff(dat, train)
  
  print(paste0("Cutoff date=", cutoff, " Training%=", nrow(train)/nrow(dat)), quote=F)
  #################################################################################################
  rf.model <- randomForest(model, data=train, mtry=m, do.trace=ntrace, ntree=no.tree)  
  #################################################################################################
  
  # Predict test dataset and calculate mse
  test.pred <- cbind(test[,c(1,36)],Pred=predict(rf.model, newdata=test))  # Uwi, Target, Pred
  mse <- sum((test.pred[,2]-test.pred[,3])^2)/nrow(test.pred)  # mean square of residuals
  
  # Prediction results (combine observed with predict)
  #pred <- rbind(cbind(train[,c(1,36)], Pred=train[,36]), test.pred)
  
  # Accuracy on test dataset at given cutoff date
  sol <- data.frame(cutoff=cutoff, train.pct=nrow(train)/nrow(dat), mse=mse, rmse=sqrt(mse), m=m, n.Tree=no.tree)
  
  return(list(sol, test.pred)) # return out of sample prediction
}



####################################################################################################
runRFRegCV <- function(dat, model, m, no.tree, k, rev=FALSE, ntrace=500, default=FALSE){
  # Run RF regression model with Cross Validation
  #
  # Args:
  #   dat: data set for rf, will be split to train and test set
  #   model: formula for RF model
  #   m: Number of variables randomly sampled as candidates at each split; mtry
  #   no.tree: Number of trees
  #   k:   K-fold CV 
  #   rev: reverse CV or not
  #   ntrace: running output is printed for every ntrace trees.
  #   default: using default setting?
  # Returns:
  #   1. Prediction mse based on CV  2. Prediction results for each fold
  
  folds <- cvFolds(nrow(dat), K=k)
  mse <- NULL;  pred <- NULL; sol <- NULL;
  
  for(i in 1:k){  
    # Split data into train/test set
    if(rev==TRUE){
      train <- dat[folds$subsets[folds$which==i],]
      test  <- dplyr::setdiff(dat, train)
    } else {
      test  <- dat[folds$subsets[folds$which==i],]
      train <- dplyr::setdiff(dat, test)
    }
    
    #####################################################################################################
    if(default==TRUE){
    	rf.model <- randomForest(model, data=train, importance=T, do.trace=ntrace, ntree=no.tree) 
    } else {
       rf.model <- randomForest(model, data=train, importance=T, mtry=m, do.trace=ntrace, ntree=no.tree)  
    }
    #####################################################################################################
    
    # Predict test dataset and calculate mse
    test.pred <- cbind(test[,c(1,33)], Pred=predict(rf.model,newdata=test), test[,c(38,39)])  # Uwi, Target, Pred, Latitude, Longitude
    #mse <- mse + sum((test.pred[,2]-test.pred[,3])^2)/nrow(test.pred)  # mean square of residuals
    mse <- c(mse, sum((test.pred[,2]-test.pred[,3])^2)/nrow(test.pred))
    pred <- rbind(pred, test.pred)  # save prediction results for fold i
    
    print(paste0("K=", k, " i=", i, " Rev=", rev), quote=F)
  }
  # CV results
  m <- rf.model$mtry  # get default value of mtry
  rmse <- sqrt(mse)
  #sol <- data.frame(K=k, Rev=rev, mse=mse/k, m=m, n.Tree=no.tree)
  sol <- data.frame(K=k, Rev=rev, mse=mean(mse), mse.sd=sd(mse), rmse=mean(rmse), rmse.sd=sd(rmse), m=m, n.Tree=no.tree)
  
  if(rev==TRUE){
    pred <- pred %>% group_by(Uwi) %>% summarise(Pred=mean(Pred))  # average reverse CV prediction
    target <- select(dat, Uwi, Target, Latitude, Longitude)
    pred <- left_join(target, pred, by="Uwi")
    pred <- select(pred, Uwi, Target, Pred, Latitude, Longitude)
  } 
  
  return(list(sol, pred))
  
}







####################################################################################################
runBartRegCV <- function(dat, no.tree, no.burn, no.after.burn, k.fold=5, k=2, q=0.9, nu=3, rev=FALSE, seed=999){
  # Run RF regression model with Cross Validation
  #
  # Args:
  #   dat: data set for bart, will be split to train and test set
  #   model: formula for RF model
  #   m: Number of variables randomly sampled as candidates at each split; mtry
  #   no.tree: Number of trees
  #   k:   K-fold CV 
  #   rev: reverse CV or not
  #   ntrace: running output is printed for every ntrace trees.
  #   default: using default setting?
  # Returns:
  #   1. Prediction mse based on CV  2. Prediction results for each fold
  
  folds <- cvFolds(nrow(dat), K=k.fold)
  mse <- 0;
  pred <- NULL; sol <- NULL;
  
  for(i in 1:k.fold){  
    # Split data into train/test set
    if(rev==TRUE){
      train <- dat[folds$subsets[folds$which==i],]
      test  <- dplyr::setdiff(dat, train)
    } else {
      test  <- dat[folds$subsets[folds$which==i],]
      train <- dplyr::setdiff(dat, test)
    }
    
  ##########################################################################################################################################################################
    bart <- bartMachine(train[,2:32], train[,33], num_trees=no.tree, num_burn_in=no.burn, num_iterations_after_burn_in=no.after.burn, k=k, q=q, nu=nu, seed=666)  # for regression, more trees seems more accurate
  ##########################################################################################################################################################################
    
    # Predict test dataset and calculate mse
    test.pred <- cbind(test[,c(1,33)], Pred=predict(bart,test[,2:32]), test[,c(36,37)])  # Uwi, Target, Pred, Latitude, Longitude
  
    mse <- mse + sum((test.pred[,2]-test.pred[,3])^2)/nrow(test.pred)  # mean square of residuals
    pred <- rbind(pred, test.pred)  # save prediction results for fold i
    
    print(paste0("K=", k.fold, " i=", i, " Rev=", rev), quote=F)
  }
  # CV results
  sol <- data.frame(K=k.fold, Rev=rev, mse=mse/k.fold, k=k, q=q, nu=nu, n.Tree=no.tree, seed=seed)
  
  return(list(sol, pred))
  
}


#######################################################################################################
# function that computes the fraction of correctly identified entries in the top n%-quantile for all n
#
#usage
# qrecoverycurve(x)
#
#arguments
#  x               a data frame with at least 2 columns that are all at least ordinal
#                   column 1 contains true values, the other columns contains predictions that are rank proxies
#
#value
#  a data frame; 
#   the i-th entry in the j-th column is the fraction of rows that are 
#    correctly (baseline = true rank) identified by the proxy in column j to rank amongst the first i rows   
#     rows with NAs in the first two columns are omitted from this computation
#   the 1st column contains as baseline the fraction which is achieved in expectation by optimal uninformed guessing
#
# example
#  x <- cbind(1:11,c(2:4,5,7,9:4))
#  recfreq <- qrecoverycurve(x)
#  plot(recfreq,type = "l",xlab = "baseline recovery", ylab = "method recovery")
#  plot(1:10,recfreq[,2],type = "l",xlab = "size of top quantile", ylab = "recovery rate")


qRecCurv <- function(x) {
  
  x <- as.data.frame(na.omit(x))
  
  n.row.x <- nrow(x)  
  n.col.x <- ncol(x)  
  
  ranks <- x %>% mutate_each(funs(row_number)) %>% arrange(desc(Target))  # ranks for each col and then ordered by 1st col(true value)
  
  rec.q <- data.frame(matrix(-1, nrow = n.row.x , ncol = n.col.x))  # recover quantiles
  rec.q[1,] <- (ranks[1,] == n.row.x)
  for (i in 2:n.row.x)
  {
    #rec.q[i,] <- ranks %>% slice(1:i) %>% summarise_each (funs(sum(.<=i)/i))
    rec.q[i,] <- ranks %>% slice(1:i) %>% summarise_each (funs(sum(.>=(n.row.x-i+1))/i))
  }
  names(rec.q)[1]<- "True"
  rec.q[,1]<-1:n.row.x/n.row.x
  
  #row.names(rec.q) <- sapply(100*(1:n.row.x)/n.row.x,  FUN = function(x) paste("P",round(x,digits = 0),sep = ""))
  
  return(rec.q)
  
}




#######################################################################################################
# function that computes the fraction of correctly identified entries in the top n%-quantile for all n
#
#usage
# qrecoverycurve(x)
#
#arguments
#  x               a data frame with at least 2 columns that are all at least ordinal
#                   column 1 contains true values, the other columns contains predictions that are rank proxies
#
#value
#  a data frame; 
#   the i-th entry in the j-th column is the fraction of rows that are 
#    correctly (baseline = true rank) identified by the proxy in column j to rank amongst the first i rows   
#     rows with NAs in the first two columns are omitted from this computation
#   the 1st column contains as baseline the fraction which is achieved in expectation by optimal uninformed guessing
#
# example
#  x <- cbind(1:11,c(2:4,5,7,9:4))
#  recfreq <- qrecoverycurve(x)
#  plot(recfreq,type = "l",xlab = "baseline recovery", ylab = "method recovery")
#  plot(1:10,recfreq[,2],type = "l",xlab = "size of top quantile", ylab = "recovery rate")


qRecCurvBottom <- function(x) {
  
  x <- as.data.frame(na.omit(x))
  
  n.row.x <- nrow(x)  
  n.col.x <- ncol(x)  
  
  ranks <- x %>% mutate_each(funs(row_number)) %>% arrange(Target)  # ranks for each col and then ordered by 1st col(true value)
  
  rec.q <- data.frame(matrix(-1, nrow = n.row.x , ncol = n.col.x))  # recover quantiles
  rec.q[1,] <- (ranks[1,] == 1)
  for (i in 2:n.row.x)
  {
    rec.q[i,] <- ranks %>% slice(1:i) %>% summarise_each (funs(sum(.<=i)/i))
  }
  names(rec.q)[1]<- "True"
  rec.q[,1]<-1:n.row.x/n.row.x
  
  #row.names(rec.q) <- sapply(100*(1:n.row.x)/n.row.x,  FUN = function(x) paste("P",round(x,digits = 0),sep = ""))
  
  return(rec.q)
  
}

