
library(dplyr)
library(bartMachine)

#######################################################################################################
runBart <- function(dat, train.pct, no.tree=50, no.burnin=250, no.afterburnin=1000, nrep=1, seed=666){
  # Run BART model n times with given pars based on full dataset
  #
  # Args:
  #   dat: data set for BART, will be split to train and test set
  #   train.pct: training data percentage 
  #   no.tree: Number of trees
  #   no.burnin: Number of burnin steps
  #   no.afterburnin: Number of steps after burnin 
  #   nrep: Number of replication to run BART model under the same set of pars
  #   seed: running seed
  #
  # Returns:
  #  1. Prediction accuracy on test data; 2. Last BART model object; 3. Prediction object
  sol <- NULL
  # Prediction data frame (Accy: classification accuracy; TP/TN: True+/True-)
  pred.df <- data.frame(Index=1:nrep,Accy=NA,TP=NA,TN=NA) 
  set_bart_machine_num_cores(1)
  
  set.seed(seed)  
  for(i in 1:nrep){
    print(paste0("i=",i, " Train%=", train.pct), quote=F)
    
    # Split data into train/test set
    train <- sample_frac(dat, train.pct, replace=F)
    test  <- dplyr::setdiff(dat, train)
    
    #####################################################################################################################################################
    bart = bartMachine(train[ ,2:32], train$Target.Q4, num_trees=no.tree, num_burn_in=no.burnin, num_iterations_after_burn_in=no.afterburnin, seed=seed)
    pred = bart_predict_for_test_data(bart, test[,2:32], test$Target.Q4, )
    #####################################################################################################################################################
        
    # Predict test dataset and calculate error rate
    pred.df$Accy[i] <- 1 - pred$confusion_matrix[3,3]  # Classification accuracy
    pred.df$TP[i] <- 1 - pred$confusion_matrix[2,3]  # TP
    pred.df$TN[i] <- 1 - pred$confusion_matrix[1,3]  # TN
  }
  
  # Avg accuracy over nrep on test dataset at fixed training %
  sol <- rbind(sol, c(Train.Perc=train.pct, Accy=mean(pred.df[,2]), TP=mean(pred.df[,3]), TN=mean(pred.df[,4]), n.Tree=no.tree, n.burnin=no.burnin, n.iterafterburnin=no.afterburnin, seed=seed) )
  
  return(list(sol, bart, pred)) 
}


#######################################################################################################
runBart2 <- function(train, test, no.tree=50, no.burnin=250, no.afterburnin=1000, seed=666){
  # Run BART model based on training and testing datesets
  #
  # Args:
  #   train: training data 
  #   test: testing data 
  #   no.tree: Number of trees
  #   no.burnin: Number of burnin steps
  #   no.afterburnin: Number of steps after burnin 
  #   seed: running seed
  #
  # Returns:
  #  1. Prediction accuracy on test data; 2. Last BART model object; 3. Prediction object
  
  #####################################################################################################################################################
  bart = bartMachine(train[ ,2:32], train$Target.Q4, num_trees=no.tree, num_burn_in=no.burnin, num_iterations_after_burn_in=no.afterburnin, seed=seed)
  pred = bart_predict_for_test_data(bart, test[,2:32], test$Target.Q4)
  #####################################################################################################################################################
  
  # Predict test dataset and calculate error rate
  pred.df$Accy <- 1 - pred$confusion_matrix[3,3]  # Classification accuracy
  pred.df$TP <- 1 - pred$confusion_matrix[2,3]  # TP
  pred.df$TN <- 1 - pred$confusion_matrix[1,3]  # TN
  
  sol <- c(Accy=pred.Accy, TP=pred.TP, TN=pred.TN, n.Tree=no.tree, n.burnin=no.burnin, n.iterafterburnin=no.afterburnin, seed=seed)
  
  return(list(sol, bart, pred)) 
}

