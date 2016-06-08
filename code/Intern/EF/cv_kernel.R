######## cross validation for kernel #############

#loc_variables is the variables names for longitude and latitude
#varname considered variables
#index_log has the same length of varname, which indicates whether a log transformation is done: 0 means log transform, 1 means original scale
# method is the method in variog, default is "cressie" here
#K number of folds, default is 10
#seed with default value 2043
#lamv and sigv are tuning parameter values for selection
#---------------------------------------------------------------

library(dplyr)
source("code/kernel_est.R")
toolbox_path = "Z:/XinWang/project/cokriging/code/kernel/Tools - validation toolbox"
cv_kernel <- function(dat, loc_variables,
                      varname, index_log, lamv = 2^(-5:5), sigv = 2^(-5:5),
                      K = 10, seed = 2043)
{
  source(file.path(toolbox_path, "/validation_toolbox.R"))
  gaussParamGrid <- tuning.grid(gridType = "expand", 
                                paramNames = c("lambda","kernel"), 
                                gridParams = list(lambda = lamv, kernel = lapply(sigv,function(x){rbfdot(sigma = x)}))
  )
  
  
  nvar <- length(varname)

  dat <- dat[,c(loc_variables,varname)]
  num0 <- colSums(!is.na(dat[,varname]))  # number of nonmissing values for each variable
  
  
  test_error <- rep(0,nvar) #output
  names(test_error) <- varname
  
  # construct K-folds data 
  set.seed(seed)
  cv_group <- cvFolds(n = nrow(dat),K = K,type = "random")

  for(k in 1:K)
  {
    test_dat  <- dat[cv_group$subsets[cv_group$which==k],]
    train_dat <- dplyr::setdiff(dat, test_dat)
    
	# log transformation
    for(idi in 1:nvar)
    {
      if(index_log[idi] ==0){train_dat[,varname[idi]] <- log(train_dat[,varname[idi]])
	                         train_dat[is.infinite(train_dat[,varname[idi]]),varname[idi]] <- NA  # set Inf to NA
	                         }
      
    }
   
    
    sd0 <- apply(train_dat[,varname],2,function(x)(sd(x,na.rm = TRUE)))  # standard deviation 
    mean0 <- colMeans(train_dat[,varname],na.rm = TRUE) # mean 
    
    train_dat[,varname] <- scale(train_dat[,varname])   
    
    #cross validation to select lambda and sigma 
    foldSetup <- initFolds(foldType = "kfold", seed = seed, foldParams = list(k = K))
    dataFolds <- sampleFolds(nrow(train_dat), foldSetup)
    
    # select best lambda and sigma
    tuningErrors <- tuning.computeErrors(train_dat, dataFolds, predictor.gaussreg, gaussParamGrid ,
                                         target = varname, covariates = loc_variables,use_foreach = FALSE)
    
    paraselect <- tuning.getBestParams(tuningErrors, gaussParamGrid)
    
  #for each variable fit model and predict on test data # 
    for(i in 1:nvar)
    {
      varnamei <- varname[i]
      index_obsi <- !is.na(train_dat[,varnamei])
      fitmodel <- fit_kernel(x = train_dat[index_obsi,loc_variables],y = train_dat[index_obsi,varnamei],params = paraselect[[varnamei]])
      pred <- predict_kernel(modelFit = fitmodel,newdata = test_dat[,loc_variables],params = paraselect[[varnamei]])
	  
	  # transform back
      if(index_log[i]==0)
      {
        pred <- exp(pred*sd0[i] + mean0[i])
      }else{
        pred <- pred*sd0[i] + mean0[i]
      }
      
      test_error[i] <-  test_error[i] + sum((test_dat[,varnamei] - pred[,1])^2,na.rm = TRUE)/num0[i]
    }
    
  }
  
  test_error <- sqrt(test_error)
  return(test_error)
}