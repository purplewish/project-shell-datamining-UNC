########################## in variogram based and mle based univariate kriging method, the smoothness parameter is fixed, which can be viewed as tuning parameter #######
# tuning_kriging function can select tuning parameter based on cross validation ####
# cv_kriging function to calculate test error based on selected tuning parameter in train data
#----------------------------------------------------------------------------------------

source("code/kriging_old.R")
source("code/kriging_new.R")

tuning_kriging <- function(dat, algorithm = "old",kappav = seq(0.2,3,0.3),
                           loc_variables,varname, varname_2nd, varname_ord, index_log,method = "cressie",
                           K, seed = 2043)
{

  tab0 <- data.frame(kappa = rep(kappav,each = length(varname)),variable = rep(varname,length(kappav)),rmse = 0)
           
  if(algorithm  == "old")
  {  
  for(m in 1:length(kappav))
  {
    tab0[tab0$kappa == kappav[m],3] <-  kriging_old(dat = dat,loc_variables = loc_variables,
                                                    varname = varname,varname_2nd = varname_2nd, 
                                                    varname_ord = varname_ord, method = method,
                                                    kappa_vec = rep(kappav[m],length(varname)), 
                                                    K = K, seed = seed)
  }
  }
  
  if(algorithm  == "new")
  {  
    for(m in 1:length(kappav))
    {
      tab0[tab0$kappa == kappav[m],3] <-  kriging_new(dat = dat,loc_variables = loc_variables,
                                                      varname = varname,index_log = index_log,
                                                      method = method,
                                                      kappa_vec = rep(kappav[m],length(varname)), 
                                                      K = K,seed = seed)
 
    }
  }
  
  tab0 <- group_by(tab0,variable)
  tabnew <- summarize(tab0, kappa = kappav[which.min(rmse)]) # kappa with smallest RMSE 
  kappa_select <- as.matrix(tabnew[,2])[,1]
  names(kappa_select) <- as.data.frame(tabnew[,1])[,1]
  kappa_select <- kappa_select[varname] # reorder to match original order 
  return(kappa_select)
}

#cross validation, in each training data, tuning parameters are selcted and then do prediction in test data 
cv_kriging <- function(dat, loc_variables,
                       varname, varname_2nd, varname_ord, index_log,
                       method = "cressie", kappav = seq(0.2,2.5,0.3),
                       K = 10, seed = 2043)
{
  nvar <- length(varname)
  
  dat <- dat[,c("UWI",loc_variables,varname)]
  num0_all <- colSums(!is.na(dat[,varname]))   # number of nonmissing values for each varaible
  
  # 0 for old 1 for new 
  test_error0 <- test_error1 <- rep(0,nvar) #output
  names(test_error0) <- names(test_error1) <- varname
  
  # construct K-folds data 
  set.seed(seed)
  cv_group <- cvFolds(n = nrow(dat),K = K,type = "random")
  
  
  for(k in 1:K)
  {
    test_dat  <- dat[cv_group$subsets[cv_group$which==k],]
    train_dat <- dplyr::setdiff(dat, test_dat)
    
    # select tuning parameter based on training data
    kappa_old <- tuning_kriging(train_dat, algorithm = "old",kappav = kappav,
                    loc_variables = loc_variables,
                    varname = varname, varname_2nd, varname_ord = varname_ord, 
					          index_log = index_log,method = method,
                    K = K, seed = seed)
    kappa_new <- tuning_kriging(train_dat, algorithm = "new",kappav = kappav,
                    loc_variables = loc_variables,
                    varname = varname, varname_2nd, varname_ord = varname_ord, 
					          index_log = index_log,method = method,
                    K = K, seed = seed)
    
    # prediction on test data
    pred0k <- fitpred_old(train_dat,test_dat,loc_variables = loc_variables,
                         varname = varname,varname_2nd = varname_2nd, varname_ord = varname_ord,
                         method = method,kappa_vec = kappa_old)
   
    pred1k <- fitpred_new(train_dat, test_dat, loc_variables = loc_variables, 
                          varname = varname, index_log = index_log,
                          method = method, kappa_vec = kappa_new)
    
    test_error0 <- test_error0 + colSums((pred0k - test_dat[,varname])^2, na.rm = TRUE)/num0_all
    test_error1 <- test_error1 + colSums((pred1k - test_dat[,varname])^2, na.rm = TRUE)/num0_all
    
  }
  
  test_error0 <- sqrt(test_error0)
  test_error1 <- sqrt(test_error1)
  tab <- cbind(test_error0,test_error1)
  colnames(tab) <- c("old","new")
  return(tab)
}
  
  print(c("tuning_kriging","cv_kriging"))