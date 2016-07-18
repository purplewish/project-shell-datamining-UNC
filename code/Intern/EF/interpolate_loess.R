#################### local polynomial ############################
# fitpred_loess: fit and predict for given span value 
#error_fun: give test error for given span value
#span_fun: a function to select tuning parameters called span
#cv_loess: cross validation
library(cvTools)
library(dplyr)
fitpred_loess <- function(train_dat, test_dat, loc_variables, varname, span0)
{
  predmat <- matrix(0,nrow(test_dat),length(varname))
  colnames(predmat) <- varname
  nvar <- length(varname)
  names(span0) <- varname
  for(j in 1:nvar)
  {
    indexj <- !is.na(train_dat[,varname[j]])
    formuj <- as.formula(paste(varname[j],"~",paste(loc_variables,collapse = "+")))
    resj <- loess(formuj, data = train_dat[indexj,], span = span0[varname[j]], control=loess.control(surface = "direct")) 
    predj <- predict(resj,newdata = test_dat[,loc_variables])
    predmat[,j] <- predj
  }
  return(predmat)
}


error_fun <- function(dat, loc_variables, varname, span0, K = 10, seed = 2043)
{
  set.seed(seed)
  # construct K-folds data 
  cv_group <- cvFolds(n = nrow(dat),K = K,type = "random")
  
  nvar <- length(varname)
  miss_inf <- colSums(!is.na(dat[,varname]))   # number of nonmissing values for each varaible
  
  test_error0  <- rep(0,nvar) #output
  names(test_error0)  <- varname
  names(span0) <- varname
  
  
  for(k in 1:K)
  {
    test_dat  <- dat[cv_group$subsets[cv_group$which==k],]
    train_dat <- dplyr::setdiff(dat, test_dat)
    
    predk <- fitpred_loess(train_dat, test_dat, loc_variables = loc_variables, 
                                   varname = varname, span0 = span0 ) # prediction on test
    
    test_error0 <- test_error0 + colSums((predk - test_dat[,varname])^2,na.rm = TRUE)/miss_inf
  }
  return(test_error0)
}


span_fun <- function(dat, spanv = seq(0.5,1,by =0.1),
                     loc_variables, varname, K=10 ,seed = 2043)
{
  tab0 <- data.frame(span = rep(spanv,each = length(varname)),variable = rep(varname,length(spanv)),rmse = 0)
  for(m in 1:length(spanv))
  {
    tab0[tab0$span == spanv[m],3] <-  error_fun(dat = dat,loc_variables = loc_variables,
                                                              varname = varname,
                                                              span0= rep(spanv[m],length(varname)), 
                                                              K = K, seed = seed)
  }
  
  
  tab0 <- group_by(tab0,variable)
  tabnew <- summarize(tab0, span = spanv[which.min(rmse)]) # kappa with smallest RMSE 
  span_select<- as.matrix(tabnew[,2])[,1]
  names(span_select) <- as.data.frame(tabnew[,1])[,1]
  span_select <- span_select[varname] # reorder to match original order 
  return(span_select)
  
}

  
cv_loess <- function(dat, loc_variables,
                     varname, spanv = seq(0.5,1,0.1),
                     K = 10, seed = 2043)
{
  nvar <- length(varname)
  
  dat <- dat[,c(loc_variables,varname)]
  num0_all <- colSums(!is.na(dat[,varname]))   # number of nonmissing values for each varaible
  
  test_error0 <- rep(0,nvar) #output
  names(test_error0)  <- varname
  
  # construct K-folds data 
  set.seed(seed)
  cv_group <- cvFolds(n = nrow(dat),K = K,type = "random")
  
  
  for(k in 1:K)
  {
    test_dat  <- dat[cv_group$subsets[cv_group$which==k],]
    train_dat <- dplyr::setdiff(dat, test_dat)
    
    # select tuning parameter based on training data
    span_select <- span_fun(train_dat,loc_variables = loc_variables,
                         varname = varname, spanv = spanv,
                         K = K, seed = seed)
    
    # prediction on test data
    predk <- fitpred_loess(train_dat,test_dat,loc_variables = loc_variables,
                                    varname = varname,span0 = span_select)
    
    test_error0 <- test_error0 + colSums((predk - test_dat[,varname])^2, na.rm = TRUE)/num0_all
    
  }
  
  test_error0 <- sqrt(test_error0)
  return(test_error0)
}
  
