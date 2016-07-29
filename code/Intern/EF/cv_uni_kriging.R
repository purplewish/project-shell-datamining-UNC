########################## in variogram based and mle based univariate kriging method, the smoothness parameter is fixed, which can be viewed as tuning parameter #######
# tuning_kriging function selects tuning parameter based on cross validation (old + new)####
# cv_kriging function  calculate cross validation error (old + new)
# tuning_kriging_old() selects tuning parameter based on cross validation based on kriging_old()
# cv_kriging_old() calculates cross validation error based on kriging_old()
# tuning_kriging_new() selects tuning parameter based on cross validation based on kriging_new()
# cv_kriging_new() calculates cross validation error based on kriging_new()
# index_log_old is the log transformation indicator for last year method
# index_log_new is the log transfromation indicator for this year
# fix.nugget: whether nugget is fixed when estimating parameters via MLE
#----------------------------------------------------------------------------------------

### last year and this year were written together in tuning_kriging and cv_kriging
# they can be separately used in tuning_kriging_old, cv_kriging_old, tuning_kriging_new and cv_kriging_new
source("code/kriging_old.R")
source("code/kriging_new.R")
library(cvTools)

# select kappa
tuning_kriging <- function(dat, algorithm = "old",kappav = seq(0.2,3,0.3),
                           loc_variables,varname, varname_2nd, varname_ord, fix.nugget = TRUE,
                           index_log_old, index_log_new, method = "cressie",
                           K, seed = 2043, const = 0.64)
{

  tab0 <- data.frame(kappa = rep(kappav,each = length(varname)),variable = rep(varname,length(kappav)),rmse = 0)


  # For each kappav value, RMSE is calucalted
  if(algorithm  == "old")
  {
  for(m in 1:length(kappav))
  {
    tab0[tab0$kappa == kappav[m],3] <-  kriging_old(dat = dat,loc_variables = loc_variables,
                                                    varname = varname,varname_2nd = varname_2nd,
                                                    varname_ord = varname_ord,
                                                    index_log = index_log_old,  method = method,
                                                    kappa_vec = rep(kappav[m],length(varname)),
                                                    K = K, seed = seed, const = const)
    print(kappav[m])
  }
  }

  if(algorithm  == "new")
  {
    for(m in 1:length(kappav))
    {
      tab0[tab0$kappa == kappav[m],3] <-  kriging_new(dat = dat,loc_variables = loc_variables,
                                                      varname = varname,index_log = index_log_new,fix.nugget = fix.nugget,
                                                      method = method,
                                                      kappa_vec = rep(kappav[m],length(varname)),
                                                      K = K,seed = seed, const = const)
    print(kappav[m])
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
                       varname, varname_2nd, varname_ord, index_log_old, index_log_new, fix.nugget = TRUE,
                       method = "cressie", kappav = seq(0.2,2.5,0.3),
                       K = 10, seed = 2043, const = 0.64)
{
  nvar <- length(varname)

  dat <- dat[,c(loc_variables,varname)]
  num0_all <- colSums(!is.na(dat[,varname, drop = FALSE]))   # number of nonmissing values for each varaible

  # 0 for last year 1 (old) for this year (new)
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
                    varname = varname, varname_2nd = varname_2nd, varname_ord = varname_ord,
					          index_log_old = index_log_old, index_log_new = index_log_new,
					          method = method,
                    K = K, seed = seed, const = const)
    kappa_new <- tuning_kriging(train_dat, algorithm = "new",kappav = kappav,
                    loc_variables = loc_variables,
                    varname = varname, varname_2nd = varname_2nd, varname_ord = varname_ord,
					          index_log_old = index_log_old, index_log_new = index_log_new, fix.nugget = fix.nugget,method = method,
                    K = K, seed = seed, const = const)

    # prediction on test data
    pred0k <- fitpred_old(train_dat,test_dat,loc_variables = loc_variables,
                         varname = varname,varname_2nd = varname_2nd, varname_ord = varname_ord,
                         index_log = index_log_old,method = method,kappa_vec = kappa_old, const = const)

    pred1k <- fitpred_new(train_dat, test_dat, loc_variables = loc_variables,
                          varname = varname, index_log = index_log_new, fix.nugget = fix.nugget,
                          method = method, kappa_vec = kappa_new, const = const)

    test_error0 <- test_error0 + colSums((pred0k - test_dat[,varname])^2, na.rm = TRUE)/num0_all
    test_error1 <- test_error1 + colSums((pred1k - test_dat[,varname])^2, na.rm = TRUE)/num0_all

  }

  test_error0 <- sqrt(test_error0)
  test_error1 <- sqrt(test_error1)
  tab <- cbind(test_error0,test_error1)
  colnames(tab) <- c("old","new")
  return(tab)
}




##### for last year only
# select kappa
tuning_kriging_old <- function(dat,kappav = seq(0.2,3,0.3),
                           loc_variables,varname, varname_2nd, varname_ord,
                           index_log_old, method = "cressie",
                           K, seed = 2043, const = 0.64)
{

  tab0 <- data.frame(kappa = rep(kappav,each = length(varname)),variable = rep(varname,length(kappav)),rmse = 0)


  # For each kappav value, RMSE is calucalted
  for(m in 1:length(kappav))
  {
    tab0[tab0$kappa == kappav[m],3] <-  kriging_old(dat = dat,loc_variables = loc_variables,
                                                    varname = varname,varname_2nd = varname_2nd,
                                                    varname_ord = varname_ord,
                                                    index_log = index_log_old,  method = method,
                                                    kappa_vec = rep(kappav[m],length(varname)),
                                                    K = K, seed = seed, const = const)
   # print(kappav[m])
  }

  tab0 <- group_by(tab0,variable)
  tabnew <- summarize(tab0, kappa = kappav[which.min(rmse)]) # kappa with smallest RMSE
  kappa_select <- as.matrix(tabnew[,2])[,1]
  names(kappa_select) <- as.data.frame(tabnew[,1])[,1]
  kappa_select <- kappa_select[varname] # reorder to match original order
  return(kappa_select)
}

#cross validation, in each training data, tuning parameters are selcted and then do prediction in test data
cv_kriging_old <- function(dat, loc_variables,
                       varname, varname_2nd, varname_ord, index_log_old,
                       method = "cressie", kappav = seq(0.2,2.5,0.3),
                       K = 10, seed = 2043, const = 0.64)
{
  nvar <- length(varname)

  dat <- dat[,c(loc_variables,varname)]
  num0_all <- colSums(!is.na(dat[,varname, drop = FALSE]))   # number of nonmissing values for each varaible

  # 0 for last year 1 (old) for this year (new)
  test_error0  <- rep(0,nvar) #output
  names(test_error0)  <- varname

  # construct K-folds data
  set.seed(seed)
  cv_group <- cvFolds(n = nrow(dat),K = K,type = "random")


  for(k in 1:K)
  {
    test_dat  <- dat[cv_group$subsets[cv_group$which==k],]
    train_dat <- dplyr::setdiff(dat, test_dat)

    # select tuning parameter based on training data
    kappa_old <- tuning_kriging_old(train_dat ,kappav = kappav,
                                loc_variables = loc_variables,
                                varname = varname, varname_2nd = varname_2nd, varname_ord = varname_ord,
                                index_log_old = index_log_old,
                                method = method,
                                K = K, seed = seed, const = const)

    # prediction on test data
    pred0k <- fitpred_old(train_dat,test_dat,loc_variables = loc_variables,
                          varname = varname,varname_2nd = varname_2nd, varname_ord = varname_ord,
                          index_log = index_log_old,method = method,kappa_vec = kappa_old, const = const)


    test_error0 <- test_error0 + colSums((pred0k - test_dat[,varname])^2, na.rm = TRUE)/num0_all

  }

  test_error0 <- sqrt(test_error0)
  return(test_error0)
}



###### for this year only
tuning_kriging_new <- function(dat, kappav = seq(0.2,3,0.3),
                           loc_variables,varname, index_log_new, fix.nugget = TRUE, method = "cressie",
                           K, seed = 2043, const = 0.64)
{

  tab0 <- data.frame(kappa = rep(kappav,each = length(varname)),variable = rep(varname,length(kappav)),rmse = 0)


  # For each kappav value, RMSE is calucalted

  for(m in 1:length(kappav))
  {
     resm <- rep(0,length(varname))
     tryCatch({resm  <-  kriging_new(dat = dat,loc_variables = loc_variables,
                                                    varname = varname,index_log = index_log_new, fix.nugget = fix.nugget,
                                                    method = method,
                                                    kappa_vec = rep(kappav[m],length(varname)),
                                                      K = K,seed = seed, const = const)},
                        error=function(e){cat("ERROR :",conditionMessage(e), "\n")})

     tab0[tab0$kappa == kappav[m],3] <- resm
  }

  tab0 <- tab0[tab0[,"rmse"]!=0,]
  tab0 <- group_by(tab0,variable)
  tabnew <- summarize(tab0, kappa = kappav[which.min(rmse)]) # kappa with smallest RMSE
  kappa_select <- as.matrix(tabnew[,2])[,1]
  names(kappa_select) <- as.data.frame(tabnew[,1])[,1]
  kappa_select <- kappa_select[varname] # reorder to match original order
  return(kappa_select)
}

#cross validation, in each training data, tuning parameters are selcted and then do prediction in test data
cv_kriging_new <- function(dat, loc_variables,
                       varname, index_log_new, fix.nugget = TRUE,
                       method = "cressie", kappav = seq(0.2,2.5,0.3),
                       K = 10, seed = 2043, const = 0.64)
{
  nvar <- length(varname)

  dat <- dat[,c(loc_variables,varname)]
  num0_all <- colSums(!is.na(dat[,varname, drop = FALSE]))   # number of nonmissing values for each varaible

  # 0 for last year 1 (old) for this year (new)
  test_error1 <- rep(0,nvar) #output
  names(test_error1) <- varname

  # construct K-folds data
  set.seed(seed)
  cv_group <- cvFolds(n = nrow(dat),K = K,type = "random")


  for(k in 1:K)
  {
    test_dat  <- dat[cv_group$subsets[cv_group$which==k],]
    train_dat <- dplyr::setdiff(dat, test_dat)

    # select tuning parameter based on training data
    kappa_new <- tuning_kriging_new(train_dat,kappav = kappav,
                                loc_variables = loc_variables,
                                varname = varname,index_log_new = index_log_new,  fix.nugget = fix.nugget,
                                method = method,
                                K = K, seed = seed, const = const)

    # prediction on test data
    pred1k <- fitpred_new(train_dat, test_dat, loc_variables = loc_variables, fix.nugget = fix.nugget,
                          varname = varname, index_log = index_log_new,
                          method = method, kappa_vec = kappa_new, const = const)


    test_error1 <- test_error1 + colSums((pred1k - test_dat[,varname])^2, na.rm = TRUE)/num0_all

  }

  test_error1 <- sqrt(test_error1)
  return(test_error1)
}


  print(c("tuning_kriging","cv_kriging","tuning_kriging_old","cv_kriging_old","tuning_kriging_new","cv_kriging_new"))
