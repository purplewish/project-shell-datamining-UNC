#####kriging based on variogram estimation (revised of last year)################
# fix.nugget = TRUE, nugget is fixed which is from variogram, fix.nugget = FALSE, nugget is estimated via MLE
# -------------------------------------
#fitpred_old_fixnugget gives predicted values in test_dat based on train_dat
#kriging_old_fixnugget gives test error in cross validation based on given kappa value, which is used in tuning function to select tuning parameter
#tuning_kriging_old_fixnugget: select tuning parameter
#cv_kriging_old_fixnugget: test error based on selected tuning parameter
# kappa is tuning parameter(smoothness parameter)

#--------------------------------------------------------------------
#loc_variables is the variables names for longitude and latitude
#varname: considered variables
#varname_2nd : variables for trend = "2nd" in last year
#varname_ord : variables for trend ="cte" in last year
#index_log has the same length of varname, which indicates whether a log transformation is done: 0 means log transform, 1 means original scale, shold have the same order as varname
# method is the method in variog, default is "cressie" here
#kappa_vec has the same length with varname, smoothness parameter in Matern model
#const: maximum distance in variogram is const*maximum distance in data

#K is the number of folds, default is K = 10
#the default seed is 2403
#------------------------------------------------------------------------

fitpred_old_fixnugget <- function(train_dat, test_dat, loc_variables, fix.nugget = TRUE,
                        varname, varname_2nd, varname_ord, index_log,
                        method = "cressie", kappa_vec = rep(0.5,length(varname)), const = 0.64)

{

  predmat <- matrix(0,nrow(test_dat),length(varname))
  colnames(predmat) <- varname
  names(kappa_vec) <- varname
  names(index_log) <- varname

  # 2nd order polynomial
  for (i in 1:length(varname_2nd))
  {
    varnamei <- varname_2nd[i]
    num0 <- sum(train_dat[,varnamei] == 0,na.rm = TRUE)
    # exclude the zero values for variables with log transformation in training data
    if(index_log[varnamei] == 0 & num0 >0)
    {
      index_obsi <- (!is.na(train_dat[,varnamei])) & ( train_dat[,varnamei] >0)
    } else{
      index_obsi <- !is.na(train_dat[,varnamei])
    }

    maxdist <- max(dist(train_dat[index_obsi,loc_variables]))
    d <-  maxdist*const
    #empirical variogram
    vgi <-variog(coords=train_dat[index_obsi,loc_variables],data=train_dat[index_obsi,varnamei],trend='2nd',max.dist=d,lambda = index_log[varnamei]) # 2nd means quadratic term
 # fit variogram with fixed kappa
  fit.vg <-variofit(vgi,cov.model = "matern", fix.nugget = FALSE,nugget = vgi$v[1],weights = method,fix.kappa = TRUE,kappa = kappa_vec[varnamei])
    # by default, the model is exponential

  if(fit.vg$nugget >0)
  {
    nugget0 <- fit.vg$nugget
  }else{
    nugget0 <- 0
  }

  #  if range is greater than maxdist, set it to maxdist in initial value
  #estimate parameters with fixed kappa g, kappa =0.5 means exponential model
  # if fix.nugget = TRUE, the value of nugget effect is fixed and is from the fit.vg
  # if fix.nugget = FALSE, then the nugget is estimated by maximizing log lik, and the initial value is the value from fit.vg
  if(fit.vg$cov.pars[2] > maxdist)
  {
    # the estimates of  cov.pars based on fit.vg are too large
    res_liki <- likfit(coords=train_dat[index_obsi,loc_variables],data = train_dat[index_obsi,varnamei],trend='2nd',fix.nugget = fix.nugget,nugget = nugget0, fix.kappa=fix.nugget,kappa =  kappa_vec[varnamei], lambda = index_log[varnamei],ini.cov.pars = c(vgi$v[length(vgi$v)] - vgi$v[1], maxdist),cov.model = "matern")
  }else{
    res_liki <- likfit(coords=train_dat[index_obsi,loc_variables],data = train_dat[index_obsi,varnamei],trend='2nd',fix.nugget = fix.nugget,nugget = nugget0, fix.kappa=fix.nugget,kappa =  kappa_vec[varnamei], lambda = index_log[varnamei],ini.cov.pars = fit.vg$cov.pars  ,cov.model = "matern")
  }

  # Kriging on test
    kgi = krige.conv(coords=train_dat[index_obsi,loc_variables],data = train_dat[index_obsi,varnamei],loc = test_dat[,loc_variables], krige = krige.control(obj.m = res_liki,lambda=index_log[varnamei], trend.d = "2nd", trend.l = "2nd", type.krige = "SK"))

    predmat[,varnamei] <-  kgi$predict

  }

  # Ordinary
  for (i in 1:length(varname_ord))
  {
    varnamei <- varname_ord[i]
    num0 <- sum(train_dat[,varnamei] == 0,na.rm = TRUE)
    # exclude the zero values for variables with log transformation in training data
    if(index_log[varnamei] == 0 & num0 >0)
    {
      index_obsi <- (!is.na(train_dat[,varnamei])) & ( train_dat[,varnamei] >0)
    } else{
      index_obsi <- !is.na(train_dat[,varnamei])
    }

    maxdist <- max(dist(train_dat[index_obsi,loc_variables]))
    d <-  maxdist*const

    vgi <- variog(coords=train_dat[index_obsi,loc_variables],data=train_dat[index_obsi,varnamei],trend='cte',max.dist=d)  # cte=constant
    fit.vg <-variofit(vgi,cov.model = "matern", fix.nugget = FALSE,nugget = vgi$v[1],weights = method,fix.kappa = TRUE,kappa = kappa_vec[varnamei])

    if(fit.vg$nugget >0)
    {
      nugget0 <- fit.vg$nugget
    }else{
      nugget0 <- 0
    }

    if(fit.vg$cov.pars[2] > maxdist)
    {
       #the estimates of  cov.pars based on fit.vg are too large
      res_liki <- likfit(coords=train_dat[index_obsi,loc_variables],data = train_dat[index_obsi,varnamei],trend='cte',fix.nugget = fix.nugget,nugget = nugget0, fix.kappa=fix.nugget,kappa =  kappa_vec[varnamei], lambda = index_log[varnamei],ini.cov.pars = c(vgi$v[length(vgi$v)] - vgi$v[1], maxdist),cov.model = "matern")


    }else{
      res_liki <- likfit(coords=train_dat[index_obsi,loc_variables],data = train_dat[index_obsi,varnamei],trend='cte',fix.nugget = fix.nugget,nugget = nugget0, fix.kappa=fix.nugget,kappa =  kappa_vec[varnamei], lambda = index_log[varnamei],ini.cov.pars = fit.vg$cov.pars,cov.model = "matern")
    }

      kgi = krige.conv(coords=train_dat[index_obsi,loc_variables],data = train_dat[index_obsi,varnamei],loc = test_dat[,loc_variables], krige = krige.control(obj.m = res_liki,lambda=index_log[varnamei], type.krige = "SK"))

    predmat[,varnamei] <-  kgi$predict
  }

  return(predmat)
}

# For given kappa values, gives corresponding test error based on cross validation
kriging_old_fixnugget <- function(dat, loc_variables, fix.nugget = TRUE,
                        varname, varname_2nd, varname_ord, index_log,
                        method = "cressie", kappa_vec = rep(0.5,length(varname)),
                        K, seed = 2043, const = 0.64)
{
  set.seed(seed)
  # construct K-folds data
  cv_group <- cvFolds(n = nrow(dat),K = K,type = "random")

  nvar <- length(varname)
  miss_inf <- colSums(!is.na(dat[,varname]))   # number of nonmissing values for each varaible

  test_error0  <- rep(0,nvar) #output
  names(test_error0)  <- varname
  names(kappa_vec) <- varname

  for(k in 1:K)
  {
    test_dat  <- dat[cv_group$subsets[cv_group$which==k],]
    train_dat <- dplyr::setdiff(dat, test_dat)

    predk <- fitpred_old_fixnugget(train_dat, test_dat, loc_variables = loc_variables, fix.nugget = fix.nugget,
                                   varname = varname,varname_2nd = varname_2nd,
                                   varname_ord = varname_ord, index_log = index_log,
                         method = method, kappa_vec = kappa_vec, const = const) # prediction on test

    test_error0 <- test_error0 + colSums((predk - test_dat[,varname])^2,na.rm = TRUE)/miss_inf
  }

  test_error0 <- sqrt(test_error0) #rmse
  return(test_error0)
}

print(c("fitpred_old_fixnugget","kriging_old_fixnugget"))



# select tuning parameter
tuning_kriging_old_fixnugget <- function(dat,kappav = seq(0.2,3,0.3),
                           loc_variables, fix.nugget = TRUE,varname, varname_2nd, varname_ord,
                           index_log, method = "cressie",
                           K, seed = 2043, const = 0.64)
{

  tab0 <- data.frame(kappa = rep(kappav,each = length(varname)),variable = rep(varname,length(kappav)),rmse = 0)

  # for each kappa value, calculate test RMSE
  for(m in 1:length(kappav))
  {
    resm <- rep(0,length(varname))
    # tryCatch can skip the error, the loop will not stop because of the error
    tryCatch({    resm <-  kriging_old_fixnugget(dat = dat,loc_variables = loc_variables, fix.nugget = fix.nugget,
                                                 varname = varname,varname_2nd = varname_2nd,
                                                 varname_ord = varname_ord,
                                                 index_log = index_log,  method = method,
                                                 kappa_vec = rep(kappav[m],length(varname)),
                                                 K = K, seed = seed, const = const)},
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
cv_kriging_old_fixnugget <- function(dat, loc_variables, fix.nugget = TRUE,
                       varname, varname_2nd, varname_ord, index_log,
                       method = "cressie", kappav = seq(0.2,2.5,0.3),
                       K = 10, seed = 2043, const = 0.64)
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
    kappa_old <- tuning_kriging_old_fixnugget(train_dat,kappav = kappav,
                                loc_variables = loc_variables, fix.nugget = fix.nugget,
                                varname = varname, varname_2nd, varname_ord = varname_ord,
                                index_log = index_log,
                                method = method,
                                K = K, seed = seed, const = const)

    # prediction on test data
    pred0k <- fitpred_old_fixnugget(train_dat,test_dat,loc_variables = loc_variables, fix.nugget = fix.nugget,
                          varname = varname,varname_2nd = varname_2nd, varname_ord = varname_ord,
                          index_log = index_log, method = method,kappa_vec = kappa_old, const = const)

    test_error0 <- test_error0 + colSums((pred0k - test_dat[,varname])^2, na.rm = TRUE)/num0_all

  }

  test_error0 <- sqrt(test_error0)

  return(test_error0)
}

print(c("tuning_kriging_old_fixnugget","cv_kriging_old_fixnugget"))
