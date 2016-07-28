###########kriging based on variogram estimation (revised of last year)################
# range parameter is fixed, which is from the estimation of variogram
#-----------------------------------------------------------------
#fitpred_old gives predicted values in test_dat based on train_dat
#kriging_old gives test error (RMSE) in cross validation based on given kappa value, which is used in tuning function to select tuning parameter
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
library(geoR)
library(fields)
fitpred_old <- function(train_dat, test_dat, loc_variables,
                        varname, varname_2nd, varname_ord, index_log,
                        method = "cressie", kappa_vec = rep(0.5,length(varname)),const = 0.64)

{

  predmat <- matrix(0,nrow(test_dat),length(varname))
  colnames(predmat) <- varname
  names(kappa_vec) <- varname
  names(index_log) <- varname

  # 2nd order polynomial
  for (i in 1:length(varname_2nd))
  {
    varnamei <- varname_2nd[i]
    num0 <- sum(train_dat[,varnamei] == 0,na.rm = TRUE)  #number of zeroes
    # exclude the zero values for variables with log transformation in training data
    if(index_log[varnamei] == 0 & num0 >0)
    {
      index_obsi <- (!is.na(train_dat[,varnamei])) & ( train_dat[,varnamei] >0)
    } else{
      index_obsi <- !is.na(train_dat[,varnamei])
    }

    if(index_log[varnamei] ==0){train_dat[,varnamei] <- log(train_dat[,varnamei])} # log transformation


    d <-  max(dist(train_dat[index_obsi,loc_variables]))*const

    # fit empirical variogram
    vgi <-variog(coords=train_dat[index_obsi,loc_variables],data=train_dat[index_obsi,varnamei],trend='2nd',max.dist=d) # 2nd means quadratic term
    # fit variogram curve based on eempirical variogram
    covpar<-variofit(vgi,cov.model = "matern", fix.nugget = FALSE,nugget = vgi$v[1],weights = method,fix.kappa = TRUE,kappa = kappa_vec[varnamei]) # by default, the model is exponential


    if(covpar$cov.pars[2]==0) # range parameter if=0 means no dependence
    {covpar$cov.pars[2]=0.001}

    #Kriging, range (theta) parameter is from covpar object
    model <- Krig(x=train_dat[index_obsi,loc_variables], Y=train_dat[index_obsi,varnamei],theta=covpar$cov.pars[2],m=3,Covariance="Matern",smoothness=kappa_vec[varnamei]) # range parameter is fixed, other parameters are estimated
    pred <- predict(model,as.matrix(test_dat[,loc_variables]))
    if(index_log[varnamei] ==0){pred <- exp(pred)}
    predmat[,varnamei] <- pred
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

    if(index_log[varnamei] ==0){train_dat[,varnamei] <- log(train_dat[,varnamei])}


    d <-  max(dist(train_dat[index_obsi,loc_variables]))*const
    vgi <-variog(coords=train_dat[index_obsi,loc_variables],data=train_dat[index_obsi,varnamei],trend='cte',max.dist=d)  # cte=constant
    covpar<-variofit(vgi,cov.model = "matern", fix.nugget = FALSE,nugget = vgi$v[1],weights = method,fix.kappa = TRUE,kappa = kappa_vec[varnamei])
    if(covpar$cov.pars[2]==0)
    {covpar$cov.pars[2]=0.001}

    model <- Krig(x=train_dat[index_obsi,loc_variables], Y=train_dat[index_obsi,varnamei],theta=covpar$cov.pars[2],m=1,Covariance="Matern",smoothness=kappa_vec[varnamei])

    pred <- predict(model,as.matrix(test_dat[,loc_variables]))
    if(index_log[varnamei] ==0){pred <- exp(pred)}
    predmat[,varnamei] <- pred
  }

  return(predmat)
}



# For given kappa values, gives corresponding test error based on cross validation
kriging_old <- function(dat, loc_variables,
                        varname, varname_2nd, varname_ord, index_log,
                        method = "cressie", kappa_vec = rep(0.5,length(varname)),
                        K, seed = 2043, const = 0.64)
{
  set.seed(seed)
  # construct K-folds data
  cv_group <- cvFolds(n = nrow(dat),K = K,type = "random")

  nvar <- length(varname)
  miss_inf <- colSums(!is.na(dat[,varname, drop = FALSE]))   # number of nonmissing values for each varaible

  test_error0  <- rep(0,nvar) #output
  names(test_error0)  <- varname
  names(kappa_vec) <- varname

    for(k in 1:K)
    {
      test_dat  <- dat[cv_group$subsets[cv_group$which==k],]
      train_dat <- dplyr::setdiff(dat, test_dat)

      # prediction on test data
      predk <- fitpred_old(train_dat, test_dat, loc_variables = loc_variables, varname = varname,
                            varname_2nd = varname_2nd, varname_ord = varname_ord, index_log = index_log,
                            method = method, kappa_vec = kappa_vec, const = const)

      test_error0 <- test_error0 + colSums((predk - test_dat[,varname])^2,na.rm = TRUE)/miss_inf
    }

    test_error0 <- sqrt(test_error0) #rmse
  return(test_error0)
}

print(c("fitpred_old","kriging_old"))
