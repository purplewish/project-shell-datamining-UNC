
#######kriging this year ################
# nugget effect is from varigram estimation, constant mean for all variables 

#--------------------------------------------------------------------
##fitpred_new gives predicted values in test_dat based on train_dat
#kriging_new gives test error (RMSE) in cross validation based on given kappa value, which is used in tuning function to select tuning parameter
#---------------------------------------------------------------------

#loc_variables is the variables names for longitude and latitude
#varname considered variables
#index_log has the same length of varname, which indicates whether a log transformation is done: 0 means log transform, 1 means original scale
# method is the method in variog, default is "cressie" here
#kappa_vec has the same length with varname, smoothness parameter in Matern model
#const: maximum distance in variogram is const*maximum distance in data 


#K is the number of folds, default is K = 10 
#the default seed is 2403 
#--------------------------------------------------------------------------

library(geoR)

fitpred_new <- function(train_dat, test_dat, loc_variables, varname, index_log,
                        method = "cressie", kappa_vec = rep(0.5,length(varname)), const = 0.64)
{
  
  nvar <- length(varname) 
  predmat <- matrix(0,nrow(test_dat),nvar)
  colnames(predmat) <- varname
  
  
  for(i in 1:nvar)
  {
    varnamei <- varname[i]
    num0 <- sum(train_dat[,varnamei] == 0,na.rm = TRUE)
    # exclude the zero values for variables with log transformation in training data 
    if(index_log[i] == 0 & num0 >0)
    {
      index_obsi <- (!is.na(train_dat[,varnamei])) & ( train_dat[,varnamei] >0)
    } else{
      index_obsi <- !is.na(train_dat[,varnamei])
    }
    
    maxdist <-  max(dist(train_dat[index_obsi,loc_variables]))
 
    
    # fit variogram
    vgi <- variog(coords=train_dat[index_obsi,loc_variables],data=train_dat[index_obsi,varnamei],max.dist=maxdist*const, lambda = index_log[i])
    fit.vg <- variofit(vgi, cov.model = "matern", fix.nugget = FALSE,nugget = vgi$v[1],fix.kappa = TRUE,kappa = kappa_vec[i],weights = method) 
    
    # nugget should be non-negative
    if(fit.vg$nugget < 0)
    {
      fit.vg$nugget <- 0
    }
    res_liki <- likfit(coords=train_dat[index_obsi,loc_variables],data = train_dat[index_obsi,varnamei],fix.nugget = TRUE,nugget = fit.vg$nugget, fix.kappa=TRUE,kappa =  kappa_vec[i], lambda = index_log[i],ini.cov.pars = fit.vg$cov.pars,cov.model = "matern") #estimate parameters with fixed nugget effect from empirical variogram fitting, kappa =0.5 means exponential model
    
    kgi = krige.conv(coords=train_dat[index_obsi,loc_variables],data = train_dat[index_obsi,varnamei],loc = test_dat[,loc_variables], krige = krige.control(obj.m = res_liki,lambda=index_log[i])) #kriging based on estimates
    
    predmat[,i] <- kgi$predict

  }
  
  return(predmat)
}

 # For given kappa values, gives corresponding test error based on cross validation 
kriging_new <- function(dat, loc_variables,
                        varname, index_log,
                        method = "cressie", kappa_vec = rep(0.5,length(varname)),
                        K, seed = 2043, const = 0.64)
{
  set.seed(seed)
   # construct K-folds data 
  cv_group <- cvFolds(n = nrow(dat),K = K,type = "random")
  
  nvar <- length(varname)
  # number of nonmissing values for each varaible
  miss_inf <- colSums(!is.na(dat[,varname]))
  
  test_error1 <- rep(0,nvar) #output
  names(test_error1) <- varname
  
  for(k in 1:K)
  {
    test_dat  <- dat[cv_group$subsets[cv_group$which==k],]
    train_dat <- dplyr::setdiff(dat, test_dat)
  

    ## predicted results on test  data
    predk <- fitpred_new(train_dat, test_dat, loc_variables = loc_variables,
                         varname = varname, index_log = index_log,
                         method = method, kappa_vec = kappa_vec, const = const)
    
    test_error1 <- test_error1 + colSums((predk - test_dat[,varname])^2,na.rm = TRUE)/miss_inf
  }
  
  test_error1 <- sqrt(test_error1)
  
  return(test_error1)
}

print(c("fitpred_new","kriging_new"))