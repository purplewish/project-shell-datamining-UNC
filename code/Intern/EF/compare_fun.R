###################### comparisons #####################################
#kriging based on variogram estimation, mle estimation and nonparametric estimation ################-----------------------------------------------------------------

# method is the method in variog 
#K is the number of folds, if K = n, that is leave one out 
#the default seed is 2403 
#index_obs is the index of considered variables 
#index_2nd is the index for trend = "2nd" in last year
#index_ord is the index for trend ="cte" in last year
#index_log has the same length of index_obs, which indicates whether a log transfromation is done: 0 means log transform, 1 means original scale
# both TRUE, run both last year and this year, otherwise, run this year
#--------------------------------------------------------------------------
compare_fun <- function(dat, index_obs, index_2nd, index_ord, index_log,
                        method = "cressie", kappa = 0.5,
                        both = TRUE, K, seed = 2043)
{
  set.seed(seed)
  varname <- colnames(dat)[-(1:3)]
  
  # number of nonmissing values for each varaible
  miss_inf <- colSums(!is.na(dat[,-c(1:3)]))
  
  # construct K-folds data 
  cv_group <- cvFolds(n = nrow(dat),K = K,type = "random")
  
  #output of RMSE, 0 is for last year, 1 is for this year 
  test_error0 <- test_error1 <- rep(0,length(index_obs))
  names(test_error0) <- names(test_error1) <- names(index_obs)
  
  if(both)
  {
    for(k in 1:K)
    {
      test_dat  <- dat[cv_group$subsets[cv_group$which==k],]
      train_dat <- dplyr::setdiff(dat, test_dat)
      
      #---------------------------------------------------------------------------------------------------
      ##### beginning of the method in last year #####    
      # 2nd order polynomial
      for (i in index_2nd)
      {
        varnamei<- varname[i]  
        index_obsi <-!is.na(train_dat[,varnamei])
        d <-  max(dist(train_dat[index_obsi,2:3]))*0.64
        vgi <-variog(coords=train_dat[index_obsi,c(3,2)],data=train_dat[index_obsi,varnamei],trend='2nd',max.dist=d) # 2nd means quadratic term
        
        covpar<-variofit(vgi,fix.nugget = FALSE,nugget = vgi$v[1],weights = method,fix.kappa = TRUE,kappa = kappa)#,cov.model='matern',fix.kappa = FALSE) # use variog to fit exponential(or matern function) cov function parameters 
        # by default, the model is exponential 
        if(covpar$cov.pars[2]==0) # range parameter if=0 means no dependence
        {covpar$cov.pars[2]=0.001}
        model <- Krig(x=train_dat[,c(3,2)], Y=train_dat[,varnamei],theta=covpar$cov.pars[2],m=3,Covariance="Matern",smoothness=kappa) # range parameter is fixed, other parameters are estimated 
        
        test_error0[varnamei] <-  test_error0[varnamei] + sum((test_dat[,varnamei] - predict(model,as.matrix(test_dat[,c(3,2)])))^2,na.rm = TRUE)/ miss_inf[varnamei]
      }
      
      # Ordinary
      for (i in index_ord)
      {
        varnamei<- varname[i]  
        index_obsi <-!is.na(train_dat[,varnamei])
        d <-  max(dist(train_dat[index_obsi,2:3]))*0.64
        vgi <-variog(coords=train_dat[index_obsi,c(3,2)],data=train_dat[index_obsi,varnamei],trend='cte',max.dist=d)  # cte=constant
        covpar<-variofit(vgi,fix.nugget = FALSE,nugget = vgi$v[1],weights = method,fix.kappa = TRUE,kappa = kappa)#,cov.model='matern',fix.kappa = FALSE)
        if(covpar$cov.pars[2]==0)  
        {covpar$cov.pars[2]=0.001}
        
        model <- Krig(x=train_dat[,c(3,2)], Y=train_dat[,varnamei],theta=covpar$cov.pars[2],m=1,Covariance="Matern",smoothness=kappa)
        test_error0[varnamei] <-  test_error0[varnamei] + sum((test_dat[,varnamei] - predict(model,as.matrix(test_dat[,c(3,2)])))^2,na.rm = TRUE)/ miss_inf[varnamei]
      }
      #### end of method for last year #  
      
      #------------------------------------------------------------------------------------------------------
      
      ####begining of method for this year ####  
      
      for(i in 1:length(index_obs))
      {
        varnamei <- varname[index_obs[i]]
        num0 <- sum(train_dat[,varnamei] == 0,na.rm = TRUE)
        # exclude the zero values for variables with log transformation in training data 
        if(index_log[i] == 0 & num0 >0)
        {
          index_obsi <- (!is.na(train_dat[,varnamei])) & ( train_dat[,varnamei] >0)
        } else{
          index_obsi <- !is.na(train_dat[,varnamei])
        }
        
        maxdist <-  max(dist(train_dat[index_obsi,2:3]))
        
        vgi <- variog(coords=train_dat[index_obsi,c(3,2)],data=train_dat[index_obsi,varnamei],max.dist=maxdist*0.64, lambda = index_log[i])
        
        fit.vg <- variofit(vgi, fix.nugget = FALSE,nugget = vgi$v[1],fix.kappa = TRUE,kappa = kappa,
                           weights = method) 
        if(fit.vg$nugget < 0)
        {
          fit.vg$nugget <- 0
        }
        res_liki <- likfit(coords=train_dat[index_obsi,c(3,2)],data = train_dat[index_obsi,varnamei],fix.nugget = TRUE,nugget = fit.vg$nugget, fix.kappa=TRUE,kappa =  kappa, lambda = index_log[i],ini.cov.pars = fit.vg$cov.pars) #estimate parameters with fixed nugget effect from empirical variogram fitting, kappa =0.5 means exponential model
        
        kgi = krige.conv(coords=train_dat[index_obsi,c(3,2)],data = train_dat[index_obsi,varnamei],loc = test_dat[,c(3,2)], krige = krige.control(obj.m = res_liki,lambda=index_log[i])) #kriging based on estimates
        
        test_error1[varnamei] <-  test_error1[varnamei] + sum((test_dat[,varnamei] - kgi$predict)^2,na.rm = TRUE)/ miss_inf[varnamei]
      }
    }
    
    test_error0 <- sqrt(test_error0) #rmse
    test_error1 <- sqrt(test_error1)
    
    tab <- cbind(test_error0,test_error1)
    colnames(tab) <- c("original","new")
  }else{
    for(k in 1:K)
    {
      test_dat  <- dat[cv_group$subsets[cv_group$which==k],]
      train_dat <- dplyr::setdiff(dat, test_dat)
      ####begining of method for this year ####  
      
      for(i in 1:length(index_obs))
      {
        varnamei <- varname[index_obs[i]]
        num0 <- sum(train_dat[,varnamei] == 0,na.rm = TRUE)
        # exclude the zero values for variables with log transformation in training data 
        if(index_log[i] == 0 & num0 >0)
        {
          index_obsi <- (!is.na(train_dat[,varnamei])) & ( train_dat[,varnamei] >0)
        } else{
          index_obsi <- !is.na(train_dat[,varnamei])
        }
        
        maxdist <-  max(dist(train_dat[index_obsi,2:3]))
        
        vgi <- variog(coords=train_dat[index_obsi,c(3,2)],data=train_dat[index_obsi,varnamei],max.dist=maxdist*0.64, lambda = index_log[i],fix.kappa = TRUE,kappa = kappa)
        
        fit.vg <- variofit(vgi, fix.nugget = FALSE,nugget = vgi$v[1],
                           weights = method) #variogram fit based on exponential model in default
        if(fit.vg$nugget < 0)
        {
          fit.vg$nugget <- 0
        }
        res_liki <- likfit(coords=train_dat[index_obsi,c(3,2)],data = train_dat[index_obsi,varnamei],fix.nugget = TRUE,nugget = fit.vg$nugget, fix.kappa=TRUE,kappa =  kappa, lambda = index_log[i],ini.cov.pars = fit.vg$cov.pars) #estimate parameters with fixed nugget effect from empirical variogram fitting, kappa =0.5 means exponential model
        
        kgi = krige.conv(coords=train_dat[index_obsi,c(3,2)],data = train_dat[index_obsi,varnamei],loc = test_dat[,c(3,2)], krige = krige.control(obj.m = res_liki,lambda=index_log[i])) #kriging based on estimates
        
        test_error1[varnamei] <-  test_error1[varnamei] + sum((test_dat[,varnamei] - kgi$predict)^2,na.rm = TRUE)/ miss_inf[varnamei]
      }
    }
    
   tab <- sqrt(test_error1)
  }
  

  
  return(tab)
}