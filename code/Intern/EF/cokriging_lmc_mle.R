
### ------------------------cokriging based on matern cross covariance in RandomFields package------- #####
#fitpred_lmc_mle fits model in train_dat predicts in test data
# cokriging_lmc_mle cross validation for given parameter values 
# tuning_lmc_mle gives the selected tuning parameter
# cv_lmc_mle gives the test RMSE based on selected tuning parameter

#-------------------------------------------------------------------------------------
#dat is dataframe with locations and variables 2nd and 3rd are locations
# variables: the variables in cokriging, here just two variables
#loc_variables: variables for locations 
#nu: value in matern model, default value is 0.5, which is exponential
#ncomp: number of components, it could be 1 and 2
#K number of folds
#seed with default value 2043
#scale is the scale(range parameter) in matern
#lower: the lower bound of parameters, which is needed to specify in RFfit function
#upper: the upper bound of parameters
#In this function, kappa or scale needs to be fixed one 
# the parameter order is nugget, nugget, coef,coef, nu(scale), coef,coef,nu(scale)
#-------------------------------------------------------

library(RandomFields)
library(dplyr)


fitpred_lmc_mle <- function(train_dat, test_dat,loc_variables,variables, ncomp = 2,
                            fix.nu = FALSE, fix.scale = !fix.nu,
                            nu1=0.5,nu2=0.5,scale1 = 100, scale2 = 100, index_log = c(1,1),
                            lower = NULL, upper = NULL)
{
  nv <- length(variables)
  
  # model
  nug <- RMmatrix(M = matrix(nc = 2, c(NA, 0, 0, NA)), RMnugget())
  #two components
  if(ncomp ==2 )
  {
    if(fix.nu)
    {
      pars.model  <- nug + RMmatrix(M = c(NA,NA), RMwhittle(nu = nu1, scale = NA))+
        RMmatrix(M = c(NA,NA), RMwhittle(nu = nu2, scale = NA))
    }
    
    if(fix.scale)
    {
      pars.model  <- nug + RMmatrix(M = c(NA,NA), RMwhittle(nu = NA, scale = scale1))+
        RMmatrix(M = c(NA,NA), RMwhittle(nu = NA, scale = scale2))
    }
    
    if((!fix.nu) & (!fix.scale))
    {
      pars.model  <- nug + RMmatrix(M = c(NA,NA), RMwhittle(nu = NA, scale = NA))+
        RMmatrix(M = c(NA,NA), RMwhittle(nu = NA, scale = NA))
    }
  }
  #one componet
  if(ncomp==1)
  {
    if(fix.nu)
    {
      pars.model  <- nug + RMmatrix(M = c(NA,NA), RMwhittle(nu = nu1, scale = NA))
    }
    
    if(fix.scale)
    {
      pars.model  <- nug + RMmatrix(M = c(NA,NA), RMwhittle(nu = NA, scale = scale1))
    }
    
    if((!fix.nu)&(!fix.scale))
    {
      pars.model  <- nug + RMmatrix(M = c(NA,NA), RMwhittle(nu = NA, scale = NA))
    }
  }
  
  
  #log transformation for training data
  for(idi in 1:nv)
  {
    if(index_log[idi] ==0){train_dat[,variables[idi]] <- log(train_dat[,variables[idi]])
	                      train_dat[is.infinite(train_dat[,variables[idi]]),variables[idi]] <- NA  # set Inf to NA
	                      }
  }
  
  sd0 <- apply(train_dat[,variables],2,function(x)(sd(x,na.rm = TRUE)))  # standard deviation 
  mean0 <- colMeans(train_dat[,variables],na.rm = TRUE) # mean 
  
  train_dat[,variables] <- scale(train_dat[,variables])   ## scale the target variables, since here assume
  
  indexk <- (!is.na(train_dat[,variables[1]]))&(!is.na(train_dat[,variables[2]])) # nonmissing index
  
  #fit lmc with given model structure
  resk <- RFfit(pars.model, x= train_dat[indexk,loc_variables[1]], y = train_dat[indexk,loc_variables[2]], 
                data = train_dat[indexk,variables],split = FALSE,lower = lower, upper = upper)
  
  pred_mat <- RFinterpolate(resk, x= as.numeric(test_dat[,loc_variables[1]]),y = as.numeric(test_dat[,loc_variables[2]]), data = train_dat[indexk,variables], given = cbind(x = train_dat[indexk,loc_variables[1]],y=train_dat[indexk,loc_variables[2]])) #prediction
  
  pred <- pred_mat@data
  for(idi in 1:nv)
  {
    if(index_log[idi] ==0){pred[,idi] <- exp((pred[,idi]*sd0[idi] + mean0[idi]))}else{
      pred[,idi] <- (pred[,idi]*sd0[idi] + mean0[idi])}
  }
  
  return(pred)
  
}


# cross validation for given tunning parameters
cokriging_lmc_mle <- function(dat,loc_variables,variables, ncomp = 2,
                              fix.nu = FALSE, fix.scale = !fix.nu,
                              nu1=0.5,nu2=0.5,scale1 = 100, scale2 = 100, index_log = c(1,1),
                             lower = NULL, upper = NULL,
                              K,seed = 2043)
{
  set.seed(seed)
  nv <- length(variables)
  cv_group <- cvFolds(n = nrow(dat),K = K,type = "random")
  test_error <-  rep(0,2)#output
  num_obs <- colSums(!is.na(dat[,variables])) # number of observations for these two variables
  sub_dat0 <- select_(dat,.dots =c("UWI",loc_variables,variables)) # dat used
  
  #K-fold cross validation
  for(k in 1:K)
  {
    test_dat  <- sub_dat0[cv_group$subsets[cv_group$which==k],]
    train_dat <- sub_dat0[cv_group$subsets[cv_group$which!=k],]
    
    pred <- fitpred_lmc_mle(train_dat, test_dat,loc_variables = loc_variables,variables = variables , 
                            ncomp = ncomp, fix.nu = fix.nu, fix.scale = fix.scale,
                            nu1= nu1,nu2=nu2,scale1 = scale1, scale2 = scale2, index_log = index_log,
                            lower = lower, upper = upper)
    test_error <- test_error + colSums((test_dat[,variables] - pred)^2,na.rm = TRUE)/num_obs
  }
  test_error <- sqrt(test_error) #RMSE
  return(test_error)
  
}



#### this is a function tunning for scale , which means nu parameter is estimated, default for 1 component
tuning_lmc_mle <- function(dat,loc_variables,variables, ncomp = 1,
                            fix.nu = FALSE, fix.scale = TRUE,scale_vec, index_log = c(1,1),
                            lower = NULL, upper = NULL,
                            K,seed = 2043)
{
  res_lmc_mle <- matrix(0, length(scale_vec),2)
  for(j in 1:length(scale_vec))
  {
    res_lmc_mle[j,] <- cokriging_lmc_mle(dat,loc_variables = loc_variables, ncomp = ncomp,
                                               variables  = variables, fix.nu = fix.nu,
                                                fix.scale = fix.scale,
                                               scale1 = scale_vec[j], scale2 = scale_vec[j],
                                               index_log = index_log,
                                               lower = lower, upper = upper, K = K,seed = seed)
  }
  
  mean_vec <- colMeans(dat[,variables],na.rm = TRUE)
  index_min <- which.min(rowSums(sweep(res_lmc_mle, 2, mean_vec,FUN = "/"))) # relative min
  scale_select <- scale_vec[index_min]
  return(scale_select)
  
}


cv_lmc_mle <- function(dat,loc_variables,variables, ncomp = 1,
                       fix.nu = FALSE, fix.scale = TRUE,scale_vec, index_log = c(1,1),
                       lower = c(0,0,-100,-100,0), upper = c(100,100,100,100,100),
                       K,seed = 2043)
{
  set.seed(seed)
  nv <- length(variables)
  cv_group <- cvFolds(n = nrow(dat),K = K,type = "random")
  test_error_lmc_mle <-  rep(0,nv) #output
  num_obs <- colSums(!is.na(dat[,variables])) 
  
  for(k in 1:K)
  {
    test_dat  <- dat[cv_group$subsets[cv_group$which==k],]
    train_dat <- dat[cv_group$subsets[cv_group$which!=k],]
    
    scale_select <- tuning_lmc_mle(train_dat, loc_variables = loc_variables, variables = variables, scale_vec = scale_vec, index_log = index_log, K = K ,seed = seed, lower = lower, upper = upper)
    
    scale_select <- as.numeric(scale_select)
    
    pred <- fitpred_lmc_mle(train_dat, test_dat,loc_variables,variables, ncomp = ncomp,
             fix.nu = fix.nu, fix.scale = fix.scale,scale1 = scale_select, scale2 = scale_select, 
             index_log = index_log,lower = lower, upper = upper)
    
    test_error_lmc_mle <- test_error_lmc_mle + colSums((test_dat[,variables] - pred)^2,na.rm = TRUE)/num_obs
  }
  test_error_lmc_mle <- sqrt(test_error_lmc_mle)
  return(test_error_lmc_mle)
}


print(c("fitpred_lmc_mle","cokriging_lmc_mle","tuning_lmc_mle","cv_lmc_mle"))