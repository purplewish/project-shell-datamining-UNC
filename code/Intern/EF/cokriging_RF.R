
### ------------------------cokriging based on matern cross covariance in RandomFields package------- #####
#fitpred_RF fits model in train_dat predicts in test data
# cokriging_RF: cross validation for given nu values 
# tuning_RF gives the selected tuning parameter
# cv_RF gives the test RMSE based on selected tuning parameter 
#-------------------------------------------------------------------------------------
#dat is dataframe with locations and variables 
##loc_variables, variables for location
# variables: the variables in cokriging, here just two variables
#nu: value in matern model, default value is 0.5, which is exponential
#K number of folds, default is 10
#seed with default value 2043
#index_log: 0 for log transformation, 1 for not
#-------------------------------------------------------

library(RandomFields)
library(dplyr)


fitpred_RF <- function(train_dat, test_dat, loc_variables, variables,
                       nu1 = 0.5, nu2 = 0.5, index_log = c(1,1))
{
  nv <- length(variables)

  # parsimonious multivariate matern model with fixed kappa value 
  nug <- RMmatrix(M = matrix(nc = 2, c(NA, 0, 0, NA)), RMnugget())
  pars.model <- nug + RMbiwm(nudiag = c(nu1, nu2), scale = NA,
                             cdiag = c(NA, NA), rhored = NA)
  
  #log transformation for training data
  for(idi in 1:nv)
  {
    if(index_log[idi] ==0){train_dat[,variables[idi]] <- log(train_dat[,variables[idi]])
	                      train_dat[is.infinite(train_dat[,variables[idi]]),variables[idi]] <- NA  # set Inf to NA
	                      }
 
  }
    
  sd0 <- apply(train_dat[,variables],2,function(x)(sd(x,na.rm = TRUE)))  # standard deviation 
  mean0 <- colMeans(train_dat[,variables],na.rm = TRUE) # mean 
  
  train_dat[,variables] <- scale(train_dat[,variables])   ## scale the target variables, since here assume zero mean
  
  
  indexk <- (!is.na(train_dat[,variables[1]]))&(!is.na(train_dat[,variables[2]])) # nonmissing index
  
  #fit parsimonious multivariate matern model with fixed kappa value 
  resk <- RFfit(pars.model, x= train_dat[indexk,loc_variables[1]], y = train_dat[indexk,loc_variables[2]], data = train_dat[indexk,variables],split = FALSE)
  
  pred_mat <- RFinterpolate(resk, x= as.numeric(test_dat[,loc_variables[1]]),y = as.numeric(test_dat[,loc_variables[2]]), data = train_dat[indexk,variables], given = cbind(x = train_dat[indexk,loc_variables[1]],y=train_dat[indexk,loc_variables[2]])) # prediction
  
  # transform back
  pred <- pred_mat@data
  for(idi in 1:nv)
  {
    if(index_log[idi] ==0){pred[,idi] <- exp((pred[,idi]*sd0[idi] + mean0[idi]))}else{
      pred[,idi] <- (pred[,idi]*sd0[idi] + mean0[idi])}
  }
  
  return(pred)
  
}


cokriging_RF <- function(dat,loc_variables,variables, nu1 =0.5,nu2=0.5,
                            index_log = c(1,1),
                            K,seed = 2043)
{

  set.seed(seed)
  nv <- length(variables)
  cv_group <- cvFolds(n = nrow(dat),K = K,type = "random")
  test_error_mat <-  rep(0,nv) #output
  num_obs <- colSums(!is.na(dat[,variables])) # number of observations for these two variables
  sub_dat0 <- select_(dat,.dots =c(loc_variables,variables)) # dat used
  
  #K-fold cross validation
  for(k in 1:K)
  {
    test_dat  <- sub_dat0[cv_group$subsets[cv_group$which==k],]
    train_dat <- sub_dat0[cv_group$subsets[cv_group$which!=k],]
    
    pred <- fitpred_RF(train_dat, test_dat, loc_variables = loc_variables, variables = variables,
                       nu1 = nu1, nu2 = nu2, index_log = index_log) # prediction on test data
  
    test_error_mat <- test_error_mat + colSums((test_dat[,variables] - pred)^2,na.rm = TRUE)/num_obs
  }
  test_error_mat <- sqrt(test_error_mat) #RMSE
  return(test_error_mat)
  
}


# tuning parameter selection
tuning_RF <- function(dat, loc_variables, variables, numat, index_log, K ,seed)
{
  res_mat <- matrix(0,nrow(numat),2)
  for(j in 1:nrow(numat))
  {
    res_mat[j,] <- cokriging_RF(dat = dat, loc_variables = loc_variables,
                                variables  = variables,
                                nu1 = numat[j,1],nu2 = numat[j,2], index_log = index_log,
                                K = K, seed = seed)
  }
  
  mean_vec <- colMeans(dat[,variables],na.rm = TRUE)
  
  index_min <- which.min(rowSums(sweep(res_mat, 2, mean_vec,FUN = "/"))) # relative min; RMSE/mean
  nu_select <- numat[index_min,]
  return(nu_select)
}


cv_RF <- function(dat, loc_variables, variables, numat, index_log, K ,seed)
{
  set.seed(seed)
  nv <- length(variables)
  cv_group <- cvFolds(n = nrow(dat),K = K,type = "random")
  test_error_mat <-  rep(0,nv) #output
  num_obs <- colSums(!is.na(dat[,variables])) 
  
  for(k in 1:K)
  {
    test_dat  <- dat[cv_group$subsets[cv_group$which==k],]
    train_dat <- dat[cv_group$subsets[cv_group$which!=k],]
    
    nu_select <- tuning_RF(train_dat, loc_variables = loc_variables, variables = variables, 
	                numat = numat, index_log = index_log, K = K ,seed = seed) # selected tuning parameter
    
    pred <- fitpred_RF(train_dat, test_dat, loc_variables = loc_variables, variables = variables, nu1 = as.numeric(nu_select[1]), nu2 = as.numeric(nu_select[2]), index_log = index_log) # prediction on test data
    
    test_error_mat <- test_error_mat + colSums((test_dat[,variables] - pred)^2,na.rm = TRUE)/num_obs
  }
  test_error_mat <- sqrt(test_error_mat)
  return(test_error_mat)
}

print(c("fitpred_RF","cokriging_RF","tunning_RF","cv_RF"))