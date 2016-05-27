
### ------------------------cokriging based on matern cross covariance in RandomFields package------- #####
#-------------------------------------------------------------------------------------
#dat is dataframe with locations and variables 2nd and 3rd are locations
# variables: the variables in cokriging, here just two variables
#nu: value in matern model, default value is 0.5, which is exponential
#K number of folds
#seed with default value 2043
#-------------------------------------------------------

library(RandomFields)
library(dplyr)
cokriging_RF <- function(dat,variables, nu1 =0.5,nu2=0.5,
                            index_log = c(1,1),
                            K,seed = 2043)
{

  set.seed(seed)
  nv <- length(variables)
  cv_group <- cvFolds(n = nrow(dat),K = K,type = "random")
  test_error_mat <-  rep(0,nv) #output
  num_obs <- colSums(!is.na(dat[,variables])) # number of observations for these two variables
  sub_dat0 <- select_(dat,.dots =c("UWI","longitude","latitude",variables)) # dat used
  
  # parsimonious multivariate matern model with fixed kappa value 
  nug <- RMmatrix(M = matrix(nc = 2, c(NA, 0, 0, NA)), RMnugget())
  pars.model <- nug + RMbiwm(nudiag = c(nu1, nu2), scale = NA,
                             cdiag = c(NA, NA), rhored = NA)
  
  #K-fold cross validation
  for(k in 1:K)
  {
    test_dat  <- sub_dat0[cv_group$subsets[cv_group$which==k],]
    train_dat <- sub_dat0[cv_group$subsets[cv_group$which!=k],]
    
    #log transformation for training data
    for(idi in 1:nv)
    {
      if(index_log[idi] ==0){train_dat[,variables[idi]] <- log(train_dat[,variables[idi]])}
    }
    
    sd0 <- apply(train_dat[,4:5],2,function(x)(sd(x,na.rm = TRUE)))  # standard deviation 
    mean0 <- colMeans(train_dat[,4:5],na.rm = TRUE) # mean 
  
    train_dat[,4:5] <- scale(train_dat[,4:5])   ## scale the target variables, since here assume zero mean
  
    
    indexk <- (!is.na(train_dat[,variables[1]]))&(!is.na(train_dat[,variables[2]])) # nonmissing index
    
    #fit parsimonious multivariate matern model with fixed kappa value 
    resk <- RFfit(pars.model, x= train_dat[indexk,"longitude"], y = train_dat[indexk,"latitude"], 
                  data = train_dat[indexk,4:5],split = FALSE)
    
    pred_mat <- RFinterpolate(resk, x= as.numeric(test_dat[,"longitude"]),y = as.numeric(test_dat[,"latitude"]), data = train_dat[indexk,4:5], given = cbind(x = train_dat[indexk,"longitude"],y=train_dat[indexk,"latitude"])) # prediction
    
    pred <- pred_mat@data
    for(idi in 1:nv)
    {
      if(index_log[idi] ==0){pred[,idi] <- exp((pred[,idi]*sd0[idi] + mean0[idi]))}else{
        pred[,idi] <- (pred[,idi]*sd0[idi] + mean0[idi])}
    }
    
    
    test_error_mat <- test_error_mat + colSums((test_dat[,4:5] - pred)^2,na.rm = TRUE)/num_obs
  }
  test_error_mat <- sqrt(test_error_mat) #RMSE
  return(test_error_mat)
  
}



