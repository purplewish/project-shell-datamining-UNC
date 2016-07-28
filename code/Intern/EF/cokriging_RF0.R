####################### cokriging Matern without tuning parameters #####################
#loc_variables, variables for location
#variables: the variables in cokriging, here just two variables
#K number of folds, default is 10
#seed with default value 2043
#index_log: 0 for log transformation, 1 for not
# scale.same: whether assume same scale parameter in the model
# fit_RF0: fit the bivariate matern model
# fitpred_RF0: fit and predict
# cv_RF0: cross validation 


library(RandomFields)
library(cvTools)
fit_RF0 <- function(dat, loc_variables, variables, scale.same = TRUE,index_log = c(1,1))
{
  nv <- length(variables)
  # parsimonious multivariate matern model with fixed kappa value
  nug <- RMmatrix(M = matrix(nc = 2, c(NA, 0, 0, NA)), RMnugget())

  if(scale.same)
  {
    pars.model <- nug + RMbiwm(nudiag = c(NA, NA), scale = NA,
                               cdiag = c(NA, NA), rhored = NA)
  }else{
    pars.model <- nug + RMbiwm(nudiag = c(NA, NA), nured = 1, s= rep(NA,3),
                               cdiag = c(NA, NA), rhored = NA)
  }

  for(idi in 1:nv)
  {
    if(index_log[idi] ==0){dat[,variables[idi]] <- log(dat[,variables[idi]])
    dat[is.infinite(dat[,variables[idi]]),variables[idi]] <- NA  # set Inf to NA
    }

  }

  sd0 <- apply(dat[,variables],2,function(x)(sd(x,na.rm = TRUE)))  # standard deviation
  mean0 <- colMeans(dat[,variables],na.rm = TRUE) # mean

  dat[,variables] <- scale(dat[,variables])   ## scale the target variables, since here assume zero mean


  indexk <- (!is.na(dat[,variables[1]]))&(!is.na(dat[,variables[2]])) # nonmissing index

  #fit parsimonious multivariate matern model with fixed kappa value
  resk <- RFfit(pars.model, x= dat[indexk,loc_variables[1]], y = dat[indexk,loc_variables[2]],
                data = dat[indexk,variables],split = FALSE)
  return(resk)
}





fitpred_RF0 <- function(train_dat, test_dat, loc_variables, variables,
                      scale.same = TRUE,index_log = c(1,1))
{
  nv <- length(variables)

  # parsimonious multivariate matern model with fixed kappa value
  nug <- RMmatrix(M = matrix(nc = 2, c(NA, 0, 0, NA)), RMnugget())

  if(scale.same)
  {
    pars.model <- nug + RMbiwm(nudiag = c(NA, NA), scale = NA,
                               cdiag = c(NA, NA), rhored = NA)
  }else{
    pars.model <- nug + RMbiwm(nudiag = c(NA, NA), nured = 1, s= rep(NA,3),
                               cdiag = c(NA, NA), rhored = NA)
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

  train_dat[,variables] <- scale(train_dat[,variables])   ## scale the target variables, since here assume zero mean


  indexk <- (!is.na(train_dat[,variables[1]]))&(!is.na(train_dat[,variables[2]])) # nonmissing index

  #fit parsimonious multivariate matern model with fixed kappa value
  resk <- RFfit(pars.model, x= train_dat[indexk,loc_variables[1]], y = train_dat[indexk,loc_variables[2]],
                data = train_dat[indexk,variables],split = FALSE)

  pred_mat <- RFinterpolate(resk, x= as.numeric(test_dat[,loc_variables[1]]),y = as.numeric(test_dat[,loc_variables[2]]),
                            data = train_dat[indexk,variables],
                            given = cbind(x = train_dat[indexk,loc_variables[1]],y=train_dat[indexk,loc_variables[2]])) # prediction

  # transform back
  pred <- pred_mat@data
  for(idi in 1:nv)
  {
    if(index_log[idi] ==0){pred[,idi] <- exp((pred[,idi]*sd0[idi] + mean0[idi]))}else{
      pred[,idi] <- (pred[,idi]*sd0[idi] + mean0[idi])}
  }

  return(pred)

}


cv_RF0<- function(dat, loc_variables, variables, scale.same, index_log, K ,seed)
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

    pred <- fitpred_RF0(train_dat, test_dat, loc_variables = loc_variables, variables = variables, scale.same = scale.same,index_log = index_log) # prediction on test data

    test_error_mat <- test_error_mat + colSums((test_dat[,variables] - pred)^2,na.rm = TRUE)/num_obs
  }
  test_error_mat <- sqrt(test_error_mat)
  return(test_error_mat)
}
