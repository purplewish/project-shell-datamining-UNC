
### ------------------------cokriging based on fit.lmc in gstat package------- #####
#fitpred_lmc_vg fits model in train_dat predicts in test data
# cokriging_lmc_vg cross validation for given nu values
# tuning_lmc_vg gives the selected tuning parameter
# cv_lmc_vg gives the cross validation RMSE based on selected tuning parameter

#-------------------------------------------------------------------------------------
#dat is dataframe with locations and variables
# variables: the variables in cokriging, here just two variables
#loc_variables, variables for location
#kappa: value in matern model, default value is 0.5, which is exponential
#K number of folds, default is 10
#seed with default value 2043
#index_log 0 means log transformation 1 means no transformation
# const is used in variogram fitting
# cor.index: positive correlated 1, negative correlated -1
#-------------------------------------------------------

library(gstat)
library(dplyr)

fitpred_lmc_vg <- function(train_dat, test_dat,loc_variables = c("longitude","latitude"),
                           variables, kappa =0.5, range = 100,
                           index_log = c(1,1),const = 0.64,cor.index)
{
  nv <- length(variables)

  #log transformation for training data
  for(idi in 1:nv)
  {
    if(index_log[idi] ==0){train_dat[,variables[idi]] <- log(train_dat[,variables[idi]])
	                       train_dat[is.infinite(train_dat[,variables[idi]]),variables[idi]] <- NA  # set Inf to NA
	                       }
  }

  indexk <- rowSums(is.na(train_dat[,variables])) ==0  # nonmissing index
  maxdist <-  max(dist(train_dat[indexk,loc_variables]))
  # data construction for gstat
  datg <- train_dat[indexk,]
  sp::coordinates(datg) <- loc_variables
  sp::coordinates(test_dat) <- loc_variables

  if(nv == 2)
  {
    # model
    g <- gstat(NULL, variables[1],as.formula(paste(variables[1],"~1")),datg)
    g <- gstat(g, variables[2],as.formula(paste(variables[2],"~1")),datg) # build the model for the two variables
  }
  # if(nv == 3)
  # {
  #   # model
  #   g <- gstat(NULL, variables[1],as.formula(paste(variables[1],"~1")),datg)
  #   g <- gstat(g, variables[2],as.formula(paste(variables[2],"~1")),datg) # build the model for the two variables
  #   g <- gstat(g, variables[3],as.formula(paste(variables[3],"~1")),datg)
  # }

  vario<-variogram(g,cressie = TRUE,cutoff = maxdist*const) #empirical variogram

  g.fit <-fit.lmc(vario,g,model=vgm(psill=1,model="Mat",range= range,nugget=0.5,kappa=kappa),correct.diagonal = 1.01) #fit variogram

  # adjust parameters
  parm <- sapply(1:length(g.fit$model),function(x) g.fit$model[[x]][,2]) # transfer the estimates to a matrix
  # if the two variables are positive correlation, the results in the cross term should be nonnegative
  # otherwise, they should be nonpositive
  # change nugget or partial sill which are not reasonable to 0
  parm[,2] <- parm[,2]*cor.index # the second column is the cross term
  if(sum(parm < 0)!=0)
  {
    loc <- which(parm < 0,arr.ind = TRUE)
    for(m in 1:nrow(loc))
    {
      g.fit$model[[loc[m,2]]][loc[m,1],2] <- 0
    }
  }

  pred <- as.data.frame(predict(g.fit,test_dat)) #prediction
  pred_name <- paste(variables,"pred",sep=".")
  pred <- pred[,pred_name]
  # transform back
  for(idi in 1:nv)
  {
    if(index_log[idi] ==0){pred[,pred_name[idi]] <- exp(pred[,pred_name[idi]])}
  }
  return(pred)
}

cokriging_lmc_vg <- function(dat,loc_variables = c("longitude","latitude"),
                            variables, kappa =0.5, range = 100,
                            index_log = c(1,1),
                            K,seed = 2043, const = 0.64,cor.index)
{
  set.seed(seed)
  nv <- length(variables)
  cv_group <- cvFolds(n = nrow(dat),K = K,type = "random")
  test_error <-  rep(0,nv)#output
  num_obs <- colSums(!is.na(dat[,variables])) # number of observations for these two variables
  sub_dat0 <- select_(dat,.dots =c(loc_variables,variables)) # dat used

  #K-fold cross validation
  for(k in 1:K)
  {
    test_dat  <- sub_dat0[cv_group$subsets[cv_group$which==k],]
    train_dat <- sub_dat0[cv_group$subsets[cv_group$which!=k],]


    pred <- fitpred_lmc_vg(train_dat, test_dat,loc_variables = loc_variables,
                           variables = variables, kappa = kappa, range = range,
                           index_log = index_log,const = const,cor.index = cor.index) # prediction on test data

    test_error <- test_error + colSums((pred - test_dat[,variables])^2,na.rm = TRUE)/num_obs

  }
  test_error <- sqrt(test_error)
  return(test_error)
}

#parm[,1] for kappa, parm[,2] for range
tuning_lmc_vg <- function(dat,loc_variables = c("longitude","latitude"),
                   variables, param,
                   index_log = c(1,1),
                   K,seed = 2043,const = 0.64,cor.index)
{
  nv <- length(variables)
  res_lmc_vg <- matrix(0,nrow(param),nv)
  # for each combination of parameters, obtain test RMSE
  for(j in 1:nrow(param))
  {
    res_lmc_vg[j,] <- cokriging_lmc_vg(dat = dat,loc_variables = loc_variables,
                                   variables = variables,
                                   index_log = index_log,
                                   kappa = param[j,1],range = param[j,2],
                                   K = K,seed =seed, const = const,cor.index = cor.index)
  }

  mean_vec <- colMeans(dat[,variables],na.rm = TRUE)
  index_min <- which.min(rowSums(sweep(res_lmc_vg, 2, mean_vec,FUN = "/"))) # relative min: RMSE/mean
  param_select <- param[index_min,]
  return(param_select)
}

# kapppa and range are two tuning parameters
cv_lmc_vg <- function(dat, loc_variables, variables, param, index_log, K ,seed, const = 0.64, cor.index)
{
  set.seed(seed)
  nv <- length(variables)
  cv_group <- cvFolds(n = nrow(dat),K = K,type = "random")
  test_error_lmc_vg <-  rep(0,nv) #output
  num_obs <- colSums(!is.na(dat[,variables])) # number of observations for these two variables

  for(k in 1:K)
  {
    test_dat  <- dat[cv_group$subsets[cv_group$which==k],]
    train_dat <- dat[cv_group$subsets[cv_group$which!=k],]

    param_select <- tuning_lmc_vg(train_dat, loc_variables = loc_variables,
                                   variables = variables, param = param,
	               index_log = index_log, K = K ,seed = seed, const = const,cor.index=cor.index) # selected tuning parameters

    pred <- fitpred_lmc_vg(train_dat, test_dat, loc_variables = loc_variables, variables = variables,
                           kappa = as.numeric(param_select[1]),
                     	range = as.numeric(param_select[2]), index_log = index_log, const = const,cor.index = cor.index) # prediction on test data

    test_error_lmc_vg <- test_error_lmc_vg + colSums((test_dat[,variables] - pred)^2,na.rm = TRUE)/num_obs
  }
  test_error_lmc_vg <- sqrt(test_error_lmc_vg)
  return(test_error_lmc_vg)
}


print(c("fitpred_lmc_vg","cokriging_lmc_vg","tuning_lmc_vg","cv_lmc_vg"))
