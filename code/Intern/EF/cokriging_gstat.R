
### ------------------------cokriging based on fit.lmc in gstat package------- #####
#-------------------------------------------------------------------------------------
#dat is dataframe with locations and variables 2nd and 3rd are locations
# variables: the variables in cokriging, here just two variables
#kappa: value in matern model, default value is 0.5, which is exponential
#K number of folds
#seed with default value 2043
#index_log 0 means log transformation 1 means no transformation
#-------------------------------------------------------

library(gstat)
library(dplyr)

cokriging_gstat <- function(dat,variables, kappa =0.5, range = 100,
                            index_log = c(1,1),
                            K,seed = 2043)
{
  set.seed(seed)
  nv <- length(variables)
  cv_group <- cvFolds(n = nrow(dat),K = K,type = "random")
  test_error <-  rep(0,nv)#output
  num_obs <- colSums(!is.na(dat[,variables])) # number of observations for these two variables
  sub_dat0 <- select_(dat,.dots =c("UWI","longitude","latitude",variables)) # dat used
  

  
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
    
    indexk <- (!is.na(train_dat[,variables[1]]))&(!is.na(train_dat[,variables[2]])) # nonmissing index
    maxdist <-  max(dist(train_dat[indexk,2:3]))
    
    datg <- train_dat[indexk,] 
    coordinates(datg) <- c("longitude","latitude")
    coordinates(test_dat) <- c("longitude","latitude")
    
    g <- gstat(NULL, variables[1],as.formula(paste(variables[1],"~1")),datg)
    g <- gstat(g, variables[2],as.formula(paste(variables[2],"~1")),datg) # build the model for the two variables
    
    vario<-variogram(g,cressie = TRUE,cutoff = maxdist*0.64) #empirical variogram
    
    g.fit <-fit.lmc(vario,g,model=vgm(psill=1,model="Mat",range= range,nugget=0.5,kappa=kappa),correct.diagonal = 1.01) #fit variogram 
    
    # if nugget or pill is less than 0, set is to 0
    parm <- sapply(1:length(g.fit$model),function(x) g.fit$model[[x]][,2])
    nrp <- nrow(parm)
    if(sum(parm < 0)!=0)
    {
      index0 <- which(parm < 0)
      loc <- sapply(index0,function(x) {c(ceiling(index0/nrp),index0%%nrp)})
      for(m in 1:ncol(loc))
      {
        g.fit$model[[loc[1,m]]][loc[2,m],2] <- 0
      }
    }
    
    pred <- as.data.frame(predict(g.fit,test_dat)) #prediction
    test_dat <- as.data.frame(test_dat)
    
    pred_name <- paste(variables,"pred",sep=".")
    for(idi in 1:2)
    {
      if(index_log[idi] ==0){pred[,pred_name[idi]] <- exp(pred[,pred_name[idi]])}
    }
    
    test_error <- test_error + colSums((pred[,pred_name] - test_dat[,variables])^2,na.rm = TRUE)/num_obs
    

    
  } 
  test_error <- sqrt(test_error)
  return(test_error)
}




