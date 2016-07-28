########## prediction based on the results of kernel estimation ###################
# x is the independent variable: locations 
# y is the dependent variable
#modelFit is the function
# fitpred_kernel: a function fit the model and do prediction based on selected tuning parameters
library(kernlab)
library(MASS)

fit_kernel  <- function(x, y, params = NULL){
    # read parameters
    kernel <- params$kernel
    lambda <- params$lambda

    # coerce data and target matrix
    x <- as.data.frame(x)
    y <- as.data.frame(y)

    # remove all non-numeric columns of x
    numericcols <- sapply(x,is.numeric)
    x <- x[,numericcols,drop=F]

    # prepare kernel matrix
    K <- kernelMatrix(kernel, x = as.matrix(x))

    # compute regression vector
    N <- nrow(K)
    beta <- ginv(as.matrix(K) + lambda*diag(N)) %*% as.matrix(y)

    return(list(xtrain = x, beta = beta, targetnames = names(y)))

}


predict_kernel <- function(modelFit, newdata, params = NULL){

    # read in parameters
    kernel <- params$kernel

    # coerce data and target matrix
    newdata <- as.data.frame(newdata)

    # remove all non-numeric columns of x
    numericcols <- sapply(newdata,is.numeric)
    newdata <- newdata[,numericcols,drop=F]

    # read in model
    beta <- modelFit$beta
    xtrain <- modelFit$xtrain
    xtrainnames <- names(xtrain)
    targetnames <- modelFit$targetnames

    # ensure that columns of y have right order
    xtrain = as.matrix(xtrain[,xtrainnames,drop=F])

    # compute ridge regression/Gaussian process regression prediction
    kappa <- kernelMatrix(kernel,x = as.matrix(newdata),y = xtrain)
    predictions <- kappa%*%beta

    predictions <- as.data.frame(predictions)
    names(predictions) <- targetnames
    row.names(predictions) <- row.names(newdata)

    return(predictions)

}

fitpred_kernel <- function(data, newdata, varname, loc_variables, K=10, seed = 2043)
{
  sd0 <- apply(data[,varname],2,function(x)(sd(x,na.rm = TRUE)))  # standard deviation
  mean0 <- colMeans(data[,varname],na.rm = TRUE) # mean
  nvar <- length(varname)

  dat0 <- data
  dat0[,varname] <- scale(dat0[,varname])

  #cross validation to select lambda and sigma
  foldSetup <- initFolds(foldType = "kfold", seed = seed, foldParams = list(k = K))
  dataFolds <- sampleFolds(nrow(data), foldSetup)

  # select best lambda and sigma
  tuningErrors <- tuning.computeErrors(dat0, dataFolds, predictor.gaussreg, gaussParamGrid ,
                                       target = varname, covariates = loc_variables,use_foreach = FALSE)

  paraselect <- tuning.getBestParams(tuningErrors, gaussParamGrid)

  pred_kernel <- matrix(0,nrow(newdata),nvar)
  colnames(pred_kernel) <- varname

  for(j in 1:nvar)
  {
    varj <- varname[j]
    indexj <- !is.na(data[,varj])
    fitmodel <- fit_kernel(x = dat0[indexj,loc_variables],y = dat0[indexj,varj],params = paraselect[[varj]])
    predj <- predict_kernel(modelFit = fitmodel,newdata = newdata[,loc_variables],params = paraselect[[varj]])
    predj <- predj*sd0[j] + mean0[j]
    pred_kernel[,varj] <- predj[,1]
  }
  return(pred_kernel)

}

print(c("fit_kernel","predict_kernel","fitpred_kernel"))
