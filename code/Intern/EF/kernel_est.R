########## prediction based on the results of kernel estimation ###################
#modelFit is the function
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


print(c("fit_kernel","predict_kernel"))