
######################################################################################################################
# Collection of useful functions
######################################################################################################################


#---------------------------------------------------------------------------------------------------------------------

CalcAccumuProd <- function(date,prod,ndays){
  # Calculate the accumulate production of the first 180 producing date for each well
  #
  # Args:
  #   date: production date vector for a well
  #   prod: daily production vector for a well
  #   ndays: number of accumulation days
  #
  # Returns:
  #   Accumulate production of the first 180 producing date
  index <- (prod>0 & !is.na(prod))
  date <- date[index]
  prod <- prod[index]
  
  if(length(unique(date))<ndays){
    sol <- 0
  } else {
    index <- date<=(unique(date)[ndays])
    prod <- prod[index]
    sol <- sum(prod)
  }
  
  return(sol)
}

#---------------------------------------------------------------------------------------------------------------------

CalcMissPct <- function(x){
  # Calculate the missing percentage of a vector
  #
  # Args:
  #   x: a vector with missing value set as NA
  #
  # Returns:
  #   Missing percentage of the vector
  sol <- length(which(is.na(x)))/length(x)
  
  return(sol)
}

#---------------------------------------------------------------------------------------------------------------------

RmOutlierIQR <- function(x, a=1.5){
  # Remove outliers based on a*IQR, set outliers=NA 
  #
  # Args:
  #   x: a vector with missing value set as NA
  #
  # Returns:
  #   vector with a*IQR outlier removed
  if(is.factor(x)){
    return(x)
  
  } else {
    iqr <- IQR(x, na.rm=TRUE, type=7)
    upper <- quantile(x, na.rm=TRUE)[4] + a*iqr
    lower <- quantile(x, na.rm=TRUE)[2] - a*iqr
  
    index <- ( (x>upper) | (x<lower) ) & (!is.na(x))
    
    if(sum(index, na.rm=TRUE)>0){
      x[index] <- NA
    }
    
    return(x)
  }
  
}

#---------------------------------------------------------------------------------------------------------------------

ImputMissValue <- function(x, method=1){
  # Imputate missing value for a vector
  #   Categorical vector: use most frequent level to replace missing value
  #   Numerical vector: use mean(method=1)/median(method=2) to replace missing value
  #
  # Args:
  #   x: a vector with missing value set as NA
  #
  # Returns:
  #   Imputated vector
  if(is.factor(x)){  # catgorical
    
    index <- is.na(x) 
    x[index] <- names(which.max(table(x)))
    return(x)
    
  } else {  # numerical
    
    index <- is.na(x) 
    
    if(method==1){
      x[index] <- mean(x, na.rm=TRUE)
    } else {
      x[index] <- median(x, na.rm=TRUE)
    }
    return(x)
  
  }

}

#---------------------------------------------------------------------------------------------------------------------


