##### variogram, histogram and correlation figure ####
#------------------- univariate  -------------------------------#
#######comparsion of kriging results between last year and this year #######
source("code/header1.R")
dat <- read.csv("data/EF/sumnewdata.csv")
nr <- nrow(dat)

#Change UTM coordinate system
cord1.dec = SpatialPoints(cbind(dat$longitude, dat$latitude), 
                          proj4string=CRS("+proj=longlat"))
cord1.UTM <- spTransform(cord1.dec, CRS("+proj=utm +north +zone=14"))
dat$longitude <- coordinates(cord1.UTM)[,1]/1000  # change the scale of location 
dat$latitude <- coordinates(cord1.UTM)[,2]/1000
dat<-as.data.frame(dat)

######################## summary of data ##################
varname <- colnames(dat)[-(1:3)]

# number of nonmissing values for each varaible
miss_inf <-  colSums(!is.na(dat[,-c(1:3)]))

# number of observations with zero value
inf0 <- apply(dat[,-c(1:3)],2,function(x){sum(x==0,na.rm = TRUE)}) 


pdf("docs/figures/histogram.pdf")
for(i in 1:29)
{
  hist(dat[,varname[i]],xlab="value",main = varname[i]) # histgram of 29 values
}
dev.off()


## index for variables with at least 80% observed values #
index_obs <- which(miss_inf >= round(nrow(dat)*0.8)) # 16 variables 

# histogram for these 16 variales 
pdf("docs/figures/histogram_sub.pdf")
for(i in index_obs)
{
  hist(dat[,varname[i]],xlab="value",main = varname[i])
}
dev.off()


#Do log transformation for some variables #
# S2 XrdClayChlorite, XrdClaylllite, XrdClayKaolinite, S1, NormalizedOil, XrdDolomite#
log_name <- c("S2", "XrdClayChlorite", "XrdClaylllite", "XrdClayKaolinite", "S1", "NormalizedOil", "XrdDolomite")

# histogram for transformated data 
pdf("docs/figures/histogram_log.pdf")
for(i in 1:length(log_name))
{
  hist(log(dat[,log_name[i]]),xlab="log_value",main = log_name[i])
}
dev.off()


index_log <- 1-  names(index_obs) %in% log_name # 1 for no transformation, 0 for log transformation

##### variogram for these 16 variables #####
# some variables are in log scale, there are two variables(XrdClayChlorite (2),XrdClayKaolinite (5)) which are in log scale but with zero values, these values are exculded. 
vg_ls <- list()
for(j in 1:length(index_obs))
{
  varnamej <- varname[index_obs[j]]
  if(index_log[j] == 0 & inf0[varnamej] >0)
  {
    index_obsj <- (!is.na(dat[,varnamej])) & (dat[,varnamej] >0)
  } else{
    index_obsj <- !is.na(dat[,varnamej])
  }
  maxdist <-  max(dist(dat[index_obsj,2:3]))
  vgj <- variog(coords=dat[index_obsj,c(3,2)],data=dat[index_obsj,varname[index_obs[j]]],max.dist=maxdist*0.64, lambda = index_log[j])
  
  vg_ls[[j]] <- vgj
}

pdf("docs/figures/variogram_sub.pdf")
for(j in 1:length(index_obs))
{
  plot(vg_ls[[j]],main = varname[index_obs[j]])
}
dev.off()


##### variogram based on original method (last year) #####
#  different mean structure for different varaibles ####
index_linear <- index_obs[index_obs %in% c(17,23)] # no variable 
index_ord <- index_obs[index_obs %in% c(1,12,16,18)] # 1, 12, 16, 18
index_2nd <- index_obs[!(index_obs %in% c(index_linear,index_ord))] # 14 variables

vg_ls0 <- list()
for(j in 1:length(index_obs))
{
  # 2nd order polynomial
  varnamej<- varname[index_obs[j]]  
  index_obsj <-!is.na(dat[,varnamej])
  d <-  max(dist(dat[index_obsj,2:3]))*0.64
  
  if(index_obs[j] %in% index_2nd)
  {
    vgj <-variog(coords=dat[index_obsj,c(3,2)],data=dat[index_obsj,varnamej],trend='2nd',max.dist=d) # 2nd means quadratic term
  }
  # Ordinary
  if(index_obs[j] %in% index_ord)
  {
    vgj <-variog(coords=dat[index_obsj,c(3,2)],data = dat[index_obsj,varnamej],trend='cte',max.dist=d)  # cte=constant
  }
  
  vg_ls0[[j]] <- vgj
  
}

pdf("docs/figures/variogram_sub0.pdf")
for(j in 1:length(index_obs))
{
  plot(vg_ls0[[j]],main = varname[index_obs[j]])
}
dev.off()

#figure which is easy to compare the variograms based on two methods
pdf("docs/figures/variogram_comapre.pdf",width = 10,height = 7)
par(mfrow = c(1,2),oma=c(0,0,2,0))
for(j in 1:length(index_obs))
{
  plot(vg_ls0[[j]],main = "original")
  plot(vg_ls[[j]],main = "new ")
  title(paste(varname[index_obs[j]],"(",index_log[j],")",sep=""),outer = TRUE)
}
dev.off()


#--------------------------------------------------------------------------------------------------

#######fit variogram to select variables which have nonzero range parameter ###
phi1 <- phi2 <- rep(0,length(index_obs))
names(phi1) <- names(phi2) <- varname[index_obs]
# 2nd order polynomial
for (i in index_2nd)
{
  varnamei<- varname[i]  
  index_obsi <-!is.na(dat[,varnamei])
  d <-  max(dist(dat[index_obsi,2:3]))*0.64
  vgi <-variog(coords = dat[index_obsi,c(3,2)],data = dat[index_obsi,varnamei],trend='2nd',max.dist=d) # 2nd means quadratic term
  covpari<-variofit(vgi,weights = "cressie",fix.nugget = FALSE,nugget = vgi$v[1])#,cov.model='matern',fix.kappa = FALSE) # use variog to fit
  phi1[varnamei] <- covpari$cov.pars[2]
}

# Ordinary
for (i in index_ord)
{
  varnamei<- varname[i]  
  index_obsi <-!is.na(dat[,varnamei])
  d <-  max(dist(dat[index_obsi,2:3]))*0.64
  vgi <-variog(coords=dat[index_obsi,c(3,2)],data=dat[index_obsi,varnamei],trend='cte',max.dist=d)  # cte=constant
  covpari<-variofit(vgi,weights = "cressie",fix.nugget = FALSE,nugget = vgi$v[1])#,cov.model='matern',fix.kappa = FALSE)
  phi1[varnamei] <- covpari$cov.pars[2]
}
#### end of method for last year #  

#------------------------------------------------------------------------------------------------------

####begining of method for this year ####  

for(i in 1:length(index_obs))
{
  varnamei <- varname[index_obs[i]]
  num0 <- sum(dat[,varnamei] == 0,na.rm = TRUE)
  # exclude the zero values for variables with log transformation in training data 
  if(index_log[i] == 0 & num0 >0)
  {
    index_obsi <- (!is.na(dat[,varnamei])) & (dat[,varnamei] >0)
  } else{
    index_obsi <- !is.na(dat[,varnamei])
  }
  
  maxdist <-  max(dist(dat[index_obsi,2:3]))
  
  vgi <- variog(coords=dat[index_obsi,c(3,2)],data=dat[index_obsi,varnamei],max.dist=maxdist*0.64, lambda = index_log[i])
  
  fit.vg <- variofit(vgi, fix.nugget =  FALSE, nugget = vgi$v[1],
                     weights = "cressie") #variogram fit based on exponential model in default
  phi2[varnamei] <- fit.vg$cov.pars[2]

}


#### variables considered 
index_use <- index_obs[(phi1 != 0) | (phi2!=0)] # one of the two is not zero phi, 12 variables 
index_ord <- index_use[index_use %in% c(1,12,16,18)] # no variales now
index_2nd <- index_use[!(index_use %in% c(index_linear,index_ord))] # 14 variables

log_name <- c("S2", "XrdClayChlorite", "XrdClaylllite", "XrdClayKaolinite", "S1", "NormalizedOil", "XrdDolomite")
index_log <- 1-  names(index_use) %in% log_name # 1 for no transformation, 0 for log transformation

save(index_use,index_ord,index_2nd, index_log,file = "docs/index.RData")



##############correlation #################

varname_use <- varname[index_use]
datc <- dat # datc contains log transformation for some variables
for(i in 1:length(index_use))
{
  if(index_log[i] ==0)
  {
    datc[,varname[index_use[i]]] <- log(datc[,varname[index_use[i]]])
  }
}

datc[abs(datc) == Inf] <- NA

### maps

pdf("docs/figures/maps.pdf")

for(i in 1:length(index_use))
{
  print(ggplot(datc,aes_string(x= "longitude",y="latitude",colour = varname[index_use[i]]))+geom_point()+scale_color_gradient2(low ="red", high ="blue",mid = "white", midpoint = median(datc[,varname[index_use[i]]],na.rm = TRUE))+theme_bw())
}
dev.off()

#correlation figure
cormat <- cor(datc[,varname[index_use]],use = "pairwise.complete.obs") # sample correlation matrix
png("docs/figures/corrplot12.png")
corrplot(cormat)
dev.off()


# sort correlation by the absolute value of sample correlation
name_comb <- unlist(lapply(1:(length(varname_use)-1),function(x){paste(varname_use[x],varname_use[(x+1):length(varname_use)],sep="_")})) # all combination of correlation
df_cor <- data.frame(name = name_comb, cor = cormat[lower.tri(cormat)]) # combine name and correlation
df_cor <- arrange(df_cor, desc(abs(cor)))
write.csv(df_cor,"docs/tabs/df_cor.csv", row.names = FALSE)