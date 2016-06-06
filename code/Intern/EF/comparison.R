##################### comparison #####################################
source("code/header1.R")
library(gstat)
library(sp)
library(corrplot)
library(dplyr)
library(cvTools)
library(RandomFields)
load("docs/index.RData")
dat <- read.csv("data/EF/sumnewdata.csv")
nr <- nrow(dat)


cord1.dec = SpatialPoints(cbind(dat$longitude, dat$latitude), 
                          proj4string=CRS("+proj=longlat"))
cord1.UTM <- spTransform(cord1.dec, CRS("+proj=utm +north +zone=14"))
dat$longitude <- coordinates(cord1.UTM)[,1]/1000
dat$latitude <- coordinates(cord1.UTM)[,2]/1000
dat<-as.data.frame(dat)

############## univariate ##################
source("code/cv_uni_kriging.R")
source("code/cv_kernel.R")
varname <- names(index_use)
varname_2nd <- names(index_2nd)
varname_ord <- names(index_ord)



kappav <- seq(0.2,2.5,by=0.3)
time1 <- Sys.time()
tab_uni <- cv_kriging(dat, loc_variables = c("longitude","latitude"),
           varname, varname_2nd, varname_ord, index_log,
           method = "cressie", kappav = seq(0.2,2.5,0.3),
           K = 10, seed = 2043)
time2 <- Sys.time()
write.csv(tab_uni,"docs/tabs/univariate/tab_uni.csv")

## based on kernel estimation####
toolbox_path = "Z:/XinWang/project/cokriging/code/kernel/Tools - validation toolbox"
time3 <- Sys.time()
tab_kernel <- cv_kernel(dat = dat, loc_variables = c("longitude","latitude"),
                        varname,index_log, lamv = 2^(-5:5), sigv = 2^(-5:5),
                        K=10, seed = 2043)

time4 <- Sys.time()
write.csv(tab_kernel,"docs/tabs/univariate/tab_kernel.csv")

time_uni <- c(difftime(time2,time1,units = "sec"), difftime(time4,time3,units = "sec"))
save(time_uni,file = "docs/tabs/univariate/time_uni.RData")

###################multivariate  ##############

#############S2, Tmax, Romeasured ######
source("code/cokriging_lmc_vg.R")
source("code/cokriging_RF.R")
source("code/cokriging_lmc_mle.R")
loc_variables <-  c("longitude","latitude")
param <- expand.grid(c(0.2,0.5,0.8,1,1.5,2),c(50,80,100,150,200))
numat <- expand.grid(seq(0.5,3,by=0.5),seq(0.5,3,by=0.5))
scale_vec <- c(50, 80, 100, 120, 150, 180, 200)

###Consider  Tmax, Romeasured and S2 Tmax

index_log_mat <- data.frame(v1 = c("Tmax","S2","S2"),v2 = c("Romeasured","Tmax","Romeasured"),log1 = c(1,0,0),log2 = c(1,1,1),stringsAsFactors=FALSE)


for(pp in 1:3)
{
  variables <- as.character(index_log_mat[pp,1:2])
  index_log <- as.numeric(index_log_mat[pp,3:4])
  tab_lmc_vg <- cv_lmc_vg(dat = dat,loc_variables = loc_variables,variables = variables,
                          param = param,index_log = index_log ,K = 10,seed = 2043)
  
  
  tab_mat <-  cv_RF(dat = dat,loc_variables,variables = variables,numat = numat,
                    index_log = index_log,K = 10, seed = 2043)
  
  
  #tab_lmc_mle <-  cv_lmc_mle(dat,loc_variables,variables = variables, ncomp = 1,
                           ## fix.nu = FALSE, fix.scale = TRUE,scale_vec, index_log = index_log,
                            #lower = c(0,0,-100,-100,0), upper = c(100,100,100,100,100),
                           # K = 10,seed = 2043)
  
  tab_bivariate <- cbind(tab_lmc_vg, tab_mat)
  colnames(tab_bivariate) <- c("LMC_vg","Matern")
  
  pathout <- paste("docs/tabs/bivariate/","tab_bivariate_",index_log_mat[pp,1],"_",index_log_mat[pp,2],".csv",sep="")
  write.csv(tab_bivariate,pathout)
}



##### LMC based likelihood ####

scale_vec <- c(50,80,100,120,150,180)
tab_lmc_mle <- data.frame(scale = rep(0,3),S2 =0 , Tmax = 0, Romeasured = 0)
for(pp in 1:3)
{
  variables <- as.character(index_log_mat[pp,1:2])
  index_log <- as.numeric(index_log_mat[pp,3:4])
  
  if(pp ==3)
  {
    scale_vec <- c(100,120,150)
  }
  
  tab_lmc_mlek <- matrix(0, length(scale_vec),2)
  for(j in 1:length(scale_vec))
  {
    tab_lmc_mlek[j,] <- cokriging_lmc_mle(dat,loc_variables = c("longitude","latitude"),
                                         variables  = variables, fix.nu = FALSE,
                                         scale1 = scale_vec[j], index_log = index_log, ncomp=1,
                                         lower = c(0,0,-10,-10,0),
                                         upper = c(100,100,100,100,100),
                                         K = 10,seed = 2043)
  }
  
  mean_vec <- colMeans(dat[,variables],na.rm = TRUE)
  index_min <- which.min(rowSums(sweep(tab_lmc_mlek, 2, mean_vec,FUN = "/"))) # relative min
  tab_lmc_mle[pp,variables] <- tab_lmc_mlek[index_min,]
  tab_lmc_mle[pp,"scale"] <- scale_vec[index_min]
}

write.csv(tab_lmc_mle,"docs/tabs/bivariate/tab_lmc_mle_S2_Tmax_Ro.csv",row.names = FALSE)

##########  time for 10-fold cross validation for fixed parameters #####

time_p1 <- Sys.time()
cokriging_lmc_vg(dat = dat,loc_variables = loc_variables,variables = c("Tmax","Romeasured"), 
                 kappa = 0.5,range = 100,index_log = c(1,1),K = 10,seed = 2043)

time_p2 <- Sys.time()

cokriging_lmc_mle(dat,loc_variables,variables = c("Tmax","Romeasured"), ncomp = 1,
                  fix.nu = FALSE, fix.scale = TRUE,scale1 = 100, index_log = c(1,1),
                  lower = c(0,0,-10,-10,0), upper = c(100,100,100,100,100),
                  K = 10,seed = 2043)
time_p3 <- Sys.time()
cokriging_RF(dat = dat,loc_variables,variables = c("Tmax","Romeasured"),nu1 = 0.5,nu2 =0.5,
             index_log = c(1,1),K = 10, seed = 2043)
time_p4 <- Sys.time()



time_biv <- diff(c(time_p1,time_p2,time_p3,time_p4))
names(time_biv) <- c("LMC_vg","LMC_mle","Matern")
save(time_biv,file = "docs/tabs/bivariate/time_biv.RData")


######## Consider other pairs with absolute correlation greater than 0.6 ####
df_cor <- read.csv("docs/tabs/df_cor.csv",stringsAsFactors = FALSE)
var_pairs <- subset(df_cor,abs(cor) > 0.6)$name[-c(1,2,3,5)]
names(index_log) <- varname
info_pairs <- data.frame(v1 = rep(0,length(var_pairs)),v2 = 0, log1 = 0, log2 = 0)
for(k in 1:length(var_pairs))
{
  pair_name <- strsplit(var_pairs[k],split ="_")[[1]]
  info_pairs[k,1:2] <- pair_name
  info_pairs[k,3:4] <- index_log[pair_name]
}


for(pp in 1:nrow(info_pairs))
{
  variables <- as.character(info_pairs[pp,1:2])
  index_log <- as.numeric(info_pairs[pp,3:4])
  tab_lmc_vg <- cv_lmc_vg(dat = dat,loc_variables = loc_variables,variables = variables,
                          param = param,index_log = index_log ,K = 10,seed = 2043)
  
  
  tab_mat <-  cv_RF(dat = dat,loc_variables,variables = variables,numat = numat,
                    index_log = index_log,K = 10, seed = 2043)
  

  tab_bivariate <- cbind(tab_lmc_vg, tab_mat)
  colnames(tab_bivariate) <- c("LMC_vg","Matern")
  
  pathout <- paste("docs/tabs/bivariate/","tab_bivariate_",info_pairs[pp,1],"_",info_pairs[pp,2],".csv",sep="")
  write.csv(tab_bivariate,pathout)
}
