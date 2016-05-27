######## cokriging #######
library(gstat)
library(sp)
library(corrplot)
library(dplyr)
library(cvTools)
library(RandomFields)
load("docs/index_use.RData")
dat <- read.csv("data/EF/sumnewdata.csv")
nr <- nrow(dat)

#Change UTM coordinate system
cord1.dec = SpatialPoints(cbind(dat$longitude, dat$latitude), 
                          proj4string=CRS("+proj=longlat"))
cord1.UTM <- spTransform(cord1.dec, CRS("+proj=utm +north +zone=14"))
dat$longitude <- coordinates(cord1.UTM)[,1]/1000
dat$latitude <- coordinates(cord1.UTM)[,2]/1000
dat<-as.data.frame(dat)

varname <- colnames(dat)[-(1:3)]
varname_use <- varname[index_use]

log_name <- c("S2", "XrdClayChlorite", "XrdClaylllite", "XrdClayKaolinite", "S1", "NormalizedOil", "XrdDolomite")
index_log <- 1-  names(index_use) %in% log_name

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


##### consider Tmax and Romeasured ########

source("code/cokriging_gstat.R")
source("code/cokriging_RF.R")
source("code/cokriging_lmc_mle.R")
#10-fold cross validation for lmc in gstat and bivariate matern in RandomFields
#same seed means they use the same split of training data and test data

parm <- expand.grid(c(0.2,0.5,0.8,1,1.5,2),c(50,80,100,150,200))

res_lmc <- matrix(0,nrow(parm),2)
for(j in 1:nrow(parm))
{
  res_lmc[j,] <- cokriging_gstat(dat = dat,variables = c("Tmax","Romeasured"),
                                 index_log = c(1,1),
                                 kappa = parm[j,1],range = parm[j,2],
                                 K = 10,seed = 2043)
}

res_lmc <-  cbind(parm,res_lmc)
colnames(res_lmc) <- c("smoothness","range","Tmax","Romeasured")
write.csv(res_lmc,"docs/tabs/res_lmc.csv",row.names = FALSE)

res_lmc[which.min(rowSums(res_lmc)),]

num <- expand.grid(seq(0.5,3,by=0.5),seq(0.5,3,by=0.5))
res_mat <- matrix(0,nrow(num),2)
for(j in 1:nrow(num))
{
  res_mat[j,] <- cokriging_RF(dat = dat, variables = c("Tmax","Romeasured"),
                            nu1 = num[j,1],nu2 = num[j,2], index_log = c(1,1),
                              K = 10,seed = 2043)
}

num[which.min(rowSums(res_mat)),]
res_mat[which.min(rowSums(res_mat)),]
res_mat <-  cbind(num,res_mat)
colnames(res_mat) <- c("smoothness1","smoothness2","Tmax","Romeasured")
write.csv(res_mat,"docs/tabs/res_mat.csv",row.names = FALSE)


scalem <- c(50,80,100,150,200)
res_lmc_mle <- matrix(0, length(scalem),2)
for(j in 1:length(scalem))
{
  res_lmc_mle[j,] <- cokriging_lmc_mle(dat,variables  = c("Tmax","Romeasured"), fix.nu = FALSE,
                                  scale1 = scalem[j], scale2 = scalem[j], index_log = c(1,1),
                                  lower = c(0,0,-10,-10,0,-10,-10,0),
                                  upper = c(100,100,100,100,100,100,100,100),
                                  K = 10,seed = 2043)
}


rownames(res_lmc_mle) <- scalem
write.csv(res_lmc_mle,"docs/tabs/res_lmc_mle.csv")



##### consider S2 and Tmax #### 
# S2 log transformation

parm <- expand.grid(c(0.2,0.5,0.8,1,1.5,2),c(50,80,100,150,200))

res_lmc2 <- matrix(0,nrow(parm),2)
for(j in 1:nrow(parm))
{
  res_lmc2[j,] <- cokriging_gstat(dat = dat,variables = c("S2","Tmax"),
                                 index_log = c(0,1),
                                 kappa = parm[j,1],range = parm[j,2],
                                 K = 10,seed = 2043)
}

parm[which.min(rowSums(res_lmc2)),]

res_lmc2[which.min(rowSums(res_lmc2)),]

res_lmc2 <-  cbind(parm,res_lmc2)
colnames(res_lmc2) <- c("smoothness","range","S2","Tmax")
write.csv(res_lmc2,"docs/tabs/res_lmc2.csv",row.names = FALSE)

num <- expand.grid(seq(0.5,3,by=0.5),seq(0.5,3,by=0.5))
res_mat2 <- matrix(0,nrow(num),2)
for(j in 1:nrow(num))
{
  res_mat2[j,] <- cokriging_RF(dat = dat, variables = c("S2","Tmax"),
                              nu1 = num[j,1],nu2 = num[j,2], index_log = c(0,1),
                              K = 10,seed = 2043)
}

num[which.min(rowSums(res_mat2)),]
res_mat2[which.min(rowSums(res_mat2)),]

res_mat2 <-  cbind(num,res_mat2)
colnames(res_mat2) <- c("smoothness1","smoothness2","S2","Tmax")

write.csv(res_mat2,"docs/tabs/res_mat2.csv",row.names = FALSE)


scalem <- c(80,100,150)
res_lmc_mle2 <- matrix(0, length(scalem),2)
for(j in 1:length(scalem))
{
  res_lmc_mle2[j,] <- cokriging_lmc_mle(dat,variables  = c("S2","Tmax"), fix.nu = FALSE,
                                       scale1 = scalem[j], scale2 = scalem[j], index_log = c(0,1),
                                       lower = c(0,0,-10,-10,0,-10,-10,0),
                                       upper = c(100,100,100,100,100,100,100,100),
                                       K = 10,seed = 2043)
}


rownames(res_lmc_mle2) <- scalem
write.csv(res_lmc_mle2,"docs/tabs/res_lmc_mle2.csv")

