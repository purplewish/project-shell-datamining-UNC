setwd("Z:/GitHup/project-shell-datamining-UNC/code/Intern/EF")
source('header.R')
source("loadData.R")
setwd(file.path(repo_path, "DataMiningUNC/Intern/EagleFord/data"))


#================================================================================================================================
# Universal Kriging Approach ###
#================================================================================================================================

##@Kriging using mean as aggregate method

##Aggregate data by mean
Sumnewdata3<-dplyr::summarise(newdata3,
                              S2=mean(S2,na.rm=TRUE),
                              Tmax=mean(Tmax,na.rm=TRUE),
                              XrdClayChlorite=mean(XrdClayChlorite,na.rm=TRUE),
                              Romeasured=mean(Romeasured,na.rm=TRUE),
                              GriWaterFilledPorosity=mean(GriWaterFilledPorosity,na.rm=TRUE),
                              XrdClaylllite=mean(XrdClaylllite,na.rm=TRUE),
                              GscCombustibleGasContent=mean(GscCombustibleGasContent,na.rm=TRUE),
                              S3=mean(S3,na.rm=TRUE),
                              GriSaturationSo=mean(GriSaturationSo,na.rm=TRUE),
                              XrdClayKaolinite=mean(XrdClayKaolinite,na.rm=TRUE),
                              Toc=mean(Toc,na.rm=TRUE),
                              S1=mean(S1,na.rm=TRUE),
                              GriSaturationSg=mean(GriSaturationSg,na.rm=TRUE),
                              NormalizedOil=mean(NormalizedOil,na.rm=TRUE),
                              GriGrainDensity=mean(GriGrainDensity,na.rm=TRUE),
                              XrdDolomite=mean(XrdDolomite,na.rm=TRUE),
                              CsgThoriumApi=mean(CsgThoriumApi,na.rm=TRUE),
                              XrdPlagioclase=mean(XrdPlagioclase,na.rm=TRUE),
                              StaticYoungsModulus=mean(StaticYoungsModulus,na.rm=TRUE),
                              GriTotalPorosity=mean(GriTotalPorosity,na.rm=TRUE),
                              GriGasFilledPorosity=mean(GriGasFilledPorosity,na.rm=TRUE),
                              GriBulkDensity=mean(GriBulkDensity,na.rm=TRUE),
                              GriTypeParameter=mean(GriTypeParameter,na.rm=TRUE),
                              XrdMarcasite=mean(XrdMarcasite,na.rm=TRUE),
                              GriMatrixPermeabilityAbsolute=mean(GriMatrixPermeabilityAbsolute,na.rm=TRUE))

Sumnewdata4<-dplyr::summarise(newdata4,
                              ConfiningStressDynamic=mean(ConfiningStressDynamic,na.rm=TRUE),
                              PoissonRatioDynamic=mean(PoissonRatioDynamic,na.rm=TRUE),
                              BulkDensityDynamic=mean(BulkDensityDynamic,na.rm=TRUE),
                              ShearVelocityDynamic=mean(ShearVelocityDynamic,na.rm=TRUE))


sumnewdata<-dplyr::full_join(Sumnewdata3,Sumnewdata4,by=c('UWI','latitude','longitude'))

#Change NaN to NA
for (i in 2:32)
{
  index<-is.na(sumnewdata[,i])
  sumnewdata[index,i]<-NA
}

###Delete Outliers###################
# manually remove two outliers to 
#reduce error (y_hat vs y)
sumnewdata$S2[60]<-NA
sumnewdata$NormalizedOil[77]<-NA
#####Check the trend#################



########Heat map##################
## NN interpolation
idw<-function(z,distance,k,num.neighs)
{
  idw.z<-rep(0,length(distance[,1]))
  for (i in 1:length(distance[,1]))
  {
    d<-sort(distance[i,],index.return=TRUE)
    w<-1/d$x[1:num.neighs]^k
    idw.z[i]<-sum(z[d$ix[1:num.neighs]]*w)/sum(w)
  }
  return(idw.z)
}  

# find grid point boundary
aq.ch<-chull(sumnewdata$longitude,sumnewdata$latitude)
aq.ch<-c(aq.ch,aq.ch[1])
aq.border<-cbind(sumnewdata$longitude[aq.ch],sumnewdata$latitude[aq.ch])

# generate grid point
aq.bbox<-sbox(as.points(sumnewdata$longitude,sumnewdata$latitude))
aq.grid<-gridpts(aq.bbox,npts=80000)
aq.grx<-unique(aq.grid[,1])
aq.gry<-unique(aq.grid[,2])
inside<-inout(aq.grid,aq.border,bound=TRUE)
aq.Grid<-aq.grid[inside,]

# Var TOC
hhh<-!is.na(sumnewdata$Toc)
distmat<-rdist(aq.Grid,cbind(sumnewdata$longitude[hhh],sumnewdata$latitude[hhh]))
m<-sumnewdata$Toc[hhh]
TToc<-idw(m,distmat,0.5,50) # NN interpolation. power 0~2, 50=# of neighbor
M<-cbind(aq.Grid,TToc)
M<-as.data.frame(M)
names(M)<-c('longitude','latitude','Toc')
jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
all_states <- map_data("county")
#county <- subset(all_states, (long>=(-100.3))&(long<=(-96))&(lat>27.6)&(lat<31) )
county <- subset(all_states, (long>=(-101))&(long<=(-95))&(lat>27.6)&(lat<31.5) )
ggplot() + geom_polygon( data=county, aes(x=long, y=lat, group = subregion),colour="grey", fill="white",size=1 )+
geom_tile(data = M, aes(x = longitude, y = latitude, fill = Toc),  alpha = 0.8)+scale_fill_gradientn(colours = jet.colors(7))#+
#stat_contour(data = M,aes( x = longitude, y = latitude, z = GriTypeParameter ))


ggplot()+geom_point(data=sumnewdata, aes(x=longitude, y=latitude)) # scatter plots
ggplot()+geom_point(data=sumnewdata, aes(x=longitude, y=latitude, colour =GriTypeParameter, size=GriTypeParameter))+
  scale_colour_gradientn(colours = jet.colors(7)) # scatter plot with size

# compared with krigging results (need to obtain krigging first)
surface(KrigCsgThoriumApi, type="C",xlab='X',ylab='Y',main='Kriging results for S3')

# calculate mse between true and krigging interpolation
hhh<-is.na(sumnewdata$CsgThoriumApi)
cbind(sumnewdata$CsgThoriumApi[!hhh],KrigCsgThoriumApi$fitted.values)
rmse(sumnewdata$CsgThoriumApi[!hhh],KrigCsgThoriumApi$fitted.values)/sd(sumnewdata$CsgThoriumApi[!hhh])


#Change UTM coordinate system
cord1.dec = SpatialPoints(cbind(sumnewdata$longitude, sumnewdata$latitude), proj4string=CRS("+proj=longlat"))
cord1.UTM <- spTransform(cord1.dec, CRS("+proj=utm +north +zone=14"))
sumnewdata$longitude <- coordinates(cord1.UTM)[,1]
sumnewdata$latitude <- coordinates(cord1.UTM)[,2]
sumnewdata<-as.data.frame(sumnewdata)


## use NN heatmap and krigged mse to judge which method to choose

##Kriging for 29 variables (Different trends for different variables)
for (i in 1:29)
{
varr<-names(sumnewdata)[i+3]  
goodname<-paste('Krig',varr,sep='') 
hhh<-!is.na(sumnewdata[,i+3])
lookb=variog(coords=sumnewdata[hhh,c(3,2)],data=sumnewdata[hhh,i+3],trend='2nd',max.dist=max(dist(sumnewdata[hhh,2:3]))*0.64) # 2nd means quadratic term
#lookc=variog(coords=sumnewdata[hhh,c(3,2)],data=sumnewdata[hhh,i+3],trend='2nd',op='cloud',max.dist=max(dist(sumnewdata[hhh,2:3]))*0.8)
#lookbc=variog(coords=sumnewdata[hhh,c(3,2)],data=sumnewdata[hhh,i+3],trend='2nd',bin.cloud=TRUE,estimator.type = "modulus",max.dist=max(dist(sumnewdata[hhh,2:3]))*0.8)
#plot(lookc, main='Variogram cloud plot of Tmax',xlab='distance',ylab='variogram')
#plot(lookbc, bin.cloud=TRUE, main="Binned variogram plot of Tmax",ylab='variogram',ylim=c(0,1000))  
#plot(lookb)  # variogram plot
covpar<-variofit(lookb)#,cov.model='matern',fix.kappa = FALSE) # use variog to fit exponential(or matern function) cov function parameters
if(covpar$cov.pars[2]==0) # range parameter if=0 means no dependence
{covpar$cov.pars[2]=0.001}
#if(covpar$kappa>2)
#{covpar$kappa=2}

# train krigged model
# exponential only need gamma input, nu=0.5 is default
# m-1 = trend
assign(goodname,Krig(x=sumnewdata[,c(3,2)], Y=sumnewdata[,i+3],theta=covpar$cov.pars[2],m=3))#,smoothness=covpar$kappa,Covariance="Matern"))
}

## Linear krigging
for (i in c(17,23))
{
  varr<-names(sumnewdata)[i+3]  
  goodname<-paste('Krig',varr,sep='') 
  hhh<-!is.na(sumnewdata[,i+3])
  lookb=variog(coords=sumnewdata[hhh,c(3,2)],data=sumnewdata[hhh,i+3],trend='1st',max.dist=max(dist(sumnewdata[hhh,2:3]))*0.64)
  #lookc=variog(coords=sumnewdata[hhh,c(3,2)],data=sumnewdata[hhh,i+3],trend='2nd',op='cloud',max.dist=max(dist(sumnewdata[hhh,2:3]))*0.8)
  #lookbc=variog(coords=sumnewdata[hhh,c(3,2)],data=sumnewdata[hhh,i+3],trend='2nd',bin.cloud=TRUE,estimator.type = "modulus",max.dist=max(dist(sumnewdata[hhh,2:3]))*0.8)
  #plot(lookc, main='Variogram cloud plot of Tmax',xlab='distance',ylab='variogram')
  #plot(lookbc, bin.cloud=TRUE, main="Binned variogram plot of Tmax",ylab='variogram',ylim=c(0,1000))  
  #plot(lookb)
  covpar<-variofit(lookb)#,cov.model='matern',fix.kappa = FALSE)
  if(covpar$cov.pars[2]==0) 
  {covpar$cov.pars[2]=0.001}
  #if(covpar$kappa>2)
  #{covpar$kappa=2}
  assign(goodname,Krig(x=sumnewdata[,c(3,2)], Y=sumnewdata[,i+3],theta=covpar$cov.pars[2],m=2))#,smoothness=covpar$kappa,Covariance="Matern"))}
}

# Ordinary
for (i in c(1,12,16,18))
{
  varr<-names(sumnewdata)[i+3]  
  goodname<-paste('Krig',varr,sep='') 
  hhh<-!is.na(sumnewdata[,i+3])
  lookb=variog(coords=sumnewdata[hhh,c(3,2)],data=sumnewdata[hhh,i+3],trend='cte',max.dist=max(dist(sumnewdata[hhh,2:3]))*0.64) # cte=constant
  #lookc=variog(coords=sumnewdata[hhh,c(3,2)],data=sumnewdata[hhh,i+3],trend='2nd',op='cloud',max.dist=max(dist(sumnewdata[hhh,2:3]))*0.8)
  #lookbc=variog(coords=sumnewdata[hhh,c(3,2)],data=sumnewdata[hhh,i+3],trend='2nd',bin.cloud=TRUE,estimator.type = "modulus",max.dist=max(dist(sumnewdata[hhh,2:3]))*0.8)
  #plot(lookc, main='Variogram cloud plot of Tmax',xlab='distance',ylab='variogram')
  #plot(lookbc, bin.cloud=TRUE, main="Binned variogram plot of Tmax",ylab='variogram',ylim=c(0,1000))  
  #plot(lookb)
  covpar<-variofit(lookb)#,cov.model='matern',fix.kappa = FALSE)
  if(covpar$cov.pars[2]==0)  
  {covpar$cov.pars[2]=0.001}
  #if(covpar$kappa>2)
  #{covpar$kappa=2}
  assign(goodname,Krig(x=sumnewdata[,c(3,2)], Y=sumnewdata[,i+3],theta=covpar$cov.pars[2],m=1))#,smoothness=covpar$kappa,Covariance="Matern"))
}




##Prediction for production well

newabc<-abc
for (i in 1:29)
{
  varr<-names(sumnewdata)[i+3]
  prename<-paste('Krig',varr,sep='') 
  newabc<-cbind(newabc,predict(get(prename),as.matrix(abc[,c(4,3)])))
  names(newabc)[i+5]<-varr
}


##Easier way to plot####
#set.panel()
#surface(KrigS2, type="C",xlab='X',ylab='Y',main='Kriging results for S3') # look at the surface 
#points(KrigS3$x)

## change coordinate back
newabc$Longitude <- coordinates(cord2.dec)[,1]
newabc$Latitude <- coordinates(cord2.dec)[,2]

write.csv(newabc,file='Interpolation for top 29 variables for production well(no truncation).csv')

## plot krigged heatmap
aq.ch<-chull(newabc$Longitude,newabc$Latitude)
aq.ch<-c(aq.ch,aq.ch[1])
aq.border<-cbind(newabc$Longitude[aq.ch],newabc$Latitude[aq.ch])

aq.bbox<-sbox(as.points(newabc$Longitude,newabc$Latitude))
aq.grid<-gridpts(aq.bbox,npts=45000)
inside<-inout(aq.grid,aq.border,bound=TRUE)
aq.Grid<-aq.grid[inside,]

cord3.dec = SpatialPoints(aq.Grid, proj4string=CRS("+proj=longlat"))
cord3.UTM <- spTransform(cord3.dec, CRS("+proj=utm +north +zone=14"))
Aq.Grid<-coordinates(cord3.UTM)

SS2<-predict(KrigS2,Aq.Grid)
M<-cbind(aq.Grid,SS2)
M<-as.data.frame(M)
names(M)<-c('longitude','latitude','S2')
M[,1:2]<-coordinates(cord3.dec)
jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))

# plot method 1
map<-get_map(location=c(left = -100.3, bottom = 27.7, right = -96.1, top = 31.2), zoom = 7, maptype='toner') 
p<-ggmap(map, extent='normal')
p+geom_tile(data = M, aes(x = longitude, y = latitude, z=S2, 
                          fill = S2),  alpha = 0.8)+scale_fill_gradientn(colours = jet.colors(7))+
  stat_contour(data = M,aes( x = longitude, y = latitude, z = S2 ))
  
# plot method 2
all_states <- map_data("county")
county <- subset(all_states, (long>=(-100.4))&(long<=(-96))&(lat>27.7)&(lat<31))
A<-unique(county$subregion)
county <- subset(all_states, (subregion %in% A) & (region=='texas'))
p <- ggplot() + geom_polygon( data=county, aes(x=long, y=lat, group = group),colour="grey", fill="white",size=1)+
  geom_tile(data = M, aes(x = longitude, y = latitude, 
                          fill = S2),  alpha = 0.8)+scale_fill_gradientn(colours = jet.colors(7))+
  stat_contour(data = M,aes( x = longitude, y = latitude, z = S2))


##############################################
#Loess for 29 variables (Just for comparison, we don't use it)
##############################################


#Loess for 29 variables


for (i in 1:29)
{
  varr<-names(sumnewdata)[i+3]  
  goodname<-paste('Loess',varr,sep='')  
  assign(goodname,loess(get(varr)~longitude*latitude, data=sumnewdata,degree=1, span=0.8, normalize=F,control=loess.control(surface = "direct"))) # degree=1 local linear; span: bigger more average, smaller put more weight on nearest ones
}

##Prediction for production well

Newabc<-abc
for (i in 1:29)
{
  varr<-names(sumnewdata)[i+3]
  prename<-paste('Loess',varr,sep='') 
  Newabc<-cbind(Newabc,predict(get(prename),as.matrix(abc[,4:3])))
  names(Newabc)[i+5]<-varr
}



Newabc$Longitude <- coordinates(cord2.dec)[,1]
Newabc$Latitude <- coordinates(cord2.dec)[,2]
write.csv(Newabc,file='Interpolation for top 29 variables for production well(no truncation).csv')



aq.ch<-chull(Newabc$Longitude,Newabc$Latitude)
aq.ch<-c(aq.ch,aq.ch[1])
aq.border<-cbind(Newabc$Longitude[aq.ch],Newabc$Latitude[aq.ch])

aq.bbox<-sbox(as.points(Newabc$Longitude,Newabc$Latitude))
aq.grid<-gridpts(aq.bbox,npts=45000)
inside<-inout(aq.grid,aq.border,bound=TRUE)

aq.Grid<-aq.grid[inside,]

cord3.dec = SpatialPoints(aq.Grid, proj4string=CRS("+proj=longlat"))
cord3.UTM <- spTransform(cord3.dec, CRS("+proj=utm +north +zone=14"))
Aq.Grid<-coordinates(cord3.UTM)

TGriTypeParameter<-predict(LoessGriTypeParameter,Aq.Grid)
M<-cbind(aq.Grid,TGriTypeParameter)
M<-as.data.frame(M)
names(M)<-c('longitude','latitude','GriTypeParameter')
M[,1:2]<-coordinates(cord3.dec)
jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))

# plot 1
map<-get_map(location=c(left = -100.3, bottom = 27.7, right = -96.1, top = 31.2), zoom = 7, maptype='toner') 
p<-ggmap(map, extent='normal')
p+geom_tile(data = M, aes(x = longitude, y = latitude, z=GriTypeParameter, 
                          fill = GriTypeParameter),  alpha = 0.8)+scale_fill_gradientn(colours = jet.colors(7))+
  stat_contour(data = M,aes( x = longitude, y = latitude, z = GriTypeParameter ))


# plot 2
all_states <- map_data("county")
#county <- subset(all_states, (long>=(-100.3))&(long<=(-96))&(lat>27.6)&(lat<31) )
county <- subset(all_states, (long>=(-101))&(long<=(-95))&(lat>27.6)&(lat<31.5) )
p <- ggplot() + geom_polygon( data=county, aes(x=long, y=lat, group = subregion),colour="grey", fill="white",size=1 )+
  geom_tile(data = M, aes(x = longitude, y = latitude, 
                         fill = GriTypeParameter),  alpha = 0.8)+scale_fill_gradientn(colours = jet.colors(7))+
  stat_contour(data = M,aes( x = longitude, y = latitude, z = GriTypeParameter))


p <- ggplot() + 
  geom_tile(data = M, aes(x = longitude, y = latitude, 
                          fill = GriTypeParameter),  alpha = 0.8)+scale_fill_gradientn(colours = jet.colors(7))
  

######Cross validation for Kriging and LOESS


#KrigingCV<- function(dat,d)
#{
#  hhh<-!is.na(dat[,d+3])
#  dat<-dat[hhh,]
#  n<-nrow(dat)
#  folds <- cvFolds(n, K=n)
#  mse <- NULL;  pred <- NULL; sol <- NULL;
#  for(i in 1:n){  
     #Split data into train/test set
#    test  <- dat[folds$subsets[folds$which==i],]
#    train <- dplyr::setdiff(dat, test)

#  lookb=variog(coords=train[,c(3,2)],data=train[,(d+3)],trend='2nd',max.dist=max(dist(train[,2:3]))*0.64)
#   covpar<-variofit(lookb)#,cov.model='matern',fix.kappa = FALSE)
#    if(covpar$cov.pars[2]==0) 
#    {covpar$cov.pars[2]=0.01}
#   #if(covpar$kappa>2)
   #{covpar$kappa=2}
#    model <- Krig(x=train[,c(3,2)],Y=train[,(d+3)],theta=covpar$cov.pars[2],m=3)#,smoothness=covpar$kappa,Covariance="Matern") 
#    test.pred <- cbind(test[,c(1,(d+3))], Pred=predict(model,as.matrix(test[,c(3,2)]))) 
#   # Uwi, Target, Pred, Latitude, Longitude
#    mse <- c(mse, sum((test.pred[,2]-test.pred[,3])^2)/nrow(test.pred))
#    pred <- rbind(pred, test.pred)  # save prediction results for fold i
#  }
  # CV results
 # sol <- data.frame(Name=names(train)[d+3],mse=mean(mse), rmse=sqrt(mean(mse)))
 # return(list(sol, pred))
#}

#Kriging <- KrigingCV(dat=sumnewdata,  d=3)



#LoessCV<- function(dat,d,span)
#{
#  hhh<-!is.na(dat[,d+3])
#  dat<-dat[hhh,]
#  n<-nrow(dat)
#  folds <- cvFolds(n, K=n)
#  varr<-names(sumnewdata)[d+3]  
#  mse <- NULL;  pred <- NULL; sol <- NULL;
#  for(i in 1:n){  
 #    Split data into train/test set
 #   test  <- dat[folds$subsets[folds$which==i],]
#    train <- dplyr::setdiff(dat, test)
#    model <- loess(get(varr)~longitude*latitude, data=train,degree=0, span=span, normalize=F,control=loess.control(surface = "direct"))
#    test.pred <- cbind(test[,c(1,(d+3))], Pred=predict(model,as.matrix(test[,c(3,2)]))) 
     
#    mse <- c(mse, sum((test.pred[,2]-test.pred[,3])^2)/nrow(test.pred))
#    pred <- rbind(pred, test.pred)  # save prediction results for fold i
#  }
   #CV results
#  sol <- data.frame(Name=names(train)[d+3],mse=mean(mse), rmse=sqrt(mean(mse)))
#  return(list(sol, pred))
#}

#Loess<-LoessCV(dat=sumnewdata,d=1, span=0.1)

#for (t in 1:100)
#{print(t*0.1)
#print(LoessCV(dat=sumnewdata,d=1, span=t*0.1)[[1]])
#}














#================================================================================================================================
# Tree, RandomForest and Boosting Algorithm ###(my data)
#================================================================================================================================

#Introducint Target variable into newabc data set

newabcY<-dplyr::inner_join(newabc,y2,by='Uwi')
newabcY<-dplyr::arrange(newabcY,Uwi)

#ggscatmat(newabcY,columns=c(3:35))

write.csv(newabcY,file='Data preparing for Machine learning.csv')

ggscatmat(newabcY,columns=c(15:16,35))

##boosting

runboostRegCV<- function(dat, no.tree, k)
{
  folds <- cvFolds(nrow(dat), K=k)
  mse <- NULL;  pred <- NULL; sol <- NULL;
  
  for(i in 1:k){  
    # Split data into train/test set
    
    test  <- dat[folds$subsets[folds$which==i],]
    train <- dplyr::setdiff(dat, test)
    #model <- gbm(Target~., data=train[,5:35], n.trees=no.tree, shrinkage=0.01,distribution='gaussian',interaction.depth=4)  
    model <- gbm(Target~., data=train[,3:35], n.trees=no.tree, shrinkage=0.01,distribution='gaussian',interaction.depth=10) 
    #####################################################################################################
    
    # Predict test dataset and calculate mse
    #test.pred <- cbind(test[,c(2,35)], Pred=predict(model,newdata=test[,5:34],n.trees<-no.tree), test[,c(36,37)])
    test.pred <- cbind(test[,c(2,35)], Pred=predict(model,newdata=test[,3:34],n.trees<-no.tree), test[,c(36,37)])# Uwi, Target, Pred, Latitude, Longitude
    mse <- c(mse, sum((test.pred[,2]-test.pred[,3])^2)/nrow(test.pred))
    pred <- rbind(pred, test.pred)  # save prediction results for fold i
  }
  # CV results
  sol <- data.frame(K=k,mse=mean(mse), rmse=mean(sqrt(mse)),n.Tree=no.tree)
  return(list(sol, pred))
}
#@@ 5-fold CV
boost5 <- runboostRegCV(dat=newabcY,  no.tree=5000, k=5)
predboost5<-boost5[[2]]

ggplot(predboost5,aes(x=Target,y=Pred))+geom_point(color='blue',size=0.2)+
  geom_abline(intercept = 0,size=1,colour='red')+xlab('True')+ylab('Prediction')+
  ggtitle('Oil production: Prediction vs True')+
  ylim(0, 70)+ylim(0, 70)+
theme(axis.title.x = element_text(size = 30, colour = 'black',vjust=-0.5),
      axis.title.y = element_text(size = 30, colour = 'black',vjust=1,hjust=0.55),
      axis.text  = element_text(size = 20, colour = 'black'),
     plot.title=element_text(size = 40, colour = 'black',vjust=2),
     panel.background = element_rect(fill = "aliceblue"))


plot(predboost5[,2],predboost5[,3],xlab='true',ylab='predict')
abline(a=0,b=1)


mmm<-rep(0,25)
for (i in 1:25)
{
  print(i)
  M<- runboostRegCV(dat=newabcY,  no.tree=8000, k=5)
  print(M[[1]])
  mmm[i]<-M[[1]][3]
}



fitControl <- trainControl(## 5-fold CV
  method = "cv",
  number = 5,
  savePredictions=TRUE)

gbmGrid <- expand.grid(interaction.depth=c(32),n.trees = c(500,800,1000,1200,1500), shrinkage=c(0.01), n.minobsinnode=10)


gbmFit <- train(Target ~ ., data = newabcY[,3:35],
                method = "gbm",
                 trControl = fitControl,
                 tuneGrid=gbmGrid,
                 verbose = FALSE)


mmm<-rep(0,25)
for (i in 1:25)
{
  
  gbmGrid <- expand.grid(interaction.depth=c(32),n.trees = c(1000), shrinkage=0.01, n.minobsinnode=10)
  print(i)
  M<- train(Target ~ ., data = newabcY[,3:35],
            method = "gbm",
            trControl = fitControl,
            tuneGrid=gbmGrid,
            verbose = FALSE)
  print(M$results)
  mmm[i]<-M$results[5]
}



#RMSE 5.233(0.051)       10,5000,0.01
#RMSE 5.180(0.038)       32,1000,0.01






###RandomForest
runRFRegCV <- function(dat, m, no.tree, k ,ntrace=500){
  
  folds <- cvFolds(nrow(dat), K=k)
  mse <- NULL;  pred <- NULL; sol <- NULL;
  
  for(i in 1:k){  
    # Split data into train/test set
    
    test  <- dat[folds$subsets[folds$which==i],]
    train <- dplyr::setdiff(dat, test)
    #model <- randomForest(Target~., data=train[,5:35], importance=T, mtry=m, do.trace=ntrace, ntree=no.tree)  
    model <- randomForest(Target~., data=train[,3:35], importance=T, mtry=m, do.trace=ntrace, ntree=no.tree)
    #####################################################################################################
    
    # Predict test dataset and calculate mse
    #test.pred <- cbind(test[,c(2,35)], Pred=predict(model,newdata=test[,5:34]), test[,c(36,37)])  # Uwi, Target, Pred, Latitude, Longitude
    test.pred <- cbind(test[,c(2,35)], Pred=predict(model,newdata=test[,3:34]), test[,c(36,37)])  # Uwi, Target, Pred, Latitude, Longitude
    
    mse <- c(mse, sum((test.pred[,2]-test.pred[,3])^2)/nrow(test.pred))
    pred <- rbind(pred, test.pred)  # save prediction results for fold i
  }
  # CV results
  m <- model$mtry  # get default value of mtry
  sol <- data.frame(K=k,mse=mean(mse), rmse=mean(sqrt(mse)), m=m, n.Tree=no.tree)
  return(list(sol, pred))
}

#@@ 5-fold CV 
set.seed(666)
rf <- runRFRegCV(dat=newabcY,  m=12, no.tree=1000, k=5)
predRF<- rf[[2]] 



ggplot(predRF,aes(x=Target,y=Pred))+geom_point(color='blue',size=0.2)+
  geom_abline(intercept = 0,size=1,colour='red')+xlab('True')+ylab('Prediction')+
  ggtitle('Oil production: Prediction vs True')+
  ylim(0, 70)+ylim(0, 70)+
  theme(axis.title.x = element_text(size = 30, colour = 'black',vjust=-0.5),
        axis.title.y = element_text(size = 30, colour = 'black',vjust=1,hjust=0.55),
        axis.text  = element_text(size = 20, colour = 'black'),
        plot.title=element_text(size = 40, colour = 'black',vjust=2),
        panel.background = element_rect(fill = "aliceblue"))



mmm<-rep(0,25)
for (i in 1:25)
{
  print(i)
  M<-runRFRegCV(dat=newabcY,  m=12, no.tree=1000, k=5)
  print(M[[1]])
  mmm[i]<-M[[1]][3]
}


#RMSE 5.157(0.046)


fitControl <- trainControl(## 5-fold CV
  method = "cv",
  number = 5,
  savePredictions=TRUE)

rfGrid <- expand.grid(mtry=c(12))


rfFit <- train(Target ~ ., data = newabcY[,3:35],
                 method = "rf",
                 trControl = fitControl,
                 tuneGrid=rfGrid,
                 ntree=1000,
                 verbose = FALSE)


mmm<-matrix(0,25,3)
mmm<-as.data.frame(mmm)
for (i in 1:25)
{ 
  print(i)
  L<- train(Target ~ ., data = newabcY[,3:35],
            method = "rf",
            trControl = fitControl,
            tuneGrid=rfGrid,
            ntree=1000,
            verbose = FALSE)
  print(L$results)
  mmm[i,1]<-L$results[2]

}


#RMSE 5.157(0.046)



#####Direct Kriging##################


runKriCV <- function(dat, k){
  
  folds <- cvFolds(nrow(dat), K=k)
  mse <- NULL;  pred <- NULL; sol <- NULL;
  
cord1.dec = SpatialPoints(cbind(dat$Longitude, dat$Latitude), proj4string=CRS("+proj=longlat"))
cord1.UTM <- spTransform(cord1.dec, CRS("+proj=utm +north +zone=14"))
dat$Longitude <- coordinates(cord1.UTM)[,1]
dat$Latitude <- coordinates(cord1.UTM)[,2]
  
  for(i in 1:k){  
    # Split data into train/test set
    
    test  <- dat[folds$subsets[folds$which==i],]
    train <- dplyr::setdiff(dat, test)
    
    
    #####################################################################################################
    
    # Predict test dataset and calculate mse
    
    lookb=variog(coords=train[,c(4,3)],data=train[,35],trend='cte')

    #lookbc=variog(coords=train[,c(4,3)],data=train[,35],trend='2nd',bin.cloud=TRUE,estimator.type = "modulus")
    #par(mfrow=c(2,2))
    #plot(lookb, main="binned variogram") 
    #plot(lookbc, bin.cloud=TRUE, main="clouds for binned variogram")  
    
    covpar<-variofit(lookb)#,cov.model='matern',fix.kappa = FALSE)
    if(covpar$cov.pars[2]==0) 
    {covpar$cov.pars[2]=0.01}
    #if(covpar$kappa>2)
    #{covpar$kappa=2}
    model <- Krig(x=train[,c(4,3)],Y=train[,35],theta=covpar$cov.pars[2],m=1)#,smoothness=covpar$kappa,Covariance="Matern") 
    test.pred <- cbind(test[,c(2,35)], Pred=predict(model,as.matrix(test[,c(4,3)])), test[,c(36,37)]) 
    
     # Uwi, Target, Pred, Latitude, Longitude
    mse <- c(mse, sum((test.pred[,2]-test.pred[,3])^2)/nrow(test.pred))
    pred <- rbind(pred, test.pred)  # save prediction results for fold i
  }
  # CV results
  sol <- data.frame(K=k,mse=mean(mse), rmse=mean(sqrt(mse)))
  return(list(sol, pred))
  
}



idw<-function(z,distance,k,num.neighs)
{
  idw.z<-rep(0,length(distance[,1]))
  for (i in 1:length(distance[,1]))
  {
    d<-sort(distance[i,],index.return=TRUE)
    w<-1/d$x[1:num.neighs]^k
    idw.z[i]<-sum(z[d$ix[1:num.neighs]]*w)/sum(w)
  }
  return(idw.z)
}  

aq.ch<-chull(newabcY$Longitude,newabcY$Latitude)
aq.ch<-c(aq.ch,aq.ch[1])
aq.border<-cbind(newabcY$Longitude[aq.ch],newabcY$Latitude[aq.ch])

aq.bbox<-sbox(as.points(newabcY$Longitude,newabcY$Latitude))
aq.grid<-gridpts(aq.bbox,npts=45000)
aq.grx<-unique(aq.grid[,1])
aq.gry<-unique(aq.grid[,2])
inside<-inout(aq.grid,aq.border,bound=TRUE)
aq.Grid<-aq.grid[inside,]
distmat<-rdist(aq.Grid,cbind(newabcY$Longitude,newabcY$Latitude))

TTarget<-idw(m,distmat,2,30)
M<-cbind(aq.Grid,TTarget)
M<-as.data.frame(M)
names(M)<-c('longitude','latitude','Target')
jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
all_states <- map_data("county")
county <- subset(all_states, (long>=(-100.3))&(long<=(-96))&(lat>27.6)&(lat<31) )
county <- subset(all_states, (long>=(-101))&(long<=(-95))&(lat>27.6)&(lat<31.5) )
p <- ggplot() + geom_polygon( data=county, aes(x=long, y=lat, group = subregion),colour="grey", fill="white",size=1 )+
  geom_tile(data = M, aes(x = longitude, y = latitude, fill = Target),  alpha = 0.8)+scale_fill_gradientn(colours = jet.colors(7))+
  stat_contour(data = M,aes( x = longitude, y = latitude, z = Target))


ggplot()+geom_point(data=newabcY, aes(x=Longitude, y=Latitude, colour = Target,size=Target))+scale_colour_gradientn(colours = jet.colors(7))


lookb=variog(coords=newabcY[,c(4,3)],data=newabcY[,35],trend='2nd',max.dist=max(dist(newabcY[,4:3]))*0.64)
plot(lookb, main="binned variogram") 
covpar<-variofit(lookb)
if(covpar$cov.pars[2]==0) 
{covpar$cov.pars[2]=0.01}
model <- Krig(x=newabcY[,c(4,3)],Y=newabcY[,35],theta=covpar$cov.pars[2],m=1)
set.panel()
surface(model, type="C",xlab='X',ylab='Y',main='Kriging results for Toc') # look at the surface 
cbind(newabcY$Target,model$fitted.values)
rmse(newabcY$Target,model$fitted.values)




predKri<- Kri[[2]] 

mmm<-rep(0,25)
for (i in 1:25)
{
  print(i)
  M<-runKriCV(dat=newabcY, k=5)
  print(M[[1]])
  mmm[i]<-M[[1]][3]
}


linear trend
#RMSE 5.278(0.039)

no trend
#RMSE 5.288(0.037)



###support vector regression########



runRegSVMCV <- function(dat, k, gamma, cost,nu){
  
  folds <- cvFolds(nrow(dat), K=k)
  mse <- NULL;  pred <- NULL; sol <- NULL;
  
  for(i in 1:k){  
    # Split data into train/test set
    
    test  <- dat[folds$subsets[folds$which==i],]
    train <- dplyr::setdiff(dat, test)
    model <- svm(Target~., train[,3:35],cost=cost,gamma=gamma,type='nu-regression',nu=nu)  
    
    #####################################################################################################
    
    # Predict test dataset and calculate mse
    test.pred <- cbind(test[,c(2,35)], Pred=predict(model,newdata=test[,3:34]), test[,c(36,37)])  # Uwi, Target, Pred, Latitude, Longitude
    mse <- c(mse, sum((test.pred[,2]-test.pred[,3])^2)/nrow(test.pred))
    pred <- rbind(pred, test.pred)  # save prediction results for fold i
  }
  # CV results
  sol <- data.frame(K=k,mse=mean(mse), rmse=mean(sqrt(mse)))
  return(list(sol, pred))
}
set.seed(897)
Svm <- runRegSVMCV(dat=newabcY, k=5, cost=20, gamma=0.2,nu=0.3)
predsvm<- Svm[[2]] 



mmm<-rep(0,50)

for (i in 1:50)
{
  Svm <- runRegSVMCV(dat=newabcY, k=5, cost=20, gamma=0.2,nu=0.3) 
  print(i)
  print(Svm[[1]])
  mmm[i]<-Svm[[1]][3]
}



A<-tune(svm, Target~., data=newabcY[,3:35], ranges = list(gamma =c(0.3), cost =10, nu=0.4,            
            type='nu-regression'),tunecontrol = tune.control(sampling = "cross",cross=5))
#print(i)
print(A$performance)

for (k in 2:34)
{
  print(k)
  for (i in 1:25)
  {
    
    A<-tune(svm, Target~., data=newabcY[,setdiff(3:35,k)], ranges = list(gamma =0.2, cost =20, nu=c(0.3),            
                           type='nu-regression'),tunecontrol = tune.control(sampling = "cross",cross=5))
    
    print(A$best.performance)
    M[k,i]<-as.vector(A$best.performance)
    
}
}



#RMSE 
#5.513(0.034)


########Neural network#########################
########(failure)


runRegneuralCV <- function(dat, k, hidden, epochs){
  
  folds <- cvFolds(nrow(dat), K=k)
  mse <- NULL;  pred <- NULL; sol <- NULL;
  localH2O <- h2o.init()
  
  for(i in 1:k){  
    # Split data into train/test set
    
    test  <- dat[folds$subsets[folds$which==i],]
    true  <- test[,35]
    train <- dplyr::setdiff(dat, test)
    train <- as.h2o(train[,3:35],conn=localH2O)
    test  <- as.h2o(test[,3:34],conn=localH2O)
    model <- h2o.deeplearning(x=1:32, y=33, training_frame=train,activation = "Tanh",# or 'Tanh'
                              #input_dropout_ratio = 0.2, # % of inputs dropout
                              #hidden_dropout_ratios = c(0.4), # % for nodes dropout
                              #balance_classes = TRUE, 
                              hidden = hidden,
                              reproducible=T,
                              epochs = epochs) # max. no. of epochs)  
    
    #####################################################################################################
    
    # Predict test dataset and calculate mse
    test.pred <- cbind(true, Pred=as.matrix(predict(model,newdata=test)))  # Uwi, Target, Pred, Latitude, Longitude

    mse <- c(mse, sum((test.pred[,1]-test.pred[,2])^2)/nrow(test.pred))
    pred <- rbind(pred, test.pred)  # save prediction results for fold i
  }
  # CV results
  sol <- data.frame(K=k,mse=mean(mse), rmse=sqrt(mean(mse)))
  return(list(sol, pred))
}

set.seed(897)
Neural <- runRegneuralCV(dat=newabcY, k=5, hidden=50, epochs=5000)
predneural<- Neural[[2]] 
predneural<-as.data.frame(predneural)
names(predneural)<-c('Target','Pred')


mmm<-rep(0,25)
mmmm<-rep(0,25)
for (i in 1:25)
{
  Neural <- runRegneuralCV(dat=newabcY, k=5, hidden=50, epochs=500)
  print(i)
  print(Neural[[1]])
  mmm[i]<-Neural[[1]][2]
  mmmm[i]<-Neural[[1]][3]
}


#MSE 43.69(0.51)
#RMSE 6.609(0.038)



#-------------------------------------------------------------------------------------------------------------------------
### Recover Curve
#-------------------------------------------------------------------------------------------------------------------------

qRecCurv <- function(x) {
  
  x <- as.data.frame(na.omit(x))
  
  n.row.x <- nrow(x)  
  n.col.x <- ncol(x)  
  
  ranks <- x %>% dplyr::mutate_each(funs(row_number)) %>% dplyr::arrange(desc(Target))  # ranks for each col and then ordered by 1st col(true value)
  
  rec.q <- data.frame(matrix(-1, nrow = n.row.x , ncol = n.col.x))  # recover quantiles
  rec.q[1,] <- (ranks[1,] == n.row.x)
  for (i in 2:n.row.x)
  {
    #rec.q[i,] <- ranks %>% slice(1:i) %>% summarise_each (funs(sum(.<=i)/i))
    rec.q[i,] <- ranks %>% dplyr::slice(1:i) %>% dplyr::summarise_each (funs(sum(.>=(n.row.x-i+1))/i))
  }
  names(rec.q)[1]<- "True"
  rec.q[,1]<-1:n.row.x/n.row.x
  
  #row.names(rec.q) <- sapply(100*(1:n.row.x)/n.row.x,  FUN = function(x) paste("P",round(x,digits = 0),sep = ""))
  
  return(rec.q)
}  


#@@ Comparison of different model
# Prediction of  models (30 vars)
pred.boost<-dplyr::select(predboost5,Uwi, Target,boost=Pred)
pred.RF<-dplyr::select(predRF, Uwi, RF=Pred)
pred.Svm<-dplyr::select(predsvm,Uwi,Svm=Pred)
pred.neural<-dplyr::select(predneural,neural=Pred,Target)


jo <- dplyr::left_join(pred.boost, pred.RF, by="Uwi")
jo <- dplyr::left_join(jo,pred.Svm,by='Uwi')
jo <- dplyr::left_join(jo,pred.neural,by='Target')
jo <- jo[,-1]  # rm Uwi

q.rec <- qRecCurv(jo) * 100

# Round to integer percentage
index <- ceiling(nrow(q.rec)*seq(0.3,100,0.3)/100)
q.rec <- q.rec[index, ]

q.rec1 <- q.rec %>% dplyr::select(True) %>% dplyr::mutate(RecRate=True, Method="Baseline")
q.rec2 <- q.rec %>% dplyr::select(True, X2) %>% dplyr::rename(RecRate=X2) %>% dplyr::mutate(Method="boost")
q.rec3 <- q.rec %>% dplyr::select(True, X3) %>% dplyr::rename(RecRate=X3) %>% dplyr::mutate(Method="RandomForest")
q.rec4 <- q.rec %>% dplyr::select(True, X4) %>% dplyr::rename(RecRate=X4) %>% dplyr::mutate(Method="SVM")
q.rec5 <- q.rec %>% dplyr::select(True, X5) %>% dplyr::rename(RecRate=X5) %>% dplyr::mutate(Method="Neural")

q.rec <- dplyr::union(q.rec1, q.rec2)
q.rec <- dplyr::union(q.rec, q.rec3)
q.rec <- dplyr::union(q.rec, q.rec4)
q.rec <- dplyr::union(q.rec, q.rec5)

ggplot(q.rec, aes(x=True, y=RecRate, colour=Method, group=Method)) + 
  geom_line(lwd=1.2) +
  scale_color_manual(values=c("#fe506e", "black", "#228b22","#0099cc",'brown')) +
  xlab("Top Quantile Percentage") + ylab("Recovery Rate") + 
  theme(#legend.position="none",
    axis.title.x = element_text(size=24),
    axis.title.y = element_text(size=24),
    axis.text.x = element_text(colour="grey20",size=15),
    axis.text.y = element_text(colour="grey20",size=15),
    legend.title=element_blank(),
    legend.text = element_text(size = 20),
    legend.justification=c(1,0), legend.position=c(1,0),
    legend.background = element_rect(fill="gray90", size=.5, linetype="dotted")
  )
# plot(q.rec, type="l", xlab="Top Quantile Percentage", ylab="Recover Rate")
# lines(q.rec[,1],q.rec[,1], col="red")







#================================================================================================================================
# Tree, RandomForest and Adaboosting Algorithm ###(Kaggle data)
#================================================================================================================================

#Introducing Target variable into newabc data set

newbbcY<-dplyr::inner_join(bbc,y2,by='Uwi')

#compare three methods

##boosting

runboostRegCV<- function(dat, no.tree, k)
{
  folds <- cvFolds(nrow(dat), K=k)
  mse <- NULL;  pred <- NULL; sol <- NULL;
  
  for(i in 1:k){  
    # Split data into train/test set
    
    test  <- dat[folds$subsets[folds$which==i],]
    train <- dplyr::setdiff(dat, test)
    model <- gbm(Target~., data=train[,2:32], n.trees=no.tree, shrinkage=0.01,distribution='gaussian',interaction.depth=4)  
    
    #####################################################################################################
    
    # Predict test dataset and calculate mse
    test.pred <- cbind(test[,c(1,32)], Pred=predict(model,newdata=test[,2:31],n.trees<-no.tree), test[,c(33,34)])  # Uwi, Target, Pred, Latitude, Longitude
    mse <- c(mse, sum((test.pred[,2]-test.pred[,3])^2)/nrow(test.pred))
    pred <- rbind(pred, test.pred)  # save prediction results for fold i
  }
  # CV results
  sol <- data.frame(K=k,mse=mean(mse), rmse=mean(sqrt(mse)),n.Tree=no.tree)
  return(list(sol, pred))
}
#@@ 5-fold CV
set.seed(666)
boost5 <- runboostRegCV(dat=newbbcY,  no.tree=500, k=5)


predboost5<-boost5[[2]]



mmm<-rep(0,25)
for (i in 1:25)
{
  print(i)
  M<- runboostRegCV(dat=newbbcY,  no.tree=5500, k=5)
  print(M[[1]])
  mmm[i]<-M[[1]][3]
}


#MSE 27.80(0.42)
#RMSE 5.273(0.040)



fitControl <- trainControl(## 5-fold CV
  method = "cv",
  number = 5,
  savePredictions=TRUE)

gbmGrid <- expand.grid(interaction.depth=c(5,12,30),n.trees = c(5000), shrinkage=c(0.01), n.minobsinnode=10)


gbmFit <- train(Target ~ ., data = newbbcY[,2:32],
                method = "gbm",
                trControl = fitControl,
                tuneGrid=gbmGrid,
                verbose = FALSE)





###RandomForest
runRFRegCV <- function(dat, m, no.tree, k ,ntrace=500){

  folds <- cvFolds(nrow(dat), K=k)
  mse <- NULL;  pred <- NULL; sol <- NULL;
  
  for(i in 1:k){  
    # Split data into train/test set

      test  <- dat[folds$subsets[folds$which==i],]
      train <- dplyr::setdiff(dat, test)
      model <- randomForest(Target~., data=train[,2:32], importance=T, mtry=m, do.trace=ntrace, ntree=no.tree)  
    
    #####################################################################################################
    
    # Predict test dataset and calculate mse
    test.pred <- cbind(test[,c(1,32)], Pred=predict(model,newdata=test[,2:31]), test[,c(33,34)])  # Uwi, Target, Pred, Latitude, Longitude
    mse <- c(mse, sum((test.pred[,2]-test.pred[,3])^2)/nrow(test.pred))
    pred <- rbind(pred, test.pred)  # save prediction results for fold i
  }
  # CV results
  m <- model$mtry  # get default value of mtry
  sol <- data.frame(K=k,mse=mean(mse), rmse=mean(sqrt(mse)), m=m, n.Tree=no.tree)
  return(list(sol, pred))
}

#@@ 5-fold CV 
set.seed(666)
rf <- runRFRegCV(dat=newbbcY,  m=12, no.tree=1000, k=5)
predRF<- rf[[2]] 


mmm<-rep(0,25)
for (i in 1:25)
{
  print(i)
  M<-runRFRegCV(dat=newbbcY,  m=12, no.tree=1000, k=5)
  print(M[[1]])
  mmm[i]<-M[[1]][3]
}


#MSE 26.53(0.520)
#RMSE 5.150(0.050)


#-------------------------------------------------------------------------------------------------------------------------
### Recover Curve
#-------------------------------------------------------------------------------------------------------------------------
  
  
  #@@ Comparison of different model
  # Prediction of  models (30 vars)
  


  
  pred.boost<-dplyr::select(predboost5,Uwi, Target,boost=Pred)
  pred.RF<-dplyr::select(predRF, Uwi, RF=Pred)
  pred.kaggle <- dplyr::select(newbbcY, Uwi, Rules.Prediction, Kaggle.Prediction)
  
  
  
  jo <- dplyr::left_join(pred.boost, pred.RF, by="Uwi")
  jo <- dplyr::left_join(jo, pred.kaggle, by="Uwi")
  jo <- jo[,-1]  # rm Uwi
  
  q.rec <- qRecCurv(jo) * 100
  
  # Round to integer percentage
  index <- ceiling(nrow(q.rec)*seq(0.3,100,0.3)/100)
  q.rec <- q.rec[index, ]
  
  q.rec1 <- q.rec %>% dplyr::select(True) %>% dplyr::mutate(RecRate=True, Method="Baseline")
  q.rec2 <- q.rec %>% dplyr::select(True, X2) %>% dplyr::rename(RecRate=X2) %>% dplyr::mutate(Method="boost")
  q.rec3 <- q.rec %>% dplyr::select(True, X3) %>% dplyr::rename(RecRate=X3) %>% dplyr::mutate(Method="RandomForest")
  q.rec4 <- q.rec %>% dplyr::select(True, X4) %>% dplyr::rename(RecRate=X4) %>% dplyr::mutate(Method="Rule Based")
  q.rec5 <- q.rec %>% dplyr::select(True, X5) %>% dplyr::rename(RecRate=X5) %>% dplyr::mutate(Method="Kaggle")
  
  q.rec <- dplyr::union(q.rec1, q.rec2)
  q.rec <- dplyr::union(q.rec, q.rec3)
  q.rec <- dplyr::union(q.rec, q.rec4)
  q.rec <- dplyr::union(q.rec, q.rec5)



  ggplot(q.rec, aes(x=True, y=RecRate, colour=Method, group=Method)) + 
    geom_line(lwd=1.2) +
    scale_color_manual(values=c("#fe506e", "black", "#228b22", "#0099cc", "#e95d3c")) +
    xlab("Top Quantile Percentage") + ylab("Recover Rate") + 
    theme(#legend.position="none",
      axis.title.x = element_text(size=24),
      axis.title.y = element_text(size=24),
      axis.text.x = element_text(colour="grey20",size=15),
      axis.text.y = element_text(colour="grey20",size=15),
      legend.title=element_blank(),
      legend.text = element_text(size = 20),
      legend.justification=c(1,0), legend.position=c(1,0),
      legend.background = element_rect(fill="gray90", size=.5, linetype="dotted")
    )
  # plot(q.rec, type="l", xlab="Top Quantile Percentage", ylab="Recover Rate")
  # lines(q.rec[,1],q.rec[,1], col="red")
  






