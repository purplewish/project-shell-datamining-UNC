######################################################################################################################
# Load Data
######################################################################################################################

#---------------------------------------------------------------------------------------------------------------------
### Data path
#---------------------------------------------------------------------------------------------------------------------
wd <- getwd()
setwd(file.path(repo_path, "DataMiningUNC/Intern/EagleFord/data"))

#---------------------------------------------------------------------------------------------------------------------
### Load data
#---------------------------------------------------------------------------------------------------------------------
#@@corewell ShaleGas
data3<-read.csv('033_Core_ShaleGas.csv',header=TRUE,as.is=TRUE)
data3<-data3[-1,]



newdata3<-dplyr::select(data3,UWI=Unique.Private.Well.ID, latitude=Well.Latitude, longitude=Well.Longitude,
                 S2=Hydrocarbon...S2..mg.g., 
                 Tmax=Tmax..degrees.C., 
                 XrdClayChlorite=XRD.Clay.Fraction.Chlorite..weight.percent., 
                 Romeasured=Ro.Measured..percent.,
                 GriWaterFilledPorosity=GRI.Water.Filled.Porosity..percent.,
                 XrdClaylllite=XRD.Clay.Fraction.Illite..weight.percent.,
                 GscCombustibleGasContent=GSC.Combustible.Gas.Content,
                 S3=CO2...S3..mg.g.,
                 GriSaturationSo=GRI.Saturations...So..percent.Vp.,
                 XrdClayKaolinite=XRD.Clay.Fraction.Kaolinite..weight.percent.,
                 Toc=Leco.TOC..wt.percent.,
                 S1=Hydrocarbon...S1..mg.g.,
                 GriSaturationSg=GRI.Saturations...Sg..percent.,
                 NormalizedOil=Normalized.Oil.Content,
                 GriGrainDensity=GRI.Grain.Density..gm.cm.3.,
                 XrdDolomite=XRD.Bulk.Rock.Dolomite..weight.percent.,
                 CsgThoriumApi=CSG...Thorium..API.Units.,
                 XrdPlagioclase=XRD.Bulk.Rock.Plagioclase..weight.percent.,
                 StaticYoungsModulus=Static.Youngs.Modulus..10.6.psi.,
                 GriTotalPorosity=GRI.Total.Porosity..percent.,
                 GriGasFilledPorosity=GRI.Gas.Filled.Porosity..percent.,
                 GriBulkDensity=GRI.Bulk.Density..gm.cm.3.,
                 GriTypeParameter=GRI.Corey.Type.Parameter,
                 XrdMarcasite=XRD.Bulk.Rock.Marcasite..weight.percent.,
                 GriMatrixPermeabilityAbsolute=GRI.Matrix.Permeability...Absolute..md.
                 )
newdata3<-dplyr::arrange(newdata3,UWI)
newdata3[,2:28]<-sapply(newdata3[,2:28],FUN=as.numeric)



for (i in c(6,13,20))
{
  newdata3[,i][newdata3[,i]<0&!is.na(newdata3[,i])] <- NA # manually fix error, if <0 set to NA
}
newdata3<-dplyr::group_by(newdata3,UWI,latitude,longitude)


#@@corewell SCAL
data4<-read.csv('034_Core_SCAL.csv',header=TRUE,as.is=TRUE)
data4<-data4[-1,]
newdata4<-dplyr::select(data4,UWI=Unique.Private.Well.ID, latitude=Well.Latitude, longitude=Well.Longitude,
                 ConfiningStressDynamic=Confining.Stress...Dynamic,
                 PoissonRatioDynamic=Poisson.s.Ratio...Dynamic,
                 BulkDensityDynamic=Bulk.Density...Dynamic,
                 ShearVelocityDynamic=Shear.Velocity...Dynamic
)
newdata4<-dplyr::arrange(newdata4,UWI)
newdata4[,2:7]<-sapply(newdata4[,2:7],FUN=as.numeric)
newdata4<-dplyr::group_by(newdata4,UWI,latitude,longitude)









#@@Production well location + Prod start time + true depth
a <- read.csv("012_Prod_Well.csv", as.is=T)
a <- a %>% dplyr::select(Entity, API, Surface.Latitude, Surface.Longitude) %>% dplyr::filter(!is.na(API),!is.na(Surface.Latitude),!is.na(Surface.Longitude))  # Location

b <- read.csv("013_Prod_Header.csv", as.is=T)
b <- b %>% dplyr::select(Entity, Date.Production.Start) %>% dplyr::filter(!is.na(Date.Production.Start))  # Prod start date

ab <- dplyr::left_join(a, b, by="Entity")
ab <- ab %>% dplyr::distinct(API) %>% dplyr::rename(Uwi=API, Latitude=Surface.Latitude, Longitude=Surface.Longitude)
ab$Date.Production.Start <- as.Date(ab$Date.Production.Start, format="%Y-%m-%d")
prod.date.loc <- ab

c<-read.csv("020_Well_Header.csv",as.is=T)
c<-c %>% dplyr::distinct(UWI) %>% dplyr::select (Uwi=UWI,Depth.True.Vertical)

abc<- dplyr::inner_join(ab,c,by='Uwi')
abc<-abc[,-5]

setwd(file.path(repo_path, "DataMiningUNC/Kaggle/Final/RulesBasedApproach Oct 8/RulesBasedApproach Oct 8"))
y <- read.csv("Rules features using recent and Jan 2014 data.csv")

y1 <- dplyr::select(y,Uwi)  # 2631 x 1


abc <- dplyr::inner_join(abc,y1,by='Uwi')

# change to UTM coordinate for production well
cord2.dec = SpatialPoints(cbind(abc$Longitude, abc$Latitude), proj4string=CRS("+proj=longlat"))
cord2.UTM <- spTransform(cord2.dec, CRS("+proj=utm +north +zone=14"))
abc$Longitude <- coordinates(cord2.UTM)[,1]
abc$Latitude <- coordinates(cord2.UTM)[,2]




y2 <- dplyr::select(y,Uwi,Target,Rules.Prediction, Kaggle.Prediction) #2632*4




####Kaggle Kriging of the 29 variables 
setwd(file.path(repo_path, "DataMiningUNC/Kaggle/Final/Documentation and Input Files September 22 2014/Documentation and Input Files"))

x <- read.csv("EagleFordOilInput.csv")
x <- dplyr::distinct(x, Uwi)  # rm duplicate records (5222 x 35)
x <- dplyr::select(x, -Latitude, -Longitude, -Producer.EstimatedLength.Joined,-Core.RoCalculated.Kriged)
x.vars <- names(x)
x.vars <- x.vars[-1]  # rm Uwi (Uwi isn't a predictor)
bbc <- dplyr::inner_join(x, y1, by="Uwi")










#---------------------------------------------------------------------------------------------------------------------
### Reset working dir
#---------------------------------------------------------------------------------------------------------------------
setwd(wd)

