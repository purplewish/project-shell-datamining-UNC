######################################################################################################################
# conbien corewell ShaleGas data and corewell SCAL data
######################################################################################################################

#@@corewell ShaleGas
data3<-read.csv('data/EF/033_Core_ShaleGas.csv',header=TRUE,as.is=TRUE)
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
data4<-read.csv('data/EF/034_Core_SCAL.csv',header=TRUE,as.is=TRUE)
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

sumnewdata[is.na(sumnewdata)] <- NA

sumnewdata <- group_by(sumnewdata,latitude, longitude)
summarize(sumnewdata,count = n())

#### 72 75 same 
sumnewdata <- sumnewdata[-72,]
# sumnewdata <- sumnewdata[!duplicated(sumnewdata[,-(1:3)]),]

write.csv(sumnewdata,"data/EF/sumnewdata.csv",row.names = FALSE)
