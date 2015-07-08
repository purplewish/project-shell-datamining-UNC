######################################################################################################################
# Load Data
######################################################################################################################

#---------------------------------------------------------------------------------------------------------------------
### Data path
#---------------------------------------------------------------------------------------------------------------------
wd <- getwd()
setwd(file.path(repo_path, "/Data/ayata"))



#---------------------------------------------------------------------------------------------------------------------
### Load data
#---------------------------------------------------------------------------------------------------------------------
#@@ Load response variable (production data)
prod.2011 <- read.csv('04A OFM Daily Vols up to 2011 64937 Reformatted v1 05072015.csv')
prod.2012 <- read.csv('04B OFM Daily Vols 2012 59844 Reformatted v1 05072015.csv')
prod.2013 <- read.csv('04C OFM Daily Vols 2013 82297 Reformatted v2 05122015.csv')  # v1 has duplicate records
prod.2014 <- read.csv('04D OFM Daily Vols 2014 110862 Reformatted v1 05072015.csv')
prod.2015 <- read.csv('04E OFM Daily Vols 2015 43937 Reformatted v1 05092015.csv')

## Combine productions
prod <- rbind(prod.2011,prod.2012,prod.2013,prod.2014,prod.2015)[,c('well','Date','Condensate')]

## Group by Well and order by Date
prod$Date <- as.Date(prod$Date, "%m/%d/%Y")
prod <- prod %>% group_by(well) %>% arrange(Date)

## Accumulate production of first 180 producing date
accumu.prod <- prod %>% summarise(Condensate=CalcAccumuProd(Date, Condensate, 180)) %>% filter(Condensate>0)


#@@ Load predictors
x <- read.csv('01a Initial Calibration PMDB 515 Reformatted 05072015.csv',na.strings='')
x <- x %>% filter(Producing_Formation=='BONE SPRING')
names(x)[4]<-'well'

bs <- inner_join(accumu.prod,x,by='well')  # Bonespring data (Xs+Y)
#write.csv(BS,file='comb_data_without_cleaning.csv')



#---------------------------------------------------------------------------------------------------------------------
### Clean data
#---------------------------------------------------------------------------------------------------------------------
bs <-as.data.frame( bs[,c(1,7:8,9:16,18:29,32:41,46:90,94:98,2)] )

#@@ Remove unrepresented vars with large pct of missing value
miss.cutoff <- 0.1
miss.pct <- bs %>% summarise_each(funs(CalcMissPct))
bs <- bs[, which(miss.pct<miss.cutoff)]

#@@ Remove var with only one value
sd <- bs %>% summarise_each(funs(sd(.,na.rm=TRUE)))
bs <- bs[, which(sd>0)]
  
#@@ Remove outlier based on IQR for numeric var
bs <- bs %>% mutate_each(funs(RmOutlierIQR(.,1.5)), -well, -Condensate)
  
#@@ Imputate missing value
bs <- bs %>% mutate_each(funs(ImputMissValue(.,method=2)), -well, -Condensate)

#write.csv(bs,file='bs_clean.csv')
  
#---------------------------------------------------------------------------------------------------------------------
### Reset working dir
#---------------------------------------------------------------------------------------------------------------------
setwd(wd)



