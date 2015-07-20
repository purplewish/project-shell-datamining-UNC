######################################################################################################################
# Data check
######################################################################################################################

#---------------------------------------------------------------------------------------------------------------------
### Setup
#---------------------------------------------------------------------------------------------------------------------
# code path
setwd("Z:/GitHup/project-shell-datamining-UNC/code/RF")
#source("header.R")


#---------------------------------------------------------------------------------------------------------------------
### Load data
#---------------------------------------------------------------------------------------------------------------------
#@@ Load Kaggle's data
source("loadData.R")  # all (dat)


#@@ Load production
wd <- getwd()
setwd(file.path(repo_path, "Data"))

## Normalized oil production = 3494 or 1678(Norm.1.yr.Num.of.Mos==12)
#  Kaggle may use monthly production data to calculate 1st 12 month production by themself
# a <- read.csv("011_Prod_IPCumNorm.csv", header=T, as.is=T)
# names(a) <- gsub("\\.+", ".", names(a))  # rm multiple dots in name
# a <- a[, c(1, 4, 7, 33, 36)]
# prod.norm <- a %>% filter(Norm.1.yr.Oil!='', Primary.Product=='O') %>% distinct(API)
#prod.norm <- a %>% filter(Norm.1.yr.Oil!='', Norm.1.yr.Num.of.Mos==12, Primary.Product=='O') %>% distinct(API)  

## Oil and horizontal well = 3596
a <- read.csv("012_Prod_Well.csv", header=T, as.is=T)
a <- a[,c(1,4,7,9,15,16)]
oil <- a %>% filter(Primary.Product=="O", Hole.Direction=='HORIZONTAL') %>% distinct(API)
#oil <- a %>% filter(Primary.Product=="O") %>% distinct(API)

## Oil well production date = 3924
a <- read.csv("013_Prod_Header.csv", header=T, as.is=T)
prod.date <- a %>% 
              select(Entity, Primary.Product, Date.Production.Start, Date.Production.Stop) %>% 
              filter(!is.na(Date.Production.Start), Primary.Product=='O')  # 3924 wells

## Oil well with production date = 3596 
oil.prod <- inner_join(oil, prod.date, by='Entity')
oil.prod <- oil.prod %>% 
                  mutate(Uwi=API, Primary.Product=Primary.Product.x) %>% 
                  select(-Primary.Product.x, -Primary.Product.y, -API)

raw.uwi <- select(oil.prod, Uwi)
new.uwi <- select(all, Uwi)

overlap1 <- intersect(raw.uwi, new.uwi)  # 1545
diff <- as.data.frame(setdiff(new.uwi, overlap1))  # 1086= 1085(Gas well) + 1(vertical Oil/gas well)


a <- read.csv("012_Prod_Well.csv", header=T, as.is=T)
a <- a[,c(1,4,7,9,15,16)]
oil.gas <- a %>% filter(Primary.Product=="G") %>% distinct(API) %>% mutate(Uwi=API)
diff.gas.prod <- inner_join(oil.gas, diff, by='Uwi')  # 1086
oil <- a %>% filter(Primary.Product=="O") %>% distinct(API) %>% mutate(Uwi=API) # 3628 wells
overlap.oil.prod <- inner_join(oil, overlap1, by='Uwi')  # 1545

wd <- getwd()
setwd(file.path(repo_path, "Code/RF/results"))
# 1085 gas well and 1 vertical oil well in Kaggle oil input
#write.csv(diff.gas.prod, "./gaswell_1085_vertoilwell_1_kaggle_oilinput.csv", row.names=F)
write.csv(diff.gas.prod, "./gaswell_1086_kaggle_oilinput.csv", row.names=F)
# 1545 horizontal oil well in kaggle's dataset (one well has two recored with primary product = O/G)
write.csv(overlap.oil.prod, "./hori_1545_oilwell_1545_kaggle_oilinput.csv", row.names=F)
setwd(wd)

  
# a <- read.csv("Master_Aug2013.csv", header=T, as.is=T)
# mas <- a[, c(1,7,11:13,22,24,27,30,36,37,51,52,70)]
# mas <- mas %>% rename(Uwi=UWI)
# mas.uwi <- select(mas, Uwi)
# overlap <- intersect(mas.uwi, new.uwi)

## Oil well with 1.yr.norm.prod = 3415
#oil.prod.date <- inner_join(oil.prod.date, prod, by='Entity')

  

  
  

