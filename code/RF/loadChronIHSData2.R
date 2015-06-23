######################################################################################################################
# Load Chronological IHS Data
######################################################################################################################

#---------------------------------------------------------------------------------------------------------------------
### Data path
#---------------------------------------------------------------------------------------------------------------------
wd <- getwd()
setwd(file.path(repo_path, "Data/compBigRulesKaggle"))



#---------------------------------------------------------------------------------------------------------------------
### Load data
#---------------------------------------------------------------------------------------------------------------------
#@@ Cov + Location + prod start date
cut.info <- read.csv("KrigedCoreDataOverview.csv", as.is=T)
keep <- c(2,8,10:42,48:50,54)

chro.dat <- NULL
for (i in 1:nrow(cut.info) ) {
  a <- read.csv(paste0("KrigedCoreData-", cut.info$cut.off.date[i],".csv"), as.is=T)
  a <- cbind(a[,keep], pct.IHS=cut.info$percentage.IHS.used[i], pct.Core=cut.info$precentage.CoreLabs.used[i])
  chro.dat <- rbind(chro.dat, a)
}


x.vars <- names(a)[-c(1,3,36:41)]  # covars with location
loc <- names(a)[c(4,5)]
x.vars.noloc <- names(a)[-c(1,3:5,36:41)] # covars w/o location

y <- names(a)[36]  # response



#@@ formula
formula.reg.chro <- formula(paste(paste(y,"~"), paste(x.vars,collapse="+")))
formula.reg.chro.loc <- formula(paste(paste(y,"~"), paste(loc,collapse="+")))
formula.reg.chro.noloc <- formula(paste(paste(y,"~"), paste(x.vars.noloc,collapse="+")))

kaggle.dat <- read.csv("KrigedCoreData-2013-08-05.csv", as.is=T)

#---------------------------------------------------------------------------------------------------------------------
### Reset working dir
#---------------------------------------------------------------------------------------------------------------------
setwd(wd)

