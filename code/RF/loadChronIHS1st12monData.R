######################################################################################################################
# Load Chronological IHS Data first 12 month
######################################################################################################################

#---------------------------------------------------------------------------------------------------------------------
### Data path
#---------------------------------------------------------------------------------------------------------------------
wd <- getwd()
setwd(file.path(repo_path, "Data/compFirst12Months"))


#---------------------------------------------------------------------------------------------------------------------
### Load data
#---------------------------------------------------------------------------------------------------------------------
#@@ Cov + Location + prod start date
cut.info <- read.csv("KrigedCoreDataOverview.csv", as.is=T)
keep <- c(2,9:49,51:53,58)

chro.dat <- NULL
for (i in 1:nrow(cut.info) ) {
  a <- read.csv(paste0("KrigedCoreData-", cut.info$cut.off.date[i],".csv"), as.is=T)
  a <- cbind(a[,keep], pct.IHS=cut.info$percentage.IHS.used[i], pct.Core=cut.info$precentage.CoreLabs.used[i])
  chro.dat <- rbind(chro.dat, a)
}

#@@ X and Y
x.vars <- names(a)[-c(1,2,44:48)]  # covars with location
x.vars.noloc <- names(a)[-c(1:4,44:48)] # covars w/o location
x.vars.geo <- x.vars.noloc[1:30]
loc <- names(a)[c(3,4)]

y <- names(a)[44]  # response


#@@ formula
formula.reg.chro <- formula(paste(paste(y,"~"), paste(x.vars,collapse="+")))
formula.reg.chro.loc <- formula(paste(paste(y,"~"), paste(loc,collapse="+")))
formula.reg.chro.geo <- formula(paste(paste(y,"~"), paste(x.vars.geo,collapse="+")))
formula.reg.chro.noloc <- formula(paste(paste(y,"~"), paste(x.vars.noloc,collapse="+")))



#---------------------------------------------------------------------------------------------------------------------
### Load grid data
#---------------------------------------------------------------------------------------------------------------------
setwd(file.path(repo_path, "Data/compFirst12MonthsCustomGrid"))

grid.dat <- NULL
for (i in 1:nrow(cut.info) ) {
  a <- read.csv(paste0("KrigedCoreDataGrid-", cut.info$cut.off.date[i],".csv"), as.is=T)
  a <- cbind(a[,x.vars.noloc], pct.IHS=cut.info$percentage.IHS.used[i], pct.Core=cut.info$precentage.CoreLabs.used[i])
  grid.dat <- rbind(grid.dat, a)
}



#---------------------------------------------------------------------------------------------------------------------
### Reset working dir
#---------------------------------------------------------------------------------------------------------------------
setwd(wd)


