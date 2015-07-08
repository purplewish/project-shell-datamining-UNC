


#---------------------------------------------------------------------------------------------------------------------
wd <- getwd()
#setwd(file.path(repo_path, "Data/compBigRulesKaggle"))
setwd(file.path(repo_path, "Data/compBigRulesKaggleMoreCov"))




#---------------------------------------------------------------------------------------------------------------------
### Load data
#---------------------------------------------------------------------------------------------------------------------
#@@ Cov + Location + prod start date
cut.info <- read.csv("KrigedCoreDataOverview.csv", as.is=T)

keep <- c(2,8,10:42,53:55,59)
i=1
a <- read.csv(paste0("KrigedCoreData-", cut.info$cut.off.date[i],".csv"), as.is=T)
a <- cbind(a[,keep], pct.IHS=cut.info$percentage.IHS.used[i], pct.Core=cut.info$precentage.CoreLabs.used[i])

# cut <- as.POSIXlt(as.Date(cut.info$cut.off.date[i]))
# cut$year <- cut$year-1
# cut <- as.Date(cut)

cut <- "2010-10-01"

d <- a[a$Date.Production.Start>cut, ]  # out of sample data (test data)

c <- d[,c(36,39)] 
c <- cbind(id=1:nrow(c),c)

n <- ceiling(nrow(c)*0.25) #1148

e <- c %>% arrange(desc(Norm.Lat.12.Month.Liquid))
true <- e[1:n,1]

e <- c %>% arrange(desc(Kriged.Production))
krig <- e[1:n,1]

overlap <- intersect(true, krig)
length(overlap)/n





  





