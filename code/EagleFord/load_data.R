
options(scipen=999)  # remove scientific notation

# Set data directory -----------------------------------------------------------------------------------
path.old <- getwd()
path <- "C:/Apps/projects/DataMiningUNC/data/"  # data path
setwd(path)




# Load production data ---------------------------------------------------------------------------------
a <- read.csv("011_Prod_IPCumNorm.csv", header=T, as.is=T)
#length(unique(a$API))  # 5814
#b <- a[a$Norm.1.yr...Num.of.Mos==12,]  # 3012 wells
keep <- c(4, 7, 33, 34)
prod <- a[, keep]
names(prod)[3] <- "oil.12m"
names(prod)[4] <- "gas.12m"


# prod.gas <- prod[prod[,1]=='G', c(1,2,4)] # 2373
# prod.gas <- prod.gas[prod.gas[,3]>0,]     # 2366 
# length(unique(prod.gas[,2]))              # 2223 unique

prod.oil <- prod[prod[,1]=='O', 2:3]        # 3924
prod.oil.uniq <- prod.oil[!duplicated(prod.oil[,1]),]  #3670 unique ID
prod.oil <- prod.oil.uniq[complete.cases(prod.oil.uniq),]  # 3475 remove NA
write.table(prod.oil, file="./analysis/EagleFordRawOilUniq.csv", row.names=F, sep=",")



# Load well data ---------------------------------------------------------------------------------------
# TestTreatment
a <- read.csv("021_Well_TestTreatment.csv", header=T, as.is=T)
keep <- c(1, 7:14, 16:19)
well.testtrt <- a[, keep]


# Treatment summary
a <- read.csv("022_Well_TreatmentSum.csv", header=T, as.is=T)
keep <- c(1, 3, 5)
well.trtsum <- a[, keep]
names(well.trtsum)[3] <- "proppant.sand.lbs"


# Test perforation
a <- read.csv("023_Well_TestPerf.csv", header=T, as.is=T)
keep <- c(1, 6:8)
well.testperf <- a[, keep]
names(well.testperf)[4] <- "lat.length"


# Mud
a <- read.csv("024_Well_Mud.csv", header=T, as.is=T)
keep <- c(1, 3, 5)
well.mud <- a[, keep]

# Survey point
a <- read.csv("025_Well_SurveyPt.csv", header=T, as.is=T)
keep <- c(1, 4:11)
well.surveypt <- a[, keep]



# Load core data ---------------------------------------------------------------------------------------
# GeoChem
col.class <- read.csv("031_Core_GeoChem.csv", header=F, skip=1, nrow=1, as.is=T)
unique(unlist(col.class, use.names=F))  # check unique class type
col.class[col.class %in% c("INTEGER", "DECIMAL")] <- "numeric"
col.class[col.class == "STRING"] <- "factor"

header <- read.csv("031_Core_GeoChem.csv", header=T, as.is=T)
header <- colnames(header)

core.geochem <- read.csv("031_Core_GeoChem.csv", header=F, skip=2, colClasses=col.class)
colnames(core.geochem) <- header


# RCA
col.class <- read.csv("032_Core_RCA.csv", header=F, skip=1, nrow=1, as.is=T)
unique(unlist(col.class, use.names=F))
col.class[col.class %in% c("INTEGER", "DECIMAL")] <- "numeric"
col.class[col.class == "STRING"] <- "factor"

header <- read.csv("032_Core_RCA.csv", header=T, as.is=T)
header <- colnames(header)

core.rca <- read.csv("032_Core_RCA.csv", header=F, skip=2, colClasses=col.class)
colnames(core.rca) <- header


# ShaleGas
col.class <- read.csv("033_Core_ShaleGas.csv", header=F, skip=1, nrow=1, as.is=T)
unique(unlist(col.class, use.names=F))
col.class[col.class %in% c("INTEGER", "DECIMAL")] <- "numeric"
col.class[col.class == "STRING"] <- "factor"

header <- read.csv("033_Core_ShaleGas.csv", header=T, as.is=T)
header <- colnames(header)

core.shalegas <- read.csv("033_Core_ShaleGas.csv", header=F, skip=2, colClasses=col.class)
colnames(core.shalegas) <- header


# SCAL
col.class <- read.csv("034_Core_SCAL.csv", header=F, skip=1, nrow=1, as.is=T)
unique(unlist(col.class, use.names=F))
col.class[col.class %in% c("INTEGER", "DECIMAL")] <- "numeric"
col.class[col.class == "STRING"] <- "factor"

header <- read.csv("034_Core_SCAL.csv", header=T, as.is=T)
header <- colnames(header)

core.scal <- read.csv("034_Core_SCAL.csv", header=F, skip=2, colClasses=col.class)
colnames(core.scal) <- header



# Load master data ---------------------------------------------------------------------------------------
a <- read.csv("Master_Aug2013.csv", header=F, skip=3, as.is=T)
a <- a[, c(1,27,36,37)]
prod.mast.dup <- a[duplicated(a[,1]), ]
prod.mast.uniq <- a[!duplicated(a[,1]),]
names(prod.mast.uniq)[1] <- "uwi"
names(prod.mast.uniq)[2] <- "fluid.type"
names(prod.mast.uniq)[3] <- "oil.12m"
names(prod.mast.uniq)[4] <- "gas.12m"

oil.type <- c("Black Oil", "Volatile Oil", "Condensate")
prod.mast.uniq.oil <- prod.mast.uniq[prod.mast.uniq$fluid.type %in% oil.type, 1:3]
prod.mast.uniq.oil <- prod.mast.uniq.oil[complete.cases(prod.mast.uniq.oil),]

write.table(prod.mast.uniq.oil, file="./analysis/EagleFordMasterOilUniq.csv", row.names=F, sep=",")



# Restore working directory -------------------------------------------------------------------------
setwd(path.old)
rm(path.old)
