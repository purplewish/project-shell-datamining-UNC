
path.old <- getwd() # save path

# Load core data =======================================================================================
path <- "Z:/Mingqi.Wu/project/Datamining/data/"
# path <- "V:/project/DataMining/data/"
setwd(path)


# GeoChem
col.class <- read.csv("EF_Core_data_Aug052013_GeoChem_TOC.csv", header=F, skip=1, nrow=1, as.is=T)
unique(unlist(col.class, use.names=F))  # check unique class type
col.class[col.class %in% c("INTEGER", "DECIMAL")] <- "numeric"
col.class[col.class == "STRING"] <- "factor"

header <- read.csv("EF_Core_data_Aug052013_GeoChem_TOC.csv", header=T, as.is=T)
header <- colnames(header)

geochem.dat <- read.csv("EF_Core_data_Aug052013_GeoChem_TOC.csv", header=F, skip=2, colClasses=col.class)
colnames(geochem.dat) <- header


# SCAL
col.class <- read.csv("EF_Core_data_Aug052013_SCAL.csv", header=F, skip=1, nrow=1, as.is=T)
unique(unlist(col.class, use.names=F))
col.class[col.class %in% c("INTEGER", "DECIMAL")] <- "numeric"
col.class[col.class == "STRING"] <- "factor"

header <- read.csv("EF_Core_data_Aug052013_SCAL.csv", header=T, as.is=T)
header <- colnames(header)

scal.dat <- read.csv("EF_Core_data_Aug052013_SCAL.csv", header=F, skip=2, colClasses=col.class)
colnames(scal.dat) <- header


# ShaleGas
col.class <- read.csv("EF_Core_data_Aug052013_ShaleGas.csv", header=F, skip=1, nrow=1, as.is=T)
unique(unlist(col.class, use.names=F))
col.class[col.class %in% c("INTEGER", "DECIMAL")] <- "numeric"
col.class[col.class == "STRING"] <- "factor"

header <- read.csv("EF_Core_data_Aug052013_ShaleGas.csv", header=T, as.is=T)
header <- colnames(header)

shalegas.dat <- read.csv("EF_Core_data_Aug052013_ShaleGas.csv", header=F, skip=2, colClasses=col.class)
colnames(shalegas.dat) <- header


# Production
col.class <- read.csv("EF_Core_data_Aug052013_Prod.csv", header=F, skip=1, nrow=1, as.is=T)
unique(unlist(col.class, use.names=F))
col.class[col.class %in% c("INTEGER", "DECIMAL")] <- "numeric"
col.class[col.class == "STRING"] <- "factor"
col.class[col.class == "DATE"] <- "Date"

header <- read.csv("EF_Core_data_Aug052013_Prod.csv", header=T, as.is=T)
header <- colnames(header)

prod.dat <- read.csv("EF_Core_data_Aug052013_Prod.csv", header=F, skip=2) #, colClasses=col.class)
colnames(prod.dat) <- header
# ======================================================================================================

setwd(path.old) # restore path
rm(path.old)
