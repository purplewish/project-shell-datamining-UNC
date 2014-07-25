
path <- "Z:/Mingqi.Wu/project/GitHub/project-shell-datamining-UNC/code/EagleFord/"
setwd(path)

### Load data ========================================================================
source("read_core_data.R")


### Explore ==========================================================================
# Dimension explore
length(unique(core.dat$Unique.Private.Well.ID)) # 84
length(unique(core.dat$Unique.Well.ID)) # 83
length(unique(core.dat$Unique.Sample.ID)) # 2599
dim(core.dat) # 2600 x78

names(core.dat)
