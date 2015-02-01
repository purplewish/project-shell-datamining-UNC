wd <- getwd()
##################################################################################################################################

## Load data file & prepare dataset for RF model
# Load predictors Xs
setwd("C:/Apps/projects/DataMiningUNC/Kaggle/Final/Documentation and Input Files September 22 2014/Documentation and Input Files")

x <- read.csv("EagleFordOilInput.csv")
x <- distinct(x, Uwi)  # rm duplicate records (5222 x 35)
x <- select(x, -Latitude, -Longitude, -Producer.EstimatedLength.Joined)
x.vars <- names(x)
x.vars <- x.vars[-1]  # rm Uwi (Uwi isn't a predictor)

# Load response Y
setwd("C:/Apps/projects/DataMiningUNC/Kaggle/Final/RulesBasedApproach Oct 8/RulesBasedApproach Oct 8")

y <- read.csv("Rules features using recent and Jan 2014 data.csv")
y <- select(y, Uwi, Target)  # 2631 x 2

# Merge X and Y
all <- left_join(x, y, by="Uwi")
all <- filter(all, !is.na(Target))  # rm missing production record

# Classify quantiles based on production
all$Target.Q <- cut(all$Target,breaks=quantile(all$Target),labels=paste0("Q",1:4),include.lowest=T)  # class: Q1 Q2 Q3 Q4
all <- mutate(all, Target.Q4=factor(Target.Q=="Q4"))  # class: Q4 ~Q4 TRUE FALSE

# Classification formula
formula.class1 <- formula(paste("Target.Q~", paste(x.vars,collapse="+")))  # class:Q1 Q2 Q3 Q4
formula.class2 <- formula(paste("Target.Q4~", paste(x.vars,collapse="+"))) # class:Q4 ~Q4, topQ vs. ~topQ

##################################################################################################################################
setwd(wd)
