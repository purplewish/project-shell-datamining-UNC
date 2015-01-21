
options(scipen=999)  

# Validation of EF Oil Data Analysis Results -----------------------------------------------------------------------------------
# 1. Kaggle claims their model correctly identified 1.63x the number of top quartile wells 
#    when compared to a rules-based approach using the data they received by Jan 2014.
# 2. This claim is stated in the email communications between Shell and Kaggle in Sept. 2014, 


# Rule-based Approach
path <- "C:/Apps/projects/DataMiningUNC/kaggle/Final/RulesBasedApproach Oct 8/RulesBasedApproach Oct 8/"
setwd(path)

ef.oil.dat <- read.csv("Rules features as of Jan 2014.csv", header=T, as.is=T)
a <- ef.oil.dat 
a[is.na(a)] <- 999999  # Missing data will be assigned with lowest score based on Kaggle's inputs

a$TOC.S <- 100
a[a[,2]<=8, 13] <- 50
a[a[,2]<=6, 13] <- 0
a[a[,2]<=4, 13] <- -50
a[a[,2]<=3, 13] <- -100
a[a[,2]==999999, 13] <- -100

a$RoMeasured.S    <- 100
a[a[,3]<=1, 14]   <- 50
a[a[,3]<=0.9, 14] <- 0
a[a[,3]<=0.8, 14] <- -50
a[a[,3]<=0.75, 14]<- -100
a[a[,3]==999999, 14] <- -100

a$oilApiGravity.S<- 100
a[a[,4]<=40, 15] <- 50
a[a[,4]<=35, 15] <- 0
a[a[,4]<=30, 15] <- -50
a[a[,4]<=25, 15] <- -100
a[a[,4]==999999, 15] <- -100

a$GriSaturationsSo.S<- 100
a[(a[,5]/100)<=0.7, 16] <- 50
a[(a[,5]/100)<=0.6, 16] <- 0
a[(a[,5]/100)<=0.5, 16] <- -50
a[(a[,5]/100)<=0.4, 16] <- -100
a[a[,5]==999999, 16] <- -100

a$ln_formationThickness.S<- 100
a[(exp(a[,6]))<=150, 17] <- 50
a[(exp(a[,6]))<=100, 17] <- 0
a[(exp(a[,6]))<=50, 17] <- -50
a[(exp(a[,6]))<=15, 17] <- -100
a[a[,6]==999999, 17] <- -100

a$GriTotalPorosity.S<- 100
a[a[,7]<=10,18] <- 50
a[a[,7]<=8, 18] <- 0
a[a[,7]<=6, 18] <- -50
a[a[,7]<=4, 18] <- -100
a[a[,7]==999999, 18] <- -100

a$XrdQuartz.S    <- 100
a[a[,8]<=70,19]  <- 50
a[a[,8]<=60, 19] <- 0
a[a[,8]<=50, 19] <- -50
a[a[,8]<=30, 19] <- -100
a[a[,8]==999999, 19] <- -100

a$mudPressure.S    <- 100
a[a[,9]<=0.7, 20]   <- 50
a[a[,9]<=0.5, 20]  <- 0
a[a[,9]<=0.475,20]<- -50
a[a[,9]<=0.45, 20] <- -100
a[a[,9]==999999, 20] <- -100


rule.s <- a[, 13:20]
a$Score <- rowSums(rule.s)
sum(a[,10] - a[,21])  # 0 => Kaggle's calculated rule-based score is consistent with the specified rules. 


c <- a[,c(1, 11, 12, 21)]  # Uwi, score (Kaggle's prediction), True production, score (rule-based) 
n.q1 <- round(nrow(a) * 0.25)  # number of wells in top quartile = 658

# Top quartile wells based on true production
q1.true <- c[order(c[,3], decreasing=T),]
q1.true <- q1.true[1:n.q1, ]
true.id <- as.data.frame(q1.true[,1])
names(true.id)[1] <- "Top.Q.True"
true.id <- sort(true.id[,1])

# Top quartile wells based on Kaggle's predictive model
q1.pred <- c[order(c[,2], decreasing=T),]
q1.pred <- q1.pred[1:n.q1, ]
pred.id <- as.data.frame(q1.pred[,1])
names(pred.id)[1] <- "Top.Q.Pred"
pred.id <- sort(pred.id[,1])

# overlap between Kaggle's model and true
ov.true.pred <- intersect(true.id, pred.id)  
n.ov.true.pred <- length(ov.true.pred)  # 366


n.ov.true.rule <- NULL  # number of overlap between rule and true
for(i in c(1:1000))  # shuffle 1000 times to get an average estimation due to many tied scores in the rule based approach
{
  
  b <- c[sample(1:nrow(a)), ]
  
  # Top quartile wells based on rule-based approach
  q1.rule <- b[order(b[,4], decreasing=T),]
  q1.rule <- q1.rule[1:n.q1, ]
  rule.id <- as.data.frame(q1.rule[,1])
  names(rule.id)[1] <- "Top.Q.Rule"
  rule.id <- sort(rule.id[,1])
    
  # overlap between rule and true
  ov.true.rule <- intersect(true.id, rule.id)
  n.ov.true.rule <- c(n.ov.true.rule, length(ov.true.rule))
  
}

# Kaggle’s model correctly identified 1.63x the number of  top quartile wells when compared to a rules-based approach.
n.ov.true.pred/mean(n.ov.true.rule)  # ~1.63

