
#---------------------------------------------------------------------------------------------------------------------
### Setup
#---------------------------------------------------------------------------------------------------------------------
# code path
setwd("Z:/GitHup/project-shell-datamining-UNC/code/RF")

source("header.R")
source("runRF.R")
source("plotFuns.R")
source("loadData4Crossplot.R")

# Results directory
setwd(file.path(repo_path, "Code/RF/results"))


#---------------------------------------------------------------------------------------------------------------------
### Cross-plot
#---------------------------------------------------------------------------------------------------------------------

plotPredvsAct(a)
fit <- lm(predict ~ true, data=a)
summary(fit) 

fit.res = resid(fit)
plot(a$true, fit.res, ylab="Residuals", xlab="Obs") 
abline(0, 0)   

plotPredvsAct(b)
fit <- lm(predict ~ true, data=b)
summary(fit)


bot.q1 <-quantile(a$true)[2]
bot.q2 <-quantile(a$true)[3]
top.q1 <-quantile(a$true)[4]

d.bot.q1 <- a[a[,2]<=bot.q1,]
plotPredvsAct(d.bot.q1)
fit <- lm(predict ~ true, data=d.bot.q1)
summary(fit) 

d.bot.q2 <- a[a[,2]<=bot.q2,]
plotPredvsAct(d.bot.q2)
fit <- lm(predict ~ true, data=d.bot.q2)
summary(fit) 

d.top.q1 <- a[a[,2]>=top.q1,]
plotPredvsAct(d.top.q1)
fit <- lm(predict ~ true, data=d.top.q1)
summary(fit) 


#---------------------------------------------------------------------------------------------------------------------
### 5-fold CV
#---------------------------------------------------------------------------------------------------------------------
set.seed(777)
rf <- runRFRegCV(dat=d, model=model1, m=8, no.tree=10, k=5)
sol <- rf[[1]]  # CV pred accuracy 
pred <- rf[[2]]  # CV pred results
saveRDS(sol, "shelldat_5CV_rf_m11_IHS33var_mse.rds")
saveRDS(pred, "shelldat_5CV_rf_m11_IHS33var_pred.rds")

#sol <- readRDS("shelldat_5CV_rf_m6_IHS33var_mse.rds")
#pred <- readRDS("shelldat_5CV_rf_m6_IHS33var_pred.rds")

a <- pred %>% select(true=Norm.Lat.12.Month.Liquid, pred=Pred)
plotPredvsAct(a)
fit <- lm(pred ~ true, data=a)
summary(fit) 

