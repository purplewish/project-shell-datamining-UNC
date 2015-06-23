######################################################################################################################
# RF Regression
######################################################################################################################

#---------------------------------------------------------------------------------------------------------------------
### Setup
#---------------------------------------------------------------------------------------------------------------------
# code path
setwd("Z:/GitHup/project-shell-datamining-UNC/code/RF")

source("header.R")
source("loadData.R")
source("runRF.R")
source("plotFuns.R")
source("loadChronIHSData.R")
source("loadChronIHSData2.R")

# Results directory
setwd(file.path(repo_path, "Code/RF/results"))



#---------------------------------------------------------------------------------------------------------------------
### Parameter tuning
#---------------------------------------------------------------------------------------------------------------------
#@@ Effect of number of trees (select no.tree=1000)
set.seed(123)
num.tree <- 3000
rf <- randomForest(formula.reg, data=all, mtry=20, ntree=num.tree)  # mtry=10, 5, 20  
#saveRDS(rf, "rf_reg_ntree_m10.rds")
#rf <- readRDS("rf_reg_ntree_m10.rds")
dat <- data.frame(ntree=1:num.tree, mse=rf$mse)
plotLine(dat, "Number of trees", "MSE")


#@@ Effect of mtry (select m=12)
m.seq <- seq(2,20,2)  # mtry=31/3=10
num.tree <- 1000

set.seed(789)
sol.all <- NULL
for(m in m.seq){
  ###############################################################################
  rf <- runRFReg(dat=all, train.pct=1.0, model=formula.reg, m=m, no.tree=num.tree)
  ###############################################################################
  rf.mod <- rf[[2]]
  sol.all <- rbind(sol.all, c(m=m, mse=rf.mod$mse[num.tree]))
}
sol.all
#saveRDS(sol.all, "mtry_train_100pct_reg.rds")
#sol.all <- readRDS("mtry_train_100pct_reg.rds")



#----------------------------------------------------------------------------------------------------------------------
### RF Reg using all data without CV (model + partial plot)
#----------------------------------------------------------------------------------------------------------------------
set.seed(123)
num.tree <- 1000
m <- 12
rf <- randomForest(formula.reg, data=all, mtry=m, ntree=num.tree)
# saveRDS(rf, "rf_reg_ntree1k_m12.rds")
# rf <- readRDS("rf_reg_ntree1k_m12.rds")


#@@ RF model:  partial plot
m <-12; num.tree <- 1000;
set.seed(789)
rf <- runRFReg(dat=all, train.pct=1.0, model=formula.reg, m=m, no.tree=num.tree)
rf.mod <- rf[[2]]
sol <- rf[[1]]

par(mfrow=c(3,1))
partialPlot(rf.mod, all, Core.Tmax.Kriged, xlab="Tmax", main="")
partialPlot(rf.mod, all, Producer.DepthTrueVertical.Joined, xlab="True Vertical Depth", main="")
partialPlot(rf.mod, all, Core.S2.Kriged, xlab="S2", main="")



#-------------------------------------------------------------------------------------------------------------------------
### RF model: Cross Validation on all data
#-------------------------------------------------------------------------------------------------------------------------
#@@ Effect of different train % using CV
K <- c(10, 5, 3, 2, 3, 5, 10)  # K-fold CV
Rev <- c(T, T, T, T, F, F, F)  # Reverse CV or not
Train.Perc <- c(0.1, 0.2, 0.333, 0.5, 0.66, 0.8, 0.9)  # Training pct <=> K-fold CV

set.seed(999)
sol.all <- NULL; 
pred.all <- NULL;
for(i in 1:length(K)){
  #sol <- runRFRegCV(dat=all, model=formula.reg, m=12, no.tree=1000, k=K[i], rev=Rev[i])
  sol <- runRFRegCV(dat=all, model=formula.reg.imp7, m=12, no.tree=1000, k=K[i], rev=Rev[i])  # parsimonial model line 262
           
  sol.all <- rbind(sol.all, sol[[1]])
  pred.all <- rbind(pred.all, data.frame(sol[[2]], Train.Perc=Train.Perc[i]))
}
#sol.all
#pred.all
#saveRDS(sol.all, "train_pct_reg_mse_cv.rds")
#saveRDS(pred.all, "train_pct_reg_pred_cv.rds")
saveRDS(pred.all, "train_pct_reg_pred_cv_6vmodel.rds")
sol.all <- readRDS("train_pct_reg_mse_cv.rds")
sol <- as.data.frame(cbind(Train.Perc=Train.Perc, rmse=sol.all[,5]))
plotLine(sol, "Percentage of Training Data", "RMSE") 


#@@ 5-fold CV 
set.seed(777)
rf <- runRFRegCV(dat=all, model=formula.reg, m=12, no.tree=1000, k=5)
sol <- rf[[1]]  # CV pred accuracy 
pred <- rf[[2]]  # CV pred results
# saveRDS(sol, "sweetspot_5CV_reg_mse.rds")
sol <- readRDS("sweetspot_5CV_reg_mse.rds")
# saveRDS(pred, "sweetspot_5CV_reg_pred.rds")
pred <- readRDS("sweetspot_5CV_reg_pred.rds")


#@@ Plot Sweetspot
n.q4 <- ceiling(nrow(pred)*0.25)  # 658
q4.true <- pred %>% top_n(n.q4, Target)
q4.pred <- pred %>% top_n(n.q4, Pred)

true.p  <- intersect(q4.true, q4.pred)  #true +  484
false.n <- setdiff(q4.true, true.p)     #false - 174
false.p <- setdiff(q4.pred, true.p)     #false + 174

dat <- pred %>% 
  mutate(Target.Q4=ifelse(Uwi %in% q4.true$Uwi, TRUE, FALSE)) %>%
  mutate(Target.Q4.Pred=ifelse(Uwi %in% q4.pred$Uwi, TRUE, FALSE)) %>%
  select(Uwi, Target.Q4, Latitude, Longitude, Target.Q4.Pred)
plotSweetspot(dat)



#-------------------------------------------------------------------------------------------------------------------------
### RF model on vars importance (avg importance score with rep runs)
#-------------------------------------------------------------------------------------------------------------------------
rep <- 50; train.pct <- 1;
mse <- NULL; purity <-NULL;

set.seed(777)
for (i in 1:rep){
  rf <- runRFReg(dat=all, train.pct=train.pct, model=formula.reg, m=12, no.tree=1000)
  rf.mod <- rf[[2]]   # Last rf model obj
  
  Predictor <- rownames(importance(rf.mod))
  #Predictor <- gsub("(Core.|.Kriged|.Joined)",'',rownames(importance(rf.mod)))
  dat <- data.frame(importance(rf.mod))
  mse<- cbind(mse, dat[,1]) # mean decrease in mse
  purity <- cbind(purity, dat[,2]) # mean decrease in purity(gini)
  
  print(paste0("rep=",i), quote=F)
}

dat <- data.frame(Predictor, mse=rowMeans(mse), purity=rowMeans(purity))  # average nrep runs importance scores
imp.mse <- dat %>% select(Predictor, mse) %>% arrange(desc(mse))   # mean decrease in mse
imp.purity <- dat %>% select(Predictor, purity) %>% arrange(desc(purity))   # mean decrease in gini
# saveRDS(imp.mse, "imp_vars_mse_50rep.rds")
# saveRDS(imp.purity, "imp_vars_purity_50rep.rds")
# imp.mse <- readRDS("imp_vars_mse_50rep.rds")
# imp.purity <- readRDS("imp_vars_purity_50rep.rds")

imp.mse$Predictor <- gsub("(Core.|.Kriged|.Joined)",'',imp.mse$Predictor)
imp.purity$Predictor <- gsub("(Core.|.Kriged|.Joined)",'',imp.purity$Predictor)
plotRFVarImp3(imp.mse)
plotRFVarImp3(imp.purity)



#-------------------------------------------------------------------------------------------------------------------------
### RF regression model 5-fold CV using default mtry setting with top K importance vars
#-------------------------------------------------------------------------------------------------------------------------
#@@ 5-fold CV with top K important vars
imp.mse <- readRDS("imp_vars_mse_50rep.rds")
sols <- NULL

set.seed(777)
for (top.n in nrow(imp.mse):1){
  sel.vars.mse <- imp.mse$Predictor[1:top.n]  
  formula.reg.imp.mse <- formula(paste("Target~", paste(sel.vars.mse, collapse="+")))
  
  rf <- runRFRegCV(dat=all, model=formula.reg.imp.mse, m=12, no.tree=1000, k=5, default=TRUE)  # default mtry setting (m!=12)
  
  sol.rmse <- data.frame(method="rmse", topn=top.n, rf[[1]])
  sols <- rbind(sols, sol.rmse)
  
  print(paste0("top.n=",top.n), quote=F)
}
# saveRDS(sols, "reg_topn_impvars_5CV_default_mtry.rds")
# sols <- readRDS("reg_topn_impvars_5CV_default_mtry.rds")

mse <- sols %>% select(topn, rmse) %>% arrange(topn)
plotLine(mse, "Top K Important Variables", "RMSE") 


#@@ Correlation
imp.mse <- readRDS("imp_vars_mse_50rep.rds")
                 
top.n <- 11
top.n.x <- imp.mse$Predictor[1:top.n]
top.n.x <- gsub("(Core.|.Kriged|.Joined)",'',top.n.x)  # strip unneccsary characters

X <- select(all, -Uwi, -Target, -Target.Q, -Target.Q4, -Latitude, -Longitude, -Date.Production.Start)
names(X) <- gsub("(Core.|.Kriged|.Joined)",'',names(X))  # strip unneccsary characters
X.top.n <- select(X, match(top.n.x, names(X)))
corr <- as.matrix(cor(X.top.n))
X.top.n.rmRoCal <- select(X.top.n, -RoCalculated)

plotCorr(X.top.n.rmRoCal)
corr <- as.matrix(cor(X.top.n.rmRoCal))


#@@ 5-fold CV with top 11 important vars but remove high correlation ones
imp.mse <- readRDS("imp_vars_mse_50rep.rds")

top.n <- 11
top.n.x <- imp.mse$Predictor[1:top.n]
top.n.x.rmRoCal <-  top.n.x[-4]
top.n.x.rmRoCalGriSo <-  top.n.x.rmRoCal[-10]
top.n.x.rmRo <-  top.n.x[-c(4,6)]
top.n.x.rmRoGriSo <- top.n.x.rmRo[-9]

formula.reg.imp <- formula(paste("Target~", paste(top.n.x.rmRoGriSo,collapse="+")))
formula.reg.imp1 <- formula(paste("Target~", paste(top.n.x.rmRoCal,collapse="+")))
formula.reg.imp0 <- formula(paste("Target~", paste(top.n.x.rmRoCalGriSo,collapse="+")))

set.seed(77)
rf.mse.imp <- runRFRegCV(dat=all, model=formula.reg.imp, m=12, no.tree=1000, k=5, default=TRUE)  # default mtry setting (m!=12)

set.seed(77)
rf.mse.imp1 <- runRFRegCV(dat=all, model=formula.reg.imp1, m=12, no.tree=1000, k=5, default=TRUE)  # default mtry setting (m!=12)

set.seed(77)
rf.mse.imp0 <- runRFRegCV(dat=all, model=formula.reg.imp0, m=12, no.tree=1000, k=5, default=TRUE)  # default mtry setting (m!=12)


top.n <- 11
top.n.x <- imp.mse$Predictor[1:top.n]
top.n.x.sel <- top.n.x[c(1,2,3,10)]

formula.reg.imp4 <- formula(paste("Target~", paste(top.n.x.sel,collapse="+")))

set.seed(77)
rf.mse.imp4 <- runRFRegCV(dat=all, model=formula.reg.imp4, m=12, no.tree=1000, k=5, default=TRUE)  # default mtry setting (m!=12)


#@@ 5-fold CV with top 6 important vars but remove high correlation ones
imp.mse <- readRDS("imp_vars_mse_50rep.rds")

top.n <- 6
top.n.x <- imp.mse$Predictor[1:top.n]
top.n.x.rmRoCal <-  top.n.x[-4]
top.n.x.rmRo <-  top.n.x.rmRoCal[-5]

formula.reg.imp2 <- formula(paste("Target~", paste(top.n.x.rmRoCal,collapse="+")))
formula.reg.imp3 <- formula(paste("Target~", paste(top.n.x.rmRo,collapse="+")))

set.seed(77)
rf.mse.imp2 <- runRFRegCV(dat=all, model=formula.reg.imp2, m=12, no.tree=1000, k=5, default=TRUE)  # default mtry setting (m!=12)

set.seed(77)
rf.mse.imp3 <- runRFRegCV(dat=all, model=formula.reg.imp3, m=12, no.tree=1000, k=5, default=TRUE)  # default mtry setting (m!=12)


#@@ 5-fold CV with top 7 important vars but remove high correlation ones (6 Vars model)
imp.mse <- readRDS("imp_vars_mse_50rep.rds")

top.n <- 7
top.n.x <- imp.mse$Predictor[1:top.n]
top.n.x.rmRoCal <-  top.n.x[-4]

formula.reg.imp7 <- formula(paste("Target~", paste(top.n.x.rmRoCal,collapse="+")))

set.seed(77)
rf.mse.imp7 <- runRFRegCV(dat=all, model=formula.reg.imp7, m=12, no.tree=1000, k=5, default=TRUE)  # default mtry setting (m!=12)

sol <- rf.mse.imp7[[1]]
pred <- rf.mse.imp7[[2]]
# saveRDS(sol, "sweetspot_5CV_reg_6impvars.rds")
# sol <- readRDS("sweetspot_5CV_reg_6impvars.rds")
# saveRDS(pred, "sweetspot_5CV_reg_6impvars_pred.rds")


n.q4 <- ceiling(nrow(pred)*0.25)

q4.true <- pred %>% top_n(n.q4, Target)
q4.pred <- pred %>% top_n(n.q4, Pred)

true.p  <- intersect(q4.true, q4.pred)  #true +  490
false.n <- setdiff(q4.true, true.p)     #false - 168
false.p <- setdiff(q4.pred, true.p)     #false + 168
true.p



#-------------------------------------------------------------------------------------------------------------------------
### Recover Curve
#-------------------------------------------------------------------------------------------------------------------------
#@@ Comparison of different model
# RF full model (31 vars)
pred <- readRDS("sweetspot_5CV_reg_pred.rds")
pred <- select(pred, Uwi, Target, RF=Pred)

# RF model (6 vars)
pred.6v <- readRDS("sweetspot_5CV_reg_6impvars_pred.rds")
pred.6v <- select(pred.6v, Uwi, RF6v=Pred)

pred.kaggle <- select(all, Uwi, Rules.Prediction, Kaggle.Prediction)

x <- left_join(pred, pred.6v, by="Uwi")
x <- left_join(x, pred.kaggle, by="Uwi")
x <- x[,-1]  # rm Uwi

q.rec0 <- qRecCurv(x) * 100
#q.rec <- q.rec0[-c(1:7),]
# Round to integer percentage
index <- ceiling(nrow(q.rec0)*seq(0.3,100,0.3)/100)
q.rec <- q.rec0[index, ]

q.rec1 <- q.rec %>% select(True) %>% mutate(RecRate=True, Method="1 Baseline")
q.rec2 <- q.rec %>% select(True, X2) %>% rename(RecRate=X2) %>% mutate(Method="4 Random Forest")
q.rec3 <- q.rec %>% select(True, X3) %>% rename(RecRate=X3) %>% mutate(Method="5 Random Forest (6 Vars)")
q.rec4 <- q.rec %>% select(True, X4) %>% rename(RecRate=X4) %>% mutate(Method="2 Rule Based")
q.rec5 <- q.rec %>% select(True, X5) %>% rename(RecRate=X5) %>% mutate(Method="3 Kaggle")

q.rec <- union(q.rec1, q.rec2)
q.rec <- union(q.rec, q.rec3)
q.rec <- union(q.rec, q.rec4)
q.rec <- union(q.rec, q.rec5)

ggplot(q.rec, aes(x=True, y=RecRate, colour=Method, group=Method)) + 
  geom_line(lwd=1.2) +
  scale_color_manual(values=c("#fe506e", "black", "#228b22", "#0099cc", "#e95d3c")) +
  xlab("Top Quantile Percentage") + ylab("Recover Rate") + 
  theme(#legend.position="none",
    axis.title.x = element_text(size=24),
    axis.title.y = element_text(size=24),
    axis.text.x = element_text(colour="grey20",size=15),
    axis.text.y = element_text(colour="grey20",size=15),
    legend.title=element_blank(),
    legend.text = element_text(size = 20),
    legend.justification=c(1,0), legend.position=c(1,0),
    legend.background = element_rect(fill="gray90", size=.5, linetype="dotted")
  )
# plot(q.rec, type="l", xlab="Top Quantile Percentage", ylab="Recover Rate")
# lines(q.rec[,1],q.rec[,1], col="red")


#@@ Comparison of different training percentage for the full model (31 vars)
#pred <- readRDS("train_pct_reg_pred_cv.rds")  # full model
pred <- readRDS("train_pct_reg_pred_cv_6vmodel.rds")  # parsimonial model
pred.10pct <- pred %>% filter(Train.Perc==0.1) %>% select(Uwi, Target, Pred)
pred.20pct <- pred %>% filter(Train.Perc==0.2) %>% select(Uwi, Pred)
pred.33pct <- pred %>% filter(Train.Perc==0.333) %>% select(Uwi, Pred)
pred.50pct <- pred %>% filter(Train.Perc==0.5) %>% select(Uwi, Pred)
pred.66pct <- pred %>% filter(Train.Perc==0.66) %>% select(Uwi, Pred)
pred.80pct <- pred %>% filter(Train.Perc==0.8) %>% select(Uwi, Pred)
pred.90pct <- pred %>% filter(Train.Perc==0.9) %>% select(Uwi, Pred)

x <- pred.10pct %>%
      left_join(pred.20pct, by="Uwi") %>% rename(pct10=Pred.x) %>% rename(pct20=Pred.y) %>%
      left_join(pred.33pct, by="Uwi") %>% rename(pct33=Pred) %>%
      left_join(pred.50pct, by="Uwi") %>% rename(pct50=Pred) %>%
      left_join(pred.66pct, by="Uwi") %>% rename(pct66=Pred) %>%
      left_join(pred.80pct, by="Uwi") %>% rename(pct80=Pred) %>%
      left_join(pred.90pct, by="Uwi") %>% rename(pct90=Pred) %>%
      select(-Uwi)

# Top quartile
#q.rec0 <- qRecCurv(x) * 100
#q.rec <- q.rec0[-c(1:10),]
# Round to integer percentage

# Bottom quartile
q.rec0 <- qRecCurvBottom(x) * 100

index <- c(ceiling(nrow(q.rec0)*seq(0.2, 100, 0.2)/100))
q.rec <- q.rec0[index, ]

q.rec1 <- q.rec %>% select(True) %>% mutate(RecRate=True, Method="1 Baseline")
q.rec2 <- q.rec %>% select(True, X2) %>% rename(RecRate=X2) %>% mutate(Method="2 10% Training")
q.rec3 <- q.rec %>% select(True, X3) %>% rename(RecRate=X3) %>% mutate(Method="3 20% Training")
q.rec4 <- q.rec %>% select(True, X4) %>% rename(RecRate=X4) %>% mutate(Method="4 33% Training")
q.rec5 <- q.rec %>% select(True, X5) %>% rename(RecRate=X5) %>% mutate(Method="5 50% Training")
q.rec6 <- q.rec %>% select(True, X6) %>% rename(RecRate=X6) %>% mutate(Method="6 66% Training")
q.rec7 <- q.rec %>% select(True, X7) %>% rename(RecRate=X7) %>% mutate(Method="7 80% Training")
q.rec8 <- q.rec %>% select(True, X8) %>% rename(RecRate=X8) %>% mutate(Method="8 90% Training")

q.rec <- q.rec1 %>% union(q.rec2) %>%
                    union(q.rec3) %>%
                    union(q.rec4) %>%
                    union(q.rec5) %>%
                    union(q.rec6) %>%
                    union(q.rec7) %>%
                    union(q.rec8)
                    


ggplot(q.rec, aes(x=True, y=RecRate, colour=Method, group=Method)) + 
  geom_line(lwd=1.2) +
  scale_color_manual(values=c("#787777", "#FBDB0C", "#5E9F37", "#007cd2", "#333333", "#FF6600", "#FF1CAE", "#ff0000")) +
  #xlab("Top Quantile Percentage") + ylab("Recover Rate") + 
  xlab("Bottom Quantile Percentage") + ylab("Recover Rate") + 
  #scale_y_continuous(limits=c(50, 90)) +
  #scale_x_continuous(limits=c(0, 50)) +
  theme(#legend.position="none",
    axis.title.x = element_text(size=24),
    axis.title.y = element_text(size=24),
    axis.text.x = element_text(colour="grey20",size=15),
    axis.text.y = element_text(colour="grey20",size=15),
    legend.title=element_blank(),
    legend.text = element_text(size = 20),
    legend.justification=c(1,0), legend.position=c(1,0),
    legend.background = element_rect(fill="gray90", size=.5, linetype="dotted")
  )



#@@ Kaggle dataset Comparison of different cutoff date for the full model (31 vars)
pred <- readRDS("cutoff_reg_pred.rds")  # full model
c <- levels(pred$cutoff)
pred.cut1 <- pred %>% filter(cutoff==c[1]) %>% select(Uwi, Target, Pred)
pred.cut2 <- pred %>% filter(cutoff==c[2]) %>% select(Uwi, Pred)
pred.cut3 <- pred %>% filter(cutoff==c[3]) %>% select(Uwi, Pred)

x <- pred.cut1 %>%
  left_join(pred.cut2, by="Uwi") %>% rename(cut1=Pred.x) %>% rename(cut2=Pred.y) %>%
  left_join(pred.cut3, by="Uwi") %>% rename(cut3=Pred) %>%
  select(-Uwi)

# top quartile
# q.rec0 <- qRecCurv(x) * 100
# q.rec <- q.rec0[-c(1:12),]  # start from 0.5%
# Round to integer percentage
# index <- ceiling(nrow(q.rec)*seq(0.1, 100, 0.1)/100)
# q.rec <- q.rec[index, ]

# bottom quartile
q.rec0 <- qRecCurvBottom(x) * 100
q.rec <- q.rec0[-c(1:12),]

q.rec1 <- q.rec %>% select(True) %>% mutate(RecRate=True, Method="1 Baseline")
q.rec2 <- q.rec %>% select(True, X2) %>% rename(RecRate=X2) %>% mutate(Method="2 11% Training ~ 297 producers, 83 cored wells by 2011-11-01")
q.rec3 <- q.rec %>% select(True, X3) %>% rename(RecRate=X3) %>% mutate(Method="3 16% Training ~ 418 producers, 83 cored wells by 2012-01-01")
q.rec4 <- q.rec %>% select(True, X4) %>% rename(RecRate=X4) %>% mutate(Method="4 21% Training ~ 542 producers, 83 cored wells by 2012-03-01")

q.rec <- q.rec1 %>% union(q.rec2) %>% union(q.rec3) %>% union(q.rec4)


ggplot(q.rec, aes(x=True, y=RecRate, colour=Method, group=Method)) + 
  geom_line(lwd=1.2) +
  scale_color_manual(values=c("#787777", "red", "#5E9F37", "#007cd2")) +
  #xlab("Top Quantile Percentage") + ylab("Recover Rate") + 
  xlab("Bottom Quantile Percentage") + ylab("Recover Rate") + 
  #scale_y_continuous(limits=c(50, 90)) +
  #scale_x_continuous(limits=c(0, 5)) +
  theme(#legend.position="none",
    axis.title.x = element_text(size=24),
    axis.title.y = element_text(size=24),
    axis.text.x = element_text(colour="grey20",size=15),
    axis.text.y = element_text(colour="grey20",size=15),
    legend.title=element_blank(),
    legend.text = element_text(size = 20),
    legend.justification=c(1,0), legend.position=c(1,0),
    legend.background = element_rect(fill="gray90", size=.5, linetype="dotted")
  )



#@@ IHS dataset Comparison of different cutoff date for the full model
pred <- readRDS("cutoff_reg_pred_IHS_m12.rds")  # full model
c <- levels(pred$cutoff)
pred.cut1 <- pred %>% filter(cutoff==c[1]) %>% select(Entity, First.12.Month.Liquid.norm.perforation, Pred)
pred.cut2 <- pred %>% filter(cutoff==c[2]) %>% select(Entity, Pred)
pred.cut3 <- pred %>% filter(cutoff==c[3]) %>% select(Entity, Pred)
pred.cut4 <- pred %>% filter(cutoff==c[4]) %>% select(Entity, Pred)

x <- pred.cut1 %>%
  left_join(pred.cut2, by="Entity") %>% rename(cut1=Pred.x) %>% rename(cut2=Pred.y) %>%
  left_join(pred.cut3, by="Entity") %>% rename(cut3=Pred) %>%
  left_join(pred.cut4, by="Entity") %>% rename(cut4=Pred) %>%
  select(-Entity) %>% rename(Target=First.12.Month.Liquid.norm.perforation)

# top quartile
q.rec0 <- qRecCurv(x) * 100
q.rec <- q.rec0[-c(1:12),]  # start from 0.5%
# Round to integer percentage
# index <- ceiling(nrow(q.rec)*seq(0.1, 100, 0.1)/100)
# q.rec <- q.rec[index, ]

# bottom quartile
# q.rec0 <- qRecCurvBottom(x) * 100
# q.rec <- q.rec0[-c(1:12),]

q.rec1 <- q.rec %>% select(True) %>% mutate(RecRate=True, Method="1 Baseline")
q.rec2 <- q.rec %>% select(True, X2) %>% rename(RecRate=X2) %>% mutate(Method="2 5%   Training ~ 233 producers, 39 cored wells by 2012-02-01")
q.rec3 <- q.rec %>% select(True, X3) %>% rename(RecRate=X3) %>% mutate(Method="3 10% Training ~ 477 producers, 50 cored wells by 2012-06-01")
q.rec4 <- q.rec %>% select(True, X4) %>% rename(RecRate=X4) %>% mutate(Method="4 16% Training ~ 741 producers, 55 cored wells by 2012-09-01")
q.rec5 <- q.rec %>% select(True, X5) %>% rename(RecRate=X5) %>% mutate(Method="5 21% Training ~ 952 producers, 57 cored wells by 2012-11-01")

q.rec <- q.rec1 %>% union(q.rec2) %>% union(q.rec3) %>% union(q.rec4) %>% union(q.rec5)
  

ggplot(q.rec, aes(x=True, y=RecRate, colour=Method, group=Method)) + 
  geom_line(lwd=1.2) +
  scale_color_manual(values=c("#787777", "orange", "#5E9F37", "#007cd2", "red")) +
  xlab("Top Quantile Percentage") + ylab("Recover Rate") + 
  #xlab("Bottom Quantile Percentage") + ylab("Recover Rate") + 
  #scale_y_continuous(limits=c(50, 90)) +
  #scale_x_continuous(limits=c(0, 5)) +
  theme(#legend.position="none",
    axis.title.x = element_text(size=24),
    axis.title.y = element_text(size=24),
    axis.text.x = element_text(colour="grey20",size=15),
    axis.text.y = element_text(colour="grey20",size=15),
    legend.title=element_blank(),
    legend.text = element_text(size = 20),
    legend.justification=c(1,0), legend.position=c(1,0),
    legend.background = element_rect(fill="gray90", size=.5, linetype="dotted")
  )




#-------------------------------------------------------------------------------------------------------------------------
### Time dependent study
#-------------------------------------------------------------------------------------------------------------------------

#@@ Kaggle's Data: Effect of different cut off date <=> different training %

# Cut off date
cutoff <- c("2010-11-01", "2011-01-01", "2011-03-01")
# sum(all$Date.Production.Start<="2010-11-01")  # 297 11.3%   1118 42.5% > cutoff+12month
# sum(all$Date.Production.Start<="2011-01-01")  # 418 15.9%   825  31.3%
# sum(all$Date.Production.Start<="2011-03-01")  # 542 20.6%   502  19.1%

set.seed(999); sol.all <- NULL; pred.all <- NULL;
for(i in 1:length(cutoff)){  
  sol <- runRFReg2(all, cutoff[i], model=formula.reg, m=12, no.tree=1000, ntrace=500)
  sol.all <- rbind(sol.all, sol[[1]])
  pred.all <- rbind(pred.all, data.frame(sol[[2]], cutoff=cutoff[i]))
}
#saveRDS(sol.all, "cutoff_reg_mse.rds")
#saveRDS(pred.all, "cutoff_reg_pred.rds")



#@@ IHS & Core Lab data: Effect of different cut off date <=> different training %

# Cut off date
cutoff <- c("2011-02-01", "2011-06-01", "2011-09-01", "2011-11-01")

set.seed(999);  sol.all <- NULL; pred.all <- NULL;
sol1 <- runRFReg3(d5p,  cutoff=cutoff[1], model=formula.reg.IHS.chron, m=12, no.tree=1000, ntrace=500)
sol2 <- runRFReg3(d10p, cutoff=cutoff[2], model=formula.reg.IHS.chron, m=12, no.tree=1000, ntrace=500)
sol3 <- runRFReg3(d15p, cutoff=cutoff[3], model=formula.reg.IHS.chron, m=12, no.tree=1000, ntrace=500)
sol4 <- runRFReg3(d20p, cutoff=cutoff[4], model=formula.reg.IHS.chron, m=12, no.tree=1000, ntrace=500)
sol.all <- rbind(sol.all, sol1[[1]], sol2[[1]], sol3[[1]], sol4[[1]])
pred.all <- rbind(pred.all, 
                  data.frame(sol1[[2]], cutoff=cutoff[1]),
                  data.frame(sol2[[2]], cutoff=cutoff[2]),
                  data.frame(sol3[[2]], cutoff=cutoff[3]),
                  data.frame(sol4[[2]], cutoff=cutoff[4])
                  )

saveRDS(sol.all, "cutoff_reg_mse_IHS_m12_noloc.rds")
saveRDS(pred.all, "cutoff_reg_pred_IHS_m12_noloc.rds")



#@@ IHS & Core Lab data: Effect of different cut off date <=> different training %
# Top Q1 by different training percentage

# RF using Loc + Cov (m=13)
set.seed(777)
q1.rec <- NULL
for (i in 1:nrow(cut.info)){
    cut <- as.POSIXlt(as.Date(cut.info$cut.off.date[i]))
    cut$year <- cut$year-1
    cut <- as.Date(cut)
    
    d <- chro.dat %>% filter(pct.IHS==cut.info$percentage.IHS.used[i])      
    
    sol <- runRFReg4(d, cutoff=cut, model=formula.reg.chro, m=15, no.tree=1000, ntrace=500)
    
    x <- sol[[2]] %>% select(-UWI) %>% rename(Target=Norm.Lat.12.Month.Liquid, RF=Pred)
    
    q.rec0 <- qRecCurv(x) * 100
    q1.rec <- c(q1.rec, q.rec0[ceiling(nrow(x) * 0.25),2])
}

q1rec.rf <- data.frame(pct.IHS=cut.info$percentage.IHS.used, pct.Core=cut.info$precentage.CoreLabs.used, q1rec.rf=q1.rec)
#saveRDS(q1rec.rf, "q1_rec_rf_reg_diff_cutoff_m15.rds")
q1rec.rf.m13 <- readRDS("q1_rec_rf_reg_diff_cutoff_m13.rds")  # full model



# RF using Loc only
set.seed(777)
q1.rec <- NULL
for (i in 1:nrow(cut.info)){
  cut <- as.POSIXlt(as.Date(cut.info$cut.off.date[i]))
  cut$year <- cut$year-1
  cut <- as.Date(cut)
  
  d <- chro.dat %>% filter(pct.IHS==cut.info$percentage.IHS.used[i])      
  
  sol <- runRFReg4(d, cutoff=cut, model=formula.reg.chro.loc, m=2, no.tree=1000, ntrace=500)
  
  x <- sol[[2]] %>% select(-UWI) %>% rename(Target=Norm.Lat.12.Month.Liquid, RF=Pred)
  
  q.rec0 <- qRecCurv(x) * 100
  q1.rec <- c(q1.rec, q.rec0[ceiling(nrow(x) * 0.25),2])
}

q1rec.rf.loc <- data.frame(pct.IHS=cut.info$percentage.IHS.used, pct.Core=cut.info$precentage.CoreLabs.used, q1rec.rf.loc=q1.rec)
#saveRDS(q1rec.rf.loc, "q1_rec_rf_reg_loconly_diff_cutoff.rds")


# RF w/o Loc 
set.seed(777)
q1.rec <- NULL
for (i in 1:nrow(cut.info)){
  cut <- as.POSIXlt(as.Date(cut.info$cut.off.date[i]))
  cut$year <- cut$year-1
  cut <- as.Date(cut)
  
  d <- chro.dat %>% filter(pct.IHS==cut.info$percentage.IHS.used[i])      
  
  sol <- runRFReg4(d, cutoff=cut, model=formula.reg.chro.noloc, m=12, no.tree=1000, ntrace=500)
  
  x <- sol[[2]] %>% select(-UWI) %>% rename(Target=Norm.Lat.12.Month.Liquid, RF=Pred)
  
  q.rec0 <- qRecCurv(x) * 100
  q1.rec <- c(q1.rec, q.rec0[ceiling(nrow(x) * 0.25),2])
}

q1rec.rf.noloc <- data.frame(pct.IHS=cut.info$percentage.IHS.used, pct.Core=cut.info$precentage.CoreLabs.used, q1rec.rf.noloc=q1.rec)
saveRDS(q1rec.rf.noloc, "q1_rec_rf_reg_noloc_diff_cutoff.rds")



# Kriging
q1.rec <- NULL
for (i in 1:nrow(cut.info)){
  cut <- as.POSIXlt(as.Date(cut.info$cut.off.date[i]))
  cut$year <- cut$year-1
  cut <- as.Date(cut)
  
  d <- chro.dat %>% filter(pct.IHS==cut.info$percentage.IHS.used[i])   
  
  d <- d[d$Date.Production.Start>cut, ]  # out of sample data (test data)
  
  x <- d %>% select(Norm.Lat.12.Month.Liquid, Kriged.Production) %>% rename(Target=Norm.Lat.12.Month.Liquid, kriged=Kriged.Production)
    
  q.rec0 <- qRecCurv(x) * 100
  q1.rec <- c(q1.rec, q.rec0[ceiling(nrow(x) * 0.25),2])
}

q1rec.kriged <- data.frame(pct.IHS=cut.info$percentage.IHS.used, pct.Core=cut.info$precentage.CoreLabs.used, q1rec.kriged=q1.rec)
#saveRDS(q1rec.kriged, "q1_rec_kriged_diff_cutoff.rds")



# Rule-based
q1.rec <- NULL
for (i in 1:nrow(cut.info)){
  cut <- as.POSIXlt(as.Date(cut.info$cut.off.date[i]))
  cut$year <- cut$year-1
  cut <- as.Date(cut)
  
  d <- chro.dat %>% filter(pct.IHS==cut.info$percentage.IHS.used[i])   
  
  d <- d[d$Date.Production.Start>cut, ]  # out of sample data (test data)
  
  x <- d %>% select(Norm.Lat.12.Month.Liquid, Big.Rules.Score) %>% rename(Target=Norm.Lat.12.Month.Liquid, Rule=Big.Rules.Score)
  
  q.rec0 <- qRecCurv(x) * 100
  q1.rec <- c(q1.rec, q.rec0[ceiling(nrow(x) * 0.25),2])
}

q1rec.rule <- data.frame(pct.IHS=cut.info$percentage.IHS.used, pct.Core=cut.info$precentage.CoreLabs.used, q1rec.rule=q1.rec)
#saveRDS(q1rec.rule, "q1_rec_rule_diff_cutoff.rds")



# Kaggle data
d <- kaggle.dat
d <- d[d$Date.Production.Start>'2012-08-05', ]
x <- d %>% select(Norm.Lat.12.Month.Liquid, Kaggle.Prediction.Oil) %>% rename(Target=Norm.Lat.12.Month.Liquid, kaggle=Kaggle.Prediction.Oil)

q.rec0 <- qRecCurv(x) * 100
q1.rec <- q.rec0[ceiling(nrow(x) * 0.25),2]
q1rec.kaggle <- data.frame(pct.IHS=0.504604840436924, pct.Core=0.987951807228916, q1rec.kaggle=q1.rec)
saveRDS(q1rec.kaggle, "q1_rec_kaggle_diff_cutoff.rds")


# Plots
q1rec.rf.m13 <- readRDS("q1_rec_rf_reg_diff_cutoff_m13.rds")
q1rec.rf.loc <- readRDS("q1_rec_rf_reg_loconly_diff_cutoff.rds")
q1rec.kriged <- readRDS("q1_rec_kriged_diff_cutoff.rds")
q1rec.rule <- readRDS("q1_rec_rule_diff_cutoff.rds")
q1rec.kaggle <- readRDS("q1_rec_kaggle_diff_cutoff.rds")

dat <- rbind(data.frame(pct=q1rec.rf.m13[,1], rec=q1rec.rf.m13[,3], method="Random Forest"), 
             data.frame(pct=q1rec.kriged[,1], rec=q1rec.kriged[,3], method="Kriged"),
             data.frame(pct=q1rec.rf.loc[,1], rec=q1rec.rf.loc[,3], method="Random Forest (Location Only)"), 
             data.frame(pct=q1rec.rule[,1], rec=q1rec.rule[,3], method="Rule based"),
             data.frame(pct=q1rec.kaggle[,1], rec=q1rec.kaggle[,3], method="Kaggle")
              )


plotMLine(dat,"Percentage of Available Producers for Training", "Out of Sample Top Quartile Recovery Rate")
