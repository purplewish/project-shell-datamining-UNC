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
# sol <- readRDS("sweetspot_5CV_reg_mse.rds")
# saveRDS(pred, "sweetspot_5CV_reg_pred.rds")


#@@ Plot Sweetspot
n.q4 <- ceiling(nrow(pred)*0.25)
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

q.rec <- qRecCurv(x) * 100
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

q.rec <- qRecCurv(x) * 100
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
  xlab("Top Quantile Percentage") + ylab("Recover Rate") + 
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



