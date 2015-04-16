# code path
setwd("Z:/GitHup/project-shell-datamining-UNC/code/RF")
source("header.R")
source("loadData.R")
source("runBart.R")
source("runRF.R")
source("plotFuns.R")
# Results directory
setwd(file.path(repo_path, "Code/RF/results"))


#------------------------------------------------------------------------------------------------------------------------------------------------------------------
## bartMachine model on one set of pars 
set_bart_machine_num_cores(1)  # seed only work for 1 core
bart <- bartMachine(all[ ,2:32], all$Target.Q4, num_trees=50, num_burn_in=500, num_iterations_after_burn_in=1000, seed=666)  # 50 trees give best results many times

print(bart)
#bart$confusion_matrix
rf.mod <- readRDS("rfMod.rds")


#-------------------------------------------------------------------------------------------------------------------------------
# Variables importance ranking
var.imp <- var_selection_by_permute(bart)

saveRDS(var.imp, "BartVarImp.rds")
#var.imp <- readRDS("BartVarImp.rds")

print(var.imp$important_vars_local_names)
print(var.imp$var_true_props_avg)


#-------------------------------------------------------------------------------------------------------------------------------
# BART model: effect of different train % 
train.pct.seq <- seq(0.1,0.9,0.1)

sol.all <- NULL
for(train.pct in train.pct.seq){
  bart <- runBart(all, train.pct=train.pct)
  sol.all <- rbind(sol.all, bart[[1]])
}
sol.all
#saveRDS(sol.all, "BartTrainPct.rds")
sol.all <- readRDS("BartTrainPct.rds")


# Partial dependence plot
pd_plot(bart, "Producer.DepthTrueVertical.Joined")



#-------------------------------------------------------------------------------------------------------------------------------
# BART model: regression
set_bart_machine_num_cores(1)  # seed only work for 1 core
bart <- bartMachine(all[,2:32], all$Target, num_trees=150, num_burn_in=500, num_iterations_after_burn_in=1000, k=3, q=0.9, nu=3,seed=666)  # for regression, more trees seems more accurate
print(bart)
#saveRDS(bart, "bart_fulldata.rds")

var_selection_by_permute_cv(bart)  # check which selection rule is the best use CV
rmse_by_num_trees(bart, num_replicates=50)  # ntree=150
bart_machine_cv <- bartMachineCV(all[,2:32], all$Target, num_tree_cvs=150)  # grid-searching of parameters
#saveRDS(bart_machine_cv, "bart_machine_cv.rds")


var.sel <- var_selection_by_permute(bart)
#saveRDS(var.sel, "bart_varsel.rds")
vimp <- var.sel$var_true_props_avg
names(vimp) <- gsub("(Core.|.Kriged|.Joined)",'',names(vimp))
vimp <- data.frame(predictor=names(vimp), score=as.numeric(vimp))
plotRFVarImp3(vimp)

#k_fold_cv(all[,2:32], all$Target, k_folds=5)


# 5-fold CV with 100% data
bart <- runBartRegCV(dat=all, no.tree=150, no.burn=500, no.after.burn=1000, k.fold=5, k=3, rev=FALSE, seed=666)
sol <- bart[[1]]  # CV pred accuracy 
pred <- bart[[2]]  # CV pred results

#saveRDS(sol, "bart_sweetspot_5CV_reg_mse.rds")
# sol <- readRDS("bart_sweetspot_5CV_reg_mse.rds")
#saveRDS(pred, "bart_sweetspot_5CV_reg_pred.rds")
# pred <- readRDS("bart_sweetspot_5CV_reg_pred.rds")
n.q4 <- ceiling(nrow(pred)*0.25)

q4.true <- pred %>% top_n(n.q4, Target)
q4.pred <- pred %>% top_n(n.q4, Pred)

true.p  <- intersect(q4.true, q4.pred)  #true +  465
false.n <- setdiff(q4.true, true.p)     #false - 193
false.p <- setdiff(q4.pred, true.p)     #false + 193



#-------------------------------------------------------------------------------------------------------------------------------
# RF model: effect of different train % using Cross Validation
K <- c(10, 5, 3, 2, 3, 5, 10)  # K-fold CV
Rev <- c(T, T, T, T, F, F, F)  # Reverse CV or not
Train.Perc <- c(0.1, 0.2, 0.333, 0.5, 0.66, 0.8, 0.9)  # Training pct <-> K-fold CV

sol.all <- NULL
for(i in 1:length(K)){
  sol <- runBartRegCV(dat=all, no.tree=150, no.burn=500, no.after.burn=1000, k.fold=K[i], k=3, rev=Rev[i], seed=666)
  sol.all <- rbind(sol.all, sol[[1]])
}
sol.all
saveRDS(sol.all, "bart_train_pct_reg_mse_cv.rds")
#sol.all <- readRDS("bart_train_pct_reg_mse_cv.rds")
sol <- as.data.frame(cbind(Train.Perc=Train.Perc, mse=sol.all[,3]))
plotLine(sol, "Percentage of Training Data", "MSE") 




#-------------------------------------------------------------------------------------------------------------------------
## RF regression model 5-fold CV using default mtry setting with top n importance vars
var.sel <- readRDS("bart_varsel.rds")
vimp <- var.sel$var_true_props_avg
vimp <- data.frame(predictor=names(vimp), score=as.numeric(vimp))
sols <- NULL

set.seed(777)
for (top.n in nrow(vimp):1){
  sel.vars.mse <- vimp$predictor[1:top.n]  
  formula.reg.imp.mse <- formula(paste("Target~", paste(sel.vars.mse, collapse="+")))
  
  rf.mse <- runRFRegCV(dat=all, model=formula.reg.imp.mse, m=12, no.tree=1000, k=5, default=TRUE)  # default mtry setting (m!=12)
  
  sol.mse <- data.frame(method="mse", topn=top.n, rf.mse[[1]])
  sols <- rbind(sols, sol.mse)
  
  print(paste0("top.n=",top.n), quote=F)
}

#saveRDS(sols, "reg_topn_impvars_5CV_default_mtry.rds")
sols <- readRDS("reg_topn_impvars_5CV_default_mtry.rds")

mse <- sols %>% select(topn, mse) %>% arrange(topn)
plotLine(mse, "Top K Important Variables", "MSE") 



# y.hat <- predict(bart, all[,2:32])
# cred.int <- calc_credible_intervals(bart, all[,2:32], ci_conf=0.8)  # P10, P90
# pred <- cbind(all$Target, y.hat, cred.int)
# print(head(cred.int))


