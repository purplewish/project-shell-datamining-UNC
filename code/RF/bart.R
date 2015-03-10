# code path
setwd("Z:/GitHup/project-shell-datamining-UNC/code/RF")
source("header.R")
source("loadData.R")
source("runBart.R")
source("plotFuns.R")
# Results directory
setwd(file.path(repo_path, "Code/RF/results"))


#------------------------------------------------------------------------------------------------------------------------------------------------------------------
## bartMachine model on one set of pars 
set_bart_machine_num_cores(5)  # seed only work for 1 core
bart <- bartMachine(all[ ,2:32], all$Target.Q4, num_trees=50, num_burn_in=500, num_iterations_after_burn_in=1000, seed=666)  # 50 trees give best results many times

print(bart)
#bart$confusion_matrix
#rf.mod <- readRDS("rfMod.rds")


#-------------------------------------------------------------------------------------------------------------------------------
# Variables importance ranking
var.imp <- var_selection_by_permute(bart)

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
#sol.all <- readRDS("BartTrainPct.rds")


# Partial dependence plot
pd_plot(bart, "Producer.DepthTrueVertical.Joined")



#-------------------------------------------------------------------------------------------------------------------------------
# BART model: regression
set_bart_machine_num_cores(1)  # seed only work for 1 core
bart <- bartMachine(all[,2:32], all$Target, num_trees=500, num_burn_in=500, num_iterations_after_burn_in=1000, seed=666)  # for regression, more trees seems more accurate
print(bart)

y.hat <- predict(bart, all[,2:32])
cred.int <- calc_credible_intervals(bart, all[,2:32], ci_conf=0.8)  # P10, P90
pred <- cbind(all$Target, y.hat, cred.int)
print(head(cred.int))


