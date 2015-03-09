# code path
setwd("Z:/GitHup/project-shell-datamining-UNC/code/RF")
source("header.R")
source("loadData.R")
source("runBart.R")
source("plotFuns.R")
# Results directory
setwd(file.path(repo_path, "Code/RF/results"))


#-------------------------------------------------------------------------------------------------------------------------
## bartMachine model on one set of pars 
set_bart_machine_num_cores(5)  # seed only work for 1 core

#bart = build_bart_machine(all[ ,2:32], all$Target.Q4, num_trees=50, num_burn_in=250, num_iterations_after_burn_in=1000)
bart = bartMachine(all[ ,2:32], all$Target.Q4, num_trees=50, num_burn_in=500, num_iterations_after_burn_in=1000, seed=666)  # 50 trees give best results many times

bart$confusion_matrix
#rf.mod <- readRDS("rfMod.rds")


#-------------------------------------------------------------------------------------------------------------------------------
# RF model: effect of different train % 
train.pct.seq <- seq(0.1,0.9,0.1)

sol.all <- NULL
for(train.pct in train.pct.seq){
  bart <- runBart(all, train.pct=train.pct)
  sol.all <- rbind(sol.all, bart[[1]])
}
sol.all


# Partial dependence plot
pd_plot(bt, "Producer.DepthTrueVertical.Joined")
