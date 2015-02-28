# code path
setwd("Z:/Mingqi.Wu/GitHup/project-shell-datamining-UNC/code/RF")
source("header.R")
source("loadData.R")
source("runRF.R")
source("plotFuns.R")
# Results directory
setwd(file.path(repo_path, "Code/RF/results"))


#-------------------------------------------------------------------------------------------------------------------------
## bartMachine model on one set of pars 
set.seed(777)
bart = build_bart_machine(all[ ,2:32], all$Target.Q4, num_trees = 1000)

##look at in-sample confusion matrix
bart$confusion_matrix

#rf.mod <- readRDS("rf.rds")
