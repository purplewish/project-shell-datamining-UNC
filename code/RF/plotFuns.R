## Plot functions

#------------------------------------------------------------------------------------------------------
# Plot vars importance for RF model (two functions)
#------------------------------------------------------------------------------------------------------
plotRFVarImp <- function(rf.mod){
# Plot variable importance of a RF model 
#
# Args:
#   rf.mod: a rf model obj
#
# Returns:
#   Two Plots: 1. based on pred accuracy 2. based on gini index 

  # Importance data
  dat <- data.frame(rownames(importance(rf.mod)),round(importance(rf.mod),2))
  names(dat)[c(1,ncol(dat)-1,ncol(dat))] <- c("Predictor","mda","mdg")
  rownames(dat) <- NULL

  pred.acc  <- select(dat, Predictor, mda)  # mean decrease in accuracy
  pred.gini <- select(dat, Predictor, mdg)  # mean decrease in gini

  # Var importance plot function
  importancePlot <- function(d,ylb,fontsize){
    fontsize <- as.numeric(fontsize)
    d <- d[order(d[,2],decreasing=T),]
    d$Predictor <- factor(as.character(d$Predictor),levels=rev(as.character(d$Predictor)))
    rownames(d) <- NULL
  
    d[,2] <- d[,2]/abs(max(d[,2])) * 100  # normalize relative to the variable with maximum score
    abs.min <- abs(min(d[,2]))
		
    g1 <- ggplot(data=d,aes_string(x="Predictor",y=ylb,group="Predictor")) + 
          geom_bar(stat="identity", colour="#639f89", fill="#639f89") + theme_grey(base_size=fontsize)

    if(ylb=="mda")      g1 <- g1 + labs(y="Mean decrease in accuracy") 
    else if(ylb=="mdg") g1 <- g1 + labs(y="Mean decrease in Gini")

    g1 <- g1 + theme(axis.title=element_text(size=25,face="bold"), 
                     axis.text.x=element_text(angle=0,hjust=1,vjust=0.4,colour='black'),
                     axis.text.y= element_text(colour='black', size=25)) + 
          geom_hline(yintercept=abs.min,linetype="dashed",colour="black") + coord_flip()
    print(g1)
	}

  importancePlot(d=pred.acc, ylb="mda", 20)
  importancePlot(d=pred.gini, ylb="mdg", 20)
}


plotRFVarImp2 <- function(rf.mod){
# Plot variable importance of a RF model 
#
# Args:
#   rf.mod: a rf model obj
#
# Returns:
#   Two side-by-side Plots: 1. based on pred accuracy 2. based on gini index
  
  varImpPlot(rf.mod)
  #varImp <- sort(rf.mod$importance[,4], dec=T)
  #plot(varImp, type="h")
}



#------------------------------------------------------------------------------------------------------
# Plot OOB error rate at different number of trees
#------------------------------------------------------------------------------------------------------
plotRFOOBErr <- function(rf.mod){
# Plot variable importance of a RF model 
#
# Args:
#   rf.mod: a rf model obj
#
# Returns:
#   OOB error rate at different number of trees

  # OOB error rate at differnt number of trees
  oob.err <- as.data.frame(cbind(No.Trees=1:nrow(rf.mod$err.rate), rf.mod$err.rate))
  names(oob.err)[2] <- "Classification error rate"
  names(oob.err)[3] <- "~Q1 error rate"
  names(oob.err)[4] <- "Q1 error rate"
  dat <- melt(oob.err, id="No.Trees")

  g1 <- ggplot(data=dat, aes(x=No.Trees, y=value, colour=variable)) + #+ xlim(0,1000) 
        geom_line(size=1.1) + geom_point(size=1) + xlab("Number of trees") + ylab("Error rate") +
        theme(axis.title.x = element_text(size=24),
              axis.title.y = element_text(size=24),
              axis.text.x = element_text(colour="grey20",size=15),
              axis.text.y = element_text(colour="grey20",size=15),
              legend.title=element_blank(),
              legend.text = element_text(size = 20),
              legend.justification=c(1,0), legend.position=c(1,0.9),
              legend.background = element_rect(fill="gray90", size=.5, linetype="dotted"))
  print(g1)
}  
  


#------------------------------------------------------------------------------------------------------
# Plot classification accuracy at different training % (avg over nrep runs)
#------------------------------------------------------------------------------------------------------
plotRFAcc <- function(sol.all){
  # Plot classification accuracy at different training % of a RF model 
  #
  # Args:
  #   sol: a matrix of summarized results of rf model at different training %
  #        each row represent an avg results over nrep for a fix training %
  #        vars for each row: m, n.Tree, Train.Perc, Accy, TP, TN
  # Returns:
  #   OOB error rate at different number of trees
  
  # Plot accuracy at different training % (Avg 50 runs results)
  err <- as.data.frame(sol.all[,seq(3,6)])
  names(err)[2] <- "Classification Accuracy"
  names(err)[3] <- "True Positive Rate"
  names(err)[4] <- "True Negative Rate"
  dat <- melt(err, id="Train.Perc")
  dat.accy <- dat[1:9,]  # classificication acc
  
  g1 <- ggplot(data=dat, aes(x=Train.Perc, y=value, colour=variable)) + 
        geom_line(size=1.1) + geom_point(size=4) + 
        ylim(0.5,1) + xlim(0.1,0.9) + scale_x_continuous(breaks=seq(0.1,0.9,0.1)) +
        xlab("Percentage of Training Data") + ylab("Test Classification Accuracy") + 
        theme(axis.title.x = element_text(size=24),
              axis.title.y = element_text(size=24),
              axis.text.x = element_text(colour="grey20",size=15),
              axis.text.y = element_text(colour="grey20",size=15),
              legend.title=element_blank(),
              legend.text = element_text(size = 20),
              legend.justification=c(1,0), legend.position=c(1,0),
              legend.background = element_rect(fill="gray90", size=.5, linetype="dotted"))
  
  print(g1) 
}



# Plot accuracy at a fix training % (e.g. 70%)
# dat.train.70pct <- read.csv("./train_perc_0.7_errrate.csv", header=T)
#   dat2 <- dat.train.70pct[,-2]
#   names(dat2)[2] <- "Classification Accuracy"
#   names(dat2)[3] <- "True Positive Rate"
#   names(dat2)[4] <- "Ture Negative Rate"
#   dat2 <- melt(dat2, id="Index")
# 
#   ggplot(data=dat2, aes(x=Index, y=value, colour=variable)) + 
#   geom_line(size=1.1) + geom_point(size=4) + 
#   ylim(0.3,1) + xlim(1,30) + scale_x_continuous(breaks=seq(1,30,1)) +
#   xlab("Runs for 70% Training Dataset") + ylab("Test Classification Accuracy") + 
#   theme(
#     axis.title.x = element_text(size=24),
#     axis.title.y = element_text(size=24),
#     axis.text.x = element_text(colour="grey20",size=15),
#     axis.text.y = element_text(colour="grey20",size=15),
#     legend.title=element_blank(),
#     legend.text = element_text(size = 20),
#     legend.justification=c(1,0), legend.position=c(1,0),
#     legend.background = element_rect(fill="gray90", size=.5, linetype="dotted")
#   )


#------------------------------------------------------------------------------------------------------
# For saving plots
#------------------------------------------------------------------------------------------------------
#pdf(file = file, width=11, height=8.5)
#dev.off()

#png(file = file, width=11, height=8.5)
#dev.off()

