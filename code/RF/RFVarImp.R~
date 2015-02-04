## Plot vars importance for RF model

plotRFVarImp <- function(rf.mod){
# Plot variable importance of a RF model 
#
# Args:
#   rf.mod: a rf model obj
#
# Returns:
#   2 Plots: 1. based on pred accuracy 2. based on gini index 

	# Importance data
	dat <- data.frame(rownames(importance(rf.mod)),round(importance(rf.mod),2))
	# measure 1:mean decrease in accuracy  2:mean decrease in gini
	names(dat)[c(1,ncol(dat)-1,ncol(dat))] <- c("Predictor","mda","mdg")  
	rownames(dat) <- NULL

	pred.acc  <- select(dat, Predictor, mda)
	pred.gini <- select(dat, Predictor, mdg)

	# Var importance plot function
	importancePlot <- function(d,ylb,fontsize){
		fontsize <- as.numeric(fontsize)
		d <- d[order(d[,2],decreasing=T),]
		d$Predictor <- factor(as.character(d$Predictor),levels=rev(as.character(d$Predictor)))
		rownames(d) <- NULL
  
    d[,2] <- d[,2]/abs(max(d[,2])) * 100  # normalize relative to the maximum
		abs.min <- abs(min(d[,2]))
		
    g1 <- ggplot(data=d,aes_string(x="Predictor",y=ylb,group="Predictor")) + geom_bar(stat="identity", colour="#639f89", fill="#639f89") + theme_grey(base_size=fontsize)
		if(ylb=="mda") g1 <- g1 + labs(y="Mean decrease in accuracy") else if(ylb=="mdg") g1 <- g1 + labs(y="Mean decrease in Gini")
		g1 <- g1 + theme(axis.title=element_text(size=25,face="bold"), 
                     axis.text.x=element_text(angle=0,hjust=1,vjust=0.4,colour='black'),
                     axis.text.y= element_text(colour='black', size=25)) + 
          geom_hline(yintercept=abs.min,linetype="dashed",colour="black") + coord_flip()
		print(g1)
	}

	importancePlot(d=pred.acc, ylb="mda", 20)
	importancePlot(d=pred.gini, ylb="mdg", 20)
}
