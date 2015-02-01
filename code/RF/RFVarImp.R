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
	d <- data.frame(rownames(importance(rf.mod)),round(importance(rf.mod),2))
	# measure 1:mean decrease in accuracy  2:mean decrease in gini
	names(d)[c(1,ncol(d)-1,ncol(d))] <- c("Predictor","mda","mdg")  
	rownames(d) <- NULL

	pred.acc  <- select(d, Predictor, mda)
	pred.gini <- select(d, Predictor, mdg)

	# Var importance plot function
	importancePlot <- function(d,ylb,fontsize){
		fontsize <- as.numeric(fontsize)
		d <- d[order(d[,2],decreasing=T),]
		d$Predictor <- factor(as.character(d$Predictor),levels=rev(as.character(d$Predictor)))
		rownames(d) <- NULL
		abs.min <- abs(min(d[,2]))
		g1 <- ggplot(data=d,aes_string(x="Predictor",y=ylb,group="Predictor",colour="Predictor",fill="Predictor")) + geom_bar(stat="identity") + theme_grey(base_size=fontsize)
		#g1 <- ggplot(data=d,aes_string(x="Predictor",y=ylb,group="Predictor",color="black",fill="black")) + geom_bar(stat="identity") + theme_grey(base_size=fontsize)
		if(ylb=="mda") g1 <- g1 + labs(y="Mean decrease in accuracy") else if(ylb=="mdg") g1 <- g1 + labs(y="Mean decrease in Gini")
		g1 <- g1 + theme(axis.text.x = element_text(angle=90,hjust=1,vjust=0.4)) + geom_hline(yintercept=abs.min,linetype="dashed",colour="black") + coord_flip()
		print(g1)
	}

	importancePlot(d=pred.acc, ylb="mda", 20)
	importancePlot(d=pred.gini, ylb="mdg", 20)
}
