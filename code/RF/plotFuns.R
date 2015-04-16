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
  dat[,1] <- gsub("(Core.|.Kriged|.Joined)",'',dat[,1])  # strip unneccsary characters
    
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
          geom_bar(stat="identity", colour="#d62d20", fill="#d62d20") + theme_grey(base_size=fontsize)
          #geom_bar(stat="identity", colour="#639f89", fill="#639f89") + theme_grey(base_size=fontsize)

    #if(ylb=="mda")      g1 <- g1 + labs(y="Mean decrease in accuracy") 
    #else if(ylb=="mdg") g1 <- g1 + labs(y="Mean decrease in Gini")
    g1 <- g1 + labs(y="Variable Importance") # Simplify for presentation purpose
    

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



plotRFVarImp3 <- function(dat){
  # Plot variable importance of a RF model 
  #
  # Args:
  #   dat: data.frame with predictors and importance scores
  # Returns:
  #   Normalized variable importance plot

  fontsize <- 20
  d <- dat[order(dat[,2],decreasing=T),]
  d[,1] <- factor(as.character(d[,1]),levels=rev(as.character(d[,1])))
  rownames(d) <- NULL
    
  d[,2] <- d[,2]/abs(max(d[,2])) * 100  # normalize relative to the variable with maximum score
  abs.min <- abs(min(d[,2]))
  
  xlb <- names(d)[1]; ylb <- names(d)[2];
  g1 <- ggplot(data=d,aes_string(x=xlb,y=ylb,group=xlb)) + 
        geom_bar(stat="identity", colour="#d62d20", fill="#d62d20") + theme_grey(base_size=fontsize)
        #geom_bar(stat="identity", colour="#639f89", fill="#639f89") + theme_grey(base_size=fontsize)
    
  #if(ylb=="mda")      g1 <- g1 + labs(y="Mean decrease in accuracy") 
  #else if(ylb=="mdg") g1 <- g1 + labs(y="Mean decrease in Gini")
  g1 <- g1 + labs(y="Variable Importance") # Simplify for presentation purpose
    
  g1 <- g1 + theme(axis.title=element_text(size=25,face="bold"), 
                   axis.text.x=element_text(angle=0,hjust=1,vjust=0.4,colour='black'),
                   axis.text.y= element_text(colour='black', size=25)) + 
        geom_hline(yintercept=abs.min,linetype="dashed",colour="black") + coord_flip()
  
  print(g1)
  
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
              legend.justification=c(1,0), legend.position=c(1,0.86),
              legend.background = element_rect(fill="gray90", size=.5, linetype="dotted"))
  g1 <- g1 + geom_vline(xintercept = 1000, linetype="longdash")
  print(g1)
}  
  


#------------------------------------------------------------------------------------------------------
# Plot classification accuracy at different training % (avg over nrep runs)
#------------------------------------------------------------------------------------------------------
plotRFAcc <- function(err, xlb, ylb, xtick){
  # Plot classification accuracy at different training % of a RF model 
  #
  # Args:
  #   err: a matrix of summarized accuracy results of rf model at different x
  #        each row represent an avg results over nrep for a fix x
  #        vars for each row: x, Accy, TP, TN
  #   xlb: x axis label
  #   ylb: y axis label
  #   xtick: x ticks 
  # Returns:
  #   Classification Accuracy plot
  
  # Plot accuracy at different training %
  names(err)[2] <- "Classification Accuracy"
  names(err)[3] <- "True Positive Rate"
  names(err)[4] <- "True Negative Rate"
  dat <- melt(err, id=names(err)[1])
  dat.accy <- dat[1:9,]  # classificication acc
  
  g1 <- ggplot(data=dat, aes_string(x=names(err)[1], y="value", colour="variable")) + 
        geom_line(size=1.1) + geom_point(size=4) + 
        ylim(0.4,1) + scale_x_continuous(breaks=xtick) +
        xlab(xlb) + ylab(ylb) + 
        theme(axis.title.x = element_text(size=24),
              axis.title.y = element_text(size=24),
              axis.text.x = element_text(colour="grey20",size=15),
              axis.text.y = element_text(colour="grey20",size=15),
              legend.title = element_blank(),
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
# Plot sweet-spots
#------------------------------------------------------------------------------------------------------
# plotSweetspot <- function(rf.mod, target){
#   # Plot Sweetspot
#   #
#   # Args:
#   #   rf.mod: a rf object that trained through data all
#   #   target: data set with 1) true top quartile indicator/target 2) Longitude and Latitude
#   #
#   # Returns:
#   #   Sweetspot plot with true+, false-, false+ 
#   
#   # Predicted results
#   a <- rf.mod$predicted
#   b <- as.numeric(names(a))
#   c <- data.frame(b, a)
#   d <- arrange(c, b)  # order predicted value from RF model
#   dat <- cbind(target, pred=d[,2])  # Uwi Target.Q4 Latitude Longitude pred
#   
#   overlap <- dat %>% filter(Target.Q4==pred & pred==TRUE) %>% select(Uwi)  # true+
#   top.q <- dat %>% 
#             filter(Target.Q4==TRUE | pred==TRUE) %>% 
#             mutate(type=ifelse(Target.Q4==pred, "Overlap", ifelse(pred==TRUE,"Pred","True")))
#   
#   g <- ggplot(top.q, aes(x=Longitude, y=Latitude, colour=type)) + geom_point(size=4) +  
#         scale_colour_manual(values=c("#f05032", "#00539e", "#fec325")) +
#         theme(
#               #legend.title=element_blank(),
#               #legend.text = element_text(size = 20),
#               #legend.position=c(.1, .9),
#               #legend.position="none",
#               #panel.background = element_blank(),
#               axis.title.x = element_text(size=28),
#               axis.title.y = element_text(size=28)
#         )
#   
#   print(g) 
# }


#------------------------------------------------------------------------------------------------------
# Plot sweet-spots
#------------------------------------------------------------------------------------------------------
plotSweetspot <- function(dat){
  # Plot Sweetspot
  #
  # Args:
  #   dat: dataframe with cols: Uwi Target.Q4 Latitude Longitude pred
  #
  # Returns:
  #   Sweetspot plot with true+, false-, false+ 
  
  overlap <- dat %>% filter(Target.Q4==Target.Q4.Pred & Target.Q4.Pred==TRUE) %>% select(Uwi)  # true+
  top.q <- dat %>% 
    filter(Target.Q4==TRUE | Target.Q4.Pred==TRUE) %>% 
    mutate(type=ifelse(Target.Q4==Target.Q4.Pred, "True Positive", ifelse(Target.Q4.Pred==TRUE,"False Positive","False Negative")))  # Pred:false+  True:fasle-
  
  g <- ggplot(top.q, aes(x=Longitude, y=Latitude, colour=type)) + geom_point(size=4) +  
    #scale_colour_manual(values=c("#f05032", "#00539e", "#")) +
    scale_colour_manual(values=c("#fec325", "#00539e", "#f05032")) +
    theme(
      legend.title=element_blank(),
      legend.text = element_text(size = 32),
      legend.key.height=unit(3,"line"),
      legend.position=c(.12, .6),
      #legend.position="none",
      #panel.background = element_blank(),
      axis.title.x = element_text(size=28),
      axis.title.y = element_text(size=28)
    )
  
  print(g) 
}





#------------------------------------------------------------------------------------------------------
# Line Plot 
#------------------------------------------------------------------------------------------------------
plotLine <- function(dat, xlab, ylab){
  # Plot x-y line 
  #
  # Args:
  #   dat: a data frame, x=1st col, y=2nd col
  #
  # Returns:
  #   x-y line plot
  a <- names(dat); 
  g <- ggplot(data=dat, aes_string(x=a[1],y=a[2])) + 
        geom_line(size=1, colour="#000000") + geom_point(size=4) + scale_x_continuous(breaks=seq(1,31,1)) + #ylim(0.5,1) +
        xlab(xlab) + ylab(ylab) + 
        theme(legend.position="none",
              axis.title.x = element_text(size=24),
              axis.title.y = element_text(size=24),
              axis.text.x = element_text(colour="grey20",size=15),
              axis.text.y = element_text(colour="grey20",size=15),
              legend.title=element_blank(),
              legend.text = element_text(size = 20),
              legend.justification=c(1,0), legend.position=c(1,0),
              legend.background = element_rect(fill="gray90", size=.5, linetype="dotted")
              )
  
  print(g)
}



#------------------------------------------------------------------------------------------------------
# Plot correlation heat map
#------------------------------------------------------------------------------------------------------
plotCorr <- function(x){
  # Plot x-y line 
  #
  # Args:
  #   x: variables data frame: x1 x2 x3...
  #
  # Returns:
  #   x-x pearson correlation heat map
  
  g <- qplot(x=Var1, y=Var2, data=melt(cor(x, use="p")), fill=value, geom="tile") + 
        scale_fill_gradient2(limits=c(-1,1)) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
        theme(axis.text.x=element_text(colour="black")) +
        theme(axis.text.y=element_text(colour="black"))
  
  print(g)
}


#------------------------------------------------------------------------------------------------------
# For saving plots
#------------------------------------------------------------------------------------------------------
#pdf(file = file, width=11, height=8.5)
#dev.off()

#png(file = file, width=11, height=8.5)
#dev.off()


#------------------------------------------------------------------------------------------------------
# Plot Well Production
#------------------------------------------------------------------------------------------------------
plotWellProd <- function(dat){
  # Plot Well production
  #
  # Args:
  #   dat: dataframe with cols:  Longitude Latitude Production
  #
  # Returns:
  #   Production well plot with production heatmap
  g <- ggplot(dat=dat, aes(x=Longitude, y=Latitude, colour=Production)) + geom_point(size=3) + 
    scale_colour_gradient(low="#56B4E9",  high="red") + 
    theme(
      legend.text = element_text(size = 16),
      legend.title = element_text(size = 20),
      axis.title.x = element_text(size=28),
      axis.title.y = element_text(size=28)
    ) 
  
  print(g) 
}


#------------------------------------------------------------------------------------------------------
# Plot Core + Production Well
#------------------------------------------------------------------------------------------------------
plotCoreProd <- function(dat){
  # Plot Well production
  #
  # Args:
  #   dat: dataframe with cols:  Longitude Latitude ID (core/prod)
  #
  # Returns:
  #    Plot Core + Production Well
  g <- ggplot(dat=dat, aes(x=Longitude, y=Latitude, colour=ID, size=ID)) + geom_point() + 
    scale_color_manual(name="Dataset", labels = c("Core wells","Production wells"), values=c("red", "grey30")) +
    scale_size_manual(name="Dataset", labels= c("Core wells","Production wells"), values=c(5, 3)) + 
    theme(
      legend.text = element_text(size = 16),
      legend.title = element_text(size = 20),
      axis.title.x = element_text(size=28),
      axis.title.y = element_text(size=28)
    ) 
  
  print(g) 
}

