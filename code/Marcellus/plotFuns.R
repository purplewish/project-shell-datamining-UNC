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
# Multiple Line Plot 
#------------------------------------------------------------------------------------------------------
plotMLine <- function(dat, xlab, ylab, xlim, ylim){
  # Plot multiple lines 
  #
  # Args:
  #   dat: a data frame, x=1st col, y=2nd col
  #
  # Returns:
  #   x-y line plot
  a <- names(dat); 
  g <- ggplot(data=dat, aes_string(x=a[1],y=a[2], colour=a[3], group=a[3])) 
  
  if(length(xlim)!=0){
      g <- g + xlim(xlim[1], xlim[2]) 
  }
  if(length(ylim)!=0) {
      g <- g + ylim(ylim[1], ylim[2])
  }
  
    #geom_line(size=1, colour="#000000") + geom_point(size=4) + scale_x_continuous(breaks=seq(1,31,1)) + #ylim(0.5,1) +
  g <- g + geom_line(lwd=1.8) + #geom_point(size=4) + 
    #scale_color_manual(values=c("red", "#007cd2","orange", "#5E9F37", "black")) +
    #scale_color_manual(values=c("red", "#007cd2", "#5E9F37", "black")) +
    scale_color_manual(values=c("#cd0000", "#ffcc00")) +
    xlab(xlab) + ylab(ylab) + 
    theme(#legend.position="none",
      axis.title.x = element_text(size=26, face='bold'),
      axis.title.y = element_text(size=26, face='bold'),
      #axis.text.x = element_text(colour="grey20",size=15),
      axis.text.x = element_text(angle=0, hjust=0.5, colour='black', size=25, face='bold'),
      axis.text.y = element_text(colour="black",size=25, face='bold'),
      legend.title=element_blank(),
      legend.text = element_text(size = 26, face='bold'),
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
  myPaletter <- colorRampPalette(rev(brewer.pal(11, "Spectral")), space="Lab")
  
  g <- ggplot(dat=dat, aes(x=Longitude, y=Latitude, colour=Production)) + geom_point(size=3, alpha=0.7) + 
       scale_colour_gradientn(colours = myPaletter(100), limits=c(0, 60)) +
       #ylim(27.5,31.4) + xlim(-100.5,-96) +
       theme(
          legend.text = element_text(size = 16),
          legend.title = element_text(size = 20),
          axis.title.x = element_text(size=28),
          axis.title.y = element_text(size=28),
          legend.justification=c(1,0), legend.position=c(1,0)
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
  g <- ggplot(dat=dat, aes(x=Longitude, y=Latitude, colour=ID, size=ID)) + geom_point(alpha=0.7) + 
    #scale_color_manual(name="", labels = c("Cored wells","Production wells"), values=c("red", "grey30")) +
    #scale_size_manual(name="", labels= c("Cored wells","Production wells"), values=c(5, 3)) + 
    scale_color_manual(labels = c("83 Cored wells","7200 Producers"), values=c("red", "grey30")) +
    scale_size_manual(labels= c("83 Cored wells","7200 Producers"), values=c(5, 3)) + 
    ylim(27.5,31.4) + xlim(-100.5,-96) +
    theme(
      legend.text = element_text(size = 34, face='bold'),
      #legend.title = element_text(size = 28),
      legend.title = element_blank(),
      axis.title.x = element_text(size=34),
      axis.title.y = element_text(size=34),
      axis.text.x = element_text(size=26),
      axis.text.y = element_text(size=26),
      legend.key.height=unit(3,"line"),
      #legend.justification=c(1,0), legend.position=c(0.95,0)
      legend.justification=c(1,0), legend.position=c(0.33,0.73)
    ) 
  
  print(g) 
}


#------------------------------------------------------------------------------------------------------
# Heatmap of production with county boundary as backgroud
#------------------------------------------------------------------------------------------------------
grid4HeatmapProd <- function(d){
  # Generate grid points for heatmap of well production within the boundary
  #
  # Args:
  #   d: dataframe use to define heatmap boundary with 2 cols: Longitude Latitude
  #
  # Returns:
  #    Grid points dataframe (Longitude, Latitude) within the boundary 

  # Define boundary through d
  aq.ch <- chull(d$Longitude, d$Latitude)
  aq.ch <- c(aq.ch, aq.ch[1])
  aq.border <- cbind(d$Longitude[aq.ch], d$Latitude[aq.ch])
  
  # Generate grid points inside boundary
  aq.bbox <- sbox(as.points(d$Longitude, d$Latitude))
  aq.grid <- gridpts(aq.bbox, npts=100000)  # number of point in larger box
  inside <- inout(aq.grid, aq.border, bound=TRUE)
  aq.Grid <- aq.grid[inside,]
  
  grid <- as.data.frame(aq.Grid)
  names(grid)<-c('Longitude', 'Latitude')
  
  return(grid)
}

# Nearest Neibor function
idw <- function(z,distance,k,num.neighs)
{
  idw.z<-rep(0,length(distance[,1]))
  for (i in 1:length(distance[,1]))
  {
    d<-sort(distance[i,],index.return=TRUE)
    w<-1/d$x[1:num.neighs]^k
    idw.z[i]<-sum(z[d$ix[1:num.neighs]]*w)/sum(w)
  }
  return(idw.z)
}  




plotHeatmapProd <- function(dat, train_well_loc, long.range, lat.range, time){
  # Plot heatmap of well production
  #
  # Args:
  #   d: dataframe use to define heatmap boundary with 2 cols: Longitude Latitude
  #   dat: dataframe with 3 cols:  Longitude Latitude Production(prediction)
  #   long.range: vector with long.min and long.max
  #   lat.range: vector with lat.min and lat.max
  #
  # Returns:
  #    Heatmap of well production with county boundary as backgroud
  
  
  jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
  myPaletter <- colorRampPalette(rev(brewer.pal(11, "Spectral")), space="Lab")
  
  # Define county boundary based on long.range and lat.range
  all_states <- map_data("county")
  #county <- all_states %>% filter(long>long.range[1] & long<long.range[2] & lat>lat.range[1] & lat<lat.range[2])
  #a<-county$subregion
  #county <- subset(all_states,(subregion %in% a)&(region=='texas'))
  b <- c('pennsylvania', 'ohio', 'west virginia', 'maryland', 'virginia', 'kentucky')
  #county <- subset(all_states,(subregion %in% a)&(region %in% b))
  county <- subset(all_states,(region %in% b))

  
  # Heatmap
  g <- ggplot() + geom_polygon(data=county, aes(x=long, y=lat, group=group), colour="#696969", fill="white", size=1) +
       geom_tile(data=dat, aes(x=Longitude, y=Latitude, z=Production, fill=Production), alpha = 0.8) + 
       scale_fill_gradientn(colours = myPaletter(100), limits=c(0, 1000)) + #scale_alpha_continuous(range=c(0,0.5)) + #scale_y_continuous(expand=c(0,0)) + scale_x_continuous(expand=c(0,0))+
       annotate(geom="text", x=long.range[1]+2, y=lat.range[2]+0.5, label=time, color= "black", size= 18) +
       #geom_point(data=train_well_loc, aes(x=Longitude, y=Latitude), shape=1, color="#191919", size=2.5) + # add well locations
       theme_bw() +  theme(line = element_blank(),
                           panel.border = element_blank(),
                           axis.text.x = element_blank(),
                           axis.text.y = element_blank(),
                           #text = element_blank(),
                           legend.key.height = unit(2.8, "cm"),
                           legend.text = element_text(size=26,face = 'bold'),
                           line = element_blank(),
                           title = element_blank()) #+
       #stat_contour(data=dat, aes(x=Longitude, y=Latitude, z=Production))
  
  print(g) 
}



plotPredvsAct <- function(dat, xlim=NULL, ylim=NULL, title=NULL){
  ## Plot cross-plot (Predict vs Actual)
  #
  # Args:
  #   dat: a data frame: predict, actual
  #   xlim: x axis limit e.g. c(0,20)
  #   ylim: y axis limit e.g. c(0,50)
  #   title: a string, plot title
  #
  # Returns:
  #   Cross-plot (Predict vs Actual)
  a <- names(dat) 
  g <- ggplot(data=dat, aes_string(x=a[2],y=a[1])) + geom_point(alpha=0.7, size=2.5) +
    geom_abline(intercept=0, size=1, colour='red') + coord_fixed(ratio=1) +
    xlab('Observed') + ylab('Predicted') +
    theme(
      axis.title.x = element_text(size=28),
      axis.title.y = element_text(size=28),
      axis.text.x = element_text(colour="grey20",size=24),
      axis.text.y = element_text(colour="grey20",size=24),
      plot.title = element_text(lineheight=.8, face="bold", size=32, vjust=2)
    ) 
  
  if(length(xlim)!=0) { g <- g + xlim(xlim[1], xlim[2]) }
  if(length(ylim)!=0) { g <- g + ylim(ylim[1], ylim[2]) }
  if(length(title)!=0) { g <- g + ggtitle(title) }
  
  print(g)
}


plotBar <- function(dat, ylim=NULL, xlab=NULL, ylab=NULL, title=NULL){
  ## Plot bar chart
  #
  # Args:
  #   dat: a data frame: category(factor), value(numerical)
  #   ylim: y axis limit e.g. c(0,50)
  #   xlab: a string, x axis lab
  #   ylab: a string, y axis lab
  #   title: a string, plot title
  #
  # Returns:
  #   Bar chart
  googPalette <- c("#008744", "#0057e7", "#d62d20", "#ffa700", "#ffa700")
  
  a <- names(dat) 
  g <- ggplot(data=dat, aes_string(x=a[1],y=a[2],fill=a[1])) +
       geom_bar(stat = "identity") +
       scale_fill_manual(values=googPalette) + 
       theme(
              legend.position="none",
              axis.title.x = element_text(size=28, face="bold", vjust=-1),
              axis.title.y = element_text(size=28, face="bold", vjust=1),
              axis.text.x = element_blank(), #element_text(colour="black",size=32, angle=0),
              axis.text.y = element_text(colour="black",size=32),
              plot.title = element_text(lineheight=.8, face="bold", size=32, vjust=2)
      ) 

  if(length(ylim)!=0) { g <- g + ylim(ylim[1], ylim[2]) }
  if(length(xlab)!=0) { g <- g + xlab(xlab) }
  if(length(ylab)!=0) { g <- g + ylab(ylab) }
  
  print(g)
}
