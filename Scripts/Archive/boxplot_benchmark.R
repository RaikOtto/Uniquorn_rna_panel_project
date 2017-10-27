### boxplot benchmark

source("~/Dropbox/PhD/Uniquorn_project/Scripts/utility.R")
require(ggplot2)
require(grid)
require(gridExtra)
library("cowplot")
library("scales")

# True Positives

benchmark_TP = data.frame(
  "Weight" = factor(
    c( "1.0","0.5","0.25","0.0","1.0","0.5","0.25","0.0" )
  ),
  "Regularization" = factor(
    c( rep("Yes",4), rep("No",4) )
  ),
  "TP" = c( 3023,3462,3451,3104,3368,3505,3512,3471 )
)

g_TP = ggplot( benchmark_TP, aes( x = Weight, y=TP, group=Regularization))
g_TP = g_TP + geom_line(aes(colour = Regularization),size=2)
g_TP = g_TP + theme(legend.position="top", legend.background = element_rect(fill="gray90", size =5 )  )
#y_pos= min( benchmark_TP[ , 3] ) + ( max( benchmark_TP[ , 3] ) - min( benchmark_TP[ , 3] ) ) / 2
#g_TP = g_TP + annotate("text", x = 2.5, y = y_pos, label = "A", size = 10)
g_TP = g_TP 

# False Positives

benchmark_FP = data.frame(
  "Weight" = factor(
    c( "1.0","0.5","0.25","0.0","1.0","0.5","0.25","0.0" )
  ),
  "Regularization" = factor(
    c( rep("Yes",4), rep("No",4) )
  ),
  "FP" = c(  29,53,73,451,25,114,175,7709 )
)


g_FP = ggplot( benchmark_FP, aes( x = Weight, y=FP, group=Regularization))
g_FP = g_FP + geom_line(aes(colour = Regularization),size=2)
g_FP = g_FP + coord_cartesian(ylim = c(0, 600))
g_FP = g_FP + theme(legend.position="top", legend.background = element_rect(fill="gray90", size =5 )  )
#y_pos= mean( g_FP$coordinates$limits$y )
#g_FP = g_FP + annotate("text", x = 2.5, y = y_pos, label = "B", size = 10)
g_FP

# True Negatives

benchmark_TN = data.frame(
  "Weight" = factor(
    c( "1.0","0.5","0.25","0.0","1.0","0.5","0.25","0.0" )
  ),
  "Regularization" = factor(
    c( rep("Yes",4), rep("No",4) )
  ),
  "TN" = 1989**2 - benchmark_FP$FP
)

g_TN = ggplot( benchmark_TN, aes( x = Weight, y=TN, group=Regularization))
g_TN = g_TN + geom_line(aes(colour = Regularization),size=2)
g_TN = g_TN + theme(legend.position="top", legend.background = element_rect(fill="gray90", size =5 )  )
#y_pos= min( benchmark_TN[ , 3] ) + ( max( benchmark_TN[ , 3] ) - min( benchmark_TN[ , 3] ) ) / 2
#g_TN = g_TN + annotate("text", x = 2.5, y = y_pos, label = "F", size = 10)
#g_TN = g_TN + scale_y_continuous(labels = scientific )
g_TN

# False Negatives

benchmark_FN = data.frame(
    "Weight" = factor(
        c( "1.0","0.5","0.25","0.0","1.0","0.5","0.25","0.0" )
    ),
    "Regularization" = factor(
        c( rep("Yes",4), rep("No",4) )
    ),
    "FN" = c(  532,93,104,451,187,50,43,84)
)

g_FN = ggplot( benchmark_FN, aes( x = Weight, y=FN, group=Regularization))
g_FN = g_FN + geom_line(aes(colour = Regularization),size=2)
g_FN = g_FN + theme(legend.position="none" )  #+ theme(legend.position="top", legend.background = element_rect(fill="gray90", size =5 )  )
#y_pos= min( benchmark_FN[ , 3] ) + ( max( benchmark_FN[ , 3] ) - min( benchmark_FN[ , 3] ) ) / 2
#g_FN = g_FN + annotate("text", x = 2.5, y = y_pos, label = "E", size = 10)
g_FN

# sensitivity

benchmark_sens = data.frame(
  "Weight" = factor(
    c( "1.0","0.5","0.25","0.0","1.0","0.5","0.25","0.0" )
  ),
  "Regularization" = factor(
    c( rep("Yes",4), rep("No",4) )
  ),
  "Sensitivity" = c( .8504,.9738,.9707,.8731,.9474,.9859,.9879,.9764 )
)

g_sens = ggplot( benchmark_sens, aes( x = Weight, y=Sensitivity, group=Regularization))
g_sens = g_sens+ geom_line(aes(colour = Regularization),size=2)
g_sens = g_sens + theme(legend.position="none")#+ theme(legend.position="top", legend.background = element_rect(fill="gray90", size =5 )  )
#y_pos= min( benchmark_sens[ , 3] ) + ( max( benchmark_sens[ , 3] ) - min( benchmark_sens[ , 3] ) ) / 2
#g_sens = g_sens + annotate("text", x = 2.5, y = y_pos, label = "C", size = 10)
g_sens

# specificity

benchmark_spec = data.frame(
  "Weight" = factor(
    c( "1.0","0.5","0.25","0.0","1.0","0.5","0.25","0.0" )
  ),
  "Regularization" = factor(
    c( rep("Yes",4), rep("No",4) )
  ),
  "Specificity" = rep(round( as.integer(99),0 ), 8 )
)

g_spec = ggplot( benchmark_spec, aes( x = Weight, y=Specificity, group=Regularization))
g_spec = g_spec+ geom_line(aes(colour = Regularization),size=2)
g_spec = g_spec+ theme(legend.position="none")# + theme(legend.position="top", legend.background = element_rect(fill="gray90", size =5 )  )
#y_pos= min( benchmark_spec[ , 3] ) + ( max( benchmark_spec[ , 3] ) - min( benchmark_spec[ , 3] ) ) / 2
#g_spec = g_spec + annotate("text", x = 2.5, y = y_pos, label = "D", size = 10)
g_spec + scale_y_continuous( )

# F1

benchmark_F1 = data.frame(
  "Weight" = factor(
    c( "1.0","0.5","0.25","0.0","1.0","0.5","0.25","0.0" )
  ),
  "Regularization" = factor(
    c( rep("Yes",4), rep("No",4) )
  ),
  "F1" = c( .92,.98,.97,.55,.97,.98,.97,.47 )
)

g_F1 = ggplot( benchmark_F1, aes( x = Weight, y=F1, group=Regularization))
g_F1 = g_F1+ geom_line(aes(colour = Regularization),size=2)
g_F1 = g_F1 + theme(legend.position="none")#+ theme(legend.position="top", legend.background = element_rect(fill="gray90", size =5 )  )
#y_pos= min( benchmark_F1[ , 3] ) + ( max( benchmark_F1[ , 3] ) - min( benchmark_F1[ , 3] ) ) / 2
#g_F1 = g_F1 + annotate("text", x = 2.5, y = y_pos, label = "G", size = 10)
g_F1

# PPV

benchmark_PPV = data.frame(
  "Weight" = factor(
    c( "1.0","0.5","0.25","0.0","1.0","0.5","0.25","0.0" )
  ),
  "Regularization" = factor(
    c( rep("Yes",4), rep("No",4) )
  ),
  "PPV" = c( .99,.98,.98,.4,.99,.97,.95,.31 )
)

g_PPV = ggplot( benchmark_PPV, aes( x = Weight, y=PPV, group=Regularization))
g_PPV = g_PPV+ geom_line(aes(colour = Regularization),size=2)
g_PPV = g_PPV + theme(legend.position="none") #+ theme(legend.position="top", legend.background = element_rect(fill="gray90", size =5 )  )
#y_pos= min( benchmark_PPV[ , 3] ) + ( max( benchmark_PPV[ , 3] ) - min( benchmark_PPV[ , 3] ) ) / 2
#g_PPV = g_PPV + annotate("text", x = 2.5, y = y_pos, label = "H", size = 10)
g_PPV 
+ theme(axis.title.x = element_text(face="bold", colour="#990000", size=20),
        axis.text.x  = element_text(angle=90, vjust=0.5, size=16))


prow = plot_grid(
    g_TP   + theme(legend.position="none") + theme( axis.text.y = element_text(angle=45) ),
    g_FN   + theme(legend.position="none") + theme( axis.text.y = element_text(angle=45) ),
    g_FP   + theme(legend.position="none") + theme( axis.text.y = element_text(angle=45) ),
    g_TN   + theme(legend.position="none")  + theme( axis.text.y = element_text(angle=45) ) + scale_y_continuous(labels = c(3.7,3.8,3.9,4.0,4.1) ),
    g_sens + theme(legend.position="none") + theme( axis.text.y = element_text(angle=45) ),
    g_F1   + theme(legend.position="none") + theme( axis.text.y = element_text(angle=45) ),
    g_spec + theme(legend.position="none") + theme( axis.text.y = element_text(angle=45) ) + scale_y_continuous(labels = c(97:100,"") ),
    g_PPV  + theme(legend.position="none") + theme( axis.text.y = element_text(angle=45) ), 
    labels=c("A", "B","C","D","E","F","G","H"), 
    ncol = 4,
    nrow = 2
)

legend_b = get_legend( g_TP )
p = plot_grid( 
    legend_b,
    prow,
    ncol = 1, rel_heights = c( .05, 1)
)
p


#pdf("~/Dropbox/PhD/Uniquorn_project/Pub/Images/5_performance.pdf")
jpeg("~/Dropbox/PhD/Uniquorn_project/Pub/Images/5_performance.jpg")
    p
dev.off()


"

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}


plots   = list( g_TP, g_FN, g_FP, g_TN, g_sens, g_F1, g_spec, g_PPV)
#plots_poster = list( g_TP, g_sens, g_spec, g_PPV)
plots_poster = plots
layout_mat_poster = matrix(
  seq( 1, 8 ),
  ncol = 4, nrow = 2
)

g = ggplotGrob(
    plots[[ 1 ]] + 
    theme(
      legend.position='top'
    )
)$grobs

legend  = g[[ which( sapply( g, function( x ) x$name ) == 'guide-box' ) ]]
lheight = sum( legend$height )
  
prepared_plots =  lapply(
  plots_poster,
  function( x ){
        
    x = x + theme( legend.position = 'none' )
    return( x )
  }
)
layout_mat = matrix(
  seq( 1, 8 ),
  ncol = 2, nrow = 4
)

AG = arrangeGrob( 
  grobs = prepared_plots,
  layout_matrix = layout_mat_poster,
  nrow = 2,
  ncol = 4
)

grid.arrange(
  
  legend,
  AG,
  ncol = 1,
  #nrow = 4,
  heights = unit.c(
    lheight,
    unit(
      1,
      'npc'
    ) - lheight
  )
)
"