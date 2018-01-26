### boxplot benchmark

require(ggplot2)
require(grid)
require(gridExtra)
library("cowplot")
library("scales")
source("~/Uniquorn_rna_panel_project/Scripts/calculations_vis.R")

# Sensitivity

g_sens = ggplot( data = benchmark_Sensitivity, aes( x = Seq_Type, y = Sensitivity))
g_sens = g_sens + geom_boxplot(aes(colour = Weight),size=1)
g_sens = g_sens + theme(legend.position="top", legend.background = element_rect(fill="gray90", size =5 )  )
g_sens = g_sens + geom_segment(aes( x = .8,xend = 1,
                                y = benchmark_Sensitivity[ which( (benchmark_Sensitivity[,1] == "0.5") & (benchmark_Sensitivity[,3] == "AVG") )[4],2], 
                                yend = benchmark_Sensitivity[ which( (benchmark_Sensitivity[,1] == "0.5") & (benchmark_Sensitivity[,3] == "AVG") )[4],2]), size = 1.5)
g_sens = g_sens + geom_segment(aes( x = 1.8,xend = 2,
                                y = benchmark_Sensitivity[ which( (benchmark_Sensitivity[,1] == "0.5") & (benchmark_Sensitivity[,3] == "RNA") )[4],2], 
                                yend = benchmark_Sensitivity[ which( (benchmark_Sensitivity[,1] == "0.5") & (benchmark_Sensitivity[,3] == "RNA") )[4],2]), size = 1.5)
g_sens = g_sens + geom_segment(aes( x = 2.8,xend = 3,
                                y = benchmark_Sensitivity[ which( (benchmark_Sensitivity[,1] == "0.5") & (benchmark_Sensitivity[,3] == "ClearSight") )[4],2], 
                                yend = benchmark_Sensitivity[ which( (benchmark_Sensitivity[,1] == "0.5") & (benchmark_Sensitivity[,3] == "ClearSight") )[4],2]), size = 1.5)
g_sens = g_sens + geom_segment(aes( x = 3.8,xend = 4,
                                y = benchmark_Sensitivity[ which( (benchmark_Sensitivity[,1] == "0.5") & (benchmark_Sensitivity[,3] == "TruSight") )[4],2], 
                                yend = benchmark_Sensitivity[ which( (benchmark_Sensitivity[,1] == "0.5") & (benchmark_Sensitivity[,3] == "TruSight") )[4],2]), size = 1.5)
g_sens = g_sens + geom_segment(aes( x = 4.8,xend = 5,
                                y = benchmark_Sensitivity[ which( (benchmark_Sensitivity[,1] == "0.5") & (benchmark_Sensitivity[,3] == "Hotspot_v2") )[4],2], 
                                yend = benchmark_Sensitivity[ which( (benchmark_Sensitivity[,1] == "0.5") & (benchmark_Sensitivity[,3] == "Hotspot_v2") )[4],2]), size = 1.5)
g_sens

# F1

g_F1 = ggplot( data = benchmark_F1, aes( x = Seq_Type, y = F1))
g_F1 = g_F1 + geom_boxplot(aes(colour = Weight),size=1)
g_F1 = g_F1 + theme(legend.position="top", legend.background = element_rect(fill="gray90", size =5 )  )
g_F1 = g_F1 + geom_segment(aes( x = .8,xend = 1,
                                  y = benchmark_F1[ which( (benchmark_F1[,1] == "0.5") & (benchmark_F1[,3] == "AVG") )[4],2], 
                                  yend = benchmark_F1[ which( (benchmark_F1[,1] == "0.5") & (benchmark_F1[,3] == "AVG") )[4],2]), size = 1.5)
g_F1 = g_F1 + geom_segment(aes( x = 1.8,xend = 2,
                                  y = benchmark_F1[ which( (benchmark_F1[,1] == "0.5") & (benchmark_F1[,3] == "RNA") )[4],2], 
                                  yend = benchmark_F1[ which( (benchmark_F1[,1] == "0.5") & (benchmark_F1[,3] == "RNA") )[4],2]), size = 1.5)
g_F1 = g_F1 + geom_segment(aes( x = 2.8,xend = 3,
                                  y = benchmark_F1[ which( (benchmark_F1[,1] == "0.5") & (benchmark_F1[,3] == "ClearSight") )[4],2], 
                                  yend = benchmark_F1[ which( (benchmark_F1[,1] == "0.5") & (benchmark_F1[,3] == "ClearSight") )[4],2]), size = 1.5)
g_F1 = g_F1 + geom_segment(aes( x = 3.8,xend = 4,
                                  y = benchmark_F1[ which( (benchmark_F1[,1] == "0.5") & (benchmark_F1[,3] == "TruSight") )[4],2], 
                                  yend = benchmark_F1[ which( (benchmark_F1[,1] == "0.5") & (benchmark_F1[,3] == "TruSight") )[4],2]), size = 1.5)
g_F1 = g_F1 + geom_segment(aes( x = 4.8,xend = 5,
                                  y = benchmark_F1[ which( (benchmark_F1[,1] == "0.5") & (benchmark_F1[,3] == "Hotspot_v2") )[4],2], 
                                  yend = benchmark_F1[ which( (benchmark_F1[,1] == "0.5") & (benchmark_F1[,3] == "Hotspot_v2") )[4],2]), size = 1.5)
g_F1

# PPV

g_PPV = ggplot( data = benchmark_PPV, aes( x = Seq_Type, y = PPV))
g_PPV = g_PPV + geom_boxplot(aes(colour = Weight),size=1)
g_PPV = g_PPV + theme(legend.position="top", legend.background = element_rect(fill="gray90", size =5 )  )
g_PPV = g_PPV + geom_segment(aes( x = .8,xend = 1,
  y = benchmark_PPV[ which( (benchmark_PPV[,1] == "0.5") & (benchmark_PPV[,3] == "AVG") )[4],2], 
  yend = benchmark_PPV[ which( (benchmark_PPV[,1] == "0.5") & (benchmark_PPV[,3] == "AVG") )[4],2]), size = 1.5)
g_PPV = g_PPV + geom_segment(aes( x = 1.8,xend = 2,
  y = benchmark_PPV[ which( (benchmark_PPV[,1] == "0.5") & (benchmark_PPV[,3] == "RNA") )[4],2], 
  yend = benchmark_PPV[ which( (benchmark_PPV[,1] == "0.5") & (benchmark_PPV[,3] == "RNA") )[4],2]), size = 1.5)
g_PPV = g_PPV + geom_segment(aes( x = 2.8,xend = 3,
  y = benchmark_PPV[ which( (benchmark_PPV[,1] == "0.5") & (benchmark_PPV[,3] == "ClearSight") )[4],2], 
  yend = benchmark_PPV[ which( (benchmark_PPV[,1] == "0.5") & (benchmark_PPV[,3] == "ClearSight") )[4],2]), size = 1.5)
g_PPV = g_PPV + geom_segment(aes( x = 3.8,xend = 4,
  y = benchmark_PPV[ which( (benchmark_PPV[,1] == "0.5") & (benchmark_PPV[,3] == "TruSight") )[4],2], 
  yend = benchmark_PPV[ which( (benchmark_PPV[,1] == "0.5") & (benchmark_PPV[,3] == "TruSight") )[4],2]), size = 1.5)
g_PPV = g_PPV + geom_segment(aes( x = 4.8,xend = 5,
  y = benchmark_PPV[ which( (benchmark_PPV[,1] == "0.5") & (benchmark_PPV[,3] == "Hotspot_v2") )[4],2], 
  yend = benchmark_PPV[ which( (benchmark_PPV[,1] == "0.5") & (benchmark_PPV[,3] == "Hotspot_v2") )[4],2]), size = 1.5)
g_PPV

# prow

prow = plot_grid(
    g_sens + theme(legend.position="none") + theme( axis.text.x = element_text(angle=45, vjust = .5), axis.title.x = element_blank() ),
    g_F1   + theme(legend.position="none") + theme( axis.text.x = element_text(angle=45, vjust = .5), axis.title.x = element_blank() ),
    g_PPV  + theme(legend.position="none") + theme( axis.text.x = element_text(angle=45, vjust = .5), axis.title.x = element_blank() ), 
    labels=c("A", "B","C"), 
    ncol = 3,
    nrow = 1
)

legend_b = get_legend( g_F1 )
p = plot_grid( 
    legend_b,
    prow,
    ncol = 1, rel_heights = c( .05, 1)
)
p


jpeg("~/Dropbox/Uniquorn_project/Figures/performance.jpg", width = 1024,height = 512)
    p
dev.off()

#

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


### Confusion matrices