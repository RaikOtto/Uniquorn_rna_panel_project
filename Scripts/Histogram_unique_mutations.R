source("~/Uniquorn_rna_panel_project//Scripts/calculations_vis.R")
library("tidyr")

# absolute_plot

abs_plot = ggplot( data = abs_mat[abs_mat$Weight=="W0",], aes( Library, log2( Count ) ) )
abs_plot = abs_plot + geom_boxplot(aes(fill = Library))
abs_plot = abs_plot + ylab("Log 2 Absolute Variant Count") + theme(legend.position="bottom")
abs_plot = abs_plot + xlab("") + theme( axis.text.x = element_blank() )
abs_plot = abs_plot + guides(fill=guide_legend(title="Libraries"))
#abs_plot

# mean counts

mean_plot = ggplot( data = mean_mat[mean_mat$Weight == "W0",], aes( x = Library, y = log2(as.double(Count)), fill = Library) ) 
#mean_plot = mean_plot + geom_bar( stat= "identity" )
mean_plot = mean_plot + geom_bar( stat= "identity" )
mean_plot = mean_plot + coord_cartesian(ylim = c(6, 11))
mean_plot = mean_plot + ylab("Log2 Mean Variant Count")
mean_plot = mean_plot + xlab("") #+ theme(legend.position="none")
mean_plot = mean_plot + theme( axis.text.x = element_blank() )
#mean_plot

# combine plots

grid.arrange( mean_plot, tidy_dist, 
  ncol = 1,
  nrow = 2,
  layout_matrix = rbind(c(1,2)),
  padding = 3
)

# prow
library("cowplot")

legend_b = get_legend( abs_plot )
prow = plot_grid(
  abs_plot + theme(legend.position="none") + theme( axis.text.x = element_blank() ),
  mean_plot + theme(legend.position="none") + theme( axis.text.x = element_blank() ), 
  labels=c("A", "B"), 
  ncol = 2,
  nrow = 1
)
p = plot_grid( 
  legend_b,
  prow,
  ncol = 1, rel_heights = c( .05, .95)
)
p

jpeg("~/Dropbox/Uniquorn_project//Figures/2_Library_composition.jpg", width = 1024,height = 512)
p
dev.off()
