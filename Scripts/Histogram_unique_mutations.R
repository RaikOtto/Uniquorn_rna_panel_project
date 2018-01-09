#source("~/Dropbox/PhD/Uniquorn_project/Scripts/calculations_vis.R")
library("tidyr")

##

abs_mat = variant_t
abs_mat$Count = log(abs_mat$Count)

# histogram

aggregate( abs_mat$Count, by = list(abs_mat$Library), FUN = mean )
aggregate( abs_mat$Count, by = list( paste(abs_mat$Library,abs_mat$Weight )) , FUN = mean )

# relative counts

interval_mat$Weight = as.character(interval_mat$Weight)

relative_plot = ggplot( data = interval_mat, aes( Group) ) 
relative_plot = relative_plot + geom_bar( aes(fill=Weight),position=position_fill() )
relative_plot = relative_plot + ylab("Variant weight")
relative_plot = relative_plot + xlab("")
relative_plot = relative_plot + theme( axis.text.x = element_text(angle = 330, vjust = 0.5) )
relative_plot

# absolute_plot

tidy_dist = ggplot( data = abs_mat, aes( Weight, Count ) )
tidy_dist = tidy_dist + geom_boxplot(aes(fill = Library))
tidy_dist = tidy_dist + ylab("Log 2")
tidy_dist = tidy_dist + xlab("Variant weight")

tidy_dist = tidy_dist + guides(fill=guide_legend(title="Libraries"))
tidy_dist

# combine plots

grid.arrange( vio_all, relative_plot, tidy_dist, 
  ncol = 2,
  nrow = 2, 
  layout_matrix = rbind(c(1,1),c(2,3)),
  padding = 3
)

library("reshape2")
library("gridExtra")
library("cowplot")

ggdraw() +
    draw_plot( vio_all, 0, .5, 1, .5) +
    draw_plot( relative_plot, 0, 0, .5, .5) +
    draw_plot( tidy_dist, .5, 0, .5, .5) +
    draw_plot_label(c("A", "B", "C"), c(0, 0, 0.5), c(1, 0.5, 0.5),
    size = 15
)

### end printing