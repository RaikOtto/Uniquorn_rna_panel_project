#source("~/Dropbox/PhD/Uniquorn_project/Scripts/calculations_vis.R")
library("tidyr")


# relative counts

interval_mat$Weight = as.character(interval_mat$Weight)

relative_plot = ggplot( data = interval_mat, aes( Group) ) 
relative_plot = relative_plot + geom_bar( aes(fill=Weight),position=position_fill() )
relative_plot = relative_plot + ylab("Variant weight")
relative_plot = relative_plot + xlab("")
relative_plot = relative_plot + theme( axis.text.x = element_text(angle = 330, vjust = 0.5) )
relative_plot

# histogram

cl_mat_merge$Group = factor( cl_mat_merge$Group, levels = c("CellMiner raw","CellMiner filtered","CoSMIC raw","CoSMIC filtered","CCLE raw","CCLE filtered") )
cl_mat_merge$Panel = factor( cl_mat_merge$Panel, levels = c("CellMiner","CoSMIC","CCLE") )

vio_all = ggplot( cl_mat_merge, aes( x = Group, y = log2( as.double(Count) )  ) ) 
vio_all = vio_all + geom_violin(alpha=0.5, aes(fill = Panel ) )
vio_all = vio_all + geom_jitter(alpha=0.25, aes( ), position = position_jitter(width = 0.05))
vio_all = vio_all + theme(legend.position=c(.5, .1),legend.direction="horizontal")
vio_all = vio_all + scale_x_discrete(
    "Amount of variants depending on their weight",
    labels = c(
        "CoSMIC raw" = "All","CoSMIC filtered" = "Unique","CCLE raw" = "All",
        "CCLE filtered" = "Unique","CellMiner raw" = "All","CellMiner filtered" = "Unique"
    )) 
vio_all = vio_all + labs(y = "Log2 amount variants", x = "")
vio_all = vio_all + guides(fill=guide_legend(title=NULL))
vio_all = vio_all + scale_fill_manual( values = c("Green","#0072B2","Red") )
vio_all
# absolute_plot

cl_count_mat = cl_count_mat %>% unite( Weight_Group, Group, Weight, sep = "_", remove = F )

if ( ! ( "Count" %in% colnames( cl_abs_count_mat ) ) )
    cl_abs_count_mat$Count = log2(as.double(cl_abs_count_mat$Count))

interval_mat$Weight = factor( interval_mat$Weight )
interval_mat$Group = factor( interval_mat$Group )
cl_abs_count_mat$Weight = factor(cl_abs_count_mat$Weight)
cl_abs_count_mat$Group = factor( cl_abs_count_mat$Group )
cl_abs_count_mat$Count = as.double(cl_abs_count_mat$Count  )

tidy_dist = ggplot( data = cl_abs_count_mat, aes( Weight, log2( Count ) ) )
tidy_dist = tidy_dist + geom_boxplot(aes(fill = Group))
tidy_dist = tidy_dist + ylab("Log 2")
tidy_dist = tidy_dist + xlab("Variant weight")

tidy_dist = tidy_dist + guides(fill=guide_legend(title="Panels"))
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