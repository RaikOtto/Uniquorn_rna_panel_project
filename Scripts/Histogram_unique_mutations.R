source("~/Uniquorn_rna_panel_project//Scripts/calculations_vis.R")
library("tidyr")

g_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
n = 5
cols = g_color_hue(n)


# absolute_plot

abs_plot = ggplot( data = abs_mat[abs_mat$Weight=="W0",], aes( Library, log2( Count ) ) )
abs_plot = abs_plot + geom_boxplot(aes(fill = Library))
abs_plot = abs_plot + ylab("Log 2 Absolute Variant Count") + theme(legend.position="bottom")
abs_plot = abs_plot + xlab("") + theme( axis.text.x = element_blank() )
abs_plot = abs_plot + guides(fill=guide_legend(title="Libraries"))
abs_plot = abs_plot + scale_fill_manual(
  values=c(cols), 
  name="Technologies",
  labels = c("Exome-Genome Analyzer  ","Hybrid capture  ","Exome-HiSeq 2500  ","RNA-HiSeq 2000  ","RNA-HiSeq 2500"))
abs_plot = abs_plot + theme(legend.text=element_text(size=13))
# mean counts

sds = aggregate( log(as.double(ref_counts$Count)), by = list(ref_counts$Library), FUN = sd)[c(-3,-7,-8),2]
means = log2(as.double( mean_mat[mean_mat$Weight == "W0",3] ))

sds_min = means - sds / 2
sds_max = means + sds / 2

mean_plot = ggplot( data = mean_mat[mean_mat$Weight == "W0",], aes( x = Library, y = log2(as.double(Count)), fill = Library) ) 
mean_plot = mean_plot + geom_bar( stat= "identity" )
mean_plot = mean_plot + coord_cartesian(ylim = c(6, 11))
mean_plot = mean_plot + ylab("Log2 Mean Variant Count")
mean_plot = mean_plot + xlab("") #+ theme(legend.position="none")
mean_plot = mean_plot + theme( axis.text.x = element_blank() )
mean_plot = mean_plot + geom_errorbar(aes( ymin = sds_min, ymax = sds_max), width=.2,position=position_dodge(.9))
mean_plot = mean_plot + coord_cartesian(ylim = c(6, 12)) 
mean_plot

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
  ncol = 1, rel_heights = c( .1, .9)
)
p

#jpeg("~/Dropbox/Uniquorn_project//Figures/2_Library_composition.jpg", width = 1024,height = 512)
p
#dev.off()
