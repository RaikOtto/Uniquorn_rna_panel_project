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


#jpeg("~/Dropbox/Uniquorn_project/Figures/performance.jpg", width = 1024,height = 512)
    p
#dev.off()

### Confusion matrices

i_table = read.table("~/Uniquorn_data/benchmark_vcf_files/0_5_Benchmark_identification_result.tab",sep="\t", header = T)
i_table[1:5,1:5]

confusion_matrix <<- matrix(integer(), ncol= 5, nrow = 0)

sapply(library_names, FUN = function(lib){
    
    index = grep(i_table$Query, pattern = lib)
    false_neg = sapply( library_names , function(vec){
        
        return(
          sum( str_count(i_table$False_negative[index], pattern = vec) )
        )
    })
    
    confusion_matrix <<- rbind( confusion_matrix, false_neg)
})
rownames(confusion_matrix) = library_names

# gold mat


gold_std = read.table("~/Dropbox/Uniquorn_project/Misc/Goldstandard.tsv",sep ="\t", header = T)
gold_matrix <<- matrix(integer(), ncol= 5, nrow = 0)

sapply(library_names, FUN = function(lib){
  
  index = grep(gold_std$CL, pattern = lib)
  pos_TP = sapply( library_names , function(vec){
    
    return(
      sum( str_count(gold_std$Merged[index], pattern = vec) )
    )
  })
  
  gold_matrix <<- rbind( gold_matrix, pos_TP)
})
rownames(gold_matrix) = library_names
new_mat = confusion_matrix / gold_matrix

#library_names[library_names=="EGA"] = "Klijn"
#library_names[library_names=="COSMIC"] = "CGP"

melted_cormat <- melt(confusion_matrix, na.rm = TRUE, as.is = T)
melted_cormat$Var1[ melted_cormat$Var1 == "COSMIC" ] = "CGP" 
melted_cormat$Var2[ melted_cormat$Var2 == "COSMIC" ] = "CGP"
melted_cormat$Var1[ melted_cormat$Var1 == "EGA" ] = "Klijn"
melted_cormat$Var2[ melted_cormat$Var2 == "EGA" ] = "Klijn"
colnames(confusion_matrix)[colnames(confusion_matrix) == "COSMIC"] = "CGP"
colnames(confusion_matrix)[colnames(confusion_matrix) == "EGA"] = "Klijn"
rownames(confusion_matrix)[rownames(confusion_matrix) == "COSMIC"] = "CGP"
rownames(confusion_matrix)[rownames(confusion_matrix) == "EGA"] = "Klijn"

# reorder
#pheatmap::pheatmap(confusion_matrix) + geom_text(aes(label = round(value, 1)))

library("ggcorrplot")
library(ggplot2)
library(scales) # for muted function

rel_plot = ggcorrplot(new_mat , method = "circle", hc.order = TRUE, lab = TRUE)+ 
    scale_fill_gradient2(
    low = ("green"), 
    mid = "white", 
    high = "red", midpoint = mean(new_mat)
)+ guides(fill=guide_legend(title="% FN"))+ theme(legend.position ="top", text=element_text(size=14,face="bold"))

abs_plot = ggplot(melted_cormat, aes(Var1, Var2)) + # x and y axes => Var1 and Var2
  geom_tile(aes(fill = value)) + # background colours are mapped according to the value column
  geom_text(aes( label = round(melted_cormat$value, 2)), col = "black") + # write the values
  scale_fill_gradient2(low = ("green"), 
                       mid = "white", 
                       high = "red", midpoint = 30) +
  theme(panel.grid.major.x=element_blank(), #no gridlines
        panel.grid.minor.x=element_blank(), 
        panel.grid.major.y=element_blank(), 
        panel.grid.minor.y=element_blank(),
        panel.background=element_rect(fill="white"), # background=white
        axis.text.x = element_text(angle=90, hjust = 1,vjust=1,size = 12,face = "bold"),
        plot.title = element_text(size=20,face="bold"),
        axis.text.y = element_text(size = 12,face = "bold"),
        legend.position = "top") + 
  ggtitle("") + 
  theme(text=element_text(face="bold", size=14)) + 
  scale_x_discrete(name="") +
  scale_y_discrete(name="") +
  labs(fill="FNs")+ theme(axis.text.x = element_text(angle = 45, hjust = 1))

prow2 = plot_grid(
  abs_plot,
  rel_plot,
  labels=c("A", "B"),
  ncol = 2,
  nrow = 1,
  scale = 1.04
)+ theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))

jpeg("~/Dropbox/Uniquorn_project/Figures/4_FN_sources.jpg", width = 2048,height = 1536)
    prow2
dev.off()

p2 = plot_grid(
  prow,
  ncol = 1, rel_heights = c( .05, 1)
)
p