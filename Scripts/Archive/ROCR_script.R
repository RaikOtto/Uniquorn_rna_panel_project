### ROCR
source("~/Uniquorn_rna_panel_project/Scripts/utility.R")
library("ROCR")
library("grid")
library("ggplot2")

perf_vec <<- c()
anchor = "~/Uniquorn_data/benchmark_vcf_files/AUC/"

run_rocr = function( anchor ){
  
  col_type = c(double(),integer(),logical() )
  
    auc_t_0    = read.table( 
        paste0( anchor, "auc_0.tsv" ), 
        header = T, stringsAsFactors = F, sep = "\t", col.names = c("P_values","Should_be_found")
    )
    auc_t_0_25    = read.table( 
        paste0( anchor, "auc_0_25.tsv" ), 
        header = T, stringsAsFactors = F, sep = "\t", col.names = c("P_values","Should_be_found")
    )
    auc_t_0_5    = read.table( 
        paste0( anchor, "auc_0_5.tsv" ), 
        header = T, stringsAsFactors = F, sep = "\t", col.names = c("P_values","Should_be_found")
    )
    auc_t_1    = read.table( 
        paste0( anchor, "auc_1.tsv" ), 
        header = T, stringsAsFactors = F, sep = "\t", col.names = c("P_values","Should_be_found")
    )

  for (input_t in c( "auc_t_0", "auc_t_0_25", "auc_t_0_5", "auc_t_1" )){
  
      print(input_t)
      auc_t = eval(
        parse(
          text = input_t
        )
      )
      
      many_labs = as.character( auc_t$Should_be_found )
      many_labs[ many_labs != "TRUE" ] = 0.0
      many_labs[ many_labs == "TRUE" ] = 1.0
      
      pred_obj = ROCR::prediction(
        predictions = as.double( auc_t$P_values  ),
        labels = many_labs 
      )
      pred_obj_2 <<- pred_obj
      
      rocr_auc = as.character(
          round(
              as.double(
                  unclass(
                      performance(
                          pred_obj,
                          "auc"
                      )
                  )@"y.values"
              ),
          2 )
      )
      
      print( paste0( c( "AUC:", rocr_auc) ), collapse = " " )
      perf = ROCR::performance(
          prediction.obj = pred_obj,
          measure = "tpr",
          x.measure = "fpr"
      );
  
      perf_vec <<- c(perf_vec,perf)
      }
}

run_rocr( anchor )

# graphics

perf_1 = unlist(perf_vec)[[1]]
perf_2 = unlist(perf_vec)[[2]]
perf_3 = unlist(perf_vec)[[3]]
perf_4 = unlist(perf_vec)[[4]]

x_val_1 = 100 - ( as.double(unlist(perf_1@x.values ) ) ) * 100
x_val_2 = 100 - ( as.double(unlist(perf_2@x.values ) ) ) * 100
x_val_3 = 100 - ( as.double(unlist(perf_3@x.values ) ) ) * 100
x_val_4 = 100 - ( as.double(unlist(perf_4@x.values ) ) ) * 100
#x_val_1[1] = 0
#x_val_1[length(x_val_1)] = 100
#x_val_2[1] = 0
#x_val_2[length(x_val_2)] = 100
#x_val_3[1] = 0
#x_val_3[length(x_val_3)] = 100
#x_val_4[1] = 0
#x_val_4[length(x_val_4)] = 100

weights = as.factor( 
  c( 
    rep("0",length(x_val_1)),
    rep("0.25",length(x_val_2)),
    rep("0.5",length(x_val_3)),
    rep("1.0",length(x_val_4) )
  )
)

y_val_1 = 100 - ( as.double( unlist(perf_1@y.values) ) ) * 100
y_val_2 = 100 - ( as.double( unlist(perf_2@y.values) ) ) * 100
y_val_3 = 100 - ( as.double( unlist(perf_3@y.values) ) ) * 100
y_val_4 = 100 - ( as.double( unlist(perf_4@y.values) ) ) * 100

plot( x= (x_val_2), y = (y_val_2), xlim = c(0,.1))

plots_frame = data.frame( 
  "X" = c(
      x_val_1,
      x_val_2, 
      x_val_3,
      x_val_4
  ),
  "Y" = c(
      y_val_1,
      y_val_2,
      y_val_3,
      y_val_4
  ),
  "Weight" = weights
)

q_bird = ggplot(
  data = plots_frame, 
  aes( x = X, y = Y )
) + geom_line( size = 2, aes( linetype = Weight, color = Weight) )
q_bird  = q_bird + xlim( 0, 0.0025)
q_bird  = q_bird + xlab( "" ) + ylab( "" )
q_bird = q_bird + theme( 
  panel.background = element_blank(),
  legend.position  = "none",
  panel.border = element_rect(colour = "black", fill=NA, size = 2 ),
  plot.margin = unit(
    c( .0,.0,.01,.01), "cm"
  )
)

q_bird


q = ggplot( 
  data = plots_frame_2,
  aes( x = X, y = Y )
) + geom_line( size = 2, aes(linetype=Weight, color = Weight) )
q = q + xlab( "FPR" ) + ylab( "Sensitivity %" )
q = q + coord_cartesian(xlim = c(0, 0.000012))
q = q + theme(
    legend.position="top",
    legend.background = element_rect(
        fill="gray90",
        size = 5
    ),
    legend.key.width = unit( 2,"cm" ),
    plot.margin = unit(
        c( .5,.5,.5,.5), "cm"
    )
)
q = q + guides(fill = guide_legend(keywidth = 100, keyheight = 10))

q = q + geom_vline( xintercept = x_val_1, lwd = 2 )

#pdf("~/Uniquorn_rna_panel_project/Misc/Images/6_ROC.pdf")

    q
    print( q_bird, vp = grid::viewport(.8, .27, .4, .4 ) )
dev.off()

opt.cut = function( perf, pred ){
  cut.ind = mapply( FUN = function( x, y, p ){
    d = (x - 0)^2 + ( y - 1 )^2
    ind = which(d == min(d))
    c(
      sensitivity = y[[ ind ]],
      specificity = 1 - x[[ ind ]], 
      cutoff = p[[ ind ]])
  }, perf@x.values, perf@y.values, pred@cutoffs)
}


auc_t = read.table( 
    paste0( anchor, "auc_1_relaxed.tab" ), 
    header = T,
    stringsAsFactors = F,
    sep = "\t",
    col.names = c("Conf_Score","Mut_Found","Should_be_found")
)

many_labs = as.character( auc_t$Should_be_found )
many_labs[ many_labs != "TRUE" ] = 0.0
many_labs[ many_labs == "TRUE" ] = 1.0

pred_obj = ROCR::prediction(
  predictions = as.double( auc_t$Conf_Score  ),
  labels = many_labs 
)

perf = ROCR::performance( prediction.obj = pred_obj, measure = "tpr", x.measure = "fpr")

print(opt.cut( perf, pred_obj))

h = cbind(
    as.character(
        unlist( pred_obj@cutoffs )
    ),
    as.character(
        unlist( pred_obj@tp )
    ),
    as.character(
        unlist( pred_obj@fp )
    )
)