
### ccl_mat 

ccl_mat = data.frame(
  "Library" = libraries,
  "CCLs" = as.integer(c(904, 60,1020,675,937)),
  stringsAsFactors = F
)


  
  ### Per_ccl_success
  
  avg_file = read.table("~/Uniquorn_data/benchmark_vcf_files/0_5_Benchmark_identification_result.tab", sep = "\t", header = T, stringsAsFactors = F)
  #c = show_contained_variants_for_ccl("105KC", library_name = "EGA", mutational_weight_inclusion_threshold = 0.5)
  
  ref_ccls_cellminer = show_contained_variants_in_library( library_name = "CELLMINER", mutational_weight_inclusion_threshold = 0.5)$Member_CCLs
  ref_ccls_cellminer = as.character(unlist( str_split( ref_ccls_cellminer, pattern = "," )))
  ref_ccls_cellminer = str_replace( ref_ccls_cellminer, pattern = "_CELLMINER", "")
  ref_ccls_cellminer = str_replace_all( ref_ccls_cellminer, pattern = "_", "")
  ref_ccls_cellminer = as.character(unlist( str_split( ref_ccls_cellminer, pattern = "," )))
  ref_ccls_cellminer = paste(ref_ccls_cellminer, "CELLMINER", sep = "_")
  
  ref_ccls_cosmic = show_contained_variants_in_library( library_name = "COSMIC", mutational_weight_inclusion_threshold = 0.5)$Member_CCLs
  ref_ccls_cosmic = str_to_upper( ref_ccls_cosmic )
  ref_ccls_cosmic = str_replace_all( ref_ccls_cosmic, pattern = "-", "")
  ref_ccls_cosmic = str_replace_all( ref_ccls_cosmic, pattern = "\\.", "")
  ref_ccls_cosmic = as.character(unlist( str_split( ref_ccls_cosmic, pattern = "," )))
  ref_ccls_cosmic = paste(ref_ccls_cosmic, "COSMIC", sep = "_")
  
  ref_ccls_ccle = show_contained_variants_in_library( library_name = "CCLE", mutational_weight_inclusion_threshold = 0.5)$Member_CCLs
  ref_ccls_ccle = as.character(unlist( str_split( ref_ccls_ccle, pattern = "," )))
  ref_ccls_ccle = paste(ref_ccls_ccle, "CCLE", sep = "_")
  
  ref_ccls_ega = show_contained_variants_in_library( library_name = "EGA", mutational_weight_inclusion_threshold = 0.5)$Member_CCLs
  ref_ccls_ega = str_replace_all( ref_ccls_ega, pattern = "\\.", "")
  ref_ccls_ega = str_replace_all( ref_ccls_ega, pattern = "-", "")
  ref_ccls_ega = as.character(unlist( str_split( ref_ccls_ega, pattern = "," )))
  
  ref_ccls_gdc = show_contained_variants_in_library( library_name = "GDC", mutational_weight_inclusion_threshold = 0.5)$Member_CCLs
  ref_ccls_gdc = str_replace_all( ref_ccls_gdc, pattern = "\\.", "")
  ref_ccls_gdc = str_replace_all( ref_ccls_gdc, pattern = "-", "")
  ref_ccls_gdc = as.character(unlist( str_split( ref_ccls_gdc, pattern = "," )))
  
  ref_counts <<- data.frame(
    "CCL" = as.character(),
    "Count" = as.character(),
    "Library" = as.character(),
    "Sensitivity" = as.character(),
    "F1" = as.character(),
    "TP" = as.character(),
    "FP" = as.character(),
    "FN" = as.character()
  )
  
  #all_ccls = c(ref_ccls_cellminer,ref_ccls_cosmic,ref_ccls_ccle,ref_ccls_ega,ref_ccls_gdc)
  
  all_ccls = avg_file$Query
  all_ccls_uni = unique( as.character(unlist( str_split( all_ccls, pattern = "," ))))
  
  i <<- 0
  
  count_per_ccl = sapply( all_ccls_uni, FUN= function(ref_ccl){
    
    ref_ccl = as.character(ref_ccl)
    ref_ccl = str_replace_all( ref_ccl, pattern = "\\)", "")
    ref_ccl = str_replace_all( ref_ccl, pattern = "\\(", "")
    ref_ccl = str_replace_all( ref_ccl, pattern = ":", "")
    ref_ccl = str_to_upper( ref_ccl)
    
    library = as.character( tail( as.character( unlist( str_split(ref_ccl, pattern = "_") ) ), 1 ) )
    #library_panel = "ClearSeq"
    #library_panel = "TruSight"
    library_panel = "Hotspot_v2"
    
    ref_ccl_pure = str_replace( ref_ccl, pattern = paste("",library, sep = "_"),""  )
    ref_ccl = str_replace_all( ref_ccl, pattern = "_", "")
    ref_ccl = str_replace_all( ref_ccl, pattern = library, "")
    
    ref_ccl = paste(ref_ccl, library,sep = "_") #
    
    count = 49 #as.character( sum ( str_detect( all_ccls, pattern = ref_ccl ) ))
    
    i <<- i + 1
    print(i)
    
    ref_match = match( ref_ccl, str_to_upper(avg_file$Query))
    
    tp = as.character( unlist( str_split( avg_file$True_positive[ref_match], pattern = ",") ) )
    tp = tp[tp != ""]
    tp = as.integer( length(tp) )
    fp = as.character( unlist( str_split( avg_file$False_positive[ref_match], pattern = ",") ) )
    fp = fp[fp != ""]
    fp = as.integer( length(fp) )
    fn = as.character( unlist( str_split( avg_file$False_negative[ref_match], pattern = ",") ) )
    fn = fn[fn != ""]
    fn = as.integer( length(fn) )
    Sensitivity = round( tp / (tp + fn) * 100, 0 )
    Specificity = 99.9
    F1 = round( 2 * (Sensitivity * Specificity) / (Sensitivity + Specificity), 0)
    PPV = round( tp / (tp + fp)  * 100,0)
    
    if(F1 == "NaN")
      F1 = "0"
    if(Sensitivity == "NaN")
      Sensitivity = "0"
    if(PPV == "NaN")
      PPV = "0"
    
    
    ref_counts <<- data.frame(
      "CCL" = c(  as.character( ref_counts$CCL) , ref_ccl_pure ), ####
      "Count" = c( as.character( ref_counts$Count), count),
      "Library" = c( as.character( ref_counts$Library), library_panel),
      "Sensitivity" = c( as.character( ref_counts$Sensitivity), as.character( Sensitivity)),
      "F1" = c( as.character( ref_counts$F1), as.character( F1)),
      "PPV" = c( as.character( ref_counts$PPV), as.character( PPV)),
      "TP" = c( as.character( ref_counts$TP), as.character( tp)),
      "FP" = c( as.character( ref_counts$FP), as.character( fp)),
      "FN" = c( as.character( ref_counts$FN), as.character( fn))
    )
  }
  )
  
  table(ref_counts$Library)
  ref_counts[(ncol(ref_counts)-5) : ncol(ref_counts),]
  ref_counts[ref_counts$Sensitivity == "NaN",]
  ref_counts$Sensitivity[ref_counts$Sensitivity == "NaN"] = "0"
  ref_counts$PPV[ref_counts$PPV == "NaN"] = "0"
  ref_counts$PPV[is.na(ref_counts$PPV)] = "0"
  
  write.table(ref_counts,"~/Uniquorn_data/benchmark_vcf_files/Benchmark_Per_CCL_stat.tsv",sep ="\t", quote = F, row.names = F)
  
  #### integrate panel seq
  
  #avg_file = read.table("~/Uniquorn_data/benchmark_vcf_files/Finished/Results/ClearSight/Conf_3_top_1/0_5_Benchmark_identification_result.tab", sep = "\t", header = T, stringsAsFactors = F)
  #avg_file = read.table("~/Uniquorn_data/benchmark_vcf_files/Finished/Results/TruSight/panel_Conf_3_top1/0_5_Benchmark_identification_result.tab", sep = "\t", header = T, stringsAsFactors = F)
  avg_file = read.table("~/Uniquorn_data/benchmark_vcf_files/Finished/Results/Hotspot/Conf_3_top_1/0_5_Benchmark_identification_result.tab", sep = "\t", header = T, stringsAsFactors = F)
  
  ####
  
  ref_counts = read.table("~/Uniquorn_data/benchmark_vcf_files/Benchmark_Per_CCL_stat.tsv", sep ="\t", header = T,stringsAsFactors = F)
  
  perf_mat = reshape2::melt(ref_counts, value.name = "CCL")
  
  # barplot
  
  ref_counts$TP = as.integer(as.character(ref_counts$TP))
  ref_counts$FN = as.integer(as.character(ref_counts$FN))
  ref_counts$FP = as.integer(as.character(ref_counts$FP))
  ref_counts$Sensitivity = as.integer(as.character(ref_counts$TP)) / (as.integer(as.character(ref_counts$TP)) + as.integer(as.character(ref_counts$FN)))  * 100
  ref_counts$Specificity = 99.9
  
  ref_counts$F1 = 2 * ( as.integer(as.character(ref_counts$Sensitivity) ) * as.integer(as.character(ref_counts$Specificity) ) ) /
    ( as.integer(as.character(ref_counts$Sensitivity) ) + as.integer(as.character(ref_counts$Specificity) ) )
  ref_counts$PPV = as.integer(as.character(ref_counts$TP)) /
    ( as.integer(as.character(ref_counts$TP) ) + as.integer(as.character(ref_counts$FP) ) ) * 100
  
  perf_mat = aggregate(ref_counts$Sensitivity, by = list(ref_counts$Library), FUN = mean)
  perf_mat$Sens_SD = aggregate(ref_counts$Sensitivity, by = list(ref_counts$Library), FUN = sd)[,2]
  perf_mat$F1 = aggregate(ref_counts$F1, by = list(ref_counts$Library), FUN = mean)[,2]
  perf_mat$SF1_SD = aggregate(ref_counts$F1, by = list(ref_counts$Library), FUN = sd)[,2]
  
  ref_counts$PPV[ ( is.na(ref_counts$PPV) ) ] = 0
  
  perf_mat$PPV = aggregate(ref_counts$PPV, by = list(ref_counts$Library), FUN = mean)[,2]
  perf_mat$PPV_SD = aggregate(ref_counts$PPV, by = list(ref_counts$Library), FUN = sd)[,2]
  perf_mat = as.data.frame(perf_mat)
  
  colnames(perf_mat) = c("Technology", "Sensitivity","Sensitivity_SD","F1","F1_SD","PPV","PPV_SD")
  vis_mat = reshape2::melt(perf_mat)
  colnames(vis_mat) = c("Technology","Parameter","Value")
  vis_mat$Technology = as.character(vis_mat$Technology)
  
  vis_mat$Technology[vis_mat$Technology == "COSMIC"] = "Exome_HiSeq_2500"
  vis_mat$Technology[vis_mat$Technology == "CCLE"] = "Hybrid_capture"
  vis_mat$Technology[vis_mat$Technology == "EGA"] = "RNA_HiSeq_2000"
  vis_mat$Technology[vis_mat$Technology == "GDC"] = "RNA_HiSeq_2500"
  vis_mat$Technology[vis_mat$Technology == "CELLMINER"] = "Exome_Genome_Analyzer"
  vis_mat$Technology = factor(vis_mat$Technology, levels = c("RNA_HiSeq_2500","RNA_HiSeq_2000","Exome_HiSeq_2500","Hybrid_capture","Exome_Genome_Analyzer","ClearSeq","TruSight","Hotspot_v2"))
  
  means = c(perf_mat$Sensitivity, perf_mat$F1, perf_mat$PPV)
  sds = c(perf_mat$Sensitivity_SD, perf_mat$F1_SD, perf_mat$PPV_SD)
  sds_min = means - sds / 2
  sds_max = means + sds / 2
  
  g_bench = ggplot(
    subset(vis_mat,Parameter %in% c("Sensitivity","F1","PPV")),
    aes( x = Parameter, y = Value, fill = Technology))
  g_bench = g_bench + geom_bar(stat="identity", position=position_dodge())
  g_bench = g_bench +  geom_errorbar(aes( ymin = sds_min, ymax = sds_max), width=.2,position=position_dodge(.9))
  g_bench = g_bench + theme(legend.position="top", legend.background = element_rect(fill="gray90", size =5 )  )
  g_bench + xlab("")

#

counts = aggregate(ref_counts$Count, by = list(ref_counts$Library), FUN = mean)[,2]
count_vs_sens_sd = count_vs_sens = as.data.frame( cbind(log(counts),perf_mat$Technology,perf_mat$Sensitivity,perf_mat$Sensitivity_SD) )
colnames(count_vs_sens) = colnames(count_vs_sens_sd) = c("Log_counts","Technology","Sensitivity","Sensitivity_SD")
count_vs_ppv_sd = count_vs_ppv = as.data.frame( cbind(log(counts),perf_mat$Technology,perf_mat$PPV,perf_mat$PPV_SD) )
colnames(count_vs_ppv_sd) = colnames(count_vs_ppv) = c("Log_counts","Technology","PPV","PPV_SD")

# sensitivity

count_vs_sens$Log_counts = as.double(as.character(count_vs_sens$Log_counts))
count_vs_sens$Sensitivity = as.double(as.character(count_vs_sens$Sensitivity))
count_vs_sens$Sensitivity_SD = as.double(as.character(count_vs_sens$Sensitivity_SD))

estimator_sens = lm(formula = Sensitivity ~ Log_counts, data = count_vs_sens)
summary( estimator_sens )
cor(count_vs_sens$Log_counts, count_vs_sens$Sensitivity)

temp_var <- predict(estimator_sens, interval="prediction")
count_vs_sens = cbind( count_vs_sens, temp_var )

reg_plot_sens = ggplot(
  count_vs_sens,
  aes(x=Log_counts, y = Sensitivity)
)
reg_plot_sens = reg_plot_sens + geom_point()
reg_plot_sens = reg_plot_sens + geom_line(aes( y = lwr ), color = "red", linetype = "dashed")
reg_plot_sens = reg_plot_sens + geom_line(aes( y = upr ), color = "red", linetype = "dashed")
reg_plot_sens = reg_plot_sens + geom_smooth( method = lm , se=TRUE )
reg_plot_sens = reg_plot_sens + xlab("Log variants per CCL") + ylab("Sensitivity")
reg_plot_sens

# sensitivity sd

count_vs_sens_sd$Log_counts = as.double(as.character(count_vs_sens_sd$Log_counts))
count_vs_sens_sd$Sensitivity = as.double(as.character(count_vs_sens_sd$Sensitivity))
count_vs_sens_sd$Sensitivity_SD = as.double(as.character(count_vs_sens_sd$Sensitivity_SD))

estimator_sens_sd = lm(formula = Sensitivity_SD ~ Log_counts, data = count_vs_sens_sd)
summary( estimator_sens_sd )
cor(count_vs_sens$Log_counts, count_vs_sens$Sensitivity_SD)

temp_var <- predict(estimator_sens_sd, interval="prediction")
count_vs_sens_sd = cbind( count_vs_sens_sd, temp_var )

reg_plot_sens_sd = ggplot(
  count_vs_sens_sd,
  aes(x=Log_counts, y = Sensitivity_SD)
)
reg_plot_sens_sd = reg_plot_sens_sd + geom_point()
reg_plot_sens_sd = reg_plot_sens_sd + geom_line(aes( y = lwr ), color = "red", linetype = "dashed")
reg_plot_sens_sd = reg_plot_sens_sd + geom_line(aes( y = upr ), color = "red", linetype = "dashed")
reg_plot_sens_sd = reg_plot_sens_sd + geom_smooth( method = lm , se=TRUE)
reg_plot_sens_sd = reg_plot_sens_sd + xlab("Log variants per CCL") + ylab("Standard Deviation Sensitivity")
reg_plot_sens_sd

# PPV

estimator_ppv = lm(formula = PPV ~ Log_counts, data = count_vs_ppv)
summary( estimator_ppv )

temp_var <- predict(estimator_ppv, interval="prediction")
count_vs_ppv = cbind( count_vs_ppv, temp_var )

reg_plot_ppv = ggplot(
  count_vs_ppv,
  aes(x=Log_counts, y = PPV)
)
reg_plot_ppv = reg_plot_ppv + geom_point()
reg_plot_ppv = reg_plot_ppv + geom_line(aes( y = lwr ), color = "red", linetype = "dashed")
reg_plot_ppv = reg_plot_ppv + geom_line(aes( y = upr ), color = "red", linetype = "dashed")
reg_plot_ppv = reg_plot_ppv + geom_smooth( method = lm , se = TRUE )
reg_plot_ppv = reg_plot_ppv + xlab("Log variants per CCL") + ylab("Standard Deviation Sensitivity")
reg_plot_ppv + theme( legend() )

library("cowplot")

legend_b = get_legend( reg_plot_sens_sd )
prow = plot_grid(
  reg_plot_sens + theme(legend.position="none") ,
  reg_plot_sens_sd + theme(legend.position="none") , 
  labels=c("A", "B"),
  ncol = 2,
  nrow = 1
)
p = plot_grid( 
  #legend_b,
  prow,
  ncol = 1, rel_heights = c( .05, .95)
) + ggtitle("Test")
p
