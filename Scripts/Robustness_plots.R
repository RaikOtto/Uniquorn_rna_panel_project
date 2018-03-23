### Per_ccl_success

avg_file = read.table("~/Uniquorn_data/benchmark_vcf_files/0_5_Benchmark_identification_result.tab", sep = "\t", header = T, stringsAsFactors = F)

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
  "FN" = as.character(),
  "TP_ref" = as.character(),
  "FP_ref" = as.character(),
  "FN_ref" = as.character()
)

all_ccls     = avg_file$Query
all_ccls_uni = unique( avg_file$Query )

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
  
  count = 49#as.character( sum ( str_detect( all_ccls, pattern = ref_ccl ) ))
  
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

ref_counts$Sensitivity[ref_counts$Sensitivity == "NaN"] = "0"
ref_counts$PPV[ref_counts$PPV == "NaN"] = "0"
ref_counts$PPV[is.na(ref_counts$PPV)] = "0"

ref_counts$F1 = 2 * ( as.integer(as.character(ref_counts$Sensitivity) ) * as.integer(as.character(ref_counts$Specificity) ) ) /
  ( as.integer(as.character(ref_counts$Sensitivity) ) + as.integer(as.character(ref_counts$Specificity) ) )
ref_counts$PPV = as.integer(as.character(ref_counts$TP)) /
  ( as.integer(as.character(ref_counts$TP) ) + as.integer(as.character(ref_counts$FP) ) ) * 100

perf_mat = aggregate( as.double(ref_counts$Sensitivity), by = list(ref_counts$Library), FUN = mean)
perf_mat$Sens_SD = aggregate( as.double( ref_counts$Sensitivity), by = list(ref_counts$Library), FUN = sd)[,2]
perf_mat$F1 = aggregate( as.double( ref_counts$F1), by = list(ref_counts$Library), FUN = mean)[,2]
perf_mat$F1_SD = aggregate( as.double( ref_counts$F1 ), by = list(ref_counts$Library), FUN = sd)[,2]

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

# capacity to identify

avg_file = read.table("~/Uniquorn_data/benchmark_vcf_files/0_5_Benchmark_identification_result.tab", sep = "\t", header = T, stringsAsFactors = F)

expected = str_trim(as.character(unlist(str_split(avg_file$Expected, pattern = ","))))
expected = expected[expected!= ""]
lib_expected = as.character(sapply( expected, function(vec){
  return( as.character( tail( as.character( unlist( str_split(vec, pattern = "_") ) ), 1 ) ) )
}))
sum_expected = table(lib_expected)

TPs_avg = str_trim(as.character(unlist(str_split(avg_file$True_positive, pattern = ","))))
TPs_avg = TPs_avg[TPs_avg!= ""]
lib_TPs = as.character(sapply( TPs_avg, function(vec){
    return( as.character( tail( as.character( unlist( str_split(vec, pattern = "_") ) ), 1 ) ) )
}))
sum_TPs = table(lib_TPs)

FNs_avg = str_trim(as.character(unlist(str_split(avg_file$False_negative, pattern = ","))))
FNs_avg = FNs_avg[FNs_avg!= ""]
lib_FNs = as.character(sapply( FNs_avg, function(vec){
  return( as.character( tail( as.character( unlist( str_split(vec, pattern = "_") ) ), 1 ) ) )
}))
sum_FNs = table(lib_FNs)

FPs_avg = str_trim(as.character(unlist(str_split(avg_file$False_positive, pattern = ","))))
FPs_avg = FPs_avg[FPs_avg!= ""]
lib_FPs = as.character(sapply( sum_FPs, function(vec){
  return( as.character( tail( as.character( unlist( str_split(vec, pattern = "_") ) ), 1 ) ) )
}))
sum_FPs = table(lib_FPs)

perf_mat <<- data.frame(
    "CCL" = as.character(),
    "Library" = as.character(),
    "Sensitivity" = as.character(),
    "F1" = as.character(),
    "PPV" = as.character(),
    "TPs" = as.character(),
    "FPs" = as.character(),
    "FNs" = as.character()
)

i <<- 0
sapply(all_ccls_uni, FUN = function(CCL){
  
  lib = tail( as.character(unlist(str_split(CCL,pattern ="_"))),1 )
  CCL = str_to_upper(CCL)
  CCL = str_replace( CCL, pattern = paste0("_", lib), "")
  CCL = str_replace_all( CCL, pattern = "_", "")
  CCL = str_replace_all( CCL, pattern = ":", "")
  CCL = str_replace_all( CCL, pattern = "\\(", "")
  CCL = str_replace_all( CCL, pattern = "\\)", "")
  CCL = str_replace_all( CCL, pattern = "-", "")
  CCL = paste(CCL, lib, sep = "_")
  
  i <<- i + 1
  print(i)
  
  TPs = sum( str_detect( TPs_avg, pattern = CCL) )
  FPs = sum( str_detect( FPs_avg, pattern = CCL) )
  FNs = sum( str_detect( FNs_avg, pattern = CCL) )
  
  Sensitivity = round(( TPs / (TPs + FNs) * 100), 1)
  F1 = round( 2 * (sensitivity * (3058 - FPs)) / (sensitivity + 3058 - FPs), 1)
  PPV = round( TPs / ( TPs + FPs)  * 100,1)

  
  perf_mat <<- data.frame(
    "CCL" = c( as.character( perf_mat$CCL), as.character(CCL)),
    "Library" = c( as.character( perf_mat$Library), as.character( lib)),
    "Sensitivity" = c( as.character( perf_mat$Sensitivity), as.character(Sensitivity)),
    "F1" = c( as.character( perf_mat$F1), as.character(F1)),
    "PPV" = c( as.character( perf_mat$PPV), as.character(PPV)),
    "TPs" = c( as.character( perf_mat$TPs), as.character(TPs)),
    "FPs" = c( as.character( perf_mat$FPs), as.character(FPs)),
    "FNs" = c( as.character( perf_mat$FNs), as.character(FNs))
  )
})
perf_mat = perf_mat[ perf_mat$Sensitivity != "NaN", ]

sensitivity = round( sum_TPs / (sum_TPs + sum_FNs) * 100, 1 )
F1 = round( 2 * (sensitivity * 99 ) / (sensitivity + 99), 1)
PPV = round( sum_TPs / ( sum_TPs + sum_FPs)  * 100,1)

F1_sd = aggregate( as.double(perf_mat$F1), by = list(perf_mat$Library), FUN = sd)
PPV_sd = aggregate( as.double(perf_mat$PPV), by = list(perf_mat$Library), FUN = sd)

vis_mat = data.frame(
    "Library" = names(table(lib_FPs)),
    "Sensitivity" = as.double( as.character( sensitivity ) ),
    "Sensitivity_SD" = as.double( as.character( sensitivity_sd[,2] ) ),
    "F1" = as.double( as.character( F1) ),
    "F1_SD" = as.double( as.character( F1_sd[,2]) ),
    "PPV" = as.double( as.character( PPV ) ),
    "PPV_SD" = as.double( as.character( PPV_sd[,2] ) )
)

vis_mat_melt = reshape2::melt(vis_mat)
colnames( vis_mat_melt ) = c("Library","Parameter", "Value")
vis_mat_melt$Library = as.character(vis_mat_melt$Library)
vis_mat_melt$Library[as.character( vis_mat_melt$Library) == "COSMIC"] = "CGP"
vis_mat_melt$Library[as.character(vis_mat_melt$Library) == "EGA"] = "Klijn"
vis_mat_melt$Library = factor( vis_mat_melt$Library, levels = c("CGP","GDC","CCLE","Klijn","CELLMINER") )

means = c(vis_mat_melt$Value[vis_mat_melt$Parameter %in% c("Sensitivity", "F1","PPV")])
sds = c(vis_mat_melt$Value[!( vis_mat_melt$Parameter %in% c("Sensitivity", "F1","PPV"))])
sds_min = means - sds / 2
sds_max = means + sds / 2

g_bench = ggplot(
  subset(vis_mat_melt,Parameter %in% c("Sensitivity","F1","PPV")),
  aes( x = Parameter, y = Value, fill = Library))
g_bench = g_bench + geom_bar(stat="identity", position=position_dodge())
g_bench = g_bench +  geom_errorbar(aes( ymin = sds_min, ymax = sds_max), width=.2,position=position_dodge(.9))
g_bench = g_bench + theme(legend.position="top", legend.background = element_rect(fill="gray90", size =5 )  )
g_bench + xlab("")
