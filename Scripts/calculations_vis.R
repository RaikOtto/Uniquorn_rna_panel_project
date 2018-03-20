library("stringr")
source("~/Uniquorn_rna_panel_project//Scripts/utility.R")
library("readODS")
library("devtools")
setwd("~/Uniquorn_rna_panel_project/Uniquorn/")
load_all()


### load ods benchmark

avg_file   = read_ods("~/Dropbox/Uniquorn_project/Pub/Benchmarks.ods", sheet = 1)
rna_file   = read_ods("~/Dropbox/Uniquorn_project/Pub/Benchmarks.ods", sheet = 2)
TruSight   = read_ods("~/Dropbox/Uniquorn_project/Pub/Benchmarks.ods", sheet = 3)
ClearSight = read_ods("~/Dropbox/Uniquorn_project/Pub/Benchmarks.ods", sheet = 4)
Hotspot    = read_ods("~/Dropbox/Uniquorn_project/Pub/Benchmarks.ods", sheet = 5)

seq_type = factor( 
  c(rep("AVG",16),rep("RNA",16),rep("TruSight",16),rep("ClearSight",16),rep("Hotspot_v2",16)),
  levels = c("AVG","RNA","ClearSight","TruSight","Hotspot_v2"))
Weights = factor(
  as.character( rep( c( rep(1.0,4),rep(0.5,4),rep(0.25,4),rep(0.0,4) ), 5 ) ),
  levels = c("1","0.5","0.25","0")
)
TPs = as.double( c( avg_file[3,-1], rna_file[3,-1], TruSight[3,-1],ClearSight[3,-1],Hotspot[3,-1] ) )

# Sensitivity

Sensitivity = as.double( c( avg_file[6,-1], rna_file[6,-1], TruSight[6,-1],ClearSight[6,-1],Hotspot[6,-1] ) )
seq_type = factor( 
  c(rep("AVG",16),rep("RNA",16),rep("TruSight",16),rep("ClearSight",16),rep("Hotspot_v2",16)),
  levels = c("AVG","RNA","ClearSight","TruSight","Hotspot_v2"))


benchmark_Sensitivity = data.frame(
  "Weight" = Weights,
  "Sensitivity" = Sensitivity,
  "Seq_Type" = seq_type
)

# F1

F1s = as.double( c( avg_file[8,-1], rna_file[8,-1], TruSight[8,-1],ClearSight[8,-1],Hotspot[8,-1] ) )

benchmark_F1 = data.frame(
  "Weight" = Weights,
  "F1" = F1s,
  "Seq_Type" = seq_type
)

# PPV

PPVs = as.double( c( avg_file[7,-1], rna_file[7,-1], TruSight[7,-1],ClearSight[7,-1],Hotspot[7,-1] ) )

benchmark_PPV = data.frame(
  "Weight" = Weights,
  "PPV" = PPVs,
  "Seq_Type" = seq_type
)

### abs_mat

library_names = Uniquorn::read_library_names(ref_gen = "GRCH37")
mat_mat = show_contained_ccls()

abs_mat <<- data.frame(
  "Library" = as.character(),
  "Weight" = as.character(),
  "Count" = as.integer(),
  stringsAsFactors = F
)

for ( var_weight in c("W0","W25","W05","W1")){
  for (library_name in library_names){
    
    lib_mat = mat_mat[ mat_mat$Library == library_name,]
    data_lib = as.integer( lib_mat[[var_weight]] )

    new_mat = data.frame(
      "Library" = rep(library_name, length(data_lib)),
      "Weight"  = rep(var_weight,   length(data_lib)),
      "Count"   = data_lib
    )
    abs_mat   = rbind(abs_mat, new_mat)
  }
}
abs_mat$Library = as.character(abs_mat$Library)
abs_mat$Library[abs_mat$Library=="COSMIC"] = "CPG"
abs_mat$Library[abs_mat$Library=="EGA"] = "Klijn"
abs_mat$Library = factor(abs_mat$Library, levels = c("CCLE","CELLMINER","CPG","Klijn","GDC"))

### mean mat

mean_mat <<- data.frame(
  "Library" = as.character(),
  "Weight" = as.character(),
  "Count" = as.double(),
  stringsAsFactors = F
)
libraries = as.character( unique(as.character(abs_mat$Library ) ))

for ( var_weight in c("W0","W25","W05","W1")){
  for (library_name in libraries){
    
    lib_vec = abs_mat[ abs_mat$Library == library_name,]
    lib_vec = lib_vec[ lib_vec$Weight == var_weight, ]
    
    mean_mat <<- data.frame(
      "Library" = c( as.character( mean_mat$Library ), library_name ),
      "Weight"  = c( as.character( mean_mat$Weight  ), var_weight   ),
      "Count"   = as.double( c( as.double(mean_mat$Count), mean( as.double( lib_vec$Count)) ) )
    )
  }
}

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

all_ccls = ref_ccls_cellminer

all_ccls_uni = unique( as.character(unlist( str_split( all_ccls, pattern = "," ))))

i <<- 0

count_per_ccl = sapply( all_ccls_uni, FUN= function(ref_ccl){
    
    ref_ccl = as.character(ref_ccl)
    ref_ccl = str_replace_all( ref_ccl, pattern = "\\)", "")
    ref_ccl = str_replace_all( ref_ccl, pattern = "\\(", "")
    
    library = as.character( tail( as.character( unlist( str_split(ref_ccl, pattern = "_") ) ), 1 ) )
    ref_ccl_pure = str_replace( ref_ccl, pattern = paste("",library, sep = "_"),""  )
    ref_ccl = str_replace_all( ref_ccl, pattern = "_", "")
    ref_ccl = str_replace_all( ref_ccl, pattern = library, "")
    ref_ccl = paste(ref_ccl, library,sep = "_")
    
    count = as.character( sum ( str_detect( all_ccls, pattern = ref_ccl ) ))

    i <<- i + 1
    print(i)
    
    perf_t = avg_file[ match( ref_ccl, str_to_upper(avg_file$Query), nomatch = 0), ]
    expected = as.character( unlist( str_split( perf_t$Expected, pattern = ",") ) )
    expected = expected[expected != ""]
    expected = length(expected)
    tp = as.character( unlist( str_split( perf_t$True_positive, pattern = ",") ) )
    tp = tp[tp != ""]
    tp = length(tp)
    fp = as.character( unlist( str_split( perf_t$False_positive, pattern = ",") ) )
    fp = fp[fp != ""]
    fp = length(fp)
    fn = as.character( unlist( str_split( perf_t$False_negative, pattern = ",") ) )
    fn = fn[fn != ""]
    fn = length(fn)
    Sensitivity = round( tp / expected * 100, 0 )
    Specificity = round( tp / (tp + fp) * 100, 0 )
    F1 = round( 2 * (Sensitivity * Specificity) / (Sensitivity + Specificity), 0)
    if(F1 == "NaN")
      F1 = "0"
    PPV = round( tp / (tp + fp)  * 100,0)
    
    ref_counts <<- data.frame(
        "CCL" = c(  as.character( ref_counts$CCL) , ref_ccl_pure ),
        "Count" = c( as.character( ref_counts$Count), count),
        "Library" = c( as.character( ref_counts$Library), library),
        "Sensitivity" = c( as.character( ref_counts$Sensitivity), Sensitivity),
        "F1" = c( as.character( ref_counts$F1), F1),
        "PPV" = c( as.character( ref_counts$PPV), PPV),
        "TP" = c( as.character( ref_counts$TP), tp),
        "FP" = c( as.character( ref_counts$FP), fp),
        "FN" = c( as.character( ref_counts$FN), fn)
    )
  }
)
table(ref_counts$Library)
ref_counts[(ncol(ref_counts)-5) : ncol(ref_counts),]
ref_counts[ref_counts$Sensitivity == "NaN",]

write.table(ref_counts,"~/Uniquorn_data/benchmark_vcf_files/Benchmark_Per_CCL_stat.tsv",sep ="\t", quote = F, row.names = F)

####

ref_counts = read.table("~/Uniquorn_data/benchmark_vcf_files/Benchmark_Per_CCL_stat.tsv", sep ="\t", header = T)

ref_counts$F1 = as.character(ref_counts$F1)
ref_counts$F1[ is.na(ref_counts$F1)] = "0"
ref_counts$F1 = as.integer(ref_counts$F1)

ref_counts$Sensitivity = as.integer(as.character(ref_counts$Sensitivity))
ref_counts$F1 = as.integer(as.character(ref_counts$F1))
ref_counts$PPV = as.integer(as.character(ref_counts$PPV))

perf_mat = ref_counts[c("Count","Sensitivity","F1","PPV","TP","FN","FP")]
perf_mat = perf_mat[ref_counts$Sensitivity <100,]
perf_mat$Count = perf_mat$Count / max(perf_mat$Count) * 100

meta_data = as.data.frame(ref_counts[c("CCL","Library")])
meta_data = meta_data[ref_counts$Sensitivity <100,]
rownames(meta_data) = rownames(perf_mat)

cor(
  y = log(perf_mat$Count),
  x = perf_mat$TP
)
summary(perf_mat$Count)

quantile(ref_counts$Count,seq(0,1,by=.1))

chisq.test(
  log(as.integer(ref_counts$Sensitivity)),
  ref_counts$Library
)

cor.test(
  as.integer(perf_mat$F1),
  log(perf_mat$Count),
  method = "pearson"
)

perf_mat_norm = perf_mat
perf_mat_norm$Count = perf_mat_norm$Count / max(perf_mat_norm$Count) * 100
per_mat_norm = perf_mat_norm[order(as.integer(perf_mat_norm$Sensitivity), decreasing = T),]
#per_mat_norm$Count = log(per_mat_norm$Count)
pheatmap::pheatmap(
  t(perf_mat_norm),
  show_colnames = F,
  annotation_col = meta_data[c("Library")],
  cluster_rows = F,
  cluster_cols = F
)

sen_p = ggplot( data = ref_counts, aes(x = Library, y = Sensitivity))
sen_p = sen_p + geom_boxplot()
sen_p

