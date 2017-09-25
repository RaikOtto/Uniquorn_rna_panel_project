library(stringr)
library(devtools)
library(Uniquorn)

t = read.table(
    "~/Uniquorn_data/benchmark_vcf_files/Benchmark_results_regularized/1_Benchmark_identification_result.tab",
    header = T,
    fill = T,
    sep ="\t",
    stringsAsFactors = F
)


weight = 1.0
ref_gen = "GRCH37"

sim_list = initiate_db_and_load_data( 
    ref_gen = ref_gen, 
    request_table = "sim_list"
)

FP = as.character( unlist( str_split( t$False_positive , "," ) ) )
FP = str_trim( FP )
FP = FP[FP!=""]

FP2 = as.character( unlist( str_replace(FP, "_COSMIC","") ) )
FP2 = as.character( unlist( str_replace(FP2, "_CCLE","") ) )
FP2 = as.character( unlist( str_replace(FP2, "_CELLMINER","") ) )
FP2 = as.character( unlist( str_replace(FP2, "_","") ) )


FP_unique = unique(FP)
c = show_contained_cls()
Not_FP = c$CL[ !( c$CL %in% FP )  ]

nr_var_fp <<- c()
for( fp in FP_unique  ){
    
    vars = sim_list[ sim_list$CL == fp, ]
    var_weight = vars$Fingerprint[  vars$Weight >= weight ]
    nr_var_fp = c(
        nr_var_fp,
        length( var_weight )
    )
    print( c( length( var_weight ), length(nr_var_fp) ) )
}

nr_var_not_fp <<- c()
for( not_fp in Not_FP  ){
    
    vars = sim_list[ sim_list$CL == not_fp, ]
    var_weight = vars$Fingerprint[  vars$Weight >= weight ]
    nr_var_not_fp = c(
        nr_var_not_fp,
        length( var_weight )
    )
    print( c( length( var_weight ), length(nr_var_not_fp) ) )
}

mean(nr_var_fp)
mean(nr_var_not_fp)
t.test(nr_var_fp, nr_var_not_fp)

ccle_nr_fp = length( grep( FP_unique, pattern = "_CCLE" ) )
cosmic_nr_fp = length( grep( FP_unique, pattern = "_COSMIC" ) )
cellminer_nr_fp = length( grep( FP_unique, pattern = "_CELLMINER" ) )

sort(nr_var_fp,decreasing = F)
sort(nr_var_not_fp,decreasing = F)

FP_unique[ which( nr_var_fp <= 1 ) ]
