library(stringr)
library(devtools)
load_all()

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

FN = as.character( unlist( str_split( t$False_negative, "," ) ) )
FN = str_trim( FN )
FN = FN[FN!=""]
FN_unique = unique(FN)
c = show_contained_cls()
Not_FN = c$CL[ !( c$CL %in% FN )  ]

nr_var_fn <<- c()
for( fn in FN_unique  ){
  
  vars = sim_list[ sim_list$CL == fn, ]
  var_weight = vars$Fingerprint[  vars$Weight >= weight ]
  nr_var_fn = c(
    nr_var_fn,
    length( var_weight )
  )
  print( c( length( var_weight ), length(nr_var_fn) ) )
}


nr_var_not_fn <<- c()
for( not_fn in Not_FN  ){
    
    vars = sim_list[ sim_list$CL == not_fn, ]
    var_weight = vars$Fingerprint[  vars$Weight >= weight ]
    nr_var_not_fn = c(
        nr_var_not_fn,
        length( var_weight )
    )
    print( c( length( var_weight ), length(nr_var_not_fn) ) )
}

nr_var_not_fn

mean(nr_var_fn)
mean(nr_var_not_fn)
t.test(nr_var_fn, nr_var_not_fn)

ccle_nr_fn = length( grep( FN_unique, pattern = "_CCLE" ) )
cosmic_nr_fn = length( grep( FN_unique, pattern = "_COSMIC" ) )
cellminer_nr_fn = length( grep( FN_unique, pattern = "_CELLMINER" ) )

sort(nr_var_fn,decreasing = F)
sort(nr_var_not_fn,decreasing = F)

FN_unique[ which( nr_var_fn <= 1 ) ]
