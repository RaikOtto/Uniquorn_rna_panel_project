library("stringr")
library("Uniquorn")

inclusion_weight      = 0.5
only_first            = FALSE
exclude_self          = FALSE
run_identification    = F
cellminer_v2          = F
distinct_mode         = TRUE
identical_mode        = FALSE
distuinguished_panels = T
panel_mode = T
#type_benchmark        = 'non_regularized'
type_benchmark        = 'regularized'

run_benchmark_endgegner = function( 
  rigid_similarity = FALSE,
  inclusion_weight = inclusion_weight,
  only_first = only_first,
  exclude_self = exclude_self,
  run_identification = run_identification,
  cellminer = F,
  distuinguished_panels = distuinguished_panels,
  type_benchmark = type_benchmark
){
  
    source("~//Uniquorn_rna_panel_project//Scripts/utility.R")
    build_path_variables( 
        inclusion_weight = inclusion_weight,
        only_first = only_first,
        exclude_self = exclude_self,
        run_identification = run_identification,
        cellminer = F,
        type_benchmark = type_benchmark
    )
    gold_t = read.table( file = "~//Uniquorn_rna_panel_project//Misc//Goldstandard.tsv",sep="\t", header = TRUE)

    
    if (panel){
        ident_result_files_path = str_replace(ident_result_files_path,pattern = "ident_files","panel_ident_files")
        benchmark_ident_file_path = str_replace(benchmark_ident_file_path,pattern = "Benchmark_results","panel_Benchmark_results")
        benchmark_res_file_path = str_replace(benchmark_res_file_path,pattern = "Benchmark_results","panel_Benchmark_results")
    }
    ## prep phase
  
    build_tables()
    
    res_table = read.table( benchmark_res_file_path,   sep ="\t", header = T)
  
    to_be_found_cls = c( gold_t$Name_identical, gold_t$Merged )
    
    # identification table
  
    to_be_found_cls = to_be_found_cls[ to_be_found_cls != "" ]
  
    found = b_table$CCL[ b_table$Identification_sig ]
    found = paste( found, b_table$CL_source[ b_table$Passed_threshold ], sep ="_" )
  
    true_pos_ident  = to_be_found_cls[ which( to_be_found_cls %in% found )  ]
    false_neg_ident = to_be_found_cls[ which( !(to_be_found_cls %in% found) )  ]
    false_pos_ident = found[ which( !(found %in% to_be_found_cls) )  ]
  
  if( length( to_be_found_cls ) == 0 ){
    
    if ( length( false_neg_ident ) == 0 ) 
      
      identification_successful = TRUE
    
    else
      
      identification_successful = FALSE
    
    }else if ( length(true_pos_ident) > 0 ){
      
      identification_successful = TRUE
      
    } else {
      
      identification_successful = FALSE
    }
  
  to_be_found_cls = concat_me( to_be_found_cls  )
  true_pos_ident  = concat_me( true_pos_ident  )
  false_pos_ident = concat_me( false_pos_ident )
  false_neg_ident = concat_me( false_neg_ident )
  found           = concat_me( found  )
  
  res_ident_table <<- data.frame( 
    "Name_query"                = c( as.character( res_ident_table$Name_query ),               as.character( b_name ) ),
    "Source_query"              = c( as.character( res_ident_table$Source_query  ),            as.character( source_file[1]   ) ),
    "Expected"                  = c( as.character( res_ident_table$Expected ),                 as.character( to_be_found_cls  ) ),
    "Found"                     = c( as.character( res_ident_table$Found ),                    as.character( found ) ),
    "True_positive"             = c( as.character( res_ident_table$True_positive ),            as.character( true_pos_ident  ) ),
    "False_negative"            = c( as.character( res_ident_table$False_negative ),           as.character( false_neg_ident  ) ),
    "False_positive"            = c( as.character( res_ident_table$False_positive ),           as.character( false_pos_ident   ) ),
    "Identification_successful" = c( as.character( res_ident_table$Identification_successful ),as.character( identification_successful ) )
  )
  
  # output

  res_table = res_table[ order( as.double( as.character(res_table$Found_muts_weighted_rel) ), decreasing = T),  ]

  print("Writing xlsx")
}

run_small_statistics = function( 
  input_path_comparison = input_path_comparison,
  rigid_similarity = rigid_similarity, 
  inclusion_weight = inclusion_weight, 
  only_first = only_first, 
  exclude_self = exclude_self,
  run_identification = run_identification, 
  panel_mode,
  distuinguished_panels = distuinguished_panels,
  type_benchmark = type_benchmark
){
  
    source("~//Uniquorn_rna_panel_project//Scripts/utility.R")
    source("~/Uniquorn_rna_panel_project//Scripts/Benchmark_Merger_Skript_Relaxed.R")
    build_path_variables( 
        inclusion_weight = inclusion_weight, 
        only_first = only_first, 
        exclude_self = exclude_self, 
        run_identification = run_identification, 
        cellminer = F,
        type_benchmark = type_benchmark
    )
    
    print( c("Inclusion weight", inclusion_weight ) )
  
    if (type_benchmark == "regularized"){
        
        input_path_comparison = paste(c(
            "~/Uniquorn_data/benchmark_vcf_files/Benchmark_results_regularized/",
            inclusion_weight
            ,"_Benchmark_comparisons_result.tab"), sep ="", collapse = "")
        input_path_identification = paste( c(
          "~/Uniquorn_data/benchmark_vcf_files/Benchmark_results_regularized/",
          inclusion_weight,
          "_Benchmark_identification_result.tab"), sep ="", collapse = ""
        )
    }
    
    if (type_benchmark == "non_regularized"){
      
        input_path_comparison = paste(c(
          "~/Uniquorn_data/benchmark_vcf_files/Benchmark_results_non_regularized/",
          inclusion_weight
          ,"_Benchmark_comparisons_result.tab"), sep ="", collapse = "")
        input_path_identification = paste( c(
          "~/Uniquorn_data/benchmark_vcf_files/Benchmark_results_non_regularized/",
          inclusion_weight,
          "_Benchmark_identification_result.tab"), sep ="", collapse = ""
        )
    }
    
    output_table_path = str_replace( 
      input_path_comparison, 
      pattern = "_Benchmark_comparisons_result.tab",
      "_Benchmark_identification_result_aggregated.tab"
    )
    
    if (panel_mode){
        ident_result_files_path = str_replace(ident_result_files_path,pattern = "ident_files","panel_ident_files")
        benchmark_ident_file_path = str_replace(benchmark_ident_file_path,pattern = "Benchmark_results","panel_Benchmark_results")
        benchmark_res_file_path = str_replace(benchmark_res_file_path,pattern = "Benchmark_results","panel_Benchmark_results")
        output_table_path = str_replace(output_table_path,pattern = "Benchmark_results","panel_Benchmark_results")
        input_path_identification = str_replace(input_path_identification,pattern = "Benchmark_results","panel_Benchmark_results")
        
    }
    
    run_relaxed_merging( 
        input_path_identification,
        output_table_path 
    )

    t_rel = read.table( output_table_path, sep ="\t", header = T)
  
    expected = t_rel$Expected[t_rel$Expected!= ""]
    nr_expected = length((as.character(unlist(str_split(expected,", ")))))
    print(c("Expected:",nr_expected) )
    
    true_positives = as.character(unlist(str_split(t_rel$True_positive,", ")))
    true_positives = true_positives[ true_positives != ""  ]
    true_positives = true_positives[ !is.na(true_positives) ]
    nr_true_positives = length(true_positives)
    print( c("TP:",as.integer(nr_true_positives), round(nr_true_positives / nr_expected * 100,1)))
    
    false_negatives = as.character(unlist(str_split(t_rel$False_negative,", ")))
    false_negatives = false_negatives[ false_negatives!= ""]
    false_negatives = false_negatives[ !is.na(false_negatives) ]
    nr_false_negatives = length( false_negatives )
    print(c("FN:",as.integer(nr_false_negatives), round(nr_false_negatives / nr_expected * 100,1)))
    
    false_positives = t_rel$False_positive[t_rel$False_positive!= ""]
    nr_false_positive = length((as.character(unlist(str_split(false_positives,",")))))
    print(c("FP",nr_false_positive))
  
}
run_small_statistics(
    inclusion_weight = inclusion_weight,
    only_first = only_first,
    exclude_self = exclude_self,
    run_identification = run_identification,
    panel_mode = panel_mode,
    type_benchmark = type_benchmark
)
