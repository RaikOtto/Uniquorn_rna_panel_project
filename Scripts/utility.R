build_path_variables = function( 
  inclusion_weight = inclusion_weight, 
  only_first = only_first, 
  exclude_self = exclude_self, 
  run_identification = run_identification, 
  cellminer = F,
  type_benchmark
){
  
  anchor                        <<- "~/Uniquorn_data/benchmark_vcf_files/"
  

  if ( type_benchmark == "regularized"){
      
      ident_result_files_path       <<- paste( anchor, "ident_files_regularized", sep ="/" )
      benchmark_res_file_path       <<- paste( anchor, "Benchmark_results_regularized/", sep ="/" )
  }
  if ( type_benchmark == "non_regularized"){
      
      ident_result_files_path       <<- paste( anchor, "ident_files_non_regularized", sep ="/" )
      benchmark_res_file_path       <<- paste( anchor, "Benchmark_results_non_regularized/", sep ="/" )
  }

    benchmark_res_file_path <<- paste(
        benchmark_res_file_path,
        paste(
            as.character( inclusion_weight ),
            "Benchmark_comparisons_result.tab"
            , sep = "_"
        ),
        sep = "/"
    )
    ident_result_files_path <<- paste(
        ident_result_files_path,
        paste0(
            as.character( inclusion_weight ),
            "/"
        ),
        sep = "/"
    )
    
    if ( ! dir.exists(ident_result_files_path  ) )
        dir.create( ident_result_files_path )
        
  
  raw_files_path                <<- paste( anchor, "raw_files", sep ="/" )
  
  known_id_file_path            <<- "~/Uniquorn_rna_panel_project/Misc/known_relationships.tsv"
  known_id_table                <<- read.table( known_id_file_path, sep ="\t", header = TRUE, fill = TRUE)
  id_pairs                      <<- c( paste( known_id_table$Cl1, known_id_table$Cl2, sep ="_" ), paste( known_id_table$Cl2, known_id_table$Cl1, sep ="_" ) )

  benchmark_ident_file_path <<- str_replace( benchmark_res_file_path, pattern = "Benchmark_comparisons_result", "Benchmark_identification_result")
  input_path_comparison <<- str_replace( benchmark_res_file_path, pattern = "comparisons", "identification")
  
}

tell_same = function( vec, symbol ){
  
  input = str_replace_all( vec, symbol, ""  )
  
  for (i in seq(1,length(vec))){
    for (j in seq(1,length(vec))){
      if (i != j){
        
        if ( input[i] == input[j]  ){
          print(c(input[i],input[j], as.character(i),as.character(j)))
        }
      }
    }
  }
}

parse_string = function( replace_name, source_file, query_name = FALSE ){
  
  replace_name = as.character( unlist( sapply( replace_name, FUN = str_replace, "T-T.vcf", "TT2.vcf" ) ) )
  replace_name = as.character( unlist( sapply( replace_name, FUN = str_replace, "T-T_", "TT2_" ) ) )
  replace_name = as.character( unlist( sapply( replace_name, FUN = str_replace, ".vcf_uniquorn_ident.tab", "" ) ) )
  replace_name = as.character( unlist( sapply( replace_name, FUN = str_replace_all, c("\\("), "" ) ) )
  replace_name = as.character( unlist( sapply( replace_name, FUN = str_replace_all, c("\\)"), "" ) ) )
  replace_name = as.character( unlist( sapply( replace_name, FUN = str_trim ) ) )
  replace_name = as.character( unlist( sapply( replace_name, FUN = str_replace_all, c("\\."), "" ) ) )
  replace_name = as.character( unlist( sapply( replace_name, FUN = str_replace_all, c("-"), "" ) ) )
  replace_name = as.character( unlist( sapply( replace_name, FUN = str_replace_all, c("_"), "" ) ) )
  replace_name = as.character( unlist( sapply( replace_name, FUN = str_replace_all, c(":"), "" ) ) )
  replace_name = as.character( unlist( sapply( replace_name, FUN = str_replace_all, c("/"), "" ) ) )
  
  replace_name = as.character( unlist( sapply( replace_name, FUN = str_to_upper ) ) )
  
  if ( source_file != "" )
    replace_name = as.character( unlist( sapply( replace_name, FUN = str_replace, source_file, "" ) ) )
  
  if ( query_name )
    replace_name = tail(replace_name, length(replace_name))
  
  return( as.character( replace_name ) )
}

build_tables = function(){
  
  res_table <<- data.frame( 
    "Query_CCL"             = as.character(),
    "Library_CCL"          = as.character(),
    "CCL_to_be_found"           = as.character(),
    "CCL_passed_threshold"      = as.character(),
    "Found_muts_abs"            = as.character(),
    "P_value"                   = as.character()
  )
  
  res_ident_table <<- data.frame( 
    "Query"                     = as.character( ),
    "Expected"                  = as.character( ),
    "Found"                     = as.character( ),
    "True_positive"             = as.character( ),
    "False_negative"            = as.character( ),
    "False_positive"            = as.character( )
  )
  
  auc_table <<- data.frame(
    "P_values" = as.double(),
    "Should_be_found" = as.logical()
  )
}

concat_me = function(to_be_concated){
  
  if (length(to_be_concated) == 0 ){
    to_be_concated = ""
  } else {
    to_be_concated = paste0( to_be_concated, collapse = ", ")
  }
  
  return(to_be_concated)
}


turn_into_ones_zero_pred = function( vec ){
  
  vec = as.integer( as.vector( unclass( vec ) ) )
  print(vec)
  
  estimate_true_neg_tmp = estimate_true_neg - as.integer( vec[2] ) -  as.integer( vec[3] ) - as.integer( vec[4] )
  
  return(
    c(
      rep( 1, as.integer( vec[2] ) ),
      rep( 1, as.integer( vec[3] ) ),
      rep( 0, as.integer( vec[4] ) ),
      rep( 0, estimate_true_neg )
    )
  )
}

turn_into_ones_zero_labs = function( vec ){
  
  vec = as.integer( as.vector( unclass( vec ) ) )
  print(vec)
  
  estimate_true_neg_tmp = estimate_true_neg - as.integer( vec[2] ) -  as.integer( vec[3] ) - as.integer( vec[4] )
  
  return(
    c(
      rep( 1, as.integer( vec[2] ) ),
      rep( 0, as.integer( vec[3] ) ),
      rep( 0, as.integer( vec[4] ) ),
      rep( 0, estimate_true_neg )
    )
  )
}