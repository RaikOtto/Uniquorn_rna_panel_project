

build_path_variables = function( 
  inclusion_weight = inclusion_weight, 
  only_first = only_first, 
  exclude_self = exclude_self, 
  similarity_threshold = similarity_threshold, 
  run_identification = run_identification, 
  cellminer = F, 
  distuinguished_panels = T ){
  
  #anchor                        <<- "/Users/raikotto/Dropbox/PhD/Uniquorn_project/benchmark_vcf_files//"
  anchor                        <<- "/local/ottoraik/benchmark_ident"
  benchmark_res_file_path       <<- paste( anchor, "Benchmark_comparisons_result", sep ="/" )
  raw_files_path                <<- paste( anchor, "raw_files", sep ="/" )
  
  known_id_file_path            <<- paste( anchor, "known_relationships.tab", sep ="/" )
  known_id_table                <<- read.table( known_id_file_path, sep ="\t", header = TRUE, fill = TRUE)
  id_pairs                      <<- c( paste( known_id_table$Cl1, known_id_table$Cl2, sep ="_" ), paste( known_id_table$Cl2, known_id_table$Cl1, sep ="_" ) )
  
  ident_result_files_path       <<- paste( anchor, "ident_files", sep ="/" )
  
  if(!dir.exists(ident_result_files_path))
    dir.create(ident_result_files_path)
  
  ident_result_files_path       <<- paste( ident_result_files_path,as.character( inclusion_weight ), sep = "/")
  benchmark_res_file_path       <<- paste( benchmark_res_file_path, as.character( inclusion_weight ) ,sep = "_")
  
  if(!dir.exists(ident_result_files_path))
    dir.create(ident_result_files_path)
  
  if (only_first){
    benchmark_res_file_path       <<- paste( benchmark_res_file_path, "only_first" ,sep = "_")
  } else {
    benchmark_res_file_path       <<- paste( benchmark_res_file_path, "not_only_first" ,sep = "_")
  }
  
  if(!dir.exists(ident_result_files_path))
    dir.create(ident_result_files_path)
  
  if (exclude_self){
    benchmark_res_file_path       <<- paste( benchmark_res_file_path, "exclude_self" ,sep = "_")
  } else {
    benchmark_res_file_path       <<- paste( benchmark_res_file_path, "not_exclude_self" ,sep = "_")
  }
  
  if (distuinguished_panels){
    benchmark_res_file_path       <<- paste( benchmark_res_file_path, "distinguished_panels" ,sep = "_")
    ident_result_files_path       <<- paste( ident_result_files_path,"distinguished_panels", sep = "/")
    
  } else {
    benchmark_res_file_path       <<- paste( benchmark_res_file_path, "not_distinguished_panels" ,sep = "_")
    ident_result_files_path       <<- paste( ident_result_files_path,"dont_distinguished_panels", sep = "/")
    
  }
  
  if(!dir.exists(ident_result_files_path))
    dir.create(ident_result_files_path)
  
  if(!dir.exists(ident_result_files_path))
    dir.create(ident_result_files_path)
  
  
  if ( identical_mode )
    benchmark_res_file_path <<- paste( benchmark_res_file_path, "_Rigid_results.tab", sep ="_" )
  else 
    benchmark_res_file_path <<- paste( benchmark_res_file_path, "_Relaxed_results.tab", sep ="_" )
  
  benchmark_ident_file_path <<- str_replace(benchmark_res_file_path, "Benchmark_comparisons_result", "Benchmark_identification_result")
  #benchmark_res_file_path_xls  <<- str_replace(benchmark_res_file_path, ".tab",".xlsx") 
  
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
  
  # replace_name = as.character( unlist( sapply( replace_name_input, FUN = str_split, "/" ) ) )
  #replace_name[ which( replace_name == "TT_CCLE.vcf_uniquorn_ident.tab"  ) ] = "TTalt1"
  #replace_name[ which( replace_name == "TT_COSMIC.vcf_uniquorn_ident.tab"  ) ] = "TTalt1"
  #replace_name[ which( replace_name == "T-T_COSMIC.vcf_uniquorn_ident.tab"  ) ] = "TTalt2"
  #replace_name[ which( replace_name == "KMH-2_COSMIC.vcf_uniquorn_ident.tab"  ) ] = "KMH2alt1"
  #replace_name[ which( replace_name == "KM-H2_COSMIC.vcf_uniquorn_ident.tab"  ) ] = "KMH2alt2"
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
    "CL_name_query"             = as.character(),
    "CL_name_training"          = as.character(),
    "Source_query"              = as.character(),
    "Source_training"           = as.character(),
    "Found_muts_abs"            = as.character(),
    #"Count_mutations_abs"      = as.character(),
    "Found_muts_rel"            = as.character(),
    #"Found_muts_weighted"      = as.character(),
    #"Count_mutations_weighted" = as.character(),
    "Found_muts_weighted_rel"   = as.character(),
    "Passed_threshold"          = as.character(),
    "Same_identity"             = as.character(),
    "Same_identity_found"       = as.character(),
    #"Same_identity_found_once"   = as.character(),
    "Related_identity"          = as.character(),
    "Related_identity_found"    = as.character()#,
    #"Related_identity_found_once"= as.character()
  )
  
  res_ident_table <<- data.frame( 
    "Name_query"                = as.character( ),
    "Source_query"              = as.character( ),
    "Expected"                  = as.character( ),
    "Found"                     = as.character( ),
    "True_positive"             = as.character( ),
    "False_negative"            = as.character( ),
    "False_positive"            = as.character( ),
    "Identification_successful" = as.character( )
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