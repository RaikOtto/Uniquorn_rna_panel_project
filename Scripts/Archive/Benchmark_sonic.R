args = commandArgs(trailingOnly = TRUE)
library("stringr")
library("Uniquorn")

inclusion_weight     = 1.0
only_first           = FALSE
exclude_self         = FALSE
similarity_threshold = as.integer(args[1])
minimum_matching_mutations = 3
run_identification   = T
cellminer_v2         = F
distinct_mode        = TRUE
identical_mode       = TRUE
auc_mode             = TRUE

if (identical_mode){
  
  auc_file_path = "/local/ottoraik/auc_folder/rigid/"
  dir.create( auc_file_path, showWarnings = FALSE)
  auc_file_path = paste0( auc_file_path, as.character(minimum_matching_mutations), paste = "/" )
  dir.create( auc_file_path, showWarnings = FALSE)
  auc_file_path = paste0( auc_file_path, "auc_", sep = "" )
  auc_file_path = paste0( c(auc_file_path, as.character(similarity_threshold), ".tab") , collapse = "" )
  
} else {
  
  auc_file_path = "/local/ottoraik/auc_folder/relaxed/"
  dir.create( auc_file_path, showWarnings = FALSE)
  auc_file_path = paste0( auc_file_path, as.character(minimum_matching_mutations), paste = "/" )
  dir.create( auc_file_path, showWarnings = FALSE)
  auc_file_path = paste0( auc_file_path, "auc_", sep = "" )
  auc_file_path = paste0( c(auc_file_path, as.character(similarity_threshold), ".tab") , collapse = "" )
}

run_benchmark = function( inclusion_weight = inclusion_weight, only_first = only_first, exclude_self = exclude_self, similarity_threshold = similarity_threshold, run_identification = run_identification, cellminer = F, distuinguished_panels = T, minimum_matching_mutations = 4 ){
  
  source("/local/ottoraik/utility.R")
  #build_path_variables( inclusion_weight = inclusion_weight, only_first = only_first, exclude_self = exclude_self, similarity_threshold = similarity_threshold, run_identification = run_identification, cellminer = F, distuinguished_panels = T )
  gold_t = read.table( file = "/local/ottoraik/Goldstandard.tab",sep="\t", header = TRUE)
  
  # identical only
  if( identical_mode )
    gold_t$Related = rep("",dim(gold_t)[1])
  
    # pre process
  
  b_files =  list.files( "/local/ottoraik/benchmark_ident" , pattern = ".vcf_uniquorn_ident.tab", full.names = T )
  
  ## benchmark results positive predictions
  
  build_tables()
  
  if (! file.exists(auc_file_path) ){
    dummy_tab = data.frame( 
      "Sim_Thresh" = as.character(   ),
      "True_pos"   = as.character(   ),
      #"True_neg"   = as.character(  ),
      "False_pos"   = as.character(  ),
      "False_neg"   = as.character(  )
    )
    write.table( dummy_tab, auc_file_path, sep ="\t", quote = F, row.names = F  )
  }
  
  auc_table = read.table( auc_file_path, sep ="\t", header = T  )
  if ( ! ( as.character(similarity_threshold  ) %in% auc_table$Sim_Thresh  ) ){
    tmp = matrix( c( as.character(similarity_threshold) ,rep("0",3) ), ncol = 4 )
    colnames(tmp) = colnames( auc_table )
    auc_table = rbind(  auc_table, tmp  )
  }
  write.table( auc_table, auc_file_path, sep ="\t", quote = F, row.names = F  )
  
  parse_identification_data = function( b_file ){
    
    print(c(b_file, as.character(which(b_files %in% b_file))))
    
    b_table = read.table(b_file, sep ="\t", header = T, stringsAsFactors = FALSE)
    
    if (auc_mode){
      
      b_table$Passed_threshold =  
        ( as.integer( b_table$Found_muts_abs )         >= minimum_matching_mutations ) &
        ( as.double( b_table$Found_muts_weighted_rel ) >= as.double( similarity_threshold ) )
      
    }
    
    if ( grepl( "COSMIC", b_file )  )
      source_file = "COSMIC"
    
    if ( grepl( "CCLE", b_file )  )
      source_file = "CCLE"
    if ( grepl( "CELLMINER", b_file )  )
      source_file = "CELLMINER"
    if ( grepl( "CUSTOM", b_file )  )
      source_file = "CUSTOM"
    
    b_name = tail(unlist(str_split(b_file, pattern = "/")  ),1)
    b_name_source = str_replace( b_name, ".vcf_uniquorn_ident.tab", "" )
    b_name = parse_string( b_name, source_file, query_name = T)
    
    identical_cls = as.character( gold_t$Name_identical[ gold_t[,1] == b_name_source ] )
    identical_cls = as.character( unlist( str_split( identical_cls, "," ) ) )
    
    related_cls   = as.character( gold_t$Related       [ gold_t[,1] == b_name_source ] )
    related_cls   = as.character( unlist( str_split( related_cls, "," ) ) )
    
    if ( exclude_self  ){
      b_table       = b_table[ as.integer(round(b_table$Found_muts_weighted_rel,0)) != 100, ]
      identical_cls = identical_cls[ identical_cls != b_name_source]
      related_cls   = related_cls[ related_cls != b_name_source]
    }
    
    if ( only_first ){
      b_table$Passed_threshold[ seq(2, length(b_table$Passed_threshold)) ] = F
    }
    
    match_name_query = rep( b_name, dim(b_table)[1] )
    
    match_name_training = parse_string( as.character( b_table$CL ), "" )
    
    same_identity       = paste( as.character( b_table$CL ), b_table$CL_source, sep ="_" ) %in% identical_cls
    same_identity_found = same_identity == b_table$Passed_threshold
    same_identity_found_once = rep("",dim(b_table)[1])
    same_identity_found_once[ match(TRUE, (same_identity==TRUE) & b_table$Passed_threshold, nomatch = 0) ] = TRUE
    
    if ( 
      ( !( TRUE %in% same_identity_found_once) ) &
      ( !( TRUE %in% same_identity) )
    ){
      
      max_entry = which( b_table$Found_muts_weighted_rel == max( b_table$Found_muts_weighted_rel )  )[1]
      same_identity_found_once[max_entry] = FALSE # think about it
    }
    if ( sum( identical_cls != "") == 0 )
      same_identity_found_once = rep("",dim(b_table)[1])
    
    related_identity = paste( as.character( b_table$CL ), b_table$CL_source, sep ="_" ) %in% related_cls | same_identity
    related_identity_found = related_identity == b_table$Passed_threshold
    related_identity_found_once = rep("",dim(b_table)[1])
    related_identity_found_once[ match(TRUE, (related_identity==TRUE) & b_table$Passed_threshold, nomatch = 0) ] = TRUE
    
    if ( 
      ( !( TRUE %in% related_identity_found_once) ) &
      ( !( TRUE %in% related_identity) )
    ){
      
      max_entry = which( b_table$Found_muts_weighted_rel == max( b_table$Found_muts_weighted_rel )  )[1]
      related_identity_found_once[ max_entry ] = FALSE # think about it
    }
    
    if ( sum(related_cls != "") == 0 )
      related_identity_found_once = rep("",dim(b_table)[1])
    
    source_file = rep( source_file, dim( b_table )[1] )
    
    res_table = data.frame(
      "CL_name_query"            = c( as.character( res_table$CL_name_query ) , as.character( match_name_query ) ),
      "CL_name_training"         = c( as.character( res_table$CL_name_training ) , as.character( match_name_training ) ),
      "Source_query"             = c( as.character( res_table$Source_query ) , as.character( source_file ) ),
      "Source_training"          = c( as.character( res_table$Source_training), as.character( b_table$CL_source ) ),
      "Found_muts_abs"           = c( as.character( res_table$Found_muts_abs), as.character( b_table$Found_muts_abs ) ),
      "Found_muts_weighted"      = c( as.character( res_table$Found_muts_weighted ), b_table$Found_muts_weighted ),
      "Found_muts_weighted_rel"  = c( as.character( res_table$Found_muts_weighted_rel ), round( b_table$Found_muts_weighted_rel, 2 ) ),
      "Passed_threshold"         = c( as.character( res_table$Passed_threshold), as.character( b_table$Passed_threshold ) ),
      "Same_identity"            = c( as.character( res_table$Same_identity )             , as.character( same_identity ) ),
      "Same_identity_found"      = c( as.character( res_table$Same_identity_found )       , as.character( same_identity_found ) ),
      "Related_identity"         = c( as.character( res_table$Related_identity )          , as.character( related_identity ) ),
      "Related_identity_found"   = c( as.character( res_table$Related_identity_found )    , as.character( related_identity_found ) )
    )
    res_table = res_table[ as.character( res_table$Found_muts_weighted_rel ) != "0",]
    
    # identification table
    
    b_table_save = b_table
    
    if (exclude_self)
      b_table = b_table[ paste( b_table$CL, b_table$CL_source, sep = "_" ) != b_name_source, ]
    
    to_be_found_cls = unique( c( identical_cls, related_cls ) ) # in the future: alternate here
    to_be_found_cls = to_be_found_cls[ to_be_found_cls != "" ]

    found = b_table$CL[ b_table$Passed_threshold  ]
    found = paste( found, b_table$CL_source[ b_table$Passed_threshold ], sep ="_" )
    
    true_pos_ident  = to_be_found_cls[ which( to_be_found_cls %in% found )  ]
    false_neg_ident = to_be_found_cls[ which( !(to_be_found_cls %in% found) )  ]
    false_pos_ident = found[ which( !(found %in% to_be_found_cls) )  ]
    
    true_pos_ident  = concat_me( true_pos_ident  )
    false_pos_ident = concat_me( false_pos_ident )
    false_neg_ident = concat_me( false_neg_ident )
    found           = concat_me( found  )
    
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
      "Identification_successful" = c( as.character( res_ident_table$Identification_successful ),as.character( identification_successful ) ),
      stringsAsFactors = F)
    
    auc_table = read.table( auc_file_path, sep ="\t", header = T  )
    
    index = which( as.character(similarity_threshold ) == auc_table$Sim_Thresh  )
    
    tmp_true_pos = as.character( unlist( str_split(  as.character( true_pos_ident  ), ",") ) )
    tmp_true_pos = tmp_true_pos[ tmp_true_pos != ""  ]
    auc_table$True_pos[ index ] = as.character( as.integer( auc_table$True_pos[ index ] ) + length( tmp_true_pos ) )
    
    tmp_false_pos = as.character( unlist( str_split( as.character( false_pos_ident  ), ",") ) )
    tmp_false_pos = tmp_false_pos[ tmp_false_pos != ""  ]
    auc_table$False_pos[ index ] = as.character( as.integer( auc_table$False_pos[ index ] ) + length( tmp_false_pos ) )
    
    tmp_false_neg = as.character( unlist( str_split(  as.character( false_neg_ident  ), ",") ) )
    tmp_false_neg = tmp_false_neg[ tmp_false_neg != ""  ]
    auc_table$False_neg[ index ] = as.character( as.integer( auc_table$False_neg[ index ] ) + length( tmp_false_neg ) )
    # output
    
    write.table( auc_table, auc_file_path, sep ="\t", quote = F, row.names = F  )
    #write.table( res_ident_table, benchmark_ident_file_path, sep ="\t", quote = F, row.names = F  )
    #write.table( res_table ,       benchmark_res_file_path, sep ="\t", quote = F, row.names = F  )
  }
  sapply( b_files, FUN = parse_identification_data)
  
  res_table = res_table[ order( as.double( as.character(res_table$Found_muts_weighted_rel) ), decreasing = T),  ]
  #write.table( res_table , benchmark_res_file_path, sep ="\t", quote = F, row.names = F  )
  
}

run_benchmark(
  inclusion_weight,
  only_first,
  exclude_self,
  similarity_threshold = similarity_threshold,
  run_identification = FALSE,
  distuinguished_panels = distuinguished_panels,
  minimum_matching_mutations = minimum_matching_mutations
)

