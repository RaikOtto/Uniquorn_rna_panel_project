library("stringr")

inclusion_weight     = .5
only_first           = FALSE
exclude_self         = FALSE
p_value = .05
q_value = .05
minimum_matching_mutations = 0
run_identification   = T
cellminer_v2         = F
auc_mode             = FALSE
#confidence_score     = 10.0
#type_benchmark       = "non_regularized"
type_benchmark       = "regularized"

run_benchmark = function( 
  inclusion_weight = inclusion_weight, 
  only_first = only_first, 
  exclude_self = exclude_self,
  run_identification = run_identification, 
  cellminer = F, 
  #confidence_score = confidence_score,
  minimum_matching_mutations = minimum_matching_mutations,
  type_benchmark
  ){
  
    source("~/Dropbox/PhD/Uniquorn_project/Scripts/utility.R")
    
    build_path_variables( 
        inclusion_weight = inclusion_weight,
        only_first = only_first,
        exclude_self = exclude_self,
        run_identification = run_identification,
        cellminer = F,
        type_benchmark = type_benchmark
    )
  gold_t = read.table( file = "~/Dropbox/PhD/Uniquorn_project/Pub/Goldstandard.tab",sep="\t", header = TRUE)
  
  # identical only

  ### !!! ###
  raw_files_path = "~/Uniquorn_data/vcfHG19/"
  i_files = list.files( raw_files_path, pattern = ".vcf$", full.names = T )
  
  if (run_identification){
    for ( i_file in i_files){
      
      file_name = paste( tail(as.character( unlist(str_split(i_file,"/"))),1) , "uniquorn_ident.tsv", sep ="_")
      ### !!! ###
      ident_result_files_path = "~/Uniquorn_data/vcfHG19/"
      output_file = paste( ident_result_files_path, file_name, sep ="/" )
      
      if( ! file.exists(output_file) ){
        Uniquorn::identify_vcf_file(
          i_file,
          output_file,
          mutational_weight_inclusion_threshold = inclusion_weight,
          only_first = F,
          minimum_matching_mutations = minimum_matching_mutations,
          #confidence_score = confidence_score,
          verbose = TRUE,
          n_threads = 6
        )
      }
    }
  }
  
  # pre process
  
  ### !!! ###
  ident_result_files_path = "/Users/raikotto/Uniquorn_data/vcfHG19/"
  b_files =  list.files(ident_result_files_path , pattern = ".vcf_uniquorn_ident.tsv", full.names = T )
  
  ## benchmark results positive predictions
  
  build_tables()
  
  parse_identification_data = function( b_file ){
    
    print(c(b_file, as.character(which(b_files %in% b_file))))

    b_table = read.table(b_file, sep ="\t", header = T, stringsAsFactors = FALSE)
    
    if ( grepl( "COSMIC", b_file )  )
      source_file = "COSMIC"
    
    if ( grepl( "CCLE", b_file )  )
      source_file = "CCLE"
    
    if ( grepl( "CELLMINER", b_file )  )
      source_file = "CELLMINER"
    
    if ( grepl( "CUSTOM", b_file )  )
        source_file = "CUSTOM"
    
    if (!exists("source_file"))
        source_file = "CUSTOM"

    b_name = tail(unlist(str_split(b_file, pattern = "/")  ),1)
    b_name_source = str_replace( b_name, ".hg19.vcf_uniquorn_ident.tsv", "" )
    b_name = parse_string( b_name, source_file, query_name = T)
    
    ### new
    
    identifier = str_replace(gold_t$CL,"_COSMIC","")
    identifier = str_replace(identifier,"_CCLE","")
    identifier = str_replace(identifier,"_CELLMINER","")
    identifier = str_replace(identifier,"_CUSTOM","")
    identifier = str_replace(identifier,"-","")
    identifier = str_replace(identifier,"_","")
    
    #identical_cls = as.character( gold_t$Name_identical[ gold_t$CL == b_name_source ] )
    identical_cls = as.character(gold_t$CL[ identifier == b_name ])
    
    # end new
    identical_cls = as.character( unlist( str_split( identical_cls, "," ) ) )
    
    ### new
    
    #related_cls   = as.character( gold_t$Related[ gold_t[,1] == b_name_source ] )
    related_cls   = as.character( gold_t$Related[ identifier == b_name ] )
    ### ned new
    
    related_cls   = as.character( unlist( str_split( related_cls, "," ) ) )
    
    match_name_query = rep( b_name, dim(b_table)[1] )
      
    match_name_training = parse_string( as.character( b_table$CL ), "" )

    same_identity       = paste( as.character( b_table$CL ), b_table$CL_source, sep ="_" ) %in% identical_cls
    same_identity_found = same_identity == b_table$Conf_score_sig
    same_identity_found_once = rep("",dim(b_table)[1])
    same_identity_found_once[ 
      match(TRUE, (same_identity==TRUE) & 
      as.logical( b_table$Conf_score_sig ), nomatch = 0) 
    ] = TRUE
      
    if ( 
      ( !( TRUE %in% same_identity_found_once) ) &
      ( !( TRUE %in% same_identity) )
    ){
        
      max_entry = which( b_table$Conf_score == min( b_table$Conf_score )  )[1]
      same_identity_found_once[max_entry] = FALSE # think about it
    }
      if ( sum( identical_cls != "") == 0 )
        same_identity_found_once = rep("",dim(b_table)[1])
      
      related_identity = paste( as.character( b_table$CL ), b_table$CL_source, sep ="_" ) %in% related_cls | same_identity
      related_identity_found = related_identity == b_table$Conf_score_sig
      related_identity_found_once = rep("",dim(b_table)[1])
      related_identity_found_once[ match(TRUE, (related_identity==TRUE) & as.logical( b_table$Conf_score_sig ), nomatch = 0) ] = TRUE
      
      if ( 
        ( !( TRUE %in% related_identity_found_once) ) &
        ( !( TRUE %in% related_identity) )
      ){
        
        max_entry = which( b_table$Found_muts_weighted_rel == min( b_table$Conf_score )  )[1]
        related_identity_found_once[ max_entry ] = FALSE # think about it
      }
      
      if ( sum(related_cls != "") == 0 )
        related_identity_found_once = rep("",dim(b_table)[1])

      source_file = rep( source_file, dim( b_table )[1] )

      res_table <<- data.frame(
        "CL_name_query"            = c( as.character( res_table$CL_name_query ) , as.character( match_name_query ) ),
        "CL_name_training"         = c( as.character( res_table$CL_name_training ) , as.character( match_name_training ) ),
        "Source_query"             = c( as.character( res_table$Source_query ) , as.character( source_file ) ),
        "Source_training"          = c( as.character( res_table$Source_training), as.character( b_table$CL_source ) ),
        "Found_muts_abs"           = c( as.character( res_table$Found_muts_abs), as.character( b_table$Found_muts ) ),
        "P_value"                  = c( as.character( res_table$P_value), as.character( b_table$P_values ) ),
        "Q_value"                  = c( as.character( res_table$Q_value), as.character( b_table$Q_values ) ),
        "Conf_score"               = c( as.character( res_table$Conf_score), as.character( b_table$Conf_score ) ),
        "Passed_threshold"         = c( as.character( res_table$Passed_threshold), as.character( b_table$Conf_score_sig ) ),
        "Same_identity"            = c( as.character( res_table$Same_identity )             , as.character( same_identity ) ),
        "Same_identity_found"      = c( as.character( res_table$Same_identity_found )       , as.character( same_identity_found ) ),
        "Related_identity"         = c( as.character( res_table$Related_identity )          , as.character( related_identity ) ),
        "Related_identity_found"   = c( as.character( res_table$Related_identity_found )    , as.character( related_identity_found ) )
      )
      res_table = res_table[ as.character( res_table$Found_muts_abs ) != "0",]
      
      # identification table
      
      b_table_save = b_table
      
      if (exclude_self)
        b_table = b_table[ paste( b_table$CL, b_table$CL_source, sep = "_" ) != b_name_source, ]
      
      to_be_found_cls = unique( c( identical_cls, related_cls ) ) # in the future: alternate here
      to_be_found_cls = to_be_found_cls[ to_be_found_cls != "" ]

      found = b_table$CL[ b_table$Conf_score_sig ]
      found = paste( found, b_table$CL_source[ b_table$Conf_score_sig ], sep ="_" )
      
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
          "Identification_successful" = c( as.character( res_ident_table$Identification_successful ),as.character( identification_successful ) ),
          stringsAsFactors = F
        )
      
      found = b_table$CL[ b_table$Conf_score_sig  ]
      found = paste( found, b_table$CL_source[ b_table$Conf_score_sig ], sep ="_" )
      
      true_pos_ident  = to_be_found_cls[ which( to_be_found_cls %in% found )  ]
      false_neg_ident = to_be_found_cls[ which( !(to_be_found_cls %in% found) )  ]
      false_pos_ident = found[ which( !(found %in% to_be_found_cls) )  ]
      
      true_pos_ident  = concat_me( true_pos_ident  )
      false_pos_ident = concat_me( false_pos_ident )
      false_neg_ident = concat_me( false_neg_ident )
      found           = concat_me( found  )
      
      index = which( as.character( inclusion_weight ) == auc_table$Sim_Thresh  )
      
      tmp_true_pos = as.character( unlist( str_split(  as.character( true_pos_ident  ), ",") ) )
      tmp_true_pos = tmp_true_pos[ tmp_true_pos != ""  ]
      
      tmp_false_pos = as.character( unlist( str_split( as.character( false_pos_ident  ), ",") ) )
      tmp_false_pos = tmp_false_pos[ tmp_false_pos != ""  ]

      tmp_false_neg = as.character( unlist( str_split(  as.character( false_neg_ident  ), ",") ) )
      tmp_false_neg = tmp_false_neg[ tmp_false_neg != ""  ]
      
      auc_table$False_neg[ index ] = as.character( as.integer( auc_table$False_neg[ index ] ) + length( tmp_false_neg ) )
      write.table( res_ident_table, benchmark_ident_file_path, sep ="\t", quote = F, row.names = F  )
      
      }
      # output
      
      write.table( res_ident_table, benchmark_ident_file_path, sep ="\t", quote = F, row.names = F  )
  sapply( b_files, FUN = parse_identification_data)
  
  res_table = res_table[ order( as.double( as.character(res_table$Found_muts_abs) ), decreasing = T),  ]

}

run_benchmark(
  inclusion_weight,
  only_first,
  exclude_self,
  run_identification = run_identification,
  #confidence_score = confidence_score,
  minimum_matching_mutations = minimum_matching_mutations,
  type_benchmark = type_benchmark
)
