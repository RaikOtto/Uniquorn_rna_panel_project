library("stringr")
library("BiocParallel")
library("doParallel")
library("foreach")
library("devtools")
setwd("~/Uniquorn_rna_panel_project/Uniquorn/")
load_all()

inclusion_weight     = .5
only_first           = FALSE
exclude_self         = FALSE
p_value = .05
minimum_matching_mutations = 0
run_identification   = F
panel_mode = F
auc_mode             = FALSE
ref_gen              = "GRCH37"
#type_benchmark       = "non_regularized"
panel_mode = TRUE
type_benchmark       = "regularized"
n_threads = 4

run_benchmark = function(
    inclusion_weight = inclusion_weight, 
    only_first = only_first,
    exclude_self = exclude_self,
    run_identification = run_identification, 
    #confidence_score = confidence_score,
    minimum_matching_mutations = minimum_matching_mutations,
    type_benchmark,
    ref_gen
){
  
    source("~/Uniquorn_rna_panel_project/Scripts/utility.R")
    
    build_path_variables( 
        inclusion_weight = inclusion_weight,
        only_first = only_first,
        exclude_self = exclude_self,
        run_identification = run_identification,
        cellminer = F,
        type_benchmark = type_benchmark
    )
    out_path_ident = str_replace(
        benchmark_ident_file_path,
        pattern = "/Benchmark_results_regularized//0.5_Benchmark_identification_result.tab",
        replacement = paste(c("ident_files_regularized",inclusion_weight,""),sep="",collapse= "/")
    )
    gold_t = read.table( file = "~//Uniquorn_rna_panel_project/Misc/Goldstandard.tsv",sep="\t", header = TRUE)
    
    if (panel_mode){
        ident_result_files_path = str_replace(ident_result_files_path,pattern = "ident_files","panel_ident_files")
        benchmark_ident_file_path = str_replace(benchmark_ident_file_path,pattern = "Benchmark_results","panel_Benchmark_results")
        benchmark_res_file_path = str_replace(benchmark_res_file_path,pattern = "Benchmark_results","panel_Benchmark_results")
    }

    b_files =  list.files(ident_result_files_path , pattern = ".ident.tsv", full.names = T )
    
    ## benchmark results positive predictions
    
    build_tables()
  
    for (b_file in b_files)
        parse_identification_data(b_file)

    #sapply( b_files, FUN = parse_identification_data)
  
    res_table = res_table[ order( as.double( as.character(res_table$Found_muts_abs) ), decreasing = T),  ]
    write.table( res_ident_table, benchmark_ident_file_path, sep ="\t", quote = F, row.names = F  )
    write.table( res_table, benchmark_res_file_path, sep ="\t", quote =F, row.names = F )
}

parse_identification_data = function( b_file ){
  
    print(c(b_file, as.character(which(b_files %in% b_file))))
    
    b_table = read.table(b_file, sep ="\t", header = T, stringsAsFactors = FALSE)
    library_names = Uniquorn::read_library_names(ref_gen = ref_gen)
    
    for (library_name in library_names)
      if ( str_detect(pattern = library_name, b_file) )
        source_file <<- library_name
    
    b_name = tail(unlist(str_split(b_file, pattern = "/")  ),1)
    b_name = str_replace( b_name, pattern = paste( c(".",source_file,".ident.tsv"), sep ="", collapse = ""),"")
    b_name_with_library = paste( b_name[1], source_file[1], sep ="_")
    b_name_plane = str_to_upper( str_replace_all(b_name, pattern = "[\\(\\*\\-\\_\\+\\)]",""))
    
    b_name_ccls = as.character(b_table$CCL)
    b_name_ccls_plane = str_to_upper( as.character(sapply( b_name_ccls, FUN = str_replace_all, pattern = "[\\(\\*\\-\\_\\+\\)]","")) )

    ### new
    
    identifier <<- gold_t$CL
    for (library_name in library_names)
        identifier <<- str_replace( identifier, pattern = paste( "_",library_name,sep=""), "")
    
    identifier_plane = str_to_upper( str_replace_all(identifier, pattern = "[\\(\\*\\-\\_\\+\\)]","") )
    identical_cls = unique( identifier_plane[ which(identifier_plane == b_name_plane) ] )
    same_identity = identical_cls == b_name_ccls_plane
    
    related_cls = as.character( gold_t$Related[ which(identifier_plane == b_name_plane) ] )
    related_cls = as.character( unlist( str_split( related_cls, "," ) ) )
    related_cls = related_cls[ related_cls != ""]
    
    match_name_query = rep( b_name, dim(b_table)[1] )
    match_name_training = identifier
    
    same_identity_found = same_identity == b_table$Identification_sig
    same_identity_found_once = rep("",dim(b_table)[1])
    same_identity_found_once[ match(TRUE, same_identity & b_table$Identification_sig, nomatch = 0) ] = TRUE
    
    #if ( 
    #  ( !( TRUE %in% same_identity_found_once) ) &
    #  ( !( TRUE %in% same_identity) )
    #){
    #  
    #  min_entry = which.min( b_table$Q_values )
    #  same_identity_found_once[min_entry] = FALSE # think about it
    #}
    if ( sum( identical_cls != "") == 0 )
      same_identity_found_once = rep("",dim(b_table)[1])
    
    related_identity = b_name_plane %in% related_cls | same_identity
    related_identity_found = related_identity == b_table$Identification_sig
    related_identity_found_once = rep("",dim(b_table)[1])
    related_identity_found_once[ match(TRUE, related_identity & b_table$Identification_sig, nomatch = 0) ] = TRUE

    if ( sum(related_cls != "") == 0 )
      related_identity_found_once = rep("",dim(b_table)[1])
    
    source_file = rep( source_file, dim( b_table )[1] )
    
    res_table <<- data.frame(
      "CL_name_query"            = c( as.character( res_table$CL_name_query ) , as.character( match_name_query ) ),
      "CL_name_training"         = c( as.character( res_table$CL_name_training ) , as.character( b_table$CCL ) ),
      "Source_query"             = c( as.character( res_table$Source_query ) , as.character( source_file ) ),
      "Source_training"          = c( as.character( res_table$Source_training), as.character( b_table$Library ) ),
      "Found_muts_abs"           = c( as.character( res_table$Found_muts_abs), as.character( b_table$Matches ) ),
      "P_value"                  = c( as.character( res_table$P_value), as.character( b_table$P_values ) ),
      "Passed_threshold"         = c( as.character( res_table$Passed_threshold), as.character( b_table$Identification_sig ) ),
      "Same_identity"            = c( as.character( res_table$Same_identity ) , as.character( same_identity ) ),
      "Same_identity_found"      = c( as.character( res_table$Same_identity_found ) , as.character( same_identity_found ) ),
      "Related_identity"         = c( as.character( res_table$Related_identity ), as.character( related_identity ) ),
      "Related_identity_found"   = c( as.character( res_table$Related_identity_found ), as.character( related_identity_found ) )
    )
    res_table = res_table[ as.character( res_table$Found_muts_abs ) != "0",]
    
    # identification table
    
    if (exclude_self)
      b_table = b_table[ paste( b_table$CCL, b_table$Library, sep = "_" ) != b_name_with_library, ]
    
    to_be_found_cls = unique( c( identical_cls, related_cls ) ) # in the future: alternate here
    to_be_found_cls = to_be_found_cls[ to_be_found_cls != "" ]
    
    found = unique( b_name_ccls_plane[ which( b_table$Identification_sig ) ] )

    true_pos_ident  = to_be_found_cls[ which( to_be_found_cls %in% found )  ]
    false_neg_ident = to_be_found_cls[ which( !(to_be_found_cls %in%  found) )  ]
    
    false_pos_ident = found[ which( !(found %in% to_be_found_cls) )  ]
    #false_pos_ident = paste( false_pos_ident, , sep ="_" )
    
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
    
    found = b_table$CCL[ which( b_table$Identification_sig )  ]
    found = paste( found, b_table$Library[ which( b_table$Identification_sig) ], sep ="_" )
    
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
    write.table( res_table, benchmark_res_file_path, sep ="\t", quote =F, row.names = F )
  
}

run_benchmark(
    inclusion_weight,
    only_first,
    exclude_self,
    run_identification = run_identification,
    #confidence_score = confidence_score,
    minimum_matching_mutations = minimum_matching_mutations,
    type_benchmark = type_benchmark,
    ref_gen = ref_gen
)
