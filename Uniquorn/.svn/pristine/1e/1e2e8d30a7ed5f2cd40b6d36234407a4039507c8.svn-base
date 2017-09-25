#' calculate_similarity_results
#' 
#' @param sim_list Contains reference mutation data
#' @param sim_list_stats Contains global reference mutation stats
#' @param found_mut_mapping Mapping to mutations from query to reference mutation set
#' @param minimum_matching_mutations Minimal amount of required matching mutations
#' @param p_value Required maximal p-value
#' @param q_value Required maximal q-value
#' @param confidence_score Threshold above which a positive prediction occurs
#' default 3.0
#' @return Results table
calculate_similarity_results = function(
    sim_list,
    sim_list_stats,
    found_mut_mapping,
    minimum_matching_mutations,
    p_value,
    q_value,
    confidence_score
){
    
    list_of_cls       = unique( sim_list$CL )
    nr_cls            = length( list_of_cls  ) # amount cls
    
    candidate_hits_abs_all = rep(0, nr_cls)
    names(candidate_hits_abs_all) = list_of_cls
    
    candidate_hits_abs      = rep( 0.0, nr_cls )
    candidate_hits_abs_all  = rep( 0.0, nr_cls )
    candidate_hits_rel      = rep( 0.0, nr_cls )
    
    cl_weight               = rep( 0.0, nr_cls )
    cl_weight_rel           = rep( 0.0, nr_cls )
    all_weighted            = rep( 0.0, nr_cls )
    
    res_cl_weighted         = rep( 0.0, nr_cls )
    res_res_cl_weighted     = rep( 0.0, nr_cls )
    stats_all_weight        = rep( 0.0, nr_cls )
    
    cl_absolute_mutation_hits = rep( 0.0, nr_cls )
    
    if ( length( found_mut_mapping ) != 0){
        
        ### scores
        
        candidate_hits_abs = stats::aggregate( 
            rep(1, length(found_mut_mapping)),
            by = list(sim_list$CL[ found_mut_mapping ]), 
            FUN = sum
        )
        
        candidate_hits_abs_all[ which( list_of_cls %in% candidate_hits_abs$Group.1) ] = 
            as.integer( candidate_hits_abs$x[ match( 
                list_of_cls, 
                candidate_hits_abs$Group.1, 
                nomatch = 0 ) ] 
            )
        
        cl_match_stats    = match( list_of_cls, sim_list_stats$CL, nomatch = 0 ) # mapping
        candidate_hits_rel = round( 
            candidate_hits_abs_all / sim_list_stats$Count[ cl_match_stats ] * 100,
            1
        ) 
        
        cl_absolute_mutation_hits = sim_list_stats$Count[ cl_match_stats ]
    }
    
    # p and q-value calculation 
    
    p_values = calculate_p_and_q_values(
        candidate_hits_abs_all,
        cl_absolute_mutation_hits,
        sim_list,
        minimum_matching_mutations,
        list_of_cls,
        p_value,
        q_value
    )
    q_values       = stats::p.adjust( p_values, "BH")
    conf_score_vec = -log( q_values ) 
    conf_score_vec[ is.infinite( conf_score_vec ) | as.double( conf_score_vec ) >= 100 ] = 100

    # treshold
    
    # threshold non weighted
    
    passed_threshold_vec_p_value = rep(FALSE,length(p_values))
    passed_threshold_vec_p_value[ p_values <= p_value] = TRUE
    passed_threshold_vec_q_value = rep(FALSE,length(q_values))
    passed_threshold_vec_q_value[ q_values <= q_value] = TRUE
    
    output_cl_names = stringr::str_replace( list_of_cls, pattern = "_CCLE|_COSMIC|_CELLMINER|_CUSTOM", replacement = "" )
    panel_vec = rep("", length( output_cl_names ))
    panel_vec[ stringr::str_detect( list_of_cls, "_CCLE" ) ] = "CCLE"
    panel_vec[ stringr::str_detect( list_of_cls, "_COSMIC" ) ] = "COSMIC"
    panel_vec[ stringr::str_detect( list_of_cls, "_CELLMINER" ) ] = "CELLMINER"
    panel_vec[ stringr::str_detect( list_of_cls, "_CUSTOM" ) ] = "CUSTOM"
    
    passed_threshold_vec_p_value[ 
        as.integer( candidate_hits_abs_all     )  < 
        as.integer( minimum_matching_mutations ) 
    ] = FALSE
    
    passed_threshold_vec_q_value[ 
        as.integer( candidate_hits_abs_all     )  < 
        as.integer( minimum_matching_mutations ) 
    ] = FALSE
    
    passed_threshold_vec_con_score =
        as.double( conf_score_vec   ) >=
        as.double( confidence_score )
    
    passed_threshold_vec_con_score[
        as.integer( candidate_hits_abs_all     )  < 
        as.integer( minimum_matching_mutations ) 
    ] = FALSE
    
    res_table = data.frame(
        "CL"                        = output_cl_names,
        "CL_source"                 = panel_vec,
        "Found_muts"                = as.character( candidate_hits_abs_all ),
        "Count_mutations"           = as.character( cl_absolute_mutation_hits ),
        "P_values"                  = as.character( p_values ),
        "Q_values"                  = as.character( q_values ),
        "Conf_score"                = as.character( conf_score_vec ),
        "P_value_sig"               = as.character( passed_threshold_vec_p_value ),
        "Q_value_sig"               = as.character( passed_threshold_vec_q_value ),
        "Conf_score_sig"            = as.character( passed_threshold_vec_con_score )
    )
    
    res_table = res_table[ order( as.integer( as.character( res_table$Found_muts) ), decreasing = TRUE),  ]
    #res_table = res_table[ order( as.double( as.character( res_table$Conf_score) ), decreasing = TRUE),  ]
    
    return(res_table)
}