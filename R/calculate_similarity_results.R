#' calculate_similarity_results
#' 
#' @param sim_list Contains reference mutation data
#' @param sim_list_stats Contains global reference mutation stats
#' @param found_mut_mapping Mapping to mutations from query to reference mutation set
#' @param minimum_matching_mutations Minimal amount of required matching mutations
#' @param p_value Required maximal p-value
#' @param q_value Required maximal q-value
#' @return Results table
calculate_similarity_results = function(
    sim_list,
    sim_list_stats,
    found_mut_mapping,
    minimum_matching_mutations,
    p_value,
    q_value
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
        
        ### unweighted scores
        
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
        
        ### weighted scores
        
        # aggregate over weights & CL identifier
        
        aggregation = stats::aggregate(
            x  = as.double( sim_list$Weight[ found_mut_mapping  ] ),
            by = list( as.character( sim_list$CL[ found_mut_mapping ] )  ),
            FUN = sum
        )
        
        weight_all = sim_list_stats$All_weights[ 
            match( aggregation$Group.1, sim_list_stats$CL )
            ]
        mapping_to_cls = match( 
            list_of_cls,
            aggregation$Group.1,
            nomatch = 0
        )
        
        names( res_cl_weighted ) = list_of_cls
        res_cl_weighted[ names( res_cl_weighted ) %in% aggregation$Group.1  ] = 
            aggregation$x[ mapping_to_cls ]
        stats_all_weight = sim_list_stats$All_weights[ 
            match( list_of_cls, sim_list_stats$CL  )
            ]
        
        res_res_cl_weighted = round( as.double(res_cl_weighted  ) / 
                                         stats_all_weight * 100, 1 )
        res_cl_weighted = round(res_cl_weighted, 0)
        
        cl_absolute_mutation_hits = sim_list_stats$Count[ cl_match_stats ]
    }
    
    # p and q-value calculation 
    
    q = as.integer( candidate_hits_abs_all ) - 1
    q[ q < 0 ] = 0
    m = as.integer( cl_absolute_mutation_hits )
    #k = sum( as.integer( candidate_hits_abs_all ) )
    #n = k - q
    p = m / sum(as.integer( cl_absolute_mutation_hits ))
    
    #p_values = phyper( q = q, m = m, n = n, k = k, lower.tail = F)
    p_values = as.double( stats::pbinom( q = q, size = sum(q), p = p, lower.tail = FALSE) )
    q_values = as.double( stats::p.adjust( p_values, "BH" ) )

    # treshold
    
    # threshold non weighted
    
    passed_threshold_vec_p_value = rep(FALSE,length(p_values))
    passed_threshold_vec_p_value[ p_values <= p_value] = TRUE
    passed_threshold_vec_q_value = rep(FALSE,length(q_values))
    passed_threshold_vec_q_value[ q_values <= q_value] = TRUE
    
    #passed_threshold_weighted = rep( "", nr_cls )
    
    #passed_threshold_weighted = rep( FALSE, nr_cls )
    #passed_threshold_weighted[ 
    #    ( candidate_hits_abs_all >= minimum_matching_mutations ) & 
    #        ( candidate_hits_rel     >= similarity_threshold ) & 
    #        ( res_res_cl_weighted    >= similarity_threshold ) 
    #    ] = TRUE
    
    output_cl_names = stringr::str_replace( list_of_cls, pattern = "_CCLE|_COSMIC|_CELLMINER|_CUSTOM", replacement = "" )
    panel_vec = rep("", length( output_cl_names ))
    panel_vec[ stringr::str_detect( list_of_cls, "_CCLE" ) ] = "CCLE"
    panel_vec[ stringr::str_detect( list_of_cls, "_COSMIC" ) ] = "COSMIC"
    panel_vec[ stringr::str_detect( list_of_cls, "_CELLMINER" ) ] = "CELLMINER"
    panel_vec[ stringr::str_detect( list_of_cls, "_CUSTOM" ) ] = "CUSTOM"
    
    res_table = data.frame(
        "CL"                        = output_cl_names,
        "CL_source"                 = panel_vec,
        "Found_muts"                = as.character( candidate_hits_abs_all ),
        "Count_mutations"           = as.character(  cl_absolute_mutation_hits ),
        "P_values"                  = as.character( round( p_values, 3 ) ),
        "Q_values"                  = as.character( round( q_values, 3 ) ),
        "P_value_sig"               = as.character( passed_threshold_vec_p_value ),
        "Q_value_sig"               = as.character( passed_threshold_vec_q_value )
    )
    
    #"Found_muts_rel"           = as.character(  candidate_hits_rel ),
    #"Found_muts_weighted"      = as.character( res_cl_weighted ),
    #"Count_mutations_weighted" = as.character( round( stats_all_weight, 0 ) ),
    #"Found_muts_weighted_rel"  = as.character( res_res_cl_weighted )
    
    res_table = res_table[ order( as.double( as.character( res_table$P_values) ), decreasing = FALSE),  ]

}