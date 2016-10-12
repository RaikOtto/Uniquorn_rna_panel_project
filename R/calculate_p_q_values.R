#' calculate_p_and_q_values
#' 
#' @param candidate_hits_abs_all Maximally possible found variants
#' @param cl_absolute_mutation_hits Matching variants
#' @param sim_list Contains reference mutation data
#' @param sim_list_stats Contains global reference mutation stats
#' @param minimum_matching_mutations Minimal amount of required 
#' matching mutations
#' @param list_of_cls List of CLs
#' @param p_value Required maximal p-value
#' @param q_value Required maximal q-value
#' @param vcf_fingerprint The start and end positions of variants in the query
#' @param panels The reference libraries
#' @import stats
#' @return Results table
calculate_p_and_q_values = function(
    candidate_hits_abs_all,
    cl_absolute_mutation_hits,
    sim_list,
    sim_list_stats,
    minimum_matching_mutations,
    list_of_cls,
    p_value,
    q_value,
    vcf_fingerprint,
    panels
){
    
    p_values <<- rep( 1.0, length(list_of_cls))
    
    for ( panel in panels ){
    
        index_panel_sim      = grep( sim_list$CL, pattern = panel )
        index_panel_sim_stats= grep( sim_list_stats$CL, pattern = panel )
        index_panel_list     = grep( list_of_cls, pattern = panel )
        
        #white_balls_possible = sim_list_stats$Count[ index_panel_sim_stats ]
        white_balls_possible   = subset( 
            x = sim_list_stats$Count, 
            grepl( panel, sim_list_stats$CL)
        )
        
        if ( ( panel == "_CUSTOM") && 
             ( length( p_values_panel ) == 1 ) && 
             ( white_balls_possible == 0 ) 
        )
            
        message("There is only one CL in the customt set. 
              A confidence scoere calculation is only 
              possible with more than one sample!")
            
        white_balls_found    = candidate_hits_abs_all[ index_panel_list ]
        white_balls_found_least_one = sum( white_balls_found > 0 )
        black_balls          = sum( white_balls_possible ) - white_balls_possible
        
        #background_cls_traces = sum( white_balls_found_least_one >= mean( white_balls_found_least_one ) )
        #sig_p_values = sum( p_values >= mean( white_balls_found_least_one ) )
        background_cls_traces = sum( white_balls_found >= mean( white_balls_found ) )
        
        x = seq(0,1, length = 100)
        penalty = max( ( stats::pbeta( x, 1, background_cls_traces ) - x ) )
        
        #likelihood = 1 / white_balls_possible
        likelihood = white_balls_possible / sum( white_balls_possible )
        
        q                    = white_balls_found - 1
        q[ q < 0 ]           = 0

        p_values_panel       = as.double( stats::pbinom( 
            q = q,
            size = sum( white_balls_found ),
            p = likelihood,
            lower.tail = FALSE 
        ) )
        
        p_values_panel[ white_balls_found == white_balls_possible ] = 0.0
        p_values_panel[ p_values_panel >= 1.0 ] = 1.0
        p_values_panel[ white_balls_found == 0 ] = 1.0
        p_values_panel[ white_balls_found < minimum_matching_mutations ] = 1.0
        
        p_values[ index_panel_list ] = p_values_panel
        
    }
    
    return( p_values )       
}