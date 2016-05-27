#' calculate_p_and_q_values
#' 
#' @param candidate_hits_abs_all Maximally possible found variants
#' @param cl_absolute_mutation_hits Matching variants
#' @param sim_list Contains reference mutation data
#' @param minimum_matching_mutations Minimal amount of required matching mutations
#' @param list_of_cls List of CLs
#' @param p_value Required maximal p-value
#' @param q_value Required maximal q-value
#' @return Results table
calculate_p_and_q_values = function(
    candidate_hits_abs_all,
    cl_absolute_mutation_hits,
    sim_list,
    minimum_matching_mutations,
    list_of_cls,
    p_value,
    q_value
){
    
    p_values <<- rep( 1.0, length(list_of_cls))

    panels = unique( sapply( list_of_cls, FUN = function( cl_name ){
        return( 
            paste0( "_",
                utils::tail( unlist( stringr::str_split( cl_name, pattern = "_") ), 1 ) 
            )
        )
    } ) )
    
    for ( panel in panels ){
    
        index_panel_sim      = grep( sim_list$CL, pattern = panel )
        index_panel_list     = grep( list_of_cls, pattern = panel )
        
        var_per_panel        = length( sim_list$Fingerprint[ index_panel_sim ] )
        muts_per_cl          = cl_absolute_mutation_hits[ index_panel_list ]
        white_balls_possible = var_per_panel - muts_per_cl
        white_balls_found    = candidate_hits_abs_all[ index_panel_list ]
        black_balls          = var_per_panel - white_balls_possible
        likelihood           = cl_absolute_mutation_hits[ index_panel_list ] / var_per_panel

        q                    = white_balls_found - 1
        q[ q < 0 ]           = 0

        p_values_panel       = as.double( stats::pbinom( 
            q = q, 
            size = sum( white_balls_found ), 
            p = likelihood, 
            lower.tail = FALSE 
        ) )
        p_values_panel[ white_balls_found < minimum_matching_mutations ] = 1.0

        p_values[ index_panel_list ] = p_values_panel
        
    }
    
    return( p_values )       
}