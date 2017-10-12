add_p_q_values_statistics = function( 
    g_query,
    match_t,
    q_value,
    ref_gen,
    minimum_matching_mutations
){
  
    library_names = read_library_names(ref_gen = ref_gen)
    p_values = rep(1.0, nrow(match_t))
    balls_in_query = length(g_query$Member_CCLs)
  
    for (library_name in library_names){
        
        index_library = which( match_t$Library == library_name  )
        white_balls_possible = as.integer(as.character(match_t$All_variants))[index_library]
        
        if( length(white_balls_possible) == 0){
            stop("There is only one CL in the custom set. 
                  A confidence scoere calculation is only 
                  possible with more than one sample!")
        }
        
        white_balls_found = as.integer(as.character(match_t$Matches))[index_library]
        white_balls_found_least_one = sum(white_balls_found > 0)
        black_balls = sum(white_balls_possible) - white_balls_possible
        
        background_cls_traces = sum(white_balls_found >= mean(white_balls_found))
        
        #x = seq(0,1, length = 100)
        #likelihood = ( 1 / balls_in_query ) ** (sum(white_balls_found) / balls_in_query )
        #likelihood = ( balls_in_query / white_balls_possible ) ** (white_balls_found / sum(white_balls_found))
        likelihood = ( sum(white_balls_found) / (sum(white_balls_found) + white_balls_found ) ) ** 
            (white_balls_found / sum(white_balls_found))
        
        q = white_balls_found - 1
        q[q < 0]   = 0
        
        p_values_panel = as.double(stats::pbinom( 
            q = q,
            size = sum(white_balls_found),
            p = likelihood,
            lower.tail = FALSE 
        ))
        
        p_values_panel[white_balls_found == white_balls_possible] = 0.0
        p_values_panel[p_values_panel >= 1.0] = 1.0
        p_values_panel[white_balls_found == 0] = 1.0
        p_values_panel[white_balls_found < minimum_matching_mutations] = 1.0
        
        p_values[index_library] = p_values_panel
      
    }
    match_t$P_values = p_values
    match_t$Q_values = p.adjust(p_values, method = "BH")
    match_t$Q_value_sig = match_t$Q_values <= q_value
  
    return(match_t)
}

add_penality_statistics = function( match_t, minimum_matching_mutations  ){

    matching_variants = as.integer(match_t$Matches)
    minimum_matching_mutations = min(matching_variants)
    matching_variants = matching_variants[ matching_variants > 0]
  
    if ( length(matching_variants) == 0 ){
      
        penalty = 0
        penalty_mutations = 0
      
    } else{ 
      
        mean_match = mean( matching_variants )
        max_match  = max( matching_variants )
        
        penalty = integrate(
            f = pbeta,
            0,
            1,
            max_match,
            max_match/ mean_match,
            stop.on.error = FALSE
            )$value
          
        penalty_mutations = 
            ceiling(
            mean_match + 
              ( max_match * penalty ) / 
              ( 1 - penalty )
        )
    }
    
    
    if ( ( minimum_matching_mutations == 0 ) & ( penalty != 0.0) ){
      
        message( paste0( collapse = "", c( 
            "Correcting the background due to traces of random, scale-freeness amounts of matches, 
            requiring at least ", 
            as.character( penalty_mutations ), " variants to match." ) )
        )
    }
    match_t$Above_Penality = as.integer(as.character(match_t$Matches)) > 
        penalty
    
    return(match_t)
}