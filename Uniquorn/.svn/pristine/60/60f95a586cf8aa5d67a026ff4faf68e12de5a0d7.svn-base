#' filter_for_weights
#' 
#' Filter the reference set
#' 
#' \code{filter_for_weights} parses vcf file and output 
#' basic information
#' @param ref_gen Reference genome version. All training sets are 
#' associated with a reference genome version. Default: GRCH37
#' @param verbose Print additional information
#' @param sim_list Contains the mutations
#' @param sim_list_stats Contains the overal mutation statistics
#' @param mutational_weight_inclusion_threshold Lower bound for mutational weight to be 
#' included
#' @usage 
#' filter_for_weights( 
#' mutational_weight_inclusion_threshold,
#' ref_gen,
#' verbose,
#' sim_list,
#' sim_list_stats)
#' @return Filtered reference sets
filter_for_weights = function( 
    mutational_weight_inclusion_threshold,
    ref_gen,
    verbose,
    sim_list,
    sim_list_stats
    ){
    
    # filter for weights
    if ( mutational_weight_inclusion_threshold != 0.0  ){
        
        if (verbose)        
            print( base::paste0( c("Adjusted mutational inclusion weight, 
                                   only using mutations that are have a weight higher than ", 
                                   as.character(mutational_weight_inclusion_threshold)), collapse="") )
        
        sim_list = sim_list[ base::as.double(sim_list$Weight) >= 
                                 as.double( mutational_weight_inclusion_threshold) ,  ]
        sum_vec = rep(1, dim(sim_list)[1])
        
        sim_list_stats = stats::aggregate( as.double( sum_vec ), 
                                           by = list( sim_list$CL), FUN = sum  )
        
        all_weights = stats::aggregate( sim_list$Weight, by = list(sim_list$CL), FUN = sum)
        mapping = match( sim_list_stats$Group.1, all_weights$Group.1 )
        sim_list_stats = cbind( sim_list_stats, all_weights$x[mapping] )
        
        colnames( sim_list_stats ) = c("CL","Count", "All_weights")
        sim_list_stats$All_weights = round(sim_list_stats$All_weights,1)
        
        if (verbose)
            print( paste0( c("Found ", as.character( dim(sim_list)[1] ), 
                             " many mutations with mutational weight of at least ", 
                             mutational_weight_inclusion_threshold), collapse="")  )
    }
    res_list = list( "sim_list" = sim_list, "sim_list_stats" = sim_list_stats)
    
    return( res_list )
}